#include <cmath>
#include <cstddef>
#include <vector>
#include "eigencore_common.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"

// Native complex Hermitian Sylvester inertia (dense).
//
// A complex Hermitian M = A - sigma B is factored by the Hermitian
// Bunch-Kaufman LDL^H (LAPACK zhetrf, lower triangle): M = P L D L^H P^T with
// D block diagonal, its 1x1 blocks real and its 2x2 blocks Hermitian
// [a b; conj(b) c] (a, c real). M is congruent to D, so the eigenvalue counts
// of A - sigma B below / at / above zero are the sign counts of D's
// eigenvalues: a 1x1 block contributes its sign, a 2x2 block the signs of its
// two eigenvalues (mean +/- hypot((a - c) / 2, |b|); the small one from the
// determinant a c - |b|^2 to avoid cancellation). This replaces dsytrf on the
// real 2n x 2n embedding [Re -Im; Im Re] (every count doubled), which costs
// (2n)^3 / 3 = 8 n^3 / 3 real flops against 4 n^3 / 3 for zhetrf (n^3 / 3
// complex multiply-adds), and four times the memory.
//
// The returned tally has the fields and meaning of the real dense helper
// (src/small_dense.cpp, eigencore_dense_symmetric_inertia): neg, zero, pos,
// min_abs_pivot, max_abs_pivot, max_abs_multiplier (largest |L_ij| outside the
// 2x2 blocks), growth (max_abs_pivot * max(1, max_abs_multiplier)^2, a proxy
// of || |L| |D| |L^H| ||), norm1 (||M||_1), info (zhetrf info) and two_by_two.
// The R layer (R/inertia.R, R/complex_inertia.R) applies the same
// backward-error reliability test as for the real count.

namespace {

struct ComplexInertiaTally {
  double neg = 0.0;
  double zero = 0.0;
  double pos = 0.0;
  double min_abs = R_PosInf;
  double max_abs = 0.0;
  double max_mult = 0.0;
  double growth = 0.0;
  double norm1 = 0.0;
  double info = 0.0;
  double two_by_two = 0.0;

  void add_pivot(double value) {
    if (value < 0.0) {
      neg += 1.0;
    } else if (value > 0.0) {
      pos += 1.0;
    } else {
      zero += 1.0;
    }
    const double a = std::fabs(value);
    if (a < min_abs) min_abs = a;
    if (a > max_abs) max_abs = a;
  }
};

SEXP complex_inertia_tally_sexp(const ComplexInertiaTally& t) {
  const char* names[] = {"neg", "zero", "pos", "min_abs_pivot", "max_abs_pivot",
                         "max_abs_multiplier", "growth", "norm1", "info",
                         "two_by_two"};
  const double values[] = {t.neg, t.zero, t.pos,
                           R_FINITE(t.min_abs) ? t.min_abs : 0.0, t.max_abs,
                           t.max_mult, t.growth, t.norm1, t.info, t.two_by_two};
  const int len = 10;
  SEXP out_ = PROTECT(allocVector(REALSXP, len));
  SEXP names_ = PROTECT(allocVector(STRSXP, len));
  for (int i = 0; i < len; ++i) {
    REAL(out_)[i] = values[i];
    SET_STRING_ELT(names_, i, mkChar(names[i]));
  }
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(2);
  return out_;
}

inline double complex_modulus(const Rcomplex& z) {
  return std::hypot(z.r, z.i);
}

}  // namespace

// Dense complex Hermitian inertia of A - sigma B (B = NULL: identity; B may
// be complex Hermitian or real symmetric). Only the lower triangles of A and
// B are read; the imaginary parts of their diagonals are ignored.
extern "C" SEXP eigencore_dense_hermitian_inertia(SEXP A_, SEXP sigma_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_)) {
    error("A must be a complex matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue || LENGTH(dimA) != 2 || INTEGER(dimA)[0] != INTEGER(dimA)[1]) {
    error("A must be a square matrix");
  }
  const int n = INTEGER(dimA)[0];
  const double sigma = asReal(sigma_);
  if (!R_FINITE(sigma)) {
    error("sigma must be finite");
  }
  const bool has_B = B_ != R_NilValue;
  const bool complex_B = has_B && isComplex(B_);
  if (has_B) {
    SEXP dimB = getAttrib(B_, R_DimSymbol);
    if ((!isReal(B_) && !isComplex(B_)) || dimB == R_NilValue || LENGTH(dimB) != 2 ||
        INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
      error("B must be a double or complex matrix with the dimensions of A");
    }
  }
  ComplexInertiaTally tally;
  if (n == 0) {
    return complex_inertia_tally_sexp(tally);
  }
  const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<Rcomplex> M(nn);
  const Rcomplex* A = COMPLEX(A_);
  const Rcomplex* Bc = complex_B ? COMPLEX(B_) : nullptr;
  const double* Br = (has_B && !complex_B) ? REAL(B_) : nullptr;
  std::vector<double> colsum(static_cast<size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      const size_t ij = static_cast<size_t>(i) + static_cast<size_t>(j) * n;
      Rcomplex v = A[ij];
      if (Bc != nullptr) {
        v.r -= sigma * Bc[ij].r;
        v.i -= sigma * Bc[ij].i;
      } else if (Br != nullptr) {
        v.r -= sigma * Br[ij];
      } else if (i == j) {
        v.r -= sigma;
      }
      if (i == j) {
        v.i = 0.0;
      }
      if (!R_FINITE(v.r) || !R_FINITE(v.i)) {
        error("A - sigma B has non-finite entries");
      }
      M[ij] = v;
      const double m = complex_modulus(v);
      colsum[static_cast<size_t>(j)] += m;
      if (i != j) {
        colsum[static_cast<size_t>(i)] += m;
      }
    }
  }
  for (int j = 0; j < n; ++j) {
    if (colsum[static_cast<size_t>(j)] > tally.norm1) {
      tally.norm1 = colsum[static_cast<size_t>(j)];
    }
  }
  std::vector<int> ipiv(static_cast<size_t>(n), 0);
  char uplo = 'L';
  int info = 0;
  int lwork = -1;
  Rcomplex work_query;
  work_query.r = 0.0;
  work_query.i = 0.0;
  F77_CALL(zhetrf)(&uplo, &n, M.data(), &n, ipiv.data(), &work_query, &lwork,
                   &info FCONE);
  if (info < 0) {
    error("LAPACK zhetrf workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query.r);
  if (lwork < n) lwork = n;
  std::vector<Rcomplex> work(static_cast<size_t>(lwork));
  F77_CALL(zhetrf)(&uplo, &n, M.data(), &n, ipiv.data(), work.data(), &lwork,
                   &info FCONE);
  if (info < 0) {
    error("LAPACK zhetrf failed with info=%d", info);
  }
  tally.info = static_cast<double>(info);
  int k = 0;
  while (k < n) {
    const size_t kk = static_cast<size_t>(k) + static_cast<size_t>(k) * n;
    if (ipiv[static_cast<size_t>(k)] > 0) {
      tally.add_pivot(M[kk].r);
      for (int i = k + 1; i < n; ++i) {
        const double l = complex_modulus(M[static_cast<size_t>(i) + static_cast<size_t>(k) * n]);
        if (l > tally.max_mult) tally.max_mult = l;
      }
      k += 1;
    } else {
      if (k + 1 >= n) {
        error("LAPACK zhetrf returned an incomplete 2x2 pivot block");
      }
      // 2x2 Hermitian pivot [a conj(b); b c] (lower storage holds b).
      const double a = M[kk].r;
      const double b = complex_modulus(M[kk + 1]);
      const double c = M[kk + 1 + static_cast<size_t>(n)].r;
      const double mean = 0.5 * (a + c);
      const double rad = std::hypot(0.5 * (a - c), b);
      const double big = mean >= 0.0 ? mean + rad : mean - rad;
      const double det = a * c - b * b;
      const double small = big != 0.0 ? det / big : 0.0;
      tally.add_pivot(big);
      tally.add_pivot(small);
      tally.two_by_two += 1.0;
      for (int col = k; col <= k + 1; ++col) {
        for (int i = k + 2; i < n; ++i) {
          const double l = complex_modulus(M[static_cast<size_t>(i) + static_cast<size_t>(col) * n]);
          if (l > tally.max_mult) tally.max_mult = l;
        }
      }
      k += 2;
    }
  }
  const double m1 = tally.max_mult > 1.0 ? tally.max_mult : 1.0;
  tally.growth = tally.max_abs * m1 * m1;
  return complex_inertia_tally_sexp(tally);
  EIGENCORE_ENTRY_END
}
