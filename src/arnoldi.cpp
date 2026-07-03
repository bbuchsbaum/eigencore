#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_common.h"
#include "native_operators.h"

static Rcomplex complex_conj(Rcomplex z) {
  Rcomplex out;
  out.r = z.r;
  out.i = -z.i;
  return out;
}

static SEXP native_arnoldi_cycle_impl(void* impl,
                                      EigencoreApplyFn apply,
                                      int64_t n64,
                                      const double* start,
                                      int max_subspace) {
  if (n64 < 1 || max_subspace < 1) {
    error("native Arnoldi requires positive dimensions");
  }
  if (!eigencore_int_indexable(n64) || !eigencore_int_indexable(max_subspace + 1)) {
    error("native Arnoldi dimensions exceed LP64 BLAS/R integer range");
  }
  const int n = static_cast<int>(n64);
  const int m_budget = std::min(max_subspace, n);
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, m_budget + 1));
  SEXP H_ = PROTECT(allocMatrix(REALSXP, m_budget + 1, m_budget));
  double* V = REAL(V_);
  double* H = REAL(H_);
  std::memset(V, 0, sizeof(double) * static_cast<size_t>(n) *
              static_cast<size_t>(m_budget + 1));
  std::memset(H, 0, sizeof(double) * static_cast<size_t>(m_budget + 1) *
              static_cast<size_t>(m_budget));

  double start_norm = 0.0;
  for (int i = 0; i < n; ++i) {
    start_norm += start[i] * start[i];
  }
  start_norm = std::sqrt(start_norm);
  if (!std::isfinite(start_norm) || start_norm <= 100.0 * DBL_EPSILON) {
    error("native Arnoldi start vector is numerically zero");
  }
  for (int i = 0; i < n; ++i) {
    V[i] = start[i] / start_norm;
  }

  std::vector<double> w(static_cast<size_t>(n));
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int iterations = 0;
  int matvecs = 0;
  int reorthogonalization_passes = 0;

  for (int j = 0; j < m_budget; ++j) {
    std::fill(w.begin(), w.end(), 0.0);
    const int status = apply(
      impl, EIGENCORE_TRANSPOSE_NONE, 1,
      V + static_cast<int64_t>(j) * n, n,
      1.0, 0.0, w.data(), n, &workspace
    );
    if (status != 0) {
      eigencore_apply_status_error("native Arnoldi operator apply", status);
    }
    ++matvecs;

    for (int pass = 0; pass < 2; ++pass) {
      for (int i = 0; i <= j; ++i) {
        const double* vi = V + static_cast<int64_t>(i) * n;
        double hij = 0.0;
        for (int row = 0; row < n; ++row) {
          hij += vi[row] * w[row];
        }
        H[i + static_cast<int64_t>(j) * (m_budget + 1)] += hij;
        for (int row = 0; row < n; ++row) {
          w[row] -= hij * vi[row];
        }
      }
      ++reorthogonalization_passes;
    }

    double beta = 0.0;
    for (int row = 0; row < n; ++row) {
      beta += w[row] * w[row];
    }
    beta = std::sqrt(beta);
    H[(j + 1) + static_cast<int64_t>(j) * (m_budget + 1)] = beta;
    iterations = j + 1;
    if (!std::isfinite(beta) || beta <= 100.0 * DBL_EPSILON ||
        iterations == m_budget) {
      break;
    }
    double* vnext = V + static_cast<int64_t>(j + 1) * n;
    for (int row = 0; row < n; ++row) {
      vnext[row] = w[row] / beta;
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out_, 0, V_);
  SET_VECTOR_ELT(out_, 1, H_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(static_cast<int>(workspace.bytes_allocated)));
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names_, 0, mkChar("V"));
  SET_STRING_ELT(names_, 1, mkChar("H"));
  SET_STRING_ELT(names_, 2, mkChar("iterations"));
  SET_STRING_ELT(names_, 3, mkChar("matvecs"));
  SET_STRING_ELT(names_, 4, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 5, mkChar("native_workspace_bytes"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_arnoldi_refined_ritz(SEXP V_, SEXP H_,
                                                SEXP iterations_,
                                                SEXP values_) {
  if (!isReal(V_) || !isReal(H_) || !isComplex(values_)) {
    error("native Arnoldi refined extraction requires real V/H and complex values");
  }
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  SEXP dimH = getAttrib(H_, R_DimSymbol);
  if (dimV == R_NilValue || dimH == R_NilValue ||
      LENGTH(dimV) != 2 || LENGTH(dimH) != 2) {
    error("native Arnoldi refined extraction requires matrix inputs");
  }

  const int n = INTEGER(dimV)[0];
  const int v_cols = INTEGER(dimV)[1];
  const int h_rows = INTEGER(dimH)[0];
  const int h_cols = INTEGER(dimH)[1];
  const int m = asInteger(iterations_);
  const int k = LENGTH(values_);
  if (n < 1 || m < 1 || k < 1 ||
      v_cols < m || h_rows < m + 1 || h_cols < m) {
    error("invalid native Arnoldi refined extraction dimensions");
  }
  if (!eigencore_int_indexable(static_cast<int64_t>(n) * k) ||
      !eigencore_int_indexable(static_cast<int64_t>(m + 1) * m)) {
    error("native Arnoldi refined extraction dimensions exceed LP64 BLAS/R integer range");
  }

  SEXP vectors_ = PROTECT(allocMatrix(CPLXSXP, n, k));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  Rcomplex* vectors = COMPLEX(vectors_);
  double* residuals = REAL(residuals_);
  const double* V = REAL(V_);
  const double* H = REAL(H_);
  const Rcomplex* values = COMPLEX(values_);
  const int rows = m + 1;

  // All workspace sizes depend only on the loop-invariant m and rows, so
  // allocate once and run the zgesvd workspace query once, outside the loop.
  std::vector<Rcomplex> z(static_cast<size_t>(rows) * static_cast<size_t>(m));
  std::vector<double> s(static_cast<size_t>(m));
  std::vector<Rcomplex> vt(static_cast<size_t>(m) * static_cast<size_t>(m));
  Rcomplex u_dummy;
  int ldu = 1;
  int ldvt = m;
  int info = 0;
  int lwork = -1;
  Rcomplex work_query;
  const int rwork_len = std::max(1, 5 * m);
  std::vector<double> rwork(static_cast<size_t>(rwork_len));
  char jobu = 'N';
  char jobvt = 'A';
  F77_CALL(zgesvd)(&jobu, &jobvt, &rows, &m, z.data(), &rows, s.data(),
                   &u_dummy, &ldu, vt.data(), &ldvt,
                   &work_query, &lwork, rwork.data(), &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK zgesvd workspace query failed for native Arnoldi refined extraction with info=%d", info);
  }
  lwork = std::max(1, static_cast<int>(work_query.r));
  std::vector<Rcomplex> work(static_cast<size_t>(lwork));
  // Real and imaginary coefficient parts split out so the Ritz vector can be
  // formed with two real dgemv calls against the real basis V.
  std::vector<double> coeff_r(static_cast<size_t>(m));
  std::vector<double> coeff_i(static_cast<size_t>(m));
  std::vector<double> ritz_r(static_cast<size_t>(n));
  std::vector<double> ritz_i(static_cast<size_t>(n));
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const int inc_one = 1;

  for (int col = 0; col < k; ++col) {
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < rows; ++i) {
        Rcomplex entry;
        entry.r = H[i + static_cast<int64_t>(j) * h_rows];
        entry.i = 0.0;
        if (i == j) {
          entry.r -= values[col].r;
          entry.i -= values[col].i;
        }
        z[i + static_cast<int64_t>(j) * rows] = entry;
      }
    }

    F77_CALL(zgesvd)(&jobu, &jobvt, &rows, &m, z.data(), &rows, s.data(),
                     &u_dummy, &ldu, vt.data(), &ldvt,
                     work.data(), &lwork, rwork.data(), &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgesvd failed for native Arnoldi refined extraction with info=%d", info);
    }
    residuals[col] = s[static_cast<size_t>(m - 1)];

    for (int j = 0; j < m; ++j) {
      const Rcomplex c = complex_conj(vt[(m - 1) + static_cast<int64_t>(j) * m]);
      coeff_r[static_cast<size_t>(j)] = c.r;
      coeff_i[static_cast<size_t>(j)] = c.i;
    }

    F77_CALL(dgemv)(&notrans, &n, &m, &one, V, &n,
                    coeff_r.data(), &inc_one,
                    &zero, ritz_r.data(), &inc_one FCONE);
    F77_CALL(dgemv)(&notrans, &n, &m, &one, V, &n,
                    coeff_i.data(), &inc_one,
                    &zero, ritz_i.data(), &inc_one FCONE);

    double norm2 = 0.0;
    for (int row = 0; row < n; ++row) {
      const double zr = ritz_r[static_cast<size_t>(row)];
      const double zi = ritz_i[static_cast<size_t>(row)];
      Rcomplex out;
      out.r = zr;
      out.i = zi;
      vectors[row + static_cast<int64_t>(col) * n] = out;
      norm2 += zr * zr + zi * zi;
    }

    const double norm = std::sqrt(norm2);
    if (std::isfinite(norm) && norm > DBL_EPSILON) {
      for (int row = 0; row < n; ++row) {
        Rcomplex* out = vectors + row + static_cast<int64_t>(col) * n;
        out->r /= norm;
        out->i /= norm;
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, vectors_);
  SET_VECTOR_ELT(out_, 1, residuals_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("vectors"));
  SET_STRING_ELT(names_, 1, mkChar("refined_residual_estimates"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_arnoldi_dense_cycle(SEXP A_, SEXP start_,
                                               SEXP max_subspace_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue || LENGTH(dimA) != 2 ||
      INTEGER(dimA)[0] != INTEGER(dimA)[1]) {
    error("A must be a square double matrix");
  }
  const int n = INTEGER(dimA)[0];
  if (LENGTH(start_) != n) {
    error("start length must equal matrix dimension");
  }
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  return native_arnoldi_cycle_impl(
    &impl, eigencore_dense_apply, n, REAL(start_), asInteger(max_subspace_)
  );
}

extern "C" SEXP eigencore_arnoldi_csc_cycle(SEXP i_, SEXP p_, SEXP x_,
                                             SEXP dim_, SEXP start_,
                                             SEXP max_subspace_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      LENGTH(dim_) != 2 || !isReal(start_)) {
    error("invalid CSC Arnoldi inputs");
  }
  const int nrow = INTEGER(dim_)[0];
  const int ncol = INTEGER(dim_)[1];
  if (nrow != ncol) {
    error("A must be a square dgCMatrix");
  }
  if (LENGTH(start_) != nrow) {
    error("start length must equal matrix dimension");
  }
  CSCOperator impl = {nrow, ncol, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_arnoldi_cycle_impl(
    &impl, eigencore_csc_apply, nrow, REAL(start_), asInteger(max_subspace_)
  );
}

extern "C" SEXP eigencore_arnoldi_r_operator_cycle(SEXP dim_, SEXP apply_,
                                                    SEXP start_,
                                                    SEXP max_subspace_) {
  if (!isInteger(dim_) || LENGTH(dim_) != 2 || TYPEOF(apply_) != CLOSXP ||
      !isReal(start_)) {
    error("invalid matrix-free Arnoldi inputs");
  }
  const int nrow = INTEGER(dim_)[0];
  const int ncol = INTEGER(dim_)[1];
  if (nrow != ncol) {
    error("A must be a square matrix-free operator");
  }
  if (LENGTH(start_) != nrow) {
    error("start length must equal operator dimension");
  }
  RApplyOperator impl = {nrow, ncol, apply_, R_NilValue};
  return native_arnoldi_cycle_impl(
    &impl, eigencore_r_operator_apply, nrow, REAL(start_), asInteger(max_subspace_)
  );
}

extern "C" SEXP eigencore_arnoldi_ritz(SEXP V_, SEXP H_, SEXP iterations_) {
  if (!isReal(V_) || !isReal(H_)) {
    error("native Arnoldi Ritz extraction requires real V and H matrices");
  }
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  SEXP dimH = getAttrib(H_, R_DimSymbol);
  if (dimV == R_NilValue || dimH == R_NilValue ||
      LENGTH(dimV) != 2 || LENGTH(dimH) != 2) {
    error("native Arnoldi Ritz extraction requires matrix inputs");
  }
  const int n = INTEGER(dimV)[0];
  const int v_cols = INTEGER(dimV)[1];
  const int h_rows = INTEGER(dimH)[0];
  const int h_cols = INTEGER(dimH)[1];
  const int m = asInteger(iterations_);
  if (n < 1 || m < 1 || v_cols < m || h_rows < m || h_cols < m) {
    error("invalid native Arnoldi Ritz dimensions");
  }
  if (!eigencore_int_indexable(static_cast<int64_t>(n) * m)) {
    error("native Arnoldi Ritz dimensions exceed LP64 BLAS/R integer range");
  }

  std::vector<double> Hm(static_cast<size_t>(m) * static_cast<size_t>(m));
  const double* H = REAL(H_);
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < m; ++row) {
      Hm[row + static_cast<int64_t>(col) * m] =
        H[row + static_cast<int64_t>(col) * h_rows];
    }
  }

  std::vector<double> wr(m), wi(m), vr(static_cast<size_t>(m) * static_cast<size_t>(m));
  double vl_dummy = 0.0;
  int ldvl = 1;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  char jobvl = 'N';
  char jobvr = 'V';
  F77_CALL(dgeev)(&jobvl, &jobvr, &m, Hm.data(), &m, wr.data(), wi.data(),
                  &vl_dummy, &ldvl, vr.data(), &m,
                  &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgeev workspace query failed for native Arnoldi Ritz extraction with info=%d", info);
  }
  lwork = std::max(1, static_cast<int>(work_query));
  std::vector<double> work(lwork);
  F77_CALL(dgeev)(&jobvl, &jobvr, &m, Hm.data(), &m, wr.data(), wi.data(),
                  &vl_dummy, &ldvl, vr.data(), &m,
                  work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgeev failed for native Arnoldi Ritz extraction with info=%d", info);
  }

  SEXP values_ = PROTECT(allocVector(CPLXSXP, m));
  SEXP vectors_ = PROTECT(allocMatrix(CPLXSXP, n, m));
  Rcomplex* values = COMPLEX(values_);
  Rcomplex* vectors = COMPLEX(vectors_);
  const double* V = REAL(V_);
  for (int col = 0; col < m; ++col) {
    values[col].r = wr[col];
    values[col].i = wi[col];
  }
  for (int ritz_col = 0; ritz_col < m; ++ritz_col) {
    std::vector<double> coeff_re(m, 0.0), coeff_im(m, 0.0);
    if (wi[ritz_col] == 0.0) {
      for (int row = 0; row < m; ++row) {
        coeff_re[row] = vr[row + static_cast<int64_t>(ritz_col) * m];
      }
    } else if (wi[ritz_col] > 0.0 && ritz_col + 1 < m) {
      for (int row = 0; row < m; ++row) {
        coeff_re[row] = vr[row + static_cast<int64_t>(ritz_col) * m];
        coeff_im[row] = vr[row + static_cast<int64_t>(ritz_col + 1) * m];
      }
    } else if (wi[ritz_col] < 0.0 && ritz_col > 0) {
      for (int row = 0; row < m; ++row) {
        coeff_re[row] = vr[row + static_cast<int64_t>(ritz_col - 1) * m];
        coeff_im[row] = -vr[row + static_cast<int64_t>(ritz_col) * m];
      }
    }

    double norm2 = 0.0;
    for (int row = 0; row < n; ++row) {
      double zr = 0.0;
      double zi = 0.0;
      for (int basis = 0; basis < m; ++basis) {
        const double vb = V[row + static_cast<int64_t>(basis) * n];
        zr += vb * coeff_re[basis];
        zi += vb * coeff_im[basis];
      }
      Rcomplex z;
      z.r = zr;
      z.i = zi;
      vectors[row + static_cast<int64_t>(ritz_col) * n] = z;
      norm2 += zr * zr + zi * zi;
    }
    const double norm = std::sqrt(norm2);
    if (std::isfinite(norm) && norm > DBL_EPSILON) {
      for (int row = 0; row < n; ++row) {
        Rcomplex* z = vectors + row + static_cast<int64_t>(ritz_col) * n;
        z->r /= norm;
        z->i /= norm;
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}
