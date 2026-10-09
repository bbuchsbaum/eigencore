#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "eigencore_common.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"
#include "native_operators.h"

static bool ritz_value_better(double candidate, double incumbent, int target_kind) {
  switch (target_kind) {
    case 2:
      return candidate < incumbent;
    case 3:
      return fabs(candidate) > fabs(incumbent);
    case 4:
      return fabs(candidate) < fabs(incumbent);
    case 1:
    default:
      return candidate > incumbent;
  }
}

static int selected_ritz_indices(const double* values,
                                 int n,
                                 int k,
                                 int target_kind,
                                 int* selected) {
  const int count = (k < n) ? k : n;
  std::vector<bool> taken(static_cast<size_t>(n), false);
  for (int i = 0; i < count; ++i) {
    int best = -1;
    for (int j = 0; j < n; ++j) {
      if (taken[static_cast<size_t>(j)]) {
        continue;
      }
      if (best < 0 || ritz_value_better(values[j], values[best], target_kind)) {
        best = j;
      }
    }
    selected[i] = best;
    if (best >= 0) {
      taken[static_cast<size_t>(best)] = true;
    }
  }
  return count;
}

// ---------------------------------------------------------------------------
// Block B-orthonormalisation and projection helpers (BLAS-3, double precision).
//
// Conventions: blocks are column-major n x m with leading dimension n. For the
// standard problem (no metric B) every "B-image" pointer is nullptr and the
// Euclidean inner product is used; for the generalized problem B-images are
// carried alongside the blocks and updated by the same linear combinations, so
// B is applied to each new block exactly once.
// ---------------------------------------------------------------------------

// Columns whose Euclidean norm shrinks below this fraction of their norm
// before projection are treated as numerically inside the span of the
// projected-out bases and dropped. The test is relative, so it is invariant
// to scaling of A, B or the block.
static const double kLobpcgDropRelative = 1e-12;
// SVQB fallback keeps Gram eigenvalues above this fraction of the largest one
// (singular values above ~3e-7 of the largest column direction).
static const double kLobpcgSvqbRelative = 1e-13;
// Cholesky-QR is used only when the scaled Gram factor is this well
// conditioned (1-norm reciprocal condition estimate of R); CholQR2 is stable
// for cond(V) well below eps^{-1/2}.
static const double kLobpcgCholRcond = 1e-6;
// Explicitly recompute A X and B X (and re-B-orthonormalise / re-project X)
// at least this often; in between they are updated by the Rayleigh-Ritz
// coefficients.
static const int kLobpcgRefreshInterval = 16;

static inline double lobpcg_nrm2(const double* x, int n) {
  const int inc = 1;
  return F77_CALL(dnrm2)(&n, x, &inc);
}

static inline double lobpcg_dot(const double* x, const double* y, int n) {
  const int inc = 1;
  return F77_CALL(ddot)(&n, x, &inc, y, &inc);
}

static inline double* lobpcg_col(double* base, int n, int col) {
  return base + static_cast<int64_t>(col) * n;
}

static inline const double* lobpcg_col(const double* base, int n, int col) {
  return base + static_cast<int64_t>(col) * n;
}

static inline void lobpcg_copy_col(double* dst, const double* src, int n) {
  if (dst != src) {
    std::memmove(dst, src, sizeof(double) * static_cast<size_t>(n));
  }
}

struct LobpcgOperator {
  void* impl;
  EigencoreApplyFn apply;
  EigencoreWorkspace* workspace;

  int run(int cols, const double* X, int n, double* Y) const {
    if (cols <= 0) {
      return 0;
    }
    return apply(impl, EIGENCORE_TRANSPOSE_NONE, cols, X, n,
                 1.0, 0.0, Y, n, workspace);
  }
};

struct LobpcgScratch {
  std::vector<double> gram;
  std::vector<double> gram_work;
  std::vector<double> evals;
  std::vector<double> coeff;
  std::vector<double> block;
  std::vector<double> norms;
  std::vector<int> iwork;

  void reserve_gram(int m) {
    const size_t mm = static_cast<size_t>(m) * static_cast<size_t>(m);
    if (gram.size() < mm) gram.resize(mm);
    if (evals.size() < static_cast<size_t>(m)) evals.resize(static_cast<size_t>(m));
    if (norms.size() < static_cast<size_t>(m)) norms.resize(static_cast<size_t>(m));
    if (iwork.size() < static_cast<size_t>(m)) iwork.resize(static_cast<size_t>(m));
    const size_t lw = static_cast<size_t>(m > 0 ? 3 * m + 64 * m : 1);
    if (gram_work.size() < lw) gram_work.resize(lw);
  }
  void reserve_block(int n, int m) {
    const size_t sz = static_cast<size_t>(n) * static_cast<size_t>(m > 0 ? m : 1);
    if (block.size() < sz) block.resize(sz);
  }
  void reserve_coeff(int rows, int cols) {
    const size_t sz = static_cast<size_t>(rows > 0 ? rows : 1) *
                      static_cast<size_t>(cols > 0 ? cols : 1);
    if (coeff.size() < sz) coeff.resize(sz);
  }
};

// V <- V - U (BU' V); BV <- BV - BU (BU' V) when BV is supplied. BU is the
// B-image of the B-orthonormal basis U (U itself for the standard problem),
// so no operator application is needed.
static void lobpcg_project_out(const double* U, const double* BU, int u_cols,
                               int n, int m, double* V, double* BV,
                               LobpcgScratch& scratch) {
  if (u_cols <= 0 || m <= 0) {
    return;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  scratch.reserve_coeff(u_cols, m);
  double* coeff = scratch.coeff.data();
  F77_CALL(dgemm)(&trans_T, &trans_N, &u_cols, &m, &n,
                  &one, BU, &n, V, &n,
                  &zero, coeff, &u_cols FCONE FCONE);
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &m, &u_cols,
                  &minus_one, U, &n, coeff, &u_cols,
                  &one, V, &n FCONE FCONE);
  if (BV != nullptr) {
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &m, &u_cols,
                    &minus_one, BU, &n, coeff, &u_cols,
                    &one, BV, &n FCONE FCONE);
  }
}

// G = V' (B V) (m x m), symmetrised.
static void lobpcg_gram(const double* V, const double* BV, int n, int m,
                        double* G) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  if (BV == nullptr) {
    const char uplo = 'U';
    F77_CALL(dsyrk)(&uplo, &trans_T, &m, &n, &one, V, &n, &zero, G, &m
                    FCONE FCONE);
    for (int j = 0; j < m; ++j) {
      for (int i = j + 1; i < m; ++i) {
        G[i + static_cast<int64_t>(j) * m] = G[j + static_cast<int64_t>(i) * m];
      }
    }
    return;
  }
  F77_CALL(dgemm)(&trans_T, &trans_N, &m, &m, &n,
                  &one, V, &n, BV, &n,
                  &zero, G, &m FCONE FCONE);
  for (int j = 0; j < m; ++j) {
    for (int i = j + 1; i < m; ++i) {
      const double avg = 0.5 * (G[i + static_cast<int64_t>(j) * m] +
                                G[j + static_cast<int64_t>(i) * m]);
      G[i + static_cast<int64_t>(j) * m] = avg;
      G[j + static_cast<int64_t>(i) * m] = avg;
    }
  }
}

// One B-orthonormalisation pass of the m columns of V (and B V): Cholesky QR
// of the diagonally scaled Gram matrix when it is well conditioned, otherwise
// SVQB (eigendecomposition of the scaled Gram matrix) which drops numerically
// dependent directions. The surviving r columns are left at the front of V
// and BV. Returns r >= 0, or a negative status.
static int lobpcg_cholqr_pass(double* V, double* BV, int n, int m,
                              LobpcgScratch& scratch) {
  if (m <= 0) {
    return 0;
  }
  scratch.reserve_gram(m);
  double* G = scratch.gram.data();
  double* d = scratch.norms.data();
  lobpcg_gram(V, BV, n, m, G);
  bool positive_diag = true;
  for (int j = 0; j < m; ++j) {
    const double gjj = G[j + static_cast<int64_t>(j) * m];
    if (!R_FINITE(gjj)) {
      return -10;
    }
    if (gjj > 0.0) {
      d[j] = sqrt(gjj);
    } else {
      d[j] = 1.0;
      positive_diag = false;
    }
  }
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      G[i + static_cast<int64_t>(j) * m] /= d[i] * d[j];
    }
  }

  if (positive_diag) {
    std::vector<double> Rf(G, G + static_cast<size_t>(m) * m);
    const char uplo = 'U';
    int info = 0;
    F77_CALL(dpotrf)(&uplo, &m, Rf.data(), &m, &info FCONE);
    if (info == 0) {
      const char norm1 = '1';
      const char diag_n = 'N';
      double rcond = 0.0;
      scratch.reserve_gram(m);
      F77_CALL(dtrcon)(&norm1, &uplo, &diag_n, &m, Rf.data(), &m, &rcond,
                       scratch.gram_work.data(), scratch.iwork.data(), &info
                       FCONE FCONE FCONE);
      if (info == 0 && rcond >= kLobpcgCholRcond) {
        // V <- V D^{-1} R^{-1}, and the same for B V.
        for (int j = 0; j < m; ++j) {
          const double inv = 1.0 / d[j];
          double* v = lobpcg_col(V, n, j);
          for (int row = 0; row < n; ++row) v[row] *= inv;
          if (BV != nullptr) {
            double* bv = lobpcg_col(BV, n, j);
            for (int row = 0; row < n; ++row) bv[row] *= inv;
          }
        }
        const char side = 'R';
        const char trans_N = 'N';
        const double one = 1.0;
        F77_CALL(dtrsm)(&side, &uplo, &trans_N, &diag_n, &n, &m, &one,
                        Rf.data(), &m, V, &n FCONE FCONE FCONE FCONE);
        if (BV != nullptr) {
          F77_CALL(dtrsm)(&side, &uplo, &trans_N, &diag_n, &n, &m, &one,
                          Rf.data(), &m, BV, &n FCONE FCONE FCONE FCONE);
        }
        return m;
      }
    }
  }

  // SVQB fallback: G_scaled = U diag(lambda) U'.
  {
    const char jobz = 'V';
    const char uplo = 'U';
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    F77_CALL(dsyev)(&jobz, &uplo, &m, G, &m, scratch.evals.data(),
                    &work_query, &lwork, &info FCONE FCONE);
    lwork = info == 0 && work_query > 0.0 ? static_cast<int>(work_query) : 3 * m;
    if (scratch.gram_work.size() < static_cast<size_t>(lwork)) {
      scratch.gram_work.resize(static_cast<size_t>(lwork));
    }
    F77_CALL(dsyev)(&jobz, &uplo, &m, G, &m, scratch.evals.data(),
                    scratch.gram_work.data(), &lwork, &info FCONE FCONE);
    if (info != 0) {
      return -3;
    }
  }
  const double* lambda = scratch.evals.data();
  const double lambda_max = lambda[m - 1];
  if (!(lambda_max > 0.0) || !R_FINITE(lambda_max)) {
    return 0;
  }
  // Columns of T = D^{-1} U_keep diag(lambda_keep)^{-1/2}, largest first,
  // written over the leading columns of G (eigenvectors are processed in
  // descending order so no column is overwritten before it is read).
  int r = 0;
  for (int idx = m - 1; idx >= 0; --idx) {
    if (lambda[idx] > kLobpcgSvqbRelative * lambda_max) {
      ++r;
    }
  }
  std::vector<double> T(static_cast<size_t>(m) * (r > 0 ? r : 1), 0.0);
  int out = 0;
  for (int idx = m - 1; idx >= 0 && out < r; --idx) {
    if (!(lambda[idx] > kLobpcgSvqbRelative * lambda_max)) {
      continue;
    }
    const double s = 1.0 / sqrt(lambda[idx]);
    for (int i = 0; i < m; ++i) {
      T[i + static_cast<size_t>(out) * m] =
        G[i + static_cast<int64_t>(idx) * m] * s / d[i];
    }
    ++out;
  }
  if (r == 0) {
    return 0;
  }
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  scratch.reserve_block(n, r);
  double* tmp = scratch.block.data();
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &r, &m, &one, V, &n, T.data(), &m,
                  &zero, tmp, &n FCONE FCONE);
  std::memcpy(V, tmp, sizeof(double) * static_cast<size_t>(n) * r);
  if (BV != nullptr) {
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &r, &m, &one, BV, &n, T.data(), &m,
                    &zero, tmp, &n FCONE FCONE);
    std::memcpy(BV, tmp, sizeof(double) * static_cast<size_t>(n) * r);
  }
  return r;
}

// B-orthonormalise the m columns of V against the B-orthonormal constraint
// basis (Y, BY) and the current block (X, BX), then among themselves. The
// projection is two-pass block CGS in the B inner product; the block is
// re-projected after the first normalisation and normalised again (CholQR2),
// so columns that came out of an ill-conditioned normalisation cannot carry
// constraint components back in. B is applied to the block exactly once.
// For the standard problem pass metric == nullptr and BV/BY/BX == nullptr.
static int lobpcg_b_orthonormalize_block(const LobpcgOperator* metric,
                                         int n, int m,
                                         double* V, double* BV,
                                         const double* Y, const double* BY,
                                         int y_cols,
                                         const double* X, const double* BX,
                                         int x_cols,
                                         LobpcgScratch& scratch) {
  if (m <= 0) {
    return 0;
  }
  const double* BYp = BY != nullptr ? BY : Y;
  const double* BXp = BX != nullptr ? BX : X;
  scratch.reserve_gram(m);
  std::vector<double> before(static_cast<size_t>(m), 0.0);
  for (int j = 0; j < m; ++j) {
    before[static_cast<size_t>(j)] = lobpcg_nrm2(lobpcg_col(V, n, j), n);
    if (!R_FINITE(before[static_cast<size_t>(j)])) {
      return -10;
    }
  }
  for (int pass = 0; pass < 2; ++pass) {
    lobpcg_project_out(Y, BYp, y_cols, n, m, V, nullptr, scratch);
    lobpcg_project_out(X, BXp, x_cols, n, m, V, nullptr, scratch);
  }
  int kept = 0;
  for (int j = 0; j < m; ++j) {
    const double after = lobpcg_nrm2(lobpcg_col(V, n, j), n);
    if (before[static_cast<size_t>(j)] > 0.0 &&
        after > kLobpcgDropRelative * before[static_cast<size_t>(j)]) {
      if (kept != j) {
        lobpcg_copy_col(lobpcg_col(V, n, kept), lobpcg_col(V, n, j), n);
      }
      ++kept;
    }
  }
  if (kept == 0) {
    return 0;
  }
  if (metric != nullptr) {
    const int status = metric->run(kept, V, n, BV);
    if (status != 0) {
      return status;
    }
  }
  int r = lobpcg_cholqr_pass(V, metric != nullptr ? BV : nullptr, n, kept, scratch);
  if (r <= 0) {
    return r;
  }
  // Re-project after normalisation, then normalise again.
  lobpcg_project_out(Y, BYp, y_cols, n, r, V,
                     metric != nullptr ? BV : nullptr, scratch);
  lobpcg_project_out(X, BXp, x_cols, n, r, V,
                     metric != nullptr ? BV : nullptr, scratch);
  r = lobpcg_cholqr_pass(V, metric != nullptr ? BV : nullptr, n, r, scratch);
  return r;
}

// Shifted tridiagonal preconditioner: LU with partial pivoting (dgttrf),
// factored once per solve and applied with dgttrs.
struct LobpcgTridiagonalFactor {
  int n = 0;
  std::vector<double> dl;
  std::vector<double> d;
  std::vector<double> du;
  std::vector<double> du2;
  std::vector<int> ipiv;

  int factor(const double* lower, const double* diag, const double* upper,
             int n_) {
    n = n_;
    if (n < 1 || diag == nullptr) {
      return -5;
    }
    const size_t off = static_cast<size_t>(n > 1 ? n - 1 : 1);
    dl.assign(off, 0.0);
    du.assign(off, 0.0);
    du2.assign(static_cast<size_t>(n > 2 ? n - 2 : 1), 0.0);
    d.assign(diag, diag + n);
    ipiv.assign(static_cast<size_t>(n), 0);
    double scale = 0.0;
    for (int i = 0; i < n; ++i) {
      if (!R_FINITE(d[static_cast<size_t>(i)])) return -5;
      scale = fmax(scale, fabs(d[static_cast<size_t>(i)]));
    }
    for (int i = 0; i + 1 < n; ++i) {
      dl[static_cast<size_t>(i)] = lower != nullptr ? lower[i] : 0.0;
      du[static_cast<size_t>(i)] = upper != nullptr ? upper[i] : 0.0;
      if (!R_FINITE(dl[static_cast<size_t>(i)]) || !R_FINITE(du[static_cast<size_t>(i)])) {
        return -5;
      }
      scale = fmax(scale, fmax(fabs(dl[static_cast<size_t>(i)]),
                               fabs(du[static_cast<size_t>(i)])));
    }
    int info = 0;
    F77_CALL(dgttrf)(&n, dl.data(), d.data(), du.data(), du2.data(),
                     ipiv.data(), &info);
    if (info != 0) {
      return -5;
    }
    // Reject numerically singular factors relative to the matrix scale.
    for (int i = 0; i < n; ++i) {
      if (fabs(d[static_cast<size_t>(i)]) <= DBL_EPSILON * scale) {
        return -5;
      }
    }
    return 0;
  }

  int solve(double* B, int nrhs) {
    if (nrhs <= 0) {
      return 0;
    }
    const char trans = 'N';
    int info = 0;
    F77_CALL(dgttrs)(&trans, &n, &nrhs, dl.data(), d.data(), du.data(),
                     du2.data(), ipiv.data(), B, &n, &info FCONE);
    return info == 0 ? 0 : -5;
  }
};

static int extract_shifted_symmetric_tridiagonal_from_csc(
    const int* i,
    const int* p,
    const double* x,
    int n,
    double shift,
    double* lower,
    double* diag,
    double* upper) {
  if (n < 1) {
    return -1;
  }
  std::memset(diag, 0, sizeof(double) * static_cast<size_t>(n));
  if (n > 1) {
    std::memset(lower, 0, sizeof(double) * static_cast<size_t>(n - 1));
    std::memset(upper, 0, sizeof(double) * static_cast<size_t>(n - 1));
  }

  for (int col = 0; col < n; ++col) {
    for (int pos = p[col]; pos < p[col + 1]; ++pos) {
      const int row = i[pos];
      const double value = x[pos];
      const int distance = row > col ? row - col : col - row;
      if (distance > 1) {
        return -6;
      }
      if (row == col) {
        diag[col] += value;
      } else if (row == col + 1) {
        lower[col] += value;
      } else {
        upper[row] += value;
      }
    }
  }

  for (int row = 0; row < n; ++row) {
    diag[row] += shift;
  }
  for (int row = 0; row < n - 1; ++row) {
    const double scale = fmax(fmax(fabs(lower[row]), fabs(upper[row])), 1.0);
    if (fabs(lower[row] - upper[row]) > 1e-12 * scale) {
      return -7;
    }
  }
  return 0;
}

static int lobpcg_constraint_matrix(SEXP constraints_,
                                    int n,
                                    const double** constraints,
                                    int* constraint_cols) {
  *constraints = nullptr;
  *constraint_cols = 0;
  if (!isReal(constraints_)) {
    return -1;
  }
  SEXP dimC = getAttrib(constraints_, R_DimSymbol);
  if (dimC == R_NilValue || INTEGER(dimC)[0] != n) {
    return -1;
  }
  const int cols = INTEGER(dimC)[1];
  if (cols < 0) {
    return -1;
  }
  const double* values = REAL(constraints_);
  for (int64_t pos = 0; pos < static_cast<int64_t>(n) * cols; ++pos) {
    if (!R_FINITE(values[pos])) {
      return -1;
    }
  }
  *constraints = values;
  *constraint_cols = cols;
  return 0;
}

static SEXP lobpcg_pack_result(int n, int k, const double* X,
                               const double* values,
                               const double* residuals,
                               const int* converged,
                               const double* hist_max_residual,
                               const int* hist_nconv,
                               int iterations, int matvecs,
                               int preconditioner_calls,
                               int q_rank_final,
                               int constraints_rank) {
  SEXP values_ = PROTECT(allocVector(REALSXP, k));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));
  SEXP hist_res_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP hist_nconv_ = PROTECT(allocVector(INTSXP, iterations));
  std::memcpy(REAL(values_), values, sizeof(double) * static_cast<size_t>(k));
  std::memcpy(REAL(vectors_), X, sizeof(double) * static_cast<size_t>(n) * k);
  std::memcpy(REAL(residuals_), residuals, sizeof(double) * static_cast<size_t>(k));
  for (int i = 0; i < k; ++i) {
    LOGICAL(converged_)[i] = converged[i] ? TRUE : FALSE;
  }
  for (int i = 0; i < iterations; ++i) {
    REAL(hist_res_)[i] = hist_max_residual[i];
    INTEGER(hist_nconv_)[i] = hist_nconv[i];
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 11));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, residuals_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, hist_res_);
  SET_VECTOR_ELT(out_, 5, hist_nconv_);
  SET_VECTOR_ELT(out_, 6, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(preconditioner_calls));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(q_rank_final));
  SET_VECTOR_ELT(out_, 10, ScalarInteger(constraints_rank));
  SEXP names_ = PROTECT(allocVector(STRSXP, 11));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("residuals"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("history_max_relative_residual"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("iterations"));
  SET_STRING_ELT(names_, 7, mkChar("matvecs"));
  SET_STRING_ELT(names_, 8, mkChar("preconditioner_calls"));
  SET_STRING_ELT(names_, 9, mkChar("q_rank_final"));
  SET_STRING_ELT(names_, 10, mkChar("constraints_rank"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(8);
  return out_;
}

// Symmetric eigensolve of an m x m matrix H in place (eigenvectors overwrite
// H, ascending eigenvalues in w).
static int lobpcg_dsyev(double* H, int m, double* w, std::vector<double>& work) {
  const char jobz = 'V';
  const char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsyev)(&jobz, &uplo, &m, H, &m, w, &work_query, &lwork, &info
                  FCONE FCONE);
  lwork = info == 0 && work_query > 0.0 ? static_cast<int>(work_query) : 3 * m;
  if (lwork < 1) lwork = 1;
  if (work.size() < static_cast<size_t>(lwork)) {
    work.resize(static_cast<size_t>(lwork));
  }
  F77_CALL(dsyev)(&jobz, &uplo, &m, H, &m, w, work.data(), &lwork, &info
                  FCONE FCONE);
  return info == 0 ? 0 : -3;
}

// LOBPCG (Knyazev 2001) in the orthogonal-basis form of Hetmaniuk & Lehoucq
// (2006): the Rayleigh-Ritz basis is [X Z] with X the current B-orthonormal
// Ritz block and Z a B-orthonormal basis of the active preconditioned
// residuals W and search directions P, B-orthogonalised against X and the
// constraints. The new search direction is the Z-part of the selected Ritz
// vectors, P = Z * Y_z, with A P and B P formed by the same coefficients, so P
// is consistent with the sign and rotation of the new Ritz vectors (C19).
// Converged pairs are soft-locked: their residuals and directions leave the
// trial basis while the pairs stay in the Rayleigh-Ritz block.
static int native_lobpcg_run(void* impl,
                             EigencoreApplyFn apply,
                             void* b_impl,
                             EigencoreApplyFn b_apply,
                             int n,
                             int k,
                             int maxit,
                             int target_kind,
                             double tol,
                             const double* start,
                             int use_tridiagonal_preconditioner,
                             const double* lower,
                             const double* diag,
                             const double* upper,
                             const double* constraints,
                             int constraint_cols,
                             double* X_out,
                             double* values_out,
                             double* residuals_out,
                             int* converged_out,
                             double* hist_max_residual,
                             int* hist_nconv,
                             int* iterations_out,
                             int* matvecs_out,
                             int* preconditioner_calls_out,
                             int* q_rank_final_out,
                             int* constraints_rank_out) {
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const bool generalized = b_apply != nullptr;
  const LobpcgOperator A_op = {impl, apply, &workspace};
  const LobpcgOperator B_op = {b_impl, b_apply, &workspace};
  const LobpcgOperator* metric = generalized ? &B_op : nullptr;

  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;

  *iterations_out = 0;
  *matvecs_out = 0;
  *preconditioner_calls_out = 0;
  *q_rank_final_out = k;
  *constraints_rank_out = 0;

  const size_t nk = static_cast<size_t>(n) * static_cast<size_t>(k);
  const size_t n2k = 2 * nk;
  LobpcgScratch scratch;

  LobpcgTridiagonalFactor preconditioner;
  if (use_tridiagonal_preconditioner) {
    const int status = preconditioner.factor(lower, diag, upper, n);
    if (status != 0) {
      return status;
    }
  }

  // Constraint basis: B-orthonormalised once; B Y recomputed explicitly so
  // the B-inner-product projections use an accurate image.
  std::vector<double> Y;
  std::vector<double> BY;
  int constraint_rank = 0;
  if (constraint_cols > 0) {
    Y.assign(constraints, constraints + static_cast<size_t>(n) * constraint_cols);
    if (generalized) {
      BY.assign(static_cast<size_t>(n) * constraint_cols, 0.0);
    }
    constraint_rank = lobpcg_b_orthonormalize_block(
      metric, n, constraint_cols, Y.data(),
      generalized ? BY.data() : nullptr,
      nullptr, nullptr, 0, nullptr, nullptr, 0, scratch);
    if (constraint_rank < 0) {
      return constraint_rank;
    }
    if (constraint_rank + k > n) {
      return -9;
    }
    if (generalized && constraint_rank > 0) {
      const int status = B_op.run(constraint_rank, Y.data(), n, BY.data());
      if (status != 0) {
        return status;
      }
    }
    *constraints_rank_out = constraint_rank;
  }
  const double* Yp = constraint_rank > 0 ? Y.data() : nullptr;
  const double* BYp = constraint_rank > 0 ? (generalized ? BY.data() : Y.data()) : nullptr;

  std::vector<double> X(start, start + nk);
  std::vector<double> AX(nk, 0.0);
  std::vector<double> BXs(generalized ? nk : 0, 0.0);
  std::vector<double> P(nk, 0.0);
  std::vector<double> R(nk, 0.0);
  std::vector<double> V(n2k, 0.0);
  std::vector<double> AV(n2k, 0.0);
  std::vector<double> BVs(generalized ? n2k : 0, 0.0);
  std::vector<double> Xt(nk, 0.0);
  std::vector<double> H(static_cast<size_t>(9) * k * k, 0.0);
  std::vector<double> theta(static_cast<size_t>(3) * k, 0.0);
  std::vector<double> Ysel(static_cast<size_t>(3) * k * k, 0.0);
  std::vector<double> eig_work;
  std::vector<int> selected(static_cast<size_t>(3) * k, 0);
  std::vector<int> active(static_cast<size_t>(k), 0);
  std::vector<double> rel(static_cast<size_t>(k), 0.0);

  // Lower bounds on ||A||_2 and ||B||_2 from every block the operators are
  // applied to (max_j ||A v_j|| / ||v_j||). They never exceed the Frobenius
  // norms the certificate uses, so the native stopping rule is never looser
  // than the certificate, and they scale with A and B.
  double norm_A_lb = 0.0;
  double norm_B_lb = generalized ? 0.0 : 1.0;
  auto update_norm_bounds = [&](const double* Vb, const double* AVb,
                                const double* BVb, int cols) {
    for (int j = 0; j < cols; ++j) {
      const double vn = lobpcg_nrm2(lobpcg_col(Vb, n, j), n);
      if (!(vn > 0.0)) continue;
      const double an = lobpcg_nrm2(lobpcg_col(AVb, n, j), n) / vn;
      if (R_FINITE(an) && an > norm_A_lb) norm_A_lb = an;
      if (BVb != nullptr) {
        const double bn = lobpcg_nrm2(lobpcg_col(BVb, n, j), n) / vn;
        if (R_FINITE(bn) && bn > norm_B_lb) norm_B_lb = bn;
      }
    }
  };

  // In-place B = B * U for an n x k block and k x k matrix U.
  auto right_multiply = [&](double* Bk, const double* U) {
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &k, &one, Bk, &n,
                    const_cast<double*>(U), &k, &zero, Xt.data(), &n
                    FCONE FCONE);
    std::memcpy(Bk, Xt.data(), sizeof(double) * nk);
  };

  // Rayleigh-Ritz on X alone: restores Ritz ordering after an explicit
  // refresh (X, AX, BX rotated by the same k x k orthogonal matrix).
  auto rayleigh_ritz_x = [&]() -> int {
    F77_CALL(dgemm)(&trans_T, &trans_N, &k, &k, &n, &one, X.data(), &n,
                    AX.data(), &n, &zero, H.data(), &k FCONE FCONE);
    for (int j = 0; j < k; ++j) {
      for (int i = j + 1; i < k; ++i) {
        const double avg = 0.5 * (H[i + static_cast<size_t>(j) * k] +
                                  H[j + static_cast<size_t>(i) * k]);
        H[i + static_cast<size_t>(j) * k] = avg;
        H[j + static_cast<size_t>(i) * k] = avg;
      }
    }
    const int status = lobpcg_dsyev(H.data(), k, theta.data(), eig_work);
    if (status != 0) return status;
    selected_ritz_indices(theta.data(), k, k, target_kind, selected.data());
    for (int col = 0; col < k; ++col) {
      std::memcpy(&Ysel[static_cast<size_t>(col) * k],
                  &H[static_cast<size_t>(selected[static_cast<size_t>(col)]) * k],
                  sizeof(double) * static_cast<size_t>(k));
    }
    right_multiply(X.data(), Ysel.data());
    right_multiply(AX.data(), Ysel.data());
    if (generalized) right_multiply(BXs.data(), Ysel.data());
    return 0;
  };

  // Explicit refresh: re-project X against the constraints, recompute B X,
  // B-orthonormalise X (Cholesky QR), recompute A X, Rayleigh-Ritz on X.
  bool fresh = false;
  int since_refresh = 0;
  auto refresh = [&]() -> int {
    for (int pass = 0; pass < 2; ++pass) {
      lobpcg_project_out(Yp, BYp, constraint_rank, n, k, X.data(), nullptr, scratch);
    }
    if (generalized) {
      const int status = B_op.run(k, X.data(), n, BXs.data());
      if (status != 0) return status;
    }
    const int r = lobpcg_cholqr_pass(X.data(), generalized ? BXs.data() : nullptr,
                                     n, k, scratch);
    if (r < 0) return r;
    if (r < k) return -4;
    const int status = A_op.run(k, X.data(), n, AX.data());
    if (status != 0) return status;
    ++(*matvecs_out);
    update_norm_bounds(X.data(), AX.data(), generalized ? BXs.data() : nullptr, k);
    fresh = true;
    since_refresh = 0;
    return rayleigh_ritz_x();
  };

  // Initial block: B-orthonormal basis of the (constraint-projected) start.
  {
    const int r = lobpcg_b_orthonormalize_block(
      metric, n, k, X.data(), generalized ? BXs.data() : nullptr,
      Yp, BYp, constraint_rank, nullptr, nullptr, 0, scratch);
    if (r < 0) return r;
    if (r < k) return -4;
    const int status = A_op.run(k, X.data(), n, AX.data());
    if (status != 0) return status;
    ++(*matvecs_out);
    update_norm_bounds(X.data(), AX.data(), generalized ? BXs.data() : nullptr, k);
    const int rr = rayleigh_ritz_x();
    if (rr != 0) return rr;
    fresh = true;
  }

  int have_p = 0;
  bool stalled = false;
  for (int iter = 0; iter < maxit; ++iter) {
    eigencore_check_interrupt();
    *iterations_out = iter + 1;
    if (!fresh && since_refresh >= kLobpcgRefreshInterval) {
      const int status = refresh();
      if (status != 0) return status;
    }

    int nconv = 0;
    double max_relative = 0.0;
    for (;;) {
      const double* BXp = generalized ? BXs.data() : X.data();
      nconv = 0;
      max_relative = 0.0;
      for (int col = 0; col < k; ++col) {
        const double* x = lobpcg_col(X.data(), n, col);
        const double* ax = lobpcg_col(AX.data(), n, col);
        const double* bx = lobpcg_col(BXp, n, col);
        const double xbx = lobpcg_dot(x, bx, n);
        const double xax = lobpcg_dot(x, ax, n);
        const double lambda = xbx > 0.0 ? xax / xbx : xax;
        values_out[col] = lambda;
        double* r = lobpcg_col(R.data(), n, col);
        for (int row = 0; row < n; ++row) {
          r[row] = ax[row] - lambda * bx[row];
        }
        residuals_out[col] = lobpcg_nrm2(r, n);
        const double scale = (norm_A_lb + fabs(lambda) * norm_B_lb) * lobpcg_nrm2(x, n);
        const double relative = scale > 0.0 ? residuals_out[col] / scale
                                            : (residuals_out[col] > 0.0 ? R_PosInf : 0.0);
        rel[static_cast<size_t>(col)] = relative;
        converged_out[col] = relative <= tol ? 1 : 0;
        if (converged_out[col]) ++nconv;
        if (!(relative <= max_relative)) max_relative = relative;
      }
      const bool stopping = nconv >= k || iter + 1 >= maxit || stalled;
      if (stopping && !fresh) {
        // Never report residuals of recurrence-updated A X / B X: recompute
        // them explicitly before the result leaves the solver.
        const int status = refresh();
        if (status != 0) return status;
        continue;
      }
      break;
    }
    hist_max_residual[iter] = max_relative;
    hist_nconv[iter] = nconv;
    if (nconv >= k || iter + 1 >= maxit || stalled) {
      break;
    }

    // Trial block [W_active, P_active].
    int na = 0;
    for (int col = 0; col < k; ++col) {
      if (!converged_out[col]) active[static_cast<size_t>(na++)] = col;
    }
    for (int a = 0; a < na; ++a) {
      std::memcpy(lobpcg_col(V.data(), n, a),
                  lobpcg_col(R.data(), n, active[static_cast<size_t>(a)]),
                  sizeof(double) * static_cast<size_t>(n));
    }
    if (constraint_rank > 0) {
      // Residual of the constrained problem: remove the Lagrange-multiplier
      // part along B Y, R <- (I - B Y Y') R. Without this the metric image of
      // the constraints leaks into W through (I - Y Y' B) when B != I and
      // the constraints are not invariant, which stalls convergence.
      lobpcg_project_out(BYp, Yp, constraint_rank, n, na, V.data(), nullptr, scratch);
    }
    if (use_tridiagonal_preconditioner) {
      const int status = preconditioner.solve(V.data(), na);
      if (status != 0) return status;
      ++(*preconditioner_calls_out);
    }
    int m = na;
    if (have_p) {
      for (int a = 0; a < na; ++a) {
        std::memcpy(lobpcg_col(V.data(), n, na + a),
                    lobpcg_col(P.data(), n, active[static_cast<size_t>(a)]),
                    sizeof(double) * static_cast<size_t>(n));
      }
      m += na;
    }
    const double* BXp = generalized ? BXs.data() : X.data();
    const int r = lobpcg_b_orthonormalize_block(
      metric, n, m, V.data(), generalized ? BVs.data() : nullptr,
      Yp, BYp, constraint_rank, X.data(), BXp, k, scratch);
    if (r < 0) return r;
    if (r == 0) {
      // No new direction survives: refresh once and stop at the next check.
      stalled = true;
      have_p = 0;
      ++since_refresh;
      fresh = false;
      continue;
    }
    {
      const int status = A_op.run(r, V.data(), n, AV.data());
      if (status != 0) return status;
      ++(*matvecs_out);
    }
    update_norm_bounds(V.data(), AV.data(), generalized ? BVs.data() : nullptr, r);

    // H = [X V]' A [X V] (q x q); the basis is B-orthonormal so the Gram
    // matrix is the identity.
    const int q = k + r;
    *q_rank_final_out = q;
    {
      std::vector<double> H11(static_cast<size_t>(k) * k);
      std::vector<double> H12(static_cast<size_t>(k) * r);
      std::vector<double> H22(static_cast<size_t>(r) * r);
      F77_CALL(dgemm)(&trans_T, &trans_N, &k, &k, &n, &one, X.data(), &n,
                      AX.data(), &n, &zero, H11.data(), &k FCONE FCONE);
      F77_CALL(dgemm)(&trans_T, &trans_N, &k, &r, &n, &one, X.data(), &n,
                      AV.data(), &n, &zero, H12.data(), &k FCONE FCONE);
      F77_CALL(dgemm)(&trans_T, &trans_N, &r, &r, &n, &one, V.data(), &n,
                      AV.data(), &n, &zero, H22.data(), &r FCONE FCONE);
      for (int j = 0; j < q; ++j) {
        for (int i = 0; i <= j; ++i) {
          double h;
          if (j < k) {
            h = 0.5 * (H11[i + static_cast<size_t>(j) * k] + H11[j + static_cast<size_t>(i) * k]);
          } else if (i < k) {
            h = H12[i + static_cast<size_t>(j - k) * k];
          } else {
            const int ii = i - k;
            const int jj = j - k;
            h = 0.5 * (H22[ii + static_cast<size_t>(jj) * r] + H22[jj + static_cast<size_t>(ii) * r]);
          }
          H[i + static_cast<size_t>(j) * q] = h;
          H[j + static_cast<size_t>(i) * q] = h;
        }
      }
    }
    {
      const int status = lobpcg_dsyev(H.data(), q, theta.data(), eig_work);
      if (status != 0) return status;
    }
    selected_ritz_indices(theta.data(), q, k, target_kind, selected.data());
    for (int col = 0; col < k; ++col) {
      std::memcpy(&Ysel[static_cast<size_t>(col) * q],
                  &H[static_cast<size_t>(selected[static_cast<size_t>(col)]) * q],
                  sizeof(double) * static_cast<size_t>(q));
    }
    const double* Yx = Ysel.data();
    const double* Yz = Ysel.data() + k;

    // P = V Y_z; X <- X Y_x + P, and A X, B X by the same coefficients
    // (A X_next = [A X, A V] y, so A is not re-applied to X). A P and B P are
    // not kept: P re-enters through the trial block, whose B image is formed
    // once there and whose A image is applied explicitly.
    auto advance = [&](double* Xb, const double* Vb, double* Pb) {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &r, &one,
                      const_cast<double*>(Vb), &n, const_cast<double*>(Yz), &q,
                      &zero, Xt.data(), &n FCONE FCONE);
      if (Pb != nullptr) {
        std::memcpy(Pb, Xt.data(), sizeof(double) * nk);
      }
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &k, &one,
                      Xb, &n, const_cast<double*>(Yx), &q,
                      &one, Xt.data(), &n FCONE FCONE);
      std::memcpy(Xb, Xt.data(), sizeof(double) * nk);
    };
    advance(X.data(), V.data(), P.data());
    advance(AX.data(), AV.data(), nullptr);
    if (generalized) {
      advance(BXs.data(), BVs.data(), nullptr);
    }
    have_p = 1;
    fresh = false;
    ++since_refresh;
  }

  std::memcpy(X_out, X.data(), sizeof(double) * nk);
  return 0;
}

static SEXP lobpcg_run_native_checked(void* impl,
                                      EigencoreApplyFn apply,
                                      void* b_impl,
                                      EigencoreApplyFn b_apply,
                                      int n,
                                      int k,
                                      int maxit,
                                      int target_kind,
                                      double tol,
                                      const double* start,
                                      const double* lower,
                                      const double* diag,
                                      const double* upper,
                                      const double* constraints,
                                      int constraint_cols,
                                      const char* error_label) {
  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0;
  int matvecs = 0;
  int preconditioner_calls = 0;
  int q_rank = 0;
  int constraints_rank = 0;
  const int status = native_lobpcg_run(
    impl, apply, b_impl, b_apply,
    n, k, maxit, target_kind, tol, start,
    diag != nullptr, lower, diag, upper,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native %s LOBPCG failed with status=%d", error_label, status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_dense(SEXP A_, SEXP k_, SEXP maxit_,
                                       SEXP target_kind_, SEXP tol_,
                                       SEXP start_, SEXP lower_, SEXP diag_,
                                       SEXP upper_, SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  return lobpcg_run_native_checked(
    &impl, eigencore_dense_apply, nullptr, nullptr,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "dense");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_dense_dense_b(SEXP A_, SEXP B_, SEXP k_,
                                               SEXP maxit_, SEXP target_kind_,
                                               SEXP tol_, SEXP start_,
                                               SEXP lower_, SEXP diag_,
                                               SEXP upper_,
                                               SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(B_) || !isReal(start_)) {
    error("A, B, and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue || dimS == R_NilValue) {
    error("A, B, and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n ||
      INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  DenseColumnMajorOperator b_impl = {n, n, REAL(B_)};
  return lobpcg_run_native_checked(
    &impl, eigencore_dense_apply, &b_impl, eigencore_dense_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "dense generalized");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_dense_diagonal_b(SEXP A_, SEXP bdiag_,
                                                  SEXP bunit_, SEXP k_,
                                                  SEXP maxit_,
                                                  SEXP target_kind_,
                                                  SEXP tol_, SEXP start_,
                                                  SEXP lower_, SEXP diag_,
                                                  SEXP upper_,
                                                  SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(bdiag_) || !isReal(start_)) {
    error("A, B diagonal, and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool unit = asLogical(bunit_) == TRUE;
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!unit && LENGTH(bdiag_) != n)) {
    error("non-conformable dense/diagonal generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  DiagonalOperator b_impl = {n, REAL(bdiag_), unit};
  return lobpcg_run_native_checked(
    &impl, eigencore_dense_apply, &b_impl, eigencore_diagonal_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "dense/diagonal generalized");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_dense_csc_b(SEXP A_, SEXP bi_, SEXP bp_,
                                             SEXP bx_, SEXP bdim_, SEXP k_,
                                             SEXP maxit_, SEXP target_kind_,
                                             SEXP tol_, SEXP start_,
                                             SEXP lower_, SEXP diag_,
                                             SEXP upper_,
                                             SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isInteger(bi_) || !isInteger(bp_) ||
      !isReal(bx_) || !isInteger(bdim_) || !isReal(start_)) {
    error("invalid dense/CSC generalized LOBPCG inputs");
  }
  eigencore_validate_csc_structure(bi_, bp_, bx_, bdim_, "LOBPCG B");
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || LENGTH(bdim_) != 2) {
    error("A and start must be matrices and B dim must have length 2");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(bdim_)[0] != n || INTEGER(bdim_)[1] != n ||
      INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense/CSC generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  CSCOperator b_impl = {
    INTEGER(bdim_)[0], INTEGER(bdim_)[1], INTEGER(bi_), INTEGER(bp_), REAL(bx_)
  };
  return lobpcg_run_native_checked(
    &impl, eigencore_dense_apply, &b_impl, eigencore_csc_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "dense/CSC generalized");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_csc_diagonal_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                                SEXP adim_, SEXP bdiag_,
                                                SEXP bunit_, SEXP k_,
                                                SEXP maxit_,
                                                SEXP target_kind_,
                                                SEXP tol_, SEXP start_,
                                                SEXP lower_, SEXP diag_,
                                                SEXP upper_,
                                                SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isReal(bdiag_) || !isReal(start_)) {
    error("invalid CSC/diagonal generalized LOBPCG inputs");
  }
  eigencore_validate_csc_structure(ai_, ap_, ax_, adim_, "LOBPCG A");
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool unit = asLogical(bunit_) == TRUE;
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!unit && LENGTH(bdiag_) != n)) {
    error("non-conformable CSC/diagonal generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  CSCOperator impl = {
    INTEGER(adim_)[0], INTEGER(adim_)[1], INTEGER(ai_), INTEGER(ap_), REAL(ax_)
  };
  DiagonalOperator b_impl = {n, REAL(bdiag_), unit};
  return lobpcg_run_native_checked(
    &impl, eigencore_csc_apply, &b_impl, eigencore_diagonal_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "CSC/diagonal generalized");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_csc_csc_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                           SEXP adim_, SEXP bi_, SEXP bp_,
                                           SEXP bx_, SEXP bdim_, SEXP k_,
                                           SEXP maxit_, SEXP target_kind_,
                                           SEXP tol_, SEXP start_,
                                           SEXP lower_, SEXP diag_,
                                           SEXP upper_,
                                           SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isInteger(bi_) || !isInteger(bp_) ||
      !isReal(bx_) || !isInteger(bdim_) || !isReal(start_)) {
    error("invalid CSC/CSC generalized LOBPCG inputs");
  }
  eigencore_validate_csc_structure(ai_, ap_, ax_, adim_, "LOBPCG A");
  eigencore_validate_csc_structure(bi_, bp_, bx_, bdim_, "LOBPCG B");
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2 || LENGTH(bdim_) != 2) {
    error("start must be a matrix and A/B dims must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(adim_)[1] != n || INTEGER(bdim_)[0] != n || INTEGER(bdim_)[1] != n ||
      INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC/CSC generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  CSCOperator impl = {
    INTEGER(adim_)[0], INTEGER(adim_)[1], INTEGER(ai_), INTEGER(ap_), REAL(ax_)
  };
  CSCOperator b_impl = {
    INTEGER(bdim_)[0], INTEGER(bdim_)[1], INTEGER(bi_), INTEGER(bp_), REAL(bx_)
  };
  return lobpcg_run_native_checked(
    &impl, eigencore_csc_apply, &b_impl, eigencore_csc_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "CSC/CSC generalized");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_diagonal_diagonal_b(SEXP adiag_, SEXP aunit_,
                                                     SEXP adim_, SEXP bdiag_,
                                                     SEXP bunit_, SEXP k_,
                                                     SEXP maxit_,
                                                     SEXP target_kind_,
                                                     SEXP tol_, SEXP start_,
                                                     SEXP lower_, SEXP diag_,
                                                     SEXP upper_,
                                                     SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(adiag_) || !isInteger(adim_) || !isReal(bdiag_) ||
      !isReal(start_)) {
    error("invalid diagonal/diagonal generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool a_unit = asLogical(aunit_) == TRUE;
  const bool b_unit = asLogical(bunit_) == TRUE;
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!a_unit && LENGTH(adiag_) != n) || (!b_unit && LENGTH(bdiag_) != n)) {
    error("non-conformable diagonal/diagonal generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  DiagonalOperator impl = {n, REAL(adiag_), a_unit};
  DiagonalOperator b_impl = {n, REAL(bdiag_), b_unit};
  return lobpcg_run_native_checked(
    &impl, eigencore_diagonal_apply, &b_impl, eigencore_diagonal_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "diagonal/diagonal generalized");
  EIGENCORE_ENTRY_END
}

static SEXP lobpcg_run_matrix_free_b(void* impl,
                                     EigencoreApplyFn apply,
                                     int n,
                                     int k,
                                     int maxit,
                                     int target_kind,
                                     double tol,
                                     const double* start,
                                     SEXP B_apply_,
                                     SEXP lower_,
                                     SEXP diag_,
                                     SEXP upper_,
                                     SEXP constraints_,
                                     const char* label) {
  if (TYPEOF(B_apply_) != CLOSXP) {
    error("matrix-free B operator apply must be an R closure");
  }
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  RApplyOperator b_impl = {n, n, B_apply_, R_NilValue};
  char error_label[128];
  std::snprintf(error_label, sizeof(error_label), "%s matrix-free-B", label);
  return lobpcg_run_native_checked(
    impl, apply, &b_impl, eigencore_r_operator_apply,
    n, k, maxit, target_kind, tol, start,
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, error_label);
}

extern "C" SEXP eigencore_lobpcg_dense_operator_b(SEXP A_, SEXP B_apply_,
                                                  SEXP k_, SEXP maxit_,
                                                  SEXP target_kind_,
                                                  SEXP tol_, SEXP start_,
                                                  SEXP lower_, SEXP diag_,
                                                  SEXP upper_,
                                                  SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense/matrix-free generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  return lobpcg_run_matrix_free_b(
    &impl, eigencore_dense_apply, n, k, maxit, target_kind, tol, REAL(start_),
    B_apply_, lower_, diag_, upper_, constraints_, "dense");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_csc_operator_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                                SEXP adim_, SEXP B_apply_,
                                                SEXP k_, SEXP maxit_,
                                                SEXP target_kind_,
                                                SEXP tol_, SEXP start_,
                                                SEXP lower_, SEXP diag_,
                                                SEXP upper_,
                                                SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isReal(start_)) {
    error("invalid CSC/matrix-free generalized LOBPCG inputs");
  }
  eigencore_validate_csc_structure(ai_, ap_, ax_, adim_, "LOBPCG A");
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC/matrix-free generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  CSCOperator impl = {
    INTEGER(adim_)[0], INTEGER(adim_)[1], INTEGER(ai_), INTEGER(ap_), REAL(ax_)
  };
  return lobpcg_run_matrix_free_b(
    &impl, eigencore_csc_apply, n, k, maxit, target_kind, tol, REAL(start_),
    B_apply_, lower_, diag_, upper_, constraints_, "CSC");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_diagonal_operator_b(SEXP adiag_,
                                                     SEXP aunit_,
                                                     SEXP adim_,
                                                     SEXP B_apply_,
                                                     SEXP k_,
                                                     SEXP maxit_,
                                                     SEXP target_kind_,
                                                     SEXP tol_,
                                                     SEXP start_,
                                                     SEXP lower_,
                                                     SEXP diag_,
                                                     SEXP upper_,
                                                     SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(adiag_) || !isInteger(adim_) || !isReal(start_)) {
    error("invalid diagonal/matrix-free generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool a_unit = asLogical(aunit_) == TRUE;
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!a_unit && LENGTH(adiag_) != n)) {
    error("non-conformable diagonal/matrix-free generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  DiagonalOperator impl = {n, REAL(adiag_), a_unit};
  return lobpcg_run_matrix_free_b(
    &impl, eigencore_diagonal_apply, n, k, maxit, target_kind, tol, REAL(start_),
    B_apply_, lower_, diag_, upper_, constraints_, "diagonal");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                     SEXP k_, SEXP maxit_, SEXP target_kind_,
                                     SEXP tol_, SEXP start_, SEXP lower_,
                                     SEXP diag_, SEXP upper_,
                                     SEXP constraints_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC LOBPCG inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "LOBPCG");
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return lobpcg_run_native_checked(
    &impl, eigencore_csc_apply, nullptr, nullptr,
    n, k, maxit, target_kind, tol, REAL(start_),
    LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols, "CSC");
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_lobpcg_csc_shifted_tridiagonal(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP k_, SEXP maxit_,
    SEXP target_kind_, SEXP tol_, SEXP start_, SEXP shift_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC shifted-tridiagonal LOBPCG inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "shifted-tridiagonal LOBPCG");
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC shifted-tridiagonal LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const double shift = asReal(shift_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  if (!R_FINITE(shift) || shift < 0.0) {
    error("shift must be a finite non-negative scalar");
  }

  std::vector<double> lower(static_cast<size_t>(n > 1 ? n - 1 : 1), 0.0);
  std::vector<double> diag(static_cast<size_t>(n), 0.0);
  std::vector<double> upper(static_cast<size_t>(n > 1 ? n - 1 : 1), 0.0);
  const int tri_status = extract_shifted_symmetric_tridiagonal_from_csc(
    INTEGER(i_), INTEGER(p_), REAL(x_), n, shift,
    lower.data(), diag.data(), upper.data()
  );
  if (tri_status == -6) {
    error("CSC matrix is not tridiagonal");
  }
  if (tri_status == -7) {
    error("CSC matrix is not symmetric tridiagonal");
  }
  if (tri_status != 0) {
    error("failed to extract shifted tridiagonal preconditioner");
  }

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return lobpcg_run_native_checked(
    &impl, eigencore_csc_apply, nullptr, nullptr,
    n, k, maxit, target_kind, tol, REAL(start_),
    lower.data(), diag.data(), upper.data(),
    nullptr, 0, "CSC shifted-tridiagonal");
  EIGENCORE_ENTRY_END
}
