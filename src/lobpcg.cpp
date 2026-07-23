#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"
#include "eigencore_common.h"
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

static int trl_orthogonalise(const double* V_locked, int n_locked,
                             const double* V_active, int m_active,
                             double* z, double* tmp, int n,
                             int passes = 2) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  int incx = 1;
  for (int pass = 0; pass < passes; ++pass) {
    if (n_locked > 0) {
      F77_CALL(dgemv)(&trans_T, &n, &n_locked, &one,
                      V_locked, &n, z, &incx,
                      &zero, tmp, &incx FCONE);
      F77_CALL(dgemv)(&trans_N, &n, &n_locked, &minus_one,
                      V_locked, &n, tmp, &incx,
                      &one, z, &incx FCONE);
    }
    if (m_active > 0) {
      F77_CALL(dgemv)(&trans_T, &n, &m_active, &one,
                      V_active, &n, z, &incx,
                      &zero, tmp, &incx FCONE);
      F77_CALL(dgemv)(&trans_N, &n, &m_active, &minus_one,
                      V_active, &n, tmp, &incx,
                      &one, z, &incx FCONE);
    }
  }
  return 0;
}

static double trl_norm2(const double* x, int n) {
  long double sum = 0.0L;
  for (int i = 0; i < n; ++i) {
    sum += static_cast<long double>(x[i]) * x[i];
  }
  return sqrt(static_cast<double>(sum));
}

static int trl_dsyev_query(int m_max) {
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork_query = -1;
  double work_query = 0.0;
  double fake = 0.0;
  double fake_w = 0.0;
  int m = m_max;
  F77_CALL(dsyev)(&jobz, &uplo, &m, &fake, &m, &fake_w,
                  &work_query, &lwork_query, &info FCONE FCONE);
  if (info != 0) {
    return 3 * m_max;
  }
  return static_cast<int>(work_query);
}

static int lobpcg_orthonormalize(const double* S, int n, int cols, double tol,
                                 double* Q, double* tmp) {
  int rank = 0;
  for (int col = 0; col < cols; ++col) {
    double* q_col = Q + static_cast<int64_t>(rank) * n;
    std::memcpy(q_col, S + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
    trl_orthogonalise(nullptr, 0, Q, rank, q_col, tmp, n);
    const double q_norm = trl_norm2(q_col, n);
    if (q_norm <= tol) {
      continue;
    }
    const double inv = 1.0 / q_norm;
    for (int row = 0; row < n; ++row) {
      q_col[row] *= inv;
    }
    ++rank;
  }
  return rank;
}

static int lobpcg_b_orthonormalize_apply(void* b_impl,
                                         EigencoreApplyFn b_apply,
                                         const double* S,
                                         int n,
                                         int cols,
                                         double tol,
                                         double* Q,
                                         double* BQ,
                                         double* coeff,
                                         double* bq_col,
                                         EigencoreWorkspace* workspace) {
  int rank = 0;
  for (int col = 0; col < cols; ++col) {
    double* q_col = Q + static_cast<int64_t>(rank) * n;
    std::memcpy(q_col, S + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
    int status = b_apply(b_impl, EIGENCORE_TRANSPOSE_NONE, 1, q_col, n,
                         1.0, 0.0, bq_col, n, workspace);
    if (status != 0) {
      return status;
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int prev = 0; prev < rank; ++prev) {
        const double* q_prev = Q + static_cast<int64_t>(prev) * n;
        long double dot = 0.0L;
        for (int row = 0; row < n; ++row) {
          dot += static_cast<long double>(q_prev[row]) * bq_col[row];
        }
        coeff[prev] = static_cast<double>(dot);
      }
      for (int prev = 0; prev < rank; ++prev) {
        const double c = coeff[prev];
        if (c == 0.0) {
          continue;
        }
        const double* q_prev = Q + static_cast<int64_t>(prev) * n;
        const double* bq_prev = BQ + static_cast<int64_t>(prev) * n;
        for (int row = 0; row < n; ++row) {
          q_col[row] -= c * q_prev[row];
          bq_col[row] -= c * bq_prev[row];
        }
      }
    }

    long double norm_sq = 0.0L;
    for (int row = 0; row < n; ++row) {
      norm_sq += static_cast<long double>(q_col[row]) * bq_col[row];
    }
    if (norm_sq <= 0.0L) {
      continue;
    }
    const double q_norm = sqrt(static_cast<double>(norm_sq));
    if (q_norm <= tol) {
      continue;
    }
    const double inv = 1.0 / q_norm;
    double* bq_dst = BQ + static_cast<int64_t>(rank) * n;
    for (int row = 0; row < n; ++row) {
      q_col[row] *= inv;
      bq_dst[row] = bq_col[row] * inv;
    }
    ++rank;
  }
  return rank;
}

static int lobpcg_apply_tridiagonal_preconditioner(
    const double* lower,
    const double* diag,
    const double* upper,
    int n,
    int cols,
    const double* R,
    double* W,
    double* cprime,
    double* dprime) {
  if (n < 1) {
    return 0;
  }
  if (fabs(diag[0]) <= DBL_EPSILON) {
    return -5;
  }
  dprime[0] = diag[0];
  if (n > 1) {
    cprime[0] = upper[0] / dprime[0];
  }
  for (int i = 1; i < n; ++i) {
    dprime[i] = diag[i] - lower[i - 1] * cprime[i - 1];
    if (fabs(dprime[i]) <= DBL_EPSILON) {
      return -5;
    }
    if (i < n - 1) {
      cprime[i] = upper[i] / dprime[i];
    }
  }

  for (int col = 0; col < cols; ++col) {
    const double* rhs = R + static_cast<int64_t>(col) * n;
    double* out = W + static_cast<int64_t>(col) * n;
    out[0] = rhs[0] / dprime[0];
    for (int i = 1; i < n; ++i) {
      out[i] = (rhs[i] - lower[i - 1] * out[i - 1]) / dprime[i];
    }
    for (int i = n - 2; i >= 0; --i) {
      out[i] -= cprime[i] * out[i + 1];
    }
  }
  return 0;
}

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

static int lobpcg_project_constraints_apply(void* b_impl,
                                            EigencoreApplyFn b_apply,
                                            const double* Qc,
                                            int constraint_rank,
                                            int n,
                                            int cols,
                                            double* X,
                                            double* coeff,
                                            double* BX_work,
                                            EigencoreWorkspace* workspace) {
  if (constraint_rank <= 0 || cols <= 0) {
    return 0;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  const double* metric_X = X;
  if (b_apply != nullptr) {
    const int status = b_apply(b_impl, EIGENCORE_TRANSPOSE_NONE, cols, X, n,
                               1.0, 0.0, BX_work, n, workspace);
    if (status != 0) {
      return status;
    }
    metric_X = BX_work;
  }
  F77_CALL(dgemm)(&trans_T, &trans_N, &constraint_rank, &cols, &n,
                  &one, const_cast<double*>(Qc), &n,
                  const_cast<double*>(metric_X), &n,
                  &zero, coeff, &constraint_rank FCONE FCONE);
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &cols, &constraint_rank,
                  &minus_one, const_cast<double*>(Qc), &n,
                  coeff, &constraint_rank,
                  &one, X, &n FCONE FCONE);
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
  const int max_trial_cols = 3 * k;
  const int tmp_cols = constraint_cols > max_trial_cols ? constraint_cols : max_trial_cols;
  const size_t nk = static_cast<size_t>(n) * static_cast<size_t>(k);
  const size_t nt = static_cast<size_t>(n) * static_cast<size_t>(max_trial_cols);
  double* X = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* Xnext = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* P = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* AX = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* BX = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* BXnext = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* R = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* W = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* S = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* Q = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* AQ = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* BQ = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* H = static_cast<double*>(std::calloc(static_cast<size_t>(max_trial_cols) * max_trial_cols, sizeof(double)));
  double* selected_vectors = static_cast<double*>(std::calloc(static_cast<size_t>(max_trial_cols) * k, sizeof(double)));
  double* theta = static_cast<double*>(std::calloc(static_cast<size_t>(max_trial_cols), sizeof(double)));
  double* tmp = static_cast<double*>(std::calloc(static_cast<size_t>(tmp_cols), sizeof(double)));
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>(max_trial_cols), sizeof(int)));
  double* cprime = static_cast<double*>(std::calloc(static_cast<size_t>(n > 1 ? n - 1 : 1), sizeof(double)));
  double* dprime = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  const int dsyev_lwork_query = trl_dsyev_query(max_trial_cols);
  int dsyev_lwork = dsyev_lwork_query > 0 ? dsyev_lwork_query : 3 * max_trial_cols;
  double* dsyev_work = static_cast<double*>(std::calloc(static_cast<size_t>(dsyev_lwork), sizeof(double)));
  if (X == nullptr || Xnext == nullptr || P == nullptr || AX == nullptr ||
      BX == nullptr || BXnext == nullptr || R == nullptr || W == nullptr ||
      S == nullptr || Q == nullptr || AQ == nullptr || BQ == nullptr ||
      H == nullptr || selected_vectors == nullptr ||
      theta == nullptr || tmp == nullptr || selected == nullptr ||
      cprime == nullptr || dprime == nullptr || dsyev_work == nullptr) {
    std::free(X); std::free(Xnext); std::free(P); std::free(AX); std::free(R);
    std::free(BX); std::free(BXnext); std::free(W); std::free(S);
    std::free(Q); std::free(AQ); std::free(BQ); std::free(H);
    std::free(selected_vectors); std::free(theta); std::free(tmp);
    std::free(selected); std::free(cprime); std::free(dprime); std::free(dsyev_work);
    return -2;
  }
  auto cleanup = [&]() {
    std::free(X); std::free(Xnext); std::free(P); std::free(AX); std::free(R);
    std::free(BX); std::free(BXnext); std::free(W); std::free(S);
    std::free(Q); std::free(AQ); std::free(BQ); std::free(H);
    std::free(selected_vectors); std::free(theta); std::free(tmp);
    std::free(selected); std::free(cprime); std::free(dprime); std::free(dsyev_work);
  };

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const bool generalized = b_apply != nullptr;
  int constraint_rank = 0;
  *constraints_rank_out = 0;
  std::vector<double> constraint_q;
  std::vector<double> constraint_bq;
  std::vector<double> constraint_bscratch;
  std::vector<double> constraint_coeff;
  std::vector<double> constraint_bx;
  if (constraint_cols > 0) {
    constraint_q.assign(static_cast<size_t>(n) * constraint_cols, 0.0);
    std::memcpy(constraint_q.data(), constraints,
                sizeof(double) * static_cast<size_t>(n) * constraint_cols);
    constraint_coeff.assign(static_cast<size_t>(constraint_cols) * max_trial_cols, 0.0);
    if (generalized) {
      constraint_bq.assign(static_cast<size_t>(n) * constraint_cols, 0.0);
      constraint_bscratch.assign(static_cast<size_t>(n), 0.0);
      constraint_bx.assign(static_cast<size_t>(n) * max_trial_cols, 0.0);
      constraint_rank = lobpcg_b_orthonormalize_apply(
        b_impl, b_apply, constraint_q.data(), n, constraint_cols,
        100.0 * DBL_EPSILON, constraint_q.data(), constraint_bq.data(), tmp,
        constraint_bscratch.data(), &workspace);
    } else {
      constraint_rank = lobpcg_orthonormalize(
        constraint_q.data(), n, constraint_cols, 100.0 * DBL_EPSILON,
        constraint_q.data(), tmp);
    }
    if (constraint_rank < 0) {
      cleanup();
      return constraint_rank;
    }
    if (constraint_rank + k > n) {
      cleanup();
      return -9;
    }
    *constraints_rank_out = constraint_rank;
  }
  std::memcpy(S, start, sizeof(double) * nk);
  if (constraint_rank > 0) {
    const int status = lobpcg_project_constraints_apply(
      b_impl, b_apply, constraint_q.data(), constraint_rank, n, k, S,
      constraint_coeff.data(),
      generalized ? constraint_bx.data() : nullptr,
      &workspace);
    if (status != 0) {
      cleanup();
      return status;
    }
  }
  int q_rank = generalized
    ? lobpcg_b_orthonormalize_apply(b_impl, b_apply, S, n, k,
                                    100.0 * DBL_EPSILON, X, BX, tmp,
                                    BQ, &workspace)
    : lobpcg_orthonormalize(S, n, k, 100.0 * DBL_EPSILON, X, tmp);
  if (q_rank < k) {
    cleanup();
    return -4;
  }
  if (!generalized) {
    std::memcpy(BX, X, sizeof(double) * nk);
  }

  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int have_p = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *preconditioner_calls_out = 0;
  *q_rank_final_out = k;

  for (int iter = 0; iter < maxit; ++iter) {
    *iterations_out = iter + 1;
    int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, k, X, n,
                       1.0, 0.0, AX, n, &workspace);
    if (status != 0) {
      cleanup();
      return status;
    }
    ++(*matvecs_out);

    int nconv = 0;
    double max_relative = 0.0;
    for (int col = 0; col < k; ++col) {
      long double lambda = 0.0L;
      for (int row = 0; row < n; ++row) {
        lambda += static_cast<long double>(X[row + static_cast<int64_t>(col) * n]) *
                  AX[row + static_cast<int64_t>(col) * n];
      }
      values_out[col] = static_cast<double>(lambda);
      long double rn = 0.0L;
      for (int row = 0; row < n; ++row) {
        const double res = AX[row + static_cast<int64_t>(col) * n] -
                           values_out[col] * BX[row + static_cast<int64_t>(col) * n];
        R[row + static_cast<int64_t>(col) * n] = res;
        rn += static_cast<long double>(res) * res;
      }
      residuals_out[col] = sqrt(static_cast<double>(rn));
      const double scale = fabs(values_out[col]) > 1.0 ? fabs(values_out[col]) : 1.0;
      const double rel = residuals_out[col] / scale;
      converged_out[col] = rel <= tol ? 1 : 0;
      if (converged_out[col]) ++nconv;
      if (rel > max_relative) max_relative = rel;
    }
    hist_max_residual[iter] = max_relative;
    hist_nconv[iter] = nconv;
    if (nconv >= k || iter + 1 >= maxit) {
      break;
    }

    if (use_tridiagonal_preconditioner) {
      status = lobpcg_apply_tridiagonal_preconditioner(
        lower, diag, upper, n, k, R, W, cprime, dprime
      );
      if (status != 0) {
        cleanup();
        return status;
      }
      ++(*preconditioner_calls_out);
    } else {
      std::memcpy(W, R, sizeof(double) * nk);
    }

    int trial_cols = 0;
    std::memcpy(S + static_cast<int64_t>(trial_cols) * n, X, sizeof(double) * nk);
    trial_cols += k;
    std::memcpy(S + static_cast<int64_t>(trial_cols) * n, W, sizeof(double) * nk);
    trial_cols += k;
    if (have_p) {
      std::memcpy(S + static_cast<int64_t>(trial_cols) * n, P, sizeof(double) * nk);
      trial_cols += k;
    }
    if (constraint_rank > 0) {
      status = lobpcg_project_constraints_apply(
        b_impl, b_apply, constraint_q.data(), constraint_rank, n, trial_cols, S,
        constraint_coeff.data(),
        generalized ? constraint_bx.data() : nullptr,
        &workspace);
      if (status != 0) {
        cleanup();
        return status;
      }
    }
    std::memset(Q, 0, sizeof(double) * nt);
    q_rank = generalized
      ? lobpcg_b_orthonormalize_apply(b_impl, b_apply, S, n, trial_cols,
                                      100.0 * DBL_EPSILON, Q, BQ, tmp,
                                      BXnext, &workspace)
      : lobpcg_orthonormalize(S, n, trial_cols, 100.0 * DBL_EPSILON, Q, tmp);
    if (q_rank < k) {
      cleanup();
      return -4;
    }
    *q_rank_final_out = q_rank;

    status = apply(impl, EIGENCORE_TRANSPOSE_NONE, q_rank, Q, n,
                   1.0, 0.0, AQ, n, &workspace);
    if (status != 0) {
      cleanup();
      return status;
    }
    ++(*matvecs_out);

    F77_CALL(dgemm)(&trans_T, &trans_N, &q_rank, &q_rank, &n,
                    &one, Q, &n, AQ, &n,
                    &zero, H, &q_rank FCONE FCONE);
    for (int i = 0; i < q_rank; ++i) {
      for (int j = i + 1; j < q_rank; ++j) {
        const double avg = 0.5 * (H[i + j * q_rank] + H[j + i * q_rank]);
        H[i + j * q_rank] = avg;
        H[j + i * q_rank] = avg;
      }
    }
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    int lwork = dsyev_lwork;
    F77_CALL(dsyev)(&jobz, &uplo, &q_rank, H, &q_rank, theta,
                    dsyev_work, &lwork, &info FCONE FCONE);
    if (info != 0) {
      cleanup();
      return -3;
    }
    selected_ritz_indices(theta, q_rank, k, target_kind, selected);
    for (int col = 0; col < k; ++col) {
      const int idx = selected[col];
      for (int row = 0; row < q_rank; ++row) {
        selected_vectors[row + static_cast<int64_t>(col) * q_rank] =
          H[row + static_cast<int64_t>(idx) * q_rank];
      }
    }
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &q_rank,
                    &one, Q, &n, selected_vectors, &q_rank,
                    &zero, Xnext, &n FCONE FCONE);
    if (generalized) {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &q_rank,
                      &one, BQ, &n, selected_vectors, &q_rank,
                      &zero, BXnext, &n FCONE FCONE);
    }
    for (size_t pos = 0; pos < nk; ++pos) {
      P[pos] = Xnext[pos] - X[pos];
      X[pos] = Xnext[pos];
      if (generalized) {
        BX[pos] = BXnext[pos];
      } else {
        BX[pos] = Xnext[pos];
      }
    }
    have_p = 1;
  }

  std::memcpy(X_out, X, sizeof(double) * nk);
  cleanup();
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
}

extern "C" SEXP eigencore_lobpcg_dense_dense_b(SEXP A_, SEXP B_, SEXP k_,
                                               SEXP maxit_, SEXP target_kind_,
                                               SEXP tol_, SEXP start_,
                                               SEXP lower_, SEXP diag_,
                                               SEXP upper_,
                                               SEXP constraints_) {
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
}

extern "C" SEXP eigencore_lobpcg_dense_diagonal_b(SEXP A_, SEXP bdiag_,
                                                  SEXP bunit_, SEXP k_,
                                                  SEXP maxit_,
                                                  SEXP target_kind_,
                                                  SEXP tol_, SEXP start_,
                                                  SEXP lower_, SEXP diag_,
                                                  SEXP upper_,
                                                  SEXP constraints_) {
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
}

extern "C" SEXP eigencore_lobpcg_dense_csc_b(SEXP A_, SEXP bi_, SEXP bp_,
                                             SEXP bx_, SEXP bdim_, SEXP k_,
                                             SEXP maxit_, SEXP target_kind_,
                                             SEXP tol_, SEXP start_,
                                             SEXP lower_, SEXP diag_,
                                             SEXP upper_,
                                             SEXP constraints_) {
  if (!isReal(A_) || !isInteger(bi_) || !isInteger(bp_) ||
      !isReal(bx_) || !isInteger(bdim_) || !isReal(start_)) {
    error("invalid dense/CSC generalized LOBPCG inputs");
  }
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
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isReal(bdiag_) || !isReal(start_)) {
    error("invalid CSC/diagonal generalized LOBPCG inputs");
  }
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
}

extern "C" SEXP eigencore_lobpcg_csc_csc_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                           SEXP adim_, SEXP bi_, SEXP bp_,
                                           SEXP bx_, SEXP bdim_, SEXP k_,
                                           SEXP maxit_, SEXP target_kind_,
                                           SEXP tol_, SEXP start_,
                                           SEXP lower_, SEXP diag_,
                                           SEXP upper_,
                                           SEXP constraints_) {
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isInteger(bi_) || !isInteger(bp_) ||
      !isReal(bx_) || !isInteger(bdim_) || !isReal(start_)) {
    error("invalid CSC/CSC generalized LOBPCG inputs");
  }
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
}

extern "C" SEXP eigencore_lobpcg_csc_operator_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                                SEXP adim_, SEXP B_apply_,
                                                SEXP k_, SEXP maxit_,
                                                SEXP target_kind_,
                                                SEXP tol_, SEXP start_,
                                                SEXP lower_, SEXP diag_,
                                                SEXP upper_,
                                                SEXP constraints_) {
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isReal(start_)) {
    error("invalid CSC/matrix-free generalized LOBPCG inputs");
  }
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
}

extern "C" SEXP eigencore_lobpcg_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                     SEXP k_, SEXP maxit_, SEXP target_kind_,
                                     SEXP tol_, SEXP start_, SEXP lower_,
                                     SEXP diag_, SEXP upper_,
                                     SEXP constraints_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC LOBPCG inputs");
  }
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
}

extern "C" SEXP eigencore_lobpcg_csc_shifted_tridiagonal(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP k_, SEXP maxit_,
    SEXP target_kind_, SEXP tol_, SEXP start_, SEXP shift_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC shifted-tridiagonal LOBPCG inputs");
  }
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
}
