#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <vector>
#include "eigencore_common.h"
#include "native_operators.h"
#include "certificates.h"
#include "projection/golub_kahan_ritz.h"
#include "scalar_krylov.h"
#include "block_golub_kahan_basis.h"

static double max_orthogonality_loss(const double* gram, int k) {
  double loss = 0.0;
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row < k; ++row) {
      const double target = (row == col) ? 1.0 : 0.0;
      const double err = fabs(gram[row + col * k] - target);
      if (err > loss) {
        loss = err;
      }
    }
  }
  return loss;
}

static void symmetrize_packed_square(double* A, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double avg = 0.5 * (A[i + j * n] + A[j + i * n]);
      A[i + j * n] = avg;
      A[j + i * n] = avg;
    }
  }
}

static SEXP block_golub_kahan_fit_pack(int n,
                                       int m,
                                       double* V,
                                       double* AV,
                                       int active_v,
                                       int active_u,
                                       int iterations,
                                       int matvecs,
                                       int ortho_passes,
                                       int cached_start_used,
                                       int rank,
                                       int target_kind,
                                       double stage_native_iteration_seconds) {
  auto ritz_timer = native_timer_now();
  SEXP ritz_ = PROTECT(eigencore_block_golub_kahan_ritz_from_ptr(
    V, n, AV, m, active_v, rank, target_kind
  ));
  const double stage_ritz_seconds = native_timer_elapsed(ritz_timer);
  SEXP stage_ = PROTECT(allocVector(REALSXP, 2));
  REAL(stage_)[0] = stage_native_iteration_seconds;
  REAL(stage_)[1] = stage_ritz_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(stage_names_, 0, mkChar("native_iteration"));
  SET_STRING_ELT(stage_names_, 1, mkChar("ritz"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 12));
  for (int i = 0; i < 5; ++i) {
    SET_VECTOR_ELT(out_, i, VECTOR_ELT(ritz_, i));
  }
  SET_VECTOR_ELT(out_, 5, ScalarInteger(active_v));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(active_u));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(ortho_passes));
  SET_VECTOR_ELT(out_, 10, ScalarLogical(cached_start_used != 0));
  SET_VECTOR_ELT(out_, 11, stage_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 12));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("Avectors"));
  SET_STRING_ELT(names_, 4, mkChar("coefficients"));
  SET_STRING_ELT(names_, 5, mkChar("active_cols"));
  SET_STRING_ELT(names_, 6, mkChar("active_left_cols"));
  SET_STRING_ELT(names_, 7, mkChar("iterations"));
  SET_STRING_ELT(names_, 8, mkChar("matvecs"));
  SET_STRING_ELT(names_, 9, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 10, mkChar("cached_start_used"));
  SET_STRING_ELT(names_, 11, mkChar("stage_seconds"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(5);
  return out_;
}

static SEXP block_golub_kahan_retained_cycle_pack(int n,
                                                  int m,
                                                  void* impl,
                                                  EigencoreApplyFn apply,
                                                  double norm_A,
                                                  double tol,
                                                  double* V,
                                                  double* AV,
                                                  int active_v,
                                                  int active_u,
                                                  int iterations,
                                                  int matvecs,
                                                  int ortho_passes,
                                                  int cached_start_used,
                                                  int requested_rank,
                                                  int active_rank,
                                                  int target_kind,
                                                  double stage_native_iteration_seconds,
                                                  double stage_restart_seconds,
                                                  SEXP attempt_history_,
                                                  double native_workspace_bytes,
                                                  const double* locked_d,
                                                  const double* locked_u,
                                                  const double* locked_v,
                                                  const double* locked_av,
                                                  int locked_count) {
  (void)requested_rank;
  if (active_rank < 0) {
    active_rank = 0;
  }
  SEXP fit_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, V, AV, active_v, active_u, iterations, matvecs, ortho_passes,
    cached_start_used, active_rank, target_kind, stage_native_iteration_seconds
  ));
  SEXP d_ = VECTOR_ELT(fit_, 0);
  SEXP u_ = VECTOR_ELT(fit_, 1);
  SEXP v_ = VECTOR_ELT(fit_, 2);
  SEXP avectors_ = VECTOR_ELT(fit_, 3);
  SEXP coeff_ = VECTOR_ELT(fit_, 4);
  SEXP combined_d_ = d_;
  SEXP combined_u_ = u_;
  SEXP combined_v_ = v_;
  SEXP combined_av_ = avectors_;
  SEXP combined_coeff_ = coeff_;
  const int active_count = LENGTH(d_);
  const int total_count = locked_count + active_count;
  if (locked_count > 0) {
    combined_d_ = PROTECT(allocVector(REALSXP, total_count));
    combined_u_ = PROTECT(allocMatrix(REALSXP, m, total_count));
    combined_v_ = PROTECT(allocMatrix(REALSXP, n, total_count));
    combined_av_ = PROTECT(allocMatrix(REALSXP, m, total_count));
    combined_coeff_ = PROTECT(allocMatrix(REALSXP, active_v, total_count));
    std::memset(REAL(combined_coeff_), 0,
                sizeof(double) * static_cast<size_t>(active_v) * total_count);
    for (int col = 0; col < locked_count; ++col) {
      REAL(combined_d_)[col] = locked_d[col];
      std::memcpy(REAL(combined_u_) + static_cast<int64_t>(col) * m,
                  locked_u + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
      std::memcpy(REAL(combined_v_) + static_cast<int64_t>(col) * n,
                  locked_v + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
      std::memcpy(REAL(combined_av_) + static_cast<int64_t>(col) * m,
                  locked_av + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
    }
    for (int col = 0; col < active_count; ++col) {
      const int dst_col = locked_count + col;
      REAL(combined_d_)[dst_col] = REAL(d_)[col];
      std::memcpy(REAL(combined_u_) + static_cast<int64_t>(dst_col) * m,
                  REAL(u_) + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
      std::memcpy(REAL(combined_v_) + static_cast<int64_t>(dst_col) * n,
                  REAL(v_) + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
      std::memcpy(REAL(combined_av_) + static_cast<int64_t>(dst_col) * m,
                  REAL(avectors_) + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
      for (int row = 0; row < active_v; ++row) {
        REAL(combined_coeff_)[row + static_cast<int64_t>(dst_col) * active_v] =
          REAL(coeff_)[row + static_cast<int64_t>(col) * active_v];
      }
    }
  }
  SEXP tol_ = PROTECT(ScalarReal(tol));
  SEXP cert_diag_ = PROTECT(native_operator_svd_certificate_cached_av(
    impl, apply, m, n, norm_A,
    combined_d_, combined_u_, combined_v_, combined_av_, tol_
  ));
  SEXP old_stage_ = VECTOR_ELT(fit_, 11);
  SEXP stage_ = PROTECT(allocVector(REALSXP, 3));
  REAL(stage_)[0] = REAL(old_stage_)[0];
  REAL(stage_)[1] = REAL(old_stage_)[1];
  REAL(stage_)[2] = stage_restart_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(stage_names_, 0, mkChar("native_iteration"));
  SET_STRING_ELT(stage_names_, 1, mkChar("ritz"));
  SET_STRING_ELT(stage_names_, 2, mkChar("restart"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 17));
  SET_VECTOR_ELT(out_, 0, combined_d_);
  SET_VECTOR_ELT(out_, 1, combined_u_);
  SET_VECTOR_ELT(out_, 2, combined_v_);
  SET_VECTOR_ELT(out_, 3, combined_av_);
  SET_VECTOR_ELT(out_, 4, combined_coeff_);
  for (int i = 5; i < 11; ++i) {
    SET_VECTOR_ELT(out_, i, VECTOR_ELT(fit_, i));
  }
  SET_VECTOR_ELT(out_, 11, stage_);
  SET_VECTOR_ELT(out_, 12, attempt_history_);
  SET_VECTOR_ELT(out_, 13, ScalarLogical(1));
  SET_VECTOR_ELT(out_, 14, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 15, cert_diag_);
  SET_VECTOR_ELT(out_, 16, ScalarInteger(locked_count));
  SEXP names_ = PROTECT(allocVector(STRSXP, 17));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("Avectors"));
  SET_STRING_ELT(names_, 4, mkChar("coefficients"));
  SET_STRING_ELT(names_, 5, mkChar("active_cols"));
  SET_STRING_ELT(names_, 6, mkChar("active_left_cols"));
  SET_STRING_ELT(names_, 7, mkChar("iterations"));
  SET_STRING_ELT(names_, 8, mkChar("matvecs"));
  SET_STRING_ELT(names_, 9, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 10, mkChar("cached_start_used"));
  SET_STRING_ELT(names_, 11, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("attempt_history"));
  SET_STRING_ELT(names_, 13, mkChar("retained_restart_native"));
  SET_STRING_ELT(names_, 14, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 15, mkChar("certificate_diagnostics"));
  SET_STRING_ELT(names_, 16, mkChar("retained_locked_count"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(locked_count > 0 ? 12 : 7);
  return out_;
}

struct BlockGolubKahanFitArrays {
  double* V = nullptr;
  double* AV = nullptr;
  double* U = nullptr;
  bool transient = false;
};

static void block_golub_kahan_fit_arrays_free(BlockGolubKahanFitArrays* arrays) {
  if (!arrays->transient) {
    std::free(arrays->V);
    std::free(arrays->AV);
    std::free(arrays->U);
  }
  arrays->V = nullptr;
  arrays->AV = nullptr;
  arrays->U = nullptr;
  arrays->transient = false;
}

static int block_golub_kahan_fit_arrays_alloc(BlockGolubKahanFitArrays* arrays,
                                              int n,
                                              int m,
                                              int max_subspace) {
  const size_t nv = static_cast<size_t>(n) * static_cast<size_t>(max_subspace);
  const size_t mv = static_cast<size_t>(m) * static_cast<size_t>(max_subspace);
  arrays->V = reinterpret_cast<double*>(R_alloc(nv > 0 ? nv : 1, sizeof(double)));
  arrays->AV = reinterpret_cast<double*>(R_alloc(mv > 0 ? mv : 1, sizeof(double)));
  arrays->U = reinterpret_cast<double*>(R_alloc(mv > 0 ? mv : 1, sizeof(double)));
  arrays->transient = true;
  return 0;
}

static int retained_subspace_sequence(int n,
                                      int rank,
                                      int block_size,
                                      int initial_max_subspace,
                                      int max_attempts,
                                      std::vector<int>* subspaces) {
  if (max_attempts < 1 || initial_max_subspace < rank || initial_max_subspace > n) {
    return -1;
  }
  subspaces->clear();
  subspaces->reserve(static_cast<size_t>(max_attempts));
  int current = initial_max_subspace;
  subspaces->push_back(current);
  for (int attempt = 1; attempt < max_attempts && current < n; ++attempt) {
    const int additive = current + ((2 * rank > 4 * block_size) ? 2 * rank : 4 * block_size);
    const int additive2 = current + 10;
    int next = additive > additive2 ? additive : additive2;
    const int scaled = static_cast<int>(ceil(1.5 * static_cast<double>(current)));
    if (scaled > next) {
      next = scaled;
    }
    if (next > n) {
      next = n;
    }
    if (next <= current) {
      break;
    }
    subspaces->push_back(next);
    current = next;
  }
  return 0;
}

static SEXP retained_attempt_history_pack(const std::vector<int>& subspaces,
                                          const std::vector<int>& active_cols,
                                          const std::vector<int>& start_cols,
                                          const std::vector<int>& iterations,
                                          const std::vector<int>& matvecs,
                                          const std::vector<int>& ortho_passes,
                                          const std::vector<int>& cached_start_used,
                                          const std::vector<int>& converged_count,
                                          const std::vector<int>& leading_converged_count,
                                          const std::vector<int>& certificate_passed,
                                          const std::vector<double>& max_backward_error,
                                          const std::vector<double>& max_residual) {
  const int rows = static_cast<int>(subspaces.size());
  SEXP out_ = PROTECT(allocVector(VECSXP, 14));
  SEXP attempt_ = PROTECT(allocVector(INTSXP, rows));
  SEXP subspace_ = PROTECT(allocVector(INTSXP, rows));
  SEXP active_ = PROTECT(allocVector(INTSXP, rows));
  SEXP start_ = PROTECT(allocVector(INTSXP, rows));
  SEXP cached_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP warm_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP iter_ = PROTECT(allocVector(INTSXP, rows));
  SEXP matvec_ = PROTECT(allocVector(INTSXP, rows));
  SEXP ortho_ = PROTECT(allocVector(INTSXP, rows));
  SEXP conv_count_ = PROTECT(allocVector(INTSXP, rows));
  SEXP leading_count_ = PROTECT(allocVector(INTSXP, rows));
  SEXP passed_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP backward_ = PROTECT(allocVector(REALSXP, rows));
  SEXP residual_ = PROTECT(allocVector(REALSXP, rows));
  for (int row = 0; row < rows; ++row) {
    INTEGER(attempt_)[row] = row + 1;
    INTEGER(subspace_)[row] = subspaces[static_cast<size_t>(row)];
    INTEGER(active_)[row] = active_cols[static_cast<size_t>(row)];
    INTEGER(start_)[row] = start_cols[static_cast<size_t>(row)];
    LOGICAL(cached_)[row] = cached_start_used[static_cast<size_t>(row)] ? TRUE : FALSE;
    LOGICAL(warm_)[row] = row > 0 ? TRUE : FALSE;
    INTEGER(iter_)[row] = iterations[static_cast<size_t>(row)];
    INTEGER(matvec_)[row] = matvecs[static_cast<size_t>(row)];
    INTEGER(ortho_)[row] = ortho_passes[static_cast<size_t>(row)];
    INTEGER(conv_count_)[row] = converged_count[static_cast<size_t>(row)];
    INTEGER(leading_count_)[row] = leading_converged_count[static_cast<size_t>(row)];
    LOGICAL(passed_)[row] = certificate_passed[static_cast<size_t>(row)] ? TRUE : FALSE;
    REAL(backward_)[row] = max_backward_error[static_cast<size_t>(row)];
    REAL(residual_)[row] = max_residual[static_cast<size_t>(row)];
  }
  SET_VECTOR_ELT(out_, 0, attempt_);
  SET_VECTOR_ELT(out_, 1, subspace_);
  SET_VECTOR_ELT(out_, 2, active_);
  SET_VECTOR_ELT(out_, 3, start_);
  SET_VECTOR_ELT(out_, 4, cached_);
  SET_VECTOR_ELT(out_, 5, warm_);
  SET_VECTOR_ELT(out_, 6, iter_);
  SET_VECTOR_ELT(out_, 7, matvec_);
  SET_VECTOR_ELT(out_, 8, ortho_);
  SET_VECTOR_ELT(out_, 9, conv_count_);
  SET_VECTOR_ELT(out_, 10, leading_count_);
  SET_VECTOR_ELT(out_, 11, passed_);
  SET_VECTOR_ELT(out_, 12, backward_);
  SET_VECTOR_ELT(out_, 13, residual_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 14));
  SET_STRING_ELT(names_, 0, mkChar("attempt"));
  SET_STRING_ELT(names_, 1, mkChar("max_subspace"));
  SET_STRING_ELT(names_, 2, mkChar("active_cols"));
  SET_STRING_ELT(names_, 3, mkChar("start_cols"));
  SET_STRING_ELT(names_, 4, mkChar("cached_start_used"));
  SET_STRING_ELT(names_, 5, mkChar("warm_started"));
  SET_STRING_ELT(names_, 6, mkChar("iterations"));
  SET_STRING_ELT(names_, 7, mkChar("matvecs"));
  SET_STRING_ELT(names_, 8, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 9, mkChar("converged_count"));
  SET_STRING_ELT(names_, 10, mkChar("leading_converged_count"));
  SET_STRING_ELT(names_, 11, mkChar("certificate_passed"));
  SET_STRING_ELT(names_, 12, mkChar("max_backward_error"));
  SET_STRING_ELT(names_, 13, mkChar("max_residual"));
  setAttrib(out_, R_NamesSymbol, names_);

  SEXP row_names_ = PROTECT(allocVector(INTSXP, rows));
  for (int row = 0; row < rows; ++row) {
    INTEGER(row_names_)[row] = row + 1;
  }
  setAttrib(out_, R_RowNamesSymbol, row_names_);
  SEXP class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(class_, 0, mkChar("data.frame"));
  setAttrib(out_, R_ClassSymbol, class_);
  UNPROTECT(18);
  return out_;
}

static int retained_cached_av_certificate_passed(void* impl,
                                                 EigencoreApplyFn apply,
                                                 int m,
                                                 int n,
                                                 double norm_A,
                                                 const double* d,
                                                 const double* u,
                                                 const double* v,
                                                 const double* av,
                                                 int k,
                                                 double tol,
                                                 double* max_backward_error,
                                                 double* max_residual,
                                                 int* converged_count,
                                                 int* leading_converged_count) {
  const double eps = DBL_EPSILON;
  const double scale_value = fmax(norm_A, eps);
  std::vector<double> right(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const int status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, k,
                           u, m, 1.0, 0.0, right.data(), n, &workspace);
  if (status != 0) {
    return status < 0 ? status : -status;
  }

  int passed = 1;
  *max_backward_error = 0.0;
  *max_residual = 0.0;
  *converged_count = 0;
  *leading_converged_count = 0;
  bool still_leading = true;
  for (int col = 0; col < k; ++col) {
    const double sigma = d[col];
    double left_sum = 0.0;
    double right_sum = 0.0;
    const int64_t left_offset = static_cast<int64_t>(col) * m;
    const int64_t right_offset = static_cast<int64_t>(col) * n;
    for (int row = 0; row < m; ++row) {
      const double residual = av[left_offset + row] -
        sigma * u[left_offset + row];
      left_sum += residual * residual;
    }
    for (int row = 0; row < n; ++row) {
      const double residual = right[right_offset + row] -
        sigma * v[right_offset + row];
      right_sum += residual * residual;
    }
    const double combined = sqrt(left_sum + right_sum);
    const double backward = combined / scale_value;
    if (backward > *max_backward_error || col == 0) {
      *max_backward_error = backward;
    }
    if (combined > *max_residual || col == 0) {
      *max_residual = combined;
    }
    const bool col_converged = R_FINITE(backward) && backward <= tol;
    if (col_converged) {
      ++(*converged_count);
      if (still_leading) {
        ++(*leading_converged_count);
      }
    } else {
      still_leading = false;
      passed = 0;
    }
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  std::vector<double> gram_u(static_cast<size_t>(k) * static_cast<size_t>(k), 0.0);
  std::vector<double> gram_v(static_cast<size_t>(k) * static_cast<size_t>(k), 0.0);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &m,
                  &one, const_cast<double*>(u), &m, const_cast<double*>(u), &m,
                  &zero, gram_u.data(), &k FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, const_cast<double*>(v), &n, const_cast<double*>(v), &n,
                  &zero, gram_v.data(), &k FCONE FCONE);
  const double orth = fmax(
    max_orthogonality_loss(gram_u.data(), k),
    max_orthogonality_loss(gram_v.data(), k)
  );
  const double orth_tol = (tol > sqrt(DBL_EPSILON)) ? tol : sqrt(DBL_EPSILON);
  if (orth > orth_tol) {
    passed = 0;
    *converged_count = 0;
    *leading_converged_count = 0;
  }
  return passed;
}

static int normalize_retained_cached_prefix(double* V,
                                            int n,
                                            double* AV,
                                            int m,
                                            int cols) {
  if (cols <= 0) {
    return 0;
  }
  std::vector<double> gram(static_cast<size_t>(cols) * static_cast<size_t>(cols), 0.0);
  const char trans = 'T';
  const char notrans = 'N';
  const char right = 'R';
  const char uplo = 'U';
  const char diag = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  F77_CALL(dgemm)(&trans, &notrans, &cols, &cols, &n,
                  &one, V, &n, V, &n,
                  &zero, gram.data(), &cols FCONE FCONE);
  symmetrize_packed_square(gram.data(), cols);
  F77_CALL(dpotrf)(&uplo, &cols, gram.data(), &cols, &info FCONE);
  if (info != 0) {
    return -1;
  }
  for (int col = 0; col < cols; ++col) {
    if (gram[col + static_cast<int64_t>(col) * cols] <= 100.0 * DBL_EPSILON) {
      return -2;
    }
  }

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &n, &cols, &one,
                  gram.data(), &cols, V, &n FCONE FCONE FCONE FCONE);
  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &m, &cols, &one,
                  gram.data(), &cols, AV, &m FCONE FCONE FCONE FCONE);
  return 0;
}

template <typename ConfigureOperator>
static SEXP block_golub_kahan_retained_cycle_impl(ConfigureOperator configure_operator,
                                                  int m,
                                                  int n,
                                                  int initial_max_subspace,
                                                  int block_size,
                                                  const double* initial_start,
                                                  const double* random_tails,
                                                  int random_tail_cols,
                                                  int max_attempts,
                                                  int rank,
                                                  int target_kind,
                                                  double norm_A,
                                                  double tol,
                                                  int use_retained_av_cache,
                                                  int use_deflation) {
  if (block_size < 1 || rank < 1 || initial_max_subspace < rank ||
      initial_max_subspace > n || max_attempts < 1) {
    error("invalid retained block Golub-Kahan cycle controls");
  }
  if (random_tail_cols < block_size * (max_attempts - 1)) {
    error("random_tails must have at least block * (max_attempts - 1) columns");
  }

  std::vector<int> subspaces;
  if (retained_subspace_sequence(
        n, rank, block_size, initial_max_subspace, max_attempts, &subspaces
      ) != 0 || subspaces.empty()) {
    error("failed to construct retained block Golub-Kahan subspace sequence");
  }
  const int final_max_subspace = subspaces.back();
  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, final_max_subspace) != 0) {
    error("failed to allocate native retained block Golub-Kahan workspace");
  }
  BlockGolubKahanBasisScratch basis_scratch;
  const int retained_start_capacity = rank + block_size;
  if (block_golub_kahan_basis_scratch_alloc(
        &basis_scratch, m, n, final_max_subspace, retained_start_capacity
      ) != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("failed to allocate native retained block Golub-Kahan basis scratch");
  }
  const double native_workspace_bytes = static_cast<double>(
    (static_cast<size_t>(n) * static_cast<size_t>(final_max_subspace) +
     2 * static_cast<size_t>(m) * static_cast<size_t>(final_max_subspace)) *
      sizeof(double) + basis_scratch.bytes
  );

  std::vector<double> current_start(static_cast<size_t>(n) *
                                    static_cast<size_t>(rank + block_size), 0.0);
  std::vector<double> current_start_av(static_cast<size_t>(m) *
                                       static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_d(static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_u(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_v(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_av(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  std::memcpy(
    current_start.data(), initial_start,
    sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(block_size)
  );
  int current_start_cols = block_size;
  int current_start_av_cols = 0;

  void* impl = nullptr;
  EigencoreApplyFn apply = nullptr;
  configure_operator(&impl, &apply);

  std::vector<int> history_subspaces;
  std::vector<int> history_active_cols;
  std::vector<int> history_start_cols;
  std::vector<int> history_iterations;
  std::vector<int> history_matvecs;
  std::vector<int> history_ortho_passes;
  std::vector<int> history_cached_start_used;
  std::vector<int> history_converged_count;
  std::vector<int> history_leading_converged_count;
  std::vector<int> history_certificate_passed;
  std::vector<double> history_max_backward_error;
  std::vector<double> history_max_residual;

  int final_active_v = 0;
  int final_active_u = 0;
  int final_iterations = 0;
  int final_matvecs = 0;
  int final_ortho_passes = 0;
  int final_cached_start_used = 0;
  int locked_count = 0;
  double total_stage_native_iteration_seconds = 0.0;
  double total_stage_restart_seconds = 0.0;

  for (size_t attempt = 0; attempt < subspaces.size(); ++attempt) {
    const int active_rank = rank - locked_count;
    if (active_rank <= 0) {
      break;
    }
    const int max_subspace = subspaces[attempt];
    int active_max_subspace = max_subspace;
    if (use_deflation && locked_count > 0) {
      active_max_subspace -= locked_count;
      if (active_max_subspace < active_rank) {
        active_max_subspace = active_rank;
      }
      if (active_max_subspace > n - locked_count) {
        active_max_subspace = n - locked_count;
      }
    }
    int active_v = 0;
    int active_u = 0;
    int iterations = 0;
    int matvecs = 0;
    int ortho_passes = 0;
    int cached_start_used = 0;

    auto stage_timer = native_timer_now();
    const int status = native_block_golub_kahan_basis_run_with_scratch(
      impl, apply, m, n, active_max_subspace, current_start_cols,
      current_start.data(),
      (use_retained_av_cache && current_start_av_cols > 0) ? current_start_av.data() : nullptr,
      use_retained_av_cache ? current_start_av_cols : 0,
      (use_deflation && locked_count > 0) ? locked_v.data() : nullptr,
      use_deflation ? locked_count : 0,
      (use_deflation && locked_count > 0) ? locked_u.data() : nullptr,
      use_deflation ? locked_count : 0,
      arrays.V, arrays.AV, arrays.U, &basis_scratch,
      &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
      &cached_start_used
    );
    total_stage_native_iteration_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      block_golub_kahan_basis_scratch_free(&basis_scratch);
      block_golub_kahan_fit_arrays_free(&arrays);
      error("native retained block Golub-Kahan cycle failed with status=%d", status);
    }

    history_subspaces.push_back(max_subspace);
    history_active_cols.push_back(active_v);
    history_start_cols.push_back(current_start_cols);
    history_iterations.push_back(iterations);
    history_matvecs.push_back(matvecs);
    history_ortho_passes.push_back(ortho_passes);
    history_cached_start_used.push_back(cached_start_used);

    final_active_v = active_v;
    final_active_u = active_u;
    final_iterations += iterations;
    final_matvecs += matvecs;
    final_ortho_passes += ortho_passes;
    final_cached_start_used = final_cached_start_used || cached_start_used;

    SEXP ritz_ = PROTECT(eigencore_block_golub_kahan_ritz_from_ptr(
      arrays.V, n, arrays.AV, m, active_v, active_rank, target_kind
    ));
    double max_backward_error = R_PosInf;
    double max_residual = R_PosInf;
    int converged_count = 0;
    int leading_converged_count = 0;
    const int certificate_passed = retained_cached_av_certificate_passed(
      impl, apply, m, n, norm_A,
      REAL(VECTOR_ELT(ritz_, 0)),
      REAL(VECTOR_ELT(ritz_, 1)),
      REAL(VECTOR_ELT(ritz_, 2)),
      REAL(VECTOR_ELT(ritz_, 3)),
      LENGTH(VECTOR_ELT(ritz_, 0)),
      tol,
      &max_backward_error,
      &max_residual,
      &converged_count,
      &leading_converged_count
    );
    if (certificate_passed < 0) {
      UNPROTECT(1);
      block_golub_kahan_basis_scratch_free(&basis_scratch);
      block_golub_kahan_fit_arrays_free(&arrays);
      error("native retained block Golub-Kahan attempt certificate failed with status=%d",
            certificate_passed);
    }
    history_converged_count.push_back(locked_count + converged_count);
    history_leading_converged_count.push_back(locked_count + leading_converged_count);
    history_certificate_passed.push_back(certificate_passed);
    history_max_backward_error.push_back(max_backward_error);
    history_max_residual.push_back(max_residual);

    if (certificate_passed) {
      UNPROTECT(1);
      break;
    }

    if (attempt + 1 >= subspaces.size()) {
      UNPROTECT(1);
      break;
    }

    auto restart_timer = native_timer_now();
    SEXP v_ = VECTOR_ELT(ritz_, 2);
    SEXP av_ = VECTOR_ELT(ritz_, 3);
    SEXP u_ = VECTOR_ELT(ritz_, 1);
    SEXP d_ = VECTOR_ELT(ritz_, 0);
    const int ritz_cols = INTEGER(getAttrib(v_, R_DimSymbol))[1];
    int newly_locked = 0;
    if (use_deflation && leading_converged_count > 0) {
      newly_locked = leading_converged_count;
      if (newly_locked > rank - locked_count) {
        newly_locked = rank - locked_count;
      }
      for (int col = 0; col < newly_locked; ++col) {
        const int dst_col = locked_count + col;
        locked_d[static_cast<size_t>(dst_col)] = REAL(d_)[col];
        std::memcpy(locked_u.data() + static_cast<int64_t>(dst_col) * m,
                    REAL(u_) + static_cast<int64_t>(col) * m,
                    sizeof(double) * static_cast<size_t>(m));
        std::memcpy(locked_v.data() + static_cast<int64_t>(dst_col) * n,
                    REAL(v_) + static_cast<int64_t>(col) * n,
                    sizeof(double) * static_cast<size_t>(n));
        std::memcpy(locked_av.data() + static_cast<int64_t>(dst_col) * m,
                    REAL(av_) + static_cast<int64_t>(col) * m,
                    sizeof(double) * static_cast<size_t>(m));
      }
      locked_count += newly_locked;
    }
    const int keep_cols = ritz_cols - newly_locked;
    current_start.assign(current_start.size(), 0.0);
    current_start_av.assign(current_start_av.size(), 0.0);
    for (int col = 0; col < keep_cols; ++col) {
      std::memcpy(
        current_start.data() + static_cast<int64_t>(col) * n,
        REAL(v_) + static_cast<int64_t>(newly_locked + col) * n,
        sizeof(double) * static_cast<size_t>(n)
      );
      if (use_retained_av_cache) {
        std::memcpy(
          current_start_av.data() + static_cast<int64_t>(col) * m,
          REAL(av_) + static_cast<int64_t>(newly_locked + col) * m,
          sizeof(double) * static_cast<size_t>(m)
        );
      }
    }
    if (use_retained_av_cache) {
      const int normalize_status = normalize_retained_cached_prefix(
        current_start.data(), n, current_start_av.data(), m, keep_cols
      );
      if (normalize_status != 0) {
        UNPROTECT(1);
        block_golub_kahan_basis_scratch_free(&basis_scratch);
        block_golub_kahan_fit_arrays_free(&arrays);
        error("native retained block Golub-Kahan cached prefix normalization failed with status=%d",
              normalize_status);
      }
    }
    const double* tail = random_tails +
      static_cast<int64_t>(attempt) * static_cast<int64_t>(n) *
        static_cast<int64_t>(block_size);
    for (int col = 0; col < block_size; ++col) {
      std::memcpy(
        current_start.data() + static_cast<int64_t>(keep_cols + col) * n,
        tail + static_cast<int64_t>(col) * n,
        sizeof(double) * static_cast<size_t>(n)
      );
    }
    current_start_cols = keep_cols + block_size;
    current_start_av_cols = use_retained_av_cache ? keep_cols : 0;
    UNPROTECT(1);
    total_stage_restart_seconds += native_timer_elapsed(restart_timer);
  }

  SEXP history_ = PROTECT(retained_attempt_history_pack(
    history_subspaces, history_active_cols, history_start_cols,
    history_iterations, history_matvecs, history_ortho_passes,
    history_cached_start_used,
    history_converged_count, history_leading_converged_count,
    history_certificate_passed,
    history_max_backward_error, history_max_residual
  ));
  const int final_active_rank = rank - locked_count;
  SEXP out_ = PROTECT(block_golub_kahan_retained_cycle_pack(
    n, m, impl, apply, norm_A, tol,
    arrays.V, arrays.AV, final_active_v, final_active_u,
    final_iterations, final_matvecs, final_ortho_passes,
    final_cached_start_used, rank, final_active_rank, target_kind,
    total_stage_native_iteration_seconds, total_stage_restart_seconds,
    history_, native_workspace_bytes,
    locked_d.data(), locked_u.data(), locked_v.data(), locked_av.data(),
    locked_count
  ));
  block_golub_kahan_basis_scratch_free(&basis_scratch);
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(2);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_fit(SEXP A_,
                                                      SEXP max_subspace_,
                                                      SEXP start_,
                                                      SEXP rank_,
                                                      SEXP target_kind_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native dense block Golub-Kahan fit workspace");
  }
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native dense block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_fit_cached(SEXP A_,
                                                             SEXP max_subspace_,
                                                             SEXP start_,
                                                             SEXP rank_,
                                                             SEXP target_kind_,
                                                             SEXP start_av_) {
  if (!isReal(A_) || !isReal(start_) || !isReal(start_av_)) {
    error("A, start, and start_av must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || dimAV == R_NilValue) {
    error("A, start, and start_av must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  const int start_av_cols = INTEGER(dimAV)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native dense cached block Golub-Kahan fit workspace");
  }
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native dense cached block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_retained_cycle(SEXP A_,
                                                                 SEXP initial_max_subspace_,
                                                                 SEXP initial_start_,
                                                                 SEXP random_tails_,
                                                                 SEXP max_attempts_,
                                                                 SEXP rank_,
                                                                 SEXP target_kind_,
                                                                 SEXP norm_A_,
                                                                 SEXP tol_,
                                                                 SEXP use_retained_av_cache_,
                                                                 SEXP use_deflation_) {
  if (!isReal(A_) || !isReal(initial_start_) || !isReal(random_tails_)) {
    error("A, initial_start, and random_tails must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(initial_start_, R_DimSymbol);
  SEXP dimT = getAttrib(random_tails_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || dimT == R_NilValue) {
    error("A, initial_start, and random_tails must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  if (INTEGER(dimS)[0] != n || block_size < 1 || INTEGER(dimT)[0] != n) {
    error("non-conformable dense retained block Golub-Kahan inputs");
  }
  DenseColumnMajorOperator impl_holder = {m, n, REAL(A_)};
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_dense_apply;
  };
  return block_golub_kahan_retained_cycle_impl(
    configure, m, n, asInteger(initial_max_subspace_), block_size,
    REAL(initial_start_), REAL(random_tails_), INTEGER(dimT)[1],
    asInteger(max_attempts_), asInteger(rank_), asInteger(target_kind_),
    asReal(norm_A_), asReal(tol_),
    asLogical(use_retained_av_cache_) == TRUE,
    asLogical(use_deflation_) == TRUE
  );
}

extern "C" SEXP eigencore_block_golub_kahan_csc_fit(SEXP i_, SEXP p_,
                                                    SEXP x_, SEXP dim_,
                                                    SEXP max_subspace_,
                                                    SEXP start_,
                                                    SEXP rank_,
                                                    SEXP target_kind_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_)) {
    error("invalid CSC block Golub-Kahan fit inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(dim_) != 2) {
    error("start must be a matrix and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native CSC block Golub-Kahan fit workspace");
  }
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native CSC block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_csc_fit_cached(SEXP i_, SEXP p_,
                                                           SEXP x_, SEXP dim_,
                                                           SEXP max_subspace_,
                                                           SEXP start_,
                                                           SEXP rank_,
                                                           SEXP target_kind_,
                                                           SEXP start_av_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_) || !isReal(start_av_)) {
    error("invalid cached CSC block Golub-Kahan fit inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimS == R_NilValue || dimAV == R_NilValue || LENGTH(dim_) != 2) {
    error("start and start_av must be matrices and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  const int start_av_cols = INTEGER(dimAV)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native CSC cached block Golub-Kahan fit workspace");
  }
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native CSC cached block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_csc_retained_cycle(SEXP i_, SEXP p_,
                                                               SEXP x_, SEXP dim_,
                                                               SEXP initial_max_subspace_,
                                                               SEXP initial_start_,
                                                               SEXP random_tails_,
                                                               SEXP max_attempts_,
                                                               SEXP rank_,
                                                               SEXP target_kind_,
                                                               SEXP norm_A_,
                                                               SEXP tol_,
                                                               SEXP use_retained_av_cache_,
                                                               SEXP use_deflation_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(initial_start_) || !isReal(random_tails_)) {
    error("invalid CSC retained block Golub-Kahan inputs");
  }
  SEXP dimS = getAttrib(initial_start_, R_DimSymbol);
  SEXP dimT = getAttrib(random_tails_, R_DimSymbol);
  if (dimS == R_NilValue || dimT == R_NilValue || LENGTH(dim_) != 2) {
    error("initial_start and random_tails must be matrices and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  if (INTEGER(dimS)[0] != n || block_size < 1 || INTEGER(dimT)[0] != n) {
    error("non-conformable CSC retained block Golub-Kahan inputs");
  }
  CSCOperator impl_holder = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_csc_apply;
  };
  return block_golub_kahan_retained_cycle_impl(
    configure, m, n, asInteger(initial_max_subspace_), block_size,
    REAL(initial_start_), REAL(random_tails_), INTEGER(dimT)[1],
    asInteger(max_attempts_), asInteger(rank_), asInteger(target_kind_),
    asReal(norm_A_), asReal(tol_),
    asLogical(use_retained_av_cache_) == TRUE,
    asLogical(use_deflation_) == TRUE
  );
}


extern "C" SEXP eigencore_golub_kahan_dense_fit(SEXP A_, SEXP maxit_, SEXP start_,
                                                SEXP rank_, SEXP target_kind_,
                                                SEXP tol_, SEXP projected_stop_,
                                                SEXP reorth_u_, SEXP reorth_v_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  if (LENGTH(start_) != n) {
    error("non-conformable dense Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  const int reorthogonalize_u = asLogical(reorth_u_) == TRUE;
  const int reorthogonalize_v = asLogical(reorth_v_) == TRUE;
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }
  if (rank < 1 || rank > maxit) {
    error("rank must be between 1 and maxit");
  }

  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(maxit), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(maxit), 0.0);
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(maxit) +
      static_cast<int64_t>(2 * m + 2 * n) +
      static_cast<int64_t>(maxit)
    );

  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int iterations = 0;
  int matvecs = 0;
  int projected_stop = 0;
  int projected_nconv = 0;
  double projected_max_residual = R_PosInf;
  int projected_checks = 0;
  double projected_seconds = 0.0;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  int reorthogonalization_passes = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_dense_apply, m, n,
                                            maxit, rank, target_kind, tol,
                                            enable_projected_stop, 1,
                                            REAL(start_),
                                            U_work.data(), V_work.data(),
                                            alpha_work.data(), beta_work.data(),
                                            &iterations, &matvecs,
                                            &projected_stop, &projected_nconv,
                                            &projected_max_residual,
                                            &projected_checks,
                                            &projected_seconds,
                                            &stage_apply_seconds,
                                            &stage_recurrence_seconds,
                                            &stage_reorthogonalization_seconds,
                                            &reorthogonalization_passes,
                                            reorthogonalize_u,
                                            reorthogonalize_v);
  if (status != 0) {
    error("native dense Golub-Kahan failed with status=%d", status);
  }

  SEXP ritz_ = PROTECT(eigencore_golub_kahan_ritz_from_ptr(
    U_work.data(), V_work.data(), m, n, iterations, alpha_work.data(),
    beta_work.data(), rank, target_kind
  ));
  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 7, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 9, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 11, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 15, ScalarLogical(FALSE));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 6, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 7, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 8, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 9, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 11, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 15, mkChar("basis_returned"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(3);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_csc_fit(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                              SEXP maxit_, SEXP start_,
                                              SEXP rank_, SEXP target_kind_,
                                              SEXP tol_, SEXP projected_stop_,
                                              SEXP reorth_u_, SEXP reorth_v_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC Golub-Kahan inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  if (LENGTH(start_) != n) {
    error("non-conformable CSC Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  const int reorthogonalize_u = asLogical(reorth_u_) == TRUE;
  const int reorthogonalize_v = asLogical(reorth_v_) == TRUE;
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }
  if (rank < 1 || rank > maxit) {
    error("rank must be between 1 and maxit");
  }

  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(maxit), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(maxit), 0.0);
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(maxit) +
      static_cast<int64_t>(2 * m + 2 * n)
    );

  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int iterations = 0;
  int matvecs = 0;
  int projected_stop = 0;
  int projected_nconv = 0;
  double projected_max_residual = R_PosInf;
  int projected_checks = 0;
  double projected_seconds = 0.0;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  int reorthogonalization_passes = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_csc_apply, m, n,
                                            maxit, rank, target_kind, tol,
                                            enable_projected_stop, 0,
                                            REAL(start_),
                                            U_work.data(), V_work.data(),
                                            alpha_work.data(), beta_work.data(),
                                            &iterations, &matvecs,
                                            &projected_stop, &projected_nconv,
                                            &projected_max_residual,
                                            &projected_checks,
                                            &projected_seconds,
                                            &stage_apply_seconds,
                                            &stage_recurrence_seconds,
                                            &stage_reorthogonalization_seconds,
                                            &reorthogonalization_passes,
                                            reorthogonalize_u,
                                            reorthogonalize_v);
  if (status != 0) {
    error("native CSC Golub-Kahan failed with status=%d", status);
  }

  SEXP ritz_ = PROTECT(eigencore_golub_kahan_ritz_from_ptr(
    U_work.data(), V_work.data(), m, n, iterations, alpha_work.data(),
    beta_work.data(), rank, target_kind
  ));
  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 7, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 9, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 11, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 15, ScalarLogical(FALSE));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 6, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 7, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 8, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 9, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 11, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 15, mkChar("basis_returned"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(3);
  return out_;
}

static void validate_real_vector_length(SEXP x, int n, const char* name) {
  if (!isReal(x) || LENGTH(x) != n) {
    error("%s must be a double vector of length %d", name, n);
  }
}

static void validate_real_matrix_dim(SEXP x, int nrow, int ncol, const char* name) {
  if (!isReal(x)) {
    error("%s must be a double matrix", name);
  }
  SEXP dim = getAttrib(x, R_DimSymbol);
  if (dim == R_NilValue || LENGTH(dim) != 2 ||
      INTEGER(dim)[0] != nrow || INTEGER(dim)[1] != ncol) {
    error("%s must have dimensions %d x %d", name, nrow, ncol);
  }
}

static void validate_irlba_lbd_retained_contract(int m, int n,
                                                 SEXP initial_start_,
                                                 SEXP retained_right_,
                                                 SEXP retained_left_,
                                                 SEXP alpha_, SEXP beta_,
                                                 SEXP random_tails_,
                                                 SEXP work_, SEXP retained_,
                                                 SEXP max_restarts_,
                                                 SEXP rank_,
                                                 SEXP target_kind_,
                                                 SEXP tol_,
                                                 SEXP reorth_policy_) {
  const int work = asInteger(work_);
  const int retained = asInteger(retained_);
  const int max_restarts = asInteger(max_restarts_);
  const int rank = asInteger(rank_);
  const int target_kind = asInteger(target_kind_);
  const double tol = asReal(tol_);
  const int reorth_policy = asInteger(reorth_policy_);
  if (m < 1 || n < 1) {
    error("active operator dimensions must be positive");
  }
  if (rank < 1 || rank > retained) {
    error("rank must be between 1 and retained");
  }
  if (retained < 1 || retained >= work) {
    error("retained must be positive and smaller than work");
  }
  if (work > ((m < n) ? m : n)) {
    error("work must not exceed min(active dimensions)");
  }
  if (max_restarts < 0) {
    error("max_restarts must be non-negative");
  }
  if (target_kind != 1) {
    error("retained IRLBA/LBD currently supports only largest singular values");
  }
  if (!R_FINITE(tol) || tol <= 0.0) {
    error("tol must be a positive finite scalar");
  }
  if (reorth_policy < 1 || reorth_policy > 3) {
    error("reorth_policy must identify a known retained IRLBA/LBD policy");
  }
  validate_real_vector_length(initial_start_, n, "initial_start");
  validate_real_matrix_dim(retained_right_, n, retained, "retained_right");
  validate_real_matrix_dim(retained_left_, m, retained, "retained_left");
  validate_real_vector_length(alpha_, work, "alpha");
  validate_real_vector_length(beta_, work, "beta");
  validate_real_matrix_dim(random_tails_, n, work - retained, "random_tails");
}

static void irlba_lbd_retained_seed(int n,
                                    int retained,
                                    const double* initial_start,
                                    const double* retained_right,
                                    const double* random_tails,
                                    int random_tail_cols,
                                    int attempt,
                                    double* start) {
  std::memset(start, 0, sizeof(double) * static_cast<size_t>(n));
  const int retained_col = (attempt < retained) ? attempt : 0;
  const double* base = retained_right +
    static_cast<int64_t>(retained_col) * static_cast<int64_t>(n);
  long double norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    start[row] = base[row];
    norm2 += static_cast<long double>(start[row]) * start[row];
  }
  if (norm2 <= 100.0L * DBL_EPSILON) {
    norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      start[row] = initial_start[row];
      norm2 += static_cast<long double>(start[row]) * start[row];
    }
  }
  if (attempt > 0 && random_tail_cols > 0) {
    const int tail_col = (attempt - 1) % random_tail_cols;
    const double* tail = random_tails +
      static_cast<int64_t>(tail_col) * static_cast<int64_t>(n);
    for (int row = 0; row < n; ++row) {
      start[row] += 1.0e-3 * tail[row];
    }
  }
  norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    norm2 += static_cast<long double>(start[row]) * start[row];
  }
  if (norm2 <= 100.0L * DBL_EPSILON) {
    start[0] = 1.0;
    norm2 = 1.0L;
  }
  const double norm = sqrt(static_cast<double>(norm2));
  for (int row = 0; row < n; ++row) {
    start[row] /= norm;
  }
}

static SEXP irlba_lbd_attempt_history_pack(const std::vector<int>& attempts,
                                           const std::vector<int>& iterations,
                                           const std::vector<int>& matvecs,
                                           const std::vector<int>& warm_started,
                                           const std::vector<double>* cheap_residuals = nullptr) {
  const int rows = static_cast<int>(attempts.size());
  const int cols = (cheap_residuals != nullptr) ? 6 : 5;
  SEXP attempt_ = PROTECT(allocVector(INTSXP, rows));
  SEXP max_subspace_ = PROTECT(allocVector(INTSXP, rows));
  SEXP iterations_ = PROTECT(allocVector(INTSXP, rows));
  SEXP matvecs_ = PROTECT(allocVector(INTSXP, rows));
  SEXP warm_started_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP cheap_residual_ = R_NilValue;
  if (cheap_residuals != nullptr) {
    cheap_residual_ = PROTECT(allocVector(REALSXP, rows));
  }
  for (int row = 0; row < rows; ++row) {
    INTEGER(attempt_)[row] = row + 1;
    INTEGER(max_subspace_)[row] = attempts[static_cast<size_t>(row)];
    INTEGER(iterations_)[row] = iterations[static_cast<size_t>(row)];
    INTEGER(matvecs_)[row] = matvecs[static_cast<size_t>(row)];
    LOGICAL(warm_started_)[row] = warm_started[static_cast<size_t>(row)] ? TRUE : FALSE;
    if (cheap_residuals != nullptr) {
      REAL(cheap_residual_)[row] = (*cheap_residuals)[static_cast<size_t>(row)];
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, cols));
  SET_VECTOR_ELT(out_, 0, attempt_);
  SET_VECTOR_ELT(out_, 1, max_subspace_);
  SET_VECTOR_ELT(out_, 2, iterations_);
  SET_VECTOR_ELT(out_, 3, matvecs_);
  SET_VECTOR_ELT(out_, 4, warm_started_);
  if (cheap_residuals != nullptr) {
    SET_VECTOR_ELT(out_, 5, cheap_residual_);
  }
  SEXP names_ = PROTECT(allocVector(STRSXP, cols));
  SET_STRING_ELT(names_, 0, mkChar("attempt"));
  SET_STRING_ELT(names_, 1, mkChar("max_subspace"));
  SET_STRING_ELT(names_, 2, mkChar("iterations"));
  SET_STRING_ELT(names_, 3, mkChar("matvecs"));
  SET_STRING_ELT(names_, 4, mkChar("warm_started"));
  if (cheap_residuals != nullptr) {
    SET_STRING_ELT(names_, 5, mkChar("cheap_residual"));
  }
  setAttrib(out_, R_NamesSymbol, names_);
  SEXP row_names_ = PROTECT(allocVector(INTSXP, 2));
  INTEGER(row_names_)[0] = NA_INTEGER;
  INTEGER(row_names_)[1] = -rows;
  setAttrib(out_, R_RowNamesSymbol, row_names_);
  SEXP class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(class_, 0, mkChar("data.frame"));
  setAttrib(out_, R_ClassSymbol, class_);
  UNPROTECT(cheap_residuals != nullptr ? 10 : 9);
  return out_;
}

struct BproAppendDiagnostics {
  double threshold;
  int monitored_appends;
  int threshold_reorthogonalizations;
  double max_estimated_orthogonality_loss;
  double max_post_append_orthogonality_loss;
  int escalation_recommended;
};

static double vector_norm2_sqrt(const double* x, int n) {
  long double norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    norm2 += static_cast<long double>(x[row]) * x[row];
  }
  return sqrt(static_cast<double>(norm2));
}

static double candidate_basis_correlation_loss(const double* basis,
                                               int n,
                                               int cols,
                                               const double* candidate,
                                               double candidate_norm) {
  if (cols <= 0) {
    return 0.0;
  }
  if (!R_FINITE(candidate_norm) || candidate_norm <= 0.0) {
    return R_PosInf;
  }
  double loss = 0.0;
  for (int col = 0; col < cols; ++col) {
    const double* q = basis + static_cast<int64_t>(col) * n;
    double dot = 0.0;
    for (int row = 0; row < n; ++row) {
      dot += q[row] * candidate[row];
    }
    const double corr = fabs(dot) / candidate_norm;
    if (corr > loss) {
      loss = corr;
    }
  }
  return loss;
}

static void mgs_orthogonalization_pass(double* basis,
                                       int n,
                                       int cols,
                                       double* z) {
  for (int col = 0; col < cols; ++col) {
    const double* q = basis + static_cast<int64_t>(col) * n;
    double dot = 0.0;
    for (int row = 0; row < n; ++row) {
      dot += q[row] * z[row];
    }
    for (int row = 0; row < n; ++row) {
      z[row] -= dot * q[row];
    }
  }
}

static int append_orthonormal_column(double* basis,
                                     int n,
                                     int capacity,
                                     int* cols,
                                     const double* candidate,
                                     double tol,
                                     int* orthogonalization_passes,
                                     double* accepted_norm,
                                     int requested_passes,
                                     BproAppendDiagnostics* bpro) {
  if (accepted_norm != nullptr) {
    *accepted_norm = 0.0;
  }
  if (*cols >= capacity) {
    return 0;
  }
  std::vector<double> z(static_cast<size_t>(n), 0.0);
  std::memcpy(z.data(), candidate, sizeof(double) * static_cast<size_t>(n));
  if (requested_passes < 1) {
    requested_passes = 1;
  }
  const int passes = (*cols > 0) ? requested_passes : 1;
  for (int pass = 0; pass < passes; ++pass) {
    mgs_orthogonalization_pass(basis, n, *cols, z.data());
    if (orthogonalization_passes != nullptr && *cols > 0) {
      ++(*orthogonalization_passes);
    }
  }
  if (bpro != nullptr && *cols > 0) {
    ++bpro->monitored_appends;
    double monitored_norm = vector_norm2_sqrt(z.data(), n);
    double estimated_loss = candidate_basis_correlation_loss(
      basis, n, *cols, z.data(), monitored_norm
    );
    if (estimated_loss > bpro->max_estimated_orthogonality_loss) {
      bpro->max_estimated_orthogonality_loss = estimated_loss;
    }
    if (!R_FINITE(estimated_loss) || estimated_loss > bpro->threshold) {
      mgs_orthogonalization_pass(basis, n, *cols, z.data());
      if (orthogonalization_passes != nullptr) {
        ++(*orthogonalization_passes);
      }
      ++bpro->threshold_reorthogonalizations;
      monitored_norm = vector_norm2_sqrt(z.data(), n);
      estimated_loss = candidate_basis_correlation_loss(
        basis, n, *cols, z.data(), monitored_norm
      );
    }
    if (estimated_loss > bpro->max_post_append_orthogonality_loss) {
      bpro->max_post_append_orthogonality_loss = estimated_loss;
    }
    if (!R_FINITE(estimated_loss) || estimated_loss > bpro->threshold) {
      bpro->escalation_recommended = 1;
    }
  }
  const double norm = vector_norm2_sqrt(z.data(), n);
  const double threshold = fmax(100.0 * DBL_EPSILON, tol * 1.0e-4);
  if (!R_FINITE(norm) || norm <= threshold) {
    return 0;
  }
  if (accepted_norm != nullptr) {
    *accepted_norm = norm;
  }
  double* dst = basis + static_cast<int64_t>(*cols) * n;
  for (int row = 0; row < n; ++row) {
    dst[row] = z[row] / norm;
  }
  ++(*cols);
  return 1;
}

static int append_orthonormal_block(double* basis,
                                    int n,
                                    int capacity,
                                    int* cols,
                                    const double* candidates,
                                    int candidate_cols,
                                    double tol,
                                    int* orthogonalization_passes,
                                    int requested_passes = 2,
                                    BproAppendDiagnostics* bpro = nullptr) {
  int accepted = 0;
  for (int col = 0; col < candidate_cols && *cols < capacity; ++col) {
    accepted += append_orthonormal_column(
      basis, n, capacity, cols,
      candidates + static_cast<int64_t>(col) * n,
      tol, orthogonalization_passes, nullptr, requested_passes, bpro
    );
  }
  return accepted;
}

static double basis_orthogonality_loss(const double* basis, int n, int cols) {
  if (cols <= 0) {
    return 0.0;
  }
  std::vector<double> gram(static_cast<size_t>(cols) * static_cast<size_t>(cols), 0.0);
  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans, &notrans, &cols, &cols, &n,
                  &one, const_cast<double*>(basis), &n,
                  const_cast<double*>(basis), &n,
                  &zero, gram.data(), &cols FCONE FCONE);
  return max_orthogonality_loss(gram.data(), cols);
}

static int apply_augmented_basis_columns(void* impl,
                                         EigencoreApplyFn apply,
                                         int m,
                                         int n,
                                         const double* Q,
                                         int from_col,
                                         int to_col,
                                         double* AQ,
                                         double* stage_apply_seconds,
                                         int* matvecs,
                                         EigencoreWorkspace* workspace) {
  if (from_col >= to_col) {
    return 0;
  }
  const int cols = to_col - from_col;
  auto stage_timer = native_timer_now();
  const int status = apply(
    impl, EIGENCORE_TRANSPOSE_NONE, cols,
    Q + static_cast<int64_t>(from_col) * n, n,
    1.0, 0.0,
    AQ + static_cast<int64_t>(from_col) * m, m,
    workspace
  );
  *stage_apply_seconds += native_timer_elapsed(stage_timer);
  if (status == 0 && matvecs != nullptr) {
    *matvecs += cols;
  }
  return status;
}

static SEXP irlba_lbd_augmented_retained_projection(
    void* impl,
    EigencoreApplyFn apply,
    int m,
    int n,
    const double* initial_start,
    const double* retained_right,
    const double* retained_left,
    const double* random_tails,
    int work,
    int retained,
    int max_restarts,
    int rank,
    int target_kind,
    double tol,
    int reorthogonalize_u,
    int reorthogonalize_v,
    double native_workspace_bytes,
    int bpro_policy) {
  const int tail_width = work - retained;
  const int retained_core = (rank < retained) ? rank : retained;
  const int requested_tail_steps = (tail_width > 0)
    ? (max_restarts + 1) * tail_width
    : 0;
  const int capacity = (n < retained_core + retained_core + requested_tail_steps + 1)
    ? n
    : retained_core + retained_core + requested_tail_steps + 1;
  if (capacity < rank) {
    return R_NilValue;
  }

  std::vector<double> Q(static_cast<size_t>(n) * static_cast<size_t>(capacity), 0.0);
  std::vector<double> AV_ret(static_cast<size_t>(m) * static_cast<size_t>(retained_core), 0.0);
  std::vector<double> ATU_ret(static_cast<size_t>(n) * static_cast<size_t>(retained_core), 0.0);
  std::vector<double> H(static_cast<size_t>(retained_core) * static_cast<size_t>(retained_core), 0.0);
  std::vector<double> residual(static_cast<size_t>(n) * static_cast<size_t>(retained_core), 0.0);
  std::vector<double> AQ(static_cast<size_t>(m) * static_cast<size_t>(capacity), 0.0);
  std::vector<double> u(static_cast<size_t>(m), 0.0);
  std::vector<double> u_prev(static_cast<size_t>(m), 0.0);
  std::vector<double> z(static_cast<size_t>(n), 0.0);

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int q_cols = 0;
  int aq_cols = 0;
  int orthogonalization_passes = 0;
  const int requested_orthogonalization_passes = bpro_policy ? 1 : 2;
  const double bpro_threshold = fmax(tol, sqrt(DBL_EPSILON));
  BproAppendDiagnostics bpro = {
    bpro_threshold,
    0,
    0,
    0.0,
    0.0,
    0
  };
  BproAppendDiagnostics* bpro_ptr = bpro_policy ? &bpro : nullptr;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  double stage_projected_seconds = 0.0;
  int matvecs = 0;
  std::vector<double> tail_beta_history;
  tail_beta_history.reserve(static_cast<size_t>(requested_tail_steps));

  auto stage_timer = native_timer_now();
  int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, retained_core,
                     retained_right, n, 1.0, 0.0, AV_ret.data(), m, &workspace);
  stage_apply_seconds += native_timer_elapsed(stage_timer);
  if (status != 0) {
    return R_NilValue;
  }
  matvecs += retained_core;
  stage_timer = native_timer_now();
  status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, retained_core,
                 retained_left, m, 1.0, 0.0, ATU_ret.data(), n, &workspace);
  stage_apply_seconds += native_timer_elapsed(stage_timer);
  if (status != 0) {
    return R_NilValue;
  }
  matvecs += retained_core;

  append_orthonormal_block(
    Q.data(), n, capacity, &q_cols, retained_right, retained_core, tol,
    &orthogonalization_passes, requested_orthogonalization_passes, bpro_ptr
  );
  if (q_cols == retained_core) {
    std::memcpy(AQ.data(), AV_ret.data(),
                sizeof(double) * static_cast<size_t>(m) * static_cast<size_t>(retained_core));
    aq_cols = q_cols;
  } else if (q_cols > aq_cols) {
    status = apply_augmented_basis_columns(
      impl, apply, m, n, Q.data(), aq_cols, q_cols, AQ.data(),
      &stage_apply_seconds, &matvecs, &workspace
    );
    if (status != 0) {
      return R_NilValue;
    }
    aq_cols = q_cols;
  }
  if (q_cols < rank) {
    append_orthonormal_column(
      Q.data(), n, capacity, &q_cols, initial_start, tol,
      &orthogonalization_passes, nullptr, requested_orthogonalization_passes,
      bpro_ptr
    );
    if (q_cols > aq_cols) {
      status = apply_augmented_basis_columns(
        impl, apply, m, n, Q.data(), aq_cols, q_cols, AQ.data(),
        &stage_apply_seconds, &matvecs, &workspace
      );
      if (status != 0) {
        return R_NilValue;
      }
      aq_cols = q_cols;
    }
  }

  stage_timer = native_timer_now();
  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans, &notrans, &retained_core, &retained_core, &m,
                  &one, const_cast<double*>(retained_left), &m,
                  AV_ret.data(), &m, &zero, H.data(), &retained_core FCONE FCONE);
  stage_projected_seconds += native_timer_elapsed(stage_timer);

  stage_timer = native_timer_now();
  std::memcpy(residual.data(), ATU_ret.data(),
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(retained_core));
  const double minus_one = -1.0;
  F77_CALL(dgemm)(&notrans, &trans, &n, &retained_core, &retained_core,
                  &minus_one, const_cast<double*>(retained_right), &n,
                  H.data(), &retained_core, &one, residual.data(), &n FCONE FCONE);
  stage_recurrence_seconds += native_timer_elapsed(stage_timer);
  const int residual_cols_before = q_cols;
  append_orthonormal_block(
    Q.data(), n, capacity, &q_cols, residual.data(), retained_core, tol,
    &orthogonalization_passes, requested_orthogonalization_passes, bpro_ptr
  );
  const int residual_cols = q_cols - residual_cols_before;
  if (q_cols > aq_cols) {
    status = apply_augmented_basis_columns(
      impl, apply, m, n, Q.data(), aq_cols, q_cols, AQ.data(),
      &stage_apply_seconds, &matvecs, &workspace
    );
    if (status != 0) {
      return R_NilValue;
    }
    aq_cols = q_cols;
  }

  if (q_cols == 0) {
    append_orthonormal_column(
      Q.data(), n, capacity, &q_cols, initial_start, tol,
      &orthogonalization_passes, nullptr, requested_orthogonalization_passes,
      bpro_ptr
    );
    if (q_cols > aq_cols) {
      status = apply_augmented_basis_columns(
        impl, apply, m, n, Q.data(), aq_cols, q_cols, AQ.data(),
        &stage_apply_seconds, &matvecs, &workspace
      );
      if (status != 0) {
        return R_NilValue;
      }
      aq_cols = q_cols;
    }
  }
  if (q_cols == 0) {
    return R_NilValue;
  }

  const double* seed = Q.data() + static_cast<int64_t>(q_cols - 1) * n;
  std::memcpy(z.data(), seed, sizeof(double) * static_cast<size_t>(n));
  double beta_prev = 0.0;
  int tail_steps_taken = 0;
  for (int step = 0; step < requested_tail_steps && q_cols < capacity; ++step) {
    stage_timer = native_timer_now();
    status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, z.data(), n,
                   1.0, 0.0, u.data(), m, &workspace);
    stage_apply_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      break;
    }
    ++matvecs;
    const int current_v_col = q_cols - 1;
    if (aq_cols < current_v_col) {
      status = apply_augmented_basis_columns(
        impl, apply, m, n, Q.data(), aq_cols, current_v_col, AQ.data(),
        &stage_apply_seconds, &matvecs, &workspace
      );
      if (status != 0) {
        break;
      }
      aq_cols = current_v_col;
    }
    if (aq_cols == current_v_col) {
      std::memcpy(
        AQ.data() + static_cast<int64_t>(current_v_col) * m,
        u.data(),
        sizeof(double) * static_cast<size_t>(m)
      );
      aq_cols = q_cols;
    }
    stage_timer = native_timer_now();
    if (beta_prev != 0.0) {
      for (int row = 0; row < m; ++row) {
        u[row] -= beta_prev * u_prev[row];
      }
    }
    long double alpha_norm2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      alpha_norm2 += static_cast<long double>(u[row]) * u[row];
    }
    const double alpha_step = sqrt(static_cast<double>(alpha_norm2));
    if (!R_FINITE(alpha_step) || alpha_step <= 100.0 * DBL_EPSILON) {
      stage_recurrence_seconds += native_timer_elapsed(stage_timer);
      break;
    }
    for (int row = 0; row < m; ++row) {
      u[row] /= alpha_step;
    }
    stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, 1, u.data(), m,
                   1.0, 0.0, z.data(), n, &workspace);
    stage_apply_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      break;
    }
    ++matvecs;
    stage_timer = native_timer_now();
    const double* current_v = Q.data() + static_cast<int64_t>(q_cols - 1) * n;
    for (int row = 0; row < n; ++row) {
      z[row] -= alpha_step * current_v[row];
    }
    stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    const int before = q_cols;
    stage_timer = native_timer_now();
    double beta_next = 0.0;
    append_orthonormal_column(
      Q.data(), n, capacity, &q_cols, z.data(), tol,
      &orthogonalization_passes, &beta_next, requested_orthogonalization_passes,
      bpro_ptr
    );
    stage_reorthogonalization_seconds += native_timer_elapsed(stage_timer);
    if (q_cols == before) {
      break;
    }
    std::memcpy(u_prev.data(), u.data(), sizeof(double) * static_cast<size_t>(m));
    const double* accepted_v = Q.data() + static_cast<int64_t>(q_cols - 1) * n;
    std::memcpy(z.data(), accepted_v, sizeof(double) * static_cast<size_t>(n));
    beta_prev = beta_next;
    tail_beta_history.push_back(beta_next);
    ++tail_steps_taken;
  }

  if (q_cols < rank) {
    return R_NilValue;
  }
  if (q_cols > aq_cols) {
    status = apply_augmented_basis_columns(
      impl, apply, m, n, Q.data(), aq_cols, q_cols, AQ.data(),
      &stage_apply_seconds, &matvecs, &workspace
    );
    if (status != 0) {
      return R_NilValue;
    }
    aq_cols = q_cols;
  }
  const double augmented_basis_orthogonality_loss =
    basis_orthogonality_loss(Q.data(), n, q_cols);
  if (!R_FINITE(augmented_basis_orthogonality_loss) ||
      augmented_basis_orthogonality_loss > bpro.threshold) {
    bpro.escalation_recommended = 1;
  }

  std::vector<int> attempted_subspaces;
  std::vector<int> attempt_iterations;
  std::vector<int> attempt_matvecs;
  std::vector<int> attempt_warm_started;
  std::vector<double> attempt_cheap_residuals;
  const int chunks = max_restarts + 1;
  attempted_subspaces.reserve(static_cast<size_t>(chunks));
  attempt_iterations.reserve(static_cast<size_t>(chunks));
  attempt_matvecs.reserve(static_cast<size_t>(chunks));
  attempt_warm_started.reserve(static_cast<size_t>(chunks));
  attempt_cheap_residuals.reserve(static_cast<size_t>(chunks));
  for (int attempt = 0; attempt < chunks; ++attempt) {
    int tail_for_attempt = tail_width * (attempt + 1);
    if (tail_for_attempt > tail_steps_taken) {
      tail_for_attempt = tail_steps_taken;
    }
    int subspace = retained_core + residual_cols + tail_for_attempt;
    if (subspace > q_cols) {
      subspace = q_cols;
    }
    attempted_subspaces.push_back(subspace);
    attempt_iterations.push_back(tail_for_attempt);
    attempt_matvecs.push_back(
      2 * retained_core + residual_cols + 2 * tail_for_attempt +
        (tail_for_attempt > 0 ? 1 : 0)
    );
    attempt_warm_started.push_back(attempt > 0 ? 1 : 0);
    double cheap_residual = R_NaReal;
    if (tail_for_attempt > 0 &&
        tail_for_attempt <= static_cast<int>(tail_beta_history.size())) {
      cheap_residual = tail_beta_history[static_cast<size_t>(tail_for_attempt - 1)];
    }
    attempt_cheap_residuals.push_back(cheap_residual);
  }

  SEXP ritz_ = R_NilValue;
  int small_svds = 0;
  double min_cheap_residual = R_PosInf;
  double final_cheap_residual = R_NaReal;
  for (int attempt = 0; attempt < chunks; ++attempt) {
    const int subspace = attempted_subspaces[static_cast<size_t>(attempt)];
    stage_timer = native_timer_now();
    SEXP attempt_ritz_ = PROTECT(eigencore_block_golub_kahan_ritz_from_ptr(
      Q.data(), n, AQ.data(), m, subspace, rank, target_kind
    ));
    stage_projected_seconds += native_timer_elapsed(stage_timer);
    ++small_svds;
    const double cheap_residual = attempt_cheap_residuals[static_cast<size_t>(attempt)];
    if (R_FINITE(cheap_residual) && cheap_residual < min_cheap_residual) {
      min_cheap_residual = cheap_residual;
    }
    if (attempt + 1 == chunks) {
      ritz_ = attempt_ritz_;
      final_cheap_residual = cheap_residual;
    } else {
      UNPROTECT(1);
    }
  }
  if (!R_FINITE(min_cheap_residual)) {
    min_cheap_residual = R_NaReal;
  }
  SEXP history_ = PROTECT(irlba_lbd_attempt_history_pack(
    attempted_subspaces, attempt_iterations, attempt_matvecs, attempt_warm_started,
    &attempt_cheap_residuals
  ));

  const int from_scratch_matvecs = 2 * retained_core + 2 * tail_steps_taken + q_cols;
  const int cached_matvec_savings = from_scratch_matvecs - matvecs;
  SEXP out_ = PROTECT(allocVector(VECSXP, 42));
  const double augmented_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      static_cast<int64_t>(n) * static_cast<int64_t>(capacity) +
      static_cast<int64_t>(m) * static_cast<int64_t>(capacity) +
      static_cast<int64_t>(m) * static_cast<int64_t>(retained_core) +
      static_cast<int64_t>(n) * static_cast<int64_t>(retained_core) +
      static_cast<int64_t>(retained_core) * static_cast<int64_t>(retained_core) +
      static_cast<int64_t>(m + m + n)
    );
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(tail_steps_taken));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(tail_steps_taken));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(max_restarts));
  SET_VECTOR_ELT(out_, 7, history_);
  SET_VECTOR_ELT(out_, 8, ScalarReal(
    augmented_workspace_bytes > native_workspace_bytes
      ? augmented_workspace_bytes
      : native_workspace_bytes
  ));
  SET_VECTOR_ELT(out_, 9, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_projected_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarInteger(orthogonalization_passes));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalize_u));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalize_v));
  SET_VECTOR_ELT(out_, 16, ScalarInteger(work));
  SET_VECTOR_ELT(out_, 17, ScalarInteger(retained));
  SET_VECTOR_ELT(out_, 18, mkString("residual_augmented_projection"));
  SET_VECTOR_ELT(out_, 19, ScalarLogical(1));
  SET_VECTOR_ELT(out_, 20, ScalarLogical(1));
  SET_VECTOR_ELT(out_, 21, ScalarInteger(residual_cols));
  SET_VECTOR_ELT(out_, 22, ScalarInteger(tail_steps_taken));
  SET_VECTOR_ELT(out_, 23, ScalarInteger(q_cols));
  SET_VECTOR_ELT(out_, 24, ScalarLogical(bpro_policy ? 1 : 0));
  SET_VECTOR_ELT(out_, 25, ScalarInteger(requested_orthogonalization_passes));
  SET_VECTOR_ELT(out_, 26, ScalarReal(bpro.threshold));
  SET_VECTOR_ELT(out_, 27, ScalarInteger(bpro.monitored_appends));
  SET_VECTOR_ELT(out_, 28, ScalarInteger(bpro.threshold_reorthogonalizations));
  SET_VECTOR_ELT(out_, 29, ScalarReal(bpro.max_estimated_orthogonality_loss));
  SET_VECTOR_ELT(out_, 30, ScalarReal(bpro.max_post_append_orthogonality_loss));
  SET_VECTOR_ELT(out_, 31, ScalarReal(augmented_basis_orthogonality_loss));
  SET_VECTOR_ELT(out_, 32, ScalarLogical(bpro.escalation_recommended ? 1 : 0));
  SET_VECTOR_ELT(out_, 33, ScalarInteger(chunks));
  SET_VECTOR_ELT(out_, 34, ScalarInteger(retained_core));
  SET_VECTOR_ELT(out_, 35, ScalarInteger(small_svds));
  SET_VECTOR_ELT(out_, 36, ScalarInteger(aq_cols));
  SET_VECTOR_ELT(out_, 37, ScalarInteger(from_scratch_matvecs));
  SET_VECTOR_ELT(out_, 38, ScalarInteger(cached_matvec_savings));
  SET_VECTOR_ELT(out_, 39, ScalarReal(min_cheap_residual));
  SET_VECTOR_ELT(out_, 40, ScalarReal(final_cheap_residual));
  SET_VECTOR_ELT(out_, 41, ScalarLogical(cached_matvec_savings > 0 ? 1 : 0));
  SEXP names_ = PROTECT(allocVector(STRSXP, 42));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("final_iterations"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("restart_count"));
  SET_STRING_ELT(names_, 7, mkChar("attempt_history"));
  SET_STRING_ELT(names_, 8, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 9, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_projected_solve_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalize_u"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalize_v"));
  SET_STRING_ELT(names_, 16, mkChar("work"));
  SET_STRING_ELT(names_, 17, mkChar("retained"));
  SET_STRING_ELT(names_, 18, mkChar("restart_state_kind"));
  SET_STRING_ELT(names_, 19, mkChar("recurrence_available"));
  SET_STRING_ELT(names_, 20, mkChar("augmented_recurrence"));
  SET_STRING_ELT(names_, 21, mkChar("residual_augmented_cols"));
  SET_STRING_ELT(names_, 22, mkChar("augmented_tail_steps"));
  SET_STRING_ELT(names_, 23, mkChar("augmented_basis_cols"));
  SET_STRING_ELT(names_, 24, mkChar("bpro_policy"));
  SET_STRING_ELT(names_, 25, mkChar("bpro_reorthogonalization_passes_per_append"));
  SET_STRING_ELT(names_, 26, mkChar("bpro_monitoring_threshold"));
  SET_STRING_ELT(names_, 27, mkChar("bpro_monitored_appends"));
  SET_STRING_ELT(names_, 28, mkChar("bpro_threshold_reorthogonalizations"));
  SET_STRING_ELT(names_, 29, mkChar("bpro_max_estimated_orthogonality_loss"));
  SET_STRING_ELT(names_, 30, mkChar("bpro_max_post_append_orthogonality_loss"));
  SET_STRING_ELT(names_, 31, mkChar("bpro_augmented_basis_orthogonality_loss"));
  SET_STRING_ELT(names_, 32, mkChar("bpro_escalation_recommended"));
  SET_STRING_ELT(names_, 33, mkChar("augmented_restart_cycles"));
  SET_STRING_ELT(names_, 34, mkChar("augmented_kept_vectors"));
  SET_STRING_ELT(names_, 35, mkChar("augmented_small_svds"));
  SET_STRING_ELT(names_, 36, mkChar("augmented_cached_aq_cols"));
  SET_STRING_ELT(names_, 37, mkChar("augmented_from_scratch_matvecs"));
  SET_STRING_ELT(names_, 38, mkChar("augmented_matvec_savings"));
  SET_STRING_ELT(names_, 39, mkChar("augmented_min_cheap_residual"));
  SET_STRING_ELT(names_, 40, mkChar("augmented_final_cheap_residual"));
  SET_STRING_ELT(names_, 41, mkChar("augmented_reduces_from_scratch_work"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}

template <typename ConfigureOperator>
static SEXP irlba_lbd_retained_impl(ConfigureOperator configure_operator,
                                    int m,
                                    int n,
                                    const double* initial_start,
                                    const double* retained_right,
                                    const double* retained_left,
                                    const double* alpha,
                                    const double* beta,
                                    const double* random_tails,
                                    int work,
                                    int retained,
                                    int max_restarts,
                                    int rank,
                                    int target_kind,
                                    double tol,
                                    int reorth_policy) {
  (void) alpha;
  (void) beta;
  const int attempts = (max_restarts < 0) ? 1 : max_restarts + 1;
  const int random_tail_cols = work - retained;
  const int reorthogonalize_u =
    (reorth_policy != 1) || (m <= n);
  const int reorthogonalize_v =
    (reorth_policy != 1) || (n <= m);
  const int use_blas_reorthogonalization = 0;
  const int enable_projected_stop = 0;
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(work) +
      static_cast<int64_t>(2 * m + 2 * n)
    );

  void* augmented_impl = nullptr;
  EigencoreApplyFn augmented_apply = nullptr;
  configure_operator(&augmented_impl, &augmented_apply);
  SEXP augmented_ = PROTECT(irlba_lbd_augmented_retained_projection(
    augmented_impl, augmented_apply, m, n, initial_start, retained_right,
    retained_left, random_tails, work, retained, max_restarts, rank,
    target_kind, tol, reorthogonalize_u, reorthogonalize_v,
    native_workspace_bytes, reorth_policy == 3 ? 1 : 0
  ));
  if (augmented_ != R_NilValue) {
    UNPROTECT(1);
    return augmented_;
  }
  UNPROTECT(1);

  std::vector<double> start(static_cast<size_t>(n), 0.0);
  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(work), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(work), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(work), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(work), 0.0);
  std::vector<int> attempted_subspaces;
  std::vector<int> attempt_iterations;
  std::vector<int> attempt_matvecs;
  std::vector<int> attempt_warm_started;
  attempted_subspaces.reserve(static_cast<size_t>(attempts));
  attempt_iterations.reserve(static_cast<size_t>(attempts));
  attempt_matvecs.reserve(static_cast<size_t>(attempts));
  attempt_warm_started.reserve(static_cast<size_t>(attempts));

  SEXP final_ritz_ = R_NilValue;
  int final_iterations = 0;
  int total_iterations = 0;
  int total_matvecs = 0;
  int total_reorthogonalization_passes = 0;
  double total_apply_seconds = 0.0;
  double total_recurrence_seconds = 0.0;
  double total_reorthogonalization_seconds = 0.0;
  double total_projected_seconds = 0.0;

  for (int attempt = 0; attempt < attempts; ++attempt) {
    irlba_lbd_retained_seed(
      n, retained, initial_start, retained_right, random_tails,
      random_tail_cols, attempt, start.data()
    );
    std::fill(U_work.begin(), U_work.end(), 0.0);
    std::fill(V_work.begin(), V_work.end(), 0.0);
    std::fill(alpha_work.begin(), alpha_work.end(), 0.0);
    std::fill(beta_work.begin(), beta_work.end(), 0.0);

    void* impl = nullptr;
    EigencoreApplyFn apply = nullptr;
    configure_operator(&impl, &apply);
    int iterations = 0;
    int matvecs = 0;
    int projected_stop = 0;
    int projected_nconv = 0;
    double projected_max_residual = R_PosInf;
    int projected_checks = 0;
    double projected_seconds = 0.0;
    double stage_apply_seconds = 0.0;
    double stage_recurrence_seconds = 0.0;
    double stage_reorthogonalization_seconds = 0.0;
    int reorthogonalization_passes = 0;
    const int status = native_golub_kahan_run(
      impl, apply, m, n, work, rank, target_kind, tol,
      enable_projected_stop, use_blas_reorthogonalization,
      start.data(), U_work.data(), V_work.data(), alpha_work.data(),
      beta_work.data(), &iterations, &matvecs, &projected_stop,
      &projected_nconv, &projected_max_residual, &projected_checks,
      &projected_seconds, &stage_apply_seconds, &stage_recurrence_seconds,
      &stage_reorthogonalization_seconds, &reorthogonalization_passes,
      reorthogonalize_u, reorthogonalize_v
    );
    if (status != 0) {
      error("native retained one-sided IRLBA/LBD failed with status=%d", status);
    }
    attempted_subspaces.push_back(work);
    attempt_iterations.push_back(iterations);
    attempt_matvecs.push_back(matvecs);
    attempt_warm_started.push_back(attempt > 0 ? 1 : 0);
    final_iterations = iterations;
    total_iterations += iterations;
    total_matvecs += matvecs;
    total_reorthogonalization_passes += reorthogonalization_passes;
    total_apply_seconds += stage_apply_seconds;
    total_recurrence_seconds += stage_recurrence_seconds;
    total_reorthogonalization_seconds += stage_reorthogonalization_seconds;
    total_projected_seconds += projected_seconds;

    if (attempt + 1 == attempts) {
      final_ritz_ = PROTECT(eigencore_golub_kahan_ritz_from_ptr(
        U_work.data(), V_work.data(), m, n, iterations, alpha_work.data(),
        beta_work.data(), rank, target_kind
      ));
    }
  }

  SEXP history_ = PROTECT(irlba_lbd_attempt_history_pack(
    attempted_subspaces, attempt_iterations, attempt_matvecs, attempt_warm_started
  ));
  SEXP out_ = PROTECT(allocVector(VECSXP, 18));
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(final_ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(final_ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(final_ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(final_iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(total_iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(total_matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(attempts - 1));
  SET_VECTOR_ELT(out_, 7, history_);
  SET_VECTOR_ELT(out_, 8, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 9, ScalarReal(total_apply_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(total_recurrence_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(total_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(total_projected_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarInteger(total_reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalize_u));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalize_v));
  SET_VECTOR_ELT(out_, 16, ScalarInteger(work));
  SET_VECTOR_ELT(out_, 17, ScalarInteger(retained));
  SEXP names_ = PROTECT(allocVector(STRSXP, 18));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("final_iterations"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("restart_count"));
  SET_STRING_ELT(names_, 7, mkChar("attempt_history"));
  SET_STRING_ELT(names_, 8, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 9, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_projected_solve_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalize_u"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalize_v"));
  SET_STRING_ELT(names_, 16, mkChar("work"));
  SET_STRING_ELT(names_, 17, mkChar("retained"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_irlba_lbd_dense_retained(SEXP A_, SEXP initial_start_,
                                                   SEXP retained_right_,
                                                   SEXP retained_left_,
                                                   SEXP alpha_, SEXP beta_,
                                                   SEXP random_tails_,
                                                   SEXP work_, SEXP retained_,
                                                   SEXP max_restarts_,
                                                   SEXP rank_,
                                                   SEXP target_kind_,
                                                   SEXP tol_,
                                                   SEXP reorth_policy_) {
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue || LENGTH(dimA) != 2) {
    error("A must be a double matrix");
  }
  validate_irlba_lbd_retained_contract(
    INTEGER(dimA)[0], INTEGER(dimA)[1], initial_start_, retained_right_,
    retained_left_, alpha_, beta_, random_tails_, work_, retained_,
    max_restarts_, rank_, target_kind_, tol_, reorth_policy_
  );
  DenseColumnMajorOperator impl_holder = {
    INTEGER(dimA)[0], INTEGER(dimA)[1], REAL(A_)
  };
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_dense_apply;
  };
  return irlba_lbd_retained_impl(
    configure, INTEGER(dimA)[0], INTEGER(dimA)[1],
    REAL(initial_start_), REAL(retained_right_), REAL(retained_left_),
    REAL(alpha_), REAL(beta_), REAL(random_tails_),
    asInteger(work_), asInteger(retained_), asInteger(max_restarts_),
    asInteger(rank_), asInteger(target_kind_), asReal(tol_),
    asInteger(reorth_policy_)
  );
}

extern "C" SEXP eigencore_irlba_lbd_csc_retained(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_,
                                                 SEXP initial_start_,
                                                 SEXP retained_right_,
                                                 SEXP retained_left_,
                                                 SEXP alpha_, SEXP beta_,
                                                 SEXP random_tails_,
                                                 SEXP work_, SEXP retained_,
                                                 SEXP max_restarts_,
                                                 SEXP rank_,
                                                 SEXP target_kind_,
                                                 SEXP tol_,
                                                 SEXP reorth_policy_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      LENGTH(dim_) != 2) {
    error("invalid CSC retained IRLBA/LBD inputs");
  }
  validate_irlba_lbd_retained_contract(
    INTEGER(dim_)[0], INTEGER(dim_)[1], initial_start_, retained_right_,
    retained_left_, alpha_, beta_, random_tails_, work_, retained_,
    max_restarts_, rank_, target_kind_, tol_, reorth_policy_
  );
  CSCOperator impl_holder = {
    INTEGER(dim_)[0], INTEGER(dim_)[1], INTEGER(i_), INTEGER(p_), REAL(x_)
  };
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_csc_apply;
  };
  return irlba_lbd_retained_impl(
    configure, INTEGER(dim_)[0], INTEGER(dim_)[1],
    REAL(initial_start_), REAL(retained_right_), REAL(retained_left_),
    REAL(alpha_), REAL(beta_), REAL(random_tails_),
    asInteger(work_), asInteger(retained_), asInteger(max_restarts_),
    asInteger(rank_), asInteger(target_kind_), asReal(tol_),
    asInteger(reorth_policy_)
  );
}
