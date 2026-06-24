#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "eigencore_common.h"
#include "native_operators.h"
#include "scalar_krylov.h"

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

static int lanczos_convergence_estimate(const double* alpha,
                                        const double* beta,
                                        int iter,
                                        int k,
                                        int target_kind,
                                        double tol,
                                        int* nconv,
                                        double* max_residual) {
  *nconv = 0;
  *max_residual = R_PosInf;
  if (iter <= 0 || k <= 0) {
    return 0;
  }

  double* diag = static_cast<double*>(std::calloc(static_cast<size_t>(iter), sizeof(double)));
  double* offdiag = (iter > 1)
    ? static_cast<double*>(std::calloc(static_cast<size_t>(iter - 1), sizeof(double)))
    : nullptr;
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  double* work = (iter > 1)
    ? static_cast<double*>(std::calloc(static_cast<size_t>(2 * iter - 2), sizeof(double)))
    : nullptr;
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>((k < iter) ? k : iter), sizeof(int)));
  if (diag == nullptr || z == nullptr || selected == nullptr ||
      (iter > 1 && (offdiag == nullptr || work == nullptr))) {
    std::free(diag);
    std::free(offdiag);
    std::free(z);
    std::free(work);
    std::free(selected);
    return -2;
  }

  for (int i = 0; i < iter; ++i) {
    diag[i] = alpha[i];
  }
  for (int i = 0; i < iter - 1; ++i) {
    offdiag[i] = beta[i];
  }

  if (iter == 1) {
    z[0] = 1.0;
  } else {
    char jobz = 'V';
    int info = 0;
    F77_CALL(dstev)(&jobz, &iter, diag, offdiag, z, &iter, work, &info FCONE);
    if (info != 0) {
      std::free(diag);
      std::free(offdiag);
      std::free(z);
      std::free(work);
      std::free(selected);
      return info;
    }
  }

  const int count = selected_ritz_indices(diag, iter, k, target_kind, selected);
  double max_selected = 0.0;
  int conv = 0;
  for (int i = 0; i < count; ++i) {
    const int idx = selected[i];
    const double residual = fabs(beta[iter - 1] * z[(iter - 1) + idx * iter]);
    const double threshold = tol * ((fabs(diag[idx]) > 1.0) ? fabs(diag[idx]) : 1.0);
    if (residual <= threshold) {
      ++conv;
    }
    if (residual > max_selected) {
      max_selected = residual;
    }
  }

  *nconv = conv;
  *max_residual = max_selected;
  std::free(diag);
  std::free(offdiag);
  std::free(z);
  std::free(work);
  std::free(selected);
  return 0;
}

static int golub_kahan_projected_convergence_estimate(const double* alpha,
                                                      const double* beta,
                                                      int iter,
                                                      int k,
                                                      int target_kind,
                                                      double tol,
                                                      int* nconv,
                                                      double* max_residual) {
  *nconv = 0;
  *max_residual = R_PosInf;
  if (iter <= 0 || k <= 0) {
    return 0;
  }

  double* B = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  double* d = static_cast<double*>(std::calloc(static_cast<size_t>(iter), sizeof(double)));
  double* u = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  double* vt = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>((k < iter) ? k : iter), sizeof(int)));
  if (B == nullptr || d == nullptr || u == nullptr || vt == nullptr || selected == nullptr) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(selected);
    return -2;
  }

  for (int i = 0; i < iter; ++i) {
    B[i + i * iter] = alpha[i];
  }
  for (int i = 0; i < iter - 1; ++i) {
    B[i + (i + 1) * iter] = beta[i];
  }

  char jobu = 'A';
  char jobvt = 'N';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgesvd)(&jobu, &jobvt, &iter, &iter, B, &iter,
                   d, u, &iter, vt, &iter, &work_query, &lwork,
                   &info FCONE FCONE);
  if (info != 0) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(selected);
    return info;
  }
  lwork = static_cast<int>(work_query);
  double* work = static_cast<double*>(std::calloc(static_cast<size_t>(lwork), sizeof(double)));
  if (work == nullptr) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(selected);
    return -2;
  }

  std::memset(B, 0, sizeof(double) * static_cast<size_t>(iter) * iter);
  for (int i = 0; i < iter; ++i) {
    B[i + i * iter] = alpha[i];
  }
  for (int i = 0; i < iter - 1; ++i) {
    B[i + (i + 1) * iter] = beta[i];
  }
  F77_CALL(dgesvd)(&jobu, &jobvt, &iter, &iter, B, &iter,
                   d, u, &iter, vt, &iter, work, &lwork,
                   &info FCONE FCONE);
  if (info != 0) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(work);
    std::free(selected);
    return info;
  }

  const int count = selected_ritz_indices(d, iter, k, target_kind, selected);
  double max_selected = 0.0;
  int conv = 0;
  for (int i = 0; i < count; ++i) {
    const int idx = selected[i];
    const double residual = fabs(beta[iter - 1] * u[(iter - 1) + idx * iter]);
    const double scale = (fabs(d[idx]) > 1.0) ? fabs(d[idx]) : 1.0;
    const double threshold = tol * scale;
    if (residual <= threshold) {
      ++conv;
    }
    if (residual > max_selected) {
      max_selected = residual;
    }
  }

  *nconv = conv;
  *max_residual = max_selected;
  std::free(B);
  std::free(d);
  std::free(u);
  std::free(vt);
  std::free(work);
  std::free(selected);
  return 0;
}

static int native_lanczos_run(void* impl,
                              EigencoreApplyFn apply,
                              int n,
                              int maxit,
                              int k,
                              int target_kind,
                              double tol,
                              const double* start,
                              double* Q,
                              double* alpha,
                              double* beta,
                              int* history_nconv,
                              double* history_max_residual,
                              int* iterations,
                              int* matvecs) {
  double* q = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* q_prev = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* coeff = static_cast<double*>(std::calloc(static_cast<size_t>(maxit), sizeof(double)));
  if (q == nullptr || q_prev == nullptr || z == nullptr || coeff == nullptr) {
    std::free(q);
    std::free(q_prev);
    std::free(z);
    std::free(coeff);
    return -2;
  }

  long double start_norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    start_norm2 += static_cast<long double>(start[row]) * start[row];
  }
  double q_norm = sqrt(static_cast<double>(start_norm2));
  if (q_norm == 0.0) {
    q[0] = 1.0;
  } else {
    for (int row = 0; row < n; ++row) {
      q[row] = start[row] / q_norm;
    }
  }

  *iterations = 0;
  *matvecs = 0;
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  for (int j = 0; j < maxit; ++j) {
    *iterations = j + 1;
    std::memcpy(Q + j * n, q, sizeof(double) * static_cast<size_t>(n));
    std::memset(z, 0, sizeof(double) * static_cast<size_t>(n));

    const int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, q, n,
                             1.0, 0.0, z, n, &workspace);
    if (status != 0) {
      std::free(q);
      std::free(q_prev);
      std::free(z);
      std::free(coeff);
      return status;
    }
    ++(*matvecs);

    if (j > 0) {
      const double bj = beta[j - 1];
      for (int row = 0; row < n; ++row) {
        z[row] -= bj * q_prev[row];
      }
    }

    long double aj = 0.0L;
    for (int row = 0; row < n; ++row) {
      aj += static_cast<long double>(q[row]) * z[row];
    }
    alpha[j] = static_cast<double>(aj);
    for (int row = 0; row < n; ++row) {
      z[row] -= alpha[j] * q[row];
    }

    // Block CGS full reorthogonalization against the stored basis with an
    // adaptive DGKS second pass: the second projection runs only when the
    // first cancelled a large fraction of ||z|| (post < eta * pre,
    // eta = 1/sqrt(2)), which carries the same orthogonality guarantee as
    // unconditional CGS2 ("twice is enough").
    {
      const char trans_T = 'T';
      const char trans_N = 'N';
      const int one_col = 1;
      const double one = 1.0;
      const double zero = 0.0;
      const double minus_one = -1.0;
      const double dgks_eta = 0.7071067811865475;
      const int active = j + 1;
      long double pre2 = 0.0L;
      for (int row = 0; row < n; ++row) {
        pre2 += static_cast<long double>(z[row]) * z[row];
      }
      const double pre_norm = sqrt(static_cast<double>(pre2));
      for (int pass = 0; pass < 2; ++pass) {
        F77_CALL(dgemv)(&trans_T, &n, &active, &one,
                        Q, &n, z, &one_col,
                        &zero, coeff, &one_col FCONE);
        F77_CALL(dgemv)(&trans_N, &n, &active, &minus_one,
                        Q, &n, coeff, &one_col,
                        &one, z, &one_col FCONE);
        if (pass == 0) {
          long double post2 = 0.0L;
          for (int row = 0; row < n; ++row) {
            post2 += static_cast<long double>(z[row]) * z[row];
          }
          const double post_norm = sqrt(static_cast<double>(post2));
          if (post_norm >= dgks_eta * pre_norm) {
            break;
          }
        }
      }
    }

    long double beta_norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      beta_norm2 += static_cast<long double>(z[row]) * z[row];
    }
    beta[j] = sqrt(static_cast<double>(beta_norm2));
    history_nconv[j] = 0;
    history_max_residual[j] = R_PosInf;
    if (j + 1 >= k) {
      int nconv = 0;
      double max_residual = R_PosInf;
      const int conv_status = lanczos_convergence_estimate(
        alpha, beta, j + 1, k, target_kind, tol, &nconv, &max_residual
      );
      if (conv_status != 0) {
        std::free(q);
        std::free(q_prev);
        std::free(z);
        std::free(coeff);
        return conv_status;
      }
      history_nconv[j] = nconv;
      history_max_residual[j] = max_residual;
      if (nconv >= k) {
        break;
      }
    }
    if (j + 1 == maxit || beta[j] <= 100.0 * DBL_EPSILON) {
      break;
    }

    std::memcpy(q_prev, q, sizeof(double) * static_cast<size_t>(n));
    for (int row = 0; row < n; ++row) {
      q[row] = z[row] / beta[j];
    }
  }

  std::free(q);
  std::free(q_prev);
  std::free(z);
  std::free(coeff);
  return 0;
}

int native_golub_kahan_run(void* impl,
                                  EigencoreApplyFn apply,
                                  int m,
                                  int n,
                                  int maxit,
                                  int k,
                                  int target_kind,
                                  double tol,
                                  int enable_projected_stop,
                                  int use_blas_reorthogonalization,
                                  const double* start,
                                  double* U,
                                  double* V,
                                  double* alpha,
                                  double* beta,
                                  int* iterations,
                                  int* matvecs,
                                  int* projected_stop,
                                  int* projected_nconv,
                                  double* projected_max_residual,
                                  int* projected_checks,
                                  double* projected_seconds,
                                  double* stage_apply_seconds,
                                  double* stage_recurrence_seconds,
                                  double* stage_reorthogonalization_seconds,
                                  int* reorthogonalization_passes,
                                  int reorthogonalize_u,
                                  int reorthogonalize_v) {
  double* v = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* u = static_cast<double*>(std::calloc(static_cast<size_t>(m), sizeof(double)));
  double* u_prev = static_cast<double*>(std::calloc(static_cast<size_t>(m), sizeof(double)));
  double* coeff = use_blas_reorthogonalization
    ? static_cast<double*>(std::calloc(static_cast<size_t>(maxit), sizeof(double)))
    : nullptr;
  if (v == nullptr || z == nullptr || u == nullptr || u_prev == nullptr || coeff == nullptr) {
    if (!use_blas_reorthogonalization && v != nullptr && z != nullptr &&
        u != nullptr && u_prev != nullptr) {
      // coeff is intentionally unused on the scalar sparse path.
    } else {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      std::free(coeff);
      return -2;
    }
  }
  if (v == nullptr || z == nullptr || u == nullptr || u_prev == nullptr) {
    std::free(v);
    std::free(z);
    std::free(u);
    std::free(u_prev);
    std::free(coeff);
    return -2;
  }

  long double start_norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    start_norm2 += static_cast<long double>(start[row]) * start[row];
  }
  double v_norm = sqrt(static_cast<double>(start_norm2));
  if (v_norm == 0.0) {
    v[0] = 1.0;
  } else {
    for (int row = 0; row < n; ++row) {
      v[row] = start[row] / v_norm;
    }
  }

  *iterations = 0;
  *matvecs = 0;
  *projected_stop = 0;
  *projected_nconv = 0;
  *projected_max_residual = R_PosInf;
  *projected_checks = 0;
  *projected_seconds = 0.0;
  *stage_apply_seconds = 0.0;
  *stage_recurrence_seconds = 0.0;
  *stage_reorthogonalization_seconds = 0.0;
  *reorthogonalization_passes = 0;
  double beta_prev = 0.0;
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const int check_interval = (2 * k > 10) ? 2 * k : 10;
  const int min_projected_savings = (k > 5) ? k : 5;

  for (int j = 0; j < maxit; ++j) {
    *iterations = j + 1;
    std::memcpy(V + j * n, v, sizeof(double) * static_cast<size_t>(n));
    std::memset(u, 0, sizeof(double) * static_cast<size_t>(m));

    auto stage_timer = native_timer_now();
    int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, v, n,
                       1.0, 0.0, u, m, &workspace);
    *stage_apply_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      std::free(coeff);
      return status;
    }
    ++(*matvecs);
    stage_timer = native_timer_now();
    if (j > 0) {
      for (int row = 0; row < m; ++row) {
        u[row] -= beta_prev * u_prev[row];
      }
    }
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    if (reorthogonalize_u && j > 0) {
      double norm_before2 = 0.0;
      for (int row = 0; row < m; ++row) {
        norm_before2 += u[row] * u[row];
      }
      const double norm_before = sqrt(norm_before2);
      const char trans_T = 'T';
      const char trans_N = 'N';
      const int one_col = 1;
      const double one = 1.0;
      const double zero = 0.0;
      const double minus_one = -1.0;
      int passes = 1;
      for (int pass = 0; pass < passes; ++pass) {
        if (use_blas_reorthogonalization) {
        F77_CALL(dgemv)(&trans_T, &m, &j, &one,
                        U, &m, u, &one_col,
                        &zero, coeff, &one_col FCONE);
        F77_CALL(dgemv)(&trans_N, &m, &j, &minus_one,
                        U, &m, coeff, &one_col,
                        &one, u, &one_col FCONE);
        } else {
          for (int prev = 0; prev < j; ++prev) {
            const double* uprev_basis = U + prev * m;
            double dot = 0.0;
            for (int row = 0; row < m; ++row) {
              dot += uprev_basis[row] * u[row];
            }
            const double scalar_coeff = dot;
            for (int row = 0; row < m; ++row) {
              u[row] -= scalar_coeff * uprev_basis[row];
            }
          }
        }
        ++(*reorthogonalization_passes);
        if (pass == 0 && norm_before > 0.0) {
          double norm_after2 = 0.0;
          for (int row = 0; row < m; ++row) {
            norm_after2 += u[row] * u[row];
          }
          const double norm_after = sqrt(norm_after2);
          if (norm_after < 0.717 * norm_before) {
            passes = 2;
          }
        }
      }
    }
    *stage_reorthogonalization_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    long double alpha_norm2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      alpha_norm2 += static_cast<long double>(u[row]) * u[row];
    }
    alpha[j] = sqrt(static_cast<double>(alpha_norm2));
    if (alpha[j] <= 100.0 * DBL_EPSILON) {
      break;
    }
    for (int row = 0; row < m; ++row) {
      u[row] /= alpha[j];
    }
    std::memcpy(U + j * m, u, sizeof(double) * static_cast<size_t>(m));
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    std::memset(z, 0, sizeof(double) * static_cast<size_t>(n));
    stage_timer = native_timer_now();
    status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, 1, u, m,
                   1.0, 0.0, z, n, &workspace);
    *stage_apply_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      std::free(coeff);
      return status;
    }
    ++(*matvecs);
    stage_timer = native_timer_now();
    for (int row = 0; row < n; ++row) {
      z[row] -= alpha[j] * v[row];
    }
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    if (reorthogonalize_v) {
      const int active_v = j + 1;
      double norm_before2 = 0.0;
      for (int row = 0; row < n; ++row) {
        norm_before2 += z[row] * z[row];
      }
      const double norm_before = sqrt(norm_before2);
      const char trans_T = 'T';
      const char trans_N = 'N';
      const int one_col = 1;
      const double one = 1.0;
      const double zero = 0.0;
      const double minus_one = -1.0;
      int passes = 1;
      for (int pass = 0; pass < passes; ++pass) {
        if (use_blas_reorthogonalization) {
        F77_CALL(dgemv)(&trans_T, &n, &active_v, &one,
                        V, &n, z, &one_col,
                        &zero, coeff, &one_col FCONE);
        F77_CALL(dgemv)(&trans_N, &n, &active_v, &minus_one,
                        V, &n, coeff, &one_col,
                        &one, z, &one_col FCONE);
        } else {
          for (int prev = 0; prev <= j; ++prev) {
            const double* vprev_basis = V + prev * n;
            double dot = 0.0;
            for (int row = 0; row < n; ++row) {
              dot += vprev_basis[row] * z[row];
            }
            const double scalar_coeff = dot;
            for (int row = 0; row < n; ++row) {
              z[row] -= scalar_coeff * vprev_basis[row];
            }
          }
        }
        ++(*reorthogonalization_passes);
        if (pass == 0 && norm_before > 0.0) {
          double norm_after2 = 0.0;
          for (int row = 0; row < n; ++row) {
            norm_after2 += z[row] * z[row];
          }
          const double norm_after = sqrt(norm_after2);
          if (norm_after < 0.717 * norm_before) {
            passes = 2;
          }
        }
      }
    }
    *stage_reorthogonalization_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    long double beta_norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      beta_norm2 += static_cast<long double>(z[row]) * z[row];
    }
    beta[j] = sqrt(static_cast<double>(beta_norm2));
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);
    const int remaining_iterations = maxit - (j + 1);
    if (enable_projected_stop && j + 1 >= k &&
        remaining_iterations >= min_projected_savings &&
        ((j + 1 == maxit) || ((j + 1) % check_interval == 0))) {
      int nconv = 0;
      double max_residual = R_PosInf;
      auto projected_timer = native_timer_now();
      const int conv_status = golub_kahan_projected_convergence_estimate(
        alpha, beta, j + 1, k, target_kind, tol, &nconv, &max_residual
      );
      *projected_seconds += native_timer_elapsed(projected_timer);
      ++(*projected_checks);
      if (conv_status != 0) {
        std::free(v);
        std::free(z);
        std::free(u);
        std::free(u_prev);
        std::free(coeff);
        return conv_status;
      }
      *projected_nconv = nconv;
      *projected_max_residual = max_residual;
      if (nconv >= k) {
        *projected_stop = 1;
        break;
      }
    }
    if (j + 1 == maxit || beta[j] <= 100.0 * DBL_EPSILON) {
      break;
    }

    stage_timer = native_timer_now();
    std::memcpy(u_prev, u, sizeof(double) * static_cast<size_t>(m));
    beta_prev = beta[j];
    for (int row = 0; row < n; ++row) {
      v[row] = z[row] / beta[j];
    }
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);
  }

  std::free(v);
  std::free(z);
  std::free(u);
  std::free(u_prev);
  std::free(coeff);
  return 0;
}

extern "C" SEXP eigencore_lanczos_dense(SEXP A_, SEXP maxit_, SEXP start_,
                                        SEXP k_, SEXP target_kind_, SEXP tol_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || LENGTH(start_) != n) {
    error("non-conformable dense Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and nrow(A)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(&impl, eigencore_dense_apply, n, maxit,
                                        k, target_kind, tol,
                                        REAL(start_), REAL(Q_), REAL(alpha_),
                                        REAL(beta_), INTEGER(history_nconv_),
                                        REAL(history_max_residual_),
                                        &iterations, &matvecs);
  if (status != 0) {
    error("native dense Lanczos failed with status=%d", status);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_shift_invert_lanczos_dense(SEXP A_, SEXP sigma_,
                                                     SEXP maxit_, SEXP start_,
                                                     SEXP k_, SEXP target_kind_,
                                                     SEXP tol_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || LENGTH(start_) != n) {
    error("non-conformable dense shift-invert Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and nrow(A)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const double sigma = asReal(sigma_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }
  if (!R_FINITE(sigma)) {
    error("sigma must be finite");
  }

  std::vector<double> lu(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
  std::memcpy(lu.data(), REAL(A_), sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    lu[static_cast<int64_t>(i) + static_cast<int64_t>(i) * n] -= sigma;
  }
  std::vector<int> pivots(static_cast<size_t>(n), 0);
  int info = 0;
  F77_CALL(dgetrf)(&n, &n, lu.data(), &n, pivots.data(), &info);
  if (info < 0) {
    error("LAPACK dgetrf failed for native dense shift-invert with info=%d", info);
  }
  if (info > 0) {
    error("native dense shift-invert factorization is singular at U[%d,%d]; perturb sigma", info, info);
  }

  double min_abs_u = R_PosInf;
  double max_abs_u = 0.0;
  for (int i = 0; i < n; ++i) {
    const double value = fabs(lu[static_cast<int64_t>(i) + static_cast<int64_t>(i) * n]);
    if (value < min_abs_u) {
      min_abs_u = value;
    }
    if (value > max_abs_u) {
      max_abs_u = value;
    }
  }
  const double pivot_ratio = (max_abs_u > 0.0 && R_FINITE(max_abs_u))
    ? min_abs_u / max_abs_u
    : NA_REAL;
  if (R_FINITE(pivot_ratio) && pivot_ratio <= sqrt(DBL_EPSILON)) {
    error("native dense shift-invert factorization is near-singular; perturb sigma");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  std::vector<double> apply_work(static_cast<size_t>(n), 0.0);
  DenseShiftInvertOperator impl = {n, lu.data(), pivots.data(), apply_work.data()};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(
    &impl, eigencore_dense_shift_invert_apply, n, maxit,
    k, target_kind, tol, REAL(start_), REAL(Q_), REAL(alpha_), REAL(beta_),
    INTEGER(history_nconv_), REAL(history_max_residual_),
    &iterations, &matvecs
  );
  if (status != 0) {
    error("native dense shift-invert Lanczos failed with status=%d", status);
  }

  SEXP cache_ = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(cache_, 0, mkString("LAPACK dgetrf/dgetrs"));
  SET_VECTOR_ELT(cache_, 1, ScalarLogical(TRUE));
  SET_VECTOR_ELT(cache_, 2, ScalarReal(pivot_ratio));
  SET_VECTOR_ELT(cache_, 3, ScalarReal(min_abs_u));
  SET_VECTOR_ELT(cache_, 4, ScalarReal(max_abs_u));
  SEXP cache_names_ = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(cache_names_, 0, mkChar("factorization"));
  SET_STRING_ELT(cache_names_, 1, mkChar("factorization_cached"));
  SET_STRING_ELT(cache_names_, 2, mkChar("condition_estimate"));
  SET_STRING_ELT(cache_names_, 3, mkChar("condition_estimate_min_pivot"));
  SET_STRING_ELT(cache_names_, 4, mkChar("condition_estimate_max_pivot"));
  setAttrib(cache_, R_NamesSymbol, cache_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SET_VECTOR_ELT(out_, 7, cache_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  SET_STRING_ELT(names_, 7, mkChar("factorization_cache"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(9);
  return out_;
}

extern "C" SEXP eigencore_shift_invert_lanczos_tridiagonal(
    SEXP lower_, SEXP diag_, SEXP upper_, SEXP maxit_, SEXP start_,
    SEXP k_, SEXP target_kind_, SEXP tol_) {
  if (!isReal(lower_) || !isReal(diag_) || !isReal(upper_) || !isReal(start_)) {
    error("lower, diag, upper, and start must be double");
  }
  const int n = LENGTH(diag_);
  if (n < 1 || LENGTH(start_) != n ||
      LENGTH(lower_) != n - 1 || LENGTH(upper_) != n - 1) {
    error("non-conformable tridiagonal shift-invert Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and length(diag)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }

  const double* lower = REAL(lower_);
  const double* diag = REAL(diag_);
  const double* upper = REAL(upper_);
  std::vector<double> cprime(static_cast<size_t>(n > 1 ? n - 1 : 1), 0.0);
  std::vector<double> denom(static_cast<size_t>(n), 0.0);
  if (fabs(diag[0]) <= DBL_EPSILON) {
    error("native tridiagonal shift-invert factorization has a zero first pivot; perturb sigma");
  }
  denom[0] = diag[0];
  if (n > 1) {
    cprime[0] = upper[0] / denom[0];
  }
  for (int i = 1; i < n; ++i) {
    denom[i] = diag[i] - lower[i - 1] * cprime[i - 1];
    if (fabs(denom[i]) <= DBL_EPSILON) {
      error("native tridiagonal shift-invert factorization has a zero pivot; perturb sigma");
    }
    if (i < n - 1) {
      cprime[i] = upper[i] / denom[i];
    }
  }

  double min_abs_pivot = R_PosInf;
  double max_abs_pivot = 0.0;
  for (int i = 0; i < n; ++i) {
    const double value = fabs(denom[i]);
    if (value < min_abs_pivot) {
      min_abs_pivot = value;
    }
    if (value > max_abs_pivot) {
      max_abs_pivot = value;
    }
  }
  const double pivot_ratio = (max_abs_pivot > 0.0 && R_FINITE(max_abs_pivot))
    ? min_abs_pivot / max_abs_pivot
    : NA_REAL;
  if (R_FINITE(pivot_ratio) && pivot_ratio <= sqrt(DBL_EPSILON)) {
    error("native tridiagonal shift-invert factorization is near-singular; perturb sigma");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  std::vector<double> apply_work(static_cast<size_t>(n), 0.0);
  TridiagonalShiftInvertOperator impl = {
    n, lower, cprime.data(), denom.data(), apply_work.data()
  };
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(
    &impl, eigencore_tridiagonal_shift_invert_apply, n, maxit,
    k, target_kind, tol, REAL(start_), REAL(Q_), REAL(alpha_), REAL(beta_),
    INTEGER(history_nconv_), REAL(history_max_residual_),
    &iterations, &matvecs
  );
  if (status != 0) {
    error("native tridiagonal shift-invert Lanczos failed with status=%d", status);
  }

  SEXP cache_ = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(cache_, 0, mkString("native tridiagonal Thomas"));
  SET_VECTOR_ELT(cache_, 1, ScalarLogical(TRUE));
  SET_VECTOR_ELT(cache_, 2, ScalarReal(pivot_ratio));
  SET_VECTOR_ELT(cache_, 3, ScalarReal(min_abs_pivot));
  SET_VECTOR_ELT(cache_, 4, ScalarReal(max_abs_pivot));
  SEXP cache_names_ = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(cache_names_, 0, mkChar("factorization"));
  SET_STRING_ELT(cache_names_, 1, mkChar("factorization_cached"));
  SET_STRING_ELT(cache_names_, 2, mkChar("condition_estimate"));
  SET_STRING_ELT(cache_names_, 3, mkChar("condition_estimate_min_pivot"));
  SET_STRING_ELT(cache_names_, 4, mkChar("condition_estimate_max_pivot"));
  setAttrib(cache_, R_NamesSymbol, cache_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SET_VECTOR_ELT(out_, 7, cache_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  SET_STRING_ELT(names_, 7, mkChar("factorization_cache"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(9);
  return out_;
}

extern "C" SEXP eigencore_shift_invert_lanczos_tridiagonal_generalized(
    SEXP lower_, SEXP diag_, SEXP upper_, SEXP sqrt_metric_, SEXP maxit_,
    SEXP start_, SEXP k_, SEXP target_kind_, SEXP tol_) {
  if (!isReal(lower_) || !isReal(diag_) || !isReal(upper_) ||
      !isReal(sqrt_metric_) || !isReal(start_)) {
    error("lower, diag, upper, sqrt_metric, and start must be double");
  }
  const int n = LENGTH(diag_);
  if (n < 1 || LENGTH(start_) != n || LENGTH(sqrt_metric_) != n ||
      LENGTH(lower_) != n - 1 || LENGTH(upper_) != n - 1) {
    error("non-conformable generalized tridiagonal shift-invert Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and length(diag)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }

  const double* lower = REAL(lower_);
  const double* diag = REAL(diag_);
  const double* upper = REAL(upper_);
  const double* sqrt_metric = REAL(sqrt_metric_);
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(sqrt_metric[i]) || sqrt_metric[i] <= 0.0) {
      error("generalized tridiagonal shift-invert requires positive finite diagonal B");
    }
  }

  std::vector<double> cprime(static_cast<size_t>(n > 1 ? n - 1 : 1), 0.0);
  std::vector<double> denom(static_cast<size_t>(n), 0.0);
  if (fabs(diag[0]) <= DBL_EPSILON) {
    error("native generalized tridiagonal shift-invert factorization has a zero first pivot; perturb sigma");
  }
  denom[0] = diag[0];
  if (n > 1) {
    cprime[0] = upper[0] / denom[0];
  }
  for (int i = 1; i < n; ++i) {
    denom[i] = diag[i] - lower[i - 1] * cprime[i - 1];
    if (fabs(denom[i]) <= DBL_EPSILON) {
      error("native generalized tridiagonal shift-invert factorization has a zero pivot; perturb sigma");
    }
    if (i < n - 1) {
      cprime[i] = upper[i] / denom[i];
    }
  }

  double min_abs_pivot = R_PosInf;
  double max_abs_pivot = 0.0;
  for (int i = 0; i < n; ++i) {
    const double value = fabs(denom[i]);
    if (value < min_abs_pivot) {
      min_abs_pivot = value;
    }
    if (value > max_abs_pivot) {
      max_abs_pivot = value;
    }
  }
  const double pivot_ratio = (max_abs_pivot > 0.0 && R_FINITE(max_abs_pivot))
    ? min_abs_pivot / max_abs_pivot
    : NA_REAL;
  if (R_FINITE(pivot_ratio) && pivot_ratio <= sqrt(DBL_EPSILON)) {
    error("native generalized tridiagonal shift-invert factorization is near-singular; perturb sigma");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  std::vector<double> apply_work(static_cast<size_t>(n), 0.0);
  TridiagonalGeneralizedShiftInvertOperator impl = {
    n, lower, cprime.data(), denom.data(), sqrt_metric, apply_work.data()
  };
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(
    &impl, eigencore_tridiagonal_generalized_shift_invert_apply, n, maxit,
    k, target_kind, tol, REAL(start_), REAL(Q_), REAL(alpha_), REAL(beta_),
    INTEGER(history_nconv_), REAL(history_max_residual_),
    &iterations, &matvecs
  );
  if (status != 0) {
    error("native generalized tridiagonal shift-invert Lanczos failed with status=%d", status);
  }

  SEXP cache_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(cache_, 0, mkString("native tridiagonal Thomas + diagonal sqrt(B)"));
  SET_VECTOR_ELT(cache_, 1, ScalarLogical(TRUE));
  SET_VECTOR_ELT(cache_, 2, ScalarReal(pivot_ratio));
  SET_VECTOR_ELT(cache_, 3, ScalarReal(min_abs_pivot));
  SET_VECTOR_ELT(cache_, 4, ScalarReal(max_abs_pivot));
  SET_VECTOR_ELT(cache_, 5, mkString("diagonal sqrt(B)"));
  SEXP cache_names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(cache_names_, 0, mkChar("factorization"));
  SET_STRING_ELT(cache_names_, 1, mkChar("factorization_cached"));
  SET_STRING_ELT(cache_names_, 2, mkChar("condition_estimate"));
  SET_STRING_ELT(cache_names_, 3, mkChar("condition_estimate_min_pivot"));
  SET_STRING_ELT(cache_names_, 4, mkChar("condition_estimate_max_pivot"));
  SET_STRING_ELT(cache_names_, 5, mkChar("metric_factorization"));
  setAttrib(cache_, R_NamesSymbol, cache_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SET_VECTOR_ELT(out_, 7, cache_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  SET_STRING_ELT(names_, 7, mkChar("factorization_cache"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(9);
  return out_;
}

extern "C" SEXP eigencore_shift_invert_lanczos_dense_generalized(
    SEXP A_, SEXP B_, SEXP sigma_, SEXP maxit_, SEXP start_,
    SEXP k_, SEXP target_kind_, SEXP tol_) {
  if (!isReal(A_) || !isReal(B_) || !isReal(start_)) {
    error("A, B, and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n ||
      LENGTH(start_) != n) {
    error("non-conformable dense generalized shift-invert Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and nrow(A)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const double sigma = asReal(sigma_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }
  if (!R_FINITE(sigma)) {
    error("sigma must be finite");
  }

  std::vector<double> chol(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
  std::memcpy(chol.data(), REAL(B_), sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
  char uplo = 'U';
  int info = 0;
  F77_CALL(dpotrf)(&uplo, &n, chol.data(), &n, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed for native generalized shift-invert B with info=%d", info);
  }
  for (int col = 0; col < n; ++col) {
    for (int row = col + 1; row < n; ++row) {
      chol[static_cast<int64_t>(row) + static_cast<int64_t>(col) * n] = 0.0;
    }
  }

  std::vector<double> lu(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
  for (int64_t pos = 0; pos < static_cast<int64_t>(n) * n; ++pos) {
    lu[static_cast<size_t>(pos)] = REAL(A_)[pos] - sigma * REAL(B_)[pos];
  }
  std::vector<int> pivots(static_cast<size_t>(n), 0);
  F77_CALL(dgetrf)(&n, &n, lu.data(), &n, pivots.data(), &info);
  if (info < 0) {
    error("LAPACK dgetrf failed for native dense generalized shift-invert with info=%d", info);
  }
  if (info > 0) {
    error("native dense generalized shift-invert factorization is singular at U[%d,%d]; perturb sigma", info, info);
  }

  double min_abs_u = R_PosInf;
  double max_abs_u = 0.0;
  for (int i = 0; i < n; ++i) {
    const double value = fabs(lu[static_cast<int64_t>(i) + static_cast<int64_t>(i) * n]);
    if (value < min_abs_u) {
      min_abs_u = value;
    }
    if (value > max_abs_u) {
      max_abs_u = value;
    }
  }
  const double pivot_ratio = (max_abs_u > 0.0 && R_FINITE(max_abs_u))
    ? min_abs_u / max_abs_u
    : NA_REAL;
  if (R_FINITE(pivot_ratio) && pivot_ratio <= sqrt(DBL_EPSILON)) {
    error("native dense generalized shift-invert factorization is near-singular; perturb sigma");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  std::vector<double> rhs(static_cast<size_t>(n), 0.0);
  std::vector<double> sol(static_cast<size_t>(n), 0.0);
  DenseGeneralizedShiftInvertOperator impl = {
    n, lu.data(), pivots.data(), chol.data(), rhs.data(), sol.data()
  };
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(
    &impl, eigencore_dense_generalized_shift_invert_apply, n, maxit,
    k, target_kind, tol, REAL(start_), REAL(Q_), REAL(alpha_), REAL(beta_),
    INTEGER(history_nconv_), REAL(history_max_residual_),
    &iterations, &matvecs
  );
  if (status != 0) {
    error("native dense generalized shift-invert Lanczos failed with status=%d", status);
  }

  SEXP cache_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(cache_, 0, mkString("LAPACK dpotrf(B) + dgetrf/dgetrs(A - sigma B)"));
  SET_VECTOR_ELT(cache_, 1, ScalarLogical(TRUE));
  SET_VECTOR_ELT(cache_, 2, ScalarReal(pivot_ratio));
  SET_VECTOR_ELT(cache_, 3, ScalarReal(min_abs_u));
  SET_VECTOR_ELT(cache_, 4, ScalarReal(max_abs_u));
  SET_VECTOR_ELT(cache_, 5, mkString("LAPACK dpotrf(B)"));
  SEXP cache_names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(cache_names_, 0, mkChar("factorization"));
  SET_STRING_ELT(cache_names_, 1, mkChar("factorization_cached"));
  SET_STRING_ELT(cache_names_, 2, mkChar("condition_estimate"));
  SET_STRING_ELT(cache_names_, 3, mkChar("condition_estimate_min_pivot"));
  SET_STRING_ELT(cache_names_, 4, mkChar("condition_estimate_max_pivot"));
  SET_STRING_ELT(cache_names_, 5, mkChar("metric_factorization"));
  setAttrib(cache_, R_NamesSymbol, cache_names_);

  SEXP chol_ = PROTECT(allocMatrix(REALSXP, n, n));
  std::memcpy(REAL(chol_), chol.data(), sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SET_VECTOR_ELT(out_, 7, cache_);
  SET_VECTOR_ELT(out_, 8, chol_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  SET_STRING_ELT(names_, 7, mkChar("factorization_cache"));
  SET_STRING_ELT(names_, 8, mkChar("chol_factor"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(10);
  return out_;
}

extern "C" SEXP eigencore_lanczos_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                      SEXP maxit_, SEXP start_, SEXP k_,
                                      SEXP target_kind_, SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC Lanczos inputs");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || LENGTH(start_) != n) {
    error("non-conformable CSC Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and nrow(A)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(&impl, eigencore_csc_apply, n, maxit,
                                        k, target_kind, tol,
                                        REAL(start_), REAL(Q_), REAL(alpha_),
                                        REAL(beta_), INTEGER(history_nconv_),
                                        REAL(history_max_residual_),
                                        &iterations, &matvecs);
  if (status != 0) {
    error("native CSC Lanczos failed with status=%d", status);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_dense(SEXP A_, SEXP maxit_, SEXP start_,
                                            SEXP rank_, SEXP target_kind_,
                                            SEXP tol_, SEXP projected_stop_) {
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
                                            1, 1);
  if (status != 0) {
    error("native dense Golub-Kahan failed with status=%d", status);
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, iterations));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, iterations));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP beta_ = PROTECT(allocVector(REALSXP, iterations));
  for (int col = 0; col < iterations; ++col) {
    std::memcpy(REAL(U_) + static_cast<size_t>(col) * static_cast<size_t>(m),
                U_work.data() + static_cast<size_t>(col) * static_cast<size_t>(m),
                sizeof(double) * static_cast<size_t>(m));
    std::memcpy(REAL(V_) + static_cast<size_t>(col) * static_cast<size_t>(n),
                V_work.data() + static_cast<size_t>(col) * static_cast<size_t>(n),
                sizeof(double) * static_cast<size_t>(n));
  }
  std::memcpy(REAL(alpha_), alpha_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));
  std::memcpy(REAL(beta_), beta_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 8, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 10, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalization_passes));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 7, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 8, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 9, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 10, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 12, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalization_passes"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                          SEXP maxit_, SEXP start_,
                                          SEXP rank_, SEXP target_kind_,
                                          SEXP tol_, SEXP projected_stop_) {
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
                                            1, 1);
  if (status != 0) {
    error("native CSC Golub-Kahan failed with status=%d", status);
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, iterations));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, iterations));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP beta_ = PROTECT(allocVector(REALSXP, iterations));
  for (int col = 0; col < iterations; ++col) {
    std::memcpy(REAL(U_) + static_cast<size_t>(col) * static_cast<size_t>(m),
                U_work.data() + static_cast<size_t>(col) * static_cast<size_t>(m),
                sizeof(double) * static_cast<size_t>(m));
    std::memcpy(REAL(V_) + static_cast<size_t>(col) * static_cast<size_t>(n),
                V_work.data() + static_cast<size_t>(col) * static_cast<size_t>(n),
                sizeof(double) * static_cast<size_t>(n));
  }
  std::memcpy(REAL(alpha_), alpha_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));
  std::memcpy(REAL(beta_), beta_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 8, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 10, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalization_passes));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 7, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 8, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 9, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 10, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 12, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalization_passes"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_r_operator(
    SEXP dim_, SEXP apply_, SEXP apply_adjoint_, SEXP maxit_, SEXP start_,
    SEXP rank_, SEXP target_kind_, SEXP tol_, SEXP projected_stop_,
    SEXP reorthogonalize_u_, SEXP reorthogonalize_v_) {
  if (!isInteger(dim_) || LENGTH(dim_) != 2 || TYPEOF(apply_) != CLOSXP ||
      TYPEOF(apply_adjoint_) != CLOSXP || !isReal(start_)) {
    error("invalid matrix-free Golub-Kahan inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  if (m < 1 || n < 1 || LENGTH(start_) != n) {
    error("non-conformable matrix-free Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  const int reorthogonalize_u = asLogical(reorthogonalize_u_) == TRUE;
  const int reorthogonalize_v = asLogical(reorthogonalize_v_) == TRUE;
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

  RApplyOperator impl = {m, n, apply_, apply_adjoint_};
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
    &impl, eigencore_r_operator_apply, m, n, maxit, rank, target_kind, tol,
    enable_projected_stop, 0, REAL(start_), U_work.data(), V_work.data(),
    alpha_work.data(), beta_work.data(), &iterations, &matvecs,
    &projected_stop, &projected_nconv, &projected_max_residual,
    &projected_checks, &projected_seconds, &stage_apply_seconds,
    &stage_recurrence_seconds, &stage_reorthogonalization_seconds,
    &reorthogonalization_passes, reorthogonalize_u, reorthogonalize_v
  );
  if (status != 0) {
    error("native matrix-free Golub-Kahan failed with status=%d", status);
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, iterations));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, iterations));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP beta_ = PROTECT(allocVector(REALSXP, iterations));
  for (int col = 0; col < iterations; ++col) {
    std::memcpy(REAL(U_) + static_cast<size_t>(col) * static_cast<size_t>(m),
                U_work.data() + static_cast<size_t>(col) * static_cast<size_t>(m),
                sizeof(double) * static_cast<size_t>(m));
    std::memcpy(REAL(V_) + static_cast<size_t>(col) * static_cast<size_t>(n),
                V_work.data() + static_cast<size_t>(col) * static_cast<size_t>(n),
                sizeof(double) * static_cast<size_t>(n));
  }
  std::memcpy(REAL(alpha_), alpha_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));
  std::memcpy(REAL(beta_), beta_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 8, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 10, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalization_passes));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 7, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 8, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 9, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 10, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 12, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalization_passes"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}
