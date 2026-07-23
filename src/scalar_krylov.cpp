#include <cmath>
#include <cfloat>
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

// Reusable scratch for the Lanczos convergence estimate. The estimate runs
// every iteration once the subspace holds k candidates, so the buffers are
// sized once for the sweep's maximum subspace (maxit) and the k requested
// vectors instead of being reallocated per check.
struct LanczosConvergenceScratch {
  double* w = nullptr;         // capacity: full eigenvalue list (dsterf)
  double* e_work = nullptr;    // capacity: offdiagonal copy consumed by dsterf
  double* d2 = nullptr;        // capacity: diagonal copy consumed by dstevr
  double* e2 = nullptr;        // capacity: offdiagonal copy consumed by dstevr
  double* w2 = nullptr;        // capacity: eigenvalues returned by dstevr
  double* z = nullptr;         // capacity * kcap: selected eigenvectors
  double* work = nullptr;      // 20 * capacity doubles: dstevr workspace
  int* iwork = nullptr;        // 10 * capacity ints: dstevr workspace
  int* isuppz = nullptr;       // 2 * kcap ints
  int* selected = nullptr;     // capacity ints
  int capacity = 0;
  int kcap = 0;

  int ensure(int needed, int k_needed) {
    if (needed <= capacity && k_needed <= kcap) {
      return 0;
    }
    release();
    const size_t cap = static_cast<size_t>(needed);
    const size_t kc = static_cast<size_t>(k_needed);
    w = static_cast<double*>(std::malloc(cap * sizeof(double)));
    e_work = static_cast<double*>(std::malloc(cap * sizeof(double)));
    d2 = static_cast<double*>(std::malloc(cap * sizeof(double)));
    e2 = static_cast<double*>(std::malloc(cap * sizeof(double)));
    w2 = static_cast<double*>(std::malloc(cap * sizeof(double)));
    z = static_cast<double*>(std::malloc(cap * kc * sizeof(double)));
    work = static_cast<double*>(std::malloc(20U * cap * sizeof(double)));
    iwork = static_cast<int*>(std::malloc(10U * cap * sizeof(int)));
    isuppz = static_cast<int*>(std::malloc(2U * kc * sizeof(int)));
    selected = static_cast<int*>(std::malloc(cap * sizeof(int)));
    if (w == nullptr || e_work == nullptr || d2 == nullptr || e2 == nullptr ||
        w2 == nullptr || z == nullptr || work == nullptr || iwork == nullptr ||
        isuppz == nullptr || selected == nullptr) {
      release();
      return -2;
    }
    capacity = needed;
    kcap = k_needed;
    return 0;
  }

  void release() {
    std::free(w);
    std::free(e_work);
    std::free(d2);
    std::free(e2);
    std::free(w2);
    std::free(z);
    std::free(work);
    std::free(iwork);
    std::free(isuppz);
    std::free(selected);
    w = nullptr;
    e_work = nullptr;
    d2 = nullptr;
    e2 = nullptr;
    w2 = nullptr;
    z = nullptr;
    work = nullptr;
    iwork = nullptr;
    isuppz = nullptr;
    selected = nullptr;
    capacity = 0;
    kcap = 0;
  }
};

// Ritz-residual convergence estimate for the Lanczos tridiagonal. The bound
// |beta_last * z[last, idx]| only needs the selected Ritz values and the
// LAST components of their tridiagonal eigenvectors, so instead of dstev
// with jobz='V' (all eigenvectors, O(iter^3) per check) we take eigenvalues
// from dsterf (O(iter^2)) and recover just the selected eigenvectors with
// dstevr over the contiguous index runs of the selection (the magnitude
// targets select from both spectrum ends, giving at most two runs).
static int lanczos_convergence_estimate(const double* alpha,
                                        const double* beta,
                                        int iter,
                                        int k,
                                        int target_kind,
                                        double tol,
                                        int* nconv,
                                        double* max_residual,
                                        LanczosConvergenceScratch* scratch) {
  *nconv = 0;
  *max_residual = R_PosInf;
  if (iter <= 0 || k <= 0) {
    return 0;
  }

  LanczosConvergenceScratch local;
  LanczosConvergenceScratch* s = (scratch != nullptr) ? scratch : &local;
  const int k_needed = (k < iter) ? k : iter;
  const int ensure_status = s->ensure(iter, k_needed);
  if (ensure_status != 0) {
    return ensure_status;
  }

  double max_selected = 0.0;
  int conv = 0;

  if (iter == 1) {
    const double residual = fabs(beta[0]);
    const double threshold = tol * ((fabs(alpha[0]) > 1.0) ? fabs(alpha[0]) : 1.0);
    if (residual <= threshold) {
      ++conv;
    }
    max_selected = residual;
    *nconv = conv;
    *max_residual = max_selected;
    if (s == &local) {
      local.release();
    }
    return 0;
  }

  std::memcpy(s->w, alpha, sizeof(double) * static_cast<size_t>(iter));
  std::memcpy(s->e_work, beta, sizeof(double) * static_cast<size_t>(iter - 1));
  int info = 0;
  F77_CALL(dsterf)(&iter, s->w, s->e_work, &info);
  if (info != 0) {
    if (s == &local) {
      local.release();
    }
    return info;
  }

  const int count = selected_ritz_indices(s->w, iter, k, target_kind, s->selected);
  // Sort the selected indices so contiguous runs can be recovered with a
  // single dstevr range each; the selection order itself is irrelevant to
  // the aggregate nconv / max_residual outputs.
  for (int i = 1; i < count; ++i) {
    const int key = s->selected[i];
    int pos = i - 1;
    while (pos >= 0 && s->selected[pos] > key) {
      s->selected[pos + 1] = s->selected[pos];
      --pos;
    }
    s->selected[pos + 1] = key;
  }

  int run_start = 0;
  while (run_start < count) {
    int run_end = run_start;
    while (run_end + 1 < count &&
           s->selected[run_end + 1] == s->selected[run_end] + 1) {
      ++run_end;
    }
    const int il = s->selected[run_start] + 1;
    const int iu = s->selected[run_end] + 1;
    const int run_len = iu - il + 1;

    std::memcpy(s->d2, alpha, sizeof(double) * static_cast<size_t>(iter));
    std::memcpy(s->e2, beta, sizeof(double) * static_cast<size_t>(iter - 1));
    char jobz = 'V';
    char range = 'I';
    double vl = 0.0;
    double vu = 0.0;
    double abstol = 0.0;
    int found = 0;
    int lwork = 20 * iter;
    int liwork = 10 * iter;
    int iter_local = iter;
    F77_CALL(dstevr)(&jobz, &range, &iter_local, s->d2, s->e2, &vl, &vu,
                     const_cast<int*>(&il), const_cast<int*>(&iu), &abstol,
                     &found, s->w2, s->z, &iter_local, s->isuppz,
                     s->work, &lwork, s->iwork, &liwork, &info FCONE FCONE);
    if (info != 0 || found != run_len) {
      if (s == &local) {
        local.release();
      }
      return (info != 0) ? info : -3;
    }

    for (int i = run_start; i <= run_end; ++i) {
      const int idx = s->selected[i];
      const int z_col = idx - (il - 1);
      const double residual =
        fabs(beta[iter - 1] * s->z[(iter - 1) + static_cast<int64_t>(z_col) * iter]);
      const double value = s->w[idx];
      const double threshold = tol * ((fabs(value) > 1.0) ? fabs(value) : 1.0);
      if (residual <= threshold) {
        ++conv;
      }
      if (residual > max_selected) {
        max_selected = residual;
      }
    }
    run_start = run_end + 1;
  }

  *nconv = conv;
  *max_residual = max_selected;
  if (s == &local) {
    local.release();
  }
  return 0;
}

// Reusable scratch for the projected convergence check so repeated checks
// inside one Golub-Kahan sweep do not reallocate. Capacity is the maximum
// subspace size (maxit) of the sweep.
struct GolubKahanProjectedScratch {
  double* d = nullptr;         // capacity doubles: singular values
  double* e = nullptr;         // capacity doubles: superdiagonal copy
  double* u_last = nullptr;    // capacity doubles: last row of projected U
  double* work = nullptr;      // 4 * capacity doubles: dbdsqr workspace
  int* selected = nullptr;     // capacity ints
  int capacity = 0;

  int ensure(int needed) {
    if (needed <= capacity) {
      return 0;
    }
    release();
    d = static_cast<double*>(std::malloc(static_cast<size_t>(needed) * sizeof(double)));
    e = static_cast<double*>(std::malloc(static_cast<size_t>(needed) * sizeof(double)));
    u_last = static_cast<double*>(std::malloc(static_cast<size_t>(needed) * sizeof(double)));
    work = static_cast<double*>(std::malloc(4U * static_cast<size_t>(needed) * sizeof(double)));
    selected = static_cast<int*>(std::malloc(static_cast<size_t>(needed) * sizeof(int)));
    if (d == nullptr || e == nullptr || u_last == nullptr ||
        work == nullptr || selected == nullptr) {
      release();
      return -2;
    }
    capacity = needed;
    return 0;
  }

  void release() {
    std::free(d);
    std::free(e);
    std::free(u_last);
    std::free(work);
    std::free(selected);
    d = nullptr;
    e = nullptr;
    u_last = nullptr;
    work = nullptr;
    selected = nullptr;
    capacity = 0;
  }
};

// Projected convergence estimate for the Golub-Kahan bidiagonalization.
// The residual bound only needs the singular values of the projected
// bidiagonal B and the LAST row of its left singular vectors, so instead of
// a dense dgesvd with jobu='A' (O(iter^3) plus per-check allocation churn)
// we run dbdsqr directly on the bidiagonal, rotating a single row vector
// seeded with e_iter^T. dbdsqr is the same kernel dgesvd applies after
// bidiagonal reduction, so singular values and the tracked row agree with
// the previous implementation to rounding.
static int golub_kahan_projected_convergence_estimate(const double* alpha,
                                                      const double* beta,
                                                      int iter,
                                                      int k,
                                                      int target_kind,
                                                      double tol,
                                                      int* nconv,
                                                      double* max_residual,
                                                      GolubKahanProjectedScratch* scratch) {
  *nconv = 0;
  *max_residual = R_PosInf;
  if (iter <= 0 || k <= 0) {
    return 0;
  }

  GolubKahanProjectedScratch local;
  GolubKahanProjectedScratch* s = (scratch != nullptr) ? scratch : &local;
  const int ensure_status = s->ensure(iter);
  if (ensure_status != 0) {
    return ensure_status;
  }

  std::memcpy(s->d, alpha, sizeof(double) * static_cast<size_t>(iter));
  if (iter > 1) {
    std::memcpy(s->e, beta, sizeof(double) * static_cast<size_t>(iter - 1));
  }
  // Track only the last row of the projected left singular vectors: dbdsqr
  // applies its rotations to any caller-supplied U block, so a 1 x iter row
  // seeded with e_iter^T ends up holding that row exactly.
  std::memset(s->u_last, 0, sizeof(double) * static_cast<size_t>(iter));
  s->u_last[iter - 1] = 1.0;

  char uplo = 'U';
  const int ncvt = 0;
  const int nru = 1;
  const int ncc = 0;
  const int ldu = 1;
  int info = 0;
  double dummy = 0.0;
  F77_CALL(dbdsqr)(&uplo, &iter, &ncvt, &nru, &ncc,
                   s->d, s->e, &dummy, &iter, s->u_last, &ldu,
                   &dummy, &iter, s->work, &info FCONE);
  if (info != 0) {
    if (s == &local) {
      local.release();
    }
    return info;
  }

  const int count = selected_ritz_indices(s->d, iter, k, target_kind, s->selected);
  double max_selected = 0.0;
  int conv = 0;
  for (int i = 0; i < count; ++i) {
    const int idx = s->selected[i];
    const double residual = fabs(beta[iter - 1] * s->u_last[idx]);
    const double scale = (fabs(s->d[idx]) > 1.0) ? fabs(s->d[idx]) : 1.0;
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
  if (s == &local) {
    local.release();
  }
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
  LanczosConvergenceScratch convergence_scratch;
  // Size the convergence scratch for the full sweep once so the
  // per-iteration estimate never reallocates.
  if (convergence_scratch.ensure(maxit, (k < maxit) ? k : maxit) != 0) {
    std::free(q);
    std::free(q_prev);
    std::free(z);
    std::free(coeff);
    return -2;
  }
  for (int j = 0; j < maxit; ++j) {
    *iterations = j + 1;
    std::memcpy(Q + j * n, q, sizeof(double) * static_cast<size_t>(n));
    std::memset(z, 0, sizeof(double) * static_cast<size_t>(n));

    const int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, q, n,
                             1.0, 0.0, z, n, &workspace);
    if (status != 0) {
      convergence_scratch.release();
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
        alpha, beta, j + 1, k, target_kind, tol, &nconv, &max_residual,
        &convergence_scratch
      );
      if (conv_status != 0) {
        convergence_scratch.release();
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

  convergence_scratch.release();
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
  GolubKahanProjectedScratch projected_scratch;
  // Size the projected-stop scratch for the full sweep once so repeated
  // convergence checks never reallocate.
  if (enable_projected_stop && projected_scratch.ensure(maxit) != 0) {
    std::free(v);
    std::free(z);
    std::free(u);
    std::free(u_prev);
    std::free(coeff);
    return -2;
  }
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
      projected_scratch.release();
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
      projected_scratch.release();
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
        alpha, beta, j + 1, k, target_kind, tol, &nconv, &max_residual,
        &projected_scratch
      );
      *projected_seconds += native_timer_elapsed(projected_timer);
      ++(*projected_checks);
      if (conv_status != 0) {
        projected_scratch.release();
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

  projected_scratch.release();
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
    enable_projected_stop, 1, REAL(start_), U_work.data(), V_work.data(),
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
