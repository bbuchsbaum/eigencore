#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstring>
#include <vector>
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

static int selected_sorted_ritz_indices(const double* values,
                                        int n,
                                        int k,
                                        int target_kind,
                                        int* selected) {
  const int count = (k < n) ? k : n;
  if (target_kind == 1) {
    for (int i = 0; i < count; ++i) {
      selected[i] = n - 1 - i;
    }
    return count;
  }
  if (target_kind == 2) {
    for (int i = 0; i < count; ++i) {
      selected[i] = i;
    }
    return count;
  }
  return selected_ritz_indices(values, n, k, target_kind, selected);
}

// =====================================================================
// Thick-restart Hermitian Lanczos with locking
//
// Implements a Krylov-Schur-style restarted Lanczos for the symmetric /
// Hermitian standard eigenproblem A x = lambda x.
// =====================================================================

static double trl_norm2(const double* x, int n);

// DGKS reorthogonalization with an adaptive second pass: the second
// projection runs only when the first one cancelled a large fraction of the
// vector norm (post < eta * pre, eta = 1/sqrt(2); Daniel-Gragg-Kaufman-
// Stewart). With that criterion two passes carry the same orthogonality
// guarantee as unconditional CGS2 ("twice is enough"). Returns the number of
// projection passes actually performed.
static const double kDgksEta = 0.7071067811865475;

static int trl_orthogonalise(const double* V_locked, int n_locked,
                             const double* V_active, int m_active,
                             double* z, double* tmp, int n,
                             int max_passes = 2) {
  if (n_locked <= 0 && m_active <= 0) {
    return 0;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  int incx = 1;
  int passes_done = 0;
  double pre_norm = trl_norm2(z, n);
  for (int pass = 0; pass < max_passes; ++pass) {
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
    ++passes_done;
    if (pass + 1 >= max_passes) {
      break;
    }
    const double post_norm = trl_norm2(z, n);
    if (post_norm >= kDgksEta * pre_norm) {
      break;
    }
    pre_norm = post_norm;
  }
  return passes_done;
}

static double trl_norm2(const double* x, int n) {
  long double sum = 0.0L;
  for (int i = 0; i < n; ++i) {
    sum += static_cast<long double>(x[i]) * x[i];
  }
  return sqrt(static_cast<double>(sum));
}

// Frobenius norm of a dense n x cols block stored contiguously (ld == n).
static double block_frobenius_norm(const double* X, int n, int cols) {
  long double sum = 0.0L;
  const int64_t total = static_cast<int64_t>(n) * cols;
  for (int64_t i = 0; i < total; ++i) {
    sum += static_cast<long double>(X[i]) * X[i];
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

static int trl_dsyevd_query(int m_max, int* liwork_out) {
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork_query = -1;
  int liwork_query = -1;
  double work_query = 0.0;
  int iwork_query = 0;
  double fake = 0.0;
  double fake_w = 0.0;
  int m = m_max;
  F77_CALL(dsyevd)(&jobz, &uplo, &m, &fake, &m, &fake_w,
                   &work_query, &lwork_query,
                   &iwork_query, &liwork_query, &info FCONE FCONE);
  if (info != 0) {
    if (liwork_out != nullptr) {
      *liwork_out = 3 + 5 * m_max;
    }
    return 1 + 6 * m_max + 2 * m_max * m_max;
  }
  if (liwork_out != nullptr) {
    *liwork_out = iwork_query;
  }
  return static_cast<int>(work_query);
}

static int symmetric_eigen_inplace(double* A, int n, double* values,
                                   double* work, int lwork,
                                   int* iwork = nullptr,
                                   int liwork = 0) {
  if (n <= 0) {
    return 0;
  }
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  if (n >= 96 && iwork != nullptr && liwork > 0) {
    F77_CALL(dsyevd)(&jobz, &uplo, &n, A, &n, values,
                     work, &lwork, iwork, &liwork, &info FCONE FCONE);
  } else {
    F77_CALL(dsyev)(&jobz, &uplo, &n, A, &n, values,
                    work, &lwork, &info FCONE FCONE);
  }
  return info == 0 ? 0 : -3;
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

static double standard_eigen_lock_scale(double norm_a, double theta,
                                        const double* v, int n) {
  if (!std::isfinite(norm_a) || norm_a <= 0.0) {
    norm_a = 1.0;
  }
  const double vnorm = trl_norm2(v, n);
  const double scale = (norm_a + fabs(theta)) *
    ((vnorm > DBL_EPSILON) ? vnorm : DBL_EPSILON);
  return (scale > DBL_EPSILON) ? scale : DBL_EPSILON;
}

static int vector_is_independent_from_locked(const double* V_locked,
                                             int n_locked,
                                             const double* v,
                                             int n) {
  const double dot_tol = 10.0 * sqrt(DBL_EPSILON);
  const double vnorm = trl_norm2(v, n);
  if (vnorm <= DBL_EPSILON) {
    return 0;
  }
  for (int col = 0; col < n_locked; ++col) {
    long double dot = 0.0L;
    long double locked_ss = 0.0L;
    const double* locked = V_locked + static_cast<int64_t>(col) * n;
    for (int row = 0; row < n; ++row) {
      dot += static_cast<long double>(locked[row]) * v[row];
      locked_ss += static_cast<long double>(locked[row]) * locked[row];
    }
    const double locked_norm = sqrt(static_cast<double>(locked_ss));
    if (locked_norm <= DBL_EPSILON) {
      continue;
    }
    if (fabs(static_cast<double>(dot)) > dot_tol * locked_norm * vnorm) {
      return 0;
    }
  }
  return 1;
}

static int block_accept_work_vector(const double* V_locked, int n_locked,
                                    double* V_active, int* m_active,
                                    int m_max, double* z, double* tmp,
                                    int n, int* ortho_passes) {
  if (*m_active >= m_max) {
    return 0;
  }
  const int passes_done =
    trl_orthogonalise(V_locked, n_locked, V_active, *m_active, z, tmp, n, 2);
  if (ortho_passes != nullptr) {
    *ortho_passes += passes_done;
  }
  const double nz = trl_norm2(z, n);
  if (nz <= 100.0 * DBL_EPSILON) {
    return 0;
  }
  const double inv_nz = 1.0 / nz;
  double* dst = V_active + static_cast<int64_t>(*m_active) * n;
  for (int row = 0; row < n; ++row) {
    dst[row] = z[row] * inv_nz;
  }
  ++(*m_active);
  return 1;
}

// Block variant of the adaptive DGKS scheme in trl_orthogonalise: the second
// projection runs only when the first cancelled a large fraction of the block
// Frobenius norm. Returns the number of projection passes performed.
static int block_reorthogonalise_against(const double* V_locked, int n_locked,
                                         const double* V_active, int m_active,
                                         double* X, int n, int cols,
                                         double* coeff, int max_passes) {
  if ((n_locked <= 0 && m_active <= 0) || cols <= 0) {
    return 0;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  int passes_done = 0;
  double pre_norm = block_frobenius_norm(X, n, cols);
  for (int pass = 0; pass < max_passes; ++pass) {
    if (n_locked > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &n_locked, &cols, &n,
                      &one, V_locked, &n, X, &n,
                      &zero, coeff, &n_locked FCONE FCONE);
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &cols, &n_locked,
                      &minus_one, V_locked, &n, coeff, &n_locked,
                      &one, X, &n FCONE FCONE);
    }
    if (m_active > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &m_active, &cols, &n,
                      &one, V_active, &n, X, &n,
                      &zero, coeff, &m_active FCONE FCONE);
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &cols, &m_active,
                      &minus_one, V_active, &n, coeff, &m_active,
                      &one, X, &n FCONE FCONE);
    }
    ++passes_done;
    if (pass + 1 >= max_passes) {
      break;
    }
    const double post_norm = block_frobenius_norm(X, n, cols);
    if (post_norm >= kDgksEta * pre_norm) {
      break;
    }
    pre_norm = post_norm;
  }
  return passes_done;
}

// Contract: X and Z_block MAY alias (X == Z_block && ldx == n is a supported
// caller pattern — the per-column memcpy below would otherwise be a self-copy,
// which is UB in standard C++ even when the bytes happen to overlap exactly).
// When X aliases Z_block at the same leading dimension, the loader is a no-op
// and we skip the memcpy. Distinct buffers (or X with a different stride) get
// the explicit column copy. If callers ever start passing partially-aliased
// regions (e.g. X == Z_block + offset) this guard is insufficient and the
// caller must materialize a non-aliased temporary first.
static int block_accept_columns_blas3(const double* X, int ldx, int x_cols,
                                      const double* V_locked, int n_locked,
                                      double* V_active, int* m_active,
                                      int m_max, double* Z_block,
                                      int block_capacity, double* coeff,
                                      double* tmp, int n, int max_accept,
                                      int* ortho_passes,
                                      bool reorthogonalize_active = true) {
  if (max_accept < 0) {
    max_accept = 0;
  }
  int cols = x_cols;
  if (cols > max_accept) {
    cols = max_accept;
  }
  if (cols > block_capacity) {
    cols = block_capacity;
  }
  if (cols > m_max - *m_active) {
    cols = m_max - *m_active;
  }
  if (cols <= 0) {
    return 0;
  }
  const bool x_aliases_z = (X == Z_block) && (ldx == n);
  if (!x_aliases_z) {
    for (int col = 0; col < cols; ++col) {
      std::memcpy(Z_block + static_cast<int64_t>(col) * n,
                  X + static_cast<int64_t>(col) * ldx,
                  sizeof(double) * static_cast<size_t>(n));
    }
  }

  const double* active_basis = reorthogonalize_active ? V_active : nullptr;
  const int active_cols = reorthogonalize_active ? *m_active : 0;
  const int reorth_passes_done = block_reorthogonalise_against(
    V_locked, n_locked, active_basis, active_cols,
    Z_block, n, cols, coeff, 2
  );
  if (ortho_passes != nullptr) {
    *ortho_passes += reorth_passes_done;
  }

  if (cols > 0) {
    const char trans_T = 'T';
    const char trans_N = 'N';
    const char right = 'R';
    const char uplo = 'U';
    const char diag = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    int info = 0;
    F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                    &one, Z_block, &n, Z_block, &n,
                    &zero, coeff, &cols FCONE FCONE);
    symmetrize_packed_square(coeff, cols);
    F77_CALL(dpotrf)(&uplo, &cols, coeff, &cols, &info FCONE);
    bool chol_ok = (info == 0);
    for (int col = 0; chol_ok && col < cols; ++col) {
      if (coeff[col + static_cast<int64_t>(col) * cols] <= 100.0 * DBL_EPSILON) {
        chol_ok = false;
      }
    }
    if (chol_ok) {
      F77_CALL(dtrsm)(&right, &uplo, &trans_N, &diag, &n, &cols, &one,
                      coeff, &cols, Z_block, &n FCONE FCONE FCONE FCONE);
      if (n < 64) {
        F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                        &one, Z_block, &n, Z_block, &n,
                        &zero, coeff, &cols FCONE FCONE);
        symmetrize_packed_square(coeff, cols);
        F77_CALL(dpotrf)(&uplo, &cols, coeff, &cols, &info FCONE);
        chol_ok = (info == 0);
        for (int col = 0; chol_ok && col < cols; ++col) {
          if (coeff[col + static_cast<int64_t>(col) * cols] <= 100.0 * DBL_EPSILON) {
            chol_ok = false;
          }
        }
      }
      if (chol_ok) {
        if (n < 64) {
          F77_CALL(dtrsm)(&right, &uplo, &trans_N, &diag, &n, &cols, &one,
                          coeff, &cols, Z_block, &n FCONE FCONE FCONE FCONE);
        }
        for (int col = 0; col < cols && *m_active < m_max; ++col) {
          std::memcpy(V_active + static_cast<int64_t>(*m_active) * n,
                      Z_block + static_cast<int64_t>(col) * n,
                      sizeof(double) * static_cast<size_t>(n));
          ++(*m_active);
        }
        return cols;
      }
    }
  }

  int accepted = 0;
  const int batch_start = *m_active;
  for (int col = 0; col < cols && *m_active < m_max; ++col) {
    double* z_col = Z_block + static_cast<int64_t>(col) * n;
    if (accepted > 0) {
      trl_orthogonalise(nullptr, 0,
                        V_active + static_cast<int64_t>(batch_start) * n,
                        accepted, z_col, tmp, n);
    }
    const double nz = trl_norm2(z_col, n);
    if (nz <= 100.0 * DBL_EPSILON) {
      continue;
    }
    const double inv_nz = 1.0 / nz;
    double* dst = V_active + static_cast<int64_t>(*m_active) * n;
    for (int row = 0; row < n; ++row) {
      dst[row] = z_col[row] * inv_nz;
    }
    ++(*m_active);
    ++accepted;
  }
  return accepted;
}

static int apply_active_block(void* impl, EigencoreApplyFn apply,
                              int n, int first_col, int cols,
                              double* V_active, double* AV_active,
                              EigencoreWorkspace* workspace,
                              int* matvecs_out) {
  if (cols <= 0) {
    return 0;
  }
  const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, cols,
                       V_active + static_cast<int64_t>(first_col) * n, n,
                       1.0, 0.0,
                       AV_active + static_cast<int64_t>(first_col) * n, n,
                       workspace);
  if (rc == 0 && matvecs_out != nullptr) {
    ++(*matvecs_out);
  }
  return rc;
}

static void subtract_projected_range(const double* V_active,
                                     const double* T_proj,
                                     int ldt,
                                     int n,
                                     int range_start,
                                     int range_cols,
                                     int current_start,
                                     int current_cols,
                                     double* W,
                                     double* coeff) {
  if (range_cols <= 0 || current_cols <= 0) {
    return;
  }
  if (range_cols <= 4 && current_cols <= 4) {
    for (int col = 0; col < current_cols; ++col) {
      double* w_col = W + static_cast<int64_t>(col) * n;
      for (int basis = 0; basis < range_cols; ++basis) {
        const double coeff_value =
          T_proj[(range_start + basis) +
                 static_cast<int64_t>(current_start + col) * ldt];
        const double* v_col =
          V_active + static_cast<int64_t>(range_start + basis) * n;
        for (int row = 0; row < n; ++row) {
          w_col[row] -= v_col[row] * coeff_value;
        }
      }
    }
    return;
  }
  for (int col = 0; col < current_cols; ++col) {
    for (int row = 0; row < range_cols; ++row) {
      coeff[row + static_cast<int64_t>(col) * range_cols] =
        T_proj[(range_start + row) +
               static_cast<int64_t>(current_start + col) * ldt];
    }
  }
  const char trans_N = 'N';
  const double one = 1.0;
  const double minus_one = -1.0;
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &current_cols, &range_cols,
                  &minus_one,
                  V_active + static_cast<int64_t>(range_start) * n, &n,
                  coeff, &range_cols,
                  &one, W, &n FCONE FCONE);
}

static void form_structured_projected_block_residual(const double* V_active,
                                                     const double* AV_active,
                                                     const double* T_proj,
                                                     int ldt,
                                                     int n,
                                                     int current_start,
                                                     int current_cols,
                                                     int previous_start,
                                                     int previous_cols,
                                                     double* W,
                                                     double* coeff) {
  if (current_cols <= 0) {
    return;
  }
  for (int col = 0; col < current_cols; ++col) {
    std::memcpy(W + static_cast<int64_t>(col) * n,
                AV_active + static_cast<int64_t>(current_start + col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }
  subtract_projected_range(V_active, T_proj, ldt, n,
                           previous_start, previous_cols,
                           current_start, current_cols, W, coeff);
  subtract_projected_range(V_active, T_proj, ldt, n,
                           current_start, current_cols,
                           current_start, current_cols, W, coeff);
}

static void projection_update_self_block(double* T_proj, int ldt,
                                         const double* V_active,
                                         const double* AV_active,
                                         int n, int start, int cols,
                                         double* scratch) {
  if (cols <= 0) {
    return;
  }
  if (cols <= 4) {
    for (int col = 0; col < cols; ++col) {
      const double* av_col =
        AV_active + static_cast<int64_t>(start + col) * n;
      for (int row = 0; row < cols; ++row) {
        const double* v_col =
          V_active + static_cast<int64_t>(start + row) * n;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
          sum += v_col[i] * av_col[i];
        }
        scratch[row + static_cast<int64_t>(col) * cols] = sum;
      }
    }
    symmetrize_packed_square(scratch, cols);
    for (int col = 0; col < cols; ++col) {
      for (int row = 0; row < cols; ++row) {
        T_proj[(start + row) + static_cast<int64_t>(start + col) * ldt] =
          scratch[row + static_cast<int64_t>(col) * cols];
      }
    }
    return;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                  &one,
                  V_active + static_cast<int64_t>(start) * n, &n,
                  AV_active + static_cast<int64_t>(start) * n, &n,
                  &zero, scratch, &cols FCONE FCONE);
  symmetrize_packed_square(scratch, cols);
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < cols; ++row) {
      T_proj[(start + row) + static_cast<int64_t>(start + col) * ldt] =
        scratch[row + static_cast<int64_t>(col) * cols];
    }
  }
}

static void projection_update_appended_block(double* T_proj, int ldt,
                                             const double* V_active,
                                             const double* AV_active,
                                             int n,
                                             int old_cols,
                                             int new_cols,
                                             double* scratch) {
  if (new_cols <= 0) {
    return;
  }
  if (old_cols <= 0) {
    projection_update_self_block(T_proj, ldt, V_active, AV_active,
                                 n, 0, new_cols, scratch);
    return;
  }

  const int new_start = old_cols;
  if (old_cols <= 4 && new_cols <= 4) {
    for (int col = 0; col < new_cols; ++col) {
      const int abs_col = new_start + col;
      const double* av_col = AV_active + static_cast<int64_t>(abs_col) * n;
      for (int row = 0; row < old_cols; ++row) {
        const double* v_col = V_active + static_cast<int64_t>(row) * n;
        double value = 0.0;
        for (int i = 0; i < n; ++i) {
          value += v_col[i] * av_col[i];
        }
        T_proj[row + static_cast<int64_t>(abs_col) * ldt] = value;
        T_proj[abs_col + static_cast<int64_t>(row) * ldt] = value;
      }
    }
  } else {
    const char trans_T = 'T';
    const char trans_N = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    F77_CALL(dgemm)(&trans_T, &trans_N, &old_cols, &new_cols, &n,
                    &one,
                    V_active, &n,
                    AV_active + static_cast<int64_t>(new_start) * n, &n,
                    &zero, scratch, &old_cols FCONE FCONE);
    for (int col = 0; col < new_cols; ++col) {
      const int abs_col = new_start + col;
      for (int row = 0; row < old_cols; ++row) {
        const double value = scratch[row + static_cast<int64_t>(col) * old_cols];
        T_proj[row + static_cast<int64_t>(abs_col) * ldt] = value;
        T_proj[abs_col + static_cast<int64_t>(row) * ldt] = value;
      }
    }
  }

  projection_update_self_block(T_proj, ldt, V_active, AV_active,
                               n, new_start, new_cols, scratch);
}

static void projection_copy_upper_compact(const double* T_proj, int ldt,
                                          double* compact, int m_active) {
  for (int col = 0; col < m_active; ++col) {
    for (int row = 0; row <= col; ++row) {
      compact[row + static_cast<int64_t>(col) * m_active] =
        T_proj[row + static_cast<int64_t>(col) * ldt];
    }
  }
}

struct ThickRestartBuffers {
  double* V_active;
  double* AV_active;
  double* T_proj;      // m_max x m_max structured projected problem
  // S_eig is a k x k scratch matrix reused across three distinct roles inside
  // a single solve cycle:
  //   role 1: V^T V Gram for the orthogonality probe in final_polish_block_ritz
  //   role 2: Cholesky factor for re-orthonormalization (dpotrf / dtrsm in-place)
  //   role 3: V^T A V projected eigenproblem (dsyev_inplace) for full polish
  // Each role overwrites the previous, in strict sequence within one call.
  // Maintaining three separate buffers would cost an extra 2 * k_max^2 doubles
  // per cycle for negligible runtime savings; the role transitions are flagged
  // with explicit comments where they happen.
  double* S_eig;       // m_max x m_max — see role notes above
  double* S_selected;  // selected Ritz vectors, m_max x selected_capacity
  double* theta;       // m_max
  double* B_v;         // n x selected_capacity
  double* B_av;        // n x selected_capacity
  double* Z_block;     // n x block_size
  double* coeff_block; // m_max x block_size
  double* z;           // n
  double* tmp;         // max(k_target, m_max)
  double* ritz_res;    // m_max
  int*    selected;    // m_max
  int*    is_locked;   // m_max
  double* dsyev_work;
  int     dsyev_lwork;
  int*    dsyevd_iwork;
  int     dsyevd_liwork;
  int     selected_capacity;
};

static void trl_buffers_free(ThickRestartBuffers* b) {
  std::free(b->V_active);
  std::free(b->AV_active);
  std::free(b->T_proj);
  std::free(b->S_eig);
  std::free(b->S_selected);
  std::free(b->theta);
  std::free(b->B_v);
  std::free(b->B_av);
  std::free(b->Z_block);
  std::free(b->coeff_block);
  std::free(b->z);
  std::free(b->tmp);
  std::free(b->ritz_res);
  std::free(b->selected);
  std::free(b->is_locked);
  std::free(b->dsyev_work);
  std::free(b->dsyevd_iwork);
}

static int trl_buffers_alloc(ThickRestartBuffers* b, int n, int k_target,
                             int m_max, int block_cols) {
  std::memset(b, 0, sizeof(*b));
  const size_t nm = static_cast<size_t>(n) * static_cast<size_t>(m_max);
  const size_t mm = static_cast<size_t>(m_max) * static_cast<size_t>(m_max);
  const size_t nb = static_cast<size_t>(n) * static_cast<size_t>(block_cols);
  const size_t mb = static_cast<size_t>(m_max) * static_cast<size_t>(block_cols);
  int selected_capacity = 2 * k_target;
  if (selected_capacity < k_target + 5) selected_capacity = k_target + 5;
  if (selected_capacity < block_cols) selected_capacity = block_cols;
  if (selected_capacity < 1) selected_capacity = 1;
  if (selected_capacity > m_max) selected_capacity = m_max;
  b->selected_capacity = selected_capacity;
  const size_t ms = static_cast<size_t>(m_max) * static_cast<size_t>(selected_capacity);
  const size_t ns = static_cast<size_t>(n) * static_cast<size_t>(selected_capacity);
  b->V_active  = static_cast<double*>(std::malloc(nm * sizeof(double)));
  b->AV_active = static_cast<double*>(std::malloc(nm * sizeof(double)));
  b->T_proj    = static_cast<double*>(std::calloc(mm, sizeof(double)));
  b->S_eig     = static_cast<double*>(std::malloc(mm * sizeof(double)));
  b->S_selected = static_cast<double*>(std::malloc(ms * sizeof(double)));
  b->theta     = static_cast<double*>(std::malloc(static_cast<size_t>(m_max) * sizeof(double)));
  b->B_v       = static_cast<double*>(std::malloc(ns * sizeof(double)));
  b->B_av      = static_cast<double*>(std::malloc(ns * sizeof(double)));
  b->Z_block   = static_cast<double*>(std::malloc(nb * sizeof(double)));
  b->coeff_block = static_cast<double*>(std::malloc(mb * sizeof(double)));
  b->z         = static_cast<double*>(std::malloc(static_cast<size_t>(n) * sizeof(double)));
  const int tmp_len = (k_target > m_max) ? k_target : m_max;
  b->tmp       = static_cast<double*>(std::malloc(static_cast<size_t>(tmp_len > 0 ? tmp_len : 1) * sizeof(double)));
  b->ritz_res  = static_cast<double*>(std::malloc(static_cast<size_t>(m_max) * sizeof(double)));
  b->selected  = static_cast<int*>(std::malloc(static_cast<size_t>(m_max) * sizeof(int)));
  b->is_locked = static_cast<int*>(std::malloc(static_cast<size_t>(m_max) * sizeof(int)));
  b->dsyev_lwork = trl_dsyevd_query(m_max, &b->dsyevd_liwork);
  if (b->dsyev_lwork < 26 * m_max) b->dsyev_lwork = 26 * m_max;
  if (b->dsyevd_liwork < 10 * m_max) b->dsyevd_liwork = 10 * m_max;
  if (b->dsyev_lwork < 1) b->dsyev_lwork = 1;
  if (b->dsyevd_liwork < 1) b->dsyevd_liwork = 1;
  b->dsyev_work = static_cast<double*>(std::malloc(static_cast<size_t>(b->dsyev_lwork) * sizeof(double)));
  b->dsyevd_iwork = static_cast<int*>(std::malloc(static_cast<size_t>(b->dsyevd_liwork) * sizeof(int)));
  if (b->V_active == nullptr || b->AV_active == nullptr ||
      b->T_proj == nullptr || b->S_eig == nullptr ||
      b->S_selected == nullptr ||
      b->theta == nullptr || b->B_v == nullptr || b->B_av == nullptr ||
      b->Z_block == nullptr || b->coeff_block == nullptr ||
      b->z == nullptr || b->tmp == nullptr || b->ritz_res == nullptr ||
      b->selected == nullptr || b->is_locked == nullptr ||
      b->dsyev_work == nullptr || b->dsyevd_iwork == nullptr) {
    trl_buffers_free(b);
    return -1;
  }
  return 0;
}

static int final_polish_block_ritz(void* impl,
                                   EigencoreApplyFn apply,
                                   int n,
                                   int k_target,
                                   int target_kind,
                                   double tol,
                                   double norm_a,
                                   double* V_out,
                                   double* lambda_out,
                                   double* residuals_out,
                                   int* converged_out,
                                   int* n_converged_out,
                                   ThickRestartBuffers* buf,
                                   EigencoreWorkspace* workspace,
                                   int* matvecs_out) {
  if (k_target <= 0) {
    if (n_converged_out != nullptr) {
      *n_converged_out = 0;
    }
    return 0;
  }

  const char trans_T = 'T';
  const char trans_N = 'N';
  const char side_R = 'R';
  const char uplo_U = 'U';
  const char diag_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  // S_eig role 1: V_out^T V_out Gram for orthogonality probe.
  F77_CALL(dgemm)(&trans_T, &trans_N, &k_target, &k_target, &n,
                  &one, V_out, &n, V_out, &n,
                  &zero, buf->S_eig, &k_target FCONE FCONE);
  double max_orthogonality = 0.0;
  for (int col = 0; col < k_target; ++col) {
    for (int row = 0; row < k_target; ++row) {
      const double expected = (row == col) ? 1.0 : 0.0;
      const double loss = fabs(buf->S_eig[row + static_cast<int64_t>(col) * k_target] - expected);
      if (loss > max_orthogonality) {
        max_orthogonality = loss;
      }
    }
  }
  const double orthogonality_tolerance =
    (tol > sqrt(DBL_EPSILON)) ? tol : sqrt(DBL_EPSILON);
  int prepolish_converged = 0;
  if (max_orthogonality <= orthogonality_tolerance) {
    for (int col = 0; col < k_target; ++col) {
      const double* vec = V_out + static_cast<int64_t>(col) * n;
      const double scale =
        standard_eigen_lock_scale(norm_a, lambda_out[col], vec, n);
      converged_out[col] = (residuals_out[col] <= tol * scale) ? 1 : 0;
      if (converged_out[col]) {
        ++prepolish_converged;
      }
    }
    if (n_converged_out != nullptr) {
      *n_converged_out = prepolish_converged;
    }
    if (prepolish_converged == k_target) {
      return 0;
    }
  }

  // S_eig role 2: in-place Cholesky factor of the Gram from role 1, used to
  // re-orthonormalize V_out via right-side dtrsm. Overwrites role-1 contents.
  symmetrize_packed_square(buf->S_eig, k_target);
  F77_CALL(dpotrf)(&uplo_U, &k_target, buf->S_eig, &k_target, &info FCONE);
  if (info != 0) {
    return 0;
  }
  F77_CALL(dtrsm)(&side_R, &uplo_U, &trans_N, &diag_N, &n, &k_target, &one,
                  buf->S_eig, &k_target, V_out, &n FCONE FCONE FCONE FCONE);

  int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, k_target,
                 V_out, n, 1.0, 0.0, buf->B_av, n, workspace);
  if (rc != 0) {
    return rc;
  }
  if (matvecs_out != nullptr) {
    ++(*matvecs_out);
  }

  int n_converged_simple = 0;
  for (int col = 0; col < k_target; ++col) {
    const double* vec = V_out + static_cast<int64_t>(col) * n;
    const double* av = buf->B_av + static_cast<int64_t>(col) * n;
    long double theta_sum = 0.0L;
    for (int row = 0; row < n; ++row) {
      theta_sum += static_cast<long double>(vec[row]) * av[row];
    }
    const double theta = static_cast<double>(theta_sum);
    long double ss = 0.0L;
    for (int row = 0; row < n; ++row) {
      const double diff = av[row] - theta * vec[row];
      ss += static_cast<long double>(diff) * diff;
    }
    const double res = sqrt(static_cast<double>(ss));
    lambda_out[col] = theta;
    residuals_out[col] = res;
    const double scale = standard_eigen_lock_scale(norm_a, theta, vec, n);
    converged_out[col] = (res <= tol * scale) ? 1 : 0;
    if (converged_out[col]) {
      ++n_converged_simple;
    }
  }
  if (n_converged_out != nullptr) {
    *n_converged_out = n_converged_simple;
  }
  if (n_converged_simple == k_target) {
    return 0;
  }

  // S_eig role 3: V_out^T A V_out projected eigenproblem, solved in place by
  // symmetric_eigen_inplace. Overwrites the role-2 Cholesky factor; columns of
  // S_eig now hold the projected-problem eigenvectors used to rotate V_out.
  // Reuse A * V_out from the simple residual check; V_out has not changed.
  F77_CALL(dgemm)(&trans_T, &trans_N, &k_target, &k_target, &n,
                  &one, V_out, &n, buf->B_av, &n,
                  &zero, buf->S_eig, &k_target FCONE FCONE);
  symmetrize_packed_square(buf->S_eig, k_target);
  rc = symmetric_eigen_inplace(buf->S_eig, k_target, buf->theta,
                               buf->dsyev_work, buf->dsyev_lwork,
                               buf->dsyevd_iwork, buf->dsyevd_liwork);
  if (rc != 0) {
    return rc;
  }
  selected_sorted_ritz_indices(buf->theta, k_target, k_target,
                               target_kind, buf->selected);

  for (int col = 0; col < k_target; ++col) {
    const int idx = buf->selected[col];
    for (int row = 0; row < k_target; ++row) {
      buf->S_selected[row + static_cast<int64_t>(col) * k_target] =
        buf->S_eig[row + static_cast<int64_t>(idx) * k_target];
    }
  }

  std::memcpy(buf->AV_active, buf->B_av,
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(k_target));
  if (k_target <= 32) {
    combine_basis_columns_small(V_out, n, k_target,
                                buf->S_selected, k_target,
                                k_target, buf->B_v);
    combine_basis_columns_small(buf->AV_active, n, k_target,
                                buf->S_selected, k_target,
                                k_target, buf->B_av);
  } else {
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &k_target,
                    &one, V_out, &n, buf->S_selected, &k_target,
                    &zero, buf->B_v, &n FCONE FCONE);
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &k_target,
                    &one, buf->AV_active, &n, buf->S_selected, &k_target,
                    &zero, buf->B_av, &n FCONE FCONE);
  }

  int n_converged = 0;
  for (int col = 0; col < k_target; ++col) {
    const int idx = buf->selected[col];
    const double theta = buf->theta[idx];
    long double ss = 0.0L;
    double* residual = buf->B_av + static_cast<int64_t>(col) * n;
    double* vec = buf->B_v + static_cast<int64_t>(col) * n;
    for (int row = 0; row < n; ++row) {
      residual[row] -= theta * vec[row];
      ss += static_cast<long double>(residual[row]) * residual[row];
      V_out[row + static_cast<int64_t>(col) * n] = vec[row];
    }
    const double res = sqrt(static_cast<double>(ss));
    lambda_out[col] = theta;
    residuals_out[col] = res;
    const double scale = standard_eigen_lock_scale(norm_a, theta, vec, n);
    converged_out[col] = (res <= tol * scale) ? 1 : 0;
    if (converged_out[col]) {
      ++n_converged;
    }
  }
  if (n_converged_out != nullptr) {
    *n_converged_out = n_converged;
  }
  return 0;
}

static SEXP trl_pack_result(int n, int k_target, const double* V_locked,
                            const double* lambda, const double* residuals,
                            const int* converged, int n_locked,
                            int iterations, int matvecs, int restarts,
                            int m_active_final) {
  SEXP values_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k_target));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k_target));
  std::memcpy(REAL(values_), lambda, sizeof(double) * static_cast<size_t>(k_target));
  std::memcpy(REAL(vectors_), V_locked,
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(k_target));
  std::memcpy(REAL(residuals_), residuals, sizeof(double) * static_cast<size_t>(k_target));
  for (int i = 0; i < k_target; ++i) {
    LOGICAL(converged_)[i] = converged[i] ? TRUE : FALSE;
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, residuals_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(n_locked));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(restarts));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(m_active_final));
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("residuals"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("n_locked"));
  SET_STRING_ELT(names_, 5, mkChar("iterations"));
  SET_STRING_ELT(names_, 6, mkChar("matvecs"));
  SET_STRING_ELT(names_, 7, mkChar("restarts"));
  SET_STRING_ELT(names_, 8, mkChar("m_active_final"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(6);
  return out_;
}

static int native_block_lanczos_run(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int block_size,
    int target_kind,
    double tol,
    const double* start_block,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* nconv_out,
    int* iterations_out,
    int* matvecs_out,
    int* m_active_final_out
) {
  *nconv_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *m_active_final_out = 0;
  for (int i = 0; i < k_target; ++i) {
    lambda_out[i] = 0.0;
    residuals_out[i] = R_PosInf;
    converged_out[i] = 0;
    std::memset(V_out + static_cast<int64_t>(i) * n, 0,
                sizeof(double) * static_cast<size_t>(n));
  }

  const size_t nm = static_cast<size_t>(n) * static_cast<size_t>(m_max);
  const size_t nb = static_cast<size_t>(n) * static_cast<size_t>(block_size);
  const size_t mm = static_cast<size_t>(m_max) * static_cast<size_t>(m_max);
  double* V = static_cast<double*>(std::calloc(nm, sizeof(double)));
  double* AV = static_cast<double*>(std::calloc(nm, sizeof(double)));
  double* Z = static_cast<double*>(std::calloc(nb, sizeof(double)));
  double* AZ = static_cast<double*>(std::calloc(nb, sizeof(double)));
  double* H = static_cast<double*>(std::calloc(mm, sizeof(double)));
  double* S_selected = static_cast<double*>(std::calloc(mm, sizeof(double)));
  double* theta = static_cast<double*>(std::calloc(static_cast<size_t>(m_max), sizeof(double)));
  double* B_v = static_cast<double*>(std::calloc(static_cast<size_t>(n) * k_target, sizeof(double)));
  double* B_av = static_cast<double*>(std::calloc(static_cast<size_t>(n) * k_target, sizeof(double)));
  double* tmp = static_cast<double*>(std::calloc(static_cast<size_t>(m_max), sizeof(double)));
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>(m_max), sizeof(int)));
  const int dsyev_lwork_query = trl_dsyev_query(m_max);
  int dsyev_lwork = dsyev_lwork_query > 0 ? dsyev_lwork_query : 3 * m_max;
  double* dsyev_work = static_cast<double*>(std::calloc(static_cast<size_t>(dsyev_lwork), sizeof(double)));
  if (V == nullptr || AV == nullptr || Z == nullptr || AZ == nullptr ||
      H == nullptr || S_selected == nullptr || theta == nullptr ||
      B_v == nullptr || B_av == nullptr || tmp == nullptr ||
      selected == nullptr || dsyev_work == nullptr) {
    std::free(V); std::free(AV); std::free(Z); std::free(AZ);
    std::free(H); std::free(S_selected); std::free(theta);
    std::free(B_v); std::free(B_av); std::free(tmp);
    std::free(selected); std::free(dsyev_work);
    return -2;
  }

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int m_active = 0;
  int last_block_start = 0;
  int last_block_cols = 0;
  int source_start = 0;
  int source_cols = block_size;
  for (int col = 0; col < block_size; ++col) {
    std::memcpy(Z + static_cast<int64_t>(col) * n,
                start_block + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  while (m_active < m_max && source_cols > 0) {
    int accepted_start = m_active;
    int accepted = 0;
    for (int col = 0; col < source_cols && m_active < m_max; ++col) {
      double* z_col = Z + static_cast<int64_t>(col) * n;
      trl_orthogonalise(nullptr, 0, V, m_active, z_col, tmp, n);
      const double nz = trl_norm2(z_col, n);
      if (nz <= 100.0 * DBL_EPSILON) {
        continue;
      }
      const double inv_nz = 1.0 / nz;
      for (int row = 0; row < n; ++row) {
        V[static_cast<int64_t>(m_active) * n + row] = z_col[row] * inv_nz;
      }
      ++m_active;
      ++accepted;
    }

    if (accepted == 0) {
      break;
    }

    const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, accepted,
                         V + static_cast<int64_t>(accepted_start) * n, n,
                         1.0, 0.0,
                         AV + static_cast<int64_t>(accepted_start) * n, n,
                         &workspace);
    if (rc != 0) {
      std::free(V); std::free(AV); std::free(Z); std::free(AZ);
      std::free(H); std::free(S_selected); std::free(theta);
      std::free(B_v); std::free(B_av); std::free(tmp);
      std::free(selected); std::free(dsyev_work);
      return rc;
    }
    ++(*matvecs_out);
    ++(*iterations_out);
    last_block_start = accepted_start;
    last_block_cols = accepted;

    if (m_active >= m_max) {
      break;
    }

    source_start = last_block_start;
    source_cols = last_block_cols;
    for (int col = 0; col < source_cols; ++col) {
      std::memcpy(Z + static_cast<int64_t>(col) * n,
                  AV + static_cast<int64_t>(source_start + col) * n,
                  sizeof(double) * static_cast<size_t>(n));
    }
  }

  if (m_active < k_target) {
    std::free(V); std::free(AV); std::free(Z); std::free(AZ);
    std::free(H); std::free(S_selected); std::free(theta);
    std::free(B_v); std::free(B_av); std::free(tmp);
    std::free(selected); std::free(dsyev_work);
    return -4;
  }

  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_T, &trans_N, &m_active, &m_active, &n,
                  &one, V, &n, AV, &n,
                  &zero, H, &m_active FCONE FCONE);
  for (int i = 0; i < m_active; ++i) {
    for (int j = i + 1; j < m_active; ++j) {
      const double avg = 0.5 * (H[i + j * m_active] + H[j + i * m_active]);
      H[i + j * m_active] = avg;
      H[j + i * m_active] = avg;
    }
  }

  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = dsyev_lwork;
  F77_CALL(dsyev)(&jobz, &uplo, &m_active, H, &m_active, theta,
                  dsyev_work, &lwork, &info FCONE FCONE);
  if (info != 0) {
    std::free(V); std::free(AV); std::free(Z); std::free(AZ);
    std::free(H); std::free(S_selected); std::free(theta);
    std::free(B_v); std::free(B_av); std::free(tmp);
    std::free(selected); std::free(dsyev_work);
    return -3;
  }

  selected_ritz_indices(theta, m_active, k_target, target_kind, selected);
  for (int p = 0; p < k_target; ++p) {
    const int idx = selected[p];
    lambda_out[p] = theta[idx];
    for (int row = 0; row < m_active; ++row) {
      S_selected[row + static_cast<int64_t>(p) * m_active] =
        H[row + static_cast<int64_t>(idx) * m_active];
    }
  }

  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &m_active,
                  &one, V, &n, S_selected, &m_active,
                  &zero, B_v, &n FCONE FCONE);
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &m_active,
                  &one, AV, &n, S_selected, &m_active,
                  &zero, B_av, &n FCONE FCONE);
  int nconv = 0;
  for (int col = 0; col < k_target; ++col) {
    long double s = 0.0L;
    for (int row = 0; row < n; ++row) {
      const double diff = B_av[row + static_cast<int64_t>(col) * n] -
                          lambda_out[col] * B_v[row + static_cast<int64_t>(col) * n];
      s += static_cast<long double>(diff) * diff;
    }
    residuals_out[col] = sqrt(static_cast<double>(s));
    const double scale_i = (fabs(lambda_out[col]) > 1.0) ? fabs(lambda_out[col]) : 1.0;
    converged_out[col] = residuals_out[col] <= tol * scale_i ? 1 : 0;
    if (converged_out[col]) ++nconv;
    std::memcpy(V_out + static_cast<int64_t>(col) * n,
                B_v + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  *nconv_out = nconv;
  *m_active_final_out = m_active;
  std::free(V); std::free(AV); std::free(Z); std::free(AZ);
  std::free(H); std::free(S_selected); std::free(theta);
  std::free(B_v); std::free(B_av); std::free(tmp);
  std::free(selected); std::free(dsyev_work);
  return 0;
}

struct BlockLanczosBestSnapshot {
  std::vector<double> V;
  std::vector<double> lambda;
  std::vector<double> residuals;
  std::vector<int> converged;
  std::vector<double> candidate_V;
  std::vector<int> candidate_converged;
  int filled = 0;
  int locked_prefix = 0;
  int nconv = -1;
  double max_backward_error = R_PosInf;

  BlockLanczosBestSnapshot(int n, int k_target) :
      V(static_cast<size_t>(n) * static_cast<size_t>(k_target), 0.0),
      lambda(static_cast<size_t>(k_target), 0.0),
      residuals(static_cast<size_t>(k_target), R_PosInf),
      converged(static_cast<size_t>(k_target), 0),
      candidate_V(static_cast<size_t>(n) * static_cast<size_t>(k_target), 0.0),
      candidate_converged(static_cast<size_t>(k_target), 0) {}
};

static int block_lanczos_expand_basis_to_budget(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int m_max,
    int block_size,
    const double* V_locked,
    int n_locked,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    int* m_active,
    int* previous_block_start,
    int* previous_block_cols,
    int* last_block_start,
    int* last_block_cols,
    int* iterations_out,
    int* matvecs_out,
    int* ortho_passes_out) {
  while (*m_active < m_max && *last_block_cols > 0) {
    auto timer = native_timer_now();
    form_structured_projected_block_residual(
      buf->V_active, buf->AV_active, buf->T_proj, m_max, n,
      *last_block_start, *last_block_cols,
      *previous_block_start, *previous_block_cols,
      buf->Z_block, buf->coeff_block
    );
    stages->recurrence += native_timer_elapsed(timer);

    const int accepted_start = *m_active;
    timer = native_timer_now();
    const int accepted = block_accept_columns_blas3(
      buf->Z_block, n, *last_block_cols, V_locked, n_locked,
      buf->V_active, m_active, m_max, buf->Z_block, block_size,
      buf->coeff_block, buf->tmp, n,
      block_size, ortho_passes_out,
      true
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    if (accepted == 0) {
      // Exact breakdown can occur in an invariant subspace that is not the
      // requested target. Continue deterministically instead of locking it.
      int continuation_accepted = 0;
      const int continuation_start = *m_active;
      for (int attempt = 0;
           continuation_accepted < block_size &&
             attempt < n + block_size &&
             *m_active < m_max;
           ++attempt) {
        std::memset(buf->z, 0, sizeof(double) * static_cast<size_t>(n));
        const int idx_basis = ((*iterations_out + 1) * 17 + attempt * 31) % n;
        buf->z[idx_basis < 0 ? -idx_basis : idx_basis] = 1.0;
        timer = native_timer_now();
        continuation_accepted += block_accept_work_vector(
          V_locked, n_locked, buf->V_active, m_active, m_max,
          buf->z, buf->tmp, n, ortho_passes_out
        );
        stages->reorthogonalization += native_timer_elapsed(timer);
      }
      if (continuation_accepted == 0) {
        break;
      }

      timer = native_timer_now();
      const int rc = apply_active_block(
        impl, apply, n, continuation_start, continuation_accepted,
        buf->V_active, buf->AV_active, workspace, matvecs_out
      );
      stages->apply += native_timer_elapsed(timer);
      if (rc != 0) {
        return rc;
      }

      timer = native_timer_now();
      projection_update_appended_block(buf->T_proj, m_max, buf->V_active,
                                       buf->AV_active, n, continuation_start,
                                       continuation_accepted, buf->S_eig);
      {
        const double elapsed = native_timer_elapsed(timer);
        stages->projected_solve += elapsed;
        stages->projection_update += elapsed;
      }
      ++(*iterations_out);
      *previous_block_start = *last_block_start;
      *previous_block_cols = *last_block_cols;
      *last_block_start = continuation_start;
      *last_block_cols = continuation_accepted;
      continue;
    }

    timer = native_timer_now();
    const int rc = apply_active_block(impl, apply, n, accepted_start, accepted,
                                      buf->V_active, buf->AV_active, workspace,
                                      matvecs_out);
    stages->apply += native_timer_elapsed(timer);
    if (rc != 0) {
      return rc;
    }

    timer = native_timer_now();
    projection_update_appended_block(buf->T_proj, m_max, buf->V_active,
                                     buf->AV_active, n, accepted_start,
                                     accepted, buf->S_eig);
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projection_update += elapsed;
    }
    ++(*iterations_out);
    *previous_block_start = *last_block_start;
    *previous_block_cols = *last_block_cols;
    *last_block_start = accepted_start;
    *last_block_cols = accepted;
  }
  return 0;
}

static void block_lanczos_maybe_capture_best_snapshot(
    int n,
    int k_target,
    int selected_count,
    int n_locked,
    double norm_a,
    double tol,
    const ThickRestartBuffers* buf,
    const double* V_out,
    const double* lambda_out,
    const double* residuals_out,
    BlockLanczosBestSnapshot* best) {
  int candidate_count = 0;
  int candidate_nconv = 0;
  double candidate_max_backward_error = 0.0;
  for (; candidate_count < n_locked && candidate_count < k_target; ++candidate_count) {
    const double* vec = V_out + static_cast<int64_t>(candidate_count) * n;
    std::memcpy(best->candidate_V.data() + static_cast<int64_t>(candidate_count) * n,
                vec,
                sizeof(double) * static_cast<size_t>(n));
    const double scale_i = standard_eigen_lock_scale(
      norm_a, lambda_out[candidate_count], vec, n
    );
    const double backward_error = residuals_out[candidate_count] / scale_i;
    if (backward_error > candidate_max_backward_error) {
      candidate_max_backward_error = backward_error;
    }
    best->candidate_converged[static_cast<size_t>(candidate_count)] =
      (residuals_out[candidate_count] <= tol * scale_i) ? 1 : 0;
    if (best->candidate_converged[static_cast<size_t>(candidate_count)]) {
      ++candidate_nconv;
    }
  }
  for (int p = 0; p < selected_count && candidate_count < k_target; ++p) {
    if (buf->is_locked[p]) {
      continue;
    }
    const int idx = buf->selected[p];
    const double* vec = buf->B_v + static_cast<int64_t>(p) * n;
    if (!vector_is_independent_from_locked(best->candidate_V.data(), candidate_count, vec, n)) {
      continue;
    }
    const double scale_i = standard_eigen_lock_scale(
      norm_a, buf->theta[idx], vec, n
    );
    const double backward_error = buf->ritz_res[p] / scale_i;
    if (backward_error > candidate_max_backward_error) {
      candidate_max_backward_error = backward_error;
    }
    best->candidate_converged[static_cast<size_t>(candidate_count)] =
      (buf->ritz_res[p] <= tol * scale_i) ? 1 : 0;
    if (best->candidate_converged[static_cast<size_t>(candidate_count)]) {
      ++candidate_nconv;
    }
    std::memcpy(best->candidate_V.data() + static_cast<int64_t>(candidate_count) * n,
                vec,
                sizeof(double) * static_cast<size_t>(n));
    ++candidate_count;
  }
  if (candidate_count != k_target ||
      (candidate_nconv < best->nconv ||
       (candidate_nconv == best->nconv &&
        candidate_max_backward_error >= best->max_backward_error))) {
    return;
  }

  int out_col = 0;
  std::memcpy(best->V.data(), best->candidate_V.data(),
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(k_target));
  for (; out_col < n_locked && out_col < k_target; ++out_col) {
    best->lambda[static_cast<size_t>(out_col)] = lambda_out[out_col];
    best->residuals[static_cast<size_t>(out_col)] = residuals_out[out_col];
  }
  for (int p = 0; p < selected_count && out_col < k_target; ++p) {
    if (buf->is_locked[p]) {
      continue;
    }
    const int idx = buf->selected[p];
    const double* vec = buf->B_v + static_cast<int64_t>(p) * n;
    if (!vector_is_independent_from_locked(best->V.data(), out_col, vec, n)) {
      continue;
    }
    best->lambda[static_cast<size_t>(out_col)] = buf->theta[idx];
    best->residuals[static_cast<size_t>(out_col)] = buf->ritz_res[p];
    ++out_col;
  }
  best->locked_prefix = n_locked;
  best->nconv = candidate_nconv;
  best->max_backward_error = candidate_max_backward_error;
  std::memcpy(best->converged.data(), best->candidate_converged.data(),
              sizeof(int) * static_cast<size_t>(k_target));
  best->filled = 1;
}

static int block_lanczos_restart_with_continuation_tail(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int block_size,
    int restart_idx,
    int selected_count,
    int n_locked,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    int* m_active,
    int* previous_block_start,
    int* previous_block_cols,
    int* last_block_start,
    int* last_block_cols,
    double* V_out,
    int* matvecs_out,
    int* restarts_out,
    int* ortho_passes_out) {
  auto timer = native_timer_now();
  const int remaining = k_target - n_locked;
  int keep_room = m_max - block_size;
  if (keep_room < 0) {
    keep_room = 0;
  }
  int pad = block_size > 4 ? block_size : 4;
  if (pad > k_target) {
    pad = k_target;
  }
  int k_keep = remaining + pad;
  if (k_keep < remaining) {
    k_keep = remaining;
  }
  if (k_keep > keep_room) {
    k_keep = keep_room;
  }
  int unlocked_count = 0;
  for (int p = 0; p < selected_count; ++p) {
    if (!buf->is_locked[p]) {
      ++unlocked_count;
    }
  }
  if (k_keep > unlocked_count) {
    k_keep = unlocked_count;
  }

  std::memset(buf->T_proj, 0,
              sizeof(double) * static_cast<size_t>(m_max) *
                static_cast<size_t>(m_max));
  *m_active = 0;
  int n_picked = 0;
  for (int p = 0; p < selected_count && n_picked < k_keep && *m_active < m_max; ++p) {
    if (buf->is_locked[p]) {
      continue;
    }
    std::memcpy(buf->z,
                buf->B_v + static_cast<int64_t>(p) * n,
                sizeof(double) * static_cast<size_t>(n));
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    const int accepted = block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
    n_picked += accepted;
  }

  const int tail_start = *m_active;
  int tail_accepted = 0;
  for (int p = 0; p < selected_count && tail_accepted < block_size && *m_active < m_max; ++p) {
    if (buf->is_locked[p] || buf->ritz_res[p] <= 100.0 * DBL_EPSILON) {
      continue;
    }
    const int idx = buf->selected[p];
    for (int row = 0; row < n; ++row) {
      buf->z[row] = buf->B_av[static_cast<int64_t>(p) * n + row] -
        buf->theta[idx] * buf->B_v[static_cast<int64_t>(p) * n + row];
    }
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    tail_accepted += block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
  }
  for (int attempt = 0; tail_accepted == 0 && attempt < n + block_size && *m_active < m_max; ++attempt) {
    std::memset(buf->z, 0, sizeof(double) * static_cast<size_t>(n));
    const int idx_basis = ((restart_idx + 1) * 17 + attempt * 31) % n;
    buf->z[idx_basis < 0 ? -idx_basis : idx_basis] = 1.0;
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    tail_accepted += block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
  }
  stages->restart += native_timer_elapsed(timer);
  if (tail_accepted == 0) {
    return 1;
  }

  timer = native_timer_now();
  int rc = apply_active_block(impl, apply, n, 0, *m_active,
                              buf->V_active, buf->AV_active, workspace,
                              matvecs_out);
  stages->apply += native_timer_elapsed(timer);
  if (rc != 0) {
    return rc;
  }

  timer = native_timer_now();
  projection_update_self_block(buf->T_proj, m_max, buf->V_active,
                               buf->AV_active, n, 0, *m_active,
                               buf->S_eig);
  {
    const double elapsed = native_timer_elapsed(timer);
    stages->projected_solve += elapsed;
    stages->projection_update += elapsed;
  }
  *previous_block_start = 0;
  *previous_block_cols = tail_start;
  *last_block_start = tail_start;
  *last_block_cols = tail_accepted;
  *restarts_out = restart_idx + 1;
  return 0;
}

static int block_lanczos_finalize_return(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int target_kind,
    double tol,
    double norm_a,
    int selected_count_final,
    bool have_last_rr,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    const BlockLanczosBestSnapshot& best,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* n_locked,
    int* matvecs_out) {
  int n_returned = *n_locked;
  if (best.filled) {
    std::memcpy(V_out, best.V.data(),
                sizeof(double) * static_cast<size_t>(n) *
                  static_cast<size_t>(k_target));
    std::memcpy(lambda_out, best.lambda.data(),
                sizeof(double) * static_cast<size_t>(k_target));
    std::memcpy(residuals_out, best.residuals.data(),
                sizeof(double) * static_cast<size_t>(k_target));
    std::memcpy(converged_out, best.converged.data(),
                sizeof(int) * static_cast<size_t>(k_target));
    *n_locked = best.locked_prefix;
    n_returned = k_target;
  } else if (have_last_rr && n_returned < k_target) {
    for (int p = 0; p < selected_count_final && n_returned < k_target; ++p) {
      if (buf->is_locked[p]) {
        continue;
      }
      const int idx = buf->selected[p];
      const double* vec = buf->B_v + static_cast<int64_t>(p) * n;
      if (!vector_is_independent_from_locked(V_out, n_returned, vec, n)) {
        continue;
      }
      std::memcpy(V_out + static_cast<int64_t>(n_returned) * n,
                  vec,
                  sizeof(double) * static_cast<size_t>(n));
      lambda_out[n_returned] = buf->theta[idx];
      residuals_out[n_returned] = buf->ritz_res[p];
      converged_out[n_returned] = 0;
      ++n_returned;
    }
  }
  if (n_returned != k_target) {
    return 0;
  }

  const int locked_prefix = *n_locked;
  const int full_best_snapshot = best.filled && best.nconv >= k_target;
  const int polish_offset = (!full_best_snapshot &&
                             locked_prefix > 0 && locked_prefix < k_target) ?
    locked_prefix : 0;
  const int polish_count = k_target - polish_offset;
  int polished_converged = 0;
  auto timer = native_timer_now();
  int polish_status = 0;
  if (polish_count > 0) {
    polish_status = final_polish_block_ritz(
      impl, apply, n, polish_count, target_kind, tol, norm_a,
      V_out + static_cast<int64_t>(polish_offset) * n,
      lambda_out + polish_offset,
      residuals_out + polish_offset,
      converged_out + polish_offset,
      &polished_converged, buf, workspace, matvecs_out
    );
  } else {
    polished_converged = 0;
  }
  {
    const double elapsed = native_timer_elapsed(timer);
    stages->ritz_residual += elapsed;
    stages->ritz_final_polish += elapsed;
  }
  if (polish_status != 0) {
    return polish_status;
  }
  *n_locked = polish_offset + polished_converged;
  if (best.filled && polish_offset > 0 && *n_locked >= k_target) {
    timer = native_timer_now();
    polish_status = final_polish_block_ritz(
      impl, apply, n, k_target, target_kind, tol, norm_a,
      V_out, lambda_out, residuals_out, converged_out,
      &polished_converged, buf, workspace, matvecs_out
    );
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_final_polish += elapsed;
    }
    if (polish_status != 0) {
      return polish_status;
    }
    *n_locked = polished_converged;
  }
  return 0;
}

static int native_block_thick_restart_lanczos_run(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int block_size,
    int target_kind,
    double tol,
    int max_restarts,
    double norm_a,
    int apply_ritz_vectors,
    const double* start_block,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* n_locked_out,
    int* iterations_out,
    int* matvecs_out,
    int* restarts_out,
    int* m_active_final_out,
    int* locking_events_out,
    int* ortho_passes_out,
    int64_t* operator_allocations_out,
    int64_t* operator_bytes_allocated_out,
    NativeBlockStageSeconds* stage_out,
    NativeBlockRestartHistory* history
) {
  *n_locked_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *restarts_out = 0;
  *m_active_final_out = 0;
  *locking_events_out = 0;
  *ortho_passes_out = 0;
  *operator_allocations_out = 0;
  *operator_bytes_allocated_out = 0;
  if (stage_out != nullptr) {
    *stage_out = NativeBlockStageSeconds();
  }
  if (history != nullptr) {
    history->length = 0;
  }
  NativeBlockStageSeconds stage_local;
  NativeBlockStageSeconds* stages = (stage_out != nullptr) ? stage_out : &stage_local;
  for (int i = 0; i < k_target; ++i) {
    lambda_out[i] = 0.0;
    residuals_out[i] = R_PosInf;
    converged_out[i] = 0;
    std::memset(V_out + static_cast<int64_t>(i) * n, 0,
                sizeof(double) * static_cast<size_t>(n));
  }

  ThickRestartBuffers buf;
  if (trl_buffers_alloc(&buf, n, k_target, m_max, block_size) != 0) {
    return -2;
  }
  BlockLanczosBestSnapshot best(n, k_target);

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int m_active = 0;
  int n_locked = 0;
  int previous_block_start = 0;
  int previous_block_cols = 0;
  int last_block_start = 0;
  auto timer = native_timer_now();
  int last_block_cols = block_accept_columns_blas3(
    start_block, n, block_size, V_out, n_locked,
    buf.V_active, &m_active, m_max, buf.Z_block, block_size,
    buf.coeff_block, buf.tmp, n,
    block_size, ortho_passes_out
  );
  stages->reorthogonalization += native_timer_elapsed(timer);
  if (last_block_cols == 0) {
    for (int attempt = 0; attempt < block_size && m_active < m_max; ++attempt) {
      std::memset(buf.z, 0, sizeof(double) * static_cast<size_t>(n));
      buf.z[attempt % n] = 1.0;
      timer = native_timer_now();
      last_block_cols += block_accept_work_vector(
        V_out, n_locked, buf.V_active, &m_active, m_max,
        buf.z, buf.tmp, n, ortho_passes_out
      );
      stages->reorthogonalization += native_timer_elapsed(timer);
    }
  }
  timer = native_timer_now();
  int rc = apply_active_block(impl, apply, n, 0, last_block_cols,
                              buf.V_active, buf.AV_active, &workspace,
                              matvecs_out);
  stages->apply += native_timer_elapsed(timer);
  if (rc != 0) {
    trl_buffers_free(&buf);
    return rc;
  }
  timer = native_timer_now();
  projection_update_self_block(buf.T_proj, m_max, buf.V_active, buf.AV_active,
                               n, 0, last_block_cols, buf.coeff_block);
  {
    const double elapsed = native_timer_elapsed(timer);
    stages->projected_solve += elapsed;
    stages->projection_update += elapsed;
  }

  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int selected_count_final = 0;
  int restart_idx = 0;
  bool have_last_rr = false;

  for (; restart_idx <= max_restarts; ++restart_idx) {
    rc = block_lanczos_expand_basis_to_budget(
      impl, apply, n, m_max, block_size, V_out, n_locked, &buf, &workspace,
      stages, &m_active, &previous_block_start, &previous_block_cols,
      &last_block_start, &last_block_cols, iterations_out, matvecs_out,
      ortho_passes_out
    );
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }

    if (m_active < 1) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    const int remaining_before_lock = k_target - n_locked;
    int pad = block_size > 4 ? block_size : 4;
    if (pad > k_target) {
      pad = k_target;
    }
    int selected_count = remaining_before_lock + pad;
    if (selected_count < remaining_before_lock) {
      selected_count = remaining_before_lock;
    }
    if (selected_count > m_active) {
      selected_count = m_active;
    }
    if (selected_count > buf.selected_capacity) {
      selected_count = buf.selected_capacity;
    }
    if (selected_count < 1) {
      selected_count = 1;
    }
    selected_count_final = selected_count;

    timer = native_timer_now();
    projection_copy_upper_compact(buf.T_proj, m_max, buf.S_eig, m_active);
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projection_copy += elapsed;
    }
    timer = native_timer_now();
    rc = symmetric_eigen_inplace(buf.S_eig, m_active, buf.theta,
                                 buf.dsyev_work, buf.dsyev_lwork,
                                 buf.dsyevd_iwork, buf.dsyevd_liwork);
    if (rc == 0) {
      selected_sorted_ritz_indices(buf.theta, m_active, selected_count, target_kind, buf.selected);
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projected_eigensolve += elapsed;
    }
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }
    have_last_rr = true;

    timer = native_timer_now();
    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      std::memcpy(buf.S_selected + static_cast<int64_t>(p) * m_active,
                  buf.S_eig + static_cast<int64_t>(idx) * m_active,
                  sizeof(double) * static_cast<size_t>(m_active));
      buf.ritz_res[p] = R_PosInf;
      buf.is_locked[p] = 0;
    }
    stages->selected_vector_copy += native_timer_elapsed(timer);

    timer = native_timer_now();
    if (selected_count <= 32) {
      combine_basis_columns_small(buf.V_active, n, m_active,
                                  buf.S_selected, m_active,
                                  selected_count, buf.B_v);
    } else {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &selected_count, &m_active,
                      &one, buf.V_active, &n, buf.S_selected, &m_active,
                      &zero, buf.B_v, &n FCONE FCONE);
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_vector_form += elapsed;
    }

    timer = native_timer_now();
    rc = 0;
    if (apply_ritz_vectors) {
      rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, selected_count,
                 buf.B_v, n, 1.0, 0.0, buf.B_av, n, &workspace);
      if (rc == 0 && matvecs_out != nullptr) {
        ++(*matvecs_out);
      }
    } else {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &selected_count, &m_active,
                      &one, buf.AV_active, &n, buf.S_selected, &m_active,
                      &zero, buf.B_av, &n FCONE FCONE);
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_operator_apply += elapsed;
    }
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }

    timer = native_timer_now();
    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      long double s = 0.0L;
      const double* av = buf.B_av + static_cast<int64_t>(p) * n;
      const double* vec = buf.B_v + static_cast<int64_t>(p) * n;
      for (int row = 0; row < n; ++row) {
        const double diff = av[row] - buf.theta[idx] * vec[row];
        s += static_cast<long double>(diff) * diff;
      }
      buf.ritz_res[p] = sqrt(static_cast<double>(s));
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_norm += elapsed;
    }

    const int wanted = k_target - n_locked;
    int history_slot = -1;
    if (history != nullptr && history->length < history->capacity) {
      history_slot = history->length;
      ++history->length;
      const int wanted_selected = wanted < selected_count ? wanted : selected_count;
      int nconv_wanted = 0;
      double max_residual = 0.0;
      double max_backward_error = 0.0;
      for (int p = 0; p < wanted_selected; ++p) {
        const int idx = buf.selected[p];
        const double scale_i = standard_eigen_lock_scale(
          norm_a, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
        );
        const double backward_error = buf.ritz_res[p] / scale_i;
        if (buf.ritz_res[p] > max_residual) {
          max_residual = buf.ritz_res[p];
        }
        if (backward_error > max_backward_error) {
          max_backward_error = backward_error;
        }
        if (buf.ritz_res[p] <= tol * scale_i) {
          ++nconv_wanted;
        }
      }
      history->restart[history_slot] = restart_idx;
      history->m_active[history_slot] = m_active;
      history->selected_count[history_slot] = selected_count;
      history->locked_before[history_slot] = n_locked;
      history->locked_after[history_slot] = n_locked;
      history->nconv_wanted[history_slot] = nconv_wanted;
      history->max_residual[history_slot] = max_residual;
      history->max_backward_error[history_slot] = max_backward_error;
    }

    timer = native_timer_now();
    int lock_now = 0;
    for (int p = 0; p < wanted && p < selected_count; ++p) {
      const int idx = buf.selected[p];
      const double scale_i = standard_eigen_lock_scale(
        norm_a, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
      );
      if (buf.ritz_res[p] <= tol * scale_i) {
        if (!vector_is_independent_from_locked(
              V_out, n_locked, buf.B_v + static_cast<int64_t>(p) * n, n)) {
          buf.is_locked[p] = 1;
          continue;
        }
        std::memcpy(V_out + static_cast<int64_t>(n_locked) * n,
                    buf.B_v + static_cast<int64_t>(p) * n,
                    sizeof(double) * static_cast<size_t>(n));
        lambda_out[n_locked] = buf.theta[idx];
        residuals_out[n_locked] = buf.ritz_res[p];
        converged_out[n_locked] = 1;
        buf.is_locked[p] = 1;
        ++n_locked;
        ++lock_now;
      } else {
        break;
      }
    }
    if (lock_now > 0) {
      ++(*locking_events_out);
    }
    if (history_slot >= 0) {
      history->locked_after[history_slot] = n_locked;
    }
    stages->locking += native_timer_elapsed(timer);

    block_lanczos_maybe_capture_best_snapshot(
      n, k_target, selected_count, n_locked, norm_a, tol, &buf,
      V_out, lambda_out, residuals_out, &best
    );

    if (n_locked >= k_target || restart_idx == max_restarts) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    rc = block_lanczos_restart_with_continuation_tail(
      impl, apply, n, k_target, m_max, block_size, restart_idx, selected_count,
      n_locked, &buf, &workspace, stages, &m_active, &previous_block_start,
      &previous_block_cols, &last_block_start, &last_block_cols, V_out,
      matvecs_out, restarts_out, ortho_passes_out
    );
    if (rc == 1) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }
  }

  rc = block_lanczos_finalize_return(
    impl, apply, n, k_target, target_kind, tol, norm_a, selected_count_final,
    have_last_rr, &buf, &workspace, stages, best, V_out, lambda_out,
    residuals_out, converged_out, &n_locked, matvecs_out
  );
  if (rc != 0) {
    trl_buffers_free(&buf);
    return rc;
  }
  *n_locked_out = n_locked;
  *m_active_final_out = m_active;
  *operator_allocations_out = workspace.allocation_count;
  *operator_bytes_allocated_out = workspace.bytes_allocated;
  trl_buffers_free(&buf);
  return 0;
}

static SEXP block_lanczos_pack_result(int n, int k_target, const double* V,
                                      const double* lambda,
                                      const double* residuals,
                                      const int* converged, int nconv,
                                      int iterations, int matvecs,
                                      int m_active_final) {
  return trl_pack_result(n, k_target, V, lambda, residuals, converged, nconv,
                         iterations, matvecs, 0, m_active_final);
}

static SEXP block_thick_lanczos_pack_result(int n, int k_target, const double* V,
                                            const double* lambda,
                                            const double* residuals,
                                            const int* converged, int n_locked,
                                            int iterations, int matvecs,
                                            int restarts, int m_active_final,
                                            int locking_events, int ortho_passes,
                                            int block_size,
                                            int64_t operator_allocations,
                                            int64_t operator_bytes_allocated,
                                            const NativeBlockStageSeconds* stage_seconds,
                                            const NativeBlockRestartHistory* history) {
  SEXP values_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k_target));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k_target));
  std::memcpy(REAL(values_), lambda, sizeof(double) * static_cast<size_t>(k_target));
  std::memcpy(REAL(vectors_), V,
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(k_target));
  std::memcpy(REAL(residuals_), residuals, sizeof(double) * static_cast<size_t>(k_target));
  for (int i = 0; i < k_target; ++i) {
    LOGICAL(converged_)[i] = converged[i] ? TRUE : FALSE;
  }

  SEXP stage_ = PROTECT(allocVector(REALSXP, 15));
  REAL(stage_)[0] = stage_seconds != nullptr ? stage_seconds->apply : 0.0;
  REAL(stage_)[1] = stage_seconds != nullptr ? stage_seconds->recurrence : 0.0;
  REAL(stage_)[2] = stage_seconds != nullptr ? stage_seconds->reorthogonalization : 0.0;
  REAL(stage_)[3] = stage_seconds != nullptr ? stage_seconds->projected_solve : 0.0;
  REAL(stage_)[4] = stage_seconds != nullptr ? stage_seconds->projection_update : 0.0;
  REAL(stage_)[5] = stage_seconds != nullptr ? stage_seconds->projection_copy : 0.0;
  REAL(stage_)[6] = stage_seconds != nullptr ? stage_seconds->projected_eigensolve : 0.0;
  REAL(stage_)[7] = stage_seconds != nullptr ? stage_seconds->selected_vector_copy : 0.0;
  REAL(stage_)[8] = stage_seconds != nullptr ? stage_seconds->ritz_residual : 0.0;
  REAL(stage_)[9] = stage_seconds != nullptr ? stage_seconds->ritz_vector_form : 0.0;
  REAL(stage_)[10] = stage_seconds != nullptr ? stage_seconds->ritz_operator_apply : 0.0;
  REAL(stage_)[11] = stage_seconds != nullptr ? stage_seconds->ritz_norm : 0.0;
  REAL(stage_)[12] = stage_seconds != nullptr ? stage_seconds->ritz_final_polish : 0.0;
  REAL(stage_)[13] = stage_seconds != nullptr ? stage_seconds->locking : 0.0;
  REAL(stage_)[14] = stage_seconds != nullptr ? stage_seconds->restart : 0.0;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 15));
  SET_STRING_ELT(stage_names_, 0, mkChar("apply"));
  SET_STRING_ELT(stage_names_, 1, mkChar("recurrence"));
  SET_STRING_ELT(stage_names_, 2, mkChar("reorthogonalization"));
  SET_STRING_ELT(stage_names_, 3, mkChar("projected_solve"));
  SET_STRING_ELT(stage_names_, 4, mkChar("projection_update"));
  SET_STRING_ELT(stage_names_, 5, mkChar("projection_copy"));
  SET_STRING_ELT(stage_names_, 6, mkChar("projected_eigensolve"));
  SET_STRING_ELT(stage_names_, 7, mkChar("selected_vector_copy"));
  SET_STRING_ELT(stage_names_, 8, mkChar("ritz_residual"));
  SET_STRING_ELT(stage_names_, 9, mkChar("ritz_vector_form"));
  SET_STRING_ELT(stage_names_, 10, mkChar("ritz_operator_apply"));
  SET_STRING_ELT(stage_names_, 11, mkChar("ritz_norm"));
  SET_STRING_ELT(stage_names_, 12, mkChar("ritz_final_polish"));
  SET_STRING_ELT(stage_names_, 13, mkChar("locking"));
  SET_STRING_ELT(stage_names_, 14, mkChar("restart"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  const int history_length =
    (history != nullptr && history->length > 0) ? history->length : 0;
  SEXP history_ = PROTECT(allocVector(VECSXP, 8));
  SEXP history_restart_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_m_active_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_selected_count_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_locked_before_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_locked_after_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_nconv_wanted_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, history_length));
  SEXP history_max_backward_error_ = PROTECT(allocVector(REALSXP, history_length));
  if (history_length > 0) {
    std::memcpy(INTEGER(history_restart_), history->restart,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_m_active_), history->m_active,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_selected_count_), history->selected_count,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_locked_before_), history->locked_before,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_locked_after_), history->locked_after,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_nconv_wanted_), history->nconv_wanted,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(REAL(history_max_residual_), history->max_residual,
                sizeof(double) * static_cast<size_t>(history_length));
    std::memcpy(REAL(history_max_backward_error_), history->max_backward_error,
                sizeof(double) * static_cast<size_t>(history_length));
  }
  SET_VECTOR_ELT(history_, 0, history_restart_);
  SET_VECTOR_ELT(history_, 1, history_m_active_);
  SET_VECTOR_ELT(history_, 2, history_selected_count_);
  SET_VECTOR_ELT(history_, 3, history_locked_before_);
  SET_VECTOR_ELT(history_, 4, history_locked_after_);
  SET_VECTOR_ELT(history_, 5, history_nconv_wanted_);
  SET_VECTOR_ELT(history_, 6, history_max_residual_);
  SET_VECTOR_ELT(history_, 7, history_max_backward_error_);
  SEXP history_names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(history_names_, 0, mkChar("restart"));
  SET_STRING_ELT(history_names_, 1, mkChar("m_active"));
  SET_STRING_ELT(history_names_, 2, mkChar("selected_count"));
  SET_STRING_ELT(history_names_, 3, mkChar("locked_before"));
  SET_STRING_ELT(history_names_, 4, mkChar("locked_after"));
  SET_STRING_ELT(history_names_, 5, mkChar("nconv_wanted"));
  SET_STRING_ELT(history_names_, 6, mkChar("max_residual"));
  SET_STRING_ELT(history_names_, 7, mkChar("max_backward_error"));
  setAttrib(history_, R_NamesSymbol, history_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, residuals_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(n_locked));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(restarts));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(m_active_final));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(locking_events));
  SET_VECTOR_ELT(out_, 10, ScalarInteger(ortho_passes));
  SET_VECTOR_ELT(out_, 11, ScalarInteger(block_size));
  SET_VECTOR_ELT(out_, 12, ScalarReal(static_cast<double>(operator_allocations)));
  SET_VECTOR_ELT(out_, 13, ScalarReal(static_cast<double>(operator_bytes_allocated)));
  SET_VECTOR_ELT(out_, 14, stage_);
  SET_VECTOR_ELT(out_, 15, history_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("residuals"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("n_locked"));
  SET_STRING_ELT(names_, 5, mkChar("iterations"));
  SET_STRING_ELT(names_, 6, mkChar("matvecs"));
  SET_STRING_ELT(names_, 7, mkChar("restarts"));
  SET_STRING_ELT(names_, 8, mkChar("m_active_final"));
  SET_STRING_ELT(names_, 9, mkChar("locking_events"));
  SET_STRING_ELT(names_, 10, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 11, mkChar("block"));
  SET_STRING_ELT(names_, 12, mkChar("operator_allocations"));
  SET_STRING_ELT(names_, 13, mkChar("operator_bytes_allocated"));
  SET_STRING_ELT(names_, 14, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("restart_history"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(18);
  return out_;
}

extern "C" SEXP eigencore_block_lanczos_dense(SEXP A_, SEXP k_,
                                              SEXP m_max_,
                                              SEXP block_size_,
                                              SEXP target_kind_,
                                              SEXP tol_,
                                              SEXP start_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable block Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k) error("m_max must be >= k");
  if (m_max > n) error("m_max must be <= nrow(A)");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int nconv = 0, iterations = 0, matvecs = 0, m_active = 0;

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  const int status = native_block_lanczos_run(
    &impl, eigencore_dense_apply, n, k, m_max, block_size, target_kind, tol,
    REAL(start_), V.data(), lambda.data(), residuals.data(), converged.data(),
    &nconv, &iterations, &matvecs, &m_active);
  if (status != 0) {
    error("native dense block Lanczos failed with status=%d", status);
  }
  return block_lanczos_pack_result(n, k, V.data(), lambda.data(),
                                   residuals.data(), converged.data(),
                                   nconv, iterations, matvecs, m_active);
}

extern "C" SEXP eigencore_block_lanczos_csc(SEXP i_, SEXP p_, SEXP x_,
                                            SEXP dim_, SEXP k_,
                                            SEXP m_max_,
                                            SEXP block_size_,
                                            SEXP target_kind_,
                                            SEXP tol_,
                                            SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC block Lanczos inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable CSC block Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k) error("m_max must be >= k");
  if (m_max > n) error("m_max must be <= nrow(A)");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int nconv = 0, iterations = 0, matvecs = 0, m_active = 0;

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_block_lanczos_run(
    &impl, eigencore_csc_apply, n, k, m_max, block_size, target_kind, tol,
    REAL(start_), V.data(), lambda.data(), residuals.data(), converged.data(),
    &nconv, &iterations, &matvecs, &m_active);
  if (status != 0) {
    error("native CSC block Lanczos failed with status=%d", status);
  }
  return block_lanczos_pack_result(n, k, V.data(), lambda.data(),
                                   residuals.data(), converged.data(),
                                   nconv, iterations, matvecs, m_active);
}

extern "C" SEXP eigencore_block_thick_restart_lanczos_dense(
    SEXP A_, SEXP k_, SEXP m_max_, SEXP block_size_,
    SEXP target_kind_, SEXP tol_, SEXP max_restarts_,
    SEXP norm_a_, SEXP start_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable block thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  const double norm_a = asReal(norm_a_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, restarts = 0, m_active = 0;
  int locking_events = 0, ortho_passes = 0;
  int64_t operator_allocations = 0, operator_bytes_allocated = 0;
  NativeBlockStageSeconds stage_seconds;
  const int history_capacity = max_restarts + 1;
  std::vector<int> history_restart(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_m_active(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_selected_count(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_before(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_after(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_nconv_wanted(static_cast<size_t>(history_capacity), 0);
  std::vector<double> history_max_residual(static_cast<size_t>(history_capacity), R_PosInf);
  std::vector<double> history_max_backward_error(static_cast<size_t>(history_capacity), R_PosInf);
  NativeBlockRestartHistory history;
  history.capacity = history_capacity;
  history.restart = history_restart.data();
  history.m_active = history_m_active.data();
  history.selected_count = history_selected_count.data();
  history.locked_before = history_locked_before.data();
  history.locked_after = history_locked_after.data();
  history.nconv_wanted = history_nconv_wanted.data();
  history.max_residual = history_max_residual.data();
  history.max_backward_error = history_max_backward_error.data();

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  const int status = native_block_thick_restart_lanczos_run(
    &impl, eigencore_dense_apply, n, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, 0, REAL(start_), V.data(), lambda.data(),
    residuals.data(), converged.data(), &n_locked, &iterations, &matvecs,
    &restarts, &m_active, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes_allocated, &stage_seconds, &history);
  if (status != 0) {
    error("native dense block thick-restart Lanczos failed with status=%d", status);
  }
  return block_thick_lanczos_pack_result(
    n, k, V.data(), lambda.data(), residuals.data(), converged.data(),
    n_locked, iterations, matvecs, restarts, m_active, locking_events,
    ortho_passes, block_size, operator_allocations, operator_bytes_allocated,
    &stage_seconds, &history);
}
extern "C" SEXP eigencore_block_thick_restart_lanczos_csc(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP k_,
    SEXP m_max_, SEXP block_size_, SEXP target_kind_,
    SEXP tol_, SEXP max_restarts_, SEXP norm_a_, SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC block thick-restart Lanczos inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable CSC block thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  const double norm_a = asReal(norm_a_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, restarts = 0, m_active = 0;
  int locking_events = 0, ortho_passes = 0;
  int64_t operator_allocations = 0, operator_bytes_allocated = 0;
  NativeBlockStageSeconds stage_seconds;
  const int history_capacity = max_restarts + 1;
  std::vector<int> history_restart(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_m_active(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_selected_count(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_before(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_after(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_nconv_wanted(static_cast<size_t>(history_capacity), 0);
  std::vector<double> history_max_residual(static_cast<size_t>(history_capacity), R_PosInf);
  std::vector<double> history_max_backward_error(static_cast<size_t>(history_capacity), R_PosInf);
  NativeBlockRestartHistory history;
  history.capacity = history_capacity;
  history.restart = history_restart.data();
  history.m_active = history_m_active.data();
  history.selected_count = history_selected_count.data();
  history.locked_before = history_locked_before.data();
  history.locked_after = history_locked_after.data();
  history.nconv_wanted = history_nconv_wanted.data();
  history.max_residual = history_max_residual.data();
  history.max_backward_error = history_max_backward_error.data();

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_block_thick_restart_lanczos_run(
    &impl, eigencore_csc_apply, n, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, 1, REAL(start_), V.data(), lambda.data(),
    residuals.data(), converged.data(), &n_locked, &iterations, &matvecs,
    &restarts, &m_active, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes_allocated, &stage_seconds, &history);
  if (status != 0) {
    error("native CSC block thick-restart Lanczos failed with status=%d", status);
  }
  return block_thick_lanczos_pack_result(
    n, k, V.data(), lambda.data(), residuals.data(), converged.data(),
    n_locked, iterations, matvecs, restarts, m_active, locking_events,
    ortho_passes, block_size, operator_allocations, operator_bytes_allocated,
    &stage_seconds, &history);
}
