#include <cmath>
#include <cfloat>
#include <cstdint>
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
  return ec_norm2(x, n);
}

// Frobenius norm of a dense n x cols block stored contiguously (ld == n).
static double block_frobenius_norm(const double* X, int n, int cols) {
  if (cols == 1) {
    return ec_norm2(X, n);
  }
  double scale = 0.0;
  double ssq = 1.0;
  for (int col = 0; col < cols; ++col) {
    const double c = ec_norm2(X + static_cast<int64_t>(col) * n, n);
    if (c > 0.0) {
      if (scale < c) {
        ssq = 1.0 + ssq * (scale / c) * (scale / c);
        scale = c;
      } else {
        ssq += (c / scale) * (c / scale);
      }
    } else if (ISNAN(c)) {
      return c;
    }
  }
  return scale * sqrt(ssq);
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

// Selected eigenpairs of the m x m symmetric projected matrix A (upper
// triangle referenced, ld m; destroyed; T/ldt is the source it was copied
// from, used to restart the solve) (P5). The algebraic targets only need
// the `count` eigenvectors at one end of the spectrum, so dsyevr runs with an
// index range instead of a full dsyev/dsyevd; magnitude targets need both ends
// and compute every pair (dsyevr, MRRR). Indices keep their meaning as
// positions in the ascending spectrum: theta[idx] is set for every selected
// idx (other entries are NaN and never read), selected[0..count) is the
// target ordering, and column p of S_selected (ld m) holds the eigenvector of
// selected[p].
static int projected_eigen_selected(const double* T, int ldt,
                                    double* A, int m, int count,
                                    int target_kind, double* theta,
                                    int* selected, double* S_selected,
                                    double* Z, double* w, int* isuppz,
                                    double* work, int lwork,
                                    int* iwork, int liwork) {
  if (m <= 0) {
    return 0;
  }
  if (count > m) {
    count = m;
  }
  if (count < 1) {
    count = 1;
  }
  char jobz = 'V';
  char uplo = 'U';
  char range = 'A';
  int il = 1;
  int iu = m;
  if (target_kind == 1 && count < m) {
    range = 'I';
    il = m - count + 1;
  } else if (target_kind == 2 && count < m) {
    range = 'I';
    iu = count;
  }
  const double vl = 0.0;
  const double vu = 0.0;
  const double abstol = 0.0;
  int found = 0;
  int info = 0;
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &m, A, &m, &vl, &vu, &il, &iu,
                   &abstol, &found, w, Z, &m, isuppz, work, &lwork,
                   iwork, &liwork, &info FCONE FCONE FCONE);
  if (info == 0 && range == 'I' && found != iu - il + 1) {
    // Bisection can return extra eigenvalues tied at a range boundary; the
    // index mapping is then ambiguous, so solve for the whole spectrum.
    for (int col = 0; col < m; ++col) {
      for (int row = 0; row <= col; ++row) {
        A[row + static_cast<int64_t>(col) * m] = T[row + static_cast<int64_t>(col) * ldt];
      }
    }
    range = 'A';
    il = 1;
    iu = m;
    found = 0;
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &m, A, &m, &vl, &vu, &il, &iu,
                     &abstol, &found, w, Z, &m, isuppz, work, &lwork,
                     iwork, &liwork, &info FCONE FCONE FCONE);
  }
  if (info != 0 || found != iu - il + 1) {
    return -3;
  }
  for (int i = 0; i < m; ++i) {
    theta[i] = R_NaN;
  }
  for (int i = 0; i < found; ++i) {
    theta[il - 1 + i] = w[i];
  }
  selected_sorted_ritz_indices(theta, m, count, target_kind, selected);
  for (int p = 0; p < count; ++p) {
    const int zcol = selected[p] - (il - 1);
    std::memcpy(S_selected + static_cast<int64_t>(p) * m,
                Z + static_cast<int64_t>(zcol) * m,
                sizeof(double) * static_cast<size_t>(m));
  }
  return 0;
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

// Ritz/locked vectors are unit vectors, so the norm thresholds below are
// relative (|v| is compared against 1, independent of the operator scale).
static int vector_is_independent_from_locked(const double* V_locked,
                                             int n_locked,
                                             const double* v,
                                             int n) {
  const double dot_tol = 10.0 * sqrt(DBL_EPSILON);
  const double vnorm = trl_norm2(v, n);
  if (vnorm <= DBL_EPSILON) {
    return 0;
  }
  if (n_locked <= 0) {
    return 1;
  }
  std::vector<double> dots(static_cast<size_t>(n_locked), 0.0);
  const char trans_T = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  int inc = 1;
  F77_CALL(dgemv)(&trans_T, &n, &n_locked, &one, V_locked, &n, v, &inc,
                  &zero, dots.data(), &inc FCONE);
  for (int col = 0; col < n_locked; ++col) {
    const double locked_norm =
      trl_norm2(V_locked + static_cast<int64_t>(col) * n, n);
    if (locked_norm <= DBL_EPSILON) {
      continue;
    }
    if (fabs(dots[static_cast<size_t>(col)]) > dot_tol * locked_norm * vnorm) {
      return 0;
    }
  }
  return 1;
}

// Orthogonalise z against the locked and active bases and append it when it
// survives. Breakdown is judged relative to ref_norm (the scale z had before
// orthogonalisation, or an operator-norm estimate for residual directions), so
// the decision is invariant under scaling of the operator (C16).
static int block_accept_work_vector(const double* V_locked, int n_locked,
                                    double* V_active, int* m_active,
                                    int m_max, double* z, double* tmp,
                                    int n, int* ortho_passes,
                                    double ref_norm = 1.0) {
  if (*m_active >= m_max) {
    return 0;
  }
  const int passes_done =
    trl_orthogonalise(V_locked, n_locked, V_active, *m_active, z, tmp, n, 2);
  if (ortho_passes != nullptr) {
    *ortho_passes += passes_done;
  }
  const double nz = trl_norm2(z, n);
  if (!(ref_norm > 0.0) || !R_FINITE(ref_norm)) {
    ref_norm = 1.0;
  }
  if (!(nz > 100.0 * DBL_EPSILON * ref_norm)) {
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
                                      bool reorthogonalize_active = true,
                                      double ref_norm = -1.0,
                                      bool reorthogonalize_locked = true) {
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

  // Breakdown thresholds are relative (C16): without a caller-supplied scale
  // the reference is the largest column norm of the block before
  // orthogonalisation, so a block scaled by 1e-12 behaves like one at scale 1.
  if (!(ref_norm > 0.0) || !R_FINITE(ref_norm)) {
    ref_norm = 0.0;
    for (int col = 0; col < cols; ++col) {
      const double c = trl_norm2(Z_block + static_cast<int64_t>(col) * n, n);
      if (c > ref_norm) {
        ref_norm = c;
      }
    }
    if (!(ref_norm > 0.0) || !R_FINITE(ref_norm)) {
      ref_norm = 1.0;
    }
  }
  const double breakdown_tol = 100.0 * DBL_EPSILON * ref_norm;
  // R_ii of the Cholesky factor scales like the column norms, i.e. like
  // ref_norm; the Gram diagonal (before dpotrf) like ref_norm^2.
  const double* active_basis = reorthogonalize_active ? V_active : nullptr;
  const int active_cols = reorthogonalize_active ? *m_active : 0;
  const int reorth_passes_done = block_reorthogonalise_against(
    reorthogonalize_locked ? V_locked : nullptr,
    reorthogonalize_locked ? n_locked : 0, active_basis, active_cols,
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
      if (!(coeff[col + static_cast<int64_t>(col) * cols] > breakdown_tol)) {
        chol_ok = false;
      }
    }
    if (chol_ok) {
      F77_CALL(dtrsm)(&right, &uplo, &trans_N, &diag, &n, &cols, &one,
                      coeff, &cols, Z_block, &n FCONE FCONE FCONE FCONE);
      // Second CholQR pass (CholQR2): a single pass loses orthogonality
      // like cond(Z)^2 * eps, so always repeat it regardless of n.
      F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                      &one, Z_block, &n, Z_block, &n,
                      &zero, coeff, &cols FCONE FCONE);
      symmetrize_packed_square(coeff, cols);
      F77_CALL(dpotrf)(&uplo, &cols, coeff, &cols, &info FCONE);
      // Second pass: Z_block is now (nearly) orthonormal, so R_ii ~ 1.
      chol_ok = (info == 0);
      for (int col = 0; chol_ok && col < cols; ++col) {
        if (!(coeff[col + static_cast<int64_t>(col) * cols] > 100.0 * DBL_EPSILON)) {
          chol_ok = false;
        }
      }
      if (chol_ok) {
        F77_CALL(dtrsm)(&right, &uplo, &trans_N, &diag, &n, &cols, &one,
                        coeff, &cols, Z_block, &n FCONE FCONE FCONE FCONE);
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
    if (!(nz > breakdown_tol)) {
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
                              int* matvecs_out,
                              int* operator_columns_out) {
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
  if (rc == 0 && operator_columns_out != nullptr) {
    *operator_columns_out += cols;
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
  double* Z_eig;       // m_max x m_max: selected projected eigenvectors (dsyevr)
  double* w_eig;       // m_max: projected eigenvalues (dsyevr)
  int*    isuppz;      // 2 * m_max (dsyevr)
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
  std::free(b->Z_eig);
  std::free(b->w_eig);
  std::free(b->isuppz);
  std::memset(b, 0, sizeof(*b));
}

struct TrlBuffersGuard {
  explicit TrlBuffersGuard(ThickRestartBuffers* buffers) : b(buffers) {}
  TrlBuffersGuard(const TrlBuffersGuard&) = delete;
  TrlBuffersGuard& operator=(const TrlBuffersGuard&) = delete;
  ~TrlBuffersGuard() { trl_buffers_free(b); }
  ThickRestartBuffers* b;
};

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
  b->Z_eig = static_cast<double*>(std::malloc((mm > 0 ? mm : 1) * sizeof(double)));
  b->w_eig = static_cast<double*>(std::malloc(static_cast<size_t>(m_max > 0 ? m_max : 1) * sizeof(double)));
  b->isuppz = static_cast<int*>(std::malloc(static_cast<size_t>(2 * (m_max > 0 ? m_max : 1)) * sizeof(int)));
  if (b->V_active == nullptr || b->AV_active == nullptr ||
      b->T_proj == nullptr || b->S_eig == nullptr ||
      b->S_selected == nullptr ||
      b->theta == nullptr || b->B_v == nullptr || b->B_av == nullptr ||
      b->Z_block == nullptr || b->coeff_block == nullptr ||
      b->z == nullptr || b->tmp == nullptr || b->ritz_res == nullptr ||
      b->selected == nullptr || b->is_locked == nullptr ||
      b->dsyev_work == nullptr || b->dsyevd_iwork == nullptr ||
      b->Z_eig == nullptr || b->w_eig == nullptr || b->isuppz == nullptr) {
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
                                   int* matvecs_out,
                                   int* operator_columns_out,
                                   int* certification_operator_columns_out) {
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
  if (operator_columns_out != nullptr) {
    *operator_columns_out += k_target;
  }
  if (certification_operator_columns_out != nullptr) {
    *certification_operator_columns_out += k_target;
  }

  int n_converged_simple = 0;
  for (int col = 0; col < k_target; ++col) {
    const double* vec = V_out + static_cast<int64_t>(col) * n;
    const double* av = buf->B_av + static_cast<int64_t>(col) * n;
    const double theta = ec_dot(vec, av, n);
    const double res = ec_residual_norm(av, vec, theta, n);
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
  combine_basis_columns(V_out, n, k_target, buf->S_selected, k_target,
                        k_target, buf->B_v);
  combine_basis_columns(buf->AV_active, n, k_target, buf->S_selected,
                        k_target, k_target, buf->B_av);

  int n_converged = 0;
  for (int col = 0; col < k_target; ++col) {
    const int idx = buf->selected[col];
    const double theta = buf->theta[idx];
    const double* residual = buf->B_av + static_cast<int64_t>(col) * n;
    const double* vec = buf->B_v + static_cast<int64_t>(col) * n;
    const double res = ec_residual_norm(residual, vec, theta, n);
    std::memcpy(V_out + static_cast<int64_t>(col) * n, vec,
                sizeof(double) * static_cast<size_t>(n));
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
  std::vector<double> V_storage(eigencore_buffer_size(nm));
  double* V = V_storage.data();
  std::vector<double> AV_storage(eigencore_buffer_size(nm));
  double* AV = AV_storage.data();
  std::vector<double> Z_storage(eigencore_buffer_size(nb));
  double* Z = Z_storage.data();
  std::vector<double> AZ_storage(eigencore_buffer_size(nb));
  double* AZ = AZ_storage.data();
  std::vector<double> H_storage(eigencore_buffer_size(mm));
  double* H = H_storage.data();
  std::vector<double> S_selected_storage(eigencore_buffer_size(mm));
  double* S_selected = S_selected_storage.data();
  std::vector<double> theta_storage(eigencore_buffer_size(static_cast<size_t>(m_max)));
  double* theta = theta_storage.data();
  std::vector<double> B_v_storage(eigencore_buffer_size(static_cast<size_t>(n) * k_target));
  double* B_v = B_v_storage.data();
  std::vector<double> B_av_storage(eigencore_buffer_size(static_cast<size_t>(n) * k_target));
  double* B_av = B_av_storage.data();
  std::vector<double> tmp_storage(eigencore_buffer_size(static_cast<size_t>(m_max)));
  double* tmp = tmp_storage.data();
  std::vector<int> selected_storage(eigencore_buffer_size(static_cast<size_t>(m_max)));
  int* selected = selected_storage.data();
  const int dsyev_lwork_query = trl_dsyev_query(m_max);
  int dsyev_lwork = dsyev_lwork_query > 0 ? dsyev_lwork_query : 3 * m_max;
  std::vector<double> dsyev_work_storage(eigencore_buffer_size(static_cast<size_t>(dsyev_lwork)));
  double* dsyev_work = dsyev_work_storage.data();

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
    eigencore_check_interrupt();
    int accepted_start = m_active;
    int accepted = 0;
    for (int col = 0; col < source_cols && m_active < m_max; ++col) {
      double* z_col = Z + static_cast<int64_t>(col) * n;
      // Breakdown relative to the vector's own scale before orthogonalisation
      // (C16): A*v for a matrix scaled by 1e-12 is not a breakdown.
      const double pre_nz = trl_norm2(z_col, n);
      trl_orthogonalise(nullptr, 0, V, m_active, z_col, tmp, n);
      const double nz = trl_norm2(z_col, n);
      if (!(nz > 100.0 * DBL_EPSILON * pre_nz)) {
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
  // Convergence relative to an operator-norm estimate (C16): ||T|| = max
  // |theta| over the whole projected spectrum, never an absolute floor of 1.
  double norm_est = 0.0;
  for (int i = 0; i < m_active; ++i) {
    if (fabs(theta[i]) > norm_est) {
      norm_est = fabs(theta[i]);
    }
  }
  int nconv = 0;
  for (int col = 0; col < k_target; ++col) {
    residuals_out[col] = ec_residual_norm(B_av + static_cast<int64_t>(col) * n,
                                          B_v + static_cast<int64_t>(col) * n,
                                          lambda_out[col], n);
    const double scale_i = (fabs(lambda_out[col]) > norm_est) ?
      fabs(lambda_out[col]) : norm_est;
    converged_out[col] = residuals_out[col] <= tol * scale_i ? 1 : 0;
    if (converged_out[col]) ++nconv;
    std::memcpy(V_out + static_cast<int64_t>(col) * n,
                B_v + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  *nconv_out = nconv;
  *m_active_final_out = m_active;
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

// Write the symmetric pair T(row, col) = T(col, row) = value.
static inline void projection_set_pair(double* T_proj, int ldt, int row, int col,
                                       double value) {
  T_proj[row + static_cast<int64_t>(col) * ldt] = value;
  T_proj[col + static_cast<int64_t>(row) * ldt] = value;
}

// Lanczos residual of the most recent block with the projected column built in
// the same sweep (P2/P4). On entry AV_active holds A * V_last. On exit W holds
//   W = A V_last - V_active H - V_locked (V_locked' W),
// orthogonal to the locked and active bases, and T_proj holds the last block's
// column H = V_active' A V_last (and its mirror row).
//
// When the column is not yet known (the normal case) the local Lanczos
// coefficients -- the coupling to the previous block and the self block -- are
// formed explicitly and subtracted (the three-term recurrence); the remaining
// entries of H fall out of one classical Gram-Schmidt pass of the recurrence
// residual against the whole basis, which is also the full reorthogonalisation
// (H = local + V'W_local exactly when V is orthonormal). A DGKS second pass
// runs only when the first cancelled a large fraction of the residual norm.
// Previously the column was a separate dense V'AV projection per step on top
// of the reorthogonalisation; folding it in saves one n x m pass per step.
//
// When the column is already explicit (mid-sweep checkpoint, or the explicit
// restart fallback), W = A V_last - V_active H is formed directly and
// reorthogonalised with the usual DGKS scheme; H is left untouched.
//
// Returns the number of full projection passes performed.
static int block_lanczos_projected_residual(const double* V_locked, int n_locked,
                                            ThickRestartBuffers* buf, int n,
                                            int m_max, int m_active,
                                            int prev_start, int prev_cols,
                                            int last_start, int last_cols,
                                            bool last_column_known,
                                            double* W) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  double* T = buf->T_proj;
  const int ldt = m_max;
  double* coeff = buf->coeff_block;
  const int b = last_cols;
  for (int col = 0; col < b; ++col) {
    std::memcpy(W + static_cast<int64_t>(col) * n,
                buf->AV_active + static_cast<int64_t>(last_start + col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  if (last_column_known) {
    subtract_projected_range(buf->V_active, T, ldt, n, 0, m_active,
                             last_start, b, W, coeff);
    return block_reorthogonalise_against(V_locked, n_locked, buf->V_active,
                                         m_active, W, n, b, coeff, 2);
  }

  // Local coefficients: coupling to the previous block and the self block.
  const double* V_last = buf->V_active + static_cast<int64_t>(last_start) * n;
  const double* AV_last = buf->AV_active + static_cast<int64_t>(last_start) * n;
  double* scratch = buf->S_eig;  // free during expansion; m_max x m_max
  if (prev_cols > 0) {
    F77_CALL(dgemm)(&trans_T, &trans_N, &prev_cols, &b, &n, &one,
                    buf->V_active + static_cast<int64_t>(prev_start) * n, &n,
                    AV_last, &n, &zero, scratch, &prev_cols FCONE FCONE);
    for (int col = 0; col < b; ++col) {
      for (int row = 0; row < prev_cols; ++row) {
        projection_set_pair(T, ldt, prev_start + row, last_start + col,
                            scratch[row + static_cast<int64_t>(col) * prev_cols]);
      }
    }
  }
  F77_CALL(dgemm)(&trans_T, &trans_N, &b, &b, &n, &one, V_last, &n,
                  AV_last, &n, &zero, scratch, &b FCONE FCONE);
  symmetrize_packed_square(scratch, b);
  for (int col = 0; col < b; ++col) {
    for (int row = 0; row < b; ++row) {
      T[(last_start + row) + static_cast<int64_t>(last_start + col) * ldt] =
        scratch[row + static_cast<int64_t>(col) * b];
    }
  }
  subtract_projected_range(buf->V_active, T, ldt, n, prev_start, prev_cols,
                           last_start, b, W, coeff);
  subtract_projected_range(buf->V_active, T, ldt, n, last_start, b,
                           last_start, b, W, coeff);

  // CGS pass(es) against [locked | active]; active coefficients complete H.
  double pre_norm = block_frobenius_norm(W, n, b);
  int passes_done = 0;
  for (int pass = 0; pass < 2; ++pass) {
    if (n_locked > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &n_locked, &b, &n, &one,
                      V_locked, &n, W, &n, &zero, coeff, &n_locked FCONE FCONE);
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &b, &n_locked, &minus_one,
                      V_locked, &n, coeff, &n_locked, &one, W, &n FCONE FCONE);
    }
    if (m_active > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &m_active, &b, &n, &one,
                      buf->V_active, &n, W, &n, &zero, coeff, &m_active FCONE FCONE);
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &b, &m_active, &minus_one,
                      buf->V_active, &n, coeff, &m_active, &one, W, &n FCONE FCONE);
      for (int col = 0; col < b; ++col) {
        const int abs_col = last_start + col;
        for (int row = 0; row < m_active; ++row) {
          const double c = coeff[row + static_cast<int64_t>(col) * m_active];
          if (row >= last_start && row < last_start + b) {
            // Self block: accumulate on the upper entry only, symmetrised below.
            T[row + static_cast<int64_t>(abs_col) * ldt] += c;
          } else {
            const double value = T[row + static_cast<int64_t>(abs_col) * ldt] + c;
            projection_set_pair(T, ldt, row, abs_col, value);
          }
        }
      }
    }
    ++passes_done;
    if (pass + 1 >= 2) {
      break;
    }
    const double post_norm = block_frobenius_norm(W, n, b);
    if (post_norm >= kDgksEta * pre_norm) {
      break;
    }
    pre_norm = post_norm;
  }
  for (int col = 0; col < b; ++col) {
    for (int row = col + 1; row < b; ++row) {
      const int r = last_start + row;
      const int c = last_start + col;
      const double avg = 0.5 * (T[r + static_cast<int64_t>(c) * ldt] +
                                T[c + static_cast<int64_t>(r) * ldt]);
      projection_set_pair(T, ldt, r, c, avg);
    }
  }
  return passes_done;
}

// Largest column norm of the n x cols block X (ld == n).
static double block_max_column_norm(const double* X, int n, int cols) {
  double out = 0.0;
  for (int col = 0; col < cols; ++col) {
    const double c = trl_norm2(X + static_cast<int64_t>(col) * n, n);
    if (c > out) {
      out = c;
    }
  }
  return out;
}

static int block_lanczos_expand_basis_to_budget(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int m_max,
    int m_stop,
    int block_size,
    const double* V_locked,
    int n_locked,
    double norm_scale,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    int* m_active,
    int* previous_block_start,
    int* previous_block_cols,
    int* last_block_start,
    int* last_block_cols,
    bool* last_column_known,
    int* iterations_out,
    int* matvecs_out,
    int* operator_columns_out,
    int* ortho_passes_out) {
  // m_stop bounds this expansion pass; buffers (V_active, AV_active, T_proj)
  // remain sized by m_max. With mid-sweep checks (m_stop < m_max) the driver
  // expands one chunk, evaluates convergence, then re-enters to continue the
  // same sweep. Legacy behaviour is m_stop == m_max (a single full sweep).
  while (*m_active < m_stop && *last_block_cols > 0) {
    eigencore_check_interrupt();
    auto timer = native_timer_now();
    const double av_norm = block_max_column_norm(
      buf->AV_active + static_cast<int64_t>(*last_block_start) * n, n,
      *last_block_cols);
    // Breakdown is judged against the operator scale (C16): the norm estimate
    // when known, and at least |A v| for the current block.
    double ref_norm = (norm_scale > av_norm) ? norm_scale : av_norm;
    if (!(ref_norm > 0.0) || !R_FINITE(ref_norm)) {
      ref_norm = 1.0;
    }
    const int residual_passes = block_lanczos_projected_residual(
      V_locked, n_locked, buf, n, m_max, *m_active,
      *previous_block_start, *previous_block_cols,
      *last_block_start, *last_block_cols, *last_column_known, buf->Z_block);
    *last_column_known = true;
    if (ortho_passes_out != nullptr) {
      *ortho_passes_out += residual_passes;
    }
    stages->reorthogonalization += native_timer_elapsed(timer);

    const int accepted_start = *m_active;
    timer = native_timer_now();
    const int accepted = block_accept_columns_blas3(
      buf->Z_block, n, *last_block_cols, V_locked, n_locked,
      buf->V_active, m_active, m_max, buf->Z_block, block_size,
      buf->coeff_block, buf->tmp, n,
      block_size, ortho_passes_out,
      false, ref_norm, false
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
          buf->z, buf->tmp, n, ortho_passes_out, 1.0
        );
        stages->reorthogonalization += native_timer_elapsed(timer);
      }
      if (continuation_accepted == 0) {
        break;
      }

      timer = native_timer_now();
      const int rc = apply_active_block(
        impl, apply, n, continuation_start, continuation_accepted,
        buf->V_active, buf->AV_active, workspace, matvecs_out,
        operator_columns_out
      );
      stages->apply += native_timer_elapsed(timer);
      if (rc != 0) {
        return rc;
      }
      ++(*iterations_out);
      *previous_block_start = *last_block_start;
      *previous_block_cols = *last_block_cols;
      *last_block_start = continuation_start;
      *last_block_cols = continuation_accepted;
      *last_column_known = false;
      continue;
    }

    timer = native_timer_now();
    const int rc = apply_active_block(impl, apply, n, accepted_start, accepted,
                                      buf->V_active, buf->AV_active, workspace,
                                      matvecs_out, operator_columns_out);
    stages->apply += native_timer_elapsed(timer);
    if (rc != 0) {
      return rc;
    }
    ++(*iterations_out);
    *previous_block_start = *last_block_start;
    *previous_block_cols = *last_block_cols;
    *last_block_start = accepted_start;
    *last_block_cols = accepted;
    *last_column_known = false;
  }

  // The Rayleigh-Ritz step needs the newest block's projected column; it is
  // only formed by the next expansion step, so build it explicitly here.
  if (!*last_column_known && *last_block_cols > 0) {
    auto timer = native_timer_now();
    projection_update_appended_block(buf->T_proj, m_max, buf->V_active,
                                     buf->AV_active, n, *last_block_start,
                                     *last_block_cols, buf->S_eig);
    *last_column_known = true;
    const double elapsed = native_timer_elapsed(timer);
    stages->projected_solve += elapsed;
    stages->projection_update += elapsed;
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

// Thick restart (Krylov-Schur style). The kept Ritz vectors Y = V S and their
// images A Y are already available from the Rayleigh-Ritz step (B_v, B_av), so
// they are copied straight into the new basis and the kept block of the
// projected matrix is diag(theta) (P2). Only the continuation tail -- the
// residual direction -- is applied to the operator; its coupling column
// Y' A t = (A Y)' t is formed by the next expansion step. This replaces
// re-orthogonalising every kept vector one at a time, re-applying A to all of
// them and recomputing the full V'AV.
//
// The reuse is only taken when the kept vectors are orthonormal and orthogonal
// to the locked set to working accuracy (they are by construction: V_active is
// orthonormal and deflated against the locked vectors, and S is orthonormal);
// a cheap Gram check guards it, and any deviation falls back to the explicit
// re-orthogonalise-and-apply path.
static const double kRestartReuseOrthTol = 1e-12;

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
    double norm_scale,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    int* m_active,
    int* previous_block_start,
    int* previous_block_cols,
    int* last_block_start,
    int* last_block_cols,
    bool* last_column_known,
    double* V_out,
    int* matvecs_out,
    int* operator_columns_out,
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
  std::vector<int> keep;
  keep.reserve(static_cast<size_t>(k_keep > 0 ? k_keep : 0));
  for (int p = 0; p < selected_count && static_cast<int>(keep.size()) < k_keep; ++p) {
    if (!buf->is_locked[p]) {
      keep.push_back(p);
    }
  }
  const int kk = static_cast<int>(keep.size());
  const double scale = (norm_scale > 0.0 && R_FINITE(norm_scale)) ? norm_scale : 1.0;

  std::memset(buf->T_proj, 0,
              sizeof(double) * static_cast<size_t>(m_max) *
                static_cast<size_t>(m_max));
  *m_active = 0;

  for (int i = 0; i < kk; ++i) {
    const int p = keep[static_cast<size_t>(i)];
    std::memcpy(buf->V_active + static_cast<int64_t>(i) * n,
                buf->B_v + static_cast<int64_t>(p) * n,
                sizeof(double) * static_cast<size_t>(n));
    std::memcpy(buf->AV_active + static_cast<int64_t>(i) * n,
                buf->B_av + static_cast<int64_t>(p) * n,
                sizeof(double) * static_cast<size_t>(n));
  }
  bool reuse = true;
  if (kk > 0) {
    const char trans_T = 'T';
    const char trans_N = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    double deviation = 0.0;
    F77_CALL(dgemm)(&trans_T, &trans_N, &kk, &kk, &n, &one,
                    buf->V_active, &n, buf->V_active, &n,
                    &zero, buf->S_eig, &kk FCONE FCONE);
    for (int col = 0; col < kk; ++col) {
      for (int row = 0; row < kk; ++row) {
        const double expected = (row == col) ? 1.0 : 0.0;
        const double d = fabs(buf->S_eig[row + static_cast<int64_t>(col) * kk] - expected);
        if (!(d <= deviation)) {
          deviation = d;
        }
      }
    }
    if (n_locked > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &n_locked, &kk, &n, &one,
                      V_out, &n, buf->V_active, &n,
                      &zero, buf->S_eig, &n_locked FCONE FCONE);
      const int64_t total = static_cast<int64_t>(n_locked) * kk;
      for (int64_t i = 0; i < total; ++i) {
        const double d = fabs(buf->S_eig[i]);
        if (!(d <= deviation)) {
          deviation = d;
        }
      }
    }
    reuse = (deviation <= kRestartReuseOrthTol);
  }

  if (reuse) {
    *m_active = kk;
    for (int i = 0; i < kk; ++i) {
      const int p = keep[static_cast<size_t>(i)];
      buf->T_proj[i + static_cast<int64_t>(i) * m_max] =
        buf->theta[buf->selected[p]];
    }
  } else {
    for (int i = 0; i < kk && *m_active < m_max; ++i) {
      const int p = keep[static_cast<size_t>(i)];
      std::memcpy(buf->z, buf->B_v + static_cast<int64_t>(p) * n,
                  sizeof(double) * static_cast<size_t>(n));
      stages->restart += native_timer_elapsed(timer);
      timer = native_timer_now();
      block_accept_work_vector(V_out, n_locked, buf->V_active, m_active, m_max,
                               buf->z, buf->tmp, n, ortho_passes_out, 1.0);
      stages->reorthogonalization += native_timer_elapsed(timer);
      timer = native_timer_now();
    }
  }

  const int tail_start = *m_active;
  int tail_accepted = 0;
  for (int p = 0; p < selected_count && tail_accepted < block_size && *m_active < m_max; ++p) {
    if (buf->is_locked[p] || !(buf->ritz_res[p] > 100.0 * DBL_EPSILON * scale)) {
      continue;
    }
    const int idx = buf->selected[p];
    for (int row = 0; row < n; ++row) {
      buf->z[row] = buf->B_av[static_cast<int64_t>(p) * n + row] -
        buf->theta[idx] * buf->B_v[static_cast<int64_t>(p) * n + row];
    }
    const double z_norm = trl_norm2(buf->z, n);
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    tail_accepted += block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out, z_norm
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
      buf->z, buf->tmp, n, ortho_passes_out, 1.0
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
  }
  stages->restart += native_timer_elapsed(timer);
  if (tail_accepted == 0) {
    return 1;
  }

  timer = native_timer_now();
  const int apply_start = reuse ? tail_start : 0;
  int rc = apply_active_block(impl, apply, n, apply_start, *m_active - apply_start,
                              buf->V_active, buf->AV_active, workspace,
                              matvecs_out, operator_columns_out);
  stages->apply += native_timer_elapsed(timer);
  if (rc != 0) {
    return rc;
  }

  if (reuse) {
    *last_column_known = false;
  } else {
    timer = native_timer_now();
    projection_update_self_block(buf->T_proj, m_max, buf->V_active,
                                 buf->AV_active, n, 0, *m_active,
                                 buf->S_eig);
    *last_column_known = true;
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
    int* matvecs_out,
    int* operator_columns_out,
    int* certification_operator_columns_out) {
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
      &polished_converged, buf, workspace, matvecs_out,
      operator_columns_out, certification_operator_columns_out
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
      &polished_converged, buf, workspace, matvecs_out,
      operator_columns_out, certification_operator_columns_out
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

// Project a vector orthogonal to a basis D (n x d_cols, orthonormal columns)
// with `passes` Gram-Schmidt sweeps. Used to build and maintain the deflated
// complement in block_lanczos_window_complement_clean.
static void block_lanczos_project_out(double* v, int n, const double* D,
                                      int d_cols, int passes) {
  int inc = 1;
  for (int pass = 0; pass < passes; ++pass) {
    for (int c = 0; c < d_cols; ++c) {
      const double* d = D + static_cast<int64_t>(c) * n;
      const double dotd = -ec_dot(d, v, n);
      F77_CALL(daxpy)(&n, &dotd, d, &inc, v, &inc);
    }
  }
}

// Does a complement extreme (d_min, d_max) beat the window edge for this target?
static bool block_lanczos_complement_intruder(int target_kind, double d_min,
                                              double d_max, double window_edge,
                                              double norm_a, double tol) {
  // Relative to the operator scale (C16): no absolute floor.
  const double scale = std::fmax(std::fmax(std::fabs(window_edge), std::fabs(d_max)),
                                 std::fmax(std::fabs(d_min), norm_a));
  const double margin = std::fmax(1e-6, 10.0 * tol) * scale;
  if (target_kind == 1) {                     // largest algebraic
    return d_max > window_edge + margin;
  } else if (target_kind == 2) {              // smallest algebraic
    return d_min < window_edge - margin;
  }
  const double comp_mag = std::fmax(std::fabs(d_min), std::fabs(d_max));
  return comp_mag > std::fabs(window_edge) + margin;  // largest magnitude
}

// Deflated-complement confirmation for a candidate mid-sweep window. The window
// (n_locked locked columns in V_out plus the `wanted` window Ritz vectors at the
// front of window_vecs) is projected out of A, and short deterministic-seed
// Lanczos segments on the deflated operator estimate the operator's extreme
// eigenvalues on that complement. If an extreme lies past the window edge -- a
// more-preferred value than the least-preferred window member for the active
// target -- the window has truncated a target direction it does not hold. That
// happens for an unresolved near-degenerate copy: its eigenvector is orthogonal
// to the (resolved) window, so it survives deflation and reappears here even
// though the un-deflated workspace already (weakly) contains it and cannot
// separate it as its own Ritz value.
//
// A single Lanczos run is not enough: if the seed coordinate is (near) an
// isolated eigenvector, the segment breaks down at step one, exploring only that
// one invariant coordinate and learning nothing about the rest of the complement
// -- treating that as a clean probe is exactly the failure this guards against.
// So segments that break down are not evidence: the explored directions are
// deflated out and the probe re-seeds from the next deterministic coordinate,
// until the honest step budget is spent or the complement is exhausted. An
// intruder in ANY segment defers. The verdict is "clean" only when either the
// complement was fully exhausted, or at least one segment ran a minimum mixed
// length (min(4, remaining complement dimension) steps) without breaking down --
// evidence that a generically mixing Krylov vector amplified dominant intruders
// and found none. If the whole budget is consumed by short breakdown segments
// and the complement is not exhausted, the probe is inconclusive and defers:
// for a genuinely diagonal-like operator with unexplored coordinates no cheap
// probe can certify the window, and deferring to the full-subspace sweep
// boundary is correct. Returns 1 (clean), 0 (defer), or -1 on operator failure.
// Deflated applies count in matvecs/operator_columns honestly. Targets a short
// Lanczos cannot resolve to the relevant extreme (smallest_magnitude, unknown)
// conservatively defer.
static int block_lanczos_window_complement_clean(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int target_kind,
    const double* V_out,
    int n_locked,
    const double* window_vecs,
    int wanted,
    double window_edge,
    double norm_a,
    double tol,
    int steps,
    int probe_seed,
    EigencoreWorkspace* workspace,
    int* matvecs_out,
    int* operator_columns_out,
    int* status_out) {
  *status_out = 0;
  if (target_kind != 1 && target_kind != 2 && target_kind != 3) {
    return 0;  // smallest_magnitude / unknown: cannot certify cheaply -> defer.
  }
  const int d0 = n_locked + wanted;
  const int comp_dim = n - d0;
  if (comp_dim <= 0) {
    return 1;  // no complement: nothing can hide.
  }
  int budget = steps;
  if (budget > comp_dim) budget = comp_dim;
  if (budget < 1) return 1;

  // Projection basis P = [locked | window | explored segment vectors]. Segment
  // vectors are appended as they are generated, so later seeds and the Lanczos
  // reorthogonalization see the whole explored subspace.
  const int cap = d0 + budget;
  std::vector<double> P(static_cast<size_t>(n) * static_cast<size_t>(cap), 0.0);
  for (int c = 0; c < n_locked; ++c) {
    std::memcpy(P.data() + static_cast<int64_t>(c) * n,
                V_out + static_cast<int64_t>(c) * n,
                sizeof(double) * static_cast<size_t>(n));
  }
  for (int c = 0; c < wanted; ++c) {
    std::memcpy(P.data() + static_cast<int64_t>(n_locked + c) * n,
                window_vecs + static_cast<int64_t>(c) * n,
                sizeof(double) * static_cast<size_t>(n));
  }
  int p_cols = d0;

  std::vector<double> q(static_cast<size_t>(n), 0.0);
  std::vector<double> Aq(static_cast<size_t>(n), 0.0);
  std::vector<double> alpha(static_cast<size_t>(budget), 0.0);
  std::vector<double> beta(static_cast<size_t>(budget), 0.0);

  double d_min = R_PosInf;
  double d_max = R_NegInf;
  bool had_mixed = false;
  bool exhausted = false;
  int total_steps = 0;
  int attempt = 0;

  while (total_steps < budget) {
    eigencore_check_interrupt();
    // Seed the next segment: the next deterministic unit vector, projected into
    // the complement of everything explored so far.
    double seed_norm = 0.0;
    bool got_seed = false;
    for (; attempt < n + budget; ++attempt) {
      std::memset(q.data(), 0, sizeof(double) * static_cast<size_t>(n));
      const int idx = ((probe_seed + 1) * 17 + attempt * 31) % n;
      q[idx < 0 ? -idx : idx] = 1.0;
      block_lanczos_project_out(q.data(), n, P.data(), p_cols, 2);
      seed_norm = trl_norm2(q.data(), n);
      if (seed_norm > 1e-8) {
        ++attempt;
        got_seed = true;
        break;
      }
    }
    if (!got_seed) {
      exhausted = true;  // no unit coordinate survives -> complement explored.
      break;
    }
    {
      const double inv = 1.0 / seed_norm;
      for (int r = 0; r < n; ++r) q[r] *= inv;
    }

    const int remaining = comp_dim - (p_cols - d0);
    int seg_max = budget - total_steps;
    if (seg_max > remaining) seg_max = remaining;
    const int min_mixed = remaining < 4 ? remaining : 4;

    int seg_len = 0;
    bool seg_broke = false;
    for (int j = 0; j < seg_max; ++j) {
      std::memcpy(P.data() + static_cast<int64_t>(p_cols + j) * n, q.data(),
                  sizeof(double) * static_cast<size_t>(n));
      const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, q.data(), n,
                           1.0, 0.0, Aq.data(), n, workspace);
      if (matvecs_out != nullptr) ++(*matvecs_out);
      if (operator_columns_out != nullptr) ++(*operator_columns_out);
      ++total_steps;
      if (rc != 0) {
        *status_out = rc;
        return -1;
      }
      if (j > 0) {
        const double b = beta[j - 1];
        const double* qprev = P.data() + static_cast<int64_t>(p_cols + j - 1) * n;
        for (int r = 0; r < n; ++r) Aq[r] -= b * qprev[r];
      }
      alpha[j] = ec_dot(q.data(), Aq.data(), n);
      for (int r = 0; r < n; ++r) Aq[r] -= alpha[j] * q[r];
      // Full reorthogonalization against everything explored (deflation basis
      // plus all segment vectors, including this segment so far).
      block_lanczos_project_out(Aq.data(), n, P.data(), p_cols + j + 1, 2);
      const double bn = trl_norm2(Aq.data(), n);
      seg_len = j + 1;
      beta[j] = bn;
      if (bn <= 1e-12 * ((norm_a > 0.0) ? norm_a : 1.0)) {
        seg_broke = true;
        break;
      }
      if (j + 1 < seg_max) {
        const double inv = 1.0 / bn;
        for (int r = 0; r < n; ++r) q[r] = Aq[r] * inv;
      }
    }
    p_cols += seg_len;

    // Extremes of this segment's tridiagonal feed the running d_min / d_max.
    std::vector<double> diag(alpha.begin(), alpha.begin() + seg_len);
    double seg_dmin = diag[0];
    double seg_dmax = diag[0];
    if (seg_len > 1) {
      std::vector<double> off(static_cast<size_t>(seg_len - 1), 0.0);
      for (int j = 0; j < seg_len - 1; ++j) off[static_cast<size_t>(j)] = beta[j];
      int info = 0;
      F77_CALL(dsterf)(&seg_len, diag.data(), off.data(), &info);
      if (info != 0) {
        return 0;  // eigensolve failed: be conservative and defer.
      }
      seg_dmin = diag[0];
      seg_dmax = diag[static_cast<size_t>(seg_len - 1)];
    }
    if (seg_dmin < d_min) d_min = seg_dmin;
    if (seg_dmax > d_max) d_max = seg_dmax;

    if (block_lanczos_complement_intruder(target_kind, d_min, d_max, window_edge,
                                          norm_a, tol)) {
      return 0;  // an intruder in any segment defers.
    }
    // A segment that runs the mixed minimum without breaking down had a
    // generically mixing seed: dominant intruders would have amplified.
    if (!seg_broke && seg_len >= min_mixed) {
      had_mixed = true;
    }
    if (p_cols >= n) {
      exhausted = true;
      break;
    }
  }

  // No intruder found. Clean only with positive evidence: the complement was
  // exhausted, or a mixed segment vouched for it. Otherwise inconclusive: defer.
  return (exhausted || had_mixed) ? 1 : 0;
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
    int check_stride,
    const double* start_block,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* n_locked_out,
    int* iterations_out,
    int* matvecs_out,
    int* operator_columns_out,
    int* certification_operator_columns_out,
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
  *operator_columns_out = 0;
  *certification_operator_columns_out = 0;
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
  // Frees buf on every exit path, including C++ exceptions (C10).
  TrlBuffersGuard buf_guard(&buf);
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
                              matvecs_out, operator_columns_out);
  stages->apply += native_timer_elapsed(timer);
  if (rc != 0) {
    return rc;
  }
  // The start block's projected column is formed by the first expansion step
  // (or explicitly when the sweep ends before one runs).
  bool last_column_known = false;
  // Operator-norm scale for relative breakdown/convergence decisions (C16):
  // the caller's estimate when available, raised by every Rayleigh quotient.
  double norm_scale = (R_FINITE(norm_a) && norm_a > 0.0) ? norm_a : 0.0;

  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int selected_count_final = 0;
  int restart_idx = 0;
  bool have_last_rr = false;

  // Probe-confirmed mid-sweep termination state. A mid-sweep window that looks
  // complete is not locked on sight: it must first pass a deflated-complement
  // check (nothing more preferred than its edge hides in its complement). The
  // check is cached by value set -- last_defer_window records the most recent
  // window that failed, so while the window is stuck the check is not repeated;
  // it reruns only when the window changes. None of this is reachable when
  // check_stride == 0 (every evaluation is a sweep boundary), preserving
  // bit-identity with the legacy path.
  std::vector<double> last_defer_window(static_cast<size_t>(k_target), 0.0);
  bool have_last_defer_window = false;
  int probe_count = 0;
  // Deflated Lanczos steps for the complement check. Enough to resolve a
  // well-separated complement extreme (a severe miss shows within a handful of
  // steps; a near-tied miss surfaces slowly but then the wrong window differs
  // from the right one only within the tie), yet small enough that a confirmed
  // window stays cheap.
  const int complement_steps = 8;

  while (true) {
    eigencore_check_interrupt();
    // Expand the active basis by one chunk. With mid-sweep checks enabled the
    // chunk is check_stride blocks wide; otherwise the whole sweep is a single
    // chunk up to the m_max budget (legacy behaviour). The chunk target is
    // computed in 64-bit and clamped to m_max so a large check_stride cannot
    // overflow int or push m_stop past the budget.
    int m_stop = m_max;
    if (check_stride > 0) {
      const int64_t chunk = static_cast<int64_t>(check_stride) *
                            static_cast<int64_t>(block_size);
      const int64_t target = static_cast<int64_t>(m_active) + chunk;
      m_stop = (target < static_cast<int64_t>(m_max)) ?
        static_cast<int>(target) : m_max;
    }
    rc = block_lanczos_expand_basis_to_budget(
      impl, apply, n, m_max, m_stop, block_size, V_out, n_locked, norm_scale,
      &buf, &workspace, stages, &m_active, &previous_block_start,
      &previous_block_cols, &last_block_start, &last_block_cols,
      &last_column_known, iterations_out,
      matvecs_out, operator_columns_out, ortho_passes_out
    );
    if (rc != 0) {
      return rc;
    }

    if (m_active < 1) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    // The sweep is complete once the budget is reached, or the expansion could
    // not fill this chunk (breakdown / exhausted continuation). A checkpoint
    // that stopped strictly at the chunk boundary (m_stop < m_max) is a
    // mid-sweep evaluation. With check_stride == 0 the sweep is always complete
    // here, so every downstream decision matches the legacy path bit-for-bit.
    const bool sweep_complete = (m_active >= m_max) || (m_active < m_stop);

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
    // Running operator-norm scale: ||A|| >= max |T_ij| for the orthonormal
    // projection T = V'AV.
    for (int col = 0; col < m_active; ++col) {
      for (int row = 0; row <= col; ++row) {
        const double t = fabs(buf.S_eig[row + static_cast<int64_t>(col) * m_active]);
        if (t > norm_scale && R_FINITE(t)) {
          norm_scale = t;
        }
      }
    }
    timer = native_timer_now();
    rc = projected_eigen_selected(buf.T_proj, m_max, buf.S_eig, m_active,
                                  selected_count,
                                  target_kind, buf.theta, buf.selected,
                                  buf.S_selected, buf.Z_eig, buf.w_eig,
                                  buf.isuppz, buf.dsyev_work, buf.dsyev_lwork,
                                  buf.dsyevd_iwork, buf.dsyevd_liwork);
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projected_eigensolve += elapsed;
    }
    if (rc != 0) {
      return rc;
    }
    have_last_rr = true;

    timer = native_timer_now();
    for (int p = 0; p < selected_count; ++p) {
      buf.ritz_res[p] = R_PosInf;
      buf.is_locked[p] = 0;
    }
    stages->selected_vector_copy += native_timer_elapsed(timer);

    timer = native_timer_now();
    combine_basis_columns(buf.V_active, n, m_active,
                          buf.S_selected, m_active,
                          selected_count, buf.B_v);
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_vector_form += elapsed;
    }

    // Mid-sweep checkpoints must not spend operator applications: recombine the
    // cached AV_active block (dgemm) for the residual regardless of the passed
    // apply_ritz_vectors policy. Sweep-boundary evaluations keep that policy
    // exactly as before. block_lanczos_finalize_return always re-certifies the
    // returned block with fresh operator applications, so the cached residual
    // here only steers mid-sweep locking, never the reported certificate.
    const bool eval_apply_ritz = (apply_ritz_vectors != 0) && sweep_complete;
    timer = native_timer_now();
    rc = 0;
    if (eval_apply_ritz) {
      rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, selected_count,
                 buf.B_v, n, 1.0, 0.0, buf.B_av, n, &workspace);
      if (rc == 0 && matvecs_out != nullptr) {
        ++(*matvecs_out);
      }
      if (rc == 0 && operator_columns_out != nullptr) {
        *operator_columns_out += selected_count;
      }
      if (rc == 0 && certification_operator_columns_out != nullptr) {
        *certification_operator_columns_out += selected_count;
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
      return rc;
    }

    timer = native_timer_now();
    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      buf.ritz_res[p] = ec_residual_norm(buf.B_av + static_cast<int64_t>(p) * n,
                                         buf.B_v + static_cast<int64_t>(p) * n,
                                         buf.theta[idx], n);
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
          norm_scale, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
        );
        const double backward_error = buf.ritz_res[p] / scale_i;
        if (ISNAN(buf.ritz_res[p]) || buf.ritz_res[p] > max_residual) {
          max_residual = buf.ritz_res[p];
        }
        if (ISNAN(backward_error) || backward_error > max_backward_error) {
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
    if (sweep_complete) {
      // Sweep-boundary evaluation: legacy incremental locking. The thick restart
      // deflates the locked directions out of the active basis before the next
      // Rayleigh-Ritz, so re-selection of a locked pair (the independence-guard
      // branch) cannot corrupt the target set. This is exactly the behaviour
      // predating mid-sweep checks, and with check_stride == 0 every evaluation
      // is a sweep boundary, so that path stays bit-for-bit identical. A sweep
      // boundary (full m_max subspace) supersedes any cached mid-sweep defer.
      have_last_defer_window = false;
      for (int p = 0; p < wanted && p < selected_count; ++p) {
        const int idx = buf.selected[p];
        const double scale_i = standard_eigen_lock_scale(
          norm_scale, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
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
    } else {
      // Mid-sweep evaluation: complement-confirmed all-or-nothing termination.
      //
      // First gate -- window readiness. The whole remaining wanted window from
      // THIS single Rayleigh-Ritz must be converged and independent of the
      // locked set (all-or-nothing; the active basis is not deflated mid-sweep,
      // so partial incremental locking could substitute a non-target pair), and
      // cluster-clear (no numerically unresolved adjacent Ritz pair at or past
      // the window edge -- a cheap first cut that defers a window holding two
      // resolved copies of a tight cluster to the full-subspace sweep boundary).
      bool window_ready = (selected_count >= wanted && wanted > 0);
      for (int p = 0; window_ready && p < wanted; ++p) {
        const int idx = buf.selected[p];
        const double scale_i = standard_eigen_lock_scale(
          norm_scale, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
        );
        if (buf.ritz_res[p] > tol * scale_i ||
            !vector_is_independent_from_locked(
                V_out, n_locked, buf.B_v + static_cast<int64_t>(p) * n, n)) {
          window_ready = false;
        }
      }
      const double set_tol = fmax(1e-6, 10.0 * tol);
      if (window_ready) {
        const int probe = (wanted + 1 <= selected_count) ? wanted + 1 : selected_count;
        for (int p = 0; p + 1 < probe; ++p) {
          const double hi = buf.theta[buf.selected[p]];
          const double lo = buf.theta[buf.selected[p + 1]];
          const double scale = fmax(fmax(fabs(hi), fabs(lo)), norm_scale);
          if (fabs(hi - lo) <= set_tol * scale) {
            window_ready = false;
            break;
          }
        }
      }

      // Second gate -- deflated-complement confirmation. A ready window is NOT
      // locked on sight: a small subspace can hold a target direction only
      // weakly, unresolved as its own Ritz value, so the window silently omits
      // it with no locally visible signal (the window is a genuine set of
      // converged, independent, cluster-clear eigenpairs). Deflate the window
      // and probe its complement; if nothing more preferred than the window edge
      // hides there, lock. Otherwise defer -- caching the window so the check is
      // not repeated while it is stuck -- and keep iterating; natural expansion
      // eventually resolves the omitted direction (which then either enters the
      // window, tripping the cluster-clearance cut, or is caught again here) and
      // the full-subspace sweep boundary reproduces the legacy set.
      if (window_ready) {
        bool same_as_last = have_last_defer_window;
        for (int p = 0; same_as_last && p < wanted; ++p) {
          const double cur = buf.theta[buf.selected[p]];
          const double ref = last_defer_window[static_cast<size_t>(p)];
          const double scale = fmax(fmax(fabs(cur), fabs(ref)), norm_scale);
          if (fabs(cur - ref) > set_tol * scale) {
            same_as_last = false;
          }
        }
        if (!same_as_last) {
          const double window_edge = buf.theta[buf.selected[wanted - 1]];
          stages->locking += native_timer_elapsed(timer);
          int status = 0;
          const int clean = block_lanczos_window_complement_clean(
            impl, apply, n, target_kind, V_out, n_locked, buf.B_v, wanted,
            window_edge, norm_scale, tol, complement_steps, probe_count, &workspace,
            matvecs_out, operator_columns_out, &status);
          ++probe_count;
          timer = native_timer_now();
          if (clean < 0) {
            return status != 0 ? status : -1;
          }
          if (clean == 1) {
            for (int p = 0; p < wanted; ++p) {
              const int idx = buf.selected[p];
              std::memcpy(V_out + static_cast<int64_t>(n_locked) * n,
                          buf.B_v + static_cast<int64_t>(p) * n,
                          sizeof(double) * static_cast<size_t>(n));
              lambda_out[n_locked] = buf.theta[idx];
              residuals_out[n_locked] = buf.ritz_res[p];
              converged_out[n_locked] = 1;
              buf.is_locked[p] = 1;
              ++n_locked;
              ++lock_now;
            }
          } else {
            for (int p = 0; p < wanted; ++p) {
              last_defer_window[static_cast<size_t>(p)] = buf.theta[buf.selected[p]];
            }
            have_last_defer_window = true;
          }
        }
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
      n, k_target, selected_count, n_locked, norm_scale, tol, &buf,
      V_out, lambda_out, residuals_out, &best
    );

    if (n_locked >= k_target) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }
    // Evaluate the final sweep before stopping, exactly as the legacy loop did
    // at restart_idx == max_restarts. A mid-sweep checkpoint never terminates
    // here: it must be able to keep expanding within the current sweep.
    if (sweep_complete && restart_idx == max_restarts) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    if (sweep_complete) {
      // A thick restart (and only a thick restart) advances restart_idx, so the
      // restart_idx / max_restarts accounting matches the legacy semantics.
      rc = block_lanczos_restart_with_continuation_tail(
        impl, apply, n, k_target, m_max, block_size, restart_idx, selected_count,
        n_locked, norm_scale, &buf, &workspace, stages, &m_active,
        &previous_block_start, &previous_block_cols, &last_block_start,
        &last_block_cols, &last_column_known, V_out,
        matvecs_out, operator_columns_out, restarts_out, ortho_passes_out
      );
      if (rc == 1) {
        *restarts_out = restart_idx;
        *m_active_final_out = m_active;
        break;
      }
      if (rc != 0) {
        return rc;
      }
      ++restart_idx;
    }
    // Otherwise this was a mid-sweep checkpoint: loop to expand the next chunk
    // of the current sweep without restarting.
  }

  rc = block_lanczos_finalize_return(
    impl, apply, n, k_target, target_kind, tol, norm_scale, selected_count_final,
    have_last_rr, &buf, &workspace, stages, best, V_out, lambda_out,
    residuals_out, converged_out, &n_locked, matvecs_out,
    operator_columns_out, certification_operator_columns_out
  );
  if (rc != 0) {
    return rc;
  }
  *n_locked_out = n_locked;
  *m_active_final_out = m_active;
  *operator_allocations_out = workspace.allocation_count;
  *operator_bytes_allocated_out = workspace.bytes_allocated;
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
                                            int operator_columns,
                                            int certification_operator_columns,
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

  SEXP out_ = PROTECT(allocVector(VECSXP, 19));
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
  SET_VECTOR_ELT(out_, 16, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 17, ScalarInteger(operator_columns));
  SET_VECTOR_ELT(out_, 18, ScalarInteger(certification_operator_columns));
  SEXP names_ = PROTECT(allocVector(STRSXP, 19));
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
  SET_STRING_ELT(names_, 16, mkChar("operator_block_calls"));
  SET_STRING_ELT(names_, 17, mkChar("operator_columns"));
  SET_STRING_ELT(names_, 18, mkChar("certification_operator_columns"));
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_block_lanczos_csc(SEXP i_, SEXP p_, SEXP x_,
                                            SEXP dim_, SEXP k_,
                                            SEXP m_max_,
                                            SEXP block_size_,
                                            SEXP target_kind_,
                                            SEXP tol_,
                                            SEXP start_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC block Lanczos inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "block Lanczos");
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
  EIGENCORE_ENTRY_END
}

// Shared driver for the three concrete block thick-restart Lanczos bindings
// (dense, CSC, matrix-free R operator). The per-storage wrappers validate and
// build their operator impl, then hand off here so allocation, restart-history
// sizing, the run call, and result packing stay identical across bindings.
// apply_ritz_vectors selects the sweep-boundary residual policy (0: recombine
// cached AV_active, 1: re-apply the operator). check_stride enables mid-sweep
// convergence checks (0: legacy single-sweep expansion).
static SEXP block_thick_restart_lanczos_impl(
    void* impl, EigencoreApplyFn apply, int n, int apply_ritz_vectors,
    int k, int m_max, int block_size, int target_kind, double tol,
    int max_restarts, double norm_a, const double* start, int check_stride) {
  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, operator_columns = 0;
  int certification_operator_columns = 0, restarts = 0, m_active = 0;
  int locking_events = 0, ortho_passes = 0;
  int64_t operator_allocations = 0, operator_bytes_allocated = 0;
  NativeBlockStageSeconds stage_seconds;
  // Each mid-sweep checkpoint records a history row, so budget one row per chunk
  // per restart cycle (capped) when check_stride > 0; legacy sizing otherwise.
  int history_capacity = max_restarts + 1;
  if (check_stride > 0) {
    const int chunk = (check_stride * block_size > 1) ? check_stride * block_size : 1;
    const long evals_per_sweep = (static_cast<long>(m_max) + chunk - 1) / chunk;
    long cap = static_cast<long>(max_restarts + 1) *
               (evals_per_sweep > 1 ? evals_per_sweep : 1);
    if (cap > 8192) cap = 8192;
    if (cap < 1) cap = 1;
    history_capacity = static_cast<int>(cap);
  }
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

  const int status = native_block_thick_restart_lanczos_run(
    impl, apply, n, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, apply_ritz_vectors, check_stride, start,
    V.data(), lambda.data(), residuals.data(), converged.data(), &n_locked,
    &iterations, &matvecs, &operator_columns, &certification_operator_columns,
    &restarts, &m_active, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes_allocated, &stage_seconds, &history);
  if (status != 0) {
    error("native block thick-restart Lanczos failed with status=%d", status);
  }
  return block_thick_lanczos_pack_result(
    n, k, V.data(), lambda.data(), residuals.data(), converged.data(),
    n_locked, iterations, matvecs, operator_columns,
    certification_operator_columns, restarts, m_active, locking_events,
    ortho_passes, block_size, operator_allocations, operator_bytes_allocated,
    &stage_seconds, &history);
}

extern "C" SEXP eigencore_block_thick_restart_lanczos_dense(
    SEXP A_, SEXP k_, SEXP m_max_, SEXP block_size_,
    SEXP target_kind_, SEXP tol_, SEXP max_restarts_,
    SEXP norm_a_, SEXP start_, SEXP check_stride_) {
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
  const int check_stride = static_cast<int>(asInteger(check_stride_));
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");
  if (check_stride < 0) error("check_stride must be >= 0");

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  return block_thick_restart_lanczos_impl(
    &impl, eigencore_dense_apply, n, 0, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, REAL(start_), check_stride);
  EIGENCORE_ENTRY_END
}
extern "C" SEXP eigencore_block_thick_restart_lanczos_csc(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP k_,
    SEXP m_max_, SEXP block_size_, SEXP target_kind_,
    SEXP tol_, SEXP max_restarts_, SEXP norm_a_, SEXP start_,
    SEXP check_stride_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC block thick-restart Lanczos inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "block thick-restart Lanczos");
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
  const int check_stride = static_cast<int>(asInteger(check_stride_));
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");
  if (check_stride < 0) error("check_stride must be >= 0");

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return block_thick_restart_lanczos_impl(
    &impl, eigencore_csc_apply, n, 1, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, REAL(start_), check_stride);
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_block_thick_restart_lanczos_r_operator(
    SEXP dim_, SEXP apply_, SEXP k_, SEXP m_max_, SEXP block_size_,
    SEXP target_kind_, SEXP tol_, SEXP max_restarts_, SEXP norm_a_, SEXP start_,
    SEXP check_stride_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(dim_) || LENGTH(dim_) != 2 || TYPEOF(apply_) != CLOSXP ||
      !isReal(start_)) {
    error("invalid matrix-free block thick-restart Lanczos inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n) {
    error("A must be a square matrix-free operator");
  }
  if (INTEGER(dimS)[0] != n) {
    error("non-conformable matrix-free block thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  const double norm_a = asReal(norm_a_);
  const int check_stride = static_cast<int>(asInteger(check_stride_));
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= the operator dimension");
  if (max_restarts < 0) error("max_restarts must be >= 0");
  if (check_stride < 0) error("check_stride must be >= 0");

  // The Hermitian kernel only applies A (no adjoint); apply_ritz_vectors = 1
  // re-applies the callback for sweep-boundary residuals, matching the CSC path.
  RApplyOperator impl = {n, n, apply_, R_NilValue};
  return block_thick_restart_lanczos_impl(
    &impl, eigencore_r_operator_apply, n, 1, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, REAL(start_), check_stride);
  EIGENCORE_ENTRY_END
}

// Shared driver for the implicit normal-equations (Gram) thick-restart
// Lanczos: runs the production block thick-restart Lanczos on A^T A
// (side == 0, subspace in R^cols) or A A^T (side == 1, subspace in R^rows)
// without materializing the Gram matrix. Each operator application costs one
// forward and one adjoint apply of A; matvecs are reported in base-A applies
// (2 per normal-operator application).
static SEXP normal_thick_restart_lanczos_impl(
    void* base_impl, EigencoreApplyFn base_apply,
    int64_t rows, int64_t cols, int side,
    SEXP k_, SEXP m_max_, SEXP block_size_, SEXP target_kind_,
    SEXP tol_, SEXP max_restarts_, SEXP norm_a_, SEXP start_) {
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (!isReal(start_) || dimS == R_NilValue) {
    error("start must be a double matrix");
  }
  if (side != 0 && side != 1) {
    error("side must be 0 (A^T A) or 1 (A A^T)");
  }
  const int64_t outer64 = (side == 0) ? cols : rows;
  const int64_t inner64 = (side == 0) ? rows : cols;
  if (!eigencore_int_indexable(outer64) || !eigencore_int_indexable(inner64)) {
    error("normal-equations Lanczos dimensions exceed native integer range");
  }
  const int n = static_cast<int>(outer64);
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  const double norm_a = asReal(norm_a_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[0] != n) error("start block has wrong number of rows");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= the normal-operator dimension");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, operator_columns = 0;
  int certification_operator_columns = 0, restarts = 0, m_active = 0;
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

  // The driver applies the operator to blocks as wide as the full active
  // subspace, so the intermediate A-product needs inner x m_max capacity.
  std::vector<double> scratch(static_cast<size_t>(inner64) *
                              static_cast<size_t>(m_max), 0.0);
  NormalEquationsOperator impl = {base_impl, base_apply, rows, cols, side,
                                  scratch.data(), m_max};
  // apply_ritz_vectors = 0: locking residuals come from the cached AV_active
  // combination (dgemm) instead of re-applying the normal operator, which
  // would cost 2*k base applies per restart. The R caller certifies the final
  // triplets with exact residuals in original coordinates regardless.
  // check_stride = 0: mid-sweep checks stay off on the normal-equations path
  // (legacy single-sweep expansion); enabling them there would require threading
  // the parameter through the separate normal-equations SEXP/R plumbing.
  const int status = native_block_thick_restart_lanczos_run(
    &impl, eigencore_normal_equations_apply, n, k, m_max, block_size,
    target_kind, tol, max_restarts, norm_a, 0, 0, REAL(start_), V.data(),
    lambda.data(), residuals.data(), converged.data(), &n_locked, &iterations,
    &matvecs, &operator_columns, &certification_operator_columns,
    &restarts, &m_active, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes_allocated, &stage_seconds, &history);
  if (status != 0) {
    error("native normal-equations thick-restart Lanczos failed with status=%d",
          status);
  }
  matvecs *= 2;  // each normal-operator application is two base-A applies
  operator_columns *= 2;
  certification_operator_columns *= 2;
  return block_thick_lanczos_pack_result(
    n, k, V.data(), lambda.data(), residuals.data(), converged.data(),
    n_locked, iterations, matvecs, operator_columns,
    certification_operator_columns, restarts, m_active, locking_events,
    ortho_passes, block_size, operator_allocations, operator_bytes_allocated,
    &stage_seconds, &history);
}

extern "C" SEXP eigencore_normal_thick_restart_lanczos_dense(
    SEXP A_, SEXP side_, SEXP k_, SEXP m_max_, SEXP block_size_,
    SEXP target_kind_, SEXP tol_, SEXP max_restarts_,
    SEXP norm_a_, SEXP start_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  DenseColumnMajorOperator base = {m, n, REAL(A_)};
  return normal_thick_restart_lanczos_impl(
    &base, eigencore_dense_apply, m, n,
    static_cast<int>(asInteger(side_)),
    k_, m_max_, block_size_, target_kind_, tol_, max_restarts_,
    norm_a_, start_);
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_normal_thick_restart_lanczos_csc(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP side_, SEXP k_,
    SEXP m_max_, SEXP block_size_, SEXP target_kind_,
    SEXP tol_, SEXP max_restarts_, SEXP norm_a_, SEXP start_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC normal-equations Lanczos inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "normal-equations Lanczos");
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  CSCOperator base = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return normal_thick_restart_lanczos_impl(
    &base, eigencore_csc_apply, m, n,
    static_cast<int>(asInteger(side_)),
    k_, m_max_, block_size_, target_kind_, tol_, max_restarts_,
    norm_a_, start_);
  EIGENCORE_ENTRY_END
}
