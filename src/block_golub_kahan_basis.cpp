#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "native_operators.h"
#include "block_golub_kahan_basis.h"

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

static void symmetrize_packed_square(double* A, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double avg = 0.5 * (A[i + j * n] + A[j + i * n]);
      A[i + j * n] = avg;
      A[j + i * n] = avg;
    }
  }
}

static int block_accept_work_vector(const double* V_locked, int n_locked,
                                    double* V_active, int* m_active,
                                    int m_max, double* z, double* tmp,
                                    int n, int* ortho_passes) {
  if (*m_active >= m_max) {
    return 0;
  }
  const int passes = 2;
  trl_orthogonalise(V_locked, n_locked, V_active, *m_active, z, tmp, n, passes);
  if (ortho_passes != nullptr) {
    *ortho_passes += passes;
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

static void block_reorthogonalise_against(const double* V_locked, int n_locked,
                                          const double* V_active, int m_active,
                                          double* X, int n, int cols,
                                          double* coeff, int passes) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  for (int pass = 0; pass < passes; ++pass) {
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
  }
}

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

  const int reorth_passes = 2;
  const double* active_basis = reorthogonalize_active ? V_active : nullptr;
  const int active_cols = reorthogonalize_active ? *m_active : 0;
  block_reorthogonalise_against(V_locked, n_locked, active_basis, active_cols,
                                Z_block, n, cols, coeff, reorth_passes);
  if (ortho_passes != nullptr) {
    *ortho_passes += reorth_passes;
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

static double* eigencore_r_alloc_zero_doubles(size_t n) {
  double* out = reinterpret_cast<double*>(R_alloc(n > 0 ? n : 1, sizeof(double)));
  std::memset(out, 0, sizeof(double) * n);
  return out;
}

void block_golub_kahan_basis_scratch_free(BlockGolubKahanBasisScratch* scratch) {
  if (!scratch->transient) {
    std::free(scratch->Z_v);
    std::free(scratch->Z_u);
    std::free(scratch->coeff);
    std::free(scratch->tmp);
  }
  scratch->Z_v = nullptr;
  scratch->Z_u = nullptr;
  scratch->coeff = nullptr;
  scratch->tmp = nullptr;
  scratch->bytes = 0;
  scratch->transient = false;
}

int block_golub_kahan_basis_scratch_alloc(BlockGolubKahanBasisScratch* scratch,
                                                 int m,
                                                 int n,
                                                 int max_subspace,
                                                 int block_size) {
  const size_t nb = static_cast<size_t>(n) * static_cast<size_t>(block_size);
  const size_t mb = static_cast<size_t>(m) * static_cast<size_t>(block_size);
  const int coeff_rows = (max_subspace > block_size) ? max_subspace : block_size;
  const size_t coeff_elems = static_cast<size_t>(coeff_rows) *
    static_cast<size_t>(block_size);
  const size_t tmp_len = static_cast<size_t>((m > n) ? m : n);
  scratch->Z_v = eigencore_r_alloc_zero_doubles(nb);
  scratch->Z_u = eigencore_r_alloc_zero_doubles(mb);
  scratch->coeff = eigencore_r_alloc_zero_doubles(coeff_elems);
  scratch->tmp = eigencore_r_alloc_zero_doubles(tmp_len);
  scratch->transient = true;
  scratch->bytes = (nb + mb + coeff_elems + tmp_len) * sizeof(double);
  return 0;
}

int native_block_golub_kahan_basis_run_with_scratch(
                                              void* impl,
                                              EigencoreApplyFn apply,
                                              int m,
                                              int n,
                                              int max_subspace,
                                              int block_size,
                                              const double* start_block,
                                              const double* start_av_block,
                                              int start_av_cols,
                                              const double* V_locked,
                                              int n_locked_v,
                                              const double* U_locked,
                                              int n_locked_u,
                                              double* V,
                                              double* AV,
                                              double* U,
                                              BlockGolubKahanBasisScratch* scratch,
                                              int* active_v_out,
                                              int* active_u_out,
                                              int* iterations_out,
                                              int* matvecs_out,
                                              int* ortho_passes_out,
                                              int* cached_start_used_out) {
  *active_v_out = 0;
  *active_u_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *ortho_passes_out = 0;
  if (cached_start_used_out != nullptr) {
    *cached_start_used_out = 0;
  }

  const size_t nv = static_cast<size_t>(n) * max_subspace;
  const size_t mv = static_cast<size_t>(m) * max_subspace;
  std::memset(V, 0, sizeof(double) * nv);
  std::memset(AV, 0, sizeof(double) * mv);
  std::memset(U, 0, sizeof(double) * mv);

  if (scratch == nullptr || scratch->Z_v == nullptr || scratch->Z_u == nullptr ||
      scratch->coeff == nullptr || scratch->tmp == nullptr) {
    return -2;
  }
  double* Z_v = scratch->Z_v;
  double* Z_u = scratch->Z_u;
  double* coeff = scratch->coeff;
  double* tmp = scratch->tmp;

  int active_v = 0;
  int active_u = 0;
  int last_v_start = 0;
  int last_v_cols = 0;
  if (start_av_block != nullptr && start_av_cols > 0) {
    int cached_cols = start_av_cols;
    if (cached_cols > block_size) {
      cached_cols = block_size;
    }
    if (cached_cols > max_subspace) {
      cached_cols = max_subspace;
    }
    for (int col = 0; col < cached_cols; ++col) {
      std::memcpy(V + static_cast<int64_t>(col) * n,
                  start_block + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
      std::memcpy(AV + static_cast<int64_t>(col) * m,
                  start_av_block + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
    }
    active_v = cached_cols;
    if (cached_start_used_out != nullptr) {
      *cached_start_used_out = 1;
    }
    if (cached_cols < block_size && active_v < max_subspace) {
      const int suffix_start = active_v;
      const int suffix_cols = block_accept_columns_blas3(
        start_block + static_cast<int64_t>(cached_cols) * n, n,
        block_size - cached_cols, V_locked, n_locked_v,
        V, &active_v, max_subspace, Z_v, block_size,
        coeff, tmp, n, block_size - cached_cols, ortho_passes_out
      );
      if (suffix_cols > 0) {
        EigencoreWorkspace suffix_workspace = {0, 0, nullptr, 0};
        const int rc_suffix = apply(
          impl, EIGENCORE_TRANSPOSE_NONE, suffix_cols,
          V + static_cast<int64_t>(suffix_start) * n, n,
          1.0, 0.0,
          AV + static_cast<int64_t>(suffix_start) * m, m,
          &suffix_workspace
        );
        if (rc_suffix != 0) {
          return rc_suffix;
        }
        ++(*matvecs_out);
      }
    }
    last_v_cols = active_v;
  } else {
    last_v_cols = block_accept_columns_blas3(
      start_block, n, block_size, V_locked, n_locked_v,
      V, &active_v, max_subspace, Z_v, block_size,
      coeff, tmp, n, block_size, ortho_passes_out
    );
  }
  if (last_v_cols == 0) {
    for (int col = 0; col < block_size && active_v < max_subspace && col < n; ++col) {
      std::memset(Z_v, 0, sizeof(double) * static_cast<size_t>(n));
      Z_v[col] = 1.0;
      last_v_cols += block_accept_work_vector(
        V_locked, n_locked_v, V, &active_v, max_subspace,
        Z_v, tmp, n, ortho_passes_out
      );
    }
  }
  if (last_v_cols == 0) {
    return -4;
  }

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int rc = 0;
  if (cached_start_used_out == nullptr || *cached_start_used_out == 0) {
    rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, last_v_cols,
               V + static_cast<int64_t>(last_v_start) * n, n,
               1.0, 0.0,
               AV + static_cast<int64_t>(last_v_start) * m, m,
               &workspace);
    if (rc != 0) {
      return rc;
    }
    ++(*matvecs_out);
  }

  while (active_v < max_subspace && last_v_cols > 0) {
    const int accepted_u_start = active_u;
    const int accepted_u = block_accept_columns_blas3(
      AV + static_cast<int64_t>(last_v_start) * m, m, last_v_cols,
      U_locked, n_locked_u, U, &active_u, max_subspace, Z_u, block_size,
      coeff, tmp, m, block_size, ortho_passes_out
    );
    if (accepted_u == 0) {
      break;
    }

    rc = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, accepted_u,
               U + static_cast<int64_t>(accepted_u_start) * m, m,
               1.0, 0.0, Z_v, n, &workspace);
    if (rc != 0) {
      return rc;
    }
    ++(*matvecs_out);

    const int accepted_v_start = active_v;
    const int accepted_v = block_accept_columns_blas3(
      Z_v, n, accepted_u, V_locked, n_locked_v,
      V, &active_v, max_subspace, Z_v, block_size,
      coeff, tmp, n, max_subspace - active_v, ortho_passes_out
    );
    if (accepted_v == 0) {
      break;
    }

    rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, accepted_v,
               V + static_cast<int64_t>(accepted_v_start) * n, n,
               1.0, 0.0,
               AV + static_cast<int64_t>(accepted_v_start) * m, m,
               &workspace);
    if (rc != 0) {
      return rc;
    }
    ++(*matvecs_out);
    ++(*iterations_out);
    last_v_start = accepted_v_start;
    last_v_cols = accepted_v;
  }

  *active_v_out = active_v;
  *active_u_out = active_u;
  return 0;
}

int native_block_golub_kahan_basis_run(void* impl,
                                              EigencoreApplyFn apply,
                                              int m,
                                              int n,
                                              int max_subspace,
                                              int block_size,
                                              const double* start_block,
                                              const double* start_av_block,
                                              int start_av_cols,
                                              const double* V_locked,
                                              int n_locked_v,
                                              const double* U_locked,
                                              int n_locked_u,
                                              double* V,
                                              double* AV,
                                              double* U,
                                              int* active_v_out,
                                              int* active_u_out,
                                              int* iterations_out,
                                              int* matvecs_out,
                                              int* ortho_passes_out,
                                              int* cached_start_used_out) {
  BlockGolubKahanBasisScratch scratch;
  const int scratch_status = block_golub_kahan_basis_scratch_alloc(
    &scratch, m, n, max_subspace, block_size
  );
  if (scratch_status != 0) {
    return scratch_status;
  }
  const int status = native_block_golub_kahan_basis_run_with_scratch(
    impl, apply, m, n, max_subspace, block_size, start_block,
    start_av_block, start_av_cols, V_locked, n_locked_v, U_locked, n_locked_u,
    V, AV, U, &scratch,
    active_v_out, active_u_out, iterations_out, matvecs_out,
    ortho_passes_out, cached_start_used_out
  );
  block_golub_kahan_basis_scratch_free(&scratch);
  return status;
}

static SEXP block_golub_kahan_basis_pack(int n,
                                         int m,
                                         int max_subspace,
                                         const double* V,
                                         const double* AV,
                                         int active_v,
                                         int active_u,
                                         int iterations,
                                         int matvecs,
                                         int ortho_passes,
                                         int cached_start_used) {
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, max_subspace));
  SEXP AV_ = PROTECT(allocMatrix(REALSXP, m, max_subspace));
  std::memcpy(REAL(V_), V, sizeof(double) * static_cast<size_t>(n) * max_subspace);
  std::memcpy(REAL(AV_), AV, sizeof(double) * static_cast<size_t>(m) * max_subspace);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, V_);
  SET_VECTOR_ELT(out_, 1, AV_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(active_v));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(active_u));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(ortho_passes));
  SET_VECTOR_ELT(out_, 7, ScalarLogical(cached_start_used != 0));
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("V"));
  SET_STRING_ELT(names_, 1, mkChar("AV"));
  SET_STRING_ELT(names_, 2, mkChar("active_cols"));
  SET_STRING_ELT(names_, 3, mkChar("active_left_cols"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 7, mkChar("cached_start_used"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_basis(SEXP A_,
                                                        SEXP max_subspace_,
                                                        SEXP start_) {
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

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native dense block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

extern "C" SEXP eigencore_block_golub_kahan_dense_basis_cached(SEXP A_,
                                                               SEXP max_subspace_,
                                                               SEXP start_,
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
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  const int start_av_cols = INTEGER(dimAV)[1];
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native dense cached block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

extern "C" SEXP eigencore_block_golub_kahan_csc_basis(SEXP i_, SEXP p_,
                                                      SEXP x_, SEXP dim_,
                                                      SEXP max_subspace_,
                                                      SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_)) {
    error("invalid CSC block Golub-Kahan basis inputs");
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

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native CSC block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

extern "C" SEXP eigencore_block_golub_kahan_csc_basis_cached(SEXP i_, SEXP p_,
                                                             SEXP x_, SEXP dim_,
                                                             SEXP max_subspace_,
                                                             SEXP start_,
                                                             SEXP start_av_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_) || !isReal(start_av_)) {
    error("invalid cached CSC block Golub-Kahan basis inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimS == R_NilValue || dimAV == R_NilValue || LENGTH(dim_) != 2) {
    error("start and start_av must be matrices and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  const int start_av_cols = INTEGER(dimAV)[1];
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native CSC cached block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}
