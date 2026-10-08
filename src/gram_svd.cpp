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
#include "eigencore_lapack_compat.h"
#include "eigencore_common.h"
#include "native_operators.h"

extern "C" void eigencore_validate_csc_structure(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_, const char* context);

// Orthogonality loss max |X^T X - I| of a rows x k block, from the upper
// triangle of a dsyrk Gram (half the flops of a dgemm cross-product).
static double block_orthogonality_loss_dsyrk(const double* X, int rows, int k) {
  if (k <= 0) {
    return 0.0;
  }
  std::vector<double> gram(static_cast<size_t>(k) * static_cast<size_t>(k), 0.0);
  const char uplo = 'U';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dsyrk)(&uplo, &trans, &k, &rows, &one, X, &rows,
                  &zero, gram.data(), &k FCONE FCONE);
  double loss = 0.0;
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row <= col; ++row) {
      const double target = (row == col) ? 1.0 : 0.0;
      const double err = fabs(gram[row + static_cast<int64_t>(col) * k] - target);
      if (ISNAN(err)) {
        return err;
      }
      if (err > loss) {
        loss = err;
      }
    }
  }
  return loss;
}

// Exact two-sided SVD certificate in ORIGINAL coordinates for a CSC matrix A
// and candidate triplets (d, U, V):
//   left_j  = ||A v_j   - d_j u_j||,   right_j = ||A^T u_j - d_j v_j||,
// with orthogonality measured on the returned U and V themselves. Neither
// residual is taken from the Gram matrix or assumed zero because one factor
// was formed from the other. `av_known` / `atu_known` may pass the exact
// products A V / A^T U that this same call computed from A and the FINAL V / U
// while forming the opposite factor (bitwise what a recompute would give);
// any side passed as nullptr is applied here. Never pass a cached, Gram-
// derived, or otherwise stale product.
static void csc_two_sided_svd_certificate(const int* Ai, const int* Ap,
                                          const double* Ax, int m, int n, int k,
                                          const double* d, const double* U,
                                          const double* V,
                                          const double* av_known,
                                          const double* atu_known,
                                          double scale_value,
                                          double tol, double* left,
                                          double* right, double* combined,
                                          double* backward, int* converged,
                                          double* orth_u, double* orth_v) {
  if (k <= 0) {
    *orth_u = 0.0;
    *orth_v = 0.0;
    return;
  }
  CSCOperator impl = {m, n, Ai, Ap, Ax};
  std::vector<double> av_work;
  std::vector<double> atu_work;
  const double* av = av_known;
  const double* atu = atu_known;
  if (av == nullptr) {
    av_work.assign(static_cast<size_t>(m) * static_cast<size_t>(k), 0.0);
    eigencore_csc_apply(&impl, EIGENCORE_TRANSPOSE_NONE, k, V, n, 1.0, 0.0,
                        av_work.data(), m, nullptr);
    av = av_work.data();
  }
  if (atu == nullptr) {
    atu_work.assign(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
    eigencore_csc_apply(&impl, EIGENCORE_TRANSPOSE_ADJOINT, k, U, m, 1.0, 0.0,
                        atu_work.data(), n, nullptr);
    atu = atu_work.data();
  }
  for (int col = 0; col < k; ++col) {
    const double sigma = d[col];
    const double* av_col = av + static_cast<int64_t>(col) * m;
    const double* u_col = U + static_cast<int64_t>(col) * m;
    const double* atu_col = atu + static_cast<int64_t>(col) * n;
    const double* v_col = V + static_cast<int64_t>(col) * n;
    long double left_sum = 0.0L;
    long double right_sum = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double residual = av_col[row] - sigma * u_col[row];
      left_sum += static_cast<long double>(residual) * residual;
    }
    for (int row = 0; row < n; ++row) {
      const double residual = atu_col[row] - sigma * v_col[row];
      right_sum += static_cast<long double>(residual) * residual;
    }
    const double l = sqrt(static_cast<double>(left_sum));
    const double r = sqrt(static_cast<double>(right_sum));
    const double c = sqrt(l * l + r * r);
    const double be = c / scale_value;
    left[col] = l;
    right[col] = r;
    combined[col] = c;
    backward[col] = be;
    converged[col] = (R_FINITE(be) && be <= tol) ? TRUE : FALSE;
  }
  *orth_u = block_orthogonality_loss_dsyrk(U, m, k);
  *orth_v = block_orthogonality_loss_dsyrk(V, n, k);
}

// Upper triangle of G = sum_j a_j a_j^T for `count` sparse vectors of length
// `dim` stored in compressed form (ptr/idx/val), then mirrored to the lower
// triangle. Short vectors use the scalar pair loop with writes running down a
// single Gram column (idx sorted ascending, as in a valid dgCMatrix; unsorted
// input stays correct via the min/max swap). Vectors with many nonzeros are
// scattered into a dense panel and folded in with dsyrk, which beats the
// scalar loop once a vector is a sizeable fraction of dense.
static void accumulate_sparse_gram_upper(int dim, int count, const int* ptr,
                                         const int* idx, const double* val,
                                         double* G) {
  const int64_t ld = dim;
  const int panel_capacity = 64;
  const int dense_threshold = std::max(32, dim / 8);
  std::vector<double> panel;
  int panel_cols = 0;
  const char uplo = 'U';
  const char notrans = 'N';
  const double one = 1.0;
  auto flush_panel = [&]() {
    if (panel_cols == 0) {
      return;
    }
    F77_CALL(dsyrk)(&uplo, &notrans, &dim, &panel_cols, &one, panel.data(), &dim,
                    &one, G, &dim FCONE FCONE);
    std::fill(panel.begin(),
              panel.begin() + static_cast<int64_t>(panel_cols) * ld, 0.0);
    panel_cols = 0;
  };
  for (int j = 0; j < count; ++j) {
    if ((j & 63) == 0) eigencore_check_interrupt();
    const int start = ptr[j];
    const int end = ptr[j + 1];
    const int nnz = end - start;
    if (nnz <= 0) {
      continue;
    }
    if (nnz >= dense_threshold) {
      if (panel.empty()) {
        panel.assign(static_cast<size_t>(ld) * panel_capacity, 0.0);
      }
      double* dst = panel.data() + static_cast<int64_t>(panel_cols) * ld;
      for (int pos = start; pos < end; ++pos) {
        dst[idx[pos]] += val[pos];
      }
      if (++panel_cols == panel_capacity) {
        flush_panel();
      }
      continue;
    }
    for (int bb = start; bb < end; ++bb) {
      const int rb = idx[bb];
      const double xb = val[bb];
      double* g_col = G + static_cast<int64_t>(rb) * ld;
      for (int aa = start; aa <= bb; ++aa) {
        const int ra = idx[aa];
        const double update = val[aa] * xb;
        if (ra <= rb) {
          g_col[ra] += update;
        } else {
          G[rb + static_cast<int64_t>(ra) * ld] += update;
        }
      }
    }
  }
  flush_panel();
  for (int col = 0; col < dim; ++col) {
    for (int row = col + 1; row < dim; ++row) {
      G[row + static_cast<int64_t>(col) * ld] = G[col + static_cast<int64_t>(row) * ld];
    }
  }
}

static int trl_orthogonalise_gram(const double* V_locked, int n_locked,
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

static double trl_norm2_gram(const double* x, int n) {
  long double sum = 0.0L;
  for (int i = 0; i < n; ++i) {
    sum += static_cast<long double>(x[i]) * x[i];
  }
  return sqrt(static_cast<double>(sum));
}

static int trl_dsyev_query_gram(int m_max) {
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

static void symmetrize_packed_square_gram(double* A, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double avg = 0.5 * (A[i + j * n] + A[j + i * n]);
      A[i + j * n] = avg;
      A[j + i * n] = avg;
    }
  }
}

static int orthonormalize_columns_inplace(double* Q, int n, int cols,
                                          double* tmp, double tol) {
  int rank = 0;
  for (int col = 0; col < cols; ++col) {
    double* q_col = Q + static_cast<int64_t>(rank) * n;
    if (rank != col) {
      std::memcpy(q_col, Q + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
    }
    trl_orthogonalise_gram(nullptr, 0, Q, rank, q_col, tmp, n);
    const double q_norm = trl_norm2_gram(q_col, n);
    if (q_norm <= tol) {
      continue;
    }
    const double inv_norm = 1.0 / q_norm;
    for (int row = 0; row < n; ++row) {
      q_col[row] *= inv_norm;
    }
    ++rank;
  }
  return rank;
}

static void csc_transpose_apply_vec(const int* Ai, const int* Ap,
                                    const double* Ax, int m, int n,
                                    const double* x, double* y) {
  (void) m;
  std::fill(y, y + n, 0.0);
  for (int col = 0; col < n; ++col) {
    long double acc = 0.0L;
    for (int jj = Ap[col]; jj < Ap[col + 1]; ++jj) {
      acc += static_cast<long double>(Ax[jj]) * x[Ai[jj]];
    }
    y[col] = static_cast<double>(acc);
  }
}

static void csc_forward_apply_vec(const int* Ai, const int* Ap,
                                  const double* Ax, int m, int n,
                                  const double* x, double* y) {
  std::fill(y, y + m, 0.0);
  for (int col = 0; col < n; ++col) {
    const double x_col = x[col];
    if (x_col == 0.0) {
      continue;
    }
    for (int jj = Ap[col]; jj < Ap[col + 1]; ++jj) {
      y[Ai[jj]] += Ax[jj] * x_col;
    }
  }
}

static void csc_left_normal_apply_vec(const int* Ai, const int* Ap,
                                      const double* Ax, int m, int n,
                                      const double* x, double* y,
                                      double* tmp_n) {
  csc_transpose_apply_vec(Ai, Ap, Ax, m, n, x, tmp_n);
  csc_forward_apply_vec(Ai, Ap, Ax, m, n, tmp_n, y);
}

static void csc_right_normal_apply_vec(const int* Ai, const int* Ap,
                                       const double* Ax, int m, int n,
                                       const double* x, double* y,
                                       double* tmp_m) {
  csc_forward_apply_vec(Ai, Ap, Ax, m, n, x, tmp_m);
  csc_transpose_apply_vec(Ai, Ap, Ax, m, n, tmp_m, y);
}

// Deterministic pseudo-random +-1 start vector entry. The raw LCG-style key
// has its low bit equal to the row parity, so it must be mixed (splitmix64
// finalizer) before taking a sign bit; otherwise the start vector is the
// alternating pattern, orthogonal to the constant vector for even m.
static inline double gram_lanczos_start_sign(int row, int rank) {
  uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(row + 1)) << 32) ^
    static_cast<uint64_t>(static_cast<uint32_t>(rank) * 2654435761u);
  key ^= key >> 30;
  key *= 0xbf58476d1ce4e5b9ULL;
  key ^= key >> 27;
  key *= 0x94d049bb133111ebULL;
  key ^= key >> 31;
  return (key >> 63) ? 1.0 : -1.0;
}

static int csc_implicit_left_normal_lanczos_attempt(const int* Ai,
                                                   const int* Ap,
                                                   const double* Ax,
                                                   int m,
                                                   int n,
                                                   int rank,
                                                   double tol,
                                                   double norm_A,
                                                   double* values,
                                                   double* U,
                                                   int* iterations_out,
                                                   double* max_backward_error_out) {
  if (m < 2 || rank < 1 || rank > m) {
    return 0;
  }
  int max_steps = std::max(43, 6 * rank + 13);
  if (max_steps > m) {
    max_steps = m;
  }
  if (max_steps < rank) {
    return 0;
  }

  std::vector<double> start(static_cast<size_t>(m), 0.0);
  for (int row = 0; row < m; ++row) {
    start[row] = gram_lanczos_start_sign(row, rank);
  }

  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(max_steps), 0.0);
  std::vector<double> z(static_cast<size_t>(m), 0.0);
  std::vector<double> tmp_m(static_cast<size_t>(m), 0.0);
  std::vector<double> tmp_n(static_cast<size_t>(n), 0.0);
  std::vector<double> alpha(static_cast<size_t>(max_steps), 0.0);
  std::vector<double> beta(static_cast<size_t>(max_steps), 0.0);
  std::memcpy(Q.data(), start.data(), sizeof(double) * static_cast<size_t>(m));
  double q_norm = trl_norm2_gram(Q.data(), m);
  if (q_norm <= 100.0 * DBL_EPSILON) {
    return 0;
  }
  for (int row = 0; row < m; ++row) {
    Q[row] /= q_norm;
  }

  int active = 0;
  for (int step = 0; step < max_steps; ++step) {
    eigencore_check_interrupt();
    const double* q = Q.data() + static_cast<int64_t>(step) * m;
    csc_left_normal_apply_vec(Ai, Ap, Ax, m, n, q, z.data(), tmp_n.data());

    long double dot = 0.0L;
    for (int row = 0; row < m; ++row) {
      dot += static_cast<long double>(q[row]) * z[static_cast<size_t>(row)];
    }
    alpha[static_cast<size_t>(step)] = static_cast<double>(dot);

    if (step > 0) {
      const double* q_prev = Q.data() + static_cast<int64_t>(step - 1) * m;
      const double b_prev = beta[static_cast<size_t>(step - 1)];
      for (int row = 0; row < m; ++row) {
        z[static_cast<size_t>(row)] -=
          alpha[static_cast<size_t>(step)] * q[row] + b_prev * q_prev[row];
      }
    } else {
      for (int row = 0; row < m; ++row) {
        z[static_cast<size_t>(row)] -= alpha[static_cast<size_t>(step)] * q[row];
      }
    }

    if ((step % 4) == 0 || step + rank >= max_steps) {
      trl_orthogonalise_gram(nullptr, 0, Q.data(), step + 1, z.data(), tmp_m.data(), m, 1);
    }
    const double b = trl_norm2_gram(z.data(), m);
    active = step + 1;
    if (step + 1 >= max_steps || b <= 100.0 * DBL_EPSILON * (norm_A * norm_A + 1.0)) {
      break;
    }
    beta[static_cast<size_t>(step)] = b;
    double* q_next = Q.data() + static_cast<int64_t>(step + 1) * m;
    const double inv_b = 1.0 / b;
    for (int row = 0; row < m; ++row) {
      q_next[row] = z[static_cast<size_t>(row)] * inv_b;
    }
  }
  if (active < rank) {
    return 0;
  }

  std::vector<double> T(static_cast<size_t>(active) * static_cast<size_t>(active), 0.0);
  std::vector<double> theta(static_cast<size_t>(active), 0.0);
  for (int col = 0; col < active; ++col) {
    T[col + static_cast<int64_t>(col) * active] = alpha[static_cast<size_t>(col)];
    if (col + 1 < active) {
      const double b = beta[static_cast<size_t>(col)];
      T[col + static_cast<int64_t>(col + 1) * active] = b;
      T[(col + 1) + static_cast<int64_t>(col) * active] = b;
    }
  }
  int lwork = trl_dsyev_query_gram(active);
  if (lwork < 3 * active) {
    lwork = 3 * active;
  }
  std::vector<double> work(static_cast<size_t>(lwork));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  F77_CALL(dsyev)(&jobz, &uplo, &active, T.data(), &active,
                  theta.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    return 0;
  }

  std::vector<double> Gu_exact(static_cast<size_t>(m), 0.0);
  const double scale_value_native = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double native_max_backward = 0.0;
  for (int out_col = 0; out_col < rank; ++out_col) {
    const int src_col = active - 1 - out_col;
    const double lambda_i = theta[static_cast<size_t>(src_col)] > 0.0 ?
      theta[static_cast<size_t>(src_col)] : 0.0;
    const double sigma = sqrt(lambda_i);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    double* u_col = U + static_cast<int64_t>(out_col) * m;
    std::fill(u_col, u_col + m, 0.0);
    for (int basis_col = 0; basis_col < active; ++basis_col) {
      const double coeff = T[basis_col + static_cast<int64_t>(src_col) * active];
      const double* q_col = Q.data() + static_cast<int64_t>(basis_col) * m;
      for (int row = 0; row < m; ++row) {
        u_col[row] += coeff * q_col[row];
      }
    }
    csc_left_normal_apply_vec(
      Ai, Ap, Ax, m, n, u_col, Gu_exact.data(), tmp_n.data()
    );
    long double residual2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double diff =
        (Gu_exact[static_cast<size_t>(row)] -
         lambda_i * u_col[row]) / sigma;
      residual2 += static_cast<long double>(diff) * diff;
    }
    const double backward = sqrt(static_cast<double>(residual2)) / scale_value_native;
    if (backward > native_max_backward) {
      native_max_backward = backward;
    }
    values[out_col] = lambda_i;
  }
  if (iterations_out != nullptr) {
    *iterations_out = active;
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = native_max_backward;
  }
  if (native_max_backward <= tol) {
    return 1;
  }
  return 0;
}

static int csc_implicit_right_normal_lanczos_attempt(const int* Ai,
                                                    const int* Ap,
                                                    const double* Ax,
                                                    int m,
                                                    int n,
                                                    int rank,
                                                    double tol,
                                                    double norm_A,
                                                    double* values,
                                                    double* V,
                                                    int* iterations_out,
                                                    double* max_backward_error_out) {
  if (n < 2 || rank < 1 || rank > n) {
    return 0;
  }
  int max_steps = std::max(38, 6 * rank + 8);
  if (max_steps > n) {
    max_steps = n;
  }
  if (max_steps < rank) {
    return 0;
  }

  std::vector<double> start(static_cast<size_t>(n), 0.0);
  for (int row = 0; row < n; ++row) {
    const uint32_t key = static_cast<uint32_t>((row + 1) * 1103515245u) ^
      static_cast<uint32_t>(rank * 2654435761u);
    start[row] = (key & 1u) ? 1.0 : -1.0;
  }

  std::vector<double> Q(static_cast<size_t>(n) * static_cast<size_t>(max_steps), 0.0);
  std::vector<double> z(static_cast<size_t>(n), 0.0);
  std::vector<double> tmp_n(static_cast<size_t>(n), 0.0);
  std::vector<double> tmp_m(static_cast<size_t>(m), 0.0);
  std::vector<double> alpha(static_cast<size_t>(max_steps), 0.0);
  std::vector<double> beta(static_cast<size_t>(max_steps), 0.0);
  std::memcpy(Q.data(), start.data(), sizeof(double) * static_cast<size_t>(n));
  double q_norm = trl_norm2_gram(Q.data(), n);
  if (q_norm <= 100.0 * DBL_EPSILON) {
    return 0;
  }
  for (int row = 0; row < n; ++row) {
    Q[row] /= q_norm;
  }

  int active = 0;
  double final_beta = 0.0;
  for (int step = 0; step < max_steps; ++step) {
    eigencore_check_interrupt();
    const double* q = Q.data() + static_cast<int64_t>(step) * n;
    csc_right_normal_apply_vec(Ai, Ap, Ax, m, n, q, z.data(), tmp_m.data());

    long double dot = 0.0L;
    for (int row = 0; row < n; ++row) {
      dot += static_cast<long double>(q[row]) * z[static_cast<size_t>(row)];
    }
    alpha[static_cast<size_t>(step)] = static_cast<double>(dot);

    if (step > 0) {
      const double* q_prev = Q.data() + static_cast<int64_t>(step - 1) * n;
      const double b_prev = beta[static_cast<size_t>(step - 1)];
      for (int row = 0; row < n; ++row) {
        z[static_cast<size_t>(row)] -=
          alpha[static_cast<size_t>(step)] * q[row] + b_prev * q_prev[row];
      }
    } else {
      for (int row = 0; row < n; ++row) {
        z[static_cast<size_t>(row)] -= alpha[static_cast<size_t>(step)] * q[row];
      }
    }

    if ((step % 4) == 0 || step + rank >= max_steps) {
      trl_orthogonalise_gram(nullptr, 0, Q.data(), step + 1, z.data(), tmp_n.data(), n, 1);
    }
    const double b = trl_norm2_gram(z.data(), n);
    final_beta = b;
    active = step + 1;
    if (step + 1 >= max_steps || b <= 100.0 * DBL_EPSILON * (norm_A * norm_A + 1.0)) {
      break;
    }
    beta[static_cast<size_t>(step)] = b;
    double* q_next = Q.data() + static_cast<int64_t>(step + 1) * n;
    const double inv_b = 1.0 / b;
    for (int row = 0; row < n; ++row) {
      q_next[row] = z[static_cast<size_t>(row)] * inv_b;
    }
  }
  if (active < rank) {
    return 0;
  }

  std::vector<double> T(static_cast<size_t>(active) * static_cast<size_t>(active), 0.0);
  std::vector<double> theta(static_cast<size_t>(active), 0.0);
  for (int col = 0; col < active; ++col) {
    T[col + static_cast<int64_t>(col) * active] = alpha[static_cast<size_t>(col)];
    if (col + 1 < active) {
      const double b = beta[static_cast<size_t>(col)];
      T[col + static_cast<int64_t>(col + 1) * active] = b;
      T[(col + 1) + static_cast<int64_t>(col) * active] = b;
    }
  }
  int lwork = trl_dsyev_query_gram(active);
  if (lwork < 3 * active) {
    lwork = 3 * active;
  }
  std::vector<double> work(static_cast<size_t>(lwork));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  F77_CALL(dsyev)(&jobz, &uplo, &active, T.data(), &active,
                  theta.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    return 0;
  }

  const double scale_value_native = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double native_max_backward = 0.0;
  for (int out_col = 0; out_col < rank; ++out_col) {
    const int src_col = active - 1 - out_col;
    const double lambda_i = theta[static_cast<size_t>(src_col)] > 0.0 ?
      theta[static_cast<size_t>(src_col)] : 0.0;
    const double sigma = sqrt(lambda_i);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    double* v_col = V + static_cast<int64_t>(out_col) * n;
    std::fill(v_col, v_col + n, 0.0);
    for (int basis_col = 0; basis_col < active; ++basis_col) {
      const double coeff = T[basis_col + static_cast<int64_t>(src_col) * active];
      const double* q_col = Q.data() + static_cast<int64_t>(basis_col) * n;
      for (int row = 0; row < n; ++row) {
        v_col[row] += coeff * q_col[row];
      }
    }
    const double last_component =
      T[(active - 1) + static_cast<int64_t>(src_col) * active];
    const double backward = fabs(final_beta * last_component) /
      (sigma * scale_value_native);
    if (backward > native_max_backward) {
      native_max_backward = backward;
    }
    values[out_col] = lambda_i;
  }
  if (iterations_out != nullptr) {
    *iterations_out = active;
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = native_max_backward;
  }
  return R_FINITE(native_max_backward) ? 1 : 0;
}

static int gram_krylov_left_normal_attempt(const double* gram,
                                           int m,
                                           int rank,
                                           double tol,
                                           double norm_A,
                                           double* values,
                                           double* U,
                                           int* iterations_out,
                                           double* max_backward_error_out) {
  if (m < 2 || rank < 1 || rank > m) {
    return 0;
  }
  int max_steps = std::max(45, 6 * rank + 15);
  if (max_steps > m) {
    max_steps = m;
  }
  if (max_steps < rank) {
    return 0;
  }

  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(max_steps), 0.0);
  std::vector<double> z(static_cast<size_t>(m), 0.0);
  std::vector<double> tmp(static_cast<size_t>(m), 0.0);
  std::vector<double> alpha(static_cast<size_t>(max_steps), 0.0);
  std::vector<double> beta(static_cast<size_t>(max_steps), 0.0);

  for (int row = 0; row < m; ++row) {
    Q[row] = gram_lanczos_start_sign(row, rank);
  }
  double q_norm = trl_norm2_gram(Q.data(), m);
  if (q_norm <= 100.0 * DBL_EPSILON) {
    return 0;
  }
  for (int row = 0; row < m; ++row) {
  Q[row] /= q_norm;
  }

  const char trans_N = 'N';
  const int inc_one = 1;
  const double one = 1.0;
  const double zero = 0.0;
  int active = 0;
  for (int step = 0; step < max_steps; ++step) {
    eigencore_check_interrupt();
    const double* q = Q.data() + static_cast<int64_t>(step) * m;
    F77_CALL(dgemv)(&trans_N, &m, &m, &one, gram, &m, q, &inc_one,
                    &zero, z.data(), &inc_one FCONE);
    long double dot = 0.0L;
    for (int row = 0; row < m; ++row) {
      dot += static_cast<long double>(q[row]) * z[static_cast<size_t>(row)];
    }
    alpha[static_cast<size_t>(step)] = static_cast<double>(dot);
    if (step > 0) {
      const double* q_prev = Q.data() + static_cast<int64_t>(step - 1) * m;
      const double b_prev = beta[static_cast<size_t>(step - 1)];
      for (int row = 0; row < m; ++row) {
        z[static_cast<size_t>(row)] -=
          alpha[static_cast<size_t>(step)] * q[row] + b_prev * q_prev[row];
      }
    } else {
      for (int row = 0; row < m; ++row) {
        z[static_cast<size_t>(row)] -= alpha[static_cast<size_t>(step)] * q[row];
      }
    }

    trl_orthogonalise_gram(nullptr, 0, Q.data(), step + 1, z.data(), tmp.data(), m);
    const double b = trl_norm2_gram(z.data(), m);
    active = step + 1;
    if (step + 1 >= max_steps || b <= 100.0 * DBL_EPSILON * (norm_A * norm_A + 1.0)) {
      break;
    }
    beta[static_cast<size_t>(step)] = b;
    double* q_next = Q.data() + static_cast<int64_t>(step + 1) * m;
    const double inv_b = 1.0 / b;
    for (int row = 0; row < m; ++row) {
      q_next[row] = z[static_cast<size_t>(row)] * inv_b;
    }
  }
  if (active < rank) {
    return 0;
  }

  std::vector<double> T(static_cast<size_t>(active) * static_cast<size_t>(active), 0.0);
  std::vector<double> theta(static_cast<size_t>(active), 0.0);
  for (int col = 0; col < active; ++col) {
    T[col + static_cast<int64_t>(col) * active] = alpha[static_cast<size_t>(col)];
    if (col + 1 < active) {
      const double b = beta[static_cast<size_t>(col)];
      T[col + static_cast<int64_t>(col + 1) * active] = b;
      T[(col + 1) + static_cast<int64_t>(col) * active] = b;
    }
  }
  int lwork = trl_dsyev_query_gram(active);
  if (lwork < 3 * active) {
    lwork = 3 * active;
  }
  std::vector<double> work(static_cast<size_t>(lwork));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  F77_CALL(dsyev)(&jobz, &uplo, &active, T.data(), &active,
                  theta.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    return 0;
  }

  std::vector<double> Gu(static_cast<size_t>(m), 0.0);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double max_backward = 0.0;
  for (int out_col = 0; out_col < rank; ++out_col) {
    const int src_col = active - 1 - out_col;
    const double lambda = theta[static_cast<size_t>(src_col)] > 0.0 ?
      theta[static_cast<size_t>(src_col)] : 0.0;
    const double sigma = sqrt(lambda);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    values[out_col] = lambda;
    double* u_col = U + static_cast<int64_t>(out_col) * m;
    std::fill(u_col, u_col + m, 0.0);
    for (int basis_col = 0; basis_col < active; ++basis_col) {
      const double coeff = T[basis_col + static_cast<int64_t>(src_col) * active];
      const double* q_col = Q.data() + static_cast<int64_t>(basis_col) * m;
      for (int row = 0; row < m; ++row) {
        u_col[row] += coeff * q_col[row];
      }
    }
    F77_CALL(dgemv)(&trans_N, &m, &m, &one, gram, &m, u_col, &inc_one,
                    &zero, Gu.data(), &inc_one FCONE);
    long double residual2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double residual = (Gu[static_cast<size_t>(row)] - lambda * u_col[row]) / sigma;
      residual2 += static_cast<long double>(residual) * residual;
    }
    const double backward = sqrt(static_cast<double>(residual2)) / scale_value;
    if (ISNAN(backward) || backward > max_backward) {
      max_backward = backward;
    }
  }
  if (iterations_out != nullptr) {
    *iterations_out = active;
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = max_backward;
  }
  return max_backward <= tol ? 1 : 0;
}

static int gram_top_subspace_attempt(const double* gram,
                                     int m,
                                     int rank,
                                     double tol,
                                     double norm_A,
                                     double* values,
                                     double* U,
                                     double* max_backward_error_out) {
  if (m < 2 || rank < 1 || rank > m) {
    return 0;
  }
  int subspace = rank + 5;
  if (subspace > m) subspace = m;
  if (subspace <= rank && rank < m) {
    return 0;
  }

  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(subspace), 0.0);
  std::vector<double> Z(static_cast<size_t>(m) * static_cast<size_t>(subspace), 0.0);
  std::vector<double> H(static_cast<size_t>(subspace) * static_cast<size_t>(subspace), 0.0);
  std::vector<double> theta(static_cast<size_t>(subspace), 0.0);
  std::vector<double> tmp(static_cast<size_t>(m), 0.0);

  for (int col = 0; col < subspace; ++col) {
    for (int row = 0; row < m; ++row) {
      const uint32_t key = static_cast<uint32_t>((row + 1) * 1103515245u) ^
        static_cast<uint32_t>((col + 1) * 2654435761u);
      Q[row + static_cast<int64_t>(col) * m] =
        ((key & 1u) ? 1.0 : -1.0) *
        (1.0 + static_cast<double>((key >> 8) & 15u) / 16.0);
    }
  }
  int q_rank = orthonormalize_columns_inplace(
    Q.data(), m, subspace, tmp.data(), 100.0 * DBL_EPSILON
  );
  if (q_rank < rank) {
    return 0;
  }
  subspace = q_rank;

  const char trans_N = 'N';
  const char trans_T = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  const int max_iter = 8;
  for (int iter = 0; iter < max_iter; ++iter) {
    eigencore_check_interrupt();
    F77_CALL(dgemm)(&trans_N, &trans_N, &m, &subspace, &m,
                    &one, gram, &m, Q.data(), &m,
                    &zero, Z.data(), &m FCONE FCONE);
    std::memcpy(Q.data(), Z.data(),
                sizeof(double) * static_cast<size_t>(m) *
                  static_cast<size_t>(subspace));
    q_rank = orthonormalize_columns_inplace(
      Q.data(), m, subspace, tmp.data(), 100.0 * DBL_EPSILON
    );
    if (q_rank < rank) {
      return 0;
    }
    subspace = q_rank;
  }

  F77_CALL(dgemm)(&trans_N, &trans_N, &m, &subspace, &m,
                  &one, gram, &m, Q.data(), &m,
                  &zero, Z.data(), &m FCONE FCONE);
  F77_CALL(dgemm)(&trans_T, &trans_N, &subspace, &subspace, &m,
                  &one, Q.data(), &m, Z.data(), &m,
                  &zero, H.data(), &subspace FCONE FCONE);
  symmetrize_packed_square_gram(H.data(), subspace);

  int lwork = trl_dsyev_query_gram(subspace);
  if (lwork < 3 * subspace) lwork = 3 * subspace;
  std::vector<double> work(static_cast<size_t>(lwork));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  F77_CALL(dsyev)(&jobz, &uplo, &subspace, H.data(), &subspace,
                  theta.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    return 0;
  }

  for (int out_col = 0; out_col < rank; ++out_col) {
    const int src_col = subspace - 1 - out_col;
    values[out_col] = theta[static_cast<size_t>(src_col)];
    for (int row = 0; row < m; ++row) {
      long double sum = 0.0L;
      for (int qcol = 0; qcol < subspace; ++qcol) {
        sum += static_cast<long double>(
          Q[row + static_cast<int64_t>(qcol) * m]
        ) * H[qcol + static_cast<int64_t>(src_col) * subspace];
      }
      U[row + static_cast<int64_t>(out_col) * m] = static_cast<double>(sum);
    }
  }

  std::vector<double> GU(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  F77_CALL(dgemm)(&trans_N, &trans_N, &m, &rank, &m,
                  &one, gram, &m, U, &m,
                  &zero, GU.data(), &m FCONE FCONE);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double max_backward = 0.0;
  for (int col = 0; col < rank; ++col) {
    const double lambda = values[col] > 0.0 ? values[col] : 0.0;
    const double sigma = sqrt(lambda);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    long double residual2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double residual = (GU[row + static_cast<int64_t>(col) * m] -
        lambda * U[row + static_cast<int64_t>(col) * m]) / sigma;
      residual2 += static_cast<long double>(residual) * residual;
    }
    const double backward = sqrt(static_cast<double>(residual2)) / scale_value;
    if (ISNAN(backward) || backward > max_backward) {
      max_backward = backward;
    }
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = max_backward;
  }
  return max_backward <= tol ? 1 : 0;
}

extern "C" SEXP eigencore_csc_left_gram_svd(SEXP i_, SEXP p_, SEXP x_,
                                            SEXP dim_, SEXP rank_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_gram_svd");
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  int rank = asInteger(rank_);
  if (m < 1 || n < 1 || rank < 1) {
    error("invalid CSC Gram SVD dimensions");
  }
  if (rank > m) {
    rank = m;
  }
  const double tol = asReal(tol_);
  const int* Ai = INTEGER(i_);
  const int* Ap = INTEGER(p_);
  const double* Ax = REAL(x_);

  auto stage_timer = native_timer_now();
  double stage_gram_seconds = 0.0;
  double stage_eigensolve_seconds = 0.0;
  double stage_vector_form_seconds = 0.0;
  double stage_diagnostics_seconds = 0.0;

  double frob2 = 0.0;
  for (int idx = 0; idx < Ap[n]; ++idx) {
    frob2 += Ax[idx] * Ax[idx];
  }
  const double norm_A = sqrt(frob2 > 0.0 ? frob2 : 0.0);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;

  std::vector<double> values(static_cast<size_t>(rank), 0.0);
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, rank));
  std::vector<double> gram;
  int used_implicit_lanczos = 0;
  int implicit_lanczos_iterations = 0;
  double implicit_lanczos_max_backward_error = R_PosInf;
  int used_subspace_eigensolve = 0;
  int used_gram_krylov = 0;
  int gram_krylov_iterations = 0;
  double subspace_max_backward_error = R_PosInf;
  const char* lapack_eigensolver = "lapack_dsyevr";

  SEXP implicit_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_normal_lanczos_attempt"));
  const int attempt_implicit_lanczos = asLogical(implicit_option_) == TRUE;
  if (attempt_implicit_lanczos && m <= 128 && rank <= 16) {
    stage_timer = native_timer_now();
    used_implicit_lanczos = csc_implicit_left_normal_lanczos_attempt(
      Ai, Ap, Ax, m, n, rank, tol, norm_A, values.data(), REAL(u_),
      &implicit_lanczos_iterations, &implicit_lanczos_max_backward_error
    );
    stage_eigensolve_seconds = native_timer_elapsed(stage_timer);
  }

  if (!used_implicit_lanczos) {
    stage_timer = native_timer_now();
    gram.assign(static_cast<size_t>(m) * static_cast<size_t>(m), 0.0);
    accumulate_sparse_gram_upper(m, n, Ap, Ai, Ax, gram.data());
    stage_gram_seconds = native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    SEXP krylov_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_gram_krylov_attempt"));
    const int attempt_gram_krylov = asLogical(krylov_option_) == TRUE;
    SEXP subspace_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_gram_subspace_attempt"));
    const int attempt_subspace_eigensolve = asLogical(subspace_option_) == TRUE;
    if (attempt_gram_krylov && m <= 90 && rank <= 8) {
      used_gram_krylov = gram_krylov_left_normal_attempt(
        gram.data(), m, rank, tol, norm_A, values.data(), REAL(u_),
        &gram_krylov_iterations, &subspace_max_backward_error
      );
    }
    if (!used_gram_krylov && attempt_subspace_eigensolve && m <= 128 && rank <= 16) {
      used_subspace_eigensolve = gram_top_subspace_attempt(
        gram.data(), m, rank, tol, norm_A, values.data(), REAL(u_),
        &subspace_max_backward_error
      );
    }
    if (m > 0 && !used_gram_krylov && !used_subspace_eigensolve) {
      std::vector<double> work_matrix(static_cast<size_t>(m) * static_cast<size_t>(m));
      std::vector<double> values_work(static_cast<size_t>(m), 0.0);
      char uplo = 'U';
      int info = 0;
      std::memcpy(work_matrix.data(), gram.data(),
                  sizeof(double) * static_cast<size_t>(m) * static_cast<size_t>(m));
      SEXP dsyevx_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_gram_dsyevx_attempt"));
      const int attempt_dsyevx = asLogical(dsyevx_option_) == TRUE;
      if (attempt_dsyevx && m <= 128 && rank <= 16) {
        lapack_eigensolver = "lapack_dsyevx";
        char jobz = 'V';
        char range = 'I';
        double vl = 0.0;
        double vu = 0.0;
        const double abstol = 0.0;
        int il = m - rank + 1;
        int iu = m;
        int m_found = 0;
        int lwork = 8 * m;
        if (lwork < 1) {
          lwork = 1;
        }
        std::vector<double> work(static_cast<size_t>(lwork));
        std::vector<int> iwork(static_cast<size_t>(5 * m));
        std::vector<int> ifail(static_cast<size_t>(m));
        F77_CALL(dsyevx)(&jobz, &range, &uplo, &m, work_matrix.data(), &m,
                         &vl, &vu, &il, &iu, &abstol,
                         &m_found, values_work.data(), REAL(u_), &m,
                         work.data(), &lwork, iwork.data(), ifail.data(),
                         &info FCONE FCONE FCONE);
        if (info != 0 || m_found != rank) {
          UNPROTECT(1);
          error("LAPACK dsyevx failed with info=%d, found=%d", info, m_found);
        }
        for (int left = 0, right = rank - 1; left < right; ++left, --right) {
          const double tmp_value = values_work[left];
          values_work[left] = values_work[right];
          values_work[right] = tmp_value;
          for (int row = 0; row < m; ++row) {
            const int64_t lpos = row + static_cast<int64_t>(left) * m;
            const int64_t rpos = row + static_cast<int64_t>(right) * m;
            const double tmp_vec = REAL(u_)[lpos];
            REAL(u_)[lpos] = REAL(u_)[rpos];
            REAL(u_)[rpos] = tmp_vec;
          }
        }
        for (int col = 0; col < rank; ++col) {
          values[static_cast<size_t>(col)] = values_work[static_cast<size_t>(col)];
        }
      } else if (m <= 128 && rank >= 16) {
        lapack_eigensolver = "lapack_dsyevd";
        char jobz = 'V';
        int lwork = -1;
        int liwork = -1;
        double work_query = 0.0;
        int iwork_query = 0;
        F77_CALL(dsyevd)(&jobz, &uplo, &m, work_matrix.data(), &m,
                         values_work.data(), &work_query, &lwork,
                         &iwork_query, &liwork, &info FCONE FCONE);
        if (info != 0) {
          UNPROTECT(1);
          error("LAPACK dsyevd workspace query failed with info=%d", info);
        }
        lwork = static_cast<int>(work_query);
        liwork = iwork_query;
        if (lwork < 1 + 6 * m + 2 * m * m) {
          lwork = 1 + 6 * m + 2 * m * m;
        }
        if (liwork < 3 + 5 * m) {
          liwork = 3 + 5 * m;
        }
        std::vector<double> work(static_cast<size_t>(lwork));
        std::vector<int> iwork(static_cast<size_t>(liwork));
        F77_CALL(dsyevd)(&jobz, &uplo, &m, work_matrix.data(), &m,
                         values_work.data(), work.data(), &lwork,
                         iwork.data(), &liwork, &info FCONE FCONE);
        if (info != 0) {
          UNPROTECT(1);
          error("LAPACK dsyevd failed with info=%d", info);
        }
        for (int col = 0; col < rank; ++col) {
          const int source_col = m - 1 - col;
          values[static_cast<size_t>(col)] = values_work[static_cast<size_t>(source_col)];
          std::memcpy(
            REAL(u_) + static_cast<int64_t>(col) * m,
            work_matrix.data() + static_cast<int64_t>(source_col) * m,
            sizeof(double) * static_cast<size_t>(m)
          );
        }
      } else {
        std::vector<int> isuppz(static_cast<size_t>(2 * rank), 0);
        char jobz = 'V';
        char range = 'I';
        double vl = 0.0;
        double vu = 0.0;
        const double abstol = 0.0;
        int il = m - rank + 1;
        int iu = m;
        int m_found = 0;
        int lwork = 26 * m;
        int liwork = 10 * m;
        std::vector<double> work(static_cast<size_t>(lwork));
        std::vector<int> iwork(static_cast<size_t>(liwork));
        F77_CALL(dsyevr)(&jobz, &range, &uplo, &m, work_matrix.data(), &m,
                         &vl, &vu, &il, &iu, &abstol,
                         &m_found, values_work.data(), REAL(u_), &m,
                         isuppz.data(), work.data(), &lwork,
                         iwork.data(), &liwork, &info FCONE FCONE FCONE);
        if (info != 0 || m_found != rank) {
          UNPROTECT(1);
          error("LAPACK dsyevr failed with info=%d, found=%d", info, m_found);
        }
        for (int left = 0, right = rank - 1; left < right; ++left, --right) {
          const double tmp_value = values_work[left];
          values_work[left] = values_work[right];
          values_work[right] = tmp_value;
          for (int row = 0; row < m; ++row) {
            const int64_t lpos = row + static_cast<int64_t>(left) * m;
            const int64_t rpos = row + static_cast<int64_t>(right) * m;
            const double tmp_vec = REAL(u_)[lpos];
            REAL(u_)[lpos] = REAL(u_)[rpos];
            REAL(u_)[rpos] = tmp_vec;
          }
        }
        for (int col = 0; col < rank; ++col) {
          values[static_cast<size_t>(col)] = values_work[static_cast<size_t>(col)];
        }
      }
      subspace_max_backward_error = R_PosInf;
    }
    stage_eigensolve_seconds = native_timer_elapsed(stage_timer);
  }

  stage_timer = native_timer_now();
  SEXP d_ = PROTECT(allocVector(REALSXP, rank));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, rank));

  for (int col = 0; col < rank; ++col) {
    const double sigma = sqrt(values[static_cast<size_t>(col)] > 0.0 ?
                              values[static_cast<size_t>(col)] : 0.0);
    REAL(d_)[col] = sigma;
  }

  // A^T U is formed once from A and the final U: it yields V = A^T U diag(1/d)
  // and is the exact adjoint product for the right certificate residual.
  std::vector<double> atu_form(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  {
    CSCOperator form_impl = {m, n, Ai, Ap, Ax};
    eigencore_csc_apply(&form_impl, EIGENCORE_TRANSPOSE_ADJOINT, rank, REAL(u_), m,
                        1.0, 0.0, atu_form.data(), n, nullptr);
  }
  for (int scol = 0; scol < rank; ++scol) {
    const double sigma = REAL(d_)[scol];
    const double inv_sigma =
      sigma > 100.0 * DBL_EPSILON * scale_value ? 1.0 / sigma : 0.0;
    const double* atu_col = atu_form.data() + static_cast<int64_t>(scol) * n;
    double* v_col = REAL(v_) + static_cast<int64_t>(scol) * n;
    for (int acol = 0; acol < n; ++acol) {
      v_col[acol] = atu_col[acol] * inv_sigma;
    }
  }
  stage_vector_form_seconds = native_timer_elapsed(stage_timer);

  stage_timer = native_timer_now();
  SEXP left_ = PROTECT(allocVector(REALSXP, rank));
  SEXP right_ = PROTECT(allocVector(REALSXP, rank));
  SEXP combined_ = PROTECT(allocVector(REALSXP, rank));
  SEXP backward_ = PROTECT(allocVector(REALSXP, rank));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, rank));
  SEXP scale_ = PROTECT(allocVector(REALSXP, rank));
  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));
  for (int col = 0; col < rank; ++col) {
    REAL(scale_)[col] = scale_value;
  }
  // Certify in original coordinates: both residuals from A itself and the
  // orthogonality of the returned U and V (never inferred from the Gram).
  csc_two_sided_svd_certificate(
    Ai, Ap, Ax, m, n, rank, REAL(d_), REAL(u_), REAL(v_),
    nullptr, atu_form.data(), scale_value, tol,
    REAL(left_), REAL(right_), REAL(combined_), REAL(backward_),
    LOGICAL(converged_), REAL(orth_), REAL(orth_) + 1
  );
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);
  stage_diagnostics_seconds = native_timer_elapsed(stage_timer);

  SEXP diagnostics_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(diagnostics_, 0, left_);
  SET_VECTOR_ELT(diagnostics_, 1, right_);
  SET_VECTOR_ELT(diagnostics_, 2, combined_);
  SET_VECTOR_ELT(diagnostics_, 3, backward_);
  SET_VECTOR_ELT(diagnostics_, 4, orth_);
  SET_VECTOR_ELT(diagnostics_, 5, converged_);
  SET_VECTOR_ELT(diagnostics_, 6, scale_);
  SEXP diag_names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(diag_names_, 0, mkChar("left"));
  SET_STRING_ELT(diag_names_, 1, mkChar("right"));
  SET_STRING_ELT(diag_names_, 2, mkChar("combined"));
  SET_STRING_ELT(diag_names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(diag_names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(diag_names_, 5, mkChar("converged"));
  SET_STRING_ELT(diag_names_, 6, mkChar("scale"));
  setAttrib(diagnostics_, R_NamesSymbol, diag_names_);

  SEXP stage_ = PROTECT(allocVector(REALSXP, 4));
  REAL(stage_)[0] = stage_gram_seconds;
  REAL(stage_)[1] = stage_eigensolve_seconds;
  REAL(stage_)[2] = stage_vector_form_seconds;
  REAL(stage_)[3] = stage_diagnostics_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(stage_names_, 0, mkChar("gram"));
  SET_STRING_ELT(stage_names_, 1, mkChar("eigensolve"));
  SET_STRING_ELT(stage_names_, 2, mkChar("vector_form"));
  SET_STRING_ELT(stage_names_, 3, mkChar("diagnostics"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP eigensolver_ = PROTECT(mkString(
    used_implicit_lanczos ? "implicit_normal_lanczos" :
      (used_gram_krylov ? "explicit_gram_krylov" :
        (used_subspace_eigensolve ? "subspace_iteration" : lapack_eigensolver))
  ));
  SEXP subspace_be_ = PROTECT(ScalarReal(subspace_max_backward_error));
  SEXP implicit_be_ = PROTECT(ScalarReal(implicit_lanczos_max_backward_error));
  SEXP implicit_iter_ = PROTECT(ScalarInteger(implicit_lanczos_iterations));
  SEXP gram_krylov_iter_ = PROTECT(ScalarInteger(gram_krylov_iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 10));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, diagnostics_);
  SET_VECTOR_ELT(out_, 4, stage_);
  SET_VECTOR_ELT(out_, 5, eigensolver_);
  SET_VECTOR_ELT(out_, 6, subspace_be_);
  SET_VECTOR_ELT(out_, 7, implicit_be_);
  SET_VECTOR_ELT(out_, 8, implicit_iter_);
  SET_VECTOR_ELT(out_, 9, gram_krylov_iter_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 10));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("diagnostics"));
  SET_STRING_ELT(names_, 4, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 5, mkChar("eigensolver"));
  SET_STRING_ELT(names_, 6, mkChar("subspace_max_backward_error"));
  SET_STRING_ELT(names_, 7, mkChar("implicit_lanczos_max_backward_error"));
  SET_STRING_ELT(names_, 8, mkChar("implicit_lanczos_iterations"));
  SET_STRING_ELT(names_, 9, mkChar("gram_krylov_iterations"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(22);
  return out_;
  EIGENCORE_ENTRY_END
}

static SEXP eigencore_csc_right_gram_svd_impl(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP rank_, SEXP tol_,
    int allow_implicit_lanczos,
    double prior_gram_seconds,
    double prior_eigensolve_seconds,
    double prior_vector_form_seconds,
    double prior_diagnostics_seconds,
    double prior_implicit_max_backward_error,
    int prior_implicit_iterations) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_gram_svd");
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  int rank = asInteger(rank_);
  if (m < 1 || n < 1 || rank < 1) {
    error("invalid CSC Gram SVD dimensions");
  }
  if (rank > n) {
    rank = n;
  }
  const double tol = asReal(tol_);
  const int* Ai = INTEGER(i_);
  const int* Ap = INTEGER(p_);
  const double* Ax = REAL(x_);

  auto stage_timer = native_timer_now();
  double stage_gram_seconds = prior_gram_seconds;
  double stage_eigensolve_seconds = prior_eigensolve_seconds;
  double stage_vector_form_seconds = prior_vector_form_seconds;
  double stage_diagnostics_seconds = prior_diagnostics_seconds;

  double frob2 = 0.0;
  for (int idx = 0; idx < Ap[n]; ++idx) {
    frob2 += Ax[idx] * Ax[idx];
  }
  const double norm_A = sqrt(frob2 > 0.0 ? frob2 : 0.0);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;

  std::vector<double> gram;
  std::vector<double> values_work(static_cast<size_t>(n), 0.0);
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, rank));
  int used_implicit_lanczos = 0;
  int implicit_lanczos_iterations = prior_implicit_iterations;
  double implicit_lanczos_max_backward_error =
    prior_implicit_max_backward_error;
  const char* lapack_eigensolver = "lapack_dsyevr";

  SEXP implicit_option_ = Rf_GetOption1(Rf_install("eigencore.csc_right_normal_lanczos_attempt"));
  const int attempt_implicit_lanczos = allow_implicit_lanczos &&
    (asLogical(implicit_option_) == TRUE);
  if (attempt_implicit_lanczos && n <= 128 && rank <= 16) {
    stage_timer = native_timer_now();
    used_implicit_lanczos = csc_implicit_right_normal_lanczos_attempt(
      Ai, Ap, Ax, m, n, rank, tol, norm_A, values_work.data(), REAL(v_),
      &implicit_lanczos_iterations, &implicit_lanczos_max_backward_error
    );
    stage_eigensolve_seconds += native_timer_elapsed(stage_timer);
  }

  if (!used_implicit_lanczos) {
    stage_timer = native_timer_now();
    gram.assign(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    std::vector<int> row_counts(static_cast<size_t>(m), 0);
    for (int idx = 0; idx < Ap[n]; ++idx) {
      const int row = Ai[idx];
      if (row >= 0 && row < m) {
        ++row_counts[static_cast<size_t>(row)];
      }
    }
    std::vector<int> row_ptr(static_cast<size_t>(m + 1), 0);
    for (int row = 0; row < m; ++row) {
      row_ptr[static_cast<size_t>(row + 1)] =
        row_ptr[static_cast<size_t>(row)] + row_counts[static_cast<size_t>(row)];
    }
    std::vector<int> row_next = row_ptr;
    const int nnz = Ap[n];
    std::vector<int> row_cols(static_cast<size_t>(nnz), 0);
    std::vector<double> row_vals(static_cast<size_t>(nnz), 0.0);
    for (int col = 0; col < n; ++col) {
      for (int jj = Ap[col]; jj < Ap[col + 1]; ++jj) {
        const int row = Ai[jj];
        const int pos = row_next[static_cast<size_t>(row)]++;
        row_cols[static_cast<size_t>(pos)] = col;
        row_vals[static_cast<size_t>(pos)] = Ax[jj];
      }
    }
    accumulate_sparse_gram_upper(n, m, row_ptr.data(), row_cols.data(),
                                 row_vals.data(), gram.data());
    stage_gram_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    std::vector<double> work_matrix(static_cast<size_t>(n) * static_cast<size_t>(n));
    std::memcpy(work_matrix.data(), gram.data(),
                sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
    char jobz = 'V';
    char range = 'I';
    char uplo = 'U';
    double vl = 0.0;
    double vu = 0.0;
    const double abstol = 0.0;
    int il = n - rank + 1;
    int iu = n;
    int m_found = 0;
    int info = 0;
    int lwork = 26 * n;
    int liwork = 10 * n;
    std::vector<int> isuppz(static_cast<size_t>(2 * rank), 0);
    std::vector<double> work(static_cast<size_t>(lwork));
    std::vector<int> iwork(static_cast<size_t>(liwork));
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, work_matrix.data(), &n,
                     &vl, &vu, &il, &iu, &abstol,
                     &m_found, values_work.data(), REAL(v_), &n,
                     isuppz.data(), work.data(), &lwork,
                     iwork.data(), &liwork, &info FCONE FCONE FCONE);
    if (info != 0 || m_found != rank) {
      UNPROTECT(1);
      error("LAPACK dsyevr failed with info=%d, found=%d", info, m_found);
    }
    for (int left = 0, right = rank - 1; left < right; ++left, --right) {
      const double tmp_value = values_work[left];
      values_work[left] = values_work[right];
      values_work[right] = tmp_value;
      for (int row = 0; row < n; ++row) {
        const int64_t lpos = row + static_cast<int64_t>(left) * n;
        const int64_t rpos = row + static_cast<int64_t>(right) * n;
        const double tmp_vec = REAL(v_)[lpos];
        REAL(v_)[lpos] = REAL(v_)[rpos];
        REAL(v_)[rpos] = tmp_vec;
      }
    }
    stage_eigensolve_seconds += native_timer_elapsed(stage_timer);
  }

  stage_timer = native_timer_now();
  SEXP d_ = PROTECT(allocVector(REALSXP, rank));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, rank));
  std::memset(REAL(u_), 0, sizeof(double) * static_cast<size_t>(m) * rank);
  for (int col = 0; col < rank; ++col) {
    const double sigma = sqrt(values_work[static_cast<size_t>(col)] > 0.0 ?
                              values_work[static_cast<size_t>(col)] : 0.0);
    REAL(d_)[col] = sigma;
  }
  // A V is formed once from A and the final V: it yields U = A V diag(1/d)
  // and is the exact forward product for the left certificate residual.
  std::vector<double> av_form(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  {
    CSCOperator form_impl = {m, n, Ai, Ap, Ax};
    eigencore_csc_apply(&form_impl, EIGENCORE_TRANSPOSE_NONE, rank, REAL(v_), n,
                        1.0, 0.0, av_form.data(), m, nullptr);
  }
  for (int scol = 0; scol < rank; ++scol) {
    const double sigma = REAL(d_)[scol];
    const double inv_sigma =
      sigma > 100.0 * DBL_EPSILON * scale_value ? 1.0 / sigma : 0.0;
    const double* av_col = av_form.data() + static_cast<int64_t>(scol) * m;
    double* u_col = REAL(u_) + static_cast<int64_t>(scol) * m;
    for (int row = 0; row < m; ++row) {
      u_col[row] = av_col[row] * inv_sigma;
    }
  }
  stage_vector_form_seconds += native_timer_elapsed(stage_timer);

  stage_timer = native_timer_now();
  SEXP left_ = PROTECT(allocVector(REALSXP, rank));
  SEXP right_ = PROTECT(allocVector(REALSXP, rank));
  SEXP combined_ = PROTECT(allocVector(REALSXP, rank));
  SEXP backward_ = PROTECT(allocVector(REALSXP, rank));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, rank));
  SEXP scale_ = PROTECT(allocVector(REALSXP, rank));
  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));
  for (int col = 0; col < rank; ++col) {
    REAL(scale_)[col] = scale_value;
  }
  // Certify in original coordinates: both residuals from A itself and the
  // orthogonality of the returned U and V (never inferred from the Gram).
  csc_two_sided_svd_certificate(
    Ai, Ap, Ax, m, n, rank, REAL(d_), REAL(u_), REAL(v_),
    av_form.data(), nullptr, scale_value, tol,
    REAL(left_), REAL(right_), REAL(combined_), REAL(backward_),
    LOGICAL(converged_), REAL(orth_), REAL(orth_) + 1
  );
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);
  if (used_implicit_lanczos) {
    double max_backward = 0.0;
    for (int scol = 0; scol < rank; ++scol) {
      if (ISNAN(REAL(backward_)[scol]) || REAL(backward_)[scol] > max_backward || scol == 0) {
        max_backward = REAL(backward_)[scol];
      }
    }
    implicit_lanczos_max_backward_error = max_backward;
  }
  stage_diagnostics_seconds += native_timer_elapsed(stage_timer);

  if (used_implicit_lanczos) {
    int certificate_passed = TRUE;
    for (int scol = 0; scol < rank; ++scol) {
      if (LOGICAL(converged_)[scol] != TRUE) {
        certificate_passed = FALSE;
        break;
      }
    }
    const double orth_tol = std::max(tol, sqrt(DBL_EPSILON));
    if (!(REAL(orth_)[0] <= orth_tol) || !(REAL(orth_)[1] <= orth_tol)) {
      certificate_passed = FALSE;
    }
    if (!certificate_passed) {
      const double failed_implicit_max_backward_error =
        implicit_lanczos_max_backward_error;
      const int failed_implicit_iterations = implicit_lanczos_iterations;
      UNPROTECT(11);
      return eigencore_csc_right_gram_svd_impl(
        i_, p_, x_, dim_, rank_, tol_, FALSE,
        stage_gram_seconds,
        stage_eigensolve_seconds,
        stage_vector_form_seconds,
        stage_diagnostics_seconds,
        failed_implicit_max_backward_error,
        failed_implicit_iterations
      );
    }
  }

  SEXP diagnostics_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(diagnostics_, 0, left_);
  SET_VECTOR_ELT(diagnostics_, 1, right_);
  SET_VECTOR_ELT(diagnostics_, 2, combined_);
  SET_VECTOR_ELT(diagnostics_, 3, backward_);
  SET_VECTOR_ELT(diagnostics_, 4, orth_);
  SET_VECTOR_ELT(diagnostics_, 5, converged_);
  SET_VECTOR_ELT(diagnostics_, 6, scale_);
  SEXP diag_names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(diag_names_, 0, mkChar("left"));
  SET_STRING_ELT(diag_names_, 1, mkChar("right"));
  SET_STRING_ELT(diag_names_, 2, mkChar("combined"));
  SET_STRING_ELT(diag_names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(diag_names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(diag_names_, 5, mkChar("converged"));
  SET_STRING_ELT(diag_names_, 6, mkChar("scale"));
  setAttrib(diagnostics_, R_NamesSymbol, diag_names_);

  SEXP stage_ = PROTECT(allocVector(REALSXP, 4));
  REAL(stage_)[0] = stage_gram_seconds;
  REAL(stage_)[1] = stage_eigensolve_seconds;
  REAL(stage_)[2] = stage_vector_form_seconds;
  REAL(stage_)[3] = stage_diagnostics_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(stage_names_, 0, mkChar("gram"));
  SET_STRING_ELT(stage_names_, 1, mkChar("eigensolve"));
  SET_STRING_ELT(stage_names_, 2, mkChar("vector_form"));
  SET_STRING_ELT(stage_names_, 3, mkChar("diagnostics"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP eigensolver_ = PROTECT(mkString(
    used_implicit_lanczos ? "implicit_normal_lanczos" : lapack_eigensolver
  ));
  SEXP subspace_be_ = PROTECT(ScalarReal(R_PosInf));
  SEXP implicit_be_ = PROTECT(ScalarReal(implicit_lanczos_max_backward_error));
  SEXP implicit_iter_ = PROTECT(ScalarInteger(implicit_lanczos_iterations));
  SEXP gram_krylov_iter_ = PROTECT(ScalarInteger(0));

  SEXP out_ = PROTECT(allocVector(VECSXP, 10));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, diagnostics_);
  SET_VECTOR_ELT(out_, 4, stage_);
  SET_VECTOR_ELT(out_, 5, eigensolver_);
  SET_VECTOR_ELT(out_, 6, subspace_be_);
  SET_VECTOR_ELT(out_, 7, implicit_be_);
  SET_VECTOR_ELT(out_, 8, implicit_iter_);
  SET_VECTOR_ELT(out_, 9, gram_krylov_iter_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 10));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("diagnostics"));
  SET_STRING_ELT(names_, 4, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 5, mkChar("eigensolver"));
  SET_STRING_ELT(names_, 6, mkChar("subspace_max_backward_error"));
  SET_STRING_ELT(names_, 7, mkChar("implicit_lanczos_max_backward_error"));
  SET_STRING_ELT(names_, 8, mkChar("implicit_lanczos_iterations"));
  SET_STRING_ELT(names_, 9, mkChar("gram_krylov_iterations"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(22);
  return out_;
}

extern "C" SEXP eigencore_csc_right_gram_svd(SEXP i_, SEXP p_, SEXP x_,
                                             SEXP dim_, SEXP rank_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  return eigencore_csc_right_gram_svd_impl(
    i_, p_, x_, dim_, rank_, tol_, TRUE,
    0.0, 0.0, 0.0, 0.0, R_PosInf, 0
  );
  EIGENCORE_ENTRY_END
}

static SEXP eigencore_csc_gram_svd_fast_result_from_native(SEXP native_,
                                                           int m, int n,
                                                           double tol,
                                                           const char* gram_side,
                                                           const char* native_kernel) {
  const int rank = LENGTH(VECTOR_ELT(native_, 0));
  SEXP gram_max_option_ = Rf_GetOption1(Rf_install("eigencore.gram_svd_max_dimension"));
  int gram_max_dimension = asInteger(gram_max_option_);
  if (gram_max_dimension == NA_INTEGER) {
    gram_max_dimension = 512;
  }

  SEXP d_ = VECTOR_ELT(native_, 0);
  SEXP u_ = VECTOR_ELT(native_, 1);
  SEXP v_ = VECTOR_ELT(native_, 2);
  SEXP diagnostics_ = VECTOR_ELT(native_, 3);
  SEXP stage_ = VECTOR_ELT(native_, 4);
  SEXP eigensolver_ = VECTOR_ELT(native_, 5);
  SEXP subspace_be_ = VECTOR_ELT(native_, 6);
  SEXP implicit_be_ = VECTOR_ELT(native_, 7);
  SEXP implicit_iter_ = VECTOR_ELT(native_, 8);
  SEXP gram_krylov_iter_ = VECTOR_ELT(native_, 9);

  // Scale-invariant cutoff: relative to the largest computed singular value
  // (an absolute floor of 1 would classify every value of a small-scale
  // matrix as numerically zero). An all-zero spectrum falls back to 1.
  double max_d = 0.0;
  for (int i = 0; i < rank; ++i) {
    const double di = REAL(d_)[i];
    if (R_finite(di) && di > max_d) max_d = di;
  }
  if (!(max_d > 0.0)) {
    max_d = 1.0;
  }
  const double zero_tol = std::max(
    std::max(100.0 * DBL_EPSILON * max_d, sqrt(DBL_EPSILON) * max_d),
    tol * max_d * 1e-3
  );
  for (int i = 0; i < rank; ++i) {
    if (REAL(d_)[i] <= zero_tol) {
      return R_NilValue;
    }
  }

  SEXP left_ = VECTOR_ELT(diagnostics_, 0);
  SEXP right_ = VECTOR_ELT(diagnostics_, 1);
  SEXP combined_ = VECTOR_ELT(diagnostics_, 2);
  SEXP backward_ = VECTOR_ELT(diagnostics_, 3);
  SEXP orth_ = VECTOR_ELT(diagnostics_, 4);
  SEXP converged_ = VECTOR_ELT(diagnostics_, 5);
  SEXP scale_ = VECTOR_ELT(diagnostics_, 6);

  double max_backward = 0.0;
  double max_residual = 0.0;
  for (int i = 0; i < rank; ++i) {
    if (ISNAN(REAL(backward_)[i]) || REAL(backward_)[i] > max_backward) max_backward = REAL(backward_)[i];
    if (ISNAN(REAL(combined_)[i]) || REAL(combined_)[i] > max_residual) max_residual = REAL(combined_)[i];
  }
  double max_orth = NA_REAL;
  if (LENGTH(orth_) > 0) {
    max_orth = REAL(orth_)[0];
    for (int i = 1; i < LENGTH(orth_); ++i) {
      if (ISNAN(REAL(orth_)[i]) || REAL(orth_)[i] > max_orth) max_orth = REAL(orth_)[i];
    }
  }
  const double orth_tol = std::max(tol, sqrt(DBL_EPSILON));
  int all_converged = TRUE;
  int failed_count = 0;
  for (int i = 0; i < rank; ++i) {
    if (LOGICAL(converged_)[i] != TRUE) {
      all_converged = FALSE;
      ++failed_count;
    }
  }
  const int orth_passed = (ISNA(max_orth) || max_orth <= orth_tol) ? TRUE : FALSE;
  if (!(all_converged && orth_passed)) {
    return R_NilValue;
  }

  SEXP residuals_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(residuals_, 0, left_);
  SET_VECTOR_ELT(residuals_, 1, right_);
  SET_VECTOR_ELT(residuals_, 2, combined_);
  SEXP residual_names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(residual_names_, 0, mkChar("left"));
  SET_STRING_ELT(residual_names_, 1, mkChar("right"));
  SET_STRING_ELT(residual_names_, 2, mkChar("combined"));
  setAttrib(residuals_, R_NamesSymbol, residual_names_);

  SEXP failed_ = PROTECT(allocVector(INTSXP, failed_count));
  for (int i = 0, out = 0; i < rank; ++i) {
    if (LOGICAL(converged_)[i] != TRUE) {
      INTEGER(failed_)[out++] = i + 1;
    }
  }

  SEXP cert_ = PROTECT(allocVector(VECSXP, 18));
  SET_VECTOR_ELT(cert_, 0, ScalarLogical(all_converged && orth_passed));
  SET_VECTOR_ELT(cert_, 1, ScalarReal(tol));
  SET_VECTOR_ELT(cert_, 2, ScalarReal(orth_tol));
  SET_VECTOR_ELT(cert_, 3, ScalarLogical(TRUE));
  SET_VECTOR_ELT(cert_, 4, mkString("residual_backward_error"));
  SET_VECTOR_ELT(cert_, 5, mkString("frobenius_exact"));
  SET_VECTOR_ELT(cert_, 6, ScalarLogical(FALSE));
  SET_VECTOR_ELT(cert_, 7, ScalarReal(max_backward));
  SET_VECTOR_ELT(cert_, 8, ScalarReal(max_residual));
  SET_VECTOR_ELT(cert_, 9, ScalarReal(max_orth));
  SET_VECTOR_ELT(cert_, 10, ScalarLogical(orth_passed));
  SET_VECTOR_ELT(cert_, 11, failed_);
  SET_VECTOR_ELT(cert_, 12, scale_);
  SET_VECTOR_ELT(cert_, 13, allocVector(STRSXP, 0));
  SET_VECTOR_ELT(cert_, 14, residuals_);
  SET_VECTOR_ELT(cert_, 15, backward_);
  SET_VECTOR_ELT(cert_, 16, orth_);
  SET_VECTOR_ELT(cert_, 17, converged_);
  SEXP cert_names_ = PROTECT(allocVector(STRSXP, 18));
  const char* cert_names[] = {
    "passed", "tolerance", "orthogonality_tolerance",
    "orthogonality_required", "certificate_type", "norm_bound_type",
    "scale_is_estimate", "max_backward_error", "max_residual",
    "max_orthogonality_loss", "orthogonality_passed", "failed_indices",
    "scale", "notes", "residuals", "backward_error", "orthogonality",
    "converged"
  };
  for (int i = 0; i < 18; ++i) SET_STRING_ELT(cert_names_, i, mkChar(cert_names[i]));
  setAttrib(cert_, R_NamesSymbol, cert_names_);
  SEXP cert_class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(cert_class_, 0, mkChar("eigencore_certificate"));
  setAttrib(cert_, R_ClassSymbol, cert_class_);

  SEXP controls_ = PROTECT(allocVector(VECSXP, 11));
  SET_VECTOR_ELT(controls_, 0, mkString(gram_side));
  SET_VECTOR_ELT(controls_, 1, ScalarInteger(m < n ? m : n));
  SET_VECTOR_ELT(controls_, 2, ScalarInteger(gram_max_dimension));
  SET_VECTOR_ELT(controls_, 3, ScalarReal(0.5));
  SET_VECTOR_ELT(controls_, 4, ScalarLogical(TRUE));
  SET_VECTOR_ELT(controls_, 5, mkString(
    "smaller Gram matrix by default; implicit normal operator only by explicit opt-in"
  ));
  SET_VECTOR_ELT(controls_, 6, mkString("certification-gated"));
  SET_VECTOR_ELT(controls_, 7, mkString(
    "opt-in implicit candidate retries explicit Gram; native Golub-Kahan if still uncertified"
  ));
  SET_VECTOR_ELT(controls_, 8, mkString("both"));
  SET_VECTOR_ELT(controls_, 9, ScalarLogical(TRUE));
  SET_VECTOR_ELT(controls_, 10, ScalarInteger(m > n ? m : n));
  SEXP control_names_ = PROTECT(allocVector(STRSXP, 11));
  const char* control_names[] = {
    "gram_side", "gram_dimension", "gram_max_dimension",
    "rank_fraction_limit", "certified_in_original_coordinates",
    "materializes", "fallback_policy", "runtime_fallback",
    "fallback_requires_vectors", "svd_partial_fastpath", "full_dimension"
  };
  for (int i = 0; i < 11; ++i) SET_STRING_ELT(control_names_, i, mkChar(control_names[i]));
  setAttrib(controls_, R_NamesSymbol, control_names_);

  SEXP reasons_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(reasons_, 0, mkChar("target: largest"));
  SET_STRING_ELT(reasons_, 1, mkChar("rectangular SVD problem"));
  SET_STRING_ELT(reasons_, 2, mkChar("adjoint is available"));
  SET_STRING_ELT(reasons_, 3, mkChar(
    "bounded smaller-side normal problem with exact original-coordinate certification; explicit Gram is the production default"
  ));
  SET_STRING_ELT(reasons_, 4, mkChar("built-in sparse CSC operator has native block apply"));
  SET_STRING_ELT(reasons_, 5, mkChar("direct svd_partial() fast path avoids S3 dispatch overhead"));

  SEXP plan_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(plan_, 0, mkString("svd"));
  SET_VECTOR_ELT(plan_, 1, ScalarInteger(rank));
  SET_VECTOR_ELT(plan_, 2, mkString("native certified Gram SVD special case"));
  SET_VECTOR_ELT(plan_, 3, mkString("largest"));
  SET_VECTOR_ELT(plan_, 4, reasons_);
  SET_VECTOR_ELT(plan_, 5, mkString(
    "opt-in implicit candidate retries explicit Gram; native Golub-Kahan if Gram remains uncertified"
  ));
  SET_VECTOR_ELT(plan_, 6, controls_);
  SEXP plan_names_ = PROTECT(allocVector(STRSXP, 7));
  const char* plan_names[] = {
    "problem_type", "requested", "method", "target", "reasons",
    "fallback", "controls"
  };
  for (int i = 0; i < 7; ++i) SET_STRING_ELT(plan_names_, i, mkChar(plan_names[i]));
  setAttrib(plan_, R_NamesSymbol, plan_names_);
  SEXP plan_class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(plan_class_, 0, mkChar("eigencore_plan"));
  setAttrib(plan_, R_ClassSymbol, plan_class_);

  SEXP restart_ = PROTECT(allocVector(VECSXP, 25));
  SET_VECTOR_ELT(restart_, 0, mkString("gram_svd_special_case"));
  SET_VECTOR_ELT(restart_, 1, ScalarLogical(TRUE));
  SET_VECTOR_ELT(restart_, 2, ScalarLogical(TRUE));
  SET_VECTOR_ELT(restart_, 3, mkString(gram_side));
  SET_VECTOR_ELT(restart_, 4, ScalarInteger(m < n ? m : n));
  SET_VECTOR_ELT(restart_, 5, mkString(native_kernel));
  SET_VECTOR_ELT(restart_, 6, eigensolver_);
  SET_VECTOR_ELT(restart_, 7, subspace_be_);
  SET_VECTOR_ELT(restart_, 8, implicit_be_);
  SET_VECTOR_ELT(restart_, 9, implicit_iter_);
  SET_VECTOR_ELT(restart_, 10, gram_krylov_iter_);
  SET_VECTOR_ELT(restart_, 11, ScalarLogical(strcmp(CHAR(STRING_ELT(eigensolver_, 0)), "implicit_normal_lanczos") == 0));
  SET_VECTOR_ELT(restart_, 12, ScalarLogical(strcmp(CHAR(STRING_ELT(eigensolver_, 0)), "implicit_normal_lanczos") != 0));
  SET_VECTOR_ELT(restart_, 13, stage_);
  SET_VECTOR_ELT(restart_, 14, ScalarLogical(FALSE));
  SET_VECTOR_ELT(restart_, 15, ScalarReal(zero_tol));
  SET_VECTOR_ELT(restart_, 16, ScalarLogical(FALSE));
  SET_VECTOR_ELT(restart_, 17, ScalarLogical(TRUE));
  SET_VECTOR_ELT(restart_, 18, ScalarLogical(FALSE));
  SET_VECTOR_ELT(restart_, 19, ScalarLogical(FALSE));
  SET_VECTOR_ELT(restart_, 20, ScalarString(NA_STRING));
  SET_VECTOR_ELT(restart_, 21, ScalarString(NA_STRING));
  SET_VECTOR_ELT(restart_, 22, ScalarLogical(all_converged && orth_passed));
  SET_VECTOR_ELT(restart_, 23, ScalarReal(max_backward));
  SET_VECTOR_ELT(
    restart_, 24,
    ScalarLogical(
      strcmp(CHAR(STRING_ELT(eigensolver_, 0)), "implicit_normal_lanczos") != 0 &&
      INTEGER(implicit_iter_)[0] > 0
    )
  );
  SEXP restart_names_ = PROTECT(allocVector(STRSXP, 25));
  const char* restart_names[] = {
    "kind", "implemented", "native", "gram_side", "gram_dimension",
    "native_gram_kernel", "native_gram_eigensolver",
    "native_gram_subspace_max_backward_error",
    "native_implicit_normal_lanczos_max_backward_error",
    "native_implicit_normal_lanczos_iterations",
    "native_gram_krylov_iterations", "normal_operator_implicit",
    "materialized_gram", "stage_seconds", "zero_singular_completion",
    "zero_singular_threshold", "certificate_reuses_gram_sides",
    "certified_in_original_coordinates", "fallback_attempted",
    "fallback_used", "fallback_method", "fallback_error",
    "gram_certificate_passed", "gram_max_backward_error",
    "explicit_gram_retry_used"
  };
  for (int i = 0; i < 25; ++i) SET_STRING_ELT(restart_names_, i, mkChar(restart_names[i]));
  setAttrib(restart_, R_NamesSymbol, restart_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 19));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, d_);
  SET_VECTOR_ELT(out_, 4, residuals_);
  SET_VECTOR_ELT(out_, 5, backward_);
  SET_VECTOR_ELT(out_, 6, orth_);
  SET_VECTOR_ELT(out_, 7, ScalarInteger(rank));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(rank));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(1));
  SET_VECTOR_ELT(out_, 10, ScalarInteger(1));
  SET_VECTOR_ELT(out_, 11, stage_);
  SET_VECTOR_ELT(out_, 12, mkString("native certified Gram SVD special case"));
  SET_VECTOR_ELT(out_, 13, mkString("largest"));
  SET_VECTOR_ELT(out_, 14, plan_);
  SET_VECTOR_ELT(out_, 15, cert_);
  SET_VECTOR_ELT(out_, 16, restart_);
  SET_VECTOR_ELT(out_, 17, mkString("using native certified Gram SVD special case; residuals certified in original coordinates"));
  SET_VECTOR_ELT(out_, 18, ScalarLogical(TRUE));
  SEXP out_names_ = PROTECT(allocVector(STRSXP, 19));
  const char* out_names[] = {
    "d", "u", "v", "values", "residuals", "backward_error",
    "orthogonality", "nconv", "requested", "iterations", "matvecs",
    "stage_seconds", "method", "target", "plan", "certificate",
    "restart", "warnings", "fastpath_native_result"
  };
  for (int i = 0; i < 19; ++i) SET_STRING_ELT(out_names_, i, mkChar(out_names[i]));
  setAttrib(out_, R_NamesSymbol, out_names_);
  SEXP out_class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(out_class_, 0, mkChar("eigencore_svd_result"));
  setAttrib(out_, R_ClassSymbol, out_class_);

  UNPROTECT(17);
  return out_;
}

extern "C" SEXP eigencore_csc_left_gram_svd_fast_result(SEXP i_, SEXP p_, SEXP x_,
                                                        SEXP dim_, SEXP rank_,
                                                        SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  SEXP native_ = PROTECT(eigencore_csc_left_gram_svd(i_, p_, x_, dim_, rank_, tol_));
  SEXP out_ = eigencore_csc_gram_svd_fast_result_from_native(
    native_,
    INTEGER(dim_)[0],
    INTEGER(dim_)[1],
    asReal(tol_),
    "left",
    "csc_left_gram"
  );
  UNPROTECT(1);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_right_gram_svd_fast_result(SEXP i_, SEXP p_, SEXP x_,
                                                         SEXP dim_, SEXP rank_,
                                                         SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  SEXP native_ = PROTECT(eigencore_csc_right_gram_svd(i_, p_, x_, dim_, rank_, tol_));
  SEXP out_ = eigencore_csc_gram_svd_fast_result_from_native(
    native_,
    INTEGER(dim_)[0],
    INTEGER(dim_)[1],
    asReal(tol_),
    "right",
    "csc_right_gram"
  );
  UNPROTECT(1);
  return out_;
  EIGENCORE_ENTRY_END
}
