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

static Rcomplex complex_conj(Rcomplex z) {
  Rcomplex out;
  out.r = z.r;
  out.i = -z.i;
  return out;
}

// ---------------------------------------------------------------------------
// Shared Krylov helpers
// ---------------------------------------------------------------------------

// Breakdown threshold factor: a new Arnoldi direction is treated as lying in
// the current Krylov space when its norm after orthogonalization falls below
// this multiple of eps times the running ||A v_j|| scale. The test is relative
// so it behaves the same for A, 1e-12 * A and 1e12 * A.
static const double kArnoldiBreakdownFactor = 100.0;

// Euclidean norm: unscaled sum of squares with four accumulators (fast path);
// falls back to the overflow/underflow-safe dnrm2 when the sum leaves the
// comfortable floating-point range, so huge or tiny scalings stay exact.
static inline double arnoldi_norm2(int n, const double* x) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 4 <= n; i += 4) {
    s0 += x[i] * x[i];
    s1 += x[i + 1] * x[i + 1];
    s2 += x[i + 2] * x[i + 2];
    s3 += x[i + 3] * x[i + 3];
  }
  for (; i < n; ++i) {
    s0 += x[i] * x[i];
  }
  const double sum = (s0 + s1) + (s2 + s3);
  if (sum > 1e-280 && sum < 1e280) {
    return std::sqrt(sum);
  }
  if (sum != sum) {
    return sum;
  }
  const int inc = 1;
  return F77_CALL(dnrm2)(&n, x, &inc);
}

static void arnoldi_gemv_t(int n, int cols, const double* __restrict__ V,
                           const double* __restrict__ w, double* __restrict__ out) {
  int c = 0;
  for (; c + 4 <= cols; c += 4) {
    const double* __restrict__ v0 = V + static_cast<int64_t>(c) * n;
    const double* __restrict__ v1 = v0 + n;
    const double* __restrict__ v2 = v1 + n;
    const double* __restrict__ v3 = v2 + n;
    double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
    double b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
    int i = 0;
    for (; i + 2 <= n; i += 2) {
      const double x0 = w[i];
      const double x1 = w[i + 1];
      a0 += v0[i] * x0; b0 += v0[i + 1] * x1;
      a1 += v1[i] * x0; b1 += v1[i + 1] * x1;
      a2 += v2[i] * x0; b2 += v2[i + 1] * x1;
      a3 += v3[i] * x0; b3 += v3[i + 1] * x1;
    }
    for (; i < n; ++i) {
      a0 += v0[i] * w[i];
      a1 += v1[i] * w[i];
      a2 += v2[i] * w[i];
      a3 += v3[i] * w[i];
    }
    out[c] = a0 + b0;
    out[c + 1] = a1 + b1;
    out[c + 2] = a2 + b2;
    out[c + 3] = a3 + b3;
  }
  for (; c < cols; ++c) {
    const double* __restrict__ v0 = V + static_cast<int64_t>(c) * n;
    double a0 = 0.0, b0 = 0.0, c0 = 0.0, d0 = 0.0;
    int i = 0;
    for (; i + 4 <= n; i += 4) {
      a0 += v0[i] * w[i];
      b0 += v0[i + 1] * w[i + 1];
      c0 += v0[i + 2] * w[i + 2];
      d0 += v0[i + 3] * w[i + 3];
    }
    for (; i < n; ++i) {
      a0 += v0[i] * w[i];
    }
    out[c] = (a0 + b0) + (c0 + d0);
  }
}

static void arnoldi_gemv_n_sub(int n, int cols, const double* __restrict__ V,
                               const double* __restrict__ h, double* __restrict__ w) {
  int c = 0;
  for (; c + 4 <= cols; c += 4) {
    const double* __restrict__ v0 = V + static_cast<int64_t>(c) * n;
    const double* __restrict__ v1 = v0 + n;
    const double* __restrict__ v2 = v1 + n;
    const double* __restrict__ v3 = v2 + n;
    const double h0 = h[c], h1 = h[c + 1], h2 = h[c + 2], h3 = h[c + 3];
    for (int i = 0; i < n; ++i) {
      w[i] -= (v0[i] * h0 + v1[i] * h1) + (v2[i] * h2 + v3[i] * h3);
    }
  }
  for (; c < cols; ++c) {
    const double* __restrict__ v0 = V + static_cast<int64_t>(c) * n;
    const double h0 = h[c];
    for (int i = 0; i < n; ++i) {
      w[i] -= v0[i] * h0;
    }
  }
}

// out (n x p) = V (n x m) * Q[, 1:p] (ld ldq). Row-blocked so each block of V
// is read from memory once and reused for all p output columns (a reference
// dgemm streams V once per output column). Used for the Krylov-Schur
// truncation V_p <- V_m Q_p at every restart.
static void arnoldi_combine_blocked(int n, int m, int p, const double* V,
                                    const double* Q, int ldq, double* out) {
  const int block = 512;
  for (int i0 = 0; i0 < n; i0 += block) {
    const int len = std::min(block, n - i0);
    for (int j = 0; j < p; ++j) {
      double* __restrict__ y = out + static_cast<int64_t>(j) * n + i0;
      const double* q = Q + static_cast<int64_t>(j) * ldq;
      std::fill(y, y + len, 0.0);
      int l = 0;
      for (; l + 4 <= m; l += 4) {
        const double* __restrict__ v0 = V + static_cast<int64_t>(l) * n + i0;
        const double* __restrict__ v1 = v0 + n;
        const double* __restrict__ v2 = v1 + n;
        const double* __restrict__ v3 = v2 + n;
        const double q0 = q[l], q1 = q[l + 1], q2 = q[l + 2], q3 = q[l + 3];
        for (int i = 0; i < len; ++i) {
          y[i] += (v0[i] * q0 + v1[i] * q1) + (v2[i] * q2 + v3[i] * q3);
        }
      }
      for (; l < m; ++l) {
        const double* __restrict__ v0 = V + static_cast<int64_t>(l) * n + i0;
        const double q0 = q[l];
        for (int i = 0; i < len; ++i) {
          y[i] += v0[i] * q0;
        }
      }
    }
  }
}

// Classical Gram-Schmidt against the first `cols` columns of V (leading
// dimension n) with the DGKS reorthogonalization test: a second CGS pass runs
// whenever the first pass removed more than 1 - 1/sqrt(2) of ||w|| (always
// when `norm_before` is negative). Each pass is one blocked V^T w product and
// one blocked w -= V h update (see arnoldi_gemv_t / arnoldi_gemv_n_sub), not
// scalar per-column loops. Accumulated projection coefficients are written to
// coeff; tmp is scratch of length cols. On return *norm_after = ||w||.
// Returns the number of passes performed.
static int arnoldi_cgs2(int n, int cols, const double* V, double* w,
                        double* coeff, double* tmp, double norm_before,
                        double* norm_after) {
  if (cols <= 0) {
    *norm_after = (norm_before >= 0.0) ? norm_before : arnoldi_norm2(n, w);
    return 0;
  }
  arnoldi_gemv_t(n, cols, V, w, coeff);
  arnoldi_gemv_n_sub(n, cols, V, coeff, w);
  double after = arnoldi_norm2(n, w);
  if (norm_before >= 0.0 && after > 0.7071067811865476 * norm_before) {
    *norm_after = after;
    return 1;
  }
  arnoldi_gemv_t(n, cols, V, w, tmp);
  arnoldi_gemv_n_sub(n, cols, V, tmp, w);
  for (int i = 0; i < cols; ++i) {
    coeff[i] += tmp[i];
  }
  *norm_after = arnoldi_norm2(n, w);
  return 2;
}

static double arnoldi_normalize_start(int n, const double* start, double* v) {
  const double start_norm = arnoldi_norm2(n, start);
  if (!std::isfinite(start_norm) || !(start_norm > 0.0)) {
    error("native Arnoldi start vector is zero or non-finite");
  }
  for (int i = 0; i < n; ++i) {
    v[i] = start[i] / start_norm;
  }
  return start_norm;
}

// Fill v with a random direction orthonormal to the first `cols` columns of V.
// Used after an Arnoldi breakdown (an invariant subspace was found) so the
// factorization A V = V H + f e^T stays valid with a zero subdiagonal entry.
// Returns false when no such direction exists numerically (cols == n).
static bool arnoldi_random_orthonormal(int n, int cols, const double* V,
                                       double* v, double* coeff, double* tmp) {
  if (cols >= n) {
    return false;
  }
  for (int attempt = 0; attempt < 3; ++attempt) {
    GetRNGstate();
    for (int i = 0; i < n; ++i) {
      v[i] = norm_rand();
    }
    PutRNGstate();
    const double before = arnoldi_norm2(n, v);
    double after = 0.0;
    arnoldi_cgs2(n, cols, V, v, coeff, tmp, -1.0, &after);
    if (std::isfinite(after) && after > 1e-3 * before) {
      for (int i = 0; i < n; ++i) {
        v[i] /= after;
      }
      return true;
    }
  }
  return false;
}

static void arnoldi_apply_column(void* impl, EigencoreApplyFn apply, int n,
                                 const double* x, double* y,
                                 EigencoreWorkspace* workspace,
                                 const char* context) {
  std::fill(y, y + n, 0.0);
  const int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, x, n,
                           1.0, 0.0, y, n, workspace);
  if (status != 0) {
    eigencore_apply_status_error(context, status);
  }
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) {
      error("%s returned a non-finite value (check input for NA/NaN/Inf)", context);
    }
  }
}

static SEXP native_arnoldi_cycle_impl(void* impl,
                                      EigencoreApplyFn apply,
                                      int64_t n64,
                                      const double* start,
                                      int max_subspace) {
  if (n64 < 1 || max_subspace < 1) {
    error("native Arnoldi requires positive dimensions");
  }
  if (!eigencore_int_indexable(n64) || !eigencore_int_indexable(max_subspace + 1) ||
      !eigencore_int_indexable(n64 * (static_cast<int64_t>(max_subspace) + 1))) {
    error("native Arnoldi dimensions exceed LP64 BLAS/R integer range");
  }
  const int n = static_cast<int>(n64);
  const int m_budget = std::min(max_subspace, n);
  const int ldh = m_budget + 1;
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, m_budget + 1));
  SEXP H_ = PROTECT(allocMatrix(REALSXP, ldh, m_budget));
  double* V = REAL(V_);
  double* H = REAL(H_);
  std::memset(V, 0, sizeof(double) * static_cast<size_t>(n) *
              static_cast<size_t>(m_budget + 1));
  std::memset(H, 0, sizeof(double) * static_cast<size_t>(ldh) *
              static_cast<size_t>(m_budget));

  arnoldi_normalize_start(n, start, V);

  double* w = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(n), sizeof(double)));
  double* tmp = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(ldh), sizeof(double)));
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int iterations = 0;
  int matvecs = 0;
  int reorthogonalization_passes = 0;
  double scale = 0.0;

  for (int j = 0; j < m_budget; ++j) {
    arnoldi_apply_column(impl, apply, n, V + static_cast<int64_t>(j) * n, w,
                         &workspace, "native Arnoldi operator apply");
    ++matvecs;
    const double wnorm = arnoldi_norm2(n, w);
    scale = std::max(scale, wnorm);
    double* hcol = H + static_cast<int64_t>(j) * ldh;
    double beta = 0.0;
    reorthogonalization_passes += arnoldi_cgs2(n, j + 1, V, w, hcol, tmp, wnorm, &beta);

    iterations = j + 1;
    if (!std::isfinite(beta)) {
      error("native Arnoldi produced a non-finite residual norm");
    }
    if (beta <= kArnoldiBreakdownFactor * DBL_EPSILON * scale ||
        iterations == m_budget) {
      H[(j + 1) + static_cast<int64_t>(j) * ldh] =
        (beta <= kArnoldiBreakdownFactor * DBL_EPSILON * scale) ? 0.0 : beta;
      break;
    }
    H[(j + 1) + static_cast<int64_t>(j) * ldh] = beta;
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

// ---------------------------------------------------------------------------
// Real Krylov-Schur restarted Arnoldi (Stewart 2001)
// ---------------------------------------------------------------------------
//
// Maintains A V_l = V_l S_l + v_{l+1} b^T with V orthonormal and S_l upper
// quasi-triangular after each restart. Each outer iteration:
//   1. expands the factorization to size m with CGS2 Arnoldi steps;
//   2. computes the real Schur form S_m = Q T Q^T (dgees);
//   3. moves the wanted Ritz values to the leading block with dtrsen, never
//      splitting a 2x2 block (conjugate pairs stay together);
//   4. estimates Ritz residuals |b^T Q y| from eigenvectors y of the leading
//      block (dtrevc) and stops once the k wanted values satisfy
//      res <= tol * max(|theta|, eps^(1/3) * ||T||_F), a scale-invariant
//      test (|theta|-relative like ARPACK, with a floor for theta ~ 0);
//   5. otherwise truncates to the leading p columns and continues.
// The returned decomposition (V: n x (p+1), H = [T_p; b^T]: (p+1) x p) feeds
// the existing projected/refined Ritz extraction; certificates are computed
// afterwards in R against the original operator.

enum ArnoldiTargetCode {
  ARNOLDI_TARGET_LARGEST_REAL = 0,
  ARNOLDI_TARGET_SMALLEST_REAL = 1,
  ARNOLDI_TARGET_LARGEST_MAGNITUDE = 2,
  ARNOLDI_TARGET_LARGEST_IMAGINARY = 3,
  ARNOLDI_TARGET_SMALLEST_IMAGINARY = 4,
  ARNOLDI_TARGET_SMALLEST_MAGNITUDE = 5,
  // Internal: nearest to a given set of values. Used by the adjoint
  // (left-eigenvector) solve, which must find the eigenvalues the right solve
  // already returned rather than re-deciding which ones are extremal.
  ARNOLDI_TARGET_NEAREST_SET = 6
};

struct ArnoldiTarget {
  int code;
  const Rcomplex* values;
  int count;
};

static inline double arnoldi_target_key(const ArnoldiTarget& target, double re, double im) {
  switch (target.code) {
  case ARNOLDI_TARGET_NEAREST_SET: {
    double best = R_PosInf;
    for (int j = 0; j < target.count; ++j) {
      best = std::min(best, std::hypot(re - target.values[j].r, im - target.values[j].i));
    }
    return -best;
  }
  case ARNOLDI_TARGET_SMALLEST_REAL:
    return -re;
  case ARNOLDI_TARGET_LARGEST_MAGNITUDE:
    return std::hypot(re, im);
  case ARNOLDI_TARGET_LARGEST_IMAGINARY:
    return im;
  case ARNOLDI_TARGET_SMALLEST_IMAGINARY:
    return -im;
  case ARNOLDI_TARGET_SMALLEST_MAGNITUDE:
    return -std::hypot(re, im);
  case ARNOLDI_TARGET_LARGEST_REAL:
  default:
    return re;
  }
}

// order[0..count) = indices sorted from most to least wanted (stable).
static void arnoldi_rank(const ArnoldiTarget& target, int count, const double* wr,
                         const double* wi, int* order) {
  for (int i = 0; i < count; ++i) {
    order[i] = i;
  }
  // Insertion sort: count is the (small) Krylov dimension and this avoids
  // heap allocation inside .Call code that can longjmp.
  for (int i = 1; i < count; ++i) {
    const int idx = order[i];
    const double key = arnoldi_target_key(target, wr[idx], wi[idx]);
    int j = i - 1;
    while (j >= 0 && arnoldi_target_key(target, wr[order[j]], wi[order[j]]) < key) {
      order[j + 1] = order[j];
      --j;
    }
    order[j + 1] = idx;
  }
}

// Partner of eigenvalue i inside a real Schur form (pairs are adjacent with
// wi[j] > 0, wi[j + 1] < 0), or -1 for a real eigenvalue.
static inline int arnoldi_schur_partner(int i, int count, const double* wi) {
  if (wi[i] > 0.0 && i + 1 < count) {
    return i + 1;
  }
  if (wi[i] < 0.0 && i > 0) {
    return i - 1;
  }
  return -1;
}

// Select the `want` most wanted eigenvalues (whole conjugate pairs) without
// exceeding max_keep. Returns the number of selected eigenvalues.
static int arnoldi_select(int count, const double* wi, const int* order,
                          int want, int max_keep, int* select) {
  std::fill(select, select + count, 0);
  int chosen = 0;
  // `taken` counts ranked picks; a conjugate partner added to keep a 2x2
  // block intact does not use up the budget (for imaginary-part targets the
  // partner of a wanted value is the least wanted one), but a partner that is
  // itself reached in rank order does count.
  int taken = 0;
  for (int r = 0; r < count && taken < want; ++r) {
    const int i = order[r];
    if (select[i]) {
      ++taken;
      continue;
    }
    const int partner = arnoldi_schur_partner(i, count, wi);
    const int size = (partner >= 0) ? 2 : 1;
    if (chosen + size > max_keep) {
      continue;
    }
    select[i] = 1;
    if (partner >= 0) {
      select[partner] = 1;
    }
    chosen += size;
    ++taken;
  }
  return chosen;
}

struct ArnoldiKrylovSchurResult {
  int size;
  int matvecs;
  int restarts;
  int nconv;
  int reorthogonalization_passes;
  int breakdowns;
  double apply_seconds;
  double orthogonalization_seconds;
  double schur_seconds;
  double restart_seconds;
};

static SEXP native_krylov_schur_impl(void* impl,
                                     EigencoreApplyFn apply,
                                     int64_t n64,
                                     const double* start,
                                     int k,
                                     int max_subspace,
                                     const ArnoldiTarget& target,
                                     double tol,
                                     int max_iterations) {
  if (n64 < 1 || max_subspace < 1 || k < 1) {
    error("native Krylov-Schur Arnoldi requires positive dimensions");
  }
  if (!eigencore_int_indexable(n64) ||
      !eigencore_int_indexable(n64 * (static_cast<int64_t>(max_subspace) + 1))) {
    error("native Krylov-Schur Arnoldi dimensions exceed LP64 BLAS/R integer range");
  }
  if (!(tol > 0.0) || !std::isfinite(tol)) {
    error("native Krylov-Schur Arnoldi requires a positive finite tolerance");
  }
  if (target.code < 0 || target.code > 6 ||
      (target.code == ARNOLDI_TARGET_NEAREST_SET && target.count < 1)) {
    error("native Krylov-Schur Arnoldi target code is out of range");
  }
  const int n = static_cast<int>(n64);
  const int m = std::min(max_subspace, n);
  if (k > m) {
    error("native Krylov-Schur Arnoldi requires k <= max_subspace");
  }
  max_iterations = std::max(1, max_iterations);
  const int lds = m + 1;
  const size_t n_sz = static_cast<size_t>(n);
  const size_t m_sz = static_cast<size_t>(m);

  // All workspace through R_alloc so an R error (including one raised by a
  // matrix-free callback) cannot leak it.
  double* V = reinterpret_cast<double*>(R_alloc(n_sz * (m_sz + 1), sizeof(double)));
  double* S = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(lds) * m_sz, sizeof(double)));
  double* T = reinterpret_cast<double*>(R_alloc(m_sz * m_sz, sizeof(double)));
  double* Q = reinterpret_cast<double*>(R_alloc(m_sz * m_sz, sizeof(double)));
  double* Y = reinterpret_cast<double*>(R_alloc(m_sz * m_sz, sizeof(double)));
  double* Wbuf = reinterpret_cast<double*>(R_alloc(n_sz * m_sz, sizeof(double)));
  double* w = reinterpret_cast<double*>(R_alloc(n_sz, sizeof(double)));
  double* coeff = reinterpret_cast<double*>(R_alloc(m_sz + 1, sizeof(double)));
  double* tmp = reinterpret_cast<double*>(R_alloc(m_sz + 1, sizeof(double)));
  double* wr = reinterpret_cast<double*>(R_alloc(m_sz, sizeof(double)));
  double* wi = reinterpret_cast<double*>(R_alloc(m_sz, sizeof(double)));
  double* b = reinterpret_cast<double*>(R_alloc(m_sz, sizeof(double)));
  double* bq = reinterpret_cast<double*>(R_alloc(m_sz, sizeof(double)));
  double* res = reinterpret_cast<double*>(R_alloc(m_sz, sizeof(double)));
  double* trevc_work = reinterpret_cast<double*>(R_alloc(3 * m_sz + 1, sizeof(double)));
  int* order = reinterpret_cast<int*>(R_alloc(m_sz, sizeof(int)));
  int* select = reinterpret_cast<int*>(R_alloc(m_sz, sizeof(int)));
  int* bwork = reinterpret_cast<int*>(R_alloc(m_sz, sizeof(int)));
  std::memset(V, 0, sizeof(double) * n_sz * (m_sz + 1));
  std::memset(S, 0, sizeof(double) * static_cast<size_t>(lds) * m_sz);

  // dgees workspace query (m is fixed for the whole run).
  const char jobvs = 'V';
  const char sort = 'N';
  int sdim = 0;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgees)(&jobvs, &sort, nullptr, &m, T, &m, &sdim, wr, wi, Q, &m,
                  &work_query, &lwork, bwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgees workspace query failed for native Krylov-Schur Arnoldi with info=%d", info);
  }
  lwork = std::max(std::max(1, 3 * m), static_cast<int>(work_query));
  double* gees_work = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(lwork), sizeof(double)));

  arnoldi_normalize_start(n, start, V);

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  ArnoldiKrylovSchurResult state = {0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0};
  const int max_keep = (m < n) ? m - 1 : m;
  int l = 0;
  double scale = 0.0;
  int p = 0;
  bool converged = false;

  for (int iter = 0; iter < max_iterations; ++iter) {
    R_CheckUserInterrupt();
    // 1. Expand A V_j = V_{j+1} S_j from size l to size m.
    for (int j = l; j < m; ++j) {
      auto stage_start = native_timer_now();
      arnoldi_apply_column(impl, apply, n, V + static_cast<int64_t>(j) * n, w,
                           &workspace, "native Krylov-Schur Arnoldi operator apply");
      state.apply_seconds += native_timer_elapsed(stage_start);
      stage_start = native_timer_now();
      ++state.matvecs;
      const double wnorm = arnoldi_norm2(n, w);
      scale = std::max(scale, wnorm);
      double* scol = S + static_cast<int64_t>(j) * lds;
      double beta = 0.0;
      state.reorthogonalization_passes +=
        arnoldi_cgs2(n, j + 1, V, w, scol, tmp, wnorm, &beta);
      if (!std::isfinite(beta)) {
        error("native Krylov-Schur Arnoldi produced a non-finite residual norm");
      }
      double* vnext = V + static_cast<int64_t>(j + 1) * n;
      const bool full_space = (j + 1 == n);
      if (full_space || beta <= kArnoldiBreakdownFactor * DBL_EPSILON * scale) {
        // Invariant subspace: the factorization is exact in this column.
        scol[j + 1] = 0.0;
        if (j + 1 < m) {
          ++state.breakdowns;
          if (!arnoldi_random_orthonormal(n, j + 1, V, vnext, coeff, tmp)) {
            error("native Krylov-Schur Arnoldi could not extend an exhausted Krylov space");
          }
        } else {
          std::fill(vnext, vnext + n, 0.0);
        }
      } else {
        scol[j + 1] = beta;
        for (int row = 0; row < n; ++row) {
          vnext[row] = w[row] / beta;
        }
      }
      state.orthogonalization_seconds += native_timer_elapsed(stage_start);
    }

    // 2. Real Schur form of the m x m projected matrix.
    auto schur_start = native_timer_now();
    double tnorm2 = 0.0;
    for (int col = 0; col < m; ++col) {
      for (int row = 0; row < m; ++row) {
        const double value = S[row + static_cast<int64_t>(col) * lds];
        T[row + static_cast<int64_t>(col) * m] = value;
        tnorm2 += value * value;
      }
      b[col] = S[m + static_cast<int64_t>(col) * lds];
    }
    const double tnorm = std::sqrt(tnorm2);
    F77_CALL(dgees)(&jobvs, &sort, nullptr, &m, T, &m, &sdim, wr, wi, Q, &m,
                    gees_work, &lwork, bwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgees failed for native Krylov-Schur Arnoldi with info=%d", info);
    }

    // 3. Choose how many Ritz values to keep (ARPACK/Spectra heuristic) and
    //    reorder them to the leading block.
    int want = k + std::min(state.nconv, (m - k) / 2);
    if (k == 1 && m >= 6) {
      want = std::max(want, m / 2);
    } else if (k == 1 && m > 2) {
      want = std::max(want, 2);
    }
    want = std::min(want, max_keep);
    arnoldi_rank(target, m, wr, wi, order);
    int chosen = arnoldi_select(m, wi, order, want, max_keep, select);
    if (chosen < 1) {
      chosen = arnoldi_select(m, wi, order, 1, m, select);
    }
    const char job = 'N';
    const char compq = 'V';
    int msel = 0;
    double s_cond = 0.0;
    double sep = 0.0;
    int liwork = 1;
    int iwork = 0;
    F77_CALL(dtrsen)(&job, &compq, select, &m, T, &m, Q, &m, wr, wi, &msel,
                     &s_cond, &sep, gees_work, &lwork, &iwork, &liwork,
                     &info FCONE FCONE);
    if (info < 0) {
      error("LAPACK dtrsen failed for native Krylov-Schur Arnoldi with info=%d", info);
    }
    // info == 1: reordering stopped early (ill-conditioned swap); T and Q are
    // still a valid Schur decomposition, so keep the leading block as is.
    p = std::max(1, std::min(msel, max_keep));
    if (p < m && T[p + static_cast<int64_t>(p - 1) * m] != 0.0) {
      p = (p + 1 <= max_keep) ? p + 1 : p - 1;
    }
    if (p < 1) {
      p = (m >= 2 && T[1] != 0.0) ? 2 : 1;
    }

    // 4. Residual estimates for the leading p Ritz pairs.
    const char trans = 'T';
    const double one = 1.0;
    const double zero = 0.0;
    const int inc = 1;
    F77_CALL(dgemv)(&trans, &m, &m, &one, Q, &m, b, &inc, &zero, bq, &inc FCONE);
    const char side = 'R';
    const char howmny = 'A';
    int vl_dummy_ld = 1;
    double vl_dummy = 0.0;
    int mout = 0;
    F77_CALL(dtrevc)(&side, &howmny, select, &p, T, &m, &vl_dummy, &vl_dummy_ld,
                     Y, &p, &p, &mout, trevc_work, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dtrevc failed for native Krylov-Schur Arnoldi with info=%d", info);
    }
    for (int j = 0; j < p;) {
      const double* yr = Y + static_cast<int64_t>(j) * p;
      const bool pair = (j + 1 < p) && T[(j + 1) + static_cast<int64_t>(j) * m] != 0.0;
      if (!pair) {
        const double ynorm = arnoldi_norm2(p, yr);
        const double dot = F77_CALL(ddot)(&p, bq, &inc, yr, &inc);
        res[j] = (ynorm > 0.0) ? std::fabs(dot) / ynorm : std::fabs(dot);
        j += 1;
      } else {
        const double* yi = Y + static_cast<int64_t>(j + 1) * p;
        const double nr = arnoldi_norm2(p, yr);
        const double ni = arnoldi_norm2(p, yi);
        const double ynorm = std::sqrt(nr * nr + ni * ni);
        const double d1 = F77_CALL(ddot)(&p, bq, &inc, yr, &inc);
        const double d2 = F77_CALL(ddot)(&p, bq, &inc, yi, &inc);
        const double r = (ynorm > 0.0) ? std::hypot(d1, d2) / ynorm : std::hypot(d1, d2);
        res[j] = r;
        res[j + 1] = r;
        j += 2;
      }
    }
    // Scale-invariant convergence test (ARPACK-style, relative to |theta|,
    // with a floor proportional to ||T||_F so eigenvalues at or near zero
    // can converge): res <= tol_eff * max(|theta|, eps^(1/3) * ||T||_F).
    // tol_eff is never below 100 * eps, the attainable residual level.
    arnoldi_rank(target, p, wr, wi, order);
    int nconv = 0;
    const int check = std::min(k, p);
    const double tol_eff = std::max(tol, 100.0 * DBL_EPSILON);
    const double floor_scale = std::cbrt(DBL_EPSILON) * tnorm;
    for (int r = 0; r < check; ++r) {
      const int idx = order[r];
      const double theta = std::hypot(wr[idx], wi[idx]);
      if (res[idx] <= tol_eff * std::max(theta, floor_scale)) {
        ++nconv;
      }
    }
    state.nconv = nconv;
    converged = (nconv >= k) || (m == n);

    state.schur_seconds += native_timer_elapsed(schur_start);
    auto restart_start = native_timer_now();
    // 5. Truncate to the leading p Schur vectors: V_p <- V_m Q_p, keep the
    //    residual vector, S_p <- [T_p; bq_p^T].
    arnoldi_combine_blocked(n, m, p, V, Q, m, Wbuf);
    std::memcpy(V, Wbuf, sizeof(double) * n_sz * static_cast<size_t>(p));
    std::memcpy(V + static_cast<int64_t>(p) * n, V + static_cast<int64_t>(m) * n,
                sizeof(double) * n_sz);
    std::memset(S, 0, sizeof(double) * static_cast<size_t>(lds) * m_sz);
    for (int col = 0; col < p; ++col) {
      for (int row = 0; row < p; ++row) {
        S[row + static_cast<int64_t>(col) * lds] = T[row + static_cast<int64_t>(col) * m];
      }
      S[p + static_cast<int64_t>(col) * lds] = bq[col];
    }
    l = p;
    state.restart_seconds += native_timer_elapsed(restart_start);
    if (converged) {
      break;
    }
    ++state.restarts;
  }
  state.size = p;

  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, p + 1));
  SEXP H_ = PROTECT(allocMatrix(REALSXP, p + 1, p));
  std::memcpy(REAL(V_), V, sizeof(double) * n_sz * static_cast<size_t>(p + 1));
  double* H = REAL(H_);
  for (int col = 0; col < p; ++col) {
    for (int row = 0; row <= p; ++row) {
      H[row + static_cast<int64_t>(col) * (p + 1)] = S[row + static_cast<int64_t>(col) * lds];
    }
  }
  SEXP res_ = PROTECT(allocVector(REALSXP, p));
  std::memcpy(REAL(res_), res, sizeof(double) * static_cast<size_t>(p));

  SEXP stage_ = PROTECT(allocVector(REALSXP, 4));
  REAL(stage_)[0] = state.apply_seconds;
  REAL(stage_)[1] = state.orthogonalization_seconds;
  REAL(stage_)[2] = state.schur_seconds;
  REAL(stage_)[3] = state.restart_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(stage_names_, 0, mkChar("apply"));
  SET_STRING_ELT(stage_names_, 1, mkChar("orthogonalization"));
  SET_STRING_ELT(stage_names_, 2, mkChar("projected_schur"));
  SET_STRING_ELT(stage_names_, 3, mkChar("restart"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  const int n_out = 13;
  SEXP out_ = PROTECT(allocVector(VECSXP, n_out));
  SEXP names_ = PROTECT(allocVector(STRSXP, n_out));
  int slot = 0;
  auto put = [&](const char* name, SEXP value) {
    SET_VECTOR_ELT(out_, slot, value);
    SET_STRING_ELT(names_, slot, mkChar(name));
    ++slot;
  };
  put("V", V_);
  put("H", H_);
  put("iterations", ScalarInteger(p));
  put("matvecs", ScalarInteger(state.matvecs));
  put("reorthogonalization_passes", ScalarInteger(state.reorthogonalization_passes));
  put("native_workspace_bytes", ScalarInteger(static_cast<int>(workspace.bytes_allocated)));
  put("krylov_schur_restarts", ScalarInteger(state.restarts));
  put("nconv", ScalarInteger(state.nconv));
  put("converged", ScalarLogical(converged ? TRUE : FALSE));
  put("residual_estimates", res_);
  put("max_subspace", ScalarInteger(m));
  put("breakdowns", ScalarInteger(state.breakdowns));
  put("stage_seconds", stage_);
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
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

static ArnoldiTarget arnoldi_target_or_error(SEXP target_, SEXP target_values_) {
  const int code = asInteger(target_);
  if (code == NA_INTEGER || code < 0 || code > 5) {
    error("invalid native Krylov-Schur Arnoldi target code");
  }
  ArnoldiTarget target = {code, nullptr, 0};
  if (target_values_ != R_NilValue) {
    if (!isComplex(target_values_) || LENGTH(target_values_) < 1) {
      error("native Krylov-Schur Arnoldi target values must be a non-empty complex vector");
    }
    target.code = ARNOLDI_TARGET_NEAREST_SET;
    target.values = COMPLEX(target_values_);
    target.count = LENGTH(target_values_);
  }
  return target;
}

extern "C" SEXP eigencore_arnoldi_ks_dense(SEXP A_, SEXP start_, SEXP k_,
                                            SEXP max_subspace_, SEXP target_,
                                            SEXP target_values_,
                                            SEXP tol_, SEXP maxit_) {
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
  return native_krylov_schur_impl(
    &impl, eigencore_dense_apply, n, REAL(start_), asInteger(k_),
    asInteger(max_subspace_), arnoldi_target_or_error(target_, target_values_),
    asReal(tol_), asInteger(maxit_)
  );
}

// Applies the adjoint of a CSC matrix: when the CSC storage holds B = A^T,
// y = B^T x = A x is a row-gather (one dot product per output entry) rather
// than the scatter y[row] += ... of the forward CSC product, which is
// noticeably faster for a single vector.
static int arnoldi_csc_adjoint_apply(void* impl, EigencoreTranspose op,
                                     int64_t block_cols, const double* X,
                                     int64_t ldx, double alpha, double beta,
                                     double* Y, int64_t ldy,
                                     EigencoreWorkspace* workspace) {
  const EigencoreTranspose flipped = (op == EIGENCORE_TRANSPOSE_NONE)
    ? EIGENCORE_TRANSPOSE_ADJOINT
    : EIGENCORE_TRANSPOSE_NONE;
  return eigencore_csc_apply(impl, flipped, block_cols, X, ldx, alpha, beta,
                             Y, ldy, workspace);
}

// CSC Krylov-Schur. With transposed = TRUE the CSC slots hold the transpose
// B = A^T of the operator and A x is applied as B^T x (row gather).
extern "C" SEXP eigencore_arnoldi_ks_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                          SEXP transposed_,
                                          SEXP start_, SEXP k_,
                                          SEXP max_subspace_, SEXP target_,
                                          SEXP target_values_,
                                          SEXP tol_, SEXP maxit_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      LENGTH(dim_) != 2 || !isReal(start_)) {
    error("invalid CSC Krylov-Schur Arnoldi inputs");
  }
  const int nrow = INTEGER(dim_)[0];
  const int ncol = INTEGER(dim_)[1];
  if (nrow != ncol) {
    error("A must be a square dgCMatrix");
  }
  if (LENGTH(start_) != nrow || LENGTH(p_) != ncol + 1 ||
      LENGTH(i_) != LENGTH(x_)) {
    error("start length must equal matrix dimension");
  }
  const bool transposed = asLogical(transposed_) == TRUE;
  CSCOperator impl = {nrow, ncol, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_krylov_schur_impl(
    &impl, transposed ? arnoldi_csc_adjoint_apply : eigencore_csc_apply,
    nrow, REAL(start_), asInteger(k_),
    asInteger(max_subspace_), arnoldi_target_or_error(target_, target_values_),
    asReal(tol_), asInteger(maxit_)
  );
}

extern "C" SEXP eigencore_arnoldi_ks_r_operator(SEXP dim_, SEXP apply_,
                                                 SEXP start_, SEXP k_,
                                                 SEXP max_subspace_,
                                                 SEXP target_,
                                                 SEXP target_values_,
                                                 SEXP tol_, SEXP maxit_) {
  if (!isInteger(dim_) || LENGTH(dim_) != 2 || TYPEOF(apply_) != CLOSXP ||
      !isReal(start_)) {
    error("invalid matrix-free Krylov-Schur Arnoldi inputs");
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
  return native_krylov_schur_impl(
    &impl, eigencore_r_operator_apply, nrow, REAL(start_), asInteger(k_),
    asInteger(max_subspace_), arnoldi_target_or_error(target_, target_values_),
    asReal(tol_), asInteger(maxit_)
  );
}

// Eigen-decomposition of the leading m x m block of H (leading dimension
// h_rows). Writes eigenvalues and complex Ritz coefficient vectors (unit
// Euclidean norm, dgeev convention) split into real/imaginary parts.
static void arnoldi_projected_eigen(const double* H, int h_rows, int m,
                                    double* wr, double* wi,
                                    double* coeff_re, double* coeff_im) {
  double* Hm = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * m, sizeof(double)));
  double* vr = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * m, sizeof(double)));
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < m; ++row) {
      Hm[row + static_cast<int64_t>(col) * m] = H[row + static_cast<int64_t>(col) * h_rows];
    }
  }
  double vl_dummy = 0.0;
  int ldvl = 1;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  const char jobvl = 'N';
  const char jobvr = 'V';
  F77_CALL(dgeev)(&jobvl, &jobvr, &m, Hm, &m, wr, wi, &vl_dummy, &ldvl, vr, &m,
                  &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgeev workspace query failed for native Arnoldi Ritz extraction with info=%d", info);
  }
  lwork = std::max(std::max(1, 4 * m), static_cast<int>(work_query));
  double* work = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(lwork), sizeof(double)));
  F77_CALL(dgeev)(&jobvl, &jobvr, &m, Hm, &m, wr, wi, &vl_dummy, &ldvl, vr, &m,
                  work, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgeev failed for native Arnoldi Ritz extraction with info=%d", info);
  }
  for (int col = 0; col < m; ++col) {
    double* cr = coeff_re + static_cast<int64_t>(col) * m;
    double* ci = coeff_im + static_cast<int64_t>(col) * m;
    if (wi[col] > 0.0 && col + 1 < m) {
      for (int row = 0; row < m; ++row) {
        cr[row] = vr[row + static_cast<int64_t>(col) * m];
        ci[row] = vr[row + static_cast<int64_t>(col + 1) * m];
      }
    } else if (wi[col] < 0.0 && col > 0) {
      for (int row = 0; row < m; ++row) {
        cr[row] = vr[row + static_cast<int64_t>(col - 1) * m];
        ci[row] = -vr[row + static_cast<int64_t>(col) * m];
      }
    } else {
      for (int row = 0; row < m; ++row) {
        cr[row] = vr[row + static_cast<int64_t>(col) * m];
        ci[row] = 0.0;
      }
    }
  }
}

// Form Z = V[, 1:m] %*% (C_re + i C_im) for `cols` coefficient columns with a
// single dgemm on the stacked [C_re | C_im] block, then normalize each column.
static SEXP arnoldi_form_ritz_vectors(const double* V, int n, int m,
                                      const double* coeff_re, int ld_re,
                                      const double* coeff_im, int ld_im,
                                      int cols) {
  SEXP vectors_ = PROTECT(allocMatrix(CPLXSXP, n, cols));
  if (cols == 0) {
    UNPROTECT(1);
    return vectors_;
  }
  const int two_cols = 2 * cols;
  double* B = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * two_cols, sizeof(double)));
  double* Z = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(n) * two_cols, sizeof(double)));
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < m; ++row) {
      B[row + static_cast<int64_t>(col) * m] = coeff_re[row + static_cast<int64_t>(col) * ld_re];
      B[row + static_cast<int64_t>(cols + col) * m] = coeff_im[row + static_cast<int64_t>(col) * ld_im];
    }
  }
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&notrans, &notrans, &n, &two_cols, &m, &one, V, &n, B, &m,
                  &zero, Z, &n FCONE FCONE);
  Rcomplex* vectors = COMPLEX(vectors_);
  for (int col = 0; col < cols; ++col) {
    const double* zr = Z + static_cast<int64_t>(col) * n;
    const double* zi = Z + static_cast<int64_t>(cols + col) * n;
    const double nr = arnoldi_norm2(n, zr);
    const double ni = arnoldi_norm2(n, zi);
    const double norm = std::sqrt(nr * nr + ni * ni);
    const double inv = (std::isfinite(norm) && norm > 0.0) ? 1.0 / norm : 1.0;
    Rcomplex* out = vectors + static_cast<int64_t>(col) * n;
    for (int row = 0; row < n; ++row) {
      out[row].r = zr[row] * inv;
      out[row].i = zi[row] * inv;
    }
  }
  UNPROTECT(1);
  return vectors_;
}

static void arnoldi_check_vh(SEXP V_, SEXP H_, int m, int* n_out, int* h_rows_out,
                             const char* context) {
  if (!isReal(V_) || !isReal(H_)) {
    error("%s requires real V and H matrices", context);
  }
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  SEXP dimH = getAttrib(H_, R_DimSymbol);
  if (dimV == R_NilValue || dimH == R_NilValue ||
      LENGTH(dimV) != 2 || LENGTH(dimH) != 2) {
    error("%s requires matrix inputs", context);
  }
  const int n = INTEGER(dimV)[0];
  const int v_cols = INTEGER(dimV)[1];
  const int h_rows = INTEGER(dimH)[0];
  const int h_cols = INTEGER(dimH)[1];
  if (n < 1 || m < 1 || v_cols < m || h_rows < m || h_cols < m) {
    error("invalid %s dimensions", context);
  }
  if (!eigencore_int_indexable(static_cast<int64_t>(n) * m * 2)) {
    error("%s dimensions exceed LP64 BLAS/R integer range", context);
  }
  *n_out = n;
  *h_rows_out = h_rows;
}

// Projected Ritz values and complex Ritz coefficient vectors (m x m, unit
// norm columns) of the leading m x m block of H. Cheap: no n-length work.
extern "C" SEXP eigencore_arnoldi_ritz_coefficients(SEXP H_, SEXP iterations_) {
  if (!isReal(H_)) {
    error("native Arnoldi Ritz coefficients require a real H matrix");
  }
  SEXP dimH = getAttrib(H_, R_DimSymbol);
  const int m = asInteger(iterations_);
  if (dimH == R_NilValue || LENGTH(dimH) != 2 || m == NA_INTEGER || m < 1 ||
      INTEGER(dimH)[0] < m || INTEGER(dimH)[1] < m) {
    error("invalid native Arnoldi Ritz coefficient dimensions");
  }
  const int h_rows = INTEGER(dimH)[0];
  double* wr = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m), sizeof(double)));
  double* wi = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m), sizeof(double)));
  double* cre = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * m, sizeof(double)));
  double* cim = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * m, sizeof(double)));
  arnoldi_projected_eigen(REAL(H_), h_rows, m, wr, wi, cre, cim);
  SEXP values_ = PROTECT(allocVector(CPLXSXP, m));
  SEXP coeff_ = PROTECT(allocMatrix(CPLXSXP, m, m));
  Rcomplex* values = COMPLEX(values_);
  Rcomplex* coeff = COMPLEX(coeff_);
  for (int col = 0; col < m; ++col) {
    values[col].r = wr[col];
    values[col].i = wi[col];
    for (int row = 0; row < m; ++row) {
      coeff[row + static_cast<int64_t>(col) * m].r = cre[row + static_cast<int64_t>(col) * m];
      coeff[row + static_cast<int64_t>(col) * m].i = cim[row + static_cast<int64_t>(col) * m];
    }
  }
  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, coeff_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("coefficients"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}

// Ritz vectors V[, 1:m] %*% C for a complex m x c coefficient block, formed
// with one real dgemm and normalized to unit 2-norm.
extern "C" SEXP eigencore_arnoldi_ritz_vectors(SEXP V_, SEXP iterations_,
                                                SEXP coefficients_) {
  if (!isReal(V_) || !isComplex(coefficients_)) {
    error("native Arnoldi Ritz vectors require real V and complex coefficients");
  }
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  SEXP dimC = getAttrib(coefficients_, R_DimSymbol);
  const int m = asInteger(iterations_);
  if (dimV == R_NilValue || dimC == R_NilValue || LENGTH(dimV) != 2 ||
      LENGTH(dimC) != 2 || m == NA_INTEGER || m < 1 ||
      INTEGER(dimV)[1] < m || INTEGER(dimC)[0] != m) {
    error("invalid native Arnoldi Ritz vector dimensions");
  }
  const int n = INTEGER(dimV)[0];
  const int cols = INTEGER(dimC)[1];
  if (!eigencore_int_indexable(static_cast<int64_t>(n) * (2 * static_cast<int64_t>(cols) + 1))) {
    error("native Arnoldi Ritz vector dimensions exceed LP64 BLAS/R integer range");
  }
  const size_t sz = static_cast<size_t>(m) * static_cast<size_t>(std::max(cols, 1));
  double* cre = reinterpret_cast<double*>(R_alloc(sz, sizeof(double)));
  double* cim = reinterpret_cast<double*>(R_alloc(sz, sizeof(double)));
  const Rcomplex* coeff = COMPLEX(coefficients_);
  for (size_t i = 0; i < static_cast<size_t>(m) * static_cast<size_t>(cols); ++i) {
    cre[i] = coeff[i].r;
    cim[i] = coeff[i].i;
  }
  return arnoldi_form_ritz_vectors(REAL(V_), n, m, cre, m, cim, m, cols);
}

extern "C" SEXP eigencore_arnoldi_ritz(SEXP V_, SEXP H_, SEXP iterations_) {
  const int m = asInteger(iterations_);
  int n = 0;
  int h_rows = 0;
  arnoldi_check_vh(V_, H_, m, &n, &h_rows, "native Arnoldi Ritz extraction");
  double* wr = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m), sizeof(double)));
  double* wi = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m), sizeof(double)));
  double* cre = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * m, sizeof(double)));
  double* cim = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(m) * m, sizeof(double)));
  arnoldi_projected_eigen(REAL(H_), h_rows, m, wr, wi, cre, cim);
  SEXP values_ = PROTECT(allocVector(CPLXSXP, m));
  Rcomplex* values = COMPLEX(values_);
  for (int col = 0; col < m; ++col) {
    values[col].r = wr[col];
    values[col].i = wi[col];
  }
  SEXP vectors_ = PROTECT(arnoldi_form_ritz_vectors(REAL(V_), n, m, cre, m, cim, m, m));
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
