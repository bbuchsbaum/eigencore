#ifndef EIGENCORE_COMMON_H
#define EIGENCORE_COMMON_H

#include <R.h>
#include <R_ext/BLAS.h>
#include <chrono>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cstring>

struct NativeBlockStageSeconds {
  double apply = 0.0;
  double recurrence = 0.0;
  double reorthogonalization = 0.0;
  double projected_solve = 0.0;
  double projection_update = 0.0;
  double projection_copy = 0.0;
  double projected_eigensolve = 0.0;
  double selected_vector_copy = 0.0;
  double ritz_residual = 0.0;
  double ritz_vector_form = 0.0;
  double ritz_operator_apply = 0.0;
  double ritz_norm = 0.0;
  double ritz_final_polish = 0.0;
  double locking = 0.0;
  double restart = 0.0;
};

struct NativeBlockRestartHistory {
  int capacity = 0;
  int length = 0;
  int* restart = nullptr;
  int* m_active = nullptr;
  int* selected_count = nullptr;
  int* locked_before = nullptr;
  int* locked_after = nullptr;
  int* nconv_wanted = nullptr;
  double* max_residual = nullptr;
  double* max_backward_error = nullptr;
};

static inline std::chrono::steady_clock::time_point native_timer_now() {
  return std::chrono::steady_clock::now();
}

static inline double native_timer_elapsed(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

static inline int eigencore_int_indexable(int64_t value) {
  return value >= 0 && value <= static_cast<int64_t>(INT_MAX);
}

// Status returned by the scalar Krylov kernels when the start vector or a
// recurrence coefficient is NaN/Inf.
#define EIGENCORE_STATUS_NONFINITE (-10)

static inline void eigencore_check_nonfinite_status(int status) {
  if (status == EIGENCORE_STATUS_NONFINITE) {
    error("non-finite value encountered (check input for NA/NaN/Inf)");
  }
}

static inline void eigencore_apply_status_error(const char* context, int status) {
  eigencore_check_nonfinite_status(status);
  if (status == -2) {
    error("%s failed: dimensions exceed LP64 BLAS/R integer range; LAPACK64 is not enabled",
          context);
  }
  error("%s failed with status=%d", context, status);
}

// Dot product / Euclidean norm through BLAS (no long double accumulators: those
// defeat vectorisation and are software binary128 on aarch64 Linux). The norm
// takes the fast sqrt(ddot) path when the sum of squares is safely inside the
// representable range and otherwise falls back to the scaled dnrm2, so
// operators scaled by 1e+-150 cannot overflow or underflow the norm.
static inline double ec_dot(const double* x, const double* y, int n) {
  if (n <= 0) {
    return 0.0;
  }
  int one = 1;
  return F77_CALL(ddot)(&n, x, &one, y, &one);
}

static inline double ec_norm2(const double* x, int n) {
  if (n <= 0) {
    return 0.0;
  }
  int one = 1;
  const double ss = F77_CALL(ddot)(&n, x, &one, x, &one);
  if (ss >= 1e-280 && ss <= 1e280) {
    return std::sqrt(ss);
  }
  if (ss == 0.0) {
    return 0.0;
  }
  return F77_CALL(dnrm2)(&n, x, &one);
}

// ||av - theta * v||_2 without a temporary and without long double: four
// independent double accumulators (keeps the FP pipeline busy), with a scaled
// two-pass fallback when the sum of squares leaves the safe range.
static inline double ec_residual_norm(const double* av, const double* v,
                                      double theta, int n) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double d0 = av[i] - theta * v[i];
    const double d1 = av[i + 1] - theta * v[i + 1];
    const double d2 = av[i + 2] - theta * v[i + 2];
    const double d3 = av[i + 3] - theta * v[i + 3];
    s0 += d0 * d0;
    s1 += d1 * d1;
    s2 += d2 * d2;
    s3 += d3 * d3;
  }
  for (; i < n; ++i) {
    const double d = av[i] - theta * v[i];
    s0 += d * d;
  }
  const double ss = (s0 + s1) + (s2 + s3);
  if ((ss >= 1e-280 && ss <= 1e280) || ss == 0.0 || ISNAN(ss)) {
    return std::sqrt(ss);
  }
  double amax = 0.0;
  for (i = 0; i < n; ++i) {
    const double d = std::fabs(av[i] - theta * v[i]);
    if (d > amax) {
      amax = d;
    }
  }
  if (!(amax > 0.0) || !R_FINITE(amax)) {
    return amax;
  }
  double scaled = 0.0;
  for (i = 0; i < n; ++i) {
    const double d = (av[i] - theta * v[i]) / amax;
    scaled += d * d;
  }
  return amax * std::sqrt(scaled);
}

static inline void combine_basis_columns_small(const double* basis,
                                               int n,
                                               int basis_cols,
                                               const double* coeff,
                                               int coeff_ld,
                                               int out_cols,
                                               double* out) {
  std::memset(out, 0,
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(out_cols));
  for (int p = 0; p < out_cols; ++p) {
    double* y = out + static_cast<int64_t>(p) * n;
    for (int col = 0; col < basis_cols; ++col) {
      const double a = coeff[col + static_cast<int64_t>(p) * coeff_ld];
      if (a == 0.0) {
        continue;
      }
      const double* x = basis + static_cast<int64_t>(col) * n;
      for (int row = 0; row < n; ++row) {
        y[row] += a * x[row];
      }
    }
  }
}

// out (n x out_cols) = basis (n x basis_cols) * coeff (basis_cols x out_cols).
// Uses BLAS-3 dgemm unless the whole product is tiny; the decision is driven
// by the total work n * basis_cols * out_cols, not by out_cols alone.
static inline void combine_basis_columns(const double* basis,
                                         int n,
                                         int basis_cols,
                                         const double* coeff,
                                         int coeff_ld,
                                         int out_cols,
                                         double* out) {
  if (out_cols <= 0 || n <= 0) {
    return;
  }
  const int64_t work = static_cast<int64_t>(n) * basis_cols * out_cols;
  if (basis_cols <= 0 || work <= 4096) {
    combine_basis_columns_small(basis, n, basis_cols, coeff, coeff_ld,
                                out_cols, out);
    return;
  }
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &out_cols, &basis_cols,
                  &one, basis, &n, coeff, &coeff_ld,
                  &zero, out, &n FCONE FCONE);
}

#endif
