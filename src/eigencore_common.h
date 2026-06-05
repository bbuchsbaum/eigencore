#ifndef EIGENCORE_COMMON_H
#define EIGENCORE_COMMON_H

#include <R.h>
#include <chrono>
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

static inline void eigencore_apply_status_error(const char* context, int status) {
  if (status == -2) {
    error("%s failed: dimensions exceed LP64 BLAS/R integer range; LAPACK64 is not enabled",
          context);
  }
  error("%s failed with status=%d", context, status);
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

#endif
