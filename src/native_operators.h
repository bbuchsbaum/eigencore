#ifndef EIGENCORE_NATIVE_OPERATORS_H
#define EIGENCORE_NATIVE_OPERATORS_H

#include <Rinternals.h>
#include <R_ext/Complex.h>
#include <stdint.h>
#include "eigencore_operator.h"

struct DenseColumnMajorOperator {
  int64_t rows;
  int64_t cols;
  double* values;
};

struct DenseComplexColumnMajorOperator {
  int64_t rows;
  int64_t cols;
  Rcomplex* values;
};

struct CSCOperator {
  int64_t rows;
  int64_t cols;
  const int* row_idx;
  const int* col_ptr;
  const double* values;
};

struct DiagonalOperator {
  int64_t rows;
  const double* values;
  bool unit;
};

struct RApplyOperator {
  int64_t rows;
  int64_t cols;
  SEXP apply;
  SEXP apply_adjoint;
};

// Symmetric normal-equations view of a rectangular base operator A: applies
// A^T A (side == 0, acting on R^cols) or A A^T (side == 1, acting on R^rows)
// without materializing the Gram matrix. The intermediate product lives in a
// caller-provided scratch buffer of at least inner_dim * block_cols doubles,
// where inner_dim is rows for side 0 and cols for side 1.
struct NormalEquationsOperator {
  void* base_impl;
  EigencoreApplyFn base_apply;
  int64_t rows;              // rows of A
  int64_t cols;              // cols of A
  int side;                  // 0: A^T A, 1: A A^T
  double* scratch;           // inner_dim * scratch_block_capacity doubles
  int64_t scratch_block_capacity;
};

struct DenseShiftInvertOperator {
  int n;
  double* lu;
  int* pivots;
  double* work;
};

struct DenseGeneralizedShiftInvertOperator {
  int n;
  double* lu;
  int* pivots;
  double* chol;
  double* rhs;
  double* sol;
};

struct TridiagonalShiftInvertOperator {
  int n;
  const double* lower;
  const double* cprime;
  const double* denom;
  double* work;
};

struct TridiagonalGeneralizedShiftInvertOperator {
  int n;
  const double* lower;
  const double* cprime;
  const double* denom;
  const double* sqrt_metric;
  double* work;
};


extern "C" int eigencore_dense_apply(void* impl,
                                      EigencoreTranspose op,
                                      int64_t block_cols,
                                      const double* X,
                                      int64_t ldx,
                                      double alpha,
                                      double beta,
                                      double* Y,
                                      int64_t ldy,
                                      EigencoreWorkspace* workspace);

extern "C" int eigencore_dense_complex_apply(void* impl,
                                             EigencoreTranspose op,
                                             int64_t block_cols,
                                             const Rcomplex* X,
                                             int64_t ldx,
                                             Rcomplex alpha,
                                             Rcomplex beta,
                                             Rcomplex* Y,
                                             int64_t ldy,
                                             EigencoreWorkspace* workspace);

extern "C" int eigencore_dense_shift_invert_apply(void* impl,
                                                   EigencoreTranspose op,
                                                   int64_t block_cols,
                                                   const double* X,
                                                   int64_t ldx,
                                                   double alpha,
                                                   double beta,
                                                   double* Y,
                                                   int64_t ldy,
                                                   EigencoreWorkspace* workspace);

extern "C" int eigencore_dense_generalized_shift_invert_apply(void* impl,
                                                               EigencoreTranspose op,
                                                               int64_t block_cols,
                                                               const double* X,
                                                               int64_t ldx,
                                                               double alpha,
                                                               double beta,
                                                               double* Y,
                                                               int64_t ldy,
                                                               EigencoreWorkspace* workspace);

extern "C" int eigencore_tridiagonal_shift_invert_apply(void* impl,
                                                         EigencoreTranspose op,
                                                         int64_t block_cols,
                                                         const double* X,
                                                         int64_t ldx,
                                                         double alpha,
                                                         double beta,
                                                         double* Y,
                                                         int64_t ldy,
                                                         EigencoreWorkspace* workspace);

extern "C" int eigencore_tridiagonal_generalized_shift_invert_apply(
    void* impl,
    EigencoreTranspose op,
    int64_t block_cols,
    const double* X,
    int64_t ldx,
    double alpha,
    double beta,
    double* Y,
    int64_t ldy,
    EigencoreWorkspace* workspace);

extern "C" int eigencore_csc_apply(void* impl,
                                    EigencoreTranspose op,
                                    int64_t block_cols,
                                    const double* X,
                                    int64_t ldx,
                                    double alpha,
                                    double beta,
                                    double* Y,
                                    int64_t ldy,
                                    EigencoreWorkspace* workspace);

extern "C" int eigencore_diagonal_apply(void* impl,
                                         EigencoreTranspose op,
                                         int64_t block_cols,
                                         const double* X,
                                         int64_t ldx,
                                         double alpha,
                                         double beta,
                                         double* Y,
                                         int64_t ldy,
                                         EigencoreWorkspace* workspace);

extern "C" int eigencore_r_operator_apply(void* impl,
                                           EigencoreTranspose op,
                                           int64_t block_cols,
                                           const double* X,
                                           int64_t ldx,
                                           double alpha,
                                           double beta,
                                           double* Y,
                                           int64_t ldy,
                                           EigencoreWorkspace* workspace);

extern "C" int eigencore_normal_equations_apply(void* impl,
                                                 EigencoreTranspose op,
                                                 int64_t block_cols,
                                                 const double* X,
                                                 int64_t ldx,
                                                 double alpha,
                                                 double beta,
                                                 double* Y,
                                                 int64_t ldy,
                                                 EigencoreWorkspace* workspace);

#endif
