#ifndef EIGENCORE_NATIVE_OPERATORS_H
#define EIGENCORE_NATIVE_OPERATORS_H

#include <Rinternals.h>
#include <stdint.h>
#include "eigencore_operator.h"

struct DenseColumnMajorOperator {
  int64_t rows;
  int64_t cols;
  double* values;
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
  SEXP apply;
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

#endif
