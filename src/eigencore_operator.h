#ifndef EIGENCORE_OPERATOR_H
#define EIGENCORE_OPERATOR_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EIGENCORE_TRANSPOSE_NONE = 0,
  EIGENCORE_TRANSPOSE_ADJOINT = 1
} EigencoreTranspose;

typedef enum {
  EIGENCORE_STRUCTURE_GENERAL = 0,
  EIGENCORE_STRUCTURE_HERMITIAN = 1,
  EIGENCORE_STRUCTURE_DIAGONAL = 2,
  EIGENCORE_STRUCTURE_SYMMETRIC = 3,
  EIGENCORE_STRUCTURE_TRIANGULAR = 4
} EigencoreStructure;

typedef enum {
  EIGENCORE_SCALAR_F64 = 0,
  EIGENCORE_SCALAR_C128 = 1
} EigencoreScalarType;

typedef struct EigencoreWorkspace {
  int64_t allocation_count;
  int64_t bytes_allocated;
  void* scratch;
  int64_t scratch_bytes;
} EigencoreWorkspace;

typedef int (*EigencoreApplyFn)(void* impl,
                                EigencoreTranspose op,
                                int64_t block_cols,
                                const double* X,
                                int64_t ldx,
                                double alpha,
                                double beta,
                                double* Y,
                                int64_t ldy,
                                EigencoreWorkspace* workspace);

typedef struct {
  int64_t rows;
  int64_t cols;
  EigencoreScalarType scalar;
  EigencoreStructure structure;
  int has_adjoint;
  double frobenius_upper;
  double two_norm_upper;
  void* impl;
  EigencoreApplyFn apply;
} EigencoreBlockOperator;

#ifdef __cplusplus
}
#endif

#endif
