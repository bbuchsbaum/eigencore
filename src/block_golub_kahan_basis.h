#ifndef EIGENCORE_BLOCK_GOLUB_KAHAN_BASIS_H
#define EIGENCORE_BLOCK_GOLUB_KAHAN_BASIS_H

#include <cstddef>
#include "native_operators.h"

struct BlockGolubKahanBasisScratch {
  double* Z_v = nullptr;
  double* Z_u = nullptr;
  double* coeff = nullptr;
  double* tmp = nullptr;
  size_t bytes = 0;
  bool transient = false;
};

void block_golub_kahan_basis_scratch_free(BlockGolubKahanBasisScratch* scratch);

int block_golub_kahan_basis_scratch_alloc(BlockGolubKahanBasisScratch* scratch,
                                          int m,
                                          int n,
                                          int max_subspace,
                                          int block_size);

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
  int* cached_start_used_out);

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
                                       int* cached_start_used_out);

#endif
