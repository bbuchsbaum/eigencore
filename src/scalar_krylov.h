#ifndef EIGENCORE_SCALAR_KRYLOV_H
#define EIGENCORE_SCALAR_KRYLOV_H

#include "native_operators.h"

int native_golub_kahan_run(void* impl,
                           EigencoreApplyFn apply,
                           int m,
                           int n,
                           int maxit,
                           int k,
                           int target_kind,
                           double tol,
                           int enable_projected_stop,
                           int use_blas_reorthogonalization,
                           const double* start,
                           double* U,
                           double* V,
                           double* alpha,
                           double* beta,
                           int* iterations,
                           int* matvecs,
                           int* projected_stop,
                           int* projected_nconv,
                           double* projected_max_residual,
                           int* projected_checks,
                           double* projected_seconds,
                           double* stage_apply_seconds,
                           double* stage_recurrence_seconds,
                           double* stage_reorthogonalization_seconds,
                           int* reorthogonalization_passes,
                           int reorthogonalize_u,
                           int reorthogonalize_v);

#endif
