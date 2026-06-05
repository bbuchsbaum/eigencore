#ifndef EIGENCORE_PROJECTION_GOLUB_KAHAN_RITZ_H
#define EIGENCORE_PROJECTION_GOLUB_KAHAN_RITZ_H

#include <Rinternals.h>


SEXP eigencore_block_golub_kahan_ritz_from_ptr(const double* V,
                                               int n,
                                               const double* AV,
                                               int m,
                                               int p,
                                               int rank,
                                               int target_kind);

extern "C" SEXP eigencore_block_golub_kahan_ritz(SEXP V_, SEXP AV_,
                                                 SEXP rank_, SEXP target_kind_,
                                                 SEXP active_p_);

SEXP eigencore_golub_kahan_ritz_from_ptr(const double* U,
                                         const double* V,
                                         int m,
                                         int n,
                                         int p,
                                         const double* alpha,
                                         const double* beta,
                                         int rank,
                                         int target_kind);

extern "C" SEXP eigencore_golub_kahan_ritz(SEXP U_, SEXP V_, SEXP alpha_,
                                           SEXP beta_, SEXP rank_,
                                           SEXP target_kind_, SEXP active_p_);

#endif
