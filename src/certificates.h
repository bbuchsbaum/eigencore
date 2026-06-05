#ifndef EIGENCORE_CERTIFICATES_H
#define EIGENCORE_CERTIFICATES_H

#include <Rinternals.h>
#include "eigencore_operator.h"

SEXP native_operator_svd_certificate_cached_av(void* impl,
                                               EigencoreApplyFn apply,
                                               int m,
                                               int n,
                                               double norm_A,
                                               SEXP d_,
                                               SEXP u_,
                                               SEXP v_,
                                               SEXP av_,
                                               SEXP tol_);

#endif
