#ifndef EIGENCORE_LAPACK_COMPAT_H
#define EIGENCORE_LAPACK_COMPAT_H

// Compatibility shims for building against R older than 4.4.0, whose
// R_ext/Lapack.h neither defines the La_INT / La_LGL aliases nor declares
// the complex QZ drivers (zggev, zgges) this package uses. On R >= 4.4.0
// this header is a no-op so the stock declarations apply unchanged.
//
// Include this after <R_ext/BLAS.h> and <R_ext/Lapack.h>.

#include <Rversion.h>

#if R_VERSION < R_Version(4, 4, 0)

#ifndef La_INT
typedef int La_INT;
#define La_INT La_INT
#endif

#ifndef La_LGL
typedef int La_LGL;
#define La_LGL La_LGL
#endif

#ifdef __cplusplus
extern "C" {
#endif

La_extern void
F77_NAME(zggev)(const char* jobvl, const char* jobvr, const int* n,
                Rcomplex* a, const int* lda, Rcomplex* b, const int* ldb,
                Rcomplex* alpha, Rcomplex* beta,
                Rcomplex* vl, const int* ldvl, Rcomplex* vr, const int* ldvr,
                Rcomplex* work, const int* lwork, double* rwork,
                int* info FCLEN FCLEN);

La_extern void
F77_NAME(zgges)(const char* jobvsl, const char* jobvsr, const char* sort,
                void* selctg, const int* n,
                Rcomplex* a, const int* lda, Rcomplex* b, const int* ldb,
                int* sdim, Rcomplex* alpha, Rcomplex* beta,
                Rcomplex* vsl, const int* ldvsl,
                Rcomplex* vsr, const int* ldvsr,
                Rcomplex* work, const int* lwork, double* rwork, int* bwork,
                int* info FCLEN FCLEN FCLEN);

#ifdef __cplusplus
}
#endif

#endif  // R_VERSION < 4.4.0

#endif  // EIGENCORE_LAPACK_COMPAT_H
