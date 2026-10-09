#include <cfloat>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>
#include "eigencore_common.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"

typedef La_LGL (*eigencore_dgges_select_fn)(double*, double*, double*);
typedef void (*eigencore_dgges_fn)(
  const char*, const char*, const char*, eigencore_dgges_select_fn,
  const La_INT*, double*, const La_INT*, double*, const La_INT*, La_INT*,
  double*, double*, double*, double*, const La_INT*, double*, const La_INT*,
  double*, const La_INT*, La_LGL*, La_INT* FCLEN FCLEN FCLEN
);

static double qz_real_scale(double alphar, double alphai, double beta) {
  double scale = 1.0;
  const double alpha_mod = hypot(alphar, alphai);
  const double beta_mod = fabs(beta);
  if (alpha_mod > scale) {
    scale = alpha_mod;
  }
  if (beta_mod > scale) {
    scale = beta_mod;
  }
  return scale;
}

static double qz_complex_scale(const Rcomplex* alpha, const Rcomplex* beta) {
  double scale = 1.0;
  const double alpha_mod = hypot(alpha->r, alpha->i);
  const double beta_mod = hypot(beta->r, beta->i);
  if (alpha_mod > scale) {
    scale = alpha_mod;
  }
  if (beta_mod > scale) {
    scale = beta_mod;
  }
  return scale;
}

static La_LGL qz_select_real_finite(double* alphar, double* alphai, double* beta) {
  const double tol = sqrt(DBL_EPSILON) * qz_real_scale(*alphar, *alphai, *beta);
  return fabs(*beta) > tol ? TRUE : FALSE;
}

static La_LGL qz_select_real_infinite(double* alphar, double* alphai, double* beta) {
  const double tol = sqrt(DBL_EPSILON) * qz_real_scale(*alphar, *alphai, *beta);
  const bool beta_zero = fabs(*beta) <= tol;
  const bool alpha_zero = hypot(*alphar, *alphai) <= tol;
  return (beta_zero && !alpha_zero) ? TRUE : FALSE;
}

static La_LGL qz_select_complex_finite(Rcomplex* alpha, Rcomplex* beta) {
  const double tol = sqrt(DBL_EPSILON) * qz_complex_scale(alpha, beta);
  return hypot(beta->r, beta->i) > tol ? TRUE : FALSE;
}

static La_LGL qz_select_complex_infinite(Rcomplex* alpha, Rcomplex* beta) {
  const double tol = sqrt(DBL_EPSILON) * qz_complex_scale(alpha, beta);
  const bool beta_zero = hypot(beta->r, beta->i) <= tol;
  const bool alpha_zero = hypot(alpha->r, alpha->i) <= tol;
  return (beta_zero && !alpha_zero) ? TRUE : FALSE;
}

static eigencore_dgges_select_fn qz_real_selector(int sort_code) {
  switch (sort_code) {
  case 0:
    return NULL;
  case 1:
    return qz_select_real_finite;
  case 2:
    return qz_select_real_infinite;
  default:
    error("unsupported generalized Schur sort code");
    return NULL;
  }
}

static void* qz_complex_selector(int sort_code) {
  switch (sort_code) {
  case 0:
    return NULL;
  case 1:
    return reinterpret_cast<void*>(qz_select_complex_finite);
  case 2:
    return reinterpret_cast<void*>(qz_select_complex_infinite);
  default:
    error("unsupported generalized Schur sort code");
    return NULL;
  }
}

// Run LAPACK dsyevr (MRRR, the driver base R's eigen() uses) on a private
// copy of the n x n column-major matrix A. RANGE = 'A' when il <= 0,
// otherwise RANGE = 'I' over [il, iu] (1-based). Eigenvalues land in w
// (length n), eigenvectors in z (n x count) when want_vectors.
// Returns LAPACK info; *m_found receives the number of eigenvalues found.
static int run_dsyevr(const double* A, int n, bool want_vectors, int il,
                      int iu, double* w, double* z, int* m_found) {
  const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
  double* work_matrix = reinterpret_cast<double*>(R_alloc(nn, sizeof(double)));
  std::memcpy(work_matrix, A, sizeof(double) * nn);
  const bool all = (il <= 0);
  const int count = all ? n : (iu - il + 1);
  char jobz = want_vectors ? 'V' : 'N';
  char range = all ? 'A' : 'I';
  char uplo = 'U';
  double vl = 0.0;
  double vu = 0.0;
  double abstol = 0.0;
  int il_la = all ? 1 : il;
  int iu_la = all ? n : iu;
  int ldz = want_vectors ? n : 1;
  double z_dummy = 0.0;
  double* z_ptr = want_vectors ? z : &z_dummy;
  int* isuppz = reinterpret_cast<int*>(
    R_alloc(static_cast<size_t>(2 * (count > 1 ? count : 1)), sizeof(int))
  );
  int info = 0;
  int lwork = -1;
  int liwork = -1;
  double work_query = 0.0;
  int iwork_query = 0;
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, work_matrix, &n,
                   &vl, &vu, &il_la, &iu_la, &abstol,
                   m_found, w, z_ptr, &ldz,
                   isuppz, &work_query, &lwork,
                   &iwork_query, &liwork, &info FCONE FCONE FCONE);
  if (info != 0) {
    return info;
  }
  lwork = static_cast<int>(work_query);
  liwork = iwork_query;
  if (lwork < 26 * n) {
    lwork = 26 * n;
  }
  if (liwork < 10 * n) {
    liwork = 10 * n;
  }
  double* work = reinterpret_cast<double*>(
    R_alloc(static_cast<size_t>(lwork), sizeof(double))
  );
  int* iwork = reinterpret_cast<int*>(
    R_alloc(static_cast<size_t>(liwork), sizeof(int))
  );
  F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, work_matrix, &n,
                   &vl, &vu, &il_la, &iu_la, &abstol,
                   m_found, w, z_ptr, &ldz,
                   isuppz, work, &lwork,
                   iwork, &liwork, &info FCONE FCONE FCONE);
  return info;
}

// QR-iteration fallback (dsyev), used only if dsyevr reports an internal
// failure (info > 0). Writes vectors into z (n x n) when want_vectors.
static void run_dsyev_fallback(const double* A, int n, bool want_vectors,
                               double* w, double* z) {
  const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
  double* a = want_vectors ? z :
    reinterpret_cast<double*>(R_alloc(nn, sizeof(double)));
  std::memcpy(a, A, sizeof(double) * nn);
  char jobz = want_vectors ? 'V' : 'N';
  char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsyev)(&jobz, &uplo, &n, a, &n, w,
                  &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dsyev workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  if (lwork < 3 * n) {
    lwork = 3 * n;
  }
  double* work = reinterpret_cast<double*>(
    R_alloc(static_cast<size_t>(lwork), sizeof(double))
  );
  F77_CALL(dsyev)(&jobz, &uplo, &n, a, &n, w, work, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dsyev failed with info=%d", info);
  }
}

static bool vectors_flag(SEXP flag_) {
  const int flag = asLogical(flag_);
  if (flag == NA_LOGICAL) {
    error("vectors must be TRUE or FALSE");
  }
  return flag != 0;
}

// Full dense symmetric eigendecomposition (ascending values). Uses dsyevr
// (MRRR); with vectors = FALSE it runs jobz = 'N' (tridiagonal reduction
// plus dsterf) and returns vectors = NULL. Falls back to dsyev on an
// internal dsyevr failure; `driver` records which one ran.
extern "C" SEXP eigencore_dense_symmetric_eigen(SEXP A_, SEXP vectors_flag_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  if (n != ncolA) {
    error("A must be square");
  }
  const bool want_vectors = vectors_flag(vectors_flag_);

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(want_vectors ? allocMatrix(REALSXP, n, n) : R_NilValue);
  const char* driver = "dsyevr";
  if (n > 0) {
    int m_found = 0;
    const int info = run_dsyevr(REAL(A_), n, want_vectors, 0, 0,
                                REAL(values_),
                                want_vectors ? REAL(vectors_) : NULL,
                                &m_found);
    if (info < 0) {
      error("LAPACK dsyevr failed with info=%d", info);
    }
    if (info > 0 || m_found != n) {
      run_dsyev_fallback(REAL(A_), n, want_vectors, REAL(values_),
                         want_vectors ? REAL(vectors_) : NULL);
      driver = "dsyev";
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, mkString(driver));
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("driver"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

// QR-iteration (dsyev) full symmetric eigendecomposition with vectors. Kept
// as a benchmark/diagnostic backend next to the dsyevd variant; production
// callers use eigencore_dense_symmetric_eigen (dsyevr).
extern "C" SEXP eigencore_dense_symmetric_eigen_dsyev(SEXP A_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  if (n != INTEGER(dimA)[1]) {
    error("A must be square");
  }
  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, n));
  if (n > 0) {
    run_dsyev_fallback(REAL(A_), n, true, REAL(values_), REAL(vectors_));
  }
  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_symmetric_eigen_dsyevd(SEXP A_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  if (n != ncolA) {
    error("A must be square");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(duplicate(A_));
  if (n > 0) {
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    F77_CALL(dsyevd)(&jobz, &uplo, &n, REAL(vectors_), &n, REAL(values_),
                     &work_query, &lwork, &iwork_query, &liwork,
                     &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyevd workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    const int64_t lwork_min = 1 + 6 * static_cast<int64_t>(n) +
      2 * static_cast<int64_t>(n) * static_cast<int64_t>(n);
    if (lwork_min > static_cast<int64_t>(INT_MAX)) {
      error("dense symmetric eigensolver workspace exceeds LP64 LAPACK integer range");
    }
    if (lwork < lwork_min) {
      lwork = static_cast<int>(lwork_min);
    }
    if (liwork < 3 + 5 * n) {
      liwork = 3 + 5 * n;
    }
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    SEXP iwork_ = PROTECT(allocVector(INTSXP, liwork));
    F77_CALL(dsyevd)(&jobz, &uplo, &n, REAL(vectors_), &n, REAL(values_),
                     REAL(work_), &lwork, INTEGER(iwork_), &liwork,
                     &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyevd failed with info=%d", info);
    }
    UNPROTECT(2);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

// Complex Hermitian eigendecomposition via zheev: R's LAPACK interface
// (R_ext/Lapack.h) does not declare zheevr, so the MRRR driver is not used
// here. With vectors = FALSE it runs jobz = 'N' and returns vectors = NULL.
extern "C" SEXP eigencore_dense_complex_hermitian_eigen(SEXP A_, SEXP vectors_flag_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_)) {
    error("A must be a complex matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  if (n != ncolA) {
    error("A must be square");
  }
  const bool want_vectors = vectors_flag(vectors_flag_);

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(duplicate(A_));
  if (n > 0) {
    char jobz = want_vectors ? 'V' : 'N';
    char uplo = 'U';
    int info = 0;
    int lwork = -1;
    Rcomplex work_query;
    const int lrwork = (3 * n - 2 > 1) ? (3 * n - 2) : 1;
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lrwork), sizeof(double))
    );
    F77_CALL(zheev)(&jobz, &uplo, &n, COMPLEX(vectors_), &n, REAL(values_),
                    &work_query, &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zheev workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    F77_CALL(zheev)(&jobz, &uplo, &n, COMPLEX(vectors_), &n, REAL(values_),
                    COMPLEX(work_), &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zheev failed with info=%d", info);
    }
    UNPROTECT(1);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, want_vectors ? vectors_ : R_NilValue);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_general_eigen(SEXP A_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_)) {
    error("A must be a complex matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  if (n != ncolA) {
    error("A must be square");
  }

  SEXP values_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP vectors_ = PROTECT(allocMatrix(CPLXSXP, n, n));
  SEXP work_matrix_ = PROTECT(duplicate(A_));
  if (n > 0) {
    char jobvl = 'N';
    char jobvr = 'V';
    int ldvl = 1;
    int ldvr = n;
    int info = 0;
    int lwork = -1;
    Rcomplex vl_dummy;
    Rcomplex work_query;
    const int lrwork = (2 * n > 1) ? (2 * n) : 1;
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lrwork), sizeof(double))
    );
    F77_CALL(zgeev)(&jobvl, &jobvr, &n, COMPLEX(work_matrix_), &n,
                    COMPLEX(values_), &vl_dummy, &ldvl, COMPLEX(vectors_),
                    &ldvr, &work_query, &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgeev workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    work_matrix_ = PROTECT(duplicate(A_));
    F77_CALL(zgeev)(&jobvl, &jobvr, &n, COMPLEX(work_matrix_), &n,
                    COMPLEX(values_), &vl_dummy, &ldvl, COMPLEX(vectors_),
                    &ldvr, COMPLEX(work_), &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgeev failed with info=%d", info);
    }
    UNPROTECT(2);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(5);
  return out_;
  EIGENCORE_ENTRY_END
}

// Expand LAPACK packed real eigenvector storage (conjugate pairs occupy two
// consecutive real columns) into full complex columns. Errors loudly on
// non-conforming pairing rather than emit a silently wrong real-only vector.
static void unpack_real_pencil_vectors(const double* packed,
                                       const double* alphai,
                                       int n, Rcomplex* out,
                                       const char* routine) {
  int j = 0;
  while (j < n) {
    if (alphai[j] > 0.0 && j + 1 < n) {
      for (int row = 0; row < n; ++row) {
        const double re = packed[row + static_cast<int64_t>(j) * n];
        const double im = packed[row + static_cast<int64_t>(j + 1) * n];
        out[row + static_cast<int64_t>(j) * n].r = re;
        out[row + static_cast<int64_t>(j) * n].i = im;
        out[row + static_cast<int64_t>(j + 1) * n].r = re;
        out[row + static_cast<int64_t>(j + 1) * n].i = -im;
      }
      j += 2;
    } else {
      if (alphai[j] != 0.0) {
        error("LAPACK %s returned a non-conforming complex eigenvalue at "
              "index %d without a consecutive conjugate partner", routine, j);
      }
      for (int row = 0; row < n; ++row) {
        out[row + static_cast<int64_t>(j) * n].r =
          packed[row + static_cast<int64_t>(j) * n];
        out[row + static_cast<int64_t>(j) * n].i = 0.0;
      }
      ++j;
    }
  }
}

extern "C" SEXP eigencore_dense_generalized_pencil_eigen(SEXP A_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(B_)) {
    error("A and B must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("A and B must be square matrices with the same dimension");
  }

  SEXP alphar_ = PROTECT(allocVector(REALSXP, n));
  SEXP alphai_ = PROTECT(allocVector(REALSXP, n));
  SEXP beta_real_ = PROTECT(allocVector(REALSXP, n));
  SEXP vr_real_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP vl_real_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP rconde_ = PROTECT(allocVector(REALSXP, n));
  SEXP rcondv_ = PROTECT(allocVector(REALSXP, n));
  double abnrm = 0.0;
  double bbnrm = 0.0;

  if (n > 0) {
    // DGGEVX: balance the pencil ('B' = permute and scale), compute left and
    // right eigenvectors, and both eigenvalue and eigenvector reciprocal
    // condition numbers ('B'). abnrm/bbnrm are the one-norms of the balanced
    // matrices; the R layer uses them for scale-aware alpha/beta
    // classification.
    char balanc = 'B';
    char jobvl = 'V';
    char jobvr = 'V';
    char sense = 'B';
    int ldvl = n;
    int ldvr = n;
    int ilo = 0;
    int ihi = 0;
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    double* lscale = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    double* rscale = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(n) + 6, sizeof(int))
    );
    La_LGL* bwork = reinterpret_cast<La_LGL*>(
      R_alloc(static_cast<size_t>(n), sizeof(La_LGL))
    );

    SEXP Awork_ = PROTECT(duplicate(A_));
    SEXP Bwork_ = PROTECT(duplicate(B_));
    F77_CALL(dggevx)(&balanc, &jobvl, &jobvr, &sense, &n,
                     REAL(Awork_), &n, REAL(Bwork_), &n,
                     REAL(alphar_), REAL(alphai_), REAL(beta_real_),
                     REAL(vl_real_), &ldvl, REAL(vr_real_), &ldvr,
                     &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm,
                     REAL(rconde_), REAL(rcondv_),
                     &work_query, &lwork, iwork, bwork, &info
                     FCONE FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK dggevx workspace query failed with info=%d", info);
    }
    UNPROTECT(2);

    lwork = static_cast<int>(work_query);
    if (lwork < 6 * n) {
      lwork = 6 * n;
    }
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    Awork_ = PROTECT(duplicate(A_));
    Bwork_ = PROTECT(duplicate(B_));
    F77_CALL(dggevx)(&balanc, &jobvl, &jobvr, &sense, &n,
                     REAL(Awork_), &n, REAL(Bwork_), &n,
                     REAL(alphar_), REAL(alphai_), REAL(beta_real_),
                     REAL(vl_real_), &ldvl, REAL(vr_real_), &ldvr,
                     &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm,
                     REAL(rconde_), REAL(rcondv_),
                     REAL(work_), &lwork, iwork, bwork, &info
                     FCONE FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK dggevx failed with info=%d", info);
    }
    UNPROTECT(3);
  }

  SEXP alpha_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP beta_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP vectors_ = PROTECT(allocMatrix(CPLXSXP, n, n));
  SEXP left_vectors_ = PROTECT(allocMatrix(CPLXSXP, n, n));
  for (int j = 0; j < n; ++j) {
    COMPLEX(alpha_)[j].r = REAL(alphar_)[j];
    COMPLEX(alpha_)[j].i = REAL(alphai_)[j];
    COMPLEX(beta_)[j].r = REAL(beta_real_)[j];
    COMPLEX(beta_)[j].i = 0.0;
  }
  unpack_real_pencil_vectors(REAL(vr_real_), REAL(alphai_), n,
                             COMPLEX(vectors_), "dggevx");
  unpack_real_pencil_vectors(REAL(vl_real_), REAL(alphai_), n,
                             COMPLEX(left_vectors_), "dggevx");

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, alpha_);
  SET_VECTOR_ELT(out_, 1, beta_);
  SET_VECTOR_ELT(out_, 2, vectors_);
  SET_VECTOR_ELT(out_, 3, left_vectors_);
  SET_VECTOR_ELT(out_, 4, ScalarReal(abnrm));
  SET_VECTOR_ELT(out_, 5, ScalarReal(bbnrm));
  SET_VECTOR_ELT(out_, 6, rconde_);
  SET_VECTOR_ELT(out_, 7, rcondv_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("alpha"));
  SET_STRING_ELT(names_, 1, mkChar("beta"));
  SET_STRING_ELT(names_, 2, mkChar("vectors"));
  SET_STRING_ELT(names_, 3, mkChar("left_vectors"));
  SET_STRING_ELT(names_, 4, mkChar("abnrm"));
  SET_STRING_ELT(names_, 5, mkChar("bbnrm"));
  SET_STRING_ELT(names_, 6, mkChar("rconde"));
  SET_STRING_ELT(names_, 7, mkChar("rcondv"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(13);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_generalized_hpd_eigen(SEXP A_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_) || !isComplex(B_)) {
    error("A and B must be complex matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("A and B must be square matrices with the same dimension");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(duplicate(A_));
  SEXP Bwork_ = PROTECT(duplicate(B_));
  if (n > 0) {
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    // Windows R does not export zhegv. Reduce B = U^H U to a standard
    // Hermitian problem C = U^-H A U^-1, then back-transform vectors.
    F77_CALL(zpotrf)(&uplo, &n, COMPLEX(Bwork_), &n, &info FCONE);
    if (info != 0) {
      error("LAPACK zpotrf failed for generalized Hermitian B with info=%d", info);
    }

    Rcomplex one;
    one.r = 1.0;
    one.i = 0.0;
    char side = 'L';
    char transa = 'C';
    char diag = 'N';
    F77_CALL(ztrsm)(&side, &uplo, &transa, &diag, &n, &n, &one,
                    COMPLEX(Bwork_), &n, COMPLEX(vectors_), &n
                    FCONE FCONE FCONE FCONE);
    side = 'R';
    transa = 'N';
    F77_CALL(ztrsm)(&side, &uplo, &transa, &diag, &n, &n, &one,
                    COMPLEX(Bwork_), &n, COMPLEX(vectors_), &n
                    FCONE FCONE FCONE FCONE);

    int lwork = -1;
    Rcomplex work_query;
    const int lrwork = (3 * n - 2 > 1) ? (3 * n - 2) : 1;
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lrwork), sizeof(double))
    );
    F77_CALL(zheev)(&jobz, &uplo, &n, COMPLEX(vectors_), &n, REAL(values_),
                    &work_query, &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zheev workspace query failed for generalized Hermitian "
            "transform with info=%d", info);
    }
    lwork = static_cast<int>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    F77_CALL(zheev)(&jobz, &uplo, &n, COMPLEX(vectors_), &n, REAL(values_),
                    COMPLEX(work_), &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zheev failed for generalized Hermitian transform "
            "with info=%d", info);
    }
    UNPROTECT(1);

    side = 'L';
    transa = 'N';
    F77_CALL(ztrsm)(&side, &uplo, &transa, &diag, &n, &n, &one,
                    COMPLEX(Bwork_), &n, COMPLEX(vectors_), &n
                    FCONE FCONE FCONE FCONE);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(5);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_generalized_pencil_eigen(SEXP A_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_) || !isComplex(B_)) {
    error("A and B must be complex matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("A and B must be square matrices with the same dimension");
  }

  SEXP alpha_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP beta_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP vectors_ = PROTECT(allocMatrix(CPLXSXP, n, n));
  SEXP left_vectors_ = PROTECT(allocMatrix(CPLXSXP, n, n));

  if (n > 0) {
    // ZGGEV with left and right eigenvectors. R's bundled LAPACK subset does
    // not ship ZGGEVX, so complex pencils get left vectors but no
    // rconde/rcondv conditioning diagnostics; the R layer documents that
    // boundary explicitly.
    char jobvl = 'V';
    char jobvr = 'V';
    int ldvl = n;
    int ldvr = n;
    int info = 0;
    int lwork = -1;
    Rcomplex work_query;
    const int lrwork = (8 * n > 1) ? (8 * n) : 1;
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lrwork), sizeof(double))
    );

    SEXP Awork_ = PROTECT(duplicate(A_));
    SEXP Bwork_ = PROTECT(duplicate(B_));
    F77_CALL(zggev)(&jobvl, &jobvr, &n, COMPLEX(Awork_), &n,
                    COMPLEX(Bwork_), &n, COMPLEX(alpha_), COMPLEX(beta_),
                    COMPLEX(left_vectors_), &ldvl, COMPLEX(vectors_), &ldvr,
                    &work_query, &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zggev workspace query failed with info=%d", info);
    }
    UNPROTECT(2);

    lwork = static_cast<int>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    Awork_ = PROTECT(duplicate(A_));
    Bwork_ = PROTECT(duplicate(B_));
    F77_CALL(zggev)(&jobvl, &jobvr, &n, COMPLEX(Awork_), &n,
                    COMPLEX(Bwork_), &n, COMPLEX(alpha_), COMPLEX(beta_),
                    COMPLEX(left_vectors_), &ldvl, COMPLEX(vectors_), &ldvr,
                    COMPLEX(work_), &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zggev failed with info=%d", info);
    }
    UNPROTECT(3);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out_, 0, alpha_);
  SET_VECTOR_ELT(out_, 1, beta_);
  SET_VECTOR_ELT(out_, 2, vectors_);
  SET_VECTOR_ELT(out_, 3, left_vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(names_, 0, mkChar("alpha"));
  SET_STRING_ELT(names_, 1, mkChar("beta"));
  SET_STRING_ELT(names_, 2, mkChar("vectors"));
  SET_STRING_ELT(names_, 3, mkChar("left_vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(6);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_generalized_schur(SEXP A_, SEXP B_,
                                                   SEXP vectors_, SEXP sort_code_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(B_)) {
    error("A and B must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  La_INT n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("A and B must be square matrices with the same dimension");
  }
  const bool want_vectors = asLogical(vectors_) == TRUE;
  const int sort_code = asInteger(sort_code_);
  eigencore_dgges_select_fn selector = qz_real_selector(sort_code);
  char jobvsl = want_vectors ? 'V' : 'N';
  char jobvsr = want_vectors ? 'V' : 'N';
  char sort = sort_code == 0 ? 'N' : 'S';
  La_INT ldvsl = want_vectors ? n : 1;
  La_INT ldvsr = want_vectors ? n : 1;
  if (ldvsl < 1) {
    ldvsl = 1;
  }
  if (ldvsr < 1) {
    ldvsr = 1;
  }

  SEXP S_ = PROTECT(duplicate(A_));
  SEXP T_ = PROTECT(duplicate(B_));
  SEXP alphar_ = PROTECT(allocVector(REALSXP, n));
  SEXP alphai_ = PROTECT(allocVector(REALSXP, n));
  SEXP beta_ = PROTECT(allocVector(REALSXP, n));
  SEXP Q_ = PROTECT(allocMatrix(REALSXP, ldvsl, want_vectors ? n : 1));
  SEXP Z_ = PROTECT(allocMatrix(REALSXP, ldvsr, want_vectors ? n : 1));
  La_INT sdim = 0;

  if (n > 0) {
    La_LGL* bwork = reinterpret_cast<La_LGL*>(
      R_alloc(static_cast<size_t>(n), sizeof(La_LGL))
    );
    La_INT info = 0;
    La_INT lwork = -1;
    double work_query = 0.0;
    eigencore_dgges_fn dgges =
      reinterpret_cast<eigencore_dgges_fn>(F77_CALL(dgges));
    dgges(&jobvsl, &jobvsr, &sort, selector, &n, REAL(S_), &n, REAL(T_), &n,
          &sdim, REAL(alphar_), REAL(alphai_), REAL(beta_), REAL(Q_), &ldvsl,
          REAL(Z_), &ldvsr, &work_query, &lwork, bwork, &info
          FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgges workspace query failed with info=%d", info);
    }
    lwork = static_cast<La_INT>(work_query);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    std::memcpy(REAL(S_), REAL(A_),
                sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
    std::memcpy(REAL(T_), REAL(B_),
                sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
    dgges(&jobvsl, &jobvsr, &sort, selector, &n, REAL(S_), &n, REAL(T_), &n,
          &sdim, REAL(alphar_), REAL(alphai_), REAL(beta_), REAL(Q_), &ldvsl,
          REAL(Z_), &ldvsr, REAL(work_), &lwork, bwork, &info
          FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgges failed with info=%d", info);
    }
    UNPROTECT(1);
  }

  SEXP alpha_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP beta_complex_ = PROTECT(allocVector(CPLXSXP, n));
  for (La_INT j = 0; j < n; ++j) {
    COMPLEX(alpha_)[j].r = REAL(alphar_)[j];
    COMPLEX(alpha_)[j].i = REAL(alphai_)[j];
    COMPLEX(beta_complex_)[j].r = REAL(beta_)[j];
    COMPLEX(beta_complex_)[j].i = 0.0;
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, S_);
  SET_VECTOR_ELT(out_, 1, T_);
  SET_VECTOR_ELT(out_, 2, want_vectors ? Q_ : R_NilValue);
  SET_VECTOR_ELT(out_, 3, want_vectors ? Z_ : R_NilValue);
  SET_VECTOR_ELT(out_, 4, alpha_);
  SET_VECTOR_ELT(out_, 5, beta_complex_);
  SET_VECTOR_ELT(out_, 6, alphar_);
  SET_VECTOR_ELT(out_, 7, alphai_);
  SET_VECTOR_ELT(out_, 8, ScalarInteger(sdim));
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("S"));
  SET_STRING_ELT(names_, 1, mkChar("T"));
  SET_STRING_ELT(names_, 2, mkChar("Q"));
  SET_STRING_ELT(names_, 3, mkChar("Z"));
  SET_STRING_ELT(names_, 4, mkChar("alpha"));
  SET_STRING_ELT(names_, 5, mkChar("beta"));
  SET_STRING_ELT(names_, 6, mkChar("alphar"));
  SET_STRING_ELT(names_, 7, mkChar("alphai"));
  SET_STRING_ELT(names_, 8, mkChar("sdim"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(11);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_generalized_schur(SEXP A_, SEXP B_,
                                                           SEXP vectors_,
                                                           SEXP sort_code_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_) || !isComplex(B_)) {
    error("A and B must be complex matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  La_INT n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("A and B must be square matrices with the same dimension");
  }
  const bool want_vectors = asLogical(vectors_) == TRUE;
  const int sort_code = asInteger(sort_code_);
  void* selector = qz_complex_selector(sort_code);
  char jobvsl = want_vectors ? 'V' : 'N';
  char jobvsr = want_vectors ? 'V' : 'N';
  char sort = sort_code == 0 ? 'N' : 'S';
  La_INT ldvsl = want_vectors ? n : 1;
  La_INT ldvsr = want_vectors ? n : 1;
  if (ldvsl < 1) {
    ldvsl = 1;
  }
  if (ldvsr < 1) {
    ldvsr = 1;
  }

  SEXP S_ = PROTECT(duplicate(A_));
  SEXP T_ = PROTECT(duplicate(B_));
  SEXP alpha_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP beta_ = PROTECT(allocVector(CPLXSXP, n));
  SEXP Q_ = PROTECT(allocMatrix(CPLXSXP, ldvsl, want_vectors ? n : 1));
  SEXP Z_ = PROTECT(allocMatrix(CPLXSXP, ldvsr, want_vectors ? n : 1));
  La_INT sdim = 0;

  if (n > 0) {
    La_LGL* bwork = reinterpret_cast<La_LGL*>(
      R_alloc(static_cast<size_t>(n), sizeof(La_LGL))
    );
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(8 * n), sizeof(double))
    );
    La_INT info = 0;
    La_INT lwork = -1;
    Rcomplex work_query;
    F77_CALL(zgges)(&jobvsl, &jobvsr, &sort, selector, &n, COMPLEX(S_), &n,
                    COMPLEX(T_), &n, &sdim, COMPLEX(alpha_), COMPLEX(beta_),
                    COMPLEX(Q_), &ldvsl, COMPLEX(Z_), &ldvsr, &work_query,
                    &lwork, rwork, bwork, &info FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgges workspace query failed with info=%d", info);
    }
    lwork = static_cast<La_INT>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    std::memcpy(COMPLEX(S_), COMPLEX(A_),
                sizeof(Rcomplex) * static_cast<size_t>(n) * static_cast<size_t>(n));
    std::memcpy(COMPLEX(T_), COMPLEX(B_),
                sizeof(Rcomplex) * static_cast<size_t>(n) * static_cast<size_t>(n));
    F77_CALL(zgges)(&jobvsl, &jobvsr, &sort, selector, &n, COMPLEX(S_), &n,
                    COMPLEX(T_), &n, &sdim, COMPLEX(alpha_), COMPLEX(beta_),
                    COMPLEX(Q_), &ldvsl, COMPLEX(Z_), &ldvsr, COMPLEX(work_),
                    &lwork, rwork, bwork, &info FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgges failed with info=%d", info);
    }
    UNPROTECT(1);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_, 0, S_);
  SET_VECTOR_ELT(out_, 1, T_);
  SET_VECTOR_ELT(out_, 2, want_vectors ? Q_ : R_NilValue);
  SET_VECTOR_ELT(out_, 3, want_vectors ? Z_ : R_NilValue);
  SET_VECTOR_ELT(out_, 4, alpha_);
  SET_VECTOR_ELT(out_, 5, beta_);
  SET_VECTOR_ELT(out_, 6, ScalarInteger(sdim));
  SEXP names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names_, 0, mkChar("S"));
  SET_STRING_ELT(names_, 1, mkChar("T"));
  SET_STRING_ELT(names_, 2, mkChar("Q"));
  SET_STRING_ELT(names_, 3, mkChar("Z"));
  SET_STRING_ELT(names_, 4, mkChar("alpha"));
  SET_STRING_ELT(names_, 5, mkChar("beta"));
  SET_STRING_ELT(names_, 6, mkChar("sdim"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(8);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_is_symmetric(SEXP A_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    return ScalarLogical(FALSE);
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    return ScalarLogical(FALSE);
  }
  const int n = INTEGER(dimA)[0];
  const int p = INTEGER(dimA)[1];
  if (n != p) {
    return ScalarLogical(FALSE);
  }
  const double tol = asReal(tol_);
  const double* A = REAL(A_);
  // Relative to the largest entry (no floor at 1) so a rescaled
  // nonsymmetric matrix is never classified as symmetric.
  double scale = 0.0;
  for (int64_t i = 0; i < static_cast<int64_t>(n) * n; ++i) {
    // Non-finite entries cannot certify symmetry.
    if (!R_FINITE(A[i])) {
      return ScalarLogical(FALSE);
    }
    const double ai = fabs(A[i]);
    if (ai > scale) {
      scale = ai;
    }
  }
  const double threshold = (R_FINITE(tol) && tol >= 0.0 ? tol : sqrt(DBL_EPSILON)) * scale;
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < col; ++row) {
      const double a = A[row + static_cast<int64_t>(col) * n];
      const double b = A[col + static_cast<int64_t>(row) * n];
      if (!(fabs(a - b) <= threshold)) {
        return ScalarLogical(FALSE);
      }
    }
  }
  return ScalarLogical(TRUE);
  EIGENCORE_ENTRY_END
}

// One pass over a double matrix returning c(all_finite, is_symmetric), with
// the same semantics as eigencore_dense_is_symmetric: symmetric means square,
// every entry finite, and max |a_ij - a_ji| <= tol * max |a_ij| (relative to
// the largest entry, no floor at 1). The pair (a_ij, a_ji) is read in
// cache-sized tiles so the transposed access stays local; each entry is read
// exactly once. Non-finite entries are detected via x * 0 (NaN unless x is
// finite) accumulated into one sum, which keeps the inner loop branch-free.
extern "C" SEXP eigencore_dense_finite_symmetric(SEXP A_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("eigencore_dense_finite_symmetric expects a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  const double* A = REAL(A_);
  const R_xlen_t len = XLENGTH(A_);
  const double tol0 = asReal(tol_);
  const double tol = (R_FINITE(tol0) && tol0 >= 0.0) ? tol0 : sqrt(DBL_EPSILON);
  bool square = false;
  int64_t n = 0;
  if (dimA != R_NilValue && LENGTH(dimA) == 2) {
    n = INTEGER(dimA)[0];
    square = (n == INTEGER(dimA)[1]);
  }
  double nonfinite = 0.0;
  bool symmetric = false;
  if (!square) {
    for (R_xlen_t i = 0; i < len; ++i) {
      nonfinite += A[i] * 0.0;
    }
  } else {
    const int64_t B = 64;
    double scale = 0.0;
    double maxdiff = 0.0;
    for (int64_t jb = 0; jb < n; jb += B) {
      const int64_t jend = (jb + B < n) ? jb + B : n;
      for (int64_t ib = 0; ib <= jb; ib += B) {
        const int64_t iend = (ib + B < n) ? ib + B : n;
        for (int64_t j = jb; j < jend; ++j) {
          const int64_t ilim = (ib == jb) ? j : iend;
          const double* colj = A + j * n;
          for (int64_t i = ib; i < ilim; ++i) {
            const double a = colj[i];
            const double b = A[j + i * n];
            nonfinite += a * 0.0 + b * 0.0;
            const double fa = fabs(a);
            const double fb = fabs(b);
            const double d = fabs(a - b);
            scale = fa > scale ? fa : scale;
            scale = fb > scale ? fb : scale;
            maxdiff = d > maxdiff ? d : maxdiff;
          }
          if (ib == jb) {
            const double djj = colj[j];
            nonfinite += djj * 0.0;
            const double fd = fabs(djj);
            scale = fd > scale ? fd : scale;
          }
        }
      }
    }
    symmetric = !ISNAN(nonfinite) && maxdiff <= tol * scale;
  }
  SEXP out = PROTECT(allocVector(LGLSXP, 2));
  LOGICAL(out)[0] = ISNAN(nonfinite) ? FALSE : TRUE;
  LOGICAL(out)[1] = symmetric ? TRUE : FALSE;
  UNPROTECT(1);
  return out;
  EIGENCORE_ENTRY_END
}

// Selected dense symmetric eigenpairs via dsyevr RANGE = 'I': the k largest
// (target_kind 1, returned descending) or k smallest (target_kind 2,
// ascending) algebraic eigenvalues. With vectors = FALSE no eigenvectors are
// formed and `vectors` is NULL.
extern "C" SEXP eigencore_dense_symmetric_eigen_selected(SEXP A_, SEXP k_,
                                                         SEXP target_kind_,
                                                         SEXP vectors_flag_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  if (n != ncolA) {
    error("A must be square");
  }
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  int k = static_cast<int>(asInteger(k_));
  if (k == NA_INTEGER || k < 1) {
    error("k must be >= 1");
  }
  if (k > n) {
    k = n;
  }
  if (target_kind != 1 && target_kind != 2) {
    error("selected dense symmetric eigen supports only largest/smallest algebraic targets");
  }
  const bool want_vectors = vectors_flag(vectors_flag_);

  SEXP values_ = PROTECT(allocVector(REALSXP, k));
  SEXP vectors_ = PROTECT(want_vectors ? allocMatrix(REALSXP, n, k) : R_NilValue);
  if (n > 0) {
    double* values_work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    int il = 1;
    int iu = k;
    if (target_kind == 1) {
      il = n - k + 1;
      iu = n;
    }
    int m_found = 0;
    const int info = run_dsyevr(REAL(A_), n, want_vectors, il, iu, values_work,
                                want_vectors ? REAL(vectors_) : NULL, &m_found);
    if (info != 0 || m_found != k) {
      error("LAPACK dsyevr failed with info=%d, found=%d", info, m_found);
    }
    for (int col = 0; col < k; ++col) {
      REAL(values_)[col] = values_work[col];
    }
    if (target_kind == 1) {
      for (int left = 0, right = k - 1; left < right; ++left, --right) {
        const double tmp_value = REAL(values_)[left];
        REAL(values_)[left] = REAL(values_)[right];
        REAL(values_)[right] = tmp_value;
        if (!want_vectors) {
          continue;
        }
        for (int row = 0; row < n; ++row) {
          const int64_t lpos = row + static_cast<int64_t>(left) * n;
          const int64_t rpos = row + static_cast<int64_t>(right) * n;
          const double tmp_vec = REAL(vectors_)[lpos];
          REAL(vectors_)[lpos] = REAL(vectors_)[rpos];
          REAL(vectors_)[rpos] = tmp_vec;
        }
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_symmetric_eigen_dsyevx_selected(SEXP A_, SEXP k_, SEXP target_kind_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  if (n != ncolA) {
    error("A must be square");
  }
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  int k = static_cast<int>(asInteger(k_));
  if (k < 1) {
    error("k must be >= 1");
  }
  if (k > n) {
    k = n;
  }
  if (target_kind != 1 && target_kind != 2) {
    error("selected dense symmetric eigen supports only largest/smallest algebraic targets");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, k));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k));
  if (n > 0) {
    double* work_matrix = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n) * static_cast<size_t>(n), sizeof(double))
    );
    double* values_work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(5 * n), sizeof(int))
    );
    int* ifail = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(n), sizeof(int))
    );
    std::memcpy(work_matrix, REAL(A_),
                sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
    char jobz = 'V';
    char range = 'I';
    char uplo = 'U';
    double vl = 0.0;
    double vu = 0.0;
    double abstol = 0.0;
    int il = 1;
    int iu = k;
    if (target_kind == 1) {
      il = n - k + 1;
      iu = n;
    }
    int m_found = 0;
    int info = 0;
    int lwork = 8 * n;
    if (lwork < 1) {
      lwork = 1;
    }
    double* work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lwork), sizeof(double))
    );
    F77_CALL(dsyevx)(&jobz, &range, &uplo, &n, work_matrix, &n,
                     &vl, &vu, &il, &iu, &abstol,
                     &m_found, values_work, REAL(vectors_), &n,
                     work, &lwork, iwork, ifail, &info FCONE FCONE FCONE);
    if (info != 0 || m_found != k) {
      error("LAPACK dsyevx failed with info=%d, found=%d", info, m_found);
    }
    for (int col = 0; col < k; ++col) {
      REAL(values_)[col] = values_work[col];
    }
    if (target_kind == 1) {
      for (int left = 0, right = k - 1; left < right; ++left, --right) {
        const double tmp_value = REAL(values_)[left];
        REAL(values_)[left] = REAL(values_)[right];
        REAL(values_)[right] = tmp_value;
        for (int row = 0; row < n; ++row) {
          const int64_t lpos = row + static_cast<int64_t>(left) * n;
          const int64_t rpos = row + static_cast<int64_t>(right) * n;
          const double tmp_vec = REAL(vectors_)[lpos];
          REAL(vectors_)[lpos] = REAL(vectors_)[rpos];
          REAL(vectors_)[rpos] = tmp_vec;
        }
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_generalized_spd_eigen(SEXP A_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(B_)) {
    error("A and B must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("A and B must be square matrices with the same dimension");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(duplicate(A_));
  SEXP Bwork_ = PROTECT(duplicate(B_));
  if (n > 0) {
    int itype = 1;
    char uplo = 'U';
    int info = 0;
    F77_CALL(dpotrf)(&uplo, &n, REAL(Bwork_), &n, &info FCONE);
    if (info != 0) {
      error("LAPACK dpotrf failed for generalized SPD B with info=%d", info);
    }

    F77_CALL(dsygst)(&itype, &uplo, &n, REAL(vectors_), &n,
                     REAL(Bwork_), &n, &info FCONE);
    if (info != 0) {
      error("LAPACK dsygst failed with info=%d", info);
    }

    char jobz = 'V';
    int lwork = -1;
    double work_query = 0.0;
    F77_CALL(dsyev)(&jobz, &uplo, &n, REAL(vectors_), &n, REAL(values_),
                    &work_query, &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyev workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    F77_CALL(dsyev)(&jobz, &uplo, &n, REAL(vectors_), &n, REAL(values_),
                    REAL(work_), &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyev failed with info=%d", info);
    }
    UNPROTECT(1);

    char side = 'L';
    char transa = 'N';
    char diag = 'N';
    double one = 1.0;
    F77_CALL(dtrsm)(&side, &uplo, &transa, &diag, &n, &n, &one,
                    REAL(Bwork_), &n, REAL(vectors_), &n FCONE FCONE FCONE FCONE);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(5);
  return out_;
  EIGENCORE_ENTRY_END
}

// Thin dense real SVD. Uses divide and conquer (dgesdd, jobz = 'S', what
// base R's svd() uses); if dgesdd reports a convergence failure (info > 0)
// it retries with the QR-iteration driver dgesvd on a fresh copy of A.
// `driver` records which LAPACK routine produced the factors.
extern "C" SEXP eigencore_dense_svd(SEXP A_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int r = (m < n) ? m : n;

  SEXP d_ = PROTECT(allocVector(REALSXP, r));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, r));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, r));
  const char* driver = "dgesdd";

  if (r > 0) {
    const size_t mn = static_cast<size_t>(m) * static_cast<size_t>(n);
    double* work_matrix = reinterpret_cast<double*>(R_alloc(mn, sizeof(double)));
    double* vt = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(r) * static_cast<size_t>(n), sizeof(double))
    );
    std::memcpy(work_matrix, REAL(A_), sizeof(double) * mn);
    int lda = m;
    int ldu = m;
    int ldvt = r;
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    char jobz = 'S';
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(8) * static_cast<size_t>(r), sizeof(int))
    );
    F77_CALL(dgesdd)(&jobz, &m, &n, work_matrix, &lda, REAL(d_),
                     REAL(u_), &ldu, vt, &ldvt,
                     &work_query, &lwork, iwork, &info FCONE);
    if (info != 0) {
      error("LAPACK dgesdd workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    if (lwork < 1) {
      lwork = 1;
    }
    double* work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lwork), sizeof(double))
    );
    F77_CALL(dgesdd)(&jobz, &m, &n, work_matrix, &lda, REAL(d_),
                     REAL(u_), &ldu, vt, &ldvt,
                     work, &lwork, iwork, &info FCONE);
    if (info < 0) {
      error("LAPACK dgesdd failed with info=%d", info);
    }
    if (info > 0) {
      // dbdsdc did not converge: fall back to bidiagonal QR iteration.
      driver = "dgesvd";
      std::memcpy(work_matrix, REAL(A_), sizeof(double) * mn);
      char jobu = 'S';
      char jobvt = 'S';
      lwork = -1;
      info = 0;
      F77_CALL(dgesvd)(&jobu, &jobvt, &m, &n, work_matrix, &lda,
                       REAL(d_), REAL(u_), &ldu, vt, &ldvt,
                       &work_query, &lwork, &info FCONE FCONE);
      if (info != 0) {
        error("LAPACK dgesvd workspace query failed with info=%d", info);
      }
      lwork = static_cast<int>(work_query);
      if (lwork < 1) {
        lwork = 1;
      }
      work = reinterpret_cast<double*>(
        R_alloc(static_cast<size_t>(lwork), sizeof(double))
      );
      F77_CALL(dgesvd)(&jobu, &jobvt, &m, &n, work_matrix, &lda,
                       REAL(d_), REAL(u_), &ldu, vt, &ldvt,
                       work, &lwork, &info FCONE FCONE);
      if (info != 0) {
        error("LAPACK dgesvd failed with info=%d", info);
      }
    }

    double* v = REAL(v_);
    for (int col = 0; col < r; ++col) {
      for (int row = 0; row < n; ++row) {
        v[row + static_cast<R_xlen_t>(col) * n] = vt[col + static_cast<R_xlen_t>(row) * r];
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, mkString(driver));
  SEXP names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("driver"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(5);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_svd(SEXP A_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_)) {
    error("A must be a complex matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int r = (m < n) ? m : n;

  SEXP d_ = PROTECT(allocVector(REALSXP, r));
  SEXP u_ = PROTECT(allocMatrix(CPLXSXP, m, r));
  SEXP vt_ = PROTECT(allocMatrix(CPLXSXP, r, n));
  SEXP v_ = PROTECT(allocMatrix(CPLXSXP, n, r));
  SEXP work_matrix_ = PROTECT(duplicate(A_));

  if (r > 0) {
    char jobu = 'S';
    char jobvt = 'S';
    int lda = m;
    int ldu = m;
    int ldvt = r;
    int info = 0;
    int lwork = -1;
    Rcomplex work_query;
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(5 * r), sizeof(double))
    );
    F77_CALL(zgesvd)(&jobu, &jobvt, &m, &n, COMPLEX(work_matrix_), &lda,
                     REAL(d_), COMPLEX(u_), &ldu, COMPLEX(vt_), &ldvt,
                     &work_query, &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgesvd workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    work_matrix_ = PROTECT(duplicate(A_));
    F77_CALL(zgesvd)(&jobu, &jobvt, &m, &n, COMPLEX(work_matrix_), &lda,
                     REAL(d_), COMPLEX(u_), &ldu, COMPLEX(vt_), &ldvt,
                     COMPLEX(work_), &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zgesvd failed with info=%d", info);
    }
    UNPROTECT(2);

    for (int col = 0; col < r; ++col) {
      for (int row = 0; row < n; ++row) {
        const Rcomplex z = COMPLEX(vt_)[col + static_cast<int64_t>(row) * r];
        Rcomplex z_conj;
        z_conj.r = z.r;
        z_conj.i = -z.i;
        COMPLEX(v_)[row + static_cast<int64_t>(col) * n] = z_conj;
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(7);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_generalized_svd(SEXP A_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(B_)) {
    error("A and B must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue) {
    error("A and B must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int p = INTEGER(dimB)[0];
  const int nB = INTEGER(dimB)[1];
  if (m <= 0 || n <= 0 || p <= 0) {
    error("A and B must have positive dimensions");
  }
  if (n != nB) {
    error("A and B must have the same number of columns");
  }

  SEXP Awork_ = PROTECT(duplicate(A_));
  SEXP Bwork_ = PROTECT(duplicate(B_));
  SEXP k_ = PROTECT(ScalarInteger(0));
  SEXP l_ = PROTECT(ScalarInteger(0));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, n));
  SEXP beta_ = PROTECT(allocVector(REALSXP, n));
  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, m));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, n));

  char jobu = 'U';
  char jobv = 'V';
  char jobq = 'Q';
  int k = 0;
  int l = 0;
  int lda = m;
  int ldb = p;
  int ldu = m;
  int ldv = p;
  int ldq = n;
  int info = 0;
  int lwork = 3 * n;
  if (m > lwork) {
    lwork = m;
  }
  if (p > lwork) {
    lwork = p;
  }
  lwork += n;
  if (lwork < 1) {
    lwork = 1;
  }
  SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
  int* iwork = reinterpret_cast<int*>(
    R_alloc(static_cast<size_t>(n), sizeof(int))
  );

  F77_CALL(dggsvd)(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l,
                   REAL(Awork_), &lda, REAL(Bwork_), &ldb,
                   REAL(alpha_), REAL(beta_), REAL(U_), &ldu,
                   REAL(V_), &ldv, REAL(Q_), &ldq, REAL(work_),
                   iwork, &info FCONE FCONE FCONE);
  if (info != 0) {
    error("LAPACK dggsvd failed with info=%d", info);
  }
  UNPROTECT(1);
  INTEGER(k_)[0] = k;
  INTEGER(l_)[0] = l;

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, Awork_);
  SET_VECTOR_ELT(out_, 1, Bwork_);
  SET_VECTOR_ELT(out_, 2, k_);
  SET_VECTOR_ELT(out_, 3, l_);
  SET_VECTOR_ELT(out_, 4, alpha_);
  SET_VECTOR_ELT(out_, 5, beta_);
  SET_VECTOR_ELT(out_, 6, U_);
  SET_VECTOR_ELT(out_, 7, V_);
  SET_VECTOR_ELT(out_, 8, Q_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("A_factor"));
  SET_STRING_ELT(names_, 1, mkChar("B_factor"));
  SET_STRING_ELT(names_, 2, mkChar("k"));
  SET_STRING_ELT(names_, 3, mkChar("l"));
  SET_STRING_ELT(names_, 4, mkChar("alpha"));
  SET_STRING_ELT(names_, 5, mkChar("beta"));
  SET_STRING_ELT(names_, 6, mkChar("U"));
  SET_STRING_ELT(names_, 7, mkChar("V"));
  SET_STRING_ELT(names_, 8, mkChar("Q"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(11);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_generalized_svd(SEXP A_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  error("native complex GSVD requires a complex LAPACK GSVD driver, "
        "which this R LAPACK interface does not export");
  EIGENCORE_ENTRY_END
}

// Full symmetric tridiagonal eigendecomposition (ascending values) via the
// MRRR driver dstevr (RANGE = 'A'); falls back to implicit QL/QR (dstev) if
// dstevr reports an internal failure.
extern "C" SEXP eigencore_tridiagonal_eigen(SEXP alpha_, SEXP beta_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(alpha_) || !isReal(beta_)) {
    error("alpha and beta must be double vectors");
  }
  const int n = LENGTH(alpha_);
  if (LENGTH(beta_) < ((n > 0) ? n - 1 : 0)) {
    error("beta must have length at least length(alpha) - 1");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, n));

  if (n == 1) {
    REAL(values_)[0] = REAL(alpha_)[0];
    REAL(vectors_)[0] = 1.0;
  } else if (n > 1) {
    double* diag = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    // Length n: dstevr may use E(n) as workspace.
    double* offdiag = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    int* isuppz = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(2 * n), sizeof(int))
    );
    for (int i = 0; i < n; ++i) {
      diag[i] = REAL(alpha_)[i];
    }
    for (int i = 0; i < n - 1; ++i) {
      offdiag[i] = REAL(beta_)[i];
    }
    offdiag[n - 1] = 0.0;
    char jobz = 'V';
    char range = 'A';
    double vl = 0.0;
    double vu = 0.0;
    int il = 1;
    int iu = n;
    double abstol = 0.0;
    int m_found = 0;
    int info = 0;
    int ldz = n;
    int n_la = n;
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    F77_CALL(dstevr)(&jobz, &range, &n_la, diag, offdiag,
                     &vl, &vu, &il, &iu, &abstol, &m_found,
                     REAL(values_), REAL(vectors_), &ldz, isuppz,
                     &work_query, &lwork, &iwork_query, &liwork,
                     &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dstevr workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    if (lwork < 20 * n) {
      lwork = 20 * n;
    }
    if (liwork < 10 * n) {
      liwork = 10 * n;
    }
    double* work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lwork), sizeof(double))
    );
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(liwork), sizeof(int))
    );
    F77_CALL(dstevr)(&jobz, &range, &n_la, diag, offdiag,
                     &vl, &vu, &il, &iu, &abstol, &m_found,
                     REAL(values_), REAL(vectors_), &ldz, isuppz,
                     work, &lwork, iwork, &liwork,
                     &info FCONE FCONE);
    if (info < 0) {
      error("LAPACK dstevr failed with info=%d", info);
    }
    if (info > 0 || m_found != n) {
      for (int i = 0; i < n; ++i) {
        REAL(values_)[i] = REAL(alpha_)[i];
      }
      for (int i = 0; i < n - 1; ++i) {
        offdiag[i] = REAL(beta_)[i];
      }
      double* qr_work = reinterpret_cast<double*>(
        R_alloc(static_cast<size_t>(2 * n - 2), sizeof(double))
      );
      char jobz_v = 'V';
      info = 0;
      F77_CALL(dstev)(&jobz_v, &n, REAL(values_), offdiag,
                      REAL(vectors_), &n, qr_work, &info FCONE);
      if (info != 0) {
        error("LAPACK dstev failed with info=%d", info);
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_tridiagonal_eigen_selected(SEXP alpha_, SEXP beta_,
                                                     SEXP k_, SEXP target_kind_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(alpha_) || !isReal(beta_)) {
    error("alpha and beta must be double vectors");
  }
  const int n = LENGTH(alpha_);
  if (n < 1) {
    error("alpha must have positive length");
  }
  if (LENGTH(beta_) < n - 1) {
    error("beta must have length at least length(alpha) - 1");
  }
  int k = static_cast<int>(asInteger(k_));
  if (k < 1) {
    error("k must be >= 1");
  }
  if (k > n) {
    k = n;
  }
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  if (target_kind != 1 && target_kind != 2) {
    error("selected tridiagonal eigen supports only largest/smallest algebraic targets");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, k));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k));
  if (n == 1) {
    REAL(values_)[0] = REAL(alpha_)[0];
    REAL(vectors_)[0] = 1.0;
  } else {
    double* diag = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    // Length n (not n - 1): dstevr documents E as workspace it may overwrite
    // through index n - 1 in some implementations.
    double* offdiag = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    double* values_work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    int* isuppz = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(2 * k), sizeof(int))
    );
    for (int i = 0; i < n; ++i) {
      diag[i] = REAL(alpha_)[i];
    }
    for (int i = 0; i < n - 1; ++i) {
      offdiag[i] = REAL(beta_)[i];
    }
    offdiag[n - 1] = 0.0;

    char jobz = 'V';
    char range = 'I';
    int il = 1;
    int iu = k;
    if (target_kind == 1) {
      il = n - k + 1;
      iu = n;
    }
    double vl = 0.0;
    double vu = 0.0;
    double abstol = 0.0;
    int m_found = 0;
    int info = 0;
    int ldz = n;
    int n_la = n;  // dstevr's header declaration is not const-qualified
    // MRRR driver (dstevr) instead of bisection + inverse iteration (dstevx):
    // far faster for selected eigenpairs when eigenvalues cluster, e.g. the
    // smallest eigenvalues of a Laplacian.
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    F77_CALL(dstevr)(&jobz, &range, &n_la, diag, offdiag,
                     &vl, &vu, &il, &iu, &abstol, &m_found,
                     values_work, REAL(vectors_), &ldz, isuppz,
                     &work_query, &lwork, &iwork_query, &liwork,
                     &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dstevr workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    liwork = iwork_query;
    double* work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lwork), sizeof(double))
    );
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(liwork), sizeof(int))
    );
    F77_CALL(dstevr)(&jobz, &range, &n_la, diag, offdiag,
                     &vl, &vu, &il, &iu, &abstol, &m_found,
                     values_work, REAL(vectors_), &ldz, isuppz,
                     work, &lwork, iwork, &liwork,
                     &info FCONE FCONE);
    if (info != 0 || m_found != k) {
      error("LAPACK dstevr failed with info=%d, found=%d", info, m_found);
    }
    for (int col = 0; col < k; ++col) {
      REAL(values_)[col] = values_work[col];
    }
    if (target_kind == 1) {
      for (int left = 0, right = k - 1; left < right; ++left, --right) {
        const double tmp_value = REAL(values_)[left];
        REAL(values_)[left] = REAL(values_)[right];
        REAL(values_)[right] = tmp_value;
        for (int row = 0; row < n; ++row) {
          const int64_t lpos = row + static_cast<int64_t>(left) * n;
          const int64_t rpos = row + static_cast<int64_t>(right) * n;
          const double tmp_vec = REAL(vectors_)[lpos];
          REAL(vectors_)[lpos] = REAL(vectors_)[rpos];
          REAL(vectors_)[rpos] = tmp_vec;
        }
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
  EIGENCORE_ENTRY_END
}

// SVD of the n x n upper bidiagonal matrix with diagonal alpha and
// superdiagonal beta, computed directly on (d, e) by the bidiagonal
// divide-and-conquer routine dbdsdc (COMPQ = 'I'), with the implicit
// zero-shift QR routine dbdsqr as fallback when dbdsdc fails to converge.
// Returns d (descending), u and v (both n x n) with B = U diag(d) V^T.
extern "C" SEXP eigencore_bidiagonal_svd(SEXP alpha_, SEXP beta_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(alpha_) || !isReal(beta_)) {
    error("alpha and beta must be double vectors");
  }
  const int n = LENGTH(alpha_);
  if (LENGTH(beta_) < ((n > 0) ? n - 1 : 0)) {
    error("beta must have length at least length(alpha) - 1");
  }

  SEXP d_ = PROTECT(allocVector(REALSXP, n));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, n));
  if (n > 0) {
    const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
    double* d = REAL(d_);
    double* u = REAL(u_);
    double* e = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    double* vt = reinterpret_cast<double*>(R_alloc(nn, sizeof(double)));
    for (int i = 0; i < n; ++i) {
      d[i] = REAL(alpha_)[i];
      e[i] = (i < n - 1) ? REAL(beta_)[i] : 0.0;
    }

    char uplo = 'U';
    char compq = 'I';
    int n_la = n;
    int ldu = n;
    int ldvt = n;
    double q_dummy = 0.0;
    int iq_dummy = 0;
    const size_t lwork = 3 * nn + 4 * static_cast<size_t>(n);
    double* work = reinterpret_cast<double*>(R_alloc(lwork, sizeof(double)));
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(8) * static_cast<size_t>(n), sizeof(int))
    );
    int info = 0;
    F77_CALL(dbdsdc)(&uplo, &compq, &n_la, d, e, u, &ldu, vt, &ldvt,
                     &q_dummy, &iq_dummy, work, iwork, &info FCONE FCONE);
    if (info < 0) {
      error("LAPACK dbdsdc failed with info=%d", info);
    }
    if (info > 0) {
      // Divide and conquer failed: restart from the input with dbdsqr.
      for (int i = 0; i < n; ++i) {
        d[i] = REAL(alpha_)[i];
        e[i] = (i < n - 1) ? REAL(beta_)[i] : 0.0;
      }
      std::memset(u, 0, sizeof(double) * nn);
      std::memset(vt, 0, sizeof(double) * nn);
      for (int i = 0; i < n; ++i) {
        u[i + static_cast<size_t>(i) * n] = 1.0;
        vt[i + static_cast<size_t>(i) * n] = 1.0;
      }
      int ncvt = n;
      int nru = n;
      int ncc = 0;
      int ldc = 1;
      double c_dummy = 0.0;
      info = 0;
      F77_CALL(dbdsqr)(&uplo, &n_la, &ncvt, &nru, &ncc, d, e, vt, &ldvt,
                       u, &ldu, &c_dummy, &ldc, work, &info FCONE);
      if (info != 0) {
        error("LAPACK dbdsqr failed with info=%d", info);
      }
    }

    double* v = REAL(v_);
    for (int col = 0; col < n; ++col) {
      for (int row = 0; row < n; ++row) {
        v[row + static_cast<R_xlen_t>(col) * n] = vt[col + static_cast<R_xlen_t>(row) * n];
      }
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(5);
  return out_;
  EIGENCORE_ENTRY_END
}

// Solve T X = B for a general tridiagonal T (subdiagonal `lower`, diagonal
// `diag`, superdiagonal `upper`) by LU with partial pivoting (dgttrf/dgttrs).
// Shifted tridiagonals are typically indefinite, where the unpivoted Thomas
// algorithm is unstable. A matrix that is exactly singular (dgttrf info > 0)
// or singular to working precision (1-norm reciprocal condition estimate
// from dgtcon below DBL_EPSILON) is rejected.
extern "C" SEXP eigencore_tridiagonal_solve(SEXP lower_, SEXP diag_,
                                            SEXP upper_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(lower_) || !isReal(diag_) || !isReal(upper_) || !isReal(B_)) {
    error("lower, diag, upper, and B must be double");
  }
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimB == R_NilValue) {
    error("B must be a matrix");
  }
  const int n = LENGTH(diag_);
  const int nrhs = INTEGER(dimB)[1];
  if (INTEGER(dimB)[0] != n || LENGTH(lower_) != (n > 0 ? n - 1 : 0) ||
      LENGTH(upper_) != (n > 0 ? n - 1 : 0)) {
    error("non-conformable tridiagonal solve inputs");
  }

  SEXP out_ = PROTECT(duplicate(B_));
  if (n == 0 || nrhs == 0) {
    UNPROTECT(1);
    return out_;
  }
  const size_t nm1 = static_cast<size_t>(n > 1 ? n - 1 : 1);
  double* dl = reinterpret_cast<double*>(R_alloc(nm1, sizeof(double)));
  double* du = reinterpret_cast<double*>(R_alloc(nm1, sizeof(double)));
  double* du2 = reinterpret_cast<double*>(R_alloc(nm1, sizeof(double)));
  double* d = reinterpret_cast<double*>(R_alloc(static_cast<size_t>(n), sizeof(double)));
  int* ipiv = reinterpret_cast<int*>(R_alloc(static_cast<size_t>(n), sizeof(int)));
  for (int i = 0; i < n - 1; ++i) {
    dl[i] = REAL(lower_)[i];
    du[i] = REAL(upper_)[i];
  }
  std::memcpy(d, REAL(diag_), sizeof(double) * static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(d[i]) || (i < n - 1 && (!R_FINITE(dl[i]) || !R_FINITE(du[i])))) {
      error("tridiagonal solve: matrix has non-finite entries");
    }
  }

  char norm = '1';
  const double anorm = F77_CALL(dlangt)(&norm, &n, dl, d, du FCONE);
  int info = 0;
  F77_CALL(dgttrf)(&n, dl, d, du, du2, ipiv, &info);
  if (info < 0) {
    error("LAPACK dgttrf failed with info=%d", info);
  }
  if (info > 0 || anorm == 0.0) {
    error("tridiagonal solve: matrix is exactly singular (zero pivot %d)", info);
  }
  double rcond = 0.0;
  double* con_work = reinterpret_cast<double*>(
    R_alloc(static_cast<size_t>(2 * n), sizeof(double))
  );
  int* con_iwork = reinterpret_cast<int*>(R_alloc(static_cast<size_t>(n), sizeof(int)));
  F77_CALL(dgtcon)(&norm, &n, dl, d, du, du2, ipiv, &anorm, &rcond,
                   con_work, con_iwork, &info FCONE);
  if (info != 0) {
    error("LAPACK dgtcon failed with info=%d", info);
  }
  if (!(rcond >= DBL_EPSILON)) {
    error("tridiagonal solve: matrix is singular to working precision (rcond=%g)",
          rcond);
  }
  char trans = 'N';
  int ldb = n;
  F77_CALL(dgttrs)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, REAL(out_), &ldb,
                   &info FCONE);
  if (info != 0) {
    error("LAPACK dgttrs failed with info=%d", info);
  }

  UNPROTECT(1);
  return out_;
  EIGENCORE_ENTRY_END
}

// ---------------------------------------------------------------------------
// Sylvester inertia helpers (tranche 5, capability gap 4).
//
// Each helper factors a shifted symmetric matrix M = A - sigma B and returns
// a named numeric vector
//   neg, zero, pos          inertia of the computed block-diagonal factor D
//   min_abs_pivot           smallest |eigenvalue| over the 1x1 / 2x2 pivots
//   max_abs_pivot           largest  |eigenvalue| over the pivots
//   max_abs_multiplier      largest |L_ij| (i > j, outside 2x2 pivot blocks)
//   growth                  || |L| |D| |L'| ||_inf (sparse LDL'); dense and
//                           tridiagonal report max_abs_pivot *
//                           max(1, max_abs_multiplier)^2 as a proxy
//   norm1                   ||M||_1 (= ||M||_inf, symmetric); NA when the
//                           helper does not see M (sparse factor diagnostics)
//   info                    LAPACK info (> 0: an exact zero pivot)
//   two_by_two              number of 2x2 pivot blocks
// By Sylvester's law of inertia M = P L D L' P' is congruent to D, so the
// counts are those of the eigenvalues of A - sigma B, i.e. (B positive
// definite) of the pencil's eigenvalues below / at / above sigma. They are
// exact for a backward-perturbed M; the R layer (R/inertia.R) judges
// reliability from min_abs_pivot and growth relative to the matrix scale.
// ---------------------------------------------------------------------------

namespace {

struct InertiaTally {
  double neg = 0.0;
  double zero = 0.0;
  double pos = 0.0;
  double min_abs = R_PosInf;
  double max_abs = 0.0;
  double max_mult = 0.0;
  double growth = 0.0;
  double norm1 = 0.0;
  double info = 0.0;
  double two_by_two = 0.0;

  void add_pivot(double value) {
    if (value < 0.0) {
      neg += 1.0;
    } else if (value > 0.0) {
      pos += 1.0;
    } else {
      zero += 1.0;
    }
    const double a = fabs(value);
    if (a < min_abs) min_abs = a;
    if (a > max_abs) max_abs = a;
  }
};

SEXP inertia_tally_sexp(const InertiaTally& t) {
  const char* names[] = {"neg", "zero", "pos", "min_abs_pivot", "max_abs_pivot",
                         "max_abs_multiplier", "growth", "norm1", "info",
                         "two_by_two"};
  const double values[] = {t.neg, t.zero, t.pos,
                           R_FINITE(t.min_abs) ? t.min_abs : 0.0, t.max_abs,
                           t.max_mult, t.growth, t.norm1, t.info, t.two_by_two};
  const int len = 10;
  SEXP out_ = PROTECT(allocVector(REALSXP, len));
  SEXP names_ = PROTECT(allocVector(STRSXP, len));
  for (int i = 0; i < len; ++i) {
    REAL(out_)[i] = values[i];
    SET_STRING_ELT(names_, i, mkChar(names[i]));
  }
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(2);
  return out_;
}

}  // namespace

// Dense real symmetric inertia of A - sigma B (B = NULL: identity) by
// Bunch-Kaufman LDL' (dsytrf on the lower triangle). Only the lower
// triangles of A and B are read.
extern "C" SEXP eigencore_dense_symmetric_inertia(SEXP A_, SEXP sigma_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue || INTEGER(dimA)[0] != INTEGER(dimA)[1]) {
    error("A must be a square matrix");
  }
  const int n = INTEGER(dimA)[0];
  const double sigma = asReal(sigma_);
  if (!R_FINITE(sigma)) {
    error("sigma must be finite");
  }
  const bool has_B = B_ != R_NilValue;
  if (has_B) {
    SEXP dimB = getAttrib(B_, R_DimSymbol);
    if (!isReal(B_) || dimB == R_NilValue || INTEGER(dimB)[0] != n ||
        INTEGER(dimB)[1] != n) {
      error("B must be a double matrix with the dimensions of A");
    }
  }
  InertiaTally tally;
  if (n == 0) {
    return inertia_tally_sexp(tally);
  }
  const size_t nn = static_cast<size_t>(n) * static_cast<size_t>(n);
  std::vector<double> M(nn, 0.0);
  const double* A = REAL(A_);
  const double* B = has_B ? REAL(B_) : nullptr;
  std::vector<double> colsum(static_cast<size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    for (int i = j; i < n; ++i) {
      const size_t ij = static_cast<size_t>(i) + static_cast<size_t>(j) * n;
      const double v = A[ij] - (has_B ? sigma * B[ij] : (i == j ? sigma : 0.0));
      if (!R_FINITE(v)) {
        error("A - sigma B has non-finite entries");
      }
      M[ij] = v;
      colsum[static_cast<size_t>(j)] += fabs(v);
      if (i != j) {
        colsum[static_cast<size_t>(i)] += fabs(v);
      }
    }
  }
  for (int j = 0; j < n; ++j) {
    if (colsum[static_cast<size_t>(j)] > tally.norm1) {
      tally.norm1 = colsum[static_cast<size_t>(j)];
    }
  }
  std::vector<int> ipiv(static_cast<size_t>(n), 0);
  char uplo = 'L';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsytrf)(&uplo, &n, M.data(), &n, ipiv.data(), &work_query, &lwork,
                   &info FCONE);
  if (info < 0) {
    error("LAPACK dsytrf workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  if (lwork < n) lwork = n;
  std::vector<double> work(static_cast<size_t>(lwork));
  F77_CALL(dsytrf)(&uplo, &n, M.data(), &n, ipiv.data(), work.data(), &lwork,
                   &info FCONE);
  if (info < 0) {
    error("LAPACK dsytrf failed with info=%d", info);
  }
  tally.info = static_cast<double>(info);
  int k = 0;
  while (k < n) {
    const size_t kk = static_cast<size_t>(k) + static_cast<size_t>(k) * n;
    if (ipiv[static_cast<size_t>(k)] > 0) {
      tally.add_pivot(M[kk]);
      for (int i = k + 1; i < n; ++i) {
        const double l = fabs(M[static_cast<size_t>(i) + static_cast<size_t>(k) * n]);
        if (l > tally.max_mult) tally.max_mult = l;
      }
      k += 1;
    } else {
      if (k + 1 >= n) {
        error("LAPACK dsytrf returned an incomplete 2x2 pivot block");
      }
      // 2x2 pivot [a b; b c]: its two eigenvalues, the small one from the
      // determinant to avoid cancellation.
      const double a = M[kk];
      const double b = M[kk + 1];
      const double c = M[kk + 1 + static_cast<size_t>(n)];
      const double mean = 0.5 * (a + c);
      const double rad = hypot(0.5 * (a - c), b);
      const double big = mean >= 0.0 ? mean + rad : mean - rad;
      const double det = a * c - b * b;
      const double small = big != 0.0 ? det / big : 0.0;
      tally.add_pivot(big);
      tally.add_pivot(small);
      tally.two_by_two += 1.0;
      for (int col = k; col <= k + 1; ++col) {
        for (int i = k + 2; i < n; ++i) {
          const double l = fabs(M[static_cast<size_t>(i) + static_cast<size_t>(col) * n]);
          if (l > tally.max_mult) tally.max_mult = l;
        }
      }
      k += 2;
    }
  }
  const double m1 = tally.max_mult > 1.0 ? tally.max_mult : 1.0;
  tally.growth = tally.max_abs * m1 * m1;
  return inertia_tally_sexp(tally);
  EIGENCORE_ENTRY_END
}

// Sturm count for the symmetric tridiagonal T - sigma B (B diagonal; NULL =
// identity) from the pivots of the unpivoted recurrence
//   q_1 = d_1 - sigma b_1,  q_i = (d_i - sigma b_i) - e_{i-1}^2 / q_{i-1}.
// The count is backward stable (it is exact for a componentwise relative
// perturbation of T; Kahan, Demmel). A pivot with |q| <= pivmin is replaced
// by -pivmin as in LAPACK dlaebz and reported as a zero pivot (not in neg).
extern "C" SEXP eigencore_tridiagonal_inertia(SEXP d_, SEXP e_, SEXP sigma_, SEXP b_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(d_) || !isReal(e_)) {
    error("tridiagonal inertia: d and e must be double vectors");
  }
  const R_xlen_t n = XLENGTH(d_);
  if (n > 0 && XLENGTH(e_) != n - 1) {
    error("tridiagonal inertia: e must have length n - 1");
  }
  const bool has_b = b_ != R_NilValue;
  if (has_b && (!isReal(b_) || XLENGTH(b_) != n)) {
    error("tridiagonal inertia: b must be a double vector of length n");
  }
  const double sigma = asReal(sigma_);
  if (!R_FINITE(sigma)) {
    error("sigma must be finite");
  }
  const double* d = REAL(d_);
  const double* e = REAL(e_);
  const double* b = has_b ? REAL(b_) : nullptr;
  InertiaTally tally;
  double emax2 = 0.0;
  for (R_xlen_t i = 0; i + 1 < n; ++i) {
    if (!R_FINITE(e[i])) error("tridiagonal inertia: non-finite off-diagonal");
    if (e[i] * e[i] > emax2) emax2 = e[i] * e[i];
  }
  const double pivmin = DBL_MIN * (emax2 > 1.0 ? emax2 : 1.0);
  double q = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double di = d[i] - sigma * (has_b ? b[i] : 1.0);
    if (!R_FINITE(di)) error("tridiagonal inertia: non-finite diagonal");
    const double left = i > 0 ? fabs(e[i - 1]) : 0.0;
    const double right = i + 1 < n ? fabs(e[i]) : 0.0;
    const double rowsum = fabs(di) + left + right;
    if (rowsum > tally.norm1) tally.norm1 = rowsum;
    q = i == 0 ? di : di - (e[i - 1] * e[i - 1]) / q;
    const double aq = fabs(q);
    if (aq < tally.min_abs) tally.min_abs = aq;
    if (aq > tally.max_abs) tally.max_abs = aq;
    if (aq <= pivmin) {
      tally.zero += 1.0;
      q = -pivmin;
    } else if (q < 0.0) {
      tally.neg += 1.0;
    } else {
      tally.pos += 1.0;
    }
  }
  tally.growth = tally.max_abs;
  return inertia_tally_sexp(tally);
  EIGENCORE_ENTRY_END
}

// Diagnostics of a CHOLMOD simplicial LDL' factor (Matrix dCHMsimpl slots
// p, i, x, nz): column j holds D_jj at x[p[j]] (row j) followed by the
// strictly lower entries of the unit lower triangular L. Returns the inertia
// of D, the pivot extremes, max |L_ij| and growth = || |L| |D| |L'| ||_inf,
// whose ratio to ||A - sigma B||_inf bounds the componentwise backward error
// of the unpivoted factorisation (up to a modest multiple of eps).
extern "C" SEXP eigencore_simplicial_ldl_diagnostics(SEXP p_, SEXP i_, SEXP x_,
                                                      SEXP nz_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(p_) || !isInteger(i_) || !isReal(x_) || !isInteger(nz_)) {
    error("simplicial LDL' diagnostics: unexpected factor slot types");
  }
  const R_xlen_t n = XLENGTH(nz_);
  if (XLENGTH(p_) < n + 1) {
    error("simplicial LDL' diagnostics: p must have length n + 1");
  }
  const int* p = INTEGER(p_);
  const int* ri = INTEGER(i_);
  const double* x = REAL(x_);
  const int* nz = INTEGER(nz_);
  const R_xlen_t len = XLENGTH(x_);
  const R_xlen_t ilen = XLENGTH(i_);
  InertiaTally tally;
  std::vector<double> y(static_cast<size_t>(n), 1.0);  // |L'| 1 (unit diagonal)
  for (R_xlen_t j = 0; j < n; ++j) {
    const R_xlen_t start = p[j];
    const R_xlen_t count = nz[j];
    if (count < 1 || start < 0 || start + count > len || start + count > ilen ||
        ri[start] != j) {
      error("simplicial LDL' diagnostics: column %ld does not start with its diagonal",
            static_cast<long>(j));
    }
    const double dj = x[start];
    if (!R_FINITE(dj)) {
      error("simplicial LDL' diagnostics: non-finite pivot in column %ld",
            static_cast<long>(j));
    }
    tally.add_pivot(dj);
    for (R_xlen_t t = start + 1; t < start + count; ++t) {
      const double l = fabs(x[t]);
      if (!R_FINITE(l)) {
        error("simplicial LDL' diagnostics: non-finite multiplier");
      }
      if (ri[t] <= j || ri[t] >= n) {
        error("simplicial LDL' diagnostics: malformed column %ld", static_cast<long>(j));
      }
      if (l > tally.max_mult) tally.max_mult = l;
      y[static_cast<size_t>(j)] += l;
    }
  }
  // z = |D| y, w = |L| z (unit diagonal plus the strictly lower part).
  std::vector<double> w(static_cast<size_t>(n), 0.0);
  for (R_xlen_t j = 0; j < n; ++j) {
    const R_xlen_t start = p[j];
    const double zj = fabs(x[start]) * y[static_cast<size_t>(j)];
    w[static_cast<size_t>(j)] += zj;
    for (R_xlen_t t = start + 1; t < start + nz[j]; ++t) {
      w[static_cast<size_t>(ri[t])] += fabs(x[t]) * zj;
    }
  }
  for (R_xlen_t j = 0; j < n; ++j) {
    if (w[static_cast<size_t>(j)] > tally.growth) tally.growth = w[static_cast<size_t>(j)];
  }
  tally.norm1 = NA_REAL;
  return inertia_tally_sexp(tally);
  EIGENCORE_ENTRY_END
}
