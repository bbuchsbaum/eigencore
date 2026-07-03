#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

extern "C" {
void F77_NAME(zhegv)(const La_INT* itype, const char* jobz, const char* uplo,
                     La_INT* n, Rcomplex* a, La_INT* lda, Rcomplex* b,
                     La_INT* ldb, double* w, Rcomplex* work, La_INT* lwork,
                     double* rwork, La_INT* info FCLEN FCLEN);
}

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

extern "C" SEXP eigencore_dense_symmetric_eigen(SEXP A_) {
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
    double work_query = 0.0;
    F77_CALL(dsyev)(&jobz, &uplo, &n, REAL(vectors_), &n, REAL(values_),
                    &work_query, &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyev workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    F77_CALL(dsyev)(&jobz, &uplo, &n, REAL(vectors_), &n, REAL(values_),
                    REAL(work_), &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyev failed with info=%d", info);
    }
    UNPROTECT(1);
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
}

extern "C" SEXP eigencore_dense_symmetric_eigen_dsyevd(SEXP A_) {
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
    if (lwork < 1 + 6 * n + 2 * n * n) {
      lwork = 1 + 6 * n + 2 * n * n;
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
}

extern "C" SEXP eigencore_dense_complex_hermitian_eigen(SEXP A_) {
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

  SEXP values_ = PROTECT(allocVector(REALSXP, n));
  SEXP vectors_ = PROTECT(duplicate(A_));
  if (n > 0) {
    char jobz = 'V';
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
  SET_VECTOR_ELT(out_, 1, vectors_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_dense_complex_general_eigen(SEXP A_) {
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
}

extern "C" SEXP eigencore_dense_complex_generalized_hpd_eigen(SEXP A_, SEXP B_) {
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
    int itype = 1;
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    int lwork = -1;
    Rcomplex work_query;
    const int lrwork = (3 * n - 2 > 1) ? (3 * n - 2) : 1;
    double* rwork = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(lrwork), sizeof(double))
    );
    F77_CALL(zhegv)(&itype, &jobz, &uplo, &n, COMPLEX(vectors_), &n,
                    COMPLEX(Bwork_), &n, REAL(values_), &work_query,
                    &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zhegv workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query.r);
    if (lwork < 1) {
      lwork = 1;
    }
    SEXP work_ = PROTECT(allocVector(CPLXSXP, lwork));
    F77_CALL(zhegv)(&itype, &jobz, &uplo, &n, COMPLEX(vectors_), &n,
                    COMPLEX(Bwork_), &n, REAL(values_), COMPLEX(work_),
                    &lwork, rwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK zhegv failed with info=%d", info);
    }
    UNPROTECT(1);
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
}

extern "C" SEXP eigencore_dense_complex_generalized_pencil_eigen(SEXP A_, SEXP B_) {
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
}

extern "C" SEXP eigencore_dense_generalized_schur(SEXP A_, SEXP B_,
                                                   SEXP vectors_, SEXP sort_code_) {
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
}

extern "C" SEXP eigencore_dense_complex_generalized_schur(SEXP A_, SEXP B_,
                                                           SEXP vectors_,
                                                           SEXP sort_code_) {
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
}

extern "C" SEXP eigencore_dense_is_symmetric(SEXP A_, SEXP tol_) {
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
  double scale = 1.0;
  for (int64_t i = 0; i < static_cast<int64_t>(n) * n; ++i) {
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
      if (fabs(a - b) > threshold) {
        return ScalarLogical(FALSE);
      }
    }
  }
  return ScalarLogical(TRUE);
}

extern "C" SEXP eigencore_dense_symmetric_eigen_selected(SEXP A_, SEXP k_, SEXP target_kind_) {
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
    int* isuppz = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(2 * (k > 1 ? k : 1)), sizeof(int))
    );
    std::memcpy(work_matrix, REAL(A_), sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
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
    int lwork = -1;
    int liwork = -1;
    double work_query = 0.0;
    int iwork_query = 0;
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, work_matrix, &n,
                     &vl, &vu, &il, &iu, &abstol,
                     &m_found, values_work, REAL(vectors_), &n,
                     isuppz, &work_query, &lwork,
                     &iwork_query, &liwork, &info FCONE FCONE FCONE);
    if (info != 0) {
      error("LAPACK dsyevr workspace query failed with info=%d", info);
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
    std::memcpy(work_matrix, REAL(A_), sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(n));
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n, work_matrix, &n,
                     &vl, &vu, &il, &iu, &abstol,
                     &m_found, values_work, REAL(vectors_), &n,
                     isuppz, work, &lwork,
                     iwork, &liwork, &info FCONE FCONE FCONE);
    if (info == 0 && m_found == k) {
      for (int col = 0; col < k; ++col) {
        REAL(values_)[col] = values_work[col];
      }
    }
    if (info != 0 || m_found != k) {
      error("LAPACK dsyevr failed with info=%d, found=%d", info, m_found);
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
}

extern "C" SEXP eigencore_dense_symmetric_eigen_dsyevx_selected(SEXP A_, SEXP k_, SEXP target_kind_) {
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
}

extern "C" SEXP eigencore_dense_generalized_spd_eigen(SEXP A_, SEXP B_) {
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
}

extern "C" SEXP eigencore_dense_svd(SEXP A_) {
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
  SEXP vt_ = PROTECT(allocMatrix(REALSXP, r, n));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, r));
  SEXP work_matrix_ = PROTECT(duplicate(A_));

  if (r > 0) {
    char jobu = 'S';
    char jobvt = 'S';
    int lda = m;
    int ldu = m;
    int ldvt = r;
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    F77_CALL(dgesvd)(&jobu, &jobvt, &m, &n, REAL(work_matrix_), &lda,
                     REAL(d_), REAL(u_), &ldu, REAL(vt_), &ldvt,
                     &work_query, &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgesvd workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    F77_CALL(dgesvd)(&jobu, &jobvt, &m, &n, REAL(work_matrix_), &lda,
                     REAL(d_), REAL(u_), &ldu, REAL(vt_), &ldvt,
                     REAL(work_), &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgesvd failed with info=%d", info);
    }
    UNPROTECT(1);

    for (int col = 0; col < r; ++col) {
      for (int row = 0; row < n; ++row) {
        REAL(v_)[row + col * n] = REAL(vt_)[col + row * r];
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
}

extern "C" SEXP eigencore_dense_complex_svd(SEXP A_) {
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
}

extern "C" SEXP eigencore_dense_generalized_svd(SEXP A_, SEXP B_) {
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
}

extern "C" SEXP eigencore_dense_complex_generalized_svd(SEXP A_, SEXP B_) {
  error("native complex GSVD requires a complex LAPACK GSVD driver, "
        "which this R LAPACK interface does not export");
}

extern "C" SEXP eigencore_tridiagonal_eigen(SEXP alpha_, SEXP beta_) {
  if (!isReal(alpha_) || !isReal(beta_)) {
    error("alpha and beta must be double vectors");
  }
  const int n = LENGTH(alpha_);
  if (LENGTH(beta_) < ((n > 0) ? n - 1 : 0)) {
    error("beta must have length at least length(alpha) - 1");
  }

  SEXP values_ = PROTECT(duplicate(alpha_));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, n));
  if (n > 0) {
    std::memset(REAL(vectors_), 0, sizeof(double) * static_cast<size_t>(n) * n);
  }

  if (n == 1) {
    REAL(vectors_)[0] = 1.0;
  } else if (n > 1) {
    SEXP offdiag_ = PROTECT(allocVector(REALSXP, n - 1));
    for (int i = 0; i < n - 1; ++i) {
      REAL(offdiag_)[i] = REAL(beta_)[i];
    }
    SEXP work_ = PROTECT(allocVector(REALSXP, 2 * n - 2));
    char jobz = 'V';
    int info = 0;
    F77_CALL(dstev)(&jobz, &n, REAL(values_), REAL(offdiag_),
                    REAL(vectors_), &n, REAL(work_), &info FCONE);
    if (info != 0) {
      error("LAPACK dstev failed with info=%d", info);
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
}

extern "C" SEXP eigencore_tridiagonal_eigen_selected(SEXP alpha_, SEXP beta_,
                                                     SEXP k_, SEXP target_kind_) {
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
}

extern "C" SEXP eigencore_bidiagonal_svd(SEXP alpha_, SEXP beta_) {
  if (!isReal(alpha_) || !isReal(beta_)) {
    error("alpha and beta must be double vectors");
  }
  const int n = LENGTH(alpha_);
  if (LENGTH(beta_) < ((n > 0) ? n - 1 : 0)) {
    error("beta must have length at least length(alpha) - 1");
  }

  SEXP d_ = PROTECT(allocVector(REALSXP, n));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP vt_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, n));
  SEXP B_ = PROTECT(allocMatrix(REALSXP, n, n));
  if (n > 0) {
    std::memset(REAL(B_), 0, sizeof(double) * static_cast<size_t>(n) * n);
    for (int i = 0; i < n; ++i) {
      REAL(B_)[i + i * n] = REAL(alpha_)[i];
    }
    for (int i = 0; i < n - 1; ++i) {
      REAL(B_)[i + (i + 1) * n] = REAL(beta_)[i];
    }

    char jobu = 'A';
    char jobvt = 'A';
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    F77_CALL(dgesvd)(&jobu, &jobvt, &n, &n, REAL(B_), &n,
                     REAL(d_), REAL(u_), &n, REAL(vt_), &n,
                     &work_query, &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgesvd workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
    F77_CALL(dgesvd)(&jobu, &jobvt, &n, &n, REAL(B_), &n,
                     REAL(d_), REAL(u_), &n, REAL(vt_), &n,
                     REAL(work_), &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dgesvd failed with info=%d", info);
    }
    UNPROTECT(1);

    for (int col = 0; col < n; ++col) {
      for (int row = 0; row < n; ++row) {
        REAL(v_)[row + col * n] = REAL(vt_)[col + row * n];
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
}

extern "C" SEXP eigencore_tridiagonal_solve(SEXP lower_, SEXP diag_,
                                            SEXP upper_, SEXP B_) {
  if (!isReal(lower_) || !isReal(diag_) || !isReal(upper_) || !isReal(B_)) {
    error("lower, diag, upper, and B must be double");
  }
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimB == R_NilValue) {
    error("B must be a matrix");
  }
  const int n = LENGTH(diag_);
  const int nrhs = INTEGER(dimB)[1];
  if (INTEGER(dimB)[0] != n || LENGTH(lower_) != n - 1 || LENGTH(upper_) != n - 1) {
    error("non-conformable tridiagonal solve inputs");
  }

  SEXP out_ = PROTECT(duplicate(B_));
  SEXP cprime_ = PROTECT(allocVector(REALSXP, n > 1 ? n - 1 : 1));
  SEXP denom_ = PROTECT(allocVector(REALSXP, n));
  double* out = REAL(out_);
  double* cprime = REAL(cprime_);
  double* denom = REAL(denom_);
  const double* lower = REAL(lower_);
  const double* diag = REAL(diag_);
  const double* upper = REAL(upper_);

  if (n == 0) {
    UNPROTECT(3);
    return out_;
  }
  if (fabs(diag[0]) <= DBL_EPSILON) {
    error("zero pivot in tridiagonal solve");
  }
  denom[0] = diag[0];
  if (n > 1) {
    cprime[0] = upper[0] / denom[0];
  }
  for (int i = 1; i < n; ++i) {
    denom[i] = diag[i] - lower[i - 1] * cprime[i - 1];
    if (fabs(denom[i]) <= DBL_EPSILON) {
      error("zero pivot in tridiagonal solve");
    }
    if (i < n - 1) {
      cprime[i] = upper[i] / denom[i];
    }
  }

  for (int rhs = 0; rhs < nrhs; ++rhs) {
    double* x = out + static_cast<int64_t>(rhs) * n;
    x[0] = x[0] / denom[0];
    for (int i = 1; i < n; ++i) {
      x[i] = (x[i] - lower[i - 1] * x[i - 1]) / denom[i];
    }
    for (int i = n - 2; i >= 0; --i) {
      x[i] -= cprime[i] * x[i + 1];
    }
  }

  UNPROTECT(3);
  return out_;
}
