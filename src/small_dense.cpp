#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

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
    double* offdiag = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n - 1), sizeof(double))
    );
    double* values_work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(n), sizeof(double))
    );
    double* work = reinterpret_cast<double*>(
      R_alloc(static_cast<size_t>(5 * n), sizeof(double))
    );
    int* iwork = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(5 * n), sizeof(int))
    );
    int* ifail = reinterpret_cast<int*>(
      R_alloc(static_cast<size_t>(n), sizeof(int))
    );
    for (int i = 0; i < n; ++i) {
      diag[i] = REAL(alpha_)[i];
    }
    for (int i = 0; i < n - 1; ++i) {
      offdiag[i] = REAL(beta_)[i];
    }

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
    F77_CALL(dstevx)(&jobz, &range, &n, diag, offdiag,
                     &vl, &vu, &il, &iu, &abstol, &m_found,
                     values_work, REAL(vectors_), &ldz,
                     work, iwork, ifail, &info FCONE FCONE);
    if (info != 0 || m_found != k) {
      error("LAPACK dstevx failed with info=%d, found=%d", info, m_found);
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
