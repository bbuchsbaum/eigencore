#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdint.h>

struct BasisWorkspace {
  int64_t rows;
  int64_t basis_cols;
  int64_t block_cols;
  double* coeff;
  int64_t allocation_count;
  int64_t bytes_allocated;
};

static BasisWorkspace* basis_workspace_from_xptr(SEXP workspace_) {
  if (TYPEOF(workspace_) != EXTPTRSXP) {
    error("workspace must be an external pointer");
  }
  BasisWorkspace* workspace = static_cast<BasisWorkspace*>(R_ExternalPtrAddr(workspace_));
  if (workspace == nullptr) {
    error("basis workspace pointer is null");
  }
  return workspace;
}

static void basis_workspace_finalizer(SEXP workspace_) {
  BasisWorkspace* workspace = static_cast<BasisWorkspace*>(R_ExternalPtrAddr(workspace_));
  if (workspace != nullptr) {
    std::free(workspace->coeff);
    delete workspace;
    R_ClearExternalPtr(workspace_);
  }
}

static int mgs_once_dense(const double* X,
                          int n,
                          int p,
                          double tol,
                          double* Q,
                          int ldq,
                          double* R,
                          int ldr) {
  int rank = 0;
  for (int col = 0; col < p; ++col) {
    double* v = Q + rank * ldq;
    const double* x_col = X + col * n;
    for (int row = 0; row < n; ++row) {
      v[row] = x_col[row];
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int prev = 0; prev < rank; ++prev) {
        const double* q_prev = Q + prev * ldq;
        long double dot = 0.0L;
        for (int row = 0; row < n; ++row) {
          dot += static_cast<long double>(q_prev[row]) * v[row];
        }
        const double r = static_cast<double>(dot);
        R[prev + col * ldr] += r;
        for (int row = 0; row < n; ++row) {
          v[row] -= r * q_prev[row];
        }
      }
    }

    long double norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      norm2 += static_cast<long double>(v[row]) * v[row];
    }
    const double rjj = sqrt(static_cast<double>(norm2));
    if (rjj > tol) {
      for (int row = 0; row < n; ++row) {
        v[row] /= rjj;
      }
      R[rank + col * ldr] = rjj;
      ++rank;
    }
  }
  return rank;
}

static void zero_strict_lower(double* A, int n) {
  for (int col = 0; col < n; ++col) {
    for (int row = col + 1; row < n; ++row) {
      A[row + col * n] = 0.0;
    }
  }
}

static void symmetrize_full(double* A, int n) {
  for (int col = 0; col < n; ++col) {
    for (int row = col + 1; row < n; ++row) {
      const double value = 0.5 * (A[row + col * n] + A[col + row * n]);
      A[row + col * n] = value;
      A[col + row * n] = value;
    }
  }
}

extern "C" SEXP eigencore_mgs2(SEXP X_, SEXP tol_) {
  if (!isReal(X_)) {
    error("X must be a double matrix");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  if (dimX == R_NilValue) {
    error("X must be a matrix");
  }

  const int n = INTEGER(dimX)[0];
  const int p = INTEGER(dimX)[1];
  const double tol = REAL(tol_)[0];

  SEXP Q1_ = PROTECT(allocMatrix(REALSXP, n, p));
  SEXP R1_ = PROTECT(allocMatrix(REALSXP, p, p));
  std::memset(REAL(Q1_), 0, sizeof(double) * static_cast<size_t>(n) * p);
  std::memset(REAL(R1_), 0, sizeof(double) * static_cast<size_t>(p) * p);
  const int rank1 = mgs_once_dense(REAL(X_), n, p, tol, REAL(Q1_), n, REAL(R1_), p);

  if (rank1 == 0) {
    SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, 0));
    SEXP R_ = PROTECT(allocMatrix(REALSXP, 0, p));
    SEXP out_ = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out_, 0, Q_);
    SET_VECTOR_ELT(out_, 1, R_);
    SET_VECTOR_ELT(out_, 2, ScalarInteger(0));
    SEXP names_ = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names_, 0, mkChar("Q"));
    SET_STRING_ELT(names_, 1, mkChar("R"));
    SET_STRING_ELT(names_, 2, mkChar("rank"));
    setAttrib(out_, R_NamesSymbol, names_);
    UNPROTECT(6);
    return out_;
  }

  SEXP Q1_compact_ = PROTECT(allocMatrix(REALSXP, n, rank1));
  for (int col = 0; col < rank1; ++col) {
    std::memcpy(REAL(Q1_compact_) + col * n,
                REAL(Q1_) + col * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  SEXP Q2_ = PROTECT(allocMatrix(REALSXP, n, rank1));
  SEXP R2_ = PROTECT(allocMatrix(REALSXP, rank1, rank1));
  std::memset(REAL(Q2_), 0, sizeof(double) * static_cast<size_t>(n) * rank1);
  std::memset(REAL(R2_), 0, sizeof(double) * static_cast<size_t>(rank1) * rank1);
  const int rank2 = mgs_once_dense(REAL(Q1_compact_), n, rank1, tol,
                                   REAL(Q2_), n, REAL(R2_), rank1);

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, rank2));
  for (int col = 0; col < rank2; ++col) {
    std::memcpy(REAL(Q_) + col * n,
                REAL(Q2_) + col * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  SEXP R_ = PROTECT(allocMatrix(REALSXP, rank2, p));
  std::memset(REAL(R_), 0, sizeof(double) * static_cast<size_t>(rank2) * p);
  for (int col = 0; col < p; ++col) {
    for (int row = 0; row < rank2; ++row) {
      long double value = 0.0L;
      for (int inner = 0; inner < rank1; ++inner) {
        value += static_cast<long double>(REAL(R2_)[row + inner * rank1]) *
                 REAL(R1_)[inner + col * p];
      }
      REAL(R_)[row + col * rank2] = static_cast<double>(value);
    }
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, R_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(rank2));
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("R"));
  SET_STRING_ELT(names_, 2, mkChar("rank"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(9);
  return out_;
}

extern "C" SEXP eigencore_cholqr2(SEXP X_) {
  if (!isReal(X_)) {
    error("X must be a double matrix");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  if (dimX == R_NilValue) {
    error("X must be a matrix");
  }

  const int n = INTEGER(dimX)[0];
  const int p = INTEGER(dimX)[1];
  SEXP Q_ = PROTECT(duplicate(X_));
  SEXP R1_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP R2_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP R_ = PROTECT(allocMatrix(REALSXP, p, p));
  if (p == 0) {
    SEXP out_ = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out_, 0, Q_);
    SET_VECTOR_ELT(out_, 1, R_);
    SET_VECTOR_ELT(out_, 2, ScalarInteger(0));
    SEXP names_ = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names_, 0, mkChar("Q"));
    SET_STRING_ELT(names_, 1, mkChar("R"));
    SET_STRING_ELT(names_, 2, mkChar("rank"));
    setAttrib(out_, R_NamesSymbol, names_);
    UNPROTECT(6);
    return out_;
  }

  const char trans = 'T';
  const char notrans = 'N';
  const char right = 'R';
  const char uplo = 'U';
  const char diag = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  F77_CALL(dgemm)(&trans, &notrans, &p, &p, &n,
                  &one, REAL(X_), &n, REAL(X_), &n,
                  &zero, REAL(R1_), &p FCONE FCONE);
  symmetrize_full(REAL(R1_), p);
  F77_CALL(dpotrf)(&uplo, &p, REAL(R1_), &p, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed in CholQR2 first pass with info=%d", info);
  }
  zero_strict_lower(REAL(R1_), p);

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &n, &p, &one,
                  REAL(R1_), &p, REAL(Q_), &n FCONE FCONE FCONE FCONE);

  F77_CALL(dgemm)(&trans, &notrans, &p, &p, &n,
                  &one, REAL(Q_), &n, REAL(Q_), &n,
                  &zero, REAL(R2_), &p FCONE FCONE);
  symmetrize_full(REAL(R2_), p);
  F77_CALL(dpotrf)(&uplo, &p, REAL(R2_), &p, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed in CholQR2 second pass with info=%d", info);
  }
  zero_strict_lower(REAL(R2_), p);

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &n, &p, &one,
                  REAL(R2_), &p, REAL(Q_), &n FCONE FCONE FCONE FCONE);

  F77_CALL(dgemm)(&notrans, &notrans, &p, &p, &p,
                  &one, REAL(R2_), &p, REAL(R1_), &p,
                  &zero, REAL(R_), &p FCONE FCONE);
  zero_strict_lower(REAL(R_), p);

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, R_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(p));
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("R"));
  SET_STRING_ELT(names_, 2, mkChar("rank"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}

extern "C" SEXP eigencore_b_cholqr2(SEXP X_, SEXP B_) {
  if (!isReal(X_) || !isReal(B_)) {
    error("X and B must be double matrices");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  if (dimX == R_NilValue || dimB == R_NilValue) {
    error("X and B must be matrices");
  }

  const int n = INTEGER(dimX)[0];
  const int p = INTEGER(dimX)[1];
  if (INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
    error("B must be square with dimension matching nrow(X)");
  }

  SEXP Q_ = PROTECT(duplicate(X_));
  SEXP BX_ = PROTECT(allocMatrix(REALSXP, n, p));
  SEXP R1_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP R2_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP R_ = PROTECT(allocMatrix(REALSXP, p, p));
  if (p == 0) {
    SEXP out_ = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out_, 0, Q_);
    SET_VECTOR_ELT(out_, 1, R_);
    SET_VECTOR_ELT(out_, 2, ScalarInteger(0));
    SEXP names_ = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names_, 0, mkChar("Q"));
    SET_STRING_ELT(names_, 1, mkChar("R"));
    SET_STRING_ELT(names_, 2, mkChar("rank"));
    setAttrib(out_, R_NamesSymbol, names_);
    UNPROTECT(7);
    return out_;
  }

  const char trans = 'T';
  const char notrans = 'N';
  const char right = 'R';
  const char uplo = 'U';
  const char diag = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  F77_CALL(dgemm)(&notrans, &notrans, &n, &p, &n,
                  &one, REAL(B_), &n, REAL(X_), &n,
                  &zero, REAL(BX_), &n FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &p, &p, &n,
                  &one, REAL(X_), &n, REAL(BX_), &n,
                  &zero, REAL(R1_), &p FCONE FCONE);
  symmetrize_full(REAL(R1_), p);
  F77_CALL(dpotrf)(&uplo, &p, REAL(R1_), &p, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed in B-CholQR2 first pass with info=%d", info);
  }
  zero_strict_lower(REAL(R1_), p);

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &n, &p, &one,
                  REAL(R1_), &p, REAL(Q_), &n FCONE FCONE FCONE FCONE);

  F77_CALL(dgemm)(&notrans, &notrans, &n, &p, &n,
                  &one, REAL(B_), &n, REAL(Q_), &n,
                  &zero, REAL(BX_), &n FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &p, &p, &n,
                  &one, REAL(Q_), &n, REAL(BX_), &n,
                  &zero, REAL(R2_), &p FCONE FCONE);
  symmetrize_full(REAL(R2_), p);
  F77_CALL(dpotrf)(&uplo, &p, REAL(R2_), &p, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed in B-CholQR2 second pass with info=%d", info);
  }
  zero_strict_lower(REAL(R2_), p);

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &n, &p, &one,
                  REAL(R2_), &p, REAL(Q_), &n FCONE FCONE FCONE FCONE);

  F77_CALL(dgemm)(&notrans, &notrans, &p, &p, &p,
                  &one, REAL(R2_), &p, REAL(R1_), &p,
                  &zero, REAL(R_), &p FCONE FCONE);
  zero_strict_lower(REAL(R_), p);

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, R_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(p));
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("R"));
  SET_STRING_ELT(names_, 2, mkChar("rank"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_diagonal_b_cholqr2(SEXP X_, SEXP diag_, SEXP unit_) {
  if (!isReal(X_) || !isReal(diag_)) {
    error("X and diagonal values must be double");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  if (dimX == R_NilValue) {
    error("X must be a matrix");
  }

  const int n = INTEGER(dimX)[0];
  const int p = INTEGER(dimX)[1];
  const int unit = asLogical(unit_) == TRUE;
  if (!unit && LENGTH(diag_) != n) {
    error("diagonal values must have length nrow(X)");
  }
  if (!unit) {
    for (int row = 0; row < n; ++row) {
      if (REAL(diag_)[row] <= 0.0) {
        error("B diagonal must be positive");
      }
    }
  }

  SEXP Q_ = PROTECT(duplicate(X_));
  SEXP BQ_ = PROTECT(allocMatrix(REALSXP, n, p));
  SEXP R1_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP R2_ = PROTECT(allocMatrix(REALSXP, p, p));
  SEXP R_ = PROTECT(allocMatrix(REALSXP, p, p));
  if (p == 0) {
    SEXP out_ = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(out_, 0, Q_);
    SET_VECTOR_ELT(out_, 1, R_);
    SET_VECTOR_ELT(out_, 2, ScalarInteger(0));
    SET_VECTOR_ELT(out_, 3, BQ_);
    SEXP names_ = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names_, 0, mkChar("Q"));
    SET_STRING_ELT(names_, 1, mkChar("R"));
    SET_STRING_ELT(names_, 2, mkChar("rank"));
    SET_STRING_ELT(names_, 3, mkChar("BQ"));
    setAttrib(out_, R_NamesSymbol, names_);
    UNPROTECT(7);
    return out_;
  }

  const char trans = 'T';
  const char notrans = 'N';
  const char right = 'R';
  const char uplo = 'U';
  const char tri_diag = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  for (int col = 0; col < p; ++col) {
    for (int row = 0; row < n; ++row) {
      const double weight = unit ? 1.0 : REAL(diag_)[row];
      REAL(BQ_)[row + static_cast<int64_t>(col) * n] =
        weight * REAL(X_)[row + static_cast<int64_t>(col) * n];
    }
  }
  F77_CALL(dgemm)(&trans, &notrans, &p, &p, &n,
                  &one, REAL(X_), &n, REAL(BQ_), &n,
                  &zero, REAL(R1_), &p FCONE FCONE);
  symmetrize_full(REAL(R1_), p);
  F77_CALL(dpotrf)(&uplo, &p, REAL(R1_), &p, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed in diagonal B-CholQR2 first pass with info=%d", info);
  }
  zero_strict_lower(REAL(R1_), p);

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &tri_diag, &n, &p, &one,
                  REAL(R1_), &p, REAL(Q_), &n FCONE FCONE FCONE FCONE);

  for (int col = 0; col < p; ++col) {
    for (int row = 0; row < n; ++row) {
      const double weight = unit ? 1.0 : REAL(diag_)[row];
      REAL(BQ_)[row + static_cast<int64_t>(col) * n] =
        weight * REAL(Q_)[row + static_cast<int64_t>(col) * n];
    }
  }
  F77_CALL(dgemm)(&trans, &notrans, &p, &p, &n,
                  &one, REAL(Q_), &n, REAL(BQ_), &n,
                  &zero, REAL(R2_), &p FCONE FCONE);
  symmetrize_full(REAL(R2_), p);
  F77_CALL(dpotrf)(&uplo, &p, REAL(R2_), &p, &info FCONE);
  if (info != 0) {
    error("LAPACK dpotrf failed in diagonal B-CholQR2 second pass with info=%d", info);
  }
  zero_strict_lower(REAL(R2_), p);

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &tri_diag, &n, &p, &one,
                  REAL(R2_), &p, REAL(Q_), &n FCONE FCONE FCONE FCONE);

  for (int col = 0; col < p; ++col) {
    for (int row = 0; row < n; ++row) {
      const double weight = unit ? 1.0 : REAL(diag_)[row];
      REAL(BQ_)[row + static_cast<int64_t>(col) * n] =
        weight * REAL(Q_)[row + static_cast<int64_t>(col) * n];
    }
  }

  F77_CALL(dgemm)(&notrans, &notrans, &p, &p, &p,
                  &one, REAL(R2_), &p, REAL(R1_), &p,
                  &zero, REAL(R_), &p FCONE FCONE);
  zero_strict_lower(REAL(R_), p);

  SEXP out_ = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, R_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(p));
  SET_VECTOR_ELT(out_, 3, BQ_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("R"));
  SET_STRING_ELT(names_, 2, mkChar("rank"));
  SET_STRING_ELT(names_, 3, mkChar("BQ"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_reorthogonalize_against(SEXP X_, SEXP Q_, SEXP passes_) {
  if (!isReal(X_) || !isReal(Q_) || !isInteger(passes_)) {
    error("X, Q, and passes must be double, double, and integer");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimQ = getAttrib(Q_, R_DimSymbol);
  if (dimX == R_NilValue || dimQ == R_NilValue) {
    error("X and Q must be matrices");
  }

  const int n = INTEGER(dimX)[0];
  const int p = INTEGER(dimX)[1];
  const int q_rows = INTEGER(dimQ)[0];
  const int q_cols = INTEGER(dimQ)[1];
  const int passes = INTEGER(passes_)[0];
  if (q_rows != n) {
    error("Q must have nrow(Q) == nrow(X)");
  }
  if (passes < 0) {
    error("passes must be non-negative");
  }

  SEXP out_ = PROTECT(duplicate(X_));
  if (n == 0 || p == 0 || q_cols == 0 || passes == 0) {
    UNPROTECT(1);
    return out_;
  }

  SEXP coeff_ = PROTECT(allocMatrix(REALSXP, q_cols, p));
  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;

  for (int pass = 0; pass < passes; ++pass) {
    F77_CALL(dgemm)(&trans, &notrans, &q_cols, &p, &n,
                    &one, REAL(Q_), &n, REAL(out_), &n,
                    &zero, REAL(coeff_), &q_cols FCONE FCONE);
    F77_CALL(dgemm)(&notrans, &notrans, &n, &p, &q_cols,
                    &minus_one, REAL(Q_), &n, REAL(coeff_), &q_cols,
                    &one, REAL(out_), &n FCONE FCONE);
  }

  UNPROTECT(2);
  return out_;
}

extern "C" SEXP eigencore_basis_workspace_create(SEXP rows_, SEXP basis_cols_, SEXP block_cols_) {
  const int64_t rows = static_cast<int64_t>(asReal(rows_));
  const int64_t basis_cols = static_cast<int64_t>(asReal(basis_cols_));
  const int64_t block_cols = static_cast<int64_t>(asReal(block_cols_));
  if (rows < 0 || basis_cols < 0 || block_cols < 0) {
    error("workspace dimensions must be non-negative");
  }

  BasisWorkspace* workspace = new BasisWorkspace;
  workspace->rows = rows;
  workspace->basis_cols = basis_cols;
  workspace->block_cols = block_cols;
  workspace->allocation_count = 0;
  workspace->bytes_allocated = 0;
  workspace->coeff = nullptr;

  const int64_t coeff_len = basis_cols * block_cols;
  if (coeff_len > 0) {
    workspace->coeff = static_cast<double*>(std::calloc(static_cast<size_t>(coeff_len), sizeof(double)));
    if (workspace->coeff == nullptr) {
      delete workspace;
      error("failed to allocate basis workspace coefficient buffer");
    }
    workspace->allocation_count = 1;
    workspace->bytes_allocated = coeff_len * static_cast<int64_t>(sizeof(double));
  }

  SEXP ptr_ = PROTECT(R_MakeExternalPtr(workspace, install("eigencore_basis_workspace"), R_NilValue));
  R_RegisterCFinalizerEx(ptr_, basis_workspace_finalizer, TRUE);
  UNPROTECT(1);
  return ptr_;
}

extern "C" SEXP eigencore_basis_workspace_info(SEXP workspace_) {
  BasisWorkspace* workspace = basis_workspace_from_xptr(workspace_);
  SEXP out_ = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out_, 0, ScalarReal(static_cast<double>(workspace->rows)));
  SET_VECTOR_ELT(out_, 1, ScalarReal(static_cast<double>(workspace->basis_cols)));
  SET_VECTOR_ELT(out_, 2, ScalarReal(static_cast<double>(workspace->block_cols)));
  SET_VECTOR_ELT(out_, 3, ScalarReal(static_cast<double>(workspace->allocation_count)));
  SET_VECTOR_ELT(out_, 4, ScalarReal(static_cast<double>(workspace->bytes_allocated)));
  SEXP names_ = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(names_, 0, mkChar("rows"));
  SET_STRING_ELT(names_, 1, mkChar("basis_cols"));
  SET_STRING_ELT(names_, 2, mkChar("block_cols"));
  SET_STRING_ELT(names_, 3, mkChar("allocation_count"));
  SET_STRING_ELT(names_, 4, mkChar("bytes_allocated"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(2);
  return out_;
}

extern "C" SEXP eigencore_reorthogonalize_against_workspace(SEXP X_, SEXP Q_,
                                                            SEXP passes_,
                                                            SEXP workspace_) {
  if (!isReal(X_) || !isReal(Q_) || !isInteger(passes_)) {
    error("X, Q, and passes must be double, double, and integer");
  }
  BasisWorkspace* workspace = basis_workspace_from_xptr(workspace_);
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimQ = getAttrib(Q_, R_DimSymbol);
  if (dimX == R_NilValue || dimQ == R_NilValue) {
    error("X and Q must be matrices");
  }

  const int n = INTEGER(dimX)[0];
  const int p = INTEGER(dimX)[1];
  const int q_rows = INTEGER(dimQ)[0];
  const int q_cols = INTEGER(dimQ)[1];
  const int passes = INTEGER(passes_)[0];
  if (q_rows != n) {
    error("Q must have nrow(Q) == nrow(X)");
  }
  if (passes < 0) {
    error("passes must be non-negative");
  }
  if (workspace->rows < n || workspace->basis_cols < q_cols || workspace->block_cols < p) {
    error("basis workspace is too small for requested reorthogonalization");
  }

  SEXP out_ = PROTECT(duplicate(X_));
  if (n == 0 || p == 0 || q_cols == 0 || passes == 0) {
    UNPROTECT(1);
    return out_;
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  for (int pass = 0; pass < passes; ++pass) {
    F77_CALL(dgemm)(&trans, &notrans, &q_cols, &p, &n,
                    &one, REAL(Q_), &n, REAL(out_), &n,
                    &zero, workspace->coeff, &q_cols FCONE FCONE);
    F77_CALL(dgemm)(&notrans, &notrans, &n, &p, &q_cols,
                    &minus_one, REAL(Q_), &n, workspace->coeff, &q_cols,
                    &one, REAL(out_), &n FCONE FCONE);
  }

  UNPROTECT(1);
  return out_;
}
