#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <cstring>
#include <stdint.h>
#include "eigencore_operator.h"

struct DenseColumnMajorOperator {
  int64_t rows;
  int64_t cols;
  double* values;
};

struct CSCOperator {
  int64_t rows;
  int64_t cols;
  const int* row_idx;
  const int* col_ptr;
  const double* values;
};

struct DiagonalOperator {
  int64_t rows;
  const double* values;
  bool unit;
};

static void scale_or_zero_output(double* Y, int64_t rows, int64_t cols, double beta) {
  const int64_t len = rows * cols;
  if (beta == 0.0) {
    for (int64_t pos = 0; pos < len; ++pos) {
      Y[pos] = 0.0;
    }
  } else if (beta != 1.0) {
    for (int64_t pos = 0; pos < len; ++pos) {
      Y[pos] *= beta;
    }
  }
}

static SEXP workspace_counters(EigencoreWorkspace* workspace) {
  SEXP out = PROTECT(allocVector(INTSXP, 2));
  INTEGER(out)[0] = static_cast<int>(workspace->allocation_count);
  INTEGER(out)[1] = static_cast<int>(workspace->bytes_allocated);
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("allocation_count"));
  SET_STRING_ELT(names, 1, mkChar("bytes_allocated"));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

extern "C" int eigencore_dense_apply(void* impl,
                                      EigencoreTranspose op,
                                      int64_t block_cols,
                                      const double* X,
                                      int64_t ldx,
                                      double alpha,
                                      double beta,
                                      double* Y,
                                      int64_t ldy,
                                      EigencoreWorkspace* workspace) {
  (void) workspace;
  DenseColumnMajorOperator* dense = static_cast<DenseColumnMajorOperator*>(impl);
  const int64_t out_rows64 = (op == EIGENCORE_TRANSPOSE_ADJOINT) ? dense->cols : dense->rows;
  const int64_t inner64 = (op == EIGENCORE_TRANSPOSE_ADJOINT) ? dense->rows : dense->cols;
  if (ldx < inner64 || ldy < out_rows64) {
    return -1;
  }

  const char transa = (op == EIGENCORE_TRANSPOSE_ADJOINT) ? 'T' : 'N';
  const char transb = 'N';
  const int out_rows = static_cast<int>(out_rows64);
  const int block_cols_i = static_cast<int>(block_cols);
  const int inner = static_cast<int>(inner64);
  const int lda = static_cast<int>(dense->rows);
  const int ldb = static_cast<int>(ldx);
  const int ldc = static_cast<int>(ldy);
  double beta_blas = beta;

  F77_CALL(dgemm)(&transa, &transb, &out_rows, &block_cols_i, &inner,
                  &alpha, dense->values, &lda, const_cast<double*>(X), &ldb,
                  &beta_blas, Y, &ldc FCONE FCONE);
  return 0;
}

extern "C" int eigencore_csc_apply(void* impl,
                                    EigencoreTranspose op,
                                    int64_t block_cols,
                                    const double* X,
                                    int64_t ldx,
                                    double alpha,
                                    double beta,
                                    double* Y,
                                    int64_t ldy,
                                    EigencoreWorkspace* workspace) {
  (void) workspace;
  CSCOperator* csc = static_cast<CSCOperator*>(impl);
  const int64_t out_rows = (op == EIGENCORE_TRANSPOSE_ADJOINT) ? csc->cols : csc->rows;
  const int64_t inner = (op == EIGENCORE_TRANSPOSE_ADJOINT) ? csc->rows : csc->cols;
  if (ldx < inner || ldy < out_rows) {
    return -1;
  }

  scale_or_zero_output(Y, out_rows, block_cols, beta);

  if (op == EIGENCORE_TRANSPOSE_NONE) {
    for (int64_t col = 0; col < csc->cols; ++col) {
      for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
        const int row = csc->row_idx[pos];
        const double a = alpha * csc->values[pos];
        for (int64_t block = 0; block < block_cols; ++block) {
          Y[row + block * ldy] += a * X[col + block * ldx];
        }
      }
    }
  } else {
    for (int64_t col = 0; col < csc->cols; ++col) {
      for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
        const int row = csc->row_idx[pos];
        const double a = alpha * csc->values[pos];
        for (int64_t block = 0; block < block_cols; ++block) {
          Y[col + block * ldy] += a * X[row + block * ldx];
        }
      }
    }
  }
  return 0;
}

extern "C" int eigencore_diagonal_apply(void* impl,
                                         EigencoreTranspose op,
                                         int64_t block_cols,
                                         const double* X,
                                         int64_t ldx,
                                         double alpha,
                                         double beta,
                                         double* Y,
                                         int64_t ldy,
                                         EigencoreWorkspace* workspace) {
  (void) op;
  (void) workspace;
  DiagonalOperator* diag = static_cast<DiagonalOperator*>(impl);
  if (ldx < diag->rows || ldy < diag->rows) {
    return -1;
  }
  scale_or_zero_output(Y, diag->rows, block_cols, beta);
  for (int64_t block = 0; block < block_cols; ++block) {
    for (int64_t row = 0; row < diag->rows; ++row) {
      const double d = diag->unit ? 1.0 : diag->values[row];
      Y[row + block * ldy] += alpha * d * X[row + block * ldx];
    }
  }
  return 0;
}

extern "C" SEXP eigencore_dense_block_apply(SEXP A_, SEXP X_, SEXP alpha_,
                                            SEXP beta_, SEXP Y_,
                                            SEXP transpose_) {
  if (!isReal(A_) || !isReal(X_) || !isReal(Y_)) {
    error("A, X, and Y must be double matrices");
  }

  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimA == R_NilValue || dimX == R_NilValue || dimY == R_NilValue) {
    error("A, X, and Y must be matrices");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = INTEGER(dimY)[0];
  const int yc = INTEGER(dimY)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const double alpha = REAL(alpha_)[0];
  const double beta = REAL(beta_)[0];

  const int inner = transpose ? m : n;
  const int out_rows = transpose ? n : m;
  if (xr != inner) {
    error("non-conformable X for dense block apply");
  }
  if (yr != out_rows || yc != xc) {
    error("non-conformable Y for dense block apply");
  }

  SEXP out_ = PROTECT(duplicate(Y_));
  double* A = REAL(A_);
  double* X = REAL(X_);
  double* out = REAL(out_);

  DenseColumnMajorOperator impl = {m, n, A};
  const int status = eigencore_dense_apply(
    &impl,
    transpose ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    X,
    xr,
    alpha,
    beta,
    out,
    out_rows,
    nullptr
  );
  if (status != 0) {
    error("dense block apply failed");
  }

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_csc_block_apply(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                          SEXP X_, SEXP alpha_, SEXP beta_,
                                          SEXP Y_, SEXP transpose_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(X_) || !isReal(Y_)) {
    error("invalid CSC block apply inputs");
  }

  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || dimY == R_NilValue) {
    error("X and Y must be matrices");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = INTEGER(dimY)[0];
  const int yc = INTEGER(dimY)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const double alpha = REAL(alpha_)[0];
  const double beta = REAL(beta_)[0];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;

  if (xr != inner) {
    error("non-conformable X for CSC block apply");
  }
  if (yr != out_rows || yc != xc) {
    error("non-conformable Y for CSC block apply");
  }

  SEXP out_ = PROTECT(duplicate(Y_));
  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  const double* X = REAL(X_);
  double* out = REAL(out_);

  CSCOperator impl = {m, n, row_idx, col_ptr, values};
  const int status = eigencore_csc_apply(
    &impl,
    transpose ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    X,
    xr,
    alpha,
    beta,
    out,
    out_rows,
    nullptr
  );
  if (status != 0) {
    error("CSC block apply failed");
  }

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_diagonal_block_apply(SEXP x_, SEXP dim_, SEXP unit_,
                                               SEXP X_, SEXP alpha_, SEXP beta_,
                                               SEXP Y_) {
  if (!isReal(x_) || !isInteger(dim_) || !isLogical(unit_) ||
      !isReal(X_) || !isReal(Y_)) {
    error("invalid diagonal block apply inputs");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || dimY == R_NilValue) {
    error("X and Y must be matrices");
  }

  const int n = INTEGER(dim_)[0];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = INTEGER(dimY)[0];
  const int yc = INTEGER(dimY)[1];
  if (INTEGER(dim_)[1] != n || xr != n || yr != n || yc != xc) {
    error("non-conformable diagonal block apply inputs");
  }

  SEXP out_ = PROTECT(duplicate(Y_));
  DiagonalOperator impl = {n, REAL(x_), static_cast<bool>(LOGICAL(unit_)[0])};
  const int status = eigencore_diagonal_apply(
    &impl,
    EIGENCORE_TRANSPOSE_NONE,
    xc,
    REAL(X_),
    xr,
    REAL(alpha_)[0],
    REAL(beta_)[0],
    REAL(out_),
    n,
    nullptr
  );
  if (status != 0) {
    error("diagonal block apply failed");
  }
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_native_apply_noalloc_check(SEXP kind_, SEXP A_,
                                                     SEXP X_, SEXP Y_) {
  if (!isString(kind_) || LENGTH(kind_) != 1 || !isReal(X_) || !isReal(Y_)) {
    error("invalid native no-allocation check inputs");
  }
  const char* kind = CHAR(STRING_ELT(kind_, 0));
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || dimY == R_NilValue) {
    error("X and Y must be matrices");
  }
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int status = -1;

  if (strcmp(kind, "dense") == 0) {
    if (!isReal(A_)) {
      error("dense check requires a double matrix");
    }
    SEXP dimA = getAttrib(A_, R_DimSymbol);
    if (dimA == R_NilValue) {
      error("dense A must be a matrix");
    }
    DenseColumnMajorOperator impl = {
      INTEGER(dimA)[0],
      INTEGER(dimA)[1],
      REAL(A_)
    };
    status = eigencore_dense_apply(&impl, EIGENCORE_TRANSPOSE_NONE,
                                   INTEGER(dimX)[1], REAL(X_), INTEGER(dimX)[0],
                                   1.0, 0.0, REAL(Y_), INTEGER(dimY)[0],
                                   &workspace);
  } else if (strcmp(kind, "csc") == 0) {
    CSCOperator impl = {
      INTEGER(GET_SLOT(A_, install("Dim")))[0],
      INTEGER(GET_SLOT(A_, install("Dim")))[1],
      INTEGER(GET_SLOT(A_, install("i"))),
      INTEGER(GET_SLOT(A_, install("p"))),
      REAL(GET_SLOT(A_, install("x")))
    };
    status = eigencore_csc_apply(&impl, EIGENCORE_TRANSPOSE_NONE,
                                 INTEGER(dimX)[1], REAL(X_), INTEGER(dimX)[0],
                                 1.0, 0.0, REAL(Y_), INTEGER(dimY)[0],
                                 &workspace);
  } else if (strcmp(kind, "diagonal") == 0) {
    SEXP diag_slot = GET_SLOT(A_, install("diag"));
    const bool unit = strcmp(CHAR(STRING_ELT(diag_slot, 0)), "U") == 0;
    DiagonalOperator impl = {
      INTEGER(GET_SLOT(A_, install("Dim")))[0],
      REAL(GET_SLOT(A_, install("x"))),
      unit
    };
    status = eigencore_diagonal_apply(&impl, EIGENCORE_TRANSPOSE_NONE,
                                      INTEGER(dimX)[1], REAL(X_), INTEGER(dimX)[0],
                                      1.0, 0.0, REAL(Y_), INTEGER(dimY)[0],
                                      &workspace);
  } else {
    error("unknown native no-allocation check kind");
  }

  if (status != 0) {
    error("native no-allocation check apply failed");
  }
  return workspace_counters(&workspace);
}

extern "C" SEXP eigencore_col_norms(SEXP X_) {
  if (!isReal(X_)) {
    error("X must be a double matrix");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  if (dimX == R_NilValue) {
    error("X must be a matrix");
  }

  const int rows = INTEGER(dimX)[0];
  const int cols = INTEGER(dimX)[1];
  const double* X = REAL(X_);
  SEXP out_ = PROTECT(allocVector(REALSXP, cols));
  double* out = REAL(out_);
  for (int col = 0; col < cols; ++col) {
    long double sum = 0.0L;
    const int offset = col * rows;
    for (int row = 0; row < rows; ++row) {
      const long double value = X[offset + row];
      sum += value * value;
    }
    out[col] = sqrt(static_cast<double>(sum));
  }
  UNPROTECT(1);
  return out_;
}

static double max_orthogonality_loss(const double* gram, int k) {
  double loss = 0.0;
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row < k; ++row) {
      const double target = (row == col) ? 1.0 : 0.0;
      const double err = fabs(gram[row + col * k] - target);
      if (err > loss) {
        loss = err;
      }
    }
  }
  return loss;
}

static double column_norm(const double* X, int rows, int col) {
  long double sum = 0.0L;
  const int offset = col * rows;
  for (int row = 0; row < rows; ++row) {
    const long double value = X[offset + row];
    sum += value * value;
  }
  return sqrt(static_cast<double>(sum));
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

extern "C" SEXP eigencore_orthogonality_loss(SEXP Q_, SEXP B_) {
  if (!isReal(Q_)) {
    error("Q must be a double matrix");
  }
  SEXP dimQ = getAttrib(Q_, R_DimSymbol);
  if (dimQ == R_NilValue) {
    error("Q must be a matrix");
  }

  const int n = INTEGER(dimQ)[0];
  const int k = INTEGER(dimQ)[1];
  if (k == 0) {
    return ScalarReal(0.0);
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  SEXP gram_ = PROTECT(allocMatrix(REALSXP, k, k));

  if (B_ == R_NilValue) {
    F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                    &one, REAL(Q_), &n, REAL(Q_), &n,
                    &zero, REAL(gram_), &k FCONE FCONE);
  } else {
    if (!isReal(B_)) {
      error("B must be a double matrix");
    }
    SEXP dimB = getAttrib(B_, R_DimSymbol);
    if (dimB == R_NilValue ||
        INTEGER(dimB)[0] != n ||
        INTEGER(dimB)[1] != n) {
      error("B must be square with dimension matching nrow(Q)");
    }
    SEXP BQ_ = PROTECT(allocMatrix(REALSXP, n, k));
    F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                    &one, REAL(B_), &n, REAL(Q_), &n,
                    &zero, REAL(BQ_), &n FCONE FCONE);
    F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                    &one, REAL(Q_), &n, REAL(BQ_), &n,
                    &zero, REAL(gram_), &k FCONE FCONE);
    UNPROTECT(1);
  }

  const double loss = max_orthogonality_loss(REAL(gram_), k);
  UNPROTECT(1);
  return ScalarReal(loss);
}

extern "C" SEXP eigencore_dense_eigen_residuals(SEXP A_, SEXP values_,
                                                SEXP vectors_, SEXP B_) {
  if (!isReal(A_) || !isReal(values_) || !isReal(vectors_)) {
    error("A, values, and vectors must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimV = getAttrib(vectors_, R_DimSymbol);
  if (dimA == R_NilValue || dimV == R_NilValue) {
    error("A and vectors must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  const int rowsV = INTEGER(dimV)[0];
  const int k = INTEGER(dimV)[1];
  if (n != ncolA || rowsV != n || LENGTH(values_) != k) {
    error("non-conformable dense eigen residual inputs");
  }
  if (B_ != R_NilValue) {
    if (!isReal(B_)) {
      error("B must be a double matrix");
    }
    SEXP dimB = getAttrib(B_, R_DimSymbol);
    if (dimB == R_NilValue || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n) {
      error("B must be square with dimension matching A");
    }
  }

  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  SEXP residual_ = PROTECT(allocMatrix(REALSXP, n, k));
  F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                  &one, REAL(A_), &n, REAL(vectors_), &n,
                  &zero, REAL(residual_), &n FCONE FCONE);

  if (B_ == R_NilValue) {
    for (int col = 0; col < k; ++col) {
      const double lambda = REAL(values_)[col];
      const int offset = col * n;
      for (int row = 0; row < n; ++row) {
        REAL(residual_)[offset + row] -= lambda * REAL(vectors_)[offset + row];
      }
    }
  } else {
    SEXP Bv_ = PROTECT(allocMatrix(REALSXP, n, k));
    F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                    &one, REAL(B_), &n, REAL(vectors_), &n,
                    &zero, REAL(Bv_), &n FCONE FCONE);
    for (int col = 0; col < k; ++col) {
      const double lambda = REAL(values_)[col];
      const int offset = col * n;
      for (int row = 0; row < n; ++row) {
        REAL(residual_)[offset + row] -= lambda * REAL(Bv_)[offset + row];
      }
    }
    UNPROTECT(1);
  }

  SEXP out_ = PROTECT(allocVector(REALSXP, k));
  for (int col = 0; col < k; ++col) {
    REAL(out_)[col] = column_norm(REAL(residual_), n, col);
  }
  UNPROTECT(2);
  return out_;
}

extern "C" SEXP eigencore_dense_svd_residuals(SEXP A_, SEXP d_,
                                              SEXP u_, SEXP v_) {
  if (!isReal(A_) || !isReal(d_) || !isReal(u_) || !isReal(v_)) {
    error("A, d, u, and v must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimU = getAttrib(u_, R_DimSymbol);
  SEXP dimV = getAttrib(v_, R_DimSymbol);
  if (dimA == R_NilValue || dimU == R_NilValue || dimV == R_NilValue) {
    error("A, u, and v must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int rowsU = INTEGER(dimU)[0];
  const int k = INTEGER(dimU)[1];
  const int rowsV = INTEGER(dimV)[0];
  const int colsV = INTEGER(dimV)[1];
  if (rowsU != m || rowsV != n || colsV != k || LENGTH(d_) != k) {
    error("non-conformable dense SVD residual inputs");
  }

  const char notrans = 'N';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  SEXP left_matrix_ = PROTECT(allocMatrix(REALSXP, m, k));
  SEXP right_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));

  F77_CALL(dgemm)(&notrans, &notrans, &m, &k, &n,
                  &one, REAL(A_), &m, REAL(v_), &n,
                  &zero, REAL(left_matrix_), &m FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &n, &k, &m,
                  &one, REAL(A_), &m, REAL(u_), &m,
                  &zero, REAL(right_matrix_), &n FCONE FCONE);

  for (int col = 0; col < k; ++col) {
    const double sigma = REAL(d_)[col];
    const int left_offset = col * m;
    const int right_offset = col * n;
    for (int row = 0; row < m; ++row) {
      REAL(left_matrix_)[left_offset + row] -= sigma * REAL(u_)[left_offset + row];
    }
    for (int row = 0; row < n; ++row) {
      REAL(right_matrix_)[right_offset + row] -= sigma * REAL(v_)[right_offset + row];
    }
  }

  SEXP left_ = PROTECT(allocVector(REALSXP, k));
  SEXP right_ = PROTECT(allocVector(REALSXP, k));
  SEXP combined_ = PROTECT(allocVector(REALSXP, k));
  for (int col = 0; col < k; ++col) {
    const double left = column_norm(REAL(left_matrix_), m, col);
    const double right = column_norm(REAL(right_matrix_), n, col);
    REAL(left_)[col] = left;
    REAL(right_)[col] = right;
    REAL(combined_)[col] = sqrt(left * left + right * right);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, left_);
  SET_VECTOR_ELT(out_, 1, right_);
  SET_VECTOR_ELT(out_, 2, combined_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("left"));
  SET_STRING_ELT(names_, 1, mkChar("right"));
  SET_STRING_ELT(names_, 2, mkChar("combined"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_rayleigh_ritz_symmetric(SEXP A_, SEXP Q_) {
  if (!isReal(A_) || !isReal(Q_)) {
    error("A and Q must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimQ = getAttrib(Q_, R_DimSymbol);
  if (dimA == R_NilValue || dimQ == R_NilValue) {
    error("A and Q must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int ncolA = INTEGER(dimA)[1];
  const int rowsQ = INTEGER(dimQ)[0];
  const int k = INTEGER(dimQ)[1];
  if (n != ncolA || rowsQ != n) {
    error("non-conformable Rayleigh-Ritz inputs");
  }

  SEXP values_ = PROTECT(allocVector(REALSXP, k));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP projected_ = PROTECT(allocMatrix(REALSXP, k, k));
  if (k == 0) {
    SEXP out_ = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out_, 0, values_);
    SET_VECTOR_ELT(out_, 1, vectors_);
    SET_VECTOR_ELT(out_, 2, projected_);
    SEXP names_ = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(names_, 0, mkChar("values"));
    SET_STRING_ELT(names_, 1, mkChar("vectors"));
    SET_STRING_ELT(names_, 2, mkChar("projected"));
    setAttrib(out_, R_NamesSymbol, names_);
    UNPROTECT(5);
    return out_;
  }

  const char notrans = 'N';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  SEXP AQ_ = PROTECT(allocMatrix(REALSXP, n, k));
  F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                  &one, REAL(A_), &n, REAL(Q_), &n,
                  &zero, REAL(AQ_), &n FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, REAL(Q_), &n, REAL(AQ_), &n,
                  &zero, REAL(projected_), &k FCONE FCONE);

  for (int col = 0; col < k; ++col) {
    for (int row = col + 1; row < k; ++row) {
      const double avg = 0.5 * (REAL(projected_)[row + col * k] +
                                REAL(projected_)[col + row * k]);
      REAL(projected_)[row + col * k] = avg;
      REAL(projected_)[col + row * k] = avg;
    }
  }

  SEXP eigvecs_ = PROTECT(duplicate(projected_));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dsyev)(&jobz, &uplo, &k, REAL(eigvecs_), &k, REAL(values_),
                  &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dsyev workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  SEXP work_ = PROTECT(allocVector(REALSXP, lwork));
  F77_CALL(dsyev)(&jobz, &uplo, &k, REAL(eigvecs_), &k, REAL(values_),
                  REAL(work_), &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dsyev failed with info=%d", info);
  }

  F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &k,
                  &one, REAL(Q_), &n, REAL(eigvecs_), &k,
                  &zero, REAL(vectors_), &n FCONE FCONE);

  SEXP out_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, projected_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("projected"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(8);
  return out_;
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

static const R_CallMethodDef CallEntries[] = {
  {"eigencore_dense_block_apply", (DL_FUNC) &eigencore_dense_block_apply, 6},
  {"eigencore_csc_block_apply", (DL_FUNC) &eigencore_csc_block_apply, 9},
  {"eigencore_diagonal_block_apply", (DL_FUNC) &eigencore_diagonal_block_apply, 7},
  {"eigencore_native_apply_noalloc_check", (DL_FUNC) &eigencore_native_apply_noalloc_check, 4},
  {"eigencore_col_norms", (DL_FUNC) &eigencore_col_norms, 1},
  {"eigencore_mgs2", (DL_FUNC) &eigencore_mgs2, 2},
  {"eigencore_orthogonality_loss", (DL_FUNC) &eigencore_orthogonality_loss, 2},
  {"eigencore_dense_eigen_residuals", (DL_FUNC) &eigencore_dense_eigen_residuals, 4},
  {"eigencore_dense_svd_residuals", (DL_FUNC) &eigencore_dense_svd_residuals, 4},
  {"eigencore_rayleigh_ritz_symmetric", (DL_FUNC) &eigencore_rayleigh_ritz_symmetric, 2},
  {"eigencore_dense_symmetric_eigen", (DL_FUNC) &eigencore_dense_symmetric_eigen, 1},
  {"eigencore_dense_generalized_spd_eigen", (DL_FUNC) &eigencore_dense_generalized_spd_eigen, 2},
  {"eigencore_dense_svd", (DL_FUNC) &eigencore_dense_svd, 1},
  {NULL, NULL, 0}
};

extern "C" void R_init_eigencore(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
