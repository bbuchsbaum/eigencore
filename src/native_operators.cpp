#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cstring>
#include <cmath>
#include <climits>
#include "eigencore_common.h"
#include "native_operators.h"


static SEXP native_operator_workspace_counters(EigencoreWorkspace* workspace) {
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

static void scale_or_zero_output(double* Y, int64_t rows, int64_t cols, double beta) {
  const int64_t len = rows * cols;
  if (beta == 0.0) {
    std::memset(Y, 0, sizeof(double) * static_cast<size_t>(len));
  } else if (beta != 1.0) {
    for (int64_t pos = 0; pos < len; ++pos) {
      Y[pos] *= beta;
    }
  }
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
  if (!eigencore_int_indexable(out_rows64) ||
      !eigencore_int_indexable(inner64) ||
      !eigencore_int_indexable(block_cols) ||
      !eigencore_int_indexable(dense->rows) ||
      !eigencore_int_indexable(ldx) ||
      !eigencore_int_indexable(ldy)) {
    return -2;
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

extern "C" int eigencore_dense_shift_invert_apply(void* impl,
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
  DenseShiftInvertOperator* si = static_cast<DenseShiftInvertOperator*>(impl);
  if (op != EIGENCORE_TRANSPOSE_NONE) {
    return -1;
  }
  if (block_cols != 1) {
    return -1;
  }
  if (ldx < si->n || ldy < si->n) {
    return -1;
  }

  scale_or_zero_output(Y, si->n, block_cols, beta);
  std::memcpy(si->work, X, sizeof(double) * static_cast<size_t>(si->n));
  char trans = 'N';
  const int nrhs = 1;
  int info = 0;
  F77_CALL(dgetrs)(&trans, &si->n, &nrhs, si->lu, &si->n, si->pivots,
                   si->work, &si->n, &info FCONE);
  if (info != 0) {
    return info;
  }
  for (int row = 0; row < si->n; ++row) {
    Y[row] += alpha * si->work[row];
  }
  return 0;
}

extern "C" int eigencore_dense_generalized_shift_invert_apply(void* impl,
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
  DenseGeneralizedShiftInvertOperator* si =
    static_cast<DenseGeneralizedShiftInvertOperator*>(impl);
  if (op != EIGENCORE_TRANSPOSE_NONE || block_cols != 1) {
    return -1;
  }
  if (ldx < si->n || ldy < si->n) {
    return -1;
  }

  scale_or_zero_output(Y, si->n, block_cols, beta);
  char uplo = 'U';
  char trans_T = 'T';
  char trans_N = 'N';
  char diag = 'N';
  int inc = 1;
  int nrhs = 1;
  std::memcpy(si->rhs, X, sizeof(double) * static_cast<size_t>(si->n));
  F77_CALL(dtrmv)(&uplo, &trans_T, &diag, &si->n, si->chol, &si->n,
                  si->rhs, &inc FCONE FCONE FCONE);
  std::memcpy(si->sol, si->rhs, sizeof(double) * static_cast<size_t>(si->n));
  int info = 0;
  F77_CALL(dgetrs)(&trans_N, &si->n, &nrhs, si->lu, &si->n, si->pivots,
                   si->sol, &si->n, &info FCONE);
  if (info != 0) {
    return info;
  }
  F77_CALL(dtrmv)(&uplo, &trans_N, &diag, &si->n, si->chol, &si->n,
                  si->sol, &inc FCONE FCONE FCONE);
  for (int row = 0; row < si->n; ++row) {
    Y[row] += alpha * si->sol[row];
  }
  return 0;
}

extern "C" int eigencore_tridiagonal_shift_invert_apply(void* impl,
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
  TridiagonalShiftInvertOperator* si =
    static_cast<TridiagonalShiftInvertOperator*>(impl);
  if (op != EIGENCORE_TRANSPOSE_NONE || block_cols != 1) {
    return -1;
  }
  if (ldx < si->n || ldy < si->n) {
    return -1;
  }

  const int n = si->n;
  scale_or_zero_output(Y, n, block_cols, beta);
  std::memcpy(si->work, X, sizeof(double) * static_cast<size_t>(n));
  si->work[0] /= si->denom[0];
  for (int i = 1; i < n; ++i) {
    si->work[i] = (si->work[i] - si->lower[i - 1] * si->work[i - 1]) /
      si->denom[i];
  }
  for (int i = n - 2; i >= 0; --i) {
    si->work[i] -= si->cprime[i] * si->work[i + 1];
  }
  for (int row = 0; row < n; ++row) {
    Y[row] += alpha * si->work[row];
  }
  return 0;
}

extern "C" int eigencore_tridiagonal_generalized_shift_invert_apply(
    void* impl,
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
  TridiagonalGeneralizedShiftInvertOperator* si =
    static_cast<TridiagonalGeneralizedShiftInvertOperator*>(impl);
  if (op != EIGENCORE_TRANSPOSE_NONE || block_cols != 1) {
    return -1;
  }
  if (ldx < si->n || ldy < si->n) {
    return -1;
  }

  const int n = si->n;
  scale_or_zero_output(Y, n, block_cols, beta);
  for (int row = 0; row < n; ++row) {
    si->work[row] = si->sqrt_metric[row] * X[row];
  }
  si->work[0] /= si->denom[0];
  for (int i = 1; i < n; ++i) {
    si->work[i] = (si->work[i] - si->lower[i - 1] * si->work[i - 1]) /
      si->denom[i];
  }
  for (int i = n - 2; i >= 0; --i) {
    si->work[i] -= si->cprime[i] * si->work[i + 1];
  }
  for (int row = 0; row < n; ++row) {
    Y[row] += alpha * si->sqrt_metric[row] * si->work[row];
  }
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
    if (block_cols == 1) {
      for (int64_t col = 0; col < csc->cols; ++col) {
        const double x_col = X[col];
        if (x_col == 0.0) continue;
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          Y[csc->row_idx[pos]] += alpha * csc->values[pos] * x_col;
        }
      }
    } else if (block_cols == 2) {
      double* y0 = Y;
      double* y1 = Y + ldy;
      const double* x0 = X;
      const double* x1 = X + ldx;
      for (int64_t col = 0; col < csc->cols; ++col) {
        const double x_col0 = x0[col];
        const double x_col1 = x1[col];
        if (x_col0 == 0.0 && x_col1 == 0.0) continue;
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          const int row = csc->row_idx[pos];
          const double a = alpha * csc->values[pos];
          y0[row] += a * x_col0;
          y1[row] += a * x_col1;
        }
      }
    } else if (block_cols <= 10) {
      double* yptr[10];
      const double* xptr[10];
      double xval[10];
      for (int64_t block = 0; block < block_cols; ++block) {
        yptr[block] = Y + block * ldy;
        xptr[block] = X + block * ldx;
      }
      for (int64_t col = 0; col < csc->cols; ++col) {
        bool all_zero = true;
        for (int64_t block = 0; block < block_cols; ++block) {
          xval[block] = xptr[block][col];
          all_zero = all_zero && xval[block] == 0.0;
        }
        if (all_zero) continue;
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          const int row = csc->row_idx[pos];
          const double a = alpha * csc->values[pos];
          for (int64_t block = 0; block < block_cols; ++block) {
            yptr[block][row] += a * xval[block];
          }
        }
      }
    } else {
      for (int64_t col = 0; col < csc->cols; ++col) {
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          const int row = csc->row_idx[pos];
          const double a = alpha * csc->values[pos];
          for (int64_t block = 0; block < block_cols; ++block) {
            Y[row + block * ldy] += a * X[col + block * ldx];
          }
        }
      }
    }
  } else {
    if (block_cols == 1) {
      for (int64_t col = 0; col < csc->cols; ++col) {
        double acc = 0.0;
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          acc += csc->values[pos] * X[csc->row_idx[pos]];
        }
        Y[col] += alpha * acc;
      }
    } else if (block_cols == 2) {
      double* y0 = Y;
      double* y1 = Y + ldy;
      const double* x0 = X;
      const double* x1 = X + ldx;
      for (int64_t col = 0; col < csc->cols; ++col) {
        double acc0 = 0.0;
        double acc1 = 0.0;
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          const int row = csc->row_idx[pos];
          const double a = csc->values[pos];
          acc0 += a * x0[row];
          acc1 += a * x1[row];
        }
        y0[col] += alpha * acc0;
        y1[col] += alpha * acc1;
      }
    } else if (block_cols <= 10) {
      double* yptr[10];
      const double* xptr[10];
      double acc[10];
      for (int64_t block = 0; block < block_cols; ++block) {
        yptr[block] = Y + block * ldy;
        xptr[block] = X + block * ldx;
      }
      for (int64_t col = 0; col < csc->cols; ++col) {
        for (int64_t block = 0; block < block_cols; ++block) {
          acc[block] = 0.0;
        }
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          const int row = csc->row_idx[pos];
          const double a = csc->values[pos];
          for (int64_t block = 0; block < block_cols; ++block) {
            acc[block] += a * xptr[block][row];
          }
        }
        for (int64_t block = 0; block < block_cols; ++block) {
          yptr[block][col] += alpha * acc[block];
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

extern "C" int eigencore_r_operator_apply(void* impl,
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
  if (op != EIGENCORE_TRANSPOSE_NONE) {
    return -1;
  }
  RApplyOperator* fn = static_cast<RApplyOperator*>(impl);
  if (!eigencore_int_indexable(fn->rows) ||
      !eigencore_int_indexable(block_cols) ||
      !eigencore_int_indexable(ldx) ||
      !eigencore_int_indexable(ldy)) {
    return -2;
  }
  const int n = static_cast<int>(fn->rows);
  const int cols = static_cast<int>(block_cols);
  if (ldx < n || ldy < n || cols < 1 || TYPEOF(fn->apply) != CLOSXP) {
    return -1;
  }

  SEXP X_ = PROTECT(allocMatrix(REALSXP, n, cols));
  SEXP Y_ = PROTECT(allocMatrix(REALSXP, n, cols));
  for (int col = 0; col < cols; ++col) {
    const double* x_col = X + static_cast<int64_t>(col) * ldx;
    double* x_dst = REAL(X_) + static_cast<int64_t>(col) * n;
    const double* y_col = Y + static_cast<int64_t>(col) * ldy;
    double* y_dst = REAL(Y_) + static_cast<int64_t>(col) * n;
    std::memcpy(x_dst, x_col, sizeof(double) * static_cast<size_t>(n));
    std::memcpy(y_dst, y_col, sizeof(double) * static_cast<size_t>(n));
  }
  SEXP alpha_ = PROTECT(ScalarReal(alpha));
  SEXP beta_ = PROTECT(ScalarReal(beta));
  SEXP call = PROTECT(lang5(fn->apply, X_, alpha_, beta_, Y_));
  SET_TAG(CDR(call), install("X"));
  SET_TAG(CDR(CDR(call)), install("alpha"));
  SET_TAG(CDR(CDR(CDR(call))), install("beta"));
  SET_TAG(CDR(CDR(CDR(CDR(call)))), install("Y"));

  int error_occurred = 0;
  SEXP out_ = PROTECT(R_tryEval(call, R_GlobalEnv, &error_occurred));
  if (error_occurred) {
    UNPROTECT(6);
    return -8;
  }
  SEXP dimY = getAttrib(out_, R_DimSymbol);
  if (!isReal(out_) || dimY == R_NilValue ||
      INTEGER(dimY)[0] != n || INTEGER(dimY)[1] != cols) {
    UNPROTECT(6);
    return -8;
  }
  for (int col = 0; col < cols; ++col) {
    const double* out_col = REAL(out_) + static_cast<int64_t>(col) * n;
    double* y_col = Y + static_cast<int64_t>(col) * ldy;
    std::memcpy(y_col, out_col, sizeof(double) * static_cast<size_t>(n));
  }
  UNPROTECT(6);
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
    eigencore_apply_status_error("dense block apply", status);
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
    eigencore_apply_status_error("CSC block apply", status);
  }

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_csc_centered_block_apply(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP row_means_, SEXP col_means_,
    SEXP rows_, SEXP columns_, SEXP X_, SEXP alpha_, SEXP beta_, SEXP Y_,
    SEXP transpose_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(row_means_) || !isReal(col_means_) ||
      !isLogical(rows_) || !isLogical(columns_) ||
      !isReal(X_) || !isReal(Y_) || !isLogical(transpose_)) {
    error("invalid centered CSC block apply inputs");
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
  const bool rows = LOGICAL(rows_)[0];
  const bool columns = LOGICAL(columns_)[0];
  const double alpha = REAL(alpha_)[0];
  const double beta = REAL(beta_)[0];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;

  if (xr != inner) {
    error("non-conformable X for centered CSC block apply");
  }
  if (yr != out_rows || yc != xc) {
    error("non-conformable Y for centered CSC block apply");
  }
  if (rows && LENGTH(row_means_) != m) {
    error("row_means length must equal CSC row dimension");
  }
  if (columns && LENGTH(col_means_) != n) {
    error("col_means length must equal CSC column dimension");
  }

  SEXP out_ = PROTECT(duplicate(Y_));
  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  const double* X = REAL(X_);
  const double* row_means = REAL(row_means_);
  const double* col_means = REAL(col_means_);
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
    eigencore_apply_status_error("centered CSC block apply", status);
  }

  for (int block_col = 0; block_col < xc; ++block_col) {
    const double* x_col = X + static_cast<int64_t>(block_col) * xr;
    double* out_col = out + static_cast<int64_t>(block_col) * out_rows;
    if (!transpose) {
      if (columns) {
        double correction = 0.0;
        for (int col = 0; col < n; ++col) {
          correction += col_means[col] * x_col[col];
        }
        correction *= alpha;
        for (int row = 0; row < m; ++row) {
          out_col[row] -= correction;
        }
      }
      if (rows) {
        double x_sum = 0.0;
        for (int col = 0; col < n; ++col) {
          x_sum += x_col[col];
        }
        x_sum *= alpha;
        for (int row = 0; row < m; ++row) {
          out_col[row] -= row_means[row] * x_sum;
        }
      }
    } else {
      if (columns) {
        double x_sum = 0.0;
        for (int row = 0; row < m; ++row) {
          x_sum += x_col[row];
        }
        x_sum *= alpha;
        for (int col = 0; col < n; ++col) {
          out_col[col] -= col_means[col] * x_sum;
        }
      }
      if (rows) {
        double correction = 0.0;
        for (int row = 0; row < m; ++row) {
          correction += row_means[row] * x_col[row];
        }
        correction *= alpha;
        for (int col = 0; col < n; ++col) {
          out_col[col] -= correction;
        }
      }
    }
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
    eigencore_apply_status_error("diagonal block apply", status);
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
    eigencore_apply_status_error("native no-allocation check apply", status);
  }
  return native_operator_workspace_counters(&workspace);
}

extern "C" SEXP eigencore_dense_apply_int_guard_check(void) {
  double scalar = 0.0;
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};

  DenseColumnMajorOperator too_many_rows = {
    static_cast<int64_t>(INT_MAX) + 1,
    1,
    &scalar
  };
  const int oversized_rows_status = eigencore_dense_apply(
    &too_many_rows, EIGENCORE_TRANSPOSE_NONE, 1,
    &scalar, static_cast<int64_t>(INT_MAX) + 1,
    1.0, 0.0,
    &scalar, static_cast<int64_t>(INT_MAX) + 1,
    &workspace
  );

  DenseColumnMajorOperator small = {1, 1, &scalar};
  const int oversized_block_status = eigencore_dense_apply(
    &small, EIGENCORE_TRANSPOSE_NONE, static_cast<int64_t>(INT_MAX) + 1,
    &scalar, 1,
    1.0, 0.0,
    &scalar, 1,
    &workspace
  );

  SEXP out_ = PROTECT(allocVector(INTSXP, 2));
  INTEGER(out_)[0] = oversized_rows_status;
  INTEGER(out_)[1] = oversized_block_status;
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("oversized_rows"));
  SET_STRING_ELT(names_, 1, mkChar("oversized_block_cols"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(2);
  return out_;
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
