#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <chrono>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <vector>
#include "eigencore_operator.h"

struct NativeBlockStageSeconds {
  double apply = 0.0;
  double recurrence = 0.0;
  double reorthogonalization = 0.0;
  double projected_solve = 0.0;
  double projection_update = 0.0;
  double projection_copy = 0.0;
  double projected_eigensolve = 0.0;
  double selected_vector_copy = 0.0;
  double ritz_residual = 0.0;
  double ritz_vector_form = 0.0;
  double ritz_operator_apply = 0.0;
  double ritz_norm = 0.0;
  double ritz_final_polish = 0.0;
  double locking = 0.0;
  double restart = 0.0;
};

struct NativeBlockRestartHistory {
  int capacity = 0;
  int length = 0;
  int* restart = nullptr;
  int* m_active = nullptr;
  int* selected_count = nullptr;
  int* locked_before = nullptr;
  int* locked_after = nullptr;
  int* nconv_wanted = nullptr;
  double* max_residual = nullptr;
  double* max_backward_error = nullptr;
};

static inline std::chrono::steady_clock::time_point native_timer_now() {
  return std::chrono::steady_clock::now();
}

static inline double native_timer_elapsed(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

static inline int eigencore_int_indexable(int64_t value) {
  return value >= 0 && value <= static_cast<int64_t>(INT_MAX);
}

static void eigencore_apply_status_error(const char* context, int status) {
  if (status == -2) {
    error("%s failed: dimensions exceed LP64 BLAS/R integer range; LAPACK64 is not enabled",
          context);
  }
  error("%s failed with status=%d", context, status);
}

static void combine_basis_columns_small(const double* basis,
                                        int n,
                                        int basis_cols,
                                        const double* coeff,
                                        int coeff_ld,
                                        int out_cols,
                                        double* out) {
  std::memset(out, 0,
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(out_cols));
  for (int p = 0; p < out_cols; ++p) {
    double* y = out + static_cast<int64_t>(p) * n;
    for (int col = 0; col < basis_cols; ++col) {
      const double a = coeff[col + static_cast<int64_t>(p) * coeff_ld];
      if (a == 0.0) {
        continue;
      }
      const double* x = basis + static_cast<int64_t>(col) * n;
      for (int row = 0; row < n; ++row) {
        y[row] += a * x[row];
      }
    }
  }
}

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

struct RApplyOperator {
  int64_t rows;
  SEXP apply;
};

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

static int eigencore_r_operator_apply(void* impl,
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
  return workspace_counters(&workspace);
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

static double frobenius_norm_dense(const double* X, int len) {
  long double sum = 0.0L;
  for (int i = 0; i < len; ++i) {
    const long double value = X[i];
    sum += value * value;
  }
  return sqrt(static_cast<double>(sum));
}

static bool ritz_value_better(double candidate, double incumbent, int target_kind) {
  switch (target_kind) {
    case 2:
      return candidate < incumbent;
    case 3:
      return fabs(candidate) > fabs(incumbent);
    case 4:
      return fabs(candidate) < fabs(incumbent);
    case 1:
    default:
      return candidate > incumbent;
  }
}

static int selected_ritz_indices(const double* values,
                                 int n,
                                 int k,
                                 int target_kind,
                                 int* selected) {
  const int count = (k < n) ? k : n;
  std::vector<bool> taken(static_cast<size_t>(n), false);
  for (int i = 0; i < count; ++i) {
    int best = -1;
    for (int j = 0; j < n; ++j) {
      if (taken[static_cast<size_t>(j)]) {
        continue;
      }
      if (best < 0 || ritz_value_better(values[j], values[best], target_kind)) {
        best = j;
      }
    }
    selected[i] = best;
    if (best >= 0) {
      taken[static_cast<size_t>(best)] = true;
    }
  }
  return count;
}

static int selected_sorted_ritz_indices(const double* values,
                                        int n,
                                        int k,
                                        int target_kind,
                                        int* selected) {
  const int count = (k < n) ? k : n;
  if (target_kind == 1) {
    for (int i = 0; i < count; ++i) {
      selected[i] = n - 1 - i;
    }
    return count;
  }
  if (target_kind == 2) {
    for (int i = 0; i < count; ++i) {
      selected[i] = i;
    }
    return count;
  }
  return selected_ritz_indices(values, n, k, target_kind, selected);
}

static int lanczos_convergence_estimate(const double* alpha,
                                        const double* beta,
                                        int iter,
                                        int k,
                                        int target_kind,
                                        double tol,
                                        int* nconv,
                                        double* max_residual) {
  *nconv = 0;
  *max_residual = R_PosInf;
  if (iter <= 0 || k <= 0) {
    return 0;
  }

  double* diag = static_cast<double*>(std::calloc(static_cast<size_t>(iter), sizeof(double)));
  double* offdiag = (iter > 1)
    ? static_cast<double*>(std::calloc(static_cast<size_t>(iter - 1), sizeof(double)))
    : nullptr;
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  double* work = (iter > 1)
    ? static_cast<double*>(std::calloc(static_cast<size_t>(2 * iter - 2), sizeof(double)))
    : nullptr;
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>((k < iter) ? k : iter), sizeof(int)));
  if (diag == nullptr || z == nullptr || selected == nullptr ||
      (iter > 1 && (offdiag == nullptr || work == nullptr))) {
    std::free(diag);
    std::free(offdiag);
    std::free(z);
    std::free(work);
    std::free(selected);
    return -2;
  }

  for (int i = 0; i < iter; ++i) {
    diag[i] = alpha[i];
  }
  for (int i = 0; i < iter - 1; ++i) {
    offdiag[i] = beta[i];
  }

  if (iter == 1) {
    z[0] = 1.0;
  } else {
    char jobz = 'V';
    int info = 0;
    F77_CALL(dstev)(&jobz, &iter, diag, offdiag, z, &iter, work, &info FCONE);
    if (info != 0) {
      std::free(diag);
      std::free(offdiag);
      std::free(z);
      std::free(work);
      std::free(selected);
      return info;
    }
  }

  const int count = selected_ritz_indices(diag, iter, k, target_kind, selected);
  double max_selected = 0.0;
  int conv = 0;
  for (int i = 0; i < count; ++i) {
    const int idx = selected[i];
    const double residual = fabs(beta[iter - 1] * z[(iter - 1) + idx * iter]);
    const double threshold = tol * ((fabs(diag[idx]) > 1.0) ? fabs(diag[idx]) : 1.0);
    if (residual <= threshold) {
      ++conv;
    }
    if (residual > max_selected) {
      max_selected = residual;
    }
  }

  *nconv = conv;
  *max_residual = max_selected;
  std::free(diag);
  std::free(offdiag);
  std::free(z);
  std::free(work);
  std::free(selected);
  return 0;
}

static int golub_kahan_projected_convergence_estimate(const double* alpha,
                                                      const double* beta,
                                                      int iter,
                                                      int k,
                                                      int target_kind,
                                                      double tol,
                                                      int* nconv,
                                                      double* max_residual) {
  *nconv = 0;
  *max_residual = R_PosInf;
  if (iter <= 0 || k <= 0) {
    return 0;
  }

  double* B = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  double* d = static_cast<double*>(std::calloc(static_cast<size_t>(iter), sizeof(double)));
  double* u = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  double* vt = static_cast<double*>(std::calloc(static_cast<size_t>(iter) * iter, sizeof(double)));
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>((k < iter) ? k : iter), sizeof(int)));
  if (B == nullptr || d == nullptr || u == nullptr || vt == nullptr || selected == nullptr) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(selected);
    return -2;
  }

  for (int i = 0; i < iter; ++i) {
    B[i + i * iter] = alpha[i];
  }
  for (int i = 0; i < iter - 1; ++i) {
    B[i + (i + 1) * iter] = beta[i];
  }

  char jobu = 'A';
  char jobvt = 'N';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgesvd)(&jobu, &jobvt, &iter, &iter, B, &iter,
                   d, u, &iter, vt, &iter, &work_query, &lwork,
                   &info FCONE FCONE);
  if (info != 0) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(selected);
    return info;
  }
  lwork = static_cast<int>(work_query);
  double* work = static_cast<double*>(std::calloc(static_cast<size_t>(lwork), sizeof(double)));
  if (work == nullptr) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(selected);
    return -2;
  }

  std::memset(B, 0, sizeof(double) * static_cast<size_t>(iter) * iter);
  for (int i = 0; i < iter; ++i) {
    B[i + i * iter] = alpha[i];
  }
  for (int i = 0; i < iter - 1; ++i) {
    B[i + (i + 1) * iter] = beta[i];
  }
  F77_CALL(dgesvd)(&jobu, &jobvt, &iter, &iter, B, &iter,
                   d, u, &iter, vt, &iter, work, &lwork,
                   &info FCONE FCONE);
  if (info != 0) {
    std::free(B);
    std::free(d);
    std::free(u);
    std::free(vt);
    std::free(work);
    std::free(selected);
    return info;
  }

  const int count = selected_ritz_indices(d, iter, k, target_kind, selected);
  double max_selected = 0.0;
  int conv = 0;
  for (int i = 0; i < count; ++i) {
    const int idx = selected[i];
    const double residual = fabs(beta[iter - 1] * u[(iter - 1) + idx * iter]);
    const double scale = (fabs(d[idx]) > 1.0) ? fabs(d[idx]) : 1.0;
    const double threshold = tol * scale;
    if (residual <= threshold) {
      ++conv;
    }
    if (residual > max_selected) {
      max_selected = residual;
    }
  }

  *nconv = conv;
  *max_residual = max_selected;
  std::free(B);
  std::free(d);
  std::free(u);
  std::free(vt);
  std::free(work);
  std::free(selected);
  return 0;
}

static int native_lanczos_run(void* impl,
                              EigencoreApplyFn apply,
                              int n,
                              int maxit,
                              int k,
                              int target_kind,
                              double tol,
                              const double* start,
                              double* Q,
                              double* alpha,
                              double* beta,
                              int* history_nconv,
                              double* history_max_residual,
                              int* iterations,
                              int* matvecs) {
  double* q = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* q_prev = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  if (q == nullptr || q_prev == nullptr || z == nullptr) {
    std::free(q);
    std::free(q_prev);
    std::free(z);
    return -2;
  }

  long double start_norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    start_norm2 += static_cast<long double>(start[row]) * start[row];
  }
  double q_norm = sqrt(static_cast<double>(start_norm2));
  if (q_norm == 0.0) {
    q[0] = 1.0;
  } else {
    for (int row = 0; row < n; ++row) {
      q[row] = start[row] / q_norm;
    }
  }

  *iterations = 0;
  *matvecs = 0;
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  for (int j = 0; j < maxit; ++j) {
    *iterations = j + 1;
    std::memcpy(Q + j * n, q, sizeof(double) * static_cast<size_t>(n));
    std::memset(z, 0, sizeof(double) * static_cast<size_t>(n));

    const int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, q, n,
                             1.0, 0.0, z, n, &workspace);
    if (status != 0) {
      std::free(q);
      std::free(q_prev);
      std::free(z);
      return status;
    }
    ++(*matvecs);

    if (j > 0) {
      const double bj = beta[j - 1];
      for (int row = 0; row < n; ++row) {
        z[row] -= bj * q_prev[row];
      }
    }

    long double aj = 0.0L;
    for (int row = 0; row < n; ++row) {
      aj += static_cast<long double>(q[row]) * z[row];
    }
    alpha[j] = static_cast<double>(aj);
    for (int row = 0; row < n; ++row) {
      z[row] -= alpha[j] * q[row];
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int prev = 0; prev <= j; ++prev) {
        const double* qprev_basis = Q + prev * n;
        long double dot = 0.0L;
        for (int row = 0; row < n; ++row) {
          dot += static_cast<long double>(qprev_basis[row]) * z[row];
        }
        const double coeff = static_cast<double>(dot);
        for (int row = 0; row < n; ++row) {
          z[row] -= coeff * qprev_basis[row];
        }
      }
    }

    long double beta_norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      beta_norm2 += static_cast<long double>(z[row]) * z[row];
    }
    beta[j] = sqrt(static_cast<double>(beta_norm2));
    history_nconv[j] = 0;
    history_max_residual[j] = R_PosInf;
    if (j + 1 >= k) {
      int nconv = 0;
      double max_residual = R_PosInf;
      const int conv_status = lanczos_convergence_estimate(
        alpha, beta, j + 1, k, target_kind, tol, &nconv, &max_residual
      );
      if (conv_status != 0) {
        std::free(q);
        std::free(q_prev);
        std::free(z);
        return conv_status;
      }
      history_nconv[j] = nconv;
      history_max_residual[j] = max_residual;
      if (nconv >= k) {
        break;
      }
    }
    if (j + 1 == maxit || beta[j] <= 100.0 * DBL_EPSILON) {
      break;
    }

    std::memcpy(q_prev, q, sizeof(double) * static_cast<size_t>(n));
    for (int row = 0; row < n; ++row) {
      q[row] = z[row] / beta[j];
    }
  }

  std::free(q);
  std::free(q_prev);
  std::free(z);
  return 0;
}

static int native_golub_kahan_run(void* impl,
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
                                  int reorthogonalize_v) {
  double* v = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* u = static_cast<double*>(std::calloc(static_cast<size_t>(m), sizeof(double)));
  double* u_prev = static_cast<double*>(std::calloc(static_cast<size_t>(m), sizeof(double)));
  double* coeff = use_blas_reorthogonalization
    ? static_cast<double*>(std::calloc(static_cast<size_t>(maxit), sizeof(double)))
    : nullptr;
  if (v == nullptr || z == nullptr || u == nullptr || u_prev == nullptr || coeff == nullptr) {
    if (!use_blas_reorthogonalization && v != nullptr && z != nullptr &&
        u != nullptr && u_prev != nullptr) {
      // coeff is intentionally unused on the scalar sparse path.
    } else {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      std::free(coeff);
      return -2;
    }
  }
  if (v == nullptr || z == nullptr || u == nullptr || u_prev == nullptr) {
    std::free(v);
    std::free(z);
    std::free(u);
    std::free(u_prev);
    std::free(coeff);
    return -2;
  }

  long double start_norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    start_norm2 += static_cast<long double>(start[row]) * start[row];
  }
  double v_norm = sqrt(static_cast<double>(start_norm2));
  if (v_norm == 0.0) {
    v[0] = 1.0;
  } else {
    for (int row = 0; row < n; ++row) {
      v[row] = start[row] / v_norm;
    }
  }

  *iterations = 0;
  *matvecs = 0;
  *projected_stop = 0;
  *projected_nconv = 0;
  *projected_max_residual = R_PosInf;
  *projected_checks = 0;
  *projected_seconds = 0.0;
  *stage_apply_seconds = 0.0;
  *stage_recurrence_seconds = 0.0;
  *stage_reorthogonalization_seconds = 0.0;
  *reorthogonalization_passes = 0;
  double beta_prev = 0.0;
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const int check_interval = (2 * k > 10) ? 2 * k : 10;
  const int min_projected_savings = (k > 5) ? k : 5;

  for (int j = 0; j < maxit; ++j) {
    *iterations = j + 1;
    std::memcpy(V + j * n, v, sizeof(double) * static_cast<size_t>(n));
    std::memset(u, 0, sizeof(double) * static_cast<size_t>(m));

    auto stage_timer = native_timer_now();
    int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, v, n,
                       1.0, 0.0, u, m, &workspace);
    *stage_apply_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      std::free(coeff);
      return status;
    }
    ++(*matvecs);
    stage_timer = native_timer_now();
    if (j > 0) {
      for (int row = 0; row < m; ++row) {
        u[row] -= beta_prev * u_prev[row];
      }
    }
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    if (reorthogonalize_u && j > 0) {
      double norm_before2 = 0.0;
      for (int row = 0; row < m; ++row) {
        norm_before2 += u[row] * u[row];
      }
      const double norm_before = sqrt(norm_before2);
      const char trans_T = 'T';
      const char trans_N = 'N';
      const int one_col = 1;
      const double one = 1.0;
      const double zero = 0.0;
      const double minus_one = -1.0;
      int passes = 1;
      for (int pass = 0; pass < passes; ++pass) {
        if (use_blas_reorthogonalization) {
        F77_CALL(dgemv)(&trans_T, &m, &j, &one,
                        U, &m, u, &one_col,
                        &zero, coeff, &one_col FCONE);
        F77_CALL(dgemv)(&trans_N, &m, &j, &minus_one,
                        U, &m, coeff, &one_col,
                        &one, u, &one_col FCONE);
        } else {
          for (int prev = 0; prev < j; ++prev) {
            const double* uprev_basis = U + prev * m;
            double dot = 0.0;
            for (int row = 0; row < m; ++row) {
              dot += uprev_basis[row] * u[row];
            }
            const double scalar_coeff = dot;
            for (int row = 0; row < m; ++row) {
              u[row] -= scalar_coeff * uprev_basis[row];
            }
          }
        }
        ++(*reorthogonalization_passes);
        if (pass == 0 && norm_before > 0.0) {
          double norm_after2 = 0.0;
          for (int row = 0; row < m; ++row) {
            norm_after2 += u[row] * u[row];
          }
          const double norm_after = sqrt(norm_after2);
          if (norm_after < 0.717 * norm_before) {
            passes = 2;
          }
        }
      }
    }
    *stage_reorthogonalization_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    long double alpha_norm2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      alpha_norm2 += static_cast<long double>(u[row]) * u[row];
    }
    alpha[j] = sqrt(static_cast<double>(alpha_norm2));
    if (alpha[j] <= 100.0 * DBL_EPSILON) {
      break;
    }
    for (int row = 0; row < m; ++row) {
      u[row] /= alpha[j];
    }
    std::memcpy(U + j * m, u, sizeof(double) * static_cast<size_t>(m));
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    std::memset(z, 0, sizeof(double) * static_cast<size_t>(n));
    stage_timer = native_timer_now();
    status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, 1, u, m,
                   1.0, 0.0, z, n, &workspace);
    *stage_apply_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      std::free(coeff);
      return status;
    }
    ++(*matvecs);
    stage_timer = native_timer_now();
    for (int row = 0; row < n; ++row) {
      z[row] -= alpha[j] * v[row];
    }
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    if (reorthogonalize_v) {
      const int active_v = j + 1;
      double norm_before2 = 0.0;
      for (int row = 0; row < n; ++row) {
        norm_before2 += z[row] * z[row];
      }
      const double norm_before = sqrt(norm_before2);
      const char trans_T = 'T';
      const char trans_N = 'N';
      const int one_col = 1;
      const double one = 1.0;
      const double zero = 0.0;
      const double minus_one = -1.0;
      int passes = 1;
      for (int pass = 0; pass < passes; ++pass) {
        if (use_blas_reorthogonalization) {
        F77_CALL(dgemv)(&trans_T, &n, &active_v, &one,
                        V, &n, z, &one_col,
                        &zero, coeff, &one_col FCONE);
        F77_CALL(dgemv)(&trans_N, &n, &active_v, &minus_one,
                        V, &n, coeff, &one_col,
                        &one, z, &one_col FCONE);
        } else {
          for (int prev = 0; prev <= j; ++prev) {
            const double* vprev_basis = V + prev * n;
            double dot = 0.0;
            for (int row = 0; row < n; ++row) {
              dot += vprev_basis[row] * z[row];
            }
            const double scalar_coeff = dot;
            for (int row = 0; row < n; ++row) {
              z[row] -= scalar_coeff * vprev_basis[row];
            }
          }
        }
        ++(*reorthogonalization_passes);
        if (pass == 0 && norm_before > 0.0) {
          double norm_after2 = 0.0;
          for (int row = 0; row < n; ++row) {
            norm_after2 += z[row] * z[row];
          }
          const double norm_after = sqrt(norm_after2);
          if (norm_after < 0.717 * norm_before) {
            passes = 2;
          }
        }
      }
    }
    *stage_reorthogonalization_seconds += native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    long double beta_norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      beta_norm2 += static_cast<long double>(z[row]) * z[row];
    }
    beta[j] = sqrt(static_cast<double>(beta_norm2));
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);
    const int remaining_iterations = maxit - (j + 1);
    if (enable_projected_stop && j + 1 >= k &&
        remaining_iterations >= min_projected_savings &&
        ((j + 1 == maxit) || ((j + 1) % check_interval == 0))) {
      int nconv = 0;
      double max_residual = R_PosInf;
      auto projected_timer = native_timer_now();
      const int conv_status = golub_kahan_projected_convergence_estimate(
        alpha, beta, j + 1, k, target_kind, tol, &nconv, &max_residual
      );
      *projected_seconds += native_timer_elapsed(projected_timer);
      ++(*projected_checks);
      if (conv_status != 0) {
        std::free(v);
        std::free(z);
        std::free(u);
        std::free(u_prev);
        std::free(coeff);
        return conv_status;
      }
      *projected_nconv = nconv;
      *projected_max_residual = max_residual;
      if (nconv >= k) {
        *projected_stop = 1;
        break;
      }
    }
    if (j + 1 == maxit || beta[j] <= 100.0 * DBL_EPSILON) {
      break;
    }

    stage_timer = native_timer_now();
    std::memcpy(u_prev, u, sizeof(double) * static_cast<size_t>(m));
    beta_prev = beta[j];
    for (int row = 0; row < n; ++row) {
      v[row] = z[row] / beta[j];
    }
    *stage_recurrence_seconds += native_timer_elapsed(stage_timer);
  }

  std::free(v);
  std::free(z);
  std::free(u);
  std::free(u_prev);
  std::free(coeff);
  return 0;
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

extern "C" SEXP eigencore_lanczos_dense(SEXP A_, SEXP maxit_, SEXP start_,
                                        SEXP k_, SEXP target_kind_, SEXP tol_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || LENGTH(start_) != n) {
    error("non-conformable dense Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and nrow(A)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(&impl, eigencore_dense_apply, n, maxit,
                                        k, target_kind, tol,
                                        REAL(start_), REAL(Q_), REAL(alpha_),
                                        REAL(beta_), INTEGER(history_nconv_),
                                        REAL(history_max_residual_),
                                        &iterations, &matvecs);
  if (status != 0) {
    error("native dense Lanczos failed with status=%d", status);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_lanczos_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                      SEXP maxit_, SEXP start_, SEXP k_,
                                      SEXP target_kind_, SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC Lanczos inputs");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || LENGTH(start_) != n) {
    error("non-conformable CSC Lanczos inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  if (maxit < 1 || maxit > n) {
    error("maxit must be between 1 and nrow(A)");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || k > maxit) {
    error("k must be between 1 and maxit");
  }

  SEXP Q_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP history_nconv_ = PROTECT(allocVector(INTSXP, maxit));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(Q_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(INTEGER(history_nconv_), 0, sizeof(int) * static_cast<size_t>(maxit));
  for (int i = 0; i < maxit; ++i) {
    REAL(history_max_residual_)[i] = R_PosInf;
  }

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_lanczos_run(&impl, eigencore_csc_apply, n, maxit,
                                        k, target_kind, tol,
                                        REAL(start_), REAL(Q_), REAL(alpha_),
                                        REAL(beta_), INTEGER(history_nconv_),
                                        REAL(history_max_residual_),
                                        &iterations, &matvecs);
  if (status != 0) {
    error("native CSC Lanczos failed with status=%d", status);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out_, 0, Q_);
  SET_VECTOR_ELT(out_, 1, alpha_);
  SET_VECTOR_ELT(out_, 2, beta_);
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, history_nconv_);
  SET_VECTOR_ELT(out_, 6, history_max_residual_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(names_, 0, mkChar("Q"));
  SET_STRING_ELT(names_, 1, mkChar("alpha"));
  SET_STRING_ELT(names_, 2, mkChar("beta"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("history_max_residual"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_dense(SEXP A_, SEXP maxit_, SEXP start_,
                                            SEXP rank_, SEXP target_kind_,
                                            SEXP tol_, SEXP projected_stop_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  if (LENGTH(start_) != n) {
    error("non-conformable dense Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }
  if (rank < 1 || rank > maxit) {
    error("rank must be between 1 and maxit");
  }

  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(maxit), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(maxit), 0.0);
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(maxit) +
      static_cast<int64_t>(2 * m + 2 * n) +
      static_cast<int64_t>(maxit)
    );

  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int iterations = 0;
  int matvecs = 0;
  int projected_stop = 0;
  int projected_nconv = 0;
  double projected_max_residual = R_PosInf;
  int projected_checks = 0;
  double projected_seconds = 0.0;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  int reorthogonalization_passes = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_dense_apply, m, n,
                                            maxit, rank, target_kind, tol,
                                            enable_projected_stop, 1,
                                            REAL(start_),
                                            U_work.data(), V_work.data(),
                                            alpha_work.data(), beta_work.data(),
                                            &iterations, &matvecs,
                                            &projected_stop, &projected_nconv,
                                            &projected_max_residual,
                                            &projected_checks,
                                            &projected_seconds,
                                            &stage_apply_seconds,
                                            &stage_recurrence_seconds,
                                            &stage_reorthogonalization_seconds,
                                            &reorthogonalization_passes,
                                            1, 1);
  if (status != 0) {
    error("native dense Golub-Kahan failed with status=%d", status);
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, iterations));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, iterations));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP beta_ = PROTECT(allocVector(REALSXP, iterations));
  for (int col = 0; col < iterations; ++col) {
    std::memcpy(REAL(U_) + static_cast<size_t>(col) * static_cast<size_t>(m),
                U_work.data() + static_cast<size_t>(col) * static_cast<size_t>(m),
                sizeof(double) * static_cast<size_t>(m));
    std::memcpy(REAL(V_) + static_cast<size_t>(col) * static_cast<size_t>(n),
                V_work.data() + static_cast<size_t>(col) * static_cast<size_t>(n),
                sizeof(double) * static_cast<size_t>(n));
  }
  std::memcpy(REAL(alpha_), alpha_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));
  std::memcpy(REAL(beta_), beta_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 8, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 10, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalization_passes));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 7, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 8, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 9, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 10, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 12, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalization_passes"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                          SEXP maxit_, SEXP start_,
                                          SEXP rank_, SEXP target_kind_,
                                          SEXP tol_, SEXP projected_stop_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC Golub-Kahan inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  if (LENGTH(start_) != n) {
    error("non-conformable CSC Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }
  if (rank < 1 || rank > maxit) {
    error("rank must be between 1 and maxit");
  }

  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(maxit), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(maxit), 0.0);
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(maxit) +
      static_cast<int64_t>(2 * m + 2 * n)
    );

  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int iterations = 0;
  int matvecs = 0;
  int projected_stop = 0;
  int projected_nconv = 0;
  double projected_max_residual = R_PosInf;
  int projected_checks = 0;
  double projected_seconds = 0.0;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  int reorthogonalization_passes = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_csc_apply, m, n,
                                            maxit, rank, target_kind, tol,
                                            enable_projected_stop, 0,
                                            REAL(start_),
                                            U_work.data(), V_work.data(),
                                            alpha_work.data(), beta_work.data(),
                                            &iterations, &matvecs,
                                            &projected_stop, &projected_nconv,
                                            &projected_max_residual,
                                            &projected_checks,
                                            &projected_seconds,
                                            &stage_apply_seconds,
                                            &stage_recurrence_seconds,
                                            &stage_reorthogonalization_seconds,
                                            &reorthogonalization_passes,
                                            1, 1);
  if (status != 0) {
    error("native CSC Golub-Kahan failed with status=%d", status);
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, iterations));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, iterations));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP beta_ = PROTECT(allocVector(REALSXP, iterations));
  for (int col = 0; col < iterations; ++col) {
    std::memcpy(REAL(U_) + static_cast<size_t>(col) * static_cast<size_t>(m),
                U_work.data() + static_cast<size_t>(col) * static_cast<size_t>(m),
                sizeof(double) * static_cast<size_t>(m));
    std::memcpy(REAL(V_) + static_cast<size_t>(col) * static_cast<size_t>(n),
                V_work.data() + static_cast<size_t>(col) * static_cast<size_t>(n),
                sizeof(double) * static_cast<size_t>(n));
  }
  std::memcpy(REAL(alpha_), alpha_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));
  std::memcpy(REAL(beta_), beta_work.data(),
              sizeof(double) * static_cast<size_t>(iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 8, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 10, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalization_passes));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 7, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 8, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 9, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 10, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 12, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalization_passes"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}

// =====================================================================
// Thick-restart Hermitian Lanczos with locking
//
// Implements a Krylov-Schur-style restarted Lanczos for the symmetric /
// Hermitian standard eigenproblem A x = lambda x. Block size is one
// (single-vector Lanczos); BLAS-3 block extension is a follow-on.
//
// Algorithm outline:
//   - Build an active Krylov basis V (n x m_max) with full
//     reorthogonalisation (DGKS x2) against both the locked basis and
//     the active basis.
//   - Cache AV = A * V column by column so Rayleigh-Ritz can run
//     entirely on stored buffers.
//   - At the end of each inner sweep, project H = V^T A V (m x m),
//     symmetrise, dsyev to get Ritz pairs (theta, S).
//   - Compute Ritz vectors B_v = V S, B_av = AV S; per-Ritz residual
//     norms ||B_av_i - theta_i B_v_i||.
//   - Lock the wanted Ritz pairs whose residual <= tol * max(|theta|, 1).
//   - If still unconverged, contract V to the top k_keep unlocked Ritz
  //     vectors (k_keep <= m_max - 2), append a Ritz-residual continuation
  //     tail orthogonalised against locked + retained, and run the next
  //     inner sweep. A deterministic random tail is only a fallback when all
  //     usable residuals have collapsed.
// =====================================================================

static int trl_orthogonalise(const double* V_locked, int n_locked,
                             const double* V_active, int m_active,
                             double* z, double* tmp, int n, int passes = 2) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  int incx = 1;
  for (int pass = 0; pass < passes; ++pass) {
    if (n_locked > 0) {
      F77_CALL(dgemv)(&trans_T, &n, &n_locked, &one,
                      V_locked, &n, z, &incx,
                      &zero, tmp, &incx FCONE);
      F77_CALL(dgemv)(&trans_N, &n, &n_locked, &minus_one,
                      V_locked, &n, tmp, &incx,
                      &one, z, &incx FCONE);
    }
    if (m_active > 0) {
      F77_CALL(dgemv)(&trans_T, &n, &m_active, &one,
                      V_active, &n, z, &incx,
                      &zero, tmp, &incx FCONE);
      F77_CALL(dgemv)(&trans_N, &n, &m_active, &minus_one,
                      V_active, &n, tmp, &incx,
                      &one, z, &incx FCONE);
    }
  }
  return 0;
}

static double trl_norm2(const double* x, int n) {
  long double sum = 0.0L;
  for (int i = 0; i < n; ++i) {
    sum += static_cast<long double>(x[i]) * x[i];
  }
  return sqrt(static_cast<double>(sum));
}

static int trl_dsyev_query(int m_max) {
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork_query = -1;
  double work_query = 0.0;
  double fake = 0.0;
  double fake_w = 0.0;
  int m = m_max;
  F77_CALL(dsyev)(&jobz, &uplo, &m, &fake, &m, &fake_w,
                  &work_query, &lwork_query, &info FCONE FCONE);
  if (info != 0) {
    return 3 * m_max;
  }
  return static_cast<int>(work_query);
}

static int trl_dsyevd_query(int m_max, int* liwork_out) {
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork_query = -1;
  int liwork_query = -1;
  double work_query = 0.0;
  int iwork_query = 0;
  double fake = 0.0;
  double fake_w = 0.0;
  int m = m_max;
  F77_CALL(dsyevd)(&jobz, &uplo, &m, &fake, &m, &fake_w,
                   &work_query, &lwork_query,
                   &iwork_query, &liwork_query, &info FCONE FCONE);
  if (info != 0) {
    if (liwork_out != nullptr) {
      *liwork_out = 3 + 5 * m_max;
    }
    return 1 + 6 * m_max + 2 * m_max * m_max;
  }
  if (liwork_out != nullptr) {
    *liwork_out = iwork_query;
  }
  return static_cast<int>(work_query);
}

static int symmetric_eigen_inplace(double* A, int n, double* values,
                                   double* work, int lwork,
                                   int* iwork = nullptr,
                                   int liwork = 0) {
  if (n <= 0) {
    return 0;
  }
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  if (n >= 96 && iwork != nullptr && liwork > 0) {
    F77_CALL(dsyevd)(&jobz, &uplo, &n, A, &n, values,
                     work, &lwork, iwork, &liwork, &info FCONE FCONE);
  } else {
    F77_CALL(dsyev)(&jobz, &uplo, &n, A, &n, values,
                    work, &lwork, &info FCONE FCONE);
  }
  return info == 0 ? 0 : -3;
}

static void symmetrize_packed_square(double* A, int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double avg = 0.5 * (A[i + j * n] + A[j + i * n]);
      A[i + j * n] = avg;
      A[j + i * n] = avg;
    }
  }
}

static double standard_eigen_lock_scale(double norm_a, double theta,
                                        const double* v, int n) {
  if (!std::isfinite(norm_a) || norm_a <= 0.0) {
    norm_a = 1.0;
  }
  const double vnorm = trl_norm2(v, n);
  const double scale = (norm_a + fabs(theta)) *
    ((vnorm > DBL_EPSILON) ? vnorm : DBL_EPSILON);
  return (scale > DBL_EPSILON) ? scale : DBL_EPSILON;
}

static int vector_is_independent_from_locked(const double* V_locked,
                                             int n_locked,
                                             const double* v,
                                             int n) {
  const double dot_tol = 10.0 * sqrt(DBL_EPSILON);
  const double vnorm = trl_norm2(v, n);
  if (vnorm <= DBL_EPSILON) {
    return 0;
  }
  for (int col = 0; col < n_locked; ++col) {
    long double dot = 0.0L;
    long double locked_ss = 0.0L;
    const double* locked = V_locked + static_cast<int64_t>(col) * n;
    for (int row = 0; row < n; ++row) {
      dot += static_cast<long double>(locked[row]) * v[row];
      locked_ss += static_cast<long double>(locked[row]) * locked[row];
    }
    const double locked_norm = sqrt(static_cast<double>(locked_ss));
    if (locked_norm <= DBL_EPSILON) {
      continue;
    }
    if (fabs(static_cast<double>(dot)) > dot_tol * locked_norm * vnorm) {
      return 0;
    }
  }
  return 1;
}

static int block_accept_work_vector(const double* V_locked, int n_locked,
                                    double* V_active, int* m_active,
                                    int m_max, double* z, double* tmp,
                                    int n, int* ortho_passes) {
  if (*m_active >= m_max) {
    return 0;
  }
  const int passes = 2;
  trl_orthogonalise(V_locked, n_locked, V_active, *m_active, z, tmp, n, passes);
  if (ortho_passes != nullptr) {
    *ortho_passes += passes;
  }
  const double nz = trl_norm2(z, n);
  if (nz <= 100.0 * DBL_EPSILON) {
    return 0;
  }
  const double inv_nz = 1.0 / nz;
  double* dst = V_active + static_cast<int64_t>(*m_active) * n;
  for (int row = 0; row < n; ++row) {
    dst[row] = z[row] * inv_nz;
  }
  ++(*m_active);
  return 1;
}

static void block_reorthogonalise_against(const double* V_locked, int n_locked,
                                          const double* V_active, int m_active,
                                          double* X, int n, int cols,
                                          double* coeff, int passes) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  for (int pass = 0; pass < passes; ++pass) {
    if (n_locked > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &n_locked, &cols, &n,
                      &one, V_locked, &n, X, &n,
                      &zero, coeff, &n_locked FCONE FCONE);
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &cols, &n_locked,
                      &minus_one, V_locked, &n, coeff, &n_locked,
                      &one, X, &n FCONE FCONE);
    }
    if (m_active > 0) {
      F77_CALL(dgemm)(&trans_T, &trans_N, &m_active, &cols, &n,
                      &one, V_active, &n, X, &n,
                      &zero, coeff, &m_active FCONE FCONE);
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &cols, &m_active,
                      &minus_one, V_active, &n, coeff, &m_active,
                      &one, X, &n FCONE FCONE);
    }
  }
}

// Contract: X and Z_block MAY alias (X == Z_block && ldx == n is a supported
// caller pattern — the per-column memcpy below would otherwise be a self-copy,
// which is UB in standard C++ even when the bytes happen to overlap exactly).
// When X aliases Z_block at the same leading dimension, the loader is a no-op
// and we skip the memcpy. Distinct buffers (or X with a different stride) get
// the explicit column copy. If callers ever start passing partially-aliased
// regions (e.g. X == Z_block + offset) this guard is insufficient and the
// caller must materialize a non-aliased temporary first.
static int block_accept_columns_blas3(const double* X, int ldx, int x_cols,
                                      const double* V_locked, int n_locked,
                                      double* V_active, int* m_active,
                                      int m_max, double* Z_block,
                                      int block_capacity, double* coeff,
                                      double* tmp, int n, int max_accept,
                                      int* ortho_passes,
                                      bool reorthogonalize_active = true) {
  if (max_accept < 0) {
    max_accept = 0;
  }
  int cols = x_cols;
  if (cols > max_accept) {
    cols = max_accept;
  }
  if (cols > block_capacity) {
    cols = block_capacity;
  }
  if (cols > m_max - *m_active) {
    cols = m_max - *m_active;
  }
  if (cols <= 0) {
    return 0;
  }
  const bool x_aliases_z = (X == Z_block) && (ldx == n);
  if (!x_aliases_z) {
    for (int col = 0; col < cols; ++col) {
      std::memcpy(Z_block + static_cast<int64_t>(col) * n,
                  X + static_cast<int64_t>(col) * ldx,
                  sizeof(double) * static_cast<size_t>(n));
    }
  }

  const int reorth_passes = 2;
  const double* active_basis = reorthogonalize_active ? V_active : nullptr;
  const int active_cols = reorthogonalize_active ? *m_active : 0;
  block_reorthogonalise_against(V_locked, n_locked, active_basis, active_cols,
                                Z_block, n, cols, coeff, reorth_passes);
  if (ortho_passes != nullptr) {
    *ortho_passes += reorth_passes;
  }

  if (cols > 0) {
    const char trans_T = 'T';
    const char trans_N = 'N';
    const char right = 'R';
    const char uplo = 'U';
    const char diag = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    int info = 0;
    F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                    &one, Z_block, &n, Z_block, &n,
                    &zero, coeff, &cols FCONE FCONE);
    symmetrize_packed_square(coeff, cols);
    F77_CALL(dpotrf)(&uplo, &cols, coeff, &cols, &info FCONE);
    bool chol_ok = (info == 0);
    for (int col = 0; chol_ok && col < cols; ++col) {
      if (coeff[col + static_cast<int64_t>(col) * cols] <= 100.0 * DBL_EPSILON) {
        chol_ok = false;
      }
    }
    if (chol_ok) {
      F77_CALL(dtrsm)(&right, &uplo, &trans_N, &diag, &n, &cols, &one,
                      coeff, &cols, Z_block, &n FCONE FCONE FCONE FCONE);
      if (n < 64) {
        F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                        &one, Z_block, &n, Z_block, &n,
                        &zero, coeff, &cols FCONE FCONE);
        symmetrize_packed_square(coeff, cols);
        F77_CALL(dpotrf)(&uplo, &cols, coeff, &cols, &info FCONE);
        chol_ok = (info == 0);
        for (int col = 0; chol_ok && col < cols; ++col) {
          if (coeff[col + static_cast<int64_t>(col) * cols] <= 100.0 * DBL_EPSILON) {
            chol_ok = false;
          }
        }
      }
      if (chol_ok) {
        if (n < 64) {
          F77_CALL(dtrsm)(&right, &uplo, &trans_N, &diag, &n, &cols, &one,
                          coeff, &cols, Z_block, &n FCONE FCONE FCONE FCONE);
        }
        for (int col = 0; col < cols && *m_active < m_max; ++col) {
          std::memcpy(V_active + static_cast<int64_t>(*m_active) * n,
                      Z_block + static_cast<int64_t>(col) * n,
                      sizeof(double) * static_cast<size_t>(n));
          ++(*m_active);
        }
        return cols;
      }
    }
  }

  int accepted = 0;
  const int batch_start = *m_active;
  for (int col = 0; col < cols && *m_active < m_max; ++col) {
    double* z_col = Z_block + static_cast<int64_t>(col) * n;
    if (accepted > 0) {
      trl_orthogonalise(nullptr, 0,
                        V_active + static_cast<int64_t>(batch_start) * n,
                        accepted, z_col, tmp, n);
    }
    const double nz = trl_norm2(z_col, n);
    if (nz <= 100.0 * DBL_EPSILON) {
      continue;
    }
    const double inv_nz = 1.0 / nz;
    double* dst = V_active + static_cast<int64_t>(*m_active) * n;
    for (int row = 0; row < n; ++row) {
      dst[row] = z_col[row] * inv_nz;
    }
    ++(*m_active);
    ++accepted;
  }
  return accepted;
}

static int apply_active_block(void* impl, EigencoreApplyFn apply,
                              int n, int first_col, int cols,
                              double* V_active, double* AV_active,
                              EigencoreWorkspace* workspace,
                              int* matvecs_out) {
  if (cols <= 0) {
    return 0;
  }
  const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, cols,
                       V_active + static_cast<int64_t>(first_col) * n, n,
                       1.0, 0.0,
                       AV_active + static_cast<int64_t>(first_col) * n, n,
                       workspace);
  if (rc == 0 && matvecs_out != nullptr) {
    ++(*matvecs_out);
  }
  return rc;
}

struct BlockGolubKahanBasisScratch {
  double* Z_v = nullptr;
  double* Z_u = nullptr;
  double* coeff = nullptr;
  double* tmp = nullptr;
  size_t bytes = 0;
  bool transient = false;
};

static double* eigencore_r_alloc_zero_doubles(size_t n) {
  double* out = reinterpret_cast<double*>(R_alloc(n > 0 ? n : 1, sizeof(double)));
  std::memset(out, 0, sizeof(double) * n);
  return out;
}

static void block_golub_kahan_basis_scratch_free(BlockGolubKahanBasisScratch* scratch) {
  if (!scratch->transient) {
    std::free(scratch->Z_v);
    std::free(scratch->Z_u);
    std::free(scratch->coeff);
    std::free(scratch->tmp);
  }
  scratch->Z_v = nullptr;
  scratch->Z_u = nullptr;
  scratch->coeff = nullptr;
  scratch->tmp = nullptr;
  scratch->bytes = 0;
  scratch->transient = false;
}

static int block_golub_kahan_basis_scratch_alloc(BlockGolubKahanBasisScratch* scratch,
                                                 int m,
                                                 int n,
                                                 int max_subspace,
                                                 int block_size) {
  const size_t nb = static_cast<size_t>(n) * static_cast<size_t>(block_size);
  const size_t mb = static_cast<size_t>(m) * static_cast<size_t>(block_size);
  const int coeff_rows = (max_subspace > block_size) ? max_subspace : block_size;
  const size_t coeff_elems = static_cast<size_t>(coeff_rows) *
    static_cast<size_t>(block_size);
  const size_t tmp_len = static_cast<size_t>((m > n) ? m : n);
  scratch->Z_v = eigencore_r_alloc_zero_doubles(nb);
  scratch->Z_u = eigencore_r_alloc_zero_doubles(mb);
  scratch->coeff = eigencore_r_alloc_zero_doubles(coeff_elems);
  scratch->tmp = eigencore_r_alloc_zero_doubles(tmp_len);
  scratch->transient = true;
  scratch->bytes = (nb + mb + coeff_elems + tmp_len) * sizeof(double);
  return 0;
}

static int native_block_golub_kahan_basis_run_with_scratch(
                                              void* impl,
                                              EigencoreApplyFn apply,
                                              int m,
                                              int n,
                                              int max_subspace,
                                              int block_size,
                                              const double* start_block,
                                              const double* start_av_block,
                                              int start_av_cols,
                                              const double* V_locked,
                                              int n_locked_v,
                                              const double* U_locked,
                                              int n_locked_u,
                                              double* V,
                                              double* AV,
                                              double* U,
                                              BlockGolubKahanBasisScratch* scratch,
                                              int* active_v_out,
                                              int* active_u_out,
                                              int* iterations_out,
                                              int* matvecs_out,
                                              int* ortho_passes_out,
                                              int* cached_start_used_out) {
  *active_v_out = 0;
  *active_u_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *ortho_passes_out = 0;
  if (cached_start_used_out != nullptr) {
    *cached_start_used_out = 0;
  }

  const size_t nv = static_cast<size_t>(n) * max_subspace;
  const size_t mv = static_cast<size_t>(m) * max_subspace;
  std::memset(V, 0, sizeof(double) * nv);
  std::memset(AV, 0, sizeof(double) * mv);
  std::memset(U, 0, sizeof(double) * mv);

  if (scratch == nullptr || scratch->Z_v == nullptr || scratch->Z_u == nullptr ||
      scratch->coeff == nullptr || scratch->tmp == nullptr) {
    return -2;
  }
  double* Z_v = scratch->Z_v;
  double* Z_u = scratch->Z_u;
  double* coeff = scratch->coeff;
  double* tmp = scratch->tmp;

  int active_v = 0;
  int active_u = 0;
  int last_v_start = 0;
  int last_v_cols = 0;
  if (start_av_block != nullptr && start_av_cols > 0) {
    int cached_cols = start_av_cols;
    if (cached_cols > block_size) {
      cached_cols = block_size;
    }
    if (cached_cols > max_subspace) {
      cached_cols = max_subspace;
    }
    for (int col = 0; col < cached_cols; ++col) {
      std::memcpy(V + static_cast<int64_t>(col) * n,
                  start_block + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
      std::memcpy(AV + static_cast<int64_t>(col) * m,
                  start_av_block + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
    }
    active_v = cached_cols;
    if (cached_start_used_out != nullptr) {
      *cached_start_used_out = 1;
    }
    if (cached_cols < block_size && active_v < max_subspace) {
      const int suffix_start = active_v;
      const int suffix_cols = block_accept_columns_blas3(
        start_block + static_cast<int64_t>(cached_cols) * n, n,
        block_size - cached_cols, V_locked, n_locked_v,
        V, &active_v, max_subspace, Z_v, block_size,
        coeff, tmp, n, block_size - cached_cols, ortho_passes_out
      );
      if (suffix_cols > 0) {
        EigencoreWorkspace suffix_workspace = {0, 0, nullptr, 0};
        const int rc_suffix = apply(
          impl, EIGENCORE_TRANSPOSE_NONE, suffix_cols,
          V + static_cast<int64_t>(suffix_start) * n, n,
          1.0, 0.0,
          AV + static_cast<int64_t>(suffix_start) * m, m,
          &suffix_workspace
        );
        if (rc_suffix != 0) {
          return rc_suffix;
        }
        ++(*matvecs_out);
      }
    }
    last_v_cols = active_v;
  } else {
    last_v_cols = block_accept_columns_blas3(
      start_block, n, block_size, V_locked, n_locked_v,
      V, &active_v, max_subspace, Z_v, block_size,
      coeff, tmp, n, block_size, ortho_passes_out
    );
  }
  if (last_v_cols == 0) {
    for (int col = 0; col < block_size && active_v < max_subspace && col < n; ++col) {
      std::memset(Z_v, 0, sizeof(double) * static_cast<size_t>(n));
      Z_v[col] = 1.0;
      last_v_cols += block_accept_work_vector(
        V_locked, n_locked_v, V, &active_v, max_subspace,
        Z_v, tmp, n, ortho_passes_out
      );
    }
  }
  if (last_v_cols == 0) {
    return -4;
  }

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int rc = 0;
  if (cached_start_used_out == nullptr || *cached_start_used_out == 0) {
    rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, last_v_cols,
               V + static_cast<int64_t>(last_v_start) * n, n,
               1.0, 0.0,
               AV + static_cast<int64_t>(last_v_start) * m, m,
               &workspace);
    if (rc != 0) {
      return rc;
    }
    ++(*matvecs_out);
  }

  while (active_v < max_subspace && last_v_cols > 0) {
    const int accepted_u_start = active_u;
    const int accepted_u = block_accept_columns_blas3(
      AV + static_cast<int64_t>(last_v_start) * m, m, last_v_cols,
      U_locked, n_locked_u, U, &active_u, max_subspace, Z_u, block_size,
      coeff, tmp, m, block_size, ortho_passes_out
    );
    if (accepted_u == 0) {
      break;
    }

    rc = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, accepted_u,
               U + static_cast<int64_t>(accepted_u_start) * m, m,
               1.0, 0.0, Z_v, n, &workspace);
    if (rc != 0) {
      return rc;
    }
    ++(*matvecs_out);

    const int accepted_v_start = active_v;
    const int accepted_v = block_accept_columns_blas3(
      Z_v, n, accepted_u, V_locked, n_locked_v,
      V, &active_v, max_subspace, Z_v, block_size,
      coeff, tmp, n, max_subspace - active_v, ortho_passes_out
    );
    if (accepted_v == 0) {
      break;
    }

    rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, accepted_v,
               V + static_cast<int64_t>(accepted_v_start) * n, n,
               1.0, 0.0,
               AV + static_cast<int64_t>(accepted_v_start) * m, m,
               &workspace);
    if (rc != 0) {
      return rc;
    }
    ++(*matvecs_out);
    ++(*iterations_out);
    last_v_start = accepted_v_start;
    last_v_cols = accepted_v;
  }

  *active_v_out = active_v;
  *active_u_out = active_u;
  return 0;
}

static int native_block_golub_kahan_basis_run(void* impl,
                                              EigencoreApplyFn apply,
                                              int m,
                                              int n,
                                              int max_subspace,
                                              int block_size,
                                              const double* start_block,
                                              const double* start_av_block,
                                              int start_av_cols,
                                              const double* V_locked,
                                              int n_locked_v,
                                              const double* U_locked,
                                              int n_locked_u,
                                              double* V,
                                              double* AV,
                                              double* U,
                                              int* active_v_out,
                                              int* active_u_out,
                                              int* iterations_out,
                                              int* matvecs_out,
                                              int* ortho_passes_out,
                                              int* cached_start_used_out) {
  BlockGolubKahanBasisScratch scratch;
  const int scratch_status = block_golub_kahan_basis_scratch_alloc(
    &scratch, m, n, max_subspace, block_size
  );
  if (scratch_status != 0) {
    return scratch_status;
  }
  const int status = native_block_golub_kahan_basis_run_with_scratch(
    impl, apply, m, n, max_subspace, block_size, start_block,
    start_av_block, start_av_cols, V_locked, n_locked_v, U_locked, n_locked_u,
    V, AV, U, &scratch,
    active_v_out, active_u_out, iterations_out, matvecs_out,
    ortho_passes_out, cached_start_used_out
  );
  block_golub_kahan_basis_scratch_free(&scratch);
  return status;
}

static SEXP block_golub_kahan_basis_pack(int n,
                                         int m,
                                         int max_subspace,
                                         const double* V,
                                         const double* AV,
                                         int active_v,
                                         int active_u,
                                         int iterations,
                                         int matvecs,
                                         int ortho_passes,
                                         int cached_start_used) {
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, max_subspace));
  SEXP AV_ = PROTECT(allocMatrix(REALSXP, m, max_subspace));
  std::memcpy(REAL(V_), V, sizeof(double) * static_cast<size_t>(n) * max_subspace);
  std::memcpy(REAL(AV_), AV, sizeof(double) * static_cast<size_t>(m) * max_subspace);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, V_);
  SET_VECTOR_ELT(out_, 1, AV_);
  SET_VECTOR_ELT(out_, 2, ScalarInteger(active_v));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(active_u));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(ortho_passes));
  SET_VECTOR_ELT(out_, 7, ScalarLogical(cached_start_used != 0));
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("V"));
  SET_STRING_ELT(names_, 1, mkChar("AV"));
  SET_STRING_ELT(names_, 2, mkChar("active_cols"));
  SET_STRING_ELT(names_, 3, mkChar("active_left_cols"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 7, mkChar("cached_start_used"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_basis(SEXP A_,
                                                        SEXP max_subspace_,
                                                        SEXP start_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native dense block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

extern "C" SEXP eigencore_block_golub_kahan_dense_basis_cached(SEXP A_,
                                                               SEXP max_subspace_,
                                                               SEXP start_,
                                                               SEXP start_av_) {
  if (!isReal(A_) || !isReal(start_) || !isReal(start_av_)) {
    error("A, start, and start_av must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || dimAV == R_NilValue) {
    error("A, start, and start_av must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  const int start_av_cols = INTEGER(dimAV)[1];
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native dense cached block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

extern "C" SEXP eigencore_block_golub_kahan_csc_basis(SEXP i_, SEXP p_,
                                                      SEXP x_, SEXP dim_,
                                                      SEXP max_subspace_,
                                                      SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_)) {
    error("invalid CSC block Golub-Kahan basis inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(dim_) != 2) {
    error("start must be a matrix and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native CSC block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

extern "C" SEXP eigencore_block_golub_kahan_csc_basis_cached(SEXP i_, SEXP p_,
                                                             SEXP x_, SEXP dim_,
                                                             SEXP max_subspace_,
                                                             SEXP start_,
                                                             SEXP start_av_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_) || !isReal(start_av_)) {
    error("invalid cached CSC block Golub-Kahan basis inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimS == R_NilValue || dimAV == R_NilValue || LENGTH(dim_) != 2) {
    error("start and start_av must be matrices and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  const int start_av_cols = INTEGER(dimAV)[1];
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  std::vector<double> V(static_cast<size_t>(n) * max_subspace, 0.0);
  std::vector<double> AV(static_cast<size_t>(m) * max_subspace, 0.0);
  std::vector<double> U(static_cast<size_t>(m) * max_subspace, 0.0);
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    V.data(), AV.data(), U.data(),
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  if (status != 0) {
    error("native CSC cached block Golub-Kahan basis failed with status=%d", status);
  }
  return block_golub_kahan_basis_pack(
    n, m, max_subspace, V.data(), AV.data(),
    active_v, active_u, iterations, matvecs, ortho_passes, cached_start_used
  );
}

static void subtract_projected_range(const double* V_active,
                                     const double* T_proj,
                                     int ldt,
                                     int n,
                                     int range_start,
                                     int range_cols,
                                     int current_start,
                                     int current_cols,
                                     double* W,
                                     double* coeff) {
  if (range_cols <= 0 || current_cols <= 0) {
    return;
  }
  if (range_cols <= 4 && current_cols <= 4) {
    for (int col = 0; col < current_cols; ++col) {
      double* w_col = W + static_cast<int64_t>(col) * n;
      for (int basis = 0; basis < range_cols; ++basis) {
        const double coeff_value =
          T_proj[(range_start + basis) +
                 static_cast<int64_t>(current_start + col) * ldt];
        const double* v_col =
          V_active + static_cast<int64_t>(range_start + basis) * n;
        for (int row = 0; row < n; ++row) {
          w_col[row] -= v_col[row] * coeff_value;
        }
      }
    }
    return;
  }
  for (int col = 0; col < current_cols; ++col) {
    for (int row = 0; row < range_cols; ++row) {
      coeff[row + static_cast<int64_t>(col) * range_cols] =
        T_proj[(range_start + row) +
               static_cast<int64_t>(current_start + col) * ldt];
    }
  }
  const char trans_N = 'N';
  const double one = 1.0;
  const double minus_one = -1.0;
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &current_cols, &range_cols,
                  &minus_one,
                  V_active + static_cast<int64_t>(range_start) * n, &n,
                  coeff, &range_cols,
                  &one, W, &n FCONE FCONE);
}

static void form_structured_projected_block_residual(const double* V_active,
                                                     const double* AV_active,
                                                     const double* T_proj,
                                                     int ldt,
                                                     int n,
                                                     int current_start,
                                                     int current_cols,
                                                     int previous_start,
                                                     int previous_cols,
                                                     double* W,
                                                     double* coeff) {
  if (current_cols <= 0) {
    return;
  }
  for (int col = 0; col < current_cols; ++col) {
    std::memcpy(W + static_cast<int64_t>(col) * n,
                AV_active + static_cast<int64_t>(current_start + col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }
  subtract_projected_range(V_active, T_proj, ldt, n,
                           previous_start, previous_cols,
                           current_start, current_cols, W, coeff);
  subtract_projected_range(V_active, T_proj, ldt, n,
                           current_start, current_cols,
                           current_start, current_cols, W, coeff);
}

static void projection_update_self_block(double* T_proj, int ldt,
                                         const double* V_active,
                                         const double* AV_active,
                                         int n, int start, int cols,
                                         double* scratch) {
  if (cols <= 0) {
    return;
  }
  if (cols <= 4) {
    for (int col = 0; col < cols; ++col) {
      const double* av_col =
        AV_active + static_cast<int64_t>(start + col) * n;
      for (int row = 0; row < cols; ++row) {
        const double* v_col =
          V_active + static_cast<int64_t>(start + row) * n;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
          sum += v_col[i] * av_col[i];
        }
        scratch[row + static_cast<int64_t>(col) * cols] = sum;
      }
    }
    symmetrize_packed_square(scratch, cols);
    for (int col = 0; col < cols; ++col) {
      for (int row = 0; row < cols; ++row) {
        T_proj[(start + row) + static_cast<int64_t>(start + col) * ldt] =
          scratch[row + static_cast<int64_t>(col) * cols];
      }
    }
    return;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_T, &trans_N, &cols, &cols, &n,
                  &one,
                  V_active + static_cast<int64_t>(start) * n, &n,
                  AV_active + static_cast<int64_t>(start) * n, &n,
                  &zero, scratch, &cols FCONE FCONE);
  symmetrize_packed_square(scratch, cols);
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < cols; ++row) {
      T_proj[(start + row) + static_cast<int64_t>(start + col) * ldt] =
        scratch[row + static_cast<int64_t>(col) * cols];
    }
  }
}

static void projection_update_cross_block(double* T_proj, int ldt,
                                          const double* V_active,
                                          const double* AV_active,
                                          int n,
                                          int left_start, int left_cols,
                                          int right_start, int right_cols,
                                          double* scratch) {
  if (left_cols <= 0 || right_cols <= 0) {
    return;
  }
  if (left_cols <= 4 && right_cols <= 4) {
    for (int col = 0; col < right_cols; ++col) {
      const double* av_col =
        AV_active + static_cast<int64_t>(right_start + col) * n;
      for (int row = 0; row < left_cols; ++row) {
        const double* v_col =
          V_active + static_cast<int64_t>(left_start + row) * n;
        double value = 0.0;
        for (int i = 0; i < n; ++i) {
          value += v_col[i] * av_col[i];
        }
        T_proj[(left_start + row) + static_cast<int64_t>(right_start + col) * ldt] = value;
        T_proj[(right_start + col) + static_cast<int64_t>(left_start + row) * ldt] = value;
      }
    }
    return;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_T, &trans_N, &left_cols, &right_cols, &n,
                  &one,
                  V_active + static_cast<int64_t>(left_start) * n, &n,
                  AV_active + static_cast<int64_t>(right_start) * n, &n,
                  &zero, scratch, &left_cols FCONE FCONE);
  for (int col = 0; col < right_cols; ++col) {
    for (int row = 0; row < left_cols; ++row) {
      const double value = scratch[row + static_cast<int64_t>(col) * left_cols];
      T_proj[(left_start + row) + static_cast<int64_t>(right_start + col) * ldt] = value;
      T_proj[(right_start + col) + static_cast<int64_t>(left_start + row) * ldt] = value;
    }
  }
}

static void projection_update_small_self_and_cross(double* T_proj, int ldt,
                                                   const double* V_active,
                                                   const double* AV_active,
                                                   int n,
                                                   int left_start,
                                                   int left_cols,
                                                   int right_start,
                                                   int right_cols,
                                                   double* scratch) {
  if (left_cols <= 0 || right_cols <= 0) {
    projection_update_self_block(T_proj, ldt, V_active, AV_active,
                                 n, right_start, right_cols, scratch);
    return;
  }
  if (left_cols > 4 || right_cols > 4) {
    projection_update_self_block(T_proj, ldt, V_active, AV_active,
                                 n, right_start, right_cols, scratch);
    projection_update_cross_block(T_proj, ldt, V_active, AV_active,
                                  n, left_start, left_cols,
                                  right_start, right_cols, scratch);
    return;
  }

  double cross[16];
  const int cross_len = left_cols * right_cols;
  const int self_len = right_cols * right_cols;
  for (int pos = 0; pos < cross_len; ++pos) {
    cross[pos] = 0.0;
  }
  for (int pos = 0; pos < self_len; ++pos) {
    scratch[pos] = 0.0;
  }

  const double* right_vecs[4];
  const double* left_vecs[4];
  for (int row = 0; row < right_cols; ++row) {
    right_vecs[row] = V_active + static_cast<int64_t>(right_start + row) * n;
  }
  for (int row = 0; row < left_cols; ++row) {
    left_vecs[row] = V_active + static_cast<int64_t>(left_start + row) * n;
  }

  for (int col = 0; col < right_cols; ++col) {
    const double* av_col =
      AV_active + static_cast<int64_t>(right_start + col) * n;
    double self_acc[4] = {0.0, 0.0, 0.0, 0.0};
    double cross_acc[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i) {
      const double av = av_col[i];
      for (int row = 0; row < right_cols; ++row) {
        self_acc[row] += right_vecs[row][i] * av;
      }
      for (int row = 0; row < left_cols; ++row) {
        cross_acc[row] += left_vecs[row][i] * av;
      }
    }
    for (int row = 0; row < right_cols; ++row) {
      scratch[row + static_cast<int64_t>(col) * right_cols] = self_acc[row];
    }
    for (int row = 0; row < left_cols; ++row) {
      cross[row + static_cast<int64_t>(col) * left_cols] = cross_acc[row];
    }
  }

  symmetrize_packed_square(scratch, right_cols);
  for (int col = 0; col < right_cols; ++col) {
    for (int row = 0; row < right_cols; ++row) {
      T_proj[(right_start + row) +
             static_cast<int64_t>(right_start + col) * ldt] =
        scratch[row + static_cast<int64_t>(col) * right_cols];
    }
    for (int row = 0; row < left_cols; ++row) {
      const double value = cross[row + static_cast<int64_t>(col) * left_cols];
      T_proj[(left_start + row) +
             static_cast<int64_t>(right_start + col) * ldt] = value;
      T_proj[(right_start + col) +
             static_cast<int64_t>(left_start + row) * ldt] = value;
    }
  }
}

static void projection_copy_upper_compact(const double* T_proj, int ldt,
                                          double* compact, int m_active) {
  for (int col = 0; col < m_active; ++col) {
    for (int row = 0; row <= col; ++row) {
      compact[row + static_cast<int64_t>(col) * m_active] =
        T_proj[row + static_cast<int64_t>(col) * ldt];
    }
  }
}

struct ThickRestartBuffers {
  double* V_active;
  double* AV_active;
  double* T_proj;      // m_max x m_max structured projected problem
  // S_eig is a k x k scratch matrix reused across three distinct roles inside
  // a single solve cycle:
  //   role 1: V^T V Gram for the orthogonality probe in final_polish_block_ritz
  //   role 2: Cholesky factor for re-orthonormalization (dpotrf / dtrsm in-place)
  //   role 3: V^T A V projected eigenproblem (dsyev_inplace) for full polish
  // Each role overwrites the previous, in strict sequence within one call.
  // Maintaining three separate buffers would cost an extra 2 * k_max^2 doubles
  // per cycle for negligible runtime savings; the role transitions are flagged
  // with explicit comments where they happen.
  double* S_eig;       // m_max x m_max — see role notes above
  double* S_selected;  // selected Ritz vectors, m_max x selected_capacity
  double* theta;       // m_max
  double* B_v;         // n x selected_capacity
  double* B_av;        // n x selected_capacity
  double* Z_block;     // n x block_size
  double* coeff_block; // m_max x block_size
  double* z;           // n
  double* tmp;         // max(k_target, m_max)
  double* ritz_res;    // m_max
  int*    selected;    // m_max
  int*    is_locked;   // m_max
  double* dsyev_work;
  int     dsyev_lwork;
  int*    dsyevd_iwork;
  int     dsyevd_liwork;
  int     selected_capacity;
};

static void trl_buffers_free(ThickRestartBuffers* b) {
  std::free(b->V_active);
  std::free(b->AV_active);
  std::free(b->T_proj);
  std::free(b->S_eig);
  std::free(b->S_selected);
  std::free(b->theta);
  std::free(b->B_v);
  std::free(b->B_av);
  std::free(b->Z_block);
  std::free(b->coeff_block);
  std::free(b->z);
  std::free(b->tmp);
  std::free(b->ritz_res);
  std::free(b->selected);
  std::free(b->is_locked);
  std::free(b->dsyev_work);
  std::free(b->dsyevd_iwork);
}

static int trl_buffers_alloc(ThickRestartBuffers* b, int n, int k_target,
                             int m_max, int block_cols) {
  std::memset(b, 0, sizeof(*b));
  const size_t nm = static_cast<size_t>(n) * static_cast<size_t>(m_max);
  const size_t mm = static_cast<size_t>(m_max) * static_cast<size_t>(m_max);
  const size_t nb = static_cast<size_t>(n) * static_cast<size_t>(block_cols);
  const size_t mb = static_cast<size_t>(m_max) * static_cast<size_t>(block_cols);
  int selected_capacity = 2 * k_target;
  if (selected_capacity < k_target + 5) selected_capacity = k_target + 5;
  if (selected_capacity < block_cols) selected_capacity = block_cols;
  if (selected_capacity < 1) selected_capacity = 1;
  if (selected_capacity > m_max) selected_capacity = m_max;
  b->selected_capacity = selected_capacity;
  const size_t ms = static_cast<size_t>(m_max) * static_cast<size_t>(selected_capacity);
  const size_t ns = static_cast<size_t>(n) * static_cast<size_t>(selected_capacity);
  b->V_active  = static_cast<double*>(std::malloc(nm * sizeof(double)));
  b->AV_active = static_cast<double*>(std::malloc(nm * sizeof(double)));
  b->T_proj    = static_cast<double*>(std::calloc(mm, sizeof(double)));
  b->S_eig     = static_cast<double*>(std::malloc(mm * sizeof(double)));
  b->S_selected = static_cast<double*>(std::malloc(ms * sizeof(double)));
  b->theta     = static_cast<double*>(std::malloc(static_cast<size_t>(m_max) * sizeof(double)));
  b->B_v       = static_cast<double*>(std::malloc(ns * sizeof(double)));
  b->B_av      = static_cast<double*>(std::malloc(ns * sizeof(double)));
  b->Z_block   = static_cast<double*>(std::malloc(nb * sizeof(double)));
  b->coeff_block = static_cast<double*>(std::malloc(mb * sizeof(double)));
  b->z         = static_cast<double*>(std::malloc(static_cast<size_t>(n) * sizeof(double)));
  const int tmp_len = (k_target > m_max) ? k_target : m_max;
  b->tmp       = static_cast<double*>(std::malloc(static_cast<size_t>(tmp_len > 0 ? tmp_len : 1) * sizeof(double)));
  b->ritz_res  = static_cast<double*>(std::malloc(static_cast<size_t>(m_max) * sizeof(double)));
  b->selected  = static_cast<int*>(std::malloc(static_cast<size_t>(m_max) * sizeof(int)));
  b->is_locked = static_cast<int*>(std::malloc(static_cast<size_t>(m_max) * sizeof(int)));
  b->dsyev_lwork = trl_dsyevd_query(m_max, &b->dsyevd_liwork);
  if (b->dsyev_lwork < 26 * m_max) b->dsyev_lwork = 26 * m_max;
  if (b->dsyevd_liwork < 10 * m_max) b->dsyevd_liwork = 10 * m_max;
  if (b->dsyev_lwork < 1) b->dsyev_lwork = 1;
  if (b->dsyevd_liwork < 1) b->dsyevd_liwork = 1;
  b->dsyev_work = static_cast<double*>(std::malloc(static_cast<size_t>(b->dsyev_lwork) * sizeof(double)));
  b->dsyevd_iwork = static_cast<int*>(std::malloc(static_cast<size_t>(b->dsyevd_liwork) * sizeof(int)));
  if (b->V_active == nullptr || b->AV_active == nullptr ||
      b->T_proj == nullptr || b->S_eig == nullptr ||
      b->S_selected == nullptr ||
      b->theta == nullptr || b->B_v == nullptr || b->B_av == nullptr ||
      b->Z_block == nullptr || b->coeff_block == nullptr ||
      b->z == nullptr || b->tmp == nullptr || b->ritz_res == nullptr ||
      b->selected == nullptr || b->is_locked == nullptr ||
      b->dsyev_work == nullptr || b->dsyevd_iwork == nullptr) {
    trl_buffers_free(b);
    return -1;
  }
  return 0;
}

static int final_polish_block_ritz(void* impl,
                                   EigencoreApplyFn apply,
                                   int n,
                                   int k_target,
                                   int target_kind,
                                   double tol,
                                   double norm_a,
                                   double* V_out,
                                   double* lambda_out,
                                   double* residuals_out,
                                   int* converged_out,
                                   int* n_converged_out,
                                   ThickRestartBuffers* buf,
                                   EigencoreWorkspace* workspace,
                                   int* matvecs_out) {
  if (k_target <= 0) {
    if (n_converged_out != nullptr) {
      *n_converged_out = 0;
    }
    return 0;
  }

  const char trans_T = 'T';
  const char trans_N = 'N';
  const char side_R = 'R';
  const char uplo_U = 'U';
  const char diag_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  // S_eig role 1: V_out^T V_out Gram for orthogonality probe.
  F77_CALL(dgemm)(&trans_T, &trans_N, &k_target, &k_target, &n,
                  &one, V_out, &n, V_out, &n,
                  &zero, buf->S_eig, &k_target FCONE FCONE);
  double max_orthogonality = 0.0;
  for (int col = 0; col < k_target; ++col) {
    for (int row = 0; row < k_target; ++row) {
      const double expected = (row == col) ? 1.0 : 0.0;
      const double loss = fabs(buf->S_eig[row + static_cast<int64_t>(col) * k_target] - expected);
      if (loss > max_orthogonality) {
        max_orthogonality = loss;
      }
    }
  }
  const double orthogonality_tolerance =
    (tol > sqrt(DBL_EPSILON)) ? tol : sqrt(DBL_EPSILON);
  int prepolish_converged = 0;
  if (max_orthogonality <= orthogonality_tolerance) {
    for (int col = 0; col < k_target; ++col) {
      const double* vec = V_out + static_cast<int64_t>(col) * n;
      const double scale =
        standard_eigen_lock_scale(norm_a, lambda_out[col], vec, n);
      converged_out[col] = (residuals_out[col] <= tol * scale) ? 1 : 0;
      if (converged_out[col]) {
        ++prepolish_converged;
      }
    }
    if (n_converged_out != nullptr) {
      *n_converged_out = prepolish_converged;
    }
    if (prepolish_converged == k_target) {
      return 0;
    }
  }

  // S_eig role 2: in-place Cholesky factor of the Gram from role 1, used to
  // re-orthonormalize V_out via right-side dtrsm. Overwrites role-1 contents.
  symmetrize_packed_square(buf->S_eig, k_target);
  F77_CALL(dpotrf)(&uplo_U, &k_target, buf->S_eig, &k_target, &info FCONE);
  if (info != 0) {
    return 0;
  }
  F77_CALL(dtrsm)(&side_R, &uplo_U, &trans_N, &diag_N, &n, &k_target, &one,
                  buf->S_eig, &k_target, V_out, &n FCONE FCONE FCONE FCONE);

  int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, k_target,
                 V_out, n, 1.0, 0.0, buf->B_av, n, workspace);
  if (rc != 0) {
    return rc;
  }
  if (matvecs_out != nullptr) {
    ++(*matvecs_out);
  }

  int n_converged_simple = 0;
  for (int col = 0; col < k_target; ++col) {
    const double* vec = V_out + static_cast<int64_t>(col) * n;
    const double* av = buf->B_av + static_cast<int64_t>(col) * n;
    long double theta_sum = 0.0L;
    for (int row = 0; row < n; ++row) {
      theta_sum += static_cast<long double>(vec[row]) * av[row];
    }
    const double theta = static_cast<double>(theta_sum);
    long double ss = 0.0L;
    for (int row = 0; row < n; ++row) {
      const double diff = av[row] - theta * vec[row];
      ss += static_cast<long double>(diff) * diff;
    }
    const double res = sqrt(static_cast<double>(ss));
    lambda_out[col] = theta;
    residuals_out[col] = res;
    const double scale = standard_eigen_lock_scale(norm_a, theta, vec, n);
    converged_out[col] = (res <= tol * scale) ? 1 : 0;
    if (converged_out[col]) {
      ++n_converged_simple;
    }
  }
  if (n_converged_out != nullptr) {
    *n_converged_out = n_converged_simple;
  }
  if (n_converged_simple == k_target) {
    return 0;
  }

  // S_eig role 3: V_out^T A V_out projected eigenproblem, solved in place by
  // symmetric_eigen_inplace. Overwrites the role-2 Cholesky factor; columns of
  // S_eig now hold the projected-problem eigenvectors used to rotate V_out.
  // Reuse A * V_out from the simple residual check; V_out has not changed.
  F77_CALL(dgemm)(&trans_T, &trans_N, &k_target, &k_target, &n,
                  &one, V_out, &n, buf->B_av, &n,
                  &zero, buf->S_eig, &k_target FCONE FCONE);
  symmetrize_packed_square(buf->S_eig, k_target);
  rc = symmetric_eigen_inplace(buf->S_eig, k_target, buf->theta,
                               buf->dsyev_work, buf->dsyev_lwork,
                               buf->dsyevd_iwork, buf->dsyevd_liwork);
  if (rc != 0) {
    return rc;
  }
  selected_sorted_ritz_indices(buf->theta, k_target, k_target,
                               target_kind, buf->selected);

  for (int col = 0; col < k_target; ++col) {
    const int idx = buf->selected[col];
    for (int row = 0; row < k_target; ++row) {
      buf->S_selected[row + static_cast<int64_t>(col) * k_target] =
        buf->S_eig[row + static_cast<int64_t>(idx) * k_target];
    }
  }

  std::memcpy(buf->AV_active, buf->B_av,
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(k_target));
  if (k_target <= 32) {
    combine_basis_columns_small(V_out, n, k_target,
                                buf->S_selected, k_target,
                                k_target, buf->B_v);
    combine_basis_columns_small(buf->AV_active, n, k_target,
                                buf->S_selected, k_target,
                                k_target, buf->B_av);
  } else {
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &k_target,
                    &one, V_out, &n, buf->S_selected, &k_target,
                    &zero, buf->B_v, &n FCONE FCONE);
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &k_target,
                    &one, buf->AV_active, &n, buf->S_selected, &k_target,
                    &zero, buf->B_av, &n FCONE FCONE);
  }

  int n_converged = 0;
  for (int col = 0; col < k_target; ++col) {
    const int idx = buf->selected[col];
    const double theta = buf->theta[idx];
    long double ss = 0.0L;
    double* residual = buf->B_av + static_cast<int64_t>(col) * n;
    double* vec = buf->B_v + static_cast<int64_t>(col) * n;
    for (int row = 0; row < n; ++row) {
      residual[row] -= theta * vec[row];
      ss += static_cast<long double>(residual[row]) * residual[row];
      V_out[row + static_cast<int64_t>(col) * n] = vec[row];
    }
    const double res = sqrt(static_cast<double>(ss));
    lambda_out[col] = theta;
    residuals_out[col] = res;
    const double scale = standard_eigen_lock_scale(norm_a, theta, vec, n);
    converged_out[col] = (res <= tol * scale) ? 1 : 0;
    if (converged_out[col]) {
      ++n_converged;
    }
  }
  if (n_converged_out != nullptr) {
    *n_converged_out = n_converged;
  }
  return 0;
}

static SEXP trl_pack_result(int n, int k_target, const double* V_locked,
                            const double* lambda, const double* residuals,
                            const int* converged, int n_locked,
                            int iterations, int matvecs, int restarts,
                            int m_active_final) {
  SEXP values_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k_target));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k_target));
  std::memcpy(REAL(values_), lambda, sizeof(double) * static_cast<size_t>(k_target));
  std::memcpy(REAL(vectors_), V_locked,
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(k_target));
  std::memcpy(REAL(residuals_), residuals, sizeof(double) * static_cast<size_t>(k_target));
  for (int i = 0; i < k_target; ++i) {
    LOGICAL(converged_)[i] = converged[i] ? TRUE : FALSE;
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, residuals_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(n_locked));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(restarts));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(m_active_final));
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("residuals"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("n_locked"));
  SET_STRING_ELT(names_, 5, mkChar("iterations"));
  SET_STRING_ELT(names_, 6, mkChar("matvecs"));
  SET_STRING_ELT(names_, 7, mkChar("restarts"));
  SET_STRING_ELT(names_, 8, mkChar("m_active_final"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(6);
  return out_;
}

static int native_block_lanczos_run(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int block_size,
    int target_kind,
    double tol,
    const double* start_block,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* nconv_out,
    int* iterations_out,
    int* matvecs_out,
    int* m_active_final_out
) {
  *nconv_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *m_active_final_out = 0;
  for (int i = 0; i < k_target; ++i) {
    lambda_out[i] = 0.0;
    residuals_out[i] = R_PosInf;
    converged_out[i] = 0;
    std::memset(V_out + static_cast<int64_t>(i) * n, 0,
                sizeof(double) * static_cast<size_t>(n));
  }

  const size_t nm = static_cast<size_t>(n) * static_cast<size_t>(m_max);
  const size_t nb = static_cast<size_t>(n) * static_cast<size_t>(block_size);
  const size_t mm = static_cast<size_t>(m_max) * static_cast<size_t>(m_max);
  double* V = static_cast<double*>(std::calloc(nm, sizeof(double)));
  double* AV = static_cast<double*>(std::calloc(nm, sizeof(double)));
  double* Z = static_cast<double*>(std::calloc(nb, sizeof(double)));
  double* AZ = static_cast<double*>(std::calloc(nb, sizeof(double)));
  double* H = static_cast<double*>(std::calloc(mm, sizeof(double)));
  double* S_selected = static_cast<double*>(std::calloc(mm, sizeof(double)));
  double* theta = static_cast<double*>(std::calloc(static_cast<size_t>(m_max), sizeof(double)));
  double* B_v = static_cast<double*>(std::calloc(static_cast<size_t>(n) * k_target, sizeof(double)));
  double* B_av = static_cast<double*>(std::calloc(static_cast<size_t>(n) * k_target, sizeof(double)));
  double* tmp = static_cast<double*>(std::calloc(static_cast<size_t>(m_max), sizeof(double)));
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>(m_max), sizeof(int)));
  const int dsyev_lwork_query = trl_dsyev_query(m_max);
  int dsyev_lwork = dsyev_lwork_query > 0 ? dsyev_lwork_query : 3 * m_max;
  double* dsyev_work = static_cast<double*>(std::calloc(static_cast<size_t>(dsyev_lwork), sizeof(double)));
  if (V == nullptr || AV == nullptr || Z == nullptr || AZ == nullptr ||
      H == nullptr || S_selected == nullptr || theta == nullptr ||
      B_v == nullptr || B_av == nullptr || tmp == nullptr ||
      selected == nullptr || dsyev_work == nullptr) {
    std::free(V); std::free(AV); std::free(Z); std::free(AZ);
    std::free(H); std::free(S_selected); std::free(theta);
    std::free(B_v); std::free(B_av); std::free(tmp);
    std::free(selected); std::free(dsyev_work);
    return -2;
  }

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int m_active = 0;
  int last_block_start = 0;
  int last_block_cols = 0;
  int source_start = 0;
  int source_cols = block_size;
  for (int col = 0; col < block_size; ++col) {
    std::memcpy(Z + static_cast<int64_t>(col) * n,
                start_block + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  while (m_active < m_max && source_cols > 0) {
    int accepted_start = m_active;
    int accepted = 0;
    for (int col = 0; col < source_cols && m_active < m_max; ++col) {
      double* z_col = Z + static_cast<int64_t>(col) * n;
      trl_orthogonalise(nullptr, 0, V, m_active, z_col, tmp, n);
      const double nz = trl_norm2(z_col, n);
      if (nz <= 100.0 * DBL_EPSILON) {
        continue;
      }
      const double inv_nz = 1.0 / nz;
      for (int row = 0; row < n; ++row) {
        V[static_cast<int64_t>(m_active) * n + row] = z_col[row] * inv_nz;
      }
      ++m_active;
      ++accepted;
    }

    if (accepted == 0) {
      break;
    }

    const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, accepted,
                         V + static_cast<int64_t>(accepted_start) * n, n,
                         1.0, 0.0,
                         AV + static_cast<int64_t>(accepted_start) * n, n,
                         &workspace);
    if (rc != 0) {
      std::free(V); std::free(AV); std::free(Z); std::free(AZ);
      std::free(H); std::free(S_selected); std::free(theta);
      std::free(B_v); std::free(B_av); std::free(tmp);
      std::free(selected); std::free(dsyev_work);
      return rc;
    }
    ++(*matvecs_out);
    ++(*iterations_out);
    last_block_start = accepted_start;
    last_block_cols = accepted;

    if (m_active >= m_max) {
      break;
    }

    source_start = last_block_start;
    source_cols = last_block_cols;
    for (int col = 0; col < source_cols; ++col) {
      std::memcpy(Z + static_cast<int64_t>(col) * n,
                  AV + static_cast<int64_t>(source_start + col) * n,
                  sizeof(double) * static_cast<size_t>(n));
    }
  }

  if (m_active < k_target) {
    std::free(V); std::free(AV); std::free(Z); std::free(AZ);
    std::free(H); std::free(S_selected); std::free(theta);
    std::free(B_v); std::free(B_av); std::free(tmp);
    std::free(selected); std::free(dsyev_work);
    return -4;
  }

  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_T, &trans_N, &m_active, &m_active, &n,
                  &one, V, &n, AV, &n,
                  &zero, H, &m_active FCONE FCONE);
  for (int i = 0; i < m_active; ++i) {
    for (int j = i + 1; j < m_active; ++j) {
      const double avg = 0.5 * (H[i + j * m_active] + H[j + i * m_active]);
      H[i + j * m_active] = avg;
      H[j + i * m_active] = avg;
    }
  }

  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  int lwork = dsyev_lwork;
  F77_CALL(dsyev)(&jobz, &uplo, &m_active, H, &m_active, theta,
                  dsyev_work, &lwork, &info FCONE FCONE);
  if (info != 0) {
    std::free(V); std::free(AV); std::free(Z); std::free(AZ);
    std::free(H); std::free(S_selected); std::free(theta);
    std::free(B_v); std::free(B_av); std::free(tmp);
    std::free(selected); std::free(dsyev_work);
    return -3;
  }

  selected_ritz_indices(theta, m_active, k_target, target_kind, selected);
  for (int p = 0; p < k_target; ++p) {
    const int idx = selected[p];
    lambda_out[p] = theta[idx];
    for (int row = 0; row < m_active; ++row) {
      S_selected[row + static_cast<int64_t>(p) * m_active] =
        H[row + static_cast<int64_t>(idx) * m_active];
    }
  }

  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &m_active,
                  &one, V, &n, S_selected, &m_active,
                  &zero, B_v, &n FCONE FCONE);
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k_target, &m_active,
                  &one, AV, &n, S_selected, &m_active,
                  &zero, B_av, &n FCONE FCONE);
  int nconv = 0;
  for (int col = 0; col < k_target; ++col) {
    long double s = 0.0L;
    for (int row = 0; row < n; ++row) {
      const double diff = B_av[row + static_cast<int64_t>(col) * n] -
                          lambda_out[col] * B_v[row + static_cast<int64_t>(col) * n];
      s += static_cast<long double>(diff) * diff;
    }
    residuals_out[col] = sqrt(static_cast<double>(s));
    const double scale_i = (fabs(lambda_out[col]) > 1.0) ? fabs(lambda_out[col]) : 1.0;
    converged_out[col] = residuals_out[col] <= tol * scale_i ? 1 : 0;
    if (converged_out[col]) ++nconv;
    std::memcpy(V_out + static_cast<int64_t>(col) * n,
                B_v + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
  }

  *nconv_out = nconv;
  *m_active_final_out = m_active;
  std::free(V); std::free(AV); std::free(Z); std::free(AZ);
  std::free(H); std::free(S_selected); std::free(theta);
  std::free(B_v); std::free(B_av); std::free(tmp);
  std::free(selected); std::free(dsyev_work);
  return 0;
}

struct BlockLanczosBestSnapshot {
  std::vector<double> V;
  std::vector<double> lambda;
  std::vector<double> residuals;
  std::vector<int> converged;
  std::vector<double> candidate_V;
  std::vector<int> candidate_converged;
  int filled = 0;
  int locked_prefix = 0;
  int nconv = -1;
  double max_backward_error = R_PosInf;

  BlockLanczosBestSnapshot(int n, int k_target) :
      V(static_cast<size_t>(n) * static_cast<size_t>(k_target), 0.0),
      lambda(static_cast<size_t>(k_target), 0.0),
      residuals(static_cast<size_t>(k_target), R_PosInf),
      converged(static_cast<size_t>(k_target), 0),
      candidate_V(static_cast<size_t>(n) * static_cast<size_t>(k_target), 0.0),
      candidate_converged(static_cast<size_t>(k_target), 0) {}
};

static int block_lanczos_expand_basis_to_budget(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int m_max,
    int block_size,
    const double* V_locked,
    int n_locked,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    int* m_active,
    int* previous_block_start,
    int* previous_block_cols,
    int* last_block_start,
    int* last_block_cols,
    int* iterations_out,
    int* matvecs_out,
    int* ortho_passes_out) {
  while (*m_active < m_max && *last_block_cols > 0) {
    auto timer = native_timer_now();
    form_structured_projected_block_residual(
      buf->V_active, buf->AV_active, buf->T_proj, m_max, n,
      *last_block_start, *last_block_cols,
      *previous_block_start, *previous_block_cols,
      buf->Z_block, buf->coeff_block
    );
    stages->recurrence += native_timer_elapsed(timer);

    const int accepted_start = *m_active;
    timer = native_timer_now();
    const int accepted = block_accept_columns_blas3(
      buf->Z_block, n, *last_block_cols, V_locked, n_locked,
      buf->V_active, m_active, m_max, buf->Z_block, block_size,
      buf->coeff_block, buf->tmp, n,
      block_size, ortho_passes_out,
      true
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    if (accepted == 0) {
      // Exact breakdown can occur in an invariant subspace that is not the
      // requested target. Continue deterministically instead of locking it.
      int continuation_accepted = 0;
      const int continuation_start = *m_active;
      for (int attempt = 0;
           continuation_accepted < block_size &&
             attempt < n + block_size &&
             *m_active < m_max;
           ++attempt) {
        std::memset(buf->z, 0, sizeof(double) * static_cast<size_t>(n));
        const int idx_basis = ((*iterations_out + 1) * 17 + attempt * 31) % n;
        buf->z[idx_basis < 0 ? -idx_basis : idx_basis] = 1.0;
        timer = native_timer_now();
        continuation_accepted += block_accept_work_vector(
          V_locked, n_locked, buf->V_active, m_active, m_max,
          buf->z, buf->tmp, n, ortho_passes_out
        );
        stages->reorthogonalization += native_timer_elapsed(timer);
      }
      if (continuation_accepted == 0) {
        break;
      }

      timer = native_timer_now();
      const int rc = apply_active_block(
        impl, apply, n, continuation_start, continuation_accepted,
        buf->V_active, buf->AV_active, workspace, matvecs_out
      );
      stages->apply += native_timer_elapsed(timer);
      if (rc != 0) {
        return rc;
      }

      timer = native_timer_now();
      projection_update_self_block(buf->T_proj, m_max, buf->V_active,
                                   buf->AV_active, n, 0, *m_active,
                                   buf->S_eig);
      {
        const double elapsed = native_timer_elapsed(timer);
        stages->projected_solve += elapsed;
        stages->projection_update += elapsed;
      }
      ++(*iterations_out);
      *previous_block_start = *last_block_start;
      *previous_block_cols = *last_block_cols;
      *last_block_start = continuation_start;
      *last_block_cols = continuation_accepted;
      continue;
    }

    timer = native_timer_now();
    const int rc = apply_active_block(impl, apply, n, accepted_start, accepted,
                                      buf->V_active, buf->AV_active, workspace,
                                      matvecs_out);
    stages->apply += native_timer_elapsed(timer);
    if (rc != 0) {
      return rc;
    }

    timer = native_timer_now();
    projection_update_self_block(buf->T_proj, m_max, buf->V_active,
                                 buf->AV_active, n, 0, *m_active,
                                 buf->S_eig);
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projection_update += elapsed;
    }
    ++(*iterations_out);
    *previous_block_start = *last_block_start;
    *previous_block_cols = *last_block_cols;
    *last_block_start = accepted_start;
    *last_block_cols = accepted;
  }
  return 0;
}

static void block_lanczos_maybe_capture_best_snapshot(
    int n,
    int k_target,
    int selected_count,
    int n_locked,
    double norm_a,
    double tol,
    const ThickRestartBuffers* buf,
    const double* V_out,
    const double* lambda_out,
    const double* residuals_out,
    BlockLanczosBestSnapshot* best) {
  int candidate_count = 0;
  int candidate_nconv = 0;
  double candidate_max_backward_error = 0.0;
  for (; candidate_count < n_locked && candidate_count < k_target; ++candidate_count) {
    const double* vec = V_out + static_cast<int64_t>(candidate_count) * n;
    std::memcpy(best->candidate_V.data() + static_cast<int64_t>(candidate_count) * n,
                vec,
                sizeof(double) * static_cast<size_t>(n));
    const double scale_i = standard_eigen_lock_scale(
      norm_a, lambda_out[candidate_count], vec, n
    );
    const double backward_error = residuals_out[candidate_count] / scale_i;
    if (backward_error > candidate_max_backward_error) {
      candidate_max_backward_error = backward_error;
    }
    best->candidate_converged[static_cast<size_t>(candidate_count)] =
      (residuals_out[candidate_count] <= tol * scale_i) ? 1 : 0;
    if (best->candidate_converged[static_cast<size_t>(candidate_count)]) {
      ++candidate_nconv;
    }
  }
  for (int p = 0; p < selected_count && candidate_count < k_target; ++p) {
    if (buf->is_locked[p]) {
      continue;
    }
    const int idx = buf->selected[p];
    const double* vec = buf->B_v + static_cast<int64_t>(p) * n;
    if (!vector_is_independent_from_locked(best->candidate_V.data(), candidate_count, vec, n)) {
      continue;
    }
    const double scale_i = standard_eigen_lock_scale(
      norm_a, buf->theta[idx], vec, n
    );
    const double backward_error = buf->ritz_res[p] / scale_i;
    if (backward_error > candidate_max_backward_error) {
      candidate_max_backward_error = backward_error;
    }
    best->candidate_converged[static_cast<size_t>(candidate_count)] =
      (buf->ritz_res[p] <= tol * scale_i) ? 1 : 0;
    if (best->candidate_converged[static_cast<size_t>(candidate_count)]) {
      ++candidate_nconv;
    }
    std::memcpy(best->candidate_V.data() + static_cast<int64_t>(candidate_count) * n,
                vec,
                sizeof(double) * static_cast<size_t>(n));
    ++candidate_count;
  }
  if (candidate_count != k_target ||
      (candidate_nconv < best->nconv ||
       (candidate_nconv == best->nconv &&
        candidate_max_backward_error >= best->max_backward_error))) {
    return;
  }

  int out_col = 0;
  std::memcpy(best->V.data(), best->candidate_V.data(),
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(k_target));
  for (; out_col < n_locked && out_col < k_target; ++out_col) {
    best->lambda[static_cast<size_t>(out_col)] = lambda_out[out_col];
    best->residuals[static_cast<size_t>(out_col)] = residuals_out[out_col];
  }
  for (int p = 0; p < selected_count && out_col < k_target; ++p) {
    if (buf->is_locked[p]) {
      continue;
    }
    const int idx = buf->selected[p];
    const double* vec = buf->B_v + static_cast<int64_t>(p) * n;
    if (!vector_is_independent_from_locked(best->V.data(), out_col, vec, n)) {
      continue;
    }
    best->lambda[static_cast<size_t>(out_col)] = buf->theta[idx];
    best->residuals[static_cast<size_t>(out_col)] = buf->ritz_res[p];
    ++out_col;
  }
  best->locked_prefix = n_locked;
  best->nconv = candidate_nconv;
  best->max_backward_error = candidate_max_backward_error;
  std::memcpy(best->converged.data(), best->candidate_converged.data(),
              sizeof(int) * static_cast<size_t>(k_target));
  best->filled = 1;
}

static int block_lanczos_restart_with_continuation_tail(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int block_size,
    int restart_idx,
    int selected_count,
    int n_locked,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    int* m_active,
    int* previous_block_start,
    int* previous_block_cols,
    int* last_block_start,
    int* last_block_cols,
    double* V_out,
    int* matvecs_out,
    int* restarts_out,
    int* ortho_passes_out) {
  auto timer = native_timer_now();
  const int remaining = k_target - n_locked;
  int keep_room = m_max - block_size;
  if (keep_room < 0) {
    keep_room = 0;
  }
  int pad = block_size > 4 ? block_size : 4;
  if (pad > k_target) {
    pad = k_target;
  }
  int k_keep = remaining + pad;
  if (k_keep < remaining) {
    k_keep = remaining;
  }
  if (k_keep > keep_room) {
    k_keep = keep_room;
  }
  int unlocked_count = 0;
  for (int p = 0; p < selected_count; ++p) {
    if (!buf->is_locked[p]) {
      ++unlocked_count;
    }
  }
  if (k_keep > unlocked_count) {
    k_keep = unlocked_count;
  }

  std::memset(buf->T_proj, 0,
              sizeof(double) * static_cast<size_t>(m_max) *
                static_cast<size_t>(m_max));
  *m_active = 0;
  int n_picked = 0;
  for (int p = 0; p < selected_count && n_picked < k_keep && *m_active < m_max; ++p) {
    if (buf->is_locked[p]) {
      continue;
    }
    std::memcpy(buf->z,
                buf->B_v + static_cast<int64_t>(p) * n,
                sizeof(double) * static_cast<size_t>(n));
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    const int accepted = block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
    n_picked += accepted;
  }

  const int tail_start = *m_active;
  int tail_accepted = 0;
  for (int p = 0; p < selected_count && tail_accepted < block_size && *m_active < m_max; ++p) {
    if (buf->is_locked[p] || buf->ritz_res[p] <= 100.0 * DBL_EPSILON) {
      continue;
    }
    const int idx = buf->selected[p];
    for (int row = 0; row < n; ++row) {
      buf->z[row] = buf->B_av[static_cast<int64_t>(p) * n + row] -
        buf->theta[idx] * buf->B_v[static_cast<int64_t>(p) * n + row];
    }
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    tail_accepted += block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
  }
  for (int attempt = 0; tail_accepted == 0 && attempt < n + block_size && *m_active < m_max; ++attempt) {
    std::memset(buf->z, 0, sizeof(double) * static_cast<size_t>(n));
    const int idx_basis = ((restart_idx + 1) * 17 + attempt * 31) % n;
    buf->z[idx_basis < 0 ? -idx_basis : idx_basis] = 1.0;
    stages->restart += native_timer_elapsed(timer);
    timer = native_timer_now();
    tail_accepted += block_accept_work_vector(
      V_out, n_locked, buf->V_active, m_active, m_max,
      buf->z, buf->tmp, n, ortho_passes_out
    );
    stages->reorthogonalization += native_timer_elapsed(timer);
    timer = native_timer_now();
  }
  stages->restart += native_timer_elapsed(timer);
  if (tail_accepted == 0) {
    return 1;
  }

  timer = native_timer_now();
  int rc = apply_active_block(impl, apply, n, 0, *m_active,
                              buf->V_active, buf->AV_active, workspace,
                              matvecs_out);
  stages->apply += native_timer_elapsed(timer);
  if (rc != 0) {
    return rc;
  }

  timer = native_timer_now();
  projection_update_self_block(buf->T_proj, m_max, buf->V_active,
                               buf->AV_active, n, 0, *m_active,
                               buf->S_eig);
  {
    const double elapsed = native_timer_elapsed(timer);
    stages->projected_solve += elapsed;
    stages->projection_update += elapsed;
  }
  *previous_block_start = 0;
  *previous_block_cols = tail_start;
  *last_block_start = tail_start;
  *last_block_cols = tail_accepted;
  *restarts_out = restart_idx + 1;
  return 0;
}

static int block_lanczos_finalize_return(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int target_kind,
    double tol,
    double norm_a,
    int selected_count_final,
    bool have_last_rr,
    ThickRestartBuffers* buf,
    EigencoreWorkspace* workspace,
    NativeBlockStageSeconds* stages,
    const BlockLanczosBestSnapshot& best,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* n_locked,
    int* matvecs_out) {
  int n_returned = *n_locked;
  if (best.filled) {
    std::memcpy(V_out, best.V.data(),
                sizeof(double) * static_cast<size_t>(n) *
                  static_cast<size_t>(k_target));
    std::memcpy(lambda_out, best.lambda.data(),
                sizeof(double) * static_cast<size_t>(k_target));
    std::memcpy(residuals_out, best.residuals.data(),
                sizeof(double) * static_cast<size_t>(k_target));
    std::memcpy(converged_out, best.converged.data(),
                sizeof(int) * static_cast<size_t>(k_target));
    *n_locked = best.locked_prefix;
    n_returned = k_target;
  } else if (have_last_rr && n_returned < k_target) {
    for (int p = 0; p < selected_count_final && n_returned < k_target; ++p) {
      if (buf->is_locked[p]) {
        continue;
      }
      const int idx = buf->selected[p];
      const double* vec = buf->B_v + static_cast<int64_t>(p) * n;
      if (!vector_is_independent_from_locked(V_out, n_returned, vec, n)) {
        continue;
      }
      std::memcpy(V_out + static_cast<int64_t>(n_returned) * n,
                  vec,
                  sizeof(double) * static_cast<size_t>(n));
      lambda_out[n_returned] = buf->theta[idx];
      residuals_out[n_returned] = buf->ritz_res[p];
      converged_out[n_returned] = 0;
      ++n_returned;
    }
  }
  if (n_returned != k_target) {
    return 0;
  }

  const int locked_prefix = *n_locked;
  const int full_best_snapshot = best.filled && best.nconv >= k_target;
  const int polish_offset = (!full_best_snapshot &&
                             locked_prefix > 0 && locked_prefix < k_target) ?
    locked_prefix : 0;
  const int polish_count = k_target - polish_offset;
  int polished_converged = 0;
  auto timer = native_timer_now();
  int polish_status = 0;
  if (polish_count > 0) {
    polish_status = final_polish_block_ritz(
      impl, apply, n, polish_count, target_kind, tol, norm_a,
      V_out + static_cast<int64_t>(polish_offset) * n,
      lambda_out + polish_offset,
      residuals_out + polish_offset,
      converged_out + polish_offset,
      &polished_converged, buf, workspace, matvecs_out
    );
  } else {
    polished_converged = 0;
  }
  {
    const double elapsed = native_timer_elapsed(timer);
    stages->ritz_residual += elapsed;
    stages->ritz_final_polish += elapsed;
  }
  if (polish_status != 0) {
    return polish_status;
  }
  *n_locked = polish_offset + polished_converged;
  if (best.filled && polish_offset > 0 && *n_locked >= k_target) {
    timer = native_timer_now();
    polish_status = final_polish_block_ritz(
      impl, apply, n, k_target, target_kind, tol, norm_a,
      V_out, lambda_out, residuals_out, converged_out,
      &polished_converged, buf, workspace, matvecs_out
    );
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_final_polish += elapsed;
    }
    if (polish_status != 0) {
      return polish_status;
    }
    *n_locked = polished_converged;
  }
  return 0;
}

static int native_block_thick_restart_lanczos_run(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int block_size,
    int target_kind,
    double tol,
    int max_restarts,
    double norm_a,
    int apply_ritz_vectors,
    const double* start_block,
    double* V_out,
    double* lambda_out,
    double* residuals_out,
    int* converged_out,
    int* n_locked_out,
    int* iterations_out,
    int* matvecs_out,
    int* restarts_out,
    int* m_active_final_out,
    int* locking_events_out,
    int* ortho_passes_out,
    int64_t* operator_allocations_out,
    int64_t* operator_bytes_allocated_out,
    NativeBlockStageSeconds* stage_out,
    NativeBlockRestartHistory* history
) {
  *n_locked_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *restarts_out = 0;
  *m_active_final_out = 0;
  *locking_events_out = 0;
  *ortho_passes_out = 0;
  *operator_allocations_out = 0;
  *operator_bytes_allocated_out = 0;
  if (stage_out != nullptr) {
    *stage_out = NativeBlockStageSeconds();
  }
  if (history != nullptr) {
    history->length = 0;
  }
  NativeBlockStageSeconds stage_local;
  NativeBlockStageSeconds* stages = (stage_out != nullptr) ? stage_out : &stage_local;
  for (int i = 0; i < k_target; ++i) {
    lambda_out[i] = 0.0;
    residuals_out[i] = R_PosInf;
    converged_out[i] = 0;
    std::memset(V_out + static_cast<int64_t>(i) * n, 0,
                sizeof(double) * static_cast<size_t>(n));
  }

  ThickRestartBuffers buf;
  if (trl_buffers_alloc(&buf, n, k_target, m_max, block_size) != 0) {
    return -2;
  }
  BlockLanczosBestSnapshot best(n, k_target);

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  int m_active = 0;
  int n_locked = 0;
  int previous_block_start = 0;
  int previous_block_cols = 0;
  int last_block_start = 0;
  auto timer = native_timer_now();
  int last_block_cols = block_accept_columns_blas3(
    start_block, n, block_size, V_out, n_locked,
    buf.V_active, &m_active, m_max, buf.Z_block, block_size,
    buf.coeff_block, buf.tmp, n,
    block_size, ortho_passes_out
  );
  stages->reorthogonalization += native_timer_elapsed(timer);
  if (last_block_cols == 0) {
    for (int attempt = 0; attempt < block_size && m_active < m_max; ++attempt) {
      std::memset(buf.z, 0, sizeof(double) * static_cast<size_t>(n));
      buf.z[attempt % n] = 1.0;
      timer = native_timer_now();
      last_block_cols += block_accept_work_vector(
        V_out, n_locked, buf.V_active, &m_active, m_max,
        buf.z, buf.tmp, n, ortho_passes_out
      );
      stages->reorthogonalization += native_timer_elapsed(timer);
    }
  }
  timer = native_timer_now();
  int rc = apply_active_block(impl, apply, n, 0, last_block_cols,
                              buf.V_active, buf.AV_active, &workspace,
                              matvecs_out);
  stages->apply += native_timer_elapsed(timer);
  if (rc != 0) {
    trl_buffers_free(&buf);
    return rc;
  }
  timer = native_timer_now();
  projection_update_self_block(buf.T_proj, m_max, buf.V_active, buf.AV_active,
                               n, 0, last_block_cols, buf.coeff_block);
  {
    const double elapsed = native_timer_elapsed(timer);
    stages->projected_solve += elapsed;
    stages->projection_update += elapsed;
  }

  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int selected_count_final = 0;
  int restart_idx = 0;
  bool have_last_rr = false;

  for (; restart_idx <= max_restarts; ++restart_idx) {
    rc = block_lanczos_expand_basis_to_budget(
      impl, apply, n, m_max, block_size, V_out, n_locked, &buf, &workspace,
      stages, &m_active, &previous_block_start, &previous_block_cols,
      &last_block_start, &last_block_cols, iterations_out, matvecs_out,
      ortho_passes_out
    );
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }

    if (m_active < 1) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    const int remaining_before_lock = k_target - n_locked;
    int pad = block_size > 4 ? block_size : 4;
    if (pad > k_target) {
      pad = k_target;
    }
    int selected_count = remaining_before_lock + pad;
    if (selected_count < remaining_before_lock) {
      selected_count = remaining_before_lock;
    }
    if (selected_count > m_active) {
      selected_count = m_active;
    }
    if (selected_count > buf.selected_capacity) {
      selected_count = buf.selected_capacity;
    }
    if (selected_count < 1) {
      selected_count = 1;
    }
    selected_count_final = selected_count;

    timer = native_timer_now();
    projection_copy_upper_compact(buf.T_proj, m_max, buf.S_eig, m_active);
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projection_copy += elapsed;
    }
    timer = native_timer_now();
    rc = symmetric_eigen_inplace(buf.S_eig, m_active, buf.theta,
                                 buf.dsyev_work, buf.dsyev_lwork,
                                 buf.dsyevd_iwork, buf.dsyevd_liwork);
    if (rc == 0) {
      selected_sorted_ritz_indices(buf.theta, m_active, selected_count, target_kind, buf.selected);
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->projected_solve += elapsed;
      stages->projected_eigensolve += elapsed;
    }
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }
    have_last_rr = true;

    timer = native_timer_now();
    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      std::memcpy(buf.S_selected + static_cast<int64_t>(p) * m_active,
                  buf.S_eig + static_cast<int64_t>(idx) * m_active,
                  sizeof(double) * static_cast<size_t>(m_active));
      buf.ritz_res[p] = R_PosInf;
      buf.is_locked[p] = 0;
    }
    stages->selected_vector_copy += native_timer_elapsed(timer);

    timer = native_timer_now();
    if (selected_count <= 32) {
      combine_basis_columns_small(buf.V_active, n, m_active,
                                  buf.S_selected, m_active,
                                  selected_count, buf.B_v);
    } else {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &selected_count, &m_active,
                      &one, buf.V_active, &n, buf.S_selected, &m_active,
                      &zero, buf.B_v, &n FCONE FCONE);
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_vector_form += elapsed;
    }

    timer = native_timer_now();
    rc = 0;
    if (apply_ritz_vectors) {
      rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, selected_count,
                 buf.B_v, n, 1.0, 0.0, buf.B_av, n, &workspace);
      if (rc == 0 && matvecs_out != nullptr) {
        ++(*matvecs_out);
      }
    } else {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &selected_count, &m_active,
                      &one, buf.AV_active, &n, buf.S_selected, &m_active,
                      &zero, buf.B_av, &n FCONE FCONE);
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_operator_apply += elapsed;
    }
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }

    timer = native_timer_now();
    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      long double s = 0.0L;
      const double* av = buf.B_av + static_cast<int64_t>(p) * n;
      const double* vec = buf.B_v + static_cast<int64_t>(p) * n;
      for (int row = 0; row < n; ++row) {
        const double diff = av[row] - buf.theta[idx] * vec[row];
        s += static_cast<long double>(diff) * diff;
      }
      buf.ritz_res[p] = sqrt(static_cast<double>(s));
    }
    {
      const double elapsed = native_timer_elapsed(timer);
      stages->ritz_residual += elapsed;
      stages->ritz_norm += elapsed;
    }

    const int wanted = k_target - n_locked;
    int history_slot = -1;
    if (history != nullptr && history->length < history->capacity) {
      history_slot = history->length;
      ++history->length;
      const int wanted_selected = wanted < selected_count ? wanted : selected_count;
      int nconv_wanted = 0;
      double max_residual = 0.0;
      double max_backward_error = 0.0;
      for (int p = 0; p < wanted_selected; ++p) {
        const int idx = buf.selected[p];
        const double scale_i = standard_eigen_lock_scale(
          norm_a, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
        );
        const double backward_error = buf.ritz_res[p] / scale_i;
        if (buf.ritz_res[p] > max_residual) {
          max_residual = buf.ritz_res[p];
        }
        if (backward_error > max_backward_error) {
          max_backward_error = backward_error;
        }
        if (buf.ritz_res[p] <= tol * scale_i) {
          ++nconv_wanted;
        }
      }
      history->restart[history_slot] = restart_idx;
      history->m_active[history_slot] = m_active;
      history->selected_count[history_slot] = selected_count;
      history->locked_before[history_slot] = n_locked;
      history->locked_after[history_slot] = n_locked;
      history->nconv_wanted[history_slot] = nconv_wanted;
      history->max_residual[history_slot] = max_residual;
      history->max_backward_error[history_slot] = max_backward_error;
    }

    timer = native_timer_now();
    int lock_now = 0;
    for (int p = 0; p < wanted && p < selected_count; ++p) {
      const int idx = buf.selected[p];
      const double scale_i = standard_eigen_lock_scale(
        norm_a, buf.theta[idx], buf.B_v + static_cast<int64_t>(p) * n, n
      );
      if (buf.ritz_res[p] <= tol * scale_i) {
        if (!vector_is_independent_from_locked(
              V_out, n_locked, buf.B_v + static_cast<int64_t>(p) * n, n)) {
          buf.is_locked[p] = 1;
          continue;
        }
        std::memcpy(V_out + static_cast<int64_t>(n_locked) * n,
                    buf.B_v + static_cast<int64_t>(p) * n,
                    sizeof(double) * static_cast<size_t>(n));
        lambda_out[n_locked] = buf.theta[idx];
        residuals_out[n_locked] = buf.ritz_res[p];
        converged_out[n_locked] = 1;
        buf.is_locked[p] = 1;
        ++n_locked;
        ++lock_now;
      } else {
        break;
      }
    }
    if (lock_now > 0) {
      ++(*locking_events_out);
    }
    if (history_slot >= 0) {
      history->locked_after[history_slot] = n_locked;
    }
    stages->locking += native_timer_elapsed(timer);

    block_lanczos_maybe_capture_best_snapshot(
      n, k_target, selected_count, n_locked, norm_a, tol, &buf,
      V_out, lambda_out, residuals_out, &best
    );

    if (n_locked >= k_target || restart_idx == max_restarts) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    rc = block_lanczos_restart_with_continuation_tail(
      impl, apply, n, k_target, m_max, block_size, restart_idx, selected_count,
      n_locked, &buf, &workspace, stages, &m_active, &previous_block_start,
      &previous_block_cols, &last_block_start, &last_block_cols, V_out,
      matvecs_out, restarts_out, ortho_passes_out
    );
    if (rc == 1) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }
    if (rc != 0) {
      trl_buffers_free(&buf);
      return rc;
    }
  }

  rc = block_lanczos_finalize_return(
    impl, apply, n, k_target, target_kind, tol, norm_a, selected_count_final,
    have_last_rr, &buf, &workspace, stages, best, V_out, lambda_out,
    residuals_out, converged_out, &n_locked, matvecs_out
  );
  if (rc != 0) {
    trl_buffers_free(&buf);
    return rc;
  }
  *n_locked_out = n_locked;
  *m_active_final_out = m_active;
  *operator_allocations_out = workspace.allocation_count;
  *operator_bytes_allocated_out = workspace.bytes_allocated;
  trl_buffers_free(&buf);
  return 0;
}

static SEXP block_lanczos_pack_result(int n, int k_target, const double* V,
                                      const double* lambda,
                                      const double* residuals,
                                      const int* converged, int nconv,
                                      int iterations, int matvecs,
                                      int m_active_final) {
  return trl_pack_result(n, k_target, V, lambda, residuals, converged, nconv,
                         iterations, matvecs, 0, m_active_final);
}

static SEXP block_thick_lanczos_pack_result(int n, int k_target, const double* V,
                                            const double* lambda,
                                            const double* residuals,
                                            const int* converged, int n_locked,
                                            int iterations, int matvecs,
                                            int restarts, int m_active_final,
                                            int locking_events, int ortho_passes,
                                            int block_size,
                                            int64_t operator_allocations,
                                            int64_t operator_bytes_allocated,
                                            const NativeBlockStageSeconds* stage_seconds,
                                            const NativeBlockRestartHistory* history) {
  SEXP values_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k_target));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k_target));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k_target));
  std::memcpy(REAL(values_), lambda, sizeof(double) * static_cast<size_t>(k_target));
  std::memcpy(REAL(vectors_), V,
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(k_target));
  std::memcpy(REAL(residuals_), residuals, sizeof(double) * static_cast<size_t>(k_target));
  for (int i = 0; i < k_target; ++i) {
    LOGICAL(converged_)[i] = converged[i] ? TRUE : FALSE;
  }

  SEXP stage_ = PROTECT(allocVector(REALSXP, 15));
  REAL(stage_)[0] = stage_seconds != nullptr ? stage_seconds->apply : 0.0;
  REAL(stage_)[1] = stage_seconds != nullptr ? stage_seconds->recurrence : 0.0;
  REAL(stage_)[2] = stage_seconds != nullptr ? stage_seconds->reorthogonalization : 0.0;
  REAL(stage_)[3] = stage_seconds != nullptr ? stage_seconds->projected_solve : 0.0;
  REAL(stage_)[4] = stage_seconds != nullptr ? stage_seconds->projection_update : 0.0;
  REAL(stage_)[5] = stage_seconds != nullptr ? stage_seconds->projection_copy : 0.0;
  REAL(stage_)[6] = stage_seconds != nullptr ? stage_seconds->projected_eigensolve : 0.0;
  REAL(stage_)[7] = stage_seconds != nullptr ? stage_seconds->selected_vector_copy : 0.0;
  REAL(stage_)[8] = stage_seconds != nullptr ? stage_seconds->ritz_residual : 0.0;
  REAL(stage_)[9] = stage_seconds != nullptr ? stage_seconds->ritz_vector_form : 0.0;
  REAL(stage_)[10] = stage_seconds != nullptr ? stage_seconds->ritz_operator_apply : 0.0;
  REAL(stage_)[11] = stage_seconds != nullptr ? stage_seconds->ritz_norm : 0.0;
  REAL(stage_)[12] = stage_seconds != nullptr ? stage_seconds->ritz_final_polish : 0.0;
  REAL(stage_)[13] = stage_seconds != nullptr ? stage_seconds->locking : 0.0;
  REAL(stage_)[14] = stage_seconds != nullptr ? stage_seconds->restart : 0.0;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 15));
  SET_STRING_ELT(stage_names_, 0, mkChar("apply"));
  SET_STRING_ELT(stage_names_, 1, mkChar("recurrence"));
  SET_STRING_ELT(stage_names_, 2, mkChar("reorthogonalization"));
  SET_STRING_ELT(stage_names_, 3, mkChar("projected_solve"));
  SET_STRING_ELT(stage_names_, 4, mkChar("projection_update"));
  SET_STRING_ELT(stage_names_, 5, mkChar("projection_copy"));
  SET_STRING_ELT(stage_names_, 6, mkChar("projected_eigensolve"));
  SET_STRING_ELT(stage_names_, 7, mkChar("selected_vector_copy"));
  SET_STRING_ELT(stage_names_, 8, mkChar("ritz_residual"));
  SET_STRING_ELT(stage_names_, 9, mkChar("ritz_vector_form"));
  SET_STRING_ELT(stage_names_, 10, mkChar("ritz_operator_apply"));
  SET_STRING_ELT(stage_names_, 11, mkChar("ritz_norm"));
  SET_STRING_ELT(stage_names_, 12, mkChar("ritz_final_polish"));
  SET_STRING_ELT(stage_names_, 13, mkChar("locking"));
  SET_STRING_ELT(stage_names_, 14, mkChar("restart"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  const int history_length =
    (history != nullptr && history->length > 0) ? history->length : 0;
  SEXP history_ = PROTECT(allocVector(VECSXP, 8));
  SEXP history_restart_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_m_active_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_selected_count_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_locked_before_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_locked_after_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_nconv_wanted_ = PROTECT(allocVector(INTSXP, history_length));
  SEXP history_max_residual_ = PROTECT(allocVector(REALSXP, history_length));
  SEXP history_max_backward_error_ = PROTECT(allocVector(REALSXP, history_length));
  if (history_length > 0) {
    std::memcpy(INTEGER(history_restart_), history->restart,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_m_active_), history->m_active,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_selected_count_), history->selected_count,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_locked_before_), history->locked_before,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_locked_after_), history->locked_after,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(INTEGER(history_nconv_wanted_), history->nconv_wanted,
                sizeof(int) * static_cast<size_t>(history_length));
    std::memcpy(REAL(history_max_residual_), history->max_residual,
                sizeof(double) * static_cast<size_t>(history_length));
    std::memcpy(REAL(history_max_backward_error_), history->max_backward_error,
                sizeof(double) * static_cast<size_t>(history_length));
  }
  SET_VECTOR_ELT(history_, 0, history_restart_);
  SET_VECTOR_ELT(history_, 1, history_m_active_);
  SET_VECTOR_ELT(history_, 2, history_selected_count_);
  SET_VECTOR_ELT(history_, 3, history_locked_before_);
  SET_VECTOR_ELT(history_, 4, history_locked_after_);
  SET_VECTOR_ELT(history_, 5, history_nconv_wanted_);
  SET_VECTOR_ELT(history_, 6, history_max_residual_);
  SET_VECTOR_ELT(history_, 7, history_max_backward_error_);
  SEXP history_names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(history_names_, 0, mkChar("restart"));
  SET_STRING_ELT(history_names_, 1, mkChar("m_active"));
  SET_STRING_ELT(history_names_, 2, mkChar("selected_count"));
  SET_STRING_ELT(history_names_, 3, mkChar("locked_before"));
  SET_STRING_ELT(history_names_, 4, mkChar("locked_after"));
  SET_STRING_ELT(history_names_, 5, mkChar("nconv_wanted"));
  SET_STRING_ELT(history_names_, 6, mkChar("max_residual"));
  SET_STRING_ELT(history_names_, 7, mkChar("max_backward_error"));
  setAttrib(history_, R_NamesSymbol, history_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, residuals_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(n_locked));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(restarts));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(m_active_final));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(locking_events));
  SET_VECTOR_ELT(out_, 10, ScalarInteger(ortho_passes));
  SET_VECTOR_ELT(out_, 11, ScalarInteger(block_size));
  SET_VECTOR_ELT(out_, 12, ScalarReal(static_cast<double>(operator_allocations)));
  SET_VECTOR_ELT(out_, 13, ScalarReal(static_cast<double>(operator_bytes_allocated)));
  SET_VECTOR_ELT(out_, 14, stage_);
  SET_VECTOR_ELT(out_, 15, history_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("residuals"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("n_locked"));
  SET_STRING_ELT(names_, 5, mkChar("iterations"));
  SET_STRING_ELT(names_, 6, mkChar("matvecs"));
  SET_STRING_ELT(names_, 7, mkChar("restarts"));
  SET_STRING_ELT(names_, 8, mkChar("m_active_final"));
  SET_STRING_ELT(names_, 9, mkChar("locking_events"));
  SET_STRING_ELT(names_, 10, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 11, mkChar("block"));
  SET_STRING_ELT(names_, 12, mkChar("operator_allocations"));
  SET_STRING_ELT(names_, 13, mkChar("operator_bytes_allocated"));
  SET_STRING_ELT(names_, 14, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 15, mkChar("restart_history"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(18);
  return out_;
}

static int lobpcg_orthonormalize(const double* S, int n, int cols, double tol,
                                 double* Q, double* tmp) {
  int rank = 0;
  for (int col = 0; col < cols; ++col) {
    double* q_col = Q + static_cast<int64_t>(rank) * n;
    std::memcpy(q_col, S + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
    trl_orthogonalise(nullptr, 0, Q, rank, q_col, tmp, n);
    const double q_norm = trl_norm2(q_col, n);
    if (q_norm <= tol) {
      continue;
    }
    const double inv = 1.0 / q_norm;
    for (int row = 0; row < n; ++row) {
      q_col[row] *= inv;
    }
    ++rank;
  }
  return rank;
}

static int lobpcg_b_orthonormalize_apply(void* b_impl,
                                         EigencoreApplyFn b_apply,
                                         const double* S,
                                         int n,
                                         int cols,
                                         double tol,
                                         double* Q,
                                         double* BQ,
                                         double* coeff,
                                         double* bq_col,
                                         EigencoreWorkspace* workspace) {
  int rank = 0;
  for (int col = 0; col < cols; ++col) {
    double* q_col = Q + static_cast<int64_t>(rank) * n;
    std::memcpy(q_col, S + static_cast<int64_t>(col) * n,
                sizeof(double) * static_cast<size_t>(n));
    int status = b_apply(b_impl, EIGENCORE_TRANSPOSE_NONE, 1, q_col, n,
                         1.0, 0.0, bq_col, n, workspace);
    if (status != 0) {
      return status;
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int prev = 0; prev < rank; ++prev) {
        const double* q_prev = Q + static_cast<int64_t>(prev) * n;
        long double dot = 0.0L;
        for (int row = 0; row < n; ++row) {
          dot += static_cast<long double>(q_prev[row]) * bq_col[row];
        }
        coeff[prev] = static_cast<double>(dot);
      }
      for (int prev = 0; prev < rank; ++prev) {
        const double c = coeff[prev];
        if (c == 0.0) {
          continue;
        }
        const double* q_prev = Q + static_cast<int64_t>(prev) * n;
        const double* bq_prev = BQ + static_cast<int64_t>(prev) * n;
        for (int row = 0; row < n; ++row) {
          q_col[row] -= c * q_prev[row];
          bq_col[row] -= c * bq_prev[row];
        }
      }
    }

    long double norm_sq = 0.0L;
    for (int row = 0; row < n; ++row) {
      norm_sq += static_cast<long double>(q_col[row]) * bq_col[row];
    }
    if (norm_sq <= 0.0L) {
      continue;
    }
    const double q_norm = sqrt(static_cast<double>(norm_sq));
    if (q_norm <= tol) {
      continue;
    }
    const double inv = 1.0 / q_norm;
    double* bq_dst = BQ + static_cast<int64_t>(rank) * n;
    for (int row = 0; row < n; ++row) {
      q_col[row] *= inv;
      bq_dst[row] = bq_col[row] * inv;
    }
    ++rank;
  }
  return rank;
}

static int lobpcg_apply_tridiagonal_preconditioner(
    const double* lower,
    const double* diag,
    const double* upper,
    int n,
    int cols,
    const double* R,
    double* W,
    double* cprime,
    double* dprime) {
  if (n < 1) {
    return 0;
  }
  if (fabs(diag[0]) <= DBL_EPSILON) {
    return -5;
  }
  dprime[0] = diag[0];
  if (n > 1) {
    cprime[0] = upper[0] / dprime[0];
  }
  for (int i = 1; i < n; ++i) {
    dprime[i] = diag[i] - lower[i - 1] * cprime[i - 1];
    if (fabs(dprime[i]) <= DBL_EPSILON) {
      return -5;
    }
    if (i < n - 1) {
      cprime[i] = upper[i] / dprime[i];
    }
  }

  for (int col = 0; col < cols; ++col) {
    const double* rhs = R + static_cast<int64_t>(col) * n;
    double* out = W + static_cast<int64_t>(col) * n;
    out[0] = rhs[0] / dprime[0];
    for (int i = 1; i < n; ++i) {
      out[i] = (rhs[i] - lower[i - 1] * out[i - 1]) / dprime[i];
    }
    for (int i = n - 2; i >= 0; --i) {
      out[i] -= cprime[i] * out[i + 1];
    }
  }
  return 0;
}

static int extract_shifted_symmetric_tridiagonal_from_csc(
    const int* i,
    const int* p,
    const double* x,
    int n,
    double shift,
    double* lower,
    double* diag,
    double* upper) {
  if (n < 1) {
    return -1;
  }
  std::memset(diag, 0, sizeof(double) * static_cast<size_t>(n));
  if (n > 1) {
    std::memset(lower, 0, sizeof(double) * static_cast<size_t>(n - 1));
    std::memset(upper, 0, sizeof(double) * static_cast<size_t>(n - 1));
  }

  for (int col = 0; col < n; ++col) {
    for (int pos = p[col]; pos < p[col + 1]; ++pos) {
      const int row = i[pos];
      const double value = x[pos];
      const int distance = row > col ? row - col : col - row;
      if (distance > 1) {
        return -6;
      }
      if (row == col) {
        diag[col] += value;
      } else if (row == col + 1) {
        lower[col] += value;
      } else {
        upper[row] += value;
      }
    }
  }

  for (int row = 0; row < n; ++row) {
    diag[row] += shift;
  }
  for (int row = 0; row < n - 1; ++row) {
    const double scale = fmax(fmax(fabs(lower[row]), fabs(upper[row])), 1.0);
    if (fabs(lower[row] - upper[row]) > 1e-12 * scale) {
      return -7;
    }
  }
  return 0;
}

static int lobpcg_project_constraints_apply(void* b_impl,
                                            EigencoreApplyFn b_apply,
                                            const double* Qc,
                                            int constraint_rank,
                                            int n,
                                            int cols,
                                            double* X,
                                            double* coeff,
                                            double* BX_work,
                                            EigencoreWorkspace* workspace) {
  if (constraint_rank <= 0 || cols <= 0) {
    return 0;
  }
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  const double* metric_X = X;
  if (b_apply != nullptr) {
    const int status = b_apply(b_impl, EIGENCORE_TRANSPOSE_NONE, cols, X, n,
                               1.0, 0.0, BX_work, n, workspace);
    if (status != 0) {
      return status;
    }
    metric_X = BX_work;
  }
  F77_CALL(dgemm)(&trans_T, &trans_N, &constraint_rank, &cols, &n,
                  &one, const_cast<double*>(Qc), &n,
                  const_cast<double*>(metric_X), &n,
                  &zero, coeff, &constraint_rank FCONE FCONE);
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &cols, &constraint_rank,
                  &minus_one, const_cast<double*>(Qc), &n,
                  coeff, &constraint_rank,
                  &one, X, &n FCONE FCONE);
  return 0;
}

static int lobpcg_constraint_matrix(SEXP constraints_,
                                    int n,
                                    const double** constraints,
                                    int* constraint_cols) {
  *constraints = nullptr;
  *constraint_cols = 0;
  if (!isReal(constraints_)) {
    return -1;
  }
  SEXP dimC = getAttrib(constraints_, R_DimSymbol);
  if (dimC == R_NilValue || INTEGER(dimC)[0] != n) {
    return -1;
  }
  const int cols = INTEGER(dimC)[1];
  if (cols < 0) {
    return -1;
  }
  const double* values = REAL(constraints_);
  for (int64_t pos = 0; pos < static_cast<int64_t>(n) * cols; ++pos) {
    if (!R_FINITE(values[pos])) {
      return -1;
    }
  }
  *constraints = values;
  *constraint_cols = cols;
  return 0;
}

static SEXP lobpcg_pack_result(int n, int k, const double* X,
                               const double* values,
                               const double* residuals,
                               const int* converged,
                               const double* hist_max_residual,
                               const int* hist_nconv,
                               int iterations, int matvecs,
                               int preconditioner_calls,
                               int q_rank_final,
                               int constraints_rank) {
  SEXP values_ = PROTECT(allocVector(REALSXP, k));
  SEXP vectors_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));
  SEXP hist_res_ = PROTECT(allocVector(REALSXP, iterations));
  SEXP hist_nconv_ = PROTECT(allocVector(INTSXP, iterations));
  std::memcpy(REAL(values_), values, sizeof(double) * static_cast<size_t>(k));
  std::memcpy(REAL(vectors_), X, sizeof(double) * static_cast<size_t>(n) * k);
  std::memcpy(REAL(residuals_), residuals, sizeof(double) * static_cast<size_t>(k));
  for (int i = 0; i < k; ++i) {
    LOGICAL(converged_)[i] = converged[i] ? TRUE : FALSE;
  }
  for (int i = 0; i < iterations; ++i) {
    REAL(hist_res_)[i] = hist_max_residual[i];
    INTEGER(hist_nconv_)[i] = hist_nconv[i];
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 11));
  SET_VECTOR_ELT(out_, 0, values_);
  SET_VECTOR_ELT(out_, 1, vectors_);
  SET_VECTOR_ELT(out_, 2, residuals_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, hist_res_);
  SET_VECTOR_ELT(out_, 5, hist_nconv_);
  SET_VECTOR_ELT(out_, 6, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(preconditioner_calls));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(q_rank_final));
  SET_VECTOR_ELT(out_, 10, ScalarInteger(constraints_rank));
  SEXP names_ = PROTECT(allocVector(STRSXP, 11));
  SET_STRING_ELT(names_, 0, mkChar("values"));
  SET_STRING_ELT(names_, 1, mkChar("vectors"));
  SET_STRING_ELT(names_, 2, mkChar("residuals"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("history_max_relative_residual"));
  SET_STRING_ELT(names_, 5, mkChar("history_nconv"));
  SET_STRING_ELT(names_, 6, mkChar("iterations"));
  SET_STRING_ELT(names_, 7, mkChar("matvecs"));
  SET_STRING_ELT(names_, 8, mkChar("preconditioner_calls"));
  SET_STRING_ELT(names_, 9, mkChar("q_rank_final"));
  SET_STRING_ELT(names_, 10, mkChar("constraints_rank"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(8);
  return out_;
}

static int native_lobpcg_run(void* impl,
                             EigencoreApplyFn apply,
                             void* b_impl,
                             EigencoreApplyFn b_apply,
                             int n,
                             int k,
                             int maxit,
                             int target_kind,
                             double tol,
                             const double* start,
                             int use_tridiagonal_preconditioner,
                             const double* lower,
                             const double* diag,
                             const double* upper,
                             const double* constraints,
                             int constraint_cols,
                             double* X_out,
                             double* values_out,
                             double* residuals_out,
                             int* converged_out,
                             double* hist_max_residual,
                             int* hist_nconv,
                             int* iterations_out,
                             int* matvecs_out,
                             int* preconditioner_calls_out,
                             int* q_rank_final_out,
                             int* constraints_rank_out) {
  const int max_trial_cols = 3 * k;
  const int tmp_cols = constraint_cols > max_trial_cols ? constraint_cols : max_trial_cols;
  const size_t nk = static_cast<size_t>(n) * static_cast<size_t>(k);
  const size_t nt = static_cast<size_t>(n) * static_cast<size_t>(max_trial_cols);
  double* X = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* Xnext = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* P = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* AX = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* BX = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* BXnext = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* R = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* W = static_cast<double*>(std::calloc(nk, sizeof(double)));
  double* S = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* Q = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* AQ = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* BQ = static_cast<double*>(std::calloc(nt, sizeof(double)));
  double* H = static_cast<double*>(std::calloc(static_cast<size_t>(max_trial_cols) * max_trial_cols, sizeof(double)));
  double* selected_vectors = static_cast<double*>(std::calloc(static_cast<size_t>(max_trial_cols) * k, sizeof(double)));
  double* theta = static_cast<double*>(std::calloc(static_cast<size_t>(max_trial_cols), sizeof(double)));
  double* tmp = static_cast<double*>(std::calloc(static_cast<size_t>(tmp_cols), sizeof(double)));
  int* selected = static_cast<int*>(std::calloc(static_cast<size_t>(max_trial_cols), sizeof(int)));
  double* cprime = static_cast<double*>(std::calloc(static_cast<size_t>(n > 1 ? n - 1 : 1), sizeof(double)));
  double* dprime = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  const int dsyev_lwork_query = trl_dsyev_query(max_trial_cols);
  int dsyev_lwork = dsyev_lwork_query > 0 ? dsyev_lwork_query : 3 * max_trial_cols;
  double* dsyev_work = static_cast<double*>(std::calloc(static_cast<size_t>(dsyev_lwork), sizeof(double)));
  if (X == nullptr || Xnext == nullptr || P == nullptr || AX == nullptr ||
      BX == nullptr || BXnext == nullptr || R == nullptr || W == nullptr ||
      S == nullptr || Q == nullptr || AQ == nullptr || BQ == nullptr ||
      H == nullptr || selected_vectors == nullptr ||
      theta == nullptr || tmp == nullptr || selected == nullptr ||
      cprime == nullptr || dprime == nullptr || dsyev_work == nullptr) {
    std::free(X); std::free(Xnext); std::free(P); std::free(AX); std::free(R);
    std::free(BX); std::free(BXnext); std::free(W); std::free(S);
    std::free(Q); std::free(AQ); std::free(BQ); std::free(H);
    std::free(selected_vectors); std::free(theta); std::free(tmp);
    std::free(selected); std::free(cprime); std::free(dprime); std::free(dsyev_work);
    return -2;
  }
  auto cleanup = [&]() {
    std::free(X); std::free(Xnext); std::free(P); std::free(AX); std::free(R);
    std::free(BX); std::free(BXnext); std::free(W); std::free(S);
    std::free(Q); std::free(AQ); std::free(BQ); std::free(H);
    std::free(selected_vectors); std::free(theta); std::free(tmp);
    std::free(selected); std::free(cprime); std::free(dprime); std::free(dsyev_work);
  };

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const bool generalized = b_apply != nullptr;
  int constraint_rank = 0;
  *constraints_rank_out = 0;
  std::vector<double> constraint_q;
  std::vector<double> constraint_bq;
  std::vector<double> constraint_bscratch;
  std::vector<double> constraint_coeff;
  std::vector<double> constraint_bx;
  if (constraint_cols > 0) {
    constraint_q.assign(static_cast<size_t>(n) * constraint_cols, 0.0);
    std::memcpy(constraint_q.data(), constraints,
                sizeof(double) * static_cast<size_t>(n) * constraint_cols);
    constraint_coeff.assign(static_cast<size_t>(constraint_cols) * max_trial_cols, 0.0);
    if (generalized) {
      constraint_bq.assign(static_cast<size_t>(n) * constraint_cols, 0.0);
      constraint_bscratch.assign(static_cast<size_t>(n), 0.0);
      constraint_bx.assign(static_cast<size_t>(n) * max_trial_cols, 0.0);
      constraint_rank = lobpcg_b_orthonormalize_apply(
        b_impl, b_apply, constraint_q.data(), n, constraint_cols,
        100.0 * DBL_EPSILON, constraint_q.data(), constraint_bq.data(), tmp,
        constraint_bscratch.data(), &workspace);
    } else {
      constraint_rank = lobpcg_orthonormalize(
        constraint_q.data(), n, constraint_cols, 100.0 * DBL_EPSILON,
        constraint_q.data(), tmp);
    }
    if (constraint_rank < 0) {
      cleanup();
      return constraint_rank;
    }
    if (constraint_rank + k > n) {
      cleanup();
      return -9;
    }
    *constraints_rank_out = constraint_rank;
  }
  std::memcpy(S, start, sizeof(double) * nk);
  if (constraint_rank > 0) {
    const int status = lobpcg_project_constraints_apply(
      b_impl, b_apply, constraint_q.data(), constraint_rank, n, k, S,
      constraint_coeff.data(),
      generalized ? constraint_bx.data() : nullptr,
      &workspace);
    if (status != 0) {
      cleanup();
      return status;
    }
  }
  int q_rank = generalized
    ? lobpcg_b_orthonormalize_apply(b_impl, b_apply, S, n, k,
                                    100.0 * DBL_EPSILON, X, BX, tmp,
                                    BQ, &workspace)
    : lobpcg_orthonormalize(S, n, k, 100.0 * DBL_EPSILON, X, tmp);
  if (q_rank < k) {
    cleanup();
    return -4;
  }
  if (!generalized) {
    std::memcpy(BX, X, sizeof(double) * nk);
  }

  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int have_p = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *preconditioner_calls_out = 0;
  *q_rank_final_out = k;

  for (int iter = 0; iter < maxit; ++iter) {
    *iterations_out = iter + 1;
    int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, k, X, n,
                       1.0, 0.0, AX, n, &workspace);
    if (status != 0) {
      cleanup();
      return status;
    }
    ++(*matvecs_out);

    int nconv = 0;
    double max_relative = 0.0;
    for (int col = 0; col < k; ++col) {
      long double lambda = 0.0L;
      for (int row = 0; row < n; ++row) {
        lambda += static_cast<long double>(X[row + static_cast<int64_t>(col) * n]) *
                  AX[row + static_cast<int64_t>(col) * n];
      }
      values_out[col] = static_cast<double>(lambda);
      long double rn = 0.0L;
      for (int row = 0; row < n; ++row) {
        const double res = AX[row + static_cast<int64_t>(col) * n] -
                           values_out[col] * BX[row + static_cast<int64_t>(col) * n];
        R[row + static_cast<int64_t>(col) * n] = res;
        rn += static_cast<long double>(res) * res;
      }
      residuals_out[col] = sqrt(static_cast<double>(rn));
      const double scale = fabs(values_out[col]) > 1.0 ? fabs(values_out[col]) : 1.0;
      const double rel = residuals_out[col] / scale;
      converged_out[col] = rel <= tol ? 1 : 0;
      if (converged_out[col]) ++nconv;
      if (rel > max_relative) max_relative = rel;
    }
    hist_max_residual[iter] = max_relative;
    hist_nconv[iter] = nconv;
    if (nconv >= k || iter + 1 >= maxit) {
      break;
    }

    if (use_tridiagonal_preconditioner) {
      status = lobpcg_apply_tridiagonal_preconditioner(
        lower, diag, upper, n, k, R, W, cprime, dprime
      );
      if (status != 0) {
        cleanup();
        return status;
      }
      ++(*preconditioner_calls_out);
    } else {
      std::memcpy(W, R, sizeof(double) * nk);
    }

    int trial_cols = 0;
    std::memcpy(S + static_cast<int64_t>(trial_cols) * n, X, sizeof(double) * nk);
    trial_cols += k;
    std::memcpy(S + static_cast<int64_t>(trial_cols) * n, W, sizeof(double) * nk);
    trial_cols += k;
    if (have_p) {
      std::memcpy(S + static_cast<int64_t>(trial_cols) * n, P, sizeof(double) * nk);
      trial_cols += k;
    }
    if (constraint_rank > 0) {
      status = lobpcg_project_constraints_apply(
        b_impl, b_apply, constraint_q.data(), constraint_rank, n, trial_cols, S,
        constraint_coeff.data(),
        generalized ? constraint_bx.data() : nullptr,
        &workspace);
      if (status != 0) {
        cleanup();
        return status;
      }
    }
    std::memset(Q, 0, sizeof(double) * nt);
    q_rank = generalized
      ? lobpcg_b_orthonormalize_apply(b_impl, b_apply, S, n, trial_cols,
                                      100.0 * DBL_EPSILON, Q, BQ, tmp,
                                      BXnext, &workspace)
      : lobpcg_orthonormalize(S, n, trial_cols, 100.0 * DBL_EPSILON, Q, tmp);
    if (q_rank < k) {
      cleanup();
      return -4;
    }
    *q_rank_final_out = q_rank;

    status = apply(impl, EIGENCORE_TRANSPOSE_NONE, q_rank, Q, n,
                   1.0, 0.0, AQ, n, &workspace);
    if (status != 0) {
      cleanup();
      return status;
    }
    ++(*matvecs_out);

    F77_CALL(dgemm)(&trans_T, &trans_N, &q_rank, &q_rank, &n,
                    &one, Q, &n, AQ, &n,
                    &zero, H, &q_rank FCONE FCONE);
    for (int i = 0; i < q_rank; ++i) {
      for (int j = i + 1; j < q_rank; ++j) {
        const double avg = 0.5 * (H[i + j * q_rank] + H[j + i * q_rank]);
        H[i + j * q_rank] = avg;
        H[j + i * q_rank] = avg;
      }
    }
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    int lwork = dsyev_lwork;
    F77_CALL(dsyev)(&jobz, &uplo, &q_rank, H, &q_rank, theta,
                    dsyev_work, &lwork, &info FCONE FCONE);
    if (info != 0) {
      cleanup();
      return -3;
    }
    selected_ritz_indices(theta, q_rank, k, target_kind, selected);
    for (int col = 0; col < k; ++col) {
      const int idx = selected[col];
      for (int row = 0; row < q_rank; ++row) {
        selected_vectors[row + static_cast<int64_t>(col) * q_rank] =
          H[row + static_cast<int64_t>(idx) * q_rank];
      }
    }
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &q_rank,
                    &one, Q, &n, selected_vectors, &q_rank,
                    &zero, Xnext, &n FCONE FCONE);
    if (generalized) {
      F77_CALL(dgemm)(&trans_N, &trans_N, &n, &k, &q_rank,
                      &one, BQ, &n, selected_vectors, &q_rank,
                      &zero, BXnext, &n FCONE FCONE);
    }
    for (size_t pos = 0; pos < nk; ++pos) {
      P[pos] = Xnext[pos] - X[pos];
      X[pos] = Xnext[pos];
      if (generalized) {
        BX[pos] = BXnext[pos];
      } else {
        BX[pos] = Xnext[pos];
      }
    }
    have_p = 1;
  }

  std::memcpy(X_out, X, sizeof(double) * nk);
  cleanup();
  return 0;
}

extern "C" SEXP eigencore_lobpcg_dense(SEXP A_, SEXP k_, SEXP maxit_,
                                       SEXP target_kind_, SEXP tol_,
                                       SEXP start_, SEXP lower_, SEXP diag_,
                                       SEXP upper_, SEXP constraints_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  const int status = native_lobpcg_run(
    &impl, eigencore_dense_apply, nullptr, nullptr,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native dense LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_dense_dense_b(SEXP A_, SEXP B_, SEXP k_,
                                               SEXP maxit_, SEXP target_kind_,
                                               SEXP tol_, SEXP start_,
                                               SEXP lower_, SEXP diag_,
                                               SEXP upper_,
                                               SEXP constraints_) {
  if (!isReal(A_) || !isReal(B_) || !isReal(start_)) {
    error("A, B, and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimB == R_NilValue || dimS == R_NilValue) {
    error("A, B, and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(dimB)[0] != n || INTEGER(dimB)[1] != n ||
      INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  DenseColumnMajorOperator b_impl = {n, n, REAL(B_)};
  const int status = native_lobpcg_run(
    &impl, eigencore_dense_apply, &b_impl, eigencore_dense_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native dense generalized LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_dense_diagonal_b(SEXP A_, SEXP bdiag_,
                                                  SEXP bunit_, SEXP k_,
                                                  SEXP maxit_,
                                                  SEXP target_kind_,
                                                  SEXP tol_, SEXP start_,
                                                  SEXP lower_, SEXP diag_,
                                                  SEXP upper_,
                                                  SEXP constraints_) {
  if (!isReal(A_) || !isReal(bdiag_) || !isReal(start_)) {
    error("A, B diagonal, and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool unit = asLogical(bunit_) == TRUE;
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!unit && LENGTH(bdiag_) != n)) {
    error("non-conformable dense/diagonal generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  DiagonalOperator b_impl = {n, REAL(bdiag_), unit};
  const int status = native_lobpcg_run(
    &impl, eigencore_dense_apply, &b_impl, eigencore_diagonal_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native dense/diagonal generalized LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_dense_csc_b(SEXP A_, SEXP bi_, SEXP bp_,
                                             SEXP bx_, SEXP bdim_, SEXP k_,
                                             SEXP maxit_, SEXP target_kind_,
                                             SEXP tol_, SEXP start_,
                                             SEXP lower_, SEXP diag_,
                                             SEXP upper_,
                                             SEXP constraints_) {
  if (!isReal(A_) || !isInteger(bi_) || !isInteger(bp_) ||
      !isReal(bx_) || !isInteger(bdim_) || !isReal(start_)) {
    error("invalid dense/CSC generalized LOBPCG inputs");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || LENGTH(bdim_) != 2) {
    error("A and start must be matrices and B dim must have length 2");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(bdim_)[0] != n || INTEGER(bdim_)[1] != n ||
      INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense/CSC generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  CSCOperator b_impl = {
    INTEGER(bdim_)[0], INTEGER(bdim_)[1], INTEGER(bi_), INTEGER(bp_), REAL(bx_)
  };
  const int status = native_lobpcg_run(
    &impl, eigencore_dense_apply, &b_impl, eigencore_csc_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native dense/CSC generalized LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_csc_diagonal_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                                SEXP adim_, SEXP bdiag_,
                                                SEXP bunit_, SEXP k_,
                                                SEXP maxit_,
                                                SEXP target_kind_,
                                                SEXP tol_, SEXP start_,
                                                SEXP lower_, SEXP diag_,
                                                SEXP upper_,
                                                SEXP constraints_) {
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isReal(bdiag_) || !isReal(start_)) {
    error("invalid CSC/diagonal generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool unit = asLogical(bunit_) == TRUE;
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!unit && LENGTH(bdiag_) != n)) {
    error("non-conformable CSC/diagonal generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  CSCOperator impl = {
    INTEGER(adim_)[0], INTEGER(adim_)[1], INTEGER(ai_), INTEGER(ap_), REAL(ax_)
  };
  DiagonalOperator b_impl = {n, REAL(bdiag_), unit};
  const int status = native_lobpcg_run(
    &impl, eigencore_csc_apply, &b_impl, eigencore_diagonal_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native CSC/diagonal generalized LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_csc_csc_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                           SEXP adim_, SEXP bi_, SEXP bp_,
                                           SEXP bx_, SEXP bdim_, SEXP k_,
                                           SEXP maxit_, SEXP target_kind_,
                                           SEXP tol_, SEXP start_,
                                           SEXP lower_, SEXP diag_,
                                           SEXP upper_,
                                           SEXP constraints_) {
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isInteger(bi_) || !isInteger(bp_) ||
      !isReal(bx_) || !isInteger(bdim_) || !isReal(start_)) {
    error("invalid CSC/CSC generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2 || LENGTH(bdim_) != 2) {
    error("start must be a matrix and A/B dims must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(adim_)[1] != n || INTEGER(bdim_)[0] != n || INTEGER(bdim_)[1] != n ||
      INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC/CSC generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  CSCOperator impl = {
    INTEGER(adim_)[0], INTEGER(adim_)[1], INTEGER(ai_), INTEGER(ap_), REAL(ax_)
  };
  CSCOperator b_impl = {
    INTEGER(bdim_)[0], INTEGER(bdim_)[1], INTEGER(bi_), INTEGER(bp_), REAL(bx_)
  };
  const int status = native_lobpcg_run(
    &impl, eigencore_csc_apply, &b_impl, eigencore_csc_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native CSC/CSC generalized LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_diagonal_diagonal_b(SEXP adiag_, SEXP aunit_,
                                                     SEXP adim_, SEXP bdiag_,
                                                     SEXP bunit_, SEXP k_,
                                                     SEXP maxit_,
                                                     SEXP target_kind_,
                                                     SEXP tol_, SEXP start_,
                                                     SEXP lower_, SEXP diag_,
                                                     SEXP upper_,
                                                     SEXP constraints_) {
  if (!isReal(adiag_) || !isInteger(adim_) || !isReal(bdiag_) ||
      !isReal(start_)) {
    error("invalid diagonal/diagonal generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool a_unit = asLogical(aunit_) == TRUE;
  const bool b_unit = asLogical(bunit_) == TRUE;
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!a_unit && LENGTH(adiag_) != n) || (!b_unit && LENGTH(bdiag_) != n)) {
    error("non-conformable diagonal/diagonal generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  DiagonalOperator impl = {n, REAL(adiag_), a_unit};
  DiagonalOperator b_impl = {n, REAL(bdiag_), b_unit};
  const int status = native_lobpcg_run(
    &impl, eigencore_diagonal_apply, &b_impl, eigencore_diagonal_apply,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native diagonal/diagonal generalized LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

static SEXP lobpcg_run_matrix_free_b(void* impl,
                                     EigencoreApplyFn apply,
                                     int n,
                                     int k,
                                     int maxit,
                                     int target_kind,
                                     double tol,
                                     const double* start,
                                     SEXP B_apply_,
                                     SEXP lower_,
                                     SEXP diag_,
                                     SEXP upper_,
                                     SEXP constraints_,
                                     const char* label) {
  if (TYPEOF(B_apply_) != CLOSXP) {
    error("matrix-free B operator apply must be an R closure");
  }
  const int use_tridiag = LENGTH(diag_) > 0;
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  RApplyOperator b_impl = {n, B_apply_};
  const int status = native_lobpcg_run(
    impl, apply, &b_impl, eigencore_r_operator_apply,
    n, k, maxit, target_kind, tol, start,
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native %s matrix-free-B LOBPCG failed with status=%d", label, status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_dense_operator_b(SEXP A_, SEXP B_apply_,
                                                  SEXP k_, SEXP maxit_,
                                                  SEXP target_kind_,
                                                  SEXP tol_, SEXP start_,
                                                  SEXP lower_, SEXP diag_,
                                                  SEXP upper_,
                                                  SEXP constraints_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable dense/matrix-free generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  return lobpcg_run_matrix_free_b(
    &impl, eigencore_dense_apply, n, k, maxit, target_kind, tol, REAL(start_),
    B_apply_, lower_, diag_, upper_, constraints_, "dense");
}

extern "C" SEXP eigencore_lobpcg_csc_operator_b(SEXP ai_, SEXP ap_, SEXP ax_,
                                                SEXP adim_, SEXP B_apply_,
                                                SEXP k_, SEXP maxit_,
                                                SEXP target_kind_,
                                                SEXP tol_, SEXP start_,
                                                SEXP lower_, SEXP diag_,
                                                SEXP upper_,
                                                SEXP constraints_) {
  if (!isInteger(ai_) || !isInteger(ap_) || !isReal(ax_) ||
      !isInteger(adim_) || !isReal(start_)) {
    error("invalid CSC/matrix-free generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC/matrix-free generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  CSCOperator impl = {
    INTEGER(adim_)[0], INTEGER(adim_)[1], INTEGER(ai_), INTEGER(ap_), REAL(ax_)
  };
  return lobpcg_run_matrix_free_b(
    &impl, eigencore_csc_apply, n, k, maxit, target_kind, tol, REAL(start_),
    B_apply_, lower_, diag_, upper_, constraints_, "CSC");
}

extern "C" SEXP eigencore_lobpcg_diagonal_operator_b(SEXP adiag_,
                                                     SEXP aunit_,
                                                     SEXP adim_,
                                                     SEXP B_apply_,
                                                     SEXP k_,
                                                     SEXP maxit_,
                                                     SEXP target_kind_,
                                                     SEXP tol_,
                                                     SEXP start_,
                                                     SEXP lower_,
                                                     SEXP diag_,
                                                     SEXP upper_,
                                                     SEXP constraints_) {
  if (!isReal(adiag_) || !isInteger(adim_) || !isReal(start_)) {
    error("invalid diagonal/matrix-free generalized LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(adim_) != 2) {
    error("start must be a matrix and A dim must have length 2");
  }
  const int n = INTEGER(adim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  const bool a_unit = asLogical(aunit_) == TRUE;
  if (INTEGER(adim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k ||
      (!a_unit && LENGTH(adiag_) != n)) {
    error("non-conformable diagonal/matrix-free generalized LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  DiagonalOperator impl = {n, REAL(adiag_), a_unit};
  return lobpcg_run_matrix_free_b(
    &impl, eigencore_diagonal_apply, n, k, maxit, target_kind, tol, REAL(start_),
    B_apply_, lower_, diag_, upper_, constraints_, "diagonal");
}

extern "C" SEXP eigencore_lobpcg_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                     SEXP k_, SEXP maxit_, SEXP target_kind_,
                                     SEXP tol_, SEXP start_, SEXP lower_,
                                     SEXP diag_, SEXP upper_,
                                     SEXP constraints_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int use_tridiag = LENGTH(diag_) > 0;
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  const double* constraints = nullptr;
  int constraint_cols = 0;
  if (lobpcg_constraint_matrix(constraints_, n, &constraints, &constraint_cols) != 0) {
    error("constraints must be a double matrix with n rows");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_lobpcg_run(
    &impl, eigencore_csc_apply, nullptr, nullptr,
    n, k, maxit, target_kind, tol, REAL(start_),
    use_tridiag, LENGTH(lower_) ? REAL(lower_) : nullptr,
    LENGTH(diag_) ? REAL(diag_) : nullptr,
    LENGTH(upper_) ? REAL(upper_) : nullptr,
    constraints, constraint_cols,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native CSC LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_lobpcg_csc_shifted_tridiagonal(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP k_, SEXP maxit_,
    SEXP target_kind_, SEXP tol_, SEXP start_, SEXP shift_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC shifted-tridiagonal LOBPCG inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  const int k = static_cast<int>(asInteger(k_));
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n || INTEGER(dimS)[1] != k) {
    error("non-conformable CSC shifted-tridiagonal LOBPCG inputs");
  }
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const double shift = asReal(shift_);
  if (k < 1 || maxit < 1) error("k and maxit must be positive");
  if (!R_FINITE(shift) || shift < 0.0) {
    error("shift must be a finite non-negative scalar");
  }

  std::vector<double> lower(static_cast<size_t>(n > 1 ? n - 1 : 1), 0.0);
  std::vector<double> diag(static_cast<size_t>(n), 0.0);
  std::vector<double> upper(static_cast<size_t>(n > 1 ? n - 1 : 1), 0.0);
  const int tri_status = extract_shifted_symmetric_tridiagonal_from_csc(
    INTEGER(i_), INTEGER(p_), REAL(x_), n, shift,
    lower.data(), diag.data(), upper.data()
  );
  if (tri_status == -6) {
    error("CSC matrix is not tridiagonal");
  }
  if (tri_status == -7) {
    error("CSC matrix is not symmetric tridiagonal");
  }
  if (tri_status != 0) {
    error("failed to extract shifted tridiagonal preconditioner");
  }

  std::vector<double> X(static_cast<size_t>(n) * k, 0.0);
  std::vector<double> values(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  std::vector<double> hist_res(static_cast<size_t>(maxit), R_PosInf);
  std::vector<int> hist_nconv(static_cast<size_t>(maxit), 0);
  int iterations = 0, matvecs = 0, preconditioner_calls = 0, q_rank = 0, constraints_rank = 0;
  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_lobpcg_run(
    &impl, eigencore_csc_apply, nullptr, nullptr,
    n, k, maxit, target_kind, tol, REAL(start_),
    1, lower.data(), diag.data(), upper.data(),
    nullptr, 0,
    X.data(), values.data(), residuals.data(), converged.data(),
    hist_res.data(), hist_nconv.data(), &iterations, &matvecs,
    &preconditioner_calls, &q_rank, &constraints_rank);
  if (status != 0) {
    error("native CSC shifted-tridiagonal LOBPCG failed with status=%d", status);
  }
  return lobpcg_pack_result(n, k, X.data(), values.data(), residuals.data(),
                            converged.data(), hist_res.data(), hist_nconv.data(),
                            iterations, matvecs, preconditioner_calls, q_rank,
                            constraints_rank);
}

extern "C" SEXP eigencore_block_lanczos_dense(SEXP A_, SEXP k_,
                                              SEXP m_max_,
                                              SEXP block_size_,
                                              SEXP target_kind_,
                                              SEXP tol_,
                                              SEXP start_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable block Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k) error("m_max must be >= k");
  if (m_max > n) error("m_max must be <= nrow(A)");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int nconv = 0, iterations = 0, matvecs = 0, m_active = 0;

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  const int status = native_block_lanczos_run(
    &impl, eigencore_dense_apply, n, k, m_max, block_size, target_kind, tol,
    REAL(start_), V.data(), lambda.data(), residuals.data(), converged.data(),
    &nconv, &iterations, &matvecs, &m_active);
  if (status != 0) {
    error("native dense block Lanczos failed with status=%d", status);
  }
  return block_lanczos_pack_result(n, k, V.data(), lambda.data(),
                                   residuals.data(), converged.data(),
                                   nconv, iterations, matvecs, m_active);
}

extern "C" SEXP eigencore_block_lanczos_csc(SEXP i_, SEXP p_, SEXP x_,
                                            SEXP dim_, SEXP k_,
                                            SEXP m_max_,
                                            SEXP block_size_,
                                            SEXP target_kind_,
                                            SEXP tol_,
                                            SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC block Lanczos inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable CSC block Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k) error("m_max must be >= k");
  if (m_max > n) error("m_max must be <= nrow(A)");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int nconv = 0, iterations = 0, matvecs = 0, m_active = 0;

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_block_lanczos_run(
    &impl, eigencore_csc_apply, n, k, m_max, block_size, target_kind, tol,
    REAL(start_), V.data(), lambda.data(), residuals.data(), converged.data(),
    &nconv, &iterations, &matvecs, &m_active);
  if (status != 0) {
    error("native CSC block Lanczos failed with status=%d", status);
  }
  return block_lanczos_pack_result(n, k, V.data(), lambda.data(),
                                   residuals.data(), converged.data(),
                                   nconv, iterations, matvecs, m_active);
}

extern "C" SEXP eigencore_block_thick_restart_lanczos_dense(
    SEXP A_, SEXP k_, SEXP m_max_, SEXP block_size_,
    SEXP target_kind_, SEXP tol_, SEXP max_restarts_,
    SEXP norm_a_, SEXP start_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable block thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  const double norm_a = asReal(norm_a_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, restarts = 0, m_active = 0;
  int locking_events = 0, ortho_passes = 0;
  int64_t operator_allocations = 0, operator_bytes_allocated = 0;
  NativeBlockStageSeconds stage_seconds;
  const int history_capacity = max_restarts + 1;
  std::vector<int> history_restart(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_m_active(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_selected_count(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_before(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_after(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_nconv_wanted(static_cast<size_t>(history_capacity), 0);
  std::vector<double> history_max_residual(static_cast<size_t>(history_capacity), R_PosInf);
  std::vector<double> history_max_backward_error(static_cast<size_t>(history_capacity), R_PosInf);
  NativeBlockRestartHistory history;
  history.capacity = history_capacity;
  history.restart = history_restart.data();
  history.m_active = history_m_active.data();
  history.selected_count = history_selected_count.data();
  history.locked_before = history_locked_before.data();
  history.locked_after = history_locked_after.data();
  history.nconv_wanted = history_nconv_wanted.data();
  history.max_residual = history_max_residual.data();
  history.max_backward_error = history_max_backward_error.data();

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  const int status = native_block_thick_restart_lanczos_run(
    &impl, eigencore_dense_apply, n, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, 0, REAL(start_), V.data(), lambda.data(),
    residuals.data(), converged.data(), &n_locked, &iterations, &matvecs,
    &restarts, &m_active, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes_allocated, &stage_seconds, &history);
  if (status != 0) {
    error("native dense block thick-restart Lanczos failed with status=%d", status);
  }
  return block_thick_lanczos_pack_result(
    n, k, V.data(), lambda.data(), residuals.data(), converged.data(),
    n_locked, iterations, matvecs, restarts, m_active, locking_events,
    ortho_passes, block_size, operator_allocations, operator_bytes_allocated,
    &stage_seconds, &history);
}

extern "C" SEXP eigencore_block_thick_restart_lanczos_csc(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP k_,
    SEXP m_max_, SEXP block_size_, SEXP target_kind_,
    SEXP tol_, SEXP max_restarts_, SEXP norm_a_, SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC block thick-restart Lanczos inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue) {
    error("start must be a matrix");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || INTEGER(dimS)[0] != n) {
    error("non-conformable CSC block thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int block_size = static_cast<int>(asInteger(block_size_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  const double norm_a = asReal(norm_a_);
  if (k < 1) error("k must be >= 1");
  if (block_size < 1) error("block_size must be >= 1");
  if (INTEGER(dimS)[1] != block_size) error("start block has wrong number of columns");
  if (m_max < k + block_size) error("m_max must be >= k + block_size");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, restarts = 0, m_active = 0;
  int locking_events = 0, ortho_passes = 0;
  int64_t operator_allocations = 0, operator_bytes_allocated = 0;
  NativeBlockStageSeconds stage_seconds;
  const int history_capacity = max_restarts + 1;
  std::vector<int> history_restart(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_m_active(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_selected_count(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_before(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_locked_after(static_cast<size_t>(history_capacity), 0);
  std::vector<int> history_nconv_wanted(static_cast<size_t>(history_capacity), 0);
  std::vector<double> history_max_residual(static_cast<size_t>(history_capacity), R_PosInf);
  std::vector<double> history_max_backward_error(static_cast<size_t>(history_capacity), R_PosInf);
  NativeBlockRestartHistory history;
  history.capacity = history_capacity;
  history.restart = history_restart.data();
  history.m_active = history_m_active.data();
  history.selected_count = history_selected_count.data();
  history.locked_before = history_locked_before.data();
  history.locked_after = history_locked_after.data();
  history.nconv_wanted = history_nconv_wanted.data();
  history.max_residual = history_max_residual.data();
  history.max_backward_error = history_max_backward_error.data();

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_block_thick_restart_lanczos_run(
    &impl, eigencore_csc_apply, n, k, m_max, block_size, target_kind,
    tol, max_restarts, norm_a, 1, REAL(start_), V.data(), lambda.data(),
    residuals.data(), converged.data(), &n_locked, &iterations, &matvecs,
    &restarts, &m_active, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes_allocated, &stage_seconds, &history);
  if (status != 0) {
    error("native CSC block thick-restart Lanczos failed with status=%d", status);
  }
  return block_thick_lanczos_pack_result(
    n, k, V.data(), lambda.data(), residuals.data(), converged.data(),
    n_locked, iterations, matvecs, restarts, m_active, locking_events,
    ortho_passes, block_size, operator_allocations, operator_bytes_allocated,
    &stage_seconds, &history);
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

extern "C" SEXP eigencore_dense_eigen_certificate(SEXP A_, SEXP values_,
                                                  SEXP vectors_, SEXP B_,
                                                  SEXP tol_) {
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
    error("non-conformable dense eigen certificate inputs");
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

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double eps = DBL_EPSILON;
  const double tol = asReal(tol_);
  const double norm_A = frobenius_norm_dense(REAL(A_), n * n);
  const double norm_B = (B_ == R_NilValue) ? 1.0 : frobenius_norm_dense(REAL(B_), n * n);

  int protect_count = 0;
  SEXP residual_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));
  ++protect_count;
  SEXP Bv_ = R_NilValue;
  SEXP gram_ = PROTECT(allocMatrix(REALSXP, k, k));
  ++protect_count;
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  ++protect_count;
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  ++protect_count;
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  ++protect_count;
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));
  ++protect_count;

  F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                  &one, REAL(A_), &n, REAL(vectors_), &n,
                  &zero, REAL(residual_matrix_), &n FCONE FCONE);

  const double* bv = REAL(vectors_);
  if (B_ != R_NilValue) {
    Bv_ = PROTECT(allocMatrix(REALSXP, n, k));
    ++protect_count;
    F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                    &one, REAL(B_), &n, REAL(vectors_), &n,
                    &zero, REAL(Bv_), &n FCONE FCONE);
    bv = REAL(Bv_);
  }

  for (int col = 0; col < k; ++col) {
    const double lambda = REAL(values_)[col];
    const int offset = col * n;
    for (int row = 0; row < n; ++row) {
      REAL(residual_matrix_)[offset + row] -= lambda * bv[offset + row];
    }
    const double residual = column_norm(REAL(residual_matrix_), n, col);
    const double vector_norm = column_norm(REAL(vectors_), n, col);
    const double scale = fmax((norm_A + fabs(lambda) * norm_B) * fmax(vector_norm, eps), eps);
    const double backward = residual / scale;
    REAL(residuals_)[col] = residual;
    REAL(scale_)[col] = scale;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, REAL(vectors_), &n, bv, &n,
                  &zero, REAL(gram_), &k FCONE FCONE);
  const double orth = max_orthogonality_loss(REAL(gram_), k);

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  ++protect_count;
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, ScalarReal(norm_A));
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  ++protect_count;
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("scale"));
  SET_STRING_ELT(names_, 4, mkChar("converged"));
  SET_STRING_ELT(names_, 5, mkChar("norm_A"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(protect_count);
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

extern "C" SEXP eigencore_dense_svd_certificate(SEXP A_, SEXP d_,
                                                SEXP u_, SEXP v_, SEXP tol_) {
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
    error("non-conformable dense SVD certificate inputs");
  }

  const char notrans = 'N';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  const double eps = DBL_EPSILON;
  const double tol = asReal(tol_);
  const double norm_A = frobenius_norm_dense(REAL(A_), m * n);
  const double scale_value = fmax(norm_A, eps);

  SEXP left_matrix_ = PROTECT(allocMatrix(REALSXP, m, k));
  SEXP right_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP left_ = PROTECT(allocVector(REALSXP, k));
  SEXP right_ = PROTECT(allocVector(REALSXP, k));
  SEXP combined_ = PROTECT(allocVector(REALSXP, k));
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));

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
    const double left = column_norm(REAL(left_matrix_), m, col);
    const double right = column_norm(REAL(right_matrix_), n, col);
    const double combined = sqrt(left * left + right * right);
    const double backward = combined / scale_value;
    REAL(left_)[col] = left;
    REAL(right_)[col] = right;
    REAL(combined_)[col] = combined;
    REAL(scale_)[col] = scale_value;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  SEXP gram_u_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP gram_v_ = PROTECT(allocMatrix(REALSXP, k, k));
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &m,
                  &one, REAL(u_), &m, REAL(u_), &m,
                  &zero, REAL(gram_u_), &k FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, REAL(v_), &n, REAL(v_), &n,
                  &zero, REAL(gram_v_), &k FCONE FCONE);
  REAL(orth_)[0] = max_orthogonality_loss(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss(REAL(gram_v_), k);
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, left_);
  SET_VECTOR_ELT(out_, 1, right_);
  SET_VECTOR_ELT(out_, 2, combined_);
  SET_VECTOR_ELT(out_, 3, backward_);
  SET_VECTOR_ELT(out_, 4, orth_);
  SET_VECTOR_ELT(out_, 5, scale_);
  SET_VECTOR_ELT(out_, 6, converged_);
  SET_VECTOR_ELT(out_, 7, ScalarReal(norm_A));
  SET_VECTOR_ELT(out_, 8, ScalarReal(scale_value));
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("left"));
  SET_STRING_ELT(names_, 1, mkChar("right"));
  SET_STRING_ELT(names_, 2, mkChar("combined"));
  SET_STRING_ELT(names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 5, mkChar("scale"));
  SET_STRING_ELT(names_, 6, mkChar("converged"));
  SET_STRING_ELT(names_, 7, mkChar("norm_A"));
  SET_STRING_ELT(names_, 8, mkChar("scale_value"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(14);
  return out_;
}

static SEXP native_operator_svd_certificate_cached_av(void* impl,
                                                      EigencoreApplyFn apply,
                                                      int m,
                                                      int n,
                                                      double norm_A,
                                                      SEXP d_,
                                                      SEXP u_,
                                                      SEXP v_,
                                                      SEXP av_,
                                                      SEXP tol_);

extern "C" SEXP eigencore_dense_svd_certificate_cached_av(SEXP A_, SEXP d_,
                                                          SEXP u_, SEXP v_,
                                                          SEXP av_, SEXP tol_) {
  if (!isReal(A_) || !isReal(d_) || !isReal(u_) || !isReal(v_) || !isReal(av_)) {
    error("A, d, u, v, and Av must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  const double norm_A = frobenius_norm_dense(REAL(A_), m * n);
  return native_operator_svd_certificate_cached_av(
    &impl, eigencore_dense_apply, m, n, norm_A, d_, u_, v_, av_, tol_
  );
}

static SEXP native_operator_eigen_certificate(void* impl,
                                              EigencoreApplyFn apply,
                                              int n,
                                              double norm_A,
                                              SEXP values_,
                                              SEXP vectors_,
                                              SEXP tol_) {
  SEXP dimV = getAttrib(vectors_, R_DimSymbol);
  if (dimV == R_NilValue) {
    error("vectors must be a matrix");
  }
  const int rowsV = INTEGER(dimV)[0];
  const int k = INTEGER(dimV)[1];
  if (rowsV != n || LENGTH(values_) != k) {
    error("non-conformable native operator eigen certificate inputs");
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double eps = DBL_EPSILON;
  const double tol = asReal(tol_);
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};

  SEXP residual_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP gram_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));

  const int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, k,
                           REAL(vectors_), n, 1.0, 0.0,
                           REAL(residual_matrix_), n, &workspace);
  if (status != 0) {
    error("native operator eigen certificate apply failed with status=%d", status);
  }

  for (int col = 0; col < k; ++col) {
    const double lambda = REAL(values_)[col];
    const int offset = col * n;
    for (int row = 0; row < n; ++row) {
      REAL(residual_matrix_)[offset + row] -= lambda * REAL(vectors_)[offset + row];
    }
    const double residual = column_norm(REAL(residual_matrix_), n, col);
    const double vector_norm = column_norm(REAL(vectors_), n, col);
    const double scale = fmax((norm_A + fabs(lambda)) * fmax(vector_norm, eps), eps);
    const double backward = residual / scale;
    REAL(residuals_)[col] = residual;
    REAL(scale_)[col] = scale;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, REAL(vectors_), &n, REAL(vectors_), &n,
                  &zero, REAL(gram_), &k FCONE FCONE);
  const double orth = max_orthogonality_loss(REAL(gram_), k);

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, workspace_counters(&workspace));
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("scale"));
  SET_STRING_ELT(names_, 4, mkChar("converged"));
  SET_STRING_ELT(names_, 5, mkChar("workspace"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(8);
  return out_;
}

static SEXP native_operator_svd_certificate(void* impl,
                                            EigencoreApplyFn apply,
                                            int m,
                                            int n,
                                            double norm_A,
                                            SEXP d_,
                                            SEXP u_,
                                            SEXP v_,
                                            SEXP tol_) {
  SEXP dimU = getAttrib(u_, R_DimSymbol);
  SEXP dimV = getAttrib(v_, R_DimSymbol);
  if (dimU == R_NilValue || dimV == R_NilValue) {
    error("u and v must be matrices");
  }
  const int rowsU = INTEGER(dimU)[0];
  const int k = INTEGER(dimU)[1];
  const int rowsV = INTEGER(dimV)[0];
  const int colsV = INTEGER(dimV)[1];
  if (rowsU != m || rowsV != n || colsV != k || LENGTH(d_) != k) {
    error("non-conformable native operator SVD certificate inputs");
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double eps = DBL_EPSILON;
  const double tol = asReal(tol_);
  const double scale_value = fmax(norm_A, eps);
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};

  SEXP left_matrix_ = PROTECT(allocMatrix(REALSXP, m, k));
  SEXP right_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP left_ = PROTECT(allocVector(REALSXP, k));
  SEXP right_ = PROTECT(allocVector(REALSXP, k));
  SEXP combined_ = PROTECT(allocVector(REALSXP, k));
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));

  int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, k,
                     REAL(v_), n, 1.0, 0.0,
                     REAL(left_matrix_), m, &workspace);
  if (status != 0) {
    error("native operator SVD certificate apply failed with status=%d", status);
  }
  status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, k,
                 REAL(u_), m, 1.0, 0.0,
                 REAL(right_matrix_), n, &workspace);
  if (status != 0) {
    error("native operator SVD certificate adjoint apply failed with status=%d", status);
  }

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
    const double left = column_norm(REAL(left_matrix_), m, col);
    const double right = column_norm(REAL(right_matrix_), n, col);
    const double combined = sqrt(left * left + right * right);
    const double backward = combined / scale_value;
    REAL(left_)[col] = left;
    REAL(right_)[col] = right;
    REAL(combined_)[col] = combined;
    REAL(scale_)[col] = scale_value;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  SEXP gram_u_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP gram_v_ = PROTECT(allocMatrix(REALSXP, k, k));
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &m,
                  &one, REAL(u_), &m, REAL(u_), &m,
                  &zero, REAL(gram_u_), &k FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, REAL(v_), &n, REAL(v_), &n,
                  &zero, REAL(gram_v_), &k FCONE FCONE);
  REAL(orth_)[0] = max_orthogonality_loss(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss(REAL(gram_v_), k);
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, left_);
  SET_VECTOR_ELT(out_, 1, right_);
  SET_VECTOR_ELT(out_, 2, combined_);
  SET_VECTOR_ELT(out_, 3, backward_);
  SET_VECTOR_ELT(out_, 4, orth_);
  SET_VECTOR_ELT(out_, 5, scale_);
  SET_VECTOR_ELT(out_, 6, converged_);
  SET_VECTOR_ELT(out_, 7, ScalarReal(scale_value));
  SET_VECTOR_ELT(out_, 8, workspace_counters(&workspace));
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("left"));
  SET_STRING_ELT(names_, 1, mkChar("right"));
  SET_STRING_ELT(names_, 2, mkChar("combined"));
  SET_STRING_ELT(names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 5, mkChar("scale"));
  SET_STRING_ELT(names_, 6, mkChar("converged"));
  SET_STRING_ELT(names_, 7, mkChar("scale_value"));
  SET_STRING_ELT(names_, 8, mkChar("workspace"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(14);
  return out_;
}

static SEXP native_operator_svd_certificate_cached_av(void* impl,
                                                      EigencoreApplyFn apply,
                                                      int m,
                                                      int n,
                                                      double norm_A,
                                                      SEXP d_,
                                                      SEXP u_,
                                                      SEXP v_,
                                                      SEXP av_,
                                                      SEXP tol_) {
  SEXP dimU = getAttrib(u_, R_DimSymbol);
  SEXP dimV = getAttrib(v_, R_DimSymbol);
  SEXP dimAV = getAttrib(av_, R_DimSymbol);
  if (dimU == R_NilValue || dimV == R_NilValue || dimAV == R_NilValue) {
    error("u, v, and Av must be matrices");
  }
  const int rowsU = INTEGER(dimU)[0];
  const int k = INTEGER(dimU)[1];
  const int rowsV = INTEGER(dimV)[0];
  const int colsV = INTEGER(dimV)[1];
  const int rowsAV = INTEGER(dimAV)[0];
  const int colsAV = INTEGER(dimAV)[1];
  if (rowsU != m || rowsV != n || colsV != k ||
      rowsAV != m || colsAV != k || LENGTH(d_) != k) {
    error("non-conformable cached-Av native operator SVD certificate inputs");
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double eps = DBL_EPSILON;
  const double tol = asReal(tol_);
  const double scale_value = fmax(norm_A, eps);
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};

  SEXP right_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP left_ = PROTECT(allocVector(REALSXP, k));
  SEXP right_ = PROTECT(allocVector(REALSXP, k));
  SEXP combined_ = PROTECT(allocVector(REALSXP, k));
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));

  const int status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, k,
                           REAL(u_), m, 1.0, 0.0,
                           REAL(right_matrix_), n, &workspace);
  if (status != 0) {
    error("cached-Av native operator SVD certificate adjoint apply failed with status=%d", status);
  }

  for (int col = 0; col < k; ++col) {
    const double sigma = REAL(d_)[col];
    const int left_offset = col * m;
    const int right_offset = col * n;
    double left_sum = 0.0;
    for (int row = 0; row < m; ++row) {
      const double residual = REAL(av_)[left_offset + row] -
        sigma * REAL(u_)[left_offset + row];
      left_sum += residual * residual;
    }
    for (int row = 0; row < n; ++row) {
      REAL(right_matrix_)[right_offset + row] -= sigma * REAL(v_)[right_offset + row];
    }
    const double left = sqrt(left_sum);
    const double right = column_norm(REAL(right_matrix_), n, col);
    const double combined = sqrt(left * left + right * right);
    const double backward = combined / scale_value;
    REAL(left_)[col] = left;
    REAL(right_)[col] = right;
    REAL(combined_)[col] = combined;
    REAL(scale_)[col] = scale_value;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  SEXP gram_u_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP gram_v_ = PROTECT(allocMatrix(REALSXP, k, k));
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &m,
                  &one, REAL(u_), &m, REAL(u_), &m,
                  &zero, REAL(gram_u_), &k FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, REAL(v_), &n, REAL(v_), &n,
                  &zero, REAL(gram_v_), &k FCONE FCONE);
  REAL(orth_)[0] = max_orthogonality_loss(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss(REAL(gram_v_), k);
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 9));
  SET_VECTOR_ELT(out_, 0, left_);
  SET_VECTOR_ELT(out_, 1, right_);
  SET_VECTOR_ELT(out_, 2, combined_);
  SET_VECTOR_ELT(out_, 3, backward_);
  SET_VECTOR_ELT(out_, 4, orth_);
  SET_VECTOR_ELT(out_, 5, scale_);
  SET_VECTOR_ELT(out_, 6, converged_);
  SET_VECTOR_ELT(out_, 7, ScalarReal(scale_value));
  SET_VECTOR_ELT(out_, 8, workspace_counters(&workspace));
  SEXP names_ = PROTECT(allocVector(STRSXP, 9));
  SET_STRING_ELT(names_, 0, mkChar("left"));
  SET_STRING_ELT(names_, 1, mkChar("right"));
  SET_STRING_ELT(names_, 2, mkChar("combined"));
  SET_STRING_ELT(names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 5, mkChar("scale"));
  SET_STRING_ELT(names_, 6, mkChar("converged"));
  SET_STRING_ELT(names_, 7, mkChar("scale_value"));
  SET_STRING_ELT(names_, 8, mkChar("workspace"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(13);
  return out_;
}

extern "C" SEXP eigencore_csc_eigen_certificate(SEXP i_, SEXP p_, SEXP x_,
                                                SEXP dim_, SEXP values_,
                                                SEXP vectors_, SEXP norm_A_,
                                                SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(values_) || !isReal(vectors_)) {
    error("invalid CSC eigen certificate inputs");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n) {
    error("CSC eigen certificate requires a square operator");
  }
  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_operator_eigen_certificate(&impl, eigencore_csc_apply, n,
                                           asReal(norm_A_), values_, vectors_, tol_);
}

extern "C" SEXP eigencore_diagonal_eigen_certificate(SEXP x_, SEXP dim_,
                                                     SEXP unit_, SEXP values_,
                                                     SEXP vectors_, SEXP norm_A_,
                                                     SEXP tol_) {
  if (!isReal(x_) || !isInteger(dim_) || !isLogical(unit_) ||
      !isReal(values_) || !isReal(vectors_)) {
    error("invalid diagonal eigen certificate inputs");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n) {
    error("diagonal eigen certificate requires a square operator");
  }
  DiagonalOperator impl = {n, REAL(x_), LOGICAL(unit_)[0] == TRUE};
  return native_operator_eigen_certificate(&impl, eigencore_diagonal_apply, n,
                                           asReal(norm_A_), values_, vectors_, tol_);
}

extern "C" SEXP eigencore_csc_svd_certificate(SEXP i_, SEXP p_, SEXP x_,
                                              SEXP dim_, SEXP d_,
                                              SEXP u_, SEXP v_,
                                              SEXP norm_A_, SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_)) {
    error("invalid CSC SVD certificate inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_operator_svd_certificate(&impl, eigencore_csc_apply, m, n,
                                         asReal(norm_A_), d_, u_, v_, tol_);
}

extern "C" SEXP eigencore_csc_svd_certificate_cached_av(SEXP i_, SEXP p_, SEXP x_,
                                                        SEXP dim_, SEXP d_,
                                                        SEXP u_, SEXP v_,
                                                        SEXP av_, SEXP norm_A_,
                                                        SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_) || !isReal(av_)) {
    error("invalid cached-Av CSC SVD certificate inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_operator_svd_certificate_cached_av(
    &impl, eigencore_csc_apply, m, n,
    asReal(norm_A_), d_, u_, v_, av_, tol_
  );
}

extern "C" SEXP eigencore_diagonal_svd_certificate(SEXP x_, SEXP dim_,
                                                   SEXP unit_, SEXP d_,
                                                   SEXP u_, SEXP v_,
                                                   SEXP norm_A_, SEXP tol_) {
  if (!isReal(x_) || !isInteger(dim_) || !isLogical(unit_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_)) {
    error("invalid diagonal SVD certificate inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  DiagonalOperator impl = {m, REAL(x_), LOGICAL(unit_)[0] == TRUE};
  return native_operator_svd_certificate(&impl, eigencore_diagonal_apply, m, n,
                                         asReal(norm_A_), d_, u_, v_, tol_);
}

extern "C" SEXP eigencore_diagonal_svd_certificate_cached_av(SEXP x_, SEXP dim_,
                                                             SEXP unit_, SEXP d_,
                                                             SEXP u_, SEXP v_,
                                                             SEXP av_,
                                                             SEXP norm_A_,
                                                             SEXP tol_) {
  if (!isReal(x_) || !isInteger(dim_) || !isLogical(unit_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_) || !isReal(av_)) {
    error("invalid cached-Av diagonal SVD certificate inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  DiagonalOperator impl = {m, REAL(x_), LOGICAL(unit_)[0] == TRUE};
  return native_operator_svd_certificate_cached_av(
    &impl, eigencore_diagonal_apply, m, n,
    asReal(norm_A_), d_, u_, v_, av_, tol_
  );
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

static SEXP block_golub_kahan_ritz_pack(const double* V, int n,
                                        const double* AV, int m,
                                        int p, int rank,
                                        int target_kind) {
  if (rank < 1) {
    error("rank must be positive");
  }
  const int s = (m < p) ? m : p;
  if (rank > s) {
    rank = s;
  }

  std::vector<double> B(static_cast<size_t>(m) * p, 0.0);
  for (int col = 0; col < p; ++col) {
    std::memcpy(B.data() + static_cast<size_t>(col) * m,
                AV + static_cast<size_t>(col) * m,
                sizeof(double) * static_cast<size_t>(m));
  }
  std::vector<double> d_all(static_cast<size_t>(s), 0.0);
  std::vector<double> u_all(static_cast<size_t>(m) * s, 0.0);
  std::vector<double> vt_all(static_cast<size_t>(s) * p, 0.0);

  char jobu = 'S';
  char jobvt = 'S';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgesvd)(&jobu, &jobvt, &m, &p, B.data(), &m,
                   d_all.data(), u_all.data(), &m, vt_all.data(), &s,
                   &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgesvd workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  std::vector<double> work(static_cast<size_t>(lwork), 0.0);
  F77_CALL(dgesvd)(&jobu, &jobvt, &m, &p, B.data(), &m,
                   d_all.data(), u_all.data(), &m, vt_all.data(), &s,
                   work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgesvd failed with info=%d", info);
  }

  std::vector<int> selected(static_cast<size_t>(rank));
  const int count = selected_ritz_indices(
    d_all.data(), s, rank, target_kind, selected.data()
  );

  SEXP d_ = PROTECT(allocVector(REALSXP, count));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, count));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, count));
  SEXP avectors_ = PROTECT(allocMatrix(REALSXP, m, count));
  SEXP coeff_ = PROTECT(allocMatrix(REALSXP, p, count));

  for (int col = 0; col < count; ++col) {
    const int idx = selected[static_cast<size_t>(col)];
    REAL(d_)[col] = d_all[static_cast<size_t>(idx)];
    for (int row = 0; row < m; ++row) {
      REAL(u_)[row + static_cast<size_t>(col) * m] =
        u_all[row + static_cast<size_t>(idx) * m];
    }
    for (int row = 0; row < p; ++row) {
      REAL(coeff_)[row + static_cast<size_t>(col) * p] =
        vt_all[idx + static_cast<size_t>(row) * s];
    }
  }

  if (count > 0) {
    const char notrans = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    F77_CALL(dgemm)(&notrans, &notrans, &n, &count, &p,
                    &one, V, &n, REAL(coeff_), &p,
                    &zero, REAL(v_), &n FCONE FCONE);
    F77_CALL(dgemm)(&notrans, &notrans, &m, &count, &p,
                    &one, AV, &m, REAL(coeff_), &p,
                    &zero, REAL(avectors_), &m FCONE FCONE);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, avectors_);
  SET_VECTOR_ELT(out_, 4, coeff_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("Avectors"));
  SET_STRING_ELT(names_, 4, mkChar("coefficients"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(7);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_ritz(SEXP V_, SEXP AV_,
                                                 SEXP rank_, SEXP target_kind_,
                                                 SEXP active_p_) {
  if (!isReal(V_) || !isReal(AV_)) {
    error("V and AV must be double matrices");
  }
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  SEXP dimAV = getAttrib(AV_, R_DimSymbol);
  if (dimV == R_NilValue || dimAV == R_NilValue) {
    error("V and AV must be matrices");
  }
  const int n = INTEGER(dimV)[0];
  const int stored_p = INTEGER(dimV)[1];
  const int m = INTEGER(dimAV)[0];
  const int av_p = INTEGER(dimAV)[1];
  int p = asInteger(active_p_);
  if (p < 1 || p > stored_p || p > av_p) {
    error("active block Golub-Kahan Ritz columns must be between 1 and ncol(V/AV)");
  }
  return block_golub_kahan_ritz_pack(
    REAL(V_), n, REAL(AV_), m, p,
    asInteger(rank_), asInteger(target_kind_)
  );
}

static SEXP block_golub_kahan_fit_pack(int n,
                                       int m,
                                       double* V,
                                       double* AV,
                                       int active_v,
                                       int active_u,
                                       int iterations,
                                       int matvecs,
                                       int ortho_passes,
                                       int cached_start_used,
                                       int rank,
                                       int target_kind,
                                       double stage_native_iteration_seconds) {
  auto ritz_timer = native_timer_now();
  SEXP ritz_ = PROTECT(block_golub_kahan_ritz_pack(
    V, n, AV, m, active_v, rank, target_kind
  ));
  const double stage_ritz_seconds = native_timer_elapsed(ritz_timer);
  SEXP stage_ = PROTECT(allocVector(REALSXP, 2));
  REAL(stage_)[0] = stage_native_iteration_seconds;
  REAL(stage_)[1] = stage_ritz_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(stage_names_, 0, mkChar("native_iteration"));
  SET_STRING_ELT(stage_names_, 1, mkChar("ritz"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 12));
  for (int i = 0; i < 5; ++i) {
    SET_VECTOR_ELT(out_, i, VECTOR_ELT(ritz_, i));
  }
  SET_VECTOR_ELT(out_, 5, ScalarInteger(active_v));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(active_u));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(ortho_passes));
  SET_VECTOR_ELT(out_, 10, ScalarLogical(cached_start_used != 0));
  SET_VECTOR_ELT(out_, 11, stage_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 12));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("Avectors"));
  SET_STRING_ELT(names_, 4, mkChar("coefficients"));
  SET_STRING_ELT(names_, 5, mkChar("active_cols"));
  SET_STRING_ELT(names_, 6, mkChar("active_left_cols"));
  SET_STRING_ELT(names_, 7, mkChar("iterations"));
  SET_STRING_ELT(names_, 8, mkChar("matvecs"));
  SET_STRING_ELT(names_, 9, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 10, mkChar("cached_start_used"));
  SET_STRING_ELT(names_, 11, mkChar("stage_seconds"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(5);
  return out_;
}

static SEXP block_golub_kahan_retained_cycle_pack(int n,
                                                  int m,
                                                  void* impl,
                                                  EigencoreApplyFn apply,
                                                  double norm_A,
                                                  double tol,
                                                  double* V,
                                                  double* AV,
                                                  int active_v,
                                                  int active_u,
                                                  int iterations,
                                                  int matvecs,
                                                  int ortho_passes,
                                                  int cached_start_used,
                                                  int requested_rank,
                                                  int active_rank,
                                                  int target_kind,
                                                  double stage_native_iteration_seconds,
                                                  double stage_restart_seconds,
                                                  SEXP attempt_history_,
                                                  double native_workspace_bytes,
                                                  const double* locked_d,
                                                  const double* locked_u,
                                                  const double* locked_v,
                                                  const double* locked_av,
                                                  int locked_count) {
  (void)requested_rank;
  if (active_rank < 0) {
    active_rank = 0;
  }
  SEXP fit_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, V, AV, active_v, active_u, iterations, matvecs, ortho_passes,
    cached_start_used, active_rank, target_kind, stage_native_iteration_seconds
  ));
  SEXP d_ = VECTOR_ELT(fit_, 0);
  SEXP u_ = VECTOR_ELT(fit_, 1);
  SEXP v_ = VECTOR_ELT(fit_, 2);
  SEXP avectors_ = VECTOR_ELT(fit_, 3);
  SEXP coeff_ = VECTOR_ELT(fit_, 4);
  SEXP combined_d_ = d_;
  SEXP combined_u_ = u_;
  SEXP combined_v_ = v_;
  SEXP combined_av_ = avectors_;
  SEXP combined_coeff_ = coeff_;
  const int active_count = LENGTH(d_);
  const int total_count = locked_count + active_count;
  if (locked_count > 0) {
    combined_d_ = PROTECT(allocVector(REALSXP, total_count));
    combined_u_ = PROTECT(allocMatrix(REALSXP, m, total_count));
    combined_v_ = PROTECT(allocMatrix(REALSXP, n, total_count));
    combined_av_ = PROTECT(allocMatrix(REALSXP, m, total_count));
    combined_coeff_ = PROTECT(allocMatrix(REALSXP, active_v, total_count));
    std::memset(REAL(combined_coeff_), 0,
                sizeof(double) * static_cast<size_t>(active_v) * total_count);
    for (int col = 0; col < locked_count; ++col) {
      REAL(combined_d_)[col] = locked_d[col];
      std::memcpy(REAL(combined_u_) + static_cast<int64_t>(col) * m,
                  locked_u + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
      std::memcpy(REAL(combined_v_) + static_cast<int64_t>(col) * n,
                  locked_v + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
      std::memcpy(REAL(combined_av_) + static_cast<int64_t>(col) * m,
                  locked_av + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
    }
    for (int col = 0; col < active_count; ++col) {
      const int dst_col = locked_count + col;
      REAL(combined_d_)[dst_col] = REAL(d_)[col];
      std::memcpy(REAL(combined_u_) + static_cast<int64_t>(dst_col) * m,
                  REAL(u_) + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
      std::memcpy(REAL(combined_v_) + static_cast<int64_t>(dst_col) * n,
                  REAL(v_) + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
      std::memcpy(REAL(combined_av_) + static_cast<int64_t>(dst_col) * m,
                  REAL(avectors_) + static_cast<int64_t>(col) * m,
                  sizeof(double) * static_cast<size_t>(m));
      for (int row = 0; row < active_v; ++row) {
        REAL(combined_coeff_)[row + static_cast<int64_t>(dst_col) * active_v] =
          REAL(coeff_)[row + static_cast<int64_t>(col) * active_v];
      }
    }
  }
  SEXP tol_ = PROTECT(ScalarReal(tol));
  SEXP cert_diag_ = PROTECT(native_operator_svd_certificate_cached_av(
    impl, apply, m, n, norm_A,
    combined_d_, combined_u_, combined_v_, combined_av_, tol_
  ));
  SEXP old_stage_ = VECTOR_ELT(fit_, 11);
  SEXP stage_ = PROTECT(allocVector(REALSXP, 3));
  REAL(stage_)[0] = REAL(old_stage_)[0];
  REAL(stage_)[1] = REAL(old_stage_)[1];
  REAL(stage_)[2] = stage_restart_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(stage_names_, 0, mkChar("native_iteration"));
  SET_STRING_ELT(stage_names_, 1, mkChar("ritz"));
  SET_STRING_ELT(stage_names_, 2, mkChar("restart"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 17));
  SET_VECTOR_ELT(out_, 0, combined_d_);
  SET_VECTOR_ELT(out_, 1, combined_u_);
  SET_VECTOR_ELT(out_, 2, combined_v_);
  SET_VECTOR_ELT(out_, 3, combined_av_);
  SET_VECTOR_ELT(out_, 4, combined_coeff_);
  for (int i = 5; i < 11; ++i) {
    SET_VECTOR_ELT(out_, i, VECTOR_ELT(fit_, i));
  }
  SET_VECTOR_ELT(out_, 11, stage_);
  SET_VECTOR_ELT(out_, 12, attempt_history_);
  SET_VECTOR_ELT(out_, 13, ScalarLogical(1));
  SET_VECTOR_ELT(out_, 14, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 15, cert_diag_);
  SET_VECTOR_ELT(out_, 16, ScalarInteger(locked_count));
  SEXP names_ = PROTECT(allocVector(STRSXP, 17));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("Avectors"));
  SET_STRING_ELT(names_, 4, mkChar("coefficients"));
  SET_STRING_ELT(names_, 5, mkChar("active_cols"));
  SET_STRING_ELT(names_, 6, mkChar("active_left_cols"));
  SET_STRING_ELT(names_, 7, mkChar("iterations"));
  SET_STRING_ELT(names_, 8, mkChar("matvecs"));
  SET_STRING_ELT(names_, 9, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 10, mkChar("cached_start_used"));
  SET_STRING_ELT(names_, 11, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("attempt_history"));
  SET_STRING_ELT(names_, 13, mkChar("retained_restart_native"));
  SET_STRING_ELT(names_, 14, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 15, mkChar("certificate_diagnostics"));
  SET_STRING_ELT(names_, 16, mkChar("retained_locked_count"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(locked_count > 0 ? 12 : 7);
  return out_;
}

struct BlockGolubKahanFitArrays {
  double* V = nullptr;
  double* AV = nullptr;
  double* U = nullptr;
  bool transient = false;
};

static void block_golub_kahan_fit_arrays_free(BlockGolubKahanFitArrays* arrays) {
  if (!arrays->transient) {
    std::free(arrays->V);
    std::free(arrays->AV);
    std::free(arrays->U);
  }
  arrays->V = nullptr;
  arrays->AV = nullptr;
  arrays->U = nullptr;
  arrays->transient = false;
}

static int block_golub_kahan_fit_arrays_alloc(BlockGolubKahanFitArrays* arrays,
                                              int n,
                                              int m,
                                              int max_subspace) {
  const size_t nv = static_cast<size_t>(n) * static_cast<size_t>(max_subspace);
  const size_t mv = static_cast<size_t>(m) * static_cast<size_t>(max_subspace);
  arrays->V = reinterpret_cast<double*>(R_alloc(nv > 0 ? nv : 1, sizeof(double)));
  arrays->AV = reinterpret_cast<double*>(R_alloc(mv > 0 ? mv : 1, sizeof(double)));
  arrays->U = reinterpret_cast<double*>(R_alloc(mv > 0 ? mv : 1, sizeof(double)));
  arrays->transient = true;
  return 0;
}

static int retained_subspace_sequence(int n,
                                      int rank,
                                      int block_size,
                                      int initial_max_subspace,
                                      int max_attempts,
                                      std::vector<int>* subspaces) {
  if (max_attempts < 1 || initial_max_subspace < rank || initial_max_subspace > n) {
    return -1;
  }
  subspaces->clear();
  subspaces->reserve(static_cast<size_t>(max_attempts));
  int current = initial_max_subspace;
  subspaces->push_back(current);
  for (int attempt = 1; attempt < max_attempts && current < n; ++attempt) {
    const int additive = current + ((2 * rank > 4 * block_size) ? 2 * rank : 4 * block_size);
    const int additive2 = current + 10;
    int next = additive > additive2 ? additive : additive2;
    const int scaled = static_cast<int>(ceil(1.5 * static_cast<double>(current)));
    if (scaled > next) {
      next = scaled;
    }
    if (next > n) {
      next = n;
    }
    if (next <= current) {
      break;
    }
    subspaces->push_back(next);
    current = next;
  }
  return 0;
}

static SEXP retained_attempt_history_pack(const std::vector<int>& subspaces,
                                          const std::vector<int>& active_cols,
                                          const std::vector<int>& start_cols,
                                          const std::vector<int>& iterations,
                                          const std::vector<int>& matvecs,
                                          const std::vector<int>& ortho_passes,
                                          const std::vector<int>& cached_start_used,
                                          const std::vector<int>& converged_count,
                                          const std::vector<int>& leading_converged_count,
                                          const std::vector<int>& certificate_passed,
                                          const std::vector<double>& max_backward_error,
                                          const std::vector<double>& max_residual) {
  const int rows = static_cast<int>(subspaces.size());
  SEXP out_ = PROTECT(allocVector(VECSXP, 14));
  SEXP attempt_ = PROTECT(allocVector(INTSXP, rows));
  SEXP subspace_ = PROTECT(allocVector(INTSXP, rows));
  SEXP active_ = PROTECT(allocVector(INTSXP, rows));
  SEXP start_ = PROTECT(allocVector(INTSXP, rows));
  SEXP cached_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP warm_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP iter_ = PROTECT(allocVector(INTSXP, rows));
  SEXP matvec_ = PROTECT(allocVector(INTSXP, rows));
  SEXP ortho_ = PROTECT(allocVector(INTSXP, rows));
  SEXP conv_count_ = PROTECT(allocVector(INTSXP, rows));
  SEXP leading_count_ = PROTECT(allocVector(INTSXP, rows));
  SEXP passed_ = PROTECT(allocVector(LGLSXP, rows));
  SEXP backward_ = PROTECT(allocVector(REALSXP, rows));
  SEXP residual_ = PROTECT(allocVector(REALSXP, rows));
  for (int row = 0; row < rows; ++row) {
    INTEGER(attempt_)[row] = row + 1;
    INTEGER(subspace_)[row] = subspaces[static_cast<size_t>(row)];
    INTEGER(active_)[row] = active_cols[static_cast<size_t>(row)];
    INTEGER(start_)[row] = start_cols[static_cast<size_t>(row)];
    LOGICAL(cached_)[row] = cached_start_used[static_cast<size_t>(row)] ? TRUE : FALSE;
    LOGICAL(warm_)[row] = row > 0 ? TRUE : FALSE;
    INTEGER(iter_)[row] = iterations[static_cast<size_t>(row)];
    INTEGER(matvec_)[row] = matvecs[static_cast<size_t>(row)];
    INTEGER(ortho_)[row] = ortho_passes[static_cast<size_t>(row)];
    INTEGER(conv_count_)[row] = converged_count[static_cast<size_t>(row)];
    INTEGER(leading_count_)[row] = leading_converged_count[static_cast<size_t>(row)];
    LOGICAL(passed_)[row] = certificate_passed[static_cast<size_t>(row)] ? TRUE : FALSE;
    REAL(backward_)[row] = max_backward_error[static_cast<size_t>(row)];
    REAL(residual_)[row] = max_residual[static_cast<size_t>(row)];
  }
  SET_VECTOR_ELT(out_, 0, attempt_);
  SET_VECTOR_ELT(out_, 1, subspace_);
  SET_VECTOR_ELT(out_, 2, active_);
  SET_VECTOR_ELT(out_, 3, start_);
  SET_VECTOR_ELT(out_, 4, cached_);
  SET_VECTOR_ELT(out_, 5, warm_);
  SET_VECTOR_ELT(out_, 6, iter_);
  SET_VECTOR_ELT(out_, 7, matvec_);
  SET_VECTOR_ELT(out_, 8, ortho_);
  SET_VECTOR_ELT(out_, 9, conv_count_);
  SET_VECTOR_ELT(out_, 10, leading_count_);
  SET_VECTOR_ELT(out_, 11, passed_);
  SET_VECTOR_ELT(out_, 12, backward_);
  SET_VECTOR_ELT(out_, 13, residual_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 14));
  SET_STRING_ELT(names_, 0, mkChar("attempt"));
  SET_STRING_ELT(names_, 1, mkChar("max_subspace"));
  SET_STRING_ELT(names_, 2, mkChar("active_cols"));
  SET_STRING_ELT(names_, 3, mkChar("start_cols"));
  SET_STRING_ELT(names_, 4, mkChar("cached_start_used"));
  SET_STRING_ELT(names_, 5, mkChar("warm_started"));
  SET_STRING_ELT(names_, 6, mkChar("iterations"));
  SET_STRING_ELT(names_, 7, mkChar("matvecs"));
  SET_STRING_ELT(names_, 8, mkChar("ortho_passes"));
  SET_STRING_ELT(names_, 9, mkChar("converged_count"));
  SET_STRING_ELT(names_, 10, mkChar("leading_converged_count"));
  SET_STRING_ELT(names_, 11, mkChar("certificate_passed"));
  SET_STRING_ELT(names_, 12, mkChar("max_backward_error"));
  SET_STRING_ELT(names_, 13, mkChar("max_residual"));
  setAttrib(out_, R_NamesSymbol, names_);

  SEXP row_names_ = PROTECT(allocVector(INTSXP, rows));
  for (int row = 0; row < rows; ++row) {
    INTEGER(row_names_)[row] = row + 1;
  }
  setAttrib(out_, R_RowNamesSymbol, row_names_);
  SEXP class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(class_, 0, mkChar("data.frame"));
  setAttrib(out_, R_ClassSymbol, class_);
  UNPROTECT(18);
  return out_;
}

static int retained_cached_av_certificate_passed(void* impl,
                                                 EigencoreApplyFn apply,
                                                 int m,
                                                 int n,
                                                 double norm_A,
                                                 const double* d,
                                                 const double* u,
                                                 const double* v,
                                                 const double* av,
                                                 int k,
                                                 double tol,
                                                 double* max_backward_error,
                                                 double* max_residual,
                                                 int* converged_count,
                                                 int* leading_converged_count) {
  const double eps = DBL_EPSILON;
  const double scale_value = fmax(norm_A, eps);
  std::vector<double> right(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};
  const int status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, k,
                           u, m, 1.0, 0.0, right.data(), n, &workspace);
  if (status != 0) {
    return status < 0 ? status : -status;
  }

  int passed = 1;
  *max_backward_error = 0.0;
  *max_residual = 0.0;
  *converged_count = 0;
  *leading_converged_count = 0;
  bool still_leading = true;
  for (int col = 0; col < k; ++col) {
    const double sigma = d[col];
    double left_sum = 0.0;
    double right_sum = 0.0;
    const int64_t left_offset = static_cast<int64_t>(col) * m;
    const int64_t right_offset = static_cast<int64_t>(col) * n;
    for (int row = 0; row < m; ++row) {
      const double residual = av[left_offset + row] -
        sigma * u[left_offset + row];
      left_sum += residual * residual;
    }
    for (int row = 0; row < n; ++row) {
      const double residual = right[right_offset + row] -
        sigma * v[right_offset + row];
      right_sum += residual * residual;
    }
    const double combined = sqrt(left_sum + right_sum);
    const double backward = combined / scale_value;
    if (backward > *max_backward_error || col == 0) {
      *max_backward_error = backward;
    }
    if (combined > *max_residual || col == 0) {
      *max_residual = combined;
    }
    const bool col_converged = R_FINITE(backward) && backward <= tol;
    if (col_converged) {
      ++(*converged_count);
      if (still_leading) {
        ++(*leading_converged_count);
      }
    } else {
      still_leading = false;
      passed = 0;
    }
  }

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  std::vector<double> gram_u(static_cast<size_t>(k) * static_cast<size_t>(k), 0.0);
  std::vector<double> gram_v(static_cast<size_t>(k) * static_cast<size_t>(k), 0.0);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &m,
                  &one, const_cast<double*>(u), &m, const_cast<double*>(u), &m,
                  &zero, gram_u.data(), &k FCONE FCONE);
  F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                  &one, const_cast<double*>(v), &n, const_cast<double*>(v), &n,
                  &zero, gram_v.data(), &k FCONE FCONE);
  const double orth = fmax(
    max_orthogonality_loss(gram_u.data(), k),
    max_orthogonality_loss(gram_v.data(), k)
  );
  const double orth_tol = (tol > sqrt(DBL_EPSILON)) ? tol : sqrt(DBL_EPSILON);
  if (orth > orth_tol) {
    passed = 0;
    *converged_count = 0;
    *leading_converged_count = 0;
  }
  return passed;
}

static int normalize_retained_cached_prefix(double* V,
                                            int n,
                                            double* AV,
                                            int m,
                                            int cols) {
  if (cols <= 0) {
    return 0;
  }
  std::vector<double> gram(static_cast<size_t>(cols) * static_cast<size_t>(cols), 0.0);
  const char trans = 'T';
  const char notrans = 'N';
  const char right = 'R';
  const char uplo = 'U';
  const char diag = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  int info = 0;

  F77_CALL(dgemm)(&trans, &notrans, &cols, &cols, &n,
                  &one, V, &n, V, &n,
                  &zero, gram.data(), &cols FCONE FCONE);
  symmetrize_packed_square(gram.data(), cols);
  F77_CALL(dpotrf)(&uplo, &cols, gram.data(), &cols, &info FCONE);
  if (info != 0) {
    return -1;
  }
  for (int col = 0; col < cols; ++col) {
    if (gram[col + static_cast<int64_t>(col) * cols] <= 100.0 * DBL_EPSILON) {
      return -2;
    }
  }

  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &n, &cols, &one,
                  gram.data(), &cols, V, &n FCONE FCONE FCONE FCONE);
  F77_CALL(dtrsm)(&right, &uplo, &notrans, &diag, &m, &cols, &one,
                  gram.data(), &cols, AV, &m FCONE FCONE FCONE FCONE);
  return 0;
}

template <typename ConfigureOperator>
static SEXP block_golub_kahan_retained_cycle_impl(ConfigureOperator configure_operator,
                                                  int m,
                                                  int n,
                                                  int initial_max_subspace,
                                                  int block_size,
                                                  const double* initial_start,
                                                  const double* random_tails,
                                                  int random_tail_cols,
                                                  int max_attempts,
                                                  int rank,
                                                  int target_kind,
                                                  double norm_A,
                                                  double tol,
                                                  int use_retained_av_cache,
                                                  int use_deflation) {
  if (block_size < 1 || rank < 1 || initial_max_subspace < rank ||
      initial_max_subspace > n || max_attempts < 1) {
    error("invalid retained block Golub-Kahan cycle controls");
  }
  if (random_tail_cols < block_size * (max_attempts - 1)) {
    error("random_tails must have at least block * (max_attempts - 1) columns");
  }

  std::vector<int> subspaces;
  if (retained_subspace_sequence(
        n, rank, block_size, initial_max_subspace, max_attempts, &subspaces
      ) != 0 || subspaces.empty()) {
    error("failed to construct retained block Golub-Kahan subspace sequence");
  }
  const int final_max_subspace = subspaces.back();
  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, final_max_subspace) != 0) {
    error("failed to allocate native retained block Golub-Kahan workspace");
  }
  BlockGolubKahanBasisScratch basis_scratch;
  const int retained_start_capacity = rank + block_size;
  if (block_golub_kahan_basis_scratch_alloc(
        &basis_scratch, m, n, final_max_subspace, retained_start_capacity
      ) != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("failed to allocate native retained block Golub-Kahan basis scratch");
  }
  const double native_workspace_bytes = static_cast<double>(
    (static_cast<size_t>(n) * static_cast<size_t>(final_max_subspace) +
     2 * static_cast<size_t>(m) * static_cast<size_t>(final_max_subspace)) *
      sizeof(double) + basis_scratch.bytes
  );

  std::vector<double> current_start(static_cast<size_t>(n) *
                                    static_cast<size_t>(rank + block_size), 0.0);
  std::vector<double> current_start_av(static_cast<size_t>(m) *
                                       static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_d(static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_u(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_v(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  std::vector<double> locked_av(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  std::memcpy(
    current_start.data(), initial_start,
    sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(block_size)
  );
  int current_start_cols = block_size;
  int current_start_av_cols = 0;

  void* impl = nullptr;
  EigencoreApplyFn apply = nullptr;
  configure_operator(&impl, &apply);

  std::vector<int> history_subspaces;
  std::vector<int> history_active_cols;
  std::vector<int> history_start_cols;
  std::vector<int> history_iterations;
  std::vector<int> history_matvecs;
  std::vector<int> history_ortho_passes;
  std::vector<int> history_cached_start_used;
  std::vector<int> history_converged_count;
  std::vector<int> history_leading_converged_count;
  std::vector<int> history_certificate_passed;
  std::vector<double> history_max_backward_error;
  std::vector<double> history_max_residual;

  int final_active_v = 0;
  int final_active_u = 0;
  int final_iterations = 0;
  int final_matvecs = 0;
  int final_ortho_passes = 0;
  int final_cached_start_used = 0;
  int locked_count = 0;
  double total_stage_native_iteration_seconds = 0.0;
  double total_stage_restart_seconds = 0.0;

  for (size_t attempt = 0; attempt < subspaces.size(); ++attempt) {
    const int active_rank = rank - locked_count;
    if (active_rank <= 0) {
      break;
    }
    const int max_subspace = subspaces[attempt];
    int active_max_subspace = max_subspace;
    if (use_deflation && locked_count > 0) {
      active_max_subspace -= locked_count;
      if (active_max_subspace < active_rank) {
        active_max_subspace = active_rank;
      }
      if (active_max_subspace > n - locked_count) {
        active_max_subspace = n - locked_count;
      }
    }
    int active_v = 0;
    int active_u = 0;
    int iterations = 0;
    int matvecs = 0;
    int ortho_passes = 0;
    int cached_start_used = 0;

    auto stage_timer = native_timer_now();
    const int status = native_block_golub_kahan_basis_run_with_scratch(
      impl, apply, m, n, active_max_subspace, current_start_cols,
      current_start.data(),
      (use_retained_av_cache && current_start_av_cols > 0) ? current_start_av.data() : nullptr,
      use_retained_av_cache ? current_start_av_cols : 0,
      (use_deflation && locked_count > 0) ? locked_v.data() : nullptr,
      use_deflation ? locked_count : 0,
      (use_deflation && locked_count > 0) ? locked_u.data() : nullptr,
      use_deflation ? locked_count : 0,
      arrays.V, arrays.AV, arrays.U, &basis_scratch,
      &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
      &cached_start_used
    );
    total_stage_native_iteration_seconds += native_timer_elapsed(stage_timer);
    if (status != 0) {
      block_golub_kahan_basis_scratch_free(&basis_scratch);
      block_golub_kahan_fit_arrays_free(&arrays);
      error("native retained block Golub-Kahan cycle failed with status=%d", status);
    }

    history_subspaces.push_back(max_subspace);
    history_active_cols.push_back(active_v);
    history_start_cols.push_back(current_start_cols);
    history_iterations.push_back(iterations);
    history_matvecs.push_back(matvecs);
    history_ortho_passes.push_back(ortho_passes);
    history_cached_start_used.push_back(cached_start_used);

    final_active_v = active_v;
    final_active_u = active_u;
    final_iterations += iterations;
    final_matvecs += matvecs;
    final_ortho_passes += ortho_passes;
    final_cached_start_used = final_cached_start_used || cached_start_used;

    SEXP ritz_ = PROTECT(block_golub_kahan_ritz_pack(
      arrays.V, n, arrays.AV, m, active_v, active_rank, target_kind
    ));
    double max_backward_error = R_PosInf;
    double max_residual = R_PosInf;
    int converged_count = 0;
    int leading_converged_count = 0;
    const int certificate_passed = retained_cached_av_certificate_passed(
      impl, apply, m, n, norm_A,
      REAL(VECTOR_ELT(ritz_, 0)),
      REAL(VECTOR_ELT(ritz_, 1)),
      REAL(VECTOR_ELT(ritz_, 2)),
      REAL(VECTOR_ELT(ritz_, 3)),
      LENGTH(VECTOR_ELT(ritz_, 0)),
      tol,
      &max_backward_error,
      &max_residual,
      &converged_count,
      &leading_converged_count
    );
    if (certificate_passed < 0) {
      UNPROTECT(1);
      block_golub_kahan_basis_scratch_free(&basis_scratch);
      block_golub_kahan_fit_arrays_free(&arrays);
      error("native retained block Golub-Kahan attempt certificate failed with status=%d",
            certificate_passed);
    }
    history_converged_count.push_back(locked_count + converged_count);
    history_leading_converged_count.push_back(locked_count + leading_converged_count);
    history_certificate_passed.push_back(certificate_passed);
    history_max_backward_error.push_back(max_backward_error);
    history_max_residual.push_back(max_residual);

    if (certificate_passed) {
      UNPROTECT(1);
      break;
    }

    if (attempt + 1 >= subspaces.size()) {
      UNPROTECT(1);
      break;
    }

    auto restart_timer = native_timer_now();
    SEXP v_ = VECTOR_ELT(ritz_, 2);
    SEXP av_ = VECTOR_ELT(ritz_, 3);
    SEXP u_ = VECTOR_ELT(ritz_, 1);
    SEXP d_ = VECTOR_ELT(ritz_, 0);
    const int ritz_cols = INTEGER(getAttrib(v_, R_DimSymbol))[1];
    int newly_locked = 0;
    if (use_deflation && leading_converged_count > 0) {
      newly_locked = leading_converged_count;
      if (newly_locked > rank - locked_count) {
        newly_locked = rank - locked_count;
      }
      for (int col = 0; col < newly_locked; ++col) {
        const int dst_col = locked_count + col;
        locked_d[static_cast<size_t>(dst_col)] = REAL(d_)[col];
        std::memcpy(locked_u.data() + static_cast<int64_t>(dst_col) * m,
                    REAL(u_) + static_cast<int64_t>(col) * m,
                    sizeof(double) * static_cast<size_t>(m));
        std::memcpy(locked_v.data() + static_cast<int64_t>(dst_col) * n,
                    REAL(v_) + static_cast<int64_t>(col) * n,
                    sizeof(double) * static_cast<size_t>(n));
        std::memcpy(locked_av.data() + static_cast<int64_t>(dst_col) * m,
                    REAL(av_) + static_cast<int64_t>(col) * m,
                    sizeof(double) * static_cast<size_t>(m));
      }
      locked_count += newly_locked;
    }
    const int keep_cols = ritz_cols - newly_locked;
    current_start.assign(current_start.size(), 0.0);
    current_start_av.assign(current_start_av.size(), 0.0);
    for (int col = 0; col < keep_cols; ++col) {
      std::memcpy(
        current_start.data() + static_cast<int64_t>(col) * n,
        REAL(v_) + static_cast<int64_t>(newly_locked + col) * n,
        sizeof(double) * static_cast<size_t>(n)
      );
      if (use_retained_av_cache) {
        std::memcpy(
          current_start_av.data() + static_cast<int64_t>(col) * m,
          REAL(av_) + static_cast<int64_t>(newly_locked + col) * m,
          sizeof(double) * static_cast<size_t>(m)
        );
      }
    }
    if (use_retained_av_cache) {
      const int normalize_status = normalize_retained_cached_prefix(
        current_start.data(), n, current_start_av.data(), m, keep_cols
      );
      if (normalize_status != 0) {
        UNPROTECT(1);
        block_golub_kahan_basis_scratch_free(&basis_scratch);
        block_golub_kahan_fit_arrays_free(&arrays);
        error("native retained block Golub-Kahan cached prefix normalization failed with status=%d",
              normalize_status);
      }
    }
    const double* tail = random_tails +
      static_cast<int64_t>(attempt) * static_cast<int64_t>(n) *
        static_cast<int64_t>(block_size);
    for (int col = 0; col < block_size; ++col) {
      std::memcpy(
        current_start.data() + static_cast<int64_t>(keep_cols + col) * n,
        tail + static_cast<int64_t>(col) * n,
        sizeof(double) * static_cast<size_t>(n)
      );
    }
    current_start_cols = keep_cols + block_size;
    current_start_av_cols = use_retained_av_cache ? keep_cols : 0;
    UNPROTECT(1);
    total_stage_restart_seconds += native_timer_elapsed(restart_timer);
  }

  SEXP history_ = PROTECT(retained_attempt_history_pack(
    history_subspaces, history_active_cols, history_start_cols,
    history_iterations, history_matvecs, history_ortho_passes,
    history_cached_start_used,
    history_converged_count, history_leading_converged_count,
    history_certificate_passed,
    history_max_backward_error, history_max_residual
  ));
  const int final_active_rank = rank - locked_count;
  SEXP out_ = PROTECT(block_golub_kahan_retained_cycle_pack(
    n, m, impl, apply, norm_A, tol,
    arrays.V, arrays.AV, final_active_v, final_active_u,
    final_iterations, final_matvecs, final_ortho_passes,
    final_cached_start_used, rank, final_active_rank, target_kind,
    total_stage_native_iteration_seconds, total_stage_restart_seconds,
    history_, native_workspace_bytes,
    locked_d.data(), locked_u.data(), locked_v.data(), locked_av.data(),
    locked_count
  ));
  block_golub_kahan_basis_scratch_free(&basis_scratch);
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(2);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_fit(SEXP A_,
                                                      SEXP max_subspace_,
                                                      SEXP start_,
                                                      SEXP rank_,
                                                      SEXP target_kind_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue) {
    error("A and start must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native dense block Golub-Kahan fit workspace");
  }
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native dense block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_fit_cached(SEXP A_,
                                                             SEXP max_subspace_,
                                                             SEXP start_,
                                                             SEXP rank_,
                                                             SEXP target_kind_,
                                                             SEXP start_av_) {
  if (!isReal(A_) || !isReal(start_) || !isReal(start_av_)) {
    error("A, start, and start_av must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || dimAV == R_NilValue) {
    error("A, start, and start_av must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  const int start_av_cols = INTEGER(dimAV)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native dense cached block Golub-Kahan fit workspace");
  }
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_dense_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native dense cached block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_dense_retained_cycle(SEXP A_,
                                                                 SEXP initial_max_subspace_,
                                                                 SEXP initial_start_,
                                                                 SEXP random_tails_,
                                                                 SEXP max_attempts_,
                                                                 SEXP rank_,
                                                                 SEXP target_kind_,
                                                                 SEXP norm_A_,
                                                                 SEXP tol_,
                                                                 SEXP use_retained_av_cache_,
                                                                 SEXP use_deflation_) {
  if (!isReal(A_) || !isReal(initial_start_) || !isReal(random_tails_)) {
    error("A, initial_start, and random_tails must be double matrices");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimS = getAttrib(initial_start_, R_DimSymbol);
  SEXP dimT = getAttrib(random_tails_, R_DimSymbol);
  if (dimA == R_NilValue || dimS == R_NilValue || dimT == R_NilValue) {
    error("A, initial_start, and random_tails must be matrices");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int block_size = INTEGER(dimS)[1];
  if (INTEGER(dimS)[0] != n || block_size < 1 || INTEGER(dimT)[0] != n) {
    error("non-conformable dense retained block Golub-Kahan inputs");
  }
  DenseColumnMajorOperator impl_holder = {m, n, REAL(A_)};
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_dense_apply;
  };
  return block_golub_kahan_retained_cycle_impl(
    configure, m, n, asInteger(initial_max_subspace_), block_size,
    REAL(initial_start_), REAL(random_tails_), INTEGER(dimT)[1],
    asInteger(max_attempts_), asInteger(rank_), asInteger(target_kind_),
    asReal(norm_A_), asReal(tol_),
    asLogical(use_retained_av_cache_) == TRUE,
    asLogical(use_deflation_) == TRUE
  );
}

extern "C" SEXP eigencore_block_golub_kahan_csc_fit(SEXP i_, SEXP p_,
                                                    SEXP x_, SEXP dim_,
                                                    SEXP max_subspace_,
                                                    SEXP start_,
                                                    SEXP rank_,
                                                    SEXP target_kind_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_)) {
    error("invalid CSC block Golub-Kahan fit inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  if (dimS == R_NilValue || LENGTH(dim_) != 2) {
    error("start must be a matrix and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native CSC block Golub-Kahan fit workspace");
  }
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    nullptr, 0, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native CSC block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_csc_fit_cached(SEXP i_, SEXP p_,
                                                           SEXP x_, SEXP dim_,
                                                           SEXP max_subspace_,
                                                           SEXP start_,
                                                           SEXP rank_,
                                                           SEXP target_kind_,
                                                           SEXP start_av_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(start_) || !isReal(start_av_)) {
    error("invalid cached CSC block Golub-Kahan fit inputs");
  }
  SEXP dimS = getAttrib(start_, R_DimSymbol);
  SEXP dimAV = getAttrib(start_av_, R_DimSymbol);
  if (dimS == R_NilValue || dimAV == R_NilValue || LENGTH(dim_) != 2) {
    error("start and start_av must be matrices and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  const int start_av_cols = INTEGER(dimAV)[1];
  int max_subspace = asInteger(max_subspace_);
  if (INTEGER(dimS)[0] != n || block_size < 1) {
    error("start must have nrow equal to ncol(A) and at least one column");
  }
  if (INTEGER(dimAV)[0] != m || start_av_cols < 1 || start_av_cols > block_size) {
    error("start_av must have nrow equal to nrow(A) and between 1 and ncol(start) columns");
  }
  if (max_subspace < 1 || max_subspace > n) {
    error("max_subspace must be between 1 and ncol(A)");
  }

  BlockGolubKahanFitArrays arrays;
  if (block_golub_kahan_fit_arrays_alloc(&arrays, n, m, max_subspace) != 0) {
    error("failed to allocate native CSC cached block Golub-Kahan fit workspace");
  }
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int active_v = 0;
  int active_u = 0;
  int iterations = 0;
  int matvecs = 0;
  int ortho_passes = 0;
  int cached_start_used = 0;
  auto stage_timer = native_timer_now();
  const int status = native_block_golub_kahan_basis_run(
    &impl, eigencore_csc_apply, m, n, max_subspace, block_size, REAL(start_),
    REAL(start_av_), start_av_cols, nullptr, 0, nullptr, 0,
    arrays.V, arrays.AV, arrays.U,
    &active_v, &active_u, &iterations, &matvecs, &ortho_passes,
    &cached_start_used
  );
  const double stage_native_iteration_seconds = native_timer_elapsed(stage_timer);
  if (status != 0) {
    block_golub_kahan_fit_arrays_free(&arrays);
    error("native CSC cached block Golub-Kahan fit failed with status=%d", status);
  }
  SEXP out_ = PROTECT(block_golub_kahan_fit_pack(
    n, m, arrays.V, arrays.AV, active_v, active_u, iterations, matvecs,
    ortho_passes, cached_start_used, asInteger(rank_), asInteger(target_kind_),
    stage_native_iteration_seconds
  ));
  block_golub_kahan_fit_arrays_free(&arrays);
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_block_golub_kahan_csc_retained_cycle(SEXP i_, SEXP p_,
                                                               SEXP x_, SEXP dim_,
                                                               SEXP initial_max_subspace_,
                                                               SEXP initial_start_,
                                                               SEXP random_tails_,
                                                               SEXP max_attempts_,
                                                               SEXP rank_,
                                                               SEXP target_kind_,
                                                               SEXP norm_A_,
                                                               SEXP tol_,
                                                               SEXP use_retained_av_cache_,
                                                               SEXP use_deflation_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) ||
      !isInteger(dim_) || !isReal(initial_start_) || !isReal(random_tails_)) {
    error("invalid CSC retained block Golub-Kahan inputs");
  }
  SEXP dimS = getAttrib(initial_start_, R_DimSymbol);
  SEXP dimT = getAttrib(random_tails_, R_DimSymbol);
  if (dimS == R_NilValue || dimT == R_NilValue || LENGTH(dim_) != 2) {
    error("initial_start and random_tails must be matrices and dim must have length 2");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int block_size = INTEGER(dimS)[1];
  if (INTEGER(dimS)[0] != n || block_size < 1 || INTEGER(dimT)[0] != n) {
    error("non-conformable CSC retained block Golub-Kahan inputs");
  }
  CSCOperator impl_holder = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_csc_apply;
  };
  return block_golub_kahan_retained_cycle_impl(
    configure, m, n, asInteger(initial_max_subspace_), block_size,
    REAL(initial_start_), REAL(random_tails_), INTEGER(dimT)[1],
    asInteger(max_attempts_), asInteger(rank_), asInteger(target_kind_),
    asReal(norm_A_), asReal(tol_),
    asLogical(use_retained_av_cache_) == TRUE,
    asLogical(use_deflation_) == TRUE
  );
}

#include "projection/golub_kahan_ritz.hpp"

extern "C" SEXP eigencore_golub_kahan_dense_fit(SEXP A_, SEXP maxit_, SEXP start_,
                                                SEXP rank_, SEXP target_kind_,
                                                SEXP tol_, SEXP projected_stop_,
                                                SEXP reorth_u_, SEXP reorth_v_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  if (LENGTH(start_) != n) {
    error("non-conformable dense Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  const int reorthogonalize_u = asLogical(reorth_u_) == TRUE;
  const int reorthogonalize_v = asLogical(reorth_v_) == TRUE;
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }
  if (rank < 1 || rank > maxit) {
    error("rank must be between 1 and maxit");
  }

  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(maxit), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(maxit), 0.0);
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(maxit) +
      static_cast<int64_t>(2 * m + 2 * n) +
      static_cast<int64_t>(maxit)
    );

  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int iterations = 0;
  int matvecs = 0;
  int projected_stop = 0;
  int projected_nconv = 0;
  double projected_max_residual = R_PosInf;
  int projected_checks = 0;
  double projected_seconds = 0.0;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  int reorthogonalization_passes = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_dense_apply, m, n,
                                            maxit, rank, target_kind, tol,
                                            enable_projected_stop, 1,
                                            REAL(start_),
                                            U_work.data(), V_work.data(),
                                            alpha_work.data(), beta_work.data(),
                                            &iterations, &matvecs,
                                            &projected_stop, &projected_nconv,
                                            &projected_max_residual,
                                            &projected_checks,
                                            &projected_seconds,
                                            &stage_apply_seconds,
                                            &stage_recurrence_seconds,
                                            &stage_reorthogonalization_seconds,
                                            &reorthogonalization_passes,
                                            reorthogonalize_u,
                                            reorthogonalize_v);
  if (status != 0) {
    error("native dense Golub-Kahan failed with status=%d", status);
  }

  SEXP ritz_ = PROTECT(eigencore_golub_kahan_ritz_from_ptr(
    U_work.data(), V_work.data(), m, n, iterations, alpha_work.data(),
    beta_work.data(), rank, target_kind
  ));
  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 7, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 9, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 11, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 15, ScalarLogical(FALSE));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 6, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 7, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 8, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 9, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 11, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 15, mkChar("basis_returned"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(3);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_csc_fit(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                              SEXP maxit_, SEXP start_,
                                              SEXP rank_, SEXP target_kind_,
                                              SEXP tol_, SEXP projected_stop_,
                                              SEXP reorth_u_, SEXP reorth_v_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC Golub-Kahan inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  if (LENGTH(start_) != n) {
    error("non-conformable CSC Golub-Kahan inputs");
  }
  const int limit = (m < n) ? m : n;
  const int maxit = static_cast<int>(asInteger(maxit_));
  const int rank = static_cast<int>(asInteger(rank_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int enable_projected_stop = asLogical(projected_stop_) == TRUE;
  const int reorthogonalize_u = asLogical(reorth_u_) == TRUE;
  const int reorthogonalize_v = asLogical(reorth_v_) == TRUE;
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }
  if (rank < 1 || rank > maxit) {
    error("rank must be between 1 and maxit");
  }

  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(maxit), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(maxit), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(maxit), 0.0);
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(maxit) +
      static_cast<int64_t>(2 * m + 2 * n)
    );

  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int iterations = 0;
  int matvecs = 0;
  int projected_stop = 0;
  int projected_nconv = 0;
  double projected_max_residual = R_PosInf;
  int projected_checks = 0;
  double projected_seconds = 0.0;
  double stage_apply_seconds = 0.0;
  double stage_recurrence_seconds = 0.0;
  double stage_reorthogonalization_seconds = 0.0;
  int reorthogonalization_passes = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_csc_apply, m, n,
                                            maxit, rank, target_kind, tol,
                                            enable_projected_stop, 0,
                                            REAL(start_),
                                            U_work.data(), V_work.data(),
                                            alpha_work.data(), beta_work.data(),
                                            &iterations, &matvecs,
                                            &projected_stop, &projected_nconv,
                                            &projected_max_residual,
                                            &projected_checks,
                                            &projected_seconds,
                                            &stage_apply_seconds,
                                            &stage_recurrence_seconds,
                                            &stage_reorthogonalization_seconds,
                                            &reorthogonalization_passes,
                                            reorthogonalize_u,
                                            reorthogonalize_v);
  if (status != 0) {
    error("native CSC Golub-Kahan failed with status=%d", status);
  }

  SEXP ritz_ = PROTECT(eigencore_golub_kahan_ritz_from_ptr(
    U_work.data(), V_work.data(), m, n, iterations, alpha_work.data(),
    beta_work.data(), rank, target_kind
  ));
  SEXP out_ = PROTECT(allocVector(VECSXP, 16));
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 5, ScalarLogical(projected_stop));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(projected_nconv));
  SET_VECTOR_ELT(out_, 7, ScalarReal(projected_max_residual));
  SET_VECTOR_ELT(out_, 8, ScalarInteger(projected_checks));
  SET_VECTOR_ELT(out_, 9, ScalarReal(projected_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 11, ScalarReal(stage_apply_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(stage_recurrence_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarReal(stage_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 15, ScalarLogical(FALSE));
  SEXP names_ = PROTECT(allocVector(STRSXP, 16));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("iterations"));
  SET_STRING_ELT(names_, 4, mkChar("matvecs"));
  SET_STRING_ELT(names_, 5, mkChar("projected_stop"));
  SET_STRING_ELT(names_, 6, mkChar("projected_nconv"));
  SET_STRING_ELT(names_, 7, mkChar("projected_max_residual"));
  SET_STRING_ELT(names_, 8, mkChar("projected_checks"));
  SET_STRING_ELT(names_, 9, mkChar("projected_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 11, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 15, mkChar("basis_returned"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(3);
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
    const double abstol = 0.0;
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
    const double abstol = 0.0;
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

static int orthonormalize_columns_inplace(double* Q, int n, int cols,
                                          double* tmp, double tol) {
  int rank = 0;
  for (int col = 0; col < cols; ++col) {
    double* q_col = Q + static_cast<int64_t>(rank) * n;
    if (rank != col) {
      std::memcpy(q_col, Q + static_cast<int64_t>(col) * n,
                  sizeof(double) * static_cast<size_t>(n));
    }
    trl_orthogonalise(nullptr, 0, Q, rank, q_col, tmp, n);
    const double q_norm = trl_norm2(q_col, n);
    if (q_norm <= tol) {
      continue;
    }
    const double inv_norm = 1.0 / q_norm;
    for (int row = 0; row < n; ++row) {
      q_col[row] *= inv_norm;
    }
    ++rank;
  }
  return rank;
}

static void csc_transpose_apply_vec(const int* Ai, const int* Ap,
                                    const double* Ax, int m, int n,
                                    const double* x, double* y) {
  (void) m;
  std::fill(y, y + n, 0.0);
  for (int col = 0; col < n; ++col) {
    long double acc = 0.0L;
    for (int jj = Ap[col]; jj < Ap[col + 1]; ++jj) {
      acc += static_cast<long double>(Ax[jj]) * x[Ai[jj]];
    }
    y[col] = static_cast<double>(acc);
  }
}

static void csc_forward_apply_vec(const int* Ai, const int* Ap,
                                  const double* Ax, int m, int n,
                                  const double* x, double* y) {
  std::fill(y, y + m, 0.0);
  for (int col = 0; col < n; ++col) {
    const double x_col = x[col];
    if (x_col == 0.0) {
      continue;
    }
    for (int jj = Ap[col]; jj < Ap[col + 1]; ++jj) {
      y[Ai[jj]] += Ax[jj] * x_col;
    }
  }
}

static void csc_left_normal_apply_vec(const int* Ai, const int* Ap,
                                      const double* Ax, int m, int n,
                                      const double* x, double* y,
                                      double* tmp_n) {
  csc_transpose_apply_vec(Ai, Ap, Ax, m, n, x, tmp_n);
  csc_forward_apply_vec(Ai, Ap, Ax, m, n, tmp_n, y);
}

struct CSCLeftNormalOperator {
  int m;
  int n;
  const int* Ai;
  const int* Ap;
  const double* Ax;
  std::vector<double> tmp_n;
  std::vector<double> tmp_m;
};

static int csc_left_normal_apply_block(void* impl,
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
  CSCLeftNormalOperator* normal = static_cast<CSCLeftNormalOperator*>(impl);
  if (ldx < normal->m || ldy < normal->m || block_cols < 1) {
    return -1;
  }
  const size_t need_n = static_cast<size_t>(normal->n) * static_cast<size_t>(block_cols);
  const size_t need_m = static_cast<size_t>(normal->m) * static_cast<size_t>(block_cols);
  if (normal->tmp_n.size() < need_n || normal->tmp_m.size() < need_m) {
    return -2;
  }
  scale_or_zero_output(Y, normal->m, block_cols, beta);
  for (int64_t col = 0; col < block_cols; ++col) {
    const double* x_col = X + col * ldx;
    double* tmp_n_col = normal->tmp_n.data() + col * normal->n;
    double* tmp_m_col = normal->tmp_m.data() + col * normal->m;
    csc_left_normal_apply_vec(
      normal->Ai, normal->Ap, normal->Ax, normal->m, normal->n,
      x_col, tmp_m_col, tmp_n_col
    );
    double* y_col = Y + col * ldy;
    for (int row = 0; row < normal->m; ++row) {
      y_col[row] += alpha * tmp_m_col[row];
    }
  }
  return 0;
}

static int csc_implicit_left_normal_lanczos_attempt(const int* Ai,
                                                   const int* Ap,
                                                   const double* Ax,
                                                   int m,
                                                   int n,
                                                   int rank,
                                                   double tol,
                                                   double norm_A,
                                                   double* values,
                                                   double* U,
                                                   int* iterations_out,
                                                   double* max_backward_error_out) {
  if (m < 2 || rank < 1 || rank > m) {
    return 0;
  }
  int max_steps = 2 * rank + 1;
  if (max_steps < 20) {
    max_steps = 20;
  }
  if (max_steps > m) {
    max_steps = m;
  }
  if (max_steps < rank) {
    return 0;
  }

  std::vector<double> start(static_cast<size_t>(m), 0.0);
  for (int row = 0; row < m; ++row) {
    const uint32_t key = static_cast<uint32_t>((row + 1) * 1103515245u) ^
      static_cast<uint32_t>(rank * 2654435761u);
    start[row] = (key & 1u) ? 1.0 : -1.0;
  }
  CSCLeftNormalOperator normal = {
    m, n, Ai, Ap, Ax,
    std::vector<double>(static_cast<size_t>(n) * static_cast<size_t>(max_steps), 0.0),
    std::vector<double>(static_cast<size_t>(m) * static_cast<size_t>(max_steps), 0.0)
  };
  std::vector<double> lambda(static_cast<size_t>(rank), 0.0);
  std::vector<double> residuals(static_cast<size_t>(rank), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(rank), 0);
  int n_locked = 0;
  int iterations = 0;
  int matvecs = 0;
  int restarts = 0;
  int m_active_final = 0;
  int locking_events = 0;
  int ortho_passes = 0;
  int64_t operator_allocations = 0;
  int64_t operator_bytes = 0;
  NativeBlockStageSeconds stage;
  const int rc = native_block_thick_restart_lanczos_run(
    &normal, csc_left_normal_apply_block, m, rank, max_steps, 1,
    1, tol * 0.1, 1000, norm_A * norm_A, 0, start.data(), U, lambda.data(),
    residuals.data(), converged.data(), &n_locked, &iterations, &matvecs,
    &restarts, &m_active_final, &locking_events, &ortho_passes,
    &operator_allocations, &operator_bytes, &stage, nullptr
  );
  if (rc != 0) {
    return 0;
  }

  std::vector<double> Gu_exact(static_cast<size_t>(m), 0.0);
  std::vector<double> tmp_n_exact(static_cast<size_t>(n), 0.0);
  const double scale_value_native = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double native_max_backward = 0.0;
  for (int out_col = 0; out_col < rank; ++out_col) {
    const double lambda_i = lambda[static_cast<size_t>(out_col)] > 0.0 ?
      lambda[static_cast<size_t>(out_col)] : 0.0;
    const double sigma = sqrt(lambda_i);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    csc_left_normal_apply_vec(
      Ai, Ap, Ax, m, n, U + static_cast<int64_t>(out_col) * m,
      Gu_exact.data(), tmp_n_exact.data()
    );
    long double residual2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double diff =
        (Gu_exact[static_cast<size_t>(row)] -
         lambda_i * U[row + static_cast<int64_t>(out_col) * m]) / sigma;
      residual2 += static_cast<long double>(diff) * diff;
    }
    const double backward = sqrt(static_cast<double>(residual2)) / scale_value_native;
    if (backward > native_max_backward) {
      native_max_backward = backward;
    }
    values[out_col] = lambda_i;
  }
  if (iterations_out != nullptr) {
    *iterations_out = matvecs;
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = native_max_backward;
  }
  if (native_max_backward <= tol) {
    return 1;
  }
  return 0;
}

static int gram_krylov_left_normal_attempt(const double* gram,
                                           int m,
                                           int rank,
                                           double tol,
                                           double norm_A,
                                           double* values,
                                           double* U,
                                           int* iterations_out,
                                           double* max_backward_error_out) {
  if (m < 2 || rank < 1 || rank > m) {
    return 0;
  }
  int max_steps = std::max(45, 6 * rank + 15);
  if (max_steps > m) {
    max_steps = m;
  }
  if (max_steps < rank) {
    return 0;
  }

  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(max_steps), 0.0);
  std::vector<double> z(static_cast<size_t>(m), 0.0);
  std::vector<double> tmp(static_cast<size_t>(m), 0.0);
  std::vector<double> alpha(static_cast<size_t>(max_steps), 0.0);
  std::vector<double> beta(static_cast<size_t>(max_steps), 0.0);

  for (int row = 0; row < m; ++row) {
    const uint32_t key = static_cast<uint32_t>((row + 1) * 1103515245u) ^
      static_cast<uint32_t>(rank * 2654435761u);
    Q[row] = (key & 1u) ? 1.0 : -1.0;
  }
  double q_norm = trl_norm2(Q.data(), m);
  if (q_norm <= 100.0 * DBL_EPSILON) {
    return 0;
  }
  for (int row = 0; row < m; ++row) {
  Q[row] /= q_norm;
  }

  const char trans_N = 'N';
  const int inc_one = 1;
  const double one = 1.0;
  const double zero = 0.0;
  int active = 0;
  for (int step = 0; step < max_steps; ++step) {
    const double* q = Q.data() + static_cast<int64_t>(step) * m;
    F77_CALL(dgemv)(&trans_N, &m, &m, &one, gram, &m, q, &inc_one,
                    &zero, z.data(), &inc_one FCONE);
    long double dot = 0.0L;
    for (int row = 0; row < m; ++row) {
      dot += static_cast<long double>(q[row]) * z[static_cast<size_t>(row)];
    }
    alpha[static_cast<size_t>(step)] = static_cast<double>(dot);
    if (step > 0) {
      const double* q_prev = Q.data() + static_cast<int64_t>(step - 1) * m;
      const double b_prev = beta[static_cast<size_t>(step - 1)];
      for (int row = 0; row < m; ++row) {
        z[static_cast<size_t>(row)] -=
          alpha[static_cast<size_t>(step)] * q[row] + b_prev * q_prev[row];
      }
    } else {
      for (int row = 0; row < m; ++row) {
        z[static_cast<size_t>(row)] -= alpha[static_cast<size_t>(step)] * q[row];
      }
    }

    trl_orthogonalise(nullptr, 0, Q.data(), step + 1, z.data(), tmp.data(), m);
    const double b = trl_norm2(z.data(), m);
    active = step + 1;
    if (step + 1 >= max_steps || b <= 100.0 * DBL_EPSILON * (norm_A * norm_A + 1.0)) {
      break;
    }
    beta[static_cast<size_t>(step)] = b;
    double* q_next = Q.data() + static_cast<int64_t>(step + 1) * m;
    const double inv_b = 1.0 / b;
    for (int row = 0; row < m; ++row) {
      q_next[row] = z[static_cast<size_t>(row)] * inv_b;
    }
  }
  if (active < rank) {
    return 0;
  }

  std::vector<double> T(static_cast<size_t>(active) * static_cast<size_t>(active), 0.0);
  std::vector<double> theta(static_cast<size_t>(active), 0.0);
  for (int col = 0; col < active; ++col) {
    T[col + static_cast<int64_t>(col) * active] = alpha[static_cast<size_t>(col)];
    if (col + 1 < active) {
      const double b = beta[static_cast<size_t>(col)];
      T[col + static_cast<int64_t>(col + 1) * active] = b;
      T[(col + 1) + static_cast<int64_t>(col) * active] = b;
    }
  }
  int lwork = trl_dsyev_query(active);
  if (lwork < 3 * active) {
    lwork = 3 * active;
  }
  std::vector<double> work(static_cast<size_t>(lwork));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  F77_CALL(dsyev)(&jobz, &uplo, &active, T.data(), &active,
                  theta.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    return 0;
  }

  std::vector<double> Gu(static_cast<size_t>(m), 0.0);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double max_backward = 0.0;
  for (int out_col = 0; out_col < rank; ++out_col) {
    const int src_col = active - 1 - out_col;
    const double lambda = theta[static_cast<size_t>(src_col)] > 0.0 ?
      theta[static_cast<size_t>(src_col)] : 0.0;
    const double sigma = sqrt(lambda);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    values[out_col] = lambda;
    double* u_col = U + static_cast<int64_t>(out_col) * m;
    std::fill(u_col, u_col + m, 0.0);
    for (int basis_col = 0; basis_col < active; ++basis_col) {
      const double coeff = T[basis_col + static_cast<int64_t>(src_col) * active];
      const double* q_col = Q.data() + static_cast<int64_t>(basis_col) * m;
      for (int row = 0; row < m; ++row) {
        u_col[row] += coeff * q_col[row];
      }
    }
    F77_CALL(dgemv)(&trans_N, &m, &m, &one, gram, &m, u_col, &inc_one,
                    &zero, Gu.data(), &inc_one FCONE);
    long double residual2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double residual = (Gu[static_cast<size_t>(row)] - lambda * u_col[row]) / sigma;
      residual2 += static_cast<long double>(residual) * residual;
    }
    const double backward = sqrt(static_cast<double>(residual2)) / scale_value;
    if (backward > max_backward) {
      max_backward = backward;
    }
  }
  if (iterations_out != nullptr) {
    *iterations_out = active;
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = max_backward;
  }
  return max_backward <= tol ? 1 : 0;
}

static int gram_top_subspace_attempt(const double* gram,
                                     int m,
                                     int rank,
                                     double tol,
                                     double norm_A,
                                     double* values,
                                     double* U,
                                     double* max_backward_error_out) {
  if (m < 2 || rank < 1 || rank > m) {
    return 0;
  }
  int subspace = rank + 5;
  if (subspace > m) subspace = m;
  if (subspace <= rank && rank < m) {
    return 0;
  }

  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(subspace), 0.0);
  std::vector<double> Z(static_cast<size_t>(m) * static_cast<size_t>(subspace), 0.0);
  std::vector<double> H(static_cast<size_t>(subspace) * static_cast<size_t>(subspace), 0.0);
  std::vector<double> theta(static_cast<size_t>(subspace), 0.0);
  std::vector<double> tmp(static_cast<size_t>(m), 0.0);

  for (int col = 0; col < subspace; ++col) {
    for (int row = 0; row < m; ++row) {
      const uint32_t key = static_cast<uint32_t>((row + 1) * 1103515245u) ^
        static_cast<uint32_t>((col + 1) * 2654435761u);
      Q[row + static_cast<int64_t>(col) * m] =
        ((key & 1u) ? 1.0 : -1.0) *
        (1.0 + static_cast<double>((key >> 8) & 15u) / 16.0);
    }
  }
  int q_rank = orthonormalize_columns_inplace(
    Q.data(), m, subspace, tmp.data(), 100.0 * DBL_EPSILON
  );
  if (q_rank < rank) {
    return 0;
  }
  subspace = q_rank;

  const char trans_N = 'N';
  const char trans_T = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  const int max_iter = 8;
  for (int iter = 0; iter < max_iter; ++iter) {
    F77_CALL(dgemm)(&trans_N, &trans_N, &m, &subspace, &m,
                    &one, gram, &m, Q.data(), &m,
                    &zero, Z.data(), &m FCONE FCONE);
    std::memcpy(Q.data(), Z.data(),
                sizeof(double) * static_cast<size_t>(m) *
                  static_cast<size_t>(subspace));
    q_rank = orthonormalize_columns_inplace(
      Q.data(), m, subspace, tmp.data(), 100.0 * DBL_EPSILON
    );
    if (q_rank < rank) {
      return 0;
    }
    subspace = q_rank;
  }

  F77_CALL(dgemm)(&trans_N, &trans_N, &m, &subspace, &m,
                  &one, gram, &m, Q.data(), &m,
                  &zero, Z.data(), &m FCONE FCONE);
  F77_CALL(dgemm)(&trans_T, &trans_N, &subspace, &subspace, &m,
                  &one, Q.data(), &m, Z.data(), &m,
                  &zero, H.data(), &subspace FCONE FCONE);
  symmetrize_packed_square(H.data(), subspace);

  int lwork = trl_dsyev_query(subspace);
  if (lwork < 3 * subspace) lwork = 3 * subspace;
  std::vector<double> work(static_cast<size_t>(lwork));
  char jobz = 'V';
  char uplo = 'U';
  int info = 0;
  F77_CALL(dsyev)(&jobz, &uplo, &subspace, H.data(), &subspace,
                  theta.data(), work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    return 0;
  }

  for (int out_col = 0; out_col < rank; ++out_col) {
    const int src_col = subspace - 1 - out_col;
    values[out_col] = theta[static_cast<size_t>(src_col)];
    for (int row = 0; row < m; ++row) {
      long double sum = 0.0L;
      for (int qcol = 0; qcol < subspace; ++qcol) {
        sum += static_cast<long double>(
          Q[row + static_cast<int64_t>(qcol) * m]
        ) * H[qcol + static_cast<int64_t>(src_col) * subspace];
      }
      U[row + static_cast<int64_t>(out_col) * m] = static_cast<double>(sum);
    }
  }

  std::vector<double> GU(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  F77_CALL(dgemm)(&trans_N, &trans_N, &m, &rank, &m,
                  &one, gram, &m, U, &m,
                  &zero, GU.data(), &m FCONE FCONE);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;
  double max_backward = 0.0;
  for (int col = 0; col < rank; ++col) {
    const double lambda = values[col] > 0.0 ? values[col] : 0.0;
    const double sigma = sqrt(lambda);
    if (sigma <= 100.0 * DBL_EPSILON) {
      return 0;
    }
    long double residual2 = 0.0L;
    for (int row = 0; row < m; ++row) {
      const double residual = (GU[row + static_cast<int64_t>(col) * m] -
        lambda * U[row + static_cast<int64_t>(col) * m]) / sigma;
      residual2 += static_cast<long double>(residual) * residual;
    }
    const double backward = sqrt(static_cast<double>(residual2)) / scale_value;
    if (backward > max_backward) {
      max_backward = backward;
    }
  }
  if (max_backward_error_out != nullptr) {
    *max_backward_error_out = max_backward;
  }
  return max_backward <= tol ? 1 : 0;
}

extern "C" SEXP eigencore_csc_left_gram_svd(SEXP i_, SEXP p_, SEXP x_,
                                            SEXP dim_, SEXP rank_, SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  int rank = asInteger(rank_);
  if (m < 1 || n < 1 || rank < 1) {
    error("invalid CSC Gram SVD dimensions");
  }
  if (rank > m) {
    rank = m;
  }
  const double tol = asReal(tol_);
  const int* Ai = INTEGER(i_);
  const int* Ap = INTEGER(p_);
  const double* Ax = REAL(x_);

  auto stage_timer = native_timer_now();
  double stage_gram_seconds = 0.0;
  double stage_eigensolve_seconds = 0.0;
  double stage_vector_form_seconds = 0.0;
  double stage_diagnostics_seconds = 0.0;

  double frob2 = 0.0;
  for (int idx = 0; idx < Ap[n]; ++idx) {
    frob2 += Ax[idx] * Ax[idx];
  }
  const double norm_A = sqrt(frob2 > 0.0 ? frob2 : 0.0);
  const double scale_value = norm_A > DBL_EPSILON ? norm_A : DBL_EPSILON;

  std::vector<double> values(static_cast<size_t>(rank), 0.0);
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, rank));
  std::vector<double> gram;
  int used_implicit_lanczos = 0;
  int implicit_lanczos_iterations = 0;
  double implicit_lanczos_max_backward_error = R_PosInf;
  int used_subspace_eigensolve = 0;
  int used_gram_krylov = 0;
  int gram_krylov_iterations = 0;
  double subspace_max_backward_error = R_PosInf;
  const char* lapack_eigensolver = "lapack_dsyevr";

  SEXP implicit_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_normal_lanczos_attempt"));
  const int attempt_implicit_lanczos = asLogical(implicit_option_) == TRUE;
  if (attempt_implicit_lanczos && m <= 128 && rank <= 16) {
    stage_timer = native_timer_now();
    used_implicit_lanczos = csc_implicit_left_normal_lanczos_attempt(
      Ai, Ap, Ax, m, n, rank, tol, norm_A, values.data(), REAL(u_),
      &implicit_lanczos_iterations, &implicit_lanczos_max_backward_error
    );
    stage_eigensolve_seconds = native_timer_elapsed(stage_timer);
  }

  if (!used_implicit_lanczos) {
    stage_timer = native_timer_now();
    gram.assign(static_cast<size_t>(m) * static_cast<size_t>(m), 0.0);
    for (int col = 0; col < n; ++col) {
      const int start = Ap[col];
      const int end = Ap[col + 1];
      for (int aa = start; aa < end; ++aa) {
        const int row_a = Ai[aa];
        const double x_a = Ax[aa];
        gram[row_a + static_cast<int64_t>(row_a) * m] += x_a * x_a;
        for (int bb = aa + 1; bb < end; ++bb) {
          const int row_b = Ai[bb];
          const double update = x_a * Ax[bb];
          gram[row_a + static_cast<int64_t>(row_b) * m] += update;
          gram[row_b + static_cast<int64_t>(row_a) * m] += update;
        }
      }
    }
    stage_gram_seconds = native_timer_elapsed(stage_timer);

    stage_timer = native_timer_now();
    SEXP krylov_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_gram_krylov_attempt"));
    const int attempt_gram_krylov = asLogical(krylov_option_) == TRUE;
    SEXP subspace_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_gram_subspace_attempt"));
    const int attempt_subspace_eigensolve = asLogical(subspace_option_) == TRUE;
    if (attempt_gram_krylov && m <= 90 && rank <= 8) {
      used_gram_krylov = gram_krylov_left_normal_attempt(
        gram.data(), m, rank, tol, norm_A, values.data(), REAL(u_),
        &gram_krylov_iterations, &subspace_max_backward_error
      );
    }
    if (!used_gram_krylov && attempt_subspace_eigensolve && m <= 128 && rank <= 16) {
      used_subspace_eigensolve = gram_top_subspace_attempt(
        gram.data(), m, rank, tol, norm_A, values.data(), REAL(u_),
        &subspace_max_backward_error
      );
    }
    if (m > 0 && !used_gram_krylov && !used_subspace_eigensolve) {
      std::vector<double> work_matrix(static_cast<size_t>(m) * static_cast<size_t>(m));
      std::vector<double> values_work(static_cast<size_t>(m), 0.0);
      char uplo = 'U';
      int info = 0;
      std::memcpy(work_matrix.data(), gram.data(),
                  sizeof(double) * static_cast<size_t>(m) * static_cast<size_t>(m));
      SEXP dsyevx_option_ = Rf_GetOption1(Rf_install("eigencore.csc_left_gram_dsyevx_attempt"));
      const int attempt_dsyevx = asLogical(dsyevx_option_) == TRUE;
      if (attempt_dsyevx && m <= 128 && rank <= 16) {
        lapack_eigensolver = "lapack_dsyevx";
        char jobz = 'V';
        char range = 'I';
        double vl = 0.0;
        double vu = 0.0;
        const double abstol = 0.0;
        int il = m - rank + 1;
        int iu = m;
        int m_found = 0;
        int lwork = 8 * m;
        if (lwork < 1) {
          lwork = 1;
        }
        std::vector<double> work(static_cast<size_t>(lwork));
        std::vector<int> iwork(static_cast<size_t>(5 * m));
        std::vector<int> ifail(static_cast<size_t>(m));
        F77_CALL(dsyevx)(&jobz, &range, &uplo, &m, work_matrix.data(), &m,
                         &vl, &vu, &il, &iu, &abstol,
                         &m_found, values_work.data(), REAL(u_), &m,
                         work.data(), &lwork, iwork.data(), ifail.data(),
                         &info FCONE FCONE FCONE);
        if (info != 0 || m_found != rank) {
          UNPROTECT(1);
          error("LAPACK dsyevx failed with info=%d, found=%d", info, m_found);
        }
        for (int left = 0, right = rank - 1; left < right; ++left, --right) {
          const double tmp_value = values_work[left];
          values_work[left] = values_work[right];
          values_work[right] = tmp_value;
          for (int row = 0; row < m; ++row) {
            const int64_t lpos = row + static_cast<int64_t>(left) * m;
            const int64_t rpos = row + static_cast<int64_t>(right) * m;
            const double tmp_vec = REAL(u_)[lpos];
            REAL(u_)[lpos] = REAL(u_)[rpos];
            REAL(u_)[rpos] = tmp_vec;
          }
        }
        for (int col = 0; col < rank; ++col) {
          values[static_cast<size_t>(col)] = values_work[static_cast<size_t>(col)];
        }
      } else if (m <= 128 && rank >= 16) {
        lapack_eigensolver = "lapack_dsyevd";
        char jobz = 'V';
        int lwork = -1;
        int liwork = -1;
        double work_query = 0.0;
        int iwork_query = 0;
        F77_CALL(dsyevd)(&jobz, &uplo, &m, work_matrix.data(), &m,
                         values_work.data(), &work_query, &lwork,
                         &iwork_query, &liwork, &info FCONE FCONE);
        if (info != 0) {
          UNPROTECT(1);
          error("LAPACK dsyevd workspace query failed with info=%d", info);
        }
        lwork = static_cast<int>(work_query);
        liwork = iwork_query;
        if (lwork < 1 + 6 * m + 2 * m * m) {
          lwork = 1 + 6 * m + 2 * m * m;
        }
        if (liwork < 3 + 5 * m) {
          liwork = 3 + 5 * m;
        }
        std::vector<double> work(static_cast<size_t>(lwork));
        std::vector<int> iwork(static_cast<size_t>(liwork));
        F77_CALL(dsyevd)(&jobz, &uplo, &m, work_matrix.data(), &m,
                         values_work.data(), work.data(), &lwork,
                         iwork.data(), &liwork, &info FCONE FCONE);
        if (info != 0) {
          UNPROTECT(1);
          error("LAPACK dsyevd failed with info=%d", info);
        }
        for (int col = 0; col < rank; ++col) {
          const int source_col = m - 1 - col;
          values[static_cast<size_t>(col)] = values_work[static_cast<size_t>(source_col)];
          std::memcpy(
            REAL(u_) + static_cast<int64_t>(col) * m,
            work_matrix.data() + static_cast<int64_t>(source_col) * m,
            sizeof(double) * static_cast<size_t>(m)
          );
        }
      } else {
        std::vector<int> isuppz(static_cast<size_t>(2 * rank), 0);
        char jobz = 'V';
        char range = 'I';
        double vl = 0.0;
        double vu = 0.0;
        const double abstol = 0.0;
        int il = m - rank + 1;
        int iu = m;
        int m_found = 0;
        int lwork = 26 * m;
        int liwork = 10 * m;
        std::vector<double> work(static_cast<size_t>(lwork));
        std::vector<int> iwork(static_cast<size_t>(liwork));
        F77_CALL(dsyevr)(&jobz, &range, &uplo, &m, work_matrix.data(), &m,
                         &vl, &vu, &il, &iu, &abstol,
                         &m_found, values_work.data(), REAL(u_), &m,
                         isuppz.data(), work.data(), &lwork,
                         iwork.data(), &liwork, &info FCONE FCONE FCONE);
        if (info != 0 || m_found != rank) {
          UNPROTECT(1);
          error("LAPACK dsyevr failed with info=%d, found=%d", info, m_found);
        }
        for (int left = 0, right = rank - 1; left < right; ++left, --right) {
          const double tmp_value = values_work[left];
          values_work[left] = values_work[right];
          values_work[right] = tmp_value;
          for (int row = 0; row < m; ++row) {
            const int64_t lpos = row + static_cast<int64_t>(left) * m;
            const int64_t rpos = row + static_cast<int64_t>(right) * m;
            const double tmp_vec = REAL(u_)[lpos];
            REAL(u_)[lpos] = REAL(u_)[rpos];
            REAL(u_)[rpos] = tmp_vec;
          }
        }
        for (int col = 0; col < rank; ++col) {
          values[static_cast<size_t>(col)] = values_work[static_cast<size_t>(col)];
        }
      }
      subspace_max_backward_error = R_PosInf;
    }
    stage_eigensolve_seconds = native_timer_elapsed(stage_timer);
  }

  stage_timer = native_timer_now();
  SEXP d_ = PROTECT(allocVector(REALSXP, rank));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, rank));
  std::memset(REAL(v_), 0, sizeof(double) * static_cast<size_t>(n) * rank);

  for (int col = 0; col < rank; ++col) {
    const double sigma = sqrt(values[static_cast<size_t>(col)] > 0.0 ?
                              values[static_cast<size_t>(col)] : 0.0);
    REAL(d_)[col] = sigma;
  }

  for (int acol = 0; acol < n; ++acol) {
    for (int jj = Ap[acol]; jj < Ap[acol + 1]; ++jj) {
      const int row = Ai[jj];
      const double val = Ax[jj];
      for (int scol = 0; scol < rank; ++scol) {
        REAL(v_)[acol + static_cast<int64_t>(scol) * n] +=
          val * REAL(u_)[row + static_cast<int64_t>(scol) * m];
      }
    }
  }
  for (int scol = 0; scol < rank; ++scol) {
    const double sigma = REAL(d_)[scol];
    const double inv_sigma = sigma > 100.0 * DBL_EPSILON ? 1.0 / sigma : 0.0;
    for (int col = 0; col < n; ++col) {
      REAL(v_)[col + static_cast<int64_t>(scol) * n] *= inv_sigma;
    }
  }
  stage_vector_form_seconds = native_timer_elapsed(stage_timer);

  stage_timer = native_timer_now();
  SEXP left_ = PROTECT(allocVector(REALSXP, rank));
  SEXP right_ = PROTECT(allocVector(REALSXP, rank));
  SEXP combined_ = PROTECT(allocVector(REALSXP, rank));
  SEXP backward_ = PROTECT(allocVector(REALSXP, rank));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, rank));
  SEXP scale_ = PROTECT(allocVector(REALSXP, rank));
  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));

  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  std::vector<double> gram_u_small(static_cast<size_t>(rank) * rank, 0.0);
  std::vector<double> gram_v_small(static_cast<size_t>(rank) * rank, 0.0);
  std::vector<double> gu;
  std::vector<double> gu_block;
  F77_CALL(dgemm)(&trans, &notrans, &rank, &rank, &m,
                  &one, REAL(u_), &m, REAL(u_), &m,
                  &zero, gram_u_small.data(), &rank FCONE FCONE);
  if (used_implicit_lanczos) {
    F77_CALL(dgemm)(&trans, &notrans, &rank, &rank, &n,
                    &one, REAL(v_), &n, REAL(v_), &n,
                    &zero, gram_v_small.data(), &rank FCONE FCONE);
    gu.assign(static_cast<size_t>(m), 0.0);
  } else {
    gu_block.assign(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
    F77_CALL(dgemm)(&notrans, &notrans, &m, &rank, &m,
                    &one, gram.data(), &m, REAL(u_), &m,
                    &zero, gu_block.data(), &m FCONE FCONE);
    for (int col = 0; col < rank; ++col) {
      const double sigma_col = REAL(d_)[col];
      const double inv_col = sigma_col > 100.0 * DBL_EPSILON ? 1.0 / sigma_col : 0.0;
      const double* gu_col = gu_block.data() + static_cast<int64_t>(col) * m;
      for (int row_col = 0; row_col < rank; ++row_col) {
        const double sigma_row = REAL(d_)[row_col];
        const double inv_row = sigma_row > 100.0 * DBL_EPSILON ? 1.0 / sigma_row : 0.0;
        const double* u_row = REAL(u_) + static_cast<int64_t>(row_col) * m;
        long double dot = 0.0L;
        for (int row = 0; row < m; ++row) {
          dot += static_cast<long double>(u_row[row]) * gu_col[row];
        }
        gram_v_small[row_col + static_cast<int64_t>(col) * rank] =
          static_cast<double>(dot) * inv_row * inv_col;
      }
    }
  }
  REAL(orth_)[0] = max_orthogonality_loss(gram_u_small.data(), rank);
  REAL(orth_)[1] = max_orthogonality_loss(gram_v_small.data(), rank);
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  std::vector<double> atu_check(static_cast<size_t>(n), 0.0);
  for (int scol = 0; scol < rank; ++scol) {
    const double sigma = REAL(d_)[scol];
    const double lambda = sigma * sigma;
    const double* gu_col = nullptr;
    if (used_implicit_lanczos) {
      csc_forward_apply_vec(
        Ai, Ap, Ax, m, n,
        REAL(v_) + static_cast<int64_t>(scol) * n,
        gu.data()
      );
      csc_transpose_apply_vec(
        Ai, Ap, Ax, m, n,
        REAL(u_) + static_cast<int64_t>(scol) * m,
        atu_check.data()
      );
      gu_col = gu.data();
    } else {
      gu_col = gu_block.data() + static_cast<int64_t>(scol) * m;
    }
    long double left_sum = 0.0L;
    long double right_sum = 0.0L;
    const double inv_sigma = sigma > 100.0 * DBL_EPSILON ? 1.0 / sigma : 0.0;
    for (int row = 0; row < m; ++row) {
      const double residual = used_implicit_lanczos
        ? (gu_col[row] -
            sigma * REAL(u_)[row + static_cast<int64_t>(scol) * m])
        : (gu_col[row] -
            lambda * REAL(u_)[row + static_cast<int64_t>(scol) * m]) * inv_sigma;
      left_sum += static_cast<long double>(residual) * residual;
    }
    if (used_implicit_lanczos) {
      for (int row = 0; row < n; ++row) {
        const double residual = atu_check[static_cast<size_t>(row)] -
          sigma * REAL(v_)[row + static_cast<int64_t>(scol) * n];
        right_sum += static_cast<long double>(residual) * residual;
      }
    }
    const double left = sqrt(static_cast<double>(left_sum));
    const double right = sqrt(static_cast<double>(right_sum));
    REAL(left_)[scol] = left;
    REAL(right_)[scol] = right;
    REAL(combined_)[scol] = sqrt(left * left + right * right);
    REAL(scale_)[scol] = scale_value;
    REAL(backward_)[scol] = REAL(combined_)[scol] / scale_value;
    LOGICAL(converged_)[scol] = (REAL(backward_)[scol] <= tol) ? TRUE : FALSE;
  }
  stage_diagnostics_seconds = native_timer_elapsed(stage_timer);

  SEXP diagnostics_ = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(diagnostics_, 0, left_);
  SET_VECTOR_ELT(diagnostics_, 1, right_);
  SET_VECTOR_ELT(diagnostics_, 2, combined_);
  SET_VECTOR_ELT(diagnostics_, 3, backward_);
  SET_VECTOR_ELT(diagnostics_, 4, orth_);
  SET_VECTOR_ELT(diagnostics_, 5, converged_);
  SET_VECTOR_ELT(diagnostics_, 6, scale_);
  SEXP diag_names_ = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(diag_names_, 0, mkChar("left"));
  SET_STRING_ELT(diag_names_, 1, mkChar("right"));
  SET_STRING_ELT(diag_names_, 2, mkChar("combined"));
  SET_STRING_ELT(diag_names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(diag_names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(diag_names_, 5, mkChar("converged"));
  SET_STRING_ELT(diag_names_, 6, mkChar("scale"));
  setAttrib(diagnostics_, R_NamesSymbol, diag_names_);

  SEXP stage_ = PROTECT(allocVector(REALSXP, 4));
  REAL(stage_)[0] = stage_gram_seconds;
  REAL(stage_)[1] = stage_eigensolve_seconds;
  REAL(stage_)[2] = stage_vector_form_seconds;
  REAL(stage_)[3] = stage_diagnostics_seconds;
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(stage_names_, 0, mkChar("gram"));
  SET_STRING_ELT(stage_names_, 1, mkChar("eigensolve"));
  SET_STRING_ELT(stage_names_, 2, mkChar("vector_form"));
  SET_STRING_ELT(stage_names_, 3, mkChar("diagnostics"));
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP eigensolver_ = PROTECT(mkString(
    used_implicit_lanczos ? "implicit_normal_lanczos" :
      (used_gram_krylov ? "explicit_gram_krylov" :
        (used_subspace_eigensolve ? "subspace_iteration" : lapack_eigensolver))
  ));
  SEXP subspace_be_ = PROTECT(ScalarReal(subspace_max_backward_error));
  SEXP implicit_be_ = PROTECT(ScalarReal(implicit_lanczos_max_backward_error));
  SEXP implicit_iter_ = PROTECT(ScalarInteger(implicit_lanczos_iterations));
  SEXP gram_krylov_iter_ = PROTECT(ScalarInteger(gram_krylov_iterations));

  SEXP out_ = PROTECT(allocVector(VECSXP, 10));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, diagnostics_);
  SET_VECTOR_ELT(out_, 4, stage_);
  SET_VECTOR_ELT(out_, 5, eigensolver_);
  SET_VECTOR_ELT(out_, 6, subspace_be_);
  SET_VECTOR_ELT(out_, 7, implicit_be_);
  SET_VECTOR_ELT(out_, 8, implicit_iter_);
  SET_VECTOR_ELT(out_, 9, gram_krylov_iter_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 10));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("diagnostics"));
  SET_STRING_ELT(names_, 4, mkChar("stage_seconds"));
  SET_STRING_ELT(names_, 5, mkChar("eigensolver"));
  SET_STRING_ELT(names_, 6, mkChar("subspace_max_backward_error"));
  SET_STRING_ELT(names_, 7, mkChar("implicit_lanczos_max_backward_error"));
  SET_STRING_ELT(names_, 8, mkChar("implicit_lanczos_iterations"));
  SET_STRING_ELT(names_, 9, mkChar("gram_krylov_iterations"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(22);
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

static void validate_real_vector_length(SEXP x, int n, const char* name) {
  if (!isReal(x) || LENGTH(x) != n) {
    error("%s must be a double vector of length %d", name, n);
  }
}

static void validate_real_matrix_dim(SEXP x, int nrow, int ncol, const char* name) {
  if (!isReal(x)) {
    error("%s must be a double matrix", name);
  }
  SEXP dim = getAttrib(x, R_DimSymbol);
  if (dim == R_NilValue || LENGTH(dim) != 2 ||
      INTEGER(dim)[0] != nrow || INTEGER(dim)[1] != ncol) {
    error("%s must have dimensions %d x %d", name, nrow, ncol);
  }
}

static void validate_irlba_lbd_retained_contract(int m, int n,
                                                 SEXP initial_start_,
                                                 SEXP retained_right_,
                                                 SEXP retained_left_,
                                                 SEXP alpha_, SEXP beta_,
                                                 SEXP random_tails_,
                                                 SEXP work_, SEXP retained_,
                                                 SEXP max_restarts_,
                                                 SEXP rank_,
                                                 SEXP target_kind_,
                                                 SEXP tol_,
                                                 SEXP reorth_policy_) {
  const int work = asInteger(work_);
  const int retained = asInteger(retained_);
  const int max_restarts = asInteger(max_restarts_);
  const int rank = asInteger(rank_);
  const int target_kind = asInteger(target_kind_);
  const double tol = asReal(tol_);
  const int reorth_policy = asInteger(reorth_policy_);
  if (m < 1 || n < 1) {
    error("active operator dimensions must be positive");
  }
  if (rank < 1 || rank > retained) {
    error("rank must be between 1 and retained");
  }
  if (retained < 1 || retained >= work) {
    error("retained must be positive and smaller than work");
  }
  if (work > ((m < n) ? m : n)) {
    error("work must not exceed min(active dimensions)");
  }
  if (max_restarts < 0) {
    error("max_restarts must be non-negative");
  }
  if (target_kind != 1) {
    error("retained IRLBA/LBD currently supports only largest singular values");
  }
  if (!R_FINITE(tol) || tol <= 0.0) {
    error("tol must be a positive finite scalar");
  }
  if (reorth_policy < 1 || reorth_policy > 3) {
    error("reorth_policy must identify a known retained IRLBA/LBD policy");
  }
  validate_real_vector_length(initial_start_, n, "initial_start");
  validate_real_matrix_dim(retained_right_, n, retained, "retained_right");
  validate_real_matrix_dim(retained_left_, m, retained, "retained_left");
  validate_real_vector_length(alpha_, work, "alpha");
  validate_real_vector_length(beta_, work, "beta");
  validate_real_matrix_dim(random_tails_, n, work - retained, "random_tails");
}

static void irlba_lbd_retained_seed(int n,
                                    int retained,
                                    const double* initial_start,
                                    const double* retained_right,
                                    const double* random_tails,
                                    int random_tail_cols,
                                    int attempt,
                                    double* start) {
  std::memset(start, 0, sizeof(double) * static_cast<size_t>(n));
  const int retained_col = (attempt < retained) ? attempt : 0;
  const double* base = retained_right +
    static_cast<int64_t>(retained_col) * static_cast<int64_t>(n);
  long double norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    start[row] = base[row];
    norm2 += static_cast<long double>(start[row]) * start[row];
  }
  if (norm2 <= 100.0L * DBL_EPSILON) {
    norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      start[row] = initial_start[row];
      norm2 += static_cast<long double>(start[row]) * start[row];
    }
  }
  if (attempt > 0 && random_tail_cols > 0) {
    const int tail_col = (attempt - 1) % random_tail_cols;
    const double* tail = random_tails +
      static_cast<int64_t>(tail_col) * static_cast<int64_t>(n);
    for (int row = 0; row < n; ++row) {
      start[row] += 1.0e-3 * tail[row];
    }
  }
  norm2 = 0.0L;
  for (int row = 0; row < n; ++row) {
    norm2 += static_cast<long double>(start[row]) * start[row];
  }
  if (norm2 <= 100.0L * DBL_EPSILON) {
    start[0] = 1.0;
    norm2 = 1.0L;
  }
  const double norm = sqrt(static_cast<double>(norm2));
  for (int row = 0; row < n; ++row) {
    start[row] /= norm;
  }
}

static SEXP irlba_lbd_attempt_history_pack(const std::vector<int>& attempts,
                                           const std::vector<int>& iterations,
                                           const std::vector<int>& matvecs,
                                           const std::vector<int>& warm_started) {
  const int rows = static_cast<int>(attempts.size());
  SEXP attempt_ = PROTECT(allocVector(INTSXP, rows));
  SEXP max_subspace_ = PROTECT(allocVector(INTSXP, rows));
  SEXP iterations_ = PROTECT(allocVector(INTSXP, rows));
  SEXP matvecs_ = PROTECT(allocVector(INTSXP, rows));
  SEXP warm_started_ = PROTECT(allocVector(LGLSXP, rows));
  for (int row = 0; row < rows; ++row) {
    INTEGER(attempt_)[row] = row + 1;
    INTEGER(max_subspace_)[row] = attempts[static_cast<size_t>(row)];
    INTEGER(iterations_)[row] = iterations[static_cast<size_t>(row)];
    INTEGER(matvecs_)[row] = matvecs[static_cast<size_t>(row)];
    LOGICAL(warm_started_)[row] = warm_started[static_cast<size_t>(row)] ? TRUE : FALSE;
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out_, 0, attempt_);
  SET_VECTOR_ELT(out_, 1, max_subspace_);
  SET_VECTOR_ELT(out_, 2, iterations_);
  SET_VECTOR_ELT(out_, 3, matvecs_);
  SET_VECTOR_ELT(out_, 4, warm_started_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(names_, 0, mkChar("attempt"));
  SET_STRING_ELT(names_, 1, mkChar("max_subspace"));
  SET_STRING_ELT(names_, 2, mkChar("iterations"));
  SET_STRING_ELT(names_, 3, mkChar("matvecs"));
  SET_STRING_ELT(names_, 4, mkChar("warm_started"));
  setAttrib(out_, R_NamesSymbol, names_);
  SEXP row_names_ = PROTECT(allocVector(INTSXP, 2));
  INTEGER(row_names_)[0] = NA_INTEGER;
  INTEGER(row_names_)[1] = -rows;
  setAttrib(out_, R_RowNamesSymbol, row_names_);
  SEXP class_ = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(class_, 0, mkChar("data.frame"));
  setAttrib(out_, R_ClassSymbol, class_);
  UNPROTECT(9);
  return out_;
}

template <typename ConfigureOperator>
static SEXP irlba_lbd_retained_impl(ConfigureOperator configure_operator,
                                    int m,
                                    int n,
                                    const double* initial_start,
                                    const double* retained_right,
                                    const double* retained_left,
                                    const double* alpha,
                                    const double* beta,
                                    const double* random_tails,
                                    int work,
                                    int retained,
                                    int max_restarts,
                                    int rank,
                                    int target_kind,
                                    double tol,
                                    int reorth_policy) {
  (void) retained_left;
  (void) alpha;
  (void) beta;
  const int attempts = (max_restarts < 0) ? 1 : max_restarts + 1;
  const int random_tail_cols = work - retained;
  const int reorthogonalize_u =
    (reorth_policy != 1) || (m <= n);
  const int reorthogonalize_v =
    (reorth_policy != 1) || (n <= m);
  const int use_blas_reorthogonalization = 0;
  const int enable_projected_stop = 0;
  const double native_workspace_bytes =
    static_cast<double>(sizeof(double)) *
    static_cast<double>(
      (static_cast<int64_t>(m) + static_cast<int64_t>(n) + 2) *
      static_cast<int64_t>(work) +
      static_cast<int64_t>(2 * m + 2 * n)
    );

  std::vector<double> start(static_cast<size_t>(n), 0.0);
  std::vector<double> U_work(static_cast<size_t>(m) * static_cast<size_t>(work), 0.0);
  std::vector<double> V_work(static_cast<size_t>(n) * static_cast<size_t>(work), 0.0);
  std::vector<double> alpha_work(static_cast<size_t>(work), 0.0);
  std::vector<double> beta_work(static_cast<size_t>(work), 0.0);
  std::vector<int> attempted_subspaces;
  std::vector<int> attempt_iterations;
  std::vector<int> attempt_matvecs;
  std::vector<int> attempt_warm_started;
  attempted_subspaces.reserve(static_cast<size_t>(attempts));
  attempt_iterations.reserve(static_cast<size_t>(attempts));
  attempt_matvecs.reserve(static_cast<size_t>(attempts));
  attempt_warm_started.reserve(static_cast<size_t>(attempts));

  SEXP final_ritz_ = R_NilValue;
  int final_iterations = 0;
  int total_iterations = 0;
  int total_matvecs = 0;
  int total_reorthogonalization_passes = 0;
  double total_apply_seconds = 0.0;
  double total_recurrence_seconds = 0.0;
  double total_reorthogonalization_seconds = 0.0;
  double total_projected_seconds = 0.0;

  for (int attempt = 0; attempt < attempts; ++attempt) {
    irlba_lbd_retained_seed(
      n, retained, initial_start, retained_right, random_tails,
      random_tail_cols, attempt, start.data()
    );
    std::fill(U_work.begin(), U_work.end(), 0.0);
    std::fill(V_work.begin(), V_work.end(), 0.0);
    std::fill(alpha_work.begin(), alpha_work.end(), 0.0);
    std::fill(beta_work.begin(), beta_work.end(), 0.0);

    void* impl = nullptr;
    EigencoreApplyFn apply = nullptr;
    configure_operator(&impl, &apply);
    int iterations = 0;
    int matvecs = 0;
    int projected_stop = 0;
    int projected_nconv = 0;
    double projected_max_residual = R_PosInf;
    int projected_checks = 0;
    double projected_seconds = 0.0;
    double stage_apply_seconds = 0.0;
    double stage_recurrence_seconds = 0.0;
    double stage_reorthogonalization_seconds = 0.0;
    int reorthogonalization_passes = 0;
    const int status = native_golub_kahan_run(
      impl, apply, m, n, work, rank, target_kind, tol,
      enable_projected_stop, use_blas_reorthogonalization,
      start.data(), U_work.data(), V_work.data(), alpha_work.data(),
      beta_work.data(), &iterations, &matvecs, &projected_stop,
      &projected_nconv, &projected_max_residual, &projected_checks,
      &projected_seconds, &stage_apply_seconds, &stage_recurrence_seconds,
      &stage_reorthogonalization_seconds, &reorthogonalization_passes,
      reorthogonalize_u, reorthogonalize_v
    );
    if (status != 0) {
      error("native retained one-sided IRLBA/LBD failed with status=%d", status);
    }
    attempted_subspaces.push_back(work);
    attempt_iterations.push_back(iterations);
    attempt_matvecs.push_back(matvecs);
    attempt_warm_started.push_back(attempt > 0 ? 1 : 0);
    final_iterations = iterations;
    total_iterations += iterations;
    total_matvecs += matvecs;
    total_reorthogonalization_passes += reorthogonalization_passes;
    total_apply_seconds += stage_apply_seconds;
    total_recurrence_seconds += stage_recurrence_seconds;
    total_reorthogonalization_seconds += stage_reorthogonalization_seconds;
    total_projected_seconds += projected_seconds;

    if (attempt + 1 == attempts) {
      final_ritz_ = PROTECT(eigencore_golub_kahan_ritz_from_ptr(
        U_work.data(), V_work.data(), m, n, iterations, alpha_work.data(),
        beta_work.data(), rank, target_kind
      ));
    }
  }

  SEXP history_ = PROTECT(irlba_lbd_attempt_history_pack(
    attempted_subspaces, attempt_iterations, attempt_matvecs, attempt_warm_started
  ));
  SEXP out_ = PROTECT(allocVector(VECSXP, 18));
  SET_VECTOR_ELT(out_, 0, VECTOR_ELT(final_ritz_, 0));
  SET_VECTOR_ELT(out_, 1, VECTOR_ELT(final_ritz_, 1));
  SET_VECTOR_ELT(out_, 2, VECTOR_ELT(final_ritz_, 2));
  SET_VECTOR_ELT(out_, 3, ScalarInteger(final_iterations));
  SET_VECTOR_ELT(out_, 4, ScalarInteger(total_iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(total_matvecs));
  SET_VECTOR_ELT(out_, 6, ScalarInteger(attempts - 1));
  SET_VECTOR_ELT(out_, 7, history_);
  SET_VECTOR_ELT(out_, 8, ScalarReal(native_workspace_bytes));
  SET_VECTOR_ELT(out_, 9, ScalarReal(total_apply_seconds));
  SET_VECTOR_ELT(out_, 10, ScalarReal(total_recurrence_seconds));
  SET_VECTOR_ELT(out_, 11, ScalarReal(total_reorthogonalization_seconds));
  SET_VECTOR_ELT(out_, 12, ScalarReal(total_projected_seconds));
  SET_VECTOR_ELT(out_, 13, ScalarInteger(total_reorthogonalization_passes));
  SET_VECTOR_ELT(out_, 14, ScalarInteger(reorthogonalize_u));
  SET_VECTOR_ELT(out_, 15, ScalarInteger(reorthogonalize_v));
  SET_VECTOR_ELT(out_, 16, ScalarInteger(work));
  SET_VECTOR_ELT(out_, 17, ScalarInteger(retained));
  SEXP names_ = PROTECT(allocVector(STRSXP, 18));
  SET_STRING_ELT(names_, 0, mkChar("d"));
  SET_STRING_ELT(names_, 1, mkChar("u"));
  SET_STRING_ELT(names_, 2, mkChar("v"));
  SET_STRING_ELT(names_, 3, mkChar("final_iterations"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  SET_STRING_ELT(names_, 6, mkChar("restart_count"));
  SET_STRING_ELT(names_, 7, mkChar("attempt_history"));
  SET_STRING_ELT(names_, 8, mkChar("native_workspace_bytes"));
  SET_STRING_ELT(names_, 9, mkChar("stage_apply_seconds"));
  SET_STRING_ELT(names_, 10, mkChar("stage_recurrence_seconds"));
  SET_STRING_ELT(names_, 11, mkChar("stage_reorthogonalization_seconds"));
  SET_STRING_ELT(names_, 12, mkChar("stage_projected_solve_seconds"));
  SET_STRING_ELT(names_, 13, mkChar("reorthogonalization_passes"));
  SET_STRING_ELT(names_, 14, mkChar("reorthogonalize_u"));
  SET_STRING_ELT(names_, 15, mkChar("reorthogonalize_v"));
  SET_STRING_ELT(names_, 16, mkChar("work"));
  SET_STRING_ELT(names_, 17, mkChar("retained"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return out_;
}

extern "C" SEXP eigencore_irlba_lbd_dense_retained(SEXP A_, SEXP initial_start_,
                                                   SEXP retained_right_,
                                                   SEXP retained_left_,
                                                   SEXP alpha_, SEXP beta_,
                                                   SEXP random_tails_,
                                                   SEXP work_, SEXP retained_,
                                                   SEXP max_restarts_,
                                                   SEXP rank_,
                                                   SEXP target_kind_,
                                                   SEXP tol_,
                                                   SEXP reorth_policy_) {
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue || LENGTH(dimA) != 2) {
    error("A must be a double matrix");
  }
  validate_irlba_lbd_retained_contract(
    INTEGER(dimA)[0], INTEGER(dimA)[1], initial_start_, retained_right_,
    retained_left_, alpha_, beta_, random_tails_, work_, retained_,
    max_restarts_, rank_, target_kind_, tol_, reorth_policy_
  );
  DenseColumnMajorOperator impl_holder = {
    INTEGER(dimA)[0], INTEGER(dimA)[1], REAL(A_)
  };
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_dense_apply;
  };
  return irlba_lbd_retained_impl(
    configure, INTEGER(dimA)[0], INTEGER(dimA)[1],
    REAL(initial_start_), REAL(retained_right_), REAL(retained_left_),
    REAL(alpha_), REAL(beta_), REAL(random_tails_),
    asInteger(work_), asInteger(retained_), asInteger(max_restarts_),
    asInteger(rank_), asInteger(target_kind_), asReal(tol_),
    asInteger(reorth_policy_)
  );
}

extern "C" SEXP eigencore_irlba_lbd_csc_retained(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_,
                                                 SEXP initial_start_,
                                                 SEXP retained_right_,
                                                 SEXP retained_left_,
                                                 SEXP alpha_, SEXP beta_,
                                                 SEXP random_tails_,
                                                 SEXP work_, SEXP retained_,
                                                 SEXP max_restarts_,
                                                 SEXP rank_,
                                                 SEXP target_kind_,
                                                 SEXP tol_,
                                                 SEXP reorth_policy_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      LENGTH(dim_) != 2) {
    error("invalid CSC retained IRLBA/LBD inputs");
  }
  validate_irlba_lbd_retained_contract(
    INTEGER(dim_)[0], INTEGER(dim_)[1], initial_start_, retained_right_,
    retained_left_, alpha_, beta_, random_tails_, work_, retained_,
    max_restarts_, rank_, target_kind_, tol_, reorth_policy_
  );
  CSCOperator impl_holder = {
    INTEGER(dim_)[0], INTEGER(dim_)[1], INTEGER(i_), INTEGER(p_), REAL(x_)
  };
  auto configure = [&impl_holder](void** impl, EigencoreApplyFn* apply) {
    *impl = &impl_holder;
    *apply = eigencore_csc_apply;
  };
  return irlba_lbd_retained_impl(
    configure, INTEGER(dim_)[0], INTEGER(dim_)[1],
    REAL(initial_start_), REAL(retained_right_), REAL(retained_left_),
    REAL(alpha_), REAL(beta_), REAL(random_tails_),
    asInteger(work_), asInteger(retained_), asInteger(max_restarts_),
    asInteger(rank_), asInteger(target_kind_), asReal(tol_),
    asInteger(reorth_policy_)
  );
}

static const R_CallMethodDef CallEntries[] = {
  {"eigencore_dense_block_apply", (DL_FUNC) &eigencore_dense_block_apply, 6},
  {"eigencore_csc_block_apply", (DL_FUNC) &eigencore_csc_block_apply, 9},
  {"eigencore_diagonal_block_apply", (DL_FUNC) &eigencore_diagonal_block_apply, 7},
  {"eigencore_native_apply_noalloc_check", (DL_FUNC) &eigencore_native_apply_noalloc_check, 4},
  {"eigencore_dense_apply_int_guard_check", (DL_FUNC) &eigencore_dense_apply_int_guard_check, 0},
  {"eigencore_col_norms", (DL_FUNC) &eigencore_col_norms, 1},
  {"eigencore_mgs2", (DL_FUNC) &eigencore_mgs2, 2},
  {"eigencore_cholqr2", (DL_FUNC) &eigencore_cholqr2, 1},
  {"eigencore_b_cholqr2", (DL_FUNC) &eigencore_b_cholqr2, 2},
  {"eigencore_diagonal_b_cholqr2", (DL_FUNC) &eigencore_diagonal_b_cholqr2, 3},
  {"eigencore_reorthogonalize_against", (DL_FUNC) &eigencore_reorthogonalize_against, 3},
  {"eigencore_basis_workspace_create", (DL_FUNC) &eigencore_basis_workspace_create, 3},
  {"eigencore_basis_workspace_info", (DL_FUNC) &eigencore_basis_workspace_info, 1},
  {"eigencore_reorthogonalize_against_workspace", (DL_FUNC) &eigencore_reorthogonalize_against_workspace, 4},
  {"eigencore_lanczos_dense", (DL_FUNC) &eigencore_lanczos_dense, 6},
  {"eigencore_lanczos_csc", (DL_FUNC) &eigencore_lanczos_csc, 9},
  {"eigencore_golub_kahan_dense", (DL_FUNC) &eigencore_golub_kahan_dense, 7},
  {"eigencore_golub_kahan_csc", (DL_FUNC) &eigencore_golub_kahan_csc, 10},
  {"eigencore_golub_kahan_dense_fit", (DL_FUNC) &eigencore_golub_kahan_dense_fit, 9},
  {"eigencore_golub_kahan_csc_fit", (DL_FUNC) &eigencore_golub_kahan_csc_fit, 12},
  {"eigencore_irlba_lbd_dense_retained", (DL_FUNC) &eigencore_irlba_lbd_dense_retained, 14},
  {"eigencore_irlba_lbd_csc_retained", (DL_FUNC) &eigencore_irlba_lbd_csc_retained, 17},
  {"eigencore_block_golub_kahan_dense_basis", (DL_FUNC) &eigencore_block_golub_kahan_dense_basis, 3},
  {"eigencore_block_golub_kahan_dense_basis_cached", (DL_FUNC) &eigencore_block_golub_kahan_dense_basis_cached, 4},
  {"eigencore_block_golub_kahan_csc_basis", (DL_FUNC) &eigencore_block_golub_kahan_csc_basis, 6},
  {"eigencore_block_golub_kahan_csc_basis_cached", (DL_FUNC) &eigencore_block_golub_kahan_csc_basis_cached, 7},
  {"eigencore_block_golub_kahan_dense_fit", (DL_FUNC) &eigencore_block_golub_kahan_dense_fit, 5},
  {"eigencore_block_golub_kahan_dense_fit_cached", (DL_FUNC) &eigencore_block_golub_kahan_dense_fit_cached, 6},
  {"eigencore_block_golub_kahan_csc_fit", (DL_FUNC) &eigencore_block_golub_kahan_csc_fit, 8},
  {"eigencore_block_golub_kahan_csc_fit_cached", (DL_FUNC) &eigencore_block_golub_kahan_csc_fit_cached, 9},
  {"eigencore_block_golub_kahan_dense_retained_cycle", (DL_FUNC) &eigencore_block_golub_kahan_dense_retained_cycle, 11},
  {"eigencore_block_golub_kahan_csc_retained_cycle", (DL_FUNC) &eigencore_block_golub_kahan_csc_retained_cycle, 14},
  {"eigencore_block_lanczos_dense", (DL_FUNC) &eigencore_block_lanczos_dense, 7},
  {"eigencore_block_lanczos_csc", (DL_FUNC) &eigencore_block_lanczos_csc, 10},
  {"eigencore_block_thick_restart_lanczos_dense", (DL_FUNC) &eigencore_block_thick_restart_lanczos_dense, 9},
  {"eigencore_block_thick_restart_lanczos_csc", (DL_FUNC) &eigencore_block_thick_restart_lanczos_csc, 12},
  {"eigencore_lobpcg_dense", (DL_FUNC) &eigencore_lobpcg_dense, 10},
  {"eigencore_lobpcg_dense_dense_b", (DL_FUNC) &eigencore_lobpcg_dense_dense_b, 11},
  {"eigencore_lobpcg_dense_diagonal_b", (DL_FUNC) &eigencore_lobpcg_dense_diagonal_b, 12},
  {"eigencore_lobpcg_dense_csc_b", (DL_FUNC) &eigencore_lobpcg_dense_csc_b, 14},
  {"eigencore_lobpcg_csc_diagonal_b", (DL_FUNC) &eigencore_lobpcg_csc_diagonal_b, 15},
  {"eigencore_lobpcg_csc_csc_b", (DL_FUNC) &eigencore_lobpcg_csc_csc_b, 17},
  {"eigencore_lobpcg_diagonal_diagonal_b", (DL_FUNC) &eigencore_lobpcg_diagonal_diagonal_b, 14},
  {"eigencore_lobpcg_dense_operator_b", (DL_FUNC) &eigencore_lobpcg_dense_operator_b, 11},
  {"eigencore_lobpcg_csc_operator_b", (DL_FUNC) &eigencore_lobpcg_csc_operator_b, 14},
  {"eigencore_lobpcg_diagonal_operator_b", (DL_FUNC) &eigencore_lobpcg_diagonal_operator_b, 13},
  {"eigencore_lobpcg_csc", (DL_FUNC) &eigencore_lobpcg_csc, 13},
  {"eigencore_lobpcg_csc_shifted_tridiagonal", (DL_FUNC) &eigencore_lobpcg_csc_shifted_tridiagonal, 10},
  {"eigencore_orthogonality_loss", (DL_FUNC) &eigencore_orthogonality_loss, 2},
  {"eigencore_dense_eigen_residuals", (DL_FUNC) &eigencore_dense_eigen_residuals, 4},
  {"eigencore_dense_eigen_certificate", (DL_FUNC) &eigencore_dense_eigen_certificate, 5},
  {"eigencore_dense_svd_residuals", (DL_FUNC) &eigencore_dense_svd_residuals, 4},
  {"eigencore_dense_svd_certificate", (DL_FUNC) &eigencore_dense_svd_certificate, 5},
  {"eigencore_dense_svd_certificate_cached_av", (DL_FUNC) &eigencore_dense_svd_certificate_cached_av, 6},
  {"eigencore_csc_eigen_certificate", (DL_FUNC) &eigencore_csc_eigen_certificate, 8},
  {"eigencore_diagonal_eigen_certificate", (DL_FUNC) &eigencore_diagonal_eigen_certificate, 7},
  {"eigencore_csc_svd_certificate", (DL_FUNC) &eigencore_csc_svd_certificate, 9},
  {"eigencore_csc_svd_certificate_cached_av", (DL_FUNC) &eigencore_csc_svd_certificate_cached_av, 10},
  {"eigencore_diagonal_svd_certificate", (DL_FUNC) &eigencore_diagonal_svd_certificate, 8},
  {"eigencore_diagonal_svd_certificate_cached_av", (DL_FUNC) &eigencore_diagonal_svd_certificate_cached_av, 9},
  {"eigencore_rayleigh_ritz_symmetric", (DL_FUNC) &eigencore_rayleigh_ritz_symmetric, 2},
  {"eigencore_tridiagonal_eigen", (DL_FUNC) &eigencore_tridiagonal_eigen, 2},
  {"eigencore_bidiagonal_svd", (DL_FUNC) &eigencore_bidiagonal_svd, 2},
  {"eigencore_block_golub_kahan_ritz", (DL_FUNC) &eigencore_block_golub_kahan_ritz, 5},
  {"eigencore_golub_kahan_ritz", (DL_FUNC) &eigencore_golub_kahan_ritz, 7},
  {"eigencore_dense_is_symmetric", (DL_FUNC) &eigencore_dense_is_symmetric, 2},
  {"eigencore_dense_symmetric_eigen", (DL_FUNC) &eigencore_dense_symmetric_eigen, 1},
  {"eigencore_dense_symmetric_eigen_dsyevd", (DL_FUNC) &eigencore_dense_symmetric_eigen_dsyevd, 1},
  {"eigencore_dense_symmetric_eigen_selected", (DL_FUNC) &eigencore_dense_symmetric_eigen_selected, 3},
  {"eigencore_dense_symmetric_eigen_dsyevx_selected", (DL_FUNC) &eigencore_dense_symmetric_eigen_dsyevx_selected, 3},
  {"eigencore_csc_left_gram_svd", (DL_FUNC) &eigencore_csc_left_gram_svd, 6},
  {"eigencore_dense_generalized_spd_eigen", (DL_FUNC) &eigencore_dense_generalized_spd_eigen, 2},
  {"eigencore_dense_svd", (DL_FUNC) &eigencore_dense_svd, 1},
  {"eigencore_tridiagonal_solve", (DL_FUNC) &eigencore_tridiagonal_solve, 4},
  {NULL, NULL, 0}
};

extern "C" void R_init_eigencore(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
