#include <cstring>
#include <cmath>
#include <climits>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"
#include <R_ext/Random.h>
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

static Rcomplex scalar_as_rcomplex(SEXP x, const char* name) {
  if (LENGTH(x) < 1) {
    error("%s must be a numeric or complex scalar", name);
  }
  Rcomplex out;
  switch (TYPEOF(x)) {
    case CPLXSXP:
      return COMPLEX(x)[0];
    case REALSXP:
      out.r = REAL(x)[0];
      out.i = 0.0;
      return out;
    case INTSXP:
      out.r = static_cast<double>(INTEGER(x)[0]);
      out.i = 0.0;
      return out;
    default:
      error("%s must be a numeric or complex scalar", name);
  }
  out.r = 0.0;
  out.i = 0.0;
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

extern "C" int eigencore_dense_complex_apply(void* impl,
                                             EigencoreTranspose op,
                                             int64_t block_cols,
                                             const Rcomplex* X,
                                             int64_t ldx,
                                             Rcomplex alpha,
                                             Rcomplex beta,
                                             Rcomplex* Y,
                                             int64_t ldy,
                                             EigencoreWorkspace* workspace) {
  (void) workspace;
  DenseComplexColumnMajorOperator* dense =
    static_cast<DenseComplexColumnMajorOperator*>(impl);
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

  const char transa = (op == EIGENCORE_TRANSPOSE_ADJOINT) ? 'C' : 'N';
  const char transb = 'N';
  const int out_rows = static_cast<int>(out_rows64);
  const int block_cols_i = static_cast<int>(block_cols);
  const int inner = static_cast<int>(inner64);
  const int lda = static_cast<int>(dense->rows);
  const int ldb = static_cast<int>(ldx);
  const int ldc = static_cast<int>(ldy);

  F77_CALL(zgemm)(&transa, &transb, &out_rows, &block_cols_i, &inner,
                  &alpha, dense->values, &lda, const_cast<Rcomplex*>(X), &ldb,
                  &beta, Y, &ldc FCONE FCONE);
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
      // Wide blocks: sweep the matrix once per chunk of at most 10 columns so
      // each pass keeps the per-column accumulators and X/Y panels hot,
      // instead of striding across the whole block per nonzero.
      for (int64_t chunk = 0; chunk < block_cols; chunk += 10) {
        const int64_t chunk_cols =
          (block_cols - chunk < 10) ? block_cols - chunk : 10;
        double* yptr[10];
        const double* xptr[10];
        double xval[10];
        for (int64_t block = 0; block < chunk_cols; ++block) {
          yptr[block] = Y + (chunk + block) * ldy;
          xptr[block] = X + (chunk + block) * ldx;
        }
        for (int64_t col = 0; col < csc->cols; ++col) {
          bool all_zero = true;
          for (int64_t block = 0; block < chunk_cols; ++block) {
            xval[block] = xptr[block][col];
            all_zero = all_zero && xval[block] == 0.0;
          }
          if (all_zero) continue;
          for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
            const int row = csc->row_idx[pos];
            const double a = alpha * csc->values[pos];
            for (int64_t block = 0; block < chunk_cols; ++block) {
              yptr[block][row] += a * xval[block];
            }
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
      // Wide blocks: same chunking as the forward path, keeping at most 10
      // running dot-product accumulators per matrix sweep.
      for (int64_t chunk = 0; chunk < block_cols; chunk += 10) {
        const int64_t chunk_cols =
          (block_cols - chunk < 10) ? block_cols - chunk : 10;
        double* yptr[10];
        const double* xptr[10];
        double acc[10];
        for (int64_t block = 0; block < chunk_cols; ++block) {
          yptr[block] = Y + (chunk + block) * ldy;
          xptr[block] = X + (chunk + block) * ldx;
        }
        for (int64_t col = 0; col < csc->cols; ++col) {
          for (int64_t block = 0; block < chunk_cols; ++block) {
            acc[block] = 0.0;
          }
          for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
            const int row = csc->row_idx[pos];
            const double a = csc->values[pos];
            for (int64_t block = 0; block < chunk_cols; ++block) {
              acc[block] += a * xptr[block][row];
            }
          }
          for (int64_t block = 0; block < chunk_cols; ++block) {
            yptr[block][col] += alpha * acc[block];
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

extern "C" int eigencore_normal_equations_apply(void* impl,
                                                 EigencoreTranspose op,
                                                 int64_t block_cols,
                                                 const double* X,
                                                 int64_t ldx,
                                                 double alpha,
                                                 double beta,
                                                 double* Y,
                                                 int64_t ldy,
                                                 EigencoreWorkspace* workspace) {
  (void) op;  // A^T A and A A^T are symmetric
  NormalEquationsOperator* normal = static_cast<NormalEquationsOperator*>(impl);
  const int64_t outer = (normal->side == 0) ? normal->cols : normal->rows;
  const int64_t inner = (normal->side == 0) ? normal->rows : normal->cols;
  if (ldx < outer || ldy < outer) {
    return -1;
  }
  if (block_cols > normal->scratch_block_capacity || normal->scratch == nullptr) {
    return -1;
  }
  const EigencoreTranspose first =
    (normal->side == 0) ? EIGENCORE_TRANSPOSE_NONE : EIGENCORE_TRANSPOSE_ADJOINT;
  const EigencoreTranspose second =
    (normal->side == 0) ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE;
  int status = normal->base_apply(normal->base_impl, first, block_cols,
                                  X, ldx, 1.0, 0.0,
                                  normal->scratch, inner, workspace);
  if (status != 0) {
    return status;
  }
  status = normal->base_apply(normal->base_impl, second, block_cols,
                              normal->scratch, inner, alpha, beta,
                              Y, ldy, workspace);
  return status;
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
  if (op != EIGENCORE_TRANSPOSE_NONE && op != EIGENCORE_TRANSPOSE_ADJOINT) {
    return -1;
  }
  RApplyOperator* fn = static_cast<RApplyOperator*>(impl);
  if (!eigencore_int_indexable(fn->rows) ||
      !eigencore_int_indexable(fn->cols) ||
      !eigencore_int_indexable(block_cols) ||
      !eigencore_int_indexable(ldx) ||
      !eigencore_int_indexable(ldy)) {
    return -2;
  }
  const bool adjoint = (op == EIGENCORE_TRANSPOSE_ADJOINT);
  SEXP closure = adjoint ? fn->apply_adjoint : fn->apply;
  const int in_rows = static_cast<int>(adjoint ? fn->rows : fn->cols);
  const int out_rows = static_cast<int>(adjoint ? fn->cols : fn->rows);
  const int cols = static_cast<int>(block_cols);
  if (ldx < in_rows || ldy < out_rows || cols < 1 || TYPEOF(closure) != CLOSXP) {
    return -1;
  }

  SEXP X_ = PROTECT(allocMatrix(REALSXP, in_rows, cols));
  SEXP Y_ = PROTECT(allocMatrix(REALSXP, out_rows, cols));
  for (int col = 0; col < cols; ++col) {
    const double* x_col = X + static_cast<int64_t>(col) * ldx;
    double* x_dst = REAL(X_) + static_cast<int64_t>(col) * in_rows;
    const double* y_col = Y + static_cast<int64_t>(col) * ldy;
    double* y_dst = REAL(Y_) + static_cast<int64_t>(col) * out_rows;
    std::memcpy(x_dst, x_col, sizeof(double) * static_cast<size_t>(in_rows));
    std::memcpy(y_dst, y_col, sizeof(double) * static_cast<size_t>(out_rows));
  }
  SEXP alpha_ = PROTECT(ScalarReal(alpha));
  SEXP beta_ = PROTECT(ScalarReal(beta));
  SEXP call = PROTECT(lang5(closure, X_, alpha_, beta_, Y_));
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
      INTEGER(dimY)[0] != out_rows || INTEGER(dimY)[1] != cols) {
    UNPROTECT(6);
    return -8;
  }
  for (int col = 0; col < cols; ++col) {
    const double* out_col = REAL(out_) + static_cast<int64_t>(col) * out_rows;
    double* y_col = Y + static_cast<int64_t>(col) * ldy;
    std::memcpy(y_col, out_col, sizeof(double) * static_cast<size_t>(out_rows));
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

extern "C" SEXP eigencore_dense_complex_block_apply(SEXP A_, SEXP X_, SEXP alpha_,
                                                    SEXP beta_, SEXP Y_,
                                                    SEXP adjoint_) {
  if (!isComplex(A_) || !isComplex(X_) || !isComplex(Y_)) {
    error("A, X, and Y must be complex matrices");
  }
  if (!isLogical(adjoint_) || LENGTH(adjoint_) != 1) {
    error("adjoint must be a logical scalar");
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
  const bool adjoint = LOGICAL(adjoint_)[0];

  const int inner = adjoint ? m : n;
  const int out_rows = adjoint ? n : m;
  if (xr != inner) {
    error("non-conformable X for dense complex block apply");
  }
  if (yr != out_rows || yc != xc) {
    error("non-conformable Y for dense complex block apply");
  }

  SEXP out_ = PROTECT(duplicate(Y_));
  DenseComplexColumnMajorOperator impl = {m, n, COMPLEX(A_)};
  const int status = eigencore_dense_complex_apply(
    &impl,
    adjoint ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    COMPLEX(X_),
    xr,
    scalar_as_rcomplex(alpha_, "alpha"),
    scalar_as_rcomplex(beta_, "beta"),
    COMPLEX(out_),
    out_rows,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error("dense complex block apply", status);
  }

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_dense_randomized_apply(SEXP A_, SEXP X_,
                                                 SEXP transpose_) {
  if (!isReal(A_) || !isReal(X_) || !isLogical(transpose_)) {
    error("invalid dense randomized apply inputs");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  if (dimA == R_NilValue || dimX == R_NilValue) {
    error("A and X must be matrices");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;
  if (xr != inner) {
    error("non-conformable X for dense randomized apply");
  }

  SEXP out_ = PROTECT(allocMatrix(REALSXP, out_rows, xc));
  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  const int status = eigencore_dense_apply(
    &impl,
    transpose ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    REAL(X_),
    xr,
    1.0,
    0.0,
    REAL(out_),
    out_rows,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error("dense randomized apply", status);
  }

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_dense_randomized_sketch(SEXP A_, SEXP cols_) {
  if (!isReal(A_)) {
    error("invalid dense randomized sketch inputs");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int sketch_cols = asInteger(cols_);
  if (sketch_cols == NA_INTEGER || sketch_cols < 0) {
    error("sketch column count must be non-negative");
  }

  SEXP out_ = PROTECT(allocMatrix(REALSXP, m, sketch_cols));
  if (sketch_cols == 0) {
    UNPROTECT(1);
    return out_;
  }

  double* omega = reinterpret_cast<double*>(
    R_alloc(static_cast<size_t>(n) * static_cast<size_t>(sketch_cols),
            sizeof(double))
  );
  GetRNGstate();
  for (int64_t pos = 0;
       pos < static_cast<int64_t>(n) * static_cast<int64_t>(sketch_cols);
       ++pos) {
    omega[pos] = norm_rand();
  }
  PutRNGstate();

  const char trans_a = 'N';
  const char trans_omega = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;
  F77_CALL(dgemm)(&trans_a, &trans_omega, &m, &sketch_cols, &n,
                  &alpha, REAL(A_), &m, omega, &n,
                  &beta, REAL(out_), &m FCONE FCONE);

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_dense_randomized_project_transposed(SEXP A_,
                                                              SEXP Q_) {
  if (!isReal(A_) || !isReal(Q_)) {
    error("invalid dense randomized projection inputs");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimQ = getAttrib(Q_, R_DimSymbol);
  if (dimA == R_NilValue || dimQ == R_NilValue) {
    error("A and Q must be matrices");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int qr = INTEGER(dimQ)[0];
  const int qcols = INTEGER(dimQ)[1];
  if (qr != m) {
    error("non-conformable Q for dense randomized projection");
  }

  SEXP out_ = PROTECT(allocMatrix(REALSXP, qcols, n));
  const char trans_q = 'T';
  const char trans_a = 'N';
  const double alpha = 1.0;
  const double beta = 0.0;
  F77_CALL(dgemm)(&trans_q, &trans_a, &qcols, &n, &m,
                  &alpha, REAL(Q_), &m, REAL(A_), &m,
                  &beta, REAL(out_), &qcols FCONE FCONE);
  SEXP transposed_ = PROTECT(ScalarLogical(TRUE));
  setAttrib(out_, install("transposed"), transposed_);
  UNPROTECT(2);
  return out_;
}

static void dense_randomized_thin_qr(std::vector<double>& X, int rows, int cols) {
  if (rows < cols) {
    error("native randomized QR requires rows >= columns");
  }
  if (cols == 0) {
    return;
  }

  std::vector<double> tau(static_cast<size_t>(cols));
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgeqrf)(&rows, &cols, X.data(), &rows, tau.data(),
                   &work_query, &lwork, &info);
  if (info != 0) {
    error("native randomized QR workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  std::vector<double> work(static_cast<size_t>(lwork));
  F77_CALL(dgeqrf)(&rows, &cols, X.data(), &rows, tau.data(),
                   work.data(), &lwork, &info);
  if (info != 0) {
    error("native randomized QR failed with info=%d", info);
  }

  lwork = -1;
  work_query = 0.0;
  F77_CALL(dorgqr)(&rows, &cols, &cols, X.data(), &rows, tau.data(),
                   &work_query, &lwork, &info);
  if (info != 0) {
    error("native randomized Q formation workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  work.assign(static_cast<size_t>(lwork), 0.0);
  F77_CALL(dorgqr)(&rows, &cols, &cols, X.data(), &rows, tau.data(),
                   work.data(), &lwork, &info);
  if (info != 0) {
    error("native randomized Q formation failed with info=%d", info);
  }
}

struct DenseRandomizedCertificate {
  std::vector<double> left;
  std::vector<double> right;
  std::vector<double> combined;
  std::vector<double> backward;
  std::vector<int> converged;
  double orth_u = 0.0;
  double orth_v = 0.0;
  double scale = 0.0;
  bool passed = false;
};

static double dense_randomized_frobenius_norm(const double* A, int len) {
  double sum = 0.0;
  for (int pos = 0; pos < len; ++pos) {
    sum += A[pos] * A[pos];
  }
  return std::sqrt(sum);
}

static double dense_randomized_column_norm(const std::vector<double>& X,
                                           int rows,
                                           int col) {
  const double* x = X.data() + static_cast<int64_t>(col) * rows;
  double sum = 0.0;
  for (int row = 0; row < rows; ++row) {
    sum += x[row] * x[row];
  }
  return std::sqrt(sum);
}

static double dense_randomized_max_orthogonality(const std::vector<double>& X,
                                                 int rows,
                                                 int cols) {
  if (cols == 0) {
    return 0.0;
  }
  std::vector<double> gram(static_cast<size_t>(cols) * static_cast<size_t>(cols), 0.0);
  // The Gram matrix is symmetric: dsyrk computes the upper triangle in half
  // the flops of the previous dgemm full-product formulation.
  const char uplo = 'U';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dsyrk)(&uplo, &trans, &cols, &rows,
                  &one, const_cast<double*>(X.data()), &rows,
                  &zero, gram.data(), &cols FCONE FCONE);
  double out = 0.0;
  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row <= col; ++row) {
      const double expected = (row == col) ? 1.0 : 0.0;
      const double loss = std::fabs(gram[row + static_cast<int64_t>(col) * cols] - expected);
      if (loss > out) {
        out = loss;
      }
    }
  }
  return out;
}

static DenseRandomizedCertificate dense_randomized_certificate(
    const double* A,
    int m,
    int n,
    const std::vector<double>& d,
    const std::vector<double>& U,
    const std::vector<double>& V,
    double tol) {
  const int rank = static_cast<int>(d.size());
  DenseRandomizedCertificate cert;
  cert.left.assign(static_cast<size_t>(rank), 0.0);
  cert.right.assign(static_cast<size_t>(rank), 0.0);
  cert.combined.assign(static_cast<size_t>(rank), 0.0);
  cert.backward.assign(static_cast<size_t>(rank), 0.0);
  cert.converged.assign(static_cast<size_t>(rank), 0);
  cert.scale = dense_randomized_frobenius_norm(A, m * n);
  if (cert.scale < 2.2204460492503131e-16) {
    cert.scale = 2.2204460492503131e-16;
  }

  std::vector<double> AV(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  std::vector<double> ATU(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  const char notrans = 'N';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  if (rank > 0) {
    F77_CALL(dgemm)(&notrans, &notrans, &m, &rank, &n,
                    &one, const_cast<double*>(A), &m,
                    const_cast<double*>(V.data()), &n,
                    &zero, AV.data(), &m FCONE FCONE);
    F77_CALL(dgemm)(&trans, &notrans, &n, &rank, &m,
                    &one, const_cast<double*>(A), &m,
                    const_cast<double*>(U.data()), &m,
                    &zero, ATU.data(), &n FCONE FCONE);
  }

  bool all_converged = true;
  for (int col = 0; col < rank; ++col) {
    double* av = AV.data() + static_cast<int64_t>(col) * m;
    double* atu = ATU.data() + static_cast<int64_t>(col) * n;
    const double* u = U.data() + static_cast<int64_t>(col) * m;
    const double* v = V.data() + static_cast<int64_t>(col) * n;
    for (int row = 0; row < m; ++row) {
      av[row] -= d[static_cast<size_t>(col)] * u[row];
    }
    for (int row = 0; row < n; ++row) {
      atu[row] -= d[static_cast<size_t>(col)] * v[row];
    }
    cert.left[static_cast<size_t>(col)] = dense_randomized_column_norm(AV, m, col);
    cert.right[static_cast<size_t>(col)] = dense_randomized_column_norm(ATU, n, col);
    cert.combined[static_cast<size_t>(col)] = std::sqrt(
      cert.left[static_cast<size_t>(col)] * cert.left[static_cast<size_t>(col)] +
      cert.right[static_cast<size_t>(col)] * cert.right[static_cast<size_t>(col)]
    );
    cert.backward[static_cast<size_t>(col)] =
      cert.combined[static_cast<size_t>(col)] / cert.scale;
    cert.converged[static_cast<size_t>(col)] =
      cert.backward[static_cast<size_t>(col)] <= tol ? 1 : 0;
    all_converged = all_converged && cert.converged[static_cast<size_t>(col)];
  }
  cert.orth_u = dense_randomized_max_orthogonality(U, m, rank);
  cert.orth_v = dense_randomized_max_orthogonality(V, n, rank);
  const double orth_tol = tol > std::sqrt(2.2204460492503131e-16)
    ? tol
    : std::sqrt(2.2204460492503131e-16);
  cert.passed = all_converged && cert.orth_u <= orth_tol && cert.orth_v <= orth_tol;
  return cert;
}

struct DenseRandomizedCandidate {
  std::vector<double> d;
  std::vector<double> U;
  std::vector<double> V;
  DenseRandomizedCertificate certificate;
};

static DenseRandomizedCandidate dense_randomized_candidate(
    const double* A,
    int m,
    int n,
    int rank,
    const std::vector<double>& Q,
    int q_cols,
    double tol,
    double* small_svd_seconds,
    double* vector_seconds,
    double* certificate_seconds) {
  DenseRandomizedCandidate candidate;
  const char trans = 'T';
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;

  std::vector<double> B(static_cast<size_t>(q_cols) * static_cast<size_t>(n), 0.0);
  F77_CALL(dgemm)(&trans, &notrans, &q_cols, &n, &m,
                  &one, const_cast<double*>(Q.data()), &m,
                  const_cast<double*>(A), &m,
                  &zero, B.data(), &q_cols FCONE FCONE);

  auto t0 = native_timer_now();
  // dgesvd destroys its input; B is not read again, so hand it over directly.
  std::vector<double>& work_B = B;
  std::vector<double> d_all(static_cast<size_t>(q_cols), 0.0);
  std::vector<double> U_small(static_cast<size_t>(q_cols) * static_cast<size_t>(q_cols), 0.0);
  std::vector<double> VT(static_cast<size_t>(q_cols) * static_cast<size_t>(n), 0.0);
  char jobu = 'S';
  char jobvt = 'S';
  int lda = q_cols;
  int ldu = q_cols;
  int ldvt = q_cols;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgesvd)(&jobu, &jobvt, &q_cols, &n, work_B.data(), &lda,
                   d_all.data(), U_small.data(), &ldu, VT.data(), &ldvt,
                   &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("native randomized projected SVD workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  std::vector<double> work(static_cast<size_t>(lwork), 0.0);
  F77_CALL(dgesvd)(&jobu, &jobvt, &q_cols, &n, work_B.data(), &lda,
                   d_all.data(), U_small.data(), &ldu, VT.data(), &ldvt,
                   work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("native randomized projected SVD failed with info=%d", info);
  }
  *small_svd_seconds += native_timer_elapsed(t0);

  t0 = native_timer_now();
  candidate.d.assign(d_all.begin(), d_all.begin() + rank);
  candidate.U.assign(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  candidate.V.assign(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  F77_CALL(dgemm)(&notrans, &notrans, &m, &rank, &q_cols,
                  &one, const_cast<double*>(Q.data()), &m,
                  U_small.data(), &q_cols,
                  &zero, candidate.U.data(), &m FCONE FCONE);
  for (int col = 0; col < rank; ++col) {
    for (int row = 0; row < n; ++row) {
      candidate.V[row + static_cast<int64_t>(col) * n] =
        VT[col + static_cast<int64_t>(row) * q_cols];
    }
  }
  *vector_seconds += native_timer_elapsed(t0);

  t0 = native_timer_now();
  candidate.certificate = dense_randomized_certificate(
    A, m, n, candidate.d, candidate.U, candidate.V, tol
  );
  *certificate_seconds += native_timer_elapsed(t0);
  return candidate;
}

static double csc_randomized_frobenius_norm(const double* values, int nnz) {
  double sum = 0.0;
  for (int pos = 0; pos < nnz; ++pos) {
    sum += values[pos] * values[pos];
  }
  return std::sqrt(sum);
}

static void csc_randomized_apply_block(const CSCOperator& impl,
                                       EigencoreTranspose transpose,
                                       int block_cols,
                                       const double* X,
                                       int ldx,
                                       double* Y,
                                       int ldy,
                                       const char* label) {
  const int status = eigencore_csc_apply(
    const_cast<CSCOperator*>(&impl),
    transpose,
    block_cols,
    X,
    ldx,
    1.0,
    0.0,
    Y,
    ldy,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error(label, status);
  }
}

static void csc_randomized_project_transposed(const CSCOperator& impl,
                                              const std::vector<double>& Q,
                                              int q_cols,
                                              std::vector<double>& B) {
  const int m = impl.rows;
  const int n = impl.cols;
  std::fill(B.begin(), B.end(), 0.0);
  // Transpose Q (m x q_cols, column-major) into row-major Qt so the per-nonzero
  // inner loop reads a contiguous q_cols-length panel instead of striding m
  // doubles per element across Q's columns.
  std::vector<double> Qt(static_cast<size_t>(m) * static_cast<size_t>(q_cols));
  for (int block = 0; block < q_cols; ++block) {
    const double* q_col = Q.data() + static_cast<int64_t>(block) * m;
    for (int row = 0; row < m; ++row) {
      Qt[static_cast<int64_t>(row) * q_cols + block] = q_col[row];
    }
  }
  for (int col = 0; col < n; ++col) {
    double* out_col = B.data() + static_cast<int64_t>(col) * q_cols;
    for (int pos = impl.col_ptr[col]; pos < impl.col_ptr[col + 1]; ++pos) {
      const int row = impl.row_idx[pos];
      const double value = impl.values[pos];
      const double* qt_row = Qt.data() + static_cast<int64_t>(row) * q_cols;
      for (int block = 0; block < q_cols; ++block) {
        out_col[block] += value * qt_row[block];
      }
    }
  }
}

static DenseRandomizedCertificate csc_randomized_certificate(
    const CSCOperator& impl,
    int nnz,
    const std::vector<double>& d,
    const std::vector<double>& U,
    const std::vector<double>& V,
    double tol) {
  const int m = impl.rows;
  const int n = impl.cols;
  const int rank = static_cast<int>(d.size());
  DenseRandomizedCertificate cert;
  cert.left.assign(static_cast<size_t>(rank), 0.0);
  cert.right.assign(static_cast<size_t>(rank), 0.0);
  cert.combined.assign(static_cast<size_t>(rank), 0.0);
  cert.backward.assign(static_cast<size_t>(rank), 0.0);
  cert.converged.assign(static_cast<size_t>(rank), 0);
  cert.scale = csc_randomized_frobenius_norm(impl.values, nnz);
  if (cert.scale < 2.2204460492503131e-16) {
    cert.scale = 2.2204460492503131e-16;
  }

  std::vector<double> AV(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  std::vector<double> ATU(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  if (rank > 0) {
    csc_randomized_apply_block(
      impl, EIGENCORE_TRANSPOSE_NONE, rank, V.data(), n, AV.data(), m,
      "CSC randomized certificate"
    );
    csc_randomized_apply_block(
      impl, EIGENCORE_TRANSPOSE_ADJOINT, rank, U.data(), m, ATU.data(), n,
      "CSC randomized certificate"
    );
  }

  bool all_converged = true;
  for (int col = 0; col < rank; ++col) {
    double* av = AV.data() + static_cast<int64_t>(col) * m;
    double* atu = ATU.data() + static_cast<int64_t>(col) * n;
    const double* u = U.data() + static_cast<int64_t>(col) * m;
    const double* v = V.data() + static_cast<int64_t>(col) * n;
    for (int row = 0; row < m; ++row) {
      av[row] -= d[static_cast<size_t>(col)] * u[row];
    }
    for (int row = 0; row < n; ++row) {
      atu[row] -= d[static_cast<size_t>(col)] * v[row];
    }
    cert.left[static_cast<size_t>(col)] = dense_randomized_column_norm(AV, m, col);
    cert.right[static_cast<size_t>(col)] = dense_randomized_column_norm(ATU, n, col);
    cert.combined[static_cast<size_t>(col)] = std::sqrt(
      cert.left[static_cast<size_t>(col)] * cert.left[static_cast<size_t>(col)] +
      cert.right[static_cast<size_t>(col)] * cert.right[static_cast<size_t>(col)]
    );
    cert.backward[static_cast<size_t>(col)] =
      cert.combined[static_cast<size_t>(col)] / cert.scale;
    cert.converged[static_cast<size_t>(col)] =
      cert.backward[static_cast<size_t>(col)] <= tol ? 1 : 0;
    all_converged = all_converged && cert.converged[static_cast<size_t>(col)];
  }
  cert.orth_u = dense_randomized_max_orthogonality(U, m, rank);
  cert.orth_v = dense_randomized_max_orthogonality(V, n, rank);
  const double orth_tol = tol > std::sqrt(2.2204460492503131e-16)
    ? tol
    : std::sqrt(2.2204460492503131e-16);
  cert.passed = all_converged && cert.orth_u <= orth_tol && cert.orth_v <= orth_tol;
  return cert;
}

static DenseRandomizedCandidate csc_randomized_candidate(
    const CSCOperator& impl,
    int nnz,
    int rank,
    const std::vector<double>& Q,
    int q_cols,
    double tol,
    double* small_svd_seconds,
    double* vector_seconds,
    double* certificate_seconds) {
  const int m = impl.rows;
  const int n = impl.cols;
  DenseRandomizedCandidate candidate;
  const char notrans = 'N';
  const double one = 1.0;
  const double zero = 0.0;

  std::vector<double> B(static_cast<size_t>(q_cols) * static_cast<size_t>(n), 0.0);
  csc_randomized_project_transposed(impl, Q, q_cols, B);

  auto t0 = native_timer_now();
  // dgesvd destroys its input; B is not read again, so hand it over directly.
  std::vector<double>& work_B = B;
  std::vector<double> d_all(static_cast<size_t>(q_cols), 0.0);
  std::vector<double> U_small(static_cast<size_t>(q_cols) * static_cast<size_t>(q_cols), 0.0);
  std::vector<double> VT(static_cast<size_t>(q_cols) * static_cast<size_t>(n), 0.0);
  char jobu = 'S';
  char jobvt = 'S';
  int lda = q_cols;
  int ldu = q_cols;
  int ldvt = q_cols;
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgesvd)(&jobu, &jobvt, &q_cols, &n, work_B.data(), &lda,
                   d_all.data(), U_small.data(), &ldu, VT.data(), &ldvt,
                   &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("native CSC randomized projected SVD workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  std::vector<double> work(static_cast<size_t>(lwork), 0.0);
  F77_CALL(dgesvd)(&jobu, &jobvt, &q_cols, &n, work_B.data(), &lda,
                   d_all.data(), U_small.data(), &ldu, VT.data(), &ldvt,
                   work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("native CSC randomized projected SVD failed with info=%d", info);
  }
  *small_svd_seconds += native_timer_elapsed(t0);

  t0 = native_timer_now();
  candidate.d.assign(d_all.begin(), d_all.begin() + rank);
  candidate.U.assign(static_cast<size_t>(m) * static_cast<size_t>(rank), 0.0);
  candidate.V.assign(static_cast<size_t>(n) * static_cast<size_t>(rank), 0.0);
  F77_CALL(dgemm)(&notrans, &notrans, &m, &rank, &q_cols,
                  &one, const_cast<double*>(Q.data()), &m,
                  U_small.data(), &q_cols,
                  &zero, candidate.U.data(), &m FCONE FCONE);
  for (int col = 0; col < rank; ++col) {
    for (int row = 0; row < n; ++row) {
      candidate.V[row + static_cast<int64_t>(col) * n] =
        VT[col + static_cast<int64_t>(row) * q_cols];
    }
  }
  *vector_seconds += native_timer_elapsed(t0);

  t0 = native_timer_now();
  candidate.certificate = csc_randomized_certificate(
    impl, nnz, candidate.d, candidate.U, candidate.V, tol
  );
  *certificate_seconds += native_timer_elapsed(t0);
  return candidate;
}

static SEXP dense_randomized_certificate_pack(const DenseRandomizedCertificate& cert) {
  const int rank = static_cast<int>(cert.backward.size());
  SEXP left_ = PROTECT(allocVector(REALSXP, rank));
  SEXP right_ = PROTECT(allocVector(REALSXP, rank));
  SEXP combined_ = PROTECT(allocVector(REALSXP, rank));
  SEXP backward_ = PROTECT(allocVector(REALSXP, rank));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, rank));
  for (int idx = 0; idx < rank; ++idx) {
    REAL(left_)[idx] = cert.left[static_cast<size_t>(idx)];
    REAL(right_)[idx] = cert.right[static_cast<size_t>(idx)];
    REAL(combined_)[idx] = cert.combined[static_cast<size_t>(idx)];
    REAL(backward_)[idx] = cert.backward[static_cast<size_t>(idx)];
    LOGICAL(converged_)[idx] = cert.converged[static_cast<size_t>(idx)];
  }

  SEXP residuals_ = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(residuals_, 0, left_);
  SET_VECTOR_ELT(residuals_, 1, right_);
  SET_VECTOR_ELT(residuals_, 2, combined_);
  SEXP residual_names_ = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(residual_names_, 0, mkChar("left"));
  SET_STRING_ELT(residual_names_, 1, mkChar("right"));
  SET_STRING_ELT(residual_names_, 2, mkChar("combined"));
  setAttrib(residuals_, R_NamesSymbol, residual_names_);

  SEXP orth_ = PROTECT(allocVector(REALSXP, 2));
  REAL(orth_)[0] = cert.orth_u;
  REAL(orth_)[1] = cert.orth_v;
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SEXP scale_ = PROTECT(ScalarReal(cert.scale));
  SEXP passed_ = PROTECT(ScalarLogical(cert.passed));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, orth_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, scale_);
  SET_VECTOR_ELT(out_, 5, passed_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("scale"));
  SET_STRING_ELT(names_, 5, mkChar("passed"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(13);
  return out_;
}

extern "C" SEXP eigencore_dense_randomized_svd_controller(
    SEXP A_, SEXP rank_, SEXP oversample_, SEXP n_iter_, SEXP normalizer_,
    SEXP tol_) {
  if (!isReal(A_)) {
    error("A must be a double matrix");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  if (!isString(normalizer_) || LENGTH(normalizer_) < 1 ||
      std::strcmp(CHAR(STRING_ELT(normalizer_, 0)), "qr") != 0) {
    error("native dense randomized controller currently supports only QR normalization");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int limit = (m < n) ? m : n;
  int rank = asInteger(rank_);
  int oversample = asInteger(oversample_);
  int n_iter = asInteger(n_iter_);
  const double tol = asReal(tol_);
  if (rank == NA_INTEGER || rank < 1) {
    error("rank must be a positive integer");
  }
  if (oversample == NA_INTEGER || oversample < 0) {
    error("oversample must be a non-negative integer");
  }
  if (n_iter == NA_INTEGER || n_iter < 0) {
    error("n_iter must be a non-negative integer");
  }
  if (rank > limit) {
    rank = limit;
  }
  const int q_cols = (rank + oversample < limit) ? rank + oversample : limit;
  const double* A = REAL(A_);
  std::vector<double> stage(7, 0.0);
  const int stage_random = 0;
  const int stage_apply = 1;
  const int stage_normalize = 2;
  const int stage_small_svd = 3;
  const int stage_vector_form = 4;
  const int stage_certificate = 5;
  const int stage_controller = 6;
  auto controller_t0 = native_timer_now();

  auto t0 = native_timer_now();
  std::vector<double> omega(static_cast<size_t>(n) * static_cast<size_t>(q_cols), 0.0);
  GetRNGstate();
  for (int64_t pos = 0;
       pos < static_cast<int64_t>(n) * static_cast<int64_t>(q_cols);
       ++pos) {
    omega[static_cast<size_t>(pos)] = norm_rand();
  }
  PutRNGstate();
  stage[stage_random] += native_timer_elapsed(t0);

  t0 = native_timer_now();
  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(q_cols), 0.0);
  const char notrans = 'N';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&notrans, &notrans, &m, &q_cols, &n,
                  &one, const_cast<double*>(A), &m, omega.data(), &n,
                  &zero, Q.data(), &m FCONE FCONE);
  stage[stage_apply] += native_timer_elapsed(t0);

  t0 = native_timer_now();
  dense_randomized_thin_qr(Q, m, q_cols);
  stage[stage_normalize] += native_timer_elapsed(t0);

  int matvecs = 1;
  DenseRandomizedCandidate candidate = dense_randomized_candidate(
    A, m, n, rank, Q, q_cols, tol,
    &stage[stage_small_svd], &stage[stage_vector_form], &stage[stage_certificate]
  );
  matvecs += 1;
  DenseRandomizedCertificate initial_cert = candidate.certificate;
  bool early_stop_used = false;
  int iterations_used = 1;

  if (n_iter > 0 && !candidate.certificate.passed) {
    std::vector<double> Z(static_cast<size_t>(n) * static_cast<size_t>(q_cols), 0.0);
    for (int iter = 0; iter < n_iter; ++iter) {
      t0 = native_timer_now();
      F77_CALL(dgemm)(&trans, &notrans, &n, &q_cols, &m,
                      &one, const_cast<double*>(A), &m, Q.data(), &m,
                      &zero, Z.data(), &n FCONE FCONE);
      stage[stage_apply] += native_timer_elapsed(t0);

      t0 = native_timer_now();
      dense_randomized_thin_qr(Z, n, q_cols);
      stage[stage_normalize] += native_timer_elapsed(t0);

      t0 = native_timer_now();
      F77_CALL(dgemm)(&notrans, &notrans, &m, &q_cols, &n,
                      &one, const_cast<double*>(A), &m, Z.data(), &n,
                      &zero, Q.data(), &m FCONE FCONE);
      stage[stage_apply] += native_timer_elapsed(t0);
      matvecs += 2;

      t0 = native_timer_now();
      dense_randomized_thin_qr(Q, m, q_cols);
      stage[stage_normalize] += native_timer_elapsed(t0);
    }
    candidate = dense_randomized_candidate(
      A, m, n, rank, Q, q_cols, tol,
      &stage[stage_small_svd], &stage[stage_vector_form], &stage[stage_certificate]
    );
    matvecs += 1;
    iterations_used = n_iter + 1;
  } else if (n_iter > 0) {
    early_stop_used = true;
  }
  stage[stage_controller] = native_timer_elapsed(controller_t0);

  SEXP d_ = PROTECT(allocVector(REALSXP, rank));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, rank));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, rank));
  std::memcpy(REAL(d_), candidate.d.data(),
              sizeof(double) * static_cast<size_t>(rank));
  std::memcpy(REAL(u_), candidate.U.data(),
              sizeof(double) * static_cast<size_t>(m) * static_cast<size_t>(rank));
  std::memcpy(REAL(v_), candidate.V.data(),
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(rank));

  SEXP stage_ = PROTECT(allocVector(REALSXP, static_cast<R_xlen_t>(stage.size())));
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, static_cast<R_xlen_t>(stage.size())));
  const char* stage_names[] = {
    "random", "apply", "normalize", "small_svd", "vector_form",
    "certificate", "native_controller"
  };
  for (R_xlen_t idx = 0; idx < static_cast<R_xlen_t>(stage.size()); ++idx) {
    REAL(stage_)[idx] = stage[static_cast<size_t>(idx)];
    SET_STRING_ELT(stage_names_, idx, mkChar(stage_names[idx]));
  }
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP cert_ = PROTECT(dense_randomized_certificate_pack(candidate.certificate));
  SEXP initial_cert_ = PROTECT(dense_randomized_certificate_pack(initial_cert));
  SEXP out_ = PROTECT(allocVector(VECSXP, 13));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, cert_);
  SET_VECTOR_ELT(out_, 4, initial_cert_);
  SET_VECTOR_ELT(out_, 5, stage_);
  SET_VECTOR_ELT(out_, 6, ScalarInteger(iterations_used));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 8, ScalarLogical(early_stop_used));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(q_cols));
  SET_VECTOR_ELT(out_, 10, mkString("native_dense_randomized_controller"));
  SET_VECTOR_ELT(out_, 11, mkString("native_dense_projected_svd"));
  SET_VECTOR_ELT(out_, 12, mkString("native_direct_qt_a"));
  SEXP names_ = PROTECT(allocVector(STRSXP, 13));
  const char* names[] = {
    "d", "u", "v", "certificate_diagnostics", "initial_certificate_diagnostics",
    "stage_seconds", "iterations", "matvecs", "adaptive_stop_used",
    "sample_dimension", "controller_kind", "core_solver", "projection_kind"
  };
  for (int idx = 0; idx < 13; ++idx) {
    SET_STRING_ELT(names_, idx, mkChar(names[idx]));
  }
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(9);
  return out_;
}

extern "C" SEXP eigencore_csc_randomized_svd_controller(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP rank_, SEXP oversample_,
    SEXP n_iter_, SEXP normalizer_, SEXP tol_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC randomized controller inputs");
  }
  if (!isString(normalizer_) || LENGTH(normalizer_) < 1 ||
      std::strcmp(CHAR(STRING_ELT(normalizer_, 0)), "qr") != 0) {
    error("native CSC randomized controller currently supports only QR normalization");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int limit = (m < n) ? m : n;
  int rank = asInteger(rank_);
  int oversample = asInteger(oversample_);
  int n_iter = asInteger(n_iter_);
  const double tol = asReal(tol_);
  if (rank == NA_INTEGER || rank < 1) {
    error("rank must be a positive integer");
  }
  if (oversample == NA_INTEGER || oversample < 0) {
    error("oversample must be a non-negative integer");
  }
  if (n_iter == NA_INTEGER || n_iter < 0) {
    error("n_iter must be a non-negative integer");
  }
  if (rank > limit) {
    rank = limit;
  }
  const int q_cols = (rank + oversample < limit) ? rank + oversample : limit;
  const int nnz = LENGTH(x_);
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  std::vector<double> stage(7, 0.0);
  const int stage_random = 0;
  const int stage_apply = 1;
  const int stage_normalize = 2;
  const int stage_small_svd = 3;
  const int stage_vector_form = 4;
  const int stage_certificate = 5;
  const int stage_controller = 6;
  auto controller_t0 = native_timer_now();

  auto t0 = native_timer_now();
  std::vector<double> omega(static_cast<size_t>(n) * static_cast<size_t>(q_cols), 0.0);
  GetRNGstate();
  for (int64_t pos = 0;
       pos < static_cast<int64_t>(n) * static_cast<int64_t>(q_cols);
       ++pos) {
    omega[static_cast<size_t>(pos)] = norm_rand();
  }
  PutRNGstate();
  stage[stage_random] += native_timer_elapsed(t0);

  t0 = native_timer_now();
  std::vector<double> Q(static_cast<size_t>(m) * static_cast<size_t>(q_cols), 0.0);
  csc_randomized_apply_block(
    impl, EIGENCORE_TRANSPOSE_NONE, q_cols, omega.data(), n, Q.data(), m,
    "CSC randomized controller"
  );
  stage[stage_apply] += native_timer_elapsed(t0);

  t0 = native_timer_now();
  dense_randomized_thin_qr(Q, m, q_cols);
  stage[stage_normalize] += native_timer_elapsed(t0);

  int matvecs = 1;
  DenseRandomizedCandidate candidate = csc_randomized_candidate(
    impl, nnz, rank, Q, q_cols, tol,
    &stage[stage_small_svd], &stage[stage_vector_form], &stage[stage_certificate]
  );
  matvecs += 1;
  DenseRandomizedCertificate initial_cert = candidate.certificate;
  bool early_stop_used = false;
  int iterations_used = 1;

  if (n_iter > 0 && !candidate.certificate.passed) {
    std::vector<double> Z(static_cast<size_t>(n) * static_cast<size_t>(q_cols), 0.0);
    for (int iter = 0; iter < n_iter; ++iter) {
      t0 = native_timer_now();
      csc_randomized_apply_block(
        impl, EIGENCORE_TRANSPOSE_ADJOINT, q_cols, Q.data(), m, Z.data(), n,
        "CSC randomized controller"
      );
      stage[stage_apply] += native_timer_elapsed(t0);

      t0 = native_timer_now();
      dense_randomized_thin_qr(Z, n, q_cols);
      stage[stage_normalize] += native_timer_elapsed(t0);

      t0 = native_timer_now();
      csc_randomized_apply_block(
        impl, EIGENCORE_TRANSPOSE_NONE, q_cols, Z.data(), n, Q.data(), m,
        "CSC randomized controller"
      );
      stage[stage_apply] += native_timer_elapsed(t0);
      matvecs += 2;

      t0 = native_timer_now();
      dense_randomized_thin_qr(Q, m, q_cols);
      stage[stage_normalize] += native_timer_elapsed(t0);
    }
    candidate = csc_randomized_candidate(
      impl, nnz, rank, Q, q_cols, tol,
      &stage[stage_small_svd], &stage[stage_vector_form], &stage[stage_certificate]
    );
    matvecs += 1;
    iterations_used = n_iter + 1;
  } else if (n_iter > 0) {
    early_stop_used = true;
  }
  stage[stage_controller] = native_timer_elapsed(controller_t0);

  SEXP d_ = PROTECT(allocVector(REALSXP, rank));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, rank));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, rank));
  std::memcpy(REAL(d_), candidate.d.data(),
              sizeof(double) * static_cast<size_t>(rank));
  std::memcpy(REAL(u_), candidate.U.data(),
              sizeof(double) * static_cast<size_t>(m) * static_cast<size_t>(rank));
  std::memcpy(REAL(v_), candidate.V.data(),
              sizeof(double) * static_cast<size_t>(n) * static_cast<size_t>(rank));

  SEXP stage_ = PROTECT(allocVector(REALSXP, static_cast<R_xlen_t>(stage.size())));
  SEXP stage_names_ = PROTECT(allocVector(STRSXP, static_cast<R_xlen_t>(stage.size())));
  const char* stage_names[] = {
    "random", "apply", "normalize", "small_svd", "vector_form",
    "certificate", "native_controller"
  };
  for (R_xlen_t idx = 0; idx < static_cast<R_xlen_t>(stage.size()); ++idx) {
    REAL(stage_)[idx] = stage[static_cast<size_t>(idx)];
    SET_STRING_ELT(stage_names_, idx, mkChar(stage_names[idx]));
  }
  setAttrib(stage_, R_NamesSymbol, stage_names_);

  SEXP cert_ = PROTECT(dense_randomized_certificate_pack(candidate.certificate));
  SEXP initial_cert_ = PROTECT(dense_randomized_certificate_pack(initial_cert));
  SEXP out_ = PROTECT(allocVector(VECSXP, 13));
  SET_VECTOR_ELT(out_, 0, d_);
  SET_VECTOR_ELT(out_, 1, u_);
  SET_VECTOR_ELT(out_, 2, v_);
  SET_VECTOR_ELT(out_, 3, cert_);
  SET_VECTOR_ELT(out_, 4, initial_cert_);
  SET_VECTOR_ELT(out_, 5, stage_);
  SET_VECTOR_ELT(out_, 6, ScalarInteger(iterations_used));
  SET_VECTOR_ELT(out_, 7, ScalarInteger(matvecs));
  SET_VECTOR_ELT(out_, 8, ScalarLogical(early_stop_used));
  SET_VECTOR_ELT(out_, 9, ScalarInteger(q_cols));
  SET_VECTOR_ELT(out_, 10, mkString("native_csc_randomized_controller"));
  SET_VECTOR_ELT(out_, 11, mkString("native_dense_projected_svd"));
  SET_VECTOR_ELT(out_, 12, mkString("native_direct_qt_a"));
  SEXP names_ = PROTECT(allocVector(STRSXP, 13));
  const char* names[] = {
    "d", "u", "v", "certificate_diagnostics", "initial_certificate_diagnostics",
    "stage_seconds", "iterations", "matvecs", "adaptive_stop_used",
    "sample_dimension", "controller_kind", "core_solver", "projection_kind"
  };
  for (int idx = 0; idx < 13; ++idx) {
    SET_STRING_ELT(names_, idx, mkChar(names[idx]));
  }
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(9);
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

extern "C" SEXP eigencore_csc_randomized_apply(SEXP i_, SEXP p_, SEXP x_,
                                               SEXP dim_, SEXP X_,
                                               SEXP transpose_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(X_) || !isLogical(transpose_)) {
    error("invalid CSC randomized apply inputs");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  if (dimX == R_NilValue) {
    error("X must be a matrix");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;
  if (xr != inner) {
    error("non-conformable X for CSC randomized apply");
  }

  SEXP out_ = PROTECT(allocMatrix(REALSXP, out_rows, xc));
  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  CSCOperator impl = {m, n, row_idx, col_ptr, values};
  const int status = eigencore_csc_apply(
    &impl,
    transpose ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    REAL(X_),
    xr,
    1.0,
    0.0,
    REAL(out_),
    out_rows,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error("CSC randomized apply", status);
  }

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_csc_randomized_sketch(SEXP i_, SEXP p_, SEXP x_,
                                                SEXP dim_, SEXP cols_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC randomized sketch inputs");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int sketch_cols = asInteger(cols_);
  if (sketch_cols == NA_INTEGER || sketch_cols < 0) {
    error("sketch column count must be non-negative");
  }

  SEXP out_ = PROTECT(allocMatrix(REALSXP, m, sketch_cols));
  double* out = REAL(out_);
  std::memset(out, 0, sizeof(double) * static_cast<size_t>(m) *
                       static_cast<size_t>(sketch_cols));
  if (sketch_cols == 0) {
    UNPROTECT(1);
    return out_;
  }

  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  GetRNGstate();
  for (int block = 0; block < sketch_cols; ++block) {
    double* out_col = out + static_cast<int64_t>(block) * m;
    for (int col = 0; col < n; ++col) {
      const double omega = norm_rand();
      for (int pos = col_ptr[col]; pos < col_ptr[col + 1]; ++pos) {
        out_col[row_idx[pos]] += values[pos] * omega;
      }
    }
  }
  PutRNGstate();

  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_csc_randomized_project_transposed(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP Q_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(Q_)) {
    error("invalid CSC randomized projection inputs");
  }
  SEXP dimQ = getAttrib(Q_, R_DimSymbol);
  if (dimQ == R_NilValue) {
    error("Q must be a matrix");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int qr = INTEGER(dimQ)[0];
  const int qcols = INTEGER(dimQ)[1];
  if (qr != m) {
    error("non-conformable Q for CSC randomized projection");
  }

  SEXP out_ = PROTECT(allocMatrix(REALSXP, qcols, n));
  double* out = REAL(out_);
  std::memset(out, 0, sizeof(double) * static_cast<size_t>(qcols) *
                       static_cast<size_t>(n));
  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  const double* Q = REAL(Q_);
  // Row-major copy of Q so each nonzero reads a contiguous qcols-length panel
  // instead of striding m doubles per element across Q's columns.
  std::vector<double> Qt(static_cast<size_t>(m) * static_cast<size_t>(qcols));
  for (int block = 0; block < qcols; ++block) {
    const double* q_col = Q + static_cast<int64_t>(block) * m;
    for (int row = 0; row < m; ++row) {
      Qt[static_cast<int64_t>(row) * qcols + block] = q_col[row];
    }
  }
  for (int col = 0; col < n; ++col) {
    double* out_col = out + static_cast<int64_t>(col) * qcols;
    for (int pos = col_ptr[col]; pos < col_ptr[col + 1]; ++pos) {
      const int row = row_idx[pos];
      const double a = values[pos];
      const double* qt_row = Qt.data() + static_cast<int64_t>(row) * qcols;
      for (int block = 0; block < qcols; ++block) {
        out_col[block] += a * qt_row[block];
      }
    }
  }
  SEXP transposed_ = PROTECT(ScalarLogical(TRUE));
  setAttrib(out_, install("transposed"), transposed_);
  UNPROTECT(2);
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
    // int64_t: col * rows overflows int once the matrix exceeds 2^31 elements.
    const int64_t offset = static_cast<int64_t>(col) * rows;
    for (int row = 0; row < rows; ++row) {
      const long double value = X[offset + row];
      sum += value * value;
    }
    out[col] = sqrt(static_cast<double>(sum));
  }
  UNPROTECT(1);
  return out_;
}
