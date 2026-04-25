#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <vector>
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
    if (block_cols == 1) {
      for (int64_t col = 0; col < csc->cols; ++col) {
        const double x_col = X[col];
        if (x_col == 0.0) continue;
        for (int pos = csc->col_ptr[col]; pos < csc->col_ptr[col + 1]; ++pos) {
          Y[csc->row_idx[pos]] += alpha * csc->values[pos] * x_col;
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
  for (int i = 0; i < count; ++i) {
    int best = -1;
    for (int j = 0; j < n; ++j) {
      bool already = false;
      for (int prev = 0; prev < i; ++prev) {
        if (selected[prev] == j) {
          already = true;
          break;
        }
      }
      if (already) {
        continue;
      }
      if (best < 0 || ritz_value_better(values[j], values[best], target_kind)) {
        best = j;
      }
    }
    selected[i] = best;
  }
  return count;
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
                                  const double* start,
                                  double* U,
                                  double* V,
                                  double* alpha,
                                  double* beta,
                                  int* iterations,
                                  int* matvecs) {
  double* v = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* z = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  double* u = static_cast<double*>(std::calloc(static_cast<size_t>(m), sizeof(double)));
  double* u_prev = static_cast<double*>(std::calloc(static_cast<size_t>(m), sizeof(double)));
  if (v == nullptr || z == nullptr || u == nullptr || u_prev == nullptr) {
    std::free(v);
    std::free(z);
    std::free(u);
    std::free(u_prev);
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
  double beta_prev = 0.0;
  EigencoreWorkspace workspace = {0, 0, nullptr, 0};

  for (int j = 0; j < maxit; ++j) {
    *iterations = j + 1;
    std::memcpy(V + j * n, v, sizeof(double) * static_cast<size_t>(n));
    std::memset(u, 0, sizeof(double) * static_cast<size_t>(m));

    int status = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1, v, n,
                       1.0, 0.0, u, m, &workspace);
    if (status != 0) {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      return status;
    }
    ++(*matvecs);
    if (j > 0) {
      for (int row = 0; row < m; ++row) {
        u[row] -= beta_prev * u_prev[row];
      }
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int prev = 0; prev < j; ++prev) {
        const double* uprev_basis = U + prev * m;
        long double dot = 0.0L;
        for (int row = 0; row < m; ++row) {
          dot += static_cast<long double>(uprev_basis[row]) * u[row];
        }
        const double coeff = static_cast<double>(dot);
        for (int row = 0; row < m; ++row) {
          u[row] -= coeff * uprev_basis[row];
        }
      }
    }

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

    std::memset(z, 0, sizeof(double) * static_cast<size_t>(n));
    status = apply(impl, EIGENCORE_TRANSPOSE_ADJOINT, 1, u, m,
                   1.0, 0.0, z, n, &workspace);
    if (status != 0) {
      std::free(v);
      std::free(z);
      std::free(u);
      std::free(u_prev);
      return status;
    }
    ++(*matvecs);
    for (int row = 0; row < n; ++row) {
      z[row] -= alpha[j] * v[row];
    }

    for (int pass = 0; pass < 2; ++pass) {
      for (int prev = 0; prev <= j; ++prev) {
        const double* vprev_basis = V + prev * n;
        long double dot = 0.0L;
        for (int row = 0; row < n; ++row) {
          dot += static_cast<long double>(vprev_basis[row]) * z[row];
        }
        const double coeff = static_cast<double>(dot);
        for (int row = 0; row < n; ++row) {
          z[row] -= coeff * vprev_basis[row];
        }
      }
    }

    long double beta_norm2 = 0.0L;
    for (int row = 0; row < n; ++row) {
      beta_norm2 += static_cast<long double>(z[row]) * z[row];
    }
    beta[j] = sqrt(static_cast<double>(beta_norm2));
    if (j + 1 == maxit || beta[j] <= 100.0 * DBL_EPSILON) {
      break;
    }

    std::memcpy(u_prev, u, sizeof(double) * static_cast<size_t>(m));
    beta_prev = beta[j];
    for (int row = 0; row < n; ++row) {
      v[row] = z[row] / beta[j];
    }
  }

  std::free(v);
  std::free(z);
  std::free(u);
  std::free(u_prev);
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

extern "C" SEXP eigencore_golub_kahan_dense(SEXP A_, SEXP maxit_, SEXP start_) {
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
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, maxit));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(U_), 0, sizeof(double) * static_cast<size_t>(m) * maxit);
  std::memset(REAL(V_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));

  DenseColumnMajorOperator impl = {m, n, REAL(A_)};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_dense_apply, m, n,
                                            maxit, REAL(start_), REAL(U_),
                                            REAL(V_), REAL(alpha_), REAL(beta_),
                                            &iterations, &matvecs);
  if (status != 0) {
    error("native dense Golub-Kahan failed with status=%d", status);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_csc(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                          SEXP maxit_, SEXP start_) {
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
  if (maxit < 1 || maxit > limit) {
    error("maxit must be between 1 and min(dim(A))");
  }

  SEXP U_ = PROTECT(allocMatrix(REALSXP, m, maxit));
  SEXP V_ = PROTECT(allocMatrix(REALSXP, n, maxit));
  SEXP alpha_ = PROTECT(allocVector(REALSXP, maxit));
  SEXP beta_ = PROTECT(allocVector(REALSXP, maxit));
  std::memset(REAL(U_), 0, sizeof(double) * static_cast<size_t>(m) * maxit);
  std::memset(REAL(V_), 0, sizeof(double) * static_cast<size_t>(n) * maxit);
  std::memset(REAL(alpha_), 0, sizeof(double) * static_cast<size_t>(maxit));
  std::memset(REAL(beta_), 0, sizeof(double) * static_cast<size_t>(maxit));

  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  int iterations = 0;
  int matvecs = 0;
  const int status = native_golub_kahan_run(&impl, eigencore_csc_apply, m, n,
                                            maxit, REAL(start_), REAL(U_),
                                            REAL(V_), REAL(alpha_), REAL(beta_),
                                            &iterations, &matvecs);
  if (status != 0) {
    error("native CSC Golub-Kahan failed with status=%d", status);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out_, 0, U_);
  SET_VECTOR_ELT(out_, 1, V_);
  SET_VECTOR_ELT(out_, 2, alpha_);
  SET_VECTOR_ELT(out_, 3, beta_);
  SET_VECTOR_ELT(out_, 4, ScalarInteger(iterations));
  SET_VECTOR_ELT(out_, 5, ScalarInteger(matvecs));
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names_, 0, mkChar("U"));
  SET_STRING_ELT(names_, 1, mkChar("V"));
  SET_STRING_ELT(names_, 2, mkChar("alpha"));
  SET_STRING_ELT(names_, 3, mkChar("beta"));
  SET_STRING_ELT(names_, 4, mkChar("iterations"));
  SET_STRING_ELT(names_, 5, mkChar("matvecs"));
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
                             double* z, double* tmp, int n) {
  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  const double minus_one = -1.0;
  int incx = 1;
  for (int pass = 0; pass < 2; ++pass) {
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

struct ThickRestartBuffers {
  double* V_active;
  double* AV_active;
  double* S_eig;       // m_max x m_max (also serves as projection scratch)
  double* S_selected;  // selected Ritz vectors, m_max x m_max
  double* theta;       // m_max
  double* B_v;         // n x m_max
  double* B_av;        // n x m_max
  double* z;           // n
  double* tmp;         // max(k_target, m_max)
  double* ritz_res;    // m_max
  int*    selected;    // m_max
  int*    is_locked;   // m_max
  double* dsyev_work;
  int     dsyev_lwork;
};

static void trl_buffers_free(ThickRestartBuffers* b) {
  std::free(b->V_active);
  std::free(b->AV_active);
  std::free(b->S_eig);
  std::free(b->S_selected);
  std::free(b->theta);
  std::free(b->B_v);
  std::free(b->B_av);
  std::free(b->z);
  std::free(b->tmp);
  std::free(b->ritz_res);
  std::free(b->selected);
  std::free(b->is_locked);
  std::free(b->dsyev_work);
}

static int trl_buffers_alloc(ThickRestartBuffers* b, int n, int k_target, int m_max) {
  std::memset(b, 0, sizeof(*b));
  const size_t nm = static_cast<size_t>(n) * static_cast<size_t>(m_max);
  const size_t mm = static_cast<size_t>(m_max) * static_cast<size_t>(m_max);
  b->V_active  = static_cast<double*>(std::calloc(nm, sizeof(double)));
  b->AV_active = static_cast<double*>(std::calloc(nm, sizeof(double)));
  b->S_eig     = static_cast<double*>(std::calloc(mm, sizeof(double)));
  b->S_selected = static_cast<double*>(std::calloc(mm, sizeof(double)));
  b->theta     = static_cast<double*>(std::calloc(static_cast<size_t>(m_max), sizeof(double)));
  b->B_v       = static_cast<double*>(std::calloc(nm, sizeof(double)));
  b->B_av      = static_cast<double*>(std::calloc(nm, sizeof(double)));
  b->z         = static_cast<double*>(std::calloc(static_cast<size_t>(n), sizeof(double)));
  const int tmp_len = (k_target > m_max) ? k_target : m_max;
  b->tmp       = static_cast<double*>(std::calloc(static_cast<size_t>(tmp_len > 0 ? tmp_len : 1), sizeof(double)));
  b->ritz_res  = static_cast<double*>(std::calloc(static_cast<size_t>(m_max), sizeof(double)));
  b->selected  = static_cast<int*>(std::calloc(static_cast<size_t>(m_max), sizeof(int)));
  b->is_locked = static_cast<int*>(std::calloc(static_cast<size_t>(m_max), sizeof(int)));
  b->dsyev_lwork = trl_dsyev_query(m_max);
  if (b->dsyev_lwork < 1) b->dsyev_lwork = 1;
  b->dsyev_work = static_cast<double*>(std::calloc(static_cast<size_t>(b->dsyev_lwork), sizeof(double)));
  if (b->V_active == nullptr || b->AV_active == nullptr || b->S_eig == nullptr ||
      b->S_selected == nullptr ||
      b->theta == nullptr || b->B_v == nullptr || b->B_av == nullptr ||
      b->z == nullptr || b->tmp == nullptr || b->ritz_res == nullptr ||
      b->selected == nullptr || b->is_locked == nullptr || b->dsyev_work == nullptr) {
    trl_buffers_free(b);
    return -1;
  }
  return 0;
}

static int native_thick_restart_lanczos_run(
    void* impl,
    EigencoreApplyFn apply,
    int n,
    int k_target,
    int m_max,
    int target_kind,
    double tol,
    int max_restarts,
    const double* start_vec,
    // Outputs
    double* V_locked_out,        // n x k_target  (column major)
    double* lambda_out,          // k_target
    double* residuals_out,       // k_target
    int* converged_out,          // k_target
    int* n_locked_out,
    int* iterations_out,
    int* matvecs_out,
    int* restarts_out,
    int* m_active_final_out
) {
  *n_locked_out = 0;
  *iterations_out = 0;
  *matvecs_out = 0;
  *restarts_out = 0;
  *m_active_final_out = 0;
  for (int i = 0; i < k_target; ++i) {
    lambda_out[i] = 0.0;
    residuals_out[i] = R_PosInf;
    converged_out[i] = 0;
    std::memset(V_locked_out + static_cast<int64_t>(i) * n, 0, sizeof(double) * static_cast<size_t>(n));
  }

  ThickRestartBuffers buf;
  if (trl_buffers_alloc(&buf, n, k_target, m_max) != 0) {
    return -2;
  }

  EigencoreWorkspace workspace = {0, 0, nullptr, 0};

  // Initialise V_active[:, 0] = start / ||start|| (or e_0 if start is zero)
  const double s_norm = trl_norm2(start_vec, n);
  if (s_norm == 0.0) {
    buf.V_active[0] = 1.0;
  } else {
    for (int i = 0; i < n; ++i) {
      buf.V_active[i] = start_vec[i] / s_norm;
    }
  }
  int m_active = 1;
  int m_av = 0;             // # of columns with AV_active up to date
  int n_locked = 0;
  int selected_count_final = 0;

  const char trans_T = 'T';
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;

  int restart_idx = 0;
  for (; restart_idx <= max_restarts; ++restart_idx) {
    // ---- Bring AV_active up to date for any column that still needs it ----
    bool lucky_breakdown = false;
    for (int jj = m_av; jj < m_active; ++jj) {
      const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1,
                           buf.V_active + static_cast<int64_t>(jj) * n, n,
                           1.0, 0.0,
                           buf.AV_active + static_cast<int64_t>(jj) * n, n,
                           &workspace);
      if (rc != 0) {
        trl_buffers_free(&buf);
        return rc;
      }
      ++(*matvecs_out);
    }
    m_av = m_active;

    // ---- Inner Lanczos sweep: extend basis to m_max ----
    while (m_active < m_max) {
      const int j = m_active - 1;
      // AV_active[:, j] is guaranteed up to date.
      std::memcpy(buf.z, buf.AV_active + static_cast<int64_t>(j) * n,
                  sizeof(double) * static_cast<size_t>(n));

      trl_orthogonalise(V_locked_out, n_locked,
                        buf.V_active, m_active,
                        buf.z, buf.tmp, n);

      const double beta = trl_norm2(buf.z, n);
      ++(*iterations_out);
      if (beta <= 100.0 * DBL_EPSILON) {
        lucky_breakdown = true;
        break;
      }
      const double inv_beta = 1.0 / beta;
      for (int i = 0; i < n; ++i) {
        buf.V_active[static_cast<int64_t>(m_active) * n + i] = buf.z[i] * inv_beta;
      }
      ++m_active;

      // Compute AV_active[:, m_active - 1] now so the next loop step (or the
      // post-loop Rayleigh-Ritz) finds the column ready.
      const int rc = apply(impl, EIGENCORE_TRANSPOSE_NONE, 1,
                           buf.V_active + static_cast<int64_t>(m_active - 1) * n, n,
                           1.0, 0.0,
                           buf.AV_active + static_cast<int64_t>(m_active - 1) * n, n,
                           &workspace);
      if (rc != 0) {
        trl_buffers_free(&buf);
        return rc;
      }
      ++(*matvecs_out);
      m_av = m_active;
    }

    // ---- Rayleigh-Ritz on V_active[:, 0..m_active] ----
    F77_CALL(dgemm)(&trans_T, &trans_N, &m_active, &m_active, &n,
                    &one, buf.V_active, &n, buf.AV_active, &n,
                    &zero, buf.S_eig, &m_active FCONE FCONE);
    // Symmetrise
    for (int i = 0; i < m_active; ++i) {
      for (int jj = i + 1; jj < m_active; ++jj) {
        const double avg = 0.5 * (buf.S_eig[i + jj * m_active] +
                                  buf.S_eig[jj + i * m_active]);
        buf.S_eig[i + jj * m_active] = avg;
        buf.S_eig[jj + i * m_active] = avg;
      }
    }
    char jobz = 'V';
    char uplo = 'U';
    int info = 0;
    int lwork = buf.dsyev_lwork;
    F77_CALL(dsyev)(&jobz, &uplo, &m_active, buf.S_eig, &m_active, buf.theta,
                    buf.dsyev_work, &lwork, &info FCONE FCONE);
    if (info != 0) {
      trl_buffers_free(&buf);
      return -3;
    }

    // ---- Sort Ritz pairs by target, then materialise only the leading
    // target-ordered subset needed for locking and restart.
    selected_ritz_indices(buf.theta, m_active, m_active, target_kind, buf.selected);

    const int k_remaining_before_lock = k_target - n_locked;
    int k_keep_budget = 2 * k_remaining_before_lock;
    if (k_keep_budget < k_remaining_before_lock + 5) {
      k_keep_budget = k_remaining_before_lock + 5;
    }
    if (k_keep_budget > m_max - 2) k_keep_budget = m_max - 2;
    if (k_keep_budget < 0) k_keep_budget = 0;
    int selected_count = k_keep_budget;
    if (selected_count < k_remaining_before_lock) {
      selected_count = k_remaining_before_lock;
    }
    if (selected_count < 1) selected_count = 1;
    if (selected_count > m_active) selected_count = m_active;
    selected_count_final = selected_count;

    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      for (int row = 0; row < m_active; ++row) {
        buf.S_selected[row + static_cast<int64_t>(p) * m_active] =
          buf.S_eig[row + static_cast<int64_t>(idx) * m_active];
      }
    }
    for (int i = 0; i < m_active; ++i) {
      buf.ritz_res[i] = R_PosInf;
    }
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &selected_count, &m_active,
                    &one, buf.V_active, &n, buf.S_selected, &m_active,
                    &zero, buf.B_v, &n FCONE FCONE);
    F77_CALL(dgemm)(&trans_N, &trans_N, &n, &selected_count, &m_active,
                    &one, buf.AV_active, &n, buf.S_selected, &m_active,
                    &zero, buf.B_av, &n FCONE FCONE);
    for (int p = 0; p < selected_count; ++p) {
      const int idx = buf.selected[p];
      long double s = 0.0L;
      for (int row = 0; row < n; ++row) {
        const double diff = buf.B_av[row + static_cast<int64_t>(p) * n] -
                            buf.theta[idx] * buf.B_v[row + static_cast<int64_t>(p) * n];
        s += static_cast<long double>(diff) * diff;
      }
      buf.ritz_res[idx] = sqrt(static_cast<double>(s));
    }

    // ---- Lock newly converged ----
    for (int i = 0; i < m_active; ++i) buf.is_locked[i] = 0;
    const int wanted = k_target - n_locked;
    for (int i = 0; i < wanted && i < selected_count; ++i) {
      const int idx = buf.selected[i];
      const double scale_i = (fabs(buf.theta[idx]) > 1.0) ? fabs(buf.theta[idx]) : 1.0;
      if (buf.ritz_res[idx] <= tol * scale_i) {
        std::memcpy(V_locked_out + static_cast<int64_t>(n_locked) * n,
                    buf.B_v + static_cast<int64_t>(i) * n,
                    sizeof(double) * static_cast<size_t>(n));
        lambda_out[n_locked] = buf.theta[idx];
        residuals_out[n_locked] = buf.ritz_res[idx];
        converged_out[n_locked] = 1;
        ++n_locked;
        buf.is_locked[idx] = 1;
      } else {
        break;  // first non-converged in target order; stop locking
      }
    }

    if (n_locked >= k_target) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }
    if (lucky_breakdown || restart_idx == max_restarts) {
      *restarts_out = restart_idx;
      *m_active_final_out = m_active;
      break;
    }

    // ---- Thick restart: keep top k_keep unlocked Ritz pairs ----
    const int k_remaining = k_target - n_locked;
    int k_keep = 2 * k_remaining;
    if (k_keep < k_remaining + 5) k_keep = k_remaining + 5;
    if (k_keep > m_max - 2) k_keep = m_max - 2;
    if (k_keep < 0) k_keep = 0;
    int unlocked_count = 0;
    for (int i = 0; i < m_active; ++i) {
      if (!buf.is_locked[i]) ++unlocked_count;
    }
    if (k_keep > unlocked_count) k_keep = unlocked_count;

    // Pick top k_keep unlocked Ritz indices in target order
    const int old_m_active = m_active;
    int n_picked = 0;
    for (int i = 0; i < selected_count && n_picked < k_keep; ++i) {
      const int idx = buf.selected[i];
      if (!buf.is_locked[idx]) {
        // Materialise into V_active / AV_active using the n_picked slot
        std::memcpy(buf.V_active + static_cast<int64_t>(n_picked) * n,
                    buf.B_v + static_cast<int64_t>(i) * n,
                    sizeof(double) * static_cast<size_t>(n));
        std::memcpy(buf.AV_active + static_cast<int64_t>(n_picked) * n,
                    buf.B_av + static_cast<int64_t>(i) * n,
                    sizeof(double) * static_cast<size_t>(n));
        ++n_picked;
      }
    }
    m_active = n_picked;
    m_av = n_picked;  // AV is up to date for all kept Ritz pairs (B_av columns)

    // ---- Generate the continuation tail from an unconverged Ritz residual ----
    bool have_residual_tail = false;
    for (int i = 0; i < selected_count && !have_residual_tail; ++i) {
      const int idx = buf.selected[i];
      if (idx < 0 || idx >= old_m_active || buf.is_locked[idx]) {
        continue;
      }
      if (buf.ritz_res[idx] <= 100.0 * DBL_EPSILON) {
        continue;
      }
      for (int row = 0; row < n; ++row) {
        buf.z[row] = buf.B_av[row + static_cast<int64_t>(i) * n] -
                     buf.theta[idx] * buf.B_v[row + static_cast<int64_t>(i) * n];
      }
      have_residual_tail = true;
    }
    if (!have_residual_tail) {
      uint32_t state = static_cast<uint32_t>((restart_idx + 1) * 2654435761u) ^ 0x9E3779B9u;
      if (state == 0u) state = 1u;
      for (int i = 0; i < n; ++i) {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        const double u = static_cast<double>(state) / 4294967295.0 - 0.5;
        buf.z[i] = u;
      }
    }
    trl_orthogonalise(V_locked_out, n_locked,
                      buf.V_active, m_active,
                      buf.z, buf.tmp, n);
    double tail_norm = trl_norm2(buf.z, n);
    if (tail_norm <= 100.0 * DBL_EPSILON) {
      // Try a sequence of canonical basis vectors as fallback tails
      bool recovered = false;
      for (int attempt = 0; attempt < 8 && !recovered; ++attempt) {
        std::memset(buf.z, 0, sizeof(double) * static_cast<size_t>(n));
        const int idx_basis = ((restart_idx + 1) * 7 + 13 * (attempt + 1)) % n;
        buf.z[idx_basis < 0 ? -idx_basis : idx_basis] = 1.0;
        trl_orthogonalise(V_locked_out, n_locked,
                          buf.V_active, m_active,
                          buf.z, buf.tmp, n);
        tail_norm = trl_norm2(buf.z, n);
        if (tail_norm > 100.0 * DBL_EPSILON) {
          recovered = true;
        }
      }
      if (!recovered) {
        // Subspace exhausted; bail out with current locks.
        *restarts_out = restart_idx;
        *m_active_final_out = m_active;
        break;
      }
    }
    const double inv_tail = 1.0 / tail_norm;
    for (int i = 0; i < n; ++i) {
      buf.V_active[static_cast<int64_t>(m_active) * n + i] = buf.z[i] * inv_tail;
    }
    ++m_active;
  }

  // ---- Pad output with top unlocked Ritz pairs from the last sweep ----
  *n_locked_out = n_locked;  // genuine locks only
  int n_returned = n_locked;
  if (n_returned < k_target) {
    for (int i = 0; i < selected_count_final && n_returned < k_target; ++i) {
      const int idx = buf.selected[i];
      if (idx < 0 || idx >= m_active) continue;
      if (buf.is_locked[idx]) continue;
      std::memcpy(V_locked_out + static_cast<int64_t>(n_returned) * n,
                  buf.B_v + static_cast<int64_t>(i) * n,
                  sizeof(double) * static_cast<size_t>(n));
      lambda_out[n_returned] = buf.theta[idx];
      residuals_out[n_returned] = buf.ritz_res[idx];
      converged_out[n_returned] = 0;
      ++n_returned;
    }
  }
  trl_buffers_free(&buf);
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

static SEXP block_lanczos_pack_result(int n, int k_target, const double* V,
                                      const double* lambda,
                                      const double* residuals,
                                      const int* converged, int nconv,
                                      int iterations, int matvecs,
                                      int m_active_final) {
  return trl_pack_result(n, k_target, V, lambda, residuals, converged, nconv,
                         iterations, matvecs, 0, m_active_final);
}

extern "C" SEXP eigencore_thick_restart_lanczos_dense(SEXP A_, SEXP k_,
                                                      SEXP m_max_,
                                                      SEXP target_kind_,
                                                      SEXP tol_,
                                                      SEXP max_restarts_,
                                                      SEXP start_) {
  if (!isReal(A_) || !isReal(start_)) {
    error("A and start must be double");
  }
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  if (dimA == R_NilValue) {
    error("A must be a matrix");
  }
  const int n = INTEGER(dimA)[0];
  if (INTEGER(dimA)[1] != n || LENGTH(start_) != n) {
    error("non-conformable thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  if (k < 1) error("k must be >= 1");
  if (m_max < k + 1) error("m_max must be >= k + 1");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V_locked(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, restarts = 0, m_active = 0;

  DenseColumnMajorOperator impl = {n, n, REAL(A_)};
  const int status = native_thick_restart_lanczos_run(
    &impl, eigencore_dense_apply,
    n, k, m_max, target_kind, tol, max_restarts,
    REAL(start_),
    V_locked.data(), lambda.data(), residuals.data(), converged.data(),
    &n_locked, &iterations, &matvecs, &restarts, &m_active);
  if (status != 0) {
    error("native thick-restart dense Lanczos failed with status=%d", status);
  }
  return trl_pack_result(n, k, V_locked.data(), lambda.data(),
                         residuals.data(), converged.data(),
                         n_locked, iterations, matvecs, restarts, m_active);
}

extern "C" SEXP eigencore_thick_restart_lanczos_csc(SEXP i_, SEXP p_, SEXP x_,
                                                    SEXP dim_, SEXP k_,
                                                    SEXP m_max_,
                                                    SEXP target_kind_,
                                                    SEXP tol_,
                                                    SEXP max_restarts_,
                                                    SEXP start_) {
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(start_)) {
    error("invalid CSC thick-restart Lanczos inputs");
  }
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n || LENGTH(start_) != n) {
    error("non-conformable CSC thick-restart Lanczos inputs");
  }
  const int k = static_cast<int>(asInteger(k_));
  const int m_max = static_cast<int>(asInteger(m_max_));
  const int target_kind = static_cast<int>(asInteger(target_kind_));
  const double tol = asReal(tol_);
  const int max_restarts = static_cast<int>(asInteger(max_restarts_));
  if (k < 1) error("k must be >= 1");
  if (m_max < k + 1) error("m_max must be >= k + 1");
  if (m_max > n) error("m_max must be <= nrow(A)");
  if (max_restarts < 0) error("max_restarts must be >= 0");

  std::vector<double> V_locked(static_cast<size_t>(n) * static_cast<size_t>(k), 0.0);
  std::vector<double> lambda(static_cast<size_t>(k), 0.0);
  std::vector<double> residuals(static_cast<size_t>(k), R_PosInf);
  std::vector<int> converged(static_cast<size_t>(k), 0);
  int n_locked = 0, iterations = 0, matvecs = 0, restarts = 0, m_active = 0;

  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  const int status = native_thick_restart_lanczos_run(
    &impl, eigencore_csc_apply,
    n, k, m_max, target_kind, tol, max_restarts,
    REAL(start_),
    V_locked.data(), lambda.data(), residuals.data(), converged.data(),
    &n_locked, &iterations, &matvecs, &restarts, &m_active);
  if (status != 0) {
    error("native thick-restart CSC Lanczos failed with status=%d", status);
  }
  return trl_pack_result(n, k, V_locked.data(), lambda.data(),
                         residuals.data(), converged.data(),
                         n_locked, iterations, matvecs, restarts, m_active);
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

  SEXP residual_matrix_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP Bv_ = PROTECT(allocMatrix(REALSXP, n, k));
  SEXP gram_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));

  F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                  &one, REAL(A_), &n, REAL(vectors_), &n,
                  &zero, REAL(residual_matrix_), &n FCONE FCONE);

  if (B_ == R_NilValue) {
    std::memcpy(REAL(Bv_), REAL(vectors_), sizeof(double) * static_cast<size_t>(n) * k);
  } else {
    F77_CALL(dgemm)(&notrans, &notrans, &n, &k, &n,
                    &one, REAL(B_), &n, REAL(vectors_), &n,
                    &zero, REAL(Bv_), &n FCONE FCONE);
  }

  for (int col = 0; col < k; ++col) {
    const double lambda = REAL(values_)[col];
    const int offset = col * n;
    for (int row = 0; row < n; ++row) {
      REAL(residual_matrix_)[offset + row] -= lambda * REAL(Bv_)[offset + row];
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
                  &one, REAL(vectors_), &n, REAL(Bv_), &n,
                  &zero, REAL(gram_), &k FCONE FCONE);
  const double orth = max_orthogonality_loss(REAL(gram_), k);

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, ScalarReal(norm_A));
  SEXP names_ = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("scale"));
  SET_STRING_ELT(names_, 4, mkChar("converged"));
  SET_STRING_ELT(names_, 5, mkChar("norm_A"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(9);
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

#include "projection/golub_kahan_ritz.hpp"

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

static const R_CallMethodDef CallEntries[] = {
  {"eigencore_dense_block_apply", (DL_FUNC) &eigencore_dense_block_apply, 6},
  {"eigencore_csc_block_apply", (DL_FUNC) &eigencore_csc_block_apply, 9},
  {"eigencore_diagonal_block_apply", (DL_FUNC) &eigencore_diagonal_block_apply, 7},
  {"eigencore_native_apply_noalloc_check", (DL_FUNC) &eigencore_native_apply_noalloc_check, 4},
  {"eigencore_col_norms", (DL_FUNC) &eigencore_col_norms, 1},
  {"eigencore_mgs2", (DL_FUNC) &eigencore_mgs2, 2},
  {"eigencore_cholqr2", (DL_FUNC) &eigencore_cholqr2, 1},
  {"eigencore_b_cholqr2", (DL_FUNC) &eigencore_b_cholqr2, 2},
  {"eigencore_reorthogonalize_against", (DL_FUNC) &eigencore_reorthogonalize_against, 3},
  {"eigencore_basis_workspace_create", (DL_FUNC) &eigencore_basis_workspace_create, 3},
  {"eigencore_basis_workspace_info", (DL_FUNC) &eigencore_basis_workspace_info, 1},
  {"eigencore_reorthogonalize_against_workspace", (DL_FUNC) &eigencore_reorthogonalize_against_workspace, 4},
  {"eigencore_lanczos_dense", (DL_FUNC) &eigencore_lanczos_dense, 6},
  {"eigencore_lanczos_csc", (DL_FUNC) &eigencore_lanczos_csc, 9},
  {"eigencore_golub_kahan_dense", (DL_FUNC) &eigencore_golub_kahan_dense, 3},
  {"eigencore_golub_kahan_csc", (DL_FUNC) &eigencore_golub_kahan_csc, 6},
  {"eigencore_thick_restart_lanczos_dense", (DL_FUNC) &eigencore_thick_restart_lanczos_dense, 7},
  {"eigencore_thick_restart_lanczos_csc", (DL_FUNC) &eigencore_thick_restart_lanczos_csc, 10},
  {"eigencore_block_lanczos_dense", (DL_FUNC) &eigencore_block_lanczos_dense, 7},
  {"eigencore_block_lanczos_csc", (DL_FUNC) &eigencore_block_lanczos_csc, 10},
  {"eigencore_orthogonality_loss", (DL_FUNC) &eigencore_orthogonality_loss, 2},
  {"eigencore_dense_eigen_residuals", (DL_FUNC) &eigencore_dense_eigen_residuals, 4},
  {"eigencore_dense_eigen_certificate", (DL_FUNC) &eigencore_dense_eigen_certificate, 5},
  {"eigencore_dense_svd_residuals", (DL_FUNC) &eigencore_dense_svd_residuals, 4},
  {"eigencore_dense_svd_certificate", (DL_FUNC) &eigencore_dense_svd_certificate, 5},
  {"eigencore_csc_eigen_certificate", (DL_FUNC) &eigencore_csc_eigen_certificate, 8},
  {"eigencore_diagonal_eigen_certificate", (DL_FUNC) &eigencore_diagonal_eigen_certificate, 7},
  {"eigencore_csc_svd_certificate", (DL_FUNC) &eigencore_csc_svd_certificate, 9},
  {"eigencore_diagonal_svd_certificate", (DL_FUNC) &eigencore_diagonal_svd_certificate, 8},
  {"eigencore_rayleigh_ritz_symmetric", (DL_FUNC) &eigencore_rayleigh_ritz_symmetric, 2},
  {"eigencore_tridiagonal_eigen", (DL_FUNC) &eigencore_tridiagonal_eigen, 2},
  {"eigencore_bidiagonal_svd", (DL_FUNC) &eigencore_bidiagonal_svd, 2},
  {"eigencore_golub_kahan_ritz", (DL_FUNC) &eigencore_golub_kahan_ritz, 6},
  {"eigencore_dense_symmetric_eigen", (DL_FUNC) &eigencore_dense_symmetric_eigen, 1},
  {"eigencore_dense_generalized_spd_eigen", (DL_FUNC) &eigencore_dense_generalized_spd_eigen, 2},
  {"eigencore_dense_svd", (DL_FUNC) &eigencore_dense_svd, 1},
  {"eigencore_tridiagonal_solve", (DL_FUNC) &eigencore_tridiagonal_solve, 4},
  {NULL, NULL, 0}
};

extern "C" void R_init_eigencore(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
