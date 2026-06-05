#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <cfloat>
#include <cmath>
#include "native_operators.h"
#include "certificates.h"

static SEXP workspace_counters_cert(EigencoreWorkspace* workspace) {
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

static double max_orthogonality_loss_cert(const double* gram, int k) {
  double loss = 0.0;
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row < k; ++row) {
      const double target = (row == col) ? 1.0 : 0.0;
      const double diff = fabs(gram[row + col * k] - target);
      if (diff > loss) {
        loss = diff;
      }
    }
  }
  return loss;
}

static double column_norm_cert(const double* X, int rows, int col) {
  long double sum = 0.0L;
  const int offset = col * rows;
  for (int row = 0; row < rows; ++row) {
    const long double value = X[offset + row];
    sum += value * value;
  }
  return sqrt(static_cast<double>(sum));
}

static double frobenius_norm_dense_cert(const double* X, int len) {
  long double sum = 0.0L;
  for (int i = 0; i < len; ++i) {
    const long double value = X[i];
    sum += value * value;
  }
  return sqrt(static_cast<double>(sum));
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

  const double loss = max_orthogonality_loss_cert(REAL(gram_), k);
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
    REAL(out_)[col] = column_norm_cert(REAL(residual_), n, col);
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
  const double norm_A = frobenius_norm_dense_cert(REAL(A_), n * n);
  const double norm_B = (B_ == R_NilValue) ? 1.0 : frobenius_norm_dense_cert(REAL(B_), n * n);

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
    const double residual = column_norm_cert(REAL(residual_matrix_), n, col);
    const double vector_norm = column_norm_cert(REAL(vectors_), n, col);
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
  const double orth = max_orthogonality_loss_cert(REAL(gram_), k);

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
    const double left = column_norm_cert(REAL(left_matrix_), m, col);
    const double right = column_norm_cert(REAL(right_matrix_), n, col);
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
  const double norm_A = frobenius_norm_dense_cert(REAL(A_), m * n);
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
    const double left = column_norm_cert(REAL(left_matrix_), m, col);
    const double right = column_norm_cert(REAL(right_matrix_), n, col);
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
  REAL(orth_)[0] = max_orthogonality_loss_cert(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss_cert(REAL(gram_v_), k);
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
  const double norm_A = frobenius_norm_dense_cert(REAL(A_), m * n);
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
    const double residual = column_norm_cert(REAL(residual_matrix_), n, col);
    const double vector_norm = column_norm_cert(REAL(vectors_), n, col);
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
  const double orth = max_orthogonality_loss_cert(REAL(gram_), k);

  SEXP out_ = PROTECT(allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, workspace_counters_cert(&workspace));
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
    const double left = column_norm_cert(REAL(left_matrix_), m, col);
    const double right = column_norm_cert(REAL(right_matrix_), n, col);
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
  REAL(orth_)[0] = max_orthogonality_loss_cert(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss_cert(REAL(gram_v_), k);
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
  SET_VECTOR_ELT(out_, 8, workspace_counters_cert(&workspace));
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

SEXP native_operator_svd_certificate_cached_av(void* impl,
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
    const double right = column_norm_cert(REAL(right_matrix_), n, col);
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
  REAL(orth_)[0] = max_orthogonality_loss_cert(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss_cert(REAL(gram_v_), k);
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
  SET_VECTOR_ELT(out_, 8, workspace_counters_cert(&workspace));
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
