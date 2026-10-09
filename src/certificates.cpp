#include <cfloat>
#include <cmath>
#include <vector>
#include "eigencore_common.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include "native_operators.h"
#include "certificates.h"

extern "C" void eigencore_validate_csc_structure(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_, const char* context);

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

// Upper-triangle variant for Gram matrices produced by dsyrk, which only
// fills the upper triangle. Equivalent to the full scan for symmetric input.
static double max_orthogonality_loss_upper_cert(const double* gram, int k) {
  double loss = 0.0;
  for (int col = 0; col < k; ++col) {
    for (int row = 0; row <= col; ++row) {
      const double target = (row == col) ? 1.0 : 0.0;
      const double diff = fabs(gram[row + col * k] - target);
      if (diff > loss) {
        loss = diff;
      }
    }
  }
  return loss;
}

// X^T X Gram product via dsyrk (upper triangle only): half the flops of the
// equivalent dgemm. Pair with max_orthogonality_loss_upper_cert.
static void gram_upper_dsyrk_cert(const double* X, int rows, int k,
                                  double* gram) {
  const char uplo = 'U';
  const char trans = 'T';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dsyrk)(&uplo, &trans, &k, &rows,
                  &one, const_cast<double*>(X), &rows,
                  &zero, gram, &k FCONE FCONE);
}

static double column_norm_cert(const double* X, int rows, int col) {
  long double sum = 0.0L;
  const R_xlen_t offset = static_cast<R_xlen_t>(col) * rows;
  for (int row = 0; row < rows; ++row) {
    const long double value = X[offset + row];
    sum += value * value;
  }
  return sqrt(static_cast<double>(sum));
}

// Largest column 2-norm of a dense column-major matrix: ||A e_j|| <= ||A||_2,
// so this is a cheap structural LOWER bound on the spectral norm (C12). The
// backward-error denominators below use lower bounds only, which makes the
// reported backward error an over-estimate and keeps `passed` sound.
static double max_column_norm_dense_cert(const double* X, int rows, int cols) {
  double best = 0.0;
  for (int col = 0; col < cols; ++col) {
    const double value = column_norm_cert(X, rows, col);
    if (value > best) {
      best = value;
    }
  }
  return best;
}

// max_j ||Y e_j|| / ||X e_j|| over columns with ||X e_j|| > 0, where Y = A X
// was computed from the operator in the same call: each ratio is <= ||A||_2.
static double applied_ratio_bound_cert(const double* applied_norms,
                                       const double* vector_norms, int k) {
  double best = 0.0;
  for (int col = 0; col < k; ++col) {
    const double denom = vector_norms[col];
    if (!(denom > 0.0) || !R_FINITE(denom) || !R_FINITE(applied_norms[col])) {
      continue;
    }
    const double ratio = applied_norms[col] / denom;
    if (ratio > best) {
      best = ratio;
    }
  }
  return best;
}

// Two-sided SVD certificate core. On entry `left_matrix` holds A V and
// `right_matrix` holds A^T U, both computed from the operator in this call.
// The spectral-norm lower bound used in the denominator is
//   L = max(norm_lower, max_j ||A v_j|| / ||v_j||, max_j ||A^T u_j|| / ||u_j||),
// each term <= ||A||_2 (C12). Returns L and the applied-vector part.
static double svd_certificate_finalize_cert(double* left_matrix,
                                            double* right_matrix,
                                            int m, int n, int k,
                                            const double* d,
                                            const double* u,
                                            const double* v,
                                            double norm_lower,
                                            double tol,
                                            double* left,
                                            double* right,
                                            double* combined,
                                            double* scale,
                                            double* backward,
                                            int* converged,
                                            double* applied_bound_out) {
  const double eps = DBL_EPSILON;
  double applied = 0.0;
  for (int col = 0; col < k; ++col) {
    const double av_norm = column_norm_cert(left_matrix, m, col);
    const double atu_norm = column_norm_cert(right_matrix, n, col);
    const double u_norm = column_norm_cert(u, m, col);
    const double v_norm = column_norm_cert(v, n, col);
    if (v_norm > 0.0 && R_FINITE(av_norm) && av_norm / v_norm > applied) {
      applied = av_norm / v_norm;
    }
    if (u_norm > 0.0 && R_FINITE(atu_norm) && atu_norm / u_norm > applied) {
      applied = atu_norm / u_norm;
    }
    const double sigma = d[col];
    const R_xlen_t left_offset = static_cast<R_xlen_t>(col) * m;
    const R_xlen_t right_offset = static_cast<R_xlen_t>(col) * n;
    for (int row = 0; row < m; ++row) {
      left_matrix[left_offset + row] -= sigma * u[left_offset + row];
    }
    for (int row = 0; row < n; ++row) {
      right_matrix[right_offset + row] -= sigma * v[right_offset + row];
    }
    left[col] = column_norm_cert(left_matrix, m, col);
    right[col] = column_norm_cert(right_matrix, n, col);
    combined[col] = sqrt(left[col] * left[col] + right[col] * right[col]);
  }
  double lower = (R_FINITE(norm_lower) && norm_lower > 0.0) ? norm_lower : 0.0;
  const double norm_A = fmax(lower, applied);
  const double scale_value = fmax(norm_A, eps);
  for (int col = 0; col < k; ++col) {
    const double be = combined[col] / scale_value;
    scale[col] = scale_value;
    backward[col] = be;
    converged[col] = (R_FINITE(be) && be <= tol) ? TRUE : FALSE;
  }
  if (applied_bound_out != nullptr) {
    *applied_bound_out = applied;
  }
  return norm_A;
}

extern "C" SEXP eigencore_orthogonality_loss(SEXP Q_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
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
    gram_upper_dsyrk_cert(REAL(Q_), n, k, REAL(gram_));
    const double loss = max_orthogonality_loss_upper_cert(REAL(gram_), k);
    UNPROTECT(1);
    return ScalarReal(loss);
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_eigen_residuals(SEXP A_, SEXP values_,
                                                SEXP vectors_, SEXP B_) {
  EIGENCORE_ENTRY_BEGIN
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
      const R_xlen_t offset = static_cast<R_xlen_t>(col) * n;
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
      const R_xlen_t offset = static_cast<R_xlen_t>(col) * n;
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_eigen_certificate(SEXP A_, SEXP values_,
                                                  SEXP vectors_, SEXP B_,
                                                  SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
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
  // Two-norm lower bounds (C12): the largest column norm of A (and B), raised
  // by ||A x_j|| / ||x_j|| for the candidate vectors themselves.
  const double col_bound_A = max_column_norm_dense_cert(REAL(A_), n, n);
  const double col_bound_B = (B_ == R_NilValue) ? 1.0 :
    max_column_norm_dense_cert(REAL(B_), n, n);

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
  SEXP vector_norms_ = PROTECT(allocVector(REALSXP, k));
  ++protect_count;
  SEXP av_norms_ = PROTECT(allocVector(REALSXP, k));
  ++protect_count;
  SEXP bv_norms_ = PROTECT(allocVector(REALSXP, k));
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
    const R_xlen_t offset = static_cast<R_xlen_t>(col) * n;
    REAL(av_norms_)[col] = column_norm_cert(REAL(residual_matrix_), n, col);
    REAL(vector_norms_)[col] = column_norm_cert(REAL(vectors_), n, col);
    REAL(bv_norms_)[col] = (B_ == R_NilValue) ? REAL(vector_norms_)[col] :
      column_norm_cert(bv, n, col);
    for (int row = 0; row < n; ++row) {
      REAL(residual_matrix_)[offset + row] -= lambda * bv[offset + row];
    }
    REAL(residuals_)[col] = column_norm_cert(REAL(residual_matrix_), n, col);
  }
  const double applied_bound_A =
    applied_ratio_bound_cert(REAL(av_norms_), REAL(vector_norms_), k);
  const double applied_bound_B = (B_ == R_NilValue) ? 1.0 :
    applied_ratio_bound_cert(REAL(bv_norms_), REAL(vector_norms_), k);
  const double norm_A = fmax(col_bound_A, applied_bound_A);
  const double norm_B = (B_ == R_NilValue) ? 1.0 : fmax(col_bound_B, applied_bound_B);

  for (int col = 0; col < k; ++col) {
    const double lambda = REAL(values_)[col];
    const double residual = REAL(residuals_)[col];
    const double vector_norm = REAL(vector_norms_)[col];
    const double scale = fmax((norm_A + fabs(lambda) * norm_B) * fmax(vector_norm, eps), eps);
    const double backward = residual / scale;
    REAL(scale_)[col] = scale;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  double orth;
  if (B_ == R_NilValue) {
    gram_upper_dsyrk_cert(REAL(vectors_), n, k, REAL(gram_));
    orth = max_orthogonality_loss_upper_cert(REAL(gram_), k);
  } else {
    F77_CALL(dgemm)(&trans, &notrans, &k, &k, &n,
                    &one, REAL(vectors_), &n, bv, &n,
                    &zero, REAL(gram_), &k FCONE FCONE);
    orth = max_orthogonality_loss_cert(REAL(gram_), k);
  }

  const int n_out = 13;
  SEXP out_ = PROTECT(allocVector(VECSXP, n_out));
  ++protect_count;
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, ScalarReal(norm_A));
  SET_VECTOR_ELT(out_, 6, ScalarReal(norm_B));
  SET_VECTOR_ELT(out_, 7, ScalarReal(col_bound_A));
  SET_VECTOR_ELT(out_, 8, ScalarReal(applied_bound_A));
  SET_VECTOR_ELT(out_, 9, ScalarReal(col_bound_B));
  SET_VECTOR_ELT(out_, 10, ScalarReal(applied_bound_B));
  SET_VECTOR_ELT(out_, 11, vector_norms_);
  SET_VECTOR_ELT(out_, 12, bv_norms_);
  SEXP names_ = PROTECT(allocVector(STRSXP, n_out));
  ++protect_count;
  const char* out_names[] = {
    "residuals", "backward_error", "orthogonality", "scale", "converged",
    "norm_A", "norm_B", "norm_A_column_bound", "norm_A_applied_bound",
    "norm_B_column_bound", "norm_B_applied_bound", "vector_norms",
    "bv_norms"
  };
  for (int i = 0; i < n_out; ++i) {
    SET_STRING_ELT(names_, i, mkChar(out_names[i]));
  }
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(protect_count);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_svd_residuals(SEXP A_, SEXP d_,
                                              SEXP u_, SEXP v_) {
  EIGENCORE_ENTRY_BEGIN
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
    const R_xlen_t left_offset = static_cast<R_xlen_t>(col) * m;
    const R_xlen_t right_offset = static_cast<R_xlen_t>(col) * n;
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_svd_certificate(SEXP A_, SEXP d_,
                                                SEXP u_, SEXP v_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
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
  const double col_bound = max_column_norm_dense_cert(REAL(A_), m, n);

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

  double applied_bound = 0.0;
  const double norm_A = svd_certificate_finalize_cert(
    REAL(left_matrix_), REAL(right_matrix_), m, n, k, REAL(d_), REAL(u_),
    REAL(v_), col_bound, tol, REAL(left_), REAL(right_), REAL(combined_),
    REAL(scale_), REAL(backward_), LOGICAL(converged_), &applied_bound
  );
  const double scale_value = fmax(norm_A, eps);

  SEXP gram_u_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP gram_v_ = PROTECT(allocMatrix(REALSXP, k, k));
  gram_upper_dsyrk_cert(REAL(u_), m, k, REAL(gram_u_));
  gram_upper_dsyrk_cert(REAL(v_), n, k, REAL(gram_v_));
  REAL(orth_)[0] = max_orthogonality_loss_upper_cert(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss_upper_cert(REAL(gram_v_), k);
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 11));
  SET_VECTOR_ELT(out_, 0, left_);
  SET_VECTOR_ELT(out_, 1, right_);
  SET_VECTOR_ELT(out_, 2, combined_);
  SET_VECTOR_ELT(out_, 3, backward_);
  SET_VECTOR_ELT(out_, 4, orth_);
  SET_VECTOR_ELT(out_, 5, scale_);
  SET_VECTOR_ELT(out_, 6, converged_);
  SET_VECTOR_ELT(out_, 7, ScalarReal(norm_A));
  SET_VECTOR_ELT(out_, 8, ScalarReal(scale_value));
  SET_VECTOR_ELT(out_, 9, ScalarReal(col_bound));
  SET_VECTOR_ELT(out_, 10, ScalarReal(applied_bound));
  SEXP names_ = PROTECT(allocVector(STRSXP, 11));
  SET_STRING_ELT(names_, 0, mkChar("left"));
  SET_STRING_ELT(names_, 1, mkChar("right"));
  SET_STRING_ELT(names_, 2, mkChar("combined"));
  SET_STRING_ELT(names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 5, mkChar("scale"));
  SET_STRING_ELT(names_, 6, mkChar("converged"));
  SET_STRING_ELT(names_, 7, mkChar("norm_A"));
  SET_STRING_ELT(names_, 8, mkChar("scale_value"));
  SET_STRING_ELT(names_, 9, mkChar("norm_A_column_bound"));
  SET_STRING_ELT(names_, 10, mkChar("norm_A_applied_bound"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(14);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_svd_certificate_cached_av(SEXP A_, SEXP d_,
                                                          SEXP u_, SEXP v_,
                                                          SEXP av_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
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
  const double norm_A = max_column_norm_dense_cert(REAL(A_), m, n);
  return native_operator_svd_certificate_cached_av(
    &impl, eigencore_dense_apply, m, n, norm_A, d_, u_, v_, av_, tol_
  );
  EIGENCORE_ENTRY_END
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

  // Two-norm lower bound (C12): the caller's structural bound raised by
  // ||A x_j|| / ||x_j|| for the candidate vectors (A x_j from this call).
  std::vector<double> vector_norms(static_cast<size_t>(k), 0.0);
  std::vector<double> av_norms(static_cast<size_t>(k), 0.0);
  for (int col = 0; col < k; ++col) {
    const double lambda = REAL(values_)[col];
    const R_xlen_t offset = static_cast<R_xlen_t>(col) * n;
    av_norms[static_cast<size_t>(col)] = column_norm_cert(REAL(residual_matrix_), n, col);
    vector_norms[static_cast<size_t>(col)] = column_norm_cert(REAL(vectors_), n, col);
    for (int row = 0; row < n; ++row) {
      REAL(residual_matrix_)[offset + row] -= lambda * REAL(vectors_)[offset + row];
    }
    REAL(residuals_)[col] = column_norm_cert(REAL(residual_matrix_), n, col);
  }
  const double applied_bound =
    applied_ratio_bound_cert(av_norms.data(), vector_norms.data(), k);
  const double norm_lower = (R_FINITE(norm_A) && norm_A > 0.0) ? norm_A : 0.0;
  const double norm_used = fmax(norm_lower, applied_bound);
  for (int col = 0; col < k; ++col) {
    const double lambda = REAL(values_)[col];
    const double residual = REAL(residuals_)[col];
    const double vector_norm = vector_norms[static_cast<size_t>(col)];
    const double scale = fmax((norm_used + fabs(lambda)) * fmax(vector_norm, eps), eps);
    const double backward = residual / scale;
    REAL(scale_)[col] = scale;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  gram_upper_dsyrk_cert(REAL(vectors_), n, k, REAL(gram_));
  const double orth = max_orthogonality_loss_upper_cert(REAL(gram_), k);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, workspace_counters_cert(&workspace));
  SET_VECTOR_ELT(out_, 6, ScalarReal(norm_used));
  SET_VECTOR_ELT(out_, 7, ScalarReal(applied_bound));
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("scale"));
  SET_STRING_ELT(names_, 4, mkChar("converged"));
  SET_STRING_ELT(names_, 5, mkChar("workspace"));
  SET_STRING_ELT(names_, 6, mkChar("norm_A"));
  SET_STRING_ELT(names_, 7, mkChar("norm_A_applied_bound"));
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

  const double eps = DBL_EPSILON;
  const double tol = asReal(tol_);
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

  double applied_bound = 0.0;
  const double norm_used = svd_certificate_finalize_cert(
    REAL(left_matrix_), REAL(right_matrix_), m, n, k, REAL(d_), REAL(u_),
    REAL(v_), norm_A, tol, REAL(left_), REAL(right_), REAL(combined_),
    REAL(scale_), REAL(backward_), LOGICAL(converged_), &applied_bound
  );
  const double scale_value = fmax(norm_used, eps);

  SEXP gram_u_ = PROTECT(allocMatrix(REALSXP, k, k));
  SEXP gram_v_ = PROTECT(allocMatrix(REALSXP, k, k));
  gram_upper_dsyrk_cert(REAL(u_), m, k, REAL(gram_u_));
  gram_upper_dsyrk_cert(REAL(v_), n, k, REAL(gram_v_));
  REAL(orth_)[0] = max_orthogonality_loss_upper_cert(REAL(gram_u_), k);
  REAL(orth_)[1] = max_orthogonality_loss_upper_cert(REAL(gram_v_), k);
  SEXP orth_names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(orth_names_, 0, mkChar("U"));
  SET_STRING_ELT(orth_names_, 1, mkChar("V"));
  setAttrib(orth_, R_NamesSymbol, orth_names_);

  SEXP out_ = PROTECT(allocVector(VECSXP, 11));
  SET_VECTOR_ELT(out_, 0, left_);
  SET_VECTOR_ELT(out_, 1, right_);
  SET_VECTOR_ELT(out_, 2, combined_);
  SET_VECTOR_ELT(out_, 3, backward_);
  SET_VECTOR_ELT(out_, 4, orth_);
  SET_VECTOR_ELT(out_, 5, scale_);
  SET_VECTOR_ELT(out_, 6, converged_);
  SET_VECTOR_ELT(out_, 7, ScalarReal(scale_value));
  SET_VECTOR_ELT(out_, 8, workspace_counters_cert(&workspace));
  SET_VECTOR_ELT(out_, 9, ScalarReal(norm_used));
  SET_VECTOR_ELT(out_, 10, ScalarReal(applied_bound));
  SEXP names_ = PROTECT(allocVector(STRSXP, 11));
  SET_STRING_ELT(names_, 0, mkChar("left"));
  SET_STRING_ELT(names_, 1, mkChar("right"));
  SET_STRING_ELT(names_, 2, mkChar("combined"));
  SET_STRING_ELT(names_, 3, mkChar("backward_error"));
  SET_STRING_ELT(names_, 4, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 5, mkChar("scale"));
  SET_STRING_ELT(names_, 6, mkChar("converged"));
  SET_STRING_ELT(names_, 7, mkChar("scale_value"));
  SET_STRING_ELT(names_, 8, mkChar("workspace"));
  SET_STRING_ELT(names_, 9, mkChar("norm_A"));
  SET_STRING_ELT(names_, 10, mkChar("norm_A_applied_bound"));
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
  // The cached Av is validated for shape only and is NOT trusted for the
  // left residual: a stale or inconsistent cache (e.g. Avectors formed as
  // U * diag(d) from the projected SVD) would make ||A v - d u|| zero by
  // construction. Both residuals are recomputed against the operator (one
  // extra forward block apply of k columns), see C13.
  return native_operator_svd_certificate(impl, apply, m, n, norm_A,
                                         d_, u_, v_, tol_);
}

extern "C" SEXP eigencore_csc_eigen_certificate(SEXP i_, SEXP p_, SEXP x_,
                                                SEXP dim_, SEXP values_,
                                                SEXP vectors_, SEXP norm_A_,
                                                SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(values_) || !isReal(vectors_)) {
    error("invalid CSC eigen certificate inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_eigen_certificate");
  const int n = INTEGER(dim_)[0];
  if (INTEGER(dim_)[1] != n) {
    error("CSC eigen certificate requires a square operator");
  }
  CSCOperator impl = {n, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_operator_eigen_certificate(&impl, eigencore_csc_apply, n,
                                           asReal(norm_A_), values_, vectors_, tol_);
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_diagonal_eigen_certificate(SEXP x_, SEXP dim_,
                                                     SEXP unit_, SEXP values_,
                                                     SEXP vectors_, SEXP norm_A_,
                                                     SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_tridiagonal_eigen_certificate(SEXP alpha_, SEXP beta_,
                                                        SEXP values_,
                                                        SEXP vectors_,
                                                        SEXP norm_A_,
                                                        SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(alpha_) || !isReal(beta_) || !isReal(values_) || !isReal(vectors_)) {
    error("invalid tridiagonal eigen certificate inputs");
  }
  const int n = LENGTH(alpha_);
  if (n < 1) {
    error("alpha must have positive length");
  }
  if (LENGTH(beta_) < n - 1) {
    error("beta must have length at least length(alpha) - 1");
  }
  SEXP dimV = getAttrib(vectors_, R_DimSymbol);
  if (dimV == R_NilValue) {
    error("vectors must be a matrix");
  }
  const int rowsV = INTEGER(dimV)[0];
  const int k = INTEGER(dimV)[1];
  if (rowsV != n || LENGTH(values_) != k) {
    error("non-conformable tridiagonal eigen certificate inputs");
  }

  const double eps = DBL_EPSILON;
  const double norm_A_in = asReal(norm_A_);
  const double tol = asReal(tol_);

  SEXP residuals_ = PROTECT(allocVector(REALSXP, k));
  SEXP scale_ = PROTECT(allocVector(REALSXP, k));
  SEXP backward_ = PROTECT(allocVector(REALSXP, k));
  SEXP converged_ = PROTECT(allocVector(LGLSXP, k));
  SEXP gram_ = PROTECT(allocMatrix(REALSXP, k, k));

  const double* alpha = REAL(alpha_);
  const double* beta = REAL(beta_);
  const double* values = REAL(values_);
  const double* vectors = REAL(vectors_);
  // Two-norm lower bounds (C12): the largest column norm of T and
  // ||T x_j|| / ||x_j|| for the candidate vectors.
  double column_bound = 0.0;
  for (int row = 0; row < n; ++row) {
    double sq = alpha[row] * alpha[row];
    if (row > 0) {
      sq += beta[row - 1] * beta[row - 1];
    }
    if (row + 1 < n) {
      sq += beta[row] * beta[row];
    }
    if (sq > column_bound) {
      column_bound = sq;
    }
  }
  column_bound = sqrt(column_bound);
  std::vector<double> vector_norms(static_cast<size_t>(k), 0.0);
  std::vector<double> av_norms(static_cast<size_t>(k), 0.0);
  for (int col = 0; col < k; ++col) {
    const double lambda = values[col];
    const R_xlen_t offset = static_cast<R_xlen_t>(col) * n;
    long double residual_sum = 0.0L;
    long double vector_sum = 0.0L;
    long double av_sum = 0.0L;
    for (int row = 0; row < n; ++row) {
      const double v = vectors[offset + row];
      double Av = alpha[row] * v;
      if (row > 0) {
        Av += beta[row - 1] * vectors[offset + row - 1];
      }
      if (row + 1 < n) {
        Av += beta[row] * vectors[offset + row + 1];
      }
      const long double residual = Av - lambda * v;
      residual_sum += residual * residual;
      vector_sum += static_cast<long double>(v) * v;
      av_sum += static_cast<long double>(Av) * Av;
    }
    REAL(residuals_)[col] = sqrt(static_cast<double>(residual_sum));
    vector_norms[static_cast<size_t>(col)] = sqrt(static_cast<double>(vector_sum));
    av_norms[static_cast<size_t>(col)] = sqrt(static_cast<double>(av_sum));
  }
  const double applied_bound =
    applied_ratio_bound_cert(av_norms.data(), vector_norms.data(), k);
  double norm_A = (R_FINITE(norm_A_in) && norm_A_in > 0.0) ? norm_A_in : 0.0;
  norm_A = fmax(norm_A, fmax(column_bound, applied_bound));
  for (int col = 0; col < k; ++col) {
    const double lambda = values[col];
    const double residual = REAL(residuals_)[col];
    const double vector_norm = vector_norms[static_cast<size_t>(col)];
    const double scale = fmax((norm_A + fabs(lambda)) * fmax(vector_norm, eps), eps);
    const double backward = residual / scale;
    REAL(scale_)[col] = scale;
    REAL(backward_)[col] = backward;
    LOGICAL(converged_)[col] = (R_FINITE(backward) && backward <= tol) ? TRUE : FALSE;
  }

  gram_upper_dsyrk_cert(REAL(vectors_), n, k, REAL(gram_));
  const double orth = max_orthogonality_loss_upper_cert(REAL(gram_), k);

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, ScalarReal(orth));
  SET_VECTOR_ELT(out_, 3, scale_);
  SET_VECTOR_ELT(out_, 4, converged_);
  SET_VECTOR_ELT(out_, 5, ScalarReal(norm_A));
  SET_VECTOR_ELT(out_, 6, ScalarReal(column_bound));
  SET_VECTOR_ELT(out_, 7, ScalarReal(applied_bound));
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("scale"));
  SET_STRING_ELT(names_, 4, mkChar("converged"));
  SET_STRING_ELT(names_, 5, mkChar("norm_A"));
  SET_STRING_ELT(names_, 6, mkChar("norm_A_column_bound"));
  SET_STRING_ELT(names_, 7, mkChar("norm_A_applied_bound"));
  setAttrib(out_, R_NamesSymbol, names_);

  UNPROTECT(7);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_svd_certificate(SEXP i_, SEXP p_, SEXP x_,
                                              SEXP dim_, SEXP d_,
                                              SEXP u_, SEXP v_,
                                              SEXP norm_A_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_)) {
    error("invalid CSC SVD certificate inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_svd_certificate");
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_operator_svd_certificate(&impl, eigencore_csc_apply, m, n,
                                         asReal(norm_A_), d_, u_, v_, tol_);
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_svd_certificate_cached_av(SEXP i_, SEXP p_, SEXP x_,
                                                        SEXP dim_, SEXP d_,
                                                        SEXP u_, SEXP v_,
                                                        SEXP av_, SEXP norm_A_,
                                                        SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_) || !isReal(av_)) {
    error("invalid cached-Av CSC SVD certificate inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_svd_certificate_cached_av");
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  CSCOperator impl = {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)};
  return native_operator_svd_certificate_cached_av(
    &impl, eigencore_csc_apply, m, n,
    asReal(norm_A_), d_, u_, v_, av_, tol_
  );
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_diagonal_svd_certificate(SEXP x_, SEXP dim_,
                                                   SEXP unit_, SEXP d_,
                                                   SEXP u_, SEXP v_,
                                                   SEXP norm_A_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(x_) || !isInteger(dim_) || !isLogical(unit_) ||
      !isReal(d_) || !isReal(u_) || !isReal(v_)) {
    error("invalid diagonal SVD certificate inputs");
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  DiagonalOperator impl = {m, REAL(x_), LOGICAL(unit_)[0] == TRUE};
  return native_operator_svd_certificate(&impl, eigencore_diagonal_apply, m, n,
                                         asReal(norm_A_), d_, u_, v_, tol_);
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_diagonal_svd_certificate_cached_av(SEXP x_, SEXP dim_,
                                                             SEXP unit_, SEXP d_,
                                                             SEXP u_, SEXP v_,
                                                             SEXP av_,
                                                             SEXP norm_A_,
                                                             SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}
