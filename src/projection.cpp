#include <cmath>
#include <cstring>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "eigencore_lapack_compat.h"
#include "eigencore_common.h"
#include "golub_kahan_ritz.h"

extern "C" SEXP eigencore_bidiagonal_svd(SEXP alpha_, SEXP beta_);

extern "C" SEXP eigencore_rayleigh_ritz_symmetric(SEXP A_, SEXP Q_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

static bool ritz_value_better_projection(double candidate,
                                         double incumbent,
                                         int target_kind) {
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

static int selected_ritz_indices_projection(const double* values,
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
      if (best < 0 ||
          ritz_value_better_projection(values[j], values[best], target_kind)) {
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

// Thin SVD of an m x p block via LAPACK dgesvd, returning the s = min(m, p)
// singular values, the full m x s left factor and the s x p right factor.
static void ritz_dense_thin_svd(double* B, int m, int p, double* d_all,
                                double* u_all, double* vt_all) {
  const int s = (m < p) ? m : p;
  char jobu = 'S';
  char jobvt = 'S';
  int info = 0;
  int lwork = -1;
  double work_query = 0.0;
  F77_CALL(dgesvd)(&jobu, &jobvt, &m, &p, B, &m,
                   d_all, u_all, &m, vt_all, &s,
                   &work_query, &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgesvd workspace query failed with info=%d", info);
  }
  lwork = static_cast<int>(work_query);
  if (lwork < 1) {
    lwork = 1;
  }
  std::vector<double> work(static_cast<size_t>(lwork), 0.0);
  F77_CALL(dgesvd)(&jobu, &jobvt, &m, &p, B, &m,
                   d_all, u_all, &m, vt_all, &s,
                   work.data(), &lwork, &info FCONE FCONE);
  if (info != 0) {
    error("LAPACK dgesvd failed with info=%d", info);
  }
}

SEXP eigencore_block_golub_kahan_ritz_from_ptr(const double* V,
                                               int n,
                                               const double* AV,
                                               int m,
                                               int p,
                                               int rank,
                                               int target_kind) {
  if (rank < 1) {
    error("rank must be positive");
  }
  const int s = (m < p) ? m : p;
  if (rank > s) {
    rank = s;
  }

  // Ritz extraction from the projected block AV = U_k S W^T. For the usual
  // tall case (m >= p) factor AV = Q R once, take the SVD of the small p x p
  // R, and apply Q only to the selected left singular vectors (dormqr on
  // count columns) instead of forming all p columns of an m x p U with a
  // full dgesvd(jobu = 'S') of AV.
  std::vector<double> B(static_cast<size_t>(m) * p, 0.0);
  std::memcpy(B.data(), AV, sizeof(double) * static_cast<size_t>(m) * p);
  std::vector<double> d_all(static_cast<size_t>(s), 0.0);
  std::vector<double> vt_all(static_cast<size_t>(s) * p, 0.0);
  std::vector<double> u_small;   // tall: p x p left factor of R
  std::vector<double> u_all;     // wide: m x s left factor of AV
  std::vector<double> tau;
  const bool tall = (m >= p);
  if (tall) {
    tau.assign(static_cast<size_t>(p), 0.0);
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    F77_CALL(dgeqrf)(&m, &p, B.data(), &m, tau.data(), &work_query, &lwork, &info);
    if (info != 0) {
      error("LAPACK dgeqrf workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    if (lwork < 1) {
      lwork = 1;
    }
    std::vector<double> work(static_cast<size_t>(lwork), 0.0);
    F77_CALL(dgeqrf)(&m, &p, B.data(), &m, tau.data(), work.data(), &lwork, &info);
    if (info != 0) {
      error("LAPACK dgeqrf failed with info=%d", info);
    }
    std::vector<double> R(static_cast<size_t>(p) * p, 0.0);
    for (int col = 0; col < p; ++col) {
      for (int row = 0; row <= col; ++row) {
        R[row + static_cast<size_t>(col) * p] = B[row + static_cast<size_t>(col) * m];
      }
    }
    u_small.assign(static_cast<size_t>(p) * p, 0.0);
    ritz_dense_thin_svd(R.data(), p, p, d_all.data(), u_small.data(), vt_all.data());
  } else {
    u_all.assign(static_cast<size_t>(m) * s, 0.0);
    ritz_dense_thin_svd(B.data(), m, p, d_all.data(), u_all.data(), vt_all.data());
  }

  std::vector<int> selected(static_cast<size_t>(rank));
  const int count = selected_ritz_indices_projection(
    d_all.data(), s, rank, target_kind, selected.data()
  );

  SEXP d_ = PROTECT(allocVector(REALSXP, count));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, count));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, count));
  SEXP avectors_ = PROTECT(allocMatrix(REALSXP, m, count));
  SEXP coeff_ = PROTECT(allocMatrix(REALSXP, p, count));

  if (tall) {
    std::memset(REAL(u_), 0, sizeof(double) * static_cast<size_t>(m) * count);
  }
  for (int col = 0; col < count; ++col) {
    const int idx = selected[static_cast<size_t>(col)];
    REAL(d_)[col] = d_all[static_cast<size_t>(idx)];
    if (tall) {
      for (int row = 0; row < p; ++row) {
        REAL(u_)[row + static_cast<size_t>(col) * m] =
          u_small[row + static_cast<size_t>(idx) * p];
      }
    } else {
      for (int row = 0; row < m; ++row) {
        REAL(u_)[row + static_cast<size_t>(col) * m] =
          u_all[row + static_cast<size_t>(idx) * m];
      }
    }
    for (int row = 0; row < p; ++row) {
      REAL(coeff_)[row + static_cast<size_t>(col) * p] =
        vt_all[idx + static_cast<size_t>(row) * s];
    }
  }

  if (tall && count > 0) {
    // u_ holds [U_R(:, selected); 0]; apply Q from the left.
    char side = 'L';
    char trans = 'N';
    int info = 0;
    int lwork = -1;
    double work_query = 0.0;
    int ncols = count;
    F77_CALL(dormqr)(&side, &trans, &m, &ncols, &p, B.data(), &m, tau.data(),
                     REAL(u_), &m, &work_query, &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dormqr workspace query failed with info=%d", info);
    }
    lwork = static_cast<int>(work_query);
    if (lwork < 1) {
      lwork = 1;
    }
    std::vector<double> work(static_cast<size_t>(lwork), 0.0);
    F77_CALL(dormqr)(&side, &trans, &m, &ncols, &p, B.data(), &m, tau.data(),
                     REAL(u_), &m, work.data(), &lwork, &info FCONE FCONE);
    if (info != 0) {
      error("LAPACK dormqr failed with info=%d", info);
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
  EIGENCORE_ENTRY_BEGIN
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
  return eigencore_block_golub_kahan_ritz_from_ptr(
    REAL(V_), n, REAL(AV_), m, p,
    asInteger(rank_), asInteger(target_kind_)
  );
  EIGENCORE_ENTRY_END
}

SEXP eigencore_golub_kahan_ritz_from_ptr(const double* U,
                                         const double* V,
                                         int m,
                                         int n,
                                         int p,
                                         const double* alpha,
                                         const double* beta,
                                         int rank,
                                         int target_kind) {
  if (rank < 1) {
    error("rank must be positive");
  }
  if (p < 1) {
    error("active Golub-Kahan Ritz columns must be positive");
  }
  if (rank > p) {
    rank = p;
  }

  SEXP alpha_active_ = PROTECT(allocVector(REALSXP, p));
  SEXP beta_active_ = PROTECT(allocVector(REALSXP, p));
  std::memcpy(REAL(alpha_active_), alpha, sizeof(double) * static_cast<size_t>(p));
  std::memcpy(REAL(beta_active_), beta, sizeof(double) * static_cast<size_t>(p));

  // Nested .Call entry point: its R error (a longjmp) is converted back into a
  // C++ exception so callers' C++ objects are destroyed (C10).
  SEXP bd_ = PROTECT(eigencore_unwind_protect([&]() {
    return eigencore_bidiagonal_svd(alpha_active_, beta_active_);
  }));
  SEXP d_all_ = VECTOR_ELT(bd_, 0);
  SEXP u_small_ = VECTOR_ELT(bd_, 1);
  SEXP v_small_ = VECTOR_ELT(bd_, 2);

  std::vector<int> selected(static_cast<size_t>(rank));
  const int count = selected_ritz_indices_projection(
    REAL(d_all_), p, rank, target_kind, selected.data()
  );

  SEXP d_ = PROTECT(allocVector(REALSXP, count));
  SEXP u_sel_ = PROTECT(allocMatrix(REALSXP, p, count));
  SEXP v_sel_ = PROTECT(allocMatrix(REALSXP, p, count));
  SEXP u_ = PROTECT(allocMatrix(REALSXP, m, count));
  SEXP v_ = PROTECT(allocMatrix(REALSXP, n, count));

  for (int col = 0; col < count; ++col) {
    const int idx = selected[static_cast<size_t>(col)];
    REAL(d_)[col] = REAL(d_all_)[idx];
    for (int row = 0; row < p; ++row) {
      REAL(u_sel_)[row + col * p] = REAL(u_small_)[row + idx * p];
      REAL(v_sel_)[row + col * p] = REAL(v_small_)[row + idx * p];
    }
  }

  if (count > 0) {
    const char notrans = 'N';
    const double one = 1.0;
    const double zero = 0.0;
    F77_CALL(dgemm)(&notrans, &notrans, &m, &count, &p,
                    &one, const_cast<double*>(U), &m, REAL(u_sel_), &p,
                    &zero, REAL(u_), &m FCONE FCONE);
    F77_CALL(dgemm)(&notrans, &notrans, &n, &count, &p,
                    &one, const_cast<double*>(V), &n, REAL(v_sel_), &p,
                    &zero, REAL(v_), &n FCONE FCONE);
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

  UNPROTECT(10);
  return out_;
}

extern "C" SEXP eigencore_golub_kahan_ritz(SEXP U_, SEXP V_, SEXP alpha_,
                                           SEXP beta_, SEXP rank_,
                                           SEXP target_kind_, SEXP active_p_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(U_) || !isReal(V_) || !isReal(alpha_) || !isReal(beta_)) {
    error("U, V, alpha, and beta must be double");
  }
  SEXP dimU = getAttrib(U_, R_DimSymbol);
  SEXP dimV = getAttrib(V_, R_DimSymbol);
  if (dimU == R_NilValue || dimV == R_NilValue) {
    error("U and V must be matrices");
  }
  const int m = INTEGER(dimU)[0];
  const int stored_p = INTEGER(dimU)[1];
  const int n = INTEGER(dimV)[0];
  int p = asInteger(active_p_);
  if (p < 1 || p > stored_p) {
    error("active Golub-Kahan Ritz columns must be between 1 and ncol(U)");
  }
  if (INTEGER(dimV)[1] != stored_p || LENGTH(alpha_) < p || LENGTH(beta_) < p) {
    error("non-conformable Golub-Kahan Ritz inputs");
  }

  return eigencore_golub_kahan_ritz_from_ptr(
    REAL(U_), REAL(V_), m, n, p, REAL(alpha_), REAL(beta_),
    asInteger(rank_), asInteger(target_kind_)
  );
  EIGENCORE_ENTRY_END
}
