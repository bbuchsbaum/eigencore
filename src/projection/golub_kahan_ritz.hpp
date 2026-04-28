#ifndef EIGENCORE_PROJECTION_GOLUB_KAHAN_RITZ_HPP
#define EIGENCORE_PROJECTION_GOLUB_KAHAN_RITZ_HPP

static SEXP eigencore_golub_kahan_ritz_from_ptr(const double* U,
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

  SEXP bd_ = PROTECT(eigencore_bidiagonal_svd(alpha_active_, beta_active_));
  SEXP d_all_ = VECTOR_ELT(bd_, 0);
  SEXP u_small_ = VECTOR_ELT(bd_, 1);
  SEXP v_small_ = VECTOR_ELT(bd_, 2);

  std::vector<int> selected(static_cast<size_t>(rank));
  const int count = selected_ritz_indices(REAL(d_all_), p, rank, target_kind, selected.data());

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
                    &one, U, &m, REAL(u_sel_), &p,
                    &zero, REAL(u_), &m FCONE FCONE);
    F77_CALL(dgemm)(&notrans, &notrans, &n, &count, &p,
                    &one, V, &n, REAL(v_sel_), &p,
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
}

#endif
