/* Symbolic CHOLMOD analysis through the Matrix package's C API (LinkingTo:
 * Matrix). eigencore factorises sparse symmetric matrices with
 * Matrix::Cholesky() at R level; this bridge only runs cholmod_analyze() to
 * predict the fill (nnz(L)) and flop count of that factorisation, which is
 * cheap (roughly linear in nnz(A)) and lets the completeness certificate and
 * the shift-invert planner decide whether an LDL' factorisation is affordable
 * before paying for it.
 *
 * The CHOLMOD structs are laid out by the headers of the Matrix version this
 * package was compiled against. eigencore_cholmod_abi() reports that version
 * so the R layer can refuse the bridge (and fall back to a heuristic) when a
 * different Matrix ABI / SuiteSparse is loaded at run time.
 */
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Matrix/Matrix.h>
#include <Matrix/stubs.c>

SEXP eigencore_cholmod_abi(void) {
  SEXP out = PROTECT(allocVector(INTSXP, 4));
  INTEGER(out)[0] = R_MATRIX_ABI_VERSION;
  INTEGER(out)[1] = R_MATRIX_SUITESPARSE_MAJOR;
  INTEGER(out)[2] = R_MATRIX_SUITESPARSE_MINOR;
  INTEGER(out)[3] = R_MATRIX_SUITESPARSE_PATCH;
  UNPROTECT(1);
  return out;
}

/* A_ must be a dsCMatrix. Returns c(lnz, flops, ordering, status). */
SEXP eigencore_cholmod_analyze(SEXP A_) {
  cholmod_common c;
  M_cholmod_start(&c);
  c.supernodal = CHOLMOD_SIMPLICIAL;
  c.final_ll = 0;
  c.print = 0;
  /* AMD only (what Matrix::Cholesky ends up using without METIS), and no
   * postordering: the fill and flop counts do not depend on it. */
  c.nmethods = 1;
  c.method[0].ordering = CHOLMOD_AMD;
  c.postorder = 0;
  cholmod_sparse tmp;
  CHM_SP A = M_sexp_as_cholmod_sparse(&tmp, A_, FALSE, FALSE);
  CHM_FR L = M_cholmod_analyze(A, &c);
  double lnz = NA_REAL, fl = NA_REAL, ordering = NA_REAL;
  double status = (double) c.status;
  if (L != NULL) {
    lnz = c.lnz;
    fl = c.fl;
    ordering = (double) L->ordering;
    M_cholmod_free_factor(&L, &c);
  }
  M_cholmod_finish(&c);
  SEXP out = PROTECT(allocVector(REALSXP, 4));
  REAL(out)[0] = lnz;
  REAL(out)[1] = fl;
  REAL(out)[2] = ordering;
  REAL(out)[3] = status;
  SEXP names = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, mkChar("lnz"));
  SET_STRING_ELT(names, 1, mkChar("flops"));
  SET_STRING_ELT(names, 2, mkChar("ordering"));
  SET_STRING_ELT(names, 3, mkChar("status"));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}
