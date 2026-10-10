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
#include <stdlib.h>
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

/* Native CHOLMOD solve for the shift-invert kernel (tranche 5, item 4).
 *
 * A solver handle wraps a Matrix CHMfactor (simplicial or supernodal, LL' or
 * LDL') without copying it, keeps its own cholmod_common and reuses the
 * cholmod_solve2() output / workspace buffers across calls, so a Lanczos step
 * costs one triangular solve pair and no allocation. The R object owning the
 * factor must outlive the handle (the composite kernel's external pointer
 * protects its spec, which holds the factor). Only built after the R layer
 * has checked the run-time Matrix ABI (cholmod_bridge_available()).
 */
typedef struct {
  cholmod_common c;
  cholmod_factor L;
  CHM_DN X, Y, E;
  int n;
} eigencore_cholmod_solver;

void* eigencore_cholmod_solver_new(SEXP factor) {
  eigencore_cholmod_solver* s =
    (eigencore_cholmod_solver*) calloc(1, sizeof(eigencore_cholmod_solver));
  if (s == NULL) {
    return NULL;
  }
  M_cholmod_start(&s->c);
  s->c.print = 0;
  s->c.error_handler = NULL;
  if (M_sexp_as_cholmod_factor(&s->L, factor) == NULL || s->L.n < 1 ||
      s->L.minor < s->L.n) {
    M_cholmod_finish(&s->c);
    free(s);
    return NULL;
  }
  s->n = (int) s->L.n;
  return s;
}

void eigencore_cholmod_solver_free(void* handle) {
  eigencore_cholmod_solver* s = (eigencore_cholmod_solver*) handle;
  if (s == NULL) {
    return;
  }
  if (s->X != NULL) M_cholmod_free_dense(&s->X, &s->c);
  if (s->Y != NULL) M_cholmod_free_dense(&s->Y, &s->c);
  if (s->E != NULL) M_cholmod_free_dense(&s->E, &s->c);
  M_cholmod_finish(&s->c);
  free(s);
}

/* out[, j] = alpha * A^{-1} X[, j] + beta * out[, j] for ncol columns.
 * Returns 0 on success, -3 on a CHOLMOD failure. */
int eigencore_cholmod_solver_apply(void* handle, int ncol, const double* X,
                                   int ldx, double alpha, double beta,
                                   double* out, int ldo) {
  eigencore_cholmod_solver* s = (eigencore_cholmod_solver*) handle;
  if (s == NULL || ncol < 1) {
    return ncol < 1 ? 0 : -1;
  }
  const int n = s->n;
  for (int col = 0; col < ncol; ++col) {
    cholmod_dense b;
    M_numeric_as_cholmod_dense(&b, (double*) (X + (size_t) col * ldx), n, 1);
    if (!M_cholmod_solve2(CHOLMOD_A, &s->L, &b, &s->X, &s->Y, &s->E, &s->c) ||
        s->X == NULL) {
      return -3;
    }
    const double* x = (const double*) s->X->x;
    double* y = out + (size_t) col * ldo;
    if (beta == 0.0) {
      for (int i = 0; i < n; ++i) y[i] = alpha * x[i];
    } else {
      for (int i = 0; i < n; ++i) y[i] = alpha * x[i] + beta * y[i];
    }
  }
  return 0;
}

/* Positive definite shift: try a supernodal LL' factorisation of
 * A - sigma I (A a dsCMatrix) and, when it succeeds (which proves
 * A - sigma I positive definite, i.e. inertia (0, 0, n)), convert it in place
 * to a simplicial LDL' factor returned as a Matrix dCHMsimpl. Supernodal
 * numeric factorisation runs on BLAS and is several times faster than the
 * simplicial LDL' that an indefinite shift needs; the simplicial LDL' form
 * keeps solves fast and lets Matrix::update() reuse its symbolic analysis
 * for later (indefinite) shifts. Returns NULL when A - sigma I is not
 * positive definite or CHOLMOD fails. */
SEXP eigencore_cholmod_spd_ldl(SEXP A_, SEXP sigma_) {
  cholmod_common c;
  M_cholmod_start(&c);
  c.print = 0;
  c.error_handler = NULL;
  c.supernodal = CHOLMOD_SUPERNODAL;
  c.final_ll = 1;
  /* Default relaxed amalgamation: it stores explicit zeros (about a third
   * more entries on a 2-D grid, so each later solve costs ~1.3x), but
   * without it the supernodal factorisation is no faster than the
   * simplicial one. */
  cholmod_sparse tmp;
  CHM_SP A = M_sexp_as_cholmod_sparse(&tmp, A_, FALSE, FALSE);
  CHM_FR L = M_cholmod_analyze(A, &c);
  if (L == NULL) {
    M_cholmod_finish(&c);
    return R_NilValue;
  }
  double beta[2];
  beta[0] = -asReal(sigma_);
  beta[1] = 0.0;
  int ok = M_cholmod_factorize_p(A, beta, NULL, 0, L, &c);
  if (!ok || c.status != CHOLMOD_OK || L->minor < L->n) {
    M_cholmod_free_factor(&L, &c);
    M_cholmod_finish(&c);
    return R_NilValue;
  }
  ok = M_cholmod_change_factor(CHOLMOD_REAL, FALSE, FALSE, TRUE, TRUE, L, &c);
  if (!ok || c.status != CHOLMOD_OK || L->is_super || L->is_ll) {
    M_cholmod_free_factor(&L, &c);
    M_cholmod_finish(&c);
    return R_NilValue;
  }
  SEXP out = PROTECT(M_cholmod_factor_as_sexp(L, 0));
  M_cholmod_free_factor(&L, &c);
  M_cholmod_finish(&c);
  UNPROTECT(1);
  return out;
}
