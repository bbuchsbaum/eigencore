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

/* C63: drop the explicit zeros that relaxed supernodal amalgamation leaves
 * in a simplicial LDL' factor converted from a supernodal one. An entry of L
 * is removed only when it is exactly zero AND outside the symbolic pattern
 * of the (unrelaxed) Cholesky factor of P A P', so the factor keeps every
 * entry a later numeric refactorisation (Matrix::update(), CHOLMOD's
 * row-oriented simplicial factorisation, which rebuilds each column's
 * pattern from the elimination tree) can write into. The symbolic pattern
 * is computed from A's pattern with Liu's elimination-tree algorithm and a
 * row-subtree traversal, O(nnz(L)). The factor is left unchanged (and the
 * function returns 0) when its layout is not the expected packed,
 * monotonic, row-sorted simplicial LDL' or the patterns disagree; it
 * returns the number of entries dropped otherwise (-1 on allocation
 * failure). */
static long eigencore_ldl_drop_relaxed_zeros(cholmod_sparse* A, cholmod_factor* L) {
  if (L->is_super || L->is_ll || !L->is_monotonic || L->xtype != CHOLMOD_REAL ||
      L->itype != CHOLMOD_INT || L->dtype != CHOLMOD_DOUBLE ||
      A->itype != CHOLMOD_INT || A->stype == 0 || L->Perm == NULL ||
      A->nrow != L->n || A->ncol != L->n) {
    return 0;
  }
  const int n = (int) L->n;
  int* Lp = (int*) L->p;
  int* Li = (int*) L->i;
  double* Lx = (double*) L->x;
  int* Lnz = (int*) L->nz;
  const int* Perm = (const int*) L->Perm;
  const int* Ap = (const int*) A->p;
  const int* Ai = (const int*) A->i;
  const int* Anz = A->packed ? NULL : (const int*) A->nz;
  if (n < 2 || Lnz == NULL) {
    return 0;
  }
  /* Columns must be packed in order with sorted row indices and the
   * diagonal first. */
  for (int j = 0; j < n; ++j) {
    if (Lp[j] + Lnz[j] != Lp[j + 1] || Lnz[j] < 1 || Li[Lp[j]] != j) {
      return 0;
    }
    for (int q = Lp[j] + 1; q < Lp[j + 1]; ++q) {
      if (Li[q] <= Li[q - 1]) {
        return 0;
      }
    }
  }
  long status = -1;
  const size_t nn = (size_t) n;
  int* pinv = (int*) malloc(nn * sizeof(int));
  int* rcount = (int*) calloc(nn + 1, sizeof(int));
  int* parent = (int*) malloc(nn * sizeof(int));
  int* ancestor = (int*) malloc(nn * sizeof(int));
  int* flag = (int*) malloc(nn * sizeof(int));
  int* cursor = (int*) malloc(nn * sizeof(int));
  int* adj = NULL;
  unsigned char* keep = NULL;
  if (pinv == NULL || rcount == NULL || parent == NULL || ancestor == NULL ||
      flag == NULL || cursor == NULL) {
    goto done;
  }
  for (int k = 0; k < n; ++k) {
    if (Perm[k] < 0 || Perm[k] >= n) {
      status = 0;
      goto done;
    }
    pinv[Perm[k]] = k;
  }
  /* Rows of the strictly lower triangle of P A P' (only one triangle of A
   * is referenced: stype > 0 upper, < 0 lower). */
  for (int c = 0; c < n; ++c) {
    const int end = Anz ? Ap[c] + Anz[c] : Ap[c + 1];
    for (int q = Ap[c]; q < end; ++q) {
      const int r = Ai[q];
      if (r == c || (A->stype > 0 ? r > c : r < c)) continue;
      const int pr = pinv[r], pc = pinv[c];
      ++rcount[(pr > pc ? pr : pc) + 1];
    }
  }
  for (int k = 0; k < n; ++k) rcount[k + 1] += rcount[k];
  adj = (int*) malloc((size_t) (rcount[n] > 0 ? rcount[n] : 1) * sizeof(int));
  keep = (unsigned char*) calloc((size_t) Lp[n] + 1, 1);
  if (adj == NULL || keep == NULL) {
    goto done;
  }
  for (int k = 0; k < n; ++k) cursor[k] = rcount[k];
  for (int c = 0; c < n; ++c) {
    const int end = Anz ? Ap[c] + Anz[c] : Ap[c + 1];
    for (int q = Ap[c]; q < end; ++q) {
      const int r = Ai[q];
      if (r == c || (A->stype > 0 ? r > c : r < c)) continue;
      const int pr = pinv[r], pc = pinv[c];
      const int k = pr > pc ? pr : pc;
      adj[cursor[k]++] = pr > pc ? pc : pr;
    }
  }
  /* Elimination tree (Liu, path compression). */
  for (int k = 0; k < n; ++k) {
    parent[k] = -1;
    ancestor[k] = -1;
    for (int q = rcount[k]; q < rcount[k + 1]; ++q) {
      int i = adj[q];
      while (ancestor[i] != -1 && ancestor[i] != k) {
        const int next = ancestor[i];
        ancestor[i] = k;
        i = next;
      }
      if (ancestor[i] == -1) {
        ancestor[i] = k;
        parent[i] = k;
      }
    }
  }
  /* Row subtrees: row k of L has entries in the columns on the etree paths
   * from each i in adj[k] up to k. Rows arrive in increasing k, so a cursor
   * per column finds each symbolic entry in L's sorted column. */
  for (int j = 0; j < n; ++j) {
    keep[Lp[j]] = 1;              /* diagonal (D) */
    cursor[j] = Lp[j] + 1;
    flag[j] = -1;
  }
  for (int k = 0; k < n; ++k) {
    flag[k] = k;
    for (int q = rcount[k]; q < rcount[k + 1]; ++q) {
      for (int j = adj[q]; flag[j] != k; j = parent[j]) {
        if (j < 0 || j >= k) {
          status = 0;
          goto done;
        }
        flag[j] = k;
        int pos = cursor[j];
        while (pos < Lp[j + 1] && Li[pos] < k) ++pos;
        if (pos >= Lp[j + 1] || Li[pos] != k) {
          status = 0;               /* symbolic entry missing from L */
          goto done;
        }
        keep[pos] = 1;
        cursor[j] = pos + 1;
      }
    }
  }
  /* Compact in place (columns are packed and in order, so writes never
   * overtake reads). */
  {
    long dropped = 0;
    int dst = 0;
    for (int j = 0; j < n; ++j) {
      const int begin = Lp[j], end = Lp[j + 1];
      Lp[j] = dst;
      for (int q = begin; q < end; ++q) {
        if (!keep[q] && Lx[q] == 0.0) {
          ++dropped;
          continue;
        }
        Li[dst] = Li[q];
        Lx[dst] = Lx[q];
        ++dst;
      }
      Lnz[j] = dst - Lp[j];
      if (L->ColCount != NULL) ((int*) L->ColCount)[j] = Lnz[j];
    }
    Lp[n] = dst;
    /* The arrays keep their allocation; nzmax only bounds what is copied
     * into the R object (and CHOLMOD's memory statistics of this local
     * Common). */
    L->nzmax = (size_t) dst;
    status = dropped;
  }
done:
  free(pinv);
  free(rcount);
  free(parent);
  free(ancestor);
  free(flag);
  free(cursor);
  free(adj);
  free(keep);
  return status;
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
  /* Default relaxed amalgamation (without it the supernodal factorisation
   * is no faster than the simplicial one). The explicit zeros it stores
   * (about a third more entries on a 2-D grid) are dropped again after the
   * conversion to simplicial LDL' (C63). */
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
  if (eigencore_ldl_drop_relaxed_zeros(A, L) < 0) {
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
