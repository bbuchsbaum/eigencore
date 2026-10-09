#include <algorithm>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <climits>
#include <memory>
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

// CSC structure validation for every native entry point that dereferences
// borrowed dgCMatrix slots. A dgCMatrix whose slots were edited after
// construction (or a hand-built list) would otherwise index out of bounds.
// Checks: dim is two non-negative ints, i/p integer and x double, length(p) ==
// n + 1, p[0] == 0, p non-decreasing, p[n] == length(x) == length(i), and
// 0 <= i < m. Pass R_NilValue for i_ to skip the row-index checks (callers
// that never read i). Signals an R error naming `context` on failure.
extern "C" void eigencore_validate_csc_structure(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_, const char* context) {
  if (!isInteger(dim_) || XLENGTH(dim_) != 2) {
    error("invalid CSC structure (%s): Dim must be an integer vector of length 2", context);
  }
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  if (m == NA_INTEGER || n == NA_INTEGER || m < 0 || n < 0) {
    error("invalid CSC structure (%s): negative or missing dimensions", context);
  }
  if (!isInteger(p_) || !isReal(x_) || (i_ != R_NilValue && !isInteger(i_))) {
    error("invalid CSC structure (%s): i and p must be integer and x double", context);
  }
  if (XLENGTH(p_) != static_cast<R_xlen_t>(n) + 1) {
    error("invalid CSC structure (%s): length(p) must equal ncol + 1", context);
  }
  const int* p = INTEGER(p_);
  if (p[0] != 0) {
    error("invalid CSC structure (%s): p[1] must be 0", context);
  }
  for (int col = 0; col < n; ++col) {
    if (p[col + 1] == NA_INTEGER || p[col + 1] < p[col]) {
      error("invalid CSC structure (%s): column pointers must be non-decreasing", context);
    }
  }
  const R_xlen_t nnz = static_cast<R_xlen_t>(p[n]);
  if (nnz != XLENGTH(x_)) {
    error("invalid CSC structure (%s): p[ncol + 1] must equal length(x)", context);
  }
  if (i_ == R_NilValue) {
    return;
  }
  if (XLENGTH(i_) != nnz) {
    error("invalid CSC structure (%s): length(i) must equal length(x)", context);
  }
  const int* idx = INTEGER(i_);
  // Unsigned comparison folds the i < 0 (including NA) and i >= m tests.
  const unsigned int bound = static_cast<unsigned int>(m);
  unsigned int bad = 0;
  for (R_xlen_t pos = 0; pos < nnz; ++pos) {
    bad |= (static_cast<unsigned int>(idx[pos]) >= bound) ? 1u : 0u;
  }
  if (bad) {
    error("invalid CSC structure (%s): row indices must lie in [0, nrow)", context);
  }
}

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

  if (block_cols_i == 1) {
    // Single-vector Krylov steps: dgemv avoids dgemm's packing and
    // threading overhead, which dominates for one column.
    const int rows = static_cast<int>(dense->rows);
    const int cols = static_cast<int>(dense->cols);
    const int inc = 1;
    F77_CALL(dgemv)(&transa, &rows, &cols, &alpha, dense->values, &lda,
                    const_cast<double*>(X), &inc, &beta_blas, Y, &inc FCONE);
    return 0;
  }
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

// ---------------------------------------------------------------------------
// Thread count (P8)
// ---------------------------------------------------------------------------
//
// The package default is computed in R at load time (R/threads.R: 1 under
// R CMD check or when _R_CHECK_LIMIT_CORES_ is set, otherwise OMP_NUM_THREADS
// or the processor count capped at 8) and stored with
// eigencore_set_default_threads(). getOption("eigencore.threads") overrides it;
// every .Call entry re-reads the option once (eigencore_refresh_thread_count,
// called from eigencore_call), so kernels below only read g_eigencore_threads.

static int g_eigencore_default_threads = 1;
static int g_eigencore_threads = 1;
// Memory cap for the cached CSR copy used by parallel forward CSC applies;
// getOption("eigencore.csr_cache_mb", 4096).
static const double kEigencoreDefaultCsrCacheMb = 4096.0;
static double g_eigencore_csr_cache_bytes = kEigencoreDefaultCsrCacheMb * 1048576.0;
static const int kEigencoreMaxThreads = 256;

static int eigencore_processor_count() {
#ifdef _OPENMP
  const int procs = omp_get_num_procs();
  return procs > 0 ? procs : 1;
#else
  return 1;
#endif
}

#ifdef _OPENMP
static int eigencore_sanitize_thread_count(double value, int fallback) {
  if (!R_FINITE(value) || value < 1.0) {
    return fallback;
  }
  if (value > static_cast<double>(kEigencoreMaxThreads)) {
    return kEigencoreMaxThreads;
  }
  return static_cast<int>(value);
}
#endif

extern "C" int eigencore_thread_count(void) {
  return g_eigencore_threads;
}

extern "C" void eigencore_refresh_thread_count(void) {
#ifdef _OPENMP
  static SEXP option_symbol = nullptr;
  if (option_symbol == nullptr) {
    option_symbol = Rf_install("eigencore.threads");
  }
  const SEXP option = Rf_GetOption1(option_symbol);
  int threads = g_eigencore_default_threads;
  if (option != R_NilValue && XLENGTH(option) >= 1) {
    if (TYPEOF(option) == REALSXP) {
      threads = eigencore_sanitize_thread_count(REAL(option)[0], threads);
    } else if (TYPEOF(option) == INTSXP && INTEGER(option)[0] != NA_INTEGER) {
      threads = eigencore_sanitize_thread_count(
        static_cast<double>(INTEGER(option)[0]), threads);
    }
  }
  g_eigencore_threads = threads;
  static SEXP csr_symbol = nullptr;
  if (csr_symbol == nullptr) {
    csr_symbol = Rf_install("eigencore.csr_cache_mb");
  }
  const SEXP csr_option = Rf_GetOption1(csr_symbol);
  double csr_mb = kEigencoreDefaultCsrCacheMb;
  if ((TYPEOF(csr_option) == REALSXP || TYPEOF(csr_option) == INTSXP) &&
      XLENGTH(csr_option) >= 1) {
    const double value = TYPEOF(csr_option) == REALSXP ? REAL(csr_option)[0] :
      (INTEGER(csr_option)[0] == NA_INTEGER ? NA_REAL :
         static_cast<double>(INTEGER(csr_option)[0]));
    if (!ISNAN(value) && value >= 0.0) {
      csr_mb = value;
    }
  }
  g_eigencore_csr_cache_bytes = csr_mb * 1048576.0;
#else
  g_eigencore_threads = 1;
#endif
}

// BLAS thread coordination. OpenBLAS (pthreads build) and FlexiBLAS keep their
// worker threads spinning for a while after every call; an OpenMP region that
// starts meanwhile runs time-sliced against them and can be several times
// slower than serial. The first multithreaded sparse kernel of a .Call
// therefore switches such a BLAS to one thread (found at run time with
// dlsym, so there is no link dependency) and eigencore_call_leave() restores
// the previous count when the outermost .Call returns or raises. While BLAS is
// single-threaded, the Lanczos reorthogonalisation runs its own OpenMP
// kernels (eigencore_reorth_threads()). An OpenMP-built OpenBLAS, MKL and BLIS
// are left alone and keep doing the dense work with their own threads.
#if defined(_OPENMP) && !defined(_WIN32)
#include <dlfcn.h>
#define EIGENCORE_HAVE_BLAS_QUIESCE 1
#endif

enum EigencoreBlasKind {
  EIGENCORE_BLAS_SERIAL = 0,        // no known threading control (reference)
  EIGENCORE_BLAS_CONTROLLABLE = 1,  // OpenBLAS pthreads / FlexiBLAS
  EIGENCORE_BLAS_THREADED = 2       // threaded, left alone (OpenMP OpenBLAS, MKL, BLIS)
};

#ifdef EIGENCORE_HAVE_BLAS_QUIESCE
typedef void (*eigencore_blas_set_threads_fn)(int);
typedef int (*eigencore_blas_get_threads_fn)(void);
static bool g_blas_probed = false;
static int g_blas_kind = EIGENCORE_BLAS_SERIAL;
static eigencore_blas_set_threads_fn g_blas_set_threads = nullptr;
static eigencore_blas_get_threads_fn g_blas_get_threads = nullptr;
static bool g_blas_quiesced = false;
static int g_blas_saved_threads = 0;

static void eigencore_blas_probe() {
  if (g_blas_probed) {
    return;
  }
  g_blas_probed = true;
  void* set = dlsym(RTLD_DEFAULT, "openblas_set_num_threads");
  void* get = dlsym(RTLD_DEFAULT, "openblas_get_num_threads");
  if (set != nullptr && get != nullptr) {
    // openblas_get_parallel(): 0 sequential, 1 pthreads, 2 OpenMP.
    void* parallel = dlsym(RTLD_DEFAULT, "openblas_get_parallel");
    const int mode = parallel != nullptr ?
      reinterpret_cast<int (*)(void)>(parallel)() : 1;
    if (mode == 2) {
      g_blas_kind = EIGENCORE_BLAS_THREADED;
    } else if (mode == 1) {
      g_blas_kind = EIGENCORE_BLAS_CONTROLLABLE;
      g_blas_set_threads = reinterpret_cast<eigencore_blas_set_threads_fn>(set);
      g_blas_get_threads = reinterpret_cast<eigencore_blas_get_threads_fn>(get);
    }
    return;
  }
  set = dlsym(RTLD_DEFAULT, "flexiblas_set_num_threads");
  get = dlsym(RTLD_DEFAULT, "flexiblas_get_num_threads");
  if (set != nullptr && get != nullptr) {
    g_blas_kind = EIGENCORE_BLAS_CONTROLLABLE;
    g_blas_set_threads = reinterpret_cast<eigencore_blas_set_threads_fn>(set);
    g_blas_get_threads = reinterpret_cast<eigencore_blas_get_threads_fn>(get);
    return;
  }
  if (dlsym(RTLD_DEFAULT, "MKL_Get_Max_Threads") != nullptr ||
      dlsym(RTLD_DEFAULT, "mkl_get_max_threads") != nullptr ||
      dlsym(RTLD_DEFAULT, "bli_thread_get_num_threads") != nullptr) {
    g_blas_kind = EIGENCORE_BLAS_THREADED;
  }
}
#endif

// Called on the main thread before a multithreaded sparse kernel.
static void eigencore_blas_quiesce() {
#ifdef EIGENCORE_HAVE_BLAS_QUIESCE
  if (g_blas_quiesced) {
    return;
  }
  eigencore_blas_probe();
  if (g_blas_set_threads == nullptr) {
    return;
  }
  const int current = g_blas_get_threads();
  if (current > 1) {
    g_blas_saved_threads = current;
    g_blas_set_threads(1);
    g_blas_quiesced = true;
  }
#endif
}

// True when a multithreaded spinning-thread BLAS is active (and not already
// quiesced by this call).
static bool eigencore_blas_busy() {
#ifdef EIGENCORE_HAVE_BLAS_QUIESCE
  if (g_blas_quiesced) {
    return false;
  }
  eigencore_blas_probe();
  return g_blas_get_threads != nullptr && g_blas_get_threads() > 1;
#else
  return false;
#endif
}

// Threads for eigencore's own OpenMP dense helpers (Lanczos
// reorthogonalisation): the thread count while BLAS runs single-threaded
// (quiesced by this call, configured to one thread, or a serial reference
// BLAS); 1 when a threaded BLAS should do the work.
extern "C" int eigencore_reorth_threads(void) {
  const int threads = g_eigencore_threads;
  if (threads <= 1) {
    return 1;
  }
#ifdef EIGENCORE_HAVE_BLAS_QUIESCE
  if (g_blas_quiesced) {
    return threads;
  }
  eigencore_blas_probe();
  if (g_blas_kind == EIGENCORE_BLAS_THREADED) {
    return 1;
  }
  if (g_blas_kind == EIGENCORE_BLAS_CONTROLLABLE) {
    return g_blas_get_threads() > 1 ? 1 : threads;
  }
  return threads;
#else
  return threads;
#endif
}

// .Call nesting depth: BLAS threads are restored only when the outermost call
// leaves, so a native solver driving an R callback operator (whose applies
// are nested .Calls) keeps BLAS quiet for the whole solve. Every exit of
// eigencore_call (return, C++ exception, interrupt, R unwind) runs
// eigencore_call_leave(); only an R longjmp from an unprotected R API call in
// a body (e.g. a failed small allocation) could skip it, which would at worst
// leave BLAS at one thread.
static int g_call_depth = 0;

extern "C" void eigencore_call_enter(void) {
  ++g_call_depth;
  eigencore_refresh_thread_count();
}

extern "C" void eigencore_call_leave(void) {
  if (g_call_depth > 0) {
    --g_call_depth;
  }
#ifdef EIGENCORE_HAVE_BLAS_QUIESCE
  if (g_call_depth == 0 && g_blas_quiesced) {
    g_blas_quiesced = false;
    g_blas_set_threads(g_blas_saved_threads);
  }
#endif
}

// c(openmp, processors, default, current, blas_kind, blas_threads):
// blas_kind 0 = no known threading control, 1 = OpenBLAS pthreads/FlexiBLAS
// (switched to one thread during multithreaded sparse solves), 2 = threaded
// BLAS left alone; blas_threads is the controllable BLAS's current thread
// count, or NA.
static SEXP eigencore_thread_info_pack() {
  SEXP out = PROTECT(allocVector(INTSXP, 6));
#ifdef _OPENMP
  INTEGER(out)[0] = 1;
#else
  INTEGER(out)[0] = 0;
#endif
  INTEGER(out)[1] = eigencore_processor_count();
  INTEGER(out)[2] = g_eigencore_default_threads;
  INTEGER(out)[3] = g_eigencore_threads;
#ifdef EIGENCORE_HAVE_BLAS_QUIESCE
  eigencore_blas_probe();
  INTEGER(out)[4] = g_blas_kind;
  INTEGER(out)[5] = g_blas_get_threads != nullptr ? g_blas_get_threads() :
    NA_INTEGER;
#else
  INTEGER(out)[4] = EIGENCORE_BLAS_SERIAL;
  INTEGER(out)[5] = NA_INTEGER;
#endif
  SEXP names = PROTECT(allocVector(STRSXP, 6));
  SET_STRING_ELT(names, 0, mkChar("openmp"));
  SET_STRING_ELT(names, 1, mkChar("processors"));
  SET_STRING_ELT(names, 2, mkChar("default"));
  SET_STRING_ELT(names, 3, mkChar("current"));
  SET_STRING_ELT(names, 4, mkChar("blas_kind"));
  SET_STRING_ELT(names, 5, mkChar("blas_threads"));
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

// eigencore_set_default_threads(n): n < 1 or NA leaves the default unchanged.
// Returns eigencore_thread_info_pack() after the update.
extern "C" SEXP eigencore_set_default_threads(SEXP n_) {
  EIGENCORE_ENTRY_BEGIN
  if ((isReal(n_) || isInteger(n_)) && XLENGTH(n_) >= 1) {
    const double value = isReal(n_) ? REAL(n_)[0] :
      (INTEGER(n_)[0] == NA_INTEGER ? NA_REAL :
         static_cast<double>(INTEGER(n_)[0]));
#ifdef _OPENMP
    g_eigencore_default_threads =
      eigencore_sanitize_thread_count(value, g_eigencore_default_threads);
#else
    (void) value;
    g_eigencore_default_threads = 1;
#endif
  }
  eigencore_refresh_thread_count();
  return eigencore_thread_info_pack();
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_thread_info(void) {
  EIGENCORE_ENTRY_BEGIN
  return eigencore_thread_info_pack();
  EIGENCORE_ENTRY_END
}

// ---------------------------------------------------------------------------
// CSC kernels (P8)
// ---------------------------------------------------------------------------
//
// Adjoint (A^T X): a gather per output column, parallel over the columns of A.
//
// Forward (A X): serially a scatter over the columns of A. In parallel it is a
// gather over the rows of a CSR copy of A (row pointers, column indices and
// values: 4 (nrow + 1) + 12 nnz bytes), built once per operator and cached on
// it. Within each row the CSR copy keeps the nonzeros in column order, and
// the gather accumulates into the output entry with the same expression and
// zero-skip rule as the serial scatter, so A X is bitwise identical for every
// thread count (no atomics, no per-thread partial sums). The copy is built on
// the second forward apply of an operator (or the first one with at least
// kCscCsrMinBlock right-hand sides), so one-off products such as the per-call
// R-level applies never pay for it, and only when it fits in
// getOption("eigencore.csr_cache_mb", 4096) MB. Without it, blocks of more
// than kCscPanelCols columns run in parallel over aligned column chunks and
// narrower blocks run serially.
//
// Multi-RHS blocks of three or more columns go through row-major panels
// (chunks of at most kCscPanelCols columns), so each nonzero touches one or
// two cache lines instead of one per right-hand side.

struct CscApplyCache {
  // Set by the per-call R-level entry points (one apply per .Call, typically
  // an R callback operator inside a solver whose BLAS work dominates): their
  // kernels run serially rather than switch a spinning multithreaded BLAS
  // off, which would slow the surrounding solver's BLAS work.
  bool per_call = false;
  int forward_calls = 0;
  bool csr_ready = false;
  bool csr_unavailable = false;
  std::vector<int> csr_ptr;            // nrow + 1
  std::vector<int> csr_col;            // nnz
  std::vector<double> csr_val;         // nnz
  std::vector<double> panel;           // max(nrow, ncol) * kCscPanelCols
  std::vector<unsigned char> skip;     // ncol zero-column flags
  // Row slabs: slab t holds the nonzeros of rows [slab_rows[t],
  // slab_rows[t + 1]) as its own CSC block; slab_ptr has (n + 1) global
  // offsets per slab.
  int slab_parts = 0;
  bool slab_unavailable = false;
  std::vector<int> slab_rows;
  std::vector<int> slab_ptr;
  std::vector<int> slab_idx;
  std::vector<double> slab_val;
};

static const int kCscPanelCols = 10;
static const int kCscCsrMinBlock = 4;
// Below this many nonzero-times-RHS updates a parallel region costs more than
// it saves.
static const int64_t kCscParallelMinWork = 32768;

static CscApplyCache* csc_apply_cache(CSCOperator* csc) {
  if (!csc->cache) {
    csc->cache = std::make_shared<CscApplyCache>();
  }
  return csc->cache.get();
}

static double* csc_panel(CSCOperator* csc, int64_t rows) {
  CscApplyCache* cache = csc_apply_cache(csc);
  const size_t need = static_cast<size_t>(rows > 0 ? rows : 1) * kCscPanelCols;
  if (cache->panel.size() < need) {
    cache->panel.resize(need);
  }
  return cache->panel.data();
}

static int csc_apply_threads(const CSCOperator* csc, int64_t block_cols) {
  const int threads = eigencore_thread_count();
  if (threads <= 1) {
    return 1;
  }
  const int64_t nnz = csc->col_ptr[csc->cols];
  if (nnz * block_cols < kCscParallelMinWork) {
    return 1;
  }
  if (csc->cache && csc->cache->per_call && eigencore_blas_busy()) {
    return 1;
  }
  return threads;
}

static void csc_mark_per_call(CSCOperator* csc) {
  csc_apply_cache(csc)->per_call = true;
}

static bool csc_build_csr(CSCOperator* csc, CscApplyCache* cache) {
  if (cache->csr_ready) {
    return true;
  }
  if (cache->csr_unavailable) {
    return false;
  }
  const int m = static_cast<int>(csc->rows);
  const int n = static_cast<int>(csc->cols);
  const int* p = csc->col_ptr;
  const int* ri = csc->row_idx;
  const int nnz = p[n];
  const double bytes = 4.0 * (static_cast<double>(m) + 1.0) +
    12.0 * static_cast<double>(nnz);
  if (bytes > g_eigencore_csr_cache_bytes) {
    cache->csr_unavailable = true;
    return false;
  }
  try {
    std::vector<int> ptr(static_cast<size_t>(m) + 1, 0);
    for (int pos = 0; pos < nnz; ++pos) {
      ++ptr[static_cast<size_t>(ri[pos]) + 1];
    }
    for (int row = 0; row < m; ++row) {
      ptr[static_cast<size_t>(row) + 1] += ptr[static_cast<size_t>(row)];
    }
    std::vector<int> next(ptr.begin(), ptr.end() - 1);
    std::vector<int> cols(eigencore_buffer_size(nnz));
    std::vector<double> vals(eigencore_buffer_size(nnz));
    // Column-major traversal keeps each row's nonzeros in column order, the
    // order in which the serial scatter accumulates them.
    for (int col = 0; col < n; ++col) {
      for (int pos = p[col]; pos < p[col + 1]; ++pos) {
        const int dst = next[static_cast<size_t>(ri[pos])]++;
        cols[static_cast<size_t>(dst)] = col;
        vals[static_cast<size_t>(dst)] = csc->values[pos];
      }
    }
    cache->csr_ptr.swap(ptr);
    cache->csr_col.swap(cols);
    cache->csr_val.swap(vals);
  } catch (const std::bad_alloc&) {
    cache->csr_unavailable = true;
    return false;
  }
  cache->csr_ready = true;
  return true;
}

static bool csc_build_slabs(CSCOperator* csc, CscApplyCache* cache, int parts) {
  if (cache->slab_parts == parts) {
    return true;
  }
  if (cache->slab_unavailable) {
    return false;
  }
  const int m = static_cast<int>(csc->rows);
  const int n = static_cast<int>(csc->cols);
  const int* p = csc->col_ptr;
  const int* ri = csc->row_idx;
  const int nnz = p[n];
  const double bytes = 4.0 * static_cast<double>(parts) * (n + 1.0) +
    12.0 * static_cast<double>(nnz);
  if (bytes > g_eigencore_csr_cache_bytes) {
    cache->slab_unavailable = true;
    return false;
  }
  try {
    std::vector<int> counts(eigencore_buffer_size(m), 0);
    for (int pos = 0; pos < nnz; ++pos) {
      ++counts[static_cast<size_t>(ri[pos])];
    }
    std::vector<int> rows(static_cast<size_t>(parts) + 1, m);
    rows[0] = 0;
    int64_t running = 0;
    int next_part = 1;
    for (int row = 0; row < m && next_part < parts; ++row) {
      running += counts[static_cast<size_t>(row)];
      while (next_part < parts &&
             running * parts >= static_cast<int64_t>(nnz) * next_part) {
        rows[static_cast<size_t>(next_part)] = row + 1;
        ++next_part;
      }
    }
    std::vector<int> slab_of_row(eigencore_buffer_size(m), 0);
    for (int t = 0; t < parts; ++t) {
      for (int row = rows[t]; row < rows[t + 1]; ++row) {
        slab_of_row[static_cast<size_t>(row)] = t;
      }
    }
    const size_t stride = static_cast<size_t>(n) + 1;
    std::vector<int> ptr(static_cast<size_t>(parts) * stride, 0);
    for (int col = 0; col < n; ++col) {
      for (int pos = p[col]; pos < p[col + 1]; ++pos) {
        ++ptr[static_cast<size_t>(slab_of_row[ri[pos]]) * stride + col + 1];
      }
    }
    int offset = 0;
    for (int t = 0; t < parts; ++t) {
      int* slab = ptr.data() + static_cast<size_t>(t) * stride;
      slab[0] = offset;
      for (int col = 0; col < n; ++col) {
        slab[col + 1] += slab[col];
      }
      offset = slab[n];
    }
    std::vector<int> next(ptr);
    std::vector<int> idx(eigencore_buffer_size(nnz));
    std::vector<double> val(eigencore_buffer_size(nnz));
    // Position order within each column is kept, so every output row still
    // accumulates its nonzeros in the serial scatter's order.
    for (int col = 0; col < n; ++col) {
      for (int pos = p[col]; pos < p[col + 1]; ++pos) {
        const int row = ri[pos];
        const int dst = next[static_cast<size_t>(slab_of_row[row]) * stride + col]++;
        idx[static_cast<size_t>(dst)] = row;
        val[static_cast<size_t>(dst)] = csc->values[pos];
      }
    }
    cache->slab_rows.swap(rows);
    cache->slab_ptr.swap(ptr);
    cache->slab_idx.swap(idx);
    cache->slab_val.swap(val);
  } catch (const std::bad_alloc&) {
    cache->slab_unavailable = true;
    return false;
  }
  cache->slab_parts = parts;
  return true;
}

// Serial forward scatter of one chunk of c <= kCscPanelCols right-hand sides.
// kScaled multiplies x by the column weights first and uses the
// coefficient-first update of the centered-scaled operator; otherwise the
// update is (alpha a_ij) x_j. Both match the historical serial expressions bit
// for bit; the parallel CSR gather below reproduces them.
template <bool kScaled>
static void csc_forward_chunk(int n, const int* cp, const int* ri,
                              const double* values, int r0, int r1, int c,
                              const double* X, int64_t ldx,
                              const double* weights, double alpha,
                              double* Y, int64_t ldy, double* panel) {
  if (c == 1) {
    for (int col = 0; col < n; ++col) {
      const double xv = kScaled ? weights[col] * X[col] : X[col];
      if (xv == 0.0) continue;
      if (kScaled) {
        const double coefficient = alpha * xv;
        for (int pos = cp[col]; pos < cp[col + 1]; ++pos) {
          Y[ri[pos]] += coefficient * values[pos];
        }
      } else {
        for (int pos = cp[col]; pos < cp[col + 1]; ++pos) {
          Y[ri[pos]] += alpha * values[pos] * xv;
        }
      }
    }
    return;
  }
  const double* xptr[kCscPanelCols];
  double* yptr[kCscPanelCols];
  for (int block = 0; block < c; ++block) {
    xptr[block] = X + block * ldx;
    yptr[block] = Y + block * ldy;
  }
  double xval[kCscPanelCols];
  if (panel != nullptr) {
    for (int row = r0; row < r1; ++row) {
      double* prow = panel + static_cast<int64_t>(row) * kCscPanelCols;
      for (int block = 0; block < c; ++block) {
        prow[block] = yptr[block][row];
      }
    }
  }
  for (int col = 0; col < n; ++col) {
    bool all_zero = true;
    for (int block = 0; block < c; ++block) {
      xval[block] = kScaled ? weights[col] * xptr[block][col] : xptr[block][col];
      all_zero = all_zero && xval[block] == 0.0;
    }
    if (all_zero) continue;
    if (kScaled) {
      for (int block = 0; block < c; ++block) {
        xval[block] *= alpha;
      }
    }
    for (int pos = cp[col]; pos < cp[col + 1]; ++pos) {
      const int row = ri[pos];
      const double a = kScaled ? values[pos] : alpha * values[pos];
      if (panel != nullptr) {
        double* prow = panel + static_cast<int64_t>(row) * kCscPanelCols;
        for (int block = 0; block < c; ++block) {
          if (kScaled) {
            prow[block] += xval[block] * a;
          } else {
            prow[block] += a * xval[block];
          }
        }
      } else {
        for (int block = 0; block < c; ++block) {
          if (kScaled) {
            yptr[block][row] += xval[block] * a;
          } else {
            yptr[block][row] += a * xval[block];
          }
        }
      }
    }
  }
  if (panel != nullptr) {
    for (int block = 0; block < c; ++block) {
      double* y = yptr[block];
      for (int row = r0; row < r1; ++row) {
        y[row] = panel[static_cast<int64_t>(row) * kCscPanelCols + block];
      }
    }
  }
}

// Rows [row_begin, row_end) of the parallel forward gather for one chunk.
// xp is the chunk's row-major column panel (ld kCscPanelCols) holding x (plain)
// or alpha * w .* x (kScaled); skip flags the columns the serial scatter
// skips. With c == 1 and !kScaled, x is read directly from X.
template <bool kScaled>
static void csc_forward_gather_rows(const CscApplyCache* cache, int row_begin,
                                    int row_end, int c, const double* X,
                                    const double* xp,
                                    const unsigned char* skip, double alpha,
                                    double* Y, int64_t ldy) {
  const int* ptr = cache->csr_ptr.data();
  const int* cols = cache->csr_col.data();
  const double* vals = cache->csr_val.data();
  if (c == 1) {
    for (int row = row_begin; row < row_end; ++row) {
      double acc = Y[row];
      for (int k = ptr[row]; k < ptr[row + 1]; ++k) {
        const int col = cols[k];
        if (kScaled) {
          if (skip[col]) continue;
          acc += xp[static_cast<int64_t>(col) * kCscPanelCols] * vals[k];
        } else {
          const double xv = X[col];
          if (xv == 0.0) continue;
          acc += alpha * vals[k] * xv;
        }
      }
      Y[row] = acc;
    }
    return;
  }
  double acc[kCscPanelCols];
  for (int row = row_begin; row < row_end; ++row) {
    for (int block = 0; block < c; ++block) {
      acc[block] = Y[row + block * ldy];
    }
    for (int k = ptr[row]; k < ptr[row + 1]; ++k) {
      const int col = cols[k];
      if (skip[col]) continue;
      const double* xrow = xp + static_cast<int64_t>(col) * kCscPanelCols;
      if (kScaled) {
        const double a = vals[k];
        for (int block = 0; block < c; ++block) {
          acc[block] += xrow[block] * a;
        }
      } else {
        const double a = alpha * vals[k];
        for (int block = 0; block < c; ++block) {
          acc[block] += a * xrow[block];
        }
      }
    }
    for (int block = 0; block < c; ++block) {
      Y[row + block * ldy] = acc[block];
    }
  }
}

template <bool kScaled>
static void csc_forward_gather(CSCOperator* csc, int threads,
                               int64_t block_cols, const double* X,
                               int64_t ldx, const double* weights,
                               double alpha, double* Y, int64_t ldy) {
  CscApplyCache* cache = csc_apply_cache(csc);
  const int m = static_cast<int>(csc->rows);
  const int n = static_cast<int>(csc->cols);
  const bool need_panel = kScaled || block_cols > 1;
  double* xp = need_panel ? csc_panel(csc, n) : nullptr;
  if (need_panel && cache->skip.size() < static_cast<size_t>(n)) {
    cache->skip.resize(static_cast<size_t>(n));
  }
  unsigned char* skip = cache->skip.data();
  const int* ptr = cache->csr_ptr.data();
  const int64_t nnz = ptr[m];
  eigencore_blas_quiesce();
  for (int64_t chunk = 0; chunk < block_cols; chunk += kCscPanelCols) {
    const int c = static_cast<int>(
      std::min<int64_t>(kCscPanelCols, block_cols - chunk));
    const double* Xc = X + chunk * ldx;
    double* Yc = Y + chunk * ldy;
    EIGENCORE_OMP(omp parallel num_threads(threads))
    {
      if (need_panel) {
        EIGENCORE_OMP(omp for schedule(static))
        for (int col = 0; col < n; ++col) {
          double* prow = xp + static_cast<int64_t>(col) * kCscPanelCols;
          bool all_zero = true;
          for (int block = 0; block < c; ++block) {
            const double x = Xc[col + block * ldx];
            prow[block] = kScaled ? weights[col] * x : x;
            all_zero = all_zero && prow[block] == 0.0;
          }
          if (kScaled) {
            for (int block = 0; block < c; ++block) {
              prow[block] *= alpha;
            }
          }
          skip[col] = all_zero ? 1 : 0;
        }
      }
      // Contiguous row ranges balanced by nonzero count.
      const int nt = eigencore_omp_num_threads();
      const int t = eigencore_omp_thread_num();
      const int64_t lo_target = nnz * t / nt;
      const int64_t hi_target = nnz * (t + 1) / nt;
      const int row_begin = (t == 0) ? 0 : static_cast<int>(
        std::lower_bound(ptr, ptr + m + 1, lo_target) - ptr);
      const int row_end = (t + 1 == nt) ? m : static_cast<int>(
        std::lower_bound(ptr, ptr + m + 1, hi_target) - ptr);
      csc_forward_gather_rows<kScaled>(cache, std::min(row_begin, m),
                                       std::min(row_end, m), c, Xc, xp, skip,
                                       alpha, Yc, ldy);
    }
  }
}

template <bool kScaled>
static void csc_forward_apply(CSCOperator* csc, int64_t block_cols,
                              const double* X, int64_t ldx,
                              const double* weights, double alpha,
                              double* Y, int64_t ldy) {
  const int m = static_cast<int>(csc->rows);
  const int n = static_cast<int>(csc->cols);
  if (block_cols <= 0 || m == 0 || n == 0) {
    return;
  }
  const int threads = csc_apply_threads(csc, block_cols);
  const int64_t chunks = (block_cols + kCscPanelCols - 1) / kCscPanelCols;
  if (threads > 1) {
    CscApplyCache* cache = csc_apply_cache(csc);
    ++cache->forward_calls;
    const bool amortized =
      cache->forward_calls >= 2 || block_cols >= kCscCsrMinBlock;
    // Long rows (n * threads >= m, i.e. at least as many nonzeros per row as
    // per column per thread) favour the CSR gather; tall matrices with short
    // rows favour per-thread row slabs, whose scatter runs over longer column
    // segments. The first choice sticks for the operator's lifetime.
    const bool use_slabs = !cache->csr_ready &&
      (cache->slab_parts > 0 ||
       static_cast<int64_t>(n) * threads < static_cast<int64_t>(m));
    if (amortized && use_slabs && csc_build_slabs(csc, cache, threads)) {
      double* panel = (block_cols >= 3 && csc->col_ptr[n] >= m) ?
        csc_panel(csc, m) : nullptr;
      const int* sp = cache->slab_ptr.data();
      const int* sidx = cache->slab_idx.data();
      const double* sval = cache->slab_val.data();
      const int* srows = cache->slab_rows.data();
      const int parts = cache->slab_parts;
      eigencore_blas_quiesce();
      EIGENCORE_OMP(omp parallel num_threads(parts))
      {
        const int stride = eigencore_omp_num_threads();
        for (int t = eigencore_omp_thread_num(); t < parts; t += stride) {
          const int* cp = sp + static_cast<int64_t>(t) * (n + 1);
          for (int64_t chunk = 0; chunk < block_cols; chunk += kCscPanelCols) {
            const int c = static_cast<int>(
              std::min<int64_t>(kCscPanelCols, block_cols - chunk));
            csc_forward_chunk<kScaled>(n, cp, sidx, sval, srows[t],
                                       srows[t + 1], c, X + chunk * ldx, ldx,
                                       weights, alpha, Y + chunk * ldy, ldy,
                                       panel);
          }
        }
      }
      return;
    }
    if (amortized && csc_build_csr(csc, cache)) {
      csc_forward_gather<kScaled>(csc, threads, block_cols, X, ldx, weights,
                                  alpha, Y, ldy);
      return;
    }
    if (chunks > 1) {
      // No CSR copy: split aligned column chunks across threads.
      const int use = static_cast<int>(std::min<int64_t>(threads, chunks));
      (void) use;  // only read by the OpenMP pragma
      eigencore_blas_quiesce();
      EIGENCORE_OMP(omp parallel for num_threads(use) schedule(static))
      for (int64_t chunk_id = 0; chunk_id < chunks; ++chunk_id) {
        const int64_t chunk = chunk_id * kCscPanelCols;
        const int c = static_cast<int>(
          std::min<int64_t>(kCscPanelCols, block_cols - chunk));
        csc_forward_chunk<kScaled>(n, csc->col_ptr, csc->row_idx,
                                   csc->values, 0, m, c, X + chunk * ldx, ldx,
                                   weights, alpha, Y + chunk * ldy, ldy,
                                   nullptr);
      }
      return;
    }
  }
  double* panel = nullptr;
  if (block_cols >= 3 && csc->col_ptr[n] >= m) {
    panel = csc_panel(csc, m);
  }
  for (int64_t chunk = 0; chunk < block_cols; chunk += kCscPanelCols) {
    const int c = static_cast<int>(
      std::min<int64_t>(kCscPanelCols, block_cols - chunk));
    csc_forward_chunk<kScaled>(n, csc->col_ptr, csc->row_idx, csc->values, 0,
                               m, c, X + chunk * ldx, ldx, weights, alpha,
                               Y + chunk * ldy, ldy, panel);
  }
}

// Adjoint gather for one column and one chunk of c right-hand sides. xt, when non-null, is the row-major panel of the chunk's X columns
// (ld kCscPanelCols). kScaled applies the centered-scaled epilogue
// alpha w_j (dot - mu_j sum(x)); otherwise alpha * dot.
template <bool kScaled>
static inline void csc_adjoint_column(const CSCOperator* csc, int col, int c,
                                      const double* const* xptr,
                                      const double* xt,
                                      const double* weights,
                                      const double* means,
                                      const double* xsum, double alpha,
                                      double* const* yptr) {
  const int* ri = csc->row_idx;
  const double* values = csc->values;
  const int begin = csc->col_ptr[col];
  const int end = csc->col_ptr[col + 1];
  if (c == 1) {
    double acc = 0.0;
    const double* x = xptr[0];
    for (int pos = begin; pos < end; ++pos) {
      acc += values[pos] * x[ri[pos]];
    }
    if (kScaled) {
      yptr[0][col] += alpha * weights[col] * (acc - means[col] * xsum[0]);
    } else {
      yptr[0][col] += alpha * acc;
    }
    return;
  }
  double acc[kCscPanelCols];
  for (int block = 0; block < c; ++block) {
    acc[block] = 0.0;
  }
  if (xt != nullptr) {
    for (int pos = begin; pos < end; ++pos) {
      const double a = values[pos];
      const double* xrow = xt + static_cast<int64_t>(ri[pos]) * kCscPanelCols;
      for (int block = 0; block < c; ++block) {
        acc[block] += a * xrow[block];
      }
    }
  } else {
    for (int pos = begin; pos < end; ++pos) {
      const int row = ri[pos];
      const double a = values[pos];
      for (int block = 0; block < c; ++block) {
        acc[block] += a * xptr[block][row];
      }
    }
  }
  for (int block = 0; block < c; ++block) {
    if (kScaled) {
      yptr[block][col] +=
        alpha * weights[col] * (acc[block] - means[col] * xsum[block]);
    } else {
      yptr[block][col] += alpha * acc[block];
    }
  }
}

template <bool kScaled>
static void csc_adjoint_apply(CSCOperator* csc, int64_t block_cols,
                              const double* X, int64_t ldx,
                              const double* weights, const double* means,
                              double alpha, double* Y, int64_t ldy) {
  const int m = static_cast<int>(csc->rows);
  const int n = static_cast<int>(csc->cols);
  if (block_cols <= 0 || n == 0) {
    return;
  }
  const int threads = csc_apply_threads(csc, block_cols);
  for (int64_t chunk = 0; chunk < block_cols; chunk += kCscPanelCols) {
    const int c = static_cast<int>(
      std::min<int64_t>(kCscPanelCols, block_cols - chunk));
    const double* xptr[kCscPanelCols];
    double* yptr[kCscPanelCols];
    double xsum[kCscPanelCols];
    for (int block = 0; block < c; ++block) {
      xptr[block] = X + (chunk + block) * ldx;
      yptr[block] = Y + (chunk + block) * ldy;
      xsum[block] = 0.0;
      if (kScaled) {
        for (int row = 0; row < m; ++row) {
          xsum[block] += xptr[block][row];
        }
      }
    }
    if (threads > 1) {
      eigencore_blas_quiesce();
    }
    const double* xt = nullptr;
    if (c >= 3 && m > 0 && csc->col_ptr[n] >= m) {
      double* panel = csc_panel(csc, m);
      EIGENCORE_OMP(omp parallel for num_threads(threads) schedule(static) if(threads > 1))
      for (int row = 0; row < m; ++row) {
        double* prow = panel + static_cast<int64_t>(row) * kCscPanelCols;
        for (int block = 0; block < c; ++block) {
          prow[block] = xptr[block][row];
        }
      }
      xt = panel;
    }
    if (threads > 1) {
      EIGENCORE_OMP(omp parallel for num_threads(threads) schedule(dynamic, 256))
      for (int col = 0; col < n; ++col) {
        csc_adjoint_column<kScaled>(csc, col, c, xptr, xt, weights, means,
                                    xsum, alpha, yptr);
      }
    } else {
      for (int col = 0; col < n; ++col) {
        csc_adjoint_column<kScaled>(csc, col, c, xptr, xt, weights, means,
                                    xsum, alpha, yptr);
      }
    }
  }
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
    csc_forward_apply<false>(csc, block_cols, X, ldx, nullptr, alpha, Y, ldy);
  } else {
    csc_adjoint_apply<false>(csc, block_cols, X, ldx, nullptr, nullptr,
                             alpha, Y, ldy);
  }
  return 0;
}

extern "C" int eigencore_centered_scaled_csc_apply(
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
  CenteredScaledCSCOperator* fused =
    static_cast<CenteredScaledCSCOperator*>(impl);
  CSCOperator* csc = &fused->base;
  if (op != EIGENCORE_TRANSPOSE_NONE &&
      op != EIGENCORE_TRANSPOSE_ADJOINT) {
    return -1;
  }
  const int64_t out_rows =
    (op == EIGENCORE_TRANSPOSE_ADJOINT) ? csc->cols : csc->rows;
  const int64_t inner =
    (op == EIGENCORE_TRANSPOSE_ADJOINT) ? csc->rows : csc->cols;
  if (block_cols < 0 || ldx < inner || ldy < out_rows ||
      fused->col_means == nullptr || fused->col_weights == nullptr) {
    return -1;
  }

  scale_or_zero_output(Y, out_rows, block_cols, beta);
  if (alpha == 0.0) {
    return 0;
  }

  if (op == EIGENCORE_TRANSPOSE_NONE) {
    // Y <- alpha (A - 1 mu^T) D X + beta Y: scatter alpha A D X, then subtract
    // the rank-one correction alpha (mu^T D x) from every row.
    csc_forward_apply<true>(csc, block_cols, X, ldx, fused->col_weights,
                            alpha, Y, ldy);
    for (int64_t block = 0; block < block_cols; ++block) {
      const double* x_col = X + block * ldx;
      double* y_col = Y + block * ldy;
      double correction = 0.0;
      for (int64_t col = 0; col < csc->cols; ++col) {
        const double scaled_x = fused->col_weights[col] * x_col[col];
        correction += fused->col_means[col] * scaled_x;
      }
      correction *= alpha;
      for (int64_t row = 0; row < csc->rows; ++row) {
        y_col[row] -= correction;
      }
    }
  } else {
    // Y <- alpha D (A^T X - mu 1^T X) + beta Y.
    csc_adjoint_apply<true>(csc, block_cols, X, ldx, fused->col_weights,
                            fused->col_means, alpha, Y, ldy);
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

// Output buffer of the R-level block applies (P11). When beta == 0 the old
// contents of Y are never read (every kernel zero-fills or uses BLAS beta = 0),
// so a fresh matrix replaces the former duplicate(Y); dimnames are kept. Y may
// be NULL, meaning a zero matrix of the output shape: callers set beta = 0
// then, so R wrappers need not allocate Y at all.
static SEXP block_apply_output(SEXP Y_, SEXPTYPE type, int rows, int cols,
                               bool beta_is_zero) {
  if (Y_ != R_NilValue && !beta_is_zero) {
    return duplicate(Y_);
  }
  SEXP out_ = PROTECT(allocMatrix(type, rows, cols));
  if (Y_ != R_NilValue) {
    SEXP dimnames_ = getAttrib(Y_, R_DimNamesSymbol);
    if (dimnames_ != R_NilValue) {
      setAttrib(out_, R_DimNamesSymbol, dimnames_);
    }
  }
  UNPROTECT(1);
  return out_;
}

extern "C" SEXP eigencore_dense_block_apply(SEXP A_, SEXP X_, SEXP alpha_,
                                            SEXP beta_, SEXP Y_,
                                            SEXP transpose_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(A_) || !isReal(X_) || !(isReal(Y_) || isNull(Y_))) {
    error("A, X, and Y must be double matrices");
  }

  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimA == R_NilValue || dimX == R_NilValue || (dimY == R_NilValue && !isNull(Y_))) {
    error("A, X, and Y must be matrices");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = isNull(Y_) ? -1 : INTEGER(dimY)[0];
  const int yc = isNull(Y_) ? -1 : INTEGER(dimY)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const double alpha = REAL(alpha_)[0];
  const double beta = isNull(Y_) ? 0.0 : REAL(beta_)[0];

  const int inner = transpose ? m : n;
  const int out_rows = transpose ? n : m;
  if (xr != inner) {
    error("non-conformable X for dense block apply");
  }
  if ((!isNull(Y_) && (yr != out_rows || yc != xc))) {
    error("non-conformable Y for dense block apply");
  }

  SEXP out_ = PROTECT(block_apply_output(Y_, REALSXP, out_rows, xc, beta == 0.0));
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_complex_block_apply(SEXP A_, SEXP X_, SEXP alpha_,
                                                    SEXP beta_, SEXP Y_,
                                                    SEXP adjoint_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isComplex(A_) || !isComplex(X_) || !(isComplex(Y_) || isNull(Y_))) {
    error("A, X, and Y must be complex matrices");
  }
  if (!isLogical(adjoint_) || LENGTH(adjoint_) != 1) {
    error("adjoint must be a logical scalar");
  }

  SEXP dimA = getAttrib(A_, R_DimSymbol);
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimA == R_NilValue || dimX == R_NilValue || (dimY == R_NilValue && !isNull(Y_))) {
    error("A, X, and Y must be matrices");
  }

  const int m = INTEGER(dimA)[0];
  const int n = INTEGER(dimA)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = isNull(Y_) ? -1 : INTEGER(dimY)[0];
  const int yc = isNull(Y_) ? -1 : INTEGER(dimY)[1];
  const bool adjoint = LOGICAL(adjoint_)[0];

  const int inner = adjoint ? m : n;
  const int out_rows = adjoint ? n : m;
  if (xr != inner) {
    error("non-conformable X for dense complex block apply");
  }
  if ((!isNull(Y_) && (yr != out_rows || yc != xc))) {
    error("non-conformable Y for dense complex block apply");
  }

  Rcomplex beta = scalar_as_rcomplex(beta_, "beta");
  if (isNull(Y_)) {
    beta.r = 0.0;
    beta.i = 0.0;
  }
  SEXP out_ = PROTECT(block_apply_output(Y_, CPLXSXP, out_rows, xc,
                                         beta.r == 0.0 && beta.i == 0.0));
  DenseComplexColumnMajorOperator impl = {m, n, COMPLEX(A_)};
  const int status = eigencore_dense_complex_apply(
    &impl,
    adjoint ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    COMPLEX(X_),
    xr,
    scalar_as_rcomplex(alpha_, "alpha"),
    beta,
    COMPLEX(out_),
    out_rows,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error("dense complex block apply", status);
  }

  UNPROTECT(1);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_randomized_apply(SEXP A_, SEXP X_,
                                                 SEXP transpose_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_randomized_sketch(SEXP A_, SEXP cols_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_randomized_project_transposed(SEXP A_,
                                                              SEXP Q_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
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
  double applied_bound = 0.0;
  bool passed = false;
};

// Largest column 2-norm: a structural lower bound on ||A||_2 (C12).
static double dense_randomized_max_column_norm(const double* A, int m, int n) {
  double best = 0.0;
  for (int col = 0; col < n; ++col) {
    const double* a = A + static_cast<int64_t>(col) * m;
    double sum = 0.0;
    for (int row = 0; row < m; ++row) {
      sum += a[row] * a[row];
    }
    if (sum > best) {
      best = sum;
    }
  }
  return std::sqrt(best);
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

// Residuals, spectral-norm lower bound, and backward errors for a randomized
// SVD certificate. On entry AV = A V and ATU = A^T U, computed from A in this
// call. The denominator is L = max(column bound, ||A v_j|| / ||v_j||,
// ||A^T u_j|| / ||u_j||) <= ||A||_2, so the backward error over-estimates the
// normwise one and `passed` stays sound (C12).
static void randomized_certificate_finalize(DenseRandomizedCertificate& cert,
                                            std::vector<double>& AV,
                                            std::vector<double>& ATU,
                                            const std::vector<double>& d,
                                            const std::vector<double>& U,
                                            const std::vector<double>& V,
                                            int m, int n, int rank,
                                            double column_bound,
                                            double tol) {
  double applied = 0.0;
  for (int col = 0; col < rank; ++col) {
    double* av = AV.data() + static_cast<int64_t>(col) * m;
    double* atu = ATU.data() + static_cast<int64_t>(col) * n;
    const double* u = U.data() + static_cast<int64_t>(col) * m;
    const double* v = V.data() + static_cast<int64_t>(col) * n;
    double av_sq = 0.0;
    double u_sq = 0.0;
    for (int row = 0; row < m; ++row) {
      av_sq += av[row] * av[row];
      u_sq += u[row] * u[row];
    }
    double atu_sq = 0.0;
    double v_sq = 0.0;
    for (int row = 0; row < n; ++row) {
      atu_sq += atu[row] * atu[row];
      v_sq += v[row] * v[row];
    }
    if (v_sq > 0.0) {
      const double ratio = std::sqrt(av_sq / v_sq);
      if (std::isfinite(ratio) && ratio > applied) applied = ratio;
    }
    if (u_sq > 0.0) {
      const double ratio = std::sqrt(atu_sq / u_sq);
      if (std::isfinite(ratio) && ratio > applied) applied = ratio;
    }
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
  }
  cert.applied_bound = applied;
  cert.scale = std::max(column_bound, applied);
  if (!(cert.scale >= 2.2204460492503131e-16)) {
    cert.scale = 2.2204460492503131e-16;
  }
  bool all_converged = true;
  for (int col = 0; col < rank; ++col) {
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
  const double column_bound = dense_randomized_max_column_norm(A, m, n);

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

  randomized_certificate_finalize(cert, AV, ATU, d, U, V, m, n, rank,
                                  column_bound, tol);
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

// Largest column 2-norm of a CSC matrix: a structural lower bound on ||A||_2.
static double csc_randomized_max_column_norm(const CSCOperator& impl) {
  double best = 0.0;
  for (int col = 0; col < impl.cols; ++col) {
    double sum = 0.0;
    for (int pos = impl.col_ptr[col]; pos < impl.col_ptr[col + 1]; ++pos) {
      sum += impl.values[pos] * impl.values[pos];
    }
    if (sum > best) {
      best = sum;
    }
  }
  return std::sqrt(best);
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

// B (q_cols x n, column-major) += (A^T Q)^T column by column: a gather per
// column of A over a row-major copy Qt of Q, parallel over columns of A
// (each output column is owned by one thread, so results do not depend on
// the thread count). Qt must hold m * q_cols doubles.
static void csc_project_transposed_kernel(const int* row_idx, const int* col_ptr,
                                          const double* values, int m, int n,
                                          const double* Q, int q_cols,
                                          double* Qt, double* B,
                                          bool per_call) {
  int threads = eigencore_thread_count();
  if (static_cast<int64_t>(col_ptr[n]) * q_cols < kCscParallelMinWork ||
      (per_call && threads > 1 && eigencore_blas_busy())) {
    threads = 1;
  }
  if (threads > 1) {
    eigencore_blas_quiesce();
  }
  EIGENCORE_OMP(omp parallel for num_threads(threads) schedule(static) if(threads > 1))
  for (int row = 0; row < m; ++row) {
    for (int block = 0; block < q_cols; ++block) {
      Qt[static_cast<int64_t>(row) * q_cols + block] =
        Q[static_cast<int64_t>(block) * m + row];
    }
  }
  EIGENCORE_OMP(omp parallel for num_threads(threads) schedule(dynamic, 256) if(threads > 1))
  for (int col = 0; col < n; ++col) {
    double* out_col = B + static_cast<int64_t>(col) * q_cols;
    for (int pos = col_ptr[col]; pos < col_ptr[col + 1]; ++pos) {
      const int row = row_idx[pos];
      const double a = values[pos];
      const double* qt_row = Qt + static_cast<int64_t>(row) * q_cols;
      for (int block = 0; block < q_cols; ++block) {
        out_col[block] += a * qt_row[block];
      }
    }
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
  csc_project_transposed_kernel(impl.row_idx, impl.col_ptr, impl.values, m, n,
                                Q.data(), q_cols, Qt.data(), B.data(), false);
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
  (void)nnz;
  const double column_bound = csc_randomized_max_column_norm(impl);

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

  randomized_certificate_finalize(cert, AV, ATU, d, U, V, m, n, rank,
                                  column_bound, tol);
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

  SEXP out_ = PROTECT(allocVector(VECSXP, 8));
  SEXP scale_ = PROTECT(ScalarReal(cert.scale));
  SEXP passed_ = PROTECT(ScalarLogical(cert.passed));
  SEXP norm_ = PROTECT(ScalarReal(cert.scale));
  SEXP applied_ = PROTECT(ScalarReal(cert.applied_bound));
  SET_VECTOR_ELT(out_, 0, residuals_);
  SET_VECTOR_ELT(out_, 1, backward_);
  SET_VECTOR_ELT(out_, 2, orth_);
  SET_VECTOR_ELT(out_, 3, converged_);
  SET_VECTOR_ELT(out_, 4, scale_);
  SET_VECTOR_ELT(out_, 5, passed_);
  SET_VECTOR_ELT(out_, 6, norm_);
  SET_VECTOR_ELT(out_, 7, applied_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(names_, 0, mkChar("residuals"));
  SET_STRING_ELT(names_, 1, mkChar("backward_error"));
  SET_STRING_ELT(names_, 2, mkChar("orthogonality"));
  SET_STRING_ELT(names_, 3, mkChar("converged"));
  SET_STRING_ELT(names_, 4, mkChar("scale"));
  SET_STRING_ELT(names_, 5, mkChar("passed"));
  SET_STRING_ELT(names_, 6, mkChar("norm_A"));
  SET_STRING_ELT(names_, 7, mkChar("norm_A_applied_bound"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(15);
  return out_;
}

extern "C" SEXP eigencore_dense_randomized_svd_controller(
    SEXP A_, SEXP rank_, SEXP oversample_, SEXP n_iter_, SEXP normalizer_,
    SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
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
      eigencore_check_interrupt();
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_randomized_svd_controller(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP rank_, SEXP oversample_,
    SEXP n_iter_, SEXP normalizer_, SEXP tol_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC randomized controller inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_randomized_svd_controller");
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
      eigencore_check_interrupt();
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_block_apply(SEXP i_, SEXP p_, SEXP x_, SEXP dim_,
                                          SEXP X_, SEXP alpha_, SEXP beta_,
                                          SEXP Y_, SEXP transpose_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(X_) || !(isReal(Y_) || isNull(Y_))) {
    error("invalid CSC block apply inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_block_apply");

  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || (dimY == R_NilValue && !isNull(Y_))) {
    error("X and Y must be matrices");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = isNull(Y_) ? -1 : INTEGER(dimY)[0];
  const int yc = isNull(Y_) ? -1 : INTEGER(dimY)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const double alpha = REAL(alpha_)[0];
  const double beta = isNull(Y_) ? 0.0 : REAL(beta_)[0];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;

  if (xr != inner) {
    error("non-conformable X for CSC block apply");
  }
  if ((!isNull(Y_) && (yr != out_rows || yc != xc))) {
    error("non-conformable Y for CSC block apply");
  }

  SEXP out_ = PROTECT(block_apply_output(Y_, REALSXP, out_rows, xc, beta == 0.0));
  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  const double* X = REAL(X_);
  double* out = REAL(out_);

  CSCOperator impl = {m, n, row_idx, col_ptr, values};
  csc_mark_per_call(&impl);
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_randomized_apply(SEXP i_, SEXP p_, SEXP x_,
                                               SEXP dim_, SEXP X_,
                                               SEXP transpose_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(X_) || !isLogical(transpose_)) {
    error("invalid CSC randomized apply inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_randomized_apply");
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
  csc_mark_per_call(&impl);
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_randomized_sketch(SEXP i_, SEXP p_, SEXP x_,
                                                SEXP dim_, SEXP cols_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_)) {
    error("invalid CSC randomized sketch inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_randomized_sketch");

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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_randomized_project_transposed(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP Q_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(Q_)) {
    error("invalid CSC randomized projection inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_randomized_project_transposed");
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
  csc_project_transposed_kernel(row_idx, col_ptr, values, m, n, Q, qcols,
                                Qt.data(), out, true);
  SEXP transposed_ = PROTECT(ScalarLogical(TRUE));
  setAttrib(out_, install("transposed"), transposed_);
  UNPROTECT(2);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_column_moments(SEXP p_, SEXP x_, SEXP dim_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      LENGTH(dim_) != 2) {
    error("invalid CSC column-moment inputs");
  }
  eigencore_validate_csc_structure(R_NilValue, p_, x_, dim_, "csc_column_moments");
  const int n = INTEGER(dim_)[1];
  const int m = INTEGER(dim_)[0];
  if (m < 0 || n < 0 || LENGTH(p_) != n + 1 || INTEGER(p_)[0] != 0 ||
      INTEGER(p_)[n] != XLENGTH(x_)) {
    error("inconsistent CSC column pointers");
  }

  SEXP sums_ = PROTECT(allocVector(REALSXP, n));
  SEXP sum_squares_ = PROTECT(allocVector(REALSXP, n));
  SEXP means_ = PROTECT(allocVector(REALSXP, n));
  SEXP centered_sum_squares_ = PROTECT(allocVector(REALSXP, n));
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  for (int col = 0; col < n; ++col) {
    long double sum = 0.0L;
    long double sum_squares = 0.0L;
    const int begin = col_ptr[col];
    const int end = col_ptr[col + 1];
    const int nonzero = end - begin;
    if (nonzero < 0 || nonzero > m) {
      error("inconsistent CSC column occupancy");
    }
    const long double shift = nonzero ?
      static_cast<long double>(values[begin]) : 0.0L;
    long double shifted_sum =
      -static_cast<long double>(m - nonzero) * shift;
    long double shifted_sum_squares =
      static_cast<long double>(m - nonzero) * shift * shift;
    for (int pos = begin; pos < end; ++pos) {
      const long double value = static_cast<long double>(values[pos]);
      sum += value;
      sum_squares += value * value;
      const long double delta = value - shift;
      shifted_sum += delta;
      shifted_sum_squares += delta * delta;
    }
    const long double mean = m > 0 ?
      shift + shifted_sum / static_cast<long double>(m) : 0.0L;
    const long double mean_correction = m > 0 ?
      shifted_sum * shifted_sum / static_cast<long double>(m) : 0.0L;
    long double centered_sum_squares =
      shifted_sum_squares - mean_correction;
    const long double central_roundoff = 64.0L * LDBL_EPSILON *
      (fabsl(shifted_sum_squares) + fabsl(mean_correction) + 1.0L);
    if (centered_sum_squares < 0.0L &&
        centered_sum_squares >= -central_roundoff) {
      centered_sum_squares = 0.0L;
    }
    REAL(sums_)[col] = static_cast<double>(sum);
    REAL(sum_squares_)[col] = static_cast<double>(sum_squares);
    REAL(means_)[col] = static_cast<double>(mean);
    REAL(centered_sum_squares_)[col] =
      static_cast<double>(centered_sum_squares);
  }

  SEXP out_ = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out_, 0, sums_);
  SET_VECTOR_ELT(out_, 1, sum_squares_);
  SET_VECTOR_ELT(out_, 2, means_);
  SET_VECTOR_ELT(out_, 3, centered_sum_squares_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(names_, 0, mkChar("sum"));
  SET_STRING_ELT(names_, 1, mkChar("sum_squares"));
  SET_STRING_ELT(names_, 2, mkChar("mean"));
  SET_STRING_ELT(names_, 3, mkChar("centered_sum_squares"));
  setAttrib(out_, R_NamesSymbol, names_);
  UNPROTECT(6);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_centered_block_apply(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP row_means_, SEXP col_means_,
    SEXP rows_, SEXP columns_, SEXP X_, SEXP alpha_, SEXP beta_, SEXP Y_,
    SEXP transpose_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(row_means_) || !isReal(col_means_) ||
      !isLogical(rows_) || !isLogical(columns_) ||
      !isReal(X_) || !(isReal(Y_) || isNull(Y_)) || !isLogical(transpose_)) {
    error("invalid centered CSC block apply inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_centered_block_apply");

  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || (dimY == R_NilValue && !isNull(Y_))) {
    error("X and Y must be matrices");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = isNull(Y_) ? -1 : INTEGER(dimY)[0];
  const int yc = isNull(Y_) ? -1 : INTEGER(dimY)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const bool rows = LOGICAL(rows_)[0];
  const bool columns = LOGICAL(columns_)[0];
  const double alpha = REAL(alpha_)[0];
  const double beta = isNull(Y_) ? 0.0 : REAL(beta_)[0];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;

  if (xr != inner) {
    error("non-conformable X for centered CSC block apply");
  }
  if ((!isNull(Y_) && (yr != out_rows || yc != xc))) {
    error("non-conformable Y for centered CSC block apply");
  }
  if (rows && LENGTH(row_means_) != m) {
    error("row_means length must equal CSC row dimension");
  }
  if (columns && LENGTH(col_means_) != n) {
    error("col_means length must equal CSC column dimension");
  }

  SEXP out_ = PROTECT(block_apply_output(Y_, REALSXP, out_rows, xc, beta == 0.0));
  const int* row_idx = INTEGER(i_);
  const int* col_ptr = INTEGER(p_);
  const double* values = REAL(x_);
  const double* X = REAL(X_);
  const double* row_means = REAL(row_means_);
  const double* col_means = REAL(col_means_);
  double* out = REAL(out_);

  CSCOperator impl = {m, n, row_idx, col_ptr, values};
  csc_mark_per_call(&impl);
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_csc_centered_scaled_block_apply(
    SEXP i_, SEXP p_, SEXP x_, SEXP dim_, SEXP col_means_, SEXP weights_,
    SEXP X_, SEXP alpha_, SEXP beta_, SEXP Y_, SEXP transpose_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isInteger(i_) || !isInteger(p_) || !isReal(x_) || !isInteger(dim_) ||
      !isReal(col_means_) || !isReal(weights_) || !isReal(X_) ||
      !isReal(alpha_) || !isReal(beta_) || !(isReal(Y_) || isNull(Y_)) ||
      !isLogical(transpose_)) {
    error("invalid centered-scaled CSC block apply inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_centered_scaled_block_apply");
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || (dimY == R_NilValue && !isNull(Y_))) {
    error("X and Y must be matrices");
  }

  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = isNull(Y_) ? -1 : INTEGER(dimY)[0];
  const int yc = isNull(Y_) ? -1 : INTEGER(dimY)[1];
  const int out_rows = transpose ? n : m;
  const int inner = transpose ? m : n;
  if (LENGTH(col_means_) != n || LENGTH(weights_) != n || xr != inner ||
      (!isNull(Y_) && (yr != out_rows || yc != xc))) {
    error("non-conformable centered-scaled CSC block apply inputs");
  }

  const double beta = isNull(Y_) ? 0.0 : REAL(beta_)[0];
  SEXP out_ = PROTECT(block_apply_output(Y_, REALSXP, out_rows, xc, beta == 0.0));
  CenteredScaledCSCOperator impl = {
    {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)},
    REAL(col_means_),
    REAL(weights_)
  };
  csc_mark_per_call(&impl.base);
  const int status = eigencore_centered_scaled_csc_apply(
    &impl,
    transpose ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE,
    xc,
    REAL(X_),
    xr,
    REAL(alpha_)[0],
    beta,
    REAL(out_),
    out_rows,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error("centered-scaled CSC block apply", status);
  }
  UNPROTECT(1);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_diagonal_block_apply(SEXP x_, SEXP dim_, SEXP unit_,
                                               SEXP X_, SEXP alpha_, SEXP beta_,
                                               SEXP Y_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(x_) || !isInteger(dim_) || !isLogical(unit_) ||
      !isReal(X_) || !(isReal(Y_) || isNull(Y_))) {
    error("invalid diagonal block apply inputs");
  }
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || (dimY == R_NilValue && !isNull(Y_))) {
    error("X and Y must be matrices");
  }

  const int n = INTEGER(dim_)[0];
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int yr = isNull(Y_) ? -1 : INTEGER(dimY)[0];
  const int yc = isNull(Y_) ? -1 : INTEGER(dimY)[1];
  if (INTEGER(dim_)[1] != n || xr != n || (!isNull(Y_) && (yr != n || yc != xc))) {
    error("non-conformable diagonal block apply inputs");
  }

  const double beta = isNull(Y_) ? 0.0 : REAL(beta_)[0];
  SEXP out_ = PROTECT(block_apply_output(Y_, REALSXP, n, xc, beta == 0.0));
  DiagonalOperator impl = {n, REAL(x_), static_cast<bool>(LOGICAL(unit_)[0])};
  const int status = eigencore_diagonal_apply(
    &impl,
    EIGENCORE_TRANSPOSE_NONE,
    xc,
    REAL(X_),
    xr,
    REAL(alpha_)[0],
    beta,
    REAL(out_),
    n,
    nullptr
  );
  if (status != 0) {
    eigencore_apply_status_error("diagonal block apply", status);
  }
  UNPROTECT(1);
  return out_;
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_native_apply_noalloc_check(SEXP kind_, SEXP A_,
                                                     SEXP X_, SEXP Y_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_dense_apply_int_guard_check(void) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

extern "C" SEXP eigencore_col_norms(SEXP X_) {
  EIGENCORE_ENTRY_BEGIN
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
  EIGENCORE_ENTRY_END
}

// Test and benchmark hook for the cached CSC kernels (P8): applies
// Y <- alpha op(A) X + beta Y0 `reps` times (reps >= 1) on ONE operator, so
// the second and later applies use the cached CSR copy / row slabs exactly
// as a native solver does. With col_means/weights non-NULL the operator is
// the centered-scaled (A - 1 mu^T) D. Returns list(Y, ms_per_apply).
extern "C" SEXP eigencore_csc_apply_repeat(SEXP i_, SEXP p_, SEXP x_,
                                           SEXP dim_, SEXP col_means_,
                                           SEXP weights_, SEXP X_,
                                           SEXP alpha_, SEXP beta_, SEXP Y_,
                                           SEXP transpose_, SEXP reps_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isReal(X_) || !isReal(Y_) || !isLogical(transpose_) ||
      !isReal(alpha_) || !isReal(beta_)) {
    error("invalid CSC repeat-apply inputs");
  }
  eigencore_validate_csc_structure(i_, p_, x_, dim_, "csc_apply_repeat");
  const int m = INTEGER(dim_)[0];
  const int n = INTEGER(dim_)[1];
  const bool transpose = LOGICAL(transpose_)[0];
  const int reps = asInteger(reps_);
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  SEXP dimY = getAttrib(Y_, R_DimSymbol);
  if (dimX == R_NilValue || dimY == R_NilValue || reps == NA_INTEGER ||
      reps < 1) {
    error("invalid CSC repeat-apply inputs");
  }
  const int xr = INTEGER(dimX)[0];
  const int xc = INTEGER(dimX)[1];
  const int out_rows = transpose ? n : m;
  if (xr != (transpose ? m : n) || INTEGER(dimY)[0] != out_rows ||
      INTEGER(dimY)[1] != xc) {
    error("non-conformable CSC repeat-apply inputs");
  }
  const bool scaled = col_means_ != R_NilValue;
  if (scaled && (!isReal(col_means_) || !isReal(weights_) ||
                 XLENGTH(col_means_) != n || XLENGTH(weights_) != n)) {
    error("col_means and weights must be double vectors of length ncol");
  }
  CenteredScaledCSCOperator impl = {
    {m, n, INTEGER(i_), INTEGER(p_), REAL(x_)},
    scaled ? REAL(col_means_) : nullptr,
    scaled ? REAL(weights_) : nullptr
  };
  const size_t len = static_cast<size_t>(out_rows) * static_cast<size_t>(xc);
  SEXP out_ = PROTECT(allocMatrix(REALSXP, out_rows, xc));
  double* out = REAL(out_);
  const EigencoreTranspose op =
    transpose ? EIGENCORE_TRANSPOSE_ADJOINT : EIGENCORE_TRANSPOSE_NONE;
  double seconds = 0.0;
  for (int rep = 0; rep < reps; ++rep) {
    if (len > 0) {
      std::memcpy(out, REAL(Y_), sizeof(double) * len);
    }
    auto start = std::chrono::steady_clock::now();
    const int status = scaled ?
      eigencore_centered_scaled_csc_apply(&impl, op, xc, REAL(X_), xr,
                                          REAL(alpha_)[0], REAL(beta_)[0],
                                          out, out_rows, nullptr) :
      eigencore_csc_apply(&impl.base, op, xc, REAL(X_), xr, REAL(alpha_)[0],
                          REAL(beta_)[0], out, out_rows, nullptr);
    if (status != 0) {
      eigencore_apply_status_error("CSC repeat apply", status);
    }
    if (rep > 0 || reps == 1) {
      seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    }
  }
  SEXP ms_ = PROTECT(ScalarReal(1e3 * seconds / (reps > 1 ? reps - 1 : 1)));
  SEXP result_ = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result_, 0, out_);
  SET_VECTOR_ELT(result_, 1, ms_);
  SEXP names_ = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names_, 0, mkChar("Y"));
  SET_STRING_ELT(names_, 1, mkChar("ms"));
  setAttrib(result_, R_NamesSymbol, names_);
  UNPROTECT(4);
  return result_;
  EIGENCORE_ENTRY_END
}
