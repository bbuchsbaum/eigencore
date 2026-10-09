#ifndef EIGENCORE_COMMON_H
#define EIGENCORE_COMMON_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <chrono>
#include <cmath>
#include <climits>
#include <csetjmp>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <new>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

// OpenMP helpers (P8). Every pragma in src/ goes through EIGENCORE_OMP so the
// package compiles warning-free and runs serially where OpenMP is unavailable
// (Apple clang, some Windows toolchains): SHLIB_OPENMP_CXXFLAGS is empty there
// and _OPENMP is undefined. Parallel regions never call the R API and never
// let a C++ exception escape.
#ifdef _OPENMP
#define EIGENCORE_OMP(directive) _Pragma(#directive)
#else
#define EIGENCORE_OMP(directive)
#endif

static inline int eigencore_omp_thread_num() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

static inline int eigencore_omp_num_threads() {
#ifdef _OPENMP
  return omp_get_num_threads();
#else
  return 1;
#endif
}

// ---------------------------------------------------------------------------
// Error, interrupt and R-unwind plumbing (review items C10 / P16)
// ---------------------------------------------------------------------------
//
// R signals errors with longjmp, which skips C++ destructors: an Rf_error()
// raised while a std::vector (or any RAII object) is alive leaks it, and a C++
// exception escaping an extern "C" .Call entry point terminates R. eigencore
// therefore uses one mechanism everywhere in src/:
//
//  * Every .Call entry point registered in init.cpp wraps its body in
//    EIGENCORE_ENTRY_BEGIN ... EIGENCORE_ENTRY_END, i.e. runs it as a lambda
//    inside eigencore_call(). eigencore_call catches every C++ exception,
//    copies the message into a stack buffer, lets the body's destructors run
//    by leaving the try scope, and only then raises the R condition
//    (Rf_error / R interrupt condition / R_ContinueUnwind).
//  * Inside eigencore sources `error(...)` (R's macro for Rf_error) is
//    remapped to eigencore_raise_error(), which formats the message exactly
//    as Rf_error would and throws eigencore::Error. Message texts are
//    therefore unchanged. Never call Rf_error() directly below an entry
//    point; use error().
//  * allocVector / allocMatrix / duplicate are remapped to unwind-protected
//    versions (small allocations such as ScalarReal/mkChar are not): an
//    R allocation failure (which longjmps) is converted into an
//    eigencore::RUnwind exception carrying the R continuation token, and the
//    entry wrapper resumes R's unwind with R_ContinueUnwind after C++ cleanup.
//    eigencore_unwind_protect(lambda) does the same for any other R API call
//    that may longjmp while C++ objects are alive.
//  * Long native loops call eigencore_check_interrupt() (throttled to ~10 Hz).
//    It probes R_CheckUserInterrupt inside R_ToplevelExec (which returns
//    instead of longjmp-ing) and throws eigencore::Interrupt; the entry
//    wrapper then signals an R "interrupt" condition (so
//    tryCatch(interrupt = ) sees it) and, if no handler exits, raises the
//    error "eigencore: computation interrupted by user".
//
// Everything in this header that throws must only run below an entry point
// (never from finalizers or R_init).

namespace eigencore {

class Error : public std::runtime_error {
 public:
  explicit Error(const char* message) : std::runtime_error(message) {}
};

// Not derived from std::exception so that no generic handler swallows it.
struct Interrupt {};

struct RUnwind {
  SEXP token;
};

}  // namespace eigencore

[[noreturn]] static inline void eigencore_raise_error(const char* format, ...)
#if defined(__GNUC__)
  __attribute__((format(printf, 1, 2)))
#endif
  ;

[[noreturn]] static inline void eigencore_raise_error(const char* format, ...) {
  // Rf_error formats into an 8192-byte buffer; match it.
  char buffer[8192];
  va_list args;
  va_start(args, format);
  std::vsnprintf(buffer, sizeof(buffer), format, args);
  va_end(args);
  throw eigencore::Error(buffer);
}

#ifdef error
#undef error
#endif
#define error(...) eigencore_raise_error(__VA_ARGS__)

static inline SEXP eigencore_unwind_token() {
  static SEXP token = nullptr;
  if (token == nullptr) {
    token = R_MakeUnwindCont();
    R_PreserveObject(token);
  }
  return token;
}

template <typename F>
static SEXP eigencore_unwind_protect_body(void* data) {
  return (*static_cast<F*>(data))();
}

static inline void eigencore_unwind_protect_cleanup(void* jmpbuf, Rboolean jump) {
  if (jump) {
    std::longjmp(*static_cast<std::jmp_buf*>(jmpbuf), 1);
  }
}

// Run `code` (a callable returning SEXP that uses the R API) so that an R
// longjmp out of it becomes an eigencore::RUnwind C++ exception. Only R frames
// are skipped by the longjmp; C++ frames are unwound by the exception.
template <typename F>
static inline SEXP eigencore_unwind_protect(F code) {
  SEXP token = eigencore_unwind_token();
  std::jmp_buf jmpbuf;
  if (setjmp(jmpbuf)) {
    throw eigencore::RUnwind{token};
  }
  return R_UnwindProtect(eigencore_unwind_protect_body<F>, &code,
                         eigencore_unwind_protect_cleanup, &jmpbuf, token);
}

static inline SEXP eigencore_alloc_vector(SEXPTYPE type, R_xlen_t length) {
  return eigencore_unwind_protect([type, length]() {
    return Rf_allocVector(type, length);
  });
}

static inline SEXP eigencore_alloc_matrix(SEXPTYPE type, int rows, int cols) {
  return eigencore_unwind_protect([type, rows, cols]() {
    return Rf_allocMatrix(type, rows, cols);
  });
}

static inline SEXP eigencore_duplicate(SEXP x) {
  return eigencore_unwind_protect([x]() {
    return Rf_duplicate(x);
  });
}

#ifdef duplicate
#undef duplicate
#endif
#define duplicate eigencore_duplicate
#ifdef allocVector
#undef allocVector
#endif
#define allocVector eigencore_alloc_vector
#ifdef allocMatrix
#undef allocMatrix
#endif
#define allocMatrix eigencore_alloc_matrix

static inline void eigencore_interrupt_probe(void*) {
  R_CheckUserInterrupt();
}

// Unthrottled check: R_ToplevelExec returns FALSE when the probe was left by
// a longjmp (a pending user interrupt), without unwinding our frames.
static inline void eigencore_check_interrupt_now() {
  if (R_ToplevelExec(eigencore_interrupt_probe, nullptr) == FALSE) {
    throw eigencore::Interrupt();
  }
}

// Throttled check for hot loops: probes at most every ~100 ms, so calling it
// once per iteration costs one steady_clock read.
static inline void eigencore_check_interrupt() {
  static std::chrono::steady_clock::time_point last =
    std::chrono::steady_clock::now();
  const std::chrono::steady_clock::time_point now =
    std::chrono::steady_clock::now();
  if (now - last < std::chrono::milliseconds(100)) {
    return;
  }
  last = now;
  eigencore_check_interrupt_now();
}

[[noreturn]] static inline void eigencore_signal_interrupt() {
  // Signal a classed "interrupt" condition so calling/exiting handlers for
  // interrupts run, then fall back to an ordinary R error.
  SEXP cond = PROTECT(Rf_allocVector(VECSXP, 0));
  SEXP klass = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(klass, 0, Rf_mkChar("interrupt"));
  SET_STRING_ELT(klass, 1, Rf_mkChar("condition"));
  Rf_setAttrib(cond, R_ClassSymbol, klass);
  SEXP call = PROTECT(Rf_lang2(Rf_install("signalCondition"), cond));
  Rf_eval(call, R_BaseEnv);
  UNPROTECT(3);
  Rf_error("eigencore: computation interrupted by user");
}

// Thread count for the OpenMP sparse kernels (P8), defined in
// native_operators.cpp. eigencore_call_enter() re-reads
// getOption("eigencore.threads") (falling back to the package default set at
// load time) once per .Call, on the main thread, so kernels only read a plain
// int; eigencore_thread_count() is always >= 1 and is 1 when the package was
// built without OpenMP. The first multithreaded sparse kernel of a call may
// switch a spinning-thread BLAS (OpenBLAS pthreads, FlexiBLAS) to one thread;
// eigencore_call_leave() restores it when the outermost .Call returns or
// raises, so BLAS and OpenMP threads never compete for cores.
extern "C" void eigencore_refresh_thread_count(void);
extern "C" int eigencore_thread_count(void);
extern "C" void eigencore_call_enter(void);
extern "C" void eigencore_call_leave(void);
// Threads for eigencore's own OpenMP dense helpers: > 1 only while BLAS runs
// single-threaded (see native_operators.cpp).
extern "C" int eigencore_reorth_threads(void);

// Body of every .Call entry point. All C++ objects created by `body` are
// destroyed before any R condition is raised.
template <typename Body>
static inline SEXP eigencore_call(Body&& body) {
  eigencore_call_enter();
  enum { kError, kInterrupt, kUnwind } kind = kError;
  char message[8192];
  message[0] = '\0';
  SEXP token = R_NilValue;
  try {
    SEXP result = body();
    eigencore_call_leave();
    return result;
  } catch (const eigencore::Interrupt&) {
    kind = kInterrupt;
  } catch (const eigencore::RUnwind& unwind) {
    kind = kUnwind;
    token = unwind.token;
  } catch (const std::bad_alloc&) {
    std::snprintf(message, sizeof(message),
                  "eigencore: memory allocation failed (std::bad_alloc)");
  } catch (const std::exception& e) {
    std::snprintf(message, sizeof(message), "%s", e.what());
  } catch (...) {
    std::snprintf(message, sizeof(message), "eigencore: unknown C++ exception");
  }
  eigencore_call_leave();
  if (kind == kUnwind) {
    R_ContinueUnwind(token);
  }
  if (kind == kInterrupt) {
    eigencore_signal_interrupt();
  }
  Rf_error("%s", message);
}

#define EIGENCORE_ENTRY_BEGIN return eigencore_call([&]() -> SEXP {
#define EIGENCORE_ENTRY_END });

// CSC slot validation (defined in native_operators.cpp); raises error() naming
// `context` when p/i/x/Dim are inconsistent. Pass R_NilValue for i_ to skip
// the row-index checks.
extern "C" void eigencore_validate_csc_structure(SEXP i_, SEXP p_, SEXP x_,
                                                 SEXP dim_, const char* context);

// Element count for an owning std::vector scratch buffer: at least one element
// so .data() is never null (kernels historically null-checked malloc results).
template <typename T>
static inline size_t eigencore_buffer_size(T count) {
  return (count > 0) ? static_cast<size_t>(count) : static_cast<size_t>(1);
}

struct NativeBlockStageSeconds {
  double apply = 0.0;
  double recurrence = 0.0;
  double reorthogonalization = 0.0;
  double projected_solve = 0.0;
  double projection_update = 0.0;
  double projection_copy = 0.0;
  double projected_eigensolve = 0.0;
  double selected_vector_copy = 0.0;
  double ritz_residual = 0.0;
  double ritz_vector_form = 0.0;
  double ritz_operator_apply = 0.0;
  double ritz_norm = 0.0;
  double ritz_final_polish = 0.0;
  double locking = 0.0;
  double restart = 0.0;
};

struct NativeBlockRestartHistory {
  int capacity = 0;
  int length = 0;
  int* restart = nullptr;
  int* m_active = nullptr;
  int* selected_count = nullptr;
  int* locked_before = nullptr;
  int* locked_after = nullptr;
  int* nconv_wanted = nullptr;
  double* max_residual = nullptr;
  double* max_backward_error = nullptr;
};

static inline std::chrono::steady_clock::time_point native_timer_now() {
  return std::chrono::steady_clock::now();
}

static inline double native_timer_elapsed(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

static inline int eigencore_int_indexable(int64_t value) {
  return value >= 0 && value <= static_cast<int64_t>(INT_MAX);
}

// Status returned by the scalar Krylov kernels when the start vector or a
// recurrence coefficient is NaN/Inf.
#define EIGENCORE_STATUS_NONFINITE (-10)

static inline void eigencore_check_nonfinite_status(int status) {
  if (status == EIGENCORE_STATUS_NONFINITE) {
    error("non-finite value encountered (check input for NA/NaN/Inf)");
  }
}

static inline void eigencore_apply_status_error(const char* context, int status) {
  eigencore_check_nonfinite_status(status);
  if (status == -2) {
    error("%s failed: dimensions exceed LP64 BLAS/R integer range; LAPACK64 is not enabled",
          context);
  }
  error("%s failed with status=%d", context, status);
}

// Dot product / Euclidean norm through BLAS (no long double accumulators: those
// defeat vectorisation and are software binary128 on aarch64 Linux). The norm
// takes the fast sqrt(ddot) path when the sum of squares is safely inside the
// representable range and otherwise falls back to the scaled dnrm2, so
// operators scaled by 1e+-150 cannot overflow or underflow the norm.
static inline double ec_dot(const double* x, const double* y, int n) {
  if (n <= 0) {
    return 0.0;
  }
  int one = 1;
  return F77_CALL(ddot)(&n, x, &one, y, &one);
}

static inline double ec_norm2(const double* x, int n) {
  if (n <= 0) {
    return 0.0;
  }
  int one = 1;
  const double ss = F77_CALL(ddot)(&n, x, &one, x, &one);
  if (ss >= 1e-280 && ss <= 1e280) {
    return std::sqrt(ss);
  }
  if (ss == 0.0) {
    return 0.0;
  }
  return F77_CALL(dnrm2)(&n, x, &one);
}

// ||av - theta * v||_2 without a temporary and without long double: four
// independent double accumulators (keeps the FP pipeline busy), with a scaled
// two-pass fallback when the sum of squares leaves the safe range.
static inline double ec_residual_norm(const double* av, const double* v,
                                      double theta, int n) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int i = 0;
  for (; i + 3 < n; i += 4) {
    const double d0 = av[i] - theta * v[i];
    const double d1 = av[i + 1] - theta * v[i + 1];
    const double d2 = av[i + 2] - theta * v[i + 2];
    const double d3 = av[i + 3] - theta * v[i + 3];
    s0 += d0 * d0;
    s1 += d1 * d1;
    s2 += d2 * d2;
    s3 += d3 * d3;
  }
  for (; i < n; ++i) {
    const double d = av[i] - theta * v[i];
    s0 += d * d;
  }
  const double ss = (s0 + s1) + (s2 + s3);
  if ((ss >= 1e-280 && ss <= 1e280) || ss == 0.0 || ISNAN(ss)) {
    return std::sqrt(ss);
  }
  double amax = 0.0;
  for (i = 0; i < n; ++i) {
    const double d = std::fabs(av[i] - theta * v[i]);
    if (d > amax) {
      amax = d;
    }
  }
  if (!(amax > 0.0) || !R_FINITE(amax)) {
    return amax;
  }
  double scaled = 0.0;
  for (i = 0; i < n; ++i) {
    const double d = (av[i] - theta * v[i]) / amax;
    scaled += d * d;
  }
  return amax * std::sqrt(scaled);
}

static inline void combine_basis_columns_small(const double* basis,
                                               int n,
                                               int basis_cols,
                                               const double* coeff,
                                               int coeff_ld,
                                               int out_cols,
                                               double* out) {
  std::memset(out, 0,
              sizeof(double) * static_cast<size_t>(n) *
                static_cast<size_t>(out_cols));
  for (int p = 0; p < out_cols; ++p) {
    double* y = out + static_cast<int64_t>(p) * n;
    for (int col = 0; col < basis_cols; ++col) {
      const double a = coeff[col + static_cast<int64_t>(p) * coeff_ld];
      if (a == 0.0) {
        continue;
      }
      const double* x = basis + static_cast<int64_t>(col) * n;
      for (int row = 0; row < n; ++row) {
        y[row] += a * x[row];
      }
    }
  }
}

// out (n x out_cols) = basis (n x basis_cols) * coeff (basis_cols x out_cols).
// Uses BLAS-3 dgemm unless the whole product is tiny; the decision is driven
// by the total work n * basis_cols * out_cols, not by out_cols alone.
static inline void combine_basis_columns(const double* basis,
                                         int n,
                                         int basis_cols,
                                         const double* coeff,
                                         int coeff_ld,
                                         int out_cols,
                                         double* out) {
  if (out_cols <= 0 || n <= 0) {
    return;
  }
  const int64_t work = static_cast<int64_t>(n) * basis_cols * out_cols;
  if (basis_cols <= 0 || work <= 4096) {
    combine_basis_columns_small(basis, n, basis_cols, coeff, coeff_ld,
                                out_cols, out);
    return;
  }
  const char trans_N = 'N';
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemm)(&trans_N, &trans_N, &n, &out_cols, &basis_cols,
                  &one, basis, &n, coeff, &coeff_ld,
                  &zero, out, &n FCONE FCONE);
}

#endif
