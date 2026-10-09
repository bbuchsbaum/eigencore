#include <cstring>
#include <new>
#include <stdexcept>
#include <vector>
#include "eigencore_common.h"
#include <R.h>
#include <Rinternals.h>

// Self-test hooks for the C10/P16 error/interrupt plumbing in
// eigencore_common.h. Each mode raises one kind of failure while a tracked C++
// object (counted in g_live_tracked) is alive; mode "live" reports how many
// tracked objects are still alive, so tests can verify that destructors ran
// before R regained control.

extern "C" {
int eigencore_unwind_trace_enabled = 0;
}

namespace {

int g_live_tracked = 0;

SEXP selftest_r_stop() {
  SEXP call = PROTECT(lang2(install("stop"),
                            mkString("eigencore selftest: R-level stop")));
  SEXP out = Rf_eval(call, R_BaseEnv);
  UNPROTECT(1);
  return out;
}

SEXP selftest_r_stop_body(void*) { return selftest_r_stop(); }
void selftest_noop_cleanup(void*, Rboolean) {}

struct Tracked {
  std::vector<double> buffer;
  Tracked() : buffer(1024, 0.0) { ++g_live_tracked; }
  Tracked(const Tracked&) = delete;
  Tracked& operator=(const Tracked&) = delete;
  ~Tracked() { --g_live_tracked; }
};

}  // namespace

extern "C" SEXP eigencore_unwind_selftest(SEXP mode_) {
  EIGENCORE_ENTRY_BEGIN
  if (!isString(mode_) || XLENGTH(mode_) != 1) {
    error("mode must be a single string");
  }
  const char* mode = CHAR(STRING_ELT(mode_, 0));
  if (std::strcmp(mode, "live") == 0) {
    return ScalarInteger(g_live_tracked);
  }
  // Stages of the r_stop path, for bisecting platform-specific failures.
  // They run before the tracked object exists: rerror_c_only and
  // rerror_continue_local deliberately leave without C++ cleanup.
  if (std::strcmp(mode, "trace_on") == 0 || std::strcmp(mode, "trace_off") == 0) {
    eigencore_unwind_trace_enabled = std::strcmp(mode, "trace_on") == 0;
    return ScalarLogical(TRUE);
  }
  if (std::strcmp(mode, "protect_ok") == 0) {
    return eigencore_unwind_protect([]() { return ScalarInteger(7); });
  }
  if (std::strcmp(mode, "rerror_c_only") == 0) {
    // R's own unwind protection with a no-op cleanup: R resumes the unwind
    // itself, no jump into C++.
    return R_UnwindProtect(selftest_r_stop_body, nullptr, selftest_noop_cleanup,
                           nullptr, eigencore_unwind_token());
  }
  if (std::strcmp(mode, "rerror_catch") == 0) {
    // Jump back into C++ and throw, then drop the R unwind (like a catch).
    try {
      eigencore_unwind_protect([]() { return selftest_r_stop(); });
    } catch (const eigencore::RUnwind&) {
      EIGENCORE_UNWIND_TRACE("rerror_catch: caught RUnwind locally");
      return ScalarLogical(TRUE);
    }
    return ScalarLogical(FALSE);
  }
  if (std::strcmp(mode, "rerror_continue_local") == 0) {
    SEXP token = R_NilValue;
    try {
      eigencore_unwind_protect([]() { return selftest_r_stop(); });
    } catch (const eigencore::RUnwind& unwind) {
      token = unwind.token;
    }
    EIGENCORE_UNWIND_TRACE("rerror_continue_local: R_ContinueUnwind");
    R_ContinueUnwind(token);
  }
  Tracked tracked;
  if (std::strcmp(mode, "error") == 0) {
    error("eigencore selftest: error with %d live C++ object(s)", g_live_tracked);
  }
  if (std::strcmp(mode, "bad_alloc") == 0) {
    throw std::bad_alloc();
  }
  if (std::strcmp(mode, "length_error") == 0) {
    std::vector<double> huge;
    huge.reserve(huge.max_size() + static_cast<size_t>(1));
  }
  if (std::strcmp(mode, "huge_vector") == 0) {
    // A size that cannot be satisfied: std::bad_alloc or std::length_error.
    std::vector<double> huge(static_cast<size_t>(1) << 60);
    return ScalarReal(huge[0]);
  }
  if (std::strcmp(mode, "interrupt") == 0) {
    throw eigencore::Interrupt();
  }
  if (std::strcmp(mode, "r_alloc_failure") == 0) {
    // R rejects an over-long vector with an error (a longjmp) before it tries
    // to allocate, on every platform; the unwind-protected allocVector turns
    // that longjmp into a C++ exception first. (Asking for R_XLEN_T_MAX
    // elements instead reaches the system allocator, whose failure mode and
    // message differ by platform.)
    SEXP big = PROTECT(allocVector(REALSXP, R_XLEN_T_MAX + static_cast<R_xlen_t>(1)));
    UNPROTECT(1);
    return big;
  }
  if (std::strcmp(mode, "r_stop") == 0) {
    return eigencore_unwind_protect([]() { return selftest_r_stop(); });
  }
  if (std::strcmp(mode, "probe") == 0) {
    // No interrupt is pending in a test run, so the probe must return.
    eigencore_check_interrupt_now();
    eigencore_check_interrupt();
    return ScalarLogical(TRUE);
  }
  if (std::strcmp(mode, "ok") == 0) {
    return ScalarInteger(g_live_tracked);
  }
  error("unknown selftest mode '%s'", mode);
  EIGENCORE_ENTRY_END
}
