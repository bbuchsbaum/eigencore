#include <cstring>
#include <new>
#include <stdexcept>
#include <vector>
#include <R.h>
#include <Rinternals.h>
#include "eigencore_common.h"

// Self-test hooks for the C10/P16 error/interrupt plumbing in
// eigencore_common.h. Each mode raises one kind of failure while a tracked C++
// object (counted in g_live_tracked) is alive; mode "live" reports how many
// tracked objects are still alive, so tests can verify that destructors ran
// before R regained control.

namespace {

int g_live_tracked = 0;

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
    // R signals "cannot allocate vector ..." by longjmp; the unwind-protected
    // allocVector turns it into a C++ exception first.
    SEXP big = PROTECT(allocVector(REALSXP, R_XLEN_T_MAX));
    UNPROTECT(1);
    return big;
  }
  if (std::strcmp(mode, "r_stop") == 0) {
    return eigencore_unwind_protect([]() {
      SEXP call = PROTECT(lang2(install("stop"),
                                mkString("eigencore selftest: R-level stop")));
      SEXP out = eval(call, R_BaseEnv);
      UNPROTECT(1);
      return out;
    });
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
