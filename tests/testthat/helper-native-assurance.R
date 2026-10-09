# Native-code assurance hooks (docs/test-assurance.md, "Native code
# assurance").
#
# EIGENCORE_TEST_THREADS=<n> runs the whole suite with
# options(eigencore.threads = n), so the OpenMP and serial branches of the
# sparse kernels can both be exercised by the full suite (CI and the thread
# determinism check run it with 1 and 4).
local({
  threads <- Sys.getenv("EIGENCORE_TEST_THREADS", "")
  if (nzchar(threads)) {
    threads <- suppressWarnings(as.integer(threads))
    if (!is.na(threads) && threads >= 1L) {
      options(eigencore.threads = threads)
    }
  }
})

# TRUE when the process runs with AddressSanitizer preloaded (R itself is not
# instrumented; the package is built with -fsanitize=address and libasan is
# LD_PRELOADed). Memory-growth (RSS) checks are meaningless there because
# ASan quarantines freed blocks; LeakSanitizer covers leaks instead.
running_under_asan <- function() {
  grepl("libasan", Sys.getenv("LD_PRELOAD"), fixed = TRUE) ||
    nzchar(Sys.getenv("ASAN_OPTIONS"))
}
