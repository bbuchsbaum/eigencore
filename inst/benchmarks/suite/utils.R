# Benchmark suite: small utilities (argument parsing, timing, memory probes).
#
# Sourced by inst/benchmarks/run-suite.R. Base R only.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

suite_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) next
    a <- substring(a, 3L)
    if (grepl("=", a, fixed = TRUE)) {
      key <- sub("=.*$", "", a)
      val <- sub("^[^=]*=", "", a)
    } else {
      key <- a
      val <- "true"
    }
    out[[gsub("-", "_", key)]] <- val
  }
  out
}

suite_csv_arg <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

suite_flag <- function(x) !is.null(x) && tolower(x) %in% c("true", "1", "yes")

suite_now <- function() as.numeric(Sys.time())

# CPU time (user + system) of this process and its finished children. On a
# shared machine it is a steadier cost measure than wall time for
# single-threaded runs (with threads > 1 it sums over threads).
suite_cpu <- function() {
  p <- proc.time()
  sum(p[c("user.self", "sys.self")], na.rm = TRUE)
}

suite_loadavg <- function() {
  if (file.exists("/proc/loadavg")) {
    x <- tryCatch(scan("/proc/loadavg", what = "", n = 3L, quiet = TRUE),
                  error = function(e) character())
    if (length(x) == 3L) return(as.numeric(x))
  }
  # macOS / BSD
  x <- tryCatch(suppressWarnings(system2("sysctl", c("-n", "vm.loadavg"),
                                         stdout = TRUE, stderr = FALSE)),
                error = function(e) character())
  if (length(x)) {
    v <- suppressWarnings(as.numeric(strsplit(gsub("[{}]", "", x[[1L]]), "\\s+")[[1L]]))
    v <- v[is.finite(v)]
    if (length(v) >= 3L) return(v[1:3])
  }
  c(NA_real_, NA_real_, NA_real_)
}

# Peak resident set size (Linux only). Writing "5" to clear_refs resets the
# high-water mark, so VmHWM after a call measures the peak *process* memory of
# that call (including C/C++ heaps that R's allocator never sees).
suite_rss_reset <- function() {
  if (!file.exists("/proc/self/clear_refs")) return(FALSE)
  isTRUE(tryCatch({
    writeLines("5", "/proc/self/clear_refs")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE))
}

suite_rss_field <- function(field) {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  x <- tryCatch(readLines("/proc/self/status"), error = function(e) character())
  line <- x[startsWith(x, paste0(field, ":"))]
  if (!length(line)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1L]])) / 1024 # MiB
}

suite_md5_string <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(x, f, useBytes = TRUE)
  unname(tools::md5sum(f))
}

suite_log <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  utils::flush.console()
}

suite_fmt_time <- function(x) {
  ifelse(!is.finite(x), "-",
         ifelse(x < 1, sprintf("%.0f ms", x * 1000), sprintf("%.2f s", x)))
}

# Deterministic RNG: every case and every method call sets this explicitly.
suite_set_seed <- function(seed) {
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  set.seed(as.integer(seed))
}
