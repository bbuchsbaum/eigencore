#' Multithreading in eigencore's sparse kernels
#'
#' eigencore parallelises its native sparse (CSC, `dgCMatrix`) operator
#' applications with OpenMP: `A^T X` products (a gather per column), `A X`
#' products (a gather over a cached CSR copy, or a scatter over cached
#' per-thread row slabs for tall matrices with short rows), the centred and
#' centred-scaled sparse operators, and the implicit normal operators
#' `A^T A` / `A A^T` behind sparse [svds()] / [svd_partial()]. Every output
#' entry is accumulated by one thread in the same order and with the same
#' arithmetic as the serial loop, so sparse products are bitwise identical for
#' any thread count. While such a solve runs, the Lanczos
#' reorthogonalisation also uses these threads (see below); its OpenMP kernels
#' do not depend on the thread count either, so results computed with one and
#' with several threads differ only by rounding.
#'
#' @section Options:
#' \describe{
#'   \item{`eigencore.threads`}{Number of threads, read at the start of every
#'     native call. When unset, a default fixed at package load is used: `1`
#'     under `R CMD check` (`_R_CHECK_PACKAGE_NAME_` set) or when
#'     `_R_CHECK_LIMIT_CORES_` is set (CRAN's core limit); otherwise the first
#'     value of `OMP_NUM_THREADS` when set; otherwise the number of processors
#'     available to OpenMP (respecting the CPU affinity mask), further limited
#'     to the physical cores (`parallel::detectCores(logical = FALSE)`) and to
#'     the cgroup CPU quota of a container, and capped at 8. Builds without
#'     OpenMP (for example Apple clang without `libomp`) always use one
#'     thread.}
#'   \item{`eigencore.adaptive_threads`}{`TRUE` (default) or `FALSE`. While
#'     enabled, eigencore watches the parallel efficiency of its own OpenMP
#'     regions and lowers the number of threads it actually uses when the
#'     machine is oversubscribed (see below). `FALSE` always uses
#'     `eigencore.threads` threads, as in eigencore 1.3.0.}
#'   \item{`eigencore.csr_cache_mb`}{Memory cap in MB (default 4096) for the
#'     per-operator CSR copy or row slabs (about 12 bytes per nonzero) that the
#'     parallel `A X` builds once per solve. When the copy would be larger,
#'     `A X` stays serial for blocks of up to 10 columns.}
#' }
#'
#' @section Busy machines:
#' OpenMP runtimes keep idle team threads spinning for a while (libgomp's
#' default `OMP_WAIT_POLICY`), which is fast on idle cores. When other
#' processes occupy the cores, every fork, join and barrier of a parallel
#' region instead waits a scheduler time slice for a team member that is not
#' running, and eigencore 1.3.0's short sparse kernels became 10-80 times
#' slower with 4 threads than with 1. eigencore therefore times each parallel
#' region: when a region takes longer than the master thread's share of the
#' work times the team size (i.e. longer than running it serially would
#' have), twice in a row, the team is halved, and a larger team is tried
#' again after a pause that grows (0.25 s to 4 s) while the machine stays
#' busy. Only the number of threads in a region changes; the algorithm and
#' the work partition that determine rounding are chosen from
#' `eigencore.threads`, so results do not depend on the machine load. The
#' wait policy itself cannot be changed from a package (R loads the OpenMP
#' runtime at start-up); on machines that are routinely oversubscribed,
#' setting `OMP_WAIT_POLICY=passive` in the environment before starting R
#' also removes the spinning. A multithreaded OpenBLAS (pthreads build)
#' spins the same way in its own dense calls (including those made between
#' eigencore's native calls and by CHOLMOD); eigencore does not change it
#' on its own because OpenBLAS results depend on its thread count, so on a
#' busy machine also limit it (`OPENBLAS_NUM_THREADS=1`, or
#' `RhpcBLASctl::blas_set_num_threads(1)`).
#'
#' @section Interaction with a multithreaded BLAS:
#' OpenBLAS (pthreads build) and FlexiBLAS keep idle worker threads spinning
#' after each call; OpenMP threads started meanwhile would compete with them
#' for cores. When a native sparse solve (for example `eigs_sym()`, `svds()` or
#' `svd_partial()` on a `dgCMatrix`) runs a multithreaded sparse kernel,
#' eigencore therefore switches such a BLAS to one thread for the rest of that
#' call and restores the previous setting when the call returns (or fails);
#' its reorthogonalisation then runs on eigencore's OpenMP threads. Sparse
#' products applied one at a time from R-level operators (for example inside
#' a matrix-free `center()` solve) stay serial while such a BLAS is
#' multithreaded, so the surrounding BLAS work keeps its threads. An
#' OpenMP-built OpenBLAS shares eigencore's thread pool; MKL, BLIS and other
#' libraries are not changed. With those, or when running several R workers
#' in parallel, consider `options(eigencore.threads = 1)` or limiting the BLAS
#' threads (for example with `RhpcBLASctl::blas_set_num_threads()`).
#'
#' @name eigencore-threads
#' @aliases eigencore.threads eigencore.adaptive_threads eigencore.csr_cache_mb
#' @examples
#' old <- options(eigencore.threads = 1)
#' getOption("eigencore.threads")
#' options(old)
NULL

# Query (and optionally set) the thread count: returns the count in effect
# with attributes "default" (load-time default) and "openmp".
eigencore_threads <- function(n = NULL) {
  if (!is.null(n)) {
    if (length(n) != 1L) {
      stop("`n` must be NULL, NA, or a single positive integer.", call. = FALSE)
    }
    if (is.na(n)) {
      options(eigencore.threads = NULL)
    } else {
      if (!is.numeric(n) || n < 1 || n != round(n)) {
        stop("`n` must be NULL, NA, or a single positive integer.", call. = FALSE)
      }
      options(eigencore.threads = as.integer(n))
    }
  }
  info <- native_thread_info()
  out <- structure(
    info[["current"]],
    default = info[["default"]],
    openmp = isTRUE(info[["openmp"]] == 1L)
  )
  if (is.null(n)) out else invisible(out)
}

native_thread_info <- function() {
  .Call("eigencore_thread_info", PACKAGE = "eigencore")
}

# Test hook: set the parallel-efficiency governor's team-size cap (0 removes
# it, NA leaves it) and return native_thread_info().
native_thread_info_governor <- function(cap = NA_real_) {
  .Call("eigencore_thread_governor", as.numeric(cap), PACKAGE = "eigencore")
}

# Load-time default thread count (see ?`eigencore-threads`).
eigencore_default_threads <- function(env = Sys.getenv(), processors = NULL,
                                      cores = NULL, quota = NULL) {
  get_env <- function(name) {
    value <- env[name]
    if (is.na(value)) "" else unname(value)
  }
  limit_cores <- get_env("_R_CHECK_LIMIT_CORES_")
  if (nzchar(limit_cores) && !identical(tolower(limit_cores), "false")) {
    return(1L)
  }
  if (nzchar(get_env("_R_CHECK_PACKAGE_NAME_"))) {
    return(1L)
  }
  omp <- get_env("OMP_NUM_THREADS")
  if (nzchar(omp)) {
    first <- suppressWarnings(as.integer(strsplit(omp, ",", fixed = TRUE)[[1L]][1L]))
    if (!is.na(first) && first >= 1L) {
      return(min(first, 256L))
    }
  }
  if (is.null(processors)) {
    processors <- native_thread_info()[["processors"]]
  }
  processors <- suppressWarnings(as.integer(processors))
  if (length(processors) != 1L || is.na(processors) || processors < 1L) {
    return(1L)
  }
  # Hyperthreads add little to these memory-bound kernels, and a container's
  # CPU quota caps what the process can use regardless of the visible cores.
  if (is.null(cores)) {
    cores <- tryCatch(parallel::detectCores(logical = FALSE),
                      error = function(e) NA_integer_)
  }
  if (is.null(quota)) {
    quota <- eigencore_cgroup_cpu_quota()
  }
  limit <- processors
  cores <- suppressWarnings(as.integer(cores))
  if (length(cores) == 1L && !is.na(cores) && cores >= 1L) {
    limit <- min(limit, cores)
  }
  quota <- suppressWarnings(as.numeric(quota))
  if (length(quota) == 1L && is.finite(quota) && quota > 0) {
    limit <- min(limit, max(1L, as.integer(floor(quota))))
  }
  min(as.integer(limit), 8L)
}

# CPU quota of this process's cgroup in CPUs (cgroup v2 cpu.max, else v1
# cpu.cfs_quota_us / cpu.cfs_period_us), or NA when unlimited or unknown.
eigencore_cgroup_cpu_quota <- function(root = "/sys/fs/cgroup",
                                       self = "/proc/self/cgroup") {
  read1 <- function(path) {
    out <- tryCatch(suppressWarnings(readLines(path, n = 1L, warn = FALSE)),
                    error = function(e) character())
    if (length(out) == 1L) trimws(out) else NA_character_
  }
  lines <- tryCatch(suppressWarnings(readLines(self, warn = FALSE)),
                    error = function(e) character())
  path_for <- function(pattern) {
    hit <- grep(pattern, lines, value = TRUE, perl = TRUE)
    if (length(hit)) sub("^[^:]*:[^:]*:", "", hit[[1L]]) else NA_character_
  }
  # The process's own cgroup first, then the namespace root (containers).
  dirs <- function(sub, rel) {
    base <- file.path(root, sub)
    unique(c(if (!is.na(rel) && nzchar(rel) && rel != "/") paste0(base, rel), base))
  }
  for (dir in dirs("", path_for("^0::"))) {
    value <- read1(file.path(dir, "cpu.max"))
    if (!is.na(value)) {
      parts <- strsplit(value, "[[:space:]]+")[[1L]]
      if (length(parts) >= 2L && parts[[1L]] != "max") {
        quota <- suppressWarnings(as.numeric(parts[[1L]]) / as.numeric(parts[[2L]]))
        if (is.finite(quota) && quota > 0) return(quota)
      }
      return(NA_real_)
    }
  }
  for (sub in c("cpu", "cpu,cpuacct", "cpuacct,cpu")) {
    for (dir in dirs(sub, path_for("^[0-9]+:[^:]*\\bcpu\\b[^:]*:"))) {
      quota <- suppressWarnings(as.numeric(read1(file.path(dir, "cpu.cfs_quota_us"))))
      period <- suppressWarnings(as.numeric(read1(file.path(dir, "cpu.cfs_period_us"))))
      if (is.finite(quota) && is.finite(period) && period > 0) {
        return(if (quota > 0) quota / period else NA_real_)
      }
    }
  }
  NA_real_
}

# Test/benchmark hook: applies op(A) (or the centered-scaled operator when
# col_means/weights are given) `reps` times on one native operator, so the
# later applies use the cached parallel structures. Returns list(Y, ms).
csc_apply_repeat <- function(A, X, alpha = 1, beta = 0, Y = NULL,
                             transpose = FALSE, reps = 2L,
                             col_means = NULL, weights = NULL) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (is.null(Y)) {
    Y <- matrix(0, if (transpose) ncol(A) else nrow(A), ncol(X))
  }
  storage.mode(Y) <- "double"
  .Call(
    "eigencore_csc_apply_repeat",
    methods::slot(A, "i"), methods::slot(A, "p"), methods::slot(A, "x"),
    methods::slot(A, "Dim"),
    if (is.null(col_means)) NULL else as.numeric(col_means),
    if (is.null(weights)) NULL else as.numeric(weights),
    X, as.numeric(alpha), as.numeric(beta), Y, isTRUE(transpose),
    as.integer(reps),
    PACKAGE = "eigencore"
  )
}

.onLoad <- function(libname, pkgname) {
  default <- tryCatch(eigencore_default_threads(), error = function(e) 1L)
  .Call("eigencore_set_default_threads", as.integer(default), PACKAGE = "eigencore")
  invisible()
}
