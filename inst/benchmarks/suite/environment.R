# Benchmark suite: environment capture.
#
# Recorded once per run (environment.json / environment.rds) and per worker
# (BLAS/OpenMP thread counts actually in effect). The hostname is never stored;
# only a short hash of hostname + CPU model identifies the machine.

suite_cpu_model <- function() {
  if (file.exists("/proc/cpuinfo")) {
    x <- readLines("/proc/cpuinfo", warn = FALSE)
    m <- x[grepl("^(model name|Hardware|Processor|cpu model)\\s*:", x)]
    if (length(m)) return(trimws(sub("^[^:]*:", "", m[[1L]])))
  }
  x <- tryCatch(suppressWarnings(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                                         stdout = TRUE, stderr = FALSE)),
                error = function(e) character())
  if (length(x) && nzchar(x[[1L]])) return(x[[1L]])
  NA_character_
}

suite_mem_total_gb <- function() {
  if (file.exists("/proc/meminfo")) {
    x <- readLines("/proc/meminfo", warn = FALSE)
    m <- x[startsWith(x, "MemTotal:")]
    if (length(m)) return(round(as.numeric(gsub("[^0-9]", "", m)) / 1024^2, 1))
  }
  x <- tryCatch(suppressWarnings(system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE)),
                error = function(e) character())
  if (length(x)) return(round(as.numeric(x[[1L]]) / 1024^3, 1))
  NA_real_
}

suite_machine_id <- function() {
  host <- Sys.info()[["nodename"]] %||% "unknown"
  substr(suite_md5_string(paste(host, suite_cpu_model(), sep = "|")), 1L, 8L)
}

suite_blas_vendor <- function(path) {
  p <- tolower(paste(path, collapse = " "))
  if (grepl("openblas", p)) "OpenBLAS"
  else if (grepl("mkl", p)) "Intel MKL"
  else if (grepl("flexiblas", p)) "FlexiBLAS"
  else if (grepl("accelerate|veclib", p)) "Apple Accelerate"
  else if (grepl("atlas", p)) "ATLAS"
  else if (grepl("blis", p)) "BLIS"
  else if (grepl("librblas|/blas/libblas|libblas\\.so", p)) "reference BLAS (or alternatives link)"
  else NA_character_
}

suite_git_info <- function(repo_dir) {
  git <- Sys.which("git")
  if (!nzchar(git) || is.null(repo_dir) || !dir.exists(repo_dir)) return(list(sha = NA_character_, dirty = NA))
  sha <- tryCatch(suppressWarnings(system2(git, c("-C", shQuote(repo_dir), "rev-parse", "--short=10", "HEAD"),
                                           stdout = TRUE, stderr = FALSE)), error = function(e) character())
  if (!length(sha) || !is.null(attr(sha, "status"))) return(list(sha = NA_character_, dirty = NA))
  st <- tryCatch(suppressWarnings(system2(git, c("-C", shQuote(repo_dir), "status", "--porcelain",
                                                 "--untracked-files=no", "--", "R", "src", "DESCRIPTION"),
                                          stdout = TRUE, stderr = FALSE)), error = function(e) character())
  list(sha = sha[[1L]], dirty = length(st) > 0L)
}

suite_pkg_version <- function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
}

suite_blas_threads <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    return(list(blas = tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_),
                omp = tryCatch(RhpcBLASctl::omp_get_max_threads(), error = function(e) NA_integer_),
                source = "RhpcBLASctl"))
  }
  list(blas = suppressWarnings(as.integer(Sys.getenv("OPENBLAS_NUM_THREADS", NA))),
       omp = suppressWarnings(as.integer(Sys.getenv("OMP_NUM_THREADS", NA))),
       source = "environment variables")
}

suite_environment <- function(repo_dir = NULL, profile = NA_character_, args = list()) {
  si <- utils::sessionInfo()
  blas_path <- si$BLAS %||% NA_character_
  lapack_path <- si$LAPACK %||% NA_character_
  la_lib <- tryCatch(La_library(), error = function(e) NA_character_)
  ext <- tryCatch(extSoftVersion(), error = function(e) character())
  ec_desc <- tryCatch(utils::packageDescription("eigencore"), error = function(e) NULL)
  git <- suite_git_info(repo_dir)
  pkgs <- c("eigencore", "Matrix", "RSpectra", "irlba", "PRIMME", "RhpcBLASctl", "bench")
  list(
    timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    profile = profile,
    machine_id = suite_machine_id(),
    r_version = R.version.string,
    platform = R.version$platform,
    os = utils::osVersion %||% Sys.info()[["sysname"]],
    cpu_model = suite_cpu_model(),
    cores_logical = parallel::detectCores(logical = TRUE),
    cores_physical = parallel::detectCores(logical = FALSE),
    mem_total_gb = suite_mem_total_gb(),
    blas_path = blas_path,
    lapack_path = lapack_path,
    la_library = la_lib,
    la_version = tryCatch(La_version(), error = function(e) NA_character_),
    ext_blas = if ("BLAS" %in% names(ext)) unname(ext[["BLAS"]]) else NA_character_,
    blas_vendor = suite_blas_vendor(c(blas_path, la_lib, if ("BLAS" %in% names(ext)) ext[["BLAS"]])),
    blas_threads_parent = suite_blas_threads(),
    env_threads = list(OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", NA),
                       OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", NA),
                       MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", NA)),
    profmem = isTRUE(capabilities("profmem")),
    eigencore_version = suite_pkg_version("eigencore"),
    eigencore_lib = if (!is.null(ec_desc)) dirname(dirname(attr(ec_desc, "file"))) else NA_character_,
    eigencore_built = if (!is.null(ec_desc)) ec_desc$Built %||% NA_character_ else NA_character_,
    eigencore_remote_sha = if (!is.null(ec_desc)) ec_desc$RemoteSha %||% NA_character_ else NA_character_,
    git_sha = git$sha,
    # fingerprint of the suite's own code (run-suite.R + suite/*.R), so a
    # result set identifies the exact harness even when it was uncommitted
    suite_md5 = if (!is.null(repo_dir)) substr(suite_md5_string(paste(unname(tools::md5sum(
      c(file.path(repo_dir, "inst", "benchmarks", "run-suite.R"),
        sort(list.files(file.path(repo_dir, "inst", "benchmarks", "suite"), full.names = TRUE)))
    )), collapse = "")), 1L, 12L) else NA_character_,
    git_dirty_package_sources = git$dirty,
    packages = as.list(stats::setNames(vapply(pkgs, suite_pkg_version, character(1)), pkgs)),
    loadavg_start = suite_loadavg(),
    args = args
  )
}

suite_write_json <- function(x, path) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
    return(invisible(TRUE))
  }
  # Minimal fallback: flatten to key: value lines.
  flat <- unlist(x)
  writeLines(c("{", paste0("  \"", names(flat), "\": \"", gsub("\"", "'", flat), "\"",
                           c(rep(",", length(flat) - 1L), "")), "}"), path)
  invisible(FALSE)
}
