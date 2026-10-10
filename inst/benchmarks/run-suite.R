#!/usr/bin/env Rscript
#
# run-suite.R -- reproducible eigencore benchmark suite.
#
# Usage (from the repository root, against an *installed* eigencore):
#
#   R CMD INSTALL --preclean --no-docs -l .rlib .
#   Rscript inst/benchmarks/run-suite.R --lib=.rlib --profile=standard
#
# Options
#   --profile=quick|standard|scaling   case set (default quick); see suite/profiles.R
#   --families=a,b                     keep only these families or sweep groups
#                                      (e.g. sym_sparse,svd_sparse,scaling_k)
#   --cases=substr,...                 keep cases whose id contains any substring
#   --methods=eigencore,RSpectra,irlba,base,PRIMME
#                                      methods to run (missing packages are skipped)
#   --threads=1,4                      thread settings; each runs in a fresh worker
#                                      process with OPENBLAS/OMP/MKL_NUM_THREADS and
#                                      options(eigencore.threads) set to the value
#   --reps=N                           timed repetitions after one warm-up call
#   --budget=SECONDS                   per method/case time budget for repetitions
#   --tol=1e-8                         requested tolerance passed to every method
#   --out=DIR                          results root (default inst/benchmarks/results)
#   --lib=DIR                          library to load eigencore from (prepended)
#   --cache=DIR                        reference cache (default: R_user_dir cache)
#   --no-cache                         recompute references
#   --suitesparse                      add SuiteSparse Matrix Collection cases
#                                      (downloaded to <cache>/suitesparse; skipped offline)
#   --inprocess                        run thread settings in this process (debugging)
#   --list                             print the selected cases and exit
#
# Output: <out>/<YYYYMMDD>-<machine-id>-<profile>/
#   results.csv       one row per case x method x threads (wall times and
#                     process CPU time per repetition: cpu_median, cpu_min)
#   cases.csv         one row per case: generator, size, reference source, ||A||_2
#   environment.json  machine, R, BLAS/LAPACK, package versions, git SHA, load
#   environment.rds   the same as an R list (plus per-worker thread settings)
#
# Summaries: source("inst/benchmarks/report.R"); see that file.

suite_script_path <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) return(normalizePath(sub("^--file=", "", f[[1L]])))
  normalizePath("inst/benchmarks/run-suite.R", mustWork = FALSE)
}
SCRIPT <- suite_script_path()
SCRIPT_DIR <- dirname(SCRIPT)
for (f in c("utils.R", "generators.R", "profiles.R", "methods.R", "accuracy.R",
            "environment.R", "suitesparse.R")) {
  source(file.path(SCRIPT_DIR, "suite", f))
}

args <- suite_parse_args()
if (suite_flag(args$help)) {
  lines <- readLines(SCRIPT)
  cat(sub("^# ?", "", lines[2:which(lines == "")[1L]]), sep = "\n")
  quit(status = 0L)
}
if (!is.null(args$lib)) .libPaths(c(normalizePath(args$lib), .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  if (!requireNamespace("eigencore", quietly = TRUE)) stop("eigencore is not installed (use --lib=)")
})

profile <- args$profile %||% "quick"
defaults <- suite_profile_defaults(profile)
reps <- as.integer(args$reps %||% defaults$reps)
budget <- as.numeric(args$budget %||% defaults$budget)
tol <- as.numeric(args$tol %||% 1e-8)
threads <- as.integer(suite_csv_arg(args$threads %||% defaults$threads))
methods <- suite_csv_arg(args$methods) %||% SUITE_METHODS
families <- suite_csv_arg(args$families)
case_ids <- suite_csv_arg(args$cases)
cache_dir <- if (suite_flag(args$no_cache)) NULL else
  (args$cache %||% file.path(tools::R_user_dir("eigencore", "cache"), "bench-suite"))
repo_dir <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(repo_dir, "DESCRIPTION"))) repo_dir <- NULL

suite_select_cases <- function() {
  cases <- suite_profile_cases(profile)
  if (suite_flag(args$suitesparse)) {
    cases <- c(cases, suite_suitesparse_cases(file.path(cache_dir %||% tempdir(), "suitesparse")))
  }
  suite_filter_cases(cases, families, case_ids)
}

# ---------------------------------------------------------------------------
# Timing harness
# ---------------------------------------------------------------------------

suite_time_method <- function(ad, reps, budget) {
  warns <- character()
  muffle <- function(w) {
    warns <<- unique(c(warns, conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
  invisible(gc(FALSE))
  rss_ok <- suite_rss_reset()
  rss_before <- suite_rss_field("VmRSS")
  g0 <- gc(reset = TRUE)
  load_start <- suite_loadavg()[1L]
  t0 <- suite_now()
  raw <- tryCatch(withCallingHandlers(ad$run(), warning = muffle),
                  error = function(e) structure(list(message = conditionMessage(e)), class = "suite_error"))
  t_first <- suite_now() - t0
  rss_peak <- if (rss_ok) suite_rss_field("VmHWM") - rss_before else NA_real_
  g1 <- gc()
  heap_peak <- sum(g1[, ncol(g1)]) - sum(g0[, 2L])
  if (inherits(raw, "suite_error")) {
    return(list(error = raw$message, t_first = t_first, warnings = warns, load_start = load_start))
  }
  n_reps <- if (t_first >= budget) 0L else as.integer(min(reps, max(1, floor(budget / max(t_first, 1e-6)))))
  times <- numeric(0)
  cpu <- numeric(0)
  for (i in seq_len(n_reps)) {
    invisible(gc(FALSE))
    c0 <- suite_cpu()
    t0 <- suite_now()
    withCallingHandlers(ad$run(), warning = muffle)
    times <- c(times, suite_now() - t0)
    cpu <- c(cpu, suite_cpu() - c0)
  }
  timed <- if (length(times)) times else t_first
  list(raw = raw, error = NULL, t_first = t_first, times = times,
       time_median = stats::median(timed), time_min = min(timed), time_max = max(timed),
       cpu_median = if (length(cpu)) stats::median(cpu) else NA_real_,
       cpu_min = if (length(cpu)) min(cpu) else NA_real_,
       reps_timed = length(times), rss_peak_mb = rss_peak, r_heap_peak_mb = heap_peak,
       warnings = warns, load_start = load_start, load_end = suite_loadavg()[1L])
}

suite_build_problem <- function(case) {
  prob <- case$gen()
  if (isTRUE(case$center)) prob$mu <- as.numeric(Matrix::colMeans(prob$A))
  prob
}

suite_case_info <- function(case, prob, ref) {
  A <- prob$A
  data.frame(
    case_id = case$id, family = case$family, group = case$group, task = case$task,
    target = case$target, k = case$k, sigma = case$sigma, center = case$center,
    vectors = case$vectors, seed = case$seed, sweep = case$sweep, sweep_value = case$sweep_value,
    m = nrow(A), n = ncol(A),
    nnz = if (inherits(A, "sparseMatrix")) length(A@x) else as.numeric(nrow(A)) * ncol(A),
    storage = if (inherits(A, "sparseMatrix")) "sparse" else "dense",
    generalized = !is.null(prob$B),
    description = prob$description %||% NA_character_,
    fingerprint = ref$fingerprint %||% suite_fingerprint(A),
    reference_source = ref$source %||% NA_character_,
    reference_crosscheck = ref$crosscheck %||% NA_real_,
    reference_backward_error = ref$reference_backward_error %||% NA_real_,
    reference_ambiguous = isTRUE(ref$ambiguous),
    reference_seconds = ref$seconds %||% NA_real_,
    norm2 = ref$norm2 %||% NA_real_, norm2_source = ref$norm2_source %||% NA_character_,
    normB = ref$normB %||% NA_real_,
    stringsAsFactors = FALSE)
}

suite_row <- function(case, method, nthreads, ad = NULL, timing = NULL, std = NULL, acc = NULL,
                      status = "ok", message = NA_character_) {
  g <- function(x, field, default = NA) {
    if (is.null(x)) return(default)
    v <- x[[field]]
    if (is.null(v) || length(v) != 1L) default else v
  }
  data.frame(
    case_id = case$id, family = case$family, group = case$group, task = case$task,
    target = case$target, k = case$k, sweep = case$sweep, sweep_value = case$sweep_value,
    method = method, threads = nthreads, status = status, message = message,
    call = g(ad, "call", NA_character_),
    method_label = g(std, "label", NA_character_),
    reps_timed = g(timing, "reps_timed", NA_integer_),
    time_median = g(timing, "time_median", NA_real_),
    time_min = g(timing, "time_min", NA_real_),
    time_max = g(timing, "time_max", NA_real_),
    time_first = g(timing, "t_first", NA_real_),
    cpu_median = g(timing, "cpu_median", NA_real_),
    cpu_min = g(timing, "cpu_min", NA_real_),
    times = if (is.null(timing) || !length(timing$times)) NA_character_ else
      paste(sprintf("%.4f", timing$times), collapse = ";"),
    matvecs = as.numeric(g(std, "matvecs", NA_real_)),
    matvec_source = g(std, "matvec_source", NA_character_),
    ec_operator_columns = as.numeric(g(std, "ec_operator_columns", NA_real_)),
    ec_adjoint_columns = as.numeric(g(std, "ec_adjoint_columns", NA_real_)),
    ec_cert_columns = as.numeric(g(std, "ec_cert_columns", NA_real_)),
    iterations = as.numeric(g(std, "iterations", NA_real_)),
    returned_k = g(acc, "returned_k", NA_integer_),
    target_ok = g(acc, "target_ok", NA),
    target_tau = g(acc, "target_tau", NA_real_),
    value_err = g(acc, "value_err", NA_real_),
    value_rel_err = g(acc, "value_rel_err", NA_real_),
    max_residual = g(acc, "max_residual", NA_real_),
    max_backward_error = g(acc, "max_backward_error", NA_real_),
    max_left_residual = g(acc, "max_left_residual", NA_real_),
    orthogonality_loss = g(acc, "orthogonality_loss", NA_real_),
    eigencore_certified = g(std, "certified", NA),
    eigencore_own_backward_error = as.numeric(g(std, "own_backward_error", NA_real_)),
    rss_peak_mb = g(timing, "rss_peak_mb", NA_real_),
    r_heap_peak_mb = g(timing, "r_heap_peak_mb", NA_real_),
    load1_start = g(timing, "load_start", NA_real_),
    load1_end = g(timing, "load_end", NA_real_),
    warnings = if (is.null(timing) || !length(timing$warnings)) NA_character_ else
      substr(paste(timing$warnings, collapse = " | "), 1L, 300L),
    stringsAsFactors = FALSE)
}

suite_set_threads <- function(nthreads) {
  options(eigencore.threads = nthreads)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(nthreads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(nthreads), silent = TRUE)
  }
}

suite_run_worker <- function(nthreads) {
  suite_set_threads(nthreads)
  avail <- methods[vapply(methods, suite_method_available, logical(1))]
  missing <- setdiff(methods, avail)
  if (length(missing)) suite_log("methods not installed, skipped: ", paste(missing, collapse = ", "))
  cases <- suite_select_cases()
  rows <- list()
  infos <- list()
  for (case in cases) {
    suite_log(sprintf("[threads=%d] case %s", nthreads, case$id))
    prob <- tryCatch(suite_build_problem(case), error = function(e) e)
    if (inherits(prob, "error")) {
      suite_log("  generator failed: ", conditionMessage(prob))
      next
    }
    ref <- tryCatch(suite_reference(case, prob, cache_dir), error = function(e) {
      suite_log("  reference failed: ", conditionMessage(e))
      list(values = NULL, norm2 = NA_real_, normB = 1, source = paste("failed:", conditionMessage(e)))
    })
    if (!isTRUE(ref$cached) && !is.null(ref$seconds)) {
      suite_log(sprintf("  reference: %s (%.1f s); ||A||_2 = %.6g [%s]", ref$source, ref$seconds,
                        ref$norm2, ref$norm2_source))
    }
    infos[[case$id]] <- suite_case_info(case, prob, ref)
    for (method in avail) {
      ad <- tryCatch(suite_adapter(method, case, prob, tol), error = function(e) NULL)
      if (is.null(ad)) next
      if (!is.null(ad$skip)) {
        rows[[length(rows) + 1L]] <- suite_row(case, method, nthreads, ad, status = "skipped", message = ad$skip)
        next
      }
      timing <- suite_time_method(ad, reps, budget)
      if (!is.null(timing$error)) {
        rows[[length(rows) + 1L]] <- suite_row(case, method, nthreads, ad, timing,
                                               status = "error", message = timing$error)
        suite_log(sprintf("  %-9s ERROR %s", method, timing$error))
        next
      }
      std <- tryCatch(ad$extract(timing$raw), error = function(e) e)
      if (inherits(std, "error")) {
        rows[[length(rows) + 1L]] <- suite_row(case, method, nthreads, ad, timing,
                                               status = "error", message = paste("extract:", conditionMessage(std)))
        next
      }
      acc <- tryCatch(suite_accuracy(case, prob, std, ref), error = function(e) {
        suite_log("  accuracy check failed: ", conditionMessage(e))
        NULL
      })
      timing$raw <- NULL
      row <- suite_row(case, method, nthreads, ad, timing, std, acc)
      rows[[length(rows) + 1L]] <- row
      suite_log(sprintf("  %-9s median %-9s (min %-9s, cpu %-9s, %d reps) matvecs %-6s bwd %.1e value_err %.1e target %s%s",
                        method, suite_fmt_time(row$time_median), suite_fmt_time(row$time_min),
                        suite_fmt_time(row$cpu_median),
                        row$reps_timed, format(row$matvecs), row$max_backward_error, row$value_err,
                        row$target_ok, if (is.na(row$eigencore_certified)) "" else
                          paste0(" certified ", row$eigencore_certified)))
    }
    rm(prob)
    invisible(gc(FALSE))
  }
  list(rows = if (length(rows)) do.call(rbind, rows) else NULL,
       cases = if (length(infos)) do.call(rbind, unname(infos)) else NULL,
       worker = list(threads = nthreads, eigencore_threads = getOption("eigencore.threads"),
                     blas_threads = suite_blas_threads(),
                     env = list(OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", NA),
                                OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", NA)),
                     methods_run = avail, methods_missing = missing,
                     loadavg_end = suite_loadavg()))
}

# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

if (suite_flag(args$list)) {
  for (cs in suite_select_cases()) cat(sprintf("%-45s %-14s %-6s %-4s k=%d\n", cs$id, cs$family, cs$task, cs$target, cs$k))
  quit(status = 0L)
}

if (suite_flag(args$worker)) {
  res <- suite_run_worker(as.integer(args$threads))
  saveRDS(res, args$part)
  quit(status = 0L)
}

t_start <- suite_now()
out_root <- args$out %||% file.path(SCRIPT_DIR, "results")
env <- suite_environment(repo_dir, profile, args)
run_id <- sprintf("%s-%s-%s", format(Sys.time(), "%Y%m%d"), env$machine_id, profile)
run_dir <- file.path(out_root, run_id)
if (dir.exists(run_dir)) {
  run_id <- paste0(run_id, "-", format(Sys.time(), "%H%M%S"))
  run_dir <- file.path(out_root, run_id)
}
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
suite_log("eigencore benchmark suite: profile ", profile, ", threads ", paste(threads, collapse = ","),
          ", reps ", reps, ", output ", run_dir)
suite_log("eigencore ", env$eigencore_version, " (git ", env$git_sha, "), ", env$blas_vendor, ", ",
          env$cpu_model, ", load ", paste(env$loadavg_start, collapse = " "))

parts <- list()
for (nt in threads) {
  if (suite_flag(args$inprocess)) {
    parts[[length(parts) + 1L]] <- suite_run_worker(nt)
    next
  }
  part <- file.path(run_dir, sprintf(".part-threads%d.rds", nt))
  pass <- grep("^--(profile|families|cases|methods|reps|budget|tol|lib|cache|no-cache|suitesparse)(=|$)",
               commandArgs(trailingOnly = TRUE), value = TRUE)
  env_vars <- sprintf(c("OPENBLAS_NUM_THREADS=%d", "OMP_NUM_THREADS=%d", "MKL_NUM_THREADS=%d",
                        "VECLIB_MAXIMUM_THREADS=%d"), nt)
  if (!is.null(args$lib)) env_vars <- c(env_vars, paste0("R_LIBS=", paste(.libPaths(), collapse = ":")))
  rscript <- file.path(R.home("bin"), "Rscript")
  status <- system2(rscript, c(shQuote(SCRIPT), pass, "--worker", sprintf("--threads=%d", nt),
                               paste0("--part=", shQuote(part))), env = env_vars)
  if (!identical(as.integer(status), 0L) || !file.exists(part)) {
    suite_log("worker for threads=", nt, " failed (status ", status, ")")
    next
  }
  parts[[length(parts) + 1L]] <- readRDS(part)
  unlink(part)
}

rows <- do.call(rbind, lapply(parts, `[[`, "rows"))
cases <- do.call(rbind, lapply(parts, `[[`, "cases"))
if (!is.null(cases)) cases <- cases[!duplicated(cases$case_id), , drop = FALSE]
if (is.null(rows)) stop("no results were produced")
rows <- cbind(run_id = run_id, profile = profile, rows, stringsAsFactors = FALSE)

env$workers <- lapply(parts, `[[`, "worker")
env$threads <- threads
env$reps <- reps
env$budget_seconds <- budget
env$tol <- tol
env$loadavg_end <- suite_loadavg()
env$wall_seconds <- round(suite_now() - t_start, 1)
env$run_id <- run_id

utils::write.csv(rows, file.path(run_dir, "results.csv"), row.names = FALSE)
utils::write.csv(cases, file.path(run_dir, "cases.csv"), row.names = FALSE)
saveRDS(env, file.path(run_dir, "environment.rds"))
suite_write_json(env, file.path(run_dir, "environment.json"))

suite_log(sprintf("done in %.0f s; wrote %d rows to %s", env$wall_seconds, nrow(rows), run_dir))
report <- file.path(SCRIPT_DIR, "report.R")
if (file.exists(report)) {
  source(report)
  res <- bench_load(run_dir)
  print(bench_time_table(res), row.names = FALSE)
}
