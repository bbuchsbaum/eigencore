# report.R -- load stored benchmark-suite results and summarise them.
#
# Base R only (jsonlite/ggplot2 are not required). Typical use:
#
#   source("inst/benchmarks/report.R")
#   res <- bench_load(bench_find_runs(profile = "standard"))
#   bench_env_table(res)
#   bench_time_table(res, threads = 1)
#   bench_ratio_table(res)
#   bench_matvec_table(res, threads = 1)
#   bench_accuracy_table(res, threads = 1)
#   bench_plot_scaling(bench_load(bench_find_runs(profile = "scaling")), "scaling_n")
#   cat(bench_readme_markdown(res), sep = "\n")
#
# or from the shell:  Rscript inst/benchmarks/report.R <run-dir> [<run-dir> ...]

bench_methods_order <- c("eigencore", "RSpectra", "irlba", "PRIMME", "base")

bench_results_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    cand <- file.path(d, "inst", "benchmarks", "results")
    if (file.exists(file.path(d, "DESCRIPTION")) && dir.exists(cand)) return(cand)
    d <- dirname(d)
  }
  NULL
}

# Run directories under `root`, newest first; optionally one profile and only
# the latest run per machine.
bench_find_runs <- function(root = bench_results_root(), profile = NULL, latest_per_machine = TRUE) {
  if (is.null(root) || !dir.exists(root)) return(character())
  dirs <- list.dirs(root, recursive = FALSE)
  dirs <- dirs[file.exists(file.path(dirs, "results.csv"))]
  if (!length(dirs)) return(character())
  b <- basename(dirs)
  parts <- strsplit(b, "-", fixed = TRUE)
  prof <- vapply(parts, function(p) if (length(p) >= 3L) p[[3L]] else NA_character_, character(1))
  mach <- vapply(parts, function(p) if (length(p) >= 2L) p[[2L]] else NA_character_, character(1))
  keep <- if (is.null(profile)) rep(TRUE, length(dirs)) else prof %in% profile
  dirs <- dirs[keep]; b <- b[keep]; mach <- mach[keep]; prof <- prof[keep]
  o <- order(b, decreasing = TRUE)
  dirs <- dirs[o]; mach <- mach[o]; prof <- prof[o]
  if (latest_per_machine) dirs <- dirs[!duplicated(paste(mach, prof))]
  dirs
}

bench_read_env <- function(dir) {
  rds <- file.path(dir, "environment.rds")
  if (file.exists(rds)) return(readRDS(rds))
  js <- file.path(dir, "environment.json")
  if (file.exists(js) && requireNamespace("jsonlite", quietly = TRUE)) return(jsonlite::read_json(js))
  list()
}

bench_load <- function(dirs) {
  dirs <- dirs[file.exists(file.path(dirs, "results.csv"))]
  if (!length(dirs)) return(NULL)
  rd <- function(f) utils::read.csv(f, stringsAsFactors = FALSE)
  results <- do.call(rbind, lapply(dirs, function(d) rd(file.path(d, "results.csv"))))
  cases <- do.call(rbind, lapply(dirs, function(d) {
    x <- rd(file.path(d, "cases.csv"))
    cbind(run_id = basename(d), x, stringsAsFactors = FALSE)
  }))
  env <- stats::setNames(lapply(dirs, bench_read_env), basename(dirs))
  results$target_ok <- as.logical(results$target_ok)
  results$eigencore_certified <- as.logical(results$eigencore_certified)
  structure(list(results = results, cases = cases, env = env, dirs = dirs), class = "bench_results")
}

bench_fmt_time <- function(x) {
  out <- ifelse(!is.finite(x), "-",
         ifelse(x < 0.01, sprintf("%.1f ms", x * 1000),
         ifelse(x < 1, sprintf("%.0f ms", x * 1000),
         ifelse(x < 10, sprintf("%.2f s", x), sprintf("%.1f s", x)))))
  out
}

bench_fmt_sci <- function(x) ifelse(is.finite(x), formatC(x, format = "e", digits = 1), "-")

bench_fmt_ratio <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "-")

# Quality flag for a result row: "" (accurate, right target set),
# "†" (wrong / incomplete target set) or "‡" (backward error > 1e-6).
bench_flag <- function(r) {
  ifelse(r$status != "ok", "",
  ifelse(!is.na(r$target_ok) & !r$target_ok, "†",
  ifelse(is.finite(r$max_backward_error) & r$max_backward_error > 1e-6, "‡", "")))
}

bench_rows <- function(res, threads = NULL, run_id = NULL) {
  r <- res$results
  if (!is.null(threads)) r <- r[r$threads %in% threads, , drop = FALSE]
  if (!is.null(run_id)) r <- r[r$run_id %in% run_id, , drop = FALSE]
  r
}

bench_case_label <- function(res, ids) {
  cs <- res$cases[!duplicated(res$cases$case_id), , drop = FALSE]
  lab <- cs$case_id
  names(lab) <- cs$case_id
  unname(lab[ids])
}

bench_wide <- function(r, value_fun, methods = NULL) {
  methods <- methods %||% intersect(bench_methods_order, unique(r$method))
  keys <- unique(r[, c("run_id", "case_id", "threads")])
  out <- keys
  for (m in methods) {
    out[[m]] <- vapply(seq_len(nrow(keys)), function(i) {
      x <- r[r$run_id == keys$run_id[i] & r$case_id == keys$case_id[i] &
               r$threads == keys$threads[i] & r$method == m, , drop = FALSE]
      if (!nrow(x)) return(NA_character_)
      value_fun(x[1L, , drop = FALSE])
    }, character(1))
  }
  out
}

if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# Median wall time per method (status/quality flags appended).
bench_time_table <- function(res, threads = NULL, methods = NULL, run_id = NULL) {
  r <- bench_rows(res, threads, run_id)
  bench_wide(r, function(x) {
    if (x$status == "skipped") return("n/a")
    if (x$status != "ok") return("error")
    paste0(bench_fmt_time(x$time_median), bench_flag(x))
  }, methods)
}

# Time of each method divided by RSpectra's (< 1: faster than RSpectra).
bench_ratio_table <- function(res, threads = NULL, methods = c("eigencore", "irlba", "PRIMME", "base"),
                              run_id = NULL) {
  r <- bench_rows(res, threads, run_id)
  rs <- r[r$method == "RSpectra" & r$status == "ok", c("run_id", "case_id", "threads", "time_median")]
  names(rs)[4] <- "rs_time"
  r <- merge(r, rs, by = c("run_id", "case_id", "threads"))
  methods <- intersect(methods, unique(r$method))
  bench_wide(r, function(x) if (x$status != "ok") "-" else bench_fmt_ratio(x$time_median / x$rs_time), methods)
}

bench_matvec_table <- function(res, threads = NULL, methods = NULL, run_id = NULL) {
  r <- bench_rows(res, threads, run_id)
  r <- r[r$method != "base", , drop = FALSE]
  bench_wide(r, function(x) if (x$status != "ok" || !is.finite(x$matvecs)) "-" else format(x$matvecs), methods)
}

bench_accuracy_table <- function(res, threads = NULL, methods = NULL, run_id = NULL) {
  r <- bench_rows(res, threads, run_id)
  bench_wide(r, function(x) {
    if (x$status != "ok") return(if (x$status == "skipped") "n/a" else "error")
    tgt <- if (is.na(x$target_ok)) "?" else if (x$target_ok) "ok" else "MISS"
    cert <- if (is.na(x$eigencore_certified)) "" else if (x$eigencore_certified) ", cert" else ", cert FAILED"
    sprintf("%s / %s / %s%s", bench_fmt_sci(x$max_backward_error), bench_fmt_sci(x$value_err), tgt, cert)
  }, methods)
}

# Long-format accuracy rows (numeric), handy for custom summaries.
bench_accuracy_long <- function(res, threads = NULL) {
  r <- bench_rows(res, threads)
  r[r$status == "ok", c("run_id", "case_id", "threads", "method", "returned_k", "target_ok",
                        "max_backward_error", "value_err", "orthogonality_loss",
                        "eigencore_certified", "eigencore_own_backward_error")]
}

bench_env_table <- function(res) {
  rows <- lapply(names(res$env), function(id) {
    e <- res$env[[id]]
    p <- e$packages %||% list()
    w <- e$workers %||% list()
    bl <- vapply(w, function(x) sprintf("%s:%s", x$threads, x$blas_threads$blas %||% NA), character(1))
    data.frame(
      run = id,
      date = e$timestamp_utc %||% NA,
      machine = e$machine_id %||% NA,
      cpu = sprintf("%s (%s logical / %s physical cores, %s GB)", e$cpu_model %||% "?",
                    e$cores_logical %||% "?", e$cores_physical %||% "?", e$mem_total_gb %||% "?"),
      R = sub("R version ", "", e$r_version %||% NA),
      platform = e$platform %||% NA,
      BLAS = sprintf("%s (%s)", e$blas_vendor %||% "?", basename(e$la_library %||% e$blas_path %||% "?")),
      eigencore = sprintf("%s @ %s%s", p$eigencore %||% "?", e$git_sha %||% "?",
                          if (isTRUE(e$git_dirty_package_sources)) " (dirty)" else ""),
      competitors = paste(na.omit(c(
        if (!is.na(p$RSpectra %||% NA)) paste("RSpectra", p$RSpectra),
        if (!is.na(p$irlba %||% NA)) paste("irlba", p$irlba),
        if (!is.na(p$PRIMME %||% NA)) paste("PRIMME", p$PRIMME))), collapse = ", "),
      threads = paste(e$threads %||% NA, collapse = ","),
      blas_threads = paste(bl, collapse = " "),
      reps = e$reps %||% NA,
      load_start_end = sprintf("%.2f / %.2f", (e$loadavg_start %||% NA)[1], (e$loadavg_end %||% NA)[1]),
      wall = sprintf("%.0f s", e$wall_seconds %||% NA),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

bench_env_line <- function(res, id = names(res$env)[1L]) {
  e <- res$env[[id]]
  p <- e$packages %||% list()
  sprintf("%s; %s, %s logical cores; R %s; %s; eigencore %s (%s); RSpectra %s, irlba %s; median of %s reps; load average %.1f at start; run `%s`.",
          sub("T.*", "", e$timestamp_utc %||% ""), e$cpu_model %||% "?", e$cores_logical %||% "?",
          sub("R version ([^ ]+).*", "\\1", e$r_version %||% "?"), e$blas_vendor %||% "?",
          p$eigencore %||% "?", e$git_sha %||% "?", p$RSpectra %||% "-", p$irlba %||% "-",
          e$reps %||% "?", (e$loadavg_start %||% NA)[1], id)
}

# Scaling curves: time (and optionally matvecs) against the sweep variable.
bench_plot_scaling <- function(res, group, threads = 1, methods = c("eigencore", "RSpectra", "irlba", "PRIMME"),
                               metric = c("time_median", "matvecs"), main = NULL) {
  metric <- match.arg(metric)
  r <- res$results
  r <- r[r$group == group & r$threads %in% threads & r$status == "ok" & r$method %in% methods, , drop = FALSE]
  if (!nrow(r)) {
    plot.new()
    title(main = paste(group, "- no data"))
    return(invisible(NULL))
  }
  r <- r[is.finite(r[[metric]]), , drop = FALSE]
  cols <- c(eigencore = "#2166ac", RSpectra = "#b2182b", irlba = "#4d9221", PRIMME = "#762a83", base = "grey40")
  pchs <- c(eigencore = 19, RSpectra = 17, irlba = 15, PRIMME = 18, base = 4)
  ms <- intersect(methods, unique(r$method))
  x <- r$sweep_value
  y <- r[[metric]]
  op <- graphics::par(mar = c(4.2, 4.4, 2.2, 1))
  on.exit(graphics::par(op))
  graphics::plot(range(x), range(y), type = "n", log = "xy", bty = "n",
                 xlab = unique(r$sweep)[1L],
                 ylab = if (metric == "time_median") "median wall time (s)" else "operator applications",
                 main = main %||% sprintf("%s (threads = %s)", group, paste(threads, collapse = ",")))
  for (m in ms) {
    for (t in threads) {
      s <- r[r$method == m & r$threads == t, , drop = FALSE]
      s <- s[order(s$sweep_value), , drop = FALSE]
      bad <- !is.na(s$target_ok) & !s$target_ok
      graphics::lines(s$sweep_value, s[[metric]], col = cols[[m]], lty = if (t == min(threads)) 1 else 2)
      graphics::points(s$sweep_value, s[[metric]], col = cols[[m]], pch = ifelse(bad, 1, pchs[[m]]))
    }
  }
  leg <- ms
  if (length(threads) > 1L) leg <- c(leg, paste("threads", threads))
  graphics::legend("topleft", legend = leg, bty = "n", cex = 0.85,
                   col = c(cols[ms], rep("grey30", length(leg) - length(ms))),
                   pch = c(pchs[ms], rep(NA, length(leg) - length(ms))),
                   lty = c(rep(1, length(ms)), if (length(threads) > 1L) seq_along(threads)))
  invisible(r)
}

# Fitted log-log slope of time vs the sweep variable per method/threads.
bench_scaling_slopes <- function(res, group) {
  r <- res$results
  r <- r[r$group == group & r$status == "ok" & is.finite(r$time_median), , drop = FALSE]
  keys <- unique(r[, c("method", "threads")])
  keys$slope <- vapply(seq_len(nrow(keys)), function(i) {
    s <- r[r$method == keys$method[i] & r$threads == keys$threads[i], , drop = FALSE]
    if (nrow(s) < 3L) return(NA_real_)
    unname(stats::coef(stats::lm(log(time_median) ~ log(sweep_value), data = s))[2L])
  }, numeric(1))
  keys
}

# Compact README table: one row per core case, eigencore at 1 and 4 threads
# against single-threaded competitors.
bench_readme_markdown <- function(res, cases = NULL, labels = NULL) {
  r <- res$results
  ids <- cases %||% unique(r$case_id)
  thr <- sort(unique(r$threads))
  get <- function(id, m, t) {
    x <- r[r$case_id == id & r$method == m & r$threads == t, , drop = FALSE]
    if (!nrow(x)) return("-")
    x <- x[1L, ]
    if (x$status == "skipped") return("n/a")
    if (x$status != "ok") return("error")
    paste0(bench_fmt_time(x$time_median), bench_flag(x))
  }
  t1 <- min(thr)
  tN <- max(thr)
  hdr <- c("Problem", sprintf("eigencore (%d thr)", t1), if (tN > t1) sprintf("eigencore (%d thr)", tN),
           "RSpectra", "irlba", "base R", "eigencore / RSpectra", "certified")
  lines <- c(paste0("| ", paste(hdr, collapse = " | "), " |"),
             paste0("|", paste(c("---", rep("---:", length(hdr) - 2L), ":---:"), collapse = "|"), "|"))
  for (id in ids) {
    ec <- r[r$case_id == id & r$method == "eigencore" & r$threads == t1, , drop = FALSE]
    rs <- r[r$case_id == id & r$method == "RSpectra" & r$threads == t1, , drop = FALSE]
    ratio <- if (nrow(ec) && nrow(rs) && ec$status[1] == "ok" && rs$status[1] == "ok")
      bench_fmt_ratio(ec$time_median[1] / rs$time_median[1]) else "-"
    cert <- if (nrow(ec) && !is.na(ec$eigencore_certified[1])) (if (ec$eigencore_certified[1]) "yes" else "**no**") else "-"
    lab <- if (!is.null(labels) && !is.null(labels[[id]])) labels[[id]] else id
    cells <- c(lab, get(id, "eigencore", t1), if (tN > t1) get(id, "eigencore", tN),
               get(id, "RSpectra", t1), get(id, "irlba", t1), get(id, "base", t1), ratio, cert)
    lines <- c(lines, paste0("| ", paste(cells, collapse = " | "), " |"))
  }
  lines
}

bench_print_all <- function(res) {
  cat("\n== Environment ==\n")
  print(bench_env_table(res), row.names = FALSE)
  for (t in sort(unique(res$results$threads))) {
    cat(sprintf("\n== Median wall time, threads = %d († wrong target set, ‡ backward error > 1e-6) ==\n", t))
    print(bench_time_table(res, threads = t), row.names = FALSE)
    cat(sprintf("\n== Time relative to RSpectra, threads = %d ==\n", t))
    print(bench_ratio_table(res, threads = t), row.names = FALSE)
  }
  cat("\n== Operator applications (threads = min) ==\n")
  print(bench_matvec_table(res, threads = min(res$results$threads)), row.names = FALSE)
  cat("\n== Accuracy: backward error / value error / target set (threads = min) ==\n")
  print(bench_accuracy_table(res, threads = min(res$results$threads)), row.names = FALSE)
}

if (sys.nframe() == 0L) {
  dirs <- commandArgs(trailingOnly = TRUE)
  if (!length(dirs)) dirs <- bench_find_runs()
  res <- bench_load(dirs)
  if (is.null(res)) stop("no results found")
  bench_print_all(res)
}
