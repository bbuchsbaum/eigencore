#!/usr/bin/env Rscript
# Oracle sweep runner (see docs/test-assurance.md).
#
#   Rscript tools/oracle_sweep.R --level extended --workers 3 --out oracle-out
#   Rscript tools/oracle_sweep.R --ids 1:200,4711 --workers 2
#
# Runs every case of the level (or the given ids) in callr subprocesses
# (crash/hang isolation, one eigencore thread and one BLAS thread per
# worker), then writes <out>/records.rds, <out>/records.csv and
# <out>/summary.md (violations with reproducers, error and soft-finding
# counts, uncertified-rate tables). Exits with status 1 when any hard
# invariant was violated. eigencore must be installed in .libPaths().

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[[i + 1L]]
}
level <- opt("level", Sys.getenv("EIGENCORE_ORACLE_LEVEL", "extended"))
workers <- as.integer(opt("workers", "2"))
out_dir <- opt("out", "oracle-out")
ids_arg <- opt("ids", "")

helper <- normalizePath(file.path("tests", "testthat", "helper-oracle.R"))
suppressPackageStartupMessages(library(eigencore))
env <- new.env()
sys.source(helper, envir = env)
Sys.setenv(EIGENCORE_ORACLE_LEVEL = level)
ids <- if (nzchar(ids_arg)) {
  unlist(lapply(strsplit(ids_arg, ",")[[1]], function(x) eval(parse(text = x))))
} else {
  env$oracle_level_ids(level)
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("oracle sweep: level=%s cases=%d workers=%d eigencore=%s\n",
            level, length(ids), workers, utils::packageVersion("eigencore")))
started <- Sys.time()
last <- 0
recs <- env$oracle_run_ids(
  ids, subprocess = TRUE, workers = workers,
  progress = function(done, total) {
    if (done - last >= 200 || done == total) {
      cat(sprintf("  %d/%d (%.0f s)\n", done, total,
                  as.numeric(difftime(Sys.time(), started, units = "secs"))))
      last <<- done
    }
  })
df <- env$oracle_record_frame(recs)
saveRDS(df, file.path(out_dir, "records.rds"))
utils::write.csv(df, file.path(out_dir, "records.csv"), row.names = FALSE)

hard <- df[nzchar(df$hard), , drop = FALSE]
md <- c(
  sprintf("# eigencore oracle sweep (%s)", level), "",
  sprintf("- cases: %d, wall time: %.0f s, eigencore %s, %s",
          nrow(df), as.numeric(difftime(Sys.time(), started, units = "secs")),
          utils::packageVersion("eigencore"), R.version.string),
  sprintf("- status: %s", paste(names(table(df$status)), table(df$status),
                                sep = " = ", collapse = ", ")),
  sprintf("- certified: %d of %d solved", sum(df$certified %in% TRUE),
          sum(df$status == "ok")),
  sprintf("- hard violations: %d", nrow(hard)), "")
if (nrow(hard)) {
  md <- c(md, "## Hard violations", "")
  for (i in seq_len(nrow(hard))) {
    md <- c(md, sprintf("- case %d: %s", hard$id[i], hard$hard[i]),
            sprintf("  - `%s`", hard$describe[i]),
            sprintf("  - reproduce: `source('tests/testthat/helper-oracle.R'); library(eigencore); str(oracle_run_case(%d))`",
                    hard$id[i]))
  }
  md <- c(md, "")
}
soft <- unlist(strsplit(df$soft[nzchar(df$soft)], " | ", fixed = TRUE))
soft <- sub(":.*", "", soft)
if (length(soft)) {
  tab <- sort(table(soft), decreasing = TRUE)
  md <- c(md, "## Soft findings", "", "| finding | cases |", "|---|---|",
          sprintf("| %s | %d |", names(tab), as.integer(tab)), "")
}
errs <- df[df$status == "error_expected", , drop = FALSE]
if (nrow(errs)) {
  tab <- sort(table(substr(gsub("-?[0-9]+([.][0-9]+)?(e[-+]?[0-9]+)?", "#",
                                errs$error), 1, 90)),
              decreasing = TRUE)
  md <- c(md, "## Expected (unsupported-input) errors", "", "| message | cases |",
          "|---|---|", sprintf("| %s | %d |", gsub("\\|", "/", names(tab)),
                               as.integer(tab)), "")
}
table_md <- function(t) {
  if (!nrow(t)) return(character())
  c(paste0("| ", paste(names(t), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(t)), collapse = "|"), "|"),
    apply(t, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
md <- c(md, "## Uncertified rate by family / method", "",
        table_md(env$oracle_summary_table(df, c("family", "method"))), "",
        "## Uncertified rate by family / storage", "",
        table_md(env$oracle_summary_table(df, c("family", "storage"))), "",
        "## Uncertified rate by family / spectrum", "",
        table_md(env$oracle_summary_table(df, c("family", "spectrum"))), "",
        "## Uncertified rate by family / method / target", "",
        table_md(env$oracle_summary_table(df, c("family", "method", "target"))))
writeLines(md, file.path(out_dir, "summary.md"))
cat(sprintf("wrote %s (hard violations: %d)\n", file.path(out_dir, "summary.md"),
            nrow(hard)))
if (nrow(hard)) {
  for (i in seq_len(min(50L, nrow(hard)))) {
    cat(sprintf("case %d: %s\n  %s\n", hard$id[i], hard$hard[i], hard$describe[i]))
  }
  quit(status = 1L)
}
