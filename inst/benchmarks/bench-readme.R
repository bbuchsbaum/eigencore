#!/usr/bin/env Rscript
#
# bench-readme.R -- retired; kept as a thin wrapper around the benchmark suite.
#
# The README "Benchmarks" table is now generated from stored results of
# inst/benchmarks/run-suite.R (profile "standard"), not from this script.
# This wrapper runs that profile and prints the README table:
#
#   R CMD INSTALL --preclean --no-docs -l .rlib .
#   Rscript inst/benchmarks/bench-readme.R --lib=.rlib [any run-suite.R option]
#
# The suite never forces R memory profiling; whether it is available here:
#   capabilities("profmem")

self <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(self)) dirname(normalizePath(sub("^--file=", "", self[[1L]]))) else "inst/benchmarks"
args <- commandArgs(trailingOnly = TRUE)
if (!any(startsWith(args, "--profile="))) args <- c("--profile=standard", args)
message("bench-readme.R now delegates to run-suite.R; R memory profiling available: ",
        isTRUE(capabilities("profmem")))
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(file.path(here, "run-suite.R")), args))
if (!identical(as.integer(status), 0L)) quit(status = 1L)

source(file.path(here, "report.R"))
out <- sub("^--out=", "", grep("^--out=", args, value = TRUE))
root <- if (length(out)) out[[1L]] else file.path(here, "results")
prof <- sub("^--profile=", "", grep("^--profile=", args, value = TRUE)[[1L]])
res <- bench_load(bench_find_runs(root, profile = prof)[1L])
cat("\n", bench_readme_markdown(res), sep = "\n")
cat("\n<sub>", bench_env_line(res), "</sub>\n", sep = "")
