#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (args$quick) 1L else 5L

rows <- benchmark_g1_candidate_baseline(
  quick = args$quick,
  iterations = iterations
)

cat("G1 candidate pre-promotion baseline rows\n")
print(rows)

if (args$save) {
  path <- "inst/benchmarks/baselines/g1_candidate_pre.csv"
  utils::write.csv(rows, path, row.names = FALSE)
  message("saved baseline: ", path)
}

if (args$strict) {
  block <- rows[rows$method == "eigencore_block_candidate", , drop = FALSE]
  required_cases <- c("path_laplacian", "dense_hermitian", "clustered", "ill_conditioned_diag")
  if (!all(required_cases %in% block$case)) {
    stop("G1 baseline is missing a required block-candidate case.", call. = FALSE)
  }
  if (any(block$certificate_type == "method_error")) {
    stop("G1 baseline block candidate has method-error rows.", call. = FALSE)
  }
  if (!all(block$certificate_passed)) {
    stop("G1 baseline block candidate has uncertified rows.", call. = FALSE)
  }
}
