#!/usr/bin/env Rscript

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
}
source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (is.na(args$iterations)) {
  if (args$quick) 5L else 30L
} else {
  args$iterations
}
dimensions <- if (args$quick) c(32L, 90L) else c(32L, 64L, 90L, 128L)
ranks <- if (args$quick) c(5L, 16L) else c(5L, 8L, 16L)

result <- benchmark_tiny_gram_eigensolvers(
  dimensions = dimensions,
  ranks = ranks,
  iterations = iterations,
  seed = 810L
)
print(result)

winners <- result[result$winner, , drop = FALSE]
print(winners[order(winners$dimension, winners$rank), , drop = FALSE])

if (args$save) {
  path <- save_benchmark_result(result, "tiny-gram-eigensolvers")
  message("saved: ", path)
}
