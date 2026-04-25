#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
sizes <- if (args$quick) 200L else c(200L, 1000L)
iterations <- 5L
k <- if (args$quick) 5L else 20L

results <- lapply(sizes, function(n) {
  A <- path_laplacian(n)
  rows <- benchmark_eigen_case(
    A,
    k = k,
    target = smallest(),
    methods = c("eigencore", "eigencore_block", "RSpectra", "PRIMME"),
    iterations = iterations,
    tol = 1e-8,
    seed = 900L
  )
  rows$case <- "path_laplacian"
  rows$n <- n
  rows$k <- k
  rows
})

rows <- do.call(rbind, results)

cat("Native block Hermitian prototype benchmark rows\n")
print(rows)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "block-hermitian-prototype-rows"))
}

