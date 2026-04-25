#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
sizes <- if (args$quick) 200L else c(1000L, 10000L)
iterations <- if (args$quick) 1L else 5L
k <- if (args$quick) 5L else 20L

results <- lapply(sizes, function(n) {
  A <- path_laplacian(n)
  out <- benchmark_eigen_case(
    A,
    k = k,
    target = smallest(),
    methods = c("eigencore", "RSpectra", "PRIMME"),
    iterations = iterations,
    seed = 100 + n
  )
  out$case <- "path_laplacian"
  out$n <- n
  out$k <- k
  out
})

result <- do.call(rbind, results)
print(result)

gates <- lapply(split(result, result$n), function(rows) {
  gate <- evaluate_native_hermitian_gate(rows, k = unique(rows$k)[[1L]])
  gate$case <- unique(rows$case)
  gate$n <- unique(rows$n)
  gate
})
gates <- do.call(rbind, gates)
row.names(gates) <- NULL
print(gates)

if (args$save) {
  message("saved: ", save_benchmark_result(result, "hermitian-sparse"))
}

if (args$strict && !all(gates$passed)) {
  stop("Hermitian sparse benchmark failed PRD release gate.", call. = FALSE)
}
