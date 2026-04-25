#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
sizes <- if (args$quick) 200L else c(1000L, 10000L)
iterations <- 5L
k <- if (args$quick) 5L else 20L

results <- lapply(sizes, function(n) {
  A <- path_laplacian(n)
  out <- benchmark_native_hermitian_gate(
    A,
    k = k,
    target = smallest(),
    iterations = iterations,
    seed = 700 + n
  )
  out$rows$case <- "path_laplacian"
  out$rows$n <- n
  out$rows$k <- k
  out$gate$case <- "path_laplacian"
  out$gate$n <- n
  out$gate$k <- k
  out
})

rows <- do.call(rbind, lapply(results, `[[`, "rows"))
gates <- do.call(rbind, lapply(results, `[[`, "gate"))

cat("Native Hermitian benchmark rows\n")
print(rows)
cat("\nNative Hermitian release gate\n")
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "native-hermitian-gate-rows"))
  message("saved gates: ", save_benchmark_result(gates, "native-hermitian-gate-summary"))
}
