#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
sizes <- if (args$quick) 200L else c(200L, 1000L, 2000L)
iterations <- 5L
k <- 5L

results <- lapply(sizes, function(n) {
  A <- path_laplacian(n)
  rows <- benchmark_eigen_case(
    A,
    k = k,
    target = smallest(),
    methods = c(
      "eigencore",
      "eigencore_lobpcg_preconditioned",
      "eigencore_lobpcg_tridiagonal",
      "RSpectra",
      "PRIMME"
    ),
    iterations = iterations,
    tol = 1e-8,
    seed = 900L
  )
  gate <- evaluate_preconditioned_lobpcg_gate(rows, k = k)
  rows$case <- "path_laplacian"
  rows$n <- n
  rows$k <- k
  gate$case <- "path_laplacian"
  gate$n <- n
  gate$k <- k
  list(rows = rows, gate = gate)
})

rows <- do.call(rbind, lapply(results, `[[`, "rows"))
gates <- do.call(rbind, lapply(results, `[[`, "gate"))

cat("Preconditioned LOBPCG benchmark rows\n")
print(rows)
cat("\nPreconditioned LOBPCG gate\n")
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "lobpcg-preconditioned-rows"))
  message("saved gates: ", save_benchmark_result(gates, "lobpcg-preconditioned-gates"))
}

if (args$strict && !all(gates$passed)) {
  stop("Preconditioned LOBPCG benchmark failed release gate.", call. = FALSE)
}
