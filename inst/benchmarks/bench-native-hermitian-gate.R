#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
sizes <- if (args$quick) 200L else c(1000L, 10000L)
iterations <- 5L
k <- if (args$quick) 5L else 20L

cases <- lapply(sizes, function(n) {
  list(
    case = "path_laplacian",
    A = path_laplacian(n),
    n = n,
    k = k,
    target = smallest(),
    seed = 700 + n
  )
})

if (args$include_dense) {
  cases <- c(cases, list(list(
    case = "dense_hermitian",
    A = dense_hermitian_with_spectrum(seq(if (args$quick) 80L else 200L, 1), seed = 1800L),
    n = if (args$quick) 80L else 200L,
    k = k,
    target = largest(),
    seed = 1800L
  )))
}

results <- lapply(cases, function(case) {
  subject <- if (args$block_candidate) "eigencore_block_candidate" else "eigencore"
  methods <- unique(c(subject, "eigencore", "RSpectra", "PRIMME"))
  out <- benchmark_native_hermitian_gate(
    case$A,
    k = case$k,
    target = case$target,
    iterations = iterations,
    seed = case$seed,
    methods = methods,
    subject = subject
  )
  out$rows$case <- case$case
  out$rows$n <- case$n
  out$rows$k <- case$k
  out$gate$case <- case$case
  out$gate$n <- case$n
  out$gate$k <- case$k
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

strict_passed <- if (args$quick) {
  all(gates$eigencore_certified)
} else {
  all(gates$passed)
}

if (args$strict && !strict_passed) {
  stop("Native Hermitian benchmark failed G1 release gate.", call. = FALSE)
}
