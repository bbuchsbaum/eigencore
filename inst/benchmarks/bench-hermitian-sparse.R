#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
sizes <- if (args$quick) 200L else c(1000L, 10000L)
iterations <- if (args$quick) 1L else 5L
k <- if (args$quick) 5L else 20L

cases <- lapply(sizes, function(n) {
  list(
    case = "path_laplacian",
    A = path_laplacian(n),
    n = n,
    k = k,
    target = smallest(),
    seed = 100 + n
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
  out <- benchmark_eigen_case(
    case$A,
    k = case$k,
    target = case$target,
    methods = methods,
    iterations = iterations,
    seed = case$seed
  )
  out$case <- case$case
  out$n <- case$n
  out$k <- case$k
  out
})

result <- do.call(rbind, results)
print(result)

gates <- lapply(split(result, interaction(result$case, result$n, drop = TRUE)), function(rows) {
  subject <- if (args$block_candidate) "eigencore_block_candidate" else "eigencore"
  gate <- evaluate_native_hermitian_gate(
    rows,
    k = unique(rows$k)[[1L]],
    subject = subject
  )
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

strict_passed <- if (args$quick) {
  all(gates$eigencore_certified)
} else {
  all(gates$passed)
}

if (args$strict && !strict_passed) {
  stop("Hermitian sparse benchmark failed PRD release gate.", call. = FALSE)
}
