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
    methods = c("eigencore", "eigencore_block_candidate", "RSpectra", "PRIMME"),
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
rows$ortho_passes_per_matvec <- with(rows, ortho_passes / matvecs)
rows$restarts_per_converged <- with(rows, restarts / pmax(nconv, 1L))
gates <- lapply(split(rows, rows$n), function(case_rows) {
  gate <- evaluate_reference_gate(
    case_rows,
    subject = "eigencore_block_candidate",
    references = intersect(c("RSpectra", "PRIMME"), case_rows$method),
    requested = unique(case_rows$k)[[1L]],
    speed_ratio_required = release_speed_gate("hermitian")
  )
  gate$case <- unique(case_rows$case)
  gate$n <- unique(case_rows$n)
  gate$k <- unique(case_rows$k)
  gate
})
gates <- do.call(rbind, gates)
row.names(gates) <- NULL

cat("Native block Hermitian prototype benchmark rows\n")
print(rows)
cat("\nNative block Hermitian prototype gate\n")
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "block-hermitian-prototype-rows"))
  message("saved gates: ", save_benchmark_result(gates, "block-hermitian-prototype-gates"))
}

if (args$strict && !all(gates$passed)) {
  stop("Native block Hermitian prototype benchmark failed promotion gate.", call. = FALSE)
}
