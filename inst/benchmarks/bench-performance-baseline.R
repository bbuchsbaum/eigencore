#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (args$quick) 1L else 5L
tol <- 1e-8

failure_rows <- function(case, kind, methods, requested, message) {
  data.frame(
    case = case,
    kind = kind,
    method = methods,
    median = NA_real_,
    min = NA_real_,
    mem_alloc = NA_real_,
    max_residual = NA_real_,
    max_backward_error = NA_real_,
    orthogonality_loss = NA_real_,
    certificate_passed = FALSE,
    certificate_type = NA_character_,
    norm_bound_type = NA_character_,
    scale_is_estimate = NA,
    nconv = 0L,
    requested = requested,
    error = message,
    stringsAsFactors = FALSE
  )
}

with_case <- function(case, kind, methods, requested, expr) {
  tryCatch({
    rows <- expr
    rows$case <- case
    rows$kind <- kind
    rows$requested <- requested
    rows$error <- ""
    rows
  }, error = function(e) {
    failure_rows(case, kind, methods, requested, conditionMessage(e))
  })
}

eigen_methods <- c("eigencore", "RSpectra", "PRIMME")
svd_methods <- c("eigencore", "eigencore_randomized", "RSpectra", "PRIMME", "irlba", "rsvd")

dense_sizes <- if (args$quick) c(180L) else c(500L, 1000L, 3000L)
dense_ks <- if (args$quick) c(5L) else c(5L, 20L)
sparse_sizes <- if (args$quick) c(200L) else c(1000L, 10000L)
sparse_ks <- if (args$quick) c(5L) else c(10L, 20L)
svd_dims <- if (args$quick) list(c(500L, 80L, 5L)) else list(c(100000L, 500L, 20L), c(100000L, 500L, 50L))
slow_svd_dims <- if (args$quick) list(c(160L, 80L, 5L)) else list(c(2000L, 500L, 20L), c(2000L, 500L, 50L))
gen_sizes <- if (args$quick) c(40L) else c(150L, 500L)

rows <- list()

for (n in dense_sizes) {
  A <- dense_low_rank_spd(n, rank = min(25L, n), seed = 1000L + n)
  for (k in dense_ks[dense_ks < n]) {
    case <- paste0("dense_partial_hermitian_n", n, "_k", k)
    rows[[length(rows) + 1L]] <- with_case(
      case, "eigen", eigen_methods, k,
      benchmark_eigen_case(A, k = k, target = largest(), methods = eigen_methods,
                           iterations = iterations, tol = tol, seed = 100L + n + k)
    )
  }
}

for (n in sparse_sizes) {
  A <- path_laplacian(n)
  for (k in sparse_ks[sparse_ks < n]) {
    case <- paste0("sparse_laplacian_smallest_n", n, "_k", k)
    rows[[length(rows) + 1L]] <- with_case(
      case, "eigen", eigen_methods, k,
      benchmark_eigen_case(A, k = k, target = smallest(), methods = eigen_methods,
                           iterations = iterations, tol = tol, seed = 200L + n + k)
    )
  }
}

for (spec in svd_dims) {
  m <- spec[[1L]]
  n <- spec[[2L]]
  rank <- spec[[3L]]
  A <- tall_skinny_sparse(m, n, density = if (args$quick) 0.03 else 0.002,
                          seed = 300L + m + n + rank)
  case <- paste0("tall_skinny_sparse_svd_", m, "x", n, "_k", rank)
  rows[[length(rows) + 1L]] <- with_case(
    case, "svd", svd_methods, rank,
    benchmark_svd_case(A, rank = rank, methods = svd_methods,
                       iterations = iterations, tol = tol, seed = 300L + rank)
  )
}

for (spec in slow_svd_dims) {
  m <- spec[[1L]]
  n <- spec[[2L]]
  rank <- spec[[3L]]
  A <- slow_decay_svd_matrix(m, n, decay = 0.35, seed = 400L + m + n + rank)
  case <- paste0("slow_decay_dense_svd_", m, "x", n, "_k", rank)
  rows[[length(rows) + 1L]] <- with_case(
    case, "svd", svd_methods, rank,
    benchmark_svd_case(A, rank = rank, methods = svd_methods,
                       iterations = iterations, tol = tol, seed = 400L + rank)
  )
}

for (n in gen_sizes) {
  pair <- generalized_spd_pair(n, rank = min(12L, n), sparse = FALSE, seed = 500L + n)
  k <- if (args$quick) 3L else min(10L, n - 1L)
  case <- paste0("dense_generalized_spd_n", n, "_k", k)
  rows[[length(rows) + 1L]] <- with_case(
    case, "generalized_eigen", c("eigencore", "base"), k,
    benchmark_generalized_eigen_case(pair$A, pair$B, k = k, target = smallest(),
                                     methods = c("eigencore", "base"),
                                     iterations = iterations, tol = tol,
                                     seed = 500L + n + k)
  )
}

result <- do.call(rbind, rows)
row.names(result) <- NULL

gate_rows <- lapply(split(result, result$case), function(case_rows) {
  requested <- unique(case_rows$requested)
  kind <- unique(case_rows$kind)
  gate <- evaluate_reference_gate(
    case_rows,
    subject = "eigencore",
    references = setdiff(unique(case_rows$method), "eigencore"),
    requested = requested[[1L]],
    speed_ratio_required = release_speed_gate(kind[[1L]]),
    memory_ratio_required = release_memory_gate(kind[[1L]])
  )
  gate$case <- unique(case_rows$case)
  gate$kind <- kind[[1L]]
  gate
})
gates <- do.call(rbind, gate_rows)
row.names(gates) <- NULL

cat("eigencore performance baseline rows\n")
print(result)
cat("\neigencore performance baseline gates\n")
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "performance-baseline-rows"))
  message("saved gates: ", save_benchmark_result(gates, "performance-baseline-gates"))
}

if (args$strict && !all(gates$passed)) {
  quit(status = 1L)
}
