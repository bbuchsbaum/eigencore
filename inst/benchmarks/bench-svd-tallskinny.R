#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
dims <- if (args$quick) {
  list(c(500L, 80L))
} else {
  list(c(100000L, 500L))
}
iterations <- if (args$quick) 1L else 5L
rank <- if (args$quick) 5L else 50L

results <- lapply(seq_along(dims), function(i) {
  m <- dims[[i]][[1L]]
  n <- dims[[i]][[2L]]
  A <- tall_skinny_sparse(m, n, density = if (args$quick) 0.03 else 0.002, seed = 200 + i)
  out <- benchmark_svd_case(
    A,
    rank = rank,
    methods = c("eigencore", "eigencore_randomized", "RSpectra", "PRIMME", "irlba", "rsvd"),
    iterations = iterations,
    seed = 200 + i
  )
  out$case <- "tall_skinny_sparse"
  out$m <- m
  out$n <- n
  out$rank <- rank
  out
})

result <- do.call(rbind, results)
print(result)

svd_gates <- lapply(split(result, paste(result$m, result$n, sep = "x")), function(rows) {
  gate <- evaluate_reference_gate(
    rows[rows$method != "eigencore_randomized", , drop = FALSE],
    requested = unique(rows$rank)[[1L]],
    speed_ratio_required = release_speed_gate("svd")
  )
  gate$case <- unique(rows$case)
  gate$m <- unique(rows$m)
  gate$n <- unique(rows$n)
  gate$gate <- "svd"
  gate
})
randomized_gates <- lapply(split(result, paste(result$m, result$n, sep = "x")), function(rows) {
  gate <- evaluate_reference_gate(
    rows,
    subject = "eigencore_randomized",
    references = intersect(c("rsvd", "irlba", "RSpectra", "PRIMME"), unique(rows$method)),
    requested = unique(rows$rank)[[1L]],
    speed_ratio_required = release_speed_gate("randomized_svd")
  )
  gate$case <- unique(rows$case)
  gate$m <- unique(rows$m)
  gate$n <- unique(rows$n)
  gate$gate <- "randomized_svd"
  gate
})
gates <- do.call(rbind, c(svd_gates, randomized_gates))
row.names(gates) <- NULL
print(gates)

if (args$save) {
  message("saved: ", save_benchmark_result(result, "svd-tallskinny"))
}

if (args$strict && !all(gates$passed)) {
  stop("SVD tall-skinny benchmark failed PRD release gate.", call. = FALSE)
}
