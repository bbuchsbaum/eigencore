#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (is.na(args$iterations)) {
  if (args$quick) 1L else 5L
} else {
  args$iterations
}

methods <- args$methods %||% c("eigencore_randomized", "rsvd")
if (!"eigencore_randomized" %in% methods) {
  stop("randomized-rsvd benchmark requires eigencore_randomized in --methods.", call. = FALSE)
}

cases <- randomized_rsvd_benchmark_cases(quick = args$quick)
cases <- filter_benchmark_cases(cases, args$cases)

rows <- lapply(cases, function(case) {
  message_benchmark_case("bench-randomized-rsvd", case)
  out <- benchmark_randomized_rsvd_case(
    case$A,
    rank = case$rank,
    methods = methods,
    iterations = iterations,
    tol = 1e-8,
    seed = case$seed
  )
  out$case <- case$case
  out$m <- nrow(case$A)
  out$n <- ncol(case$A)
  out$rank <- case$rank
  out
})
result <- do.call(rbind, rows)
row.names(result) <- NULL
print(result)

gates <- lapply(split(result, result$case), function(case_rows) {
  gate <- evaluate_randomized_rsvd_gate(
    case_rows,
    requested = unique(case_rows$rank)[[1L]]
  )
  gate$case <- unique(case_rows$case)
  gate$m <- unique(case_rows$m)
  gate$n <- unique(case_rows$n)
  gate$rank <- unique(case_rows$rank)
  gate
})
gates <- do.call(rbind, gates)
row.names(gates) <- NULL
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "randomized-rsvd-rows"))
  message("saved gates: ", save_benchmark_result(gates, "randomized-rsvd-gates"))
}

if (args$strict && !all(gates$passed)) {
  stop("randomized-rsvd benchmark failed rsvd parity/performance gate.", call. = FALSE)
}
