#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (is.na(args$iterations)) {
  if (args$quick) 1L else 5L
} else {
  args$iterations
}
methods <- c(
  "eigencore",
  "eigencore_golub_kahan",
  "eigencore_randomized",
  "RSpectra",
  "PRIMME",
  "irlba",
  "rsvd",
  "base"
)
if (isTRUE(args$svd_projected_stop)) {
  methods <- append(methods, "eigencore_golub_kahan_projected", after = 2L)
}
if (!is.null(args$methods)) {
  methods <- args$methods
}

rank_deficient_sparse <- function(m, n, intrinsic_rank, density = 0.1, seed = 1L) {
  set.seed(seed)
  L <- Matrix::rsparsematrix(m, intrinsic_rank, density = density)
  R <- Matrix::rsparsematrix(intrinsic_rank, n, density = density)
  L %*% R
}

clustered_dense_svd <- function(m, n, rank, seed = 1L) {
  set.seed(seed)
  r <- min(m, n)
  U <- qr.Q(qr(matrix(stats::rnorm(m * r), nrow = m, ncol = r)))
  V <- qr.Q(qr(matrix(stats::rnorm(n * r), nrow = n, ncol = r)))
  d <- c(seq(10, 9.99, length.out = min(rank + 2L, r)),
         if (r > rank + 2L) seq(1, 0.1, length.out = r - rank - 2L) else numeric())
  U %*% (d * t(V))
}

projected_stop_comparison <- function(result) {
  rows <- lapply(split(result, result$case), function(case_rows) {
    plain <- case_rows[case_rows$method == "eigencore_golub_kahan", , drop = FALSE]
    projected <- case_rows[
      case_rows$method == "eigencore_golub_kahan_projected",
      ,
      drop = FALSE
    ]
    if (nrow(plain) != 1L || nrow(projected) != 1L) {
      return(NULL)
    }
    data.frame(
      case = unique(case_rows$case),
      rank = unique(case_rows$rank),
      plain_certified = plain$certificate_passed,
      projected_certified = projected$certificate_passed,
      plain_median = plain$median,
      projected_median = projected$median,
      projected_speed_ratio = plain$median / projected$median,
      plain_mem_alloc = plain$mem_alloc,
      projected_mem_alloc = projected$mem_alloc,
      projected_memory_ratio = plain$mem_alloc / projected$mem_alloc,
      plain_iterations = plain$final_iterations,
      projected_iterations = projected$final_iterations,
      iteration_savings = plain$final_iterations - projected$final_iterations,
      plain_matvecs = plain$final_matvecs,
      projected_matvecs = projected$final_matvecs,
      matvec_savings = plain$final_matvecs - projected$final_matvecs,
      projected_stop_requested = projected$projected_stop_requested,
      projected_stop_enabled = projected$projected_stop_enabled,
      projected_stop_disable_reason = projected$projected_stop_disable_reason,
      projected_stop = projected$projected_stop,
      projected_nconv = projected$projected_nconv,
      projected_max_residual = projected$projected_max_residual,
      plain_projected_checks = plain$projected_checks,
      projected_checks = projected$projected_checks,
      projected_seconds = projected$projected_seconds,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

cases <- if (args$quick) {
  list(
    list(case = "tall_sparse", A = tall_skinny_sparse(600L, 90L, density = 0.03, seed = 701), rank = 5L),
    list(case = "wide_sparse", A = Matrix::t(tall_skinny_sparse(600L, 90L, density = 0.03, seed = 702)), rank = 5L),
    list(case = "rank_deficient_sparse", A = rank_deficient_sparse(160L, 50L, 4L, seed = 703), rank = 6L),
    list(case = "clustered_dense", A = clustered_dense_svd(120L, 60L, 6L, seed = 704), rank = 6L),
    list(case = "slow_decay_dense", A = slow_decay_svd_matrix(120L, 60L, decay = 0.35, seed = 705), rank = 6L),
    list(case = "low_rank_sparse", A = rank_deficient_sparse(300L, 80L, 8L, density = 0.08, seed = 706), rank = 5L)
  )
} else {
  list(
    list(case = "tall_sparse", A = tall_skinny_sparse(100000L, 500L, density = 0.002, seed = 701), rank = 20L),
    list(case = "wide_sparse", A = Matrix::t(tall_skinny_sparse(100000L, 500L, density = 0.002, seed = 702)), rank = 20L),
    list(case = "rank_deficient_sparse", A = rank_deficient_sparse(5000L, 500L, 20L, density = 0.01, seed = 703), rank = 30L),
    list(case = "clustered_dense", A = clustered_dense_svd(2000L, 500L, 20L, seed = 704), rank = 20L),
    list(case = "slow_decay_dense", A = slow_decay_svd_matrix(2000L, 500L, decay = 0.35, seed = 705), rank = 20L),
    list(case = "low_rank_sparse", A = rank_deficient_sparse(10000L, 500L, 40L, density = 0.01, seed = 706), rank = 20L)
  )
}
if (!is.null(args$cases)) {
  wanted_cases <- args$cases
  cases <- Filter(function(case) case$case %in% wanted_cases, cases)
  missing_cases <- setdiff(wanted_cases, vapply(cases, `[[`, character(1), "case"))
  if (length(missing_cases)) {
    stop("Unknown SVD benchmark case(s): ", paste(missing_cases, collapse = ", "), call. = FALSE)
  }
}

rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  active_methods <- if (args$quick && inherits(case$A, "matrix")) {
    methods
  } else {
    setdiff(methods, "base")
  }
  out <- benchmark_svd_case(
    case$A,
    rank = case$rank,
    methods = active_methods,
    iterations = iterations,
    tol = 1e-8,
    seed = 700L + i
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

can_evaluate_gates <- "eigencore" %in% result$method &&
  any(!result$method %in% c(
    "eigencore",
    "eigencore_golub_kahan",
    "eigencore_golub_kahan_projected",
    "eigencore_randomized"
  ))
gates <- if (isTRUE(can_evaluate_gates)) lapply(split(result, result$case), function(case_rows) {
  internal_methods <- c(
    "eigencore_golub_kahan",
    "eigencore_golub_kahan_projected",
    "eigencore_randomized"
  )
  gate <- evaluate_reference_gate(
    case_rows[!case_rows$method %in% internal_methods, , drop = FALSE],
    subject = "eigencore",
    references = setdiff(unique(case_rows$method), c("eigencore", internal_methods)),
    requested = unique(case_rows$rank)[[1L]],
    speed_ratio_required = release_speed_gate("svd"),
    memory_ratio_required = release_memory_gate("svd")
  )
  gate$case <- unique(case_rows$case)
  gate$m <- unique(case_rows$m)
  gate$n <- unique(case_rows$n)
  gate$rank <- unique(case_rows$rank)
  gate
}) else {
  list()
}
gates <- if (length(gates)) do.call(rbind, gates) else data.frame()
row.names(gates) <- NULL
print(gates)

projected_comparison <- NULL
if (isTRUE(args$svd_projected_stop)) {
  projected_comparison <- projected_stop_comparison(result)
  print(projected_comparison)
}

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "svd-surface-rows"))
  message("saved gates: ", save_benchmark_result(gates, "svd-surface-gates"))
  if (!is.null(projected_comparison)) {
    message(
      "saved projected stop comparison: ",
      save_benchmark_result(projected_comparison, "svd-projected-stop-comparison")
    )
  }
}

if (args$strict) {
  if (!nrow(gates)) {
    stop("SVD surface strict mode requires eigencore and at least one external reference.", call. = FALSE)
  }
  if (!all(gates$passed)) {
    stop("SVD surface benchmark failed PRD release gate.", call. = FALSE)
  }
}
