#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (is.na(args$iterations)) {
  if (args$quick) 1L else 3L
} else {
  args$iterations
}
methods <- args$methods %||% c("eigencore", "RSpectra", "irlba")
gate_subject <- args$subject %||% "eigencore"
if (!gate_subject %in% methods) {
  stop("Implicit-Gram gate subject is not in --methods.", call. = FALSE)
}

rank_deficient_dense <- function(m, n, intrinsic_rank, seed) {
  set.seed(seed)
  matrix(stats::rnorm(m * intrinsic_rank), m, intrinsic_rank) %*%
    matrix(stats::rnorm(intrinsic_rank * n), intrinsic_rank, n)
}

implicit_gram_cases <- function(quick = FALSE) {
  if (isTRUE(quick)) {
    set.seed(1102L)
    return(list(
      list(
        case = "dense_tall_implicit",
        id = "dense_tall_implicit:600x550",
        A = matrix(stats::rnorm(600L * 550L), 600L, 550L),
        rank = 6L,
        orientation = "tall",
        storage = "dense"
      ),
      list(
        case = "sparse_tall_implicit",
        id = "sparse_tall_implicit:4000x900",
        A = Matrix::rsparsematrix(4000L, 900L, density = 0.01),
        rank = 8L,
        orientation = "tall",
        storage = "dgCMatrix"
      ),
      list(
        case = "sparse_wide_implicit",
        id = "sparse_wide_implicit:1100x6000",
        A = Matrix::rsparsematrix(1100L, 6000L, density = 0.005),
        rank = 5L,
        orientation = "wide",
        storage = "dgCMatrix"
      ),
      list(
        case = "dense_rank_deficient_fallback",
        id = "dense_rank_deficient_fallback:180x100",
        A = rank_deficient_dense(180L, 100L, 3L, seed = 117L),
        rank = 6L,
        orientation = "tall",
        storage = "dense",
        fallback_probe = TRUE,
        methods = "eigencore"
      )
    ))
  }

  set.seed(2102L)
  sparse_tall <- Matrix::rsparsematrix(20000L, 5000L, density = 0.001)
  list(
    list(
      case = "dense_tall_implicit",
      id = "dense_tall_implicit:4000x1000",
      A = matrix(stats::rnorm(4000L * 1000L), 4000L, 1000L),
      rank = 10L,
      orientation = "tall",
      storage = "dense"
    ),
    list(
      case = "dense_square_implicit",
      id = "dense_square_implicit:2000x2000",
      A = matrix(stats::rnorm(2000L * 2000L), 2000L, 2000L),
      rank = 10L,
      orientation = "square",
      storage = "dense"
    ),
    list(
      case = "sparse_tall_implicit_k10",
      id = "sparse_tall_implicit_k10:20000x5000",
      A = sparse_tall,
      rank = 10L,
      orientation = "tall",
      storage = "dgCMatrix"
    ),
    list(
      case = "sparse_tall_implicit_k50",
      id = "sparse_tall_implicit_k50:20000x5000",
      A = sparse_tall,
      rank = 50L,
      orientation = "tall",
      storage = "dgCMatrix"
    ),
    list(
      case = "sparse_wide_implicit",
      id = "sparse_wide_implicit:1500x10000",
      A = Matrix::rsparsematrix(1500L, 10000L, density = 0.001),
      rank = 10L,
      orientation = "wide",
      storage = "dgCMatrix"
    ),
    list(
      case = "dense_rank_deficient_fallback",
      id = "dense_rank_deficient_fallback:180x100",
      A = rank_deficient_dense(180L, 100L, 3L, seed = 117L),
      rank = 6L,
      orientation = "tall",
      storage = "dense",
      fallback_probe = TRUE,
      methods = "eigencore"
    )
  )
}

run_implicit_gram_case <- function(case, methods, iterations, seed) {
  if (isTRUE(case$fallback_probe)) {
    old_options <- options(
      eigencore.gram_svd_max_dimension = 64L,
      eigencore.gram_svd_max_dimension_wide = 64L
    )
    on.exit(options(old_options), add = TRUE)
  }
  benchmark_svd_case(
    case$A,
    rank = case$rank,
    methods = case$methods %||% methods,
    iterations = iterations,
    tol = 1e-8,
    seed = seed
  )
}

cases <- filter_benchmark_cases(implicit_gram_cases(args$quick), args$cases)
rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  message_benchmark_case("bench-implicit-gram-svd", case)
  out <- run_implicit_gram_case(
    case,
    methods = methods,
    iterations = iterations,
    seed = 3100L + i
  )
  out$case <- case$case
  out$m <- nrow(case$A)
  out$n <- ncol(case$A)
  out$rank <- case$rank
  out$orientation <- case$orientation
  out$storage <- case$storage
  out$fallback_probe <- isTRUE(case$fallback_probe)
  out
})
result <- do.call(rbind, rows)
row.names(result) <- NULL

display_fields <- c(
  "case", "method", "solver_label", "m", "n", "rank", "median",
  "mem_alloc", "certificate_passed", "nconv", "max_left_residual",
  "max_right_residual", "max_backward_error", "normal_operator_implicit",
  "materialized_gram", "certified_in_original_coordinates",
  "fallback_attempted", "fallback_used", "fallback_method"
)
print(result[, display_fields, drop = FALSE])

contracts <- lapply(cases, function(case) {
  row <- result[
    result$case == case$case & result$method == gate_subject,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Implicit-Gram contract requires exactly one gate-subject row.", call. = FALSE)
  }
  reference_rows <- result[
    result$case == case$case & result$method != gate_subject,
    ,
    drop = FALSE
  ]
  two_sided_certificate <- isTRUE(row$certificate_passed) &&
    isTRUE(row$nconv >= row$rank) &&
    is.finite(row$max_left_residual) &&
    is.finite(row$max_right_residual) &&
    is.finite(row$max_backward_error) &&
    identical(row$certificate_type, "residual_backward_error") &&
    !isTRUE(row$scale_is_estimate)
  references_certified <- if (isTRUE(case$fallback_probe)) {
    TRUE
  } else {
    nrow(reference_rows) > 0L &&
      all(reference_rows$certificate_passed) &&
      all(reference_rows$nconv >= reference_rows$rank)
  }
  best_reference <- if (nrow(reference_rows)) {
    reference_rows[which.min(reference_rows$median), , drop = FALSE]
  } else {
    NULL
  }
  speed_ratio <- if (is.null(best_reference)) {
    NA_real_
  } else {
    best_reference$median / row$median
  }
  memory_ratio <- if (is.null(best_reference)) {
    NA_real_
  } else {
    best_reference$mem_alloc / row$mem_alloc
  }

  if (isTRUE(case$fallback_probe)) {
    provenance <- identical(
      row$solver_label,
      "native prototype Golub-Kahan fallback from implicit Gram SVD"
    ) &&
      isTRUE(row$fallback_attempted) &&
      isTRUE(row$fallback_used) &&
      identical(row$fallback_method, "native prototype Golub-Kahan") &&
      identical(row$implicit_gram_certificate_passed, FALSE) &&
      is.finite(row$implicit_gram_max_backward_error)
  } else {
    provenance <- identical(
      row$solver_label,
      "native certified implicit Gram SVD (thick-restart Lanczos)"
    ) &&
      isTRUE(row$normal_operator_implicit) &&
      !isTRUE(row$materialized_gram) &&
      isTRUE(row$certified_in_original_coordinates) &&
      !isTRUE(row$fallback_attempted) &&
      !isTRUE(row$fallback_used)
  }

  data.frame(
    case = case$case,
    storage = case$storage,
    orientation = case$orientation,
    rank = case$rank,
    fallback_probe = isTRUE(case$fallback_probe),
    two_sided_certificate = two_sided_certificate,
    provenance = provenance,
    references_certified = references_certified,
    best_reference = if (is.null(best_reference)) NA_character_ else best_reference$method,
    certified_time_ratio = speed_ratio,
    memory_ratio = memory_ratio,
    passed = two_sided_certificate && provenance && references_certified,
    stringsAsFactors = FALSE
  )
})
contracts <- do.call(rbind, contracts)
row.names(contracts) <- NULL

cat("\nImplicit-Gram SVD contracts\n")
print(contracts)

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "implicit-gram-svd-rows"))
  message(
    "saved contracts: ",
    save_benchmark_result(contracts, "implicit-gram-svd-contracts")
  )
}

if (args$strict && !all(contracts$passed)) {
  stop("Implicit-Gram SVD benchmark contract failed.", call. = FALSE)
}
