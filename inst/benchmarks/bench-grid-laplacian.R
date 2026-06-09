#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (is.na(args$iterations)) {
  if (args$quick) 1L else 3L
} else {
  args$iterations
}
methods <- args$methods %||% c("eigencore", "RSpectra")
if ("RSpectra" %in% methods && !requireNamespace("RSpectra", quietly = TRUE)) {
  methods <- setdiff(methods, "RSpectra")
}
if (!"eigencore" %in% methods) {
  stop("2D grid benchmark requires eigencore as the diagnostic subject.", call. = FALSE)
}

cases <- if (args$quick) {
  list(
    list(case = "grid_2d_path", id = "grid_2d_path:20x18", nx = 20L, ny = 18L, k = 6L)
  )
} else {
  list(
    list(case = "grid_2d_path", id = "grid_2d_path:100x100", nx = 100L, ny = 100L, k = 10L)
  )
}
cases <- filter_benchmark_cases(cases, args$cases)

rspectra_grid_eigen <- function(A, k, tol = 1e-8) {
  fit <- RSpectra::eigs_sym(A, k = k, which = "SA", opts = list(tol = tol))
  idx <- order(fit$values, decreasing = FALSE)
  values <- fit$values[idx]
  vectors <- fit$vectors[, idx, drop = FALSE]
  residuals <- eigencore:::col_norms(A %*% vectors - sweep(vectors, 2L, values, `*`))
  cert <- eigencore:::certify_eigen_operator_residuals(
    as_operator(A),
    values,
    vectors,
    residuals,
    tol = tol
  )
  list(values = values, vectors = vectors, certificate = cert)
}

benchmark_grid_case <- function(case, methods, iterations, seed = 1L) {
  op <- eigencore:::grid_laplacian_2d_operator(case$nx, case$ny)
  A <- op$metadata$matrix
  rows <- lapply(methods, function(method) {
    timed <- if (identical(method, "eigencore")) {
      run_timed(
        eig_partial(op, k = case$k, target = smallest(), tol = 1e-8),
        iterations = iterations,
        seed = seed
      )
    } else if (identical(method, "RSpectra")) {
      run_timed(
        rspectra_grid_eigen(A, k = case$k, tol = 1e-8),
        iterations = iterations,
        seed = seed
      )
    } else {
      stop("Unsupported 2D grid benchmark method: ", method, call. = FALSE)
    }
    fit <- timed$value
    cert <- fit$certificate
    data.frame(
      method = method,
      solver_label = if (identical(method, "eigencore")) {
        fit$method
      } else {
        NA_character_
      },
      median = timed$median,
      min = timed$min,
      mem_alloc = timed$mem_alloc,
      certificate_passed = cert$passed,
      nconv = sum(cert$converged),
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      scale_is_estimate = cert$scale_is_estimate,
      materialized_dense_operator = if (identical(method, "eigencore")) {
        isTRUE(fit$restart$materialized_dense_operator)
      } else {
        FALSE
      },
      case = case$case,
      nx = case$nx,
      ny = case$ny,
      n = case$nx * case$ny,
      k = case$k,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  message_benchmark_case("bench-grid-laplacian", list(
    case = case$case,
    id = case$id,
    k = case$k
  ))
  benchmark_grid_case(case, methods, iterations, seed = 9500L + i)
})
result <- do.call(rbind, rows)
row.names(result) <- NULL
print(result)

gates <- lapply(split(result, result$case), function(case_rows) {
  eig <- case_rows[case_rows$method == "eigencore", , drop = FALSE]
  ref <- case_rows[case_rows$method == "RSpectra" & case_rows$certificate_passed, , drop = FALSE]
  if (nrow(eig) != 1L || !nrow(ref)) {
    return(data.frame(
      case = unique(case_rows$case),
      subject_certified = nrow(eig) == 1L && isTRUE(eig$certificate_passed),
      speed_ratio_vs_rspectra = NA_real_,
      memory_ratio_vs_rspectra = NA_real_,
      status = "future_only_missing_reference",
      passed = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  speed_ratio <- ref$median / eig$median
  memory_ratio <- ref$mem_alloc / eig$mem_alloc
  passed <- isTRUE(eig$certificate_passed) &&
    isTRUE(eig$nconv >= eig$k) &&
    isTRUE(speed_ratio >= release_speed_gate("hermitian"))
  data.frame(
    case = eig$case,
    subject_certified = eig$certificate_passed,
    speed_ratio_vs_rspectra = speed_ratio,
    memory_ratio_vs_rspectra = memory_ratio,
    status = if (passed) "green_diagnostic_prototype" else "red_future_only",
    passed = passed,
    stringsAsFactors = FALSE
  )
})
gates <- do.call(rbind, gates)
row.names(gates) <- NULL
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "grid-laplacian-rows"))
  message("saved gates: ", save_benchmark_result(gates, "grid-laplacian-gates"))
}

if (args$strict) {
  if (!all(result[result$method == "eigencore", "certificate_passed"])) {
    stop("2D grid benchmark strict mode requires certified eigencore rows.", call. = FALSE)
  }
  if ("RSpectra" %in% methods && !any(result$method == "RSpectra")) {
    stop("2D grid benchmark strict mode requires RSpectra reference rows.", call. = FALSE)
  }
}
