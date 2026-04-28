#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (!is.na(args$iterations)) args$iterations else if (args$quick) 1L else 5L
tol <- 1e-8

failure_rows <- function(case, target, methods, requested, message) {
  data.frame(
    case = case,
    kind = "generalized_eigen",
    target = target,
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
    iterations = NA_integer_,
    matvecs = NA_integer_,
    preconditioner_kind = NA_character_,
    preconditioner_native = NA,
    preconditioner_calls = NA_integer_,
    seed = NA_integer_,
    pkg_version = as.character(utils::packageVersion("eigencore")),
    error = message,
    stringsAsFactors = FALSE
  )
}

with_case <- function(case, target, methods, requested, expr) {
  tryCatch({
    rows <- expr
    rows$case <- case
    rows$kind <- "generalized_eigen"
    rows$target <- target
    rows$requested <- requested
    rows$error <- ""
    rows
  }, error = function(e) {
    failure_rows(case, target, methods, requested, conditionMessage(e))
  })
}

benchmark_generalized_lobpcg_case <- function(A, B, k, target = smallest(),
                                              methods = c("eigencore", "base"),
                                              iterations = 3L, tol = 1e-8,
                                              seed = 1L, maxit = 200L) {
  methods <- intersect(methods, c("eigencore", "base"))
  rows <- lapply(methods, function(method) {
    timed <- run_timed({
      if (identical(method, "eigencore")) {
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          method = lobpcg(maxit = maxit),
          tol = tol,
          allow_dense_fallback = "never"
        )
      } else {
        eig <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
        idx <- eigencore:::order_indices(eig$values, target)
        idx <- idx[seq_len(k)]
        list(values = eig$values[idx], vectors = eig$vectors[, idx, drop = FALSE])
      }
    }, iterations = iterations, seed = seed)
    cert <- eigencore:::certify_eigen(
      as.matrix(A),
      eigencore:::method_values(timed$value, kind = "eigen"),
      timed$value$vectors,
      B = as.matrix(B),
      tol = tol
    )
    data.frame(
      method = method,
      median = timed$median,
      min = timed$min,
      mem_alloc = timed$mem_alloc,
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      orthogonality_loss = cert$max_orthogonality_loss,
      certificate_passed = cert$passed,
      certificate_type = cert$certificate_type,
      norm_bound_type = cert$norm_bound_type,
      scale_is_estimate = cert$scale_is_estimate,
      nconv = sum(cert$converged),
      iterations = result_iterations(timed$value),
      matvecs = result_matvecs(timed$value),
      restarts = result_restarts(timed$value),
      ortho_passes = result_ortho_passes(timed$value),
      locking_events = result_locking_events(timed$value),
      block_size = result_block_size(timed$value),
      preconditioner_kind = result_preconditioner_field(timed$value, "kind"),
      preconditioner_native = result_preconditioner_field(timed$value, "native"),
      preconditioner_calls = result_preconditioner_calls(timed$value),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

case_specs <- if (args$quick) {
  list(
    list(case = "sparse_generalized_path_smallest", n = 80L, k = 3L, target = smallest(), sparse = TRUE),
    list(case = "sparse_generalized_path_largest", n = 80L, k = 3L, target = largest(), sparse = TRUE)
  )
} else {
  list(
    list(case = "sparse_generalized_path_smallest", n = 500L, k = 10L, target = smallest(), sparse = TRUE),
    list(case = "sparse_generalized_path_largest", n = 500L, k = 10L, target = largest(), sparse = TRUE),
    list(case = "dense_generalized_partial_smallest", n = 180L, k = 8L, target = smallest(), sparse = FALSE),
    list(case = "dense_generalized_partial_largest", n = 180L, k = 8L, target = largest(), sparse = FALSE)
  )
}

methods <- c("eigencore", "base")
rows <- lapply(seq_along(case_specs), function(i) {
  spec <- case_specs[[i]]
  pair <- generalized_spd_pair(
    spec$n,
    rank = min(12L, spec$n),
    sparse = spec$sparse,
    seed = 12000L + spec$n + i
  )
  with_case(
    spec$case,
    eigencore:::target_label(spec$target),
    methods,
    spec$k,
    benchmark_generalized_lobpcg_case(
      pair$A,
      pair$B,
      k = spec$k,
      target = spec$target,
      methods = methods,
      iterations = iterations,
      tol = tol,
      seed = 12100L + spec$n + spec$k + i
    )
  )
})
rows <- do.call(rbind, rows)
row.names(rows) <- NULL

gate_rows <- lapply(split(rows, rows$case), function(case_rows) {
  requested <- unique(case_rows$requested)
  gate <- evaluate_reference_gate(
    case_rows,
    subject = "eigencore",
    references = setdiff(unique(case_rows$method), "eigencore"),
    requested = requested[[1L]],
    speed_ratio_required = if (args$quick) 0 else release_speed_gate("generalized_eigen"),
    memory_ratio_required = if (args$quick) 0 else release_memory_gate("generalized_eigen")
  )
  gate$case <- unique(case_rows$case)
  gate$target <- unique(case_rows$target)
  gate$kind <- "generalized_eigen"
  gate
})
gates <- do.call(rbind, gate_rows)
row.names(gates) <- NULL

cat("Generalized SPD LOBPCG benchmark rows\n")
print(rows)
cat("\nGeneralized SPD LOBPCG gates\n")
print(gates)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "generalized-lobpcg-rows"))
  message("saved gates: ", save_benchmark_result(gates, "generalized-lobpcg-gates"))
}

if (args$strict && !all(gates$passed)) {
  stop("Generalized SPD LOBPCG benchmark failed release gate.", call. = FALSE)
}
