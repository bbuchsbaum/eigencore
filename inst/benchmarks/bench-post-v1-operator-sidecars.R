#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (!is.na(args$iterations)) args$iterations else if (args$quick) 1L else 3L
tol <- 1e-8

sidecar_row <- function(gate_id, surface, case_id, requested, expected_label,
                        expected_native, timed, fit, extra_gate = TRUE,
                        error = "") {
  cert <- fit$certificate %||% certificate(fit)
  restart <- fit$restart %||% list()
  label <- fit$plan$method %||% fit$method %||% NA_character_
  native <- isTRUE(restart$native)
  label_gate <- identical(label, expected_label)
  native_gate <- identical(native, expected_native)
  certificate_gate <- isTRUE(cert$passed) &&
    is.finite(cert$max_backward_error) &&
    cert$max_backward_error <= tol

  data.frame(
    gate_id = gate_id,
    surface = surface,
    case = case_id,
    requested = requested,
    planner_label = label,
    expected_label = expected_label,
    native_boundary = native,
    expected_native_boundary = expected_native,
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
    iterations = result_iterations(fit),
    matvecs = result_matvecs(fit),
    restarts = result_restarts(fit),
    ortho_passes = result_ortho_passes(fit),
    block_size = result_block_size(fit),
    orthogonalization_methods = result_restart_character(fit, "orthogonalization_methods"),
    label_gate = label_gate,
    native_gate = native_gate,
    certificate_gate = certificate_gate,
    extra_gate = isTRUE(extra_gate),
    strict_pass = isTRUE(label_gate) &&
      isTRUE(native_gate) &&
      isTRUE(certificate_gate) &&
      isTRUE(extra_gate),
    seed = NA_integer_,
    pkg_version = as.character(utils::packageVersion("eigencore")),
    error = error,
    stringsAsFactors = FALSE
  )
}

run_sidecar_case <- function(gate_id, surface, case_id, requested,
                             expected_label, expected_native, seed, expr,
                             extra_gate = function(fit) TRUE) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
          "bench-post-v1-operator-sidecars: ", case_id)
  flush.console()
  tryCatch({
    timed <- run_timed(expr, iterations = iterations, seed = seed)
    fit <- timed$value
    row <- sidecar_row(
      gate_id = gate_id,
      surface = surface,
      case_id = case_id,
      requested = requested,
      expected_label = expected_label,
      expected_native = expected_native,
      timed = timed,
      fit = fit,
      extra_gate = extra_gate(fit)
    )
    row$seed <- seed
    row
  }, error = function(e) {
    data.frame(
      gate_id = gate_id,
      surface = surface,
      case = case_id,
      requested = requested,
      planner_label = NA_character_,
      expected_label = expected_label,
      native_boundary = NA,
      expected_native_boundary = expected_native,
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
      iterations = NA_integer_,
      matvecs = NA_integer_,
      restarts = NA_integer_,
      ortho_passes = NA_integer_,
      block_size = NA_integer_,
      orthogonalization_methods = NA_character_,
      label_gate = FALSE,
      native_gate = FALSE,
      certificate_gate = FALSE,
      extra_gate = FALSE,
      strict_pass = FALSE,
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

matrix_free_svd_case <- function() {
  singular <- if (args$quick) c(7, 5, 3, 1) else c(12, 8, 5, 3, 1)
  m <- if (args$quick) 8L else 80L
  n <- length(singular)
  A <- rbind(diag(singular), matrix(0, m - n, n))
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (t(A) %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    name = "post_v1_matrix_free_svd",
    metadata = list(frobenius_norm = norm(A, type = "F"))
  )
  run_sidecar_case(
    gate_id = "post_v1_matrix_free_svd_native_callback_boundary",
    surface = "matrix_free_svd",
    case_id = if (args$quick) "matrix_free_svd:8x4" else "matrix_free_svd:80x5",
    requested = 2L,
    expected_label = eigencore:::native_matrix_free_golub_kahan_label(),
    expected_native = TRUE,
    seed = 9101L,
    expr = svd_partial(op, rank = 2L, tol = tol, seed = 9101L)
  )
}

matrix_free_nonsymmetric_case <- function() {
  A <- if (args$quick) {
    rbind(
      c(5, 2, 0, 0),
      c(0, 3, 1, 0),
      c(0, 0, 1, 0.5),
      c(0, 0, 0, -2)
    )
  } else {
    A <- diag(seq(30, 1))
    A[cbind(seq_len(29L), seq_len(29L) + 1L)] <- 4
    A
  }
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = general(),
    name = "post_v1_matrix_free_nonnormal",
    metadata = list(frobenius_norm = sqrt(sum(A^2)))
  )
  run_sidecar_case(
    gate_id = "post_v1_matrix_free_nonsymmetric_native_boundary",
    surface = "matrix_free_nonsymmetric",
    case_id = if (args$quick) "matrix_free_nonnormal:4" else "matrix_free_nonnormal:30",
    requested = 2L,
    expected_label = "native matrix-free Arnoldi callback cycle + native Ritz extraction",
    expected_native = TRUE,
    seed = 9102L,
    expr = eig_partial(op, k = 2L, target = largest_real(), tol = tol, seed = 9102L)
  )
}

matrix_free_generalized_b_case <- function() {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- diag(c(1, 2, 3, 4, 5, 6))
  Bop <- linear_operator(
    dim = dim(B),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = hermitian(),
    name = "post_v1_explicit_spd_matrix_free_b",
    metadata = list(
      frobenius_norm = sqrt(sum(B^2)),
      positive_definite = TRUE
    )
  )
  run_sidecar_case(
    gate_id = "post_v1_matrix_free_b_native_generalized_contract",
    surface = "matrix_free_generalized_b",
    case_id = "matrix_free_b:diagonal6",
    requested = 2L,
    expected_label = eigencore:::native_generalized_lobpcg_label(),
    expected_native = TRUE,
    seed = 9103L,
    expr = eig_partial(
      A,
      B = Bop,
      k = 2L,
      target = smallest(),
      method = lobpcg(maxit = 80L),
      tol = tol,
      seed = 9103L,
      allow_dense_fallback = "never"
    ),
    extra_gate = function(fit) {
      identical(
        fit$restart$orthogonalization$methods,
        "native_matrix_free_b_mgs2"
      )
    }
  )
}

rows <- do.call(rbind, list(
  matrix_free_svd_case(),
  matrix_free_nonsymmetric_case(),
  matrix_free_generalized_b_case()
))
row.names(rows) <- NULL
print(rows)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "post-v1-operator-sidecars-rows"))
}

if (args$strict && !all(rows$strict_pass)) {
  failed <- rows[!rows$strict_pass, c("gate_id", "case", "planner_label", "error")]
  print(failed)
  stop("post-V1 operator sidecar strict gate failed.", call. = FALSE)
}
