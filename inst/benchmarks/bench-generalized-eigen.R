#!/usr/bin/env Rscript

# Generalized eigen replacement-surface benchmark and planner gates.
#
# Gate families (see mote bd-01KVWKVFC8QFZSTR6KFQ3MM3NH):
# - dense_full_generalized: dense SPD, dense general real, dense general
#   complex, beta-zero/infinite classification.
# - sparse_general_pencil_partial: diagonal-B transformed Arnoldi plus the
#   explicit unsupported boundary (no dense fallback, no densification).
# - qz_dense: generalized_schur real/complex, unsorted and sorted.
# - gsvd_dense: real dense GSVD reconstruction.
#
# Sparse SPD partial gates (diagonal/CSC/matrix-free B, constraints) live in
# inst/benchmarks/bench-generalized-lobpcg.R and are not duplicated here.
#
# Usage:
#   Rscript inst/benchmarks/bench-generalized-eigen.R --quick
#   Rscript inst/benchmarks/bench-generalized-eigen.R --strict --save

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (!is.na(args$iterations)) args$iterations else if (args$quick) 1L else 5L
tol <- 1e-8

dense_n <- if (args$quick) 40L else 200L
sparse_n <- if (args$quick) 60L else 400L

`%||%` <- function(x, y) if (is.null(x)) y else x

new_gate_row <- function(family, case, checks, timing = NULL) {
  failed <- names(checks)[!vapply(checks, isTRUE, logical(1))]
  median_s <- timing$median %||% NA_real_
  baseline_s <- timing$baseline_median %||% NA_real_
  data.frame(
    family = family,
    case = case,
    median = median_s,
    baseline_median = baseline_s,
    speed_ratio_vs_baseline = if (is.finite(median_s) && is.finite(baseline_s)) {
      baseline_s / median_s
    } else {
      NA_real_
    },
    checks_total = length(checks),
    checks_failed = length(failed),
    failed_checks = paste(failed, collapse = ","),
    passed = !length(failed),
    stringsAsFactors = FALSE
  )
}

gate_case <- function(family, case, expr) {
  tryCatch(
    expr,
    error = function(e) {
      new_gate_row(family, case, list(no_error = FALSE)) |>
        transform(failed_checks = paste0("error: ", conditionMessage(e)))
    }
  )
}

sorted_complex <- function(z) {
  z <- as.complex(z)
  z[order(round(Re(z), 8), round(Im(z), 8))]
}

phase_free_match <- function(actual, expected, tolerance = 1e-6) {
  max(Mod(sorted_complex(actual) - sorted_complex(expected))) <= tolerance
}

rows <- list()

## ---------------------------------------------------------------------------
## Family 1: dense_full_generalized
## ---------------------------------------------------------------------------

rows$dense_spd <- gate_case("dense_full_generalized", "spd_real", {
  set.seed(101)
  A <- crossprod(matrix(rnorm(dense_n^2), dense_n)) + diag(dense_n)
  B <- crossprod(matrix(rnorm(dense_n^2), dense_n)) + diag(dense_n)
  timed <- run_timed(eig_full(A, B = B, tol = tol), iterations = iterations)
  fit <- timed$value
  baseline <- run_timed({
    eig <- eigen(solve(B, A), only.values = TRUE)
    eig$values
  }, iterations = iterations)
  V <- vectors(fit)
  b_orth <- max(abs(crossprod(V, B %*% V) - diag(dense_n)))
  oracle <- sort(eigen(solve(B, A), only.values = TRUE)$values)
  new_gate_row(
    "dense_full_generalized", "spd_real",
    list(
      native_label = identical(fit$method, eigencore:::native_dense_generalized_spd_full_label()),
      promotion_gate = identical(fit$plan$controls$promotion_gate, "dense_full_generalized:spd"),
      certificate_passed = isTRUE(certificate(fit)$passed),
      exact_scale = !isTRUE(certificate(fit)$scale_is_estimate),
      oracle_match = max(abs(sort(values(fit)) - oracle)) <= 1e-6 * max(abs(oracle)),
      b_orthogonality = b_orth <= 1e-8,
      all_finite = all(fit$classification == "finite")
    ),
    timing = list(median = timed$median, baseline_median = baseline$median)
  )
})

rows$dense_general_real <- gate_case("dense_full_generalized", "general_pencil_real", {
  set.seed(102)
  A <- matrix(rnorm(dense_n^2), dense_n)
  B <- crossprod(matrix(rnorm(dense_n^2), dense_n)) + diag(dense_n)
  timed <- run_timed(
    eig_full(A, B = B, structure = general(), tol = tol),
    iterations = iterations
  )
  fit <- timed$value
  baseline <- run_timed(
    eigen(solve(B, A), only.values = TRUE)$values,
    iterations = iterations
  )
  oracle <- eigen(solve(B, A), only.values = TRUE)$values
  W <- left_vectors(fit)
  left_residual <- Conj(t(W)) %*% A - diag(values(fit)) %*% (Conj(t(W)) %*% B)
  new_gate_row(
    "dense_full_generalized", "general_pencil_real",
    list(
      native_label = identical(fit$method, eigencore:::native_dense_generalized_pencil_full_label()),
      promotion_gate = identical(fit$plan$controls$promotion_gate,
                                 "dense_full_generalized:general_pencil_real"),
      certificate_passed = isTRUE(certificate(fit)$passed),
      original_coordinate_backward_error = max(fit$certificate$backward_error) <= tol,
      oracle_match = phase_free_match(values(fit), oracle),
      left_vectors_certified = max(Mod(left_residual)) <= 1e-6,
      conditioning_available = isTRUE(fit$conditioning$available),
      norm_scaled_policy = identical(fit$classification_policy$policy, "pencil_norm_scaled")
    ),
    timing = list(median = timed$median, baseline_median = baseline$median)
  )
})

rows$dense_general_complex <- gate_case("dense_full_generalized", "general_pencil_complex", {
  set.seed(103)
  nc <- max(10L, dense_n %/% 4L)
  A <- matrix(complex(real = rnorm(nc^2), imaginary = rnorm(nc^2)), nc)
  B <- matrix(complex(real = rnorm(nc^2), imaginary = rnorm(nc^2)), nc) + diag(nc) * 4
  timed <- run_timed(
    eig_full(A, B = B, structure = general(), tol = tol),
    iterations = iterations
  )
  fit <- timed$value
  oracle <- eigen(solve(B, A), only.values = TRUE)$values
  W <- left_vectors(fit)
  left_residual <- Conj(t(W)) %*% A - diag(values(fit)) %*% (Conj(t(W)) %*% B)
  new_gate_row(
    "dense_full_generalized", "general_pencil_complex",
    list(
      native_label = identical(fit$method, eigencore:::native_dense_generalized_pencil_full_label()),
      promotion_gate = identical(fit$plan$controls$promotion_gate,
                                 "dense_full_generalized:general_pencil_complex"),
      certificate_passed = isTRUE(certificate(fit)$passed),
      oracle_match = phase_free_match(values(fit), oracle),
      left_vectors_certified = max(Mod(left_residual)) <= 1e-6,
      conditioning_boundary_documented = !isTRUE(fit$conditioning$available) &&
        nzchar(fit$conditioning$note %||% "")
    ),
    timing = list(median = timed$median)
  )
})

rows$dense_beta_zero <- gate_case("dense_full_generalized", "beta_zero_infinite", {
  A <- diag(c(2, 3, 5, 7))
  B <- diag(c(1, 1, 1, 0))
  fit <- eig_full(A, B = B, structure = general(), tol = tol)
  scaled <- eig_full(1e-10 * A, B = 1e-10 * B, structure = general(), tol = tol)
  new_gate_row(
    "dense_full_generalized", "beta_zero_infinite",
    list(
      infinite_labelled = identical(fit$classification,
                                    c("finite", "finite", "finite", "infinite")),
      infinite_value = is.infinite(values(fit)[[4L]]),
      finite_certified = all(certificate(fit)$converged[1:3]),
      infinite_not_certified = !certificate(fit)$converged[[4L]],
      scale_invariant = identical(scaled$classification, fit$classification),
      policy_recorded = identical(alpha_beta(fit)$classification_policy$policy,
                                  "pencil_norm_scaled")
    )
  )
})

## ---------------------------------------------------------------------------
## Family 2: sparse_general_pencil_partial
## ---------------------------------------------------------------------------

rows$sparse_diag_b <- gate_case("sparse_general_pencil_partial", "diagonal_B_transform", {
  set.seed(104)
  A <- Matrix::bandSparse(
    sparse_n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, sparse_n - 1L), seq(2, 4, length.out = sparse_n),
                     rep(0.5, sparse_n - 1L))
  )
  A <- methods::as(A, "generalMatrix")
  A <- methods::as(A, "CsparseMatrix")
  B <- Matrix::Diagonal(x = seq(1, 3, length.out = sparse_n))
  k <- 2L
  timed <- run_timed(
    eig_partial(A, B = B, k = k, target = largest_real(), tol = 1e-9,
                allow_dense_fallback = "never"),
    iterations = iterations
  )
  fit <- timed$value
  baseline <- run_timed({
    vals <- eigen(solve(as.matrix(B), as.matrix(A)), only.values = TRUE)$values
    vals[order(Re(vals), decreasing = TRUE)][seq_len(k)]
  }, iterations = iterations)
  oracle <- {
    vals <- eigen(solve(as.matrix(B), as.matrix(A)), only.values = TRUE)$values
    vals[order(Re(vals), decreasing = TRUE)][seq_len(k)]
  }
  new_gate_row(
    "sparse_general_pencil_partial", "diagonal_B_transform",
    list(
      native_label = identical(fit$method,
                               eigencore:::sparse_general_pencil_diagonal_arnoldi_label()),
      promotion_gate = identical(fit$plan$controls$promotion_gate,
                                 "sparse_general_pencil_partial:diagonal_B"),
      certificate_passed = isTRUE(certificate(fit)$passed),
      original_coordinate_residual = identical(
        fit$transform$certification$residual_formula,
        "A * x - lambda * B * x"
      ),
      requested_count_returned = length(values(fit)) == k,
      oracle_match = phase_free_match(values(fit), oracle, tolerance = 1e-6),
      no_densification = !isTRUE(fit$restart$materialized_dense_operator),
      sparse_transform = isTRUE(fit$restart$materialized_sparse_operator),
      native_provenance = isTRUE(fit$restart$native),
      alpha_beta_finite = all(fit$classification == "finite")
    ),
    timing = list(median = timed$median, baseline_median = baseline$median)
  )
})

rows$sparse_unsupported <- gate_case("sparse_general_pencil_partial", "unsupported_boundary", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3), j = c(1, 2, 2, 3), x = c(4, 1, 2, -1), dims = c(3, 3)
  )
  B_singular <- Matrix::Diagonal(x = c(1, 0, 3))
  B_general <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3), j = c(1, 2, 2, 3), x = c(1, 1, 2, 3), dims = c(3, 3)
  )
  plan_singular <- plan_solver(
    eigen_problem(A, metric = B_singular, structure = general()), k = 2L
  )
  plan_general <- plan_solver(
    eigen_problem(A, metric = B_general, structure = general()), k = 2L
  )
  solve_fails <- tryCatch({
    solve(eigen_problem(A, metric = B_general, structure = general()),
          k = 2L, allow_dense_fallback = "never")
    FALSE
  }, error = function(e) grepl("nonsingular diagonal B", conditionMessage(e)))
  new_gate_row(
    "sparse_general_pencil_partial", "unsupported_boundary",
    list(
      singular_B_unsupported_label = identical(
        plan_singular$method, eigencore:::sparse_general_pencil_unsupported_label()
      ),
      general_B_unsupported_label = identical(
        plan_general$method, eigencore:::sparse_general_pencil_unsupported_label()
      ),
      fails_before_dense_fallback = isTRUE(solve_fails),
      no_native_label_leak = !grepl("^native", plan_general$method)
    )
  )
})

## ---------------------------------------------------------------------------
## Family 3: qz_dense
## ---------------------------------------------------------------------------

rows$qz_real <- gate_case("qz_dense", "real_unsorted_and_sorted", {
  set.seed(105)
  nq <- max(12L, dense_n %/% 4L)
  A <- matrix(rnorm(nq^2), nq)
  B <- matrix(rnorm(nq^2), nq)
  timed <- run_timed(generalized_schur(A, B), iterations = iterations)
  qz <- timed$value
  recon_A <- max(abs(A - qz$Q %*% qz$S %*% t(qz$Z)))
  recon_B <- max(abs(B - qz$Q %*% qz$T %*% t(qz$Z)))
  oracle <- eigen(solve(B, A), only.values = TRUE)$values
  pencil_sorted <- generalized_schur(diag(c(2, 3, 0)), diag(c(1, 0, 0)),
                                     sort = "infinite")
  new_gate_row(
    "qz_dense", "real_unsorted_and_sorted",
    list(
      native_label = identical(qz$method, eigencore:::native_dense_generalized_schur_label()),
      promotion_gate = identical(qz$plan$controls$promotion_gate, "qz_dense:generalized_schur"),
      reconstruction_A = recon_A <= 1e-8 * max(1, max(abs(A))),
      reconstruction_B = recon_B <= 1e-8 * max(1, max(abs(B))),
      oracle_match = phase_free_match(qz$values, oracle, tolerance = 1e-6),
      norm_scaled_policy = identical(qz$classification_policy$policy, "pencil_norm_scaled"),
      sorted_sdim = identical(pencil_sorted$sdim, 1L),
      sorted_leading_infinite = identical(pencil_sorted$classification[[1L]], "infinite")
    ),
    timing = list(median = timed$median)
  )
})

rows$qz_complex <- gate_case("qz_dense", "complex", {
  set.seed(106)
  nq <- max(10L, dense_n %/% 5L)
  A <- matrix(complex(real = rnorm(nq^2), imaginary = rnorm(nq^2)), nq)
  B <- matrix(complex(real = rnorm(nq^2), imaginary = rnorm(nq^2)), nq)
  timed <- run_timed(generalized_schur(A, B), iterations = iterations)
  qz <- timed$value
  recon_A <- max(Mod(A - qz$Q %*% qz$S %*% Conj(t(qz$Z))))
  recon_B <- max(Mod(B - qz$Q %*% qz$T %*% Conj(t(qz$Z))))
  oracle <- eigen(solve(B, A), only.values = TRUE)$values
  new_gate_row(
    "qz_dense", "complex",
    list(
      native_label = identical(qz$method, eigencore:::native_dense_generalized_schur_label()),
      reconstruction_A = recon_A <= 1e-8 * max(1, max(Mod(A))),
      reconstruction_B = recon_B <= 1e-8 * max(1, max(Mod(B))),
      oracle_match = phase_free_match(qz$values, oracle, tolerance = 1e-6),
      norm_scaled_policy = identical(qz$classification_policy$policy, "pencil_norm_scaled")
    ),
    timing = list(median = timed$median)
  )
})

## ---------------------------------------------------------------------------
## Family 4: gsvd_dense
## ---------------------------------------------------------------------------

rows$gsvd_real <- gate_case("gsvd_dense", "real_reconstruction", {
  set.seed(107)
  m <- max(8L, dense_n %/% 5L)
  p <- m + 2L
  nn <- m - 1L
  A <- matrix(rnorm(m * nn), m, nn)
  B <- matrix(rnorm(p * nn), p, nn)
  timed <- run_timed(generalized_svd(A, B, tol = tol), iterations = iterations)
  fit <- timed$value
  new_gate_row(
    "gsvd_dense", "real_reconstruction",
    list(
      native_label = identical(fit$method, eigencore:::native_dense_generalized_svd_label()),
      promotion_gate = identical(fit$plan$controls$promotion_gate, "gsvd_dense:real"),
      certificate_passed = isTRUE(certificate(fit)$passed),
      reconstruction_backward_error = max(fit$backward_error) <= tol,
      orthogonality = max(fit$orthogonality) <= 1e-8,
      no_densification_flag = !isTRUE(fit$plan$controls$sparse_densified)
    ),
    timing = list(median = timed$median)
  )
})

## ---------------------------------------------------------------------------
## Label honesty gate: native labels only on native paths
## ---------------------------------------------------------------------------

rows$label_honesty <- gate_case("planner_labels", "native_label_honesty", {
  # Reference generalized SPD LOBPCG (matrix-free) must not carry a native
  # label; the promoted sparse general-pencil boundary must carry one.
  mf_diag <- seq(1, 2, length.out = 8L)
  mf_B <- linear_operator(
    dim = c(8L, 8L),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (mf_diag * X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = hermitian(),
    name = "mf_diag_B"
  )
  problem <- eigen_problem(Matrix::Diagonal(x = seq(1, 8)), metric = mf_B,
                           structure = hermitian())
  plan_reference <- plan_solver(problem, k = 2L, method = lobpcg())
  reference_is_honest <- grepl("reference", plan_reference$method, fixed = TRUE) &&
    !grepl("^native", plan_reference$method)

  A_sparse <- methods::as(Matrix::sparseMatrix(
    i = c(1, 1, 2, 3), j = c(1, 2, 2, 3), x = c(4, 1, 2, -1), dims = c(3, 3)
  ), "generalMatrix")
  plan_native <- plan_solver(
    eigen_problem(A_sparse, metric = Matrix::Diagonal(x = c(1, 2, 3)),
                  structure = general()),
    k = 2L
  )
  native_is_labelled <- identical(
    plan_native$method,
    eigencore:::sparse_general_pencil_diagonal_arnoldi_label()
  )
  new_gate_row(
    "planner_labels", "native_label_honesty",
    list(
      reference_path_not_native = reference_is_honest,
      native_path_labelled = native_is_labelled,
      promotion_gate_exposed = identical(plan_native$controls$promotion_gate,
                                         "sparse_general_pencil_partial:diagonal_B"),
      target_family_exposed = identical(plan_native$controls$target_family,
                                        "sparse_general_pencil_partial"),
      dense_fallback_policy_exposed = nzchar(plan_native$controls$dense_fallback_policy %||% ""),
      alpha_beta_semantics_exposed = nzchar(plan_native$controls$alpha_beta_semantics %||% "")
    )
  )
})

gates <- do.call(rbind, rows)
row.names(gates) <- NULL

cat("Generalized eigen replacement-surface gates\n")
print(gates)

if (args$save) {
  message("saved gates: ", save_benchmark_result(gates, "generalized-eigen-gates"))
}

if (args$strict && !all(gates$passed)) {
  failed <- gates[!gates$passed, c("family", "case", "failed_checks")]
  print(failed)
  stop("Generalized eigen benchmark failed release gate.", call. = FALSE)
}
