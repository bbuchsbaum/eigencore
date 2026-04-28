test_that("benchmark harness produces certificate-inclusive rows", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  eigen_rows <- benchmark_eigen_case(
    A,
    k = 2,
    target = smallest(),
    methods = "eigencore",
    iterations = 1,
    seed = 300
  )

  X <- tall_skinny_sparse(40, 12, density = 0.08, seed = 301)
  svd_rows <- benchmark_svd_case(
    X,
    rank = 2,
    methods = "eigencore",
    iterations = 1,
    seed = 301
  )

  required_common <- c(
    "method", "median", "mem_alloc", "max_residual", "max_backward_error",
    "solver_median", "solver_mem_alloc", "certificate_median",
    "certificate_mem_alloc", "total_median", "total_mem_alloc",
    "orthogonality_loss", "certificate_passed", "certificate_type",
    "norm_bound_type", "scale_is_estimate", "nconv", "iterations", "matvecs",
    "restarts", "ortho_passes", "locking_events", "block_size",
    "stage_apply_seconds", "stage_recurrence_seconds",
    "stage_reorthogonalization_seconds", "stage_projected_solve_seconds",
    "stage_projection_update_seconds",
    "stage_projection_copy_seconds", "stage_projected_eigensolve_seconds",
    "stage_selected_vector_copy_seconds", "stage_ritz_residual_seconds",
    "stage_ritz_vector_form_seconds", "stage_ritz_operator_apply_seconds",
    "stage_ritz_norm_seconds", "stage_ritz_final_polish_seconds",
    "preconditioner_kind", "preconditioner_native", "preconditioner_calls",
    "seed", "pkg_version"
  )
  required_svd <- c(
    required_common,
    "restart_attempts", "final_iterations", "final_matvecs",
    "projected_stop_requested", "projected_stop_enabled",
    "projected_stop_disable_reason", "projected_stop", "projected_nconv",
    "projected_max_residual", "projected_checks", "projected_seconds",
    "native_workspace_bytes", "basis_returned",
    "first_certified_prefix", "final_prefix_iteration_overshoot",
    "final_prefix_matvec_overshoot",
    "stage_native_iteration_seconds", "stage_golub_kahan_ritz_seconds",
    "stage_retry_overhead_seconds",
    "fallback_attempted", "fallback_used", "fallback_method",
    "gram_max_backward_error", "fallback_max_backward_error"
  )
  expect_true(all(required_common %in% names(eigen_rows)))
  expect_true(all(required_svd %in% names(svd_rows)))
  expect_true(all(is.finite(eigen_rows$median)))
  expect_true(all(is.finite(svd_rows$median)))
  expect_equal(eigen_rows$median, eigen_rows$total_median)
  expect_equal(eigen_rows$mem_alloc, eigen_rows$total_mem_alloc)
  expect_equal(svd_rows$median, svd_rows$total_median)
  expect_equal(svd_rows$mem_alloc, svd_rows$total_mem_alloc)
  expect_true(all(is.finite(eigen_rows$solver_median)))
  expect_true(all(is.finite(eigen_rows$certificate_median)))
  expect_true(all(is.finite(svd_rows$solver_median)))
  expect_true(all(is.finite(svd_rows$certificate_median)))
  expect_true(all(is.finite(eigen_rows$max_backward_error)))
  expect_true(all(is.finite(svd_rows$max_backward_error)))
  expect_true(all(eigen_rows$nconv >= 0))
  expect_true(all(svd_rows$nconv >= 0))
})

test_that("benchmark timing resets seed for each measured iteration", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  set.seed(410)
  expected <- stats::runif(1)
  timed <- run_timed(stats::runif(1), iterations = 3, seed = 410)

  expect_equal(timed$value, expected)
})

test_that("benchmark argument parser keeps dense diagnostics opt-in", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  expect_false(benchmark_args(character())$include_dense)
  expect_true(benchmark_args("--include-dense")$include_dense)
  expect_true(benchmark_args(c("--quick", "--block-candidate"))$quick)
  expect_false(benchmark_args(c("--quick", "--block-candidate"))$include_dense)
  expect_true(benchmark_args("--projected-stop")$svd_projected_stop)
  expect_true(benchmark_args("--h-candidate")$h_candidate)
  expect_equal(benchmark_args("--iterations=2")$iterations, 2L)
  expect_equal(benchmark_args("--subject=eigencore_golub_kahan_projected")$subject,
               "eigencore_golub_kahan_projected")
  expect_equal(
    benchmark_args("--methods=eigencore,eigencore_golub_kahan")$methods,
    c("eigencore", "eigencore_golub_kahan")
  )
  expect_equal(
    benchmark_args("--cases=tall_sparse,wide_sparse")$cases,
    c("tall_sparse", "wide_sparse")
  )
  expect_error(benchmark_args("--iterations=0"), "positive integer")
})

test_that("SVD surface benchmark script is available", {
  script <- test_path("../../inst/benchmarks/bench-svd-surface.R")
  expect_true(file.exists(script))
  lines <- readLines(script, warn = FALSE)
  expect_true(any(grepl("eigencore_golub_kahan", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_golub_kahan_projected", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_projected_stop", lines, fixed = TRUE)))
  expect_true(any(grepl("h_candidate", lines, fixed = TRUE)))
  expect_true(any(grepl("projected_stop_comparison", lines, fixed = TRUE)))
  expect_true(any(grepl("memory_diagnostics", lines, fixed = TRUE)))
  expect_true(any(grepl("evaluate_memory_diagnostics", lines, fixed = TRUE)))
  expect_true(any(grepl("gate_subject", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_surface_gate_subject", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_surface_default_methods", lines, fixed = TRUE)))
  expect_true(any(grepl("args$cases", lines, fixed = TRUE)))
  expect_true(any(grepl("rank_deficient_sparse", lines, fixed = TRUE)))
  expect_true(any(grepl("clustered_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("slow_decay_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("benchmark_svd_case", lines, fixed = TRUE)))
})

test_that("SVD surface H candidate preset selects projected GK subject", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  args <- benchmark_args("--h-candidate")
  methods <- svd_surface_default_methods(args)
  expect_equal(methods[[1L]], "eigencore_golub_kahan")
  expect_true("eigencore_golub_kahan_projected" %in% methods)
  expect_false("eigencore" %in% methods)
  expect_equal(
    svd_surface_gate_subject(args, methods),
    "eigencore_golub_kahan_projected"
  )

  args <- benchmark_args("--subject=eigencore_golub_kahan_projected")
  methods <- svd_surface_default_methods(args)
  expect_error(
    svd_surface_gate_subject(args, methods),
    "not in the selected methods"
  )

  args <- benchmark_args(c(
    "--subject=eigencore_golub_kahan_projected",
    "--projected-stop"
  ))
  methods <- svd_surface_default_methods(args)
  expect_equal(
    svd_surface_gate_subject(args, methods),
    "eigencore_golub_kahan_projected"
  )
})

test_that("SVD benchmark harness exposes Golub-Kahan candidate separately", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)
  old_options <- options(eigencore.golub_kahan_prefix_diagnostics = TRUE)
  on.exit(options(old_options), add = TRUE)

  A <- Matrix::rsparsematrix(60, 20, density = 0.08)
  rows <- benchmark_svd_case(
    A,
    rank = 3,
    methods = c("eigencore", "eigencore_golub_kahan"),
    iterations = 1,
    seed = 901
  )

  expect_true(all(c("eigencore", "eigencore_golub_kahan") %in% rows$method))
  gk <- rows[rows$method == "eigencore_golub_kahan", , drop = FALSE]
  expect_true(gk$certificate_passed)
  expect_gte(gk$nconv, 3L)
  expect_gt(gk$matvecs, 0L)
  expect_gte(gk$matvecs, gk$final_matvecs)
  expect_gte(gk$restart_attempts, 1L)
  expect_true(is.finite(gk$first_certified_prefix))
  expect_gte(gk$final_iterations, gk$first_certified_prefix)
  expect_gte(gk$final_prefix_iteration_overshoot, 0L)
  expect_true(is.finite(gk$stage_native_iteration_seconds))
  expect_true(is.finite(gk$stage_golub_kahan_ritz_seconds))
})

test_that("default Golub-Kahan candidate is not accidentally fixed-budget", {
  old_options <- options(eigencore.golub_kahan_prefix_diagnostics = TRUE)
  on.exit(options(old_options), add = TRUE)
  A <- Matrix::rsparsematrix(120, 30, density = 0.05)
  fit <- svd_partial(A, rank = 4, method = golub_kahan(), seed = 902)

  expect_false(fit$restart$fixed_max_subspace)
  expect_true(is.data.frame(fit$restart$history))
  expect_gte(fit$restart$final_max_subspace, fit$plan$controls$initial_max_subspace)
  expect_equal(fit$matvecs, sum(fit$restart$history$matvecs))
  expect_equal(fit$iterations, sum(fit$restart$history$iterations))
  expect_equal(fit$restart$total_matvecs, fit$matvecs)
  expect_equal(fit$restart$final_matvecs, utils::tail(fit$restart$history$matvecs, 1L))
  expect_true(is.data.frame(fit$restart$prefix_history))
  expect_true(all(c(
    "prefix_iterations", "prefix_matvecs", "nconv", "certificate_passed",
    "max_residual", "max_backward_error"
  ) %in% names(fit$restart$prefix_history)))
  expect_equal(
    fit$restart$first_certified_prefix,
    min(fit$restart$prefix_history$prefix_iterations[fit$restart$prefix_history$certificate_passed])
  )
  expect_true(all(c("native_iteration", "ritz", "retry_overhead") %in% names(fit$stage_seconds)))
})

test_that("wide Golub-Kahan adaptive default starts with planner budget", {
  set.seed(903)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))
  plan <- plan_solver(svd_problem(A), rank = 5, method = golub_kahan())
  fit <- svd_partial(A, rank = 5, method = golub_kahan(), seed = 903)

  expect_true(plan$controls$adaptive_subspace)
  expect_equal(plan$controls$initial_max_subspace, 60L)
  expect_false(fit$restart$fixed_max_subspace)
  expect_equal(fit$restart$history$max_subspace[[1L]], plan$controls$initial_max_subspace)
  expect_lte(fit$restart$attempts, 1L)
  expect_true(fit$certificate$passed)
})

test_that("opt-in projected Golub-Kahan stop can shorten a wide attempt", {
  old_options <- options(eigencore.golub_kahan_projected_stop = TRUE)
  on.exit(options(old_options), add = TRUE)
  A <- Matrix::t(tall_skinny_sparse(600L, 90L, density = 0.03, seed = 702))
  plan <- plan_solver(svd_problem(A), rank = 5, method = golub_kahan())
  fit <- svd_partial(A, rank = 5, method = golub_kahan(), seed = 777)

  expect_true(fit$restart$projected_stop_enabled)
  expect_true(fit$restart$projected_stop)
  expect_gte(fit$restart$projected_nconv, 5L)
  expect_lt(fit$restart$final_iterations, plan$controls$initial_max_subspace)
  expect_true(fit$certificate$passed)
})

test_that("projected Golub-Kahan auto policy skips high-aspect tall sparse checks", {
  old_options <- options(eigencore.golub_kahan_projected_stop = TRUE)
  on.exit(options(old_options), add = TRUE)
  A <- tall_skinny_sparse(600L, 90L, density = 0.03, seed = 701)
  fit <- svd_partial(A, rank = 5, method = golub_kahan(), seed = 701)

  expect_true(fit$restart$projected_stop_requested)
  expect_false(fit$restart$projected_stop_enabled)
  expect_match(fit$restart$projected_stop_disable_reason, "high-aspect tall sparse")
  expect_equal(fit$restart$projected_checks, 0L)
  expect_true(fit$certificate$passed)
})

test_that("SVD benchmark can expose projected Golub-Kahan as a separate row", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- Matrix::t(tall_skinny_sparse(600L, 90L, density = 0.03, seed = 702))
  rows <- benchmark_svd_case(
    A,
    rank = 5,
    methods = c("eigencore_golub_kahan", "eigencore_golub_kahan_projected"),
    iterations = 1,
    seed = 904
  )

  plain <- rows[rows$method == "eigencore_golub_kahan", , drop = FALSE]
  projected <- rows[rows$method == "eigencore_golub_kahan_projected", , drop = FALSE]
  expect_false(plain$projected_stop_enabled)
  expect_true(projected$projected_stop_requested)
  expect_true(projected$projected_stop_enabled)
  expect_true(projected$projected_stop)
  expect_gt(projected$projected_checks, 0L)
  expect_gte(projected$projected_seconds, 0)
  expect_lt(projected$final_iterations, plain$final_iterations)
  expect_true(projected$certificate_passed)
})

test_that("SVD reference gate can evaluate an explicit H candidate subject", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore", "eigencore_golub_kahan_projected", "RSpectra"),
    median = c(3, 2, 1),
    mem_alloc = c(300, 200, 100),
    certificate_passed = c(TRUE, TRUE, TRUE),
    nconv = c(2L, 2L, 2L)
  )
  gate <- evaluate_reference_gate(
    rows[rows$method != "eigencore", , drop = FALSE],
    subject = "eigencore_golub_kahan_projected",
    references = "RSpectra",
    requested = 2L,
    speed_ratio_required = release_speed_gate("svd")
  )

  expect_equal(gate$subject, "eigencore_golub_kahan_projected")
  expect_false(gate$passed)
  expect_equal(gate$speed_ratio_vs_best_reference, 0.5)
  expect_equal(gate$memory_ratio_vs_best_reference, 0.5)
})

test_that("SVD memory diagnostics expose subject/reference allocation gaps", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore_golub_kahan_projected", "RSpectra", "irlba"),
    mem_alloc = c(200, 80, 100),
    solver_mem_alloc = c(180, 50, 70),
    certificate_mem_alloc = c(20, 30, 30),
    certificate_passed = c(TRUE, TRUE, TRUE),
    nconv = c(2L, 2L, 2L)
  )
  diagnostics <- evaluate_memory_diagnostics(
    rows,
    subject = "eigencore_golub_kahan_projected",
    references = c("RSpectra", "irlba"),
    requested = 2L
  )

  expect_equal(diagnostics$subject, "eigencore_golub_kahan_projected")
  expect_equal(diagnostics$best_reference, "RSpectra")
  expect_equal(diagnostics$total_memory_gap_bytes, 120)
  expect_equal(diagnostics$solver_memory_gap_bytes, 130)
  expect_equal(diagnostics$certificate_memory_gap_bytes, -10)
  expect_equal(diagnostics$total_memory_ratio_vs_best_reference, 0.4)
  expect_equal(diagnostics$solver_memory_ratio_vs_best_reference, 50 / 180)
  expect_equal(diagnostics$subject_solver_memory_fraction, 0.9)
})

test_that("G1 candidate baseline covers required pre-promotion cases", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  cases <- g1_candidate_baseline_cases(quick = TRUE)
  expect_equal(
    vapply(cases, `[[`, character(1), "case"),
    c("path_laplacian", "dense_hermitian", "clustered", "ill_conditioned_diag")
  )

  rows <- benchmark_g1_candidate_baseline(
    quick = TRUE,
    iterations = 1,
    methods = "eigencore_block_candidate"
  )
  expect_equal(unique(rows$method), "eigencore_block_candidate")
  expect_true(all(vapply(cases, `[[`, character(1), "case") %in% rows$case))
  expect_false(any(rows$certificate_type == "method_error"))
  expect_true(all(rows$certificate_passed))
  expect_true(all(rows$nconv >= rows$k))
  expect_true(all(c(
    "solver_median", "certificate_median", "total_median",
    "iterations", "matvecs", "restarts", "ortho_passes",
    "locking_events", "block_size"
  ) %in% names(rows)))
})

test_that("benchmark harness records failed eigen references as uncertified rows", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  row <- benchmark_failed_eigen_row("broken_reference", 304L, simpleError("bad vectors"))
  expect_equal(row$method, "broken_reference")
  expect_false(row$certificate_passed)
  expect_equal(row$nconv, 0L)
  expect_equal(row$certificate_type, "method_error")
  expect_match(row$error, "bad vectors")
})

test_that("promoted block Lanczos exposes G1 counters", {
  A <- diag(c(8, 5, 3, 1, 0))
  fit <- eig_partial(
    A,
    k = 2,
    target = largest(),
    method = lanczos(block = 2L, max_subspace = 4L),
    seed = 302
  )
  diag <- diagnostics(fit)

  expect_equal(fit$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_gte(fit$restarts, 0L)
  expect_gte(fit$locking_events, 0L)
  expect_equal(fit$block, 2L)
  expect_gte(fit$ortho_passes, 2L)
  expect_equal(diag$restart$ortho_passes, fit$ortho_passes)
  expect_equal(diag$restart$locking_events, fit$locking_events)
  expect_true(all(c(
    "apply", "recurrence", "reorthogonalization", "projected_solve",
    "projection_update", "projection_copy", "projected_eigensolve",
    "selected_vector_copy",
    "ritz_residual", "ritz_vector_form", "ritz_operator_apply",
    "ritz_norm", "ritz_final_polish", "locking", "restart"
  ) %in% names(fit$stage_seconds)))
  expect_true(all(is.finite(fit$stage_seconds)))
  projected_parts <- unname(fit$stage_seconds[c(
    "projection_update", "projection_copy", "projected_eigensolve"
  )])
  expect_equal(unname(fit$stage_seconds[["projected_solve"]]), sum(projected_parts))
  ritz_parts <- unname(fit$stage_seconds[c(
    "ritz_vector_form", "ritz_operator_apply", "ritz_norm", "ritz_final_polish"
  )])
  expect_equal(unname(fit$stage_seconds[["ritz_residual"]]), sum(ritz_parts))
  expect_true(sum(fit$stage_seconds) > 0)
  expect_equal(diag$stage_seconds, fit$stage_seconds)
  expect_s3_class(fit$restart$history, "data.frame")
  expect_true(all(c(
    "restart", "m_active", "selected_count", "locked_before",
    "locked_after", "nconv_wanted", "max_residual", "max_backward_error"
  ) %in% names(fit$restart$history)))
  expect_equal(diag$convergence_history, fit$restart$history)
  expect_gte(nrow(fit$restart$history), 1L)
  expect_true(all(is.finite(fit$restart$history$max_backward_error)))
})

test_that("benchmark rows expose typed preconditioner diagnostics", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  rows <- benchmark_eigen_case(
    A,
    k = 2,
    target = smallest(),
    methods = "eigencore_lobpcg_tridiagonal",
    iterations = 1,
    seed = 303
  )

  expect_equal(rows$preconditioner_kind, "shifted_tridiagonal")
  expect_true(rows$preconditioner_native)
  expect_gte(rows$preconditioner_calls, 1L)
})

test_that("native shifted-tridiagonal LOBPCG wrapper is native and guarded", {
  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(40)
  fit <- eigencore:::native_lobpcg_tridiagonal_hermitian(
    as_operator(A),
    k = 2,
    target = smallest(),
    maxit = 40L,
    tol = 1e-8,
    seed = 304
  )
  expect_true(fit$certificate$passed)
  expect_equal(fit$preconditioner$kind, "shifted_tridiagonal")
  expect_true(fit$preconditioner$native)
  expect_gte(fit$preconditioner_calls, 1L)

  B <- Matrix::sparseMatrix(
    i = c(1L, 4L),
    j = c(4L, 1L),
    x = c(1, 1),
    dims = c(4L, 4L)
  )
  B <- as_operator(B, structure = hermitian())
  expect_error(
    eigencore:::native_lobpcg_tridiagonal_hermitian(B, k = 1, maxit = 5L),
    "not tridiagonal"
  )
})

test_that("native Hermitian benchmark gate reports current pass/fail state", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")
  skip_if_not_installed("RSpectra")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  out <- benchmark_native_hermitian_gate(
    A,
    k = 2,
    target = smallest(),
    methods = c("eigencore", "RSpectra"),
    iterations = 1,
    seed = 702
  )

  expect_true(all(c("rows", "gate") %in% names(out)))
  expect_true(all(c(
    "eigencore_certified", "speed_ratio_vs_best_reference",
    "parity_ratio_vs_best_reference", "parity_gate", "subject", "passed"
  ) %in% names(out$gate)))
  expect_equal(out$gate$requested, 2)
  expect_equal(out$gate$subject, "eigencore")
  expect_true(is.logical(out$gate$passed))
  expect_true(all(c("eigencore", "RSpectra") %in% out$rows$method))
})

test_that("Hermitian benchmark harness can gate the explicit block candidate", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")
  skip_if_not_installed("RSpectra")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- path_laplacian(30)
  out <- benchmark_native_hermitian_gate(
    A,
    k = 2,
    target = smallest(),
    methods = c("eigencore_block_candidate", "eigencore", "RSpectra"),
    subject = "eigencore_block_candidate",
    iterations = 1,
    seed = 703
  )

  expect_true("eigencore_block_candidate" %in% out$rows$method)
  expect_equal(out$gate$subject, "eigencore_block_candidate")
  expect_equal(out$gate$requested, 2)
  expect_true(is.logical(out$gate$passed))
  block_row <- out$rows[out$rows$method == "eigencore_block_candidate", , drop = FALSE]
  expect_equal(block_row$block_size, 2L)
  expect_gte(block_row$restarts, 0L)
  expect_true(is.finite(block_row$stage_reorthogonalization_seconds))
  expect_gt(block_row$stage_reorthogonalization_seconds, 0)
})

test_that("G1 block controls select adaptive large-sparse settings", {
  small <- eigencore:::benchmark_block_candidate_lanczos_method(
    Matrix::Diagonal(1000L),
    20L
  )
  large <- eigencore:::benchmark_block_candidate_lanczos_method(
    Matrix::Diagonal(5000L),
    20L
  )

  expect_equal(small$block, 2L)
  expect_null(small$max_subspace)
  expect_equal(large$block, 4L)
  expect_equal(large$max_subspace, 320L)
})

test_that("G1 RSpectra reference uses robust large-sparse controls", {
  small <- eigencore:::benchmark_rspectra_eigen_opts(
    Matrix::Diagonal(1000L),
    20L,
    smallest()
  )
  large <- eigencore:::benchmark_rspectra_eigen_opts(
    Matrix::sparseMatrix(i = 1:5000, j = 1:5000, x = 1, dims = c(5000L, 5000L)),
    20L,
    smallest()
  )

  expect_equal(small, list())
  expect_equal(large$ncv, 120L)
  expect_equal(large$maxitr, 20000L)
})

test_that("auto planner promotes only benchmark-proven block regimes", {
  mid <- Matrix::bandSparse(
    1000L,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, 999L), c(1, rep(2, 998L), 1), rep(-1, 999L))
  )
  mid_plan <- plan_solver(eigen_problem(mid, target = smallest()), k = 20L)
  expect_equal(mid_plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(mid_plan$controls$block, 2L)
  expect_equal(mid_plan$controls$max_subspace, 160L)

  large <- Matrix::sparseMatrix(
    i = 1:5000,
    j = 1:5000,
    x = 1,
    dims = c(5000L, 5000L)
  )
  large_plan <- plan_solver(eigen_problem(large, target = smallest()), k = 20L)
  expect_equal(large_plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(large_plan$controls$block, 4L)
  expect_equal(large_plan$controls$max_subspace, 320L)

  small_k <- plan_solver(eigen_problem(mid, target = smallest()), k = 5L)
  expect_equal(small_k$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(small_k$controls$block, 1L)

  dense <- diag(as.numeric(seq(200L, 1L)))
  dense_plan <- plan_solver(eigen_problem(dense, target = largest()), k = 20L)
  expect_equal(dense_plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(dense_plan$controls$block, 2L)
  expect_equal(dense_plan$controls$max_subspace, 200L)
})

test_that("native Hermitian gate separates RSpectra threshold from PRIMME parity", {
  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore", "RSpectra", "PRIMME"),
    median = c(1.0, 1.3, 0.9),
    mem_alloc = c(1.0, 1.1, 1.1),
    certificate_passed = c(TRUE, TRUE, TRUE),
    nconv = c(2L, 2L, 2L)
  )
  gate <- evaluate_native_hermitian_gate(
    rows,
    k = 2L,
    reference_methods = "RSpectra",
    parity_methods = "PRIMME"
  )

  expect_true(gate$speed_gate)
  expect_true(gate$memory_gate)
  expect_false(gate$parity_gate)
  expect_false(gate$passed)
  expect_equal(gate$subject, "eigencore")
  expect_equal(gate$speed_reference_methods, "RSpectra")
  expect_equal(gate$parity_reference_methods, "PRIMME")
})

test_that("native Hermitian gate records uncertified references as failed rows", {
  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore", "RSpectra", "PRIMME"),
    median = c(1.0, 0.5, 0.8),
    mem_alloc = c(1.0, 1.1, 1.1),
    certificate_passed = c(TRUE, FALSE, FALSE),
    nconv = c(2L, 0L, 0L)
  )
  gate <- evaluate_native_hermitian_gate(
    rows,
    k = 2L,
    reference_methods = "RSpectra",
    parity_methods = "PRIMME"
  )

  expect_false(gate$speed_gate)
  expect_false(gate$parity_gate)
  expect_false(gate$passed)
  expect_true(is.na(gate$speed_ratio_vs_best_reference))
  expect_match(gate$note, "no certified speed reference rows", fixed = TRUE)
  expect_match(gate$note, "parity reference present but uncertified", fixed = TRUE)
})

test_that("performance baseline helpers cover release regimes", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)
  if (!exists("dense_low_rank_spd", mode = "function")) {
    source(test_path("../../inst/benchmarks/_helpers.R"))
  }

  dense <- dense_low_rank_spd(20, rank = 4, seed = 11)
  expect_equal(dim(dense), c(20L, 20L))
  expect_true(isSymmetric(dense))

  slow <- slow_decay_svd_matrix(18, 9, seed = 12)
  expect_equal(dim(slow), c(18L, 9L))

  pair <- generalized_spd_pair(12, sparse = FALSE, seed = 13)
  expect_equal(dim(pair$A), c(12L, 12L))
  expect_equal(dim(pair$B), c(12L, 12L))

  rows <- data.frame(
    method = c("eigencore", "reference"),
    median = c(2, 1),
    mem_alloc = c(4, 2),
    certificate_passed = c(TRUE, TRUE),
    nconv = c(2L, 2L)
  )
  gate <- evaluate_reference_gate(rows, requested = 2L)
  expect_false(gate$passed)
  expect_equal(gate$speed_ratio_vs_best_reference, 0.5)
  expect_equal(gate$memory_ratio_vs_best_reference, 0.5)
  expect_equal(release_speed_gate("hermitian"), 1.25)
  expect_equal(release_speed_gate("svd"), 1.5)
  expect_equal(release_speed_gate("randomized_svd"), 2.0)
})
