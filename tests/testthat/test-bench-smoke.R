benchmark_file <- function(...) {
  installed <- system.file("benchmarks", ..., package = "eigencore")
  if (nzchar(installed)) {
    return(installed)
  }
  test_path("../../inst/benchmarks", ...)
}

validation_file <- function(...) {
  installed <- system.file("validation", ..., package = "eigencore")
  if (nzchar(installed)) {
    return(installed)
  }
  test_path("../../inst/validation", ...)
}

test_that("benchmark harness produces certificate-inclusive rows", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- benchmark_file("_helpers.R")
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
    "native_workspace_bytes", "native_workspace_allocator", "basis_returned",
    "reorthogonalization_passes",
    "retained_restart", "retained_restart_native", "retained_av_cache",
    "retained_converged_count", "retained_leading_converged_count",
    "native_attempt_certification", "native_early_stop",
    "native_stage_accounted_seconds",
    "stage_reorthogonalization_fraction",
    "reorthogonalization_seconds_per_pass",
    "reorthogonalization_passes_per_iteration",
    "native_seconds_per_matvec", "projected_seconds_per_check",
    "first_certified_prefix", "final_prefix_iteration_overshoot",
    "final_prefix_matvec_overshoot",
    "stage_native_iteration_seconds", "stage_golub_kahan_ritz_seconds",
    "stage_retry_overhead_seconds",
    "attempted_subspaces", "max_attempted_subspace", "max_start_cols",
    "warm_started_attempts", "cached_start_attempts", "certified_attempt",
    "irlba_lbd_small_work_matvecs", "irlba_lbd_fallback_matvecs",
    "irlba_lbd_scout_matvec_overhead_fraction",
    "irlba_lbd_small_work_accounted_seconds",
    "irlba_lbd_fallback_accounted_seconds",
    "irlba_lbd_retained_native_attempted",
    "irlba_lbd_retained_matvecs",
    "irlba_lbd_total_matvecs",
    "irlba_lbd_retained_native_fallback_reason",
    "irlba_lbd_locking_policy", "irlba_lbd_lock_source",
    "irlba_lbd_soft_locked_count", "irlba_lbd_hard_locked_count",
    "irlba_lbd_locked_triplets_certified",
    "irlba_lbd_locked_orthogonality_loss",
    "irlba_lbd_future_vectors_orthogonal_to_locks",
    "irlba_lbd_lock_fallback_reason",
    "irlba_lbd_restart_state_kind",
    "irlba_lbd_recurrence_available",
    "irlba_lbd_augmented_recurrence",
    "irlba_lbd_residual_augmented_cols",
    "irlba_lbd_augmented_tail_steps",
    "irlba_lbd_augmented_basis_cols",
    "irlba_lbd_augmented_restart_cycles",
    "irlba_lbd_augmented_kept_vectors",
    "irlba_lbd_augmented_small_svds",
    "irlba_lbd_augmented_cached_aq_cols",
    "irlba_lbd_augmented_from_scratch_matvecs",
    "irlba_lbd_augmented_matvec_savings",
    "irlba_lbd_augmented_min_cheap_residual",
    "irlba_lbd_augmented_final_cheap_residual",
    "irlba_lbd_augmented_reduces_from_scratch_work",
    "irlba_lbd_native_certificate_diagnostics_reused",
    "irlba_lbd_native_certificate_diagnostics_swapped",
    "irlba_lbd_bpro_policy",
    "irlba_lbd_bpro_passes_per_append",
    "irlba_lbd_bpro_monitoring_threshold",
    "irlba_lbd_bpro_monitored_appends",
    "irlba_lbd_bpro_threshold_reorthogonalizations",
    "irlba_lbd_bpro_max_estimated_orthogonality_loss",
    "irlba_lbd_bpro_max_post_append_orthogonality_loss",
    "irlba_lbd_bpro_basis_orthogonality_loss",
    "irlba_lbd_bpro_escalation_recommended",
    "irlba_lbd_reorth_mode",
    "irlba_lbd_one_sided_reorth_used",
    "irlba_lbd_bpro_block_size",
    "irlba_lbd_bpro_exact_orthogonality_loss",
    "irlba_lbd_bpro_exact_orthogonality_passed",
    "irlba_lbd_bpro_guard_fallback_reason",
    "irlba_lbd_retained_seed_strategy",
    "irlba_lbd_retained_from_scout",
    "irlba_lbd_retained_padding",
    "irlba_lbd_retained_fixed_work_attempts",
    "irlba_lbd_scout_matvecs", "irlba_lbd_scout_certificate_passed",
    "irlba_lbd_normal_scout_attempted", "irlba_lbd_normal_scout_steps",
    "irlba_lbd_normal_scout_chosen_steps", "irlba_lbd_normal_scout_count",
    "irlba_lbd_normal_scout_side", "irlba_lbd_normal_scout_materialized",
    "irlba_lbd_normal_scout_certificate_trusted",
    "irlba_lbd_normal_scout_matvecs",
    "irlba_lbd_normal_scout_operator_matvecs",
    "irlba_lbd_normal_scout_iterations",
    "irlba_lbd_normal_scout_accounted_seconds",
    "irlba_lbd_normal_scout_polish_matvecs",
    "irlba_lbd_normal_scout_polish_iterations",
    "irlba_lbd_normal_scout_polish_accounted_seconds",
    "final_attempt_matvecs", "final_attempt_ortho_passes", "total_ortho_passes",
    "fallback_attempted", "fallback_used", "fallback_method",
    "gram_max_backward_error", "gram_certificate_passed", "gram_dimension",
    "native_gram_kernel", "native_gram_eigensolver",
    "native_gram_subspace_max_backward_error",
    "normal_operator_implicit", "materialized_gram",
    "certified_in_original_coordinates", "fallback_max_backward_error"
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

test_that("Gram SVD cutoff benchmark rows expose 600-side provenance", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  old_options <- options(eigencore.gram_svd_max_dimension = 768L)
  on.exit(options(old_options), add = TRUE)

  X <- tall_skinny_sparse(1200, 600, density = 0.004, seed = 302)
  rows <- benchmark_svd_case(
    X,
    rank = 3,
    methods = "eigencore",
    iterations = 1,
    seed = 302
  )

  eig <- rows[rows$method == "eigencore", , drop = FALSE]
  expect_equal(eig$solver_label, "native certified Gram SVD special case")
  expect_true(eig$certificate_passed)
  expect_equal(eig$gram_dimension, 600L)
  expect_true(eig$materialized_gram)
  expect_false(eig$normal_operator_implicit)
  expect_true(eig$certified_in_original_coordinates)
  expect_match(eig$native_gram_kernel, "csc_.*_gram")
  expect_true(is.finite(eig$max_backward_error))
})

test_that("PRIMME is an optional certified SVD benchmark baseline", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")
  skip_if_not_installed("PRIMME")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  expect_true("PRIMME" %in% benchmark_available_svd_methods())

  X <- tall_skinny_sparse(60, 20, density = 0.08, seed = 303)
  rows <- benchmark_svd_case(
    X,
    rank = 2,
    methods = "PRIMME",
    iterations = 1,
    seed = 303
  )
  expect_true(rows$certificate_passed)
  expect_gte(rows$nconv, 2L)
})

test_that("benchmark timing resets seed for each measured iteration", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  set.seed(410)
  expected <- stats::runif(1)
  timed <- run_timed(stats::runif(1), iterations = 3, seed = 410)

  expect_equal(timed$value, expected)
})

test_that("native validation smoke script is available for sanitizer runs", {
  script <- validation_file("native-smoke.R")
  expect_true(file.exists(script))
  lines <- readLines(script, warn = FALSE)
  required <- c(
    "eigencore native smoke passed",
    "dense Hermitian eigen",
    "sparse CSC Hermitian Lanczos",
    "dense generalized SPD LOBPCG",
    "dense SVD",
    "sparse CSC SVD",
    "dense shift-invert",
    "native tridiagonal shift-invert",
    "--load-all"
  )
  for (needle in required) {
    expect_true(any(grepl(needle, lines, fixed = TRUE)), info = needle)
  }
})

test_that("generalized LOBPCG benchmark exposes native B-orthogonal diagnostics", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  pair <- generalized_spd_pair(40, sparse = TRUE, seed = 421)
  rows <- benchmark_generalized_eigen_case(
    pair$A,
    pair$B,
    k = 2,
    target = smallest(),
    methods = "eigencore",
    iterations = 1,
    seed = 421
  )

  required <- c(
    "native", "native_kernels", "generalized", "orthogonalization_native",
    "orthogonalization_methods", "q_rank_final", "constrained",
    "constraints_rank"
  )
  expect_true(all(required %in% names(rows)))

  eig <- rows[rows$method == "eigencore", , drop = FALSE]
  expect_true(eig$certificate_passed)
  expect_true(eig$native)
  expect_true(eig$native_kernels)
  expect_true(eig$generalized)
  expect_true(eig$orthogonalization_native)
  expect_match(eig$orthogonalization_methods, "native_.*_b_mgs2")
  expect_gte(eig$q_rank_final, 2L)
  expect_false(eig$constrained)
  expect_equal(eig$constraints_rank, 0L)
})

test_that("generalized LOBPCG release script gates native contract rows", {
  script <- benchmark_file("bench-generalized-lobpcg.R")
  expect_true(file.exists(script))
  lines <- readLines(script, warn = FALSE)
  expect_true(any(grepl("eigencore_auto", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_lanczos_native", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_lanczos_block_native", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_lanczos_reference", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_shifted_diagonal", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_shifted_tridiagonal", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lobpcg_tridiagonal_shift", lines, fixed = TRUE)))
  expect_true(any(grepl("as.numeric(n)", lines, fixed = TRUE)))
  expect_true(any(grepl("n = 1000L, k = 10L, target = smallest()", lines, fixed = TRUE)))
  expect_true(any(grepl("n = 1000L, k = 10L, target = largest()", lines, fixed = TRUE)))
  expect_true(any(grepl("maxit = 300L", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_constrained", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lobpcg_native_contract", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lanczos_native_contract", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lanczos_block_native_contract", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lanczos_reference_contract", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized-lanczos-native-contracts", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized-lanczos-block-native-contracts", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized-lanczos-reference-contracts", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lobpcg_adversarial_b_specs", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_lobpcg_adversarial_b_contract", lines, fixed = TRUE)))
  expect_true(any(grepl("filter_benchmark_cases(case_specs, args$cases)", lines, fixed = TRUE)))
  expect_true(any(grepl("message_benchmark_case(\"bench-generalized-lobpcg\"", lines, fixed = TRUE)))
  expect_true(any(grepl("args$methods", lines, fixed = TRUE)))
  expect_true(any(grepl("args$subject", lines, fixed = TRUE)))
  expect_true(any(grepl("preconditioner_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("constraint_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("metric_solve_label", lines, fixed = TRUE)))
  expect_true(any(grepl("metric_solve_native", lines, fixed = TRUE)))
  expect_true(any(grepl("metric_boundary_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("ill_conditioned_diagonal_b", lines, fixed = TRUE)))
  expect_true(any(grepl("diagonal_generalized_lanczos_native_smallest", lines, fixed = TRUE)))
  expect_true(any(grepl("diagonal_generalized_block_lanczos_native_smallest", lines, fixed = TRUE)))
  expect_true(any(grepl("dense_generalized_lanczos_native_largest", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_tridiagonal_b_generalized_lanczos_native_metric_smallest", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_csc_generalized_lanczos_ref_cholesky_smallest", lines, fixed = TRUE)))
  expect_true(any(grepl("symmetric_csc_generalized_lobpcg_smallest_magnitude", lines, fixed = TRUE)))
  expect_true(any(grepl("explicit_spd_matrix_free_b", lines, fixed = TRUE)))
  expect_true(any(grepl("expected_orthogonalization", lines, fixed = TRUE)))
})

test_that("generalized eigen replacement-surface script gates planner families", {
  script <- benchmark_file("bench-generalized-eigen.R")
  expect_true(file.exists(script))
  lines <- readLines(script, warn = FALSE)
  expect_true(any(grepl("dense_full_generalized", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_general_pencil_partial", lines, fixed = TRUE)))
  expect_true(any(grepl("qz_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("gsvd_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("native_dense_generalized_spd_full_label", lines, fixed = TRUE)))
  expect_true(any(grepl("native_dense_generalized_pencil_full_label", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_general_pencil_diagonal_arnoldi_label", lines, fixed = TRUE)))
  expect_true(any(grepl("pencil_norm_scaled", lines, fixed = TRUE)))
  expect_true(any(grepl("left_vectors(fit)", lines, fixed = TRUE)))
  expect_true(any(grepl("diagonal_B_transform", lines, fixed = TRUE)))
  expect_true(any(grepl("unsupported_boundary", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_schur", lines, fixed = TRUE)))
  expect_true(any(grepl("generalized_svd", lines, fixed = TRUE)))
  expect_true(any(grepl("new_gate_row", lines, fixed = TRUE)))
  expect_true(any(grepl("args$strict", lines, fixed = TRUE)))
})

test_that("benchmark argument parser keeps dense diagnostics opt-in", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- benchmark_file("_helpers.R")
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

test_that("benchmark case filtering supports stable ids and names", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  cases <- list(
    list(case = "path_laplacian", n = 1000L, k = 20L),
    list(case = "path_laplacian", n = 10000L, k = 20L),
    list(case = "dense_hermitian", n = 200L, k = 20L)
  )

  expect_equal(benchmark_case_id(cases[[2L]]), "path_laplacian:10000")
  by_id <- filter_benchmark_cases(cases, "path_laplacian:10000")
  expect_length(by_id, 1L)
  expect_equal(by_id[[1L]]$n, 10000L)

  svd_case <- list(case = "tall_sparse", id = "tall_sparse:600x90", rank = 5L)
  expect_equal(benchmark_case_id(svd_case), "tall_sparse:600x90")
  by_svd_id <- filter_benchmark_cases(list(svd_case), "tall_sparse:600x90")
  expect_length(by_svd_id, 1L)

  generalized_case <- list(
    case = "sparse_generalized_path_smallest",
    n = 80L,
    k = 3L,
    methods = c("eigencore", "eigencore_shifted_diagonal", "base")
  )
  expect_equal(
    benchmark_case_id(generalized_case),
    "sparse_generalized_path_smallest:80"
  )

  by_name <- filter_benchmark_cases(cases, "path_laplacian")
  expect_length(by_name, 2L)

  expect_error(
    filter_benchmark_cases(cases, "missing"),
    "Available cases: path_laplacian:1000, path_laplacian:10000, dense_hermitian:200",
    fixed = TRUE
  )
})

test_that("SVD surface benchmark script is available", {
  script <- benchmark_file("bench-svd-surface.R")
  helper <- benchmark_file("_helpers.R")
  expect_true(file.exists(script))
  expect_true(file.exists(helper))
  lines <- c(readLines(script, warn = FALSE), readLines(helper, warn = FALSE))
  expect_true(any(grepl("eigencore_golub_kahan", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_golub_kahan_projected", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_block_golub_kahan_cycle", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_block_golub_kahan_cycle_cached", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_block_golub_kahan_cycle_lean", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_projected_stop", lines, fixed = TRUE)))
  expect_true(any(grepl("h_candidate", lines, fixed = TRUE)))
  expect_true(any(grepl("projected_stop_comparison", lines, fixed = TRUE)))
  expect_true(any(grepl("memory_diagnostics", lines, fixed = TRUE)))
  expect_true(any(grepl("iteration_savings_fraction", lines, fixed = TRUE)))
  expect_true(any(grepl("reorthogonalization_pass_savings_fraction", lines, fixed = TRUE)))
  expect_true(any(grepl("evaluate_memory_diagnostics", lines, fixed = TRUE)))
  expect_true(any(grepl("gate_subject", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_surface_gate_subject", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_surface_default_methods", lines, fixed = TRUE)))
  expect_true(any(grepl("args$cases", lines, fixed = TRUE)))
  expect_true(any(grepl("filter_benchmark_cases(cases, args$cases)", lines, fixed = TRUE)))
  expect_true(any(grepl("message_benchmark_case(\"bench-svd-surface\"", lines, fixed = TRUE)))
  expect_true(any(grepl("tall_sparse:600x90", lines, fixed = TRUE)))
  expect_true(any(grepl("rank_deficient_sparse", lines, fixed = TRUE)))
  expect_true(any(grepl("smallest_sparse", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_smallest", lines, fixed = TRUE)))
  expect_true(any(grepl("base_smallest", lines, fixed = TRUE)))
  expect_true(any(grepl("interior_sparse", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_interior", lines, fixed = TRUE)))
  expect_true(any(grepl("base_interior", lines, fixed = TRUE)))
  expect_true(any(grepl("solver_label", lines, fixed = TRUE)))
  expect_true(any(grepl("svd_target_contract", lines, fixed = TRUE)))
  expect_true(any(grepl("svd-surface-target-contracts", lines, fixed = TRUE)))
  expect_true(any(grepl("native certified Gram SVD special case", lines, fixed = TRUE)))
  expect_true(any(grepl("quick_reference_contract_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("speed and memory ratios are diagnostics", lines, fixed = TRUE)))
  expect_true(any(grepl("clustered_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("complex_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("gate = FALSE", lines, fixed = TRUE)))
  expect_true(any(grepl("slow_decay_dense", lines, fixed = TRUE)))
  expect_true(any(grepl("benchmark_svd_case", lines, fixed = TRUE)))
})

test_that("post-V1 benchmark gate manifest covers hard promotion surfaces", {
  manifest_path <- benchmark_file("post-v1-gate-manifest.R")
  expect_true(file.exists(manifest_path))

  env <- new.env(parent = globalenv())
  sys.source(manifest_path, envir = env)
  manifest <- env$post_v1_gate_manifest

  expect_type(manifest, "list")
  expect_equal(manifest$version, 1L)
  expect_equal(manifest$generated_on, "2026-06-06")
  expect_true(grepl("YYYYMMDD-", manifest$artifact_policy$naming, fixed = TRUE))
  expect_true(all(c("smoke", "strict", "long") %in% names(manifest$tier_profile)))
  expect_equal(manifest$tier_profile$smoke$default_gate_ids, "post_v1_operator_sidecars")
  expect_equal(manifest$tier_profile$strict$default_gate_ids, "all")
  expect_equal(manifest$tier_profile$long$default_gate_ids, "all")
  expect_true(all(c(
    "time_to_certified_answer",
    "memory",
    "max_residual",
    "max_backward_error",
    "orthogonality_loss",
    "planner_label",
    "certificate_type",
    "native_or_reference_boundary",
    "nconv"
  ) %in% manifest$required_metrics))
  current_gate_owners <- c(
    "bd-01KTE8G6RYE4RD5F6CN7SNKKC6",
    "bd-01KVWKVFC8QFZSTR6KFQ3MM3NH"
  )
  expect_equal(manifest$current_gate_owner_issue_ids, current_gate_owners)

  gates <- manifest$gates
  gate_ids <- vapply(gates, `[[`, character(1), "id")
  expect_true(all(c(
    "post_v1_svd_hard_surface",
    "post_v1_operator_sidecars",
    "post_v1_randomized_svd_hard_surface",
    "post_v1_generalized_preconditioned_surface",
    "post_v1_generalized_eigen_surface",
    "post_v1_shift_invert_boundaries",
    "post_v1_nonsymmetric_matrix_free_surface"
  ) %in% gate_ids))

  required_fields <- c(
    "id", "owner_issue", "surface", "script", "command",
    "quick_smoke_command", "long_command", "cases", "baselines",
    "artifacts", "thresholds"
  )
  for (gate in gates) {
    expect_true(all(required_fields %in% names(gate)), info = gate$id)
    expect_true(gate$owner_issue %in% current_gate_owners, info = gate$id)
    expect_true(file.exists(benchmark_file(basename(gate$script))), info = gate$id)
    expect_true(grepl(gate$script, gate$command, fixed = TRUE), info = gate$id)
    expect_true(grepl("Rscript", gate$quick_smoke_command, fixed = TRUE), info = gate$id)
    expect_true(grepl("Rscript", gate$long_command, fixed = TRUE), info = gate$id)
    expect_true(grepl("--iterations=10", gate$long_command, fixed = TRUE), info = gate$id)
    expect_true(length(gate$cases) > 0L, info = gate$id)
    expect_true(length(gate$baselines) > 0L, info = gate$id)
    expect_true(length(gate$thresholds) > 0L, info = gate$id)
    expect_true(all(grepl("\\.rds$", gate$artifacts)), info = gate$id)
  }

  svd_gate <- gates[[match("post_v1_svd_hard_surface", gate_ids)]]
  expect_true(any(grepl("rank_deficient_sparse", svd_gate$cases, fixed = TRUE)))
  expect_true(any(grepl("slow_decay_dense", svd_gate$cases, fixed = TRUE)))
  expect_true("RSpectra" %in% svd_gate$baselines)
  expect_true("irlba" %in% svd_gate$baselines)
  expect_equal(svd_gate$thresholds$speed_ratio_min, 1.10)

  operator_gate <- gates[[match("post_v1_operator_sidecars", gate_ids)]]
  expect_equal(operator_gate$owner_issue, "bd-01KTE8G6RYE4RD5F6CN7SNKKC6")
  expect_true(any(grepl("matrix_free_svd", operator_gate$cases, fixed = TRUE)))
  expect_true(any(grepl("matrix_free_nonnormal", operator_gate$cases, fixed = TRUE)))
  expect_true(any(grepl("matrix_free_b", operator_gate$cases, fixed = TRUE)))
  expect_true(operator_gate$thresholds$planner_label_exact)
  expect_true(operator_gate$thresholds$native_boundary_exact)

  closed_owner_ids <- c(
    "bd-01KTE8J396W5SGK6FQBSQX7BY8",
    "bd-01KTE8JNVAJW5SHR6AHFCH9T4B",
    "bd-01KTE8K131G6EBK3QKYA3PRH0E",
    "bd-01KTE8K6PRRY1SRCFXXQR84YQW",
    "bd-01KTE8J9SF16Y1832D8HQQ9KEC",
    "bd-01KTE8JFKPA90ZJTXK496SBMK4",
    "bd-01KTEH5JM64A4CBZG7ECBWT9WB",
    "bd-01KTE8JVEPGA1EEQYERZS1V7S1",
    "bd-01KTEH6862GB19JJWX2M3FQP6T",
    "bd-01KTEH60X91VZRSW7NGV65FBDR"
  )
  expect_false(any(vapply(gates, function(gate) {
    gate$owner_issue %in% closed_owner_ids
  }, logical(1))))
})

test_that("post-V1 benchmark profile runner and workflow are wired", {
  runner <- benchmark_file("run-post-v1-gates.R")
  workflow <- test_path("../../.github/workflows/post-v1-benchmarks.yaml")
  expect_true(file.exists(runner))

  runner_lines <- readLines(runner, warn = FALSE)
  required_runner <- c(
    "validate_gate_manifest",
    "select_gate_ids",
    "gate_command",
    "--tier=",
    "--gates=",
    "--dry-run",
    "--load-all",
    "command_with_load_all",
    "smoke",
    "strict",
    "long",
    "owner_issue",
    "Post-V1 gate failed"
  )
  for (needle in required_runner) {
    expect_true(any(grepl(needle, runner_lines, fixed = TRUE)), info = needle)
  }

  skip_if_not(
    file.exists(workflow),
    "GitHub workflow is source-only and excluded from the built package"
  )
  expect_true(file.exists(workflow))
  workflow_lines <- readLines(workflow, warn = FALSE)
  required_workflow <- c(
    "workflow_dispatch",
    "schedule:",
    "post-v1-benchmarks.yaml",
    "BENCH_LIB",
    "setup-r-dependencies",
    "Run post-V1 benchmark profile",
    "run-post-v1-gates.R",
    "actions/upload-artifact",
    "inst/benchmarks/results"
  )
  for (needle in required_workflow) {
    expect_true(any(grepl(needle, workflow_lines, fixed = TRUE)), info = needle)
  }
})

test_that("post-V1 operator sidecar benchmark gates matrix-free boundaries", {
  script <- benchmark_file("bench-post-v1-operator-sidecars.R")
  expect_true(file.exists(script))
  lines <- readLines(script, warn = FALSE)

  required <- c(
    "post_v1_matrix_free_svd_native_callback_boundary",
    "post_v1_matrix_free_nonsymmetric_native_boundary",
    "post_v1_matrix_free_b_native_generalized_contract",
    "native_matrix_free_golub_kahan_label",
    "native matrix-free Arnoldi callback cycle + native Ritz extraction",
    "native_matrix_free_b_mgs2",
    "expected_native = TRUE",
    "planner_label",
    "native_boundary",
    "strict_pass",
    "post-v1-operator-sidecars-rows"
  )
  for (needle in required) {
    expect_true(any(grepl(needle, lines, fixed = TRUE)), info = needle)
  }
})

test_that("randomized-rsvd benchmark script is available", {
  script <- benchmark_file("bench-randomized-rsvd.R")
  helper <- benchmark_file("_helpers.R")
  expect_true(file.exists(script))
  expect_true(file.exists(helper))
  lines <- c(readLines(script, warn = FALSE), readLines(helper, warn = FALSE))
  expect_true(any(grepl("benchmark_randomized_rsvd_case", lines, fixed = TRUE)))
  expect_true(any(grepl("evaluate_randomized_rsvd_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("randomized_rsvd_benchmark_cases", lines, fixed = TRUE)))
  expect_true(any(grepl("filter_benchmark_cases(cases, args$cases)", lines, fixed = TRUE)))
  expect_true(any(grepl("message_benchmark_case(\"bench-randomized-rsvd\"", lines, fixed = TRUE)))
  expect_true(any(grepl("release_gate_required", lines, fixed = TRUE)))
  expect_true(any(grepl("exact_low_rank_dense:120x80", lines, fixed = TRUE)))
  expect_true(any(grepl("slow_decay_dense:2000x500", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_randomized", lines, fixed = TRUE)))
  expect_true(any(grepl("rsvd", lines, fixed = TRUE)))
})

test_that("randomized-rsvd cases mark release and diagnostic rows", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  source(benchmark_file("_helpers.R"))
  quick <- randomized_rsvd_benchmark_cases(quick = TRUE)
  full <- randomized_rsvd_benchmark_cases(quick = FALSE)

  quick_required <- vapply(quick, `[[`, logical(1), "release_gate_required")
  full_required <- vapply(full, `[[`, logical(1), "release_gate_required")

  expect_equal(
    vapply(quick[quick_required], `[[`, character(1), "case"),
    character()
  )
  expect_equal(
    vapply(full[full_required], `[[`, character(1), "case"),
    "exact_low_rank_dense"
  )
  expect_true(all(vapply(c(quick, full), function(case) {
    is.character(case$release_gate_note) && nzchar(case$release_gate_note)
  }, logical(1))))
})

test_that("shift-invert benchmark script is available", {
  script <- benchmark_file("bench-shift-invert.R")
  helper <- benchmark_file("_helpers.R")
  expect_true(file.exists(script))
  expect_true(file.exists(helper))
  lines <- c(readLines(script, warn = FALSE), readLines(helper, warn = FALSE))
  expect_true(any(grepl("shift_invert_cases", lines, fixed = TRUE)))
  expect_true(any(grepl("filter_benchmark_cases(shift_invert_cases(args$quick), args$cases)", lines, fixed = TRUE)))
  expect_true(any(grepl("message_benchmark_case(\"bench-shift-invert\"", lines, fixed = TRUE)))
  expect_true(any(grepl("dense_lu_native", lines, fixed = TRUE)))
  expect_true(any(grepl("dense_lu_generalized_native", lines, fixed = TRUE)))
  expect_true(any(grepl("tridiagonal_thomas_native", lines, fixed = TRUE)))
  expect_true(any(grepl("tridiagonal_thomas_generalized_native", lines, fixed = TRUE)))
  expect_true(any(grepl("diagonal_standard_native", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_tridiagonal_generalized_native", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_general_reference", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_general_diagonal_b_reference", lines, fixed = TRUE)))
  expect_true(any(grepl("matrix_free_user_solve_reference", lines, fixed = TRUE)))
  expect_true(any(grepl("shift_invert_matrix_free_user_solve", lines, fixed = TRUE)))
  expect_true(any(grepl("estimated_converged", lines, fixed = TRUE)))
  expect_true(any(grepl("user_solve", lines, fixed = TRUE)))
  expect_true(any(grepl("Matrix::lu", lines, fixed = TRUE)))
  expect_true(any(grepl("shift_invert_factorization_contract_v1", lines, fixed = TRUE)))
  expect_true(any(grepl("factorization_contract_provider", lines, fixed = TRUE)))
  expect_true(any(grepl("factorization_contract_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("eigencore_native_factorization", lines, fixed = TRUE)))
  expect_true(any(grepl("Matrix::lu_reference_factorization", lines, fixed = TRUE)))
  expect_true(any(grepl("user_supplied_solve", lines, fixed = TRUE)))
  expect_true(any(grepl("shift_invert_contract", lines, fixed = TRUE)))
})

test_that("nonsymmetric benchmark script is available", {
  script <- benchmark_file("bench-nonsymmetric.R")
  helper <- benchmark_file("_helpers.R")
  expect_true(file.exists(script))
  expect_true(file.exists(helper))
  lines <- c(readLines(script, warn = FALSE), readLines(helper, warn = FALSE))
  expect_true(any(grepl("nonsymmetric_cases", lines, fixed = TRUE)))
  expect_true(any(grepl("filter_benchmark_cases(nonsymmetric_cases(args$quick), args$cases)", lines, fixed = TRUE)))
  expect_true(any(grepl("message_benchmark_case(\"bench-nonsymmetric\"", lines, fixed = TRUE)))
  expect_true(any(grepl("dense LAPACK general eigen oracle", lines, fixed = TRUE)))
  expect_true(any(grepl("native Arnoldi cycle", lines, fixed = TRUE)))
  expect_true(any(grepl("native Ritz extraction", lines, fixed = TRUE)))
  expect_true(any(grepl("matrix_free_nonnormal", lines, fixed = TRUE)))
  expect_true(any(grepl("native matrix-free Arnoldi callback cycle", lines, fixed = TRUE)))
  expect_true(any(grepl("native_matrix_free_arnoldi_label", lines, fixed = TRUE)))
  expect_true(any(grepl("matrix_free_native", lines, fixed = TRUE)))
  expect_true(any(grepl("reference Arnoldi (prototype/oracle fallback)", lines, fixed = TRUE)))
  expect_true(any(grepl("right_residual_backward_error", lines, fixed = TRUE)))
  expect_true(any(grepl("arnoldi_native", lines, fixed = TRUE)))
  expect_true(any(grepl("restart_count", lines, fixed = TRUE)))
  expect_true(any(grepl("max_restarts", lines, fixed = TRUE)))
  expect_true(any(grepl("restart_gate", lines, fixed = TRUE)))
  expect_true(any(grepl("stage_arnoldi_cycle_seconds", lines, fixed = TRUE)))
  expect_true(any(grepl("stage_ritz_extraction_seconds", lines, fixed = TRUE)))
  expect_true(any(grepl("ritz_extraction_native", lines, fixed = TRUE)))
  expect_true(any(grepl("dense_native_arnoldi_lm", lines, fixed = TRUE)))
  expect_true(any(grepl("dense_eigs_native_arnoldi_li", lines, fixed = TRUE)))
  expect_true(any(grepl("dense_complex_lapack_li", lines, fixed = TRUE)))
  expect_true(any(grepl("native dense complex general LAPACK fallback", lines, fixed = TRUE)))
  expect_true(any(grepl("native_dense_complex_label", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_native_arnoldi_lr", lines, fixed = TRUE)))
  expect_true(any(grepl("sparse_native_arnoldi_li", lines, fixed = TRUE)))
  expect_true(any(grepl("nonsymmetric_contract", lines, fixed = TRUE)))
})

test_that("V2 CRAN benchmark manifest names release benchmark surfaces", {
  manifest <- test_path("../../docs/v1-benchmark-manifest.md")
  skip_if_not(file.exists(manifest), "source docs are not available in installed-package checks")
  expect_true(file.exists(manifest))
  lines <- readLines(manifest, warn = FALSE)
  required <- c(
    "bench-native-hermitian-gate.R",
    "bench-hermitian-sparse.R",
    "bench-svd-surface.R",
    "bench-randomized-rsvd.R",
    "bench-generalized-lobpcg.R",
    "bench-generalized-eigen.R",
    "bench-lobpcg-preconditioned.R",
    "bench-shift-invert.R",
    "bench-nonsymmetric.R",
    "R CMD check --no-manual"
  )
  for (needle in required) {
    expect_true(any(grepl(needle, lines, fixed = TRUE)), info = needle)
  }
  expect_true(any(grepl("not a release signoff", lines, fixed = TRUE)))
  expect_true(any(grepl("V2 CRAN release", lines, fixed = TRUE)))
})

test_that("contribution methods artifact ties claims to evidence and losses", {
  artifact <- test_path("../../docs/contribution-methods-artifact.md")
  skip_if_not(file.exists(artifact), "source docs are not available in installed-package checks")
  lines <- readLines(artifact, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  required <- c(
    "Contribution Thesis",
    "Method Design",
    "Benchmark Report",
    "Losses and boundaries",
    "Migration Package",
    "V3 Deferral Program",
    "V2 does not add solver families beyond the CRAN release surface",
    "docs/v1-benchmark-manifest.md",
    "bench-native-hermitian-gate.R",
    "bench-svd-surface.R",
    "bench-randomized-rsvd.R",
    "bench-generalized-lobpcg.R",
    "bench-generalized-eigen.R",
    "bench-shift-invert.R",
    "bench-nonsymmetric.R",
    "docs/rspectra-migration.md",
    "docs/known-limitations.md",
    "docs/post-v1-benchmark-gates.md",
    "Base complex dense matrices use native dense complex",
    "shift_invert_factorization_contract_v1",
    "bd-01KTE8J9SF16Y1832D8HQQ9KEC",
    "bd-01KTE8JFKPA90ZJTXK496SBMK4",
    "bd-01KTEH4RA4NJBQSKF9JQV3V4BS",
    "bd-01KTEH4G1QPR4RT14B4G78PF1M",
    "bd-01KTEH4ZPPBESDDPR0Z5MRSQCG",
    "bd-01KTEH57J1RR3SP1SJ27YRC0ZE",
    "bd-01KTEH5JM64A4CBZG7ECBWT9WB",
    "bd-01KTEH60X91VZRSW7NGV65FBDR",
    "bd-01KTEH6862GB19JJWX2M3FQP6T",
    "bd-01KTEH5SRWDHXBZXK5CPHBT6G2",
    "bd-01KTEH6HNN33M15YNGW7T35RQR",
    "bd-01KTE8JVEPGA1EEQYERZS1V7S1",
    "Closed No-Promotion Decisions",
    "Closed Non-Goal Decisions",
    "bd-01KTEH48HH4X8G9Q69HHAJ983B",
    "documented no-promotion decision",
    "explicit PRD non-goal",
    "owner ids stay current",
    "Jacobi-Davidson",
    "GraphBLAS/GPU/distributed/SLEPc/PRIMME plugins"
  )
  for (needle in required) {
    expect_true(grepl(needle, text, fixed = TRUE), info = needle)
  }
})

test_that("V2 CRAN completion audit maps the active goal to stop-rule evidence", {
  audit <- test_path("../../docs/v1-completion-audit.md")
  skip_if_not(file.exists(audit), "source docs are not available in installed-package checks")
  lines <- readLines(audit, warn = FALSE)
  required <- c(
    "Package checks remain clean",
    "Public solver paths either run native engine code or carry honest",
    "Production SVD and randomized SVD",
    "Generalized SPD, shift-invert, and nonsymmetric surfaces",
    "Sanitizer / valgrind-style evidence",
    "Current decision: **V2 CRAN release candidate pending final validation and mote closure**",
    "Mark V2 CRAN complete only after",
    "benchmarks/RELEASES.md",
    "mote board"
  )
  for (needle in required) {
    expect_true(any(grepl(needle, lines, fixed = TRUE)), info = needle)
  }
})

test_that("V2 CRAN documentation scope audit names required doc surfaces", {
  audit <- test_path("../../docs/v1-doc-scope-audit.md")
  skip_if_not(file.exists(audit), "source docs are not available in installed-package checks")
  lines <- readLines(audit, warn = FALSE)
  required <- c(
    "README.md",
    "README.Rmd",
    "vignettes/eigencore.Rmd",
    "vignettes/certificates.Rmd",
    "docs/rspectra-migration.md",
    "docs/known-limitations.md",
    "docs/method-selection-and-workflows.md",
    "docs/v1-readiness-audit.md",
    "docs/v1-benchmark-manifest.md",
    "benchmarks/RELEASES.md",
    "docs/native-lobpcg.md",
    "docs/native-generalized-spd-lobpcg.md",
    "docs/native-block-lanczos.md",
    "docs/hegelsvd_svd_acceleration.md",
    "V2 CRAN Documentation Scope Audit",
    "tridiagonal generalized-with-diagonal-B native labels",
    "Documentation Boundaries",
    "Stop Rule"
  )
  for (needle in required) {
    expect_true(any(grepl(needle, lines, fixed = TRUE)), info = needle)
  }
  expect_true(any(grepl("documentation-scope companion", lines, fixed = TRUE)))
})

test_that("prd defines V2 as CRAN release and moves hard solver expansion to V3", {
  prd <- test_path("../../prd.json")
  skip_if_not(file.exists(prd), "source PRD is not available in installed-package checks")
  text <- paste(readLines(prd, warn = FALSE), collapse = "\n")
  required <- c(
    "\"release_strategy\"",
    "\"v2_cran_release\"",
    "Ship eigencore V2 as the first CRAN release",
    "V2 is a release-hardening and publication boundary",
    "\"v2_scope\"",
    "\"release_name\": \"V2 CRAN release\"",
    "\"v3_scope\"",
    "Jacobi-Davidson",
    "full nonsymmetric Krylov-Schur workflows",
    "scalable sparse/matrix-free interior SVD",
    "native general sparse LU ownership",
    "native complex sparse/operator kernels",
    "GraphBLAS-style sparse kernels",
    "optional SLEPc/PETSc plugin",
    "automatic adapters for broad matrix ecosystems"
  )
  for (needle in required) {
    expect_true(grepl(needle, text, fixed = TRUE), info = needle)
  }
  v2_pos <- regexpr("\"v2_scope\"", text, fixed = TRUE)[[1]]
  v3_pos <- regexpr("\"v3_scope\"", text, fixed = TRUE)[[1]]
  jd_pos <- regexpr("Jacobi-Davidson", text, fixed = TRUE)[[1]]
  expect_gt(v2_pos, 0)
  expect_gt(v3_pos, v2_pos)
  expect_gt(jd_pos, v3_pos)
})

test_that("known limitations match current shift-invert boundary", {
  limits <- test_path("../../docs/known-limitations.md")
  skip_if_not(file.exists(limits), "source docs are not available in installed-package checks")
  lines <- readLines(limits, warn = FALSE)
  required <- c(
    "Dense standard, dense generalized SPD, diagonal standard, sparse symmetric-tridiagonal standard, and sparse/diagonal tridiagonal generalized paths are native",
    "shift_invert_factorization_contract_v1",
    "general sparse uses an honest `Matrix::lu` reference contract",
    "user solves use an external-cache contract"
  )
  for (needle in required) {
    expect_true(any(grepl(needle, lines, fixed = TRUE)), info = needle)
  }
})

test_that("known limitations keep graph Fiedler off sparse general-pencil route", {
  limits <- test_path("../../docs/known-limitations.md")
  skip_if_not(file.exists(limits), "source docs are not available in installed-package checks")
  text <- paste(readLines(limits, warn = FALSE), collapse = "\n")
  required <- c(
    "sparse general-pencil `smallest()` on graph-Laplacian-like near-null spectra is not a promoted Fiedler path",
    "`native transformed sparse general-pencil Arnoldi`; this label is a",
    "general-pencil boundary, not graph/Fiedler guidance"
  )
  for (needle in required) {
    expect_true(grepl(needle, text, fixed = TRUE), info = needle)
  }
})

test_that("known limitations document complex ABI certificate contract", {
  limits <- test_path("../../docs/known-limitations.md")
  skip_if_not(file.exists(limits), "source docs are not available in installed-package checks")
  text <- paste(readLines(limits, warn = FALSE), collapse = "\n")
  required <- c(
    "Complex ABI And Certificate Contract",
    "ScalarType::C128",
    "`dtype = \"complex\"`",
    "`apply_adjoint()` to implement the conjugate",
    "`storage = \"complex_dense_matrix\"`",
    "`native = TRUE`",
    "`native_operator_kernel = \"dense_complex_zgemm\"`",
    "`V^* V`",
    "`A v - sigma u`",
    "`A^* u - sigma v`",
    "Explicit dense complex sources use exact Frobenius scales",
    "`norm_bound_type = \"frobenius_metadata\"`",
    "`scale_is_estimate = FALSE`",
    "Complex matrix-free solver operators fail with actionable future-scope messages"
  )
  for (needle in required) {
    expect_true(grepl(needle, text, fixed = TRUE), info = needle)
  }
})

test_that("randomized-rsvd gate enforces accuracy and speed versus rsvd", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore_randomized", "rsvd"),
    median = c(0.5, 1.0),
    certificate_passed = c(TRUE, TRUE),
    nconv = c(4L, 4L),
    singular_value_relative_error = c(1e-8, 1.1e-8),
    left_subspace_error = c(2e-8, 2.2e-8),
    right_subspace_error = c(3e-8, 3.3e-8)
  )
  gate <- evaluate_randomized_rsvd_gate(rows, requested = 4L)
  expect_equal(gate$subject, "eigencore_randomized")
  expect_equal(gate$baseline, "rsvd")
  expect_true(gate$subject_certified)
  expect_true(gate$baseline_certified)
  expect_equal(gate$baseline_nconv, 4L)
  expect_true(gate$accuracy_gate)
  expect_true(gate$speed_gate)
  expect_true(gate$passed)
  expect_equal(gate$note, "")

  rows$singular_value_relative_error[[1L]] <- 1e-4
  failed <- evaluate_randomized_rsvd_gate(rows, requested = 4L)
  expect_false(failed$accuracy_gate)
  expect_false(failed$passed)

  rows$singular_value_relative_error[[1L]] <- 1e-8
  rows$certificate_passed[[2L]] <- FALSE
  uncertified_ref <- evaluate_randomized_rsvd_gate(rows, requested = 4L)
  expect_false(uncertified_ref$baseline_certified)
  expect_equal(uncertified_ref$baseline_nconv, 4L)
  expect_false(uncertified_ref$speed_gate)
  expect_false(uncertified_ref$passed)
  expect_match(uncertified_ref$note, "rsvd baseline did not satisfy")
})

test_that("randomized-rsvd controller contract enforces native dense and sparse provenance", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  source(benchmark_file("_helpers.R"))

  rows <- data.frame(
    method = c("eigencore_randomized", "eigencore_randomized"),
    case = c("exact_low_rank_dense", "low_rank_sparse"),
    rank = c(4L, 4L),
    certificate_passed = c(TRUE, TRUE),
    nconv = c(4L, 4L),
    scale_is_estimate = c(FALSE, FALSE),
    randomized_controller_native = c(TRUE, TRUE),
    randomized_controller_kind = c(
      "native_dense_randomized_controller",
      "native_csc_randomized_controller"
    ),
    randomized_dense_native_controller = c(TRUE, FALSE),
    randomized_sparse_native_controller = c(FALSE, TRUE),
    randomized_native_certificate_diagnostics = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  contracts <- randomized_controller_contract(rows)
  expect_true(all(contracts$certificate_gate))
  expect_true(all(contracts$provenance_gate))
  expect_true(all(contracts$passed))

  rows$randomized_controller_kind[[2L]] <- "randomized_range_finder"
  failed <- randomized_controller_contract(rows)
  expect_false(failed$passed[[2L]])
})

test_that("randomized-rsvd benchmark rows expose native projection diagnostics", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- benchmark_file("_helpers.R")
  source(helper_path)

  dense <- nearly_low_rank_svd_matrix(40L, 24L, rank = 5L, noise = 0, seed = 1601L)
  dense_rows <- benchmark_randomized_rsvd_case(
    dense,
    rank = 3L,
    methods = "eigencore_randomized",
    iterations = 1L,
    seed = 1601L
  )

  set.seed(1603)
  A <- Matrix::rsparsematrix(140L, 12L, density = 0.12) %*%
    Matrix::rsparsematrix(12L, 90L, density = 0.12)
  sparse_rows <- benchmark_randomized_rsvd_case(
    A,
    rank = 8L,
    methods = "eigencore_randomized",
    iterations = 1L,
    seed = 1603L
  )

  for (rows in list(dense_rows, sparse_rows)) {
    expect_true("randomized_native_sketch" %in% names(rows))
    expect_true("randomized_sketch_kind" %in% names(rows))
    expect_true("randomized_controller_native" %in% names(rows))
    expect_true("randomized_controller_kind" %in% names(rows))
    expect_true("randomized_dense_native_controller" %in% names(rows))
    expect_true("randomized_sparse_native_controller" %in% names(rows))
    expect_true("randomized_native_certificate_diagnostics" %in% names(rows))
    expect_true("randomized_adaptive_stop_used" %in% names(rows))
    expect_true("randomized_core_solver" %in% names(rows))
    expect_true("randomized_projection_kind" %in% names(rows))
    expect_true("randomized_projection_transposed" %in% names(rows))
    expect_true(rows$randomized_native_sketch)
    expect_equal(rows$randomized_sketch_kind, "native_fused_a_omega")
    expect_match(rows$randomized_core_solver, "^native_")
    expect_equal(rows$randomized_projection_kind, "native_direct_qt_a")
    expect_true(rows$randomized_projection_transposed)
  }
  expect_true(dense_rows$randomized_controller_native)
  expect_equal(dense_rows$randomized_controller_kind, "native_dense_randomized_controller")
  expect_true(dense_rows$randomized_dense_native_controller)
  expect_false(dense_rows$randomized_sparse_native_controller)
  expect_true(dense_rows$randomized_native_certificate_diagnostics)
  expect_true(sparse_rows$randomized_controller_native)
  expect_equal(sparse_rows$randomized_controller_kind, "native_csc_randomized_controller")
  expect_true(sparse_rows$randomized_sparse_native_controller)
  expect_false(sparse_rows$randomized_dense_native_controller)
  expect_true(sparse_rows$randomized_native_certificate_diagnostics)
})

test_that("SVD surface H candidate preset selects production eigencore subject", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  source(benchmark_file("_helpers.R"))

  args <- benchmark_args("--h-candidate")
  methods <- svd_surface_default_methods(args)
  expect_equal(methods[[1L]], "eigencore")
  expect_true("eigencore_golub_kahan_one_sided" %in% methods)
  expect_true("eigencore_irlba_lbd_one_sided" %in% methods)
  expect_true("eigencore_irlba_lbd_retained_native" %in% methods)
  expect_true("eigencore_irlba_lbd_retained_bpro" %in% methods)
  expect_true("eigencore_irlba_lbd_retained_bpro_one_sided_guarded" %in% methods)
  expect_true("eigencore_irlba_lbd_retained_bpro_block_guarded" %in% methods)
  expect_true("eigencore_irlba_lbd_normal_scout" %in% methods)
  expect_true("eigencore_golub_kahan_projected" %in% methods)
  expect_true("eigencore_implicit_normal_lanczos" %in% methods)
  expect_true("eigencore_gram_dsyevx" %in% methods)
  expect_true("eigencore_block_golub_kahan_cycle" %in% methods)
  expect_true("eigencore_block_golub_kahan_retained" %in% methods)
  expect_true("eigencore_block_golub_kahan_retained_cached" %in% methods)
  expect_equal(
    svd_surface_gate_subject(args, methods),
    "eigencore"
  )

  expect_true("eigencore_block_golub_kahan_retained_cached" %in% svd_internal_methods())
  expect_true("eigencore_implicit_normal_lanczos" %in% svd_internal_methods())
  expect_true("eigencore_gram_dsyevx" %in% svd_internal_methods())
  expect_true("eigencore_smallest" %in% svd_internal_methods())
  expect_true("eigencore_golub_kahan_smallest" %in% svd_internal_methods())
  expect_true("eigencore_interior" %in% svd_internal_methods())
  expect_true("eigencore_golub_kahan_interior" %in% svd_internal_methods())
  expect_true("eigencore_golub_kahan_one_sided" %in% svd_internal_methods())
  expect_true("eigencore_irlba_lbd_one_sided" %in% svd_internal_methods())
  expect_true("eigencore_irlba_lbd_retained_native" %in% svd_internal_methods())
  expect_true("eigencore_irlba_lbd_retained_bpro" %in% svd_internal_methods())
  expect_true("eigencore_irlba_lbd_retained_bpro_one_sided_guarded" %in% svd_internal_methods())
  expect_true("eigencore_irlba_lbd_retained_bpro_block_guarded" %in% svd_internal_methods())
  expect_true("eigencore_irlba_lbd_normal_scout" %in% svd_internal_methods())

  default_methods <- svd_surface_default_methods(benchmark_args(character()))
  expect_false("eigencore_irlba_lbd_retained_bpro_one_sided_guarded" %in% default_methods)
  expect_false("eigencore_irlba_lbd_retained_bpro_block_guarded" %in% default_methods)
  expect_false("eigencore_block_golub_kahan_retained" %in% default_methods)
  expect_false("eigencore_block_golub_kahan_retained_cached" %in% default_methods)

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

  requested <- c("eigencore_irlba_lbd_retained_native", "eigencore_golub_kahan")
  expect_equal(bench_methods("svd", requested), requested)
  expect_error(
    bench_methods("svd", c("eigencore", "definitely_missing_svd_method")),
    "Requested SVD benchmark method"
  )
  expect_error(
    benchmark_svd_case(
      tall_skinny_sparse(20, 8, density = 0.1, seed = 323),
      rank = 2,
      methods = "definitely_missing_svd_method",
      iterations = 1,
      seed = 323
    ),
    "Requested SVD benchmark method"
  )
})

test_that("SVD benchmark exposes guarded BPRO diagnostic rows", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  set.seed(702)
  A <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  rows <- benchmark_svd_case(
    A,
    rank = 5L,
    methods = c(
      "eigencore_irlba_lbd_retained_bpro_one_sided_guarded",
      "eigencore_irlba_lbd_retained_bpro_block_guarded"
    ),
    iterations = 1L,
    seed = 702L
  )

  one_sided <- rows[
    rows$method == "eigencore_irlba_lbd_retained_bpro_one_sided_guarded",
    ,
    drop = FALSE
  ]
  block <- rows[
    rows$method == "eigencore_irlba_lbd_retained_bpro_block_guarded",
    ,
    drop = FALSE
  ]
  expect_true(one_sided$certificate_passed)
  expect_true(block$certificate_passed)
  expect_identical(one_sided$irlba_lbd_reorth_mode, "bpro_one_sided_guarded")
  expect_identical(block$irlba_lbd_reorth_mode, "bpro_block_guarded")
  expect_true(one_sided$irlba_lbd_one_sided_reorth_used)
  expect_false(block$irlba_lbd_one_sided_reorth_used)
  expect_equal(one_sided$irlba_lbd_bpro_block_size, 1L)
  expect_equal(block$irlba_lbd_bpro_block_size, 5L)
  expect_true(one_sided$irlba_lbd_bpro_exact_orthogonality_passed)
  expect_true(block$irlba_lbd_bpro_exact_orthogonality_passed)
  expect_true(is.na(one_sided$irlba_lbd_bpro_guard_fallback_reason))
  expect_true(is.na(block$irlba_lbd_bpro_guard_fallback_reason))
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
  expect_true(is.finite(gk$stage_reorthogonalization_fraction))
  expect_true(is.finite(gk$reorthogonalization_seconds_per_pass))
  expect_true(is.finite(gk$reorthogonalization_passes_per_iteration))
  expect_true(is.finite(gk$native_seconds_per_matvec))
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
  expect_true(is.finite(projected$projected_seconds_per_check))
  expect_true(is.finite(projected$stage_reorthogonalization_fraction))
  expect_true(is.finite(projected$reorthogonalization_seconds_per_pass))
  expect_true(is.finite(projected$reorthogonalization_passes_per_iteration))
  expect_lt(projected$final_iterations, plain$final_iterations)
  expect_true(projected$certificate_passed)
})

test_that("SVD benchmark can expose native block Golub-Kahan cycle candidate", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- Matrix::rsparsematrix(80, 24, density = 0.08)
  rows <- benchmark_svd_case(
    A,
    rank = 4,
    methods = "eigencore_block_golub_kahan_cycle",
    iterations = 1,
    seed = 905
  )

  expect_equal(unique(rows$method), "eigencore_block_golub_kahan_cycle")
  expect_true(rows$certificate_passed)
  expect_gte(rows$nconv, 4L)
  expect_gt(rows$matvecs, 0L)
  expect_equal(rows$block_size, 2L)
  expect_false(rows$basis_returned)
  expect_gt(rows$stage_native_iteration_seconds, 0)
  expect_gt(rows$stage_golub_kahan_ritz_seconds, 0)
})

test_that("SVD benchmark can expose lean native block Golub-Kahan restart candidate", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  set.seed(702)
  A <- Matrix::t(Matrix::rsparsematrix(600, 90, density = 0.03))
  rows <- benchmark_svd_case(
    A,
    rank = 5,
    methods = c(
      "eigencore_block_golub_kahan_cycle",
      "eigencore_block_golub_kahan_cycle_cached",
      "eigencore_block_golub_kahan_cycle_cached_random",
      "eigencore_block_golub_kahan_cycle_residual",
      "eigencore_block_golub_kahan_cycle_lean",
      "eigencore_block_golub_kahan_retained",
      "eigencore_block_golub_kahan_retained_cached"
    ),
    iterations = 1,
    seed = 701
  )

  expect_true(all(c(
    "eigencore_block_golub_kahan_cycle",
    "eigencore_block_golub_kahan_cycle_cached",
    "eigencore_block_golub_kahan_cycle_cached_random",
    "eigencore_block_golub_kahan_cycle_residual",
    "eigencore_block_golub_kahan_cycle_lean",
    "eigencore_block_golub_kahan_retained",
    "eigencore_block_golub_kahan_retained_cached"
  ) %in% rows$method))

  regular <- rows[rows$method == "eigencore_block_golub_kahan_cycle", , drop = FALSE]
  cached <- rows[rows$method == "eigencore_block_golub_kahan_cycle_cached", , drop = FALSE]
  cached_random <- rows[rows$method == "eigencore_block_golub_kahan_cycle_cached_random", , drop = FALSE]
  residual <- rows[rows$method == "eigencore_block_golub_kahan_cycle_residual", , drop = FALSE]
  lean <- rows[rows$method == "eigencore_block_golub_kahan_cycle_lean", , drop = FALSE]
  retained <- rows[rows$method == "eigencore_block_golub_kahan_retained", , drop = FALSE]
  retained_cached <- rows[rows$method == "eigencore_block_golub_kahan_retained_cached", , drop = FALSE]
  expect_true(cached$certificate_passed)
  expect_true(cached_random$certificate_passed)
  expect_true(residual$certificate_passed)
  expect_true(lean$certificate_passed)
  expect_true(retained$certificate_passed)
  expect_true(retained_cached$certificate_passed)
  expect_gte(cached$nconv, 5L)
  expect_gte(cached_random$nconv, 5L)
  expect_gte(residual$nconv, 5L)
  expect_gte(lean$nconv, 5L)
  expect_gte(retained$nconv, 5L)
  expect_gte(retained_cached$nconv, 5L)
  expect_gt(lean$matvecs, 0L)
  expect_gte(lean$restart_attempts, 1L)
  expect_match(lean$attempted_subspaces, ",")
  expect_gte(lean$max_attempted_subspace, 5L)
  expect_gte(lean$max_start_cols, 5L)
  expect_gte(lean$warm_started_attempts, 1L)
  expect_gte(lean$certified_attempt, 1L)
  expect_equal(lean$final_attempt_matvecs, lean$final_matvecs)
  expect_gte(lean$total_ortho_passes, lean$final_attempt_ortho_passes)
  expect_false(regular$basis_returned)
  expect_false(cached$basis_returned)
  expect_false(cached_random$basis_returned)
  expect_false(residual$basis_returned)
  expect_false(lean$basis_returned)
  expect_false(retained$basis_returned)
  expect_false(retained_cached$basis_returned)
  expect_true(retained$retained_restart)
  expect_true(retained$retained_restart_native)
  expect_gt(retained$native_workspace_bytes, 0)
  expect_equal(retained$native_workspace_allocator, "native_malloc")
  expect_false(retained$retained_av_cache)
  expect_true(retained$native_attempt_certification)
  expect_false(retained$native_early_stop)
  expect_false(retained$fallback_attempted)
  expect_false(retained$fallback_used)
  expect_true(retained_cached$retained_restart)
  expect_true(retained_cached$retained_restart_native)
  expect_gt(retained_cached$native_workspace_bytes, 0)
  expect_equal(retained_cached$native_workspace_allocator, "native_malloc")
  expect_true(retained_cached$retained_av_cache)
  expect_true(retained_cached$native_attempt_certification)
  expect_false(retained_cached$fallback_attempted)
  expect_false(retained_cached$fallback_used)
  expect_true(is.na(retained_cached$fallback_method))
  expect_true(is.na(retained_cached$fallback_max_backward_error))
  expect_true(all(rows$stage_native_iteration_seconds > 0))
  expect_true(all(rows$stage_golub_kahan_ritz_seconds > 0))
  expect_gte(cached$cached_start_attempts, 1L)
  expect_gte(cached_random$cached_start_attempts, 1L)
  expect_equal(retained$cached_start_attempts, 0L)
  expect_gte(retained_cached$cached_start_attempts, 1L)
  expect_lt(cached$matvecs, lean$matvecs)
  expect_lte(cached_random$matvecs, regular$matvecs)
  expect_lte(lean$matvecs, regular$matvecs * 2L)
})

test_that("SVD reference gate can evaluate an explicit H candidate subject", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore", "eigencore_block_golub_kahan_retained", "RSpectra"),
    median = c(3, 2, 1),
    mem_alloc = c(300, 200, 100),
    certificate_passed = c(TRUE, TRUE, TRUE),
    nconv = c(2L, 2L, 2L)
  )
  gate <- evaluate_reference_gate(
    rows[rows$method != "eigencore", , drop = FALSE],
    subject = "eigencore_block_golub_kahan_retained",
    references = "RSpectra",
    requested = 2L,
    speed_ratio_required = release_speed_gate("svd")
  )

  expect_equal(gate$subject, "eigencore_block_golub_kahan_retained")
  expect_false(gate$passed)
  expect_equal(gate$speed_ratio_vs_best_reference, 0.5)
  expect_equal(gate$memory_ratio_vs_best_reference, 0.5)
})

test_that("SVD benchmark rows audit raw and eigencore-certified reference timing", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  A <- diag(c(5, 4, 3, 2, 1))
  rows <- benchmark_svd_case(
    A,
    rank = 2L,
    methods = c("eigencore", "base"),
    iterations = 1L,
    tol = 1e-8,
    seed = 991L
  )

  expect_true(all(c(
    "raw_solver_median",
    "eigencore_certificate_median",
    "eigencore_certified_total_median",
    "certificate_recomputed_by_eigencore",
    "max_left_residual",
    "max_right_residual",
    "max_cyclic_residual",
    "orthogonality_U",
    "orthogonality_V",
    "singular_values_sorted"
  ) %in% names(rows)))
  expect_equal(
    rows$certificate_recomputed_by_eigencore[rows$method == "base"],
    TRUE
  )
  expect_equal(
    rows$certificate_recomputed_by_eigencore[rows$method == "eigencore"],
    FALSE
  )
  expect_true(all(rows$singular_values_sorted))
  expect_true(all(rows$certificate_passed))
  expect_true(all(is.finite(rows$max_left_residual)))
  expect_true(all(is.finite(rows$max_right_residual)))
})

test_that("tiny Gram eigensolver benchmark compares native backends", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")
  skip_if_not_installed("bench")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- benchmark_tiny_gram_eigensolvers(
    dimensions = 16L,
    ranks = c(3L, 5L),
    iterations = 1L,
    seed = 992L
  )

  expect_equal(
    sort(unique(rows$backend)),
    sort(c(
      "lapack_dsyevr_selected",
      "lapack_dsyevx_selected",
      "lapack_dsyev_full",
      "lapack_dsyevd_full"
    ))
  )
  expect_true(all(c(
    "dimension", "rank", "backend", "median", "mem_alloc",
    "max_value_error", "values_sorted", "winner"
  ) %in% names(rows)))
  expect_true(all(rows$values_sorted))
  expect_lte(max(rows$max_value_error), 1e-8)
  expect_true(all(tapply(rows$winner, list(rows$dimension, rows$rank), sum) >= 1L))
})

test_that("SVD memory diagnostics expose subject/reference allocation gaps", {
  skip_if(identical(Sys.getenv("CRAN"), "true"), "skip benchmark smoke on CRAN")

  helper_path <- system.file("benchmarks/_helpers.R", package = "eigencore")
  if (!nzchar(helper_path)) {
    helper_path <- test_path("../../inst/benchmarks/_helpers.R")
  }
  source(helper_path)

  rows <- data.frame(
    method = c("eigencore_block_golub_kahan_retained", "RSpectra", "irlba"),
    mem_alloc = c(200, 80, 100),
    solver_mem_alloc = c(180, 50, 70),
    certificate_mem_alloc = c(20, 30, 30),
    certificate_passed = c(TRUE, TRUE, TRUE),
    nconv = c(2L, 2L, 2L)
  )
  diagnostics <- evaluate_memory_diagnostics(
    rows,
    subject = "eigencore_block_golub_kahan_retained",
    references = c("RSpectra", "irlba"),
    requested = 2L
  )

  expect_equal(diagnostics$subject, "eigencore_block_golub_kahan_retained")
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

test_that("auto planner keeps sparse block promotion diagnostic-only by default", {
  mid <- Matrix::bandSparse(
    1000L,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, 999L), c(1, rep(2, 998L), 1), rep(-1, 999L))
  )
  mid[1, 1000] <- 0.01
  mid[1000, 1] <- 0.01
  mid_plan <- plan_solver(eigen_problem(mid, target = smallest()), k = 20L)
  expect_equal(mid_plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(mid_plan$controls$block, 1L)

  large <- Matrix::sparseMatrix(
    i = c(1:5000, 1, 5000),
    j = c(1:5000, 5000, 1),
    x = c(rep(1, 5000), 0.01, 0.01),
    dims = c(5000L, 5000L)
  )
  large_plan <- plan_solver(eigen_problem(large, target = smallest()), k = 20L)
  expect_equal(large_plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(large_plan$controls$block, 1L)

  small_k <- plan_solver(eigen_problem(mid, target = smallest()), k = 5L)
  expect_equal(small_k$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(small_k$controls$block, 1L)

  dense <- diag(as.numeric(seq(200L, 1L)))
  dense_plan <- plan_solver(eigen_problem(dense, target = largest()), k = 20L)
  expect_equal(dense_plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(dense_plan$controls$block, 2L)
  expect_equal(dense_plan$controls$max_subspace, 200L)
})

test_that("sparse block auto promotion remains available as diagnostic opt-in", {
  old <- options(eigencore.promote_sparse_block_lanczos = TRUE)
  on.exit(options(old), add = TRUE)

  mid <- Matrix::bandSparse(
    1000L,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, 999L), c(1, rep(2, 998L), 1), rep(-1, 999L))
  )
  mid[1, 1000] <- 0.01
  mid[1000, 1] <- 0.01
  mid_plan <- plan_solver(eigen_problem(mid, target = smallest()), k = 20L)
  expect_equal(mid_plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(mid_plan$controls$block, 2L)
  expect_equal(mid_plan$controls$max_subspace, 160L)

  large <- Matrix::sparseMatrix(
    i = c(1:5000, 1, 5000),
    j = c(1:5000, 5000, 1),
    x = c(rep(1, 5000), 0.01, 0.01),
    dims = c(5000L, 5000L)
  )
  large_plan <- plan_solver(eigen_problem(large, target = smallest()), k = 20L)
  expect_equal(large_plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(large_plan$controls$block, 4L)
  expect_equal(large_plan$controls$max_subspace, 320L)
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
  expect_equal(release_speed_gate("svd"), 1.1)
  expect_equal(release_speed_gate("randomized_svd"), 2.0)
})
