test_that("eigen and SVD plans are executable, frozen, and immutable", {
  eigen_plan <- plan_solver(
    eigen_problem(diag(c(5, 4, 3, 2, 1))),
    k = 2L,
    tol = 1e-10,
    vectors = TRUE,
    certify = TRUE,
    allow_dense_fallback = "auto"
  )
  svd_plan <- plan_solver(
    svd_problem(matrix(seq_len(30), 6, 5)),
    rank = 2L,
    tol = 1e-9,
    vectors = "both",
    certify = TRUE,
    allow_dense_fallback = "auto"
  )

  required <- eigencore:::required_plan_fields()
  expect_true(all(required %in% names(eigen_plan)))
  expect_true(all(required %in% names(svd_plan)))
  expect_identical(eigen_plan$schema_version, 1L)
  expect_identical(svd_plan$schema_version, 1L)
  expect_s3_class(eigen_plan$planner_policy, "eigencore_planner_policy")
  expect_s3_class(eigen_plan$memory, "eigencore_memory")
  expect_identical(names(eigen_plan$planner_policy),
                   eigencore:::planner_policy_keys())

  eigen_before <- serialize(eigen_plan, NULL, version = 3L)
  svd_before <- serialize(svd_plan, NULL, version = 3L)
  eigen_fit <- solve(eigen_plan)
  svd_fit <- solve(svd_plan)

  expect_identical(serialize(eigen_plan, NULL, version = 3L), eigen_before)
  expect_identical(serialize(svd_plan, NULL, version = 3L), svd_before)
  expect_identical(eigen_fit$plan, eigen_plan)
  expect_identical(svd_fit$plan, svd_plan)
  expect_equal(values(eigen_fit), c(5, 4), tolerance = 1e-12)
  expect_equal(values(svd_fit), svd(matrix(seq_len(30), 6, 5), nu = 2, nv = 2)$d[1:2],
               tolerance = 1e-10)
  expect_identical(eigen_fit$planned_method, eigen_plan$planned_method)
  expect_identical(eigen_fit$actual_method, eigen_fit$method)
  expect_false(eigen_fit$fallback_used)
  expect_null(eigen_fit$fallback_reason)
})

test_that("solve(plan) reads frozen policy rather than changed global options", {
  policy_keys <- eigencore:::planner_policy_keys()
  old <- setNames(lapply(policy_keys, getOption), policy_keys)
  options(
    eigencore.dense_partial_lanczos_min_n = 128L,
    eigencore.dense_partial_lanczos_max_fraction = 0.25,
    eigencore.block_dense_full_subspace_max_n = 256L,
    eigencore.dense_fallback_mb = 256,
    eigencore.randomized_stage_timing = FALSE,
    eigencore.randomized_adaptive_stop = TRUE,
    eigencore.golub_kahan_projected_stop = FALSE,
    eigencore.golub_kahan_prefix_diagnostics = FALSE
  )

  A <- diag(seq(160, 1))
  eigen_plan <- plan_solver(eigen_problem(A), k = 8L)
  expect_identical(eigen_plan$method,
                   "native scalar thick-restart Hermitian Lanczos")

  X <- matrix(seq_len(48), 8, 6)
  dense_plan <- plan_solver(svd_problem(X), rank = 2L)
  set.seed(81)
  dense_before <- solve(dense_plan)

  options(
    eigencore.arnoldi_max_restarts = 0L,
    eigencore.block_dense_full_subspace_max_n = 1L,
    eigencore.dense_partial_lanczos_max_fraction = 1e-6,
    eigencore.dense_partial_lanczos_min_n = 10000L,
    eigencore.promote_sparse_block_lanczos = TRUE,
    eigencore.lobpcg_maxit = 1L,
    eigencore.dense_fallback_mb = 0,
    eigencore.max_dense_fallback_bytes = 0,
    eigencore.gram_svd_max_dimension = 1,
    eigencore.gram_svd_max_dimension_wide = 1,
    eigencore.gram_svd_memory_mb = 0,
    eigencore.gram_svd_rank_fraction_limit = .Machine$double.eps,
    eigencore.gram_svd_min_aspect_ratio = 1e6,
    eigencore.gram_svd_work_budget = 0,
    eigencore.implicit_gram_svd_min_dimension = 10000L,
    eigencore.promote_retained_golub_kahan = TRUE,
    eigencore.golub_kahan_projected_stop = TRUE,
    eigencore.golub_kahan_prefix_diagnostics = TRUE,
    eigencore.randomized_stage_timing = TRUE,
    eigencore.randomized_adaptive_stop = FALSE
  )

  eigen_fit <- solve(eigen_plan)
  set.seed(81)
  dense_after <- solve(dense_plan)
  options(old)
  expect_identical(eigen_fit$planned_method, eigen_plan$planned_method)
  expect_identical(eigen_fit$actual_method, eigen_plan$planned_method)
  expect_equal(values(dense_after), values(dense_before), tolerance = 0)
  expect_identical(dense_after$plan$planner_policy,
                   dense_plan$planner_policy)
})

test_that("replan is explicit and leaves the original plan unchanged", {
  old <- options(
    eigencore.dense_partial_lanczos_min_n = 128L,
    eigencore.dense_fallback_mb = 256
  )
  A <- diag(seq(200, 1))
  plan <- plan_solver(eigen_problem(A), k = 10L)
  before <- serialize(plan, NULL, version = 3L)

  options(eigencore.dense_partial_lanczos_min_n = 1000L)
  frozen <- solve(plan)
  replanned <- solve(plan, replan = TRUE)
  options(old)

  expect_identical(frozen$method,
                   "native scalar thick-restart Hermitian Lanczos")
  expect_identical(replanned$method,
                   "native dense Hermitian LAPACK fallback")
  expect_false(identical(replanned$plan, plan))
  expect_identical(serialize(plan, NULL, version = 3L), before)
})

test_that("execution overrides and corrupt plans fail closed with typed codes", {
  plan <- plan_solver(eigen_problem(diag(c(4, 3, 2, 1))), k = 2L)

  override <- tryCatch(solve(plan, tol = 1e-4), error = identity)
  expect_s3_class(override, "eigencore_plan_error")
  expect_identical(override$code, "execution_override")
  expect_identical(override$field, "tol")

  legacy <- structure(
    list(
      problem_type = "eigen", requested = 2L, method = plan$method,
      target = "largest", reasons = character(), fallback = "none",
      controls = list()
    ),
    class = "eigencore_plan"
  )
  legacy_error <- tryCatch(solve(legacy), error = identity)
  expect_s3_class(legacy_error, "eigencore_plan_error")
  expect_identical(legacy_error$code, "unsupported_schema")

  missing <- plan
  missing$controls <- NULL
  missing_error <- tryCatch(solve(missing), error = identity)
  expect_identical(missing_error$code, "missing_field")
  expect_identical(missing_error$field, "controls")

  corrupt <- plan
  corrupt$controls$invented <- TRUE
  corrupt_error <- tryCatch(solve(corrupt), error = identity)
  expect_identical(corrupt_error$code, "corrupt_plan")
  expect_identical(corrupt_error$field, "controls")

  unavailable <- plan
  unavailable$planned_method <- unavailable$method <- "invented solver route"
  unavailable$serialization$dispatch_token <- eigencore:::stable_raw_hash(
    list(1L, unavailable$planned_method)
  )
  unavailable_error <- tryCatch(solve(unavailable), error = identity)
  expect_identical(unavailable_error$code, "dispatch_unavailable")
})

test_that("operator identity is deterministic, revisioned, and validated before apply", {
  A <- diag(c(5, 3, 1))
  same_a <- operator_identity(as_operator(A))
  same_b <- operator_identity(as_operator(A))
  changed <- operator_identity(as_operator(A + diag(c(0, 0, 0.25))))
  expect_identical(same_a, same_b)
  expect_false(identical(same_a$revision, changed$revision))
  expect_true(same_a$portable)

  calls <- 0L
  callback <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      alpha * (A %*% X)
    },
    structure = hermitian(),
    operator_id = "identity-test",
    revision = "r1",
    portable = TRUE,
    metadata = list(frobenius_norm = sqrt(sum(A^2)))
  )
  explicit <- operator_identity(callback)
  expect_identical(explicit$operator_id, "identity-test")
  expect_identical(explicit$revision, "r1")
  expect_true(explicit$portable)

  opaque_a <- linear_operator(dim(A), function(X, ...) A %*% X,
                              structure = hermitian())
  opaque_b <- linear_operator(dim(A), function(X, ...) A %*% X,
                              structure = hermitian())
  expect_false(operator_identity(opaque_a)$portable)
  expect_false(identical(operator_identity(opaque_a)$operator_id,
                         operator_identity(opaque_b)$operator_id))
  expect_error(
    linear_operator(dim(A), function(X, ...) A %*% X,
                    structure = hermitian(), portable = TRUE),
    "requires explicit operator_id and revision"
  )

  plan <- plan_solver(
    eigen_problem(callback),
    k = 2L,
    method = lanczos(block = 2L, max_subspace = 3L)
  )
  plan$problem$A$identity$session_id <- "different-session"
  plan$problem$A$identity$portable <- FALSE
  plan$operator_identity <- operator_identity(plan$problem)
  plan$serialization$operator_identity_token <- eigencore:::stable_raw_hash(
    plan$operator_identity
  )
  session_error <- tryCatch(solve(plan), error = identity)
  expect_identical(session_error$code, "session_incompatible")
  expect_identical(calls, 0L)
})

test_that("matrix-backed plans survive RDS round trips", {
  old <- options(eigencore.dense_fallback_mb = 256)
  plan <- plan_solver(
    eigen_problem(diag(c(7, 5, 3, 1))),
    k = 2L,
    tol = 1e-10
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(plan, path)
  restored <- readRDS(path)
  restored_values <- values(solve(restored))
  options(old)
  expect_true(restored$serialization$portable)
  expect_equal(restored_values, c(7, 5), tolerance = 1e-12)
})

test_that("certificate fallback preserves planned and actual provenance", {
  old <- options(
    eigencore.promote_retained_golub_kahan = FALSE,
    eigencore.gram_svd_max_dimension = 512,
    eigencore.gram_svd_max_dimension_wide = 1024,
    eigencore.gram_svd_memory_mb = 64,
    eigencore.gram_svd_rank_fraction_limit = 0.5,
    eigencore.gram_svd_min_aspect_ratio = 2,
    eigencore.gram_svd_work_budget = Inf
  )
  set.seed(2)
  U <- Matrix::rsparsematrix(120, 2, density = 0.2)
  V <- Matrix::rsparsematrix(2, 30, density = 0.2)
  A <- U %*% Matrix::Diagonal(x = c(1, 1e-6)) %*% V
  plan <- plan_solver(svd_problem(A), rank = 2L, tol = 1e-10)

  fit <- solve(plan)
  options(old)
  expect_identical(plan$planned_method,
                   "native certified Gram SVD special case")
  expect_identical(fit$planned_method, plan$planned_method)
  expect_identical(fit$actual_method, fit$method)
  expect_true(fit$fallback_used)
  expect_false(identical(fit$actual_method, fit$planned_method))
  expect_s3_class(fit$fallback_reason, "eigencore_fallback_reason")
  expect_identical(fit$fallback_reason$code, "certificate_failed")
  expect_identical(fit$fallback_reason$planned_method, fit$planned_method)
  expect_identical(fit$fallback_reason$attempted_method, fit$actual_method)
  expect_true(certificate(fit)$passed)
  expect_identical(fit$plan, plan)
})
