restart_counted_operator <- function(values, operator_id = "restart-counted",
                                     revision = "r1", structure = hermitian(),
                                     portable = TRUE) {
  counts <- new.env(parent = emptyenv())
  counts$calls <- 0L
  counts$columns <- 0L
  apply <- function(X, alpha = 1, beta = 0, Y = NULL) {
    X <- as.matrix(X)
    counts$calls <- counts$calls + 1L
    counts$columns <- counts$columns + ncol(X)
    out <- alpha * (values * X)
    if (is.null(Y) || beta == 0) out else out + beta * Y
  }
  list(
    operator = linear_operator(
      c(length(values), length(values)),
      apply,
      structure = structure,
      operator_id = operator_id,
      revision = revision,
      portable = portable,
      metadata = list(frobenius_norm = sqrt(sum(values^2)))
    ),
    counts = counts
  )
}

restart_lanczos_plan <- function(A, k = 3L, block = k,
                                 target = largest(), tol = 1e-8) {
  plan_solver(
    eigen_problem(A, target = target),
    k = k,
    method = lanczos(
      block = block,
      max_subspace = max(16L, 6L * k),
      max_restarts = 100L
    ),
    tol = tol
  )
}

repack_restart_state <- function(state) {
  state$serialization$integrity_token <- eigencore:::stable_raw_hash(
    eigencore:::restart_state_integrity_payload(state)
  )
  state$memory <- eigencore:::new_restart_memory(state)
  state
}

expect_restart_error <- function(expr, code) {
  error <- tryCatch(expr, error = identity)
  expect_s3_class(error, "eigencore_restart_state_error")
  expect_identical(error$code, code)
  error
}

test_that("restart states have a frozen schema, immutable copies, and exact memory", {
  plan <- restart_lanczos_plan(diag(seq(50, 1)))
  fit <- solve(plan, retain_state = "same_operator")
  state <- fit$restart_state

  expect_s3_class(state, "eigencore_restart_state")
  expect_named(state, eigencore:::required_restart_state_fields())
  expect_identical(state$schema_version, 1L)
  expect_equal(crossprod(state$basis), diag(ncol(state$basis)), tolerance = 1e-10)
  expect_s3_class(state$problem_signature, "eigencore_problem_signature")
  expect_s3_class(state$method_state, "eigencore_method_state")
  expect_true(state$provenance$method_state_available)
  expect_s3_class(state$memory, "eigencore_memory")
  expect_true(state$memory$complete)
  expect_identical(state$memory$native_bytes, 0)
  expect_named(
    state$memory$components,
    c("basis", "method_state", "cached_operator_actions", "metadata")
  )
  expect_identical(state$memory$components$cached_operator_actions, 0)

  bytes <- retained_bytes(state)
  expect_identical(as.numeric(bytes), state$memory$total_bytes)
  expect_identical(attr(bytes, "memory"), state$memory)
  expect_gt(retained_bytes(plan), 0)

  before <- serialize(state, NULL, version = 3L)
  basis_only <- restart_state(state, retention = "basis")
  expect_identical(serialize(state, NULL, version = 3L), before)
  expect_null(basis_only$method_state)
  expect_lt(retained_bytes(basis_only), retained_bytes(state))

  path <- tempfile(fileext = ".rds")
  saveRDS(state, path)
  restored <- readRDS(path)
  expect_identical(restart_state(restored, retention = "same_operator"), state)
})

test_that("same-operator Lanczos reuses only a fitted start and certifies freshly", {
  A <- diag(seq(60, 1))
  plan <- restart_lanczos_plan(A, k = 4L, block = 4L)
  first <- solve(plan, retain_state = "same_operator")
  state <- restart_state(first, retention = "same_operator")
  before <- serialize(state, NULL, version = 3L)

  second <- solve(
    plan,
    restart_state = state,
    reuse = "same_operator",
    retain_state = "same_operator"
  )

  expect_identical(serialize(state, NULL, version = 3L), before)
  expect_true(certificate(second)$passed)
  expect_equal(values(second), values(first), tolerance = 1e-8)
  expect_identical(second$state_transition$relation, "same_operator")
  expect_true(second$state_transition$basis_used)
  expect_true(second$state_transition$method_state_used)
  expect_identical(
    second$state_transition$adapter$kind,
    "hermitian_lanczos_start_block"
  )
  expect_true(all(c(
    "recurrence", "locked", "cached_operator_actions", "projection",
    "residuals", "convergence", "certificate"
  ) %in% second$state_transition$invalidated))
  expect_gt(
    work(second)$operator_columns +
      work(second)$certification_operator_columns,
    0L
  )
  expect_s3_class(second$restart_state, "eigencore_restart_state")
  expect_false(identical(second$restart_state, state))
  expect_named(
    second$state_transition,
    c(
      "schema_version", "relation", "basis_used", "method_state_used",
      "invalidated", "reason", "source_operator_identity",
      "destination_operator_identity", "reuse", "adapter"
    )
  )
})

test_that("changed revisions and lineages use only the public basis", {
  source <- restart_counted_operator(seq(70, 1), revision = "rho=0")
  first_plan <- restart_lanczos_plan(source$operator, k = 3L, block = 3L)
  first <- solve(first_plan, retain_state = "same_operator")

  changed <- restart_counted_operator(
    seq(70, 1) - 0.01 * seq_len(70),
    revision = "rho=0.01"
  )
  changed_plan <- restart_lanczos_plan(changed$operator, k = 3L, block = 3L)
  continued <- solve(
    changed_plan,
    restart_state = first$restart_state,
    reuse = "auto"
  )
  continued_work <- work(continued)

  expect_identical(continued$state_transition$relation, "changed_revision")
  expect_true(continued$state_transition$basis_used)
  expect_false(continued$state_transition$method_state_used)
  expect_true(all(c(
    "method_state", "recurrence", "locked", "cached_operator_actions",
    "projection", "residuals", "convergence", "certificate"
  ) %in% continued$state_transition$invalidated))
  expect_true(certificate(continued)$passed)
  expect_identical(
    continued_work$operator_block_calls +
      continued_work$certification_operator_block_calls,
    changed$counts$calls
  )
  expect_identical(
    continued_work$operator_columns +
      continued_work$certification_operator_columns,
    changed$counts$columns
  )

  other <- restart_counted_operator(
    seq(70, 1) - 0.02 * seq_len(70),
    operator_id = "different-lineage",
    revision = "r1"
  )
  other_fit <- solve(
    restart_lanczos_plan(other$operator, k = 3L, block = 3L),
    restart_state = first$restart_state,
    reuse = "basis_only"
  )
  expect_identical(other_fit$state_transition$relation, "changed_operator")
  expect_true(other_fit$state_transition$basis_used)
  expect_false(other_fit$state_transition$method_state_used)
  expect_true(certificate(other_fit)$passed)
})

test_that("stale method state downgrades only under auto and strict reuse fails closed", {
  counted <- restart_counted_operator(seq(50, 1))
  source_plan <- restart_lanczos_plan(counted$operator, block = 2L)
  first <- solve(source_plan, retain_state = "same_operator")
  state <- first$restart_state

  counted$counts$calls <- 0L
  counted$counts$columns <- 0L
  changed_controls <- restart_lanczos_plan(counted$operator, block = 3L)
  strict <- expect_restart_error(
    solve(changed_controls, restart_state = state, reuse = "same_operator"),
    "method_incompatible"
  )
  expect_identical(strict$field, "method_state$controls_token")
  expect_identical(counted$counts$calls, 0L)
  expect_identical(counted$counts$columns, 0L)

  automatic <- solve(changed_controls, restart_state = state, reuse = "auto")
  expect_true(automatic$state_transition$basis_used)
  expect_false(automatic$state_transition$method_state_used)
  expect_identical(
    automatic$state_transition$reason$code,
    "method_incompatible"
  )
  expect_true(certificate(automatic)$passed)

  future <- state
  future$method_state$adapter_version <- 2L
  future <- repack_restart_state(future)
  auto_future <- solve(source_plan, restart_state = future, reuse = "auto")
  expect_true(auto_future$state_transition$basis_used)
  expect_false(auto_future$state_transition$method_state_used)
  expect_identical(auto_future$state_transition$reason$code,
                   "stale_method_state")
  expect_restart_error(
    solve(source_plan, restart_state = future, reuse = "same_operator"),
    "stale_method_state"
  )
})

test_that("target changes invalidate payloads while preserving safe basis reuse", {
  A <- diag(seq(45, 1))
  largest_plan <- restart_lanczos_plan(A, target = largest())
  state <- solve(largest_plan, retain_state = "same_operator")$restart_state
  smallest_plan <- restart_lanczos_plan(A, target = smallest())

  automatic <- solve(
    smallest_plan, restart_state = state, reuse = "auto"
  )
  expect_true(automatic$state_transition$basis_used)
  expect_false(automatic$state_transition$method_state_used)
  expect_identical(automatic$state_transition$reason$code,
                   "method_incompatible")
  expect_true(certificate(automatic)$passed)
  expect_equal(sort(values(automatic)), 1:3, tolerance = 1e-7)

  expect_restart_error(
    solve(smallest_plan, restart_state = state, reuse = "same_operator"),
    "method_incompatible"
  )
})

test_that("corrupt, incompatible, and unknown states fail before operator work", {
  source <- restart_counted_operator(seq(40, 1))
  plan <- restart_lanczos_plan(source$operator)
  state <- solve(plan, retain_state = "same_operator")$restart_state

  source$counts$calls <- 0L
  source$counts$columns <- 0L
  future <- state
  future$schema_version <- 2L
  expect_restart_error(
    solve(plan, restart_state = future, reuse = "auto"),
    "unsupported_schema"
  )
  expect_identical(source$counts$calls, 0L)

  mutated <- state
  mutated$basis[1, 1] <- mutated$basis[1, 1] + 0.01
  expect_restart_error(
    solve(plan, restart_state = mutated, reuse = "basis_only"),
    "corrupt_state"
  )
  expect_identical(source$counts$calls, 0L)

  nonfinite <- state
  nonfinite$basis[1, 1] <- NA_real_
  expect_restart_error(
    solve(plan, restart_state = nonfinite, reuse = "auto"),
    "corrupt_state"
  )
  expect_identical(source$counts$calls, 0L)

  wrong_size <- restart_counted_operator(seq(41, 1), operator_id = "other")
  wrong_plan <- restart_lanczos_plan(wrong_size$operator)
  expect_restart_error(
    solve(wrong_plan, restart_state = state, reuse = "auto"),
    "coordinate_incompatible"
  )
  expect_identical(wrong_size$counts$calls, 0L)
})

test_that("unsupported receiving routes reject states rather than ignoring them", {
  counted <- restart_counted_operator(seq(36, 1))
  lanczos_plan <- restart_lanczos_plan(counted$operator)
  state <- solve(lanczos_plan, retain_state = "basis")$restart_state

  counted$counts$calls <- 0L
  lobpcg_plan <- plan_solver(
    eigen_problem(counted$operator),
    k = 3L,
    method = lobpcg(maxit = 30L),
    tol = 1e-7
  )
  expect_restart_error(
    solve(lobpcg_plan, restart_state = state),
    "method_incompatible"
  )
  expect_identical(counted$counts$calls, 0L)

  fallback_plan <- plan_solver(
    eigen_problem(counted$operator),
    k = 3L,
    method = auto(),
    tol = 1e-7
  )
  if (!eigencore:::warm_start_plan_consumes_start(
    fallback_plan$problem, fallback_plan
  )) {
    expect_restart_error(
      solve(fallback_plan, restart_state = state),
      "method_incompatible"
    )
    expect_identical(counted$counts$calls, 0L)
  }

  A <- diag(seq(36, 1))
  dense_state <- solve(
    restart_lanczos_plan(A), retain_state = "basis"
  )$restart_state
  dense_fallback <- plan_solver(
    eigen_problem(A), k = 3L, method = auto(), tol = 1e-7
  )
  if (!eigencore:::warm_start_plan_consumes_start(
    dense_fallback$problem, dense_fallback
  )) {
    expect_restart_error(
      solve(dense_fallback, restart_state = dense_state),
      "method_incompatible"
    )
  }

  shift_plan <- plan_solver(
    eigen_problem(
      A,
      target = nearest(10.5),
      transform = shift_invert(10.5)
    ),
    k = 2L,
    method = shift_invert(10.5),
    tol = 1e-7
  )
  expect_restart_error(
    solve(shift_plan, restart_state = dense_state),
    "method_incompatible"
  )

  generalized_plan <- plan_solver(
    eigen_problem(A, metric = diag(seq(1, 2, length.out = nrow(A)))),
    k = 3L,
    method = lanczos(block = 3L),
    tol = 1e-7
  )
  expect_restart_error(
    solve(generalized_plan, restart_state = dense_state),
    "coordinate_incompatible"
  )
})

test_that("generalized states preserve the certified B-orthonormal basis", {
  A <- diag(c(12, 9, 5, 2, 1))
  B <- diag(c(1, 1.5, 2, 2.5, 3))
  plan <- plan_solver(
    eigen_problem(A, metric = B, target = smallest()),
    k = 2L,
    method = lobpcg(maxit = 100L),
    tol = 1e-8
  )
  fit <- solve(plan, retain_state = "same_operator")
  state <- fit$restart_state

  expect_identical(state$problem_signature$basis_metric, "B")
  expect_equal(
    crossprod(state$basis, B %*% state$basis),
    diag(2L),
    tolerance = 1e-8
  )
  expect_null(state$method_state)
  expect_false(state$provenance$method_state_available)
  expect_restart_error(
    solve(plan, restart_state = state, reuse = "basis_only"),
    "method_incompatible"
  )
})

test_that("dense, CSC, explicit, and opaque Lanczos bases follow admitted paths", {
  skip_if_not_installed("Matrix")
  n <- 70L
  diagonal <- seq(n, 1)
  A1 <- methods::as(
    methods::as(Matrix::Diagonal(x = diagonal), "generalMatrix"),
    "CsparseMatrix"
  )
  A2 <- methods::as(
    methods::as(
      Matrix::Diagonal(x = diagonal - 0.005 * seq_len(n)),
      "generalMatrix"
    ),
    "CsparseMatrix"
  )
  first <- solve(
    restart_lanczos_plan(A1), retain_state = "same_operator"
  )
  second <- solve(
    restart_lanczos_plan(A2),
    restart_state = first$restart_state,
    reuse = "auto"
  )
  expect_identical(second$state_transition$relation, "changed_operator")
  expect_true(second$state_transition$basis_used)
  expect_true(certificate(second)$passed)

  opaque <- linear_operator(
    c(n, n),
    function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (diagonal * X)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    structure = hermitian(),
    metadata = list(frobenius_norm = sqrt(sum(diagonal^2)))
  )
  opaque_plan <- restart_lanczos_plan(opaque)
  opaque_state <- solve(
    opaque_plan, retain_state = "same_operator"
  )$restart_state
  expect_false(opaque_state$serialization$portable)
  expect_false(opaque_state$method_state$portable)
  automatic <- solve(
    opaque_plan, restart_state = opaque_state, reuse = "auto"
  )
  expect_true(automatic$state_transition$basis_used)
  expect_false(automatic$state_transition$method_state_used)
  expect_identical(automatic$state_transition$reason$code,
                   "session_incompatible")
  expect_true(certificate(automatic)$passed)
  expect_restart_error(
    solve(opaque_plan, restart_state = opaque_state, reuse = "same_operator"),
    "session_incompatible"
  )
})

test_that("general Arnoldi states are inspectable but rejected before callbacks", {
  n <- 24L
  A <- diag(seq(n, 1))
  A[upper.tri(A)] <- 0.02
  calls <- 0L
  op <- linear_operator(
    c(n, n),
    function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      out <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      out <- alpha * (t(A) %*% X)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    structure = general(),
    operator_id = "restart-arnoldi",
    revision = "r1",
    portable = TRUE,
    metadata = list(frobenius_norm = sqrt(sum(A^2)))
  )
  plan <- plan_solver(eigen_problem(op), k = 2L, tol = 1e-7)
  fit <- solve(plan)
  state <- restart_state(fit)
  calls <- 0L

  expect_restart_error(
    solve(plan, restart_state = state, reuse = "basis_only"),
    "method_incompatible"
  )
  expect_identical(calls, 0L)
})

test_that("SVD states preserve original sides but have no receiving adapter", {
  set.seed(91)
  A <- matrix(rnorm(240), 20, 12)
  plan <- plan_solver(
    svd_problem(A), rank = 3L, method = golub_kahan(max_subspace = 12L),
    tol = 1e-7, vectors = "both"
  )
  fit <- solve(plan)
  state <- restart_state(fit)

  expect_s3_class(state$basis, "eigencore_svd_basis")
  expect_identical(state$problem_signature$basis_sides, c("left", "right"))
  expect_equal(crossprod(state$basis$left), diag(3), tolerance = 1e-7)
  expect_equal(crossprod(state$basis$right), diag(3), tolerance = 1e-7)
  expect_null(state$method_state)
  expect_restart_error(
    solve(plan, restart_state = state, reuse = "basis_only"),
    "method_incompatible"
  )
})

test_that("state extraction requires a passed certificate and available vectors", {
  plan <- restart_lanczos_plan(diag(seq(20, 1)), k = 2L, block = 2L)
  fit <- solve(plan)
  failed <- fit
  failed$certificate$passed <- FALSE
  expect_restart_error(restart_state(failed), "corrupt_state")

  missing <- fit
  missing$vectors <- NULL
  expect_restart_error(restart_state(missing), "corrupt_state")
})
