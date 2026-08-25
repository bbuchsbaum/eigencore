expected_work_fields <- function() {
  c(
    "schema_version", "complete",
    "operator_block_calls", "operator_columns",
    "adjoint_block_calls", "adjoint_columns",
    "metric_block_calls", "metric_columns",
    "preconditioner_calls", "preconditioner_columns",
    "certification_operator_block_calls",
    "certification_operator_columns",
    "certification_adjoint_block_calls",
    "certification_adjoint_columns",
    "iterations", "restarts", "setup_seconds", "solve_seconds",
    "certification_seconds", "total_seconds", "legacy_matvecs"
  )
}

counted_operator <- function(A, structure = hermitian(), with_adjoint = FALSE,
                             operator_id = "counted-operator") {
  counts <- new.env(parent = emptyenv())
  counts$forward_calls <- 0L
  counts$forward_columns <- 0L
  counts$adjoint_calls <- 0L
  counts$adjoint_columns <- 0L
  forward <- function(X, alpha = 1, beta = 0, Y = NULL) {
    X <- as.matrix(X)
    counts$forward_calls <- counts$forward_calls + 1L
    counts$forward_columns <- counts$forward_columns + ncol(X)
    out <- alpha * (A %*% X)
    if (is.null(Y) || beta == 0) out else out + beta * Y
  }
  adjoint_apply <- if (isTRUE(with_adjoint)) {
    function(X, alpha = 1, beta = 0, Y = NULL) {
      X <- as.matrix(X)
      counts$adjoint_calls <- counts$adjoint_calls + 1L
      counts$adjoint_columns <- counts$adjoint_columns + ncol(X)
      out <- alpha * (t(A) %*% X)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    }
  } else {
    NULL
  }
  list(
    operator = linear_operator(
      dim(A), forward, adjoint_apply,
      structure = structure,
      operator_id = operator_id,
      revision = "r1",
      metadata = list(frobenius_norm = sqrt(sum(A^2)))
    ),
    counts = counts
  )
}

test_that("typed work schema is stable and preserves legacy matvecs", {
  fit <- eig_partial(diag(c(5, 3, 1)), k = 2L)
  record <- work(fit)

  expect_s3_class(record, "eigencore_work")
  expect_named(record, expected_work_fields())
  expect_identical(record$schema_version, 1L)
  expect_identical(record$legacy_matvecs, fit$matvecs)
  expect_true(record$complete)
  expect_true(all(vapply(
    record[eigencore:::work_counter_fields()],
    function(x) is.integer(x) && length(x) == 1L && x >= 0L,
    logical(1L)
  )))
  expect_gte(record$total_seconds, record$setup_seconds)
  expect_identical(work(record), record)
  expect_identical(diagnostics(fit)$work, record)

  legacy <- structure(
    list(matvecs = 7L, iterations = 5L, requested = 2L),
    class = "eigencore_eigen_result"
  )
  legacy_work <- work(legacy)
  expect_s3_class(legacy_work, "eigencore_work")
  expect_false(legacy_work$complete)
  expect_identical(legacy_work$legacy_matvecs, 7L)
  expect_true(anyNA(unlist(
    legacy_work[eigencore:::work_counter_fields()],
    use.names = FALSE
  )))
})

test_that("Lanczos callback work separates solve and certification exactly", {
  set.seed(11)
  counted <- counted_operator(diag(seq(60, 1)), operator_id = "work-lanczos")
  fit <- eig_partial(
    counted$operator,
    k = 3L,
    method = lanczos(block = 2L, max_subspace = 30L),
    tol = 1e-8
  )
  record <- work(fit)

  expect_identical(
    record$operator_block_calls + record$certification_operator_block_calls,
    counted$counts$forward_calls
  )
  expect_identical(
    record$operator_columns + record$certification_operator_columns,
    counted$counts$forward_columns
  )
  expect_identical(record$adjoint_block_calls, 0L)
  expect_identical(record$metric_block_calls, 0L)
  expect_true(record$complete)
})

test_that("built-in dense warm starts retain native work diagnostics", {
  A0 <- diag(seq(40, 1))
  cold <- eig_partial(
    A0,
    k = 3L,
    method = lanczos(block = 3L, max_subspace = 20L),
    tol = 1e-8
  )
  warm <- eig_partial(
    A0 - 0.1 * diag(seq_len(nrow(A0))),
    k = 3L,
    method = lanczos(block = 3L, max_subspace = 20L),
    tol = 1e-8,
    initial_subspace = vectors(cold)
  )
  record <- work(warm)

  expect_true(record$complete)
  expect_identical(record$certification_operator_block_calls, 0L)
  expect_gt(record$operator_block_calls, 1L)
  expect_gt(record$operator_columns, ncol(vectors(cold)))
})

test_that("Golub-Kahan callback work keeps forward and adjoint sides distinct", {
  set.seed(12)
  counted <- counted_operator(
    matrix(rnorm(800), 40, 20),
    structure = general(),
    with_adjoint = TRUE,
    operator_id = "work-golub-kahan"
  )
  fit <- svd_partial(
    counted$operator,
    rank = 3L,
    method = golub_kahan(max_subspace = 20L),
    tol = 1e-7
  )
  record <- work(fit)

  expect_identical(
    record$operator_block_calls + record$certification_operator_block_calls,
    counted$counts$forward_calls
  )
  expect_identical(
    record$adjoint_block_calls + record$certification_adjoint_block_calls,
    counted$counts$adjoint_calls
  )
  expect_identical(
    record$operator_columns + record$certification_operator_columns,
    counted$counts$forward_columns
  )
  expect_identical(
    record$adjoint_columns + record$certification_adjoint_columns,
    counted$counts$adjoint_columns
  )
  expect_true(record$complete)
})

test_that("Arnoldi callback work includes independently observed left solves", {
  set.seed(13)
  counted <- counted_operator(
    matrix(rnorm(625), 25, 25),
    structure = general(),
    with_adjoint = TRUE,
    operator_id = "work-arnoldi"
  )
  fit <- eig_partial(counted$operator, k = 2L, tol = 1e-7)
  record <- work(fit)

  expect_identical(
    record$operator_block_calls + record$certification_operator_block_calls,
    counted$counts$forward_calls
  )
  expect_identical(
    record$adjoint_block_calls + record$certification_adjoint_block_calls,
    counted$counts$adjoint_calls
  )
  expect_identical(
    record$operator_columns + record$certification_operator_columns,
    counted$counts$forward_columns
  )
  expect_identical(
    record$adjoint_columns + record$certification_adjoint_columns,
    counted$counts$adjoint_columns
  )
  expect_true(record$complete)
})

test_that("metric and preconditioner work agrees with independent counters", {
  set.seed(14)
  n <- 20L
  counted_A <- counted_operator(diag(seq(n, 1)), operator_id = "work-metric-A")
  counted_B <- counted_operator(
    diag(seq(1, 2, length.out = n)),
    operator_id = "work-metric-B"
  )
  counted_B$operator$metadata$symmetric_positive_definite <- TRUE

  preconditioner_calls <- 0L
  preconditioner_columns <- 0L
  preconditioner <- function(R) {
    preconditioner_calls <<- preconditioner_calls + 1L
    preconditioner_columns <<- preconditioner_columns + ncol(R)
    R
  }
  fit <- eig_partial(
    counted_A$operator,
    B = counted_B$operator,
    k = 2L,
    target = smallest(),
    method = lobpcg(maxit = 30L, preconditioner = preconditioner),
    tol = 1e-7
  )
  record <- work(fit)

  expect_identical(
    record$operator_block_calls + record$certification_operator_block_calls,
    counted_A$counts$forward_calls
  )
  expect_identical(record$metric_block_calls, counted_B$counts$forward_calls)
  expect_identical(record$metric_columns, counted_B$counts$forward_columns)
  expect_identical(record$preconditioner_calls, preconditioner_calls)
  expect_identical(record$preconditioner_columns, preconditioner_columns)
  expect_true(record$complete)
})
