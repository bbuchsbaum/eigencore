test_that("auto promotes structured generalized SPD LOBPCG without dense fallback", {
  n <- 24L
  A <- Matrix::Diagonal(x = seq_len(n) + 1)
  B <- Matrix::Diagonal(x = seq(1, 2, length.out = n))

  fit <- eig_partial(
    A,
    B = B,
    k = 3,
    target = smallest(),
    seed = 311,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(3), tolerance = 1e-8)
})

test_that("matrix-free B remains a certified reference fallback", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- diag(c(1, 2, 3, 4, 5, 6))
  Bop <- linear_operator(
    dim = dim(B),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    metadata = list(frobenius_norm = sqrt(sum(B^2)))
  )

  plan <- plan_solver(
    eigen_problem(A, metric = Bop, target = smallest()),
    k = 2,
    method = lobpcg(maxit = 80L)
  )
  expect_equal(plan$method, "reference generalized SPD LOBPCG prototype")

  auto_plan <- plan_solver(
    eigen_problem(A, metric = Bop, target = smallest()),
    k = 2
  )
  expect_equal(auto_plan$method, eigencore:::reference_generalized_lobpcg_label())

  auto_fit <- eig_partial(
    A,
    B = Bop,
    k = 2,
    target = smallest(),
    seed = 312,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )
  expect_equal(auto_fit$method, eigencore:::reference_generalized_lobpcg_label())
  expect_true(auto_fit$restart$generalized)
  expect_false(auto_fit$restart$native)
  expect_true(certificate(auto_fit)$passed)

  fit <- eig_partial(
    A,
    B = Bop,
    k = 2,
    target = smallest(),
    method = lobpcg(maxit = 80L),
    seed = 312,
    tol = 1e-8
  )

  expect_equal(fit$method, "reference generalized SPD LOBPCG prototype")
  expect_true(fit$restart$generalized)
  expect_false(fit$restart$native)
  expect_true(certificate(fit)$passed)
})

test_that("explicit SPD matrix-free B runs through native generalized LOBPCG", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- diag(c(1, 2, 3, 4, 5, 6))
  Bop <- linear_operator(
    dim = dim(B),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    metadata = list(
      frobenius_norm = sqrt(sum(B^2)),
      positive_definite = TRUE
    )
  )

  auto_plan <- plan_solver(
    eigen_problem(A, metric = Bop, target = smallest()),
    k = 2
  )
  expect_equal(auto_plan$method, eigencore:::native_generalized_lobpcg_label())

  fit <- eig_partial(
    A,
    B = Bop,
    k = 2,
    target = smallest(),
    method = lobpcg(maxit = 80L),
    seed = 316,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_equal(fit$restart$orthogonalization$methods, "native_matrix_free_b_mgs2")
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
})

test_that("non-SPD explicit B is not promoted", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16))
  B <- Matrix::Diagonal(x = c(1, 2, 0, 4))
  problem <- eigen_problem(A, metric = B, target = smallest())
  plan <- plan_solver(problem, k = 2, method = lobpcg(maxit = 50L))

  expect_equal(plan$method, "reference generalized SPD LOBPCG prototype")
  expect_error(
    eig_partial(A, B = B, k = 2, target = smallest(), method = lobpcg(maxit = 50L)),
    "positive definite"
  )
})

test_that("generalized preconditioners run through honest reference fallback", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))
  preconditioner <- shifted_cholesky_preconditioner(A, shift = 1e-3)
  method <- lobpcg(maxit = 80L, preconditioner = preconditioner)
  problem <- eigen_problem(A, metric = B, target = smallest())
  plan <- plan_solver(problem, k = 2, method = method)

  expect_equal(plan$method, "reference generalized SPD LOBPCG prototype")
  expect_match(paste(plan$reasons, collapse = "\n"), "preconditioner: shifted_cholesky")

  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = method,
    seed = 315,
    tol = 1e-8
  )

  expect_equal(fit$method, "reference generalized SPD LOBPCG prototype")
  expect_true(fit$restart$generalized)
  expect_false(fit$restart$native)
  expect_false(fit$restart$preconditioner_native)
  expect_equal(fit$restart$preconditioner_kind, "shifted_cholesky")
  expect_gt(fit$restart$preconditioner_calls, 0L)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)
})

test_that("native generalized LOBPCG accepts native shifted-tridiagonal preconditioners", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))
  preconditioner <- shifted_tridiagonal_preconditioner(A, shift = 1e-3)
  method <- lobpcg(maxit = 80L, preconditioner = preconditioner)
  problem <- eigen_problem(A, metric = B, target = smallest())
  plan <- plan_solver(problem, k = 2, method = method)

  expect_equal(plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_match(paste(plan$reasons, collapse = "\n"), "preconditioner: shifted_tridiagonal")

  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = method,
    seed = 317,
    tol = 1e-8
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_true(fit$restart$preconditioner_native)
  expect_equal(fit$restart$preconditioner_kind, "shifted_tridiagonal")
  expect_gt(fit$restart$preconditioner_calls, 0L)
  expect_true(certificate(fit)$passed)
})

test_that("native generalized LOBPCG accepts native shifted-diagonal preconditioners", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- diag(c(1, 2, 3, 4, 5, 6))
  Bop <- linear_operator(
    dim = dim(B),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (B %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    metadata = list(
      frobenius_norm = sqrt(sum(B^2)),
      positive_definite = TRUE
    )
  )
  preconditioner <- shifted_diagonal_preconditioner(A, shift = 1e-3)
  info <- eigencore:::eigencore_preconditioner_info(preconditioner, include_arrays = FALSE)
  method <- lobpcg(maxit = 80L, preconditioner = preconditioner)
  problem <- eigen_problem(A, metric = Bop, target = smallest())
  plan <- plan_solver(problem, k = 2, method = method)

  expect_equal(info$kind, "shifted_diagonal")
  expect_true(info$native)
  expect_equal(plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_match(paste(plan$reasons, collapse = "\n"), "preconditioner: shifted_diagonal")

  fit <- eig_partial(
    A,
    B = Bop,
    k = 2,
    target = smallest(),
    method = method,
    seed = 318,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_true(fit$restart$preconditioner_native)
  expect_equal(fit$restart$preconditioner_kind, "shifted_diagonal")
  expect_gt(fit$restart$preconditioner_calls, 0L)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
})

test_that("generalized constraints deflate known nullspace on native path", {
  A <- Matrix::Diagonal(x = c(0, 1, 4, 9))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4))
  nullspace <- matrix(c(1, 0, 0, 0), ncol = 1)
  method <- lobpcg(maxit = 80L, constraints = nullspace)
  problem <- eigen_problem(A, metric = B, target = smallest())
  plan <- plan_solver(problem, k = 1, method = method)

  expect_equal(plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_match(paste(plan$reasons, collapse = "\n"), "constraints: deflating 1 vector")

  fit <- eig_partial(
    A,
    B = B,
    k = 1,
    target = smallest(),
    method = method,
    seed = 313,
    tol = 1e-8
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$constrained)
  expect_true(fit$restart$native)
  expect_equal(fit$restart$constraints_rank, 1L)
  expect_equal(values(fit), 0.5, tolerance = 1e-8)
  expect_equal(crossprod(nullspace, as.matrix(B) %*% vectors(fit)), matrix(0, 1, 1), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
})

test_that("standard constraints deflate graph Laplacian nullspace", {
  n <- 10L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1))
  )
  nullspace <- matrix(rep(1, n), ncol = 1)
  fit <- eig_partial(
    A,
    k = 1,
    target = smallest(),
    method = lobpcg(maxit = 100L, constraints = nullspace),
    seed = 314,
    tol = 1e-8
  )
  oracle <- eigen(as.matrix(A), symmetric = TRUE)

  expect_equal(fit$method, "native standard Hermitian LOBPCG prototype")
  expect_true(fit$restart$native)
  expect_true(fit$restart$constrained)
  expect_equal(values(fit), sort(oracle$values)[2L], tolerance = 1e-8)
  expect_lt(abs(sum(vectors(fit))), 1e-8)
  expect_true(certificate(fit)$passed)
})
