test_that("eig_partial returns certified eigenpairs", {
  A <- diag(c(4, 1, 3, 2))
  fit <- eig_partial(A, k = 2, target = largest())

  expect_equal(values(fit), c(4, 3))
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$plan$method, "native dense Hermitian LAPACK fallback")
  expect_equal(fit$method, "native dense Hermitian LAPACK fallback")
})

test_that("native dense Hermitian fallback matches base eigen ordering", {
  set.seed(24)
  A0 <- matrix(rnorm(36), nrow = 6)
  A <- crossprod(A0)
  fit <- eig_partial(A, k = 3, target = smallest())
  oracle <- eigen(A, symmetric = TRUE)

  expect_equal(values(fit), sort(oracle$values)[1:3], tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_match(fit$warnings, "native dense Hermitian LAPACK fallback")
})

test_that("auto routes large dense partial Hermitian problems to native Lanczos", {
  A <- diag(c(10, 8, 6, rep(1, 157)))
  fit <- eig_partial(A, k = 3, target = largest(), seed = 27)

  expect_equal(fit$plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(fit$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(values(fit), c(10, 8, 6), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("auto keeps dense LAPACK for near-full dense Hermitian requests", {
  A <- diag(seq_len(160))
  plan <- plan_solver(eigen_problem(A), k = 80)

  expect_equal(plan$method, "native dense Hermitian LAPACK fallback")
})

test_that("solve dispatch helpers follow plan method labels", {
  expect_true(eigencore:::plan_dispatches_lanczos(list(
    method = "native scalar thick-restart Hermitian Lanczos"
  )))
  expect_true(eigencore:::plan_dispatches_native_lanczos(list(
    method = "native block Hermitian Lanczos (thick restart, locking)"
  )))
  expect_false(eigencore:::plan_dispatches_lanczos(list(
    method = "native dense Hermitian LAPACK fallback"
  )))

  expect_true(eigencore:::plan_dispatches_lobpcg(list(
    method = eigencore:::native_generalized_lobpcg_label()
  )))
  expect_true(eigencore:::plan_dispatches_lobpcg(list(
    method = eigencore:::reference_generalized_lobpcg_label()
  )))
  expect_false(eigencore:::plan_dispatches_lobpcg(list(
    method = "dense LAPACK general eigen oracle (prototype fallback)"
  )))

  expect_true(eigencore:::plan_dispatches_golub_kahan(list(
    method = "native prototype Golub-Kahan"
  )))
  expect_true(eigencore:::plan_dispatches_golub_kahan(list(
    method = "prototype Golub-Kahan"
  )))
  expect_false(eigencore:::plan_dispatches_golub_kahan(list(
    method = "native certified Gram SVD special case"
  )))
})

test_that("generalized SPD eigenproblem is certified in original coordinates", {
  A <- diag(c(6, 4, 2))
  B <- diag(c(3, 2, 1))
  fit <- eig_partial(A, B = B, k = 2, target = largest())

  expect_equal(unname(values(fit)), c(2, 2), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$plan$method, "native dense generalized SPD LAPACK fallback")
  expect_match(fit$warnings, "native dense generalized SPD LAPACK fallback")
})

test_that("native generalized SPD fallback produces B-orthonormal vectors", {
  set.seed(26)
  A0 <- matrix(rnorm(25), nrow = 5)
  B0 <- matrix(rnorm(25), nrow = 5)
  A <- crossprod(A0) + diag(5)
  B <- crossprod(B0) + diag(5)
  fit <- eig_partial(A, B = B, k = 3, target = smallest())

  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(3), tolerance = 1e-10)
  expect_lt(max(abs(A %*% vectors(fit) - B %*% vectors(fit) %*% diag(values(fit)))), 1e-8)
})

test_that("native generalized SPD LOBPCG supports dense SPD problems", {
  A <- diag(c(1, 4, 9, 16, 25))
  B <- diag(c(1, 2, 3, 4, 5))
  fit <- eig_partial(A, B = B, k = 2, target = smallest(),
                     method = lobpcg(maxit = 50L), seed = 28, tol = 1e-8)

  expect_equal(fit$plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_equal(values(fit), c(1, 2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(2), tolerance = 1e-8)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_true(fit$restart$orthogonalization_native)
  expect_true("native_dense_b_mgs2" %in% fit$restart$orthogonalization_methods)
})

test_that("native generalized SPD LOBPCG slice supports dense A with diagonal B", {
  A <- diag(c(1, 4, 9, 16, 25))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5))
  fit <- eig_partial(A, B = B, k = 2, target = smallest(),
                     method = lobpcg(maxit = 50L), seed = 281, tol = 1e-8)

  expect_equal(fit$plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_equal(values(fit), c(1, 2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_true("native_diagonal_b_mgs2" %in% fit$restart$orthogonalization_methods)
})

test_that("native generalized SPD LOBPCG handles sparse operators without densifying", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5))
  fit <- eig_partial(A, B = B, k = 2, target = smallest(),
                     method = lobpcg(maxit = 50L), seed = 29, tol = 1e-8)

  expect_equal(fit$plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_equal(values(fit), c(1, 2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$native)
  expect_true(fit$restart$orthogonalization_native)
  expect_true("native_diagonal_b_mgs2" %in% fit$restart$orthogonalization_methods)
})

test_that("native generalized SPD LOBPCG supports sparse CSC A and sparse CSC B", {
  A <- Matrix::sparseMatrix(
    i = 1:5,
    j = 1:5,
    x = c(1, 4, 9, 16, 25),
    dims = c(5, 5)
  )
  B <- Matrix::sparseMatrix(
    i = 1:5,
    j = 1:5,
    x = c(1, 2, 3, 4, 5),
    dims = c(5, 5)
  )
  fit <- eig_partial(A, B = B, k = 2, target = smallest(),
                     method = lobpcg(maxit = 50L), seed = 291, tol = 1e-8)

  expect_equal(fit$plan$method, eigencore:::native_generalized_lobpcg_label())
  expect_equal(values(fit), c(1, 2), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)
  expect_true(fit$restart$native)
  expect_true("native_csc_b_mgs2" %in% fit$restart$orthogonalization_methods)
})

test_that("preconditioned generalized SPD LOBPCG uses native shifted-tridiagonal path", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5))
  preconditioner <- shifted_tridiagonal_preconditioner(
    Matrix::bandSparse(5, k = 0, diagonals = list(c(1, 4, 9, 16, 25))),
    shift = 1e-3
  )
  problem <- eigen_problem(A, metric = B, target = smallest())
  plan <- plan_solver(problem, k = 2, method = lobpcg(maxit = 50L, preconditioner = preconditioner))

  expect_false(eigencore:::native_lobpcg_supported(
    as_operator(A),
    target = smallest(),
    preconditioner = preconditioner,
    Bop = as_operator(B)
  ))
  expect_true(eigencore:::native_generalized_lobpcg_supported(
    as_operator(A),
    as_operator(B),
    target = smallest(),
    preconditioner = preconditioner
  ))
  expect_equal(plan$method, eigencore:::native_generalized_lobpcg_label())

  fit <- eig_partial(A, B = B, k = 2, target = smallest(),
                     method = lobpcg(maxit = 50L, preconditioner = preconditioner),
                     seed = 30, tol = 1e-8)

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$preconditioner_native)
  expect_true(certificate(fit)$passed)
})

test_that("svd_partial returns sorted singular values and certificate", {
  A <- diag(c(5, 1, 3))
  fit <- svd_partial(A, rank = 2)

  expect_equal(values(fit), c(5, 3))
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$plan$method, "native dense LAPACK SVD fallback")
  expect_equal(fit$method, "native dense LAPACK SVD fallback")
})

test_that("native dense SVD fallback matches base SVD on rectangular inputs", {
  set.seed(25)
  A <- matrix(rnorm(35), nrow = 7)
  fit <- svd_partial(A, rank = 3, target = largest())
  oracle <- svd(A, nu = 3, nv = 3)

  expect_equal(values(fit), oracle$d[1:3], tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_match(fit$warnings, "native dense LAPACK SVD fallback")
})

test_that("randomized SVD path returns honest certified results on low-rank input", {
  set.seed(31)
  U <- qr.Q(qr(matrix(rnorm(12 * 3), nrow = 12, ncol = 3)))
  V <- qr.Q(qr(matrix(rnorm(8 * 3), nrow = 8, ncol = 3)))
  A <- U %*% diag(c(9, 5, 2), nrow = 3) %*% t(V)
  fit <- svd_partial(A, rank = 2, method = randomized(oversample = 4, n_iter = 1),
                     seed = 31, tol = 1e-8)

  expect_equal(fit$plan$method, "reference randomized SVD prototype")
  expect_equal(values(fit), c(9, 5), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$restart$kind, "randomized_range_finder")
  expect_false(fit$restart$native)
})

test_that("shift_invert dense reference path returns certified original eigenpairs", {
  A <- diag(c(1, 3, 7))
  fit <- eig_partial(A, k = 1, target = nearest(2.8), method = shift_invert(2.8))

  expect_equal(values(fit), 3, tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$transform$kind, "shift_invert")
  expect_equal(fit$transform$label_kind, "dense_qr")
  expect_match(fit$warnings, "reference Hermitian Lanczos shift-invert")
})

test_that("dense fallback is memory-budgeted", {
  old <- getOption("eigencore.dense_fallback_mb")
  on.exit(options(eigencore.dense_fallback_mb = old), add = TRUE)
  options(eigencore.dense_fallback_mb = 16 / 1e6)

  expect_error(
    svd_partial(matrix(1, nrow = 4, ncol = 4), rank = 1),
    "exceeding eigencore.dense_fallback_mb"
  )
})

test_that("dense fallback can be disabled explicitly", {
  expect_error(
    svd_partial(matrix(1, nrow = 3, ncol = 3), rank = 1, allow_dense_fallback = "never"),
    "allow_dense_fallback = 'never'"
  )
})

test_that("solver paths refuse implicit sparse densification", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 2, 2, 3, 3),
    x = c(3, 1, 2, 1, 1),
    dims = c(3, 3)
  )

  expect_error(
    eig_partial(A, k = 1),
    "Refusing to densify sparse A"
  )

  explicit_fit <- eig_partial(A, k = 1, allow_dense_fallback = "always")
  expect_equal(explicit_fit$plan$method, "dense LAPACK general eigen oracle (prototype fallback)")

  fit <- eig_partial(as.matrix(A), k = 1)
  expect_equal(fit$plan$method, "dense LAPACK general eigen oracle (prototype fallback)")
})

test_that("RSpectra-compatible shims expose core fields", {
  A <- diag(c(3, 2, 1))
  ef <- eigs_sym(A, k = 1)
  sf <- svds(A, k = 1)

  expect_equal(ef$values, 3)
  expect_equal(sf$d, 3)
  expect_s3_class(ef$certificate, "eigencore_certificate")
  expect_s3_class(sf$certificate, "eigencore_certificate")
})

test_that("RSpectra SM maps to smallest magnitude, not smallest algebraic", {
  A <- diag(c(-10, 2, 5))
  ef <- eigs_sym(A, k = 1, which = "SM")

  expect_equal(ef$values, 2)
})

test_that("prototype Lanczos solves matrix-free Hermitian operators", {
  vals <- c(6, 5, 4, 3, 2, 1)
  op <- linear_operator(
    dim = c(length(vals), length(vals)),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    name = "diagonal_matrix_free",
    metadata = list(frobenius_norm = sqrt(sum(vals^2)))
  )

  fit <- eig_partial(op, k = 2, method = lanczos(max_subspace = 6), seed = 123)

  expect_equal(values(fit), c(6, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$method, "reference Hermitian Lanczos (prototype/oracle fallback)")
})

test_that("auto uses Lanczos for matrix-free Hermitian operators", {
  vals <- c(4, 3, 2, 1)
  op <- linear_operator(
    dim = c(4, 4),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (vals * X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = hermitian()
  )

  fit <- eig_partial(op, k = 1, seed = 321)
  expect_equal(values(fit), 4, tolerance = 1e-10)
  expect_equal(fit$plan$method, "reference Hermitian Lanczos (prototype/oracle fallback)")
})

test_that("auto uses native CSC-backed Lanczos for sparse Hermitian matrices", {
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(1, 2, 3, 4),
    x = c(8, 5, 3, 1),
    dims = c(4, 4)
  )
  op <- as_operator(A)
  fit <- eig_partial(A, k = 2, seed = 101)

  expect_null(eigencore:::source_or_null(op))
  expect_equal(fit$plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(fit$method, "native scalar thick-restart Hermitian Lanczos")
  expect_identical(fit$warnings, character())
  expect_equal(values(fit), c(8, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("explicit Lanczos uses the native thick-restart path for dense Hermitian matrices", {
  A <- diag(c(8, 5, 3, 1))
  fit <- eig_partial(A, k = 2, method = lanczos(max_subspace = 4), seed = 102)

  expect_equal(fit$plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(fit$method, "native scalar thick-restart Hermitian Lanczos")
  expect_identical(fit$warnings, character())
  expect_equal(values(fit), c(8, 5), tolerance = 1e-10)
  expect_lte(fit$iterations, 4L)
  expect_gte(fit$matvecs, fit$iterations)
  expect_true(all(c("restart", "iteration", "n_locked") %in% names(fit$convergence_history)))
  expect_equal(fit$restart$kind, "thick_restart")
  expect_true(fit$restart$implemented)
  expect_equal(fit$restart$locking, "in_native_loop")
  expect_equal(fit$restart$locked_count, fit$nconv)
  expect_equal(fit$locked, seq_len(fit$nconv))
  expect_true(certificate(fit)$passed)
})

test_that("native Lanczos does not silently remap unsupported targets", {
  A <- diag(c(1, 5, 10))
  fit <- eig_partial(A, k = 1, target = nearest(5),
                     method = lanczos(max_subspace = 3), seed = 103)

  expect_equal(fit$plan$method, "reference Hermitian Lanczos (prototype/oracle fallback)")
  expect_equal(values(fit), 5, tolerance = 1e-12)
  expect_true(certificate(fit)$passed)
})

test_that("explicit Lanczos on non-Hermitian inputs is honestly routed to an oracle fallback", {
  A <- matrix(c(1, 10, 0, 2), 2, 2)
  fit <- eig_partial(A, k = 1, method = lanczos(max_subspace = 2), seed = 104)

  expect_equal(fit$plan$method, "dense LAPACK eigen oracle (Lanczos requires Hermitian structure)")
  expect_equal(values(fit), 2, tolerance = 1e-12)
  expect_true(certificate(fit)$passed)
})

test_that("native Lanczos iteration matches diagonal dense spectrum", {
  A <- diag(c(7, 4, 2, 1))
  op <- as_operator(A)
  set.seed(103)
  out <- eigencore:::native_lanczos_hermitian(op, k = 2, maxit = 4, target = largest())

  expect_equal(out$values, c(7, 4), tolerance = 1e-10)
  expect_lte(out$iterations, 4L)
  expect_gte(out$matvecs, out$iterations)
  expect_equal(out$restart$kind, "thick_restart")
  expect_true(out$restart$implemented)
  expect_equal(out$restart$locked_count, length(out$locked))
  expect_true(out$certificate$passed)
})

test_that("native Lanczos convergence control can stop before maxit", {
  A <- diag(c(10, 1, 0.1, 0.01))
  op <- as_operator(A)
  set.seed(104)
  out <- eigencore:::native_lanczos_hermitian(op, k = 1, maxit = 4, target = largest(), tol = 1)

  expect_equal(out$values, 10, tolerance = 1e-10)
  expect_true(out$certificate$passed)
  expect_equal(out$restart$locked_count, 1L)
})

test_that("native Lanczos thick-restart metadata is visible in diagnostics", {
  fit <- eig_partial(diag(c(6, 3, 1)), k = 1, method = lanczos(max_subspace = 3), seed = 105)
  diag <- diagnostics(fit)

  expect_identical(fit$warnings, character())
  expect_equal(diag$restart$kind, "thick_restart")
  expect_equal(diag$restart$locking, "in_native_loop")
  expect_equal(diag$locked, fit$locked)
  expect_equal(fit$locked, 1L)
})

test_that("native thick-restart Lanczos converges with m_max << k+convergence-need", {
  set.seed(42)
  values <- c(50, 40, 30, 20, 10, 5, 4, 3, 2, 1)
  A <- symmetric_with_spectrum(values, seed = 11)
  k <- 3L
  m_max <- as.integer(k + 2L)  # well below the 10-dim full subspace

  fit <- eig_partial(A, k = k, target = largest(),
                     method = lanczos(max_subspace = m_max, max_restarts = 200L),
                     seed = 42, tol = 1e-10)

  expect_equal(fit$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(values(fit), c(50, 40, 30), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_gt(fit$restart$restarts_used, 0L)
  expect_equal(fit$restart$max_subspace, m_max)
  expect_lte(fit$restart$final_active_subspace, m_max)
  expect_gte(fit$restart$final_active_subspace, k)
  expect_equal(fit$restart$locked_count, k)
  expect_equal(fit$nconv, k)
})

test_that("native block Lanczos prototype matches dense oracle on small Hermitian problems", {
  set.seed(44)
  values0 <- c(9, 5, 2, 1, -1, -3)
  A <- symmetric_with_spectrum(values0, seed = 44)

  fit <- eig_partial(A, k = 2L, target = largest(),
                     method = lanczos(block = 2L, max_subspace = 6L),
                     seed = 44, tol = 1e-8)

  expect_equal(fit$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(values(fit), c(9, 5), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$restart$kind, "block_full_subspace_dense_lapack")
  expect_equal(fit$restart$locking, "not_required_full_subspace")
  expect_equal(fit$restart$block, 2L)
  expect_equal(fit$restart$final_active_subspace, nrow(A))
  expect_equal(fit$matvecs, 0L)
})

test_that("planner labels explicit block Lanczos with production label and stores controls", {
  A <- Matrix::sparseMatrix(i = 1:4, j = 1:4, x = c(4, 3, 2, 1))
  P <- eigen_problem(A, target = largest())
  plan <- plan_solver(P, k = 2L, method = lanczos(block = 2L, max_subspace = 4L))

  expect_equal(plan$method, "native block Hermitian Lanczos (thick restart, locking)")
  expect_equal(plan$controls$block, 2L)
  expect_equal(plan$controls$max_subspace, 4L)
  expect_equal(plan$controls$max_restarts, 100L)
})

test_that("auto planner keeps small-k sparse requests on scalar Lanczos", {
  A <- Matrix::sparseMatrix(i = 1:8, j = 1:8, x = 8:1)
  P <- eigen_problem(A, target = largest())
  plan <- plan_solver(P, k = 4L, method = auto())

  expect_equal(plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(plan$controls$block, 1L)
})

test_that("reference LOBPCG spike shows preconditioner leverage on Laplacian", {
  n <- 80L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  ridge_factor <- Matrix::Cholesky(A + Matrix::Diagonal(n) * 1e-3, LDL = FALSE)
  preconditioner <- function(R) as.matrix(Matrix::solve(ridge_factor, R))

  unpre <- eigencore:::reference_lobpcg_hermitian(
    A, k = 3L, target = smallest(), maxit = 80L, tol = 1e-8, seed = 91
  )
  pre <- eigencore:::reference_lobpcg_hermitian(
    A, k = 3L, target = smallest(), maxit = 80L, tol = 1e-8,
    preconditioner = preconditioner, seed = 91
  )

  expect_true(pre$certificate$passed)
  expect_lt(pre$iterations, unpre$iterations)
  expect_lt(pre$iterations, 20L)
})

test_that("eig_partial routes explicit lobpcg method with preconditioner", {
  n <- 80L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  preconditioner <- shifted_cholesky_preconditioner(A, shift = 1e-3)
  expect_s3_class(preconditioner, "eigencore_preconditioner")
  cholesky_info <- eigencore:::eigencore_preconditioner_info(preconditioner)
  expect_equal(cholesky_info$kind, "shifted_cholesky")
  expect_false(cholesky_info$native)
  fit <- eig_partial(
    A,
    k = 3L,
    target = smallest(),
    method = lobpcg(maxit = 80L, preconditioner = preconditioner),
    tol = 1e-8,
    seed = 91
  )

  expect_equal(fit$method, "reference LOBPCG prototype")
  expect_true(certificate(fit)$passed)
  expect_true(fit$restart$preconditioned)
  expect_equal(fit$restart$preconditioner_kind, "shifted_cholesky")
  expect_false(fit$restart$preconditioner_native)
  expect_equal(fit$restart$preconditioner_calls, fit$preconditioner_calls)
  expect_equal(fit$restart$preconditioner_calls, fit$iterations - 1L)
  expect_equal(diagnostics(fit)$preconditioner$kind, "shifted_cholesky")
  expect_match(paste(fit$plan$reasons, collapse = "\n"), "preconditioner: shifted_cholesky")
  expect_lt(fit$iterations, 20L)
  expect_equal(nrow(fit$convergence_history), fit$iterations)
  expect_equal(tail(fit$convergence_history$nconv, 1), fit$nconv)
  expect_lt(tail(fit$convergence_history$max_relative_residual, 1), 1e-8)
  expect_match(fit$warnings, "R-level reference LOBPCG")
})

test_that("shifted tridiagonal preconditioner matches shifted dense solve", {
  n <- 20L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  R <- matrix(rnorm(n * 3L), nrow = n)
  preconditioner <- shifted_tridiagonal_preconditioner(A, shift = 1e-3)
  expect_s3_class(preconditioner, "eigencore_preconditioner")
  tridiagonal_info <- eigencore:::eigencore_preconditioner_info(preconditioner)
  expect_equal(tridiagonal_info$kind, "shifted_tridiagonal")
  expect_true(tridiagonal_info$native)
  got <- preconditioner(R)
  expected <- solve(as.matrix(A + Matrix::Diagonal(n) * 1e-3), R)

  expect_equal(got, expected, tolerance = 1e-10)
})

test_that("native LOBPCG routes typed tridiagonal preconditioner through native loop", {
  n <- 80L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  preconditioner <- shifted_tridiagonal_preconditioner(A, shift = 1e-3)
  info <- eigencore:::eigencore_preconditioner_info(preconditioner)
  expect_true(all(c("lower", "diag", "upper") %in% names(info)))

  fit <- eig_partial(
    A,
    k = 3L,
    target = smallest(),
    method = lobpcg(maxit = 80L, preconditioner = preconditioner),
    tol = 1e-8,
    seed = 91
  )

  expect_equal(fit$method, "native standard Hermitian LOBPCG prototype")
  expect_equal(fit$plan$method, "native standard Hermitian LOBPCG prototype")
  expect_true(certificate(fit)$passed)
  expect_true(fit$restart$native)
  expect_true(fit$restart$native_kernels)
  expect_true(fit$restart$preconditioned)
  expect_equal(fit$restart$preconditioner_kind, "shifted_tridiagonal")
  expect_true(fit$restart$preconditioner_native)
  expect_equal(fit$restart$preconditioner_calls, fit$iterations - 1L)
  expect_equal(diagnostics(fit)$preconditioner$kind, "shifted_tridiagonal")
  expect_equal(tail(fit$convergence_history$nconv, 1), fit$nconv)
  expect_lt(tail(fit$convergence_history$max_relative_residual, 1), 1e-8)
  expect_match(fit$warnings, "native standard LOBPCG")
})

test_that("native thick-restart Lanczos locks the converged Ritz pairs in target order", {
  set.seed(7)
  values <- c(20, 19, 5, 4, 3, 2)
  A <- symmetric_with_spectrum(values, seed = 12)
  fit <- eig_partial(A, k = 2L, target = largest(),
                     method = lanczos(max_subspace = 4L, max_restarts = 100L),
                     seed = 7)

  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), c(20, 19), tolerance = 1e-8)
  # Locked pairs are precisely the converged ones
  expect_equal(fit$locked, which(fit$certificate$converged))
  # Each locked Ritz pair satisfies its own residual <= tol scale
  cert <- certificate(fit)
  for (i in fit$locked) {
    scale_i <- max(abs(values(fit)[i]), 1)
    expect_lte(cert$residuals[i], cert$tolerance * scale_i * 10)
  }
})

test_that("native thick-restart Lanczos honestly reports non-convergence under restart starvation", {
  set.seed(1)
  # Closely spaced spectrum so a tiny restart budget can't finish.
  values <- 10 + (1:30) * 1e-3
  A <- symmetric_with_spectrum(values, seed = 13)

  fit <- eig_partial(A, k = 5L, target = largest(),
                     method = lanczos(max_subspace = 7L, max_restarts = 1L),
                     seed = 1, tol = 1e-12)

  expect_equal(fit$method, "native scalar thick-restart Hermitian Lanczos")
  # Some pairs may converge, but not all five.
  expect_lt(fit$nconv, 5L)
  expect_false(certificate(fit)$passed)
  expect_match(fit$warnings, "exhausted")
  # The result still ships k_target slots, with non-converged ones flagged.
  expect_length(values(fit), 5L)
  expect_length(fit$certificate$converged, 5L)
  expect_true(any(!fit$certificate$converged))
})

test_that("native thick-restart Lanczos residuals match independently recomputed residuals", {
  set.seed(5)
  values <- c(15, 12, 8, 5, 3, 2, 1)
  A <- symmetric_with_spectrum(values, seed = 14)

  fit <- eig_partial(A, k = 3L, target = largest(),
                     method = lanczos(max_subspace = 5L, max_restarts = 100L),
                     seed = 5, tol = 1e-10)

  expect_true(certificate(fit)$passed)
  V <- vectors(fit)
  lambda <- values(fit)
  recomputed <- vapply(seq_along(lambda), function(i) {
    sqrt(sum((A %*% V[, i] - lambda[i] * V[, i])^2))
  }, numeric(1))
  expect_equal(unname(fit$certificate$residuals), recomputed, tolerance = 1e-10)
})

test_that("native thick-restart Lanczos respects sparse non-densification", {
  set.seed(9)
  n <- 60L
  # Top 4 dominant eigenvalues well-separated from a tail; restart converges.
  diag_vals <- c(50, 30, 20, 10, sort(runif(n - 4L, 0.1, 5), decreasing = TRUE))
  A <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n),
                            x = diag_vals, dims = c(n, n))
  fit <- eig_partial(A, k = 4L, target = largest(),
                     method = lanczos(max_subspace = 12L, max_restarts = 200L),
                     seed = 9, tol = 1e-8)

  expect_equal(fit$plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(values(fit), diag_vals[1:4], tolerance = 1e-6)
  expect_true(certificate(fit)$passed)
})

test_that("Lanczos Ritz extraction preserves projected eigenvector rotation", {
  A <- diag(c(5, 3, 1))
  op <- as_operator(A)
  Q <- cbind(c(1, 0, 0), c(1e-3, sqrt(1 - 1e-6), 0))
  alpha <- c(2, 1)
  beta <- c(0.2, 0)
  projected <- eigencore:::tridiagonal_matrix(alpha, beta)
  eig <- eigen(projected, symmetric = TRUE)
  idx <- eigencore:::order_indices(eig$values, largest())
  expected <- Q %*% eig$vectors[, idx, drop = FALSE]

  out <- eigencore:::reference_lanczos_ritz(op, Q, alpha, beta, k = 2, target = largest(), tol = 1e-8)

  expect_equal(out$values, eig$values[idx], tolerance = 1e-14)
  expect_equal(out$vectors, expected, tolerance = 1e-14)
})

test_that("native tridiagonal eigensolve matches dense projected oracle", {
  alpha <- c(2, 1, 3, 4)
  beta <- c(0.2, -0.5, 0.1, 0)
  projected <- eigencore:::tridiagonal_matrix(alpha, beta)
  native <- eigencore:::native_tridiagonal_eigen(alpha, beta)
  oracle <- eigen(projected, symmetric = TRUE)
  idx <- order(native$values, decreasing = TRUE)

  expect_equal(native$values[idx], oracle$values, tolerance = 1e-12)
  expect_equal(abs(crossprod(native$vectors[, idx], oracle$vectors)), diag(length(alpha)), tolerance = 1e-10)
})

test_that("prototype Golub-Kahan solves matrix-free rectangular SVD", {
  sing <- c(7, 5, 3, 1)
  m <- 6
  n <- length(sing)
  A <- rbind(diag(sing), matrix(0, m - n, n))
  op <- linear_operator(
    dim = c(m, n),
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
    name = "rectangular_matrix_free",
    metadata = list(frobenius_norm = norm(A, type = "F"))
  )

  fit <- svd_partial(op, rank = 2, method = golub_kahan(max_subspace = 4), seed = 123)

  expect_equal(values(fit), c(7, 5), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_equal(fit$nconv, 2)
  expect_equal(fit$method, "prototype Golub-Kahan")
})

test_that("auto uses Golub-Kahan for matrix-free SVD", {
  A <- rbind(diag(c(4, 2, 1)), 0)
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
    }
  )

  fit <- svd_partial(op, rank = 1, seed = 456)
  expect_equal(values(fit), 4, tolerance = 1e-10)
  expect_equal(fit$plan$method, "prototype Golub-Kahan")
})

test_that("auto uses native CSC-backed Golub-Kahan for sparse rectangular SVD", {
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4),
    j = c(1, 2, 3, 4),
    x = c(9, 6, 2, 1),
    dims = c(6, 4)
  )
  op <- as_operator(A)
  fit <- svd_partial(A, rank = 2, seed = 202)

  expect_null(eigencore:::source_or_null(op))
  expect_equal(fit$plan$method, "native retained Golub-Kahan SVD (thick restart)")
  expect_equal(fit$method, "native retained Golub-Kahan SVD (thick restart)")
  expect_match(fit$warnings, "native retained Golub-Kahan")
  expect_equal(values(fit), c(9, 6), tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
})

test_that("explicit Golub-Kahan uses native prototype for dense rectangular SVD", {
  A <- rbind(diag(c(9, 6, 2, 1)), matrix(0, 2, 4))
  fit <- svd_partial(A, rank = 2, method = golub_kahan(max_subspace = 4), seed = 203)

  expect_equal(fit$plan$method, "native prototype Golub-Kahan")
  expect_equal(fit$method, "native prototype Golub-Kahan")
  expect_match(fit$warnings, "native prototype Golub-Kahan")
  expect_equal(values(fit), c(9, 6), tolerance = 1e-10)
  expect_equal(fit$matvecs, 2L * fit$iterations)
  expect_true(certificate(fit)$passed)
})

test_that("native Golub-Kahan iteration matches dense diagonal singular values", {
  A <- rbind(diag(c(7, 4, 2, 1)), matrix(0, 2, 4))
  op <- as_operator(A)
  set.seed(204)
  out <- eigencore:::native_golub_kahan_svd(op, rank = 2, maxit = 4, target = largest())

  expect_equal(out$d, c(7, 4), tolerance = 1e-10)
  expect_equal(out$matvecs, 2L * out$iterations)
  expect_true(out$certificate$passed)
})

test_that("native Golub-Kahan exposes adaptive subspace metadata", {
  old_options <- options(
    eigencore.golub_kahan_prefix_diagnostics = TRUE,
    eigencore.golub_kahan_projected_stop = TRUE
  )
  on.exit(options(old_options), add = TRUE)
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 5),
    j = c(1, 2, 3, 4, 5),
    x = c(9, 7, 5, 2, 1),
    dims = c(8, 5)
  )
  fit <- svd_partial(A, rank = 3, tol = 1e-8, seed = 445)

  expect_equal(fit$method, "native retained Golub-Kahan SVD (thick restart)")
  expect_true(all(c(
    "kind", "implemented", "native", "thick_restart",
    "retained_restart", "retained_restart_native",
    "native_workspace_bytes", "basis_returned", "attempt_history"
  ) %in% names(fit$restart)))
  expect_equal(fit$restart$kind, "block_golub_kahan_native_retained_cycle")
  expect_true(fit$restart$thick_restart)
  expect_true(fit$restart$retained_restart)
  expect_true(fit$restart$implemented)
  expect_true(fit$restart$retained_restart_native)
  expect_true(fit$restart$native_attempt_certification)
  expect_true(is.data.frame(fit$restart$attempt_history))
  expect_true(any(fit$restart$attempt_history$certificate_passed))
  expect_gte(fit$restart$final_max_subspace, fit$restart$final_iterations)
  expect_gte(fit$nconv, 3L)
  expect_gt(fit$restart$native_workspace_bytes, 0)
  expect_false(fit$restart$basis_returned)
  expect_gt(fit$restart$total_ortho_passes, 0L)
  expect_true(all(c(
    "attempt", "max_subspace", "iterations", "matvecs",
    "certificate_passed", "max_residual", "max_backward_error",
    "ortho_passes"
  ) %in% names(fit$restart$attempt_history)))
  expect_true(all(c(
    "native_iteration", "ritz", "restart"
  ) %in% names(fit$restart$stage_seconds)))
  expect_gte(fit$restart$stage_seconds[["native_iteration"]], 0)
  expect_gte(fit$restart$stage_seconds[["ritz"]], 0)
  expect_true(utils::tail(fit$restart$attempt_history$certificate_passed, 1L))
  expect_lte(fit$restart$certified_attempt, fit$restart$attempts)
  expect_gte(fit$restart$final_attempt_matvecs, 0L)
  expect_true(certificate(fit)$passed)
})

test_that("native Golub-Kahan compact fit avoids returning Krylov basis by default", {
  old_options <- options(
    eigencore.golub_kahan_prefix_diagnostics = FALSE,
    eigencore.golub_kahan_projected_stop = TRUE
  )
  on.exit(options(old_options), add = TRUE)
  A <- matrix(rnorm(12 * 7), 12, 7)
  fit <- svd_partial(A, rank = 3, tol = 1e-8, seed = 446, method = golub_kahan())

  expect_equal(fit$method, "native prototype Golub-Kahan")
  expect_true(fit$certificate$passed)
  expect_false(fit$restart$basis_returned)
  expect_gt(fit$restart$native_workspace_bytes, 0)
  expect_true(is.data.frame(fit$restart$prefix_history))
  expect_equal(nrow(fit$restart$prefix_history), 0L)
})

test_that("Golub-Kahan Ritz extraction preserves coupled SVD rotations", {
  A <- rbind(diag(c(5, 2)), c(0, 0))
  op <- as_operator(A)
  U <- cbind(c(1, 0, 0), c(1e-3, sqrt(1 - 1e-6), 0))
  V <- cbind(c(1, 0), c(1e-3, sqrt(1 - 1e-6)))
  alpha <- c(3, 2)
  beta <- c(0.4, 0)
  projected <- eigencore:::bidiagonal_matrix(alpha, beta)
  bd <- svd(projected, nu = nrow(projected), nv = ncol(projected))
  idx <- eigencore:::order_indices(bd$d, largest())

  out <- eigencore:::reference_golub_kahan_ritz(op, U, V, alpha, beta, rank = 2, target = largest(), tol = 1e-8)
  native <- eigencore:::native_golub_kahan_ritz(op, U, V, alpha, beta, rank = 2, target = largest(), tol = 1e-8)

  expect_equal(out$d, bd$d[idx], tolerance = 1e-14)
  expect_equal(out$u, U %*% bd$u[, idx, drop = FALSE], tolerance = 1e-14)
  expect_equal(out$v, V %*% bd$v[, idx, drop = FALSE], tolerance = 1e-14)
  expect_equal(native$d, out$d, tolerance = 1e-14)
  expect_equal(native$u, out$u, tolerance = 1e-14)
  expect_equal(native$v, out$v, tolerance = 1e-14)

  U_padded <- cbind(U, c(0, 0, 1))
  V_padded <- cbind(V, c(0, 0))
  alpha_padded <- c(alpha, 999)
  beta_padded <- c(beta, 999)
  active <- eigencore:::native_golub_kahan_ritz(
    op,
    U_padded,
    V_padded,
    alpha_padded,
    beta_padded,
    rank = 2,
    target = largest(),
    tol = 1e-8,
    active_iterations = 2
  )
  expect_equal(active$d, out$d, tolerance = 1e-14)
  expect_equal(active$u, out$u, tolerance = 1e-14)
  expect_equal(active$v, out$v, tolerance = 1e-14)
})

test_that("native bidiagonal SVD matches dense projected oracle", {
  alpha <- c(3, 2, 1)
  beta <- c(0.4, -0.2, 0)
  projected <- eigencore:::bidiagonal_matrix(alpha, beta)
  native <- eigencore:::native_bidiagonal_svd(alpha, beta)
  oracle <- svd(projected, nu = length(alpha), nv = length(alpha))

  expect_equal(native$d, oracle$d, tolerance = 1e-12)
  expect_equal(abs(crossprod(native$u, oracle$u)), diag(length(alpha)), tolerance = 1e-10)
  expect_equal(abs(crossprod(native$v, oracle$v)), diag(length(alpha)), tolerance = 1e-10)
})

test_that("Golub-Kahan vector modes withhold full certificate when needed", {
  A <- rbind(diag(c(3, 2, 1)), 0)
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (A %*% X),
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (t(A) %*% X)
  )

  fit <- svd_partial(op, rank = 1, method = golub_kahan(max_subspace = 3), vectors = "right", seed = 789)
  expect_null(left_vectors(fit))
  expect_false(certificate(fit)$passed)
  expect_match(certificate(fit)$notes, "both left and right vectors")
})
