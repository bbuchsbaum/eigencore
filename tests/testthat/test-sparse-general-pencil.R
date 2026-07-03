expect_complex_set_equal <- function(actual, expected, tolerance = 1e-8) {
  actual <- as.complex(actual)
  expected <- as.complex(expected)
  expect_equal(length(actual), length(expected))
  remaining <- seq_along(actual)
  for (target in expected) {
    distances <- Mod(actual[remaining] - target)
    best <- which.min(distances)
    expect_lte(distances[[best]], tolerance)
    remaining <- remaining[-best]
  }
}

general_pencil_oracle <- function(A, B, target, k) {
  vals <- eigen(solve(as.matrix(B), as.matrix(A)), only.values = TRUE)$values
  idx <- eigencore:::order_indices(vals, target)
  vals[idx[seq_len(k)]]
}

test_that("sparse diagonal-B general pencils use transformed Arnoldi", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 2, 2, 3, 3),
    x = c(5, 1, 3, 2, 1),
    dims = c(3, 3)
  )
  B <- Matrix::Diagonal(x = c(1, 2, 3))
  target <- largest_real()
  expected <- general_pencil_oracle(A, B, target, k = 2L)

  fit <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = target,
    tol = 1e-10,
    allow_dense_fallback = "never"
  )

  expect_equal(
    fit$method,
    eigencore:::sparse_general_pencil_diagonal_arnoldi_label()
  )
  expect_equal(fit$plan$fallback, "none; sparse general-pencil boundary is explicit")
  expect_equal(sort(values(fit), decreasing = TRUE), sort(Re(expected), decreasing = TRUE),
               tolerance = 1e-10)
  expect_equal(fit$certificate$certificate_type,
               "generalized_pencil_right_residual_backward_error")
  expect_true(certificate(fit)$passed)
  expect_equal(fit$classification, rep("finite", 2L))
  expect_equal(alpha_beta(fit)$beta, rep(1, 2L))
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$right_hand_pencil)
  expect_equal(fit$restart$metric_solve, "nonsingular diagonal B row scaling")
  expect_false(fit$restart$materialized_dense_operator)
  expect_true(fit$restart$materialized_sparse_operator)
  expect_equal(fit$transform$certification$residual_formula,
               "A * x - lambda * B * x")
})

test_that("explicit general sparse pencils with diagonal B expose planner provenance", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3),
    j = c(1, 2, 2, 3),
    x = c(4, 1, 2, -1),
    dims = c(3, 3)
  )
  B <- Matrix::Diagonal(x = c(1, -2, 3))
  problem <- eigen_problem(A, metric = B, structure = general(),
                           target = largest_magnitude())
  plan <- plan_solver(problem, k = 2L)

  expect_equal(
    plan$method,
    eigencore:::sparse_general_pencil_diagonal_arnoldi_label()
  )
  expect_equal(plan$controls$transformed_operator, "B^{-1} A")
  expect_true(plan$controls$right_hand_pencil)
  expect_equal(plan$controls$metric_solve, "nonsingular diagonal B row scaling")
  expect_true(plan$controls$unsupported_sparse_qz)

  fit <- solve(problem, k = 2L, tol = 1e-10, allow_dense_fallback = "never")
  expected <- general_pencil_oracle(A, B, largest_magnitude(), k = 2L)

  expect_complex_set_equal(values(fit), expected, tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  expect_false(fit$certificate$orthogonality_required)
})

test_that("real sparse general pencils can return certified complex pairs", {
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 3),
    j = c(2, 1, 3),
    x = c(-2, 2, 3),
    dims = c(3, 3)
  )
  B <- Matrix::Diagonal(x = c(1, 2, 3))

  fit <- eig_partial(
    A,
    B = B,
    k = 1L,
    target = largest_imaginary(),
    tol = 1e-10,
    allow_dense_fallback = "never"
  )

  expect_equal(Im(values(fit))[[1L]], sqrt(2), tolerance = 1e-8)
  expect_true(is.complex(values(fit)))
  expect_true(is.complex(vectors(fit)))
  expect_true(certificate(fit)$passed)
  expect_equal(fit$classification, "finite")
  expect_false(fit$restart$materialized_dense_operator)
})

test_that("unsupported sparse general pencils fail before dense fallback", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3),
    j = c(1, 2, 2, 3),
    x = c(4, 1, 2, -1),
    dims = c(3, 3)
  )
  B_singular <- Matrix::Diagonal(x = c(1, 0, 3))
  B_general <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3),
    j = c(1, 2, 2, 3),
    x = c(1, 1, 2, 3),
    dims = c(3, 3)
  )

  singular_problem <- eigen_problem(A, metric = B_singular,
                                    structure = general())
  general_problem <- eigen_problem(A, metric = B_general,
                                   structure = general())

  expect_equal(
    plan_solver(singular_problem, k = 2L)$method,
    eigencore:::sparse_general_pencil_unsupported_label()
  )
  expect_equal(
    plan_solver(general_problem, k = 2L)$method,
    eigencore:::sparse_general_pencil_unsupported_label()
  )
  expect_error(
    solve(singular_problem, k = 2L, allow_dense_fallback = "never"),
    "nonsingular diagonal B"
  )
  expect_error(
    solve(general_problem, k = 2L, allow_dense_fallback = "never"),
    "nonsingular diagonal B"
  )
})

test_that("sparse general-pencil planner uses a larger Arnoldi subspace budget", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3),
    j = c(1, 2, 2, 3),
    x = c(4, 1, 2, -1),
    dims = c(60, 60)
  )
  B <- Matrix::Diagonal(x = seq(1, 3, length.out = 60L))
  plan <- plan_solver(
    eigen_problem(A, metric = B, structure = general(), target = largest_real()),
    k = 2L
  )
  expect_gte(plan$controls$max_subspace, 35L)
})

test_that("sparse general-pencil Arnoldi certifies reliably under random starts", {
  skip_on_cran()
  sparse_n <- 60L
  k <- 2L
  set.seed(104)
  A <- Matrix::bandSparse(
    sparse_n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-1, sparse_n - 1L),
      seq(2, 4, length.out = sparse_n),
      rep(0.5, sparse_n - 1L)
    )
  )
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  B <- Matrix::Diagonal(x = seq(1, 3, length.out = sparse_n))
  passed <- vapply(seq_len(20L), function(i) {
    fit <- eig_partial(
      A,
      B = B,
      k = k,
      target = largest_real(),
      tol = 1e-9,
      allow_dense_fallback = "never"
    )
    certificate(fit)$passed
  }, logical(1))
  expect_true(all(passed))
})

test_that("Hermitian metric guard still rejects nonsymmetric B", {
  A <- Matrix::Diagonal(3, x = c(1, 4, 9))
  B_nonsym <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3),
    j = c(1, 2, 2, 3),
    x = c(1, 1, 2, 3),
    dims = c(3, 3)
  )

  expect_error(
    eigen_problem(A, metric = B_nonsym, structure = hermitian()),
    "Hermitian metric"
  )
})
