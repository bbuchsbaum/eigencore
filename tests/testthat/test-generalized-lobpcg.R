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

test_that("explicit generalized Lanczos uses native transformed diagonal-B path", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))
  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = lanczos(max_subspace = 6L),
    seed = 321,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lanczos_label())
  expect_equal(fit$plan$method, eigencore:::native_generalized_lanczos_label())
  expect_match(fit$warnings, "native transformed generalized SPD B-orthogonal Lanczos", fixed = TRUE)
  expect_true(fit$generalized)
  expect_equal(fit$restart$kind, "native_transformed_generalized_b_orthogonal_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$native_kernels)
  expect_equal(fit$restart$metric_solve, "diagonal scaling similarity transform for B")
  expect_equal(fit$restart$transformed_operator_storage, "dgCMatrix")
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), c(1, 2), tolerance = 1e-8)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)

  lobpcg_fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = lobpcg(maxit = 80L),
    seed = 321,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )
  expect_true(certificate(lobpcg_fit)$passed)
  expect_equal(values(fit), values(lobpcg_fit), tolerance = 1e-8)
})

test_that("native generalized LOBPCG handles repeated clusters and magnitude targets", {
  A <- Matrix::Diagonal(x = c(1, 1, 4, 4, 9, 16))
  B <- Matrix::Diagonal(x = c(1, 1, 2, 2, 3, 4))
  cases <- list(
    list(target = smallest(), expected = c(1, 1)),
    list(target = largest_magnitude(), expected = c(4, 3))
  )

  for (case in cases) {
    fit <- eig_partial(
      A,
      B = B,
      k = 2,
      target = case$target,
      method = lobpcg(maxit = 120L),
      seed = 326,
      tol = 1e-8,
      allow_dense_fallback = "never"
    )

    expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
    expect_true(fit$restart$native)
    expect_true(fit$restart$native_kernels)
    expect_true(fit$restart$generalized)
    expect_true(certificate(fit)$passed)
    expect_equal(values(fit), case$expected, tolerance = 1e-8)
    expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2),
                 tolerance = 1e-8)
  }
})

test_that("explicit generalized Lanczos supports native dense SPD metric transforms", {
  A <- diag(c(2, 5, 9, 14))
  B <- diag(c(1, 2, 3, 4))
  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = largest(),
    method = lanczos(max_subspace = 4L),
    seed = 322,
    tol = 1e-8
  )

  expect_equal(fit$method, eigencore:::native_generalized_lanczos_label())
  expect_true(fit$restart$native)
  expect_equal(fit$restart$metric_solve, "dense Cholesky similarity transform for B")
  expect_equal(fit$restart$transformed_operator_storage, "dense")
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), c(3.5, 3), tolerance = 1e-8)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(2), tolerance = 1e-8)
})

test_that("explicit block generalized Lanczos uses native transformed diagonal-B path", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))
  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = lanczos(max_subspace = 6L, block = 2L),
    seed = 324,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lanczos_label())
  expect_equal(fit$restart$kind, "native_transformed_generalized_b_orthogonal_lanczos")
  expect_true(fit$restart$native)
  expect_equal(fit$restart$block, 2L)
  expect_equal(fit$restart$metric_solve, "diagonal scaling similarity transform for B")
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), c(1, 2), tolerance = 1e-8)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)
})

test_that("explicit generalized Lanczos supports sparse CSC SPD metric solves", {
  n <- 10L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), rep(2.5, n), rep(-1, n - 1L))
  )
  B <- methods::as(Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-0.05, n - 1L),
      seq(1.2, 2, length.out = n),
      rep(-0.05, n - 1L)
    )
  ), "dgCMatrix")

  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = lanczos(max_subspace = 10L),
    seed = 323,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::generalized_lanczos_label())
  expect_equal(fit$restart$metric_solve, "native sparse tridiagonal Thomas solve for B")
  expect_equal(fit$restart$metric_solve_kind, "native_sparse_tridiagonal_thomas")
  expect_equal(fit$restart$metric_factorization, "tridiagonal_thomas")
  expect_true(fit$restart$metric_solve_native)
  expect_true(fit$restart$native_metric_solve)
  expect_gt(fit$restart$metric_solves, 0L)
  expect_false(fit$restart$native)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)

  lobpcg_fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = lobpcg(maxit = 120L),
    seed = 323,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )
  expect_true(certificate(lobpcg_fit)$passed)
  expect_equal(values(fit), values(lobpcg_fit), tolerance = 1e-8)
})

test_that("dsCMatrix generalized SPD inputs keep native sparse eligibility", {
  n <- 14L
  A <- methods::as(Matrix::forceSymmetric(Matrix::bandSparse(
    n,
    k = c(0, 1),
    diagonals = list(rep(2.5, n), rep(-1, n - 1L))
  )), "CsparseMatrix")
  B <- methods::as(Matrix::forceSymmetric(Matrix::bandSparse(
    n,
    k = c(0, 1),
    diagonals = list(seq(1.5, 2.4, length.out = n), rep(0.03, n - 1L))
  )), "CsparseMatrix")
  target <- smallest_magnitude()
  Aop <- as_operator(A)
  Bop <- as_operator(B)

  expect_s4_class(A, "dsCMatrix")
  expect_s4_class(B, "dsCMatrix")
  expect_equal(Aop$metadata$storage, "dgCMatrix")
  expect_equal(Bop$metadata$storage, "dgCMatrix")
  expect_equal(Aop$metadata$input_storage, "dsCMatrix")
  expect_equal(Bop$metadata$input_storage, "dsCMatrix")
  expect_true(Aop$metadata$symmetric_storage)
  expect_true(Bop$metadata$symmetric_storage)
  expect_lt(length(methods::slot(Aop$metadata$matrix, "x")), n * n / 2)

  problem <- eigen_problem(A, metric = B, target = target)
  plan <- plan_solver(problem, k = 2L, method = lobpcg(maxit = 180L))
  expect_equal(plan$method, eigencore:::native_generalized_lobpcg_label())

  fit <- eig_partial(
    A,
    B = B,
    k = 2L,
    target = target,
    method = lobpcg(maxit = 180L),
    seed = 325,
    tol = 1e-7,
    allow_dense_fallback = "never"
  )
  oracle <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
  idx <- eigencore:::order_indices(oracle$values, target)[seq_len(2L)]
  expected <- oracle$values[idx]

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$native)
  expect_true(fit$restart$native_kernels)
  expect_equal(fit$restart$orthogonalization_methods, "native_csc_b_mgs2")
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), expected, tolerance = 1e-5)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2),
               tolerance = 1e-6)
})

test_that("general sparse CSC SPD metric solves keep reference Cholesky provenance", {
  n <- 10L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), rep(2.5, n), rep(-1, n - 1L))
  )
  B <- methods::as(
    Matrix::Diagonal(x = seq(1.5, 2.4, length.out = n)) +
      Matrix::bandSparse(
        n,
        k = c(-2, 2),
        diagonals = list(rep(0.04, n - 2L), rep(0.04, n - 2L))
      ),
    "dgCMatrix"
  )

  fit <- eig_partial(
    A,
    B = B,
    k = 2,
    target = smallest(),
    method = lanczos(max_subspace = 10L),
    seed = 324,
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::generalized_lanczos_label())
  expect_equal(fit$restart$metric_solve, "sparse Cholesky solve for B")
  expect_equal(fit$restart$metric_solve_kind, "sparse_cholesky")
  expect_equal(fit$restart$metric_factorization, "Matrix::Cholesky")
  expect_false(fit$restart$metric_solve_native)
  expect_false(fit$restart$native_metric_solve)
  expect_gt(fit$restart$metric_solves, 0L)
  expect_false(fit$restart$native)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), as.matrix(B) %*% vectors(fit)), diag(2), tolerance = 1e-8)
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

test_that("native shifted-tridiagonal generalized LOBPCG certifies largest sparse targets", {
  n <- 40L
  A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n - 1L), seq_len(n - 1L) + 1L),
    j = c(seq_len(n), seq_len(n - 1L) + 1L, seq_len(n - 1L)),
    x = c(rep(2.25, n), rep(-1, 2L * (n - 1L))),
    dims = c(n, n)
  )
  B <- Matrix::Diagonal(n, seq(0.5, 2, length.out = n))
  preconditioner <- shifted_tridiagonal_preconditioner(A, shift = 4)
  fit <- eig_partial(
    A,
    B = B,
    k = 3L,
    target = largest(),
    method = lobpcg(maxit = 160L, preconditioner = preconditioner),
    tol = 1e-8,
    seed = 77,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
  expect_equal(fit$restart$preconditioner_kind, "shifted_tridiagonal")
  expect_gt(fit$restart$preconditioner_calls, 0L)
  expect_certificate_clean(fit)
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

test_that("adversarial generalized B bank stays native for largest and smallest targets", {
  csc_n <- 10L
  csc_A <- Matrix::bandSparse(
    csc_n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, csc_n - 1L), rep(2.5, csc_n), rep(-1, csc_n - 1L))
  )
  csc_B <- methods::as(Matrix::bandSparse(
    csc_n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-0.05, csc_n - 1L),
      seq(1.2, 2, length.out = csc_n),
      rep(-0.05, csc_n - 1L)
    )
  ), "dgCMatrix")

  mf_n <- 9L
  mf_bdiag <- seq(0.75, 2.5, length.out = mf_n)
  mf_B_dense <- diag(mf_bdiag)
  mf_B <- linear_operator(
    dim = c(mf_n, mf_n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (mf_bdiag * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (mf_bdiag * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    metadata = list(
      frobenius_norm = sqrt(sum(mf_bdiag^2)),
      positive_definite = TRUE
    )
  )

  cases <- list(
    list(
      name = "ill_conditioned_diagonal_b",
      A = Matrix::Diagonal(x = c(1, 3, 7, 12, 20, 33, 50, 80)),
      B = Matrix::Diagonal(x = 10^seq(-3, 3, length.out = 8L)),
      B_dense = diag(10^seq(-3, 3, length.out = 8L)),
      orthogonalization = "native_diagonal_b_mgs2",
      seed = 901L
    ),
    list(
      name = "sparse_csc_b",
      A = csc_A,
      B = csc_B,
      B_dense = as.matrix(csc_B),
      orthogonalization = "native_csc_b_mgs2",
      seed = 902L
    ),
    list(
      name = "explicit_spd_matrix_free_b",
      A = Matrix::Diagonal(x = seq(2, 18, by = 2)),
      B = mf_B,
      B_dense = mf_B_dense,
      orthogonalization = "native_matrix_free_b_mgs2",
      seed = 903L
    )
  )

  for (case in cases) {
    for (target in list(smallest(), largest())) {
      problem <- eigen_problem(case$A, metric = case$B, target = target)
      plan <- plan_solver(problem, k = 2L, method = lobpcg(maxit = 220L))
      expect_equal(plan$method, eigencore:::native_generalized_lobpcg_label())

      fit <- eig_partial(
        case$A,
        B = case$B,
        k = 2L,
        target = target,
        method = lobpcg(maxit = 220L),
        seed = case$seed,
        tol = 1e-8,
        allow_dense_fallback = "never"
      )
      oracle <- eigencore:::dense_generalized_spd_eigen(
        as.matrix(case$A),
        case$B_dense
      )
      idx <- eigencore:::order_indices(oracle$values, target)[seq_len(2L)]
      expected <- oracle$values[idx]
      relative_value_error <- max(
        abs(values(fit) - expected) / pmax(1, abs(expected))
      )

      expect_equal(fit$method, eigencore:::native_generalized_lobpcg_label())
      expect_true(fit$restart$native)
      expect_true(fit$restart$native_kernels)
      expect_true(fit$restart$generalized)
      expect_equal(fit$restart$orthogonalization_methods, case$orthogonalization)
      expect_equal(fit$nconv, 2L)
      expect_true(certificate(fit)$passed)
      expect_lte(certificate(fit)$max_backward_error, 1e-8)
      expect_lte(relative_value_error, 1e-6)
      expect_equal(
        crossprod(vectors(fit), case$B_dense %*% vectors(fit)),
        diag(2L),
        tolerance = 1e-8
      )
    }
  }
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
