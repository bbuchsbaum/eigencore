test_that("shift-invert returns interior eigenvalues nearest sigma on dense Hermitian", {
  set.seed(42)
  n <- 20L
  vals <- c(-3, -1, 0.5, 1.5, 2.7, 4, seq(5, 30, length.out = n - 6L))
  A    <- symmetric_with_spectrum(vals, seed = 42)

  fit <- eig_partial(A, k = 3L, target = nearest(2),
                     method = shift_invert(sigma = 2))

  expect_identical(fit$method,
                   eigencore:::native_dense_shift_invert_label())
  expected <- vals[order(abs(vals - 2))][1:3]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  direct <- eigencore:::certify_eigen_operator(as_operator(A), fit$values,
                                               fit$vectors, tol = 1e-8)
  expect_equal(fit$certificate$residuals, direct$residuals, tolerance = 1e-12)
  expect_equal(fit$certificate$backward_error, direct$backward_error,
               tolerance = 1e-12)
  expect_identical(fit$transform$kind, "shift_invert")
  expect_identical(fit$transform$label_kind, "dense_lu_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_identical(fit$transform$factorization_cache$factorization,
                   "LAPACK dgetrf/dgetrs")
  expect_identical(fit$restart$kind, "native_dense_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_identical(fit$transform$certification$problem, "original")
  expect_false(fit$transform$certification$transformed_residuals_used)
  expect_true(fit$plan$controls$certified_in_original_coordinates)
  expect_identical(fit$plan$controls$transformed_operator_target,
                   "largest_magnitude")
  expect_equal(fit$sigma, 2)
})

test_that("shift-invert handles a sparse CSC source via factorized solve", {
  set.seed(7)
  vals <- seq(1, 30)
  A_csc <- methods::as(symmetric_with_spectrum(vals, seed = 7),
                       "CsparseMatrix")

  fit <- eig_partial(A_csc, k = 4L, target = nearest(15.5),
                     method = shift_invert(sigma = 15.5))

  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (sparse LU)")
  expected <- vals[order(abs(vals - 15.5))][1:4]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  # Sparse paths use a Hutchinson norm estimate — all converged but
  # passed is withheld by design.
  expect_true(all(fit$certificate$converged))
  expect_true(fit$certificate$scale_is_estimate)
  expect_lt(max(fit$certificate$backward_error), 1e-7)
  cache <- fit$transform$factorization_cache
  expect_equal(cache$label_kind, "sparse_lu")
  expect_equal(cache$factorization, "Matrix::lu")
  expect_true(cache$factorization_cached)
  expect_equal(cache$condition_estimate_type, "sparse_lu_pivot_ratio")
  expect_true(is.finite(cache$condition_estimate))
  expect_gt(cache$condition_estimate, 0)
  expect_false(isTRUE(cache$near_singular))
})

test_that("shift-invert uses native tridiagonal factorized Lanczos for sparse CSC paths", {
  n <- 50L
  A_csc <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  sigma <- 0.01

  fit <- eig_partial(A_csc, k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-8,
                     allow_dense_fallback = "never")
  expected <- eigen(as.matrix(A_csc), symmetric = TRUE, only.values = TRUE)$values
  expected <- expected[order(abs(expected - sigma))][1:3]

  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$label_kind, "tridiagonal_thomas_native")
  cache <- fit$transform$factorization_cache
  expect_true(cache$native)
  expect_equal(cache$factorization, "native tridiagonal Thomas")
  expect_equal(cache$condition_estimate_type, "tridiagonal_thomas_pivot_ratio")
  expect_true(is.finite(cache$condition_estimate))
  expect_gt(cache$condition_estimate, 0)
  expect_identical(fit$restart$kind, "native_tridiagonal_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$factorization_native)
  expect_identical(fit$transform$certification$problem, "original")
  expect_false(fit$transform$certification$transformed_residuals_used)
})

test_that("shift-invert uses native tridiagonal factorized Lanczos for diagonal sources", {
  A <- Matrix::Diagonal(x = c(1, 2, 4, 8, 16, 32))
  sigma <- 7

  fit <- eig_partial(A, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-10,
                     allow_dense_fallback = "never")

  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_equal(sort(fit$values), c(4, 8), tolerance = 1e-10)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$label_kind, "tridiagonal_thomas_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_identical(fit$restart$kind, "native_tridiagonal_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$factorization_native)
})

test_that("shift-invert accepts a user-supplied solve operator", {
  set.seed(11)
  n <- 15L
  vals <- seq(0.1, 14.5, length.out = n)
  A    <- symmetric_with_spectrum(vals, seed = 11)
  sigma <- 7.4   # not an eigenvalue
  M     <- A - sigma * diag(n)
  user_solve <- function(X) base::solve(M, X)

  fit <- eig_partial(A, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma, solve = user_solve))

  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (user solve)")
  expected <- vals[order(abs(vals - sigma))][1:2]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_equal(fit$transform$factorization_cache$factorization, "user_solve")
  expect_true(fit$transform$factorization_cache$external_cache)
})

test_that("shift-invert refuses ignored public factorization handles", {
  A <- diag(c(1, 2, 4, 8))
  expect_error(
    eig_partial(A, k = 1L, target = nearest(2.5),
                method = shift_invert(sigma = 2.5, factorization = list())),
    "factorization = .*not implemented"
  )
})

test_that("generalized shift-invert certifies original generalized residuals on dense SPD problems", {
  set.seed(9)
  n <- 10L
  A <- symmetric_with_spectrum(seq(1, n), seed = 9)
  B <- diag(seq(1, 2, length.out = n))
  sigma <- 3.2

  fit <- eig_partial(A, B = B, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-9)
  oracle <- eigencore:::dense_generalized_spd_eigen(A, B)
  expected <- oracle$values[order(abs(oracle$values - sigma))][1:2]

  expect_identical(
    fit$method,
    eigencore:::native_dense_generalized_shift_invert_label()
  )
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  direct <- eigencore:::certify_eigen_operator(
    as_operator(A), fit$values, fit$vectors, Bop = as_operator(B), tol = 1e-9
  )
  expect_equal(fit$certificate$residuals, direct$residuals, tolerance = 1e-12)
  expect_equal(fit$certificate$backward_error, direct$backward_error,
               tolerance = 1e-12)
  expect_equal(crossprod(fit$vectors, B %*% fit$vectors), diag(2),
               tolerance = 1e-8)
  expect_identical(fit$transform$certification$problem, "original")
  expect_identical(fit$transform$certification$residual_formula,
                   "A * x - lambda * B * x")
  expect_false(fit$transform$certification$transformed_residuals_used)
  expect_true(fit$plan$controls$certified_in_original_coordinates)
  expect_equal(fit$transform$factorization_cache$label_kind,
               "dense_lu_generalized_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_true(fit$transform$factorization_cache$generalized)
  expect_equal(fit$transform$factorization_cache$metric_factorization,
               "LAPACK dpotrf(B)")
  expect_identical(fit$restart$kind,
                   "native_dense_generalized_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
})

test_that("generalized shift-invert handles sparse A with diagonal B without densifying B", {
  n <- 30L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), rep(2.5, n), rep(-1, n - 1))
  )
  B <- Matrix::Diagonal(x = seq(1, 2, length.out = n))
  oracle <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
  sigma <- sort(oracle$values)[4L] + 0.01

  fit <- eig_partial(A, B = B, k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-8,
                     allow_dense_fallback = "never")
  expected <- oracle$values[order(abs(oracle$values - sigma))][1:3]

  expect_identical(
    fit$method,
    eigencore:::native_tridiagonal_generalized_shift_invert_label()
  )
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-6)
  expect_certificate_clean(fit)
  expect_equal(crossprod(fit$vectors, as.matrix(B) %*% fit$vectors), diag(3),
               tolerance = 1e-8)
  expect_equal(fit$transform$factorization_cache$label_kind,
               "tridiagonal_thomas_generalized_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_equal(fit$transform$factorization_cache$factorization,
               "native tridiagonal Thomas + diagonal sqrt(B)")
  expect_equal(fit$transform$factorization_cache$metric_factorization,
               "diagonal sqrt(B)")
  expect_identical(fit$restart$kind,
                   "native_tridiagonal_generalized_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$factorization_native)
})

test_that("generalized shift-invert keeps non-tridiagonal sparse A reference-labelled", {
  set.seed(10)
  A <- methods::as(symmetric_with_spectrum(seq(1, 12), seed = 10), "CsparseMatrix")
  B <- Matrix::Diagonal(x = seq(1.1, 1.8, length.out = 12L))
  oracle <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
  sigma <- oracle$values[[5L]] + 0.02

  fit <- eig_partial(A, B = B, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-8,
                     allow_dense_fallback = "never")
  expected <- oracle$values[order(abs(oracle$values - sigma))][1:2]

  expect_identical(
    fit$method,
    "reference generalized SPD Lanczos shift-invert (sparse LU)"
  )
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-6)
  expect_equal(fit$transform$factorization_cache$label_kind,
               "sparse_lu_generalized")
  expect_false(isTRUE(fit$transform$factorization_cache$native))
})

test_that("generalized shift-invert rejects unsupported sparse B without densifying", {
  n <- 8L
  A <- Matrix::Diagonal(x = seq_len(n) + 1)
  B <- methods::as(Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-0.05, n - 1L),
      rep(1.5, n),
      rep(-0.05, n - 1L)
    )
  ), "dgCMatrix")

  expect_error(
    eig_partial(A, B = B, k = 2L, target = nearest(3),
                method = shift_invert(sigma = 3),
                allow_dense_fallback = "never"),
    "dense B or diagonal B|avoid silent densification"
  )
})

test_that("planner emits an honest shift-invert label for the chosen path", {
  set.seed(13)
  vals <- seq(1, 12)
  A    <- symmetric_with_spectrum(vals, seed = 13)
  P    <- eigen_problem(A, target = nearest(5),
                        transform = shift_invert(sigma = 5))
  plan <- plan_solver(P, k = 3, method = shift_invert(sigma = 5))
  expect_identical(plan$method,
                   eigencore:::native_dense_shift_invert_label())
  expect_true(plan$controls$certified_in_original_coordinates)
  expect_identical(plan$controls$eigenvalue_recovery,
                   "lambda = sigma + 1 / mu")
  expect_match(plan$controls$certification_policy, "original eigenproblem")
  expect_true(any(grepl("shift-invert transform requested", plan$reasons,
                        fixed = TRUE)))
})

test_that("shift-invert factorization cache keys invalidate on A, sigma, and B changes", {
  A <- diag(c(1, 2, 4, 8))
  B <- diag(c(1, 2, 3, 4))
  A_changed <- A
  A_changed[2, 2] <- 2.5
  B_changed <- B
  B_changed[3, 3] <- 3.5

  key <- eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.25)
  expect_identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.25)
  )
  expect_false(identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(as_operator(A_changed), 0.25)
  ))
  expect_false(identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.5)
  ))
  expect_false(identical(
    eigencore:::shift_invert_factorization_cache_key(
      as_operator(A), 0.25, as_operator(B)
    ),
    eigencore:::shift_invert_factorization_cache_key(
      as_operator(A), 0.25, as_operator(B_changed)
    )
  ))
  expect_true(key$standard_problem)
  expect_null(key$B)
  expect_identical(key$structure, "hermitian")

  general_op <- as_operator(A)
  general_op$structure <- general()
  expect_false(identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(general_op, 0.25)
  ))
})

test_that("shift-invert result exposes factorization-cache provenance", {
  vals <- seq(1, 8)
  A <- symmetric_with_spectrum(vals, seed = 91)

  fit <- eig_partial(A, k = 2L, target = nearest(4.2),
                     method = shift_invert(sigma = 4.2))
  cache <- fit$transform$factorization_cache

  expect_s3_class(cache$key, "eigencore_shift_invert_cache_key")
  expect_equal(cache$key$sigma, 4.2)
  expect_equal(cache$label_kind, "dense_lu_native")
  expect_equal(cache$factorization, "LAPACK dgetrf/dgetrs")
  expect_true(cache$factorization_cached)
  expect_equal(cache$condition_estimate_type, "dense_lu_pivot_ratio")
  expect_true(is.finite(cache$condition_estimate))
  expect_true(cache$native)
  expect_true(cache$reusable_within_operator)
  expect_false(cache$external_cache)
})

test_that("shift-invert near a true eigenvalue surfaces a clear error", {
  set.seed(17)
  vals <- c(1, 2, 3, 4, 5)
  A    <- symmetric_with_spectrum(vals, seed = 17)
  # sigma equal to an eigenvalue makes (A - sigma I) singular; base::solve errors
  expect_error(
    eig_partial(A, k = 2L, target = nearest(3),
                method = shift_invert(sigma = 3)),
    "singular|computationally"
  )
})

test_that("shift-invert recovers smallest eigenvalues of a 1D Laplacian", {
  # discrete 1D Laplacian -- the canonical Milestone L test case
  n <- 50L
  L <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1))
  )
  fit <- eig_partial(L, k = 4L, target = nearest(0.01),
                     method = shift_invert(sigma = 0.01))
  oracle <- sort(eigen(as.matrix(L), symmetric = TRUE)$values)[1:4]
  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_equal(sort(fit$values), sort(oracle), tolerance = 1e-6)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$label_kind, "tridiagonal_thomas_native")
  expect_true(fit$restart$native)
  expect_true(fit$restart$factorization_native)
})
