test_that("shift-invert returns interior eigenvalues nearest sigma on dense Hermitian", {
  set.seed(42)
  n <- 20L
  vals <- c(-3, -1, 0.5, 1.5, 2.7, 4, seq(5, 30, length.out = n - 6L))
  A    <- symmetric_with_spectrum(vals, seed = 42)

  fit <- eig_partial(A, k = 3L, target = nearest(2),
                     method = shift_invert(sigma = 2))

  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (dense QR)")
  expected <- vals[order(abs(vals - 2))][1:3]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  direct <- eigencore:::certify_eigen_operator(as_operator(A), fit$values,
                                               fit$vectors, tol = 1e-8)
  expect_equal(fit$certificate$residuals, direct$residuals, tolerance = 1e-12)
  expect_equal(fit$certificate$backward_error, direct$backward_error,
               tolerance = 1e-12)
  expect_identical(fit$transform$kind, "shift_invert")
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
  expect_equal(cache$condition_estimate_type, "uncomputed_sparse_no_dense_rcond")
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

test_that("generalized shift-invert is rejected at plan time with a roadmap note", {
  set.seed(9)
  n <- 10L
  A <- symmetric_with_spectrum(seq(1, n), seed = 9)
  B <- diag(seq(1, 2, length.out = n))
  expect_error(
    eig_partial(A, B = B, k = 2L, target = nearest(3),
                method = shift_invert(sigma = 3)),
    "[Gg]eneralized SPD shift-invert"
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
                   "reference Hermitian Lanczos shift-invert (dense QR)")
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
})

test_that("shift-invert result exposes factorization-cache provenance", {
  vals <- seq(1, 8)
  A <- symmetric_with_spectrum(vals, seed = 91)

  fit <- eig_partial(A, k = 2L, target = nearest(4.2),
                     method = shift_invert(sigma = 4.2))
  cache <- fit$transform$factorization_cache

  expect_s3_class(cache$key, "eigencore_shift_invert_cache_key")
  expect_equal(cache$key$sigma, 4.2)
  expect_equal(cache$label_kind, "dense_qr")
  expect_equal(cache$factorization, "base::qr")
  expect_true(cache$factorization_cached)
  expect_equal(cache$condition_estimate_type, "dense_rcond")
  expect_true(is.finite(cache$condition_estimate))
  expect_false(cache$native)
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
                   "reference Hermitian Lanczos shift-invert (sparse LU)")
  expect_equal(sort(fit$values), sort(oracle), tolerance = 1e-6)
  expect_true(all(fit$certificate$converged))
})
