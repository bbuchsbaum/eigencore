test_that("shift-invert returns interior eigenvalues nearest sigma on dense Hermitian", {
  set.seed(42)
  n <- 20L
  vals <- c(-3, -1, 0.5, 1.5, 2.7, 4, seq(5, 30, length.out = n - 6L))
  A    <- symmetric_with_spectrum(vals, seed = 42)

  fit <- eig_partial(A, k = 3L, target = nearest(2),
                     method = shift_invert(sigma = 2))

  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (dense LU)")
  expected <- vals[order(abs(vals - 2))][1:3]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$kind, "shift_invert")
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
                   "reference Hermitian Lanczos shift-invert (dense LU)")
  expect_true(any(grepl("shift-invert transform requested", plan$reasons,
                        fixed = TRUE)))
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
