# SVD target completeness (R/completeness_svd.R): certificate$passed requires
# the returned triplets to be the requested set.

svd_completeness_diag <- function() {
  c(rep(9, 3), rep(7, 3), 5, 4, 3, rep(1, 3), rep(0.5, 3), 0.2)
}

test_that("dense LAPACK and explicit Gram routes are exact", {
  set.seed(31)
  X <- matrix(rnorm(400), 40, 10)
  fit <- svd_partial(X, 3, target = smallest())
  expect_identical(fit$certificate$target_completeness, "exact")
  expect_true(fit$certificate$passed)
  fit <- svd_partial(X, 3)
  expect_identical(fit$certificate$target_completeness, "exact")
  expect_true(fit$certificate$passed)
  S <- Matrix::rsparsematrix(200, 40, 0.1)
  fit <- svd_partial(S, 4)
  expect_identical(fit$method, "native certified Gram SVD special case")
  expect_identical(fit$certificate$target_completeness, "exact")
  expect_true(fit$certificate$passed)
})

test_that("Golub-Kahan on repeated singular values is repaired or proved complete", {
  d <- svd_completeness_diag()
  D <- Matrix::Diagonal(x = d)
  for (mode in c("auto", "probe", "inertia")) {
    withr::local_options(eigencore.target_completeness = mode)
    fit <- svd_partial(D, 5, method = golub_kahan())
    expect_true(fit$certificate$passed)
    expect_true(fit$certificate$target_completeness %in%
                  c("repaired", "probed", "inertia_verified"))
    expect_equal(fit$d, c(9, 9, 9, 7, 7), tolerance = 1e-8)
    fit <- svd_partial(D, 4, target = smallest(), method = golub_kahan())
    expect_true(fit$certificate$passed)
    expect_equal(sort(fit$d), c(0.2, 0.5, 0.5, 0.5), tolerance = 1e-8)
    expect_lte(fit$certificate$max_backward_error, 1e-8)
  }
})

test_that("the Gram probe repairs a missing copy and re-certifies it", {
  withr::local_options(eigencore.target_completeness = "probe")
  D <- Matrix::Diagonal(x = svd_completeness_diag())
  fit <- svd_partial(D, 5, method = golub_kahan())
  cert <- fit$certificate
  expect_identical(cert$target_completeness, "repaired")
  expect_true(cert$completeness$intruder_found)
  expect_identical(cert$completeness$method, "gram_deflated_complement_probe")
  expect_true(cert$passed)
  expect_true(any(grepl("repaired", fit$warnings)))
  # The repaired triplets are certified from scratch.
  U <- as.matrix(fit$u)
  V <- as.matrix(fit$v)
  expect_lt(max(abs(as.matrix(D %*% V) - sweep(U, 2L, fit$d, `*`))), 1e-8)
  expect_lt(max(abs(crossprod(V) - diag(5))), 1e-8)
})

test_that("a short Golub-Kahan result is completed from the Gram complement", {
  D <- Matrix::Diagonal(x = rep(c(3, 2, 1), each = 4))
  withr::local_options(eigencore.target_completeness = "none")
  short <- svd_partial(D, 6, method = golub_kahan())
  expect_lt(length(short$d), 6L)
  expect_false(short$certificate$passed)
  withr::local_options(eigencore.target_completeness = "auto")
  fit <- svd_partial(D, 6, method = golub_kahan())
  expect_true(fit$certificate$passed)
  expect_equal(fit$d, c(3, 3, 3, 3, 2, 2), tolerance = 1e-8)
  expect_gt(fit$certificate$completeness$filled, 0L)
  fit <- svd_partial(D, 6, target = smallest(), method = golub_kahan())
  expect_true(fit$certificate$passed)
  expect_equal(sort(fit$d), c(1, 1, 1, 1, 2, 2), tolerance = 1e-8)
  # A complement too large for the dense path uses block Lanczos.
  big <- Matrix::Diagonal(x = rep(c(5, 4, 3, 2, 1), each = 30))
  withr::local_options(eigencore.target_completeness = "probe")
  fit <- svd_partial(big, 8, method = golub_kahan())
  expect_true(fit$certificate$passed)
  expect_equal(fit$d, rep(5, 8), tolerance = 1e-8)
})

test_that("matrix-free and centred operators are probed", {
  D <- Matrix::Diagonal(x = svd_completeness_diag())
  op <- linear_operator(
    dim(D),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) alpha * as.matrix(D %*% X),
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) alpha * as.matrix(D %*% X)
  )
  fit <- svd_partial(op, 5)
  expect_true(fit$certificate$passed)
  expect_true(fit$certificate$target_completeness %in% c("probed", "repaired"))
  expect_equal(fit$d, c(9, 9, 9, 7, 7), tolerance = 1e-8)
  set.seed(4)
  S <- Matrix::rsparsematrix(150, 30, 0.2)
  fit <- svd_partial(center(S), 3)
  expect_true(fit$certificate$passed)
  expect_true(fit$certificate$target_completeness %in%
                c("probed", "repaired", "exact", "inertia_verified"))
  oracle <- svd(sweep(as.matrix(S), 2L, Matrix::colMeans(S)))$d[1:3]
  expect_equal(fit$d, oracle, tolerance = 1e-8)
})

test_that("an inertia count on the augmented matrix proves completeness", {
  set.seed(8)
  S <- Matrix::rsparsematrix(120, 40, 0.15)
  withr::local_options(eigencore.target_completeness = "inertia")
  fit <- svd_partial(S, 4, method = golub_kahan())
  cert <- fit$certificate
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_identical(cert$completeness$method, "augmented_inertia_count")
  expect_identical(cert$completeness$count_upper, 4)
  expect_true(cert$passed)
  fit <- svd_partial(S, 3, target = nearest(1))
  expect_true(fit$certificate$target_completeness %in% c("inertia_verified", "exact"))
  expect_true(fit$certificate$passed)
  oracle <- svd(as.matrix(S))$d
  expect_equal(sort(fit$d), sort(oracle[order(abs(oracle - 1))][1:3]), tolerance = 1e-8)
})

test_that("an interior SVD that misses a copy reports inertia_failed", {
  D <- Matrix::Diagonal(x = svd_completeness_diag())
  fit <- svd_partial(D, 3, target = nearest(0.9))
  if (!isTRUE(all.equal(sort(fit$d), c(1, 1, 1), tolerance = 1e-8))) {
    expect_false(fit$certificate$passed)
    expect_true(fit$certificate$residual_passed)
    expect_identical(fit$certificate$target_completeness, "inertia_failed")
  } else {
    expect_true(fit$certificate$passed)
  }
})

test_that("svds() carries the SVD completeness status, also for nu = 0 / nv = 0", {
  set.seed(9)
  S <- Matrix::rsparsematrix(80, 30, 0.2)
  for (nunv in list(c(3L, 3L), c(0L, 0L), c(0L, 3L), c(3L, 0L))) {
    fit <- svds(S, 3, nu = nunv[1], nv = nunv[2], opts = list(tol = 1e-10))
    expect_true(fit$certificate$passed)
    expect_true(fit$certificate$target_completeness %in%
                  eigencore:::verified_completeness_states())
  }
  fit <- svds(S, 3, opts = list(center = TRUE, scale = TRUE))
  expect_true(fit$certificate$passed)
})

test_that("completeness = none and require_completeness = FALSE behave as documented", {
  D <- Matrix::Diagonal(x = svd_completeness_diag())
  withr::local_options(eigencore.target_completeness = "none")
  fit <- svd_partial(D, 5, method = golub_kahan())
  expect_identical(fit$certificate$target_completeness, "not_checked")
  expect_true(fit$certificate$residual_passed)
  expect_false(fit$certificate$passed)
  withr::local_options(eigencore.require_completeness = FALSE)
  fit <- svd_partial(D, 5, method = golub_kahan())
  expect_true(fit$certificate$passed)
})

test_that("the SVD completeness probe leaves the global RNG stream alone", {
  D <- Matrix::Diagonal(x = svd_completeness_diag())
  withr::local_options(eigencore.target_completeness = "none")
  set.seed(77)
  fit <- svd_partial(D, 5, method = golub_kahan())
  unchecked <- .Random.seed
  withr::local_options(eigencore.target_completeness = "probe")
  set.seed(77)
  fit <- svd_partial(D, 5, method = golub_kahan())
  expect_identical(fit$certificate$target_completeness, "repaired")
  expect_identical(.Random.seed, unchecked)
})
