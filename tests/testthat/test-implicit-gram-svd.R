test_that("implicit Gram SVD certifies dense partial SVD beyond the explicit Gram cap", {
  set.seed(101)
  A <- matrix(rnorm(600 * 550), 600, 550)

  plan <- plan_solver(svd_problem(A), rank = 6)
  expect_identical(
    plan$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )

  fit <- svd_partial(A, rank = 6, tol = 1e-8)
  expect_identical(
    fit$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )
  expect_true(fit$certificate$passed)
  expect_false(fit$restart$materialized_gram)
  expect_true(fit$restart$normal_operator_implicit)
  expect_true(fit$restart$certified_in_original_coordinates)

  ref <- svd(A, nu = 6, nv = 6)
  expect_equal(fit$d, ref$d[1:6], tolerance = 1e-7)
  # singular vectors up to sign
  agreement_u <- abs(colSums(fit$u * ref$u))
  agreement_v <- abs(colSums(fit$v * ref$v))
  expect_true(all(agreement_u > 1 - 1e-6))
  expect_true(all(agreement_v > 1 - 1e-6))
})

test_that("implicit Gram SVD certifies sparse partial SVD beyond the explicit Gram cap", {
  set.seed(102)
  A <- Matrix::rsparsematrix(4000, 900, density = 0.01)

  plan <- plan_solver(svd_problem(A), rank = 8)
  expect_identical(
    plan$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )

  fit <- svd_partial(A, rank = 8, tol = 1e-8)
  expect_true(fit$certificate$passed)
  expect_identical(dim(fit$u), c(4000L, 8L))
  expect_identical(dim(fit$v), c(900L, 8L))

  gk_ref <- svd_partial(A, rank = 8, tol = 1e-8, method = golub_kahan())
  expect_equal(fit$d, gk_ref$d, tolerance = 1e-7)
})

test_that("implicit Gram SVD handles wide operators via the left normal side", {
  set.seed(103)
  # small side above the wide explicit-Gram cap (1024) so the implicit path owns it
  A <- Matrix::rsparsematrix(1100, 6000, density = 0.01)

  fit <- svd_partial(A, rank = 5, tol = 1e-8)
  expect_identical(
    fit$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )
  expect_true(fit$certificate$passed)
  expect_identical(fit$restart$gram_side, "left")

  gk_ref <- svd_partial(A, rank = 5, tol = 1e-8, method = golub_kahan())
  expect_equal(fit$d, gk_ref$d, tolerance = 1e-7)
})

test_that("implicit Gram SVD respects vector selection modes", {
  set.seed(104)
  A <- matrix(rnorm(400 * 200), 400, 200)

  fit_left <- svd_partial(A, rank = 3, vectors = "left")
  expect_false(is.null(fit_left$u))
  expect_null(fit_left$v)

  fit_none <- svd_partial(A, rank = 3, vectors = "none")
  expect_null(fit_none$u)
  expect_null(fit_none$v)
})

test_that("implicit Gram SVD policy leaves explicit Gram and small problems alone", {
  set.seed(105)
  # small side within the explicit Gram cap and aspect ratio satisfied:
  # explicit Gram keeps priority
  A <- Matrix::rsparsematrix(2000, 300, density = 0.01)
  plan <- plan_solver(svd_problem(A), rank = 5)
  expect_identical(plan$method, "native certified Gram SVD special case")

  # tiny problems stay on their existing paths
  B <- matrix(rnorm(40 * 30), 40, 30)
  plan_small <- plan_solver(svd_problem(B), rank = 3)
  expect_false(identical(
    plan_small$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  ))

  # smallest targets are not captured by the implicit Gram path
  plan_smallest <- plan_solver(svd_problem(A), rank = 3, method = auto())
  expect_identical(plan_smallest$method, "native certified Gram SVD special case")

  # explicit user method choices are honored
  plan_gk <- plan_solver(svd_problem(A), rank = 5, method = golub_kahan())
  expect_identical(plan_gk$method, "native prototype Golub-Kahan")
})

test_that("implicit Gram SVD result matches Golub-Kahan on a fixed spectrum", {
  set.seed(106)
  m <- 500L; n <- 400L
  d_true <- c(50, 40, 30, 20, 10, rep(1, n - 5))
  U0 <- qr.Q(qr(matrix(rnorm(m * n), m, n)))
  V0 <- qr.Q(qr(matrix(rnorm(n * n), n, n)))
  A <- U0 %*% (d_true * t(V0))

  fit <- svd_partial(A, rank = 5, tol = 1e-9)
  expect_true(fit$certificate$passed)
  expect_equal(fit$d, d_true[1:5], tolerance = 1e-8)
})
