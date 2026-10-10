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

test_that("implicitly centred sparse PCA runs on the implicit Gram kernel (P20)", {
  set.seed(107)
  A <- abs(Matrix::rsparsematrix(3000, 300, density = 0.02))
  dense <- as.matrix(A)
  centred <- sweep(dense, 2L, colMeans(dense), `-`)
  ref <- svd(centred, nu = 6, nv = 6)

  op <- center(A)
  plan <- plan_solver(svd_problem(op), rank = 6)
  expect_identical(plan$method, eigencore:::native_implicit_gram_svd_label())
  fit <- svd_partial(op, rank = 6, tol = 1e-8)
  expect_identical(fit$method, eigencore:::native_implicit_gram_svd_label())
  expect_true(fit$restart$centered_csc)
  expect_false(fit$restart$materialized_gram)
  expect_true(fit$certificate$passed)
  expect_equal(fit$d, ref$d[1:6], tolerance = 1e-8)
  expect_true(all(abs(colSums(fit$u * ref$u)) > 1 - 1e-6))
  expect_true(all(abs(colSums(fit$v * ref$v)) > 1 - 1e-6))
  # The residuals are those of the centred matrix, not of A.
  expect_lt(max(abs(centred %*% fit$v - sweep(fit$u, 2L, fit$d, `*`))),
            1e-6 * ref$d[[1L]])

  # Centred and column-scaled.
  w <- seq(0.5, 2, length.out = ncol(A))
  fit_s <- svd_partial(scale_cols(center(A), w), rank = 4, tol = 1e-8)
  expect_true(fit_s$restart$centered_csc)
  expect_true(fit_s$certificate$passed)
  expect_equal(fit_s$d, svd(sweep(centred, 2L, w, `*`), nu = 0, nv = 0)$d[1:4],
               tolerance = 1e-8)

  # Wide operator: the left (A A') Gram side.
  At <- methods::as(Matrix::t(A), "CsparseMatrix")
  wide <- t(dense)
  fit_w <- svd_partial(center(At), rank = 4, tol = 1e-8)
  expect_true(fit_w$restart$centered_csc)
  expect_identical(fit_w$restart$gram_side, "left")
  expect_true(fit_w$certificate$passed)
  expect_equal(fit_w$d,
               svd(sweep(wide, 2L, colMeans(wide), `-`), nu = 0, nv = 0)$d[1:4],
               tolerance = 1e-8)

  # Row (or double) centring is not a column rank-one correction: it keeps
  # the Golub-Kahan route.
  expect_null(eigencore:::implicit_gram_centered_csc_parts(center(A, rows = TRUE)))
})
