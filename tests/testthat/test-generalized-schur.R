schur_sort <- function(z) z[order(round(Re(z), 7), round(Im(z), 7))]

expect_real_qz_reconstruction <- function(qz, A, B, tolerance = 1e-10) {
  expect_equal(qz$Q %*% qz$S %*% t(qz$Z), A, tolerance = tolerance)
  expect_equal(qz$Q %*% qz$T %*% t(qz$Z), B, tolerance = tolerance)
  expect_equal(crossprod(qz$Q), diag(nrow(A)), tolerance = tolerance)
  expect_equal(crossprod(qz$Z), diag(nrow(A)), tolerance = tolerance)
}

expect_complex_qz_reconstruction <- function(qz, A, B, tolerance = 1e-10) {
  expect_equal(qz$Q %*% qz$S %*% Conj(t(qz$Z)), A, tolerance = tolerance)
  expect_equal(qz$Q %*% qz$T %*% Conj(t(qz$Z)), B, tolerance = tolerance)
  I <- diag(as.complex(rep(1, nrow(A))))
  expect_equal(Conj(t(qz$Q)) %*% qz$Q, I, tolerance = tolerance)
  expect_equal(Conj(t(qz$Z)) %*% qz$Z, I, tolerance = tolerance)
}

test_that("generalized_schur computes real dense QZ factors and values", {
  A <- matrix(c(0, -1, 1, 0), 2, 2)
  B <- diag(2)

  qz <- generalized_schur(A, B)
  eig <- eig_full(A, B = B, structure = general(), tol = 1e-10)

  expect_s3_class(qz, "eigencore_generalized_schur_result")
  expect_identical(qz$method, eigencore:::native_dense_generalized_schur_label())
  expect_identical(qz$plan$method, qz$method)
  expect_identical(qz$plan$controls$qz, TRUE)
  expect_equal(schur_sort(values(qz)), schur_sort(values(eig)), tolerance = 1e-10)
  expect_equal(qz$classification, rep("finite", 2L))
  expect_real_qz_reconstruction(qz, A, B)
})

test_that("generalized_schur computes complex dense QZ factors and values", {
  A <- matrix(c(1 + 1i, 2, -1i, 3 + 0.5i), 2, 2)
  B <- matrix(c(2, 0.25i, -0.5i, 1.5), 2, 2)

  qz <- generalized_schur(A, B)
  eig <- eig_full(A, B = B, structure = general(), tol = 1e-9)

  expect_equal(schur_sort(values(qz)), schur_sort(values(eig)), tolerance = 1e-8)
  expect_equal(qz$classification, rep("finite", 2L))
  expect_true(is.complex(qz$S))
  expect_true(is.complex(qz$Q))
  expect_complex_qz_reconstruction(qz, A, B, tolerance = 1e-9)
})

test_that("generalized_schur exposes beta-zero classifications", {
  A <- diag(c(2, 3, 0))
  B <- diag(c(1, 0, 0))

  qz <- generalized_schur(A, B, sort = "none")

  expect_equal(Re(qz$alpha), c(2, 3, 0), tolerance = 1e-12)
  expect_equal(Re(qz$beta), c(1, 0, 0), tolerance = 1e-12)
  expect_equal(qz$classification, c("finite", "infinite", "undefined"))
  expect_equal(values(qz)[[1L]], 2 + 0i, tolerance = 1e-12)
  expect_true(is.infinite(values(qz)[[2L]]))
  expect_true(is.na(values(qz)[[3L]]))
})

test_that("generalized_schur supports finite/infinite sort classes", {
  A <- diag(c(2, 3, 0))
  B <- diag(c(1, 0, 0))

  finite <- generalized_schur(A, B, sort = "finite")
  infinite <- generalized_schur(A, B, sort = "infinite")

  expect_equal(finite$sdim, 1L)
  expect_equal(finite$classification[[1L]], "finite")
  expect_equal(infinite$sdim, 1L)
  expect_equal(infinite$classification[[1L]], "infinite")
})

test_that("generalized_schur validates inputs and rejects unsupported sort predicates", {
  expect_error(generalized_schur(Matrix::Diagonal(2), diag(2)), "base dense matrix")
  expect_error(generalized_schur(diag(2), diag(3)), "same dimension")
  expect_error(generalized_schur(diag(2), diag(2), vectors = NA), "vectors")
  expect_error(
    generalized_schur(diag(2), diag(2), sort = function(alpha, beta) TRUE),
    "unsupported generalized_schur\\(\\) sort"
  )
  expect_error(generalized_schur(diag(2), diag(2), sort = "undefined"), "arg")
})
