# Tranche 4a miscellany: C48 (fused finite/symmetry screen), P10 (operator
# algebra materialisation policy), C49 (PSD factor hash_format guard), C51
# (matrix-free smallest/interior SVD without Frobenius metadata).

test_that("C48: fused finite/symmetry pass matches the symmetry kernel", {
  fused <- eigencore:::dense_finite_symmetric
  sym <- function(x) {
    isTRUE(.Call("eigencore_dense_is_symmetric", x, sqrt(.Machine$double.eps),
                 PACKAGE = "eigencore"))
  }
  set.seed(4801)
  for (n in c(1L, 2L, 63L, 64L, 65L, 130L)) {
    A <- crossprod(matrix(rnorm(n * n), n))
    expect_identical(fused(A), c(TRUE, TRUE))
    if (n > 1L) {
      B <- A
      B[n, 1L] <- B[n, 1L] + 1
      expect_identical(fused(B), c(TRUE, sym(B)))
      expect_false(fused(B)[[2L]])
      # Relative tolerance with no floor at 1: a tiny-scale asymmetry is
      # still an asymmetry.
      S <- A * 1e-200
      S[n, 1L] <- S[n, 1L] * (1 + 1e-6)
      expect_identical(fused(S)[[2L]], sym(S))
      expect_false(fused(S)[[2L]])
    }
    for (bad in c(NA, NaN, Inf, -Inf)) {
      C <- A
      C[sample.int(n * n, 1L)] <- bad
      expect_identical(fused(C), c(FALSE, FALSE))
    }
  }
  expect_identical(fused(matrix(c(1, 2, 3, 4, 5, 6), 2)), c(TRUE, FALSE))
  expect_identical(fused(matrix(c(1, NaN, 3, 4, 5, 6), 2)), c(FALSE, FALSE))
  expect_identical(fused(matrix(0, 3, 3)), c(TRUE, TRUE))
})

test_that("C48: as_operator keeps its non-finite errors and structure detection", {
  msg <- "Matrix input contains NA, NaN, or Inf entries."
  expect_error(as_operator(matrix(c(1, NA, 2, 3), 2)), msg, fixed = TRUE)
  expect_error(as_operator(matrix(c(1L, NA, 2L, 3L), 2)), msg, fixed = TRUE)
  expect_error(as_operator(matrix(c(1, Inf, 2, 3), 2)), msg, fixed = TRUE)
  expect_error(
    as_operator(Matrix::Matrix(matrix(c(1, NA, 2, 3), 2), sparse = FALSE)),
    msg, fixed = TRUE
  )
  expect_error(as_operator(matrix(complex(real = c(1, NA, 2, 3)), 2)), msg, fixed = TRUE)

  set.seed(4802)
  S <- crossprod(matrix(rnorm(25), 5))
  expect_identical(as_operator(S)$structure$kind, "hermitian")
  expect_identical(as_operator(S + 1e-3 * upper.tri(S))$structure$kind, "general")
  expect_identical(as_operator(matrix(1:4, 2))$structure$kind, "general")
  expect_identical(as_operator(diag(1:3))$structure$kind, "hermitian")
  dense_matrix <- Matrix::Matrix(S, sparse = FALSE)
  op <- as_operator(dense_matrix)
  expect_identical(op$structure$kind, "hermitian")
  expect_true(op$metadata$native)
  expect_equal(op$apply(diag(5)), S)
})

test_that("P10: small dense and sparse products materialise once", {
  set.seed(1001)
  A <- matrix(rnorm(15), 5)
  B <- matrix(rnorm(12), 3)
  op <- compose(A, B)
  expect_identical(op$metadata$fused, "compose")
  expect_equal(op$metadata$source, A %*% B)
  expect_true(op$metadata$native)
  cp <- crossprod_operator(A)
  expect_true(isTRUE(cp$metadata$materialized_crossprod))
  expect_equal(cp$metadata$source, crossprod(A))
})

test_that("P10: a large low-rank dense composition stays lazy", {
  set.seed(1002)
  n <- 400L
  k <- 3L
  L <- matrix(rnorm(n * k), n, k)
  R <- matrix(rnorm(k * n), k, n)
  op <- compose(L, R)
  # n x n product (160000 entries) is larger than the factors (2400): lazy.
  expect_null(op$metadata$fused)
  expect_null(eigencore:::source_or_null(op))
  expect_false(isTRUE(op$metadata$native))
  X <- matrix(rnorm(n * 2), n, 2)
  expect_equal(op$apply(X), L %*% (R %*% X))
  expect_equal(op$apply_adjoint(X), t(R) %*% (t(L) %*% X))
  expect_true(check_adjoint(op, seed = 7)$passed)

  # crossprod of a wide dense matrix (ncol > nrow) stays lazy as well.
  W <- matrix(rnorm(3 * n), 3, n)
  cpw <- crossprod_operator(W)
  expect_false(isTRUE(cpw$metadata$materialized_crossprod))
  expect_identical(cpw$structure$kind, "hermitian")
  expect_equal(cpw$apply(X), crossprod(W) %*% X)

  # a tall dense crossprod is small and materialises.
  cpt <- crossprod_operator(t(W))
  expect_true(isTRUE(cpt$metadata$materialized_crossprod))
})

test_that("P10: a sparse product that would densify stays lazy", {
  set.seed(1003)
  n <- 600L
  # One dense column in A and one dense row in B: A %*% B is fully dense.
  A <- Matrix::sparseMatrix(i = seq_len(n), j = rep(1L, n), x = rnorm(n),
                            dims = c(n, n))
  B <- Matrix::sparseMatrix(i = rep(1L, n), j = seq_len(n), x = rnorm(n),
                            dims = c(n, n))
  op <- compose(A, B)
  expect_null(op$metadata$fused)
  expect_null(op$metadata$matrix)
  X <- matrix(rnorm(n * 2), n, 2)
  expect_equal(op$apply(X), as.matrix(A %*% (B %*% X)))

  # crossprod of a matrix with one dense row densifies A^T A: lazy.
  D <- Matrix::sparseMatrix(i = rep(1L, n), j = seq_len(n), x = rnorm(n),
                            dims = c(n, n))
  cpd <- crossprod_operator(D)
  expect_false(isTRUE(cpd$metadata$materialized_crossprod))
  expect_equal(cpd$apply(X), as.matrix(Matrix::crossprod(D) %*% X))

  # A banded product stays sparse and materialises.
  T1 <- Matrix::bandSparse(n, k = c(-1, 0, 1),
                           diagonals = list(rep(-1, n - 1), rep(2, n), rep(-1, n - 1)))
  T1 <- methods::as(T1, "generalMatrix")
  banded <- compose(T1, T1)
  expect_identical(banded$metadata$fused, "compose")
  expect_s4_class(banded$metadata$matrix, "dgCMatrix")
  expect_equal(banded$apply(X), as.matrix(T1 %*% (T1 %*% X)))
})

test_that("P10: dense sums still fuse and are evaluated correctly", {
  set.seed(1004)
  A <- matrix(rnorm(20), 4)
  B <- matrix(rnorm(20), 4)
  op <- eigencore:::operator_sum(as_operator(A), as_operator(B))
  expect_identical(op$metadata$fused, "sum")
  expect_equal(op$metadata$source, A + B)
})

test_that("C49: PSD factors record and enforce the identity hash format", {
  factor <- psd_factor(c(9, 4, 0))
  expect_identical(
    factor$serialization$hash_format,
    eigencore:::identity_hash_format()
  )
  gram <- psd_gram_factor(matrix(c(1, 2, 3, 4, 5, 6), 3))
  expect_identical(
    gram$serialization$hash_format,
    eigencore:::identity_hash_format()
  )
  expect_identical(psd_rank(factor), 2L)

  old <- factor
  old$serialization$hash_format <- NULL
  err <- tryCatch(psd_rank(old), error = function(e) e)
  expect_s3_class(err, "eigencore_psd_corrupt_state")
  expect_identical(err$code, "identity_format_changed")
  expect_identical(err$field, "serialization$hash_format")
  expect_match(conditionMessage(err), "identity format changed", fixed = TRUE)
  expect_match(conditionMessage(err), "Re-factor", fixed = TRUE)

  stale <- factor
  stale$serialization$hash_format <- "eigencore-identity-hash-v1"
  err <- tryCatch(psd_apply(stale, c(1, 2, 3)), error = function(e) e)
  expect_identical(err$code, "identity_format_changed")
  expect_identical(err$actual, "eigencore-identity-hash-v1")
})

test_that("C49: the legacy FNV raw hash entry point is gone", {
  expect_false(is.loaded("eigencore_stable_raw_hash", PACKAGE = "eigencore"))
})

c51_callback_operator <- function(M) {
  linear_operator(
    dim = dim(M),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (M %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * crossprod(M, X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    }
  )
}

test_that("C51: matrix-free smallest SVD without metadata routes natively and certifies", {
  set.seed(5101)
  m <- 60L
  n <- 40L
  d <- seq(10, 0.5, length.out = n)
  U <- qr.Q(qr(matrix(rnorm(m * n), m)))
  V <- qr.Q(qr(matrix(rnorm(n * n), n)))
  M <- U %*% diag(d) %*% t(V)
  op <- c51_callback_operator(M)
  expect_null(op$metadata$frobenius_norm)
  expect_null(eigencore:::source_or_null(op))

  plan <- plan_solver(svd_problem(op, target = smallest()), rank = 3L)
  expect_identical(
    plan$method,
    eigencore:::native_matrix_free_smallest_golub_kahan_label()
  )
  expect_false(plan$controls$requires_nonestimated_norm_scale)

  fit <- svd_partial(op, rank = 3L, target = smallest(), tol = 1e-10,
                     seed = 51, allow_dense_fallback = "never")
  expect_identical(fit$method, eigencore:::native_matrix_free_smallest_golub_kahan_label())
  expect_true(fit$certificate$passed)
  expect_false(fit$certificate$scale_is_estimate)
  expect_identical(fit$certificate$norm_bound_type, "two_norm_lower_bound")
  expect_true(fit$restart$matrix_free)
  expect_equal(sort(fit$d), sort(tail(d, 3L)), tolerance = 1e-9)
  expect_certificate_clean(fit, tol = 1e-10)
})

test_that("C51: matrix-free interior SVD without metadata routes natively and certifies", {
  set.seed(5102)
  m <- 30L
  n <- 20L
  d <- seq(6, 0.3, length.out = n)
  U <- qr.Q(qr(matrix(rnorm(m * n), m)))
  V <- qr.Q(qr(matrix(rnorm(n * n), n)))
  M <- U %*% diag(d) %*% t(V)
  op <- c51_callback_operator(M)
  fit <- svd_partial(op, rank = 2L, target = nearest(2), tol = 1e-10,
                     seed = 52, allow_dense_fallback = "never")
  expect_identical(fit$method, eigencore:::native_matrix_free_interior_golub_kahan_label())
  # Matrix-free nearest() SVD: no inertia count and no sound complement probe
  # for an interior target, so only the residuals are certified.
  expect_true(fit$certificate$residual_passed)
  expect_false(fit$certificate$passed)
  oracle <- d[order(abs(d - 2))][1:2]
  expect_equal(sort(fit$d), sort(oracle), tolerance = 1e-9)
})
