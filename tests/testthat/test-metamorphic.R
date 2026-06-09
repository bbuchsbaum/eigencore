expect_native_nonfallback <- function(fit) {
  expect_match(fit$method, "native")
  expect_false(grepl("oracle|reference|fallback", fit$method))
}

test_that("sparse Hermitian eigenvalues are invariant under permutation and positive scaling", {
  A <- as(Matrix::bandSparse(
    12L,
    k = c(-1L, 0L, 1L),
    diagonals = list(rep(-0.25, 11L), seq_len(12L), rep(-0.25, 11L))
  ), "dgCMatrix")

  perm <- c(12L, 1L, 6L, 3L, 9L, 2L, 11L, 4L, 8L, 5L, 10L, 7L)

  fit <- eig_partial(
    A,
    k = 4L,
    target = largest(),
    allow_dense_fallback = "never",
    tol = 1e-8,
    seed = 41L
  )
  permuted_fit <- eig_partial(
    A[perm, perm],
    k = 4L,
    target = largest(),
    allow_dense_fallback = "never",
    tol = 1e-8,
    seed = 41L
  )
  scaled_fit <- eig_partial(
    3 * A,
    k = 4L,
    target = largest(),
    allow_dense_fallback = "never",
    tol = 1e-8,
    seed = 41L
  )

  expect_native_nonfallback(fit)
  expect_native_nonfallback(permuted_fit)
  expect_native_nonfallback(scaled_fit)
  expect_certificate_clean(fit)
  expect_certificate_clean(permuted_fit)
  expect_certificate_clean(scaled_fit)
  expect_equal(sort(values(permuted_fit), decreasing = TRUE),
               sort(values(fit), decreasing = TRUE),
               tolerance = 1e-8)
  expect_equal(sort(values(scaled_fit), decreasing = TRUE),
               3 * sort(values(fit), decreasing = TRUE),
               tolerance = 1e-8)
})

test_that("sparse SVD singular values are invariant under row and column permutations", {
  M <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 3L, 7L),
    j = c(1L, 1L, 2L, 3L, 3L, 4L, 5L, 6L, 7L, 8L, 6L, 2L),
    x = c(4, -1, 3, 2, -0.5, 1.5, 2.5, -2, 1, 0.75, 0.4, -0.3),
    dims = c(10L, 8L)
  )
  row_perm <- c(10L, 1L, 7L, 2L, 9L, 4L, 3L, 8L, 5L, 6L)
  col_perm <- c(8L, 1L, 6L, 2L, 7L, 3L, 5L, 4L)

  fit <- svd_partial(
    M,
    rank = 4L,
    target = largest(),
    allow_dense_fallback = "never",
    tol = 1e-8,
    seed = 42L
  )
  permuted_fit <- svd_partial(
    M[row_perm, col_perm],
    rank = 4L,
    target = largest(),
    allow_dense_fallback = "never",
    tol = 1e-8,
    seed = 42L
  )
  oracle <- svd(as.matrix(M), nu = 0L, nv = 0L)$d[seq_len(4L)]

  expect_native_nonfallback(fit)
  expect_native_nonfallback(permuted_fit)
  expect_certificate_clean(fit)
  expect_certificate_clean(permuted_fit)
  expect_equal(fit$d, oracle, tolerance = 1e-8)
  expect_equal(permuted_fit$d, oracle, tolerance = 1e-8)
})
