geigen_migration_bank <- data.frame(
  geigen_function = c(
    "geigen",
    "geigen",
    "geigen",
    "geigen",
    "geigen",
    "geigen",
    "gqz",
    "gevalues",
    "gsvd"
  ),
  use_class = c(
    "dense SPD",
    "dense complex Hermitian SPD",
    "dense general pencil",
    "beta-zero general pencil",
    "sparse SPD",
    "sparse transformed SPD",
    "QZ",
    "homogeneous values",
    "GSVD"
  ),
  eigencore_surface = c(
    "eig_full(A, B = ...)",
    "eig_full(A, B = ...)",
    "eig_full(A, B = ..., structure = general())",
    "eig_full(A, B = ..., structure = general())",
    "eig_partial(A, B = ..., allow_dense_fallback = 'never')",
    "eig_partial(A, B = ..., method = shift_invert(), allow_dense_fallback = 'never')",
    "generalized_schur(A, B)",
    "values(x) and alpha_beta(x)",
    "generalized_svd(A, B) when available"
  ),
  status = c(
    "covered",
    "covered",
    "covered",
    "covered",
    "covered",
    "covered",
    "covered",
    "covered",
    "deferred"
  ),
  stringsAsFactors = FALSE
)

expect_same_complex_set <- function(actual, expected, tolerance = 1e-8) {
  actual <- as.complex(actual)
  expected <- as.complex(expected)
  expect_equal(length(actual), length(expected))
  remaining <- seq_along(actual)
  for (target in expected) {
    distances <- Mod(actual[remaining] - target)
    best <- which.min(distances)
    expect_lte(distances[[best]], tolerance)
    remaining <- remaining[-best]
  }
}

nearest_expected <- function(values, sigma, k) {
  values[order(abs(values - sigma))[seq_len(k)]]
}

test_that("migration bank covers the geigen manual and reverse-dep use classes", {
  expect_setequal(
    unique(geigen_migration_bank$geigen_function),
    c("geigen", "gqz", "gevalues", "gsvd")
  )
  expect_setequal(
    geigen_migration_bank$use_class,
    c(
      "dense SPD",
      "dense complex Hermitian SPD",
      "dense general pencil",
      "beta-zero general pencil",
      "sparse SPD",
      "sparse transformed SPD",
      "QZ",
      "homogeneous values",
      "GSVD"
    )
  )
  expect_true(any(
    geigen_migration_bank$geigen_function == "gsvd" &
      geigen_migration_bank$status == "deferred"
  ))
})

test_that("geigen symmetric dense cases map to certified eig_full results", {
  A <- diag(c(2, 8, 18))
  B <- diag(c(1, 2, 3))

  fit <- eig_full(A, B = B, tol = 1e-10)

  expect_equal(fit$method, eigencore:::native_dense_generalized_spd_full_label())
  expect_equal(sort(unname(values(fit))), c(2, 4, 6), tolerance = 1e-12)
  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(3),
               tolerance = 1e-10)
})

test_that("geigen complex Hermitian SPD cases map to certified eig_full results", {
  A <- matrix(c(3, 1 - 1i, 1 + 1i, 4), 2, 2)
  B <- diag(c(2, 3))
  expected <- eigen(solve(B, A), only.values = TRUE)$values

  fit <- eig_full(A, B = B, tol = 1e-10)

  expect_equal(fit$method, eigencore:::native_dense_generalized_spd_full_label())
  expect_same_complex_set(values(fit), expected, tolerance = 1e-10)
  expect_true(certificate(fit)$passed)
  gram <- Conj(t(vectors(fit))) %*% B %*% vectors(fit)
  expect_equal(gram, diag(as.complex(rep(1, 2))),
               tolerance = 1e-10)
})

test_that("geigen general pencils map to dense eig_full with alpha/beta evidence", {
  A <- matrix(c(1, 4, 2, 3), 2, 2)
  B <- matrix(c(2, 1, 0, -1), 2, 2)
  expected <- eigen(solve(B, A), only.values = TRUE)$values

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-10)
  coords <- alpha_beta(fit)

  expect_equal(fit$method, eigencore:::native_dense_generalized_pencil_full_label())
  expect_true(all(coords$classification == "finite"))
  expect_same_complex_set(values(fit), expected, tolerance = 1e-10)
  expect_equal(coords$values, values(fit))
  expect_true(certificate(fit)$passed)
})

test_that("beta-zero geigen edge cases have explicit finite/infinite labels", {
  A <- diag(c(2, 3, 0))
  B <- diag(c(1, 0, 0))

  fit <- eig_full(A, B = B, structure = general(), tol = 1e-10)
  coords <- alpha_beta(fit)

  expect_equal(coords$classification, c("finite", "infinite", "undefined"))
  expect_equal(coords$finite, c(TRUE, FALSE, FALSE))
  expect_equal(coords$infinite, c(FALSE, TRUE, FALSE))
  expect_equal(coords$undefined, c(FALSE, FALSE, TRUE))
  expect_equal(values(fit)[1L], 2 + 0i)
  expect_true(is.infinite(values(fit)[2L]))
  expect_true(is.na(values(fit)[3L]))
  expect_false(certificate(fit)$passed)
})

test_that("sparse SPD migration paths stay native and do not densify silently", {
  A <- Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36))
  B <- Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6))
  expected <- sort(c(1, 2, 3, 4, 5, 6))[seq_len(3)]

  fit <- eig_partial(
    A,
    B = B,
    k = 3L,
    target = smallest(),
    method = lanczos(max_subspace = 6L),
    seed = 3001,
    tol = 1e-9,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$method, eigencore:::native_generalized_lanczos_label())
  expect_equal(sort(unname(values(fit))), expected, tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_true(fit$restart$native)
  expect_equal(
    fit$restart$metric_solve,
    "diagonal scaling similarity transform for B"
  )
  expect_error(
    eig_full(A, B = B, allow_dense_fallback = "always"),
    "not silently densified"
  )
})

test_that("sparse transformed migration paths use explicit shift-invert labels", {
  n <- 16L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), rep(2.5, n), rep(-1, n - 1L))
  )
  B <- Matrix::Diagonal(x = seq(1, 2, length.out = n))
  oracle <- eigen(solve(as.matrix(B), as.matrix(A)), only.values = TRUE)$values
  sigma <- sort(Re(oracle))[5L] + 0.01
  expected <- nearest_expected(oracle, sigma, 3L)

  fit <- eig_partial(
    A,
    B = B,
    k = 3L,
    target = nearest(sigma),
    method = shift_invert(sigma = sigma),
    tol = 1e-8,
    allow_dense_fallback = "never"
  )

  expect_equal(
    fit$method,
    eigencore:::native_tridiagonal_generalized_shift_invert_label()
  )
  expect_same_complex_set(values(fit), expected, tolerance = 1e-6)
  expect_true(certificate(fit)$passed)
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
})

test_that("gqz and gevalues migrate to generalized_schur and accessors", {
  A <- matrix(c(2, 1, 0, 3), 2, 2)
  B <- matrix(c(1, 0, 0, 2), 2, 2)
  eig <- eig_full(A, B = B, structure = general(), tol = 1e-10)

  qz <- generalized_schur(A, B)
  coords <- alpha_beta(qz)

  expect_equal(qz$method, eigencore:::native_dense_generalized_schur_label())
  expect_equal(qz$classification, rep("finite", 2L))
  expect_equal(values(qz), coords$values)
  expect_equal(coords$alpha / coords$beta, values(qz), tolerance = 1e-12)
  expect_same_complex_set(values(qz), values(eig), tolerance = 1e-10)
  expect_error(generalized_schur(A, B, sort = "R"), "arg")
  expect_error(
    generalized_schur(A, B, sort = function(alpha, beta) TRUE),
    "unsupported generalized_schur\\(\\) sort"
  )
})

test_that("gsvd migration is explicitly deferred, not accidentally exported", {
  exports <- getNamespaceExports("eigencore")

  expect_false("gsvd" %in% exports)
  expect_false("generalized_svd" %in% exports)
  expect_equal(
    geigen_migration_bank$status[geigen_migration_bank$geigen_function == "gsvd"],
    "deferred"
  )
  expect_match(
    geigen_migration_bank$eigencore_surface[
      geigen_migration_bank$geigen_function == "gsvd"
    ],
    "generalized_svd"
  )
})

test_that("alpha_beta rejects result types without homogeneous coordinates", {
  fit <- eig_partial(diag(c(3, 2, 1)), k = 1L, target = largest())

  expect_error(
    alpha_beta(fit),
    "alpha/beta coordinates are not available"
  )
})
