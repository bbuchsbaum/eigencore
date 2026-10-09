# Determinism invariants (docs/test-assurance.md, invariant d):
#   * same inputs + same seed => identical results, on every sweep family;
#   * results with eigencore.threads = 1 and 4 agree to ~1e-12 (relative);
#   * seed = never moves the caller's RNG stream.

oracle_determinism_ids <- function() {
  # a few cases of every family/method from the sweep's id stream
  ids <- seq_len(if (oracle_level() == "cran") 40L else 160L)
  ids[vapply(ids, function(i) {
    case <- oracle_case(i)
    !isTRUE(case$shim)
  }, NA)]
}

test_that("same inputs and seed give identical results across the sweep families", {
  old <- options(eigencore.threads = 1L)
  on.exit(options(old), add = TRUE)
  checked <- 0L
  for (id in oracle_determinism_ids()) {
    case <- oracle_case(id)
    prob <- oracle_build(case)
    truth <- oracle_truth(case, prob)
    sigma <- oracle_sigma(case, truth$values)
    a <- tryCatch(suppressWarnings(oracle_invoke(case, prob, sigma)), error = function(e) e)
    if (inherits(a, "error")) next
    prob2 <- oracle_build(case)
    b <- suppressWarnings(oracle_invoke(case, prob2, sigma))
    label <- sprintf("case %d [%s]", id, oracle_describe(case))
    va <- oracle_or(a$d, a$values)
    vb <- oracle_or(b$d, b$values)
    expect_identical(va, vb, label = label)
    expect_identical(oracle_or(a$vectors, a$u), oracle_or(b$vectors, b$u), label = label)
    expect_identical(a$certificate$passed, b$certificate$passed, label = label)
    checked <- checked + 1L
  }
  expect_gt(checked, 10L)
})

test_that("seed = leaves the global RNG stream untouched", {
  set.seed(2024)
  A <- crossprod(matrix(rnorm(60 * 40), 60, 40))
  X <- matrix(rnorm(80 * 30), 80, 30)
  S <- Matrix::rsparsematrix(200, 200, 0.03, symmetric = TRUE)
  calls <- list(
    function() eig_partial(A, 3, seed = 1),
    function() eig_partial(A, 3, method = lanczos(block = 2), seed = 1),
    function() eig_partial(A, 3, method = lobpcg(), seed = 1),
    function() eig_partial(S, 3, target = smallest(), seed = 1),
    function() eig_partial(matrix(rnorm(400), 20), 2,
                           target = largest_magnitude(), seed = 1),
    function() svd_partial(X, 3, seed = 1),
    function() svd_partial(X, 3, method = randomized(), seed = 1),
    function() svd_partial(X, 3, method = golub_kahan(), seed = 1),
    function() svd_partial(Matrix::Matrix(X, sparse = TRUE), 3, seed = 1)
  )
  for (i in seq_along(calls)) {
    set.seed(77)
    before <- .Random.seed
    suppressWarnings(calls[[i]]())
    expect_identical(.Random.seed, before, label = sprintf("call %d", i))
  }
})

test_that("results agree across eigencore.threads = 1 and 4", {
  skip_on_cran()
  old <- getOption("eigencore.threads")
  on.exit(options(eigencore.threads = old), add = TRUE)
  set.seed(31)
  n <- 2500L
  S <- Matrix::rsparsematrix(n, n, 6 / n, symmetric = TRUE) +
    Matrix::Diagonal(n, x = seq_len(n) / n)
  S <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
  N <- Matrix::rsparsematrix(n, n, 6 / n) + Matrix::Diagonal(n, x = seq_len(n) / n)
  R <- Matrix::rsparsematrix(5000L, 500L, 0.01)
  runs <- list(
    sym = function() eig_partial(S, 6, seed = 5)$values,
    sym_small = function() eig_partial(S, 4, target = smallest(), seed = 5)$values,
    block = function() eig_partial(S, 6, method = lanczos(block = 3), seed = 5)$values,
    nonsym = function() eig_partial(N, 4, target = largest_magnitude(), seed = 5,
                                    left_vectors = "none")$values,
    svd = function() svd_partial(R, 6, seed = 5)$d,
    svd_center = function() svd_partial(center(R), 5, seed = 5)$d,
    crossprod = function() eig_partial(crossprod_operator(scale_cols(R, seq_len(500) / 500)),
                                       4, seed = 5)$values
  )
  for (nm in names(runs)) {
    options(eigencore.threads = 1L)
    one <- runs[[nm]]()
    options(eigencore.threads = 4L)
    four <- runs[[nm]]()
    expect_equal(length(one), length(four), label = nm)
    rel <- max(Mod(sort(Mod(one)) - sort(Mod(four)))) / max(Mod(one))
    expect_lte(rel, 1e-12, label = sprintf("%s: threads 1 vs 4 relative difference", nm))
  }
})
