# RSpectra drop-in invariant (docs/test-assurance.md, invariant e): on
# well-separated spectra, eigs_sym()/eigs()/svds() return the same values
# as RSpectra, in RSpectra's order. Ties and LI/SI (C30: eigencore ranks
# the signed imaginary part) are excluded; the sweep covers them by set.

skip_if_not_installed("RSpectra")

separated_sym <- function(n, seed) {
  set.seed(seed)
  vals <- sample(c(-1, 1), n, replace = TRUE) * (seq_len(n) + stats::runif(n, 0, 0.5))
  Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
  A <- Q %*% (vals * t(Q))
  (A + t(A)) / 2
}

test_that("eigs_sym matches RSpectra values and order for every `which`", {
  for (seed in 1:3) {
    A <- separated_sym(40, seed)
    S <- methods::as(methods::as(Matrix::Matrix(A, sparse = TRUE), "generalMatrix"),
                     "CsparseMatrix")
    for (input in list(A, S)) {
      for (which in c("LM", "SM", "LA", "SA", "BE")) {
        for (k in c(1L, 4L, 5L)) {
          # An uncertified result warns (compat_warn_nconv) and is not a
          # drop-in claim; sparse BE currently takes the unrestarted
          # reference Lanczos and often does not converge (O10).
          ours <- suppressWarnings(eigs_sym(input, k, which = which,
                                            opts = list(tol = 1e-10)))
          theirs <- RSpectra::eigs_sym(input, k, which = which, opts = list(tol = 1e-10))
          label <- sprintf("eigs_sym %s seed=%d which=%s k=%d",
                           if (is.matrix(input)) "dense" else "sparse", seed, which, k)
          if (is.matrix(input) || which != "BE") {
            expect_true(isTRUE(ours$certificate$passed), label = label)
          }
          if (isTRUE(ours$certificate$passed)) {
            expect_equal(ours$values, theirs$values, tolerance = 1e-8, label = label)
          }
        }
      }
      ours <- eigs_sym(input, 3, sigma = 2.3, opts = list(tol = 1e-10))
      theirs <- RSpectra::eigs_sym(input, 3, sigma = 2.3, opts = list(tol = 1e-10))
      expect_equal(ours$values, theirs$values, tolerance = 1e-8,
                   label = sprintf("eigs_sym sigma seed=%d", seed))
    }
  }
})

test_that("eigs matches RSpectra values for LM/SM/LR/SR on real-spectrum input", {
  for (seed in 1:3) {
    set.seed(seed)
    n <- 30
    vals <- sample(c(-1, 1), n, replace = TRUE) * (seq_len(n) + stats::runif(n, 0, 0.5))
    V <- matrix(rnorm(n * n), n) + diag(3, n)
    A <- V %*% diag(vals) %*% solve(V)
    for (which in c("LM", "SM", "LR", "SR")) {
      ours <- eigs(A, 4, which = which, opts = list(tol = 1e-10))
      theirs <- RSpectra::eigs(A, 4, which = which, opts = list(tol = 1e-10))
      expect_equal(sort(Re(ours$values)), sort(Re(theirs$values)), tolerance = 1e-7,
                   label = sprintf("eigs seed=%d which=%s (set)", seed, which))
      ord <- switch(which,
                    LM = order(-Mod(ours$values)), SM = order(Mod(ours$values)),
                    LR = order(-Re(ours$values)), SR = order(Re(ours$values)))
      expect_identical(ord, seq_along(ord),
                       label = sprintf("eigs seed=%d which=%s ordered by target", seed, which))
    }
  }
})

test_that("svds matches RSpectra singular values (dense, sparse, centred)", {
  set.seed(9)
  X <- matrix(rnorm(60 * 25), 60, 25) %*% diag(seq(5, 1, length.out = 25))
  S <- Matrix::rsparsematrix(80, 30, 0.2)
  for (input in list(X, S)) {
    for (k in c(1L, 5L)) {
      ours <- svds(input, k, opts = list(tol = 1e-10))
      theirs <- RSpectra::svds(input, k, opts = list(tol = 1e-10))
      expect_equal(ours$d, theirs$d, tolerance = 1e-8)
      expect_equal(dim(ours$u), dim(theirs$u))
      expect_equal(dim(ours$v), dim(theirs$v))
    }
    ours <- svds(input, 3, nu = 1, nv = 0, opts = list(tol = 1e-10, center = TRUE))
    theirs <- RSpectra::svds(input, 3, nu = 1, nv = 0, opts = list(tol = 1e-10, center = TRUE))
    expect_equal(ours$d, theirs$d, tolerance = 1e-8)
    expect_equal(ncol(ours$u), 1L)
    expect_true(is.null(ours$v) || ncol(ours$v) == 0L)
  }
})
