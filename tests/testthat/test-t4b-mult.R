# C50: target completeness. Residual certificates prove each returned pair is
# an eigenpair; they cannot prove the returned SET is the requested one. A
# single-vector Krylov method can miss a copy of an exactly repeated
# eigenvalue (9, 9, 7, 7, 7 -> 9, 9, 7, 7, 5) and still certify every pair.
# The deflated-complement probe must catch and repair that.
#
# These are the probe's regression tests, so the file pins the completeness
# mode to "probe"; the default "auto" mode replaces the probe by the inertia
# certificate when the matrix is explicit (tests in test-t5-inertia.R).
withr::local_options(list(eigencore.target_completeness = "probe"))

t4b_diag_sparse <- function(d, seed) {
  n <- length(d)
  set.seed(seed)
  P <- sample(n)
  Matrix::sparseMatrix(i = P, j = P, x = d, dims = c(n, n))
}

# Sparse but not diagonal: random 2x2 rotations on disjoint index pairs.
t4b_rotated_sparse <- function(d, seed) {
  n <- length(d)
  set.seed(seed)
  P <- sample(n)
  ii <- jj <- integer()
  xx <- numeric()
  for (t in seq_len(n %/% 2L)) {
    a <- P[2L * t - 1L]
    b <- P[2L * t]
    th <- stats::runif(1L, 0, 2 * pi)
    ii <- c(ii, a, a, b, b)
    jj <- c(jj, a, b, a, b)
    xx <- c(xx, cos(th), -sin(th), sin(th), cos(th))
  }
  if (n %% 2L) {
    ii <- c(ii, P[n]); jj <- c(jj, P[n]); xx <- c(xx, 1)
  }
  R <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  S <- R %*% Matrix::Diagonal(n, d) %*% Matrix::t(R)
  S <- (S + Matrix::t(S)) / 2
  methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
}

t4b_dense <- function(d, seed) {
  set.seed(seed)
  n <- length(d)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% (d * t(Q))
  (A + t(A)) / 2
}

t4b_callback <- function(M) {
  linear_operator(
    dim(M),
    function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * as.matrix(M %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian()
  )
}

t4b_spectrum <- function(n, top = c(9, 9, 7, 7, 7)) {
  c(top, seq(5, 0.1, length.out = n - length(top)))
}

t4b_same_set <- function(fit, want, tol = 1e-6) {
  got <- sort(values(fit))
  length(got) == length(want) && max(abs(got - sort(want))) <= tol
}

with_completeness <- function(mode, code) {
  old <- options(eigencore.target_completeness = mode)
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("scalar Lanczos misses a repeated copy without the probe (the C50 hazard)", {
  d <- t4b_spectrum(60L)
  wrong_certified <- 0L
  for (s in 1:6) {
    S <- t4b_diag_sparse(d, 1000L + s)
    fit <- eig_partial(S, 5L, largest(),
                       method = lanczos(completeness = "none"), seed = s)
    if (isTRUE(certificate(fit)$passed) && !t4b_same_set(fit, c(9, 9, 7, 7, 7))) {
      wrong_certified <- wrong_certified + 1L
    }
    expect_identical(certificate(fit)$target_completeness, "not_checked")
  }
  # Documents the hazard the probe exists for; if the base solver ever stops
  # missing copies this expectation can be relaxed.
  expect_gt(wrong_certified, 0L)
})

test_that("the probe repairs missed copies: sparse diagonal, LA and SA", {
  d <- t4b_spectrum(60L)
  for (s in 1:6) {
    S <- t4b_diag_sparse(d, 1000L + s)
    fit <- eig_partial(S, 5L, largest(), method = lanczos(), seed = s)
    cert <- certificate(fit)
    expect_true(cert$passed)
    expect_true(t4b_same_set(fit, c(9, 9, 7, 7, 7)), info = paste("seed", s))
    expect_true(cert$target_completeness %in% c("probed", "repaired"))
    expect_true(cert$target_passed)
    fit <- eig_partial(-S, 5L, smallest(), method = lanczos(), seed = s)
    expect_true(certificate(fit)$passed)
    expect_true(t4b_same_set(fit, -c(9, 9, 7, 7, 7)), info = paste("SA seed", s))
  }
  # At least one of these seeds actually exercises the repair.
  S <- t4b_diag_sparse(d, 1001L)
  fit <- eig_partial(S, 5L, largest(), method = lanczos(), seed = 1L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "repaired")
  expect_true(cert$completeness$intruder_found)
  expect_gte(cert$completeness$rounds, 1L)
  expect_true(any(grepl("repaired", fit$warnings)))
  expect_true(cert$residual_passed)
})

test_that("rotated sparse and dense operators return the full repeated set", {
  d <- t4b_spectrum(60L)
  for (s in 1:5) {
    S <- t4b_rotated_sparse(d, 2000L + s)
    fit <- eig_partial(S, 5L, largest(), method = lanczos(max_subspace = 12L),
                       seed = s)
    expect_true(certificate(fit)$passed)
    expect_true(t4b_same_set(fit, c(9, 9, 7, 7, 7)), info = paste("sparse", s))
    A <- t4b_dense(d, 3000L + s)
    fit <- eig_partial(A, 5L, smallest(),
                       method = lanczos(max_subspace = 12L), seed = s)
    # smallest of d: 0.1 and the next four grid values (distinct) -- the probe
    # must not disturb a correct non-degenerate answer.
    expect_true(certificate(fit)$passed)
    expect_true(t4b_same_set(fit, sort(d)[1:5]))
    fit <- eig_partial(-A, 5L, smallest(),
                       method = lanczos(max_subspace = 12L), seed = s)
    expect_true(certificate(fit)$passed)
    expect_true(t4b_same_set(fit, -c(9, 9, 7, 7, 7)), info = paste("dense", s))
  }
})

test_that("matrix-free block Lanczos (block < multiplicity) returns the full set", {
  d <- t4b_spectrum(60L, top = c(9, 9, 9, 9, 7))
  for (s in 1:4) {
    op <- t4b_callback(t4b_diag_sparse(d, 4000L + s))
    fit <- eig_partial(op, 5L, largest(), method = lanczos(block = 2L), seed = s)
    cert <- certificate(fit)
    expect_true(cert$passed)
    expect_true(t4b_same_set(fit, c(9, 9, 9, 9, 7)))
    expect_true(cert$target_completeness %in% c("probed", "repaired"))
    expect_gt(cert$completeness$operator_columns, 0L)
  }
})

test_that("largest_magnitude with +/- repeated values returns the full set", {
  d <- c(9, -9, 9, -9, 7, seq(5, -5, length.out = 55L))
  for (s in 1:4) {
    S <- t4b_diag_sparse(d, 5000L + s)
    fit <- eig_partial(S, 4L, largest_magnitude(), method = lanczos(), seed = s)
    expect_true(certificate(fit)$passed)
    expect_true(t4b_same_set(fit, c(-9, -9, 9, 9)), info = paste("seed", s))
  }
})

test_that("near-multiplicity clusters are returned completely", {
  d <- t4b_spectrum(60L, top = c(9, 9 + 1e-9, 7, 7 + 1e-10, 7 - 1e-10))
  for (s in 1:5) {
    S <- t4b_diag_sparse(d, 6000L + s)
    fit <- eig_partial(S, 5L, largest(), method = lanczos(), seed = s)
    expect_true(certificate(fit)$passed)
    expect_true(t4b_same_set(fit, d[1:5], tol = 1e-6), info = paste("seed", s))
  }
})

test_that("generalized SPD (diagonal B) Lanczos is probed in the transformed space", {
  n <- 60L
  bdiag <- seq(1, 3, length.out = n)
  d <- t4b_spectrum(n)
  for (s in 1:4) {
    C <- t4b_diag_sparse(d, 7000L + s)
    Bh <- Matrix::Diagonal(n, sqrt(bdiag))
    A <- as.matrix(Bh %*% C %*% Bh)
    A <- (A + t(A)) / 2
    B <- Matrix::Diagonal(n, bdiag)
    fit <- eig_partial(A, 5L, largest(), B = B, method = lanczos(max_subspace = 12L),
                       seed = s)
    cert <- certificate(fit)
    expect_identical(fit$method, "native transformed generalized SPD B-orthogonal Lanczos")
    expect_true(cert$target_completeness %in% c("probed", "repaired"))
    expect_true(cert$passed)
    expect_true(t4b_same_set(fit, c(9, 9, 7, 7, 7)), info = paste("seed", s))
    expect_identical(cert$completeness$space, "transformed_standard_problem")
    if (s == 2L) {
      # Without the probe this seed returns 9, 9, 7, 7, 5 and certifies.
      expect_identical(cert$target_completeness, "repaired")
      bare <- eig_partial(A, 5L, largest(), B = B,
                          method = lanczos(max_subspace = 12L, completeness = "none"),
                          seed = s)
      expect_true(certificate(bare)$passed)
      expect_false(t4b_same_set(bare, c(9, 9, 7, 7, 7)))
    }
  }
})

test_that("the probe never consumes or perturbs the global RNG stream", {
  d <- t4b_spectrum(60L)
  S <- t4b_diag_sparse(d, 1001L)
  rng_after <- function(mode) {
    with_completeness(mode, {
      set.seed(99)
      fit <- eig_partial(S, 5L, largest(), method = lanczos())
      list(seed = .Random.seed, fit = fit)
    })
  }
  a <- rng_after("none")
  b <- rng_after("probe")
  # Whether the solve misses a copy first (and needs repair) depends on BLAS
  # rounding; this test is about the RNG stream, so accept either outcome.
  expect_true(certificate(b$fit)$target_completeness %in% c("probed", "repaired"))
  expect_identical(a$seed, b$seed)

  # Direct check with a non-default RNG kind: kind and state survive.
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(7)
  before <- .Random.seed
  # eig_partial itself draws its start vector from the stream; compare against
  # a probe-free solve from the same state.
  fit <- eig_partial(S, 5L, largest(), method = lanczos())
  expect_identical(certificate(fit)$target_completeness, "repaired")
  after_probe <- .Random.seed
  assign(".Random.seed", before, envir = .GlobalEnv)
  fit0 <- eig_partial(S, 5L, largest(), method = lanczos(completeness = "none"))
  expect_identical(.Random.seed, after_probe)
  expect_identical(RNGkind()[[1L]], "L'Ecuyer-CMRG")
  st <- .Random.seed
  invisible(eigencore:::completeness_probe_start(10L, 2L))
  expect_identical(.Random.seed, st)
})

test_that("vectors = FALSE still probes and repairs", {
  d <- t4b_spectrum(60L)
  S <- t4b_diag_sparse(d, 1001L)
  fit <- eig_partial(S, 5L, largest(), method = lanczos(), seed = 1L,
                     vectors = FALSE)
  expect_null(vectors(fit))
  expect_identical(certificate(fit)$target_completeness, "repaired")
  expect_true(t4b_same_set(fit, c(9, 9, 7, 7, 7)))
})

test_that("an unresolved intruder fails the certificate", {
  d <- t4b_spectrum(60L)
  S <- t4b_diag_sparse(d, 1001L)
  old <- options(eigencore.completeness_max_rounds = 0L)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(S, 5L, largest(), method = lanczos(), seed = 1L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "failed")
  expect_false(cert$target_passed)
  expect_true(cert$residual_passed)
  expect_false(cert$passed)
  expect_true(cert$completeness$intruder_found)
  expect_gt(cert$completeness$most_preferred_complement, cert$completeness$edge)
  expect_true(any(grepl("completeness", fit$warnings)))
})

test_that("completeness labels: exact, not_checked, option and descriptor", {
  A <- t4b_dense(t4b_spectrum(40L), 11L)
  expect_identical(certificate(eig_partial(A, 3L))$target_completeness, "exact")
  fit <- eig_partial(A, 3L, method = lanczos(), seed = 1L)
  expect_identical(certificate(fit)$target_completeness, "probed")
  expect_true(certificate(fit)$target_passed)
  rec <- certificate(fit)$completeness
  expect_identical(rec$method, "deflated_complement_probe")
  expect_gte(rec$steps, 1L)
  expect_false(rec$intruder_found)

  fit <- eig_partial(A, 3L, method = lanczos(completeness = "none"), seed = 1L)
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  expect_true(is.na(certificate(fit)$target_passed))
  fit <- with_completeness("none", eig_partial(A, 3L, method = lanczos(), seed = 1L))
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  # The descriptor wins over the option.
  fit <- with_completeness("none", eig_partial(A, 3L, method = lanczos(completeness = "probe"),
                                               seed = 1L))
  expect_identical(certificate(fit)$target_completeness, "probed")
  fit <- eig_partial(A, 3L, method = lanczos(), seed = 1L, certify = FALSE)
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  fit <- eig_partial(A, 4L, target = both_ends(2L, 2L), method = lanczos(), seed = 1L)
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  expect_error(lanczos(completeness = "always"), "completeness")
})

test_that("probe work is counted in the typed work record", {
  d <- t4b_spectrum(200L, top = c(9, 8, 7, 6, 5.5))
  S <- t4b_rotated_sparse(d, 12L)
  none <- eig_partial(S, 4L, method = lanczos(completeness = "none"), seed = 3L)
  probe <- eig_partial(S, 4L, method = lanczos(), seed = 3L)
  rec <- certificate(probe)$completeness
  expect_identical(certificate(probe)$target_completeness, "probed")
  expect_identical(values(probe), values(none))
  expect_identical(work(probe)$operator_columns - work(none)$operator_columns,
                   as.integer(rec$operator_columns))
  # Matrix-free: callback applies are observed directly.
  calls <- 0L
  op <- linear_operator(dim(S), function(X, alpha = 1, beta = 0, Y = NULL) {
    calls <<- calls + ncol(as.matrix(X))
    Z <- alpha * as.matrix(S %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  }, structure = hermitian())
  fit <- eig_partial(op, 4L, method = lanczos(block = 2L), seed = 3L)
  expect_identical(certificate(fit)$target_completeness, "probed")
  w <- work(fit)
  expect_identical(as.integer(w$operator_columns + w$certification_operator_columns),
                   as.integer(calls))
})
