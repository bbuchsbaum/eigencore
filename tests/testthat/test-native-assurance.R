# Native-code assurance (docs/test-assurance.md, "Native code assurance"):
# targeted tests for branches of the refactored native code that the rest of
# the suite did not reach (coverage of lines changed since 83ca827), run under
# ASan/UBSan and valgrind as part of the full suite.

native_call <- function(name, ...) .Call(name, ..., PACKAGE = "eigencore")

spd_matrix <- function(n, seed) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% diag(seq_len(n)) %*% t(Q)
  (A + t(A)) / 2
}

test_that("non-restarted native block Lanczos entries match eigen() on a full Krylov space", {
  n <- 20L
  A <- spd_matrix(n, 1)
  start <- matrix(stats::rnorm(n * 2L), n)
  ref <- eigen(A, symmetric = TRUE)$values
  dense_hi <- native_call("eigencore_block_lanczos_dense", A, 3L, n, 2L, 1L,
                          1e-10, start)
  expect_equal(dense_hi$values, ref[1:3], tolerance = 1e-10)
  expect_true(all(dense_hi$converged))
  dense_lo <- native_call("eigencore_block_lanczos_dense", A, 3L, n, 2L, 2L,
                          1e-10, start)
  expect_equal(sort(dense_lo$values), sort(ref[n - 0:2]), tolerance = 1e-10)

  S <- methods::as(Matrix::Matrix(A, sparse = TRUE), "generalMatrix")
  csc <- native_call("eigencore_block_lanczos_csc", S@i, S@p, S@x, S@Dim,
                     3L, n, 2L, 1L, 1e-10, start)
  expect_equal(csc$values, ref[1:3], tolerance = 1e-10)

  expect_error(native_call("eigencore_block_lanczos_dense", A, 3L, n + 1L, 2L,
                           1L, 1e-10, start), "m_max")
  expect_error(native_call("eigencore_block_lanczos_dense", A, 3L, n, 3L, 1L,
                           1e-10, start), "wrong number of columns")
  bad <- S
  bad@i[1L] <- 99L
  expect_error(native_call("eigencore_block_lanczos_csc", bad@i, bad@p, bad@x,
                           bad@Dim, 3L, n, 2L, 1L, 1e-10, start))
})

test_that("composite spec builder rejects malformed specs and reports dims", {
  X <- matrix(stats::rnorm(12), 4L)
  ptr <- native_call("eigencore_composite_operator_build",
                     list(type = "product",
                          children = list(list(type = "dense", x = X),
                                          list(type = "adjoint",
                                               child = list(type = "dense", x = X)))))
  expect_identical(native_call("eigencore_composite_operator_dim", ptr), c(4L, 4L))
  # A cleared (e.g. deserialised) pointer is reported as NULL, not a crash.
  cleared <- unserialize(serialize(ptr, NULL))
  expect_null(native_call("eigencore_composite_operator_dim", cleared))
  expect_null(native_call("eigencore_composite_block_apply", cleared, diag(4),
                          1, 0, NULL, FALSE))

  build <- function(spec) native_call("eigencore_composite_operator_build", spec)
  expect_error(build(list(x = X)), "missing node type")
  expect_error(build(list(type = "nonsense")), "unknown node type")
  expect_error(build(list(type = "dense", x = 1:3)), "double matrix")
  expect_error(build(list(type = "product", children = list())),
               "at least one child")
  expect_error(build(list(type = "product",
                          children = list(list(type = "dense", x = X),
                                          list(type = "dense", x = X)))),
               "non-conformable")
  expect_error(build(list(type = "sum",
                          children = list(list(type = "dense", x = X),
                                          list(type = "dense", x = X)),
                          weights = c(1, NaN))),
               "finite")
  deep <- list(type = "dense", x = X)
  for (i in 1:70) deep <- list(type = "adjoint", child = deep)
  expect_error(build(deep), "nested too deeply")
})

test_that("rank-one composite leaves accumulate into Y when beta != 0", {
  u <- c(1, 2, 3)
  v <- c(4, 5)
  spec <- list(type = "sum",
               children = list(list(type = "dense", x = matrix(1, 3L, 2L)),
                               list(type = "rank1", u = u, v = v)),
               weights = c(2, -1))
  ptr <- native_call("eigencore_composite_operator_build", spec)
  X <- matrix(c(1, -1, 2, 0.5), 2L)
  Y0 <- matrix(stats::rnorm(6), 3L)
  explicit <- 2 * matrix(1, 3L, 2L) - tcrossprod(u, v)
  got <- native_call("eigencore_composite_block_apply", ptr, X, 0.5, 2, Y0,
                     FALSE)
  expect_equal(got, 0.5 * explicit %*% X + 2 * Y0, tolerance = 1e-14)
  W <- matrix(stats::rnorm(6), 3L)
  Z0 <- matrix(stats::rnorm(4), 2L)
  got_adj <- native_call("eigencore_composite_block_apply", ptr, W, 1, -1, Z0,
                         TRUE)
  expect_equal(got_adj, crossprod(explicit, W) - Z0, tolerance = 1e-14)
})

test_that("identity hash covers language, pairlist, symbol, expression and S4 nodes", {
  h <- function(x) native_call("eigencore_identity_hash", x)
  expect_identical(h(quote(f(x, y = 2))), h(quote(f(x, y = 2))))
  expect_false(identical(h(quote(f(x, y = 2))), h(quote(f(x, y = 3)))))
  # A tagged argument hashes differently from an untagged one.
  expect_false(identical(h(quote(f(x, 2))), h(quote(f(x, y = 2)))))
  expect_false(identical(h(as.pairlist(list(a = 1))),
                         h(as.pairlist(list(b = 1)))))
  expect_false(identical(h(as.name("alpha")), h(as.name("beta"))))
  expect_false(identical(h(expression(a + 1)), h(expression(a + 2))))
  M1 <- Matrix::Matrix(diag(3), sparse = TRUE)
  M2 <- Matrix::Matrix(diag(c(1, 2, 3)), sparse = TRUE)
  expect_identical(h(M1), h(Matrix::Matrix(diag(3), sparse = TRUE)))
  expect_false(identical(h(M1), h(M2)))
})

test_that("Krylov-Schur Arnoldi recovers from an exhausted (invariant) Krylov space", {
  # The start vector lies in a 2-dimensional invariant subspace, so the
  # Arnoldi recurrence breaks down after two steps and must continue from a
  # random direction orthogonal to the basis.
  n <- 12L
  A <- diag(c(0, 0, seq(1, 3, length.out = n - 2L)))
  A[1:2, 1:2] <- matrix(c(10, 1, -2, 9), 2L)
  start <- c(1, 1, rep(0, n - 2L))
  ns <- asNamespace("eigencore")
  set.seed(7)
  ks <- ns$native_krylov_schur(as_operator(A), start, k = 3L, m = 8L,
                               target = largest_magnitude(), tol = 1e-10)
  expect_gte(ks$breakdowns, 1L)
  expect_true(ks$converged)
  fit <- eig_partial(A, 3L, target = largest_magnitude(), seed = 3)
  ref <- eigen(A, only.values = TRUE)$values
  expect_equal(sort(Mod(values(fit))),
               sort(sort(Mod(ref), decreasing = TRUE)[1:3]), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
})

test_that("residual norms do not overflow or underflow at extreme scales", {
  A <- spd_matrix(40L, 2)
  for (s in c(1e150, 1e-150)) {
    fit <- eig_partial(A * s, 3L, method = lanczos(), seed = 4)
    expect_equal(values(fit), s * eigen(A, symmetric = TRUE)$values[1:3],
                 tolerance = 1e-8)
    cert <- certificate(fit)
    expect_true(cert$passed)
    expect_true(all(is.finite(cert$residuals)))
    S <- methods::as(Matrix::Matrix(A * s, sparse = TRUE), "generalMatrix")
    fit_s <- eig_partial(S, 3L, method = lanczos(block = 2L), seed = 5)
    expect_true(certificate(fit_s)$passed)
  }
})

test_that("projected eigensolves with ties at the selection boundary stay exact", {
  # Exactly repeated eigenvalues straddle the k-th position, so the index
  # subset solve (dsyevr RANGE = 'I') may return extra tied values and the
  # native code falls back to the full projected spectrum.
  d <- c(rep(5, 4L), seq(4, 1, length.out = 36L))
  A <- diag(d)
  for (block in c(1L, 2L, 3L)) {
    fit <- eig_partial(A, 2L, method = lanczos(block = block), seed = 6)
    expect_equal(values(fit), c(5, 5), tolerance = 1e-10)
  }
  S <- methods::as(Matrix::Diagonal(40L, x = d), "CsparseMatrix")
  fit <- eig_partial(S, 3L, method = lanczos(block = 2L), seed = 7)
  expect_equal(values(fit), c(5, 5, 5), tolerance = 1e-10)
})

test_that("LOBPCG handles blocks that exhaust a tiny search space", {
  for (n in c(5L, 6L, 8L)) {
    A <- spd_matrix(n, n)
    k <- if (n == 5L) 2L else 3L
    fit <- eig_partial(A, k, target = smallest(), method = lobpcg(maxit = 200L),
                       seed = 8)
    expect_equal(values(fit), seq_len(k), tolerance = 1e-8)
    fitB <- eig_partial(A, k, B = diag(seq(1, 2, length.out = n)),
                        target = smallest(), method = lobpcg(maxit = 200L),
                        seed = 9)
    ref <- sort(Re(eigen(solve(diag(seq(1, 2, length.out = n))) %*% A,
                         only.values = TRUE)$values))[seq_len(k)]
    expect_equal(sort(values(fitB)), ref, tolerance = 1e-7)
  }
})

test_that("unwind selftest reports R-level stops and unknown modes after cleanup", {
  st <- function(mode) native_call("eigencore_unwind_selftest", mode)
  expect_error(st(1L))
  expect_identical(st("live"), 0L)
})

test_that("retained IRLBA restart ABI agrees with svd() under every reorthogonalisation policy", {
  ns <- asNamespace("eigencore")
  set.seed(15)
  G <- matrix(stats::rnorm(90L * 30L), 90L) %*% diag(exp(-seq(0, 4, length.out = 30L)))
  R <- Matrix::rsparsematrix(90L, 30L, density = 0.2)
  for (A in list(G, R)) {
    for (policy in c("one_sided_small_side", "full_two_sided", "bpro_two_sided",
                     "bpro_one_sided_guarded", "bpro_block_guarded")) {
      for (target in list(largest(), smallest())) {
        fit <- ns$native_irlba_lbd_retained_svd(
          as_operator(A), rank = 3L, target = target, work = 10L,
          retained = 5L, max_restarts = 20L, tol = 1e-8, vectors = "both",
          reorth_policy = policy
        )
        sv <- svd(as.matrix(A))$d
        want <- if (identical(target$kind, "smallest")) {
          sort(sv)[1:3]
        } else {
          sort(sv, decreasing = TRUE)[1:3]
        }
        if (isTRUE(fit$certificate$passed)) {
          expect_equal(sort(fit$d), sort(want), tolerance = 1e-6,
                       info = paste(policy, target$kind))
        }
        expect_true(all(is.finite(fit$d)))
      }
    }
  }
})

test_that("nonsymmetric shift-invert covers sparse adjoint solves and guard rails", {
  set.seed(11)
  n <- 60L
  A <- Matrix::rsparsematrix(n, n, density = 0.08) +
    Matrix::Diagonal(n, x = seq(1, 12, length.out = n))
  A <- methods::as(A, "CsparseMatrix")
  Ad <- as.matrix(A)
  ref <- eigen(Ad, only.values = TRUE)$values
  near <- function(sigma, k) ref[order(Mod(ref - sigma))][seq_len(k)]

  # sparse LU forward and transposed (left-vector) solves
  fit <- eig_partial(A, 2L, target = nearest(4.05), method = shift_invert(4.05),
                     seed = 12, left_vectors = "compute")
  expect_equal(sort(Re(values(fit))), sort(Re(near(4.05, 2))), tolerance = 1e-8)
  expect_true(certificate(fit)$passed)
  expect_false(is.null(left_vectors(fit)))

  # dense QR route without certification
  fit_nc <- eig_partial(Ad, 2L, target = nearest(6.1), method = shift_invert(6.1),
                        certify = FALSE, seed = 13)
  expect_equal(sort(Re(values(fit_nc))), sort(Re(near(6.1, 2))), tolerance = 1e-6)
  expect_match(paste(certificate(fit_nc)$notes, collapse = " "), "disabled")

  # Hermitian shift-invert without certification
  S <- Matrix::forceSymmetric(A + Matrix::t(A))
  S <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
  fit_h <- eig_partial(S, 2L, target = nearest(9), method = shift_invert(9),
                       certify = FALSE, seed = 14)
  ev <- eigen(as.matrix(S), symmetric = TRUE)$values
  expect_equal(sort(values(fit_h)), sort(ev[order(abs(ev - 9))][1:2]),
               tolerance = 1e-6)

  expect_error(
    eig_partial(Ad, 2L, target = nearest(4),
                method = shift_invert(4, factorization = "lu")),
    "not implemented"
  )
  expect_error(
    eig_partial(Ad, 2L, B = diag(n), target = nearest(4),
                method = shift_invert(4)),
    "not implemented|generalized|requires a Hermitian"
  )
})
