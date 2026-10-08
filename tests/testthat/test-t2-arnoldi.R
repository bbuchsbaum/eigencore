## Tranche 2 (review P1, P12, C16, C25): native Krylov-Schur Arnoldi for
## nonsymmetric partial eigenproblems.

# A real nonsymmetric matrix with a well-separated spectrum: real eigenvalues
# plus conjugate pairs with distinct real and imaginary parts, hidden by a
# random sparse similarity-like perturbation.
t2_test_matrix <- function(n, seed = 1L) {
  set.seed(seed)
  A <- Matrix::Matrix(0, n, n, sparse = TRUE)
  pairs <- 4L
  for (j in seq_len(pairs)) {
    i <- 2L * j - 1L
    a <- c(1.5, -2.5, 0.5, -0.7)[[j]]
    b <- c(4, 3, 2.2, 1.4)[[j]]
    A[i, i] <- a
    A[i + 1L, i + 1L] <- a
    A[i, i + 1L] <- b
    A[i + 1L, i] <- -b
  }
  rest <- (2L * pairs + 1L):n
  # Real eigenvalues: interior ones in [-1, 1], plus 3.5 and -3.2 at the ends.
  A[cbind(rest, rest)] <- c(3.5, seq(-1, 1, length.out = length(rest) - 2L), -3.2)
  noise <- Matrix::rsparsematrix(n, n, density = 3 / n) * 0.02
  methods::as(A + noise, "CsparseMatrix")
}

# Operator explicitly marked general. The automatic symmetry probes use
# absolute tolerances (the dense probe floors its scale at 1 in
# src/small_dense.cpp; Matrix::isSymmetric(tol =) is absolute for sparse
# input), so a tiny nonsymmetric matrix would otherwise be classified as
# Hermitian before the Arnoldi path is ever reached.
t2_general <- function(A) {
  op <- eigencore:::as_operator(A)
  op$structure <- general()
  op
}

t2_target_key <- function(values, target) {
  switch(
    target$kind,
    largest_real = Re(values),
    smallest_real = -Re(values),
    largest_magnitude = Mod(values),
    largest_imaginary = Im(values),
    smallest_imaginary = -Im(values),
    stop("unexpected target")
  )
}

t2_expect_matches_eigen <- function(fit, A, target, k, tol = 1e-7) {
  ev <- eigen(as.matrix(A), only.values = TRUE)$values
  key <- t2_target_key(ev, target)
  expected <- sort(key, decreasing = TRUE)[seq_len(k)]
  got <- sort(t2_target_key(values(fit), target), decreasing = TRUE)
  expect_equal(got, expected, tolerance = tol)
  # Every returned value is an eigenvalue of A.
  for (v in values(fit)) {
    expect_lt(min(Mod(ev - v)), tol * max(1, Mod(v)))
  }
}

test_that("Krylov-Schur matches eigen() for every target on sparse and dense inputs", {
  n <- 120L
  A <- t2_test_matrix(n, seed = 11L)
  targets <- list(
    largest_real(), smallest_real(), largest_magnitude(),
    largest_imaginary(), smallest_imaginary()
  )
  for (input in list(A, as.matrix(A))) {
    for (target in targets) {
      fit <- eig_partial(
        input, k = 4L, target = target, tol = 1e-10, seed = 3,
        allow_dense_fallback = "never"
      )
      label <- paste(class(input)[[1L]], target$kind)
      expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label(),
                   label = label)
      expect_true(fit$restart$krylov_schur, label = label)
      expect_true(fit$certificate$passed, label = label)
      t2_expect_matches_eigen(fit, input, target, 4L)
    }
  }
})

test_that("Krylov-Schur keeps conjugate pairs intact", {
  n <- 100L
  A <- t2_test_matrix(n, seed = 5L)
  # The two largest-magnitude eigenvalues form the pair ~1.5 +/- 4i.
  fit <- eig_partial(A, k = 2L, target = largest_magnitude(), tol = 1e-10,
                     seed = 9, allow_dense_fallback = "never")
  vals <- values(fit)
  expect_true(is.complex(vals))
  expect_equal(sort(Im(vals)), c(-1, 1) * abs(Im(vals[[1L]])), tolerance = 1e-10)
  expect_equal(Re(vals[[1L]]), Re(vals[[2L]]), tolerance = 1e-10)
  # Matching eigenvectors are complex conjugates (up to phase).
  V <- vectors(fit)
  overlap <- Mod(sum(Conj(V[, 1L]) * Conj(V[, 2L])))
  expect_equal(overlap, 1, tolerance = 1e-6)
  expect_true(fit$certificate$passed)

  # Imaginary-part targets: partners of the wanted values do not crowd out
  # the next wanted value.
  fit_li <- eig_partial(A, k = 2L, target = largest_imaginary(), tol = 1e-10,
                        seed = 9, allow_dense_fallback = "never")
  expect_equal(sort(Im(values(fit_li)), decreasing = TRUE), c(4, 3),
               tolerance = 2e-2)
  expect_true(all(Im(values(fit_li)) > 0))
  expect_true(fit_li$certificate$passed)
})

test_that("Krylov-Schur is scale invariant for tiny and huge matrices", {
  n <- 80L
  A <- t2_test_matrix(n, seed = 21L)
  base <- eig_partial(A, k = 4L, target = largest_magnitude(), tol = 1e-10,
                      seed = 4, allow_dense_fallback = "never")
  for (s in c(1e-12, 1e12)) {
    for (input in list(t2_general(A * s), t2_general(as.matrix(A) * s))) {
      fit <- eig_partial(input, k = 4L, target = largest_magnitude(),
                         tol = 1e-10, seed = 4, allow_dense_fallback = "never")
      expect_true(fit$certificate$passed, label = paste("scale", s))
      expect_equal(sort(Mod(values(fit))), sort(Mod(values(base))) * s,
                   tolerance = 1e-8)
      # Complex pairs stay complex at every scale (C25: the realness test is
      # relative, not an absolute sqrt(tol)).
      expect_true(is.complex(values(fit)))
      expect_equal(sum(abs(Im(values(fit))) > 0), 4L)
    }
  }
})

test_that("relative realness test keeps tiny genuine imaginary parts (C25)", {
  A <- t2_general(1e-9 * matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 0.5), 3, 3))
  fit <- eig_partial(A, k = 1L, target = largest_imaginary(), tol = 1e-8,
                     allow_dense_fallback = "never")
  expect_true(is.complex(values(fit)))
  expect_equal(Im(values(fit)), 1e-9, tolerance = 1e-6)

  expect_true(eigencore:::arnoldi_real_values(complex(real = 2, imaginary = 1e-17)))
  expect_false(eigencore:::arnoldi_real_values(complex(real = 2e-12, imaginary = 1e-13)))
  expect_true(all(eigencore:::arnoldi_real_values(c(3, 1))))
})

test_that("Arnoldi breakdown thresholds are relative to the operator scale (C16)", {
  set.seed(8)
  n <- 40L
  A <- matrix(rnorm(n * n), n, n) * 1e-20
  op <- eigencore:::as_operator(A)
  start <- rnorm(n)
  cycle <- eigencore:::native_arnoldi_cycle(op, start, 10L)
  # An absolute 100 * eps threshold would declare breakdown at step 1.
  expect_equal(cycle$iterations, 10L)
  Vm <- cycle$V[, 1:10]
  expect_equal(crossprod(Vm), diag(10), tolerance = 1e-12)

  ref <- eigencore:::reference_arnoldi_cycle(op, start, 10L)
  expect_equal(ref$iterations, 10L)

  # A genuine invariant subspace is still detected at any scale.
  D <- diag(c(3, 2, 1, rep(0, n - 3L))) * 1e-20
  cyc_d <- eigencore:::native_arnoldi_cycle(eigencore:::as_operator(D),
                                            c(1, 1, 1, rep(0, n - 3L)), 10L)
  expect_equal(cyc_d$iterations, 3L)
  expect_equal(cyc_d$H[4L, 3L], 0)
})

test_that("native Krylov-Schur returns a valid Krylov-Schur decomposition", {
  n <- 150L
  A <- t2_test_matrix(n, seed = 31L)
  op <- eigencore:::as_operator(A)
  set.seed(2)
  ks <- eigencore:::native_krylov_schur(op, rnorm(n), k = 4L, m = 20L,
                                        target = largest_real(), tol = 1e-10)
  p <- ks$iterations
  expect_true(ks$converged)
  expect_gte(p, 4L)
  expect_lt(p, 20L)
  V <- ks$V
  expect_equal(dim(V), c(n, p + 1L))
  expect_equal(crossprod(V), diag(p + 1L), tolerance = 1e-12)
  AV <- as.matrix(A %*% V[, seq_len(p)])
  expect_lt(max(abs(AV - V %*% ks$H)), 1e-12 * max(abs(as.matrix(A))) * n)
  # The leading block is quasi-triangular (real Schur form).
  Tp <- ks$H[seq_len(p), seq_len(p)]
  expect_true(all(Tp[row(Tp) > col(Tp) + 1L] == 0))
  expect_true(all(c("apply", "orthogonalization", "projected_schur", "restart") %in%
                    names(ks$stage_seconds)))
})

test_that("Ritz extraction forms only the requested vectors", {
  n <- 60L
  A <- t2_test_matrix(n, seed = 41L)
  op <- eigencore:::as_operator(A)
  set.seed(3)
  cycle <- eigencore:::native_arnoldi_cycle(op, rnorm(n), 20L)
  coef <- eigencore:::native_arnoldi_ritz_coefficients(cycle)
  expect_equal(dim(coef$coefficients), c(20L, 20L))
  all_ritz <- eigencore:::native_arnoldi_projected_ritz(cycle)
  idx <- c(3L, 1L)
  sel <- eigencore:::native_arnoldi_ritz_vectors(cycle, coef$coefficients[, idx])
  expect_equal(dim(sel), c(n, 2L))
  expect_equal(sel, all_ritz$vectors[, idx], tolerance = 1e-12)
  expect_equal(coef$values, all_ritz$values)
})

test_that("dense inputs restart in a small subspace and match eigen()", {
  set.seed(51)
  n <- 300L
  A <- matrix(rnorm(n * n), n, n) / sqrt(n)
  diag(A) <- diag(A) + c(6, 5, 4, rep(0, n - 3L))
  fit <- eig_partial(A, k = 3L, target = largest_real(), tol = 1e-10,
                     seed = 1, allow_dense_fallback = "never")
  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_lt(fit$plan$controls$max_subspace, n)
  expect_true(fit$certificate$passed)
  t2_expect_matches_eigen(fit, A, largest_real(), 3L)
})

test_that("sparse complex certificates and left vectors stay native and sparse", {
  n <- 200L
  A <- t2_test_matrix(n, seed = 61L)
  sparse_as_matrix_calls <- 0L
  invisible(trace(
    "as.matrix",
    tracer = quote({
      if (inherits(x, "sparseMatrix")) {
        sparse_as_matrix_calls <<- sparse_as_matrix_calls + 1L
      }
    }),
    print = FALSE,
    where = asNamespace("Matrix")
  ))
  on.exit(untrace("as.matrix", where = asNamespace("Matrix")), add = TRUE)
  fit <- eig_partial(A, k = 4L, target = largest_magnitude(), tol = 1e-10,
                     seed = 7, allow_dense_fallback = "never")
  expect_equal(sparse_as_matrix_calls, 0L)
  expect_true(is.complex(values(fit)))
  expect_true(fit$certificate$passed)
  expect_true(fit$left_certificate$passed)
  expect_lt(max(Mod(fit$biorthogonality - diag(4L))), 1e-6)
  # The adjoint of a dgCMatrix operator resolves to the native CSC kernel.
  adj <- eigencore:::adjoint(eigencore:::as_operator(A))
  expect_true(eigencore:::native_arnoldi_available(adj))
})

test_that("adjoint solve targets the eigenvalues the right solve returned", {
  vals <- c(3 + 1i, 3 - 1i, -2)
  ritz <- c(-2.0001, 5, 3 - 1.0001i, 0.1, 3 + 0.9999i)
  expect_setequal(
    eigencore:::arnoldi_order_indices(ritz, largest_real(), vals)[1:3],
    c(1L, 3L, 5L)
  )
  expect_equal(eigencore:::arnoldi_order_indices(ritz, largest_real()),
               order(Re(ritz), decreasing = TRUE))

  # Clustered random spectrum (circular law): the left solve must find the
  # same eigenvalues as the right solve, so the left certificate passes.
  set.seed(1)
  n <- 1500L
  A <- Matrix::rsparsematrix(n, n, density = 10 / n)
  fit <- eig_partial(A, k = 6L, target = largest_magnitude(), seed = 2,
                     allow_dense_fallback = "never")
  expect_true(fit$certificate$passed)
  expect_true(fit$left_certificate$passed)
})

test_that("nonsymmetric sparse eigs performance smoke (n = 3000, k = 6)", {
  skip_on_cran()
  set.seed(1)
  n <- 3000L
  A <- Matrix::rsparsematrix(n, n, density = 10 / n)
  elapsed <- system.time(
    fit <- eigencore::eigs(A, 6)
  )[["elapsed"]]
  expect_true(fit$certificate$passed)
  expect_equal(fit$nconv, 6L)
  # The pre-Krylov-Schur implementation took ~10 s here and failed to
  # converge; the restarted solver takes about one second.
  expect_lt(elapsed, 30)
})
