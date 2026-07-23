# Warm-start (initial_subspace) contract and adversarial coverage.
#
# Scope: planning/prd.json `initial_subspace_contract`. The public
# initial_subspace argument seeds standard real Hermitian Lanczos: the native
# paths on explicit dense double / dgCMatrix operators and the matrix-free
# reference Hermitian Lanczos path. The subspace is only a hint; every solve
# recomputes a fresh current-operator certificate. Correctness is checked
# against an independent dense oracle (base eigen()).

build_dense_sym <- function(n, seed = 3L) {
  random_symmetric_with_spectrum(n, pattern = "uniform", seed = seed)
}

build_sparse_sym <- function(n) {
  # Explicit general (dgCMatrix) storage that is numerically symmetric: supply
  # both off-diagonal triangles so no dsCMatrix coercion is involved.
  i <- c(seq_len(n), seq_len(n - 1L), 2:n)
  j <- c(seq_len(n), 2:n, seq_len(n - 1L))
  x <- c(seq(n, 1), rep(0.4, n - 1L), rep(0.4, n - 1L))
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
}

# Perturbed leading eigenvectors: a realistic "good" warm start.
warm_from_truth <- function(A, k, noise = 1e-2, seed = 1L) {
  V <- eigen(as.matrix(A), symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  set.seed(seed)
  V + matrix(stats::rnorm(nrow(V) * k, sd = noise), nrow(V), k)
}

block_method <- function(k, n) {
  # max_subspace < n keeps both cold and warm on the iterative path (no exact
  # full-subspace dsyev shortcut), so the supplied start is genuinely consumed.
  lanczos(block = k, max_subspace = min(n - 1L, 8L * k), max_restarts = 200L)
}

# --------------------------------------------------------------------------
# P1-1: NULL preserves cold behavior exactly
# --------------------------------------------------------------------------

test_that("initial_subspace = NULL is regression-identical to cold", {
  n <- 60L; k <- 4L
  A <- build_dense_sym(n)
  cold1 <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                       seed = 11)
  cold2 <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                       seed = 11, initial_subspace = NULL)
  expect_identical(cold1$start_source, "cold")
  expect_identical(cold2$start_source, "cold")
  expect_equal(values(cold1), values(cold2))
  expect_equal(cold1$matvecs, cold2$matvecs)
  expect_identical(cold1$method, cold2$method)
})

# --------------------------------------------------------------------------
# P1-1: supported native plans consume a start and certify (dense + CSC)
# --------------------------------------------------------------------------

test_that("dense block warm start certifies against the dense oracle", {
  n <- 80L; k <- 5L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  warm <- warm_from_truth(A, k)
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 7, initial_subspace = warm)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
  expect_identical(fit$start_source, "user_supplied")
  expect_identical(fit$method, "native block Hermitian Lanczos (thick restart, locking)")
  prov <- fit$initial_subspace
  expect_equal(prov$supplied, k)
  expect_equal(prov$accepted, k)
  expect_equal(prov$rejected, 0L)
  expect_equal(prov$augmented, 0L)
})

test_that("sparse dgCMatrix block warm start certifies against the dense oracle", {
  n <- 70L; k <- 4L
  sp <- build_sparse_sym(n)
  expect_s4_class(sp, "dgCMatrix")
  truth <- sort(eigen(as.matrix(sp), symmetric = TRUE)$values,
                decreasing = TRUE)[seq_len(k)]
  warm <- eigen(as.matrix(sp), symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  fit <- eig_partial(sp, k = k, target = largest(), method = block_method(k, n),
                     seed = 5, initial_subspace = warm)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
  expect_identical(fit$start_source, "user_supplied")
})

test_that("scalar Lanczos fits a wide subspace to a single-column start", {
  n <- 60L; k <- 3L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  # method = lanczos() defaults to block = 1 -> native scalar path, width 1.
  warm <- eigen(A, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  fit <- eig_partial(A, k = k, target = largest(),
                     method = lanczos(max_subspace = 40L, max_restarts = 200L),
                     seed = 9, initial_subspace = warm)
  expect_identical(fit$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
  prov <- fit$initial_subspace
  expect_equal(prov$supplied, k)
  expect_equal(prov$accepted, 1L)   # width-1 start block
  expect_equal(prov$rejected, k - 1L)
  expect_equal(prov$augmented, 0L)
  expect_true(prov$compressed)      # rotation of the full k-basis, not X[, 1]
})

test_that("compression fit mixes every supplied direction into the start", {
  # rank > width: the start block must be a rotation of the whole accepted
  # basis, so each supplied direction carries generic (non-negligible) weight.
  n <- 50L; k <- 4L
  set.seed(21)
  basis <- qr.Q(qr(matrix(stats::rnorm(n * k), n, k)))
  set.seed(33)
  prep <- eigencore:::prepare_initial_subspace(basis, n = n, width = 1L,
                                               solver_tol = 1e-8)
  expect_true(prep$compressed)
  expect_equal(prep$accepted, 1L)
  expect_equal(dim(prep$start), c(n, 1L))
  overlaps <- abs(drop(crossprod(basis, prep$start)))
  # A Gaussian rotation gives each direction weight ~ 1/sqrt(k); anything
  # above noise level proves no direction was truncated away.
  expect_true(all(overlaps > 1e-2))

  # width >= rank keeps the basis as-is (no compression).
  set.seed(33)
  wide <- eigencore:::prepare_initial_subspace(basis, n = n, width = k,
                                               solver_tol = 1e-8)
  expect_false(wide$compressed)
  expect_equal(wide$accepted, k)
})

# --------------------------------------------------------------------------
# P1-2: adversarial starts still recover the certified target spectrum
# --------------------------------------------------------------------------

test_that("badly scaled start does not corrupt the certified answer", {
  n <- 80L; k <- 4L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  warm <- warm_from_truth(A, k) * 1e8   # huge column scale
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 2, initial_subspace = warm)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
  expect_equal(fit$initial_subspace$accepted, k)
})

test_that("duplicate columns are rejected and deterministically augmented", {
  n <- 80L; k <- 4L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  V <- eigen(A, symmetric = TRUE)$vectors
  dup <- cbind(V[, 1], V[, 1], V[, 2], V[, 2])   # 4 columns, numerical rank 2
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 4, initial_subspace = dup)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
  prov <- fit$initial_subspace
  expect_equal(prov$supplied, 4L)
  expect_equal(prov$rank, 2L)
  expect_equal(prov$accepted, 2L)
  expect_equal(prov$rejected, 2L)
  expect_equal(prov$augmented, 2L)   # filled up to block width
  expect_identical(fit$start_source, "user_supplied")
})

test_that("rank-deficient (too-few) start is augmented up to block width", {
  n <- 80L; k <- 5L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  warm <- eigen(A, symmetric = TRUE)$vectors[, seq_len(2L), drop = FALSE]  # only 2 cols
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 6, initial_subspace = warm)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
  prov <- fit$initial_subspace
  expect_equal(prov$supplied, 2L)
  expect_equal(prov$accepted, 2L)
  expect_equal(prov$augmented, k - 2L)
})

test_that("zero-overlap start (non-target invariant subspace) still finds the target", {
  n <- 80L; k <- 4L
  A <- build_dense_sym(n)
  E <- eigen(A, symmetric = TRUE)
  truth <- E$values[seq_len(k)]
  # Exact eigenvectors for the NEXT block down: zero overlap with the target.
  bad <- E$vectors[, (k + 1L):(2L * k), drop = FALSE]
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 1, initial_subspace = bad)
  expect_equal(sort(values(fit), decreasing = TRUE), sort(truth, decreasing = TRUE),
               tolerance = 1e-6)
  expect_certificate_clean(fit, tol = 1e-7)
  expect_equal(fit$initial_subspace$accepted, k)  # full rank, just wrong directions
})

test_that("non-target invariant-subspace start never yields a wrong certified answer", {
  # An exact non-target invariant subspace would trap Krylov iteration and, at
  # a minimal max_subspace, certify the non-target eigenpairs it was handed. The
  # blended random escape must prevent that: the result is either the certified
  # target spectrum or an honest non-convergence, but never wrong-and-certified.
  n <- 40L; k <- 3L
  A <- build_dense_sym(n)
  E <- eigen(A, symmetric = TRUE)
  truth <- sort(E$values, decreasing = TRUE)[seq_len(k)]
  bad <- E$vectors[, (k + 1L):(2L * k), drop = FALSE]

  # Minimal subspace budget (the dangerous regime).
  tight <- eig_partial(A, k = k, target = largest(),
                       method = lanczos(block = k, max_subspace = 2L * k,
                                        max_restarts = 200L),
                       seed = 1, initial_subspace = bad)
  tight_vals <- sort(values(tight), decreasing = TRUE)
  target_match <- isTRUE(all.equal(tight_vals, truth, tolerance = 1e-5))
  wrong_certified <- certificate(tight)$passed && !target_match
  expect_false(wrong_certified)

  # Adequate budget: must fully recover and certify the target.
  ample <- eig_partial(A, k = k, target = largest(),
                       method = lanczos(block = k, max_subspace = 8L * k,
                                        max_restarts = 200L),
                       seed = 1, initial_subspace = bad)
  expect_equal(sort(values(ample), decreasing = TRUE), truth, tolerance = 1e-6)
  expect_certificate_clean(ample, tol = 1e-7)
  expect_true(tight$initial_subspace$escape_blended)
})

test_that("partially overlapping start recovers the certified target spectrum", {
  n <- 80L; k <- 4L
  A <- build_dense_sym(n)
  E <- eigen(A, symmetric = TRUE)
  truth <- E$values[seq_len(k)]
  # Half target directions, half wrong directions.
  warm <- cbind(E$vectors[, 1:2], E$vectors[, (k + 1L):(k + 2L)])
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 8, initial_subspace = warm)
  expect_equal(sort(values(fit), decreasing = TRUE), sort(truth, decreasing = TRUE),
               tolerance = 1e-6)
  expect_certificate_clean(fit, tol = 1e-7)
})

test_that("all-zero start is an honest degenerate warm start that still certifies", {
  n <- 60L; k <- 4L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 3, initial_subspace = matrix(0, n, k))
  expect_identical(fit$start_source, "user_supplied_degenerate")
  prov <- fit$initial_subspace
  expect_equal(prov$accepted, 0L)
  expect_equal(prov$rejected, k)
  expect_equal(prov$augmented, k)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  expect_certificate_clean(fit)
})

# --------------------------------------------------------------------------
# P1-2: input validation
# --------------------------------------------------------------------------

test_that("dimension mismatch is a clear error", {
  n <- 60L; k <- 4L
  A <- build_dense_sym(n)
  expect_error(
    eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                initial_subspace = matrix(1, n + 3L, k)),
    "operator domain has dimension 60"
  )
})

test_that("non-finite start entries are rejected", {
  n <- 60L; k <- 4L
  A <- build_dense_sym(n)
  warm <- warm_from_truth(A, k)
  warm[1, 1] <- NA_real_
  expect_error(
    eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                initial_subspace = warm),
    "finite numeric"
  )
})

test_that("a plain numeric vector is accepted as a single start direction", {
  n <- 50L; k <- 1L
  A <- build_dense_sym(n)
  truth <- max(eigen(A, symmetric = TRUE)$values)
  v <- eigen(A, symmetric = TRUE)$vectors[, 1]
  fit <- eig_partial(A, k = k, target = largest(),
                     method = lanczos(max_subspace = 30L, max_restarts = 100L),
                     seed = 1, initial_subspace = v)
  expect_equal(values(fit), truth, tolerance = 1e-7)
  expect_equal(fit$initial_subspace$supplied, 1L)
})

# --------------------------------------------------------------------------
# Matrix-free reference Hermitian Lanczos consumes a start
# --------------------------------------------------------------------------

matrix_free_hermitian <- function(A) {
  n <- nrow(A)
  linear_operator(
    dim = c(n, n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian(),
    name = "matrix-free test operator"
  )
}

test_that("matrix-free warm start is consumed, certified, and cheaper than cold", {
  n <- 80L; k <- 4L
  A <- build_dense_sym(n)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  op <- matrix_free_hermitian(A)

  cold <- eig_partial(op, k = k, target = largest(), method = lanczos(),
                      maxit = n, seed = 12)
  expect_identical(cold$method,
                   "reference Hermitian Lanczos (prototype/oracle fallback)")
  expect_identical(cold$start_source, "cold")
  expect_equal(sort(values(cold), decreasing = TRUE), truth, tolerance = 1e-7)

  # NULL is regression-identical to cold on the matrix-free path too.
  cold2 <- eig_partial(op, k = k, target = largest(), method = lanczos(),
                       maxit = n, seed = 12, initial_subspace = NULL)
  expect_equal(values(cold), values(cold2))
  expect_equal(cold$matvecs, cold2$matvecs)

  # Exact previous-solve eigenvectors: the spectral-continuation access pattern.
  warm_basis <- eigen(A, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  warm <- eig_partial(op, k = k, target = largest(), method = lanczos(),
                      maxit = n, seed = 12, initial_subspace = warm_basis)
  expect_identical(warm$method,
                   "reference Hermitian Lanczos (prototype/oracle fallback)")
  expect_identical(warm$start_source, "user_supplied")
  expect_equal(sort(values(warm), decreasing = TRUE), truth, tolerance = 1e-7)
  # Matrix-free certificates withhold `passed` (stochastic norm estimate);
  # assert the underlying evidence instead.
  cert <- certificate(warm)
  expect_true(all(cert$converged))
  expect_lte(cert$max_backward_error, cert$tolerance)
  expect_true(cert$orthogonality_passed)
  prov <- warm$initial_subspace
  expect_equal(prov$supplied, k)
  expect_equal(prov$accepted, 1L)   # scalar reference path, width 1
  expect_true(prov$compressed)
  expect_true(prov$escape_blended)
  expect_lt(warm$matvecs, cold$matvecs)
})

test_that("matrix-free warm start accepts a perturbed continuation basis", {
  n <- 60L; k <- 3L
  A <- build_dense_sym(n, seed = 8L)
  truth <- sort(eigen(A, symmetric = TRUE)$values, decreasing = TRUE)[seq_len(k)]
  op <- matrix_free_hermitian(A)
  warm <- warm_from_truth(A, k, noise = 1e-2)
  fit <- eig_partial(op, k = k, target = largest(), method = lanczos(),
                     maxit = n, seed = 4, initial_subspace = warm)
  expect_equal(sort(values(fit), decreasing = TRUE), truth, tolerance = 1e-7)
  cert <- certificate(fit)
  expect_true(all(cert$converged))
  expect_lte(cert$max_backward_error, cert$tolerance)
  expect_true(cert$orthogonality_passed)
  expect_identical(fit$start_source, "user_supplied")
})

# --------------------------------------------------------------------------
# P1-1: unsupported plans reject initial_subspace explicitly
# --------------------------------------------------------------------------

test_that("dense LAPACK fallback plan rejects initial_subspace", {
  A <- build_dense_sym(6L)
  expect_error(
    eig_partial(A, k = 5L, target = largest(),
                initial_subspace = matrix(stats::rnorm(30), 6L, 5L)),
    "does not consume a starting subspace"
  )
})

test_that("generalized (metric B) plan rejects initial_subspace", {
  n <- 40L; k <- 3L
  A <- build_dense_sym(n)
  B <- crossprod(matrix(stats::rnorm(n * n), n)) + diag(n)  # SPD
  expect_error(
    eig_partial(A, B = B, k = k, target = largest(), method = lanczos(block = k),
                initial_subspace = matrix(stats::rnorm(n * k), n, k)),
    "does not consume a starting subspace"
  )
})

test_that("shift-invert transform plan rejects initial_subspace", {
  n <- 40L; k <- 3L
  A <- build_dense_sym(n)
  expect_error(
    eig_partial(A, k = k, target = nearest(0), method = shift_invert(sigma = 0),
                initial_subspace = matrix(stats::rnorm(n * k), n, k)),
    "does not consume a starting subspace"
  )
})

# --------------------------------------------------------------------------
# Contract: no new exported symbols; diagnostics surface provenance
# --------------------------------------------------------------------------

test_that("warm-start helpers are internal, not exported", {
  exports <- getNamespaceExports("eigencore")
  expect_false("prepare_initial_subspace" %in% exports)
  expect_false("validate_initial_subspace_plan_support" %in% exports)
  expect_false("warm_start_plan_consumes_start" %in% exports)
  expect_false("warm_start_cold_provenance" %in% exports)
})

test_that("diagnostics() surfaces start provenance on warm eigen results only", {
  n <- 60L; k <- 4L
  A <- build_dense_sym(n)
  warm <- warm_from_truth(A, k)
  fit <- eig_partial(A, k = k, target = largest(), method = block_method(k, n),
                     seed = 1, initial_subspace = warm)
  d <- diagnostics(fit)
  expect_true(all(c("start_source", "initial_subspace") %in% names(d)))
  expect_identical(d$start_source, "user_supplied")
  expect_equal(d$initial_subspace$accepted, k)

  # SVD diagnostics schema must be unchanged (no warm-start fields).
  sv <- svd_partial(matrix(stats::rnorm(60), 10, 6), rank = 3)
  expect_false("start_source" %in% names(diagnostics(sv)))
})
