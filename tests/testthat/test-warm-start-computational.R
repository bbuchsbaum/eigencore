# Computational warm-start invariants.
#
# These tests complement the boundary/regression coverage in test-warm-start.R.
# They deliberately test mathematical relationships rather than snapshots:
# coordinate changes inside a supplied subspace, simultaneous permutation of
# an operator and its start, changed-operator residuals, clustered targets, and
# independently observed matrix-free work.

ws_comp_symmetric <- function(spectrum, seed = 1L) {
  n <- length(spectrum)
  set.seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n, n)))
  list(
    A = Q %*% (spectrum * t(Q)),
    vectors = Q,
    spectrum = spectrum
  )
}

ws_comp_method <- function(k, n) {
  lanczos(
    block = k,
    max_subspace = min(n - 1L, 8L * k),
    max_restarts = 300L
  )
}

ws_comp_projector_distance <- function(X, Y) {
  X <- qr.Q(qr(as.matrix(X)))
  Y <- qr.Q(qr(as.matrix(Y)))
  sqrt(sum((tcrossprod(X) - tcrossprod(Y))^2))
}

ws_comp_target_order <- function(values, name) {
  switch(
    name,
    largest = order(values, decreasing = TRUE),
    smallest = order(values, decreasing = FALSE),
    largest_magnitude = order(abs(values), decreasing = TRUE),
    smallest_magnitude = order(abs(values), decreasing = FALSE),
    stop("unknown target")
  )
}

test_that("warm starts are invariant to coordinates within the supplied subspace", {
  n <- 56L
  k <- 4L
  fixture <- ws_comp_symmetric(seq(9, -3, length.out = n), seed = 101L)
  oracle <- eigen(fixture$A, symmetric = TRUE)

  set.seed(102L)
  start <- oracle$vectors[, seq_len(k), drop = FALSE] +
    matrix(stats::rnorm(n * k, sd = 1e-3), n, k)
  set.seed(103L)
  rotation <- qr.Q(qr(matrix(stats::rnorm(k * k), k, k)))
  mixing <- rotation %*% diag(c(1e-4, 1e-2, 1e2, 1e4))
  mixed_start <- start %*% mixing

  fit <- eig_partial(
    fixture$A, k = k, target = largest(),
    method = ws_comp_method(k, n), tol = 1e-8, seed = 104L,
    initial_subspace = start
  )
  mixed_fit <- eig_partial(
    fixture$A, k = k, target = largest(),
    method = ws_comp_method(k, n), tol = 1e-8, seed = 104L,
    initial_subspace = mixed_start
  )

  truth <- oracle$values[seq_len(k)]
  expect_certificate_clean(fit)
  expect_certificate_clean(mixed_fit)
  expect_equal(sort(values(fit)), sort(truth), tolerance = 1e-7)
  expect_equal(sort(values(mixed_fit)), sort(truth), tolerance = 1e-7)
  expect_lt(ws_comp_projector_distance(vectors(fit), vectors(mixed_fit)), 1e-6)
  expect_equal(fit$initial_subspace$accepted, k)
  expect_equal(mixed_fit$initial_subspace$accepted, k)
  expect_equal(mixed_fit$initial_subspace$rejected, 0L)
})

test_that("changed-operator continuation is permutation equivariant and freshly certified", {
  n <- 60L
  k <- 4L
  fixture <- ws_comp_symmetric(seq(7, -4, length.out = n), seed = 111L)
  grid <- seq_len(n) / (n + 1)
  perturbation <- diag(0.5 + grid + 0.2 * sin(4 * pi * grid))
  changed <- fixture$A - 0.07 * perturbation
  method <- ws_comp_method(k, n)

  previous <- eig_partial(
    fixture$A, k = k, target = smallest(), method = method,
    tol = 1e-8, seed = 112L
  )
  cold <- eig_partial(
    changed, k = k, target = smallest(), method = method,
    tol = 1e-8, seed = 112L
  )
  warm <- eig_partial(
    changed, k = k, target = smallest(), method = method,
    tol = 1e-8, seed = 112L, initial_subspace = vectors(previous)
  )

  permutation <- c(seq(2L, n, by = 2L), seq(1L, n, by = 2L))
  permuted <- eig_partial(
    changed[permutation, permutation], k = k, target = smallest(),
    method = method, tol = 1e-8, seed = 112L,
    initial_subspace = vectors(previous)[permutation, , drop = FALSE]
  )
  unpermuted_vectors <- vectors(permuted)[order(permutation), , drop = FALSE]

  oracle <- sort(eigen(changed, symmetric = TRUE, only.values = TRUE)$values)
  truth <- oracle[seq_len(k)]
  expect_certificate_clean(previous)
  expect_certificate_clean(cold)
  expect_certificate_clean(warm)
  expect_certificate_clean(permuted)
  expect_identical(warm$start_source, "user_supplied")
  expect_identical(permuted$start_source, "user_supplied")
  expect_equal(sort(values(warm)), truth, tolerance = 1e-7)
  expect_equal(sort(values(cold)), truth, tolerance = 1e-7)
  expect_equal(sort(values(permuted)), truth, tolerance = 1e-7)
  expect_lt(
    ws_comp_projector_distance(vectors(warm), unpermuted_vectors),
    1e-6
  )

  manual_residuals <- sqrt(colSums(
    (changed %*% vectors(warm) -
       sweep(vectors(warm), 2L, values(warm), `*`))^2
  ))
  expect_lt(
    max(abs(certificate(warm)$residuals - manual_residuals)),
    1e-12
  )
})

test_that("warm starts recover clustered spectra across native target families", {
  spectrum <- c(
    13, 12.9999, 12.9998,
    seq(8, 2, length.out = 18L),
    0.004, 0.002, -0.001, -0.003,
    seq(-2, -8, length.out = 20L),
    -10.9998, -10.9999, -11
  )
  n <- length(spectrum)
  k <- 3L
  fixture <- ws_comp_symmetric(spectrum, seed = 121L)
  oracle <- eigen(fixture$A, symmetric = TRUE)
  cases <- list(
    largest = largest(),
    smallest = smallest(),
    largest_magnitude = largest_magnitude(),
    smallest_magnitude = smallest_magnitude()
  )

  for (case_name in names(cases)) {
    indices <- ws_comp_target_order(oracle$values, case_name)[seq_len(k)]
    set.seed(122L)
    start <- oracle$vectors[, indices, drop = FALSE] +
      matrix(stats::rnorm(n * k, sd = 1e-4), n, k)
    fit <- eig_partial(
      fixture$A, k = k, target = cases[[case_name]],
      method = ws_comp_method(k, n), tol = 1e-8, seed = 123L,
      initial_subspace = start
    )

    expect_certificate_clean(fit)
    expect_true(all(is.finite(values(fit))))
    expect_equal(
      sort(values(fit)),
      sort(oracle$values[indices]),
      tolerance = 1e-7,
      info = case_name
    )
    expect_lt(
      ws_comp_projector_distance(
        vectors(fit),
        oracle$vectors[, indices, drop = FALSE]
      ),
      1e-5
    )
  }
})

test_that("matrix-free work accounting matches an independently observed callback", {
  n <- 50L
  k <- 3L
  fixture <- ws_comp_symmetric(seq(8, 1, length.out = n), seed = 131L)
  observed <- new.env(parent = emptyenv())
  observed$calls <- 0L
  observed$columns <- 0L

  op <- linear_operator(
    dim = c(n, n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      observed$calls <- observed$calls + 1L
      observed$columns <- observed$columns +
        if (is.null(dim(X))) 1L else ncol(X)
      Z <- alpha * (fixture$A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian(),
    metadata = list(frobenius_norm = sqrt(sum(fixture$A^2)))
  )
  set.seed(132L)
  start <- fixture$vectors[, seq_len(k), drop = FALSE] +
    matrix(stats::rnorm(n * k, sd = 1e-3), n, k)
  fit <- eig_partial(
    op, k = k, target = largest(), method = lanczos(),
    maxit = n, tol = 1e-8, seed = 133L, initial_subspace = start
  )

  expect_certificate_clean(fit)
  expect_equal(
    sort(values(fit), decreasing = TRUE),
    fixture$spectrum[seq_len(k)],
    tolerance = 1e-7
  )
  expect_identical(fit$operator_block_calls, observed$calls)
  expect_identical(fit$operator_columns, observed$columns)
  expect_equal(
    fit$operator_columns,
    fit$matvecs +
      fit$certification_operator_columns +
      fit$initial_subspace$guard_operator_columns
  )
  expect_gt(fit$certification_operator_columns, 0L)
  expect_gt(fit$initial_subspace$guard_operator_columns, 0L)
})

test_that("rank-deficient augmentation is reproducible and preserves accepted directions", {
  n <- 40L
  width <- 4L
  set.seed(141L)
  basis <- qr.Q(qr(matrix(stats::rnorm(n * 2L), n, 2L)))
  supplied <- cbind(
    basis[, 1L],
    1e100 * basis[, 1L],
    basis[, 2L],
    1e-100 * basis[, 2L]
  )

  set.seed(142L)
  first <- eigencore:::prepare_initial_subspace(
    supplied, n = n, width = width
  )
  set.seed(142L)
  repeated <- eigencore:::prepare_initial_subspace(
    supplied, n = n, width = width
  )
  set.seed(143L)
  other_seed <- eigencore:::prepare_initial_subspace(
    supplied, n = n, width = width
  )

  expect_identical(first$start, repeated$start)
  expect_equal(first$rank, 2L)
  expect_equal(first$accepted, 2L)
  expect_equal(first$rejected, 2L)
  expect_equal(first$augmented, 2L)
  expect_equal(crossprod(first$start), diag(width), tolerance = 1e-12)
  expect_lt(
    sqrt(sum(
      (basis - first$start %*% crossprod(first$start, basis))^2
    )),
    1e-12
  )
  expect_false(isTRUE(all.equal(first$start, other_seed$start)))
  expect_equal(
    tcrossprod(first$accepted_basis),
    tcrossprod(other_seed$accepted_basis),
    tolerance = 1e-12
  )
})
