psd_seeded_spectral_fixture <- function(seed, scale = 1, n = 5L) {
  set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(n * n), n, n)))
  values <- scale * c(16, 4, 1, rep(0, n - 3L))
  list(
    Q = Q,
    values = values,
    K = Q %*% diag(values) %*% t(Q),
    retained = seq_len(3L)
  )
}

psd_relative_frobenius_error <- function(actual, expected) {
  numerator <- norm(as.matrix(actual - expected), "F")
  denominator <- norm(as.matrix(expected), "F")
  if (denominator == 0) numerator else numerator / denominator
}

test_that("seeded dense factors satisfy every canonical action oracle across scales", {
  cases <- expand.grid(
    seed = c(104729L, 130363L, 155921L),
    scale = c(1e-12, 1, 1e12)
  )
  for (case in seq_len(nrow(cases))) {
    fixture <- psd_seeded_spectral_fixture(cases$seed[[case]], cases$scale[[case]])
    factor <- psd_factor(fixture$K)
    Qr <- fixture$Q[, fixture$retained, drop = FALSE]
    lambda <- fixture$values[fixture$retained]
    form <- Qr %*% diag(lambda) %*% t(Qr)
    root <- Qr %*% diag(sqrt(lambda)) %*% t(Qr)
    inverse_root <- Qr %*% diag(1 / sqrt(lambda)) %*% t(Qr)
    pseudoinverse <- Qr %*% diag(1 / lambda) %*% t(Qr)
    projector <- tcrossprod(Qr)
    null_projector <- diag(nrow(Qr)) - projector
    X <- matrix(sin(seq_len(nrow(Qr) * 4L)), nrow(Qr), 4L)

    oracles <- list(
      form = form,
      sqrt = root,
      inverse_sqrt = inverse_root,
      pseudoinverse = pseudoinverse,
      image_projector = projector,
      null_projector = null_projector
    )
    for (action in names(oracles)) {
      actual <- psd_apply(factor, X, action)
      expected <- oracles[[action]] %*% X
      expect_lt(
        psd_relative_frobenius_error(actual, expected),
        5e-12
      )
      columnwise <- do.call(cbind, lapply(seq_len(ncol(X)), function(j) {
        psd_apply(factor, X[, j], action)
      }))
      expect_equal(actual, columnwise, tolerance = 1e-13)
    }

    reduced <- psd_reduce(factor, X)
    expect_lt(
      psd_relative_frobenius_error(crossprod(reduced), crossprod(X, form %*% X)),
      5e-12
    )
    expect_lt(
      psd_relative_frobenius_error(psd_lift(factor, reduced), projector %*% X),
      5e-12
    )
    expect_identical(psd_rank(factor), 3L)
    expect_identical(psd_nullity(factor), as.integer(nrow(Qr) - 3L))
  }
})

test_that("positive scaling has the exact action exponents and no absolute floor", {
  base <- psd_seeded_spectral_fixture(196613L, scale = 1)
  base_factor <- psd_factor(base$K)
  X <- matrix(c(1, 2, 3, 4, 5, -1, 0, 2, 1, 3), 5, 2)
  actions <- list(
    form = function(scale) scale,
    sqrt = function(scale) sqrt(scale),
    inverse_sqrt = function(scale) 1 / sqrt(scale),
    pseudoinverse = function(scale) 1 / scale,
    image_projector = function(scale) 1,
    null_projector = function(scale) 1
  )
  for (scale in c(1e-12, 1e-6, 1e6, 1e12)) {
    scaled <- psd_factor(scale * base$K)
    expect_identical(psd_rank(scaled), psd_rank(base_factor))
    for (action in names(actions)) {
      expected <- actions[[action]](scale) * psd_apply(base_factor, X, action)
      expect_lt(
        psd_relative_frobenius_error(psd_apply(scaled, X, action), expected),
        1e-10
      )
    }
  }

  tiny <- psd_factor(1e-12 * c(4, 1, 1e-9))
  expect_identical(psd_rank(tiny), 2L)
  mutant_floor <- sqrt(.Machine$double.eps) * max(1, tiny$scale)
  expect_gt(mutant_floor, max(psd_spectrum(tiny)))
  expect_lt(tiny$thresholds$rank, min(psd_spectrum(tiny)[1:2]))
})

test_that("threshold neighbors enforce the frozen inclusive inequalities", {
  above_one <- 1 + .Machine$double.eps
  below_one <- 1 - .Machine$double.eps / 2
  expect_gt(above_one, 1)
  expect_lt(below_one, 1)

  positivity <- psd_policy(
    positivity = psd_tolerance(abs = 1, rel = 0),
    rank = psd_tolerance(abs = 0, rel = 0)
  )
  admitted <- psd_factor(c(2, -below_one, -1), policy = positivity)
  expect_identical(
    admitted$spectrum$category,
    c("retained_positive", "accepted_negative", "accepted_negative")
  )
  error <- expect_psd_condition(
    psd_factor(c(2, -above_one), policy = positivity),
    "eigencore_psd_indefinite_error", "indefinite_input", "spectrum"
  )
  expect_lt(error$defect, -1)
  expect_identical(error$threshold, 1)

  rank_policy <- psd_policy(rank = psd_tolerance(abs = 1, rel = 0))
  rank_factor <- psd_factor(c(2, above_one, 1, below_one, 0), policy = rank_policy)
  expect_identical(
    rank_factor$spectrum$category,
    c("retained_positive", "retained_positive", "numerical_null", "numerical_null", "exact_zero")
  )
})

test_that("coordinate permutations and block partitions are action equivariant", {
  fixture <- psd_seeded_spectral_fixture(262147L)
  factor <- psd_factor(fixture$K)
  permutation <- c(5, 2, 4, 1, 3)
  permuted_factor <- psd_factor(fixture$K[permutation, permutation])
  X <- matrix(c(
    1, 2, 3, 4, 5,
    -2, 1, 0, 3, 2,
    4, 0, -1, 2, 1
  ), 5, 3)

  for (action in c(
    "form", "sqrt", "inverse_sqrt", "pseudoinverse",
    "image_projector", "null_projector"
  )) {
    expected <- psd_apply(factor, X, action)[permutation, , drop = FALSE]
    actual <- psd_apply(
      permuted_factor,
      X[permutation, , drop = FALSE],
      action
    )
    expect_equal(actual, expected, tolerance = 1e-11)
    partitioned <- cbind(
      psd_apply(factor, X[, 1, drop = FALSE], action),
      psd_apply(factor, X[, 2:3, drop = FALSE], action)
    )
    expect_equal(partitioned, psd_apply(factor, X, action), tolerance = 1e-13)
  }
})

test_that("rotated repair diagnostics distinguish symmetric and original sources", {
  set.seed(32452843L)
  Q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  surrogate <- Q %*% diag(c(4, 1, -5e-8)) %*% t(Q)
  skew <- matrix(0, 3, 3)
  skew[1, 2] <- 1e-7
  skew[2, 1] <- -1e-7
  source <- surrogate + skew
  repaired <- Q %*% diag(c(4, 1, 0)) %*% t(Q)
  policy <- psd_policy(
    symmetry = psd_tolerance(abs = 5e-7, rel = 0),
    positivity = psd_tolerance(abs = 1e-7, rel = 0)
  )
  factor <- psd_factor(source, policy = policy)

  expect_true(certificate(factor)$repair_applied)
  expect_identical(factor$classification$accepted_negative, 1L)
  expect_lt(
    abs(certificate(factor)$symmetry_defect - norm(source - t(source), "F")),
    1e-15
  )
  expect_equal(
    certificate(factor)$repair_defect,
    norm(surrogate - repaired, "F"),
    tolerance = 1e-8
  )
  expect_equal(
    certificate(factor)$source_action_defect,
    norm(source - repaired, "F"),
    tolerance = 1e-8
  )
  expect_lt(
    psd_relative_frobenius_error(psd_apply(factor, diag(3)), repaired),
    1e-12
  )
})

test_that("mutation sentinels reject wrong root and reduction formulas", {
  fixture <- psd_seeded_spectral_fixture(49979687L)
  factor <- psd_factor(fixture$K)
  X <- matrix(seq_len(10), 5, 2)
  storage.mode(X) <- "double"
  original_apply <- eigencore:::psd_apply_matrix
  testthat::expect_failure(testthat::with_mocked_bindings(
    expect_equal(
      psd_apply(factor, X, "sqrt"),
      fixture$Q[, 1:3] %*% diag(sqrt(fixture$values[1:3])) %*%
        t(fixture$Q[, 1:3]) %*% X,
      tolerance = 1e-11
    ),
    psd_apply_matrix = function(x, X, action) {
      if (identical(action, "sqrt")) {
        original_apply(x, X, "form")
      } else {
        original_apply(x, X, action)
      }
    },
    .package = "eigencore"
  ))

  original_reduce <- eigencore:::psd_reduce_matrix
  testthat::expect_failure(testthat::with_mocked_bindings(
    expect_equal(
      crossprod(psd_reduce(factor, X)),
      psd_gram(factor, X),
      tolerance = 1e-11
    ),
    psd_reduce_matrix = function(x, X) {
      state <- attr(x, "eigencore_psd_state", exact = TRUE)
      crossprod(state$basis, X)
    },
    .package = "eigencore"
  ))
  expect_true(is.function(original_reduce))
})

test_that("mutation sentinels reject Gram roots, zero inversion, and solve conflation", {
  L <- matrix(c(1, 2, 0, 1), 2, 2)
  K <- tcrossprod(L)
  principal <- psd_apply(psd_gram_factor(L), diag(2), "sqrt")
  expect_equal(principal %*% principal, K, tolerance = 1e-12)
  expect_equal(principal, t(principal), tolerance = 1e-14)
  expect_false(isTRUE(all.equal(L %*% L, K, tolerance = 1e-12)))
  expect_false(isTRUE(all.equal(L, t(L), tolerance = 1e-14)))

  singular <- psd_factor(c(4, 0))
  pseudoinverse <- psd_apply(singular, diag(2), "pseudoinverse")
  expect_true(all(is.finite(pseudoinverse)))
  expect_identical(pseudoinverse[2, 2], 0)
  incompatible <- c(8, 1)
  expect_equal(psd_apply(singular, incompatible, "pseudoinverse"), c(2, 0))
  expect_psd_condition(
    psd_solve(singular, incompatible),
    "eigencore_psd_incompatible_rhs", "incompatible_rhs", "B"
  )
})

test_that("evidence, memory, and SPD boundaries cannot be promoted by mutation", {
  sparse_factor <- Matrix::sparseMatrix(
    i = c(1, 2, 3), j = c(1, 1, 2), x = c(1, 1, 1), dims = c(3, 2)
  )
  sparse_factor <- methods::as(sparse_factor, "dgCMatrix")
  structural <- psd_gram_factor(sparse_factor)
  expect_identical(structural$evidence$spectrum_coverage, "structural")
  expect_false(structural$capabilities$numerical_rank$available)
  expect_psd_condition(
    psd_rank(structural),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", "numerical_rank"
  )
  state_bytes <- as.numeric(object.size(attr(structural, "eigencore_psd_state")))
  expect_identical(structural$memory$components$factor_state, state_bytes)
  expect_identical(as.numeric(retained_bytes(structural)), structural$memory$total_bytes)

  expect_error(
    eig_partial(diag(c(3, 2, 1)), B = diag(c(1, 1, 0)), k = 1),
    "dpotrf failed"
  )
  reduced <- psd_reduced_operator(psd_factor(diag(c(1, 1, 0))), diag(c(3, 2, 1)))
  expect_no_error(eig_partial(reduced, k = 1, target = largest()))
})
