psd_hadamard4 <- function() {
  matrix(c(
    1, 1, 1, 1,
    1, -1, 1, -1,
    1, 1, -1, -1,
    1, -1, -1, 1
  ), 4, 4, byrow = TRUE) / 2
}

psd_known_fixture <- function(values = c(9, 1, 0, 0)) {
  Q <- psd_hadamard4()
  list(Q = Q, values = values, K = Q %*% diag(values) %*% t(Q))
}

test_that("identity and diagonal paths are analytic and shape-complete", {
  identity <- psd_identity(4)
  X <- matrix(seq_len(12), 4, 3)
  storage.mode(X) <- "double"
  expect_identical(identity$representation, "identity")
  expect_identical(identity$method, "analytic identity PSD factor")
  expect_false(identity$materialization$dense_n_by_n)
  expect_identical(psd_rank(identity), 4L)
  expect_identical(psd_nullity(identity), 0L)
  expect_identical(psd_apply(identity, X, "form"), X)
  expect_identical(psd_apply(identity, X, "sqrt"), X)
  expect_identical(psd_apply(identity, X, "inverse_sqrt"), X)
  expect_identical(psd_apply(identity, X, "pseudoinverse"), X)
  expect_identical(psd_apply(identity, X, "image_projector"), X)
  expect_equal(psd_apply(identity, X, "null_projector"), 0 * X)
  expect_identical(psd_reduce(identity, X), X)
  expect_identical(psd_lift(identity, X), X)

  diagonal <- psd_factor(c(a = 9, b = 4, c = 1, d = 0))
  expect_identical(diagonal$representation, "diagonal")
  expect_false(diagonal$materialization$dense_n_by_n)
  expect_equal(psd_spectrum(diagonal), c(9, 4, 1, 0))
  expect_equal(psd_apply(diagonal, c(1, 2, 3, 4), "form"), c(a = 9, b = 8, c = 3, d = 0))
  expect_equal(psd_apply(diagonal, c(1, 2, 3, 4), "sqrt"), c(a = 3, b = 4, c = 3, d = 0))
  expect_equal(psd_apply(diagonal, c(1, 2, 3, 4), "pseudoinverse"), c(a = 1 / 9, b = 1 / 2, c = 3, d = 0))
  expect_identical(psd_rank(diagonal), 3L)
  expect_identical(psd_nullity(diagonal), 1L)
  expect_identical(psd_nullity(diagonal, "algebraic"), 1L)

  empty <- psd_factor(numeric())
  expect_identical(empty$dim, c(0L, 0L))
  expect_identical(psd_rank(empty), 0L)
  expect_identical(psd_nullity(empty), 0L)
  expect_identical(psd_reduce(empty, numeric()), numeric())
  expect_identical(psd_lift(empty, numeric()), numeric())
  expect_equal(dim(psd_apply(empty, matrix(numeric(), 0, 3), "form")), c(0, 3))

  zero <- psd_factor(rep(0, 5))
  expect_identical(psd_rank(zero), 0L)
  expect_identical(psd_nullity(zero), 5L)
  expect_equal(dim(psd_reduce(zero, matrix(1, 5, 2))), c(0, 2))
  expect_equal(psd_lift(zero, matrix(numeric(), 0, 2)), matrix(0, 5, 2))

  high_rank_threshold <- psd_policy(rank = psd_tolerance(abs = 2, rel = 0))
  truncated_identity <- psd_identity(3, policy = high_rank_threshold)
  expect_identical(psd_rank(truncated_identity), 0L)
  expect_equal(psd_apply(truncated_identity, diag(3), "form"), matrix(0, 3, 3))
  expect_equal(psd_apply(truncated_identity, diag(3), "null_projector"), diag(3))
  expect_true(truncated_identity$certificate$repair_applied)
  expect_equal(truncated_identity$certificate$repair_defect, sqrt(3))
  expect_equal(truncated_identity$certificate$source_action_defect, sqrt(3))
})

test_that("classification boundaries are independent, inclusive, and diagnostic", {
  policy <- psd_policy(
    positivity = psd_tolerance(abs = 1e-8, rel = 0),
    rank = psd_tolerance(abs = 1e-6, rel = 0)
  )
  values <- c(2, 1e-6 + 1e-12, 1e-6, 0, -1e-8)
  factor <- psd_factor(values, policy = policy)
  expect_identical(
    factor$spectrum$category,
    c("retained_positive", "retained_positive", "numerical_null", "exact_zero", "accepted_negative")
  )
  expect_equal(psd_spectrum(factor, repaired = TRUE), c(2, 1e-6 + 1e-12, 0, 0, 0))
  expect_identical(factor$classification$accepted_negative, 1L)
  expect_identical(factor$classification$numerical_null, 1L)
  expect_true(factor$certificate$repair_applied)
  expect_identical(factor$evidence$action_fidelity, "repaired_with_defect")
  expect_equal(factor$certificate$repair_defect, sqrt((1e-6)^2 + (1e-8)^2))

  error <- expect_psd_condition(
    psd_factor(c(2, -1e-8 - 1e-15), policy = policy),
    "eigencore_psd_indefinite_error", "indefinite_input", "spectrum"
  )
  expect_lt(error$defect, -error$threshold)
  expect_s3_class(error$details$certificate, "eigencore_psd_certificate")
  expect_false(error$details$certificate$passed)
  expect_identical(error$details$certificate$classification$materially_negative, 1L)

  reject <- psd_policy(
    positivity = psd_tolerance(abs = 1e-8, rel = 0),
    negative_repair = "reject"
  )
  expect_psd_condition(
    psd_factor(c(2, -1e-8), policy = reject),
    "eigencore_psd_indefinite_error", "indefinite_input", "spectrum"
  )
})

test_that("default classification is invariant from 1e-12 through 1e12", {
  scales <- c(1e-12, 1e-6, 1, 1e6, 1e12)
  factors <- lapply(scales, function(scale) {
    psd_factor(scale * c(4, 1, 1e-9, -1e-15))
  })
  categories <- lapply(factors, function(x) x$spectrum$category)
  for (category in categories) expect_identical(category, categories[[1L]])
  expect_identical(
    categories[[1L]],
    c("retained_positive", "retained_positive", "numerical_null", "accepted_negative")
  )
  expect_equal(
    vapply(factors, psd_rank, integer(1L)),
    rep(2L, length(scales))
  )
  base_projector <- psd_apply(factors[[3L]], diag(4), "image_projector")
  for (factor in factors) {
    expect_equal(psd_apply(factor, diag(4), "image_projector"), base_projector, tolerance = 1e-14)
  }
  expect_equal(
    vapply(factors, function(x) x$thresholds$rank / x$scale, numeric(1L)),
    rep(sqrt(.Machine$double.eps), length(scales))
  )
})

test_that("symmetry validation runs before the upper-triangle spectral kernel", {
  upper_spd_lower_corrupt <- matrix(c(
    2, 0, 0,
    0, 2, 9,
    0, 0, 2
  ), 3, 3)
  expect_false(isSymmetric(upper_spd_lower_corrupt))
  error <- expect_psd_condition(
    psd_factor(upper_spd_lower_corrupt),
    "eigencore_psd_asymmetry_error", "asymmetric_input", "x"
  )
  expect_gt(error$defect, error$threshold)

  small_skew <- diag(c(3, 2, 1))
  small_skew[2, 1] <- 5e-7
  policy <- psd_policy(symmetry = psd_tolerance(abs = 1e-6, rel = 0))
  admitted <- psd_factor(small_skew, policy = policy)
  expected <- (small_skew + t(small_skew)) / 2
  expect_true(admitted$certificate$repair_applied)
  expect_equal(
    admitted$certificate$symmetry_defect,
    sqrt(2) * 5e-7,
    tolerance = 1e-15
  )
  expect_equal(psd_apply(admitted, diag(3), "form"), expected, tolerance = 1e-12)
  expect_equal(
    admitted$certificate$source_action_defect,
    norm(small_skew - expected, "F"),
    tolerance = 1e-12
  )

  reject <- psd_policy(
    symmetry = psd_tolerance(abs = 1e-6, rel = 0),
    symmetry_repair = "reject"
  )
  expect_psd_condition(
    psd_factor(small_skew, policy = reject),
    "eigencore_psd_asymmetry_error", "asymmetric_input", "x"
  )
})

test_that("dense singular actions satisfy independent Moore-Penrose and image laws", {
  fixture <- psd_known_fixture()
  factor <- psd_factor(fixture$K)
  Qr <- fixture$Q[, 1:2, drop = FALSE]
  lambda <- fixture$values[1:2]
  K <- fixture$K
  P <- tcrossprod(Qr)
  Khalf <- Qr %*% diag(sqrt(lambda)) %*% t(Qr)
  Kihalf <- Qr %*% diag(1 / sqrt(lambda)) %*% t(Qr)
  Kplus <- Qr %*% diag(1 / lambda) %*% t(Qr)
  X <- matrix(c(1, 2, 3, 4, -1, 0, 2, 1), 4, 2)

  expect_identical(psd_rank(factor), 2L)
  expect_identical(psd_nullity(factor), 2L)
  expect_psd_condition(
    psd_nullity(factor, "algebraic"),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence",
    "algebraic_nullity"
  )
  expect_equal(psd_apply(factor, X, "form"), K %*% X, tolerance = 1e-12)
  expect_equal(psd_apply(factor, X, "sqrt"), Khalf %*% X, tolerance = 1e-12)
  expect_equal(psd_apply(factor, X, "inverse_sqrt"), Kihalf %*% X, tolerance = 1e-12)
  expect_equal(psd_apply(factor, X, "pseudoinverse"), Kplus %*% X, tolerance = 1e-12)
  expect_equal(psd_apply(factor, X, "image_projector"), P %*% X, tolerance = 1e-12)
  expect_equal(psd_apply(factor, X, "null_projector"), (diag(4) - P) %*% X, tolerance = 1e-12)

  root <- psd_apply(factor, diag(4), "sqrt")
  invroot <- psd_apply(factor, diag(4), "inverse_sqrt")
  pinv <- psd_apply(factor, diag(4), "pseudoinverse")
  projector <- psd_apply(factor, diag(4), "image_projector")
  expect_equal(root %*% root, psd_apply(factor, diag(4), "form"), tolerance = 1e-11)
  expect_equal(root %*% invroot, projector, tolerance = 1e-11)
  expect_equal(K %*% pinv, projector, tolerance = 1e-11)
  expect_equal(pinv %*% K, projector, tolerance = 1e-11)
  expect_equal(projector %*% projector, projector, tolerance = 1e-11)
  expect_equal(projector, t(projector), tolerance = 1e-14)

  reduced <- psd_reduce(factor, X)
  lifted <- psd_lift(factor, reduced)
  expect_equal(crossprod(reduced), crossprod(X, K %*% X), tolerance = 1e-11)
  expect_equal(lifted, P %*% X, tolerance = 1e-11)
  expect_equal(psd_reduce(factor, psd_lift(factor, diag(2))), diag(2), tolerance = 1e-11)
})

test_that("repeated eigenspaces are certified through projectors, not basis columns", {
  Q <- psd_hadamard4()
  values <- c(9, 9, 0, 0)
  K <- Q %*% diag(values) %*% t(Q)
  theta <- 0.37
  rotation <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  Q2 <- Q
  Q2[, 1:2] <- Q[, 1:2] %*% rotation
  Q2[, 3:4] <- Q[, 3:4] %*% rotation
  K2 <- Q2 %*% diag(values) %*% t(Q2)
  expect_equal(K2, K, tolerance = 1e-14)

  factor1 <- psd_factor(K)
  factor2 <- psd_factor(K2)
  P1 <- psd_apply(factor1, diag(4), "image_projector")
  P2 <- psd_apply(factor2, diag(4), "image_projector")
  expect_equal(P1, P2, tolerance = 1e-12)
  expect_equal(psd_apply(factor1, diag(4), "form"), psd_apply(factor2, diag(4), "form"), tolerance = 1e-12)
})

test_that("pseudoinverse application is total and strict solve is image-compatible", {
  fixture <- psd_known_fixture()
  factor <- psd_factor(fixture$K)
  P <- psd_apply(factor, diag(4), "image_projector")
  B <- cbind(P %*% c(1, 2, 3, 4), P %*% c(-1, 0, 2, 1))
  fit <- psd_solve(factor, B)
  expect_s3_class(fit, "eigencore_psd_solve_result")
  expect_true(all(fit$compatible))
  expect_true(certificate(fit)$passed)
  expect_equal(certificate(fit)$scale, factor$scale)
  expect_named(certificate(fit)$thresholds, c("compatibility", "equation"))
  expect_named(certificate(fit)$action_bounds, c("compatibility", "equation"))
  expect_true(all(
    certificate(fit)$residuals$equation <=
      certificate(fit)$action_bounds$equation
  ))
  expect_equal(psd_apply(factor, fit$solution, "form"), B, tolerance = 1e-11)
  expect_equal(fit$solution, psd_apply(factor, B, "pseudoinverse"), tolerance = 1e-12)

  bad <- B
  bad[, 2] <- bad[, 2] + psd_apply(factor, c(1, 2, 3, 4), "null_projector")
  error <- expect_psd_condition(
    psd_solve(factor, bad),
    "eigencore_psd_incompatible_rhs", "incompatible_rhs", "B"
  )
  expect_identical(error$indices, 2L)
  expect_gt(error$defect[[2L]], error$threshold[[2L]])
  expect_no_error(psd_apply(factor, bad, "pseudoinverse"))

  expect_psd_condition(
    psd_solve(factor, B, tolerance = 1e-8),
    "eigencore_psd_invalid_input", "invalid_policy", "tolerance"
  )
})

test_that("PSD Gram, block orthonormalization, and reduced operator are diagnostic", {
  factor <- psd_factor(c(9, 4, 0))
  X <- matrix(c(1, 0, 1, 0, 1, 1), 3, 2)
  Y <- matrix(c(1, 2, 3, -1, 0, 2), 3, 2)
  K <- diag(c(9, 4, 0))
  expect_equal(psd_gram(factor, X), crossprod(X, K %*% X))
  expect_equal(psd_gram(factor, X, Y), crossprod(X, K %*% Y))

  block <- psd_orthonormalize(factor, diag(3))
  expect_s3_class(block, "eigencore_psd_block_result")
  expect_identical(block$rank, 2L)
  expect_true(certificate(block)$passed)
  expect_equal(certificate(block)$scale, factor$scale)
  expect_named(certificate(block)$thresholds, c("gram_rank", "postcondition"))
  expect_equal(
    certificate(block)$action_bounds$postcondition,
    sqrt(.Machine$double.eps) * sqrt(2)
  )
  expect_equal(psd_gram(factor, block$basis), diag(2), tolerance = 1e-12)
  expect_psd_condition(
    psd_orthonormalize(factor, diag(3), required_rank = 3),
    "eigencore_psd_infeasible_block", "infeasible_block_rank", "required_rank"
  )

  A <- diag(c(2, 3, 4))
  reduced <- psd_reduced_operator(factor, A)
  expected <- diag(c(1 / 3, 1 / 2)) %*% A[1:2, 1:2] %*% diag(c(1 / 3, 1 / 2))
  expect_equal(reduced$apply(diag(2)), expected, tolerance = 1e-14)
  expect_identical(reduced$dim, c(2L, 2L))
})

test_that("snapshots serialize, identities are deterministic, and mutation fails closed", {
  factor1 <- psd_factor(c(9, 4, 0))
  factor2 <- psd_factor(c(9, 4, 0))
  changed <- psd_factor(c(9, 4, 1e-3))
  expect_identical(operator_identity(factor1), operator_identity(factor2))
  expect_false(identical(operator_identity(factor1), operator_identity(changed)))

  path <- tempfile(fileext = ".rds")
  saveRDS(factor1, path)
  restored <- readRDS(path)
  expect_identical(psd_apply(restored, c(1, 2, 3), "form"), c(9, 8, 0))
  expect_identical(certificate(restored), certificate(factor1))

  operator <- as_operator(factor1)
  expect_identical(operator_identity(operator), operator_identity(factor1))
  mutated_caller <- factor1
  mutated_caller$rank <- 99L
  expect_psd_condition(
    psd_rank(mutated_caller),
    "eigencore_psd_corrupt_state", "corrupt_factor", "integrity_token"
  )
  expect_identical(operator$apply(matrix(c(1, 2, 3), 3, 1)), matrix(c(9, 8, 0), 3, 1))

  mutated_state <- factor1
  state <- attr(mutated_state, "eigencore_psd_state")
  state$repaired[[1L]] <- 100
  attr(mutated_state, "eigencore_psd_state") <- state
  expect_psd_condition(
    psd_apply(mutated_state, c(1, 2, 3)),
    "eigencore_psd_corrupt_state", "corrupt_factor", "factor_state"
  )
})

test_that("work, retained memory, and analytic storage are honest", {
  factor <- psd_factor(c(9, 4, 0))
  expect_s3_class(work(factor), "eigencore_work")
  expect_true(work(factor)$complete)
  expect_identical(work(factor)$operator_block_calls, 0L)
  expect_identical(work(factor)$iterations, 0L)
  expect_gte(work(factor)$setup_seconds, 0)
  expect_s3_class(factor$memory, "eigencore_memory")
  expect_true(factor$memory$complete)
  expect_identical(factor$memory$native_bytes, 0)
  expect_named(
    factor$memory$components,
    c("source_snapshot", "factor_state", "cached_actions", "metadata")
  )
  expect_identical(as.numeric(retained_bytes(factor)), factor$memory$total_bytes)

  small <- psd_identity(1000)
  large <- psd_identity(2000)
  expect_lt(retained_bytes(large) / retained_bytes(small), 3)
  expect_false(any(vapply(
    attr(large, "eigencore_psd_state"),
    function(value) is.matrix(value) && all(dim(value) == c(2000, 2000)),
    logical(1L)
  )))
})

test_that("type, dimension, dimname, Matrix, signed-zero, and unsupported boundaries fail closed", {
  expect_psd_condition(
    psd_factor(c(1L, 0L)),
    "eigencore_psd_invalid_input", "invalid_dtype", "x"
  )
  expect_psd_condition(
    psd_factor(as.complex(diag(2))),
    "eigencore_psd_invalid_input", "invalid_dtype", "x"
  )
  expect_psd_condition(
    psd_factor(matrix(1, 2, 3)),
    "eigencore_psd_invalid_input", "invalid_dimension", "x"
  )
  nonfinite <- diag(2)
  nonfinite[1, 1] <- Inf
  expect_psd_condition(
    psd_factor(nonfinite),
    "eigencore_psd_invalid_input", "nonfinite_input", "x"
  )
  named <- diag(2)
  dimnames(named) <- list(c("a", "b"), c("a", "c"))
  expect_psd_condition(
    psd_factor(named),
    "eigencore_psd_invalid_input", "invalid_dimension", "dimnames"
  )

  diagonal <- Matrix::Diagonal(x = c(4, 1, 0))
  diagonal_factor <- psd_factor(diagonal)
  expect_identical(diagonal_factor$representation, "diagonal")
  expect_false(diagonal_factor$materialization$dense_n_by_n)

  dense_symmetric <- Matrix::forceSymmetric(Matrix::Matrix(matrix(c(
    4, 0.2, 0,
    0.2, 1, 0,
    0, 0, 0
  ), 3, 3), sparse = FALSE))
  dense_factor <- psd_factor(dense_symmetric)
  expect_identical(dense_factor$representation, "dense_spectral")

  sparse <- methods::as(Matrix::Diagonal(x = c(4, 1, 0)), "CsparseMatrix")
  expect_psd_condition(
    psd_factor(sparse),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", "x"
  )

  signed <- psd_factor(c(-0, 0, 1))
  expect_identical(signed$spectrum$signed_zero_count, 1L)
  expect_equal(psd_spectrum(signed, repaired = TRUE), c(1, 0, 0))

  factor <- psd_factor(c(a = 4, b = 1, c = 0))
  wrong_names <- c(a = 1, c = 2, b = 3)
  expect_psd_condition(
    psd_apply(factor, wrong_names),
    "eigencore_psd_invalid_input", "invalid_dimension", "X$rownames"
  )
  expect_psd_condition(
    psd_apply(factor, c(1L, 2L, 3L)),
    "eigencore_psd_invalid_input", "invalid_dtype", "X"
  )
  action_error <- expect_psd_condition(
    psd_apply(factor, c(1, Inf, 3)),
    "eigencore_psd_invalid_input", "nonfinite_input", "X"
  )
  expect_identical(action_error$source_identity, factor$source_identity)
  expect_identical(action_error$factor_identity, factor$operator_identity)
  expect_identical(action_error$representation, factor$representation)

  expect_psd_condition(
    psd_gram_factor(matrix(1, 2, 1)),
    "eigencore_psd_unsupported_action", "unsupported_action", "x"
  )
  lap <- Matrix::sparseMatrix(i = c(1, 1, 2, 2), j = c(1, 2, 1, 2), x = c(1, -1, -1, 1))
  expect_psd_condition(
    psd_laplacian(lap),
    "eigencore_psd_unsupported_action", "unsupported_action", "x"
  )
})

test_that("zero-column blocks are conformable and do no invented work", {
  factor <- psd_factor(c(9, 4, 0))
  X <- matrix(numeric(), 3, 0)
  for (action in c(
    "form", "sqrt", "inverse_sqrt", "pseudoinverse",
    "image_projector", "null_projector"
  )) {
    expect_equal(dim(psd_apply(factor, X, action)), c(3, 0))
  }
  expect_equal(dim(psd_reduce(factor, X)), c(2, 0))
  expect_equal(dim(psd_lift(factor, matrix(numeric(), 2, 0))), c(3, 0))
  expect_equal(dim(psd_gram(factor, X)), c(0, 0))
})

test_that("a nonorthogonal dense kernel result cannot become a certified factor", {
  testthat::local_mocked_bindings(
    native_dense_symmetric_eigen = function(A) {
      list(
        values = c(2, 1),
        vectors = matrix(c(1, 0, 1, 1), 2, 2)
      )
    },
    .package = "eigencore"
  )
  error <- expect_psd_condition(
    psd_factor(diag(c(2, 1))),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", "certificate"
  )
  expect_false(error$details$certificate$passed)
  expect_false(error$details$certificate$orthogonality_passed)
  expect_gt(
    error$details$certificate$max_orthogonality_loss,
    error$details$certificate$orthogonality_tolerance
  )
})
