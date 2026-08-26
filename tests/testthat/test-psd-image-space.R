psd_image_fixture <- function() {
  Q <- matrix(c(
    1 / sqrt(2), 1 / sqrt(2), 0,
    1 / sqrt(2), -1 / sqrt(2), 0,
    0, 0, 1
  ), 3, 3)
  values <- c(4, 1, 0)
  list(Q = Q, values = values, K = Q %*% diag(values) %*% t(Q))
}

test_that("reduction represents the PSD quotient and discards only null contamination", {
  fixture <- psd_image_fixture()
  factor <- psd_factor(fixture$K)
  X <- matrix(c(1, 2, 3, -2, 1, 4), 3, 2)
  Y <- matrix(c(0, 1, -1, 3, 2, 5), 3, 2)
  null_X <- psd_apply(factor, matrix(c(2, -1, 3, 4, 0, -2), 3, 2), "null_projector")
  null_Y <- psd_apply(factor, matrix(c(-3, 1, 2, 0, 4, 1), 3, 2), "null_projector")

  reduced_X <- psd_reduce(factor, X)
  reduced_Y <- psd_reduce(factor, Y)
  expect_equal(psd_reduce(factor, X + null_X), reduced_X, tolerance = 1e-13)
  expect_equal(psd_reduce(factor, Y + null_Y), reduced_Y, tolerance = 1e-13)
  expect_equal(crossprod(reduced_X, reduced_Y), psd_gram(factor, X, Y), tolerance = 1e-12)
  expect_equal(crossprod(reduced_X), psd_gram(factor, X), tolerance = 1e-12)
  expect_equal(
    psd_lift(factor, reduced_X),
    psd_apply(factor, X, "image_projector"),
    tolerance = 1e-12
  )
  expect_equal(
    psd_lift(factor, psd_reduce(factor, X + null_X)),
    psd_lift(factor, reduced_X),
    tolerance = 1e-12
  )
})

test_that("image and null projectors satisfy complementary projector laws", {
  factor <- psd_factor(psd_image_fixture()$K)
  P <- psd_apply(factor, diag(3), "image_projector")
  N <- psd_apply(factor, diag(3), "null_projector")
  K <- psd_apply(factor, diag(3), "form")

  expect_equal(P, t(P), tolerance = 1e-14)
  expect_equal(N, t(N), tolerance = 1e-14)
  expect_equal(P %*% P, P, tolerance = 1e-13)
  expect_equal(N %*% N, N, tolerance = 1e-13)
  expect_equal(P %*% N, matrix(0, 3, 3), tolerance = 1e-13)
  expect_equal(P + N, diag(3), tolerance = 1e-13)
  expect_equal(K %*% P, K, tolerance = 1e-12)
  expect_equal(K %*% N, matrix(0, 3, 3), tolerance = 1e-12)
})

test_that("strict solve never aliases pseudoinverse application", {
  factor <- psd_factor(c(4, 0))
  incompatible <- c(8, 1e-6)

  expect_equal(psd_apply(factor, incompatible, "pseudoinverse"), c(2, 0))
  error <- expect_psd_condition(
    psd_solve(factor, incompatible),
    "eigencore_psd_incompatible_rhs", "incompatible_rhs", "B"
  )
  expect_identical(error$indices, 1L)
  expect_equal(error$defect, 1e-6)

  admitted <- psd_solve(
    factor,
    incompatible,
    tolerance = psd_tolerance(abs = 1e-6, rel = 0)
  )
  expect_true(admitted$compatible)
  expect_true(certificate(admitted)$passed)
  expect_equal(admitted$solution, c(2, 0))
  expect_equal(certificate(admitted)$residuals$equation, 1e-6)
  expect_equal(certificate(admitted)$thresholds$compatibility, 1e-6)
})

test_that("block orthonormalization reports quotient rank and conditioning", {
  factor <- psd_factor(c(1, 1, 0))
  rank_deficient <- cbind(c(1, 0, 3), c(2, 0, -5))
  collapsed <- psd_orthonormalize(factor, rank_deficient)
  expect_identical(collapsed$rank, 1L)
  expect_length(collapsed$dropped, 1L)
  expect_equal(psd_gram(factor, collapsed$basis), matrix(1), tolerance = 1e-13)

  ill_conditioned <- cbind(c(1, 0, 0), c(1, 1e-3, 7))
  resolved <- psd_orthonormalize(factor, ill_conditioned, required_rank = 2)
  expect_identical(resolved$rank, 2L)
  expect_gt(resolved$condition, 1e6)
  expect_true(certificate(resolved)$passed)
  expect_lte(
    resolved$postcondition_error,
    certificate(resolved)$action_bounds$postcondition
  )

  zero <- psd_orthonormalize(psd_factor(rep(0, 3)), diag(3))
  expect_identical(zero$rank, 0L)
  expect_equal(dim(zero$basis), c(3, 0))
  expect_true(is.na(zero$condition))
  expect_identical(zero$dropped, 1:3)
  expect_true(certificate(zero)$passed)
  expect_psd_condition(
    psd_orthonormalize(psd_factor(rep(0, 3)), diag(3), required_rank = 1),
    "eigencore_psd_infeasible_block", "infeasible_block_rank", "required_rank"
  )
})

test_that("the reduced operator solves the finite image-space eigenproblem", {
  fixture <- psd_image_fixture()
  factor <- psd_factor(fixture$K)
  A <- fixture$Q %*% diag(c(8, 3, 5)) %*% t(fixture$Q)
  reduced <- psd_reduced_operator(factor, A)
  dense_reduced <- reduced$apply(diag(2))
  oracle <- eigen(dense_reduced, symmetric = TRUE)

  largest_fit <- eig_partial(reduced, k = 1, target = largest(), tol = 1e-12)
  smallest_fit <- eig_partial(reduced, k = 1, target = smallest(), tol = 1e-12)
  expect_equal(values(largest_fit), max(oracle$values), tolerance = 1e-10)
  expect_equal(values(smallest_fit), min(oracle$values), tolerance = 1e-10)

  for (fit in list(largest_fit, smallest_fit)) {
    value <- values(fit)[[1L]]
    original_vector <- psd_lift(factor, vectors(fit)[, 1L])
    residual <- A %*% original_vector -
      value * psd_apply(factor, original_vector, "form")
    expect_lt(sqrt(sum(residual^2)), 1e-9)
  }
})

test_that("fresh operation certificates detect backend postcondition mutations", {
  factor <- psd_factor(c(4, 1))
  original_apply <- eigencore:::psd_apply_matrix
  testthat::local_mocked_bindings(
    psd_apply_matrix = function(x, X, action) {
      out <- original_apply(x, X, action)
      if (identical(action, "form")) out + 1e-3 else out
    },
    .package = "eigencore"
  )
  solve_result <- psd_solve(factor, cbind(c(8, 2), c(4, -1)))
  expect_false(certificate(solve_result)$passed)
  expect_identical(certificate(solve_result)$failed_indices, c(1L, 2L))
  expect_true(all(
    certificate(solve_result)$residuals$equation >
      certificate(solve_result)$action_bounds$equation
  ))
})

test_that("block certificates do not borrow a permissive rank threshold", {
  factor <- psd_factor(c(4, 1))
  calls <- 0L
  original_gram <- eigencore:::psd_gram
  testthat::local_mocked_bindings(
    psd_gram = function(x, X, Y = NULL) {
      calls <<- calls + 1L
      out <- original_gram(x, X, Y)
      if (calls > 1L) out + diag(1e-3, nrow(out)) else out
    },
    .package = "eigencore"
  )
  result <- psd_orthonormalize(
    factor,
    diag(2),
    tolerance = psd_tolerance(abs = 0.5, rel = 0)
  )
  expect_false(certificate(result)$passed)
  expect_gt(result$postcondition_error, certificate(result)$action_bounds$postcondition)
  expect_lt(certificate(result)$action_bounds$postcondition, 1e-6)
})
