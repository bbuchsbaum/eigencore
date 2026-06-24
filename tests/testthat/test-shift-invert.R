test_that("shift-invert returns interior eigenvalues nearest sigma on dense Hermitian", {
  set.seed(42)
  n <- 20L
  vals <- c(-3, -1, 0.5, 1.5, 2.7, 4, seq(5, 30, length.out = n - 6L))
  A    <- symmetric_with_spectrum(vals, seed = 42)

  fit <- eig_partial(A, k = 3L, target = nearest(2),
                     method = shift_invert(sigma = 2))

  expect_identical(fit$method,
                   eigencore:::native_dense_shift_invert_label())
  expected <- vals[order(abs(vals - 2))][1:3]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  direct <- eigencore:::certify_eigen_operator(as_operator(A), fit$values,
                                               fit$vectors, tol = 1e-8)
  expect_equal(fit$certificate$residuals, direct$residuals, tolerance = 1e-12)
  expect_equal(fit$certificate$backward_error, direct$backward_error,
               tolerance = 1e-12)
  expect_identical(fit$transform$kind, "shift_invert")
  expect_identical(fit$transform$label_kind, "dense_lu_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_identical(fit$transform$factorization_cache$factorization,
                   "LAPACK dgetrf/dgetrs")
  expect_identical(fit$restart$kind, "native_dense_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_identical(fit$transform$certification$problem, "original")
  expect_false(fit$transform$certification$transformed_residuals_used)
  expect_true(fit$plan$controls$certified_in_original_coordinates)
  expect_identical(fit$plan$controls$transformed_operator_target,
                   "largest_magnitude")
  expect_equal(fit$sigma, 2)
})

test_that("auto nearest target routes through shift-invert with original-coordinate certification", {
  set.seed(420)
  vals <- c(-3, -1, 0.5, 1.5, 2.7, 4, seq(5, 30, length.out = 14L))
  A <- symmetric_with_spectrum(vals, seed = 420)

  plan <- plan_solver(eigen_problem(A, target = nearest(2)), k = 3L)
  expect_identical(plan$method, eigencore:::native_dense_shift_invert_label())
  expect_true(any(grepl("nearest target routed through shift_invert", plan$reasons,
                        fixed = TRUE)))
  expect_identical(plan$controls$transform, "shift_invert")
  expect_identical(plan$controls$transformed_operator_target,
                   "largest_magnitude")
  expect_true(plan$controls$certified_in_original_coordinates)

  fit <- eig_partial(A, k = 3L, target = nearest(2))
  expected <- vals[order(abs(vals - 2))][1:3]
  expect_identical(fit$method, eigencore:::native_dense_shift_invert_label())
  expect_equal(sort(values(fit)), sort(expected), tolerance = 1e-7)
  expect_identical(fit$transform$kind, "shift_invert")
  expect_identical(fit$target, "nearest(2)")
  expect_certificate_clean(fit)

  P <- eigen_problem(A, target = nearest(2))
  fit_direct <- solve(P, k = 3L)
  expect_identical(fit_direct$method, eigencore:::native_dense_shift_invert_label())
  expect_equal(sort(values(fit_direct)), sort(expected), tolerance = 1e-7)
  expect_identical(fit_direct$transform$kind, "shift_invert")
})

test_that("auto nearest target preserves sparse shift-invert boundary labels", {
  set.seed(421)
  vals <- seq(1, 30)
  A_csc <- methods::as(
    Matrix::Matrix(symmetric_with_spectrum(vals, seed = 421), sparse = TRUE),
    "CsparseMatrix"
  )

  fit <- eig_partial(A_csc, k = 4L, target = nearest(15.5))
  expected <- vals[order(abs(vals - 15.5))][1:4]
  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (sparse LU)")
  expect_equal(sort(values(fit)), sort(expected), tolerance = 1e-7)
  expect_identical(fit$transform$kind, "shift_invert")
  expect_equal(fit$transform$factorization_cache$contract$provider,
               "Matrix::lu_reference_factorization")
  expect_equal(fit$transform$factorization_cache$contract$promotion_status,
               "reference_boundary")
  expect_true(all(fit$certificate$converged))
  # The dsCMatrix source carries an exact Frobenius norm (Matrix::norm), so the
  # sparse-LU shift-invert certificate uses an exact original-coordinate scale
  # and passed is not withheld. See review item bd-01KVWRKQ2JMJJJ0CKM939C4NZ6.
  expect_false(fit$certificate$scale_is_estimate)
  expect_true(fit$certificate$passed)
})

test_that("auto nearest target fails loudly for matrix-free operators without a solve", {
  A <- diag(c(1, 3, 6, 10))
  op <- linear_operator(
    dim = c(4L, 4L),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- A %*% X
      if (is.null(Y)) {
        alpha * out
      } else {
        alpha * out + beta * Y
      }
    },
    structure = hermitian(),
    name = "matrix_free_hermitian"
  )

  plan <- plan_solver(eigen_problem(op, target = nearest(3)), k = 1L)
  expect_match(plan$method, "provide method\\$solve for matrix-free A")
  expect_true(any(grepl("nearest target routed through shift_invert", plan$reasons,
                        fixed = TRUE)))

  expect_error(
    eig_partial(op, k = 1L, target = nearest(3)),
    "user-supplied solve operator"
  )
})

test_that("shift-invert handles a sparse CSC source via factorized solve", {
  set.seed(7)
  vals <- seq(1, 30)
  A_csc <- methods::as(
    Matrix::Matrix(symmetric_with_spectrum(vals, seed = 7), sparse = TRUE),
    "CsparseMatrix"
  )

  fit <- eig_partial(A_csc, k = 4L, target = nearest(15.5),
                     method = shift_invert(sigma = 15.5))

  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (sparse LU)")
  expected <- vals[order(abs(vals - 15.5))][1:4]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  # The dsCMatrix source carries an exact Frobenius norm (Matrix::norm), so the
  # sparse-LU shift-invert certificate uses an exact original-coordinate scale;
  # all pairs converge and passed is no longer withheld (was previously withheld
  # when this path fell back to a Hutchinson norm estimate).
  # See review item bd-01KVWRKQ2JMJJJ0CKM939C4NZ6.
  expect_true(all(fit$certificate$converged))
  expect_false(fit$certificate$scale_is_estimate)
  expect_true(fit$certificate$passed)
  expect_lt(max(fit$certificate$backward_error), 1e-7)
  cache <- fit$transform$factorization_cache
  expect_equal(cache$label_kind, "sparse_lu")
  expect_equal(cache$factorization, "Matrix::lu")
  expect_true(cache$factorization_cached)
  expect_equal(cache$condition_estimate_type, "sparse_lu_pivot_ratio")
  expect_true(is.finite(cache$condition_estimate))
  expect_gt(cache$condition_estimate, 0)
  expect_false(isTRUE(cache$near_singular))
  contract <- cache$contract
  expect_equal(contract$contract_version, "shift_invert_factorization_contract_v1")
  expect_equal(contract$provider, "Matrix::lu_reference_factorization")
  expect_equal(contract$promotion_status, "reference_boundary")
  expect_false(contract$owned_by_eigencore)
  expect_false(contract$external_cache)
  expect_equal(contract$memory_policy, "sparse_factorization_no_dense_rcond")
  expect_equal(contract$certificate_policy, "original_coordinate_residual_required")
})

test_that("shift-invert uses native tridiagonal factorized Lanczos for sparse CSC paths", {
  n <- 50L
  A_csc <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
  sigma <- 0.01

  fit <- eig_partial(A_csc, k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-8,
                     allow_dense_fallback = "never")
  expected <- eigen(as.matrix(A_csc), symmetric = TRUE, only.values = TRUE)$values
  expected <- expected[order(abs(expected - sigma))][1:3]

  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$label_kind, "tridiagonal_thomas_native")
  cache <- fit$transform$factorization_cache
  expect_true(cache$native)
  expect_equal(cache$factorization, "native tridiagonal Thomas")
  expect_equal(cache$condition_estimate_type, "tridiagonal_thomas_pivot_ratio")
  expect_true(is.finite(cache$condition_estimate))
  expect_gt(cache$condition_estimate, 0)
  expect_identical(fit$restart$kind, "native_tridiagonal_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$factorization_native)
  expect_identical(fit$transform$certification$problem, "original")
  expect_false(fit$transform$certification$transformed_residuals_used)
})

test_that("shift-invert uses native tridiagonal factorized Lanczos for diagonal sources", {
  A <- Matrix::Diagonal(x = c(1, 2, 4, 8, 16, 32))
  sigma <- 7

  fit <- eig_partial(A, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-10,
                     allow_dense_fallback = "never")

  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_equal(sort(fit$values), c(4, 8), tolerance = 1e-10)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$label_kind, "tridiagonal_thomas_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_identical(fit$restart$kind, "native_tridiagonal_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$factorization_native)
})

test_that("shift-invert accepts a user-supplied solve operator", {
  set.seed(11)
  n <- 15L
  vals <- seq(0.1, 14.5, length.out = n)
  A    <- symmetric_with_spectrum(vals, seed = 11)
  sigma <- 7.4   # not an eigenvalue
  M     <- A - sigma * diag(n)
  user_solve <- function(X) base::solve(M, X)

  fit <- eig_partial(A, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma, solve = user_solve))

  expect_identical(fit$method,
                   "reference Hermitian Lanczos shift-invert (user solve)")
  expected <- vals[order(abs(vals - sigma))][1:2]
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_equal(fit$transform$factorization_cache$factorization, "user_solve")
  expect_true(fit$transform$factorization_cache$external_cache)
  contract <- fit$transform$factorization_cache$contract
  expect_equal(contract$contract_version, "shift_invert_factorization_contract_v1")
  expect_equal(contract$provider, "user_supplied_solve")
  expect_equal(contract$promotion_status, "reference_boundary")
  expect_false(contract$owned_by_eigencore)
  expect_true(contract$external_cache)
  expect_equal(contract$memory_policy, "external_cache_user_owned_no_dense_fallback")
  expect_equal(contract$certificate_policy, "original_coordinate_residual_required")
})

test_that("shift-invert refuses ignored public factorization handles", {
  A <- diag(c(1, 2, 4, 8))
  expect_error(
    eig_partial(A, k = 1L, target = nearest(2.5),
                method = shift_invert(sigma = 2.5, factorization = list())),
    "factorization = .*not implemented"
  )
})

test_that("generalized shift-invert certifies original generalized residuals on dense SPD problems", {
  set.seed(9)
  n <- 10L
  A <- symmetric_with_spectrum(seq(1, n), seed = 9)
  B <- diag(seq(1, 2, length.out = n))
  sigma <- 3.2

  fit <- eig_partial(A, B = B, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-9)
  oracle <- eigencore:::dense_generalized_spd_eigen(A, B)
  expected <- oracle$values[order(abs(oracle$values - sigma))][1:2]

  expect_identical(
    fit$method,
    eigencore:::native_dense_generalized_shift_invert_label()
  )
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-7)
  expect_certificate_clean(fit)
  direct <- eigencore:::certify_eigen_operator(
    as_operator(A), fit$values, fit$vectors, Bop = as_operator(B), tol = 1e-9
  )
  expect_equal(fit$certificate$residuals, direct$residuals, tolerance = 1e-12)
  expect_equal(fit$certificate$backward_error, direct$backward_error,
               tolerance = 1e-12)
  expect_equal(crossprod(fit$vectors, B %*% fit$vectors), diag(2),
               tolerance = 1e-8)
  expect_identical(fit$transform$certification$problem, "original")
  expect_identical(fit$transform$certification$residual_formula,
                   "A * x - lambda * B * x")
  expect_false(fit$transform$certification$transformed_residuals_used)
  expect_true(fit$plan$controls$certified_in_original_coordinates)
  expect_equal(fit$transform$factorization_cache$label_kind,
               "dense_lu_generalized_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_true(fit$transform$factorization_cache$generalized)
  expect_equal(fit$transform$factorization_cache$metric_factorization,
               "LAPACK dpotrf(B)")
  expect_identical(fit$restart$kind,
                   "native_dense_generalized_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
})

test_that("generalized shift-invert handles sparse A with diagonal B without densifying B", {
  n <- 30L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), rep(2.5, n), rep(-1, n - 1))
  )
  B <- Matrix::Diagonal(x = seq(1, 2, length.out = n))
  oracle <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
  sigma <- sort(oracle$values)[4L] + 0.01

  fit <- eig_partial(A, B = B, k = 3L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-8,
                     allow_dense_fallback = "never")
  expected <- oracle$values[order(abs(oracle$values - sigma))][1:3]

  expect_identical(
    fit$method,
    eigencore:::native_tridiagonal_generalized_shift_invert_label()
  )
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-6)
  expect_certificate_clean(fit)
  expect_equal(crossprod(fit$vectors, as.matrix(B) %*% fit$vectors), diag(3),
               tolerance = 1e-8)
  expect_equal(fit$transform$factorization_cache$label_kind,
               "tridiagonal_thomas_generalized_native")
  expect_true(fit$transform$factorization_cache$native)
  expect_equal(fit$transform$factorization_cache$factorization,
               "native tridiagonal Thomas + diagonal sqrt(B)")
  expect_equal(fit$transform$factorization_cache$metric_factorization,
               "diagonal sqrt(B)")
  expect_identical(fit$restart$kind,
                   "native_tridiagonal_generalized_shift_invert_lanczos")
  expect_true(fit$restart$native)
  expect_true(fit$restart$generalized)
  expect_true(fit$restart$factorization_native)
})

test_that("generalized shift-invert keeps non-tridiagonal sparse A reference-labelled", {
  set.seed(10)
  A <- methods::as(symmetric_with_spectrum(seq(1, 12), seed = 10), "CsparseMatrix")
  B <- Matrix::Diagonal(x = seq(1.1, 1.8, length.out = 12L))
  oracle <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
  sigma <- oracle$values[[5L]] + 0.02

  fit <- eig_partial(A, B = B, k = 2L, target = nearest(sigma),
                     method = shift_invert(sigma = sigma), tol = 1e-8,
                     allow_dense_fallback = "never")
  expected <- oracle$values[order(abs(oracle$values - sigma))][1:2]

  expect_identical(
    fit$method,
    "reference generalized SPD Lanczos shift-invert (sparse LU)"
  )
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-6)
  expect_equal(fit$transform$factorization_cache$label_kind,
               "sparse_lu_generalized")
  expect_false(isTRUE(fit$transform$factorization_cache$native))
  contract <- fit$transform$factorization_cache$contract
  expect_equal(contract$provider, "Matrix::lu_reference_factorization")
  expect_equal(contract$promotion_status, "reference_boundary")
  expect_false(contract$owned_by_eigencore)
  expect_true(contract$generalized)
  expect_equal(contract$memory_policy, "sparse_factorization_no_dense_rcond")
})

test_that("generalized shift-invert rejects unsupported sparse B without densifying", {
  n <- 8L
  A <- Matrix::Diagonal(x = seq_len(n) + 1)
  B <- methods::as(Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-0.05, n - 1L),
      rep(1.5, n),
      rep(-0.05, n - 1L)
    )
  ), "dgCMatrix")

  expect_error(
    eig_partial(A, B = B, k = 2L, target = nearest(3),
                method = shift_invert(sigma = 3),
                allow_dense_fallback = "never"),
    "dense B or diagonal B|avoid silent densification"
  )
})

test_that("planner emits an honest shift-invert label for the chosen path", {
  set.seed(13)
  vals <- seq(1, 12)
  A    <- symmetric_with_spectrum(vals, seed = 13)
  P    <- eigen_problem(A, target = nearest(5),
                        transform = shift_invert(sigma = 5))
  plan <- plan_solver(P, k = 3, method = shift_invert(sigma = 5))
  expect_identical(plan$method,
                   eigencore:::native_dense_shift_invert_label())
  expect_true(plan$controls$certified_in_original_coordinates)
  expect_identical(plan$controls$eigenvalue_recovery,
                   "lambda = sigma + 1 / mu")
  expect_match(plan$controls$certification_policy, "original eigenproblem")
  expect_true(any(grepl("shift-invert transform requested", plan$reasons,
                        fixed = TRUE)))
})

test_that("shift-invert factorization cache keys invalidate on A, sigma, and B changes", {
  A <- diag(c(1, 2, 4, 8))
  B <- diag(c(1, 2, 3, 4))
  A_changed <- A
  A_changed[2, 2] <- 2.5
  B_changed <- B
  B_changed[3, 3] <- 3.5

  key <- eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.25)
  expect_identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.25)
  )
  expect_false(identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(as_operator(A_changed), 0.25)
  ))
  expect_false(identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.5)
  ))
  expect_false(identical(
    eigencore:::shift_invert_factorization_cache_key(
      as_operator(A), 0.25, as_operator(B)
    ),
    eigencore:::shift_invert_factorization_cache_key(
      as_operator(A), 0.25, as_operator(B_changed)
    )
  ))
  expect_true(key$standard_problem)
  expect_null(key$B)
  expect_identical(key$structure, "hermitian")

  general_op <- as_operator(A)
  general_op$structure <- general()
  expect_false(identical(
    key,
    eigencore:::shift_invert_factorization_cache_key(general_op, 0.25)
  ))
})

test_that("shift-invert result exposes factorization-cache provenance", {
  vals <- seq(1, 8)
  A <- symmetric_with_spectrum(vals, seed = 91)

  fit <- eig_partial(A, k = 2L, target = nearest(4.2),
                     method = shift_invert(sigma = 4.2))
  cache <- fit$transform$factorization_cache

  expect_s3_class(cache$key, "eigencore_shift_invert_cache_key")
  expect_equal(cache$key$sigma, 4.2)
  expect_equal(cache$label_kind, "dense_lu_native")
  expect_equal(cache$factorization, "LAPACK dgetrf/dgetrs")
  expect_true(cache$factorization_cached)
  expect_equal(cache$condition_estimate_type, "dense_lu_pivot_ratio")
  expect_true(is.finite(cache$condition_estimate))
  expect_true(cache$native)
  expect_true(cache$reusable_within_operator)
  expect_false(cache$external_cache)
  contract <- cache$contract
  expect_equal(contract$contract_version, "shift_invert_factorization_contract_v1")
  expect_equal(contract$provider, "eigencore_native_factorization")
  expect_equal(contract$promotion_status, "promoted_native")
  expect_true(contract$owned_by_eigencore)
  expect_false(contract$external_cache)
  expect_equal(contract$cache_key_scope, "A_fingerprint+B_fingerprint+sigma+structure")
  expect_equal(contract$memory_policy, "native_factorized_apply_no_dense_fallback")
  expect_equal(contract$certificate_policy, "original_coordinate_residual_required")
  expect_true(contract$native_label_requires_owned_factorized_apply)
})

test_that("shift-invert near a true eigenvalue surfaces a clear error", {
  set.seed(17)
  vals <- c(1, 2, 3, 4, 5)
  A    <- symmetric_with_spectrum(vals, seed = 17)
  # sigma equal to an eigenvalue makes (A - sigma I) singular; base::solve errors
  expect_error(
    eig_partial(A, k = 2L, target = nearest(3),
                method = shift_invert(sigma = 3)),
    "singular|computationally"
  )
})

test_that("native tridiagonal shift-invert perturbs singular requested sigma", {
  n <- 8L
  L <- Matrix::bandSparse(
    n,
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), rep(2, n), rep(-1, n - 1L))
  )
  oracle <- eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values
  expected <- oracle[order(abs(oracle - 1))][1:2]

  fit <- eig_partial(L, k = 2L, target = nearest(1), seed = 27L,
                     allow_dense_fallback = "never", tol = 1e-10)

  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_false(isTRUE(all.equal(fit$sigma, 1)))
  expect_equal(fit$transform$requested_sigma, 1)
  expect_true(fit$transform$sigma_perturbed)
  expect_true(is.finite(fit$transform$sigma_perturbation))
  expect_match(fit$warnings, "perturbed requested sigma")
  expect_equal(sort(fit$values), sort(expected), tolerance = 1e-9)
  expect_certificate_clean(fit, tol = 1e-8)
  expect_identical(fit$restart$kind, "native_tridiagonal_shift_invert_lanczos")
  expect_true(fit$restart$factorization_native)
})

test_that("shift-invert recovers smallest eigenvalues of a 1D Laplacian", {
  # discrete 1D Laplacian -- the canonical Milestone L test case
  n <- 50L
  L <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1))
  )
  fit <- eig_partial(L, k = 4L, target = nearest(0.01),
                     method = shift_invert(sigma = 0.01))
  oracle <- sort(eigen(as.matrix(L), symmetric = TRUE)$values)[1:4]
  expect_identical(fit$method,
                   eigencore:::native_tridiagonal_shift_invert_label())
  expect_equal(sort(fit$values), sort(oracle), tolerance = 1e-6)
  expect_certificate_clean(fit)
  expect_identical(fit$transform$label_kind, "tridiagonal_thomas_native")
  expect_true(fit$restart$native)
  expect_true(fit$restart$factorization_native)
})
