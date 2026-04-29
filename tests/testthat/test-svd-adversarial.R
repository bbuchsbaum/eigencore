test_that("planner labels for SVD paths match the kernel that actually runs", {
  # Dense double + auto() routes to native dense LAPACK SVD fallback unless the
  # smaller Gram matrix is an explicit bounded special case.
  set.seed(101)
  M_dense <- matrix(rnorm(80 * 30), 80, 30)
  fit <- svd_partial(M_dense, rank = 5, target = largest())
  expect_identical(fit$plan$method, "native certified Gram SVD special case")
  expect_identical(fit$method,      "native certified Gram SVD special case")
  expect_identical(fit$plan$controls$gram_side, "right")
  expect_true(fit$plan$controls$certified_in_original_coordinates)
  expect_certificate_clean(fit)

  M_squareish <- matrix(rnorm(45 * 30), 45, 30)
  fit_squareish <- svd_partial(M_squareish, rank = 5, target = largest())
  expect_identical(fit_squareish$plan$method, "native dense LAPACK SVD fallback")
  expect_identical(fit_squareish$method,      "native dense LAPACK SVD fallback")
  expect_certificate_clean(fit_squareish)

  # Dense double + explicit golub_kahan() routes to native prototype Golub-Kahan
  fit_gk <- svd_partial(M_dense, rank = 5, target = largest(),
                        method = golub_kahan())
  expect_identical(fit_gk$plan$method, "native prototype Golub-Kahan")
  expect_identical(fit_gk$method,      "native prototype Golub-Kahan")
  expect_certificate_clean(fit_gk)

  # Sparse CSC + auto() can route to the certified Gram special case when the
  # reduced dimension is small enough to make it an explicit safe plan.
  set.seed(102)
  M_csc <- Matrix::rsparsematrix(100, 30, density = 0.05)
  fit_csc <- svd_partial(M_csc, rank = 4, target = largest())
  expect_identical(fit_csc$plan$method, "native certified Gram SVD special case")
  expect_identical(fit_csc$method,      "native certified Gram SVD special case")
  expect_identical(fit_csc$restart$kind, "gram_svd_special_case")
  expect_identical(fit_csc$plan$controls$gram_side, "right")
  expect_identical(fit_csc$plan$controls$gram_dimension, 30L)
  expect_true(fit_csc$plan$controls$certified_in_original_coordinates)
  expect_true(fit_csc$restart$certified_in_original_coordinates)
  expect_certificate_clean(fit_csc)

  fit_csc_gk <- svd_partial(M_csc, rank = 4, target = largest(),
                            method = golub_kahan())
  expect_identical(fit_csc_gk$plan$method, "native prototype Golub-Kahan")
  expect_identical(fit_csc_gk$method,      "native prototype Golub-Kahan")
  expect_false(fit_csc_gk$plan$controls$default_normal_equations)
  expect_certificate_clean(fit_csc_gk)
})

test_that("matrix-free operator with adjoint runs the reference Golub-Kahan path", {
  set.seed(103)
  m <- 60L; n <- 25L
  source_M <- matrix(rnorm(m * n), m, n)
  op <- linear_operator(
    dim = c(m, n),
    apply         = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (source_M %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (crossprod(source_M, X))
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    name = "matrix-free wrapper"
  )
  fit <- svd_partial(op, rank = 4, target = largest())
  expect_identical(fit$plan$method, "prototype Golub-Kahan")
  oracle <- svd(source_M, nu = 4, nv = 4)
  expect_equal(fit$d, oracle$d[1:4], tolerance = 1e-6)
})

test_that("tall-skinny sparse SVD certifies its top triplets", {
  set.seed(104)
  m <- 200L; n <- 40L
  M <- Matrix::rsparsematrix(m, n, density = 0.05)
  fit <- svd_partial(M, rank = 5, target = largest())
  oracle <- svd(as.matrix(M), nu = 5, nv = 5)
  expect_equal(fit$d, oracle$d[1:5], tolerance = 1e-6)
  expect_certificate_clean(fit)
})

test_that("wide-short dense SVD certifies via golub_kahan", {
  set.seed(105)
  m <- 30L; n <- 100L
  M <- matrix(rnorm(m * n), m, n)
  fit <- svd_partial(M, rank = 5, target = largest(), method = golub_kahan())
  oracle <- svd(M, nu = 5, nv = 5)
  expect_equal(fit$d, oracle$d[1:5], tolerance = 1e-6)
  expect_certificate_clean(fit)
})

test_that("near-rank-deficient rectangular SVD returns finite, sorted triplets", {
  s <- c(8, 4, 2, 1e-12, 1e-13)  # numerical near-zero tail
  A <- rectangular_with_singular_values(s, m = 12L, n = 8L, seed = 106)
  fit <- svd_partial(A, rank = 4, target = largest(), method = golub_kahan(),
                     tol = 1e-8)
  expect_true(all(is.finite(fit$d)))
  expect_true(!is.unsorted(rev(fit$d)))  # decreasing
  expect_equal(fit$d[1:3], c(8, 4, 2), tolerance = 1e-8)
})

test_that("tiny gap separating retained from discarded SVs does not break sorting", {
  # rank=3 retained values are 5, 4.999999, 4.999998; tail values are 1, 0.5
  s <- c(5, 5 - 1e-6, 5 - 2e-6, 1, 0.5)
  A <- rectangular_with_singular_values(s, m = 10L, n = 7L, seed = 107)
  fit <- svd_partial(A, rank = 3, target = largest(), tol = 1e-9)
  oracle <- svd(A, nu = 3, nv = 3)
  expect_equal(fit$d, oracle$d[1:3], tolerance = 1e-7)
  expect_lt(subspace_distance(left_vectors(fit),  oracle$u[, 1:3]), 1e-5)
  expect_lt(subspace_distance(right_vectors(fit), oracle$v[, 1:3]), 1e-5)
  expect_certificate_clean(fit)
})

test_that("sparse rank-deficient CSC SVD does not silently densify", {
  set.seed(108)
  m <- 80L; n <- 30L
  # rank-15 ground truth via low-rank construction
  L <- Matrix::rsparsematrix(m, 15L, density = 0.2)
  R <- Matrix::rsparsematrix(15L, n, density = 0.2)
  M <- L %*% R
  expect_s4_class(M, "dgCMatrix")
  fit <- svd_partial(M, rank = 5, target = largest())
  expect_true(all(is.finite(fit$d)))
  oracle <- svd(as.matrix(M), nu = 5, nv = 5)
  expect_equal(fit$d, oracle$d[1:5], tolerance = 1e-6)
})

test_that("sparse Gram SVD certifies exact zero singular triplets", {
  set.seed(514)
  m <- 80L; n <- 30L
  L <- Matrix::rsparsematrix(m, 3L, density = 0.2)
  R <- Matrix::rsparsematrix(3L, n, density = 0.2)
  M <- L %*% R

  fit <- svd_partial(M, rank = 5, target = largest(), tol = 1e-8)

  expect_identical(fit$method, "native certified Gram SVD special case")
  expect_true(fit$restart$zero_singular_completion)
  expect_true(all(is.finite(fit$d)))
  expect_equal(fit$d[4:5], c(0, 0), tolerance = 1e-10)
  expect_certificate_clean(fit)
})

test_that("native Golub-Kahan completes exact zero singular triplets", {
  set.seed(517)
  M <- Matrix::rsparsematrix(160L, 4L, density = 0.1) %*%
    Matrix::rsparsematrix(4L, 50L, density = 0.1)

  fit <- svd_partial(
    M,
    rank = 6L,
    target = largest(),
    method = golub_kahan(),
    tol = 1e-8,
    seed = 702
  )

  expect_identical(fit$method, "native prototype Golub-Kahan")
  expect_length(fit$d, 6L)
  expect_equal(fit$nconv, 6L)
  expect_true(fit$restart$zero_singular_completion)
  expect_true(all(tail(fit$d, 2L) <= fit$restart$zero_singular_threshold))
  expect_certificate_clean(fit)
})

test_that("Gram SVD uses selected dense eigensolve for top singular values", {
  set.seed(516)
  M <- Matrix::rsparsematrix(100, 30, density = 0.05)
  fit <- svd_partial(M, rank = 4, target = largest(), tol = 1e-8)
  oracle <- svd(as.matrix(M), nu = 4, nv = 4)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_equal(fit$d, oracle$d[1:4], tolerance = 1e-8)
  expect_certificate_clean(fit)
})

test_that("wide sparse Gram SVD uses native CSC left-Gram kernel", {
  set.seed(517)
  M <- Matrix::t(Matrix::rsparsematrix(100, 30, density = 0.05))
  fit <- svd_partial(M, rank = 4, target = largest(), tol = 1e-8)
  oracle <- svd(as.matrix(M), nu = 4, nv = 4)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_side, "left")
  expect_identical(fit$restart$native_gram_kernel, "csc_left_gram")
  expect_identical(fit$restart$native_gram_eigensolver, "lapack_dsyevr")
  expect_true(is.infinite(fit$restart$native_gram_subspace_max_backward_error))
  expect_false(fit$restart$normal_operator_implicit)
  expect_true(fit$restart$materialized_gram)
  expect_identical(fit$certificate$norm_bound_type, "frobenius_exact")
  expect_true(all(c("gram", "eigensolve", "vector_form", "diagnostics") %in%
                    names(fit$stage_seconds)))
  expect_true(all(is.finite(fit$stage_seconds)))
  expect_true(all(fit$stage_seconds >= 0))
  expect_equal(fit$d, oracle$d[1:4], tolerance = 1e-8)
  expect_certificate_clean(fit)
})

test_that("wide sparse Gram SVD exposes opt-in certified subspace eigensolve", {
  old_options <- options(eigencore.csc_left_gram_subspace_attempt = TRUE)
  on.exit(options(old_options), add = TRUE)
  M <- cbind(
    Matrix::Diagonal(x = c(10, 8, 6, 1, 0.5)),
    Matrix::Matrix(0, 5, 20, sparse = TRUE)
  )

  fit <- svd_partial(M, rank = 2, target = largest(), tol = 1e-8)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$native_gram_eigensolver, "subspace_iteration")
  expect_lte(fit$restart$native_gram_subspace_max_backward_error, 1e-8)
  expect_equal(fit$d, c(10, 8), tolerance = 1e-10)
  expect_certificate_clean(fit)
})

test_that("wide sparse Gram SVD exposes opt-in implicit normal Lanczos", {
  old_options <- options(
    eigencore.csc_left_normal_lanczos_attempt = TRUE,
    eigencore.csc_left_gram_subspace_attempt = FALSE
  )
  on.exit(options(old_options), add = TRUE)
  M <- cbind(
    Matrix::Diagonal(x = c(10, 8, 6, 1, 0.5)),
    Matrix::Matrix(0, 5, 20, sparse = TRUE)
  )

  fit <- svd_partial(M, rank = 2, target = largest(), tol = 1e-8)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$native_gram_eigensolver, "implicit_normal_lanczos")
  expect_true(fit$restart$normal_operator_implicit)
  expect_false(fit$restart$materialized_gram)
  expect_lte(fit$restart$native_implicit_normal_lanczos_max_backward_error, 1e-8)
  expect_lte(fit$restart$native_implicit_normal_lanczos_iterations, 5L)
  expect_equal(fit$d, c(10, 8), tolerance = 1e-10)
  expect_certificate_clean(fit)
})

test_that("sparse Gram SVD falls back when certification is weaker than Golub-Kahan", {
  set.seed(2)
  U <- Matrix::rsparsematrix(120, 2, density = 0.2)
  V <- Matrix::rsparsematrix(2, 30, density = 0.2)
  M <- U %*% Matrix::Diagonal(x = c(1, 1e-6)) %*% V

  fit <- svd_partial(M, rank = 2, target = largest(), tol = 1e-10)

  expect_identical(fit$plan$method, "native certified Gram SVD special case")
  expect_identical(fit$method, "native prototype Golub-Kahan fallback from Gram SVD")
  expect_true(fit$restart$fallback_attempted)
  expect_true(fit$restart$fallback_used)
  expect_equal(fit$restart$fallback_method, "native prototype Golub-Kahan")
  expect_false(fit$restart$gram_certificate_passed)
  expect_true(is.finite(fit$restart$fallback_max_backward_error))
  expect_true(fit$certificate$passed)
})

test_that("svd_partial honors vectors='left' / 'right' / 'none' modes", {
  set.seed(109)
  M <- matrix(rnorm(50 * 20), 50, 20)
  left  <- svd_partial(M, rank = 3, vectors = "left",  method = golub_kahan())
  right <- svd_partial(M, rank = 3, vectors = "right", method = golub_kahan())
  none  <- svd_partial(M, rank = 3, vectors = "none",  method = golub_kahan())
  expect_false(is.null(left$u));   expect_null(left$v)
  expect_null(right$u);            expect_false(is.null(right$v))
  expect_null(none$u);             expect_null(none$v)
})

test_that("randomized SVD refine option can recover a certified result", {
  set.seed(511)
  M <- Matrix::rsparsematrix(120, 30, density = 0.08)
  rough <- svd_partial(
    M,
    rank = 4,
    method = randomized(oversample = 2, n_iter = 0, refine = FALSE),
    tol = 1e-8,
    seed = 511
  )
  refined <- svd_partial(
    M,
    rank = 4,
    method = randomized(oversample = 2, n_iter = 0, refine = TRUE),
    tol = 1e-8,
    seed = 511
  )

  expect_false(rough$restart$refinement_attempted)
  expect_true(refined$restart$refinement_attempted)
  expect_true(refined$restart$refinement_passed)
  expect_equal(refined$restart$refinement_kind, "gram_svd_special_case")
  expect_true(refined$certificate$passed)
  expect_lte(refined$certificate$max_backward_error,
             rough$certificate$max_backward_error)
})

test_that("randomized SVD records and honors normalizer choices", {
  set.seed(512)
  U <- qr.Q(qr(matrix(rnorm(60), nrow = 20, ncol = 3)))
  V <- qr.Q(qr(matrix(rnorm(36), nrow = 12, ncol = 3)))
  M <- U %*% diag(c(7, 4, 1), nrow = 3) %*% t(V)

  for (normalizer in c("qr", "lu", "none")) {
    fit <- svd_partial(
      M,
      rank = 2,
      method = randomized(
        oversample = 4,
        n_iter = 1,
        normalizer = normalizer,
        refine = FALSE
      ),
      tol = 1e-8,
      seed = 512
    )
    expect_equal(fit$restart$normalizer, normalizer)
    expect_equal(fit$restart$apply_kind, "dense_direct")
    expect_true(fit$restart$certificate_reuses_projection)
    expect_equal(fit$restart$adaptive_stop, TRUE)
    expect_equal(fit$restart$adaptive_stop_used, identical(normalizer, "qr"))
    expect_true(fit$certificate$passed)
    expect_equal(fit$d, c(7, 4), tolerance = 1e-8)
  }
})

test_that("randomized SVD stops after q0 when projected certificate already passes", {
  set.seed(1904)
  A <- rectangular_with_singular_values(
    c(9, 7, 5, 3, 2, 1, rep(0, 34)),
    m = 80L,
    n = 40L,
    seed = 1904
  )
  fit <- svd_partial(
    A,
    rank = 5,
    method = randomized(oversample = 10, n_iter = 2, refine = FALSE),
    tol = 1e-8,
    seed = 1904
  )

  expect_true(fit$restart$adaptive_stop)
  expect_true(fit$restart$adaptive_stop_used)
  expect_equal(fit$restart$iterations_used, 1L)
  expect_equal(fit$iterations, 1L)
  expect_true(fit$certificate$passed)
  expect_equal(fit$d, c(9, 7, 5, 3, 2), tolerance = 1e-8)
})

test_that("randomized SVD projected certificate matches direct recomputation", {
  set.seed(1905)
  A <- matrix(rnorm(70 * 40), nrow = 70, ncol = 40)
  fit <- svd_partial(
    A,
    rank = 5,
    method = randomized(oversample = 8, n_iter = 1, refine = FALSE),
    tol = 1e-8,
    seed = 1905
  )
  direct <- eigencore:::certify_svd_operator(
    as_operator(A),
    fit$d,
    fit$u,
    fit$v,
    tol = 1e-8
  )

  expect_true(fit$restart$certificate_reuses_projection)
  expect_equal(fit$certificate$residuals$left, direct$residuals$left, tolerance = 1e-10)
  expect_equal(fit$certificate$residuals$right, direct$residuals$right, tolerance = 1e-10)
  expect_equal(fit$certificate$backward_error, direct$backward_error, tolerance = 1e-10)
  expect_equal(fit$certificate$passed, direct$passed)
})

test_that("SVD planner records inspectable method controls", {
  set.seed(513)
  M <- Matrix::rsparsematrix(120, 30, density = 0.05)

  gram_plan <- plan_solver(svd_problem(M), rank = 4)
  expect_identical(gram_plan$method, "native certified Gram SVD special case")
  expect_identical(gram_plan$controls$gram_side, "right")
  expect_identical(gram_plan$controls$gram_dimension, 30L)
  expect_identical(gram_plan$controls$materializes, "smaller Gram matrix only")
  expect_identical(gram_plan$controls$fallback_policy, "certification-gated")
  expect_match(gram_plan$controls$runtime_fallback, "Golub-Kahan")

  dense_gram_plan <- plan_solver(
    svd_problem(matrix(rnorm(80 * 30), 80, 30)),
    rank = 4
  )
  expect_identical(dense_gram_plan$method, "native certified Gram SVD special case")
  expect_identical(dense_gram_plan$controls$gram_dimension, 30L)

  gk_plan <- plan_solver(
    svd_problem(M),
    rank = 4,
    method = golub_kahan(max_subspace = 11L, reorthogonalize = FALSE)
  )
  expect_identical(gk_plan$method, "native prototype Golub-Kahan")
  expect_identical(gk_plan$controls$max_subspace, 11L)
  expect_false(gk_plan$controls$adaptive_subspace)
  expect_identical(gk_plan$controls$initial_max_subspace, 11L)
  expect_false(gk_plan$controls$reorthogonalize)
  expect_true(gk_plan$controls$requires_adjoint)
  expect_false(gk_plan$controls$default_normal_equations)

  adaptive_gk_plan <- plan_solver(
    svd_problem(M),
    rank = 4,
    method = golub_kahan()
  )
  expect_identical(adaptive_gk_plan$method, "native prototype Golub-Kahan")
  expect_true(adaptive_gk_plan$controls$adaptive_subspace)
  expect_null(adaptive_gk_plan$controls$max_subspace)
  expect_gte(adaptive_gk_plan$controls$initial_max_subspace, 4L)

  randomized_plan <- plan_solver(
    svd_problem(M),
    rank = 4,
    method = randomized(oversample = 3L, n_iter = 1L, normalizer = "lu", refine = TRUE)
  )
  expect_identical(randomized_plan$method, "reference randomized SVD prototype")
  expect_identical(randomized_plan$controls$oversample, 3L)
  expect_identical(randomized_plan$controls$n_iter, 1L)
  expect_identical(randomized_plan$controls$sample_dimension, 7L)
  expect_identical(randomized_plan$controls$normalizer, "lu")
  expect_true(randomized_plan$controls$approximate)
  expect_false(randomized_plan$controls$auto_selected)
  expect_true(randomized_plan$controls$refine)
  expect_match(randomized_plan$controls$certification_policy, "residual certificate")
  expect_match(randomized_plan$controls$certification_refinement, "native Gram SVD")
})

test_that("randomized SVD reports approximation and certificate policy", {
  set.seed(515)
  M <- Matrix::rsparsematrix(120, 30, density = 0.08)
  fit <- svd_partial(
    M,
    rank = 4,
    method = randomized(oversample = 2, n_iter = 0, refine = FALSE),
    tol = 1e-8,
    seed = 515
  )

  expect_true(fit$restart$approximate)
  expect_equal(fit$restart$apply_kind, "csc_direct")
  expect_true(fit$restart$certificate_reuses_projection)
  expect_match(fit$restart$certificate_policy, "stochastic sketch is not sufficient")
  expect_false(fit$restart$refine)
  expect_false(fit$restart$refinement_attempted)
  expect_true(is.logical(fit$restart$initial_certificate_passed))
  expect_true(is.finite(fit$restart$initial_max_backward_error))
})
