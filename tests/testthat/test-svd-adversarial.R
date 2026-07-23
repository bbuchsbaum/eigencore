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
  expect_true(isTRUE(fit_csc$fastpath_native_result))
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

test_that("matrix-free operator with adjoint runs the native callback Golub-Kahan path", {
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
  expect_identical(fit$plan$method, eigencore:::native_matrix_free_golub_kahan_label())
  expect_true(fit$restart$native_callback)
  expect_true(fit$restart$callback_boundary)
  oracle <- svd(source_M, nu = 4, nv = 4)
  expect_equal(fit$d, oracle$d[1:4], tolerance = 1e-6)
})

test_that("tall-skinny sparse SVD defaults to the explicit right Gram solve", {
  set.seed(104)
  m <- 200L; n <- 40L
  M <- Matrix::rsparsematrix(m, n, density = 0.05)
  fit <- svd_partial(M, rank = 5, target = largest())
  oracle <- svd(as.matrix(M), nu = 5, nv = 5)
  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_side, "right")
  expect_identical(fit$restart$native_gram_kernel, "csc_right_gram")
  expect_true(isTRUE(fit$fastpath_native_result))
  expect_identical(fit$restart$native_gram_eigensolver, "lapack_dsyevr")
  expect_false(fit$restart$normal_operator_implicit)
  expect_true(fit$restart$materialized_gram)
  expect_false(fit$restart$explicit_gram_retry_used)
  expect_true(all(c("gram", "eigensolve", "vector_form", "diagnostics") %in%
                    names(fit$stage_seconds)))
  expect_true(all(is.finite(fit$stage_seconds)))
  expect_true(all(fit$stage_seconds >= 0))
  expect_equal(fit$d, oracle$d[1:5], tolerance = 1e-6)
  expect_certificate_clean(fit)
})

test_that("tall-skinny sparse SVD can disable right-normal diagnostic", {
  old_options <- options(eigencore.csc_right_normal_lanczos_attempt = FALSE)
  on.exit(options(old_options), add = TRUE)
  set.seed(104)
  M <- Matrix::rsparsematrix(200L, 40L, density = 0.05)

  fit <- svd_partial(M, rank = 5, target = largest())

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_side, "right")
  expect_identical(fit$restart$native_gram_eigensolver, "lapack_dsyevr")
  expect_false(fit$restart$normal_operator_implicit)
  expect_true(fit$restart$materialized_gram)
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

test_that("sparse Gram SVD completes near-null rank-deficient triplets without fallback", {
  set.seed(703)
  L <- Matrix::rsparsematrix(5000L, 20L, density = 0.01)
  R <- Matrix::rsparsematrix(20L, 500L, density = 0.01)
  M <- L %*% R

  fit <- svd_partial(M, rank = 30L, target = largest(), tol = 1e-8, seed = 701L)

  expect_identical(fit$method, "native certified Gram SVD special case")
  expect_true(fit$restart$zero_singular_completion)
  expect_false(fit$restart$fallback_attempted)
  expect_false(fit$restart$fallback_used)
  expect_true(fit$restart$gram_certificate_passed)
  expect_equal(fit$d[21:30], rep(0, 10L), tolerance = 1e-10)
  expect_certificate_clean(fit)
})

test_that("smallest and interior dense SVD targets have exact certificate policy", {
  M <- diag(c(10, 5, 1, 0.2, 0.1))

  smallest_fit <- svd_partial(M, rank = 2L, target = smallest(), tol = 1e-10)
  expect_identical(smallest_fit$method, "native dense LAPACK SVD fallback")
  expect_equal(smallest_fit$d, c(0.1, 0.2), tolerance = 1e-12)
  expect_identical(smallest_fit$plan$controls$svd_target_family, "smallest")
  expect_identical(
    smallest_fit$plan$controls$svd_target_certificate_policy,
    "exact two-sided residual certificate in original coordinates"
  )
  expect_certificate_clean(smallest_fit)

  interior_fit <- svd_partial(M, rank = 2L, target = nearest(0.8), tol = 1e-10)
  expect_identical(interior_fit$method, "native dense LAPACK SVD fallback")
  expect_equal(interior_fit$d, c(1, 0.2), tolerance = 1e-12)
  expect_identical(interior_fit$plan$controls$svd_target_family, "interior")
  expect_identical(interior_fit$plan$controls$svd_target_boundary, "dense exact fallback")
  expect_certificate_clean(interior_fit)
})

test_that("sparse interior SVD uses native full-subspace boundary", {
  M <- Matrix::sparseMatrix(
    i = 1:5,
    j = 1:5,
    x = c(10, 5, 1, 0.2, 0.1),
    dims = c(8L, 5L)
  )

  fit <- svd_partial(
    M,
    rank = 2L,
    target = nearest(0.8),
    tol = 1e-10,
    seed = 90,
    allow_dense_fallback = "never"
  )

  expect_identical(fit$method, eigencore:::native_interior_golub_kahan_label())
  expect_identical(
    fit$plan$controls$svd_target_boundary,
    "native full-subspace interior SVD boundary"
  )
  expect_true(fit$plan$controls$full_subspace_interior)
  expect_equal(fit$restart$final_max_subspace, min(dim(M)))
  expect_false(fit$restart$projected_stop_enabled)
  expect_equal(fit$d, c(1, 0.2), tolerance = 1e-12)
  expect_certificate_clean(fit)
})

test_that("diagonal interior SVD reference boundary certifies without dense fallback", {
  M <- Matrix::Diagonal(x = c(10, 5, 1, 0.2, 0.1))

  fit <- svd_partial(
    M,
    rank = 1L,
    target = nearest(1),
    tol = 1e-10,
    allow_dense_fallback = "never"
  )

  expect_identical(fit$method, "prototype Golub-Kahan")
  expect_identical(fit$plan$controls$svd_target_family, "interior")
  expect_identical(
    fit$plan$controls$svd_target_boundary,
    "reference/prototype interior selection"
  )
  expect_equal(fit$d, 1, tolerance = 1e-12)
  expect_certificate_clean(fit)
})

test_that("smallest sparse CSC SVD uses native certified production boundary", {
  M <- Matrix::sparseMatrix(
    i = 1:6,
    j = 1:6,
    x = c(10, 5, 2, 1, 0.4, 0.1),
    dims = c(8L, 6L)
  )

  fit <- svd_partial(
    M,
    rank = 2L,
    target = smallest(),
    tol = 1e-10,
    seed = 91,
    allow_dense_fallback = "never"
  )

  expect_identical(fit$method, eigencore:::native_smallest_golub_kahan_label())
  expect_identical(
    fit$plan$controls$svd_target_boundary,
    "native certified smallest SVD production boundary"
  )
  expect_identical(
    fit$plan$controls$promotion_gate_issue,
    "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
  )
  expect_identical(
    fit$plan$controls$closed_decision_issue,
    "bd-01KTEH6862GB19JJWX2M3FQP6T"
  )
  expect_equal(fit$d, c(0.1, 0.4), tolerance = 1e-12)
  expect_true(fit$restart$native)
  expect_false(fit$restart$matrix_free)
  expect_certificate_clean(fit)
})

test_that("smallest tall sparse CSC SVD uses native Gram production boundary", {
  M <- Matrix::sparseMatrix(
    i = 1:6,
    j = 1:6,
    x = c(10, 5, 2, 1, 0.4, 0.1),
    dims = c(20L, 6L)
  )

  fit <- svd_partial(
    M,
    rank = 2L,
    target = smallest(),
    tol = 1e-10,
    seed = 92,
    allow_dense_fallback = "never"
  )

  expect_identical(fit$method, "native certified Gram SVD special case")
  expect_identical(fit$plan$controls$svd_target_family, "smallest")
  expect_identical(
    fit$plan$controls$svd_target_boundary,
    "native certified smallest SVD production boundary"
  )
  expect_identical(fit$plan$controls$promotion_status, "production_smallest_gram_certified")
  expect_identical(
    fit$plan$controls$promotion_gate_issue,
    "bd-01KTE8G6RYE4RD5F6CN7SNKKC6"
  )
  expect_identical(
    fit$plan$controls$closed_decision_issue,
    "bd-01KTEH6862GB19JJWX2M3FQP6T"
  )
  expect_equal(fit$d, c(0.1, 0.4), tolerance = 1e-12)
  expect_true(fit$restart$native)
  expect_equal(fit$restart$kind, "gram_svd_special_case")
  expect_equal(fit$restart$gram_side, "right")
  expect_equal(fit$restart$native_gram_kernel, "materialized_right_gram")
  expect_equal(fit$restart$native_gram_eigensolver, "lapack_dsyev_full")
  expect_true(fit$restart$materialized_gram)
  expect_certificate_clean(fit)
})

test_that("smallest matrix-free SVD promotion requires exact norm metadata", {
  A <- rbind(diag(c(10, 5, 2, 1, 0.4, 0.1)), matrix(0, 2, 6))
  make_op <- function(metadata = list()) {
    linear_operator(
      dim = dim(A),
      apply = function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- alpha * (A %*% X)
        if (!is.null(Y) && beta != 0) out <- out + beta * Y
        out
      },
      apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- alpha * (t(A) %*% X)
        if (!is.null(Y) && beta != 0) out <- out + beta * Y
        out
      },
      structure = general(),
      metadata = metadata
    )
  }

  no_metadata_plan <- plan_solver(
    svd_problem(make_op(), target = smallest()),
    rank = 2L
  )
  expect_identical(
    no_metadata_plan$method,
    eigencore:::native_matrix_free_golub_kahan_label()
  )

  fit <- svd_partial(
    make_op(metadata = list(frobenius_norm = norm(A, "F"))),
    rank = 2L,
    target = smallest(),
    tol = 1e-10,
    seed = 92,
    allow_dense_fallback = "never"
  )

  expect_identical(
    fit$method,
    eigencore:::native_matrix_free_smallest_golub_kahan_label()
  )
  expect_identical(fit$certificate$norm_bound_type, "frobenius_metadata")
  expect_false(fit$certificate$scale_is_estimate)
  expect_true(fit$restart$matrix_free)
  expect_true(fit$restart$native_callback)
  expect_true(fit$plan$controls$requires_nonestimated_norm_scale)
  expect_equal(fit$d, c(0.1, 0.4), tolerance = 1e-12)
  expect_certificate_clean(fit)
})

test_that("interior matrix-free SVD promotion requires exact norm metadata", {
  A <- rbind(diag(c(10, 5, 1, 0.2, 0.1)), matrix(0, 3, 5))
  make_op <- function(metadata = list()) {
    linear_operator(
      dim = dim(A),
      apply = function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- alpha * (A %*% X)
        if (!is.null(Y) && beta != 0) out <- out + beta * Y
        out
      },
      apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
        out <- alpha * (t(A) %*% X)
        if (!is.null(Y) && beta != 0) out <- out + beta * Y
        out
      },
      structure = general(),
      metadata = metadata
    )
  }

  expect_error(
    svd_partial(
      make_op(),
      rank = 2L,
      target = nearest(0.8),
      tol = 1e-10,
      seed = 93,
      allow_dense_fallback = "never"
    ),
    "Interior SVD target nearest\\(0.8\\) is not supported by native matrix-free Golub-Kahan callback"
  )

  fit <- svd_partial(
    make_op(metadata = list(frobenius_norm = norm(A, "F"))),
    rank = 2L,
    target = nearest(0.8),
    tol = 1e-10,
    seed = 94,
    allow_dense_fallback = "never"
  )

  expect_identical(
    fit$method,
    eigencore:::native_matrix_free_interior_golub_kahan_label()
  )
  expect_identical(
    fit$plan$controls$svd_target_boundary,
    "native full-subspace interior SVD boundary"
  )
  expect_true(fit$plan$controls$requires_nonestimated_norm_scale)
  expect_true(fit$plan$controls$full_subspace_interior)
  expect_identical(fit$certificate$norm_bound_type, "frobenius_metadata")
  expect_false(fit$certificate$scale_is_estimate)
  expect_equal(fit$restart$final_max_subspace, min(dim(A)))
  expect_equal(fit$d, c(1, 0.2), tolerance = 1e-12)
  expect_certificate_clean(fit)
})

test_that("complex dense SVD uses native dense complex certification", {
  A <- matrix(c(1 + 1i, 0, 0, 2), 2, 2)

  fit <- svd_partial(A, rank = 2L, tol = 1e-10)

  expect_identical(fit$method, eigencore:::native_dense_complex_svd_label())
  expect_equal(fit$d, c(2, sqrt(2)), tolerance = 1e-10)
  expect_true(is.complex(fit$u))
  expect_true(is.complex(fit$v))
  expect_true(fit$certificate$passed)
  expect_identical(fit$certificate$norm_bound_type, "frobenius_exact")
  expect_lt(fit$certificate$max_orthogonality_loss, 1e-10)
})

test_that("complex matrix-free SVD fails with actionable boundary", {
  A <- matrix(c(1 + 1i, 0, 0, 2), 2, 2)
  op <- linear_operator(
    dim = c(2, 2),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (Conj(t(A)) %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    dtype = "complex",
    structure = general()
  )

  expect_error(
    svd_partial(op, rank = 1L, tol = 1e-10),
    "Complex matrix-free SVD operators are future scope"
  )
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
  expect_true(fit$fastpath_native_result)
  expect_true(fit$plan$controls$svd_partial_fastpath)
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

test_that("wide sparse Gram SVD exposes guarded DSYEVX backend", {
  old_options <- options(eigencore.csc_left_gram_dsyevx_attempt = TRUE)
  on.exit(options(old_options), add = TRUE)
  set.seed(519)
  M <- Matrix::t(Matrix::rsparsematrix(100, 30, density = 0.05))
  fit <- svd_partial(M, rank = 4, target = largest(), tol = 1e-8)
  oracle <- svd(as.matrix(M), nu = 4, nv = 4)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_side, "left")
  expect_identical(fit$restart$native_gram_eigensolver, "lapack_dsyevx")
  expect_equal(fit$d, oracle$d[1:4], tolerance = 1e-8)
  expect_certificate_clean(fit)
})

test_that("wide sparse Gram SVD dispatches larger tiny ranks to DSYEVD", {
  set.seed(518)
  M <- Matrix::t(Matrix::rsparsematrix(90, 32, density = 0.08))
  fit <- svd_partial(M, rank = 16, target = largest(), tol = 1e-8)
  oracle <- svd(as.matrix(M), nu = 16, nv = 16)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_side, "left")
  expect_identical(fit$restart$native_gram_eigensolver, "lapack_dsyevd")
  expect_equal(fit$d, oracle$d[1:16], tolerance = 1e-8)
  expect_certificate_clean(fit)
})

test_that("native Golub-Kahan supports one-sided small-side reorthogonalization", {
  M <- cbind(
    Matrix::Diagonal(x = c(9, 6, 3, 1, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01)),
    Matrix::Matrix(0, 10, 60, sparse = TRUE)
  )
  fit <- svd_partial(
    M,
    rank = 2,
    target = largest(),
    method = golub_kahan(max_subspace = 8, reorthogonalize = FALSE),
    tol = 1e-8,
    seed = 519
  )

  expect_identical(fit$method, "native prototype Golub-Kahan")
  expect_identical(fit$restart$reorthogonalization_mode, "one_sided_small_side")
  expect_identical(fit$restart$internal_orientation, "transposed_wide_operator")
  expect_true(fit$restart$internal_transposed)
  expect_false(fit$restart$reorthogonalize_u)
  expect_true(fit$restart$reorthogonalize_v)
  expect_equal(fit$d, c(9, 6), tolerance = 1e-10)
  expect_certificate_clean(fit)

  tall <- Matrix::t(M)
  tall_fit <- svd_partial(
    tall,
    rank = 2,
    target = largest(),
    method = golub_kahan(max_subspace = 8, reorthogonalize = FALSE),
    tol = 1e-8,
    seed = 520
  )
  expect_identical(tall_fit$restart$reorthogonalization_mode, "one_sided_small_side")
  expect_identical(tall_fit$restart$internal_orientation, "as_given")
  expect_false(tall_fit$restart$internal_transposed)
  expect_false(tall_fit$restart$reorthogonalize_u)
  expect_true(tall_fit$restart$reorthogonalize_v)
  expect_equal(tall_fit$d, c(9, 6), tolerance = 1e-10)
  expect_certificate_clean(tall_fit)
})

test_that("one-sided IRLBA LBD benchmark policy falls back to certified adaptive work", {
  set.seed(702)
  M <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  fit <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_one_sided",
    M,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )

  expect_true(fit$restart$irlba_lbd_small_work_attempted)
  expect_true(fit$restart$irlba_lbd_fallback_attempted)
  expect_true(fit$restart$irlba_lbd_fallback_used)
  expect_true(fit$restart$irlba_lbd_fallback_warm_started)
  expect_true(fit$restart$warm_started)
  expect_true(fit$restart$fallback_attempted)
  expect_true(fit$restart$fallback_used)
  expect_false(fit$restart$irlba_lbd_small_work_certificate_passed)
  expect_equal(fit$restart$attempted_subspaces, c(12L, 45L))
  expect_true(is.data.frame(fit$restart$attempt_history))
  expect_equal(fit$restart$attempt_history$warm_started, c(FALSE, TRUE))
  expect_false(fit$restart$attempt_history$certificate_passed[[1L]])
  expect_true(fit$restart$attempt_history$certificate_passed[[2L]])
  expect_equal(fit$restart$irlba_lbd_small_work_matvecs, fit$restart$attempt_history$matvecs[[1L]])
  expect_equal(fit$restart$irlba_lbd_fallback_matvecs, fit$restart$attempt_history$matvecs[[2L]])
  expect_gt(fit$restart$irlba_lbd_scout_matvec_overhead_fraction, 0)
  expect_lt(fit$restart$irlba_lbd_scout_matvec_overhead_fraction, 1)
  expect_certificate_clean(fit)
})

test_that("one-sided IRLBA LBD restart ABI fixes the native implementation contract", {
  set.seed(703)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  abi <- eigencore:::native_irlba_lbd_restart_abi(
    wide,
    rank = 5L,
    work = 12L,
    retained = 7L,
    max_restarts = 4L
  )

  expect_s3_class(abi, "eigencore_irlba_lbd_restart_abi")
  expect_equal(abi$version, 1L)
  expect_true(abi$implemented)
  expect_equal(abi$native_storage, "dgCMatrix")
  expect_true(abi$internal_transposed)
  expect_identical(abi$internal_orientation, "transposed_wide_operator")
  expect_equal(abi$original_dim, c(90L, 600L))
  expect_equal(abi$active_dim, c(600L, 90L))
  expect_equal(abi$small_side_dimension, 90L)
  expect_equal(abi$work, 12L)
  expect_equal(abi$retained, 7L)
  expect_equal(abi$input_schema$initial_start, 90L)
  expect_equal(abi$input_schema$retained_right_subspace, c(90L, 7L))
  expect_equal(abi$input_schema$retained_left_subspace, c(600L, 7L))
  expect_true("attempt_history" %in% abi$output_schema)
  expect_true(any(grepl("run internally on A\\^T", abi$invariants)))
  expect_true(any(grepl("rotated together", abi$invariants, fixed = TRUE)))
  expect_true(any(grepl("not thrown away", abi$invariants, fixed = TRUE)))
  expect_equal(unname(abi$entry_points[["csc"]]), "eigencore_irlba_lbd_csc_retained")

  tall <- Matrix::rsparsematrix(600L, 90L, density = 0.03)
  tall_abi <- eigencore:::native_irlba_lbd_restart_abi(tall, rank = 5L)
  expect_false(tall_abi$internal_transposed)
  expect_identical(tall_abi$internal_orientation, "as_given")
  expect_equal(tall_abi$active_dim, c(600L, 90L))
  expect_equal(tall_abi$input_schema$initial_start, 90L)
})

test_that("retained IRLBA LBD native core certifies or falls back honestly", {
  set.seed(704)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  fit <- eigencore:::native_irlba_lbd_retained_svd(
    wide,
    rank = 5L,
    work = 12L,
    retained = 7L,
    max_restarts = 7L,
    tol = 1e-8,
    vectors = "both"
  )

  expect_true(fit$certificate$passed)
  expect_equal(length(fit$d), 5L)
  expect_true(fit$restart$retained_restart)
  expect_true(fit$restart$native_attempt_certification)
  expect_equal(fit$restart$retained_restart_abi_version, 1L)
  expect_false(fit$restart$fallback_used)
  expect_identical(fit$restart$irlba_lbd_restart_state_kind, "residual_augmented_projection")
  expect_true(fit$restart$irlba_lbd_recurrence_available)
  expect_true(fit$restart$irlba_lbd_augmented_recurrence)
  expect_identical(fit$restart$irlba_lbd_retained_seed_strategy, "ritz_residual_augmented_krylov_projection")
  expect_equal(fit$restart$irlba_lbd_retained_from_scout, 5L)
  expect_equal(fit$restart$irlba_lbd_retained_padding, 2L)
  expect_equal(fit$restart$irlba_lbd_residual_augmented_cols, 1L)
  expect_equal(fit$restart$irlba_lbd_augmented_tail_steps, 32L)
  expect_equal(fit$restart$irlba_lbd_augmented_basis_cols, 38L)
  expect_equal(fit$restart$irlba_lbd_augmented_small_svds, 3L)
  expect_equal(fit$restart$irlba_lbd_augmented_cached_aq_cols,
               fit$restart$irlba_lbd_augmented_basis_cols)
  expect_true("certificate_passed" %in% names(fit$restart$attempt_history))
  expect_true("converged_count" %in% names(fit$restart$attempt_history))
  expect_true("leading_converged_count" %in% names(fit$restart$attempt_history))
  expect_true(tail(fit$restart$attempt_history$certificate_passed, 1L))
  expect_equal(fit$restart$attempt_history$iterations, c(30L, 31L, 32L))
  expect_equal(tail(fit$restart$attempt_history$converged_count, 1L), 5L)
  expect_equal(tail(fit$restart$attempt_history$leading_converged_count, 1L), 5L)
  expect_true(any(
    fit$restart$attempt_history$leading_converged_count[
      !fit$restart$attempt_history$certificate_passed
    ] < 5L
  ))
  expect_true(fit$restart$irlba_lbd_augmented_reduces_from_scratch_work)
  expect_gt(fit$restart$irlba_lbd_augmented_matvec_savings, 0L)
  expect_identical(fit$restart$internal_orientation, "transposed_wide_operator")
  expect_true(fit$restart$internal_transposed)
  expect_equal(
    fit$matvecs,
    fit$restart$irlba_lbd_scout_matvecs +
      fit$restart$irlba_lbd_retained_matvecs
  )
  expect_equal(fit$restart$irlba_lbd_total_matvecs, fit$matvecs)
  expect_false(fit$restart$irlba_lbd_scout_certificate_passed)
  expect_certificate_clean(fit)
})

test_that("retained IRLBA benchmark candidate avoids repeated fixed-work native scouts", {
  set.seed(702)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  fit <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_retained_native",
    wide,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )

  expect_true(fit$certificate$passed)
  expect_true(fit$restart$irlba_lbd_retained_native_attempted)
  expect_false(fit$restart$fallback_attempted)
  expect_false(fit$restart$fallback_used)
  expect_identical(fit$restart$irlba_lbd_restart_state_kind, "residual_augmented_projection")
  expect_true(fit$restart$irlba_lbd_recurrence_available)
  expect_true(fit$restart$irlba_lbd_augmented_recurrence)
  expect_equal(fit$restart$irlba_lbd_retained_fixed_work_attempts, 0L)
  expect_equal(fit$restart$irlba_lbd_scout_matvecs, 24L)
  expect_lt(fit$restart$irlba_lbd_retained_matvecs, 136L)
  expect_equal(fit$restart$irlba_lbd_augmented_tail_steps, 30L)
  expect_equal(fit$restart$irlba_lbd_augmented_basis_cols, 36L)
  expect_equal(fit$restart$irlba_lbd_augmented_restart_cycles, 8L)
  expect_equal(fit$restart$irlba_lbd_augmented_kept_vectors, 5L)
  expect_equal(fit$restart$irlba_lbd_augmented_small_svds, 1L)
  expect_equal(fit$restart$irlba_lbd_augmented_cached_aq_cols, 36L)
  expect_true(fit$restart$irlba_lbd_augmented_reduces_from_scratch_work)
  expect_gt(fit$restart$irlba_lbd_augmented_matvec_savings, 0L)
  expect_true(is.finite(fit$restart$irlba_lbd_augmented_min_cheap_residual))
  expect_true("certificate_passed" %in% names(fit$restart$attempt_history))
  expect_true("converged_count" %in% names(fit$restart$attempt_history))
  expect_true("leading_converged_count" %in% names(fit$restart$attempt_history))
  expect_true(fit$restart$attempt_history$certificate_passed[[1L]])
  expect_equal(fit$restart$attempt_history$converged_count[[1L]], 5L)
  expect_equal(fit$restart$attempt_history$leading_converged_count[[1L]], 5L)
  expect_equal(
    fit$matvecs,
    fit$restart$irlba_lbd_scout_matvecs +
      fit$restart$irlba_lbd_retained_matvecs
  )
  expect_equal(fit$restart$irlba_lbd_total_matvecs, fit$matvecs)
  expect_certificate_clean(fit)
})

test_that("retained IRLBA BPRO policy certifies with monitored partial reorthogonalization", {
  set.seed(702)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  full <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_retained_native",
    wide,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )
  bpro <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_retained_bpro",
    wide,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )

  expect_true(bpro$certificate$passed)
  expect_false(bpro$restart$fallback_used)
  expect_false(bpro$restart$retained_deflation)
  expect_equal(bpro$restart$retained_locked_count, 5L)
  expect_equal(bpro$restart$irlba_lbd_soft_locked_count, 0L)
  expect_equal(bpro$restart$irlba_lbd_hard_locked_count, 5L)
  expect_identical(
    bpro$restart$irlba_lbd_lock_source,
    "exact_retained_restart_certificate"
  )
  expect_true(bpro$restart$irlba_lbd_locked_triplets_certified)
  expect_true(bpro$restart$irlba_lbd_future_vectors_orthogonal_to_locks)
  expect_lte(
    bpro$restart$irlba_lbd_locked_orthogonality_loss,
    bpro$certificate$orthogonality_tolerance
  )
  expect_true(is.na(bpro$restart$irlba_lbd_lock_fallback_reason))
  expect_true(bpro$restart$irlba_lbd_bpro_policy)
  expect_equal(bpro$restart$irlba_lbd_bpro_passes_per_append, 1L)
  expect_gte(bpro$restart$irlba_lbd_bpro_monitored_appends, 32L)
  expect_gte(bpro$restart$irlba_lbd_bpro_threshold_reorthogonalizations, 0L)
  expect_lte(
    bpro$restart$irlba_lbd_bpro_max_post_append_orthogonality_loss,
    bpro$restart$irlba_lbd_bpro_monitoring_threshold
  )
  expect_lte(
    bpro$restart$irlba_lbd_bpro_basis_orthogonality_loss,
    bpro$restart$irlba_lbd_bpro_monitoring_threshold
  )
  expect_false(bpro$restart$irlba_lbd_bpro_escalation_recommended)
  expect_equal(bpro$matvecs, full$matvecs)
  expect_lt(bpro$restart$reorthogonalization_passes, full$restart$reorthogonalization_passes)
  expect_certificate_clean(bpro)
})

test_that("guarded retained IRLBA BPRO modes expose exact orthogonality guard diagnostics", {
  set.seed(702)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  one_sided <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_retained_bpro_one_sided_guarded",
    wide,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )
  block <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_retained_bpro_block_guarded",
    wide,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )

  expect_true(one_sided$certificate$passed)
  expect_identical(
    one_sided$restart$irlba_lbd_reorth_mode,
    "bpro_one_sided_guarded"
  )
  expect_true(one_sided$restart$irlba_lbd_bpro_policy)
  expect_true(xor(
    one_sided$restart$reorthogonalize_u,
    one_sided$restart$reorthogonalize_v
  ))
  expect_true(one_sided$restart$irlba_lbd_one_sided_reorth_used)
  expect_equal(one_sided$restart$irlba_lbd_bpro_block_size, 1L)
  expect_true(one_sided$restart$irlba_lbd_bpro_exact_orthogonality_passed)
  expect_lte(
    one_sided$restart$irlba_lbd_bpro_exact_orthogonality_loss,
    one_sided$certificate$orthogonality_tolerance
  )
  expect_true(is.na(one_sided$restart$irlba_lbd_bpro_guard_fallback_reason))
  expect_certificate_clean(one_sided)

  expect_true(block$certificate$passed)
  expect_identical(
    block$restart$irlba_lbd_reorth_mode,
    "bpro_block_guarded"
  )
  expect_true(block$restart$irlba_lbd_bpro_policy)
  expect_false(block$restart$irlba_lbd_one_sided_reorth_used)
  expect_true(block$restart$reorthogonalize_u)
  expect_true(block$restart$reorthogonalize_v)
  expect_equal(block$restart$irlba_lbd_bpro_block_size, 5L)
  expect_true(block$restart$irlba_lbd_bpro_exact_orthogonality_passed)
  expect_lte(
    block$restart$irlba_lbd_bpro_exact_orthogonality_loss,
    block$certificate$orthogonality_tolerance
  )
  expect_true(is.na(block$restart$irlba_lbd_bpro_guard_fallback_reason))
  expect_certificate_clean(block)
})

test_that("retained IRLBA BPRO augmented restart covers clustered and slow-decay fixtures", {
  make_fixture <- function(m, n, values, seed) {
    set.seed(seed)
    basis_rank <- length(values)
    U <- qr.Q(qr(matrix(rnorm(m * basis_rank), m, basis_rank)))
    V <- qr.Q(qr(matrix(rnorm(n * basis_rank), n, basis_rank)))
    U %*% (diag(values, basis_rank, basis_rank) %*% t(V))
  }
  fixtures <- list(
    clustered = make_fixture(
      180L, 120L,
      c(rep(10, 8), rep(9.999, 8), rep(9.99, 8), exp(-0.03 * seq_len(30L))),
      812L
    ),
    slow_decay = make_fixture(140L, 90L, exp(-0.08 * seq_len(24L)), 811L)
  )

  for (A in fixtures) {
    fit <- eigencore:::native_irlba_lbd_retained_svd(
      A,
      rank = 5L,
      work = 12L,
      retained = 7L,
      max_restarts = 7L,
      tol = 1e-8,
      vectors = "both",
      reorth_policy = "bpro_two_sided"
    )

    expect_true(fit$certificate$passed)
    expect_false(fit$restart$fallback_used)
    expect_true(fit$restart$irlba_lbd_retained_native_attempted)
    expect_identical(fit$restart$irlba_lbd_restart_state_kind, "residual_augmented_projection")
    expect_equal(fit$restart$irlba_lbd_augmented_restart_cycles, 8L)
    expect_equal(fit$restart$irlba_lbd_augmented_kept_vectors, 5L)
    expect_lte(fit$restart$irlba_lbd_augmented_small_svds, 8L)
    expect_true(any(isTRUE(fit$restart$attempt_history$certificate_passed)))
    expect_equal(
      fit$restart$irlba_lbd_augmented_cached_aq_cols,
      fit$restart$irlba_lbd_augmented_basis_cols
    )
    expect_true(fit$restart$irlba_lbd_augmented_reduces_from_scratch_work)
    expect_gt(fit$restart$irlba_lbd_augmented_matvec_savings, 0L)
    expect_true("cheap_residual" %in% names(fit$restart$attempt_history))
    expect_true(any(is.finite(fit$restart$attempt_history$cheap_residual)))
    expect_certificate_clean(fit)
  }
})

test_that("retained IRLBA residual-augmented path certifies a larger sparse wide fixture", {
  set.seed(703)
  wide <- Matrix::t(Matrix::rsparsematrix(2000L, 200L, density = 0.02))
  fit <- eigencore:::native_irlba_lbd_retained_svd(
    wide,
    rank = 5L,
    work = 12L,
    retained = 7L,
    max_restarts = 7L,
    tol = 1e-8,
    vectors = "both"
  )

  expect_true(fit$certificate$passed)
  expect_false(fit$restart$fallback_used)
  expect_identical(fit$restart$irlba_lbd_restart_state_kind, "residual_augmented_projection")
  expect_true(fit$restart$irlba_lbd_augmented_recurrence)
  expect_true(fit$restart$irlba_lbd_native_certificate_diagnostics_reused)
  expect_true(fit$restart$irlba_lbd_native_certificate_diagnostics_swapped)
  expect_lte(fit$certificate$max_backward_error, 1e-8)
  expect_gte(fit$restart$irlba_lbd_augmented_basis_cols, 40L)
  expect_certificate_clean(fit)
})

test_that("retained IRLBA reused native certificate diagnostics match direct recomputation", {
  set.seed(702)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  fit <- eigencore:::native_irlba_lbd_retained_svd(
    wide,
    rank = 5L,
    work = 12L,
    retained = 7L,
    max_restarts = 7L,
    tol = 1e-8,
    vectors = "both",
    reorth_policy = "bpro_two_sided"
  )
  direct <- eigencore:::certify_svd_operator(
    as_operator(wide),
    fit$d,
    fit$u,
    fit$v,
    tol = 1e-8
  )

  expect_true(fit$restart$irlba_lbd_native_certificate_diagnostics_reused)
  expect_true(fit$restart$irlba_lbd_native_certificate_diagnostics_swapped)
  expect_lt(max(abs(fit$certificate$residuals$left - direct$residuals$left)), 1e-10)
  expect_lt(max(abs(fit$certificate$residuals$right - direct$residuals$right)), 1e-10)
  expect_lt(max(abs(fit$certificate$backward_error - direct$backward_error)), 1e-10)
  expect_lt(max(abs(fit$certificate$orthogonality - direct$orthogonality)), 1e-10)
  expect_equal(fit$certificate$passed, direct$passed)
})

test_that("normal-scout IRLBA benchmark candidate only trusts final SVD certificate", {
  set.seed(702)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  fit <- eigencore:::run_svd_method(
    "eigencore_irlba_lbd_normal_scout",
    wide,
    rank = 5L,
    tol = 1e-8,
    seed = 702L
  )

  expect_true(fit$certificate$passed)
  expect_true(fit$restart$irlba_lbd_normal_scout_attempted)
  expect_equal(fit$restart$irlba_lbd_normal_scout_steps, "8,12,16,20")
  expect_equal(fit$restart$irlba_lbd_normal_scout_chosen_steps, 20L)
  expect_equal(fit$restart$irlba_lbd_normal_scout_count, 4L)
  expect_equal(fit$restart$irlba_lbd_normal_scout_side, "left")
  expect_false(fit$restart$irlba_lbd_normal_scout_materialized)
  expect_false(fit$restart$irlba_lbd_normal_scout_certificate_trusted)
  expect_true(fit$restart$irlba_lbd_fallback_warm_started)
  expect_gt(fit$restart$irlba_lbd_normal_scout_matvecs, 0L)
  expect_equal(
    fit$restart$irlba_lbd_normal_scout_operator_matvecs,
    2L * fit$restart$irlba_lbd_normal_scout_matvecs
  )
  expect_equal(
    tail(fit$restart$attempt_history$certificate_passed, 1L),
    TRUE
  )
  expect_true(all(!head(
    fit$restart$attempt_history$certificate_passed,
    fit$restart$irlba_lbd_normal_scout_count
  )))
  expect_certificate_clean(fit)
})

test_that("retained IRLBA LBD native wrapper supports tall CSC and scout early return", {
  tall <- rbind(
    Matrix::Diagonal(x = c(9, 6, 3, 1, 0.5, 0.2)),
    Matrix::Matrix(0, 20, 6, sparse = TRUE)
  )
  tall <- methods::as(tall, "dgCMatrix")
  fit <- eigencore:::native_irlba_lbd_retained_svd(
    tall,
    rank = 2L,
    work = 6L,
    retained = 3L,
    max_restarts = 1L,
    tol = 1e-8,
    vectors = "both"
  )

  expect_true(fit$certificate$passed)
  expect_equal(fit$d, c(9, 6), tolerance = 1e-10)
  expect_identical(fit$restart$internal_orientation, "as_given")
  expect_false(fit$restart$internal_transposed)
  expect_false(fit$restart$irlba_lbd_retained_native_attempted)
  expect_false(fit$restart$retained_restart)
  expect_false(fit$restart$retained_deflation)
  expect_equal(fit$restart$retained_locked_count, 2L)
  expect_equal(fit$restart$irlba_lbd_hard_locked_count, 2L)
  expect_identical(fit$restart$irlba_lbd_lock_source, "exact_scout_certificate")
  expect_true(fit$restart$irlba_lbd_locked_triplets_certified)
  expect_true(fit$restart$irlba_lbd_future_vectors_orthogonal_to_locks)
  expect_false(fit$restart$fallback_attempted)
  expect_equal(fit$restart$certified_attempt, 1L)
  expect_equal(fit$restart$attempted_subspaces, 6L)
  expect_true(is.data.frame(fit$restart$attempt_history))
  expect_true(fit$restart$attempt_history$certificate_passed[[1L]])
  expect_certificate_clean(fit)
})

test_that("retained IRLBA fallback warm start matches transposed orientation", {
  set.seed(706)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  fit <- eigencore:::native_irlba_lbd_retained_svd(
    wide,
    rank = 5L,
    work = 12L,
    retained = 7L,
    max_restarts = 1L,
    tol = 1e-8,
    vectors = "both",
    reorth_policy = "full_two_sided"
  )

  expect_true(fit$certificate$passed)
  expect_true(fit$restart$internal_transposed)
  expect_true(fit$restart$retained_restart)
  expect_true(fit$restart$fallback_attempted)
  expect_identical(fit$restart$irlba_lbd_lock_source, "exact_fallback_certificate")
  expect_equal(fit$restart$retained_locked_count, 5L)
  expect_equal(fit$restart$irlba_lbd_hard_locked_count, 5L)
  expect_true(fit$restart$irlba_lbd_locked_triplets_certified)
  expect_true(fit$restart$irlba_lbd_future_vectors_orthogonal_to_locks)
  expect_certificate_clean(fit)
})

test_that("retained IRLBA LBD scout state matches the native ABI orientation", {
  set.seed(705)
  wide <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))
  scout <- svd_partial(
    wide,
    rank = 5L,
    method = golub_kahan(max_subspace = 12L, reorthogonalize = FALSE),
    tol = 1e-8,
    seed = 705L
  )
  abi <- eigencore:::native_irlba_lbd_restart_abi(
    wide,
    rank = 5L,
    work = 12L,
    retained = 7L,
    max_restarts = 4L
  )
  state <- eigencore:::native_irlba_lbd_retained_state_from_scout(
    wide,
    scout,
    abi = abi
  )

  expect_s3_class(state, "eigencore_irlba_lbd_retained_state")
  expect_true(state$internal_transposed)
  expect_identical(state$internal_orientation, "transposed_wide_operator")
  expect_length(state$initial_start, 90L)
  expect_equal(dim(state$retained_right_subspace), c(90L, 7L))
  expect_equal(dim(state$retained_left_subspace), c(600L, 7L))
  expect_equal(dim(state$restart_random_tail), c(90L, 5L))
  expect_equal(length(state$alpha), 12L)
  expect_equal(length(state$beta), 12L)
  expect_equal(state$retained_from_scout, 5L)
  expect_equal(state$retained_padding, 2L)
  expect_false(state$recurrence_available)
  expect_identical(state$restart_state_kind, "ritz_subspace_only")
  expect_equal(abs(diag(crossprod(state$retained_right_subspace[, 1:5], scout$u))), rep(1, 5), tolerance = 1e-8)
  expect_equal(abs(diag(crossprod(state$retained_left_subspace[, 1:5], scout$v))), rep(1, 5), tolerance = 1e-8)

  active_wide <- Matrix::t(wide)
  out <- .Call(
    "eigencore_irlba_lbd_csc_retained",
    active_wide@i, active_wide@p, active_wide@x, as.integer(active_wide@Dim),
    state$initial_start,
    state$retained_right_subspace,
    state$retained_left_subspace,
    state$alpha,
    state$beta,
    state$restart_random_tail,
    abi$work, abi$retained, 1L, abi$rank,
    abi$target_kind, 1e-8, 1L,
    PACKAGE = "eigencore"
  )
  expect_equal(length(out$d), 5L)
  expect_true(is.data.frame(out$attempt_history))
  expect_equal(out$restart_count, 1L)
})

test_that("retained IRLBA LBD native ABI entry points are registered", {
  dense_info <- getNativeSymbolInfo(
    "eigencore_irlba_lbd_dense_retained",
    PACKAGE = "eigencore"
  )
  csc_info <- getNativeSymbolInfo(
    "eigencore_irlba_lbd_csc_retained",
    PACKAGE = "eigencore"
  )

  expect_equal(dense_info$numParameters, 14L)
  expect_equal(csc_info$numParameters, 17L)
  dense <- diag(c(6, 4, 2, 1), nrow = 6L, ncol = 4L)
  csc <- Matrix::Matrix(dense, sparse = TRUE)
  start <- c(1, 1, 0, 0) / sqrt(2)
  right <- qr.Q(qr(cbind(start, c(1, -1, 0, 0) / sqrt(2))))
  left <- qr.Q(qr(dense %*% right))
  alpha <- numeric(4L)
  beta <- numeric(4L)
  tails <- matrix(c(0, 0, 1, 0, 0, 0, 0, 1), 4L, 2L)
  dense_out <- .Call(
    "eigencore_irlba_lbd_dense_retained",
    dense, start, right, left, alpha, beta, tails,
    4L, 2L, 1L, 2L, 1L, 1e-8, 1L,
    PACKAGE = "eigencore"
  )
  csc_out <- .Call(
    "eigencore_irlba_lbd_csc_retained",
    csc@i, csc@p, csc@x, as.integer(csc@Dim),
    start, right, left, alpha, beta, tails,
    4L, 2L, 1L, 2L, 1L, 1e-8, 1L,
    PACKAGE = "eigencore"
  )
  expect_equal(dense_out$d, c(6, 4), tolerance = 1e-10)
  expect_equal(csc_out$d, c(6, 4), tolerance = 1e-10)
  expect_true(is.data.frame(dense_out$attempt_history))
  expect_true("certificate_diagnostics" %in% names(dense_out))
  expect_true(all(c(
    "left", "right", "combined", "backward_error",
    "orthogonality", "scale", "converged"
  ) %in% names(dense_out$certificate_diagnostics)))
  expect_true(all(dense_out$certificate_diagnostics$converged))
  expect_lte(nrow(dense_out$attempt_history), 2L)
  expect_true(any(isTRUE(dense_out$attempt_history$certificate_passed)))
  expect_equal(dense_out$restart_count, 1L)
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

test_that("wide sparse Gram SVD exposes opt-in explicit Gram Krylov eigensolve", {
  old_options <- options(
    eigencore.csc_left_gram_krylov_attempt = TRUE,
    eigencore.csc_left_gram_subspace_attempt = FALSE,
    eigencore.csc_left_normal_lanczos_attempt = FALSE
  )
  on.exit(options(old_options), add = TRUE)
  M <- cbind(
    Matrix::Diagonal(x = c(10, 8, 6, 1, 0.5)),
    Matrix::Matrix(0, 5, 20, sparse = TRUE)
  )

  fit <- svd_partial(M, rank = 2, target = largest(), tol = 1e-8)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$native_gram_eigensolver, "explicit_gram_krylov")
  expect_lte(fit$restart$native_gram_subspace_max_backward_error, 1e-8)
  expect_lte(fit$restart$native_gram_krylov_iterations, 5L)
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

test_that("tall sparse Gram SVD exposes implicit right-normal Lanczos", {
  old_options <- options(eigencore.csc_right_normal_lanczos_attempt = TRUE)
  on.exit(options(old_options), add = TRUE)
  M <- rbind(
    Matrix::Diagonal(x = c(10, 8, 6, 1, 0.5)),
    Matrix::Matrix(0, 20, 5, sparse = TRUE)
  )

  fit <- svd_partial(M, rank = 2, target = largest(), tol = 1e-8)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_side, "right")
  expect_identical(fit$restart$native_gram_eigensolver, "implicit_normal_lanczos")
  expect_true(fit$restart$normal_operator_implicit)
  expect_false(fit$restart$materialized_gram)
  expect_lte(fit$restart$native_implicit_normal_lanczos_max_backward_error, 1e-8)
  expect_lte(fit$restart$native_implicit_normal_lanczos_iterations, 5L)
  expect_equal(fit$d, c(10, 8), tolerance = 1e-10)
  expect_certificate_clean(fit)
})

test_that("wide sparse implicit normal Lanczos uses bounded restarted work", {
  old_options <- options(
    eigencore.csc_left_normal_lanczos_attempt = TRUE,
    eigencore.csc_left_gram_subspace_attempt = FALSE
  )
  on.exit(options(old_options), add = TRUE)
  set.seed(1906)
  M <- Matrix::t(Matrix::rsparsematrix(600L, 90L, density = 0.03))

  fit <- svd_partial(M, rank = 5L, target = largest(), tol = 1e-8, seed = 1906)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$native_gram_eigensolver, "implicit_normal_lanczos")
  expect_true(fit$restart$normal_operator_implicit)
  expect_false(fit$restart$materialized_gram)
  expect_lte(fit$restart$native_implicit_normal_lanczos_iterations, 90L)
  expect_lte(fit$restart$native_implicit_normal_lanczos_max_backward_error, 1e-8)
  expect_certificate_clean(fit)
})

test_that("tall sparse implicit right-normal Lanczos uses bounded restarted work", {
  old_options <- options(eigencore.csc_right_normal_lanczos_attempt = TRUE)
  on.exit(options(old_options), add = TRUE)
  set.seed(701)
  M <- Matrix::rsparsematrix(600L, 90L, density = 0.03)

  fit <- svd_partial(M, rank = 5L, target = largest(), tol = 1e-8, seed = 701)

  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$native_gram_eigensolver, "implicit_normal_lanczos")
  expect_true(fit$restart$normal_operator_implicit)
  expect_false(fit$restart$materialized_gram)
  expect_lte(fit$restart$native_implicit_normal_lanczos_iterations, 90L)
  expect_lte(fit$restart$native_implicit_normal_lanczos_max_backward_error, 1e-8)
  expect_certificate_clean(fit)
})

test_that("failed tall implicit right-normal candidate retries explicit Gram", {
  old_options <- options(eigencore.csc_right_normal_lanczos_attempt = TRUE)
  on.exit(options(old_options), add = TRUE)
  set.seed(1907)
  M <- Matrix::rsparsematrix(600L, 90L, density = 0.03)

  fit <- svd_partial(M, rank = 5L, target = largest(), tol = 1e-8, seed = 1907)

  expect_identical(fit$method, "native certified Gram SVD special case")
  expect_identical(fit$restart$native_gram_eigensolver, "lapack_dsyevr")
  expect_true(fit$restart$explicit_gram_retry_used)
  expect_gt(fit$restart$native_implicit_normal_lanczos_iterations, 0L)
  expect_false(fit$restart$fallback_attempted)
  expect_false(fit$restart$fallback_used)
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
    expect_true(fit$restart$native_sketch)
    expect_equal(fit$restart$projection_kind, "native_direct_qt_a")
    expect_true(fit$restart$projection_transposed)
    if (identical(normalizer, "qr")) {
      expect_true(fit$restart$controller_native)
      expect_true(fit$restart$native_certificate_diagnostics)
      expect_false(fit$restart$certificate_reuses_projection)
    } else {
      expect_false(isTRUE(fit$restart$controller_native))
      expect_true(fit$restart$certificate_reuses_projection)
    }
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

test_that("native randomized SVD certificate matches direct recomputation", {
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

  expect_true(fit$restart$controller_native)
  expect_true(fit$restart$native_certificate_diagnostics)
  expect_false(fit$restart$certificate_reuses_projection)
  expect_equal(fit$certificate$residuals$left, direct$residuals$left, tolerance = 1e-10)
  expect_equal(fit$certificate$residuals$right, direct$residuals$right, tolerance = 1e-10)
  expect_equal(fit$certificate$backward_error, direct$backward_error, tolerance = 1e-10)
  expect_equal(fit$certificate$passed, direct$passed)
})

test_that("randomized SVD wide-core eigensolve matches dense SVD values", {
  set.seed(1906)
  core <- matrix(rnorm(14 * 80), nrow = 14, ncol = 80)
  fast <- eigencore:::randomized_svd_core_decomposition(core, rank = 6)
  dense <- svd(core, nu = 6, nv = 6)

  expect_identical(fast$solver, "native_left_gram_eigen_selected")
  expect_equal(fast$d, dense$d[seq_len(6)], tolerance = 1e-10)
  expect_equal(crossprod(fast$u), diag(6), tolerance = 1e-10)
  expect_equal(crossprod(fast$v), diag(6), tolerance = 1e-10)
  expect_equal(core %*% fast$v, sweep(fast$u, 2L, fast$d, `*`), tolerance = 1e-10)
})

test_that("SVD planner records inspectable method controls", {
  old_options <- options(
    eigencore.gram_svd_max_dimension = 512,
    eigencore.gram_svd_memory_mb = 64,
    eigencore.gram_svd_work_budget = Inf,
    eigencore.gram_svd_rank_fraction_limit = 0.5,
    eigencore.gram_svd_min_aspect_ratio = 2
  )
  on.exit(options(old_options), add = TRUE)

  set.seed(513)
  M <- Matrix::rsparsematrix(120, 30, density = 0.05)

  gram_plan <- plan_solver(svd_problem(M), rank = 4)
  expect_identical(gram_plan$method, "native certified Gram SVD special case")
  expect_identical(gram_plan$controls$gram_side, "right")
  expect_identical(gram_plan$controls$gram_dimension, 30L)
  expect_identical(gram_plan$controls$gram_max_dimension, 512L)
  expect_equal(gram_plan$controls$gram_memory_budget_mb, 64)
  expect_equal(gram_plan$controls$rank_fraction_limit, 0.5)
  expect_equal(gram_plan$controls$min_aspect_ratio, 2)
  expect_true(gram_plan$controls$gram_policy_passed)
  expect_identical(
    gram_plan$controls$materialization_policy,
    "budgeted smaller Gram materialization"
  )
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

  one_sided_wide_plan <- plan_solver(
    svd_problem(Matrix::rsparsematrix(90, 600, density = 0.01)),
    rank = 5,
    method = golub_kahan(reorthogonalize = FALSE)
  )
  expect_false(one_sided_wide_plan$controls$reorthogonalize)
  expect_identical(one_sided_wide_plan$controls$initial_max_subspace, 45L)

  one_sided_rank10_plan <- plan_solver(
    svd_problem(Matrix::rsparsematrix(3000, 800, density = 0.002)),
    rank = 10,
    method = auto()
  )
  expect_identical(
    one_sided_rank10_plan$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )

  one_sided_rank10_gk_plan <- plan_solver(
    svd_problem(Matrix::rsparsematrix(3000, 800, density = 0.002)),
    rank = 10,
    method = golub_kahan(reorthogonalize = FALSE)
  )
  expect_identical(one_sided_rank10_gk_plan$method, "native prototype Golub-Kahan")
  expect_false(one_sided_rank10_gk_plan$controls$reorthogonalize)
  expect_identical(one_sided_rank10_gk_plan$controls$initial_max_subspace, 90L)

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

test_that("memory-budgeted Gram SVD policy keeps default cutoff and supports opt-in just-over-512 route", {
  old_options <- options(
    eigencore.gram_svd_max_dimension = 512,
    eigencore.gram_svd_memory_mb = 64,
    eigencore.gram_svd_work_budget = Inf,
    eigencore.gram_svd_rank_fraction_limit = 0.5,
    eigencore.gram_svd_min_aspect_ratio = 2
  )
  on.exit(options(old_options), add = TRUE)

  set.seed(514)
  M <- Matrix::rsparsematrix(1200, 600, density = 0.004)

  default_plan <- plan_solver(svd_problem(M), rank = 3)
  expect_identical(
    default_plan$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )
  default_policy <- eigencore:::gram_svd_policy(dim(M), rank = 3, target = largest())
  expect_false(default_policy$eligible)
  expect_true("gram_dimension_exceeds_max" %in% default_policy$rejection_reasons)

  options(eigencore.gram_svd_max_dimension = 768)
  opt_in_plan <- plan_solver(svd_problem(M), rank = 3)
  expect_identical(opt_in_plan$method, "native certified Gram SVD special case")
  expect_identical(opt_in_plan$controls$gram_dimension, 600L)
  expect_identical(opt_in_plan$controls$gram_max_dimension, 768L)
  expect_lte(
    opt_in_plan$controls$estimated_total_materialization_bytes,
    opt_in_plan$controls$gram_memory_budget_bytes
  )

  fit <- svd_partial(M, rank = 3, target = largest(), tol = 1e-8, seed = 514)
  expect_identical(fit$restart$kind, "gram_svd_special_case")
  expect_identical(fit$restart$gram_dimension, 600L)
  expect_true(fit$restart$materialized_gram)
  expect_true(fit$restart$certified_in_original_coordinates)
  expect_certificate_clean(fit)
})

test_that("Gram SVD policy rejects large small sides by memory budget even when dimension cap allows them", {
  old_options <- options(
    eigencore.gram_svd_max_dimension = 2048,
    eigencore.gram_svd_memory_mb = 1,
    eigencore.gram_svd_work_budget = Inf,
    eigencore.gram_svd_rank_fraction_limit = 0.5,
    eigencore.gram_svd_min_aspect_ratio = 2
  )
  on.exit(options(old_options), add = TRUE)

  set.seed(515)
  M <- Matrix::rsparsematrix(1400, 600, density = 0.002)
  policy <- eigencore:::gram_svd_policy(dim(M), rank = 3, target = largest())

  expect_false(policy$eligible)
  expect_true("gram_memory_budget_exceeded" %in% policy$rejection_reasons)
  expect_gt(policy$estimated_total_materialization_bytes, policy$gram_memory_budget_bytes)

  plan <- plan_solver(svd_problem(M), rank = 3)
  expect_identical(
    plan$method,
    "native certified implicit Gram SVD (thick-restart Lanczos)"
  )
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
  expect_true(fit$restart$native_sketch)
  expect_true(fit$restart$controller_native)
  expect_true(fit$restart$sparse_native_controller)
  expect_true(fit$restart$native_certificate_diagnostics)
  expect_equal(fit$restart$projection_kind, "native_direct_qt_a")
  expect_true(fit$restart$projection_transposed)
  expect_false(fit$restart$certificate_reuses_projection)
  expect_match(fit$restart$certificate_policy, "stochastic sketch is not sufficient")
  expect_false(fit$restart$refine)
  expect_false(fit$restart$refinement_attempted)
  expect_true(is.logical(fit$restart$initial_certificate_passed))
  expect_true(is.finite(fit$restart$initial_max_backward_error))

  pair <- eigencore:::randomized_svd_apply_pair(as_operator(M))
  set.seed(817)
  fused_sketch <- pair$sketch(6L)
  set.seed(817)
  Omega <- matrix(rnorm(ncol(M) * 6L), ncol(M), 6L)
  expect_equal(fused_sketch, as.matrix(M %*% Omega), tolerance = 1e-12)
  sketch <- pair$apply(Omega)
  expect_equal(sketch, as.matrix(M %*% Omega), tolerance = 1e-12)
  adjoint_sketch <- pair$apply_adjoint(sketch)
  expect_equal(adjoint_sketch, as.matrix(Matrix::crossprod(M, sketch)), tolerance = 1e-12)

  Q <- qr.Q(qr(as.matrix(M %*% Omega)))
  projected <- pair$project(Q)
  adjoint <- pair$apply_adjoint(Q)
  expect_true(isTRUE(attr(projected, "transposed", exact = TRUE)))
  attr(projected, "transposed") <- NULL
  expect_equal(projected, t(adjoint), tolerance = 1e-12)
})

test_that("randomized SVD dense fused sketch preserves seeded R Gaussian contract", {
  set.seed(818)
  M <- matrix(rnorm(30 * 18), nrow = 30, ncol = 18)
  pair <- eigencore:::randomized_svd_apply_pair(as_operator(M))

  set.seed(819)
  fused_sketch <- pair$sketch(7L)
  set.seed(819)
  Omega <- matrix(rnorm(ncol(M) * 7L), ncol(M), 7L)

  expect_equal(pair$sketch_kind, "native_fused_a_omega")
  expect_equal(fused_sketch, M %*% Omega, tolerance = 1e-12)
})
