test_that("clustered Hermitian eigenvalues are certified by residuals and subspaces", {
  values <- c(10, 10 - 1e-9, 10 - 2e-9, 2, 1, 0.5)
  A <- symmetric_with_spectrum(values, seed = 31)
  fit <- eig_partial(A, k = 3, target = largest(), tol = 1e-8)
  oracle <- eigen(A, symmetric = TRUE)

  expect_equal(values(fit), oracle$values[1:3], tolerance = 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, 1:3]), 1e-6)
  expect_certificate_clean(fit)
})

test_that("nearly repeated singular values are judged by subspace accuracy", {
  s <- c(7, 7 - 1e-9, 7 - 2e-9, 2, 0.5)
  A <- rectangular_with_singular_values(s, m = 9, n = 6, seed = 32)
  fit <- svd_partial(A, rank = 3, target = largest(), tol = 1e-8)
  oracle <- svd(A, nu = 3, nv = 3)

  expect_equal(values(fit), oracle$d[1:3], tolerance = 1e-8)
  expect_lt(subspace_distance(left_vectors(fit), oracle$u[, 1:3]), 1e-6)
  expect_lt(subspace_distance(right_vectors(fit), oracle$v[, 1:3]), 1e-6)
  expect_certificate_clean(fit)
})

test_that("rank-deficient rectangular SVD returns finite certified triplets", {
  A <- rectangular_with_singular_values(c(6, 3, 0, 0), m = 8, n = 5, seed = 33)
  fit <- svd_partial(A, rank = 4, target = largest(), tol = 1e-8)

  expect_equal(values(fit)[1:2], c(6, 3), tolerance = 1e-10)
  expect_true(all(is.finite(values(fit))))
  expect_false(any(is.nan(values(fit))))
  expect_certificate_clean(fit)
})

test_that("graph Laplacian nullspace is recovered without sparse densification", {
  n <- 8L
  A <- Matrix::bandSparse(n, k = c(-1, 0, 1), diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1)))
  fit <- eig_partial(A, k = 1, target = smallest(), method = lanczos(max_subspace = n), seed = 34, tol = 1e-8)

  expect_equal(fit$plan$method, "native scalar thick-restart Hermitian Lanczos")
  expect_lt(abs(values(fit)), 1e-8)
  expect_lt(subspace_distance(vectors(fit), matrix(rep(1 / sqrt(n), n), ncol = 1)), 1e-5)
  expect_certificate_clean(fit)
})

test_that("native LOBPCG certifies clustered smallest dense eigenpairs", {
  values <- c(0, 1e-8, 2e-8, 1, 2, 4, 8, 16)
  A <- symmetric_with_spectrum(values, seed = 71)
  fit <- eig_partial(
    A,
    k = 3,
    target = smallest(),
    method = lobpcg(maxit = 120L),
    seed = 71,
    tol = 1e-8
  )
  oracle <- eigen(A, symmetric = TRUE)
  idx <- order(oracle$values)[1:3]

  expect_equal(fit$method, "native standard Hermitian LOBPCG prototype")
  expect_lt(max(abs(values(fit) - sort(oracle$values)[1:3])), 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, idx]), 1e-6)
  expect_certificate_clean(fit)
})

test_that("native preconditioned LOBPCG handles Laplacian near-nullspace clusters", {
  n <- 50L
  A <- Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1), c(1, rep(2, n - 2), 1), rep(-1, n - 1))
  )
  preconditioner <- shifted_tridiagonal_preconditioner(A, shift = 1e-4)
  fit <- eig_partial(
    A,
    k = 4,
    target = smallest(),
    method = lobpcg(maxit = 80L, preconditioner = preconditioner),
    seed = 72,
    tol = 1e-8
  )
  oracle <- eigen(as.matrix(A), symmetric = TRUE)
  idx <- order(oracle$values)[1:4]

  expect_equal(fit$method, "native standard Hermitian LOBPCG prototype")
  expect_true(fit$restart$preconditioner_native)
  expect_lt(max(abs(values(fit) - sort(oracle$values)[1:4])), 1e-8)
  expect_lt(subspace_distance(vectors(fit), oracle$vectors[, idx]), 1e-5)
  expect_certificate_clean(fit)
})

test_that("ill-conditioned generalized SPD problems remain B-orthonormal", {
  A <- diag(c(9, 6, 4, 2, 1))
  B <- diag(c(1, 1e-2, 1e-4, 1e-6, 1e-8))
  fit <- eig_partial(A, B = B, k = 3, target = smallest(), tol = 1e-8)

  expect_true(certificate(fit)$passed)
  expect_equal(crossprod(vectors(fit), B %*% vectors(fit)), diag(3), tolerance = 1e-8)
  expect_lt(max(abs(A %*% vectors(fit) - B %*% vectors(fit) %*% diag(values(fit)))), 1e-7)
})

test_that("non-normal nonsymmetric matrices with real spectra are certified", {
  V <- matrix(c(
    1, 1, 1,
    0, 1e-3, 1,
    0, 0, 1
  ), nrow = 3, byrow = TRUE)
  A <- V %*% diag(c(5, 3, 1)) %*% solve(V)
  fit <- eig_partial(A, k = 2, target = largest_magnitude(), tol = 1e-8)

  expect_equal(values(fit), c(5, 3), tolerance = 1e-8)
  expect_certificate_clean(fit)
})

test_that("clustered nonsymmetric spectra return certified biorthogonal left vectors", {
  A <- rbind(
    c(5, 1e-6, 0, 0),
    c(0, 5 - 1e-5, 2e-6, 0),
    c(0, 0, 4, 1e-6),
    c(0, 0, 0, -1)
  )

  fit <- eig_partial(A, k = 2, target = largest_real(), tol = 1e-8)

  expect_equal(values(fit), c(5, 5 - 1e-5), tolerance = 1e-8)
  expect_true(fit$certificate$passed)
  expect_false(is.null(left_vectors(fit)))
  expect_false(is.null(right_vectors(fit)))
  expect_true(fit$left_certificate$passed)
  expect_lt(fit$left_certificate$max_backward_error, 1e-8)
  expect_lt(max(abs(fit$biorthogonality - diag(2L))), 1e-8)
  expect_match(fit$warnings, "left residuals and biorthogonality certified")
})

test_that("defective nonsymmetric spectra do not report false biorthogonal success", {
  A <- matrix(c(
    3, 1, 0,
    0, 3, 1,
    0, 0, 1
  ), nrow = 3, byrow = TRUE)

  fit <- eig_partial(A, k = 2, target = largest_real(), tol = 1e-10)

  expect_equal(values(fit), c(3, 3), tolerance = 1e-8)
  expect_true(fit$certificate$passed)
  expect_true(fit$left_eigenvectors$supported)
  expect_false(fit$left_certificate$passed)
  expect_gte(fit$left_certificate$max_orthogonality_loss, 0.5)
  expect_match(fit$warnings, "left residuals or biorthogonality did not pass certificate")
  expect_false(grepl("left residuals and biorthogonality certified", fit$warnings, fixed = TRUE))
})

test_that("dense nonsymmetric native Arnoldi certifies complex right residuals", {
  A <- rbind(
    c(0, -2, 0),
    c(2,  0, 0),
    c(0,  0, 1)
  )
  fit <- eig_partial(A, k = 2, target = largest_imaginary(), tol = 1e-8)

  expect_equal(Im(values(fit))[[1L]], 2, tolerance = 1e-10)
  expect_true(is.complex(values(fit)))
  expect_true(is.complex(vectors(fit)))
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$certificate$passed)
  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_true(fit$restart$native)
  expect_true(fit$restart$ritz_extraction_native)
  expect_equal(fit$restart$extraction, "refined_ritz")
  expect_true(fit$restart$refined_extraction_native)
  expect_false(fit$restart$krylov_schur)
  expect_match(fit$warnings, "right residuals certified")
})

test_that("native matrix-free Arnoldi handles nonsymmetric real spectra", {
  A <- rbind(
    c(5, 2, 0, 0),
    c(0, 3, 1, 0),
    c(0, 0, 1, 0.5),
    c(0, 0, 0, -2)
  )
  op <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    },
    structure = general(),
    name = "matrix_free_nonnormal",
    metadata = list(frobenius_norm = sqrt(sum(A^2)))
  )

  fit <- eig_partial(op, k = 2L, target = largest_real(), tol = 1e-10, seed = 11)

  expect_equal(fit$plan$method, eigencore:::native_matrix_free_arnoldi_label())
  expect_equal(fit$plan$controls$arnoldi_extraction, "projected_ritz")
  expect_match(fit$warnings, "native matrix-free Arnoldi callback cycle", fixed = TRUE)
  expect_true(fit$certificate$passed)
  expect_equal(sort(Re(fit$values), decreasing = TRUE), c(5, 3), tolerance = 1e-8)
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$restart$implemented)
  expect_true(fit$restart$native)
  expect_true(fit$restart$matrix_free)
  expect_true(fit$restart$ritz_extraction_native)
  expect_equal(fit$restart$extraction, "projected_ritz")
  expect_false(fit$restart$refined_extraction_native)
  expect_error(
    eigencore:::native_arnoldi_general(
      op,
      k = 2L,
      target = largest_real(),
      extraction = "refined_ritz"
    ),
    "matrix-free refined Ritz extraction is future scope"
  )
  expect_certificate_clean(fit)
})

test_that("native Arnoldi handles nonsymmetric sparse real spectra without densifying", {
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3, 4),
    j = c(1, 2, 2, 3, 3, 4),
    x = c(5, 2, 3, 1, 1, -2),
    dims = c(4, 4)
  )

  fit <- eig_partial(
    A,
    k = 2L,
    target = largest_real(),
    tol = 1e-10,
    seed = 11,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_match(fit$warnings, "native Arnoldi cycle", fixed = TRUE)
  expect_equal(fit$plan$controls$max_restarts, 5L)
  expect_equal(fit$plan$controls$arnoldi_extraction, "refined_ritz")
  expect_true(fit$plan$controls$refined_extraction_native)
  expect_true(fit$certificate$passed)
  expect_equal(sort(Re(fit$values), decreasing = TRUE), c(5, 3), tolerance = 1e-8)
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$restart$implemented)
  expect_true(fit$restart$native)
  expect_equal(fit$restart$max_restarts, 5L)
  expect_true(fit$restart$ritz_extraction_native)
  expect_equal(fit$restart$extraction, "refined_ritz")
  expect_true(fit$restart$refined_extraction_native)
  expect_false(fit$restart$krylov_schur)
  expect_equal(fit$restart$v2_issue, "bd-01KTF6H41S9XDN286TR3V184P4")
  expect_true(all(c("cycle", "ritz_extraction") %in% names(fit$restart$stage_seconds)))
  expect_true(all(is.finite(fit$restart$stage_seconds)))
  expect_true(all(fit$restart$stage_seconds >= 0))
  expect_true(all(c("cycle_seconds", "ritz_extraction_seconds") %in%
                    names(fit$restart$attempt_history)))
  expect_certificate_clean(fit)
})

test_that("native Arnoldi handles nonsymmetric dense real spectra without oracle fallback", {
  A <- diag(c(5, 3, 1, -2))
  A[1, 2] <- 2
  A[2, 3] <- 1

  fit <- eig_partial(
    A,
    k = 2L,
    target = largest_real(),
    tol = 1e-10,
    seed = 11,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_match(fit$warnings, "native Arnoldi cycle", fixed = TRUE)
  expect_true(fit$certificate$passed)
  expect_equal(sort(Re(fit$values), decreasing = TRUE), c(5, 3), tolerance = 1e-8)
  expect_equal(fit$certificate$certificate_type, "right_residual_backward_error")
  expect_false(fit$certificate$orthogonality_required)
  expect_true(fit$restart$implemented)
  expect_true(fit$restart$native)
  expect_true(fit$restart$ritz_extraction_native)
  expect_equal(fit$restart$extraction, "refined_ritz")
  expect_true(fit$restart$refined_extraction_native)
  expect_certificate_clean(fit)
})

test_that("native Arnoldi uses full dense subspace for dense compatibility", {
  n <- 40L
  A <- matrix(0, n, n)
  A[1:2, 1:2] <- matrix(c(0, -2, 2, 0), 2, byrow = TRUE)
  A[3:4, 3:4] <- matrix(c(0, -1, 1, 0), 2, byrow = TRUE)
  diag(A)[5:n] <- seq(0.2, 3.7, length.out = n - 4L)

  fit <- eig_partial(
    A,
    k = 4L,
    target = largest_imaginary(),
    tol = 1e-10,
    seed = 15046,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_equal(fit$plan$controls$max_subspace, n)
  expect_equal(fit$plan$controls$arnoldi_extraction, "refined_ritz")
  expect_equal(fit$restart$max_subspace, n)
  expect_true(fit$certificate$passed)
  expect_equal(sum(fit$certificate$converged), 4L)
  expect_true(fit$restart$native)
  expect_true(fit$restart$ritz_extraction_native)
  expect_true(fit$restart$refined_extraction_native)
  expect_certificate_clean(fit)
})

test_that("native Arnoldi Ritz extraction preserves complex Ritz pairs", {
  A <- matrix(c(
    0, -2, 0, 0,
    2,  0, 0, 0,
    0,  0, 3, 1,
    0,  0, 0, 1
  ), 4, 4, byrow = TRUE)
  cycle <- eigencore:::native_arnoldi_cycle(
    eigencore:::as_operator(A),
    start = rep(0.5, 4),
    m = 4L
  )
  ritz <- eigencore:::native_arnoldi_ritz(
    eigencore:::as_operator(A),
    cycle,
    k = 2L,
    target = largest_imaginary(),
    tol = 1e-10
  )

  expect_true(is.complex(ritz$values))
  expect_true(is.complex(ritz$vectors))
  expect_equal(sort(Im(ritz$values), decreasing = TRUE)[1L], 2, tolerance = 1e-8)
  residual <- A %*% ritz$vectors - ritz$vectors %*% diag(ritz$values)
  expect_lt(max(Mod(residual)), 1e-8)
})

test_that("native refined Ritz extraction improves clustered non-normal residuals", {
  set.seed(42)
  n <- 12L
  Q <- matrix(rnorm(n * n), n)
  A <- Q %*% diag(c(5, 5 - 1e-5, 4.99, seq(2, -2, length.out = n - 3L))) %*% solve(Q)
  op <- eigencore:::as_operator(A)
  start <- rnorm(n)
  start <- start / sqrt(sum(start^2))
  cycle <- eigencore:::native_arnoldi_cycle(op, start, 10L)

  projected <- eigencore:::native_arnoldi_ritz(
    op, cycle, k = 1L, target = largest_real(), tol = 1e-12,
    extraction = "projected_ritz"
  )
  refined <- eigencore:::native_arnoldi_ritz(
    op, cycle, k = 1L, target = largest_real(), tol = 1e-12,
    extraction = "refined_ritz"
  )

  expect_equal(refined$extraction, "refined_ritz")
  expect_true(length(refined$refined_residual_estimates) == 1L)
  expect_lt(refined$certificate$max_residual, projected$certificate$max_residual / 5)
  expect_lt(refined$certificate$max_backward_error,
            projected$certificate$max_backward_error / 5)
  expect_equal(refined$certificate$max_residual,
               refined$refined_residual_estimates,
               tolerance = 1e-10)
})

test_that("native Arnoldi sparse path supports imaginary-part targets", {
  A <- Matrix::bdiag(
    matrix(c(0, -3, 3, 0), 2, 2, byrow = TRUE),
    matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE),
    Matrix::Diagonal(3, c(0.5, 0.25, -0.5))
  )
  A <- methods::as(A, "dgCMatrix")

  fit <- eig_partial(
    A,
    k = 2L,
    target = largest_imaginary(),
    tol = 1e-10,
    seed = 15091,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_equal(Im(fit$values), c(3, 1), tolerance = 1e-8)
  expect_true(is.complex(fit$vectors))
  expect_true(fit$restart$native)
  expect_true(fit$restart$ritz_extraction_native)
  expect_true(fit$restart$refined_extraction_native)
  expect_certificate_clean(fit)
})

test_that("native Arnoldi restart budget is wired and keeps best attempt", {
  old_options <- options(eigencore.arnoldi_max_restarts = 2L)
  on.exit(options(old_options), add = TRUE)
  n <- 30L
  A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n - 1L)),
    j = c(seq_len(n), seq_len(n - 1L) + 1L),
    x = c(seq(n, 1), rep(2, n - 1L)),
    dims = c(n, n)
  )

  fit <- eig_partial(
    A,
    k = 5L,
    target = largest_real(),
    tol = 1e-12,
    maxit = 8L,
    seed = 1,
    allow_dense_fallback = "never"
  )

  history <- fit$restart$attempt_history
  finite_errors <- history$max_backward_error[is.finite(history$max_backward_error)]

  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_equal(fit$plan$controls$max_restarts, 2L)
  expect_equal(fit$plan$controls$arnoldi_extraction, "refined_ritz")
  expect_equal(fit$restart$max_restarts, 2L)
  expect_equal(fit$restart$restart_count, 2L)
  expect_equal(nrow(history), 3L)
  expect_equal(
    fit$restart$selected_attempt,
    history$attempt[which.min(history$max_backward_error)]
  )
  expect_equal(fit$certificate$max_backward_error, min(finite_errors), tolerance = 1e-12)
  expect_false(fit$certificate$passed)
  expect_true(fit$restart$native)
  expect_true(fit$restart$refined_extraction_native)
})

test_that("native Arnoldi default subspace certifies sparse benchmark-sized row", {
  n <- 80L
  A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n - 1L)),
    j = c(seq_len(n), seq_len(n - 1L) + 1L),
    x = c(seq(n, 1), rep(2, n - 1L)),
    dims = c(n, n)
  )

  fit <- eig_partial(
    A,
    k = 8L,
    target = largest_real(),
    tol = 1e-10,
    seed = 15092,
    allow_dense_fallback = "never"
  )

  expect_equal(fit$plan$method, eigencore:::native_refined_arnoldi_label())
  expect_equal(fit$restart$max_subspace, 72L)
  expect_equal(fit$restart$max_restarts, 5L)
  expect_true(fit$certificate$passed)
  expect_equal(sum(fit$certificate$converged), 8L)
  expect_true(fit$restart$native)
  expect_true(fit$restart$refined_extraction_native)
  expect_certificate_clean(fit)
})

test_that("poorly scaled SVD inputs remain finite and certified", {
  A <- diag(c(1e6, 1e2, 1, 1e-2, 1e-6))
  A <- diag(c(1e-3, 1, 1e3, 1e-2, 1e2)) %*% A
  fit <- svd_partial(A, rank = 3, target = largest(), tol = 1e-8)
  oracle <- svd(A, nu = 3, nv = 3)

  expect_equal(values(fit), oracle$d[1:3], tolerance = 1e-8)
  expect_true(all(is.finite(certificate(fit)$backward_error)))
  expect_certificate_clean(fit)
})

test_that("near-singular shift-invert requests fail loudly", {
  A <- diag(c(-1e-10, 0, 1e-10))
  expect_error(
    eig_partial(A, k = 1, target = nearest(0), method = shift_invert(0)),
    "singular|near-singular"
  )
})
