# Regression tests for the native-code assurance findings F4-F7
# (docs/test-assurance.md, "Native code assurance", Findings).

assurance_path_plus_diagonal <- function(n = 60L) {
  methods::as(
    Matrix::bandSparse(n, k = c(0, 1),
                       diagonals = list(rep(2, n), rep(-1, n - 1)),
                       symmetric = TRUE) +
      Matrix::Diagonal(n, seq_len(n) / n),
    "generalMatrix"
  )
}

assurance_callback_op <- function(A) {
  linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * as.matrix(A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * as.matrix(A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian()
  )
}

expect_verified_certificate <- function(cert) {
  expect_true(isTRUE(cert$passed))
  expect_true(cert$target_completeness %in%
                eigencore:::verified_completeness_states())
}

# RSpectra "BE": k values alternating from both ends of the spectrum, the
# extra one from the high end when k is odd; returned in decreasing order.
assurance_be_reference <- function(spectrum, k) {
  s <- sort(spectrum)
  n <- length(s)
  high <- k - k %/% 2L
  low <- k %/% 2L
  sort(c(s[seq_len(low)], s[n - seq_len(high) + 1L]), decreasing = TRUE)
}

test_that("F4: eigs_sym(sparse, which = 'BE') takes the native both-ends route", {
  A <- assurance_path_plus_diagonal()
  spectrum <- eigen(as.matrix(A), symmetric = TRUE, only.values = TRUE)$values
  for (k in 1:6) {
    res <- eigs_sym(A, k = k, which = "BE")
    expect_identical(res$diagnostics$method,
                     eigencore:::native_both_ends_lanczos_label())
    expect_identical(res$nconv, k)
    expect_verified_certificate(res$certificate)
    expect_equal(res$values, assurance_be_reference(spectrum, k),
                 tolerance = 1e-8)
    dense <- eigs_sym(as.matrix(A), k = k, which = "BE")
    expect_equal(res$values, dense$values, tolerance = 1e-8)
    if (requireNamespace("RSpectra", quietly = TRUE)) {
      expect_equal(res$values, RSpectra::eigs_sym(A, k = k, which = "BE")$values,
                   tolerance = 1e-8)
    }
  }
})

test_that("F4: 'BE' with opts$initvec and with a function input stays native", {
  A <- assurance_path_plus_diagonal()
  n <- nrow(A)
  spectrum <- eigen(as.matrix(A), symmetric = TRUE, only.values = TRUE)$values
  expected <- assurance_be_reference(spectrum, 3L)
  # initvec used to force the unrestarted reference Lanczos: 0 of 3 pairs.
  res <- expect_silent(eigs_sym(A, k = 3, which = "BE",
                                opts = list(initvec = rep(1, n))))
  expect_identical(res$diagnostics$method,
                   eigencore:::native_both_ends_lanczos_label())
  expect_identical(res$diagnostics$start_source, "user_supplied")
  expect_verified_certificate(res$certificate)
  expect_equal(res$values, expected, tolerance = 1e-8)

  f <- function(x, args) as.numeric(args %*% x)
  res <- eigs_sym(f, k = 3, which = "BE", n = n, args = A,
                  opts = list(initvec = rep(1, n)))
  expect_identical(res$diagnostics$method,
                   eigencore:::native_both_ends_lanczos_label())
  expect_identical(res$nconv, 3L)
  expect_equal(res$values, expected, tolerance = 1e-8)
})

test_that("F5: explicit lanczos() on a matrix-free operator uses the native kernel", {
  A <- assurance_path_plus_diagonal()
  op <- assurance_callback_op(A)
  truth <- sort(eigen(as.matrix(A), symmetric = TRUE, only.values = TRUE)$values,
                decreasing = TRUE)
  for (fit in list(
    eig_partial(op, 3L, method = lanczos(), seed = 1),
    eig_partial(op, 3L, method = lanczos(max_subspace = 30), seed = 1),
    eig_partial(op, 3L, method = lanczos(), maxit = 500, seed = 1)
  )) {
    expect_identical(fit$method, eigencore:::native_matrix_free_block_lanczos_label())
    expect_identical(fit$nconv, 3L)
    expect_verified_certificate(certificate(fit))
    expect_equal(values(fit), truth[1:3], tolerance = 1e-8)
  }
  low <- eig_partial(op, 3L, target = smallest(), method = lanczos(), seed = 2)
  expect_identical(low$method, eigencore:::native_matrix_free_block_lanczos_label())
  expect_verified_certificate(certificate(low))
  # A warm start is consumed by the native kernel too.
  warm <- eig_partial(op, 3L, method = lanczos(), seed = 3,
                      initial_subspace = rep(1, nrow(A)))
  expect_identical(warm$method, eigencore:::native_matrix_free_block_lanczos_label())
  expect_identical(warm$start_source, "user_supplied")
  expect_verified_certificate(certificate(warm))
})

test_that("F5: maxit sets the step budget of the unrestarted reference Lanczos", {
  A <- assurance_path_plus_diagonal()
  op <- assurance_callback_op(A)
  # nearest() has no native matrix-free kernel: the reference route remains.
  P <- eigen_problem(op, target = nearest(2))
  default <- plan_solver(P, k = 3L, method = lanczos())
  expect_identical(default$method,
                   "reference Hermitian Lanczos (prototype/oracle fallback)")
  expect_identical(default$controls$max_subspace, 29L)
  raised <- plan_solver(P, k = 3L, method = lanczos(), maxit = 60L)
  expect_identical(raised$controls$max_subspace, 60L)
  expect_identical(raised$controls$iteration_limit_kind, "lanczos_steps")
  # An explicit max_subspace is still capped, never raised, by maxit.
  capped <- plan_solver(P, k = 3L, method = lanczos(max_subspace = 20L), maxit = 50L)
  expect_identical(capped$controls$max_subspace, 20L)
  fit <- eig_partial(op, 3L, target = nearest(2), method = lanczos(), maxit = 60L,
                     seed = 4)
  expect_identical(fit$nconv, 3L)
  expect_verified_certificate(certificate(fit))
})

test_that("F6: a singular shifted_tridiagonal_preconditioner fails with a clear error", {
  n <- 60L
  L <- Matrix::bandSparse(n, k = c(0, 1),
                          diagonals = list(c(1, rep(2, n - 2), 1), rep(-1, n - 1)),
                          symmetric = TRUE)
  expect_error(shifted_tridiagonal_preconditioner(L, shift = 0),
               "shifted_tridiagonal_preconditioner\\(\\).*singular")
  expect_s3_class(shifted_tridiagonal_preconditioner(L, shift = 1e-3),
                  "eigencore_preconditioner")

  # A factor that turns singular after construction is named by the native
  # LOBPCG instead of the bare "status=-5".
  P <- shifted_tridiagonal_preconditioner(L, shift = 1e-3)
  meta <- attr(P, "eigencore_preconditioner")
  meta$diag <- meta$diag - 1e-3
  attr(P, "eigencore_preconditioner") <- meta
  A <- L + Matrix::Diagonal(n, seq_len(n) / n)
  expect_error(
    eig_partial(A, 3L, target = smallest(), method = lobpcg(preconditioner = P)),
    "shifted_tridiagonal_preconditioner\\(\\).*singular"
  )
})

test_that("F7: sparse identity hashes the matrix content, not derived moments", {
  S <- Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 2, 3, 3),
                            x = c(1, 1 / 3, 2.5, 0.1), dims = c(3, 3))
  op <- as_operator(S)
  id <- operator_identity(op)$revision
  # Pinned: the digest is a function of the exact slots only, so it is the
  # same on every platform; a change here needs a new identity_hash_format().
  expect_identical(id, "651bddc5ec4b559be41d83837863c170")
  # Column moments are accumulated in long double; a platform with another
  # long double width (or valgrind) may round them differently. Simulate
  # that: one ulp in every derived moment leaves the identity unchanged.
  ulp <- function(x) x * (1 + .Machine$double.eps)
  other <- op
  for (key in c("column_sums", "column_sum_squares", "column_means",
                "column_centered_sum_squares", "frobenius_norm")) {
    other$metadata[[key]] <- ulp(other$metadata[[key]])
  }
  expect_identical(operator_identity(other)$revision, id)
  # The Matrix factorization cache is not part of the identity ...
  cached <- S
  cached@factors <- list(dummy = Matrix::Diagonal(3))
  expect_identical(operator_identity(as_operator(cached))$revision, id)
  # ... but every content slot is.
  moved <- S
  moved@x[2] <- moved@x[2] * (1 + .Machine$double.eps)
  expect_false(identical(operator_identity(as_operator(moved))$revision, id))
  named <- S
  dimnames(named) <- list(c("a", "b", "c"), NULL)
  expect_false(identical(operator_identity(as_operator(named))$revision, id))
  sym <- Matrix::forceSymmetric(S + Matrix::t(S))
  expect_false(identical(
    operator_identity(as_operator(sym))$revision,
    operator_identity(as_operator(methods::as(sym, "generalMatrix")))$revision
  ))
})

test_that("F7: plans persisted under identity format v2 fail with a typed plan error", {
  A <- assurance_path_plus_diagonal(30L)
  plan <- plan_solver(eigen_problem(A), k = 2L)
  expect_identical(plan$serialization$hash_format, "eigencore-identity-hash-v3")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(plan, path)
  restored <- readRDS(path)
  expect_true(isTRUE(certificate(solve(restored))$passed))
  # A plan saved before F7 recorded the v2 format (and a v2 operator digest).
  old <- restored
  old$serialization$hash_format <- "eigencore-identity-hash-v2"
  err <- tryCatch(solve(old), error = identity)
  expect_s3_class(err, "eigencore_plan_error")
  expect_identical(err$code, "identity_format_changed")
  expect_match(conditionMessage(err), "Re-plan")

  state <- restart_state(solve(plan))
  state$serialization$hash_format <- "eigencore-identity-hash-v2"
  err <- tryCatch(restart_state(state), error = identity)
  expect_s3_class(err, "eigencore_restart_state_error")
  expect_identical(err$code, "identity_format_changed")
})
