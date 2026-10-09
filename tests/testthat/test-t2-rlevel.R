# Tranche 2, R-level performance items (docs/review-2026-10.md P6, P7, P9,
# P13 shift-invert part, P15). Each test pins the contract the optimisation
# must preserve as well as the cheaper behaviour.

large_hash_counter <- function() {
  env <- new.env()
  env$count <- 0L
  original <- eigencore:::stable_raw_hash
  env$fn <- function(x) {
    if (as.numeric(utils::object.size(x)) > 1e5) {
      env$count <- env$count + 1L
    }
    original(x)
  }
  env
}

test_that("P6: built-in operator identity is deferred, cached, and unchanged in value", {
  set.seed(1)
  A <- crossprod(matrix(rnorm(200 * 200), 200)) / 200
  counter <- large_hash_counter()
  testthat::local_mocked_bindings(stable_raw_hash = counter$fn, .package = "eigencore")

  op <- as_operator(A)
  expect_identical(counter$count, 0L)
  expect_true(isTRUE(attr(op$identity, "deferred")))

  id1 <- operator_identity(op)
  expect_identical(counter$count, 1L)
  id2 <- operator_identity(op)
  expect_identical(counter$count, 1L)
  expect_identical(id1, id2)
  expect_null(attributes(id1)[["deferred"]])

  # Same value as an eager digest of the same parts.
  eager <- eigencore:::builtin_identity_from_parts(
    dim(A), "double", hermitian(), op$metadata
  )
  expect_identical(id1, eager)
  expect_identical(id1, operator_identity(as_operator(A)))

  # A whole solve hashes the source at most once (plan freeze); validation
  # and the plan's operator reuse the cached digest.
  counter$count <- 0L
  fit <- eig_partial(A, k = 3)
  expect_lte(counter$count, 1L)
  expect_true(fit$certificate$passed)
  expect_identical(operator_identity(fit), list(A = id1))
})

test_that("P6: an operator modified after construction is re-hashed, plans still fail closed", {
  A <- diag(c(5, 4, 3, 2, 1))
  op <- as_operator(A)
  original <- operator_identity(op)
  changed <- op
  changed$metadata$source <- diag(c(5, 4, 3, 2, 1.5))
  expect_false(identical(operator_identity(changed), original))
  expect_identical(
    operator_identity(changed),
    operator_identity(as_operator(diag(c(5, 4, 3, 2, 1.5))))
  )

  plan <- plan_solver(eigen_problem(A), k = 2L)
  plan$problem$A$metadata$source <- diag(c(9, 4, 3, 2, 1))
  err <- tryCatch(solve(plan), error = identity)
  expect_s3_class(err, "eigencore_plan_error")
  expect_identical(err$code, "operator_incompatible")
})

test_that("P6: identity survives serialisation and solving leaves the plan bytes unchanged", {
  set.seed(2)
  S <- Matrix::rsparsematrix(60, 60, density = 0.1)
  S <- S + Matrix::t(S) + Matrix::Diagonal(60, x = 1:60)
  S <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
  plan <- plan_solver(eigen_problem(S), k = 3L)
  before <- serialize(plan, NULL, version = 3L)
  fit <- solve(plan)
  expect_true(fit$certificate$passed)
  expect_identical(serialize(plan, NULL, version = 3L), before)
  restored <- unserialize(before)
  expect_identical(operator_identity(restored), plan$operator_identity)
  expect_true(solve(restored)$certificate$passed)
})

test_that("P9: non-dgC Matrix classes are converted once to native storage", {
  set.seed(3)
  S <- Matrix::rsparsematrix(80, 80, density = 0.08)
  S <- S + Matrix::t(S) + Matrix::Diagonal(80, x = seq(1, 8, length.out = 80))
  dgC <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
  reference <- values(eig_partial(dgC, k = 3))

  inputs <- list(
    dgT = methods::as(dgC, "TsparseMatrix"),
    dgR = methods::as(dgC, "RsparseMatrix"),
    dsT = methods::as(methods::as(dgC, "symmetricMatrix"), "TsparseMatrix"),
    dge = methods::as(dgC, "unpackedMatrix"),
    dsy = methods::as(methods::as(dgC, "symmetricMatrix"), "unpackedMatrix")
  )
  for (name in names(inputs)) {
    x <- inputs[[name]]
    op <- as_operator(x)
    expect_true(isTRUE(op$metadata$native), info = name)
    expect_identical(op$structure$kind, "hermitian", info = name)
    if (inherits(x, "sparseMatrix")) {
      expect_identical(op$metadata$storage, "dgCMatrix", info = name)
      expect_identical(op$metadata$input_storage, class(x)[[1L]], info = name)
    } else {
      expect_true(is.matrix(op$metadata$source), info = name)
    }
    fit <- eig_partial(x, k = 3)
    expect_false(grepl("reference", fit$method), info = name)
    expect_true(fit$certificate$passed, info = name)
    expect_equal(values(fit), reference, tolerance = 1e-8, info = name)
  }

  # Unit-triangular sparse storage materialises its implicit diagonal.
  tri <- Matrix::sparseMatrix(i = c(1, 2, 3, 2, 3), j = c(1, 2, 3, 1, 2),
                              x = c(1, 1, 1, 0.5, 0.25),
                              dims = c(3, 3), triangular = TRUE)
  tri <- Matrix::.diagN2U(tri)
  expect_true(inherits(tri, "dtCMatrix"))
  expect_identical(tri@diag, "U")
  op_tri <- as_operator(tri)
  X <- diag(3)
  expect_equal(op_tri$apply(X), as.matrix(tri) %*% X)
  expect_equal(op_tri$apply_adjoint(X), t(as.matrix(tri)) %*% X)
  expect_equal(unname(as.matrix(tri))[cbind(1:3, 1:3)], c(1, 1, 1))
})

test_that("P15: real-part targets on Hermitian problems use the algebraic native routes", {
  set.seed(4)
  S <- Matrix::rsparsematrix(300, 300, density = 0.02)
  S <- S + Matrix::t(S) + Matrix::Diagonal(300, x = seq_len(300) / 30)
  S <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
  lr <- eig_partial(S, k = 4, target = largest_real())
  la <- eig_partial(S, k = 4, target = largest())
  expect_identical(lr$method, la$method)
  expect_false(grepl("reference", lr$method))
  expect_identical(lr$target, "largest")
  expect_equal(values(lr), values(la), tolerance = 1e-8)
  expect_true(lr$certificate$passed)

  sr <- eig_partial(S, k = 4, target = smallest_real())
  sa <- eig_partial(S, k = 4, target = smallest())
  expect_identical(sr$method, sa$method)
  expect_equal(values(sr), values(sa), tolerance = 1e-8)

  # General structure keeps the real-part target.
  P <- eigen_problem(matrix(c(1, 2, 0, 3), 2), target = largest_real())
  expect_identical(P$target$kind, "largest_real")
})

test_that("P7: reference Lanczos certifies once instead of every iteration", {
  set.seed(5)
  d <- c(seq(10, 6, length.out = 4), seq(1, 0.01, length.out = 196))
  calls <- 0L
  op <- linear_operator(
    dim = c(200, 200),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      alpha * (d * X)
    },
    structure = hermitian()
  )
  fit <- eig_partial(op, k = 3, method = lanczos(), seed = 5)
  expect_match(fit$method, "reference Hermitian Lanczos")
  # The two-norm lower-bound scale (C12) lets a matrix-free certificate pass.
  expect_true(fit$certificate$passed)
  expect_true(all(fit$certificate$converged))
  expect_equal(values(fit), d[1:3], tolerance = 1e-8)
  # Lanczos steps + a single k-column certificate: the dominant Ritz values
  # already give the norm bound, so no extra probes are needed. The
  # per-iteration certificate used to add k + 8 applies for every step j >= k.
  # The target-completeness probe (C50) adds its own, separately recorded,
  # block applies after the solve.
  probe_calls <- fit$certificate$completeness$operator_block_calls
  expect_identical(fit$certificate$target_completeness, "probed")
  expect_lte(calls, fit$iterations + 2L + probe_calls)
  expect_lte(as.integer(fit$work$certification_operator_columns), 3L)
})

test_that("P7/C12: the Krylov two-norm bound is memoised per operator", {
  calls <- 0L
  op <- linear_operator(
    dim = c(50, 50),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      alpha * (seq_len(50) * X)
    },
    structure = hermitian()
  )
  # The pre-solve value needs no applies and no randomness.
  info <- eigencore:::operator_norm_for_certificate_info(op)
  expect_identical(calls, 0L)
  expect_false(info$scale_is_estimate)
  expect_identical(info$norm_bound_type, "two_norm_lower_bound")

  first <- eigencore:::two_norm_krylov_bound(op)
  used <- calls
  second <- eigencore:::two_norm_krylov_bound(op)
  expect_gt(used, 0L)
  expect_identical(calls, used)
  expect_identical(first, second)
  expect_lte(first$value, 50 * (1 + 1e-12))
  expect_gt(first$value, 45)

  # A copy whose apply closure was replaced does not reuse the memo.
  other <- op
  other$apply <- function(X, alpha = 1, beta = 0, Y = NULL) alpha * (2 * X)
  third <- eigencore:::two_norm_krylov_bound(other)
  expect_equal(third$value, 2, tolerance = 1e-12)
})

test_that("P7: sparse shift-invert runs the native thick-restart callback and certifies", {
  set.seed(6)
  n <- 400L
  S <- Matrix::rsparsematrix(n, n, density = 0.01)
  S <- S + Matrix::t(S) + Matrix::Diagonal(n, x = seq_len(n))
  S <- methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
  sigma <- 200.3
  fit <- eig_partial(S, k = 6, target = nearest(sigma))
  expect_identical(fit$method, "native thick-restart Hermitian Lanczos shift-invert (sparse LU solve callback)")
  expect_identical(fit$restart$kind, "native_thick_restart_shift_invert_callback")
  expect_true(fit$certificate$passed)
  dense <- eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values
  expected <- dense[order(abs(dense - sigma))][1:6]
  expect_equal(sort(values(fit)), sort(expected), tolerance = 1e-8)

  # A subspace that cannot hold all wanted pairs at once is restarted rather
  # than returned unconverged.
  tight <- eig_partial(S, k = 6, target = nearest(sigma),
                       method = auto(max_subspace = 8L))
  expect_equal(tight$plan$controls$max_subspace, 8L)
  expect_true(tight$certificate$passed)
  expect_equal(sort(values(tight)), sort(expected), tolerance = 1e-8)

  rs <- eigs_sym(S, 4, sigma = sigma)
  expect_equal(sort(rs$values), sort(expected[1:4]), tolerance = 1e-8)
})

test_that("P13: dense shift-invert factors once and keeps its singularity gates", {
  set.seed(7)
  A <- crossprod(matrix(rnorm(40 * 40), 40))
  prep <- eigencore:::shift_invert_solver_dense(A, sigma = 3.3)
  M <- A - 3.3 * diag(40)
  expect_equal(prep$M, M)
  X <- matrix(rnorm(80), 40, 2)
  expect_equal(prep$solve_fn(X), solve(M, X), tolerance = 1e-10)
  expect_true(is.finite(prep$cache$condition_estimate))
  expect_gt(prep$cache$condition_estimate, 0)
  expect_identical(prep$cache$condition_estimate_type, "dense_qr_triangular_rcond")

  B <- diag(seq(1, 2, length.out = 40))
  gprep <- eigencore:::shift_invert_solver_dense(A, sigma = 3.3, B = B)
  expect_equal(gprep$M, A - 3.3 * B)

  expect_error(
    eigencore:::shift_invert_solver_dense(diag(c(-1e-10, 0, 1e-10)), sigma = 0),
    "singular|near-singular|rank-deficient"
  )
})
