# Tranche 4b: native composed operators (C52) and R wrappers that let the
# native block applies allocate their own output (P11).

t4b_dense_of <- function(op) {
  op$apply(diag(op$dim[[2L]]))
}

t4b_expect_matches <- function(op, reference, tol = 1e-10) {
  set.seed(11)
  n <- op$dim[[2L]]
  m <- op$dim[[1L]]
  X <- matrix(stats::rnorm(n * 3L), n, 3L)
  Z <- matrix(stats::rnorm(m * 3L), m, 3L)
  Y <- matrix(stats::rnorm(m * 3L), m, 3L)
  scale <- max(1, max(abs(reference)))
  expect_equal(op$apply(X), reference %*% X, tolerance = tol * scale)
  expect_equal(op$apply(X[, 1L, drop = FALSE]), reference %*% X[, 1L, drop = FALSE],
               tolerance = tol * scale)
  expect_equal(op$apply(X, alpha = 2, beta = -0.5, Y = Y),
               2 * reference %*% X - 0.5 * Y, tolerance = tol * scale)
  if (!is.null(op$apply_adjoint)) {
    expect_equal(op$apply_adjoint(Z), crossprod(reference, Z),
                 tolerance = tol * scale)
    W <- matrix(stats::rnorm(n * 3L), n, 3L)
    expect_equal(op$apply_adjoint(Z, alpha = -1, beta = 3, Y = W),
                 -crossprod(reference, Z) + 3 * W, tolerance = tol * scale)
  }
}

t4b_sparse <- function(m, n, density = 0.05, seed = 1) {
  set.seed(seed)
  Matrix::rsparsematrix(m, n, density = density)
}

test_that("lazy compose of large dense factors gets a native composite kernel", {
  old <- options(eigencore.native_composite = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(1)
  # The product has more entries than both factors, so compose() keeps it
  # lazy (P10 policy).
  A <- matrix(stats::rnorm(400 * 40), 400, 40)
  B <- matrix(stats::rnorm(40 * 450), 40, 450)
  op <- compose(A, B)
  expect_null(op$metadata$source)
  expect_identical(op$metadata$storage, "native_composite")
  expect_true(eigencore:::has_native_composite_kernel(op))
  expect_false(is.null(attr(op$apply, "eigencore_native_kernel")))
  expect_false(is.null(attr(op$apply_adjoint, "eigencore_native_kernel")))
  expect_true(is.na(eigencore:::native_kernel_kind(op)))
  t4b_expect_matches(op, A %*% B)
})

test_that("composite kernels cover sparse, diagonal, scaled, summed and centered leaves", {
  S <- t4b_sparse(60, 40, seed = 2)
  D <- matrix(stats::rnorm(60 * 40), 60, 40)
  w_rows <- stats::runif(60) + 0.5
  w_cols <- stats::runif(40) + 0.5
  Sd <- as.matrix(S)

  # Sum of a lazy row scaling of a composite and a scaled composite.
  big <- matrix(stats::rnorm(40 * 300), 40, 300)
  inner <- compose(as_operator(S), compose(big, t(big)))
  expect_identical(inner$metadata$storage, "native_composite")
  ref_inner <- Sd %*% big %*% t(big)
  t4b_expect_matches(inner, ref_inner)

  rows <- scale_rows(inner, w_rows)
  expect_identical(rows$metadata$storage, "native_composite")
  t4b_expect_matches(rows, w_rows * ref_inner)

  cols <- scale_cols(inner, w_cols)
  t4b_expect_matches(cols, sweep(ref_inner, 2L, w_cols, `*`))

  scaled <- eigencore:::operator_scale(inner, -2.5)
  t4b_expect_matches(scaled, -2.5 * ref_inner)

  summed <- eigencore:::operator_sum(rows, as_operator(S), as_operator(D))
  expect_identical(summed$metadata$storage, "native_composite")
  t4b_expect_matches(summed, w_rows * ref_inner + Sd + D)

  # Centering of a composite (column, row, and double centering).
  for (cfg in list(c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))) {
    cm <- colMeans(ref_inner)
    rm <- rowMeans(ref_inner)
    centered <- center(inner, rows = cfg[[1L]], columns = cfg[[2L]],
                       row_means = rm, col_means = cm)
    expect_identical(centered$metadata$storage, "native_composite")
    ref <- ref_inner
    if (cfg[[2L]]) ref <- sweep(ref, 2L, cm, `-`)
    if (cfg[[1L]]) ref <- sweep(ref, 1L, rowMeans(ref), `-`)
    t4b_expect_matches(centered, ref)
  }

  # Diagonal leaves and adjoint views.
  Dg <- Matrix::Diagonal(x = w_rows)
  diag_compose <- compose(compose(as_operator(Dg), as_operator(S)), big)
  t4b_expect_matches(diag_compose, w_rows * (Sd %*% big))
  adj <- adjoint(inner)
  expect_identical(adj$metadata$storage, "native_composite")
  t4b_expect_matches(adj, t(ref_inner))
  t4b_expect_matches(adjoint(adj), ref_inner)
})

test_that("crossprod_operator of a sparse matrix is a native composite", {
  S <- t4b_sparse(2000, 400, density = 0.05, seed = 3)
  op <- crossprod_operator(as_operator(S))
  expect_null(op$metadata$materialized_crossprod)
  expect_identical(op$metadata$storage, "native_composite")
  expect_identical(op$structure$kind, "hermitian")
  t4b_expect_matches(op, as.matrix(Matrix::crossprod(S)))
})

test_that("centered and centered-scaled CSC leaves compose natively", {
  S <- t4b_sparse(80, 30, density = 0.2, seed = 4)
  Sd <- as.matrix(S)
  cm <- colMeans(Sd)
  w <- stats::runif(30) + 0.5
  centered <- center(as_operator(S))
  expect_identical(centered$metadata$storage, "centered_dgCMatrix")
  # The fused centered CSC operator also carries a kernel for matrix-free
  # native solvers, without changing its storage label.
  expect_false(is.null(attr(centered$apply, "eigencore_native_kernel")))
  cs <- scale_cols(centered, w)
  expect_identical(cs$metadata$storage, "centered_scaled_dgCMatrix")
  ref_cs <- sweep(sweep(Sd, 2L, cm, `-`), 2L, w, `*`)
  big <- matrix(stats::rnorm(30 * 300), 30, 300)
  op <- compose(cs, compose(big, t(big)))
  expect_identical(op$metadata$storage, "native_composite")
  t4b_expect_matches(op, ref_cs %*% big %*% t(big))
  op2 <- compose(compose(big, t(big)), adjoint(centered))
  expect_identical(op2$metadata$storage, "native_composite")
  t4b_expect_matches(op2, big %*% t(big) %*% t(sweep(Sd, 2L, cm, `-`)))
  # adjoint() views of plain CSC leaves.
  L <- matrix(stats::rnorm(80 * 2), 80, 2)
  R <- matrix(stats::rnorm(2 * 500), 2, 500)
  op3 <- compose(compose(adjoint(as_operator(S)), L), R)
  expect_identical(op3$metadata$storage, "native_composite")
  t4b_expect_matches(op3, t(Sd) %*% L %*% R)
})

test_that("callback operands keep the R composition", {
  base <- matrix(stats::rnorm(40 * 30), 40, 30)
  cb <- linear_operator(
    dim = dim(base),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (base %*% X)
      if (!is.null(Y) && beta != 0) out + beta * Y else out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * crossprod(base, X)
      if (!is.null(Y) && beta != 0) out + beta * Y else out
    }
  )
  big <- matrix(stats::rnorm(30 * 300), 30, 300)
  op <- compose(cb, compose(big, t(big)))
  expect_null(op$metadata$storage)
  expect_null(attr(op$apply, "eigencore_native_kernel"))
  t4b_expect_matches(op, base %*% big %*% t(big))
})

test_that("the composite kernel can be disabled and survives serialisation", {
  set.seed(5)
  A <- matrix(stats::rnorm(400 * 40), 400, 40)
  B <- matrix(stats::rnorm(40 * 450), 40, 450)
  old <- options(eigencore.native_composite = FALSE)
  plain <- compose(A, B)
  options(old)
  expect_null(plain$metadata$storage)
  expect_null(attr(plain$apply, "eigencore_native_kernel"))

  op <- compose(A, B)
  restored <- unserialize(serialize(op, NULL))
  expect_identical(restored$metadata$storage, "native_composite")
  X <- matrix(stats::rnorm(450 * 2), 450, 2)
  Z <- matrix(stats::rnorm(400 * 2), 400, 2)
  expect_equal(restored$apply(X), A %*% (B %*% X), tolerance = 1e-10)
  expect_equal(restored$apply_adjoint(Z), crossprod(A %*% B, Z),
               tolerance = 1e-10)
  # A restored operator also drives a native solver.
  fit <- svd_partial(restored, rank = 2, seed = 1)
  expect_equal(fit$d, svd(A %*% B, nu = 0, nv = 0)$d[1:2], tolerance = 1e-7)
})

test_that("native solvers drive composites without R callbacks and match explicit results", {
  set.seed(6)
  A <- matrix(stats::rnorm(500 * 60), 500, 60)
  B <- matrix(stats::rnorm(60 * 450), 60, 450)
  op <- compose(A, B)
  expect_identical(op$metadata$storage, "native_composite")
  fit <- svd_partial(op, rank = 4, seed = 1)
  ref <- svd(A %*% B, nu = 0, nv = 0)$d[1:4]
  expect_equal(fit$d, ref, tolerance = 1e-7)
  expect_true(isTRUE(certificate(fit)$passed))
  # Natively dispatched applies are still booked in the typed work record.
  w <- work(fit)
  expect_gt(w$operator_block_calls, 0L)
  expect_gt(w$adjoint_block_calls, 0L)

  # Wrap the R apply closures with counters but keep the kernel attribute:
  # the native solve loop dispatches on the attribute and never evaluates
  # them.
  r_calls <- 0L
  wrap <- function(f) {
    g <- function(X, alpha = 1, beta = 0, Y = NULL) {
      r_calls <<- r_calls + 1L
      f(X, alpha = alpha, beta = beta, Y = Y)
    }
    attributes(g) <- attributes(f)
    g
  }
  probe <- op
  probe$apply <- wrap(op$apply)
  probe$apply_adjoint <- wrap(op$apply_adjoint)
  fit_probe <- svd_partial(probe, rank = 4, seed = 1)
  expect_equal(fit_probe$d, fit$d, tolerance = 1e-12)
  wp <- work(fit_probe)
  expect_equal(wp$operator_block_calls, w$operator_block_calls)
  expect_lt(r_calls, wp$operator_block_calls + wp$adjoint_block_calls)
})

test_that("native composite applies inside the matrix-free native kernels", {
  S <- t4b_sparse(3000, 300, density = 0.03, seed = 7)
  op <- crossprod_operator(as_operator(S))
  expect_identical(op$metadata$storage, "native_composite")
  fit <- eig_partial(op, k = 4, method = lanczos(block = 2L), seed = 2)
  ref <- svd(as.matrix(S), nu = 0, nv = 0)$d[1:4]^2
  expect_equal(sort(values(fit), decreasing = TRUE), ref, tolerance = 1e-7)
  expect_true(isTRUE(certificate(fit)$passed))

  w_rows <- stats::runif(3000) + 0.5
  sum_op <- eigencore:::operator_sum(
    scale_rows(compose(as_operator(S), diag(300) + 0), w_rows),
    eigencore:::operator_scale(as_operator(S), 2)
  )
  expect_identical(sum_op$metadata$storage, "native_composite")
  ref_mat <- w_rows * as.matrix(S) + 2 * as.matrix(S)
  fit2 <- svd_partial(sum_op, rank = 3, seed = 3)
  expect_equal(fit2$d, svd(ref_mat, nu = 0, nv = 0)$d[1:3], tolerance = 1e-7)
})

test_that("R-level native applies accept a NULL output when beta is zero (P11)", {
  set.seed(8)
  A <- matrix(stats::rnorm(20 * 10), 20, 10)
  S <- t4b_sparse(20, 10, density = 0.3, seed = 9)
  X <- matrix(stats::rnorm(30), 10, 3)
  Y <- matrix(stats::rnorm(60), 20, 3)
  expect_null(eigencore:::block_apply_y(NULL, 1))
  expect_null(eigencore:::block_apply_y(Y, 0))
  named <- Y
  dimnames(named) <- list(NULL, c("a", "b", "c"))
  expect_identical(eigencore:::block_apply_y(named, 0), named)
  expect_identical(eigencore:::block_apply_y(Y, 2), Y)

  expect_equal(eigencore:::dense_block_apply(A, X), A %*% X)
  expect_equal(eigencore:::dense_block_apply(A, X, alpha = 2, beta = 0, Y = Y),
               2 * A %*% X)
  expect_equal(eigencore:::dense_block_apply(A, X, alpha = 2, beta = 1, Y = Y),
               2 * A %*% X + Y)
  out <- eigencore:::dense_block_apply(A, X, beta = 0, Y = named)
  expect_identical(colnames(out), c("a", "b", "c"))
  expect_equal(eigencore:::csc_block_apply(S, X), as.matrix(S %*% X))
  expect_equal(eigencore:::csc_block_apply(S, Y, transpose = TRUE, beta = 0.5,
                                           Y = X),
               as.matrix(Matrix::crossprod(S, Y)) + 0.5 * X)
  Dg <- Matrix::Diagonal(x = 1:20 + 0)
  expect_equal(eigencore:::diagonal_block_apply(Dg, Y), (1:20) * Y)
  Ac <- A + 1i * A
  expect_equal(eigencore:::complex_dense_block_apply(Ac, X), Ac %*% X)
  cm <- colMeans(as.matrix(S))
  expect_equal(
    eigencore:::csc_centered_block_apply(S, X, col_means = cm),
    sweep(as.matrix(S), 2L, cm, `-`) %*% X
  )
  expect_equal(
    eigencore:::csc_centered_scaled_block_apply(S, cm, rep(2, 10), X),
    sweep(as.matrix(S), 2L, cm, `-`) %*% (2 * X)
  )
})
