# Complex Hermitian eigenproblems on the iterative routes
# (R/complex_hermitian.R): realified native routes, reference complex
# Lanczos / LOBPCG, complex_operator() sparse input, complex pencils, and the
# completeness verdicts they reach.

complex_hermitian_fixture <- function(n, lambda = NULL, seed = 1L) {
  set.seed(seed)
  Z <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n)
  if (is.null(lambda)) {
    return((Z + Conj(t(Z))) / 2)
  }
  Q <- qr.Q(qr(Z))
  H <- Q %*% (lambda * Conj(t(Q)))
  (H + Conj(t(H))) / 2
}

# Random-phase (disordered magnetic) Laplacian on an m x m grid: a sparse
# complex Hermitian matrix given by its real and imaginary parts.
complex_hermitian_grid <- function(m, seed = 1L) {
  set.seed(seed)
  n <- m * m
  g <- expand.grid(i = seq_len(m), j = seq_len(m))
  id <- function(i, j) (j - 1L) * m + i
  h <- g[g$i < m, ]
  v <- g[g$j < m, ]
  I <- c(id(h$i, h$j), id(v$i, v$j))
  J <- c(id(h$i + 1L, h$j), id(v$i, v$j + 1L))
  th <- runif(length(I), 0, 2 * pi)
  deg <- tabulate(c(I, J), n)
  re <- Matrix::sparseMatrix(i = c(I, J, seq_len(n)), j = c(J, I, seq_len(n)),
                             x = c(-cos(th), -cos(th), deg), dims = c(n, n))
  im <- Matrix::sparseMatrix(i = c(I, J), j = c(J, I),
                             x = c(-sin(th), sin(th)), dims = c(n, n))
  list(re = re, im = im, H = as.matrix(re) + 1i * as.matrix(im))
}

expect_complex_pairs <- function(fit, H, truth, B = NULL, tol = 1e-8,
                                 completeness = c("inertia_verified", "probed", "exact")) {
  cert <- certificate(fit)
  expect_true(isTRUE(cert$passed), info = paste(cert$notes, collapse = "; "))
  expect_true(cert$target_completeness %in% completeness,
              info = cert$target_completeness)
  expect_equal(sort(fit$values), sort(truth), tolerance = 1e-9)
  V <- fit$vectors
  expect_true(is.complex(V))
  BV <- if (is.null(B)) V else B %*% V
  R <- H %*% V - sweep(BV, 2L, fit$values, `*`)
  scale <- max(abs(truth)) + 1
  expect_lt(max(sqrt(colSums(Mod(R)^2))) / scale, 10 * tol)
  G <- Conj(t(V)) %*% BV
  expect_lt(max(Mod(G - diag(ncol(V)))), 1e-6)
}

test_that("realification maps vectors, operators and spectra consistently", {
  ns <- asNamespace("eigencore")
  H <- complex_hermitian_fixture(7)
  R <- ns$complex_hermitian_realify_dense(H)
  expect_true(isSymmetric(R))
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  er <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(er, rep(ev, each = 2L), tolerance = 1e-12)
  V <- eigen(H, symmetric = TRUE)$vectors[, 1:2]
  U <- ns$complex_hermitian_realify_vectors(V)
  expect_equal(dim(U), c(14L, 4L))
  expect_equal(crossprod(U), diag(4), tolerance = 1e-12)
  expect_equal(R %*% U, U %*% diag(rep(ev[1:2], each = 2L)), tolerance = 1e-10)
  expect_equal(ns$complex_hermitian_complexify(U[, c(1, 3)], 7L), V, tolerance = 1e-14)
  # Matrix-free embedding agrees with the explicit one.
  op <- linear_operator(dim = c(7, 7), apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    alpha * (H %*% X)
  }, dtype = "complex", structure = hermitian())
  Rop <- ns$complex_hermitian_realify_operator(op)
  X <- matrix(rnorm(14 * 3), 14)
  expect_equal(ns$apply_operator(Rop, X), R %*% X, tolerance = 1e-12)
})

test_that("dense complex Hermitian: lanczos(), lobpcg() and auto() match eigen()", {
  H <- complex_hermitian_fixture(60)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  for (method in list(lanczos(), lanczos(block = 2L), lobpcg())) {
    fit <- eig_partial(H, k = 4, method = method, seed = 1)
    expect_complex_pairs(fit, H, ev[1:4])
    fit <- eig_partial(H, k = 3, target = smallest(), method = method, seed = 1)
    expect_complex_pairs(fit, H, sort(ev)[1:3])
  }
  # The lanczos() route is the realified native block kernel.
  fit <- eig_partial(H, k = 4, method = lanczos(), seed = 1)
  expect_identical(fit$method, "native realified complex Hermitian route (real 2n embedding)")
  expect_match(fit$restart$inner_method, "native block Hermitian Lanczos")
  expect_identical(fit$restart$complex_rank, 4L)
  # Small dense auto() keeps zheev, for every target.
  fit <- eig_partial(H, k = 4, target = nearest(0.3))
  expect_identical(fit$method, "native dense complex Hermitian LAPACK fallback")
  expect_identical(certificate(fit)$target_completeness, "exact")
  expect_equal(sort(fit$values), sort(ev[order(abs(ev - 0.3))][1:4]), tolerance = 1e-10)
})

test_that("auto() routes larger dense complex Hermitian problems to the realified kernel", {
  H <- complex_hermitian_fixture(300, seed = 2L)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  fit <- eig_partial(H, k = 6, seed = 3)
  expect_identical(fit$method, "native realified complex Hermitian route (real 2n embedding)")
  expect_complex_pairs(fit, H, ev[1:6])
  plan <- plan_solver(eigen_problem(H), k = 6)
  expect_identical(plan$controls$embedding, "real_2n")
  expect_identical(plan$controls$inner_k, 12L)
  expect_match(plan$controls$inner_method, "lanczos\\(block = 2\\)")
  cases <- list(
    list(target = smallest(), truth = sort(ev)[1:6]),
    list(target = largest_magnitude(), truth = ev[order(-abs(ev))][1:6]),
    list(target = smallest_magnitude(), truth = ev[order(abs(ev))][1:6]),
    list(target = nearest(1), truth = ev[order(abs(ev - 1))][1:6]),
    list(target = both_ends(3, 3), truth = c(sort(ev)[1:3], ev[1:3]))
  )
  for (case in cases) {
    fit <- eig_partial(H, k = 6, target = case$target, seed = 3)
    expect_complex_pairs(fit, H, case$truth)
  }
  # shift_invert(sigma) on the embedding (real sigma).
  fit <- eig_partial(H, k = 4, method = shift_invert(1), seed = 3)
  expect_complex_pairs(fit, H, ev[order(abs(ev - 1))][1:4])
  expect_error(eig_partial(H, k = 4, method = shift_invert(1 + 1i)), "real")
})

test_that("repeated complex eigenvalues are returned with their multiplicity", {
  lambda <- c(rep(5, 3), rep(4, 2), seq(-3, 3, length.out = 195))
  H <- complex_hermitian_fixture(200, lambda, seed = 4L)
  for (method in list(auto(), lanczos(), lobpcg())) {
    fit <- eig_partial(H, k = 5, method = method, seed = 5)
    expect_complex_pairs(fit, H, c(5, 5, 5, 4, 4))
  }
  # Straddling the edge: k = 4 cuts the double eigenvalue 4.
  fit <- eig_partial(H, k = 4, seed = 5)
  expect_complex_pairs(fit, H, c(5, 5, 5, 4))
})

test_that("reference complex Lanczos is fixed and verified on the embedding", {
  H <- complex_hermitian_fixture(150, seed = 6L)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  old <- options(eigencore.complex_hermitian_realify = FALSE)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(H, k = 4, method = lanczos(max_subspace = 120), seed = 7)
  expect_identical(fit$method, "reference Hermitian Lanczos (prototype/oracle fallback)")
  expect_complex_pairs(fit, H, ev[1:4])
  # Direct kernel call: complex basis, real tridiagonal.
  out <- asNamespace("eigencore")$reference_lanczos_hermitian(
    as_operator(H), k = 3, maxit = 100, tol = 1e-10
  )
  expect_true(is.complex(out$vectors))
  expect_equal(out$values, ev[1:3], tolerance = 1e-9)
  # A scalar complex Krylov run misses copies of a repeated eigenvalue; the
  # completeness check on the embedding repairs or flags it, never certifies
  # a wrong set.
  lambda <- c(rep(5, 3), seq(-3, 3, length.out = 97))
  Hr <- complex_hermitian_fixture(100, lambda, seed = 8L)
  fit <- eig_partial(Hr, k = 3, method = lanczos(max_subspace = 60), seed = 9)
  cert <- certificate(fit)
  if (isTRUE(cert$passed)) {
    expect_equal(fit$values, c(5, 5, 5), tolerance = 1e-8)
    expect_true(cert$target_completeness %in% c("inertia_verified", "probed", "repaired"))
  } else {
    expect_false(cert$target_completeness %in% c("inertia_verified", "probed", "exact"))
  }
})

test_that("complex matrix-free Hermitian operators solve with verified completeness", {
  H <- complex_hermitian_fixture(200, seed = 10L)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  op <- linear_operator(dim = c(200, 200), apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (H %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  }, dtype = "complex", structure = hermitian())
  fit <- eig_partial(op, k = 5, seed = 11)
  expect_complex_pairs(fit, H, ev[1:5])
  # Below the materialisation limit the count is exact.
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
  fit <- eig_partial(op, k = 5, target = smallest(), method = lanczos(), seed = 11)
  expect_complex_pairs(fit, H, sort(ev)[1:5])
  fit <- eig_partial(op, k = 4, method = lobpcg(), seed = 11)
  expect_identical(fit$method, "reference LOBPCG prototype")
  expect_complex_pairs(fit, H, ev[1:4])
  # Without the materialised count the probe on the embedding decides.
  old <- options(eigencore.completeness_materialize_limit = 10L)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(op, k = 5, seed = 11)
  expect_complex_pairs(fit, H, ev[1:5], completeness = c("probed", "repaired"))
})

test_that("complex_operator() brings sparse complex Hermitian matrices to native sparse routes", {
  g <- complex_hermitian_grid(20)
  op <- complex_operator(g$re, g$im)
  expect_identical(op$dtype, "complex")
  expect_identical(op$structure$kind, "hermitian")
  X <- matrix(complex(real = rnorm(800), imaginary = rnorm(800)), 400)
  expect_equal(asNamespace("eigencore")$apply_operator(op, X), g$H %*% X, tolerance = 1e-12)
  expect_equal(asNamespace("eigencore")$apply_adjoint_operator(op, X),
               Conj(t(g$H)) %*% X, tolerance = 1e-12)
  ev <- eigen(g$H, symmetric = TRUE, only.values = TRUE)$values
  fit <- eig_partial(op, k = 5, target = smallest(), seed = 12)
  expect_complex_pairs(fit, g$H, sort(ev)[1:5])
  inner <- fit$plan$problem
  expect_identical(fit$plan$method,
                   "native realified complex Hermitian route (real 2n embedding)")
  fit <- eig_partial(op, k = 5, target = nearest(2), seed = 12)
  expect_complex_pairs(fit, g$H, ev[order(abs(ev - 2))][1:5])
  expect_match(fit$restart$inner_method, "shift-invert")
  fit <- eig_partial(op, k = 4, method = shift_invert(3), seed = 12)
  expect_complex_pairs(fit, g$H, ev[order(abs(ev - 3))][1:4])
  # interval(): sliced and counted on the sparse embedding.
  fit <- eig_partial(op, target = interval(2, 2.5), seed = 12)
  inside <- sort(ev[ev >= 2 & ev <= 2.5])
  expect_equal(sort(fit$values), inside, tolerance = 1e-9)
  expect_true(certificate(fit)$target_completeness %in% c("exact", "inertia_verified"))
  expect_true(certificate(fit)$passed)
  # Above the dense interval limit the embedding is sliced with sparse LDL'.
  old <- options(eigencore.interval_dense_limit = 100L)
  on.exit(options(old), add = TRUE)
  fit <- eig_partial(op, target = interval(2, 2.5), seed = 12)
  expect_equal(sort(fit$values), inside, tolerance = 1e-9)
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
  expect_match(fit$restart$inner_method, "slicing")
  # A non-Hermitian split operator is general.
  expect_identical(complex_operator(g$re, g$re)$structure$kind, "general")
  expect_error(complex_operator(g$re, matrix(0, 3, 3)), "same dimensions")
})

test_that("complex Hermitian-definite pencils solve through the embedding", {
  n <- 120
  H <- complex_hermitian_fixture(n, seed = 13L)
  set.seed(14)
  W <- matrix(complex(real = rnorm(n * n), imaginary = rnorm(n * n)), n)
  B <- Conj(t(W)) %*% W / n + diag(n)
  B <- (B + Conj(t(B))) / 2
  eb <- eigen(B, symmetric = TRUE)
  Bmh <- eb$vectors %*% (Conj(t(eb$vectors)) / sqrt(eb$values))
  C <- Bmh %*% H %*% Bmh
  ev <- eigen((C + Conj(t(C))) / 2, symmetric = TRUE, only.values = TRUE)$values
  for (method in list(auto(), lanczos(), lobpcg())) {
    fit <- eig_partial(H, k = 4, B = B, method = method, seed = 15)
    expect_complex_pairs(fit, H, ev[1:4], B = B)
    fit <- eig_partial(H, k = 3, B = B, target = smallest(), method = method, seed = 15)
    expect_complex_pairs(fit, H, sort(ev)[1:3], B = B)
  }
  # A real SPD metric with a complex A.
  d <- seq(1, 2, length.out = n)
  fit <- eig_partial(H, k = 3, B = diag(d), seed = 15)
  Dmh <- diag(1 / sqrt(d))
  ed <- eigen(Dmh %*% H %*% Dmh, symmetric = TRUE, only.values = TRUE)$values
  expect_complex_pairs(fit, H, ed[1:3], B = diag(d))
})

test_that("eigs_sym() accepts complex Hermitian input", {
  H <- complex_hermitian_fixture(200, seed = 16L)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  res <- eigs_sym(H, k = 4, which = "LA")
  expect_equal(res$values, ev[1:4], tolerance = 1e-9)
  expect_true(isTRUE(res$certificate$passed))
  res <- eigs_sym(H, k = 3, which = "SA")
  expect_equal(sort(res$values), sort(ev)[1:3], tolerance = 1e-9)
})

test_that("realified results report the inner route and transfer work", {
  H <- complex_hermitian_fixture(200, seed = 17L)
  fit <- eig_partial(H, k = 3, seed = 18)
  expect_identical(fit$restart$kind, "realified_complex_hermitian")
  expect_identical(fit$restart$inner_returned, 6L)
  expect_identical(certificate(fit)$completeness$realified, TRUE)
  w <- work(fit)
  expect_gt(w$operator_columns, 0)
  # values-only and uncertified requests
  fit <- eig_partial(H, k = 3, vectors = FALSE, seed = 18)
  expect_null(fit$vectors)
  expect_true(isTRUE(certificate(fit)$passed))
})

test_that("mapped-back pairs that miss the tolerance are polished and re-certified", {
  # Integration regression: a repaired realified set (interior target,
  # squared-shift repair) mapped back with residuals ~1e-6 against tol 1e-8.
  ns <- asNamespace("eigencore")
  H <- complex_hermitian_fixture(300, seed = 2L)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  fit <- eig_partial(H, k = 6, target = smallest_magnitude(), seed = 3)
  expect_complex_pairs(fit, H, ev[order(abs(ev))][1:6])

  g <- complex_hermitian_grid(12)
  cases <- list(
    dense = list(op = H, H = H),
    split = list(op = complex_operator(Re(H), Im(H)), H = H),
    sparse = list(op = complex_operator(g$re, g$im), H = g$H),
    matrix_free = list(op = linear_operator(
      dim = c(300, 300),
      apply = function(X, alpha = 1, beta = 0, Y = NULL) alpha * (H %*% X),
      dtype = "complex", structure = hermitian()), H = H)
  )
  for (nm in names(cases)) {
    Hm <- cases[[nm]]$H
    P <- eigen_problem(cases[[nm]]$op, target = smallest_magnitude())
    e <- eigen(Hm, symmetric = TRUE)
    id <- order(abs(e$values))[1:4]
    set.seed(21)
    V <- e$vectors[, id] + 1e-6 * matrix(complex(real = rnorm(nrow(Hm) * 4),
                                                 imaginary = rnorm(nrow(Hm) * 4)), nrow(Hm))
    V <- qr.Q(qr(V))
    vals <- Re(colSums(Conj(V) * (Hm %*% V)))
    cert <- ns$complex_hermitian_certify(P, vals, V, 1e-10)
    expect_false(isTRUE(cert$passed), info = nm)
    out <- ns$complex_hermitian_polish(P, vals, V, cert, 1e-10)
    expect_false(is.null(out), info = nm)
    expect_true(isTRUE(out$certificate$passed), info = nm)
    expect_equal(sort(out$values), sort(e$values[id]), tolerance = 1e-10, info = nm)
    expect_identical(out$method, "shifted_inverse_iteration", info = nm)
  }
  # Above the materialisation limit a matrix-free operator has no solve:
  # the residual-started block Krylov expansion still improves the pairs
  # (and never returns a worse set).
  old <- options(eigencore.completeness_materialize_limit = 10L)
  on.exit(options(old), add = TRUE)
  P <- eigen_problem(cases$matrix_free$op, target = largest())
  e <- eigen(H, symmetric = TRUE)
  set.seed(22)
  V <- qr.Q(qr(e$vectors[, 1:3] + 1e-6 * matrix(complex(real = rnorm(900),
                                                         imaginary = rnorm(900)), 300)))
  vals <- Re(colSums(Conj(V) * (H %*% V)))
  cert <- ns$complex_hermitian_certify(P, vals, V, 1e-10)
  out <- ns$complex_hermitian_polish(P, vals, V, cert, 1e-10)
  expect_identical(out$method, "block_krylov_residual_expansion")
  expect_true(isTRUE(out$certificate$passed))
  expect_equal(out$values, e$values[1:3], tolerance = 1e-10)
})
