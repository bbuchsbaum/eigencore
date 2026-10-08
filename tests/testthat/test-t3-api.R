# Tranche 3 API workstream: maxit is an iteration limit (C15), left
# eigenvectors on request in the eigs() shim (C46), default subspace sizes
# (C43), target-ordered thick-restart output (C44), nonsymmetric
# smallest-magnitude and shift-invert Arnoldi (C41), and truthful
# factorization/route labels (C33/C35).

t3_sym <- function(n, seed) {
  set.seed(seed)
  crossprod(matrix(rnorm(n * n), n)) / n
}

t3_sparse_sym <- function(n, density, seed) {
  set.seed(seed)
  S <- Matrix::rsparsematrix(n, n, density = density)
  S <- S + Matrix::t(S) + Matrix::Diagonal(n, x = seq_len(n) / n)
  methods::as(methods::as(S, "generalMatrix"), "CsparseMatrix")
}

# ---------------------------------------------------------------- C15 maxit

test_that("maxit is a restart limit and never the Krylov subspace size", {
  S <- t3_sparse_sym(600, 0.01, 1)
  plain <- plan_solver(eigen_problem(S), k = 4L)
  capped <- plan_solver(eigen_problem(S), k = 4L, maxit = 7L)
  expect_identical(capped$method, "native scalar thick-restart Hermitian Lanczos")
  expect_identical(capped$controls$max_subspace, plain$controls$max_subspace)
  expect_identical(capped$controls$max_restarts, 7L)
  expect_identical(capped$controls$iteration_limit, 7L)
  expect_identical(capped$controls$iteration_limit_kind, "thick_restart_cycles")
  expect_identical(plain$controls$max_restarts, 100L)

  fit <- eig_partial(S, k = 4L, maxit = 7L, seed = 3)
  expect_identical(fit$restart$max_restarts, 7L)
  expect_identical(fit$restart$max_subspace, plain$controls$max_subspace)
  expect_lte(fit$restart$restarts_used, 7L)

  # A tiny restart budget on a hard target genuinely limits the work.
  starved <- eig_partial(S, k = 4L, target = smallest(), maxit = 1L, seed = 3,
                         tol = 1e-14)
  expect_lte(starved$restart$restarts_used, 1L)
})

test_that("max_subspace on the method descriptor sets the subspace", {
  S <- t3_sparse_sym(500, 0.01, 2)
  fit <- eig_partial(S, k = 3L, method = lanczos(max_subspace = 15L), seed = 1)
  expect_identical(fit$plan$controls$max_subspace, 15L)
  expect_identical(fit$restart$max_subspace, 15L)
  expect_true(fit$certificate$passed)

  auto_fit <- eig_partial(S, k = 3L, method = auto(max_subspace = 15L), seed = 1)
  expect_identical(auto_fit$plan$method, fit$plan$method)
  expect_identical(auto_fit$plan$controls$max_subspace, 15L)
  expect_true(auto_fit$certificate$passed)

  expect_error(auto(max_subspace = 1.5), "max_subspace")
  expect_error(lanczos(max_subspace = 0), "max_subspace")
  expect_error(shift_invert(1, max_subspace = -3), "max_subspace")
})

test_that("conflicting iteration limits are an error; agreeing ones are not", {
  A <- t3_sym(60, 4)
  expect_error(
    eig_partial(A, k = 2L, method = lanczos(max_restarts = 3L), maxit = 4L),
    "conflicts with lanczos\\(max_restarts = 3\\)"
  )
  expect_error(
    eig_partial(A, k = 2L, method = lobpcg(maxit = 30L), maxit = 40L),
    "conflicts with lobpcg\\(maxit = 30\\)"
  )
  ok <- plan_solver(eigen_problem(A), k = 2L,
                    method = lanczos(max_restarts = 9L), maxit = 9L)
  expect_identical(ok$controls$max_restarts, 9L)
})

test_that("maxit maps to each route's own iteration limit", {
  A <- t3_sym(80, 5)
  lob <- plan_solver(eigen_problem(A), k = 2L, method = lobpcg(), maxit = 33L)
  expect_identical(lob$controls$maxit, 33L)
  expect_identical(lob$controls$iteration_limit_kind, "lobpcg_iterations")
  expect_identical(plan_solver(eigen_problem(A), k = 2L, method = lobpcg())$controls$maxit,
                   200L)

  set.seed(6)
  N <- matrix(rnorm(80 * 80), 80)
  arn <- plan_solver(eigen_problem(N), k = 3L, maxit = 25L)
  expect_identical(arn$method, eigencore:::native_refined_arnoldi_label())
  expect_identical(arn$controls$krylov_schur_max_iterations, 25L)
  expect_identical(arn$controls$iteration_limit_kind, "krylov_schur_restarts")
  expect_identical(arn$controls$max_subspace,
                   eigencore:::native_krylov_schur_default_ncv(80L, 3L))
  fit <- eig_partial(N, k = 3L, target = largest_magnitude(), maxit = 25L)
  expect_identical(fit$restart$krylov_schur_maxit, 25L)

  # Unrestarted reference Lanczos: one iteration is one Lanczos step.
  op <- linear_operator(dim = dim(A), apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (A %*% X)
    if (!is.null(Y) && beta != 0) Z + beta * Y else Z
  }, structure = hermitian())
  ref <- plan_solver(eigen_problem(op), k = 2L, maxit = 12L)
  expect_identical(ref$method, "reference Hermitian Lanczos (prototype/oracle fallback)")
  expect_identical(ref$controls$iteration_limit_kind, "lanczos_steps")
  expect_identical(ref$controls$max_subspace, 12L)

  dense <- plan_solver(eigen_problem(diag(5)), k = 2L, maxit = 3L)
  expect_identical(dense$controls$iteration_limit_kind, "not_applicable")
})

test_that("shift-invert Lanczos takes maxit as a restart limit", {
  S <- t3_sparse_sym(400, 0.01, 7)
  p <- plan_solver(eigen_problem(S, target = nearest(0.4)), k = 3L, maxit = 50L)
  expect_identical(p$method,
                   "native thick-restart Hermitian Lanczos shift-invert (sparse LU solve callback)")
  expect_identical(p$controls$iteration_limit_kind, "thick_restart_cycles")
  expect_identical(p$controls$max_restarts, 50L)
  expect_identical(p$controls$max_subspace,
                   eigencore:::default_shift_invert_max_subspace(400L, 3L))
  sub <- plan_solver(eigen_problem(S, target = nearest(0.4)), k = 3L,
                     method = auto(max_subspace = 12L))
  expect_identical(sub$controls$max_subspace, 12L)
})

test_that("RSpectra opts map ncv to max_subspace and maxitr to maxit", {
  S <- t3_sparse_sym(500, 0.01, 8)
  base <- eigs_sym(S, 3, which = "LA")
  fit <- eigs_sym(S, 3, which = "LA", opts = list(ncv = 17, maxitr = 40))
  plan <- fit$diagnostics$plan
  expect_identical(plan$method, base$diagnostics$plan$method)
  expect_identical(plan$controls$max_subspace, 17L)
  expect_identical(plan$controls$max_restarts, 40L)
  expect_identical(plan$execution$maxit, 40L)
  expect_equal(fit$values, base$values, tolerance = 1e-8)

  set.seed(9)
  N <- matrix(rnorm(100 * 100), 100)
  g <- eigs(N, 3, opts = list(ncv = 30, maxitr = 77))
  expect_identical(g$diagnostics$plan$controls$max_subspace, 30L)
  expect_identical(g$diagnostics$plan$controls$krylov_schur_max_iterations, 77L)
  expect_error(eigs_sym(S, 3, opts = list(ncv = 10), method = lanczos()),
               "opts\\$ncv cannot be combined")

  # svds: ncv reaches the Golub-Kahan / implicit-Gram subspace.
  set.seed(10)
  X <- Matrix::rsparsematrix(300, 120, density = 0.05)
  sv <- svd_partial(X, 3, method = auto(max_subspace = 25L))
  expect_true(sv$certificate$passed)
  sv2 <- svds(X, 3, opts = list(ncv = 25))
  expect_equal(sv2$d, sv$d, tolerance = 1e-8)
})

# ---------------------------------------------------------- C46 left vectors

test_that("eigs() computes left eigenvectors only on request", {
  set.seed(11)
  N <- matrix(rnorm(120 * 120), 120)
  right <- eigs(N, 4)
  expect_null(right$left_vectors)
  expect_null(right$left_certificate)
  expect_equal(right$diagnostics$work$adjoint_block_calls, 0)
  expect_true(right$certificate$passed)

  both <- eigs(N, 4, left = TRUE)
  expect_equal(ncol(both$left_vectors), 4L)
  expect_true(both$left_certificate$passed)
  expect_gt(both$diagnostics$work$adjoint_block_calls, 0)
  expect_equal(both$values, right$values, tolerance = 1e-8)
  expect_error(eigs(N, 2, left = NA), "left must be TRUE or FALSE")
})

test_that("eig_partial keeps two-sided certificates by default", {
  set.seed(12)
  N <- matrix(rnorm(90 * 90), 90)
  fit <- eig_partial(N, 3, target = largest_magnitude())
  expect_identical(fit$plan$execution$left_vectors, "auto")
  expect_false(is.null(left_vectors(fit)))
  expect_true(fit$left_certificate$passed)

  none <- eig_partial(N, 3, target = largest_magnitude(), left_vectors = "none")
  expect_null(none$left_vectors)
  expect_null(none$left_certificate)
  expect_false(any(grepl("left eigenvectors", none$warnings)))

  # A matrix-free operator without an adjoint cannot return left vectors.
  op <- linear_operator(dim = dim(N), apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (N %*% X)
    if (!is.null(Y) && beta != 0) Z + beta * Y else Z
  }, structure = general())
  expect_error(
    eig_partial(op, 2, target = largest_magnitude(), left_vectors = "compute"),
    "left eigenvectors are unavailable"
  )
  replanned <- solve(none$plan, replan = TRUE)
  expect_identical(replanned$plan$execution$left_vectors, "none")
})

# ------------------------------------------------------------ C44 ordering

test_that("thick-restart results come back in target order", {
  # Near-repeated clusters lock out of order inside the kernel.
  top <- c(9, 9 - 1e-3, 7, 7 - 1e-3, 7 - 2e-3)
  d <- c(top, 5, seq(4, 0.01, length.out = 194))
  set.seed(13)
  Q <- qr.Q(qr(matrix(rnorm(200 * 200), 200)))
  A <- Q %*% diag(d) %*% t(Q)
  A <- (A + t(A)) / 2
  fit <- eig_partial(A, k = 5L, method = lanczos(max_subspace = 12L), seed = 2)
  expect_identical(fit$method, "native scalar thick-restart Hermitian Lanczos")
  expect_equal(fit$values, top, tolerance = 1e-8)
  expect_equal(sort(fit$locked), seq_len(5L))
  expect_false(is.unsorted(rev(fit$values)))
  r <- colSums((A %*% fit$vectors - sweep(fit$vectors, 2L, fit$values, `*`))^2)
  expect_lt(max(sqrt(r)), 1e-6)

  small <- eig_partial(A, k = 4L, target = smallest(),
                       method = lanczos(max_subspace = 14L), seed = 2)
  expect_false(is.unsorted(small$values))
})

# ----------------------------------------------- C41 nonsymmetric targets

test_that("nonsymmetric smallest_magnitude matches eigen() (dense and sparse)", {
  for (seed in 21:23) {
    set.seed(seed)
    A <- matrix(rnorm(120 * 120), 120) / sqrt(120)
    ev <- eigen(A, only.values = TRUE)$values
    fit <- eig_partial(A, 4, target = smallest_magnitude())
    expect_identical(fit$plan$method,
                     eigencore:::shift_invert_arnoldi_label("dense_qr"))
    expect_true(fit$certificate$passed)
    expect_equal(sort(Mod(fit$values)), sort(Mod(ev[order(Mod(ev))][1:4])),
                 tolerance = 1e-8)
    expect_false(is.unsorted(Mod(fit$values)))
  }
  set.seed(24)
  S <- Matrix::rsparsematrix(800, 800, density = 0.005) +
    Matrix::Diagonal(800, x = stats::runif(800, 0.5, 2))
  ev <- eigen(as.matrix(S), only.values = TRUE)$values
  fit <- eig_partial(S, 3, target = smallest_magnitude())
  expect_identical(fit$plan$method,
                   eigencore:::shift_invert_arnoldi_label("sparse_lu"))
  expect_true(fit$certificate$passed)
  expect_equal(sort(Mod(fit$values)), sort(Mod(ev[order(Mod(ev))][1:3])),
               tolerance = 1e-8)
})

test_that("singular A perturbs the implicit sigma = 0 shift", {
  set.seed(25)
  B <- matrix(rnorm(80 * 80), 80)
  B[, 1] <- B[, 2]
  ev <- eigen(B, only.values = TRUE)$values
  fit <- eig_partial(B, 3, target = smallest_magnitude())
  expect_true(fit$transform$sigma_perturbed)
  expect_true(fit$certificate$passed)
  expect_equal(sort(Mod(fit$values)), sort(Mod(ev[order(Mod(ev))][1:3])),
               tolerance = 1e-7)
  # An explicit sigma is never perturbed silently.
  expect_error(eig_partial(B, 2, method = shift_invert(0)), "singular")
})

test_that("matrix-free smallest_magnitude runs the native Krylov-Schur ranking", {
  set.seed(26)
  A <- matrix(rnorm(40 * 40), 40) + diag(seq(1, 8, length.out = 40))
  ev <- eigen(A, only.values = TRUE)$values
  op <- linear_operator(dim = dim(A), apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (A %*% X)
    if (!is.null(Y) && beta != 0) Z + beta * Y else Z
  }, structure = general(), metadata = list(frobenius_norm = sqrt(sum(A^2))))
  fit <- eig_partial(op, 2, target = smallest_magnitude(), left_vectors = "none")
  expect_identical(fit$plan$method, eigencore:::native_matrix_free_arnoldi_label())
  expect_true(all(fit$certificate$converged))
  expect_equal(sort(Mod(fit$values)), sort(Mod(ev[order(Mod(ev))][1:2])),
               tolerance = 1e-8)
})

test_that("nonsymmetric nearest(sigma) uses shift-invert Arnoldi", {
  set.seed(27)
  A <- matrix(rnorm(150 * 150), 150)
  ev <- eigen(A, only.values = TRUE)$values
  sigma <- 0.7
  fit <- eig_partial(A, 4, target = nearest(sigma))
  expect_identical(fit$plan$method,
                   eigencore:::shift_invert_arnoldi_label("dense_qr"))
  expect_identical(fit$plan$controls$iteration_limit_kind, "krylov_schur_restarts")
  expect_true(fit$certificate$passed)
  expect_true(fit$left_certificate$passed)
  expect_lt(max(Mod(fit$biorthogonality - diag(4))), 1e-6)
  ref <- ev[order(Mod(ev - sigma))][1:4]
  expect_equal(sort(Mod(fit$values - sigma)), sort(Mod(ref - sigma)),
               tolerance = 1e-8)
  expect_equal(fit$transform$certification$problem, "original")

  shim <- eigs(A, 4, sigma = sigma)
  expect_equal(shim$values, fit$values, tolerance = 1e-8)
  expect_null(shim$left_vectors)
  expect_error(eigs(A, 2, sigma = 1 + 2i), "real sigma")

  set.seed(28)
  S <- Matrix::rsparsematrix(1500, 1500, density = 0.004) +
    Matrix::Diagonal(1500, x = stats::rnorm(1500))
  evs <- eigen(as.matrix(S), only.values = TRUE)$values
  sfit <- eig_partial(S, 3, target = nearest(0.25), left_vectors = "compute")
  expect_identical(sfit$plan$method,
                   eigencore:::shift_invert_arnoldi_label("sparse_lu"))
  expect_true(sfit$certificate$passed)
  expect_true(sfit$left_certificate$passed)
  sref <- evs[order(Mod(evs - 0.25))][1:3]
  expect_equal(sort(Mod(sfit$values - 0.25)), sort(Mod(sref - 0.25)),
               tolerance = 1e-8)

  sub <- eig_partial(A, 3, target = nearest(sigma), method = auto(max_subspace = 12L))
  expect_identical(sub$plan$controls$max_subspace, 12L)
  expect_true(sub$certificate$passed)
})

test_that("nonsymmetric shift_invert(solve =) runs without left vectors", {
  set.seed(29)
  A <- matrix(rnorm(60 * 60), 60)
  sigma <- 0.2
  M <- A - sigma * diag(60)
  fit <- eig_partial(A, 2, target = nearest(sigma),
                     method = shift_invert(sigma, solve = function(X) solve(M, X)))
  expect_identical(fit$plan$method,
                   eigencore:::shift_invert_arnoldi_label("user_solve"))
  expect_true(fit$certificate$passed)
  expect_null(fit$left_vectors)
  expect_true(any(grepl("left eigenvectors unavailable", fit$warnings)))
})

# ------------------------------------------------------ C33/C35 labels

test_that("shift-invert and metric-solve labels name the real factorization", {
  n <- 60
  T <- Matrix::bandSparse(n, k = c(-1, 0, 1),
                          diagonals = list(rep(-1, n - 1), rep(2, n), rep(-1, n - 1)))
  T <- methods::as(methods::as(T, "generalMatrix"), "CsparseMatrix")
  fit <- eig_partial(T, 2, target = nearest(1), method = shift_invert(1))
  expect_identical(fit$transform$label_kind, "tridiagonal_lu_native")
  expect_identical(fit$transform$factorization_cache$condition_estimate_type,
                   "tridiagonal_lu_pivot_ratio")
  expect_false(grepl("thomas", fit$transform$label_kind, ignore.case = TRUE))

  S <- t3_sparse_sym(300, 0.02, 30)
  sfit <- eig_partial(S, 2, target = nearest(0.5))
  expect_identical(sfit$method,
                   "native thick-restart Hermitian Lanczos shift-invert (sparse LU solve callback)")
  expect_identical(sfit$restart$kind, "native_thick_restart_shift_invert_callback")

  pre <- shifted_tridiagonal_preconditioner(T, shift = 0.1)
  expect_identical(attr(pre, "eigencore_preconditioner")$factorization,
                   "tridiagonal_lu")
})
