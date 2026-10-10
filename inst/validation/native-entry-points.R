#!/usr/bin/env Rscript
# Exercise every native (.Call) entry point of eigencore at small sizes,
# mostly through the public API, for memory checkers:
#
#   R -d "valgrind --error-exitcode=1 --track-origins=yes" --vanilla \
#     -f inst/validation/native-entry-points.R
#
#   LD_PRELOAD=$(gcc -print-file-name=libasan.so) \
#     ASAN_OPTIONS=detect_leaks=0:abort_on_error=1 \
#     Rscript inst/validation/native-entry-points.R      # ASan/UBSan build
#
# Each section runs inside tryCatch; an R-level error in a section is reported
# and counted (the script stops with an error at the end if any section
# failed), so a memory checker sees every section even when one breaks.
# Expected-error sections (malformed input) assert that the error is raised.
# Sizes are kept tiny so the whole script finishes in a few minutes under
# valgrind.

suppressPackageStartupMessages({
  library(eigencore)
  library(Matrix)
})

options(warn = 1)
threads <- as.integer(Sys.getenv("EIGENCORE_ENTRY_THREADS", "2"))
options(eigencore.threads = threads)

failures <- character()
section <- function(name, expr) {
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({
    force(expr)
    TRUE
  }, error = function(e) {
    message(sprintf("SECTION FAILED [%s]: %s", name, conditionMessage(e)))
    FALSE
  })
  if (!ok) failures <<- c(failures, name)
  cat(sprintf("%-48s %s (%.1fs)\n", name, if (ok) "ok" else "FAILED",
              proc.time()[["elapsed"]] - t0))
  invisible(ok)
}
expect_err <- function(expr, pattern = NULL) {
  msg <- tryCatch({
    force(expr)
    NULL
  }, error = function(e) conditionMessage(e))
  if (is.null(msg)) stop("expected an error, none raised")
  if (!is.null(pattern) && !grepl(pattern, msg)) {
    stop("unexpected error message: ", msg)
  }
  invisible(msg)
}
cert_ok <- function(fit) {
  cert <- certificate(fit)
  if (!isTRUE(cert$passed)) {
    stop("certificate failed (max backward error ",
         format(cert$max_backward_error), ")")
  }
  invisible(cert)
}
ns <- asNamespace("eigencore")
call_native <- function(name, ...) .Call(name, ..., PACKAGE = "eigencore")

set.seed(20261009)
path_laplacian <- function(n) {
  Matrix::bandSparse(n, k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L)))
}
n <- 60L
S <- path_laplacian(n) + Matrix::Diagonal(n, x = seq(0, 1, length.out = n))
S <- as(as(S, "generalMatrix"), "CsparseMatrix")
Sd <- as.matrix(S)
Q <- qr.Q(qr(matrix(rnorm(n * n), n)))
D <- Q %*% diag(seq_len(n) / 4) %*% t(Q)
D <- (D + t(D)) / 2
R <- Matrix::rsparsematrix(80L, 30L, density = 0.15)
Rd <- as.matrix(R)
G <- matrix(rnorm(70L * 25L), 70L) %*% diag(exp(-seq(0, 5, length.out = 25L)))
NS <- Matrix::rsparsematrix(n, n, density = 0.08) + Matrix::Diagonal(n, x = 1:n)
NS <- as(NS, "CsparseMatrix")
NSd <- as.matrix(NS)
Bspd <- diag(seq(1, 2, length.out = n))
Bsp <- as(Matrix::Diagonal(n, x = seq(1, 2, length.out = n)) +
            0.1 * path_laplacian(n), "CsparseMatrix")
Bsp <- as(as(Bsp, "generalMatrix"), "CsparseMatrix")

cat("eigencore threads:", getOption("eigencore.threads"), "\n")

# --- runtime / threads / unwind ---------------------------------------------
section("thread info + set default", {
  call_native("eigencore_thread_info")
  call_native("eigencore_thread_governor", NA_real_)
  call_native("eigencore_set_default_threads", NA_real_)
  ns$eigencore_threads()
})
section("unwind selftest (ok/error/bad_alloc/r_alloc)", {
  st <- function(mode) call_native("eigencore_unwind_selftest", mode)
  stopifnot(identical(st("live"), 0L))
  st("ok")
  expect_err(st("error"))
  expect_err(st("bad_alloc"))
  expect_err(st("length_error"))
  expect_err(st("huge_vector"))
  expect_err(st("r_alloc_failure"))
  for (mode in c("interrupt", "nested_unwind", "r_error")) {
    try(tryCatch(st(mode), interrupt = function(i) NULL), silent = TRUE)
  }
  stopifnot(identical(st("live"), 0L))
})
section("int guards / noalloc check", {
  call_native("eigencore_dense_apply_int_guard_check")
  X <- matrix(rnorm(n * 2), n)
  ns$native_apply_noalloc_check("dense", D, X)
  ns$native_apply_noalloc_check("csc", S, X)
})

# --- identity hash ----------------------------------------------------------
section("identity hash (dense, sparse, operator)", {
  operator_identity(as_operator(D))
  operator_identity(as_operator(S))
  operator_identity(as_operator(R))
  call_native("eigencore_identity_hash", list(1:3, c(1.5, NA), "a", NULL, TRUE,
                                               list(x = 2)))
})

# --- dense / CSC / diagonal / composite applies ------------------------------
section("operator applies (dense, complex, CSC, diagonal)", {
  for (A in list(D, S, R, Diagonal(n, x = 1:n))) {
    op <- as_operator(A)
    X <- matrix(rnorm(op$dim[2L] * 3L), ncol = 3L)
    Y <- ns$apply_operator(op, X)
    stopifnot(max(abs(as.matrix(Y) - as.matrix(A %*% X))) < 1e-10)
    Z <- ns$apply_adjoint_operator(op, matrix(rnorm(op$dim[1L] * 2L), ncol = 2L))
  }
  Cz <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), 3)
  opz <- as_operator(Cz)
  ns$apply_operator(opz, diag(3))
  ns$csc_column_moments(R)
  ns$csc_apply_repeat(R, matrix(rnorm(30 * 12), 30), reps = 3L)
  ns$csc_apply_repeat(R, matrix(rnorm(80 * 12), 80), transpose = TRUE, reps = 3L)
})
section("centered / centered-scaled / composite operators", {
  ap <- ns$apply_operator
  apt <- ns$apply_adjoint_operator
  for (cols in c(1L, 3L, 12L)) {
    X <- matrix(rnorm(ncol(R) * cols), ncol = cols)
    W <- matrix(rnorm(nrow(R) * cols), ncol = cols)
    cen <- center(R)
    cs <- scale_cols(center(R), runif(ncol(R), 0.5, 2))
    stopifnot(max(abs(ap(cen, X) - scale(Rd, scale = FALSE) %*% X)) < 1e-10)
    ap(cs, X); apt(cs, W); apt(cen, W)
    both <- center(R, rows = TRUE, columns = TRUE)
    ap(both, X)
    big <- Matrix::rsparsematrix(400L, 300L, density = 0.05)
    cp <- compose(as_operator(big), as_operator(Matrix::rsparsematrix(300L, 30L, 0.1)))
    ap(cp, X); apt(cp, matrix(rnorm(400 * cols), ncol = cols))
    cpo <- crossprod_operator(big)
    ap(cpo, matrix(rnorm(300 * cols), ncol = cols))
    sr <- scale_rows(as_operator(big), runif(400))
    ap(sr, matrix(rnorm(300 * cols), ncol = cols))
    adj <- adjoint(as_operator(R))
    ap(adj, W)
  }
})

# --- orthogonalisation utilities --------------------------------------------
section("orthogonalisation kernels", {
  X <- matrix(rnorm(n * 4), n)
  ns$native_mgs2(X)
  ns$native_cholqr2(X)
  ns$native_b_cholqr2(X, Bspd)
  ns$native_diagonal_b_cholqr2(X, diag(Bspd))
  ns$native_diagonal_b_cholqr2(X, rep(1, n), unit = TRUE)
  Qb <- qr.Q(qr(matrix(rnorm(n * 3), n)))
  ns$reorthogonalize_against(X, Qb, 2L)
  ws <- ns$basis_workspace(n, 8L, 4L)
  ns$basis_workspace_info(ws)
  ns$reorthogonalize_against(X, Qb, 2L, workspace = ws)
  ns$orthogonality_loss(Qb)
  ns$orthogonality_loss(Qb, Bspd)
  ns$native_rayleigh_ritz_symmetric(D, Qb)
  ns$rayleigh_ritz(D, Qb)
  ns$cholqr2(X); ns$mgs2(X); ns$b_orthogonalize(X, Bspd)
})

# --- Hermitian eigen: dense / CSC / matrix-free ------------------------------
section("Lanczos dense + CSC (largest/smallest)", {
  cert_ok(eig_partial(D, 4L, target = largest(), method = lanczos(), seed = 1))
  eig_partial(S, 4L, target = smallest(), method = lanczos(max_subspace = 40L), seed = 2)
  cert_ok(eigs_sym(S, 3L, which = "LA"))
  eigs_sym(S, 3L, which = "BE")
  call_native("eigencore_lanczos_dense", D, 20L, rnorm(n), 3L, 1L, 1e-8)
})
section("block Lanczos dense + CSC + matrix-free", {
  cert_ok(eig_partial(D, 4L, method = lanczos(block = 3L), seed = 3))
  cert_ok(eig_partial(S, 4L, target = smallest(),
                      method = lanczos(block = 2L, max_subspace = 16L), seed = 4))
  mf <- linear_operator(c(n, n), apply = function(X, ...) Sd %*% X,
                        apply_adjoint = function(X, ...) Sd %*% X,
                        structure = hermitian())
  eig_partial(mf, 3L, method = lanczos(), seed = 5)
  cert_ok(eig_partial(mf, 3L, seed = 5))
  cert_ok(eig_partial(crossprod_operator(R), 3L, seed = 6))
})
section("completeness probe (repeated eigenvalue)", {
  Dr <- diag(c(9, 9, 7, 7, 7, seq(5, 1, length.out = 25)))
  f <- eig_partial(Dr, 5L, method = lanczos(completeness = "probe"), seed = 7)
  f2 <- eig_partial(as(Dr, "CsparseMatrix"), 5L,
                    method = lanczos(completeness = "probe"), seed = 7)
  mf <- linear_operator(dim(Dr), apply = function(X, ...) Dr %*% X,
                        apply_adjoint = function(X, ...) Dr %*% X,
                        structure = hermitian())
  f3 <- eig_partial(mf, 5L, method = lanczos(completeness = "probe"), seed = 7)
})
section("shift-invert (dense, tridiagonal, generalized)", {
  cert_ok(eig_partial(D, 2L, target = nearest(3.1), method = shift_invert(3.1)))
  cert_ok(eig_partial(S, 2L, target = nearest(0.05), method = shift_invert(0.05),
                      seed = 8, allow_dense_fallback = "never"))
  T3 <- as(as(path_laplacian(n), "generalMatrix"), "CsparseMatrix")
  eig_partial(T3, 2L, B = Diagonal(n, x = rep(2, n)), target = nearest(0.05),
              method = shift_invert(0.05), seed = 9)
  eig_partial(D, 2L, B = Bspd, target = nearest(3.1), method = shift_invert(3.1))
  eig_partial(S, 2L, target = nearest(0.7), method = shift_invert(0.7), seed = 10)
})
section("shift-invert singular shift error path", {
  Ds <- diag(c(1, 2, 3, 4))
  # sigma exactly at an eigenvalue: the factorisation is singular and the
  # driver must either perturb sigma or fail cleanly
  try(eig_partial(Ds, 1L, target = nearest(2), method = shift_invert(2)),
      silent = TRUE)
  try(eig_partial(as(Ds, "CsparseMatrix"), 1L, target = nearest(2),
                  method = shift_invert(2)), silent = TRUE)
  call_native("eigencore_tridiagonal_solve", c(-1, -1, -1), c(2, 2, 2, 2),
              c(-1, -1, -1), matrix(1, 4, 2))
  expect_err(call_native("eigencore_tridiagonal_solve", c(0, 0, 0),
                         c(0, 0, 0, 0), c(0, 0, 0), matrix(1, 4, 2)))
  call_native("eigencore_shift_invert_lanczos_dense", Ds + 0, 2.5, 4L,
              rep(1, 4), 1L, 0L, 1e-10)
  expect_err(call_native("eigencore_shift_invert_lanczos_dense", Ds + 0, 2, 4L,
                         rep(1, 4), 1L, 0L, 1e-10), "perturb sigma")
})
section("LOBPCG standard + generalized (all B storages)", {
  cert_ok(eig_partial(D, 3L, target = smallest(), method = lobpcg(maxit = 200L), seed = 11))
  cert_ok(eig_partial(S, 3L, target = smallest(), method = lobpcg(maxit = 300L), seed = 12))
  pre <- shifted_tridiagonal_preconditioner(path_laplacian(n), 0.1)
  eig_partial(S, 3L, target = smallest(), method = lobpcg(maxit = 300L, preconditioner = pre), seed = 13)
  eig_partial(D, 3L, B = Bspd, target = smallest(), method = lobpcg(maxit = 200L), seed = 14)
  eig_partial(D, 3L, B = diag(Bspd) * diag(n), target = smallest(),
              method = lobpcg(maxit = 200L), seed = 15)
  eig_partial(D, 3L, B = Diagonal(n, x = diag(Bspd)), target = smallest(),
              method = lobpcg(maxit = 200L), seed = 15)
  eig_partial(D, 3L, B = Bsp, target = smallest(), method = lobpcg(maxit = 200L), seed = 16)
  eig_partial(S, 3L, B = Diagonal(n, x = diag(Bspd)), target = smallest(),
              method = lobpcg(maxit = 300L), seed = 17)
  eig_partial(S, 3L, B = Bsp, target = smallest(), method = lobpcg(maxit = 300L), seed = 18)
  eig_partial(Diagonal(n, x = 1:n), 3L, B = Diagonal(n, x = diag(Bspd)),
              target = smallest(), method = lobpcg(maxit = 200L), seed = 19)
  Bop <- linear_operator(c(n, n), apply = function(X, ...) Bspd %*% X,
                         apply_adjoint = function(X, ...) Bspd %*% X,
                         structure = hermitian())
  eig_partial(D, 3L, B = Bop, target = smallest(), method = lobpcg(maxit = 200L), seed = 20)
  eig_partial(S, 3L, B = Bop, target = smallest(), method = lobpcg(maxit = 200L), seed = 21)
  eig_partial(Diagonal(n, x = 1:n), 3L, B = Bop, target = smallest(),
              method = lobpcg(maxit = 200L), seed = 22)
  Cn <- qr.Q(qr(matrix(rnorm(n), n)))
  # constrained Ritz pairs are not eigenpairs of D: no certificate check
  eig_partial(D, 2L, target = smallest(),
              method = lobpcg(maxit = 200L, constraints = Cn), seed = 23)
})

# --- nonsymmetric: Arnoldi / Krylov-Schur ------------------------------------
section("Arnoldi Krylov-Schur dense/CSC/matrix-free", {
  eigs(NSd, 4L, which = "LM")
  eigs(NS, 4L, which = "LR", left = TRUE)
  eigs(NS, 3L, sigma = 0.5)
  eig_partial(NS, 3L, target = smallest_magnitude(), seed = 24)
  mf <- linear_operator(c(n, n), apply = function(X, ...) NSd %*% X,
                        apply_adjoint = function(X, ...) crossprod(NSd, X))
  eig_partial(mf, 3L, target = largest_magnitude(), seed = 25)
  eig_partial(mf, 2L, target = largest_real(), seed = 25)
})
section("Arnoldi cycle/ritz helpers (reference route)", {
  mfa <- linear_operator(c(n, n), apply = function(X, ...) NSd %*% X,
                         apply_adjoint = function(X, ...) crossprod(NSd, X))
  for (A in list(NSd, NS, mfa)) {
    cyc <- ns$native_arnoldi_cycle(A, rnorm(n), 12L)
    ns$native_arnoldi_projected_ritz(cyc)
    co <- ns$native_arnoldi_ritz_coefficients(cyc)
    ns$native_arnoldi_ritz_vectors(cyc, co$coefficients[, 1:3, drop = FALSE])
    ns$native_arnoldi_refined_ritz_vectors(cyc, co$values[1:2])
  }
  ns$native_arnoldi_general(NSd, 2L, target = largest_magnitude(),
                            extraction = "refined_ritz")
})

# --- SVD --------------------------------------------------------------------
section("SVD: Gram (left/right), GK, retained, block GK", {
  cert_ok(svd_partial(R, 3L, seed = 26))
  cert_ok(svd_partial(t(R), 3L, seed = 27))
  cert_ok(svds(R, 3L))
  cert_ok(svd_partial(R, 3L, method = golub_kahan(), seed = 28))
  cert_ok(svd_partial(Rd, 3L, method = golub_kahan(), seed = 29))
  cert_ok(svd_partial(G, 3L, method = golub_kahan(), seed = 30))
  svd_partial(R, 2L, target = smallest(), seed = 31)
  svd_partial(G, 2L, target = smallest(), method = golub_kahan(), seed = 32)
  svd_partial(R, 3L, vectors = "none", seed = 33)
  svd_partial(R, 3L, vectors = "left", seed = 33)
  mf <- linear_operator(dim(Rd), apply = function(X, ...) Rd %*% X,
                        apply_adjoint = function(X, ...) crossprod(Rd, X))
  svd_partial(mf, 3L, seed = 34)
  svd_partial(mf, 2L, target = smallest(), seed = 34)
  svd_partial(center(R), 3L, seed = 35)
  svd_partial(scale_cols(center(R), runif(30, 0.5, 2)), 3L, seed = 35)
  svd_partial(compose(as_operator(R), as_operator(diag(30))), 3L, seed = 36)
})
section("SVD: randomized dense/CSC", {
  svd_partial(G, 3L, method = randomized(), seed = 37)
  svd_partial(R, 3L, method = randomized(normalizer = "lu"), seed = 38)
  svd_partial(R, 3L, method = randomized(normalizer = "none", refine = FALSE),
              seed = 39)
})
section("dense drivers (eig_full, svd, schur, gsvd)", {
  eig_full(D)
  eig_full(D, vectors = FALSE)
  eig_full(NSd[1:10, 1:10])
  eig_full(D[1:10, 1:10], B = Bspd[1:10, 1:10])
  eig_full(NSd[1:8, 1:8], B = NSd[11:18, 11:18])
  Cz <- matrix(complex(real = rnorm(16), imaginary = rnorm(16)), 4)
  Hz <- Cz + Conj(t(Cz))
  eig_full(Hz)
  eig_full(Cz)
  eig_full(Hz, B = crossprod(Cz) + diag(4))
  eig_full(Cz, B = Cz + diag(4))
  generalized_schur(NSd[1:8, 1:8], NSd[11:18, 11:18])
  generalized_schur(Cz, Cz + diag(4))
  generalized_svd(G[1:12, 1:5], G[20:29, 1:5])
  svd_partial(Cz, 2L)
  eig_partial(D, 3L, method = auto())
  ns$native_dense_symmetric_eigen_selected(D, 3L, largest())
  ns$native_dense_symmetric_eigen_selected(D, 3L, smallest(), vectors = FALSE)
})
section("restart state + warm start", {
  f <- eig_partial(S, 3L, method = lanczos(), seed = 40)
  st <- restart_state(f)
  eig_partial(S, 3L, method = lanczos(), initial_subspace = vectors(f), seed = 41)
  p <- plan_solver(eigen_problem(S), k = 3L)
})

# --- certificates on every storage ------------------------------------------
section("certificates (dense/CSC/diagonal/tridiagonal, eigen+SVD)", {
  for (A in list(D, S, Diagonal(n, x = 1:n))) {
    f <- eig_partial(A, 2L, seed = 42)
    certificate(f)
  }
  for (A in list(G, R, Diagonal(n, x = 1:n))) {
    certificate(svd_partial(A, 2L, seed = 43))
  }
})

# --- internal routes the public planner does not pick at these sizes --------
# (randomized controllers, prototype Golub-Kahan, retained IRLBA, block
# Golub-Kahan, implicit Gram, tridiagonal LOBPCG, cached-AV certificates and
# the alternative dense drivers). Called through the package's own wrappers.
section("internal SVD routes (randomized/GK/IRLBA/block GK/Gram)", {
  Gop <- as_operator(G)
  Rop <- as_operator(R)
  mfR <- linear_operator(dim(Rd), apply = function(X, ...) Rd %*% X,
                         apply_adjoint = function(X, ...) crossprod(Rd, X))
  ns$native_dense_randomized_svd(Gop, 3L)
  ns$native_dense_randomized_svd(Gop, 2L, refine = FALSE, vectors = "left")
  ns$native_csc_randomized_svd(Rop, 3L)
  ns$native_csc_randomized_svd(Rop, 3L, vectors = "none")
  pair <- ns$randomized_svd_apply_pair(Gop)
  pair2 <- ns$randomized_svd_apply_pair(Rop)
  for (op in list(Gop, Rop, mfR, center(R))) {
    ns$native_golub_kahan_svd(op, 3L)
    ns$native_golub_kahan_svd(op, 2L, target = smallest())
  }
  for (op in list(Gop, Rop)) {
    ns$native_irlba_lbd_restart_abi(op, 3L)
    ns$native_block_golub_kahan_basis(op, 12L, block = 2L)
    ns$native_block_golub_kahan_fit(op, 12L, 3L, block = 2L)
    ns$native_block_golub_kahan_retained_cycle_svd(op, 3L, block = 2L)
    ns$native_block_golub_kahan_retained_cycle_svd(op, 3L, block = 2L,
                                                   retained_av_cache = TRUE)
    ns$native_block_golub_kahan_retained_restart_abi(op, 3L, block = 2L)
    ns$native_implicit_gram_svd(op, 3L)
    ns$native_implicit_gram_svd(op, 2L, vectors = "right")
  }
  sv <- svd(G)
  ns$native_dense_svd_certificate_cached_av(G, sv$d[1:3], sv$u[, 1:3],
                                            sv$v[, 1:3], G %*% sv$v[, 1:3])
  ns$dense_svd_residuals(G, sv$d[1:3], sv$u[, 1:3], sv$v[, 1:3])
  Vb <- qr.Q(qr(matrix(rnorm(25 * 6), 25)))
  ns$native_block_golub_kahan_ritz(Vb, G %*% Vb, 3L)
  ns$native_dense_svd(G)
})
section("internal eigen routes (tridiagonal LOBPCG, drivers, Lanczos)", {
  T3 <- as(as(path_laplacian(n), "generalMatrix"), "CsparseMatrix")
  ns$native_lobpcg_tridiagonal_hermitian(as_operator(T3), 3L, maxit = 200L)
  ns$native_tridiagonal_eigen_selected(rep(2, 20), rep(-1, 19), 3L, smallest())
  ns$native_tridiagonal_eigen_selected(rep(2, 20), rep(-1, 19), 3L, largest())
  ns$native_dense_symmetric_eigen_dsyev(D)
  Cz <- matrix(complex(real = rnorm(16), imaginary = rnorm(16)), 4)
  ns$native_dense_complex_generalized_hpd_eigen(Cz + Conj(t(Cz)),
                                               Conj(t(Cz)) %*% Cz + diag(4))
  ns$native_lanczos_hermitian(as_operator(D), 3L)
  ns$native_lanczos_hermitian(as_operator(S), 3L, target = smallest())
})

# --- malformed input / error paths ------------------------------------------
section("malformed CSC rejected", {
  bad <- S
  bad@i[3] <- 1000L
  expect_err(call_native("eigencore_csc_block_apply", bad@i, bad@p, bad@x,
                         nrow(bad), ncol(bad), diag(n)[, 1:2], FALSE))
  bad2 <- S
  bad2@p[5] <- bad2@p[5] + 100L
  expect_err(call_native("eigencore_csc_block_apply", bad2@i, bad2@p, bad2@x,
                         nrow(bad2), ncol(bad2), diag(n)[, 1:2], FALSE))
})
section("non-finite input rejected", {
  Dn <- D; Dn[2, 2] <- NaN
  expect_err(eig_partial(Dn, 2L))
  expect_err(svd_partial(Dn, 2L))
})
section("R callback error inside native solver unwinds", {
  calls <- 0L
  mf <- linear_operator(c(n, n), apply = function(X, ...) {
    calls <<- calls + 1L
    if (calls > 3L) stop("callback boom")
    Sd %*% X
  }, apply_adjoint = function(X, ...) Sd %*% X, structure = hermitian())
  expect_err(eig_partial(mf, 3L, method = lanczos(), seed = 44), "boom")
  calls <- 0L
  mf2 <- linear_operator(c(n, n), apply = function(X, ...) {
    calls <<- calls + 1L
    if (calls > 3L) stop("callback boom")
    NSd %*% X
  }, apply_adjoint = function(X, ...) crossprod(NSd, X))
  expect_err(eig_partial(mf2, 3L, target = largest_magnitude(), seed = 45), "boom")
  calls <- 0L
  mf3 <- linear_operator(dim(Rd), apply = function(X, ...) {
    calls <<- calls + 1L
    if (calls > 3L) stop("callback boom")
    Rd %*% X
  }, apply_adjoint = function(X, ...) crossprod(Rd, X))
  expect_err(svd_partial(mf3, 3L, method = golub_kahan(), seed = 46), "boom")
  calls <- 0L
  expect_err(svd_partial(mf3, 3L, seed = 46), "boom")
})

if (length(failures)) {
  stop("native entry-point sections failed: ", paste(failures, collapse = ", "))
}
cat("eigencore native entry points: all sections passed\n")
