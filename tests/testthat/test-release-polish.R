# Release polish: C63 (supernodal explicit zeros), P25 (dense completeness
# count cost) and P26 (parallel work threshold).

rp_laplacian_2d <- function(m) {
  T1 <- Matrix::bandSparse(m, k = c(0, 1),
                           diagonals = list(rep(2, m), rep(-1, m - 1L)),
                           symmetric = TRUE)
  I <- Matrix::Diagonal(m)
  methods::as(kronecker(T1, I) + kronecker(I, T1), "generalMatrix")
}

test_that("C63: the supernodal SPD factor drops relaxed-amalgamation zeros", {
  skip_if_not(eigencore:::cholmod_bridge_available())
  set.seed(63)
  n_r <- 900L
  R <- Matrix::rsparsematrix(n_r, n_r, density = 4 / n_r, symmetric = TRUE)
  R <- R + Matrix::Diagonal(n_r, 2 * max(Matrix::rowSums(abs(R))) + 1)
  # Two disconnected blocks give an elimination forest.
  blocks <- Matrix::bdiag(rp_laplacian_2d(20), rp_laplacian_2d(12))
  mats <- list(
    grid_upper = Matrix::forceSymmetric(rp_laplacian_2d(40), uplo = "U"),
    grid_lower = Matrix::forceSymmetric(rp_laplacian_2d(40), uplo = "L"),
    random = Matrix::forceSymmetric(methods::as(R, "generalMatrix"), uplo = "U"),
    forest = Matrix::forceSymmetric(methods::as(blocks, "generalMatrix"), uplo = "U")
  )
  for (name in names(mats)) {
    As <- methods::as(mats[[name]], "CsparseMatrix")
    n <- nrow(As)
    sigma <- -0.25
    F <- eigencore:::ldl_try_spd_factor(As, sigma)
    expect_s4_class(F, "CHMsimpl")
    nz <- methods::slot(F, "nz")
    p <- methods::slot(F, "p")
    x <- methods::slot(F, "x")
    expect_identical(length(x), as.integer(sum(nz)), info = name)
    expect_identical(as.numeric(p[n + 1L]), sum(as.numeric(nz)), info = name)
    # Off-diagonal stored entries are the symbolic pattern: a numeric
    # refactorisation (which rebuilds the pattern from the elimination
    # tree) has exactly as many, and none of them is an explicit zero here.
    diag_pos <- p[seq_len(n)] + 1L
    expect_false(any(x[-diag_pos] == 0), info = name)
    F2 <- Matrix::update(F, As, mult = -sigma)
    expect_identical(sum(as.numeric(methods::slot(F2, "nz"))), sum(as.numeric(nz)),
                     info = name)
    # Solves and later (indefinite) shifts are unchanged.
    set.seed(1)
    b <- matrix(stats::rnorm(n * 2), n, 2)
    M <- As - sigma * Matrix::Diagonal(n)
    xs <- as.matrix(Matrix::solve(F, b, system = "A"))
    expect_lt(max(abs(as.matrix(M %*% xs) - b)), 1e-10 * max(abs(b)))
    s2 <- stats::median(Matrix::diag(As)) + 0.123
    F3 <- Matrix::update(F, As, mult = -s2)
    F4 <- Matrix::Cholesky(As, LDL = TRUE, super = FALSE, perm = TRUE, Imult = -s2)
    expect_equal(as.matrix(Matrix::solve(F3, b, system = "A")),
                 as.matrix(Matrix::solve(F4, b, system = "A")), tolerance = 1e-8,
                 info = name)
  }
  # The relaxed zeros are gone on a 2-D grid: the converted factor is no
  # larger than CHOLMOD's simplicial LDL' with its own (AMD) ordering by more
  # than the ordering difference, and the native solve kernel accepts it.
  S <- rp_laplacian_2d(60)
  As <- methods::as(Matrix::forceSymmetric(S, uplo = "U"), "CsparseMatrix")
  F <- eigencore:::ldl_try_spd_factor(As, -0.5)
  Fs <- Matrix::Cholesky(As, LDL = TRUE, super = FALSE, perm = TRUE, Imult = 0.5)
  expect_lt(sum(as.numeric(methods::slot(F, "nz"))),
            1.1 * sum(as.numeric(methods::slot(Fs, "nz"))))
  spec <- eigencore:::ldl_native_solve_spec(F)
  kernel <- eigencore:::new_native_composite_kernel(spec)
  set.seed(3)
  X <- matrix(stats::rnorm(nrow(S) * 3), nrow(S), 3)
  expect_equal(eigencore:::native_composite_block_apply(kernel, X),
               as.matrix(Matrix::solve(Fs, X, system = "A")), tolerance = 1e-10)
  # End to end: the C60 smallest route on the converted factor.
  fit <- eig_partial(S, k = 4, target = nearest(-0.5), method = shift_invert(-0.5))
  expect_true(certificate(fit)$passed)
  expect_true(isTRUE(fit$transform$factorization_cache$positive_definite_factor))
  ev <- sort(as.vector(outer(4 * sin(pi * (1:60) / 122)^2, 4 * sin(pi * (1:60) / 122)^2, `+`)))
  expect_equal(sort(values(fit)), ev[1:4], tolerance = 1e-9)
})

rp_dense_with_spectrum <- function(d, seed) {
  set.seed(seed)
  n <- length(d)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  A <- Q %*% (d * t(Q))
  list(A = (A + t(A)) / 2, Q = Q)
}

rp_local_options <- function(opts, .env = parent.frame()) {
  old <- options(opts)
  do.call(on.exit, list(call("options", old), add = TRUE), envir = .env)
  invisible(old)
}

test_that("P25: the dense extremal count is proved by one deflated Cholesky", {
  rp_local_options(list(eigencore.completeness_cholesky_min_n = 100))
  n <- 300L
  d <- c(10, 9.5, 9, 8, 7.9, seq(7, -3, length.out = n - 5L))
  A <- rp_dense_with_spectrum(d, 2501L)$A
  fit <- eig_partial(A, 5L, largest(), seed = 1L)
  cert <- certificate(fit)
  expect_true(cert$passed)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_identical(cert$completeness$inertia_method, "dense_deflated_cholesky")
  expect_identical(cert$completeness$proof, "deflated_cholesky")
  expect_identical(cert$completeness$count_upper, 5L)
  expect_lt(cert$completeness$cholesky_delta, cert$completeness$margin)
  expect_equal(sort(values(fit)), sort(d)[n - 4:0], tolerance = 1e-8)
  fit <- eig_partial(A, 4L, smallest(), seed = 2L)
  cert <- certificate(fit)
  expect_identical(cert$target_completeness, "inertia_verified")
  expect_identical(cert$completeness$inertia_method, "dense_deflated_cholesky")
  expect_equal(sort(values(fit)), sort(d)[1:4], tolerance = 1e-8)
  # With the size cut-off at Inf the Bunch-Kaufman count is used.
  rp_local_options(list(eigencore.completeness_cholesky_min_n = Inf))
  fit <- eig_partial(A, 5L, largest(), seed = 1L)
  expect_identical(certificate(fit)$completeness$inertia_method, "dense_bunch_kaufman")
  expect_identical(certificate(fit)$target_completeness, "inertia_verified")
})

test_that("P25: the deflated Cholesky never verifies a wrong or unresolved set", {
  rp_local_options(list(eigencore.completeness_cholesky_min_n = 10))
  n <- 120L
  d <- c(9, 8, 7, 6, seq(5, 0.1, length.out = n - 4L))
  sp <- rp_dense_with_spectrum(d, 2502L)
  A <- sp$A
  Q <- sp$Q
  ctx <- eigencore:::inertia_context(A)
  op <- as_operator(A)
  check <- function(idx, values = d[idx], residuals = rep(1e-12, length(idx)),
                    target = largest()) {
    eigencore:::inertia_completeness_check(op, values, residuals, 0, target,
                                           ctx = ctx, vectors = Q[, idx, drop = FALSE])
  }
  ok <- check(1:3)
  expect_identical(ok$status, "inertia_verified")
  expect_identical(ok$record$proof, "deflated_cholesky")
  # A missing eigenvalue (8): H is indefinite, the count then reports failure.
  bad <- check(c(1, 3, 4))
  expect_identical(bad$status, "inertia_failed")
  expect_null(bad$record$proof)
  expect_match(bad$record$cholesky_fallback, "not positive definite")
  # Smallest end, wrong set.
  bad <- check(c(n, n - 2L), target = smallest())
  expect_identical(bad$status, "inertia_failed")
  # Residual bound wider than the gap: never verified by the Cholesky proof.
  wide <- check(1:3, residuals = rep(0.8, 3))
  expect_null(wide$record$proof)
  # Vectors that are not eigenvectors: the proof only uses null(X'), so a
  # garbage basis cannot verify the set.
  set.seed(9)
  junk <- eigencore:::inertia_completeness_check(
    op, d[1:3], rep(1e-12, 3), 0, largest(), ctx = ctx,
    vectors = matrix(stats::rnorm(n * 3), n, 3))
  expect_null(junk$record$proof)
  # Agreement with the count on random spectra and subsets.
  for (seed in 1:6) {
    set.seed(100 + seed)
    dd <- sort(stats::rnorm(n), decreasing = TRUE)
    sp2 <- rp_dense_with_spectrum(dd, 2600L + seed)
    ctx2 <- eigencore:::inertia_context(sp2$A)
    op2 <- as_operator(sp2$A)
    for (idx in list(1:5, c(1:4, 6), 2:6)) {
      run <- function(min_n) {
        rp_local_options(list(eigencore.completeness_cholesky_min_n = min_n))
        eigencore:::inertia_completeness_check(
          op2, dd[idx], rep(1e-12, length(idx)), 0, largest(), ctx = ctx2,
          vectors = sp2$Q[, idx, drop = FALSE])$status
      }
      expect_identical(run(10), run(Inf), info = paste(seed, toString(idx)))
    }
  }
})

test_that("P25: the structure check reports exact symmetry for the dense Lanczos", {
  fused <- eigencore:::dense_finite_symmetric
  set.seed(2503)
  for (n in c(1L, 2L, 63L, 64L, 65L, 130L)) {
    A <- crossprod(matrix(stats::rnorm(n * n), n))
    expect_identical(fused(A), c(TRUE, TRUE, TRUE))
    if (n > 1L) {
      B <- A
      B[n, 1L] <- B[n, 1L] * (1 + 1e-14)
      expect_identical(fused(B), c(TRUE, TRUE, FALSE))
      op <- as_operator(B)
      expect_identical(op$structure$kind, "hermitian")
      expect_false(eigencore:::dense_exact_symmetric_flag(op))
    }
    expect_true(eigencore:::dense_exact_symmetric_flag(as_operator(A)))
  }
  expect_identical(fused(matrix(c(1, NaN, NaN, 1), 2)), c(FALSE, FALSE, FALSE))
  # A derived operator has no memoised flag: the kernel checks itself.
  op <- as_operator(crossprod(matrix(stats::rnorm(400), 20)))
  expect_true(is.na(eigencore:::dense_exact_symmetric_flag(eigencore:::operator_scale(op, 2))))
  # Same values whichever way the flag arrives.
  A <- crossprod(matrix(stats::rnorm(160 * 160), 160))
  fit1 <- eig_partial(A, 4L, largest(), seed = 3L)
  fit2 <- eig_partial(eigencore:::operator_scale(as_operator(A), 1), 4L, largest(), seed = 3L)
  expect_equal(values(fit1), values(fit2), tolerance = 1e-10)
})

test_that("P26: the parallel work threshold sets team sizes, never results", {
  skip_if_not(isTRUE(eigencore:::native_thread_info()[["openmp"]] == 1L),
              "built without OpenMP")
  set.seed(2601)
  shapes <- list(
    square = Matrix::rsparsematrix(3000, 3000, density = 0.004),  # CSR gather
    tall = Matrix::rsparsematrix(9000, 500, density = 0.006),     # row slabs
    wide = Matrix::rsparsematrix(500, 7000, density = 0.006)
  )
  run_apply <- function(min_work, threads, A, X, transpose) {
    old <- options(eigencore.threads = threads, eigencore.parallel_min_work = min_work)
    on.exit(options(old))
    eigencore:::csc_apply_repeat(A, X, transpose = transpose, reps = 3L,
                                 col_means = Matrix::colMeans(A),
                                 weights = rep(1.5, ncol(A)))$Y
  }
  for (name in names(shapes)) {
    A <- shapes[[name]]
    for (b in c(1L, 12L)) {
      for (transpose in c(FALSE, TRUE)) {
        X <- matrix(rnorm((if (transpose) nrow(A) else ncol(A)) * b), ncol = b)
        ref <- run_apply(0, 4L, A, X, transpose)
        for (mw in list(Inf, 1e5, NULL)) {
          expect_identical(run_apply(mw, 4L, A, X, transpose), ref,
                           info = paste(name, b, transpose, format(mw)))
        }
        expect_identical(run_apply(NULL, 1L, A, X, transpose), ref)
      }
    }
  }
  # Solves: CSC Lanczos (slab reorthogonalisation, n > 1024 rows), block
  # Lanczos, implicit Gram SVD and Golub-Kahan agree bitwise for every
  # threshold and every thread count > 1.
  B <- Matrix::rsparsematrix(5000, 5000, density = 8e-4)
  S <- methods::as(methods::as(Matrix::forceSymmetric(B + Matrix::t(B)),
                               "generalMatrix"), "CsparseMatrix")
  R <- shapes$tall
  solve_all <- function(threads, min_work) {
    old <- options(eigencore.threads = threads, eigencore.parallel_min_work = min_work)
    on.exit(options(old))
    list(
      e = eig_partial(S, 6L, seed = 1),
      b = eig_partial(S, 6L, method = lanczos(block = 3L), seed = 2),
      s = svd_partial(R, 5L, seed = 3),
      g = svd_partial(R, 4L, method = golub_kahan(), seed = 4)
    )
  }
  strip <- function(r) lapply(r, function(f) list(values(f), f$vectors %||% f$u))
  `%||%` <- function(a, b) if (is.null(a)) b else a
  ref <- strip(solve_all(4L, 0))
  expect_identical(strip(solve_all(4L, Inf)), ref)
  expect_identical(strip(solve_all(4L, NULL)), ref)
  expect_identical(strip(solve_all(2L, NULL)), ref)
  one <- strip(solve_all(1L, NULL))
  for (i in seq_along(ref)) {
    expect_equal(one[[i]][[1L]], ref[[i]][[1L]], tolerance = 1e-10)
  }
})

test_that("P26: invalid parallel_min_work values fall back to the default", {
  set.seed(2602)
  A <- Matrix::rsparsematrix(2000, 2000, density = 0.004)
  X <- matrix(rnorm(2000 * 3), ncol = 3)
  ref <- eigencore:::csc_apply_repeat(A, X, reps = 2L)$Y
  for (bad in list(-1, NA_real_, "x", NA_integer_, c(0, 1))) {
    old <- options(eigencore.threads = 4L, eigencore.parallel_min_work = bad)
    out <- tryCatch(eigencore:::csc_apply_repeat(A, X, reps = 2L)$Y,
                    finally = options(old))
    expect_identical(out, ref)
  }
})
