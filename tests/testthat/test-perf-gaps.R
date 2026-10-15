# Single-thread performance-gap fixes (docs/review-2026-10.md, P19-P24):
# each test pins the cheaper path and that its answer and certificate are
# unchanged.

test_that("sparse shift-invert nearest() counts with two factorisations (P19)", {
  g <- 100L
  T1 <- Matrix::bandSparse(g, k = c(-1, 0, 1),
                           diagonals = list(rep(-1, g - 1), rep(2, g), rep(-1, g - 1)))
  I1 <- Matrix::Diagonal(g)
  A <- methods::as(Matrix::kronecker(T1, I1) + Matrix::kronecker(I1, T1), "CsparseMatrix")
  sigma <- 4.01
  fit <- eig_partial(A, k = 6, target = nearest(sigma), tol = 1e-8, seed = 1L)
  cert <- certificate(fit)
  expect_true(cert$passed)
  expect_identical(cert$target_completeness, "inertia_verified")
  # The margin starts above the shift-invert factor's backward bound, so the
  # count needs no widening recount: one LDL' per side of the window.
  expect_equal(cert$completeness$factorizations, 2)
  t1 <- 2 - 2 * cos(seq_len(g) * pi / (g + 1))
  ev <- as.vector(outer(t1, t1, `+`))
  expect_equal(sort(values(fit)), sort(ev[order(abs(ev - sigma))[1:6]]),
               tolerance = 1e-10)
})

test_that("shift-invert cache keys use the native digest (P19)", {
  A <- Matrix::rsparsematrix(40, 40, 0.1, symmetric = TRUE)
  key <- eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.5)
  B <- A
  B@x[1] <- B@x[1] + 1
  expect_identical(key, eigencore:::shift_invert_factorization_cache_key(as_operator(A), 0.5))
  expect_false(identical(key, eigencore:::shift_invert_factorization_cache_key(as_operator(B), 0.5)))
})

test_that("exactly symmetric dense sources use the one-triangle apply (P21)", {
  set.seed(11)
  X <- matrix(rnorm(400 * 300), 400)
  A <- crossprod(X) / 400
  expect_true(isSymmetric(A, tol = 0))
  fit <- eig_partial(A, k = 5, target = largest(), tol = 1e-10)
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  expect_true(certificate(fit)$passed)
  expect_equal(values(fit), ev[1:5], tolerance = 1e-10)
  # A source symmetric only within the structure tolerance keeps the full
  # (dgemv) product and gives the same answer.
  B <- A
  B[2, 1] <- B[2, 1] * (1 + 1e-15)
  fitB <- eig_partial(B, k = 5, target = largest(), tol = 1e-10)
  expect_true(certificate(fitB)$passed)
  expect_equal(values(fitB), ev[1:5], tolerance = 1e-10)
  # Tiny exact-symmetry helper cases through the dense kernel.
  for (n in c(1L, 2L, 65L)) {
    S <- crossprod(matrix(rnorm(n * n), n))
    f <- eig_partial(S, k = 1, target = largest())
    expect_equal(values(f), max(eigen(S, symmetric = TRUE, only.values = TRUE)$values),
                 tolerance = 1e-8)
  }
})

test_that("the dense inertia context uses the LAPACK 1-norm (P21)", {
  set.seed(12)
  A <- crossprod(matrix(rnorm(50 * 40), 50))
  ctx <- eigencore:::inertia_context(A)
  expect_equal(ctx$normA, max(colSums(abs(A))), tolerance = 1e-14)
})

test_that("identity hashes are unchanged by the bulk double loop (P22)", {
  # Lengths around the 4-word lane boundary, attribute-carrying doubles and
  # a mixed list: hashing word by word or in whole rounds must agree, which
  # the two calls below check through the boundary-sensitive prefix sums.
  set.seed(13)
  x <- rnorm(1031)
  hs <- vapply(0:9, function(i) eigencore:::stable_raw_hash(x[seq_len(1021 + i)]), "")
  expect_identical(length(unique(hs)), 10L)
  expect_identical(eigencore:::stable_raw_hash(c(-0, NaN, NA, 1)),
                   eigencore:::stable_raw_hash(c(0, NaN, NA, 1)))
  expect_false(identical(eigencore:::stable_raw_hash(c(NaN, 1)),
                         eigencore:::stable_raw_hash(c(NA, 1))))
  m <- matrix(x[1:1024], 32)
  expect_identical(eigencore:::stable_raw_hash(m), eigencore:::stable_raw_hash(m + 0))
})

test_that("a screen-converged intruder repairs without the full-tolerance probe (P23)", {
  set.seed(3)
  A <- matrix(rnorm(120^2), 120) / sqrt(120)
  op <- eigencore:::as_operator(A)
  e <- eigen(A)
  idx <- order(Re(e$values), decreasing = TRUE)[2:5]
  cert_op <- eigencore:::arnoldi_certificate_operator(op)
  certify_fn <- function(values, vectors)
    eigencore:::certify_general_eigen_operator(cert_op, values, vectors, tol = 1e-8)
  ctl <- modifyList(eigencore:::nonsym_completeness_controls(), list(exhaust = 0L))
  fixed <- eigencore:::nonsym_completeness_check(
    op, e$values[idx], e$vectors[, idx], largest_real(), tol = 1e-8,
    certify_fn = certify_fn, controls = ctl)
  # Screening disabled (screen tolerance = solve tolerance): the old path.
  slow <- eigencore:::nonsym_completeness_check(
    op, e$values[idx], e$vectors[, idx], largest_real(), tol = 1e-8,
    certify_fn = certify_fn, controls = modifyList(ctl, list(screen_tol = 1e-8)))
  expect_identical(fixed$status, "repaired")
  expect_identical(slow$status, "repaired")
  expect_equal(sort(Re(fixed$values)), sort(Re(slow$values)), tolerance = 1e-8)
  expect_true(isTRUE(fixed$certificate$passed))
  expect_lte(fixed$record$operator_columns, slow$record$operator_columns)
})

test_that("sparse Lanczos uses the ARPACK subspace for cheap applies (P24)", {
  A <- Matrix::rsparsematrix(3000, 3000, 0.002, symmetric = TRUE)
  plan <- plan_solver(eigen_problem(A, target = largest()), k = 10)
  expect_identical(plan$controls$max_subspace, 21L)
  D <- crossprod(matrix(rnorm(300 * 200), 300))
  plan_d <- plan_solver(eigen_problem(D, target = largest()), k = 10)
  expect_identical(plan_d$controls$max_subspace, 50L)
  old <- options(eigencore.lanczos_sparse_row_nnz = 0)
  on.exit(options(old))
  plan_o <- plan_solver(eigen_problem(A, target = largest()), k = 10)
  expect_identical(plan_o$controls$max_subspace, 50L)
})
