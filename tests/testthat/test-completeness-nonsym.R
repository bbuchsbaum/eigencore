# Nonsymmetric target-completeness probe (R/completeness_nonsym.R).

top_by <- function(ev, target, k) ev[eigencore:::order_indices(ev, target)[seq_len(k)]]

test_that("a Krylov-Schur result that misses a largest-magnitude pair is repaired", {
  set.seed(1)
  A <- matrix(rnorm(300^2), 300) / sqrt(300)
  ev <- eigen(A, only.values = TRUE)$values
  fit <- eig_partial(A, k = 6, target = largest_magnitude())
  cert <- certificate(fit)
  expect_true(cert$passed)
  expect_true(cert$residual_passed)
  expect_true(cert$target_completeness %in% c("probed", "repaired"))
  expect_equal(sort(Mod(values(fit))), sort(Mod(top_by(ev, largest_magnitude(), 6))),
               tolerance = 1e-8)
  expect_identical(cert$completeness$method, "deflated_krylov_schur_probe")
})

test_that("the probe never consumes the caller's random stream", {
  set.seed(2)
  A <- matrix(rnorm(200^2), 200)
  states <- lapply(c("none", "auto"), function(mode) {
    old <- options(eigencore.target_completeness = mode)
    on.exit(options(old))
    set.seed(99)
    fit <- eig_partial(A, k = 4, target = largest_real())
    list(seed = .Random.seed, status = certificate(fit)$target_completeness)
  })
  expect_identical(states[[1]]$seed, states[[2]]$seed)
  expect_identical(states[[1]]$status, "not_checked")
  expect_true(states[[2]]$status %in% c("probed", "repaired"))
})

test_that("an intruder found by the probe and not repaired withholds passed", {
  set.seed(3)
  A <- matrix(rnorm(120^2), 120) / sqrt(120)
  op <- eigencore:::as_operator(A)
  e <- eigen(A)
  # A deliberately wrong set: the 2nd..5th largest-real eigenpairs.
  idx <- order(Re(e$values), decreasing = TRUE)[2:5]
  vals <- e$values[idx]
  vecs <- e$vectors[, idx]
  check <- eigencore:::nonsym_completeness_check(
    op, vals, vecs, largest_real(), tol = 1e-8,
    certify_fn = function(values, vectors) list(passed = FALSE),
    controls = modifyList(eigencore:::nonsym_completeness_controls(),
                          list(exhaust = 0L)))
  expect_identical(check$status, "failed")
  expect_true(check$record$intruder_found)
  expect_false(check$repaired)
  # With a real certificate the merge recovers the true set.
  cert_op <- eigencore:::arnoldi_certificate_operator(op)
  fixed <- eigencore:::nonsym_completeness_check(
    op, vals, vecs, largest_real(), tol = 1e-8,
    certify_fn = function(values, vectors)
      eigencore:::certify_general_eigen_operator(cert_op, values, vectors, tol = 1e-8))
  expect_identical(fixed$status, "repaired")
  expect_equal(sort(Re(fixed$values)), sort(Re(top_by(e$values, largest_real(), 4))),
               tolerance = 1e-8)
})

test_that("a conjugate partner more preferred than the edge is an intruder", {
  # Spectrum 3 +- 2i, 1, 0.5, ...: for largest_imaginary k = 2 the set
  # {3 - 2i, 1} is wrong although 3 - 2i's partner is in span(Re x, Im x).
  set.seed(4)
  D <- diag(c(1, 0.5, -1, -2, 0))
  D[1:2, 1:2] <- matrix(c(3, -2, 2, 3), 2)
  P <- qr.Q(qr(matrix(rnorm(25), 5)))
  A <- P %*% D %*% t(P)
  e <- eigen(A)
  pick <- c(which.min(Im(e$values)), which.min(Mod(e$values - 1)))
  check <- eigencore:::nonsym_completeness_check(
    eigencore:::as_operator(A), e$values[pick], e$vectors[, pick],
    largest_imaginary(), tol = 1e-8,
    certify_fn = function(values, vectors) list(passed = FALSE))
  expect_true(check$record$intruder_found)
  expect_identical(check$status, "failed")
})

test_that("linearly dependent returned vectors are inconclusive, never verified", {
  set.seed(5)
  A <- matrix(rnorm(60^2), 60)
  e <- eigen(A)
  i <- which.max(Re(e$values))
  check <- eigencore:::nonsym_completeness_check(
    eigencore:::as_operator(A), rep(e$values[i], 2), cbind(e$vectors[, i], e$vectors[, i]),
    largest_real(), tol = 1e-8,
    certify_fn = function(values, vectors) list(passed = FALSE))
  expect_identical(check$status, "inconclusive")
})

test_that("shift-invert and sparse general-pencil Arnoldi results are probed", {
  set.seed(6)
  A <- matrix(rnorm(150^2), 150)
  fit <- eig_partial(A, 4, target = nearest(0.3))
  expect_true(certificate(fit)$target_completeness %in% c("probed", "repaired"))
  expect_true(certificate(fit)$passed)
  ev <- eigen(A, only.values = TRUE)$values
  expect_equal(sort(Mod(values(fit) - 0.3)), sort(Mod(ev - 0.3))[1:4], tolerance = 1e-8)

  S <- Matrix::rsparsematrix(400, 400, 0.02) + Matrix::Diagonal(400, 1)
  B <- Matrix::Diagonal(x = stats::runif(400, 1, 2))
  pfit <- eig_partial(S, B = B, k = 3, target = largest_real())
  expect_true(certificate(pfit)$target_completeness %in% c("probed", "repaired"))
  pev <- eigen(solve(as.matrix(B), as.matrix(S)), only.values = TRUE)$values
  expect_equal(sort(Re(values(pfit))), sort(Re(top_by(pev, largest_real(), 3))),
               tolerance = 1e-7)
})

test_that("completeness = 'none' leaves nonsymmetric results not_checked", {
  set.seed(7)
  A <- matrix(rnorm(100^2), 100)
  old <- options(eigencore.target_completeness = "none")
  on.exit(options(old))
  fit <- eig_partial(A, 3, target = largest_real())
  expect_identical(certificate(fit)$target_completeness, "not_checked")
  expect_false(certificate(fit)$passed)
  expect_true(certificate(fit)$residual_passed)
})
