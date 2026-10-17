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

# P27: a dense ring of eigenvalue moduli at the wanted edge. Krylov-Schur
# converges ring eigenvalues in an order set by their angular isolation, not
# by moduli that differ by ~1e-5, so before P27 the main solve returned a set
# missing several pairs and the single complement probe certified it.
nonsym_ring_matrix <- function(n, npairs, width, seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  T <- matrix(0, n, n)
  vals <- complex(0)
  j <- 1L
  for (p in seq_len(npairs)) {
    r <- 1 - width * stats::runif(1)
    phi <- stats::runif(1, 0.05, pi - 0.05)
    T[j:(j + 1L), j:(j + 1L)] <- matrix(c(r * cos(phi), -r * sin(phi),
                                          r * sin(phi), r * cos(phi)), 2)
    vals <- c(vals, r * exp(1i * phi), r * exp(-1i * phi))
    j <- j + 2L
  }
  rest <- n - j + 1L
  d <- 0.95 * sqrt(stats::runif(rest)) * sample(c(-1, 1), rest, TRUE)
  T[cbind(j:n, j:n)] <- d
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  list(A = Q %*% T %*% t(Q), values = c(vals, d))
}

test_that("a dense ring at the edge is repaired or declined, never certified wrong (P27)", {
  ring <- nonsym_ring_matrix(600, 60, 2e-3, seed = 2)
  truth <- sort(Mod(ring$values), decreasing = TRUE)[1:6]
  for (s in c(1L, 3L)) {
    fit <- eig_partial(ring$A, k = 6, target = largest_magnitude(), seed = s,
                       left_vectors = "none")
    cert <- certificate(fit)
    correct <- isTRUE(all.equal(sort(Mod(values(fit)), decreasing = TRUE), truth,
                                tolerance = 1e-7))
    expect_true(correct || !isTRUE(cert$passed))
    if (!correct) {
      expect_false(identical(cert$target_completeness, "probed"))
    }
    if (identical(cert$target_completeness, "inconclusive")) {
      expect_match(cert$completeness$reason, "boundary unresolved")
    }
    expect_false(is.null(cert$completeness$crowding_gap))
  }
})

test_that("a well-separated edge runs a single probe pass (P27)", {
  set.seed(11)
  n <- 400
  d <- c(10, 9, 8, 7, 6, 5, stats::runif(n - 6, -1, 1))
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
  T <- diag(d)
  T[upper.tri(T)] <- 0.1 * stats::rnorm(n * (n - 1) / 2) / sqrt(n)
  A <- Q %*% T %*% t(Q)
  fit <- eig_partial(A, k = 3, target = largest_real(), left_vectors = "none")
  cert <- certificate(fit)
  expect_true(cert$passed)
  expect_identical(cert$target_completeness, "probed")
  expect_null(cert$completeness$boundary_passes)
  expect_equal(sort(values(fit)), c(8, 9, 10), tolerance = 1e-8)
})

test_that("boundary helpers: edge spacing and probe agreement (P27)", {
  ns <- asNamespace("eigencore")
  # Conjugate partners share a key and collapse; spacing is relative.
  expect_equal(ns$nonsym_edge_spacing(c(2, 2, 1.9, 1.9, 1.7), take = 10,
                                      collapse = 1e-12, scale = 2),
               stats::median(c(0.1, 0.2)) / 2)
  expect_true(is.na(ns$nonsym_edge_spacing(c(1, 1), 5, 1e-12, 1)))
  a <- list(key = c(1, 0.9), res = c(1e-10, 1e-10), converged_values = c(TRUE, TRUE))
  b <- list(key = c(0.95, 0.9), res = c(1e-10, 1e-10), converged_values = c(TRUE, TRUE))
  expect_true(ns$nonsym_probe_agree(a, a, resolution = 1e-4, margin = 0))
  expect_false(ns$nonsym_probe_agree(a, b, resolution = 1e-4, margin = 0))
  expect_false(ns$nonsym_probe_agree(a, list(key = 1, res = 1, converged_values = FALSE),
                                     resolution = 1e-4, margin = 0))
})
