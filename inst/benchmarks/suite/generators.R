# Benchmark suite: deterministic problem generators.
#
# Every generator takes explicit size parameters and a seed, sets the RNG
# itself (Mersenne-Twister / Inversion / Rejection), and returns a list with
#
#   A            the matrix (dgCMatrix for sparse families, base matrix for dense)
#   B            optional SPD metric (generalized families), else NULL
#   center       TRUE when the problem is the column-centred operator A - 1 mu'
#   exact        optional numeric vector holding the *complete* exact spectrum
#                (analytic families); used as the trusted reference
#   description  one-line human description (stored in cases.csv)
#
# Bump SUITE_GENERATOR_VERSION whenever a generator changes its output: it is
# part of the reference cache key.

SUITE_GENERATOR_VERSION <- 1L

as_dgc <- function(A) {
  methods::as(methods::as(methods::as(A, "CsparseMatrix"), "generalMatrix"),
              "dMatrix")
}

# Random sparse symmetric matrix, ~nnz_row nonzeros per row, N(0,1) entries.
gen_sparse_sym <- function(n, nnz_row = 5, seed) {
  suite_set_seed(seed)
  A <- Matrix::rsparsematrix(n, n, density = nnz_row / n, symmetric = TRUE)
  list(A = as_dgc(A),
       description = sprintf("random sparse symmetric, n=%d, ~%g nnz/row, N(0,1) values", n, nnz_row))
}

# Random sparse nonsymmetric square matrix.
gen_sparse_nonsym <- function(n, nnz_row = 5, seed) {
  suite_set_seed(seed)
  A <- Matrix::rsparsematrix(n, n, density = nnz_row / n)
  list(A = as_dgc(A),
       description = sprintf("random sparse nonsymmetric, n=%d, ~%g nnz/row, N(0,1) values", n, nnz_row))
}

# Random sparse rectangular matrix (SVD families).
gen_sparse_rect <- function(m, n, density, seed, positive = FALSE) {
  suite_set_seed(seed)
  rx <- if (positive) function(k) stats::rexp(k) else function(k) stats::rnorm(k)
  A <- Matrix::rsparsematrix(m, n, density = density, rand.x = rx)
  list(A = as_dgc(A),
       description = sprintf("random sparse %d x %d, density %g, %s values", m, n,
                             density, if (positive) "Exp(1)" else "N(0,1)"))
}

# 1-D Dirichlet Laplacian (banded, tridiagonal). Exact spectrum known.
gen_laplacian_1d <- function(n, seed = NULL) {
  A <- Matrix::bandSparse(n, n, k = c(-1L, 0L, 1L),
                          diagonals = list(rep(-1, n - 1L), rep(2, n), rep(-1, n - 1L)))
  j <- seq_len(n)
  list(A = as_dgc(A), exact = 2 - 2 * cos(j * pi / (n + 1)),
       description = sprintf("1-D Dirichlet Laplacian (tridiagonal), n=%d; analytic spectrum", n))
}

lap1d_parts <- function(g) {
  Matrix::bandSparse(g, g, k = c(-1L, 0L, 1L),
                     diagonals = list(rep(-1, g - 1L), rep(2, g), rep(-1, g - 1L)))
}

# 2-D 5-point Dirichlet Laplacian on a g x g grid (n = g^2). Exact spectrum
# 4 sin^2(i pi / 2(g+1)) + 4 sin^2(j pi / 2(g+1)); eigenvalues are repeated
# (pairs from the grid symmetry, large multiplicity at 4).
gen_laplacian_2d <- function(g, seed = NULL) {
  T <- lap1d_parts(g)
  I <- Matrix::Diagonal(g)
  A <- Matrix::kronecker(T, I) + Matrix::kronecker(I, T)
  s <- 4 * sin(seq_len(g) * pi / (2 * (g + 1)))^2
  list(A = as_dgc(A), exact = as.vector(outer(s, s, "+")),
       description = sprintf("2-D 5-point Dirichlet Laplacian, %d x %d grid (n=%d); analytic spectrum, repeated eigenvalues", g, g, g * g))
}

# Generalized SPD pencil from bilinear finite elements on a g x g grid:
# A = K (x) M + M (x) K (stiffness), B = M (x) M (mass), with
# K = tridiag(-1, 2, -1), M = tridiag(1/6, 2/3, 1/6). Both share the sine
# basis, so lambda_ij = kappa_i / mu_i + kappa_j / mu_j exactly.
gen_fem_pencil_2d <- function(g, seed = NULL) {
  K <- lap1d_parts(g)
  M <- Matrix::bandSparse(g, g, k = c(-1L, 0L, 1L),
                          diagonals = list(rep(1 / 6, g - 1L), rep(2 / 3, g), rep(1 / 6, g - 1L)))
  A <- Matrix::kronecker(K, M) + Matrix::kronecker(M, K)
  B <- Matrix::kronecker(M, M)
  th <- seq_len(g) * pi / (g + 1)
  r <- (2 - 2 * cos(th)) / (2 / 3 + cos(th) / 3)
  list(A = as_dgc(A), B = as_dgc(B), exact = as.vector(outer(r, r, "+")),
       description = sprintf("generalized SPD FEM pencil (Q1 stiffness/mass), %d x %d grid (n=%d); analytic spectrum", g, g, g * g))
}

# Random banded symmetric matrix (half-bandwidth bw, N(0,1) entries).
gen_banded_sym <- function(n, bw, seed) {
  suite_set_seed(seed)
  ks <- 0:bw
  diags <- lapply(ks, function(k) stats::rnorm(n - k))
  U <- Matrix::bandSparse(n, n, k = ks, diagonals = diags)
  A <- U + Matrix::t(U) - Matrix::Diagonal(n, x = diags[[1L]])
  list(A = as_dgc(A),
       description = sprintf("random banded symmetric, n=%d, half-bandwidth %d", n, bw))
}

# Dense symmetric with a prescribed spectrum: Q diag(values) Q', Q Haar-ish.
dense_with_spectrum <- function(values, seed) {
  n <- length(values)
  suite_set_seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), n, n)))
  A <- Q %*% (values * t(Q))
  (A + t(A)) / 2
}

# Dense Wishart-type symmetric matrix crossprod(X)/n (Marchenko-Pastur bulk).
gen_dense_wishart <- function(n, seed) {
  suite_set_seed(seed)
  X <- matrix(stats::rnorm(n * n), n, n)
  list(A = crossprod(X) / n,
       description = sprintf("dense symmetric crossprod(X)/n, X %d x %d N(0,1)", n, n))
}

# Power-law spectrum lambda_j = j^-alpha (random signs off for LA targets).
gen_powerlaw_dense <- function(n, alpha = 1, seed) {
  values <- seq_len(n)^(-alpha)
  list(A = dense_with_spectrum(values, seed), exact = values,
       description = sprintf("dense symmetric Q diag(j^-%g) Q', n=%d; analytic spectrum", alpha, n))
}

# Clustered / repeated top eigenvalues: 1 (x3), 1 - 1e-6 j (j=1..3),
# 0.9 (x2), 0.85, 0.8, then a bulk uniform on [0, 0.7].
gen_clustered_dense <- function(n, seed) {
  top <- c(1, 1, 1, 1 - 1e-6, 1 - 2e-6, 1 - 3e-6, 0.9, 0.9, 0.85, 0.8)
  suite_set_seed(seed + 1L)
  values <- c(top, sort(stats::runif(n - length(top), 0, 0.7), decreasing = TRUE))
  list(A = dense_with_spectrum(values, seed), exact = values,
       description = sprintf("dense symmetric, n=%d, top spectrum 1 (x3), 1-1e-6..3e-6, 0.9 (x2), 0.85, 0.8, bulk U(0,0.7); analytic", n))
}

# Low-rank plus noise (dense): U diag(s) V' + noise * N(0,1).
gen_lowrank_noise <- function(m, n, rank, noise, seed) {
  suite_set_seed(seed)
  U <- qr.Q(qr(matrix(stats::rnorm(m * rank), m, rank)))
  V <- qr.Q(qr(matrix(stats::rnorm(n * rank), n, rank)))
  s <- 10 * seq(rank, 1, length.out = rank) / rank
  A <- U %*% (s * t(V)) + noise * matrix(stats::rnorm(m * n), m, n)
  list(A = A,
       description = sprintf("dense %d x %d, rank-%d signal (s = 10..%.2g) + %g N(0,1) noise", m, n, rank, min(s), noise))
}
