# Oracle-based differential testing harness.
#
# Every case is a pure function of an integer case id: `oracle_case(id)`
# draws the problem family, storage, structure, spectrum, size, k, target,
# method and tolerance from a private RNG seeded by the id (the global RNG
# stream is saved and restored), `oracle_build()` materialises the input
# object together with its dense ground truth, and `oracle_run_case()` solves
# it and checks the result against base eigen()/svd() (or a dense Cholesky
# reduction for SPD pencils). A failing case is re-run in isolation with
#
#   source("tests/testthat/helper-oracle.R"); library(eigencore)
#   str(oracle_run_case(<id>))
#
# Invariants (see docs/test-assurance.md):
#   hard (a test failure):
#     crash / hang / unexpected error on a valid input;
#     certificate$passed but the harness's own normwise backward error
#       ||r|| / ((||A||_2 + |lambda| ||B||_2) ||x||) (SVD: ||[Av-su; A'u-sv]|| /
#       ||A||_2, with ||A||_2 from a dense SVD) exceeds tol;
#     reported backward error below the true one (C12 promises an over-estimate);
#     certified vectors not orthonormal (Hermitian, B-orthonormal for pencils,
#       U and V for SVD);
#     a certified value farther from every true eigenvalue than the residual
#       bound allows (Hermitian/SVD: ||r||/||x||; pencils: ||r||_{B^-1}/||x||_B;
#       nonsymmetric: Bauer-Fike / Jordan-block bound with the generator's V);
#     target_completeness in exact/probed/repaired/inertia_verified but the
#       returned multiset differs from the oracle's target set;
#     API: passed with fewer than k values, nconv > k, wrong dims, values not
#       ordered per target.
#   soft (counted and reported):
#     certified but the set differs from the oracle target set while the
#       certificate makes no completeness claim; uncertified results (with
#       their accuracy); expected "unsupported" errors.

oracle_or <- function(x, y) if (is.null(x)) y else x
`%||%` <- function(x, y) if (is.null(x)) y else x

# num / den with 0/0 = 0 and x/0 = Inf (zero operators).
oracle_div <- function(num, den) {
  out <- num / den
  out[den == 0 & num == 0] <- 0
  out[den == 0 & num != 0] <- Inf
  out
}

# ---------------------------------------------------------------------------
# Levels
# ---------------------------------------------------------------------------

oracle_level <- function() {
  level <- Sys.getenv("EIGENCORE_ORACLE_LEVEL", "")
  if (!nzchar(level)) {
    level <- if (identical(Sys.getenv("NOT_CRAN"), "true")) "ci" else "cran"
  }
  level <- tolower(level)
  if (!level %in% c("cran", "ci", "extended")) {
    stop("EIGENCORE_ORACLE_LEVEL must be cran, ci or extended.", call. = FALSE)
  }
  level
}

# Case ids per level. Levels are nested prefixes of one id stream, and the
# size cap depends on the id only: "cran" is ids 1..40 (n <= 40),
# "ci" 1..400 (n <= 160), "extended" 1..15000 (n <= 300). The id alone
# reproduces a case.
oracle_level_ids <- function(level = oracle_level()) {
  switch(level,
         cran = seq_len(40L),
         ci = seq_len(400L),
         extended = seq_len(15000L))
}

oracle_level_max_n <- function(level) {
  switch(level, cran = 40L, ci = 160L, extended = 300L)
}

# ---------------------------------------------------------------------------
# Private RNG
# ---------------------------------------------------------------------------

oracle_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}

oracle_pick <- function(x, prob = NULL) {
  x[[sample.int(length(x), 1L, prob = prob)]]
}

# ---------------------------------------------------------------------------
# Case descriptors
# ---------------------------------------------------------------------------

# Family schedule: a 20-slot cycle so every family appears in any 20
# consecutive ids (cran level covers all of them).
oracle_family_cycle <- c(
  "herm", "svd", "nonsym", "herm", "shim", "svd", "gen", "herm", "complex",
  "svd", "herm", "nonsym", "shim", "herm", "svd", "gen", "herm", "nonsym",
  "svd", "shim"
)

oracle_sizes <- c(5L, 6L, 8L, 12L, 20L, 30L, 50L, 80L, 120L, 200L, 300L)
oracle_size_prob <- c(4, 3, 4, 4, 4, 4, 4, 3, 2, 1, 0.6)

oracle_draw_n <- function(max_n, min_n = 5L) {
  ok <- oracle_sizes >= min_n & oracle_sizes <= max_n
  oracle_pick(oracle_sizes[ok], oracle_size_prob[ok])
}

oracle_draw_k <- function(n, kmax = n) {
  kind <- oracle_pick(c("one", "small", "near_n"), c(3, 5, 2))
  k <- switch(kind,
              one = 1L,
              small = sample.int(max(1L, min(6L, kmax - 1L)), 1L),
              near_n = max(1L, kmax - sample(0:2, 1L)))
  as.integer(min(max(k, 1L), kmax))
}

oracle_draw_tol <- function() {
  oracle_pick(c(1e-8, 1e-10, 1e-6, 1e-12), c(6, 2, 1.5, 0.5))
}

oracle_case <- function(id, level = NULL) {
  id <- as.integer(id)
  max_n <- oracle_level_max_n(oracle_or(level, if (id <= 40L) "cran" else
    if (id <= 400L) "ci" else "extended"))
  oracle_with_seed(1000003L * 7L + id, {
    family <- oracle_family_cycle[[(id - 1L) %% length(oracle_family_cycle) + 1L]]
    case <- switch(family,
                   herm = oracle_case_herm(max_n),
                   nonsym = oracle_case_nonsym(max_n),
                   complex = oracle_case_complex(max_n),
                   gen = oracle_case_gen(max_n),
                   svd = oracle_case_svd(max_n),
                   shim = oracle_case_shim(max_n))
    case$id <- id
    case$family <- family
    case$seed <- 10000L + id
    case
  })
}

oracle_herm_spectra <- c("random", "clustered", "repeated", "near_repeated",
                         "power_law", "tiny", "huge", "indefinite", "singular",
                         "very_tiny")

oracle_case_herm <- function(max_n) {
  storage <- oracle_pick(
    c("dense", "dgC", "dsC_U", "dsC_L", "dgT", "dgR", "ddi", "dge", "dsy",
      "linop", "composite"),
    c(4, 3, 2, 2, 1, 1, 1, 1, 1, 2, 2))
  n <- oracle_draw_n(max_n)
  target <- oracle_pick(
    c("largest", "smallest", "largest_magnitude", "smallest_magnitude",
      "both_ends", "nearest", "largest_real", "smallest_real"),
    c(4, 3, 2, 2, 2, 2, 0.5, 0.5))
  method <- oracle_pick(
    c("auto", "lanczos1", "lanczos2", "lanczos3", "lobpcg", "shift_invert"),
    c(5, 2, 1.5, 1, 1.5, 1))
  # shift_invert(sigma) means "nearest sigma"; other targets are an error
  if (method == "shift_invert" && stats::runif(1) < 0.85) target <- "nearest"
  list(api = "eig", structure = "herm", storage = storage, n = n,
       spectrum = oracle_pick(oracle_herm_spectra),
       composite = if (storage == "composite")
         oracle_pick(c("crossprod", "crossprod_center", "crossprod_scale",
                       "crossprod_sparse")),
       k = oracle_draw_k(n), target = target, method = method,
       sigma_pos = stats::runif(1), tol = oracle_draw_tol())
}

oracle_case_nonsym <- function(max_n) {
  storage <- oracle_pick(c("dense", "dgC", "dgT", "dgR", "dge", "linop",
                           "composite"), c(5, 3, 1, 1, 1, 2, 1))
  n <- oracle_draw_n(min(max_n, 200L))
  target <- oracle_pick(
    c("largest_magnitude", "smallest_magnitude", "largest_real",
      "smallest_real", "largest_imaginary", "smallest_imaginary", "nearest"),
    c(4, 2, 2, 2, 1, 1, 2))
  method <- if (target == "nearest") {
    oracle_pick(c("auto", "shift_invert"), c(1, 1))
  } else {
    oracle_pick(c("auto", "lanczos1"), c(10, 1))
  }
  list(api = "eig", structure = "nonsym", storage = storage, n = n,
       spectrum = oracle_pick(c("random_matrix", "real", "complex_pairs",
                                "clustered", "defective", "tiny", "huge",
                                "singular", "repeated")),
       k = oracle_draw_k(n), target = target, method = method,
       sigma_pos = stats::runif(1), tol = oracle_draw_tol())
}

oracle_case_complex <- function(max_n) {
  n <- oracle_draw_n(min(max_n, 120L))
  herm <- stats::runif(1) < 0.5
  target <- if (herm) {
    oracle_pick(c("largest", "smallest", "largest_magnitude",
                  "smallest_magnitude", "both_ends", "nearest"))
  } else {
    oracle_pick(c("largest_magnitude", "smallest_magnitude", "largest_real",
                  "smallest_real", "largest_imaginary", "smallest_imaginary",
                  "nearest"))
  }
  list(api = "eig", structure = if (herm) "cherm" else "cgen",
       storage = "dense", n = n,
       spectrum = if (herm) oracle_pick(oracle_herm_spectra) else
         oracle_pick(c("random_matrix", "complex", "repeated", "tiny", "huge")),
       k = oracle_draw_k(n), target = target, method = "auto",
       sigma_pos = stats::runif(1), tol = oracle_draw_tol())
}

oracle_case_gen <- function(max_n) {
  n <- oracle_draw_n(min(max_n, 200L))
  storage <- oracle_pick(c("dense", "dgC", "dsC_U", "B_ddi"), c(4, 3, 1, 1))
  target <- oracle_pick(c("largest", "smallest", "largest_magnitude",
                          "smallest_magnitude", "both_ends", "nearest"),
                        c(4, 4, 1, 1, 1, 1))
  method <- oracle_pick(c("auto", "lanczos1", "lobpcg", "shift_invert"),
                        c(5, 2, 2, 1))
  if (method == "shift_invert" && stats::runif(1) < 0.85) target <- "nearest"
  list(api = "eig", structure = "gen", storage = storage, n = n,
       spectrum = oracle_pick(c("random", "clustered", "repeated",
                                "power_law", "indefinite", "tiny", "huge")),
       k = oracle_draw_k(n), target = target, method = method,
       sigma_pos = stats::runif(1), tol = oracle_draw_tol())
}

oracle_case_svd <- function(max_n) {
  storage <- oracle_pick(
    c("dense", "dgC", "dgT", "dgR", "dge", "ddi", "linop", "center",
      "scale_cols", "compose", "center_sparse"),
    c(4, 3, 1, 1, 1, 0.7, 2, 1.5, 1, 1, 1))
  shape <- if (storage == "ddi") "square" else
    oracle_pick(c("tall", "wide", "square"), c(3, 2, 1))
  a <- oracle_draw_n(max_n)
  b <- oracle_draw_n(max_n)
  m <- switch(shape, tall = max(a, b), wide = min(a, b), square = a)
  n <- switch(shape, tall = min(a, b), wide = max(a, b), square = a)
  target <- oracle_pick(c("largest", "smallest", "nearest"), c(5, 3, 2))
  method <- oracle_pick(c("auto", "golub_kahan", "randomized", "lanczos1"),
                        c(6, 3, 2, 0.5))
  list(api = "svd", structure = "svd", storage = storage, m = m, n = n,
       spectrum = oracle_pick(c("random", "clustered", "repeated",
                                "near_repeated", "power_law", "tiny", "huge",
                                "singular", "very_tiny"),
                              c(3, 1, 1, 1, 1, 1, 1, 1, 0.5)),
       k = oracle_draw_k(min(m, n)), target = target, method = method,
       sigma_pos = stats::runif(1), tol = oracle_draw_tol())
}

oracle_case_shim <- function(max_n) {
  fun <- oracle_pick(c("eigs_sym", "eigs", "svds"), c(4, 3, 3))
  n <- oracle_draw_n(min(max_n, 200L), min_n = 8L)
  case <- list(api = fun, shim = TRUE, n = n,
               storage = oracle_pick(c("dense", "dgC", "dsC_U", "function"),
                                     c(4, 3, if (fun == "eigs_sym") 1 else 0, 1)),
               tol = oracle_pick(c(1e-10, 1e-8), c(2, 1)),
               sigma_pos = stats::runif(1),
               use_ncv = stats::runif(1) < 0.25,
               use_maxitr = stats::runif(1) < 0.15,
               retvec = stats::runif(1) > 0.1)
  if (fun == "eigs_sym") {
    case$structure <- "herm"
    case$spectrum <- oracle_pick(oracle_herm_spectra)
    case$which <- oracle_pick(c("LM", "SM", "LA", "SA", "BE", "sigma"),
                              c(3, 1, 3, 2, 1, 1.5))
    case$lower <- stats::runif(1) < 0.6
  } else if (fun == "eigs") {
    case$structure <- "nonsym"
    case$spectrum <- oracle_pick(c("random_matrix", "real", "complex_pairs",
                                   "clustered"))
    case$which <- oracle_pick(c("LM", "SM", "LR", "SR", "LI", "SI", "sigma"),
                              c(4, 1.5, 2, 1, 1, 1, 1))
  } else {
    case$structure <- "svd"
    case$m <- n
    case$n <- oracle_draw_n(min(max_n, 200L), min_n = 8L)
    case$spectrum <- oracle_pick(c("random", "clustered", "power_law",
                                   "repeated"))
    case$center <- stats::runif(1) < 0.25
    case$scale <- case$center && stats::runif(1) < 0.4
    case$nu_nv <- oracle_pick(c("kk", "k0", "0k", "1k"), c(5, 1, 1, 1))
  }
  case$k <- oracle_draw_k(min(case$m %||% case$n, case$n),
                          kmax = max(1L, min(oracle_or(case$m, case$n), case$n) - 2L))
  case
}


oracle_describe <- function(case) {
  drop <- c("sigma_pos")
  keep <- case[setdiff(names(case), drop)]
  paste(vapply(names(keep), function(nm) {
    v <- keep[[nm]]
    paste0(nm, "=", if (is.numeric(v)) format(v, digits = 3) else as.character(v))
  }, ""), collapse = " ")
}

# ---------------------------------------------------------------------------
# Spectra and matrix generators
# ---------------------------------------------------------------------------

oracle_spectrum <- function(kind, n) {
  v <- switch(kind,
    random = stats::rnorm(n) * 3 + 1,
    clustered = {
      centers <- c(10, 5, -2, 1)[seq_len(min(4L, n))]
      centers[sample(length(centers), n, replace = TRUE)] *
        (1 + 1e-4 * stats::rnorm(n))
    },
    repeated = sample(c(9, 7, 5, 3, 1, -1), n, replace = TRUE,
                      prob = c(2, 3, 1, 1, 1, 1)),
    near_repeated = sample(c(9, 7, 5, 3, 1), n, replace = TRUE) +
      1e-10 * stats::rnorm(n),
    power_law = 10 * seq_len(n)^(-2),
    tiny = (stats::rnorm(n) + 2) * 1e-12,
    very_tiny = (stats::rnorm(n) + 2) * 1e-20,
    huge = (stats::rnorm(n) + 2) * 1e12,
    indefinite = c(-rev(seq_len(n %/% 2)), seq_len(n - n %/% 2)) *
      stats::runif(1, 0.5, 2),
    singular = {
      z <- stats::rnorm(n) * 2
      z[sample.int(n, max(1L, n %/% 3))] <- 0
      z
    },
    stats::rnorm(n))
  v[sample.int(n)]
}

oracle_orthogonal <- function(n, sparse = FALSE) {
  if (!sparse) {
    q <- qr(matrix(stats::rnorm(n * n), n))
    Q <- qr.Q(q)
    return(Q %*% diag(sign(diag(qr.R(q))), n))
  }
  # Block-diagonal orthogonal factor (blocks of size 1..4), rows and columns
  # permuted: an exactly orthogonal, sparse Q.
  Q <- matrix(0, n, n)
  i <- 1L
  while (i <= n) {
    b <- min(n - i + 1L, sample.int(4L, 1L))
    idx <- i:(i + b - 1L)
    Q[idx, idx] <- qr.Q(qr(matrix(stats::rnorm(b * b), b)))
    i <- i + b
  }
  p <- sample.int(n)
  Q[p, sample.int(n)]
}

oracle_sym <- function(A) (A + t(A)) / 2

oracle_herm_matrix <- function(n, spectrum, sparse = FALSE, diagonal = FALSE) {
  lambda <- oracle_spectrum(spectrum, n)
  if (diagonal) {
    return(diag(lambda, n))
  }
  Q <- oracle_orthogonal(n, sparse = sparse)
  A <- oracle_sym(Q %*% (lambda * t(Q)))
  if (sparse) A[abs(A) < 1e-300] <- 0
  A
}

# Real nonsymmetric matrix A = V J V^-1 with known J (eigenvalues, Jordan
# structure) and V, used by the Bauer-Fike / Jordan-block bound.
oracle_nonsym_matrix <- function(n, spectrum, sparse = FALSE) {
  if (spectrum == "random_matrix") {
    A <- matrix(stats::rnorm(n * n), n) / sqrt(n)
    if (sparse) {
      A[matrix(stats::runif(n * n) > max(0.08, 6 / n), n)] <- 0
      diag(A) <- stats::rnorm(n)
    }
    return(list(A = A, V = NULL, jordan = 1L))
  }
  scale <- switch(spectrum, tiny = 1e-12, huge = 1e12, 1)
  J <- matrix(0, n, n)
  lambda <- complex(0)
  jordan <- 1L
  i <- 1L
  while (i <= n) {
    left <- n - i + 1L
    kind <- if (spectrum %in% c("complex_pairs", "random", "tiny", "huge") &&
                left >= 2L && stats::runif(1) < 0.5) {
      "pair"
    } else if (spectrum == "defective" && left >= 2L && stats::runif(1) < 0.4) {
      "jordan"
    } else {
      "real"
    }
    if (kind == "pair") {
      a <- stats::rnorm(1) * 2
      b <- abs(stats::rnorm(1)) * 2 + 0.1
      J[i:(i + 1L), i:(i + 1L)] <- matrix(c(a, -b, b, a), 2)
      lambda <- c(lambda, complex(real = a, imaginary = b),
                  complex(real = a, imaginary = -b))
      i <- i + 2L
    } else if (kind == "jordan") {
      m <- min(left, sample(2:3, 1L))
      mu <- round(stats::rnorm(1) * 3, 1)
      for (j in 0:(m - 1L)) {
        J[i + j, i + j] <- mu
        if (j > 0L) J[i + j - 1L, i + j] <- 1
      }
      lambda <- c(lambda, rep(complex(real = mu), m))
      jordan <- max(jordan, m)
      i <- i + m
    } else {
      mu <- switch(spectrum,
                   clustered = sample(c(4, -3, 1), 1L) * (1 + 1e-5 * stats::rnorm(1)),
                   repeated = sample(c(4, 2, -1), 1L),
                   singular = if (stats::runif(1) < 0.35) 0 else stats::rnorm(1) * 2,
                   stats::rnorm(1) * 2)
      J[i, i] <- mu
      lambda <- c(lambda, complex(real = mu))
      i <- i + 1L
    }
  }
  J <- J * scale
  lambda <- lambda * scale
  if (sparse) {
    V <- oracle_orthogonal(n, sparse = TRUE) %*%
      diag(stats::runif(n, 1, 3), n)
    p <- seq_len(n)
  } else {
    V <- oracle_orthogonal(n) %*% diag(stats::runif(n, 1, 4), n) %*%
      oracle_orthogonal(n)
  }
  A <- V %*% J %*% solve(V)
  if (sparse) A[abs(A) < 1e-14 * max(abs(A))] <- 0
  list(A = A, V = V, J = J, lambda = lambda, jordan = jordan)
}

oracle_rect_matrix <- function(m, n, spectrum, sparse = FALSE, diagonal = FALSE,
                               root = FALSE) {
  p <- min(m, n)
  sigma <- abs(oracle_spectrum(spectrum, p))
  if (spectrum == "random") sigma <- abs(stats::rnorm(p) * 3) + 0.01
  if (root) sigma <- sqrt(sigma)
  if (diagonal) return(diag(sigma, p))
  U <- oracle_orthogonal(m, sparse = sparse)[, seq_len(p), drop = FALSE]
  V <- oracle_orthogonal(n, sparse = sparse)[, seq_len(p), drop = FALSE]
  A <- U %*% (sigma * t(V))
  if (sparse) A[abs(A) < 1e-300] <- 0
  A
}

# Exact conversion (Matrix's dense -> sparse coercions may drop tiny entries).
oracle_dgC <- function(A) {
  nz <- which(A != 0, arr.ind = TRUE)
  methods::as(Matrix::sparseMatrix(i = nz[, 1L], j = nz[, 2L], x = A[nz],
                                   dims = dim(A)), "generalMatrix")
}

oracle_store <- function(A, storage, structure = "herm") {
  switch(storage,
    dense = A,
    dgC = oracle_dgC(A),
    dsC_U = Matrix::forceSymmetric(oracle_dgC(A), uplo = "U"),
    dsC_L = Matrix::forceSymmetric(oracle_dgC(A), uplo = "L"),
    dgT = methods::as(oracle_dgC(A), "TsparseMatrix"),
    dgR = methods::as(oracle_dgC(A), "RsparseMatrix"),
    ddi = Matrix::Diagonal(x = diag(A)),
    # Built from slots: Matrix::Matrix() tests symmetry with an absolute
    # tolerance and silently symmetrises tiny-scale nonsymmetric input.
    dge = methods::new("dgeMatrix", Dim = dim(A), x = as.numeric(A)),
    dsy = methods::new("dsyMatrix", Dim = dim(A), uplo = "U", x = as.numeric(A)),
    linop = oracle_linop(A, structure),
    stop("unknown storage ", storage))
}

oracle_linop <- function(A, structure = "herm") {
  force(A)
  At <- t(A)
  eigencore::linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (At %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = if (structure == "herm") eigencore::hermitian() else
      eigencore::general(),
    name = "oracle_linop"
  )
}

# ---------------------------------------------------------------------------
# Building a case: input object + dense truth
# ---------------------------------------------------------------------------

oracle_build <- function(case) {
  oracle_with_seed(case$seed + 7L, oracle_build_impl(case))
}

oracle_build_impl <- function(case) {
  st <- case$structure
  out <- list(V = NULL, jordan = 1L, B = NULL, Bobj = NULL)
  if (isTRUE(case$shim)) {
    return(oracle_build_shim(case))
  }
  if (st == "herm") {
    n <- case$n
    if (case$storage == "composite") {
      m <- n + sample(0:10, 1L)
      sparse <- case$composite == "crossprod_sparse"
      # singular values sqrt(|lambda|): X'X has the requested spectrum (abs)
      X <- oracle_rect_matrix(m, n, case$spectrum, sparse = sparse, root = TRUE)
      Xobj <- if (sparse) oracle_dgC(X) else X
      if (case$composite == "crossprod_center") {
        Xe <- sweep(X, 2L, colMeans(X))
        obj <- eigencore::crossprod_operator(eigencore::center(Xobj))
      } else if (case$composite == "crossprod_scale") {
        w <- stats::runif(n, 0.5, 2)
        Xe <- sweep(X, 2L, w, `*`)
        obj <- eigencore::crossprod_operator(eigencore::scale_cols(Xobj, w))
      } else {
        Xe <- X
        obj <- eigencore::crossprod_operator(Xobj)
      }
      out$A <- crossprod(Xe)
      out$obj <- obj
    } else {
      sparse <- case$storage %in% c("dgC", "dsC_U", "dsC_L", "dgT", "dgR")
      A <- oracle_herm_matrix(n, case$spectrum, sparse = sparse,
                              diagonal = case$storage == "ddi")
      out$A <- A
      out$obj <- oracle_store(A, case$storage, "herm")
    }
  } else if (st == "nonsym") {
    n <- case$n
    if (case$storage == "composite") {
      g1 <- oracle_nonsym_matrix(n, case$spectrum)
      D <- diag(stats::runif(n, 0.5, 2), n)
      out$A <- g1$A %*% D
      out$obj <- eigencore::compose(g1$A, D)
    } else {
      sparse <- case$storage %in% c("dgC", "dgT", "dgR")
      g <- oracle_nonsym_matrix(n, case$spectrum, sparse = sparse)
      out$A <- g$A
      out$V <- g$V
      out$jordan <- g$jordan
      out$J <- g$J
      out$obj <- oracle_store(g$A, case$storage, "nonsym")
    }
  } else if (st == "cherm") {
    n <- case$n
    lambda <- oracle_spectrum(case$spectrum, n)
    Z <- matrix(complex(real = stats::rnorm(n * n), imaginary = stats::rnorm(n * n)), n)
    Q <- qr.Q(qr(Z))
    A <- Q %*% (lambda * Conj(t(Q)))
    A <- (A + Conj(t(A))) / 2
    out$A <- A
    out$obj <- A
  } else if (st == "cgen") {
    n <- case$n
    if (case$spectrum == "random_matrix") {
      A <- matrix(complex(real = stats::rnorm(n * n),
                          imaginary = stats::rnorm(n * n)), n) / sqrt(n)
    } else {
      lambda <- switch(case$spectrum,
        repeated = sample(c(2 + 1i, -1, 3i), n, replace = TRUE),
        tiny = complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)) * 1e-12,
        huge = complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)) * 1e12,
        complex(real = stats::rnorm(n), imaginary = stats::rnorm(n)) * 2)
      Z <- matrix(complex(real = stats::rnorm(n * n), imaginary = stats::rnorm(n * n)), n)
      V <- qr.Q(qr(Z)) %*% diag(stats::runif(n, 1, 3), n)
      A <- V %*% diag(lambda, n) %*% solve(V)
      out$V <- V
    }
    out$A <- A
    out$obj <- A
  } else if (st == "gen") {
    n <- case$n
    lambda <- oracle_spectrum(case$spectrum, n)
    if (case$storage == "B_ddi") {
      d <- stats::runif(n, 0.5, 3)
      Q <- oracle_orthogonal(n)
      # B = diag(d); A = B^{1/2} Q diag(lambda) Q' B^{1/2}
      Bh <- sqrt(d)
      A <- oracle_sym((Bh * Q) %*% (lambda * t(Bh * Q)))
      B <- diag(d, n)
      out$Bobj <- Matrix::Diagonal(x = d)
      out$obj <- A
    } else {
      sparse <- case$storage %in% c("dgC", "dsC_U")
      Q <- oracle_orthogonal(n, sparse = sparse)
      s <- stats::runif(n, 0.5, 2)
      X <- Q %*% diag(s, n)           # B-orthonormal eigenvectors
      Xi <- solve(X)
      A <- oracle_sym(t(Xi) %*% (lambda * Xi))
      B <- oracle_sym(t(Xi) %*% Xi)
      if (sparse) {
        A[abs(A) < 1e-300] <- 0
        B[abs(B) < 1e-300] <- 0
      }
      out$obj <- oracle_store(A, if (case$storage == "dense") "dense" else case$storage)
      out$Bobj <- oracle_store(B, if (case$storage == "dense") "dense" else case$storage)
    }
    out$A <- A
    out$B <- B
  } else if (st == "svd") {
    m <- case$m
    n <- case$n
    stor <- case$storage
    sparse <- stor %in% c("dgC", "dgT", "dgR", "center_sparse")
    X <- oracle_rect_matrix(m, n, case$spectrum, sparse = sparse,
                            diagonal = stor == "ddi")
    if (stor == "center" || stor == "center_sparse") {
      out$A <- sweep(X, 2L, colMeans(X))
      out$obj <- eigencore::center(if (sparse) oracle_dgC(X) else X)
    } else if (stor == "scale_cols") {
      w <- stats::runif(n, 0.5, 2)
      out$A <- sweep(X, 2L, w, `*`)
      out$obj <- eigencore::scale_cols(X, w)
    } else if (stor == "compose") {
      p <- n
      Y <- oracle_orthogonal(n) %*% diag(stats::runif(n, 0.5, 2), n)
      out$A <- X %*% Y
      out$obj <- eigencore::compose(X, Y)
    } else {
      out$A <- X
      out$obj <- oracle_store(X, stor, "svd")
    }
  }
  out
}

oracle_build_shim <- function(case) {
  out <- list(V = NULL, jordan = 1L, B = NULL)
  if (case$api == "eigs_sym") {
    sparse <- case$storage %in% c("dgC", "dsC_U")
    A <- oracle_herm_matrix(case$n, case$spectrum, sparse = sparse)
    out$A <- A
    stored <- A
    if (case$storage %in% c("dense", "dgC")) {
      # RSpectra reads only one triangle: poison the other one.
      junk <- matrix(stats::rnorm(case$n^2), case$n) * (sparse == FALSE)
      if (isTRUE(case$lower)) {
        stored[upper.tri(stored)] <- stored[upper.tri(stored)] + junk[upper.tri(junk)]
      } else {
        stored[lower.tri(stored)] <- stored[lower.tri(stored)] + junk[lower.tri(junk)]
      }
    }
    out$obj <- switch(case$storage,
                      dense = stored,
                      dgC = oracle_dgC(stored),
                      dsC_U = Matrix::forceSymmetric(oracle_dgC(A), uplo = "U"),
                      "function" = {
                        force(A)
                        function(x, args) as.numeric(A %*% x)
                      })
  } else if (case$api == "eigs") {
    sparse <- case$storage == "dgC"
    g <- oracle_nonsym_matrix(case$n, case$spectrum, sparse = sparse)
    out$A <- g$A
    out$V <- g$V
    out$jordan <- g$jordan
    out$obj <- switch(case$storage,
                      dense = g$A,
                      dgC = oracle_dgC(g$A),
                      "function" = {
                        A <- g$A
                        function(x, args) as.numeric(A %*% x)
                      })
  } else {
    sparse <- case$storage == "dgC"
    X <- oracle_rect_matrix(case$m, case$n, case$spectrum, sparse = sparse)
    Xe <- X
    if (isTRUE(case$center)) Xe <- sweep(Xe, 2L, colMeans(Xe))
    if (isTRUE(case$scale)) {
      s <- sqrt(colSums(Xe^2) / max(1, nrow(Xe) - 1))
      s[s == 0] <- 1
      Xe <- sweep(Xe, 2L, s, `/`)
    }
    out$A <- Xe
    out$raw <- X
    out$obj <- switch(case$storage,
                      dense = X,
                      dgC = oracle_dgC(X),
                      "function" = {
                        force(X)
                        function(x, args) as.numeric(X %*% x)
                      })
  }
  out
}

# ---------------------------------------------------------------------------
# Targets, methods, oracle
# ---------------------------------------------------------------------------

oracle_sigma <- function(case, truth_values) {
  r <- range(Re(truth_values))
  if (case$structure %in% c("nonsym", "cgen")) {
    # a point inside the spectrum's bounding box, off the real axis only for
    # complex input (real shift-invert needs a real sigma)
    s <- r[1] + case$sigma_pos * diff(r)
    return(s + 0.0123 * max(1e-300, diff(r)))
  }
  w <- diff(r)
  r[1] - 0.1 * w + case$sigma_pos * 1.2 * w + 1e-3 * max(w, abs(r))
}

oracle_target <- function(case, sigma) {
  switch(case$target,
    largest = eigencore::largest(),
    smallest = eigencore::smallest(),
    largest_magnitude = eigencore::largest_magnitude(),
    smallest_magnitude = eigencore::smallest_magnitude(),
    largest_real = eigencore::largest_real(),
    smallest_real = eigencore::smallest_real(),
    largest_imaginary = eigencore::largest_imaginary(),
    smallest_imaginary = eigencore::smallest_imaginary(),
    nearest = eigencore::nearest(sigma),
    both_ends = {
      kl <- case$k %/% 2L
      eigencore::both_ends(kl, case$k - kl)
    })
}

oracle_method <- function(case, sigma) {
  switch(case$method,
    auto = eigencore::auto(),
    lanczos1 = eigencore::lanczos(block = 1L),
    lanczos2 = eigencore::lanczos(block = 2L),
    lanczos3 = eigencore::lanczos(block = 3L),
    lobpcg = eigencore::lobpcg(),
    shift_invert = eigencore::shift_invert(sigma),
    golub_kahan = eigencore::golub_kahan(),
    randomized = eigencore::randomized())
}

# Truth: values of the dense problem plus the norms the harness divides by.
oracle_truth <- function(case, prob) {
  A <- prob$A
  st <- case$structure
  if (st %in% c("herm", "cherm")) {
    ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    list(values = ev, normA = max(abs(ev)), normB = 1)
  } else if (st == "gen") {
    L <- chol(prob$B)
    C <- backsolve(L, t(backsolve(L, A, transpose = TRUE)), transpose = TRUE)
    ev <- eigen(oracle_sym(C), symmetric = TRUE, only.values = TRUE)$values
    list(values = ev, normA = svd(A, 0, 0)$d[1],
         normB = max(eigen(prob$B, symmetric = TRUE, only.values = TRUE)$values),
         L = L)
  } else if (st %in% c("nonsym", "cgen")) {
    e <- eigen(A)
    list(values = e$values, vectors = e$vectors, normA = svd(A, 0, 0)$d[1],
         normB = 1)
  } else {
    d <- svd(A, 0, 0)$d
    list(values = d, normA = d[1], normB = 1)
  }
}

oracle_score <- function(x, case, sigma) {
  switch(case$target,
    largest = Re(x), smallest = -Re(x),
    largest_real = Re(x), smallest_real = -Re(x),
    largest_magnitude = Mod(x), smallest_magnitude = -Mod(x),
    largest_imaginary = Im(x), smallest_imaginary = -Im(x),
    nearest = -Mod(x - sigma),
    Re(x))
}

# ---------------------------------------------------------------------------
# Running one case
# ---------------------------------------------------------------------------

oracle_expected_error_patterns <- c(
  "not supported", "unsupported", "does not support", "only supports",
  "requires", "must be", "cannot", "can only", "only available",
  "not available", "not implemented", "dispatch_unavailable",
  "Invalid eigencore plan", "is only defined", "only for",
  "not defined", "must not exceed", "at most", "less than",
  "Refusing to densify", "supports", "needs an explicit matrix-backed",
  "computes the eigenvalues nearest sigma", "require a matrix input",
  "produced zero or non-finite column scales",
  "more than the operator dimension", "is an eigensolver method"
)

oracle_unexpected_error_patterns <- c(
  "subscript out of bounds", "non-conformable", "object '.*' not found",
  "missing value where TRUE/FALSE", "argument of length 0", "argument is of length zero",
  "\\$ operator is invalid", "attempt to apply non-function",
  "invalid 'type'", "NA/NaN/Inf in foreign function", "replacement has",
  "incorrect number of dimensions", "could not find function",
  "arguments imply differing", "unused argument", "not a matrix",
  "infinite or missing values", "LAPACK routine", "non-numeric argument",
  "the condition has length", "status=-", "invalid argument", "Lapack routine",
  "max_subspace must be at least"
)

oracle_classify_error <- function(e) {
  msg <- conditionMessage(e)
  if (any(vapply(oracle_unexpected_error_patterns, grepl, NA, x = msg)))
    return("unexpected")
  if (any(vapply(oracle_expected_error_patterns, grepl, NA, x = msg,
                 ignore.case = TRUE)))
    return("expected")
  "unexpected"
}

oracle_invoke <- function(case, prob, sigma) {
  if (isTRUE(case$shim)) return(oracle_invoke_shim(case, prob, sigma))
  target <- oracle_target(case, sigma)
  method <- oracle_method(case, sigma)
  if (case$api == "eig") {
    eigencore::eig_partial(prob$obj, k = case$k, target = target, B = prob$Bobj,
                           method = method, tol = case$tol, seed = case$seed)
  } else {
    eigencore::svd_partial(prob$obj, rank = case$k, target = target,
                           method = method, tol = case$tol, seed = case$seed)
  }
}

oracle_shim_opts <- function(case) {
  opts <- list(tol = case$tol)
  if (isTRUE(case$use_ncv)) opts$ncv <- min(case$n, max(2L * case$k + 1L, 20L))
  if (isTRUE(case$use_maxitr)) opts$maxitr <- 2000L
  if (!isTRUE(case$retvec)) opts$retvec <- FALSE
  opts
}

oracle_invoke_shim <- function(case, prob, sigma) {
  opts <- oracle_shim_opts(case)
  fn_args <- if (case$storage == "function") list(n = case$n) else list()
  # The shims have no seed= argument: seed the global stream like a user would.
  set.seed(case$seed)
  if (case$api == "eigs_sym") {
    which <- if (case$which == "sigma") "LM" else case$which
    do.call(eigencore::eigs_sym, c(list(prob$obj, k = case$k, which = which,
      sigma = if (case$which == "sigma") sigma, opts = opts,
      lower = isTRUE(case$lower)), fn_args))
  } else if (case$api == "eigs") {
    which <- if (case$which == "sigma") "LM" else case$which
    do.call(eigencore::eigs, c(list(prob$obj, k = case$k, which = which,
      sigma = if (case$which == "sigma") sigma, opts = opts), fn_args))
  } else {
    nu <- switch(case$nu_nv, kk = case$k, k0 = case$k, "0k" = 0L, "1k" = 1L)
    nv <- switch(case$nu_nv, kk = case$k, k0 = 0L, "0k" = case$k, "1k" = case$k)
    o <- list(tol = case$tol)
    if (isTRUE(case$center)) o$center <- TRUE
    if (isTRUE(case$scale)) o$scale <- TRUE
    if (case$storage == "function") {
      X <- prob$raw
      eigencore::svds(prob$obj, k = case$k, nu = nu, nv = nv, opts = o,
                      Atrans = function(x, args) as.numeric(crossprod(X, x)),
                      dim = dim(X))
    } else {
      eigencore::svds(prob$obj, k = case$k, nu = nu, nv = nv, opts = o)
    }
  }
}

# Map shim targets onto the harness's target vocabulary.
oracle_shim_target <- function(case) {
  if (case$api == "svds") return("largest")
  switch(case$which,
         LM = "largest_magnitude", SM = "smallest_magnitude",
         LA = "largest", SA = "smallest", BE = "both_ends_shim",
         LR = "largest_real", SR = "smallest_real",
         LI = "largest_imaginary", SI = "smallest_imaginary",
         sigma = "nearest")
}

oracle_run_case <- function(id_or_case, check_rspectra = TRUE) {
  case <- if (is.list(id_or_case)) id_or_case else oracle_case(id_or_case)
  started <- proc.time()[["elapsed"]]
  rec <- list(id = case$id, family = case$family, api = case$api,
              structure = case$structure, storage = case$storage,
              spectrum = case$spectrum, target = case$target %||% case$which,
              method = case$method %||% "shim", k = case$k, n = case$n,
              m = case$m %||% case$n, tol = case$tol,
              status = NA_character_, error = NA_character_,
              certified = NA, completeness = NA_character_, route = NA_character_,
              true_backward = NA_real_, reported_backward = NA_real_,
              ortho_loss = NA_real_, value_error = NA_real_, set_ok = NA,
              hard = character(), soft = character(), seconds = NA_real_,
              describe = oracle_describe(case))
  prob <- tryCatch(oracle_build(case), error = function(e) e)
  if (inherits(prob, "error")) {
    rec$status <- "harness_error"
    rec$error <- conditionMessage(prob)
    rec$hard <- paste0("harness: build failed: ", conditionMessage(prob))
    return(rec)
  }
  truth <- oracle_truth(case, prob)
  sigma <- oracle_sigma(case, truth$values)
  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch(oracle_invoke(case, prob, sigma), error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  rec$seconds <- proc.time()[["elapsed"]] - started
  if (inherits(fit, "error")) {
    rec$error <- conditionMessage(fit)
    cls <- oracle_classify_error(fit)
    rec$status <- if (cls == "expected") "error_expected" else "error_unexpected"
    if (cls != "expected") {
      rec$hard <- paste0("unexpected error: ", conditionMessage(fit))
    }
    return(rec)
  }
  rec$status <- "ok"
  chk <- tryCatch(oracle_check(case, prob, truth, sigma, fit),
                  error = function(e) list(hard = paste0("harness: check failed: ",
                                                          conditionMessage(e))))
  for (nm in names(chk)) rec[[nm]] <- chk[[nm]]
  if (isTRUE(case$shim) && check_rspectra) {
    rs <- tryCatch(oracle_check_rspectra(case, prob, truth, sigma, fit, chk),
                   error = function(e) list(soft = paste0("rspectra check error: ",
                                                          conditionMessage(e))))
    rec$hard <- c(rec$hard, rs$hard)
    rec$soft <- c(rec$soft, rs$soft)
  }
  rec
}

# Normalise eigencore / shim results.
oracle_normalise <- function(case, fit) {
  if (case$api %in% c("svd", "svds")) {
    d <- oracle_or(fit$d, fit$values)
    list(kind = "svd", values = d, u = fit$u, v = fit$v,
         cert = fit$certificate, nconv = fit$nconv,
         completeness = fit$certificate$target_completeness,
         route = oracle_or(fit$actual_method, oracle_or(fit$diagnostics$actual_method, NA_character_)))
  } else {
    list(kind = "eig", values = fit$values,
         vectors = if (!is.null(fit$right_vectors)) fit$right_vectors else fit$vectors,
         cert = fit$certificate, nconv = fit$nconv,
         completeness = fit$certificate$target_completeness,
         route = oracle_or(fit$actual_method, oracle_or(fit$diagnostics$actual_method, NA_character_)))
  }
}

oracle_target_kind <- function(case) {
  if (isTRUE(case$shim)) oracle_shim_target(case) else case$target
}

oracle_check <- function(case, prob, truth, sigma, fit) {
  hard <- character()
  soft <- character()
  res <- oracle_normalise(case, fit)
  cert <- res$cert
  A <- prob$A
  eps <- .Machine$double.eps
  tol <- case$tol
  k <- case$k
  nA <- truth$normA
  dims <- dim(A)
  vals <- res$values
  nval <- length(vals)
  passed <- isTRUE(cert$passed)
  tkind <- oracle_target_kind(case)
  out <- list(certified = passed,
              completeness = oracle_or(res$completeness, NA_character_),
              route = if (is.null(res$route)) NA_character_ else as.character(res$route)[1])

  # ---- API invariants -----------------------------------------------------
  if (nval > k) hard <- c(hard, sprintf("api: %d values returned for k = %d", nval, k))
  if (!is.null(res$nconv) && length(res$nconv) == 1L && !is.na(res$nconv) &&
      res$nconv > k) {
    hard <- c(hard, sprintf("api: nconv %d > k %d", res$nconv, k))
  }
  if (passed && nval != k) {
    hard <- c(hard, sprintf("api: certified with %d of %d values", nval, k))
  }
  if (any(!is.finite(vals))) {
    if (passed) hard <- c(hard, "soundness: certified non-finite values")
    else soft <- c(soft, "non-finite values")
  }
  want_vectors <- !(isTRUE(case$shim) && !isTRUE(case$retvec))
  if (res$kind == "eig") {
    X <- res$vectors
    if (want_vectors && nval > 0L) {
      if (is.null(X)) {
        hard <- c(hard, "api: vectors missing")
      } else if (nrow(X) != dims[1] || ncol(X) != nval) {
        hard <- c(hard, sprintf("api: vectors %dx%d for n=%d, %d values",
                                nrow(X), ncol(X), dims[1], nval))
        X <- NULL
      }
    }
  } else {
    if (want_vectors && nval > 0L) {
      nu <- if (isTRUE(case$shim)) switch(case$nu_nv, kk = k, k0 = k, "0k" = 0L, "1k" = 1L) else nval
      nv <- if (isTRUE(case$shim)) switch(case$nu_nv, kk = k, k0 = 0L, "0k" = k, "1k" = k) else nval
      if (nu > 0L && (is.null(res$u) || nrow(res$u) != dims[1] ||
                      ncol(res$u) != min(nu, nval))) {
        hard <- c(hard, "api: u has wrong dims")
      }
      if (nv > 0L && (is.null(res$v) || nrow(res$v) != dims[2] ||
                      ncol(res$v) != min(nv, nval))) {
        hard <- c(hard, "api: v has wrong dims")
      }
      if (nu == 0L && !is.null(res$u) && ncol(res$u) > 0L) hard <- c(hard, "api: u returned for nu = 0")
      if (nv == 0L && !is.null(res$v) && ncol(res$v) > 0L) hard <- c(hard, "api: v returned for nv = 0")
    }
  }
  # Ordering per target.
  if (nval > 1L && all(is.finite(vals))) {
    ord_tol <- 1e-8 * max(nA, 1e-300) + 1e3 * eps * max(abs(vals))
    if (isTRUE(case$shim) && case$api == "eigs_sym") {
      if (any(diff(vals) > ord_tol)) hard <- c(hard, "api: eigs_sym values not decreasing")
    } else if (isTRUE(case$shim) && case$api == "svds" || res$kind == "svd" && tkind == "largest") {
      if (any(diff(vals) > ord_tol)) hard <- c(hard, "api: singular values not decreasing")
    } else if (tkind == "both_ends") {
      kl <- k %/% 2L
      lo <- vals[seq_len(min(kl, nval))]
      hi <- vals[-seq_len(min(kl, nval))]
      if (any(diff(Re(lo)) < -ord_tol) || any(diff(Re(hi)) > ord_tol)) {
        hard <- c(hard, "api: both_ends values not (low ascending, high descending)")
      }
    } else if (!isTRUE(case$shim) && !is.null(tkind) && tkind != "both_ends_shim") {
      s <- oracle_score(vals, list(target = tkind), sigma)
      if (any(diff(s) > ord_tol)) {
        hard <- c(hard, sprintf("api: values not ordered for target %s", tkind))
      }
    }
  }

  if (nval == 0L) {
    out$hard <- hard
    out$soft <- c(soft, "no values returned")
    return(out)
  }

  # ---- Independent residuals ---------------------------------------------
  if (res$kind == "eig") {
    X <- res$vectors
    if (is.null(X) || !want_vectors) {
      out$hard <- hard
      out$soft <- c(soft, "no vectors: soundness not checkable")
      return(oracle_check_values_only(out, case, truth, sigma, vals, passed, hard, soft))
    }
    X <- as.matrix(X)
    AX <- A %*% X
    if (case$structure == "gen") {
      BX <- prob$B %*% X
      R <- AX - sweep(BX, 2L, vals, `*`)
      xn <- sqrt(colSums(Mod(X)^2))
      eta <- oracle_div(sqrt(colSums(Mod(R)^2)), (nA + Mod(vals) * truth$normB) * xn)
      # residual bound in the B-inner product: |theta - lambda| <= ||r||_{B^-1} / ||x||_B
      Linv_r <- backsolve(truth$L, R, transpose = TRUE)
      xB <- sqrt(colSums((truth$L %*% X)^2))
      vbound <- sqrt(colSums(Linv_r^2)) / xB
      G <- crossprod(X, BX)
      ortho <- max(abs(G - diag(ncol(X))))
    } else {
      R <- AX - sweep(X, 2L, vals, `*`)
      xn <- sqrt(colSums(Mod(X)^2))
      eta <- oracle_div(sqrt(colSums(Mod(R)^2)), (nA + Mod(vals)) * xn)
      vbound <- sqrt(colSums(Mod(R)^2)) / xn
      G <- Conj(t(X)) %*% X
      ortho <- max(Mod(G - diag(ncol(X))))
    }
  } else {
    U <- res$u
    V <- res$v
    if (is.null(U) || is.null(V) || ncol(U) != nval || ncol(V) != nval) {
      return(oracle_check_values_only(out, case, truth, sigma, vals, passed, hard,
                                      c(soft, "partial singular vectors: residual not checkable")))
    }
    U <- as.matrix(U)
    V <- as.matrix(V)
    r1 <- A %*% V - sweep(U, 2L, vals, `*`)
    r2 <- crossprod(A, U) - sweep(V, 2L, vals, `*`)
    comb <- sqrt(colSums(r1^2) + colSums(r2^2))
    eta <- oracle_div(comb, rep(nA, length(comb)))
    # Jordan-Wielandt: |s - sigma_j| <= ||[r1; r2]|| / ||[u; v]||
    vbound <- comb / sqrt(colSums(U^2) + colSums(V^2))
    ortho <- max(abs(crossprod(U) - diag(nval)), abs(crossprod(V) - diag(nval)))
  }
  out$true_backward <- max(eta)
  out$reported_backward <- oracle_or(cert$max_backward_error, NA_real_)
  out$ortho_loss <- ortho
  n_eff <- max(dims)
  slack_abs <- 64 * n_eff * eps
  # Reported backward error must over-estimate the true one (C12).
  rep_be <- cert$backward_error
  if (!is.null(rep_be) && length(rep_be) == length(eta) && all(is.finite(rep_be))) {
    under <- eta > rep_be * 1.01 + slack_abs & eta > 1e3 * n_eff * eps
    if (any(under)) {
      hard <- c(hard, sprintf(
        "soundness: reported backward error %.3g under true %.3g (index %d)",
        rep_be[which(under)[1]], eta[which(under)[1]], which(under)[1]))
    }
  }
  if (passed) {
    if (any(!is.finite(eta)) || max(eta) > tol * 1.05 + slack_abs) {
      hard <- c(hard, sprintf("soundness: certified but true backward error %.3g > tol %.1g",
                              max(eta), tol))
    }
    orth_required <- !identical(cert$orthogonality_required, FALSE) &&
      case$structure %in% c("herm", "cherm", "gen", "svd")
    orth_tol <- max(tol, sqrt(eps))
    if (orth_required && ortho > orth_tol * 2 + 64 * n_eff * eps) {
      hard <- c(hard, sprintf("soundness: certified but orthogonality loss %.3g > %.3g",
                              ortho, orth_tol))
    }
  }

  # ---- Values against the oracle -----------------------------------------
  vchk <- oracle_value_checks(case, prob, truth, sigma, vals, vbound, tkind)
  out$value_error <- vchk$value_error
  out$set_ok <- vchk$set_ok
  if (passed && length(vchk$far)) {
    hard <- c(hard, paste0("soundness: certified value far from spectrum: ", vchk$far))
  }
  claims <- c("exact", "probed", "repaired", "inertia_verified")
  if (passed && isFALSE(vchk$set_ok)) {
    msg <- sprintf("certified set differs from oracle target set (completeness=%s): %s",
                   oracle_or(res$completeness, "none"), vchk$set_msg)
    if (!is.null(res$completeness) && res$completeness %in% claims) {
      hard <- c(hard, paste0("soundness: ", msg))
    } else {
      soft <- c(soft, msg)
    }
  }
  if (!passed) {
    soft <- c(soft, "uncertified")
  }
  out$hard <- hard
  out$soft <- soft
  out
}

oracle_check_values_only <- function(out, case, truth, sigma, vals, passed, hard, soft) {
  vchk <- oracle_value_checks(case, NULL, truth, sigma, vals, NULL,
                              oracle_target_kind(case))
  out$value_error <- vchk$value_error
  out$set_ok <- vchk$set_ok
  if (!passed) soft <- c(soft, "uncertified")
  out$hard <- hard
  out$soft <- soft
  out
}

# Oracle target set as a vector of values (best first by score).
oracle_target_set <- function(case, truth, sigma, tkind, k) {
  ev <- truth$values
  if (case$api %in% c("svd", "svds")) {
    s <- switch(tkind, largest = ev, smallest = -ev, nearest = -abs(ev - sigma), ev)
    return(ev[order(-s)][seq_len(k)])
  }
  if (tkind == "both_ends" || tkind == "both_ends_shim") {
    kh <- if (tkind == "both_ends") k - k %/% 2L else (k + 1L) %/% 2L
    kl <- k - kh
    srt <- sort(Re(ev))
    return(c(head(srt, kl), rev(tail(srt, kh))))
  }
  s <- oracle_score(ev, list(target = tkind), sigma)
  ev[order(-s)][seq_len(k)]
}

oracle_value_checks <- function(case, prob, truth, sigma, vals, vbound, tkind) {
  eps <- .Machine$double.eps
  ev <- truth$values
  nA <- truth$normA
  n_eff <- length(ev)
  st <- case$structure
  k <- length(vals)
  noise <- 256 * n_eff * eps * max(nA, 1e-300)
  far <- character()
  # per-value bound
  if (!is.null(vbound)) {
    cand <- ev
    if (st == "svd" && !is.null(prob) && nrow(prob$A) != ncol(prob$A)) cand <- c(ev, 0)
    if (st %in% c("nonsym", "cgen")) {
      kap <- oracle_bf_constant(prob, truth)
      for (i in seq_len(k)) {
        d <- min(Mod(vals[i] - cand))
        eps_i <- vbound[i] + noise
        b <- oracle_jordan_bound(kap * eps_i, prob$jordan %||% 1L)
        if (is.finite(b) && d > 1.5 * b + noise) {
          far <- c(far, sprintf("value %s: distance %.3g > Bauer-Fike bound %.3g",
                                format(vals[i], digits = 6), d, b))
        }
      }
    } else {
      for (i in seq_len(k)) {
        d <- min(abs(Re(vals[i]) - cand))
        if (abs(Im(vals[i])) > 0) d <- max(d, abs(Im(vals[i])))
        if (d > 1.05 * vbound[i] + noise) {
          far <- c(far, sprintf("value %.10g: distance %.3g > residual bound %.3g",
                                Re(vals[i]), d, vbound[i]))
        }
      }
    }
  }
  # set comparison
  K <- case$k
  set_ok <- NA
  set_msg <- ""
  value_error <- NA_real_
  if (k == K && K <= length(ev) && all(is.finite(vals))) {
    target_set <- oracle_target_set(case, truth, sigma, tkind, K)
    if (st %in% c("nonsym", "cgen")) {
      kap <- if (!is.null(prob)) oracle_bf_constant(prob, truth) else 1
      dmax <- if (is.null(vbound)) 1e-6 * nA else max(vbound)
      delta <- 2 * oracle_jordan_bound(kap * (dmax + noise), prob$jordan %||% 1L) +
        1e3 * noise
    } else {
      dmax <- if (is.null(vbound)) case$tol * 2 * max(nA, 1e-300) else max(vbound)
      delta <- 2 * sqrt(K) * dmax + noise
    }
    if (tkind %in% c("both_ends", "both_ends_shim")) {
      kh <- if (tkind == "both_ends") K - K %/% 2L else (K + 1L) %/% 2L
      kl <- K - kh
      srt <- sort(Re(vals))
      got <- c(head(srt, kl), rev(tail(srt, kh)))
      err <- abs(got - Re(target_set))
    } else if (case$api %in% c("svd", "svds")) {
      sk <- switch(tkind, largest = vals, smallest = -vals,
                   nearest = -abs(vals - sigma), vals)
      st_k <- switch(tkind, largest = target_set, smallest = -target_set,
                     nearest = -abs(target_set - sigma), target_set)
      err <- abs(sort(sk) - sort(st_k))
    } else {
      sk <- oracle_score(vals, list(target = tkind), sigma)
      st_k <- oracle_score(target_set, list(target = tkind), sigma)
      err <- abs(sort(sk) - sort(st_k))
    }
    value_error <- max(err) / max(nA, 1e-300)
    set_ok <- is.finite(delta) && max(err) <= delta
    if (!set_ok) {
      set_msg <- sprintf("max score mismatch %.3g > %.3g; got [%s] want [%s]",
                         max(err), delta,
                         paste(format(head(vals, 8), digits = 6), collapse = ", "),
                         paste(format(head(target_set, 8), digits = 6), collapse = ", "))
    }
  }
  list(far = far, set_ok = set_ok, set_msg = set_msg, value_error = value_error)
}

# kappa(V) for Bauer-Fike, from the generator's exact V when known.
oracle_bf_constant <- function(prob, truth) {
  V <- prob$V
  if (is.null(V)) V <- truth$vectors
  if (is.null(V)) return(Inf)
  s <- svd(V, 0, 0)$d
  if (min(s) <= 0) return(Inf)
  s[1] / s[length(s)]
}

# Largest t with t^m <= delta (1 + t)^(m - 1) (Chatelin's bound for Jordan
# blocks of size m); m = 1 gives Bauer-Fike, t = delta.
oracle_jordan_bound <- function(delta, m) {
  if (!is.finite(delta)) return(Inf)
  if (m <= 1L) return(delta)
  f <- function(t) t^m - delta * (1 + t)^(m - 1)
  hi <- max(1, 2 * delta * 2^m)
  tryCatch(stats::uniroot(f, c(0, hi), tol = 1e-14 * hi)$root, error = function(e) Inf)
}

# ---------------------------------------------------------------------------
# RSpectra comparison for shim cases
# ---------------------------------------------------------------------------

oracle_check_rspectra <- function(case, prob, truth, sigma, fit, chk) {
  if (!requireNamespace("RSpectra", quietly = TRUE)) return(list())
  hard <- character()
  soft <- character()
  opts <- list(tol = case$tol)
  set.seed(case$seed)
  rs <- tryCatch(suppressWarnings(switch(case$api,
    eigs_sym = {
      obj <- if (case$storage == "function") prob$A else prob$obj
      RSpectra::eigs_sym(obj, k = case$k,
                         which = if (case$which == "sigma") "LM" else case$which,
                         sigma = if (case$which == "sigma") sigma,
                         opts = list(tol = case$tol), lower = isTRUE(case$lower))
    },
    eigs = {
      obj <- if (case$storage == "function") prob$A else prob$obj
      RSpectra::eigs(obj, k = case$k,
                     which = if (case$which == "sigma") "LM" else case$which,
                     sigma = if (case$which == "sigma") sigma,
                     opts = list(tol = case$tol))
    },
    svds = {
      obj <- if (case$storage == "function") prob$raw else prob$obj
      o <- list(tol = case$tol)
      if (isTRUE(case$center)) o$center <- TRUE
      if (isTRUE(case$scale)) o$scale <- TRUE
      RSpectra::svds(obj, k = case$k, nu = 0, nv = 0, opts = o)
    })), error = function(e) e)
  if (inherits(rs, "error")) return(list(soft = "rspectra errored"))
  rv <- if (case$api == "svds") rs$d else rs$values
  ev <- if (case$api == "svds") fit$d else fit$values
  if (length(rv) != case$k || length(ev) != case$k) {
    return(list(soft = "rspectra/eigencore returned fewer values"))
  }
  # Is RSpectra itself right (matches the oracle set)?
  tkind <- oracle_target_kind(case)
  ts <- oracle_target_set(case, truth, sigma, tkind, case$k)
  scale <- max(truth$normA, 1e-300)
  # Compare in target-score space: ties in the score (+-lambda for LM,
  # conjugate pairs for LR) may legitimately pick different values.
  score <- function(x) {
    if (case$api == "svds") return(sort(x))
    if (tkind == "both_ends_shim") return(sort(Re(x)))
    sort(oracle_score(x, list(target = tkind), sigma))
  }
  rs_ok <- max(abs(score(rv) - score(ts))) <= 1e-6 * scale
  ec_ok <- isTRUE(chk$certified) && isTRUE(chk$set_ok)
  if (rs_ok && ec_ok) {
    if (max(abs(score(ev) - score(rv))) > 1e-6 * scale) {
      hard <- c(hard, sprintf("shim: values differ from RSpectra: eigencore [%s] RSpectra [%s]",
                              paste(format(ev, digits = 5), collapse = ", "),
                              paste(format(rv, digits = 5), collapse = ", ")))
    } else if (max(Mod(ev - rv)) > 1e-6 * scale) {
      # Same set, different order (or a different member of a score tie).
      # RSpectra's order for SR/SM/LI/SI is not monotone in the target score,
      # so this is reported, not failed.
      soft <- c(soft, "shim: same set as RSpectra, different order or tie member")
    }
  } else if (!rs_ok && ec_ok) {
    soft <- c(soft, "rspectra wrong, eigencore right")
  } else if (rs_ok && !isTRUE(chk$certified)) {
    soft <- c(soft, "rspectra right, eigencore uncertified")
  }
  list(hard = hard, soft = soft)
}

# ---------------------------------------------------------------------------
# Batch execution with crash isolation
# ---------------------------------------------------------------------------

oracle_helper_path <- function() {
  p <- tryCatch(testthat::test_path("helper-oracle.R"), error = function(e) "")
  if (nzchar(p) && file.exists(p)) return(normalizePath(p))
  for (cand in c("tests/testthat/helper-oracle.R", "helper-oracle.R")) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  stop("cannot locate helper-oracle.R", call. = FALSE)
}

oracle_subprocess_available <- function() {
  if (!requireNamespace("callr", quietly = TRUE)) return(FALSE)
  installed <- tryCatch(find.package("eigencore", quiet = TRUE), error = function(e) "")
  if (!length(installed) || !nzchar(installed[1])) return(FALSE)
  # under pkgload::load_all the installed copy may be stale
  loaded <- tryCatch(getNamespaceInfo(asNamespace("eigencore"), "path"),
                     error = function(e) installed[1])
  identical(normalizePath(loaded), normalizePath(installed[1]))
}

oracle_record_frame <- function(recs) {
  scalar <- c("id", "family", "api", "structure", "storage", "spectrum",
              "target", "method", "k", "n", "m", "tol", "status", "error",
              "certified", "completeness", "route", "true_backward",
              "reported_backward", "ortho_loss", "value_error", "set_ok",
              "seconds", "describe")
  rows <- lapply(recs, function(r) {
    vals <- lapply(scalar, function(nm) {
      v <- r[[nm]]
      if (is.null(v) || length(v) == 0L) NA else v[[1]]
    })
    names(vals) <- scalar
    vals$hard <- paste(r$hard, collapse = " | ")
    vals$soft <- paste(r$soft, collapse = " | ")
    as.data.frame(vals, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

oracle_env <- function() {
  c(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
    EIGENCORE_ORACLE_CHILD = "1")
}


oracle_child <- function(ids, helper) {
  suppressPackageStartupMessages(library(eigencore))
  options(eigencore.threads = 1L)
  env <- new.env()
  sys.source(helper, envir = env)
  lapply(ids, function(i) env$oracle_run_case(i))
}

oracle_dead_record <- function(i, msg) {
  kind <- if (grepl("timeout|timed out|time limit", msg, ignore.case = TRUE)) "hang" else "crash"
  case <- oracle_case(i)
  list(id = i, family = case$family, api = case$api,
       structure = case$structure, storage = case$storage,
       spectrum = case$spectrum, target = case$target %||% case$which,
       method = case$method %||% "shim", k = case$k, n = case$n,
       m = case$m %||% case$n, tol = case$tol, status = kind,
       error = substr(msg, 1, 300),
       hard = sprintf("%s: subprocess died: %s", kind, substr(msg, 1, 200)),
       soft = character(), describe = oracle_describe(case))
}

# Re-run each case of a failed chunk alone so a crash or hang is attributed
# to its case id.
oracle_isolate <- function(ids, helper, case_timeout) {
  lapply(ids, function(i) {
    one <- tryCatch(
      callr::r(oracle_child, args = list(ids = i, helper = helper),
               env = c(callr::rcmd_safe_env(), oracle_env()),
               timeout = case_timeout, libpath = .libPaths()),
      error = function(e) e)
    if (inherits(one, "error")) oracle_dead_record(i, conditionMessage(one)) else one[[1]]
  })
}

# Runs ids in chunks; in subprocesses (callr) when possible so a crash or
# hang is attributed to its case, with up to `workers` chunks in parallel.
# Returns a list of records.
oracle_run_ids <- function(ids, subprocess = oracle_subprocess_available(),
                           workers = 1L, chunk = 20L, timeout = 900,
                           case_timeout = 180, progress = NULL) {
  chunks <- split(ids, ceiling(seq_along(ids) / chunk))
  if (!subprocess) {
    old <- options(eigencore.threads = 1L)
    on.exit(options(old), add = TRUE)
    out <- list()
    for (ch in chunks) {
      out <- c(out, lapply(ch, oracle_run_case))
      if (is.function(progress)) progress(length(out), length(ids))
    }
    return(out)
  }
  helper <- oracle_helper_path()
  results <- vector("list", length(chunks))
  pending <- seq_along(chunks)
  running <- list()
  done <- 0L
  while (length(pending) || length(running)) {
    while (length(pending) && length(running) < workers) {
      j <- pending[[1L]]
      pending <- pending[-1L]
      running[[as.character(j)]] <- list(
        proc = callr::r_bg(oracle_child, args = list(ids = chunks[[j]], helper = helper),
                           env = c(callr::rcmd_safe_env(), oracle_env()),
                           libpath = .libPaths(), supervise = TRUE),
        started = Sys.time())
    }
    finished <- character()
    for (key in names(running)) {
      job <- running[[key]]
      j <- as.integer(key)
      elapsed <- as.numeric(difftime(Sys.time(), job$started, units = "secs"))
      if (job$proc$is_alive()) {
        if (elapsed > timeout) {
          job$proc$kill()
          results[[j]] <- oracle_isolate(chunks[[j]], helper, case_timeout)
          finished <- c(finished, key)
        }
        next
      }
      res <- tryCatch(job$proc$get_result(), error = function(e) e)
      results[[j]] <- if (inherits(res, "error")) {
        oracle_isolate(chunks[[j]], helper, case_timeout)
      } else {
        res
      }
      finished <- c(finished, key)
    }
    for (key in finished) {
      done <- done + length(chunks[[as.integer(key)]])
      running[[key]] <- NULL
      if (is.function(progress)) progress(done, length(ids))
    }
    if (!length(finished)) Sys.sleep(0.2)
  }
  do.call(c, results)
}


oracle_reproducer <- function(rec) {
  sprintf(paste0("case %d [%s]\n  reproduce: source('tests/testthat/helper-oracle.R'); ",
                 "library(eigencore); str(oracle_run_case(%d))\n  %s"),
          rec$id, rec$describe, rec$id, rec$hard)
}

# Uncertified-rate table per configuration (family x storage x target x method).
oracle_summary_table <- function(df, by = c("family", "method", "target")) {
  ok <- df[df$status == "ok", , drop = FALSE]
  if (!nrow(ok)) return(data.frame())
  key <- do.call(paste, c(ok[by], sep = " / "))
  agg <- lapply(split(ok, key), function(d) {
    data.frame(config = paste(d[1, by], collapse = " / "),
               runs = nrow(d),
               certified = sum(d$certified %in% TRUE),
               uncertified_pct = round(100 * mean(!(d$certified %in% TRUE)), 1),
               wrong_set_certified = sum(d$certified %in% TRUE & d$set_ok %in% FALSE),
               median_value_err_uncert = signif(stats::median(
                 d$value_error[!(d$certified %in% TRUE)], na.rm = TRUE), 2),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, agg)
  out[order(-out$uncertified_pct, -out$runs), , drop = FALSE]
}
