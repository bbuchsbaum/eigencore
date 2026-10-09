#!/usr/bin/env Rscript
# Thread-count determinism check for eigencore's OpenMP kernels.
#
#   Rscript inst/validation/thread-determinism.R [threads...]   (default 1 2 4)
#
# Runs the same seeded sparse solves (CSC Lanczos, block thick restart,
# implicit-Gram / Golub-Kahan SVD, centred-scaled CSC, nonsymmetric
# Krylov-Schur, LOBPCG, randomized SVD) and raw operator applies at every
# thread count and compares each against the single-threaded result.
#
# Contract (?eigencore-threads): sparse products (A X, A^T X, centred and
# centred-scaled CSC, implicit normal operators) are bitwise identical for any
# thread count; whole solves may differ by rounding only, because with more
# than one thread the Lanczos reorthogonalisation switches from BLAS dgemm to
# OpenMP kernels (themselves independent of the thread count, so results with
# 2 and 4 threads are expected to agree bitwise). The script fails if an
# apply differs at all, if two multi-threaded runs differ, or if any solve
# differs by more than 1e-10 relative (values) / 1e-6 (vector subspace; the
# solves converge to tol = 1e-8, so vectors computed along different rounding
# paths legitimately differ by about tol / gap).

suppressPackageStartupMessages({
  library(eigencore)
  library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
thread_counts <- if (length(args)) as.integer(args) else c(1L, 2L, 4L)
ns <- asNamespace("eigencore")

make_inputs <- function() {
  set.seed(4242)
  list(
    S = {
      A <- Matrix::rsparsematrix(3000L, 3000L, density = 0.002)
      A <- A + t(A) + Matrix::Diagonal(3000L, x = seq(1, 50, length.out = 3000L))
      as(as(A, "generalMatrix"), "CsparseMatrix")
    },
    R = Matrix::rsparsematrix(4000L, 600L, density = 0.01),
    W = Matrix::rsparsematrix(400L, 6000L, density = 0.01),
    N = as(Matrix::rsparsematrix(1500L, 1500L, density = 0.004) +
             Matrix::Diagonal(1500L, x = seq(1, 10, length.out = 1500L)),
           "CsparseMatrix"),
    X = matrix(stats::rnorm(600L * 13L), 600L),
    Y = matrix(stats::rnorm(4000L * 13L), 4000L),
    Z = matrix(stats::rnorm(6000L * 3L), 6000L)
  )
}
inp <- make_inputs()

run_all <- function(threads) {
  options(eigencore.threads = threads)
  on.exit(options(eigencore.threads = NULL))
  ap <- ns$apply_operator
  apt <- ns$apply_adjoint_operator
  Rop <- as_operator(inp$R)
  cs <- scale_cols(center(inp$R), seq(0.5, 2, length.out = 600L))
  out <- list()
  # raw applies: repeated so the CSR / row-slab caches are built and used
  for (rep in 1:3) {
    out$apply_fwd <- ap(Rop, inp$X)
    out$apply_fwd1 <- ap(Rop, inp$X[, 1, drop = FALSE])
    out$apply_adj <- apt(Rop, inp$Y)
    out$apply_cs_fwd <- ap(cs, inp$X)
    out$apply_cs_adj <- apt(cs, inp$Y)
  }
  out$apply_repeat <- ns$csc_apply_repeat(inp$R, inp$X, reps = 4L)$Y
  out$apply_repeat_wide <- ns$csc_apply_repeat(
    inp$W, inp$Z, reps = 4L)$Y
  out$apply_repeat_cs <- ns$csc_apply_repeat(
    inp$R, inp$X, reps = 4L, col_means = Matrix::colMeans(inp$R),
    weights = seq(0.5, 2, length.out = 600L))$Y
  out$apply_repeat_adj <- ns$csc_apply_repeat(inp$R, inp$Y, transpose = TRUE,
                                              reps = 3L)$Y
  fit <- function(f) list(values = values(f), vectors = vectors(f))
  out$lanczos_csc <- fit(eig_partial(inp$S, 8L, seed = 1))
  out$lanczos_csc_small <- fit(eig_partial(inp$S, 6L, target = smallest(), seed = 2))
  out$block_lanczos <- fit(eig_partial(inp$S, 8L, method = lanczos(block = 4L), seed = 3))
  out$svd_tall <- fit(svd_partial(inp$R, 10L, seed = 4))
  out$svd_wide <- fit(svd_partial(inp$W, 10L, seed = 5))
  out$svd_gk <- fit(svd_partial(inp$R, 6L, method = golub_kahan(), seed = 6))
  out$svd_cs <- fit(svd_partial(cs, 6L, seed = 7))
  set.seed(11)  # svds() has no seed argument; its start vector uses the RNG
  out$svds <- svds(inp$R, 5L)[c("d", "u", "v")]
  out$arnoldi <- fit(eig_partial(inp$N, 4L, target = largest_magnitude(), seed = 8))
  out$lobpcg <- fit(eig_partial(inp$S, 4L, target = smallest(),
                                method = lobpcg(maxit = 400L), seed = 9))
  out$randomized <- fit(svd_partial(inp$R, 5L, method = randomized(), seed = 10))
  out
}

subspace_gap <- function(a, b) {
  if (is.null(a) || is.null(b)) return(0)
  a <- as.matrix(a); b <- as.matrix(b)
  if (is.complex(a) || is.complex(b)) {
    qa <- qr.Q(qr(a)); qb <- qr.Q(qr(b))
    return(max(Mod(qa %*% (Conj(t(qa)) %*% qb) - qb)))
  }
  # sign- and rotation-invariant: || (I - P_a) b ||
  qa <- qr.Q(qr(a))
  max(abs(b - qa %*% crossprod(qa, b)))
}

compare <- function(ref, cur) {
  rows <- list()
  for (name in names(ref)) {
    r <- ref[[name]]; c <- cur[[name]]
    if (startsWith(name, "apply")) {
      num <- function(x) if (is.list(x)) unlist(Filter(is.numeric, x)) else as.numeric(as.matrix(x))
      d <- max(abs(num(r) - num(c)))
      rows[[name]] <- data.frame(case = name, identical = identical(r, c),
                                 value_rel = d, subspace = NA_real_)
    } else {
      vr <- r$values %||% r$d; vc <- c$values %||% c$d
      rel <- max(Mod(vr - vc) / pmax(Mod(vr), 1e-300))
      vec_r <- r$vectors %||% r$u; vec_c <- c$vectors %||% c$u
      gap <- if (is.list(vec_r) && !is.matrix(vec_r)) {
        max(subspace_gap(vec_r$u, vec_c$u), subspace_gap(vec_r$v, vec_c$v))
      } else {
        subspace_gap(vec_r, vec_c)
      }
      rows[[name]] <- data.frame(case = name, identical = identical(r, c),
                                 value_rel = rel, subspace = gap)
    }
  }
  do.call(rbind, rows)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("eigencore threads available:", ns$eigencore_threads(), "\n")
results <- lapply(thread_counts, run_all)
names(results) <- paste0("t", thread_counts)
bad <- FALSE
for (i in seq_along(thread_counts)[-1L]) {
  cmp <- compare(results[[1L]], results[[i]])
  cat(sprintf("\n== threads %d vs %d ==\n", thread_counts[[i]], thread_counts[[1L]]))
  print(cmp, row.names = FALSE, digits = 3)
  apply_rows <- startsWith(cmp$case, "apply")
  if (!all(cmp$identical[apply_rows])) {
    message("FAIL: an operator apply is not bitwise identical across thread counts")
    bad <- TRUE
  }
  solve_rows <- !apply_rows
  if (any(cmp$value_rel[solve_rows] > 1e-10) ||
      any(cmp$subspace[solve_rows] > 1e-6, na.rm = TRUE)) {
    message("FAIL: a solve differs by more than rounding across thread counts")
    bad <- TRUE
  }
}
multi <- which(thread_counts > 1L)
if (length(multi) >= 2L) {
  cmp <- compare(results[[multi[1L]]], results[[multi[2L]]])
  cat(sprintf("\n== threads %d vs %d (both OpenMP) ==\n",
              thread_counts[[multi[2L]]], thread_counts[[multi[1L]]]))
  print(cmp, row.names = FALSE, digits = 3)
  if (!all(cmp$identical)) {
    message("FAIL: multi-threaded runs are not bitwise identical")
    bad <- TRUE
  }
}
if (bad) stop("thread determinism check failed")
cat("\nthread determinism check passed\n")
