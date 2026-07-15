#!/usr/bin/env Rscript
#
# bench-readme.R — reproduces the comparison table shown in README.Rmd.
#
# This is a transparent, self-contained head-to-head against the partial-solver
# baselines a user would otherwise reach for (RSpectra and irlba). It uses
# only the *public* eigencore API. The README table reports eigencore's full
# certified call. Optional console ratios compare that call with external raw
# solver calls at matched tolerance; they are labelled raw diagnostics and are
# not time-to-certified-answer claims. Numbers are machine/BLAS dependent;
# rerun locally to get your own.
#
#   Rscript inst/benchmarks/bench-readme.R
#
# It is a documentation aid, not a release gate. The strict, installed-package
# release gates live in the other inst/benchmarks/bench-*.R scripts and in
# docs/v1-benchmark-manifest.md.

suppressMessages({
  library(eigencore)
  library(Matrix)
})
have_rspectra <- requireNamespace("RSpectra", quietly = TRUE)
have_irlba    <- requireNamespace("irlba", quietly = TRUE)

ITERS <- 5L

timeit <- function(expr, iters = ITERS) {
  expr <- substitute(expr)
  env <- parent.frame()
  mark <- tryCatch(
    bench::mark(eval(expr, env), iterations = iters, check = FALSE,
                time_unit = "s", memory = TRUE, filter_gc = FALSE),
    error = function(e) NULL
  )
  if (is.null(mark)) return(list(t = NA_real_, mem = NA_real_))
  list(t = as.numeric(stats::median(mark$time[[1L]])),
       mem = as.numeric(mark$mem_alloc[[1L]]))
}

ms     <- function(x) if (is.finite(x)) sprintf("%.0f ms", x * 1000) else "—"
speedx <- function(base, ref) if (is.finite(base) && is.finite(ref) && ref > 0) sprintf("%.1f×", base / ref) else "—"
memx   <- function(base, ref) if (is.finite(base) && is.finite(ref) && ref > 0) sprintf("%.1f×", base / ref) else "—"

rows <- list()
record <- function(label, eig_fit, ec, rs = NULL, ir = NULL) {
  rows[[length(rows) + 1L]] <<- data.frame(
    problem = label,
    method  = eig_fit$method,
    passed  = isTRUE(eig_fit$certificate$passed),
    norm_bound = eig_fit$certificate$norm_bound_type,
    eigencore_ms = ec$t * 1000,
    rspectra_raw_ratio = if (!is.null(rs)) rs$t / ec$t else NA_real_,
    irlba_raw_ratio    = if (!is.null(ir)) ir$t / ec$t else NA_real_,
    irlba_raw_mem_ratio = if (!is.null(ir)) ir$mem / ec$mem else NA_real_,
    stringsAsFactors = FALSE
  )
  cat(sprintf("\n[%s]\n  path: %s | cert passed=%s (%s)\n",
              label, eig_fit$method, isTRUE(eig_fit$certificate$passed),
              eig_fit$certificate$norm_bound_type))
  cat(sprintf("  eigencore : %s\n", ms(ec$t)))
  if (!is.null(rs)) cat(sprintf("  RSpectra raw: %s   (raw/eigencore-certified ratio %s)\n", ms(rs$t), speedx(rs$t, ec$t)))
  if (!is.null(ir)) cat(sprintf("  irlba raw   : %s   (raw/eigencore-certified ratio %s, raw memory ratio %s)\n", ms(ir$t), speedx(ir$t, ec$t), memx(ir$mem, ec$mem)))
}

set.seed(11)

## ---- Tall/wide sparse SVD: the certified Gram special case ----------------
svd_case <- function(label, M, r = 10L) {
  tol <- 1e-8
  fit <- svd_partial(M, rank = r, target = largest(), tol = tol)
  ec <- timeit(svd_partial(M, rank = r, target = largest(), tol = tol))
  rs <- if (have_rspectra) {
    timeit(RSpectra::svds(
      M,
      k = r,
      nu = r,
      nv = r,
      opts = list(tol = tol, maxitr = 1000L)
    ))
  } else NULL
  ir <- if (have_irlba) timeit(irlba::irlba(M, nv = r, nu = r, tol = tol)) else NULL
  record(label, fit, ec, rs, ir)
}
svd_case("Tall sparse SVD  100000 x 500, k=10",
         as(rsparsematrix(100000, 500, density = 0.002), "dgCMatrix"))
svd_case("Wide sparse SVD  500 x 100000, k=10",
         as(rsparsematrix(500, 100000, density = 0.002), "dgCMatrix"))

## ---- Banded/structured Hermitian smallest eigenvalues ---------------------
path_laplacian <- function(n) {
  as(bandSparse(n, n, k = c(-1, 0, 1),
                diagonals = list(rep(-1, n - 1), rep(2, n), rep(-1, n - 1))),
     "dgCMatrix")
}
herm_case <- function(label, A, k = 8L) {
  fit <- eig_partial(A, k = k, target = smallest())
  ec <- timeit(eig_partial(A, k = k, target = smallest()))
  rs <- NULL
  conv <- NA_integer_; t_sa <- NA_real_; t_si <- NA_real_
  if (have_rspectra) {
    conv <- length(suppressWarnings(RSpectra::eigs_sym(A, k = k, which = "SA"))$values)
    t_sa <- timeit(suppressWarnings(RSpectra::eigs_sym(A, k = k, which = "SA")))$t
    t_si <- timeit(RSpectra::eigs_sym(A, k = k, which = "LM", sigma = 0))$t
    # Only feed the README table a speedup when the baseline actually solved
    # the problem; a ratio against a non-converged run is not a fair claim.
    if (identical(conv, k)) rs <- list(t = t_sa, mem = NA_real_)
  }
  record(label, fit, ec, rs = rs)
  if (have_rspectra) {
    cat(sprintf('  RSpectra which="SA"      : %s (converged %d of %d)\n', ms(t_sa), conv, k))
    cat(sprintf('  RSpectra shift-invert    : %s (which="LM", sigma=0)\n', ms(t_si)))
  }
}
herm_case("Banded Hermitian smallest  n=20000, k=8", path_laplacian(20000L))

## ---- Operators without densifying: centered tall sparse SVD ---------------
cat("\n[Centered tall sparse SVD without densifying]\n")
M <- as(rsparsematrix(50000, 600, density = 0.01), "dgCMatrix")
op <- center(M, columns = TRUE)
dense_mb <- as.numeric(nrow(M)) * ncol(M) * 8 / 1e6
fitc <- svd_partial(op, rank = 10L, target = largest())
ecc <- timeit(svd_partial(op, rank = 10L, target = largest()), iters = 3L)
dmm <- timeit({ Md <- scale(as.matrix(M), center = TRUE, scale = FALSE)
                svd_partial(Md, rank = 10L, target = largest()) }, iters = 3L)
cat(sprintf("  dense centered matrix would occupy : %.0f MB (never formed)\n", dense_mb))
cat(sprintf("  sparse M object size               : %.1f MB\n", as.numeric(object.size(M)) / 1e6))
cat(sprintf("  center() operator path             : %s\n", fitc$method))
cat(sprintf("  certificate passed=%s (%s)\n", isTRUE(fitc$certificate$passed), fitc$certificate$norm_bound_type))
cat(sprintf("  operator solve allocation          : %.0f MB\n", ecc$mem / 1e6))
cat(sprintf("  densify-then-solve allocation      : %.0f MB  (%.1fx more)\n",
            dmm$mem / 1e6, dmm$mem / ecc$mem))

## ---- Markdown table for the README ----------------------------------------
tab <- do.call(rbind, rows)
cat("\n\n================= README table (markdown) =================\n\n")
cat("| Problem (all certified by eigencore) | eigencore | vs RSpectra raw | vs irlba raw |\n")
cat("|---|---|---|---|\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("| %s | %s | %s | %s |\n",
              tab$problem[i], ms(tab$eigencore_ms[i] / 1000),
              if (is.finite(tab$rspectra_raw_ratio[i])) sprintf("%.1f×", tab$rspectra_raw_ratio[i]) else "—",
              if (is.finite(tab$irlba_raw_ratio[i])) sprintf("%.1f×", tab$irlba_raw_ratio[i]) else "—"))
}
cat("\nR ", as.character(getRversion()), " | ",
    R.version$platform, " | ",
    "RSpectra ", if (have_rspectra) as.character(packageVersion("RSpectra")) else "NA",
    ", irlba ", if (have_irlba) as.character(packageVersion("irlba")) else "NA", "\n", sep = "")
