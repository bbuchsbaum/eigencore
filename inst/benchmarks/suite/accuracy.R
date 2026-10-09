# Benchmark suite: independent accuracy checks and trusted references.
#
# Every method's output (eigencore included) is checked by this file, not by
# the method itself:
#
#   residuals      eigen:  r_i = A x_i - lambda_i B x_i          (B = I if absent)
#                  SVD:    r_i = sqrt(||A v_i - s_i u_i||^2 + ||A' u_i - s_i v_i||^2)
#                          (two-sided; u_i, v_i normalised)
#                  nonsymmetric eigen: right residual; left residual
#                          ||A' y_i - conj(lambda_i) y_i|| when left vectors exist
#   backward error eigen:  ||r_i|| / ((||A||_2 + |lambda_i| ||B||_2) ||x_i||)
#                  SVD:    r_i / ||A||_2
#                  i.e. the normwise 2-norm backward error, the same definition
#                  eigencore's certificates use since the C12 change. The suite
#                  divides by a high-accuracy ||A||_2 (exact or a converged
#                  Golub-Kahan-Lanczos estimate), not by a lower bound.
#   orthogonality  max |X' X - I| (X' B X - I for generalized; U and V for SVD)
#   values         matched against a trusted reference: the analytic spectrum
#                  when known, dense LAPACK when feasible, otherwise a
#                  tight certified iterative reference (see
#                  suite_certified_reference()), cross-checked.
#                  value_err = max_i |lambda_i - lambda_ref_i| / (||A||_2 + |lambda_ref_i| ||B||_2)
#   target_ok      the method returned k values and every reference value in
#                  the target set was matched within a gap-aware tolerance
#                  tau = clamp(gap/2, 1e-10, 1e-6) (gap = smallest scaled gap
#                  between distinct reference values, including the first
#                  values outside the target set).

# Bump when the reference computation changes (part of the cache key).
SUITE_REFERENCE_VERSION <- 3L

# Operator closures for A (with optional column centring).
suite_ops <- function(prob) {
  A <- prob$A
  mu <- prob$mu
  if (is.null(mu)) {
    list(m = nrow(A), n = ncol(A),
         mv = function(X) as.matrix(A %*% X),
         tmv = function(Y) as.matrix(Matrix::crossprod(A, Y)))
  } else {
    list(m = nrow(A), n = ncol(A),
         mv = function(X) { X <- as.matrix(X); as.matrix(A %*% X) - matrix(colSums(mu * X), nrow(A), ncol(X), byrow = TRUE) },
         tmv = function(Y) { Y <- as.matrix(Y); as.matrix(Matrix::crossprod(A, Y)) - outer(mu, colSums(Y)) })
  }
}

# Largest singular value by Golub-Kahan-Lanczos bidiagonalisation with full
# reorthogonalisation, run until the estimate is stable to ~1e-15.
suite_norm2_lanczos <- function(mv, tmv, m, n, steps = 80L, seed = 1L) {
  suite_set_seed(seed)
  steps <- min(steps, m, n)
  V <- matrix(0, n, steps)
  U <- matrix(0, m, steps)
  alpha <- beta <- numeric(steps)
  v <- stats::rnorm(n)
  v <- v / sqrt(sum(v^2))
  prev <- 0
  stable <- 0L
  est <- NA_real_
  for (j in seq_len(steps)) {
    V[, j] <- v
    u <- as.numeric(mv(v))
    if (j > 1L) u <- u - beta[j - 1L] * U[, j - 1L]
    if (j > 1L) u <- u - U[, 1:(j - 1L), drop = FALSE] %*% crossprod(U[, 1:(j - 1L), drop = FALSE], u)
    alpha[j] <- sqrt(sum(u^2))
    if (alpha[j] == 0) break
    u <- as.numeric(u / alpha[j])
    U[, j] <- u
    w <- as.numeric(tmv(u)) - alpha[j] * v
    w <- w - V[, 1:j, drop = FALSE] %*% crossprod(V[, 1:j, drop = FALSE], w)
    beta[j] <- sqrt(sum(w^2))
    Bj <- matrix(0, j + 1L, j)
    for (i in seq_len(j)) {
      Bj[i, i] <- alpha[i]
      Bj[i + 1L, i] <- beta[i]
    }
    # ||B_j|| = ||A V_j|| <= ||A||_2: a monotone lower bound that converges fast
    est <- svd(Bj, nu = 0, nv = 0)$d[1L]
    stable <- if (abs(est - prev) <= 1e-15 * est) stable + 1L else 0L
    prev <- est
    if ((j >= 5L && stable >= 2L) || beta[j] == 0) break
    v <- as.numeric(w / beta[j])
  }
  list(value = est, steps = j)
}

suite_fingerprint <- function(A) {
  if (inherits(A, "sparseMatrix")) {
    x <- A@x
    nz <- length(x)
  } else {
    x <- as.vector(A)
    nz <- sum(x != 0)
  }
  sprintf("%dx%d-nnz%d-s%.12e-a%.12e", nrow(A), ncol(A), nz, sum(x), sum(abs(x)))
}

suite_target_sort <- function(values, case, extra = 0L) {
  idx <- suite_select(values, case$target, case$k + extra, case$sigma)
  values[idx]
}

# ---- references ---------------------------------------------------------

suite_reference <- function(case, prob, cache_dir = NULL, verbose = TRUE) {
  key <- paste0(case$id, "-g", SUITE_GENERATOR_VERSION, "-r", SUITE_REFERENCE_VERSION, "-",
                substr(suite_md5_string(suite_fingerprint(prob$A)), 1L, 10L))
  cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(key, ".rds")) else NULL
  if (!is.null(cache_file) && file.exists(cache_file)) {
    ref <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(ref)) {
      ref$cached <- TRUE
      return(ref)
    }
  }
  t0 <- suite_now()
  ref <- suite_compute_reference(case, prob)
  ref$seconds <- suite_now() - t0
  ref$key <- key
  ref$cached <- FALSE
  if (!is.null(cache_file)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    try(saveRDS(ref, cache_file), silent = TRUE)
  }
  ref
}

suite_compute_reference <- function(case, prob) {
  A <- prob$A
  ops <- suite_ops(prob)
  n <- ncol(A)
  extra <- 3L
  normB <- 1
  normB_source <- "identity"
  norm_lanczos <- function() {
    est <- suite_norm2_lanczos(ops$mv, ops$tmv, ops$m, ops$n, seed = case$seed + 7L)
    list(value = est$value, source = sprintf("Golub-Kahan-Lanczos (%d steps, full reorth.)", est$steps))
  }
  if (!is.null(prob$B)) {
    Bm <- prob$B
    eB <- suite_norm2_lanczos(function(x) as.matrix(Bm %*% x), function(x) as.matrix(Bm %*% x),
                              n, n, seed = case$seed + 9L)
    normB <- eB$value
    normB_source <- "Golub-Kahan-Lanczos"
  }

  vals <- NULL
  source <- NULL
  crosscheck <- NA_real_
  ref_be <- NA_real_
  ambiguous <- FALSE
  nrm <- NULL

  if (!is.null(prob$exact)) {
    vals <- if (case$task == "full") sort(prob$exact) else suite_target_sort(prob$exact, case, extra)
    source <- "analytic spectrum"
    if (is.null(prob$B) && case$task %in% c("sym", "full")) {
      nrm <- list(value = max(abs(prob$exact)), source = "exact (analytic spectrum)")
    }
  } else if (case$task %in% c("sym", "full") && n <= suite_limits$base_sym_n) {
    ev <- eigen(as.matrix(A), symmetric = TRUE, only.values = TRUE)$values
    vals <- if (case$task == "full") sort(ev) else suite_target_sort(ev, case, extra)
    source <- "dense LAPACK eigen()"
    nrm <- list(value = max(abs(ev)), source = "exact (dense eigen)")
  } else if (case$task == "nonsym" && n <= suite_limits$base_nonsym_n) {
    Ad <- as.matrix(A)
    ev <- eigen(Ad, only.values = TRUE)$values
    vals <- suite_target_sort(ev, case, extra)
    source <- "dense LAPACK eigen()"
    nrm <- list(value = svd(Ad, nu = 0, nv = 0)$d[1L], source = "exact (dense svd)")
  } else if (case$task == "gen" && n <= suite_limits$base_gen_n) {
    R <- chol(as.matrix(prob$B))
    C <- backsolve(R, t(backsolve(R, as.matrix(A), transpose = TRUE)), transpose = TRUE)
    ev <- eigen((C + t(C)) / 2, symmetric = TRUE, only.values = TRUE)$values
    vals <- suite_target_sort(ev, case, extra)
    source <- "dense LAPACK (Cholesky-reduced) eigen()"
  } else if (case$task == "svd" && min(dim(A)) <= 3000L) {
    small_cols <- ncol(A) <= nrow(A)
    if (is.null(prob$mu)) {
      G <- as.matrix(if (small_cols) Matrix::crossprod(A) else Matrix::tcrossprod(A))
    } else {
      # centred Gram (columns centred): X'X with X = A - 1 mu'
      if (!small_cols) stop("centred wide SVD reference not implemented")
      G <- as.matrix(Matrix::crossprod(A)) - nrow(A) * tcrossprod(prob$mu)
    }
    ev <- eigen((G + t(G)) / 2, symmetric = TRUE, only.values = TRUE)$values
    sv <- sqrt(pmax(ev, 0))
    vals <- sv[seq_len(min(length(sv), case$k + extra))]
    source <- "dense LAPACK eigen() of the Gram matrix (top singular values)"
    nrm <- list(value = sv[1L], source = "exact (dense Gram eigen)")
  }

  if (is.null(nrm)) nrm <- norm_lanczos()

  if (is.null(vals)) {
    cert <- suite_certified_reference(case, prob, extra, nrm$value, normB)
    vals <- cert$values
    source <- cert$source
    crosscheck <- cert$crosscheck
    ref_be <- cert$backward_error
    ambiguous <- cert$ambiguous
  }

  list(values = vals, source = source, norm2 = nrm$value, norm2_source = nrm$source,
       normB = normB, normB_source = normB_source, crosscheck = crosscheck,
       reference_backward_error = ref_be, ambiguous = ambiguous,
       fingerprint = suite_fingerprint(A))
}

# Tight iterative reference when no analytic or dense reference is feasible.
#
# Candidates: RSpectra at tol 1e-13 asking for k+10 values with a large
# subspace, RSpectra again asking for k+20 values with an even larger one
# (extremal sets with nearly equal keys, e.g. the circular-law edge of a random
# nonsymmetric matrix, are only found reliably that way), and eigencore at tol
# 1e-12. Every candidate is verified with the independent backward error;
# among the accurate ones (<= 1e-9) the one whose k wanted values are most
# extreme for the target wins: a value with a tiny backward error is a genuine
# eigenvalue, so a set that reaches further cannot be wrong where a set that
# stops short can have skipped values. Disagreement between candidates is
# recorded as `crosscheck` (scaled value distance of the top k).
suite_certified_reference <- function(case, prob, extra, norm2, normB) {
  k <- case$k
  tight <- 1e-13
  fake_ref <- list(norm2 = norm2, normB = normB, values = NULL)
  rs_run <- function(kk, ncv) {
    A <- prob$A
    o <- list(tol = tight, maxitr = 20000L, ncv = min(ncv, min(dim(A)) - 1L))
    suite_set_seed(case$seed)
    res <- switch(case$task,
      sym = RSpectra::eigs_sym(A, kk, which = if (case$target == "near") "LM" else case$target,
                               sigma = if (case$target == "near") case$sigma else NULL, opts = o),
      nonsym = RSpectra::eigs(A, kk, which = "LM", opts = o),
      svd = if (is.null(prob$mu)) RSpectra::svds(A, kk, nu = kk, nv = kk, opts = o) else {
        mu <- prob$mu
        RSpectra::svds(function(x, args) as.numeric(A %*% x) - sum(mu * x), kk, nu = kk, nv = kk,
                       Atrans = function(y, args) as.numeric(Matrix::crossprod(A, y)) - mu * sum(y),
                       dim = dim(A), opts = o)
      })
    if (case$task == "svd") suite_std(values = res$d, u = res$u, v = res$v)
    else suite_std(values = res$values, vectors = res$vectors)
  }
  cands <- list()
  if (requireNamespace("RSpectra", quietly = TRUE) && case$task != "gen") {
    k1 <- k + max(extra, 10L)
    k2 <- k + 20L
    r1 <- tryCatch(rs_run(k1, max(3L * k1, 60L)), error = function(e) NULL)
    if (!is.null(r1)) cands[[sprintf("RSpectra(k=%d, tol=%g)", k1, tight)]] <- r1
    r2 <- tryCatch(rs_run(k2, max(4L * k2, 100L)), error = function(e) NULL)
    if (!is.null(r2)) cands[[sprintf("RSpectra(k=%d, tol=%g)", k2, tight)]] <- r2
  }
  if (requireNamespace("eigencore", quietly = TRUE)) {
    tcase <- case
    tcase$k <- k + max(extra, 10L)
    ad <- suite_adapter_eigencore(tcase, prob, 1e-12)
    res <- tryCatch(ad$extract(ad$run()), error = function(e) NULL)
    if (!is.null(res)) cands[[sprintf("eigencore(k=%d, tol=1e-12)", tcase$k)]] <- res
  }
  if (!length(cands)) stop("no method available to compute a reference")
  info <- lapply(cands, function(r) {
    tc <- case
    tc$k <- length(r$values)
    acc <- suite_accuracy(tc, prob, r, fake_ref)
    vals <- suite_target_sort(r$values, tc)
    list(vals = vals, be = acc$max_backward_error,
         score = if (length(vals) >= k) sum(suite_target_key(vals[seq_len(k)], case)) else Inf)
  })
  be <- vapply(info, `[[`, numeric(1), "be")
  score <- vapply(info, `[[`, numeric(1), "score")
  ok <- is.finite(be) & be <= 1e-9 & is.finite(score)
  if (!any(ok)) ok <- is.finite(score)
  best <- names(cands)[ok][which.min(score[ok])]
  vals <- info[[best]]$vals
  others <- setdiff(names(cands), best)
  cross <- if (length(others)) max(vapply(others, function(o) {
    ov <- info[[o]]$vals
    m <- min(length(ov), k)
    if (m == 0) return(Inf)
    suite_match_values(ov[seq_len(m)], vals[seq_len(k)], norm2, normB)$err
  }, numeric(1))) else NA_real_
  list(values = vals,
       source = sprintf("certified iterative: %s, independent backward error %.1e; candidates %s; max top-k disagreement %.1e",
                        best, be[[best]],
                        paste(sprintf("%s (bwd %.1e)", names(cands), be), collapse = ", "), cross),
       crosscheck = cross, backward_error = be[[best]],
       # No candidate was accurate: there is no trusted reference.
       ambiguous = !any(is.finite(be) & be <= 1e-9))
}

# Greedy nearest matching: each computed value (in the order given) takes the
# nearest unused reference value. Returns the scaled and relative errors and
# the matched reference indices (NA when the reference set ran out).
suite_match_values <- function(computed, reference, norm2, normB = 1) {
  if (!length(computed) || !length(reference)) return(list(err = NA_real_, rel = NA_real_, idx = integer()))
  if (!is.complex(computed) && !is.complex(reference) && length(computed) > 60L) {
    m <- min(length(computed), length(reference))
    a <- sort(computed)[seq_len(m)]
    b <- sort(reference)[seq_len(m)]
    d <- abs(a - b)
    return(list(err = max(d / (norm2 + abs(b) * normB)), rel = max(d / pmax(abs(b), .Machine$double.xmin)),
                idx = seq_len(m)))
  }
  used <- rep(FALSE, length(reference))
  err <- rel <- numeric(0)
  idx <- integer(0)
  for (c in computed) {
    if (all(used)) {
      err <- c(err, Inf); rel <- c(rel, Inf); idx <- c(idx, NA_integer_)
      next
    }
    d <- abs(reference - c)
    d[used] <- Inf
    j <- which.min(d)
    used[j] <- TRUE
    idx <- c(idx, j)
    err <- c(err, d[j] / (norm2 + abs(reference[j]) * normB))
    rel <- c(rel, d[j] / max(abs(reference[j]), .Machine$double.xmin))
  }
  list(err = max(err), rel = max(rel), idx = idx)
}

# Ordering key of a target: smaller is "more wanted".
suite_target_key <- function(v, case) {
  switch(case$target, LA = -Re(v), SA = Re(v), LM = -Mod(v), near = abs(v - case$sigma),
         top = -Re(v), all = Re(v))
}

suite_gap_tau <- function(ref_values, norm2, normB = 1) {
  v <- ref_values
  if (length(v) < 2L) return(1e-6)
  d <- as.matrix(stats::dist(cbind(Re(v), Im(v))))
  d <- d[upper.tri(d)]
  scale <- norm2 + max(abs(v)) * normB
  d <- d[d > 1e-13 * scale]
  if (!length(d)) return(1e-6)
  min(1e-6, max(1e-10, 0.5 * min(d) / scale))
}

colnorms <- function(X) sqrt(colSums(Mod(X)^2))

suite_apply_complex <- function(f, X) {
  if (is.complex(X)) f(Re(X)) + 1i * f(Im(X)) else f(X)
}

suite_accuracy <- function(case, prob, std, ref) {
  out <- list(returned_k = length(std$values), max_residual = NA_real_,
              max_backward_error = NA_real_, max_left_residual = NA_real_,
              orthogonality_loss = NA_real_, value_err = NA_real_,
              value_rel_err = NA_real_, target_ok = NA, target_tau = NA_real_)
  ops <- suite_ops(prob)
  norm2 <- ref$norm2
  normB <- ref$normB %||% 1
  lam <- std$values
  if (case$task == "svd") {
    if (!is.null(std$u) && !is.null(std$v) && length(lam)) {
      kk <- length(lam)
      U <- as.matrix(std$u)[, seq_len(kk), drop = FALSE]
      V <- as.matrix(std$v)[, seq_len(kk), drop = FALSE]
      U <- sweep(U, 2L, colnorms(U), "/")
      V <- sweep(V, 2L, colnorms(V), "/")
      rr <- colnorms(ops$mv(V) - sweep(U, 2L, lam, "*"))
      rl <- colnorms(ops$tmv(U) - sweep(V, 2L, lam, "*"))
      r <- sqrt(rr^2 + rl^2)
      out$max_residual <- max(r)
      out$max_backward_error <- max(r / norm2)
      out$max_left_residual <- max(rl)
      out$orthogonality_loss <- max(abs(crossprod(as.matrix(std$u)[, seq_len(kk), drop = FALSE]) - diag(kk)),
                                    abs(crossprod(as.matrix(std$v)[, seq_len(kk), drop = FALSE]) - diag(kk)))
    }
  } else if (!is.null(std$vectors) && length(lam)) {
    X <- as.matrix(std$vectors)
    kk <- length(lam)
    X <- X[, seq_len(kk), drop = FALSE]
    AX <- suite_apply_complex(ops$mv, X)
    BX <- if (is.null(prob$B)) X else suite_apply_complex(function(Z) as.matrix(prob$B %*% Z), X)
    R <- AX - sweep(BX, 2L, lam, "*")
    rn <- colnorms(R)
    out$max_residual <- max(rn / colnorms(X))
    out$max_backward_error <- max(rn / ((norm2 + Mod(lam) * normB) * colnorms(X)))
    if (case$task != "nonsym") {
      G <- if (is.null(prob$B)) crossprod(X) else crossprod(X, as.matrix(BX))
      out$orthogonality_loss <- max(abs(G - diag(kk)))
    }
    if (!is.null(std$left) && case$task == "nonsym") {
      Y <- as.matrix(std$left)[, seq_len(kk), drop = FALSE]
      RL <- suite_apply_complex(ops$tmv, Y) - sweep(Y, 2L, Conj(lam), "*")
      out$max_left_residual <- max(colnorms(RL) / colnorms(Y))
    }
  }
  if (!is.null(ref$values) && !length(lam)) {
    out$target_ok <- FALSE
  } else if (!is.null(ref$values)) {
    if (case$task == "full") {
      k_ref <- length(ref$values)
      mv <- suite_match_values(lam, ref$values, norm2, normB)
      tau <- 1e-10
      in_target <- TRUE
    } else {
      k_ref <- case$k
      # Match against the reference target set plus a few values beyond it,
      # then require every matched reference value to lie inside the target
      # set; ties at the boundary (repeated values, conjugate pairs) count as
      # inside.
      mv <- suite_match_values(lam, ref$values, norm2, normB)
      tau <- suite_gap_tau(ref$values, norm2, normB)
      kth <- ref$values[min(k_ref, length(ref$values))]
      bound <- suite_target_key(kth, case) + tau * (norm2 + Mod(kth) * normB)
      in_target <- length(mv$idx) && !anyNA(mv$idx) &&
        all(suite_target_key(ref$values[mv$idx], case) <= bound)
    }
    out$value_err <- mv$err
    out$value_rel_err <- mv$rel
    out$target_tau <- tau
    out$target_ok <- length(lam) >= k_ref && is.finite(mv$err) && mv$err <= tau && isTRUE(in_target)
    # No trusted reference set: report the value error, but no verdict.
    if (isTRUE(ref$ambiguous)) out$target_ok <- NA
  }
  out
}
