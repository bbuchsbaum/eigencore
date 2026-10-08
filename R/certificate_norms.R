# Spectral-norm (2-norm) bounds for certificate backward errors (C12).
#
# Certificates report the normwise backward error
#
#   eigen: eta_i = ||A x_i - lambda_i B x_i|| / ((||A||_2 + |lambda_i| ||B||_2) ||x_i||)
#   SVD:   eta_i = sqrt(||A v_i - s_i u_i||^2 + ||A^H u_i - s_i v_i||^2) / ||A||_2
#
# ||A||_2 is rarely known exactly, so every denominator uses a value L with
# L <= ||A||_2 (up to rounding). A lower bound in the denominator makes the
# reported eta an OVER-estimate of the true normwise backward error, so
# `passed` stays sound: a certificate that passes with L also passes with the
# exact norm. Bounds are collected, cheapest first:
#
#   * exact values: identity (1), diagonal storage (max |d|), a full computed
#     spectrum, or a user-asserted `metadata$two_norm`;
#   * structural: the largest column norm ||A e_j|| (dense, CSC, tridiagonal),
#     or ||A||_F / sqrt(min(m, n)) when only a Frobenius norm is known;
#   * the solve's own results: ||A x_i|| / ||x_i|| for the certified vectors
#     ("applied_vectors"), or the residual-corrected Ritz value
#     |lambda_i| ||B x_i|| / ||x_i|| - r_i / ||x_i|| ("ritz");
#   * a short deterministic Krylov estimate ("lanczos"): 24 steps of Lanczos
#     on A (Hermitian) or A^H A, verified as ||A y|| / ||y|| for the extreme
#     Ritz vector y. It runs only in the "gray zone" -- when some pair fails
#     with the free bounds but could pass with a known upper bound (or no
#     upper bound is known) -- and is memoised per operator. It uses a fixed
#     quasi-random start vector and never touches the R RNG stream.

#' @keywords internal
cert_norm_bound <- function(value, exact = FALSE, source = "none") {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value < 0) {
    value <- 0
    exact <- FALSE
    source <- "none"
  }
  list(value = value, exact = isTRUE(exact), source = source)
}

#' @keywords internal
cert_norm_max <- function(...) {
  bounds <- list(...)
  best <- cert_norm_bound(0)
  for (b in bounds) {
    if (is.null(b)) {
      next
    }
    if (isTRUE(b$exact) && !isTRUE(best$exact)) {
      # An exact value dominates every valid lower bound (up to rounding).
      best <- b
      next
    }
    if (isTRUE(best$exact)) {
      next
    }
    if (b$value > best$value) {
      best <- b
    }
  }
  best
}

# A "norm target" wraps an operator or a plain matrix behind the handful of
# capabilities the bound machinery needs.
#' @keywords internal
cert_norm_target <- function(x, hermitian = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "cert_norm_target")) {
    return(x)
  }
  if (inherits(x, "eigencore_operator")) {
    dim <- as.integer(x$dim)
    herm <- if (is.null(hermitian)) {
      identical(x$structure$kind, "hermitian") && dim[[1L]] == dim[[2L]]
    } else {
      isTRUE(hermitian)
    }
    adjoint_fn <- if (!is.null(x$apply_adjoint)) {
      function(X) apply_adjoint_operator(x, X)
    } else if (herm) {
      function(X) apply_operator(x, X)
    } else {
      NULL
    }
    out <- list(
      kind = "operator",
      op = x,
      dim = dim,
      hermitian = herm,
      apply = function(X) apply_operator(x, X),
      adjoint = adjoint_fn
    )
  } else {
    dim <- as.integer(dim(x))
    complex_x <- is.complex(x)
    out <- list(
      kind = "matrix",
      matrix = x,
      dim = dim,
      hermitian = isTRUE(hermitian) && dim[[1L]] == dim[[2L]],
      apply = function(X) as.matrix(x %*% X),
      adjoint = function(X) {
        if (complex_x) {
          as.matrix(Conj(t(x)) %*% X)
        } else {
          as.matrix(Matrix::crossprod(x, X))
        }
      }
    )
  }
  class(out) <- "cert_norm_target"
  out
}

#' @keywords internal
cert_norm_memo <- function(target, key, value) {
  if (identical(target$kind, "operator")) {
    return(operator_memoised_value(target$op, key, value))
  }
  value
}

# Exact value or cheap structural lower bound for ||A||_2 (no operator
# applies). Used as the pre-solve scale handed to native kernels and as the
# floor of every certificate denominator.
#' @keywords internal
two_norm_structural_bound <- function(x) {
  if (is.null(x)) {
    return(cert_norm_bound(1, exact = TRUE, source = "identity"))
  }
  target <- cert_norm_target(x)
  if (identical(target$kind, "operator")) {
    op <- target$op
    meta <- op$metadata
    # Exact name match: $ would partially match e.g. two_norm_upper.
    asserted <- meta[["two_norm", exact = TRUE]]
    if (is.numeric(asserted) && length(asserted) == 1L &&
        is.finite(asserted) && asserted >= 0) {
      return(cert_norm_bound(asserted, exact = TRUE, source = "metadata"))
    }
    if (identical(meta[["storage", exact = TRUE]] %||% NULL, "ddiMatrix")) {
      return(matrix_two_norm_structural(meta[["matrix", exact = TRUE]]))
    }
    return(cert_norm_memo(target, "two_norm_structural_bound",
                          operator_two_norm_structural_compute(op)))
  }
  matrix_two_norm_structural(target$matrix)
}

#' @keywords internal
operator_two_norm_structural_compute <- function(op) {
  meta <- op$metadata
  src <- source_or_null(op)
  if (!is.null(src)) {
    return(matrix_two_norm_structural(src))
  }
  if (identical(meta[["storage", exact = TRUE]] %||% NULL, "dgCMatrix")) {
    css <- meta[["column_sum_squares", exact = TRUE]] %||% NULL
    if (is.numeric(css) && length(css)) {
      return(cert_norm_bound(sqrt(max(css)), source = "column_norms"))
    }
    if (!is.null(meta[["matrix", exact = TRUE]])) {
      return(matrix_two_norm_structural(meta[["matrix", exact = TRUE]]))
    }
  }
  frob <- meta[["frobenius_norm", exact = TRUE]] %||% NULL
  if (is.numeric(frob) && length(frob) == 1L && is.finite(frob) && frob >= 0) {
    # sum(sigma_i^2) <= rank * sigma_1^2 and rank <= min(m, n).
    r <- max(1, min(op$dim))
    return(cert_norm_bound(frob / sqrt(r), source = "frobenius_rank_bound"))
  }
  cert_norm_bound(0)
}

#' @keywords internal
matrix_two_norm_structural <- function(x) {
  if (inherits(x, "ddiMatrix")) {
    unit <- identical(methods::slot(x, "diag"), "U")
    vals <- if (unit) {
      if (nrow(x) > 0L) 1 else 0
    } else {
      methods::slot(x, "x")
    }
    value <- if (length(vals)) max(abs(vals)) else 0
    return(cert_norm_bound(value, exact = TRUE, source = "diagonal"))
  }
  if (!length(x)) {
    return(cert_norm_bound(0, exact = TRUE, source = "empty"))
  }
  value <- if (inherits(x, "sparseMatrix")) {
    sqrt(max(Matrix::colSums(abs(x)^2)))
  } else {
    x <- as.matrix(x)
    max(col_norms(x))
  }
  cert_norm_bound(value, source = "column_norms")
}

# A cheap upper bound on ||A||_2 (or Inf when none is known), used only to
# decide whether refining the lower bound could change a pass/fail decision.
#' @keywords internal
two_norm_upper_bound <- function(x) {
  if (is.null(x)) {
    return(1)
  }
  target <- cert_norm_target(x)
  structural <- two_norm_structural_bound(target)
  if (isTRUE(structural$exact)) {
    return(structural$value)
  }
  if (identical(target$kind, "operator")) {
    op <- target$op
    frob <- op$metadata[["frobenius_norm", exact = TRUE]] %||% NULL
    if (is.numeric(frob) && length(frob) == 1L && is.finite(frob)) {
      return(frob)
    }
    src <- source_or_null(op) %||% NULL
    if (!is.null(src)) {
      return(cert_norm_memo(target, "frobenius_norm_exact", matrix_norm(src)))
    }
    return(Inf)
  }
  matrix_norm(target$matrix)
}

# Deterministic, RNG-free start vector with a nonzero mean and generic
# components (golden-ratio sequence).
#' @keywords internal
two_norm_start_vector <- function(n) {
  x <- ((seq_len(n) * 0.6180339887498949) %% 1) - 0.25
  x / sqrt(sum(x^2))
}

#' @keywords internal
cert_inner <- function(Q, w) {
  if (is.complex(Q) || is.complex(w)) {
    return(Conj(t(Q)) %*% w)
  }
  crossprod(Q, w)
}

# Short Krylov lower bound on ||A||_2: Lanczos with full reorthogonalisation
# on A (Hermitian) or A^H A, then the extreme Ritz vector y is verified as
# ||A y|| / ||y||, a rigorous lower bound regardless of orthogonality. The
# running max of ||A q_j|| over the unit basis vectors is also a lower bound.
# Without an adjoint, a short power iteration on A is used instead.
#' @keywords internal
two_norm_krylov_bound <- function(x, steps = 24L) {
  target <- cert_norm_target(x)
  cert_norm_memo(target, "two_norm_krylov_lower_bound",
                 cert_norm_bound(two_norm_krylov_compute(target, steps),
                                 source = "lanczos"))
}

#' @keywords internal
two_norm_krylov_compute <- function(target, steps = 24L) {
  n <- target$dim[[2L]]
  if (!n || !target$dim[[1L]]) {
    return(0)
  }
  steps <- as.integer(min(steps, n))
  q <- matrix(two_norm_start_vector(n), ncol = 1L)
  best <- 0
  ratio <- function(Y, X) {
    xn <- col_norms(X)
    yn <- col_norms(Y)
    ok <- is.finite(xn) & xn > 0 & is.finite(yn)
    if (!any(ok)) 0 else max(yn[ok] / xn[ok])
  }
  if (is.null(target$adjoint) && !isTRUE(target$hermitian)) {
    # Power iteration on a square operator without an adjoint:
    # ||A q|| / ||q|| <= ||A||_2 for every iterate.
    if (target$dim[[1L]] != n) {
      return(0)
    }
    for (j in seq_len(steps)) {
      z <- as.matrix(target$apply(q))
      best <- max(best, ratio(z, q))
      zn <- sqrt(sum(Mod(z)^2))
      if (!is.finite(zn) || zn <= 0) {
        break
      }
      q <- z / zn
    }
    return(best)
  }
  normal <- !isTRUE(target$hermitian)
  Q <- NULL
  Z <- NULL
  for (j in seq_len(steps)) {
    z <- as.matrix(target$apply(q))
    Q <- cbind(Q, q)
    Z <- cbind(Z, z)
    if (j == steps) {
      break
    }
    w <- if (normal) as.matrix(target$adjoint(z)) else z
    if (!all(is.finite(Mod(w)))) {
      break
    }
    w_before <- sqrt(sum(Mod(w)^2))
    for (pass in 1:2) {
      w <- w - Q %*% cert_inner(Q, w)
    }
    wn <- sqrt(sum(Mod(w)^2))
    # Stop once the Krylov space is (numerically) invariant.
    if (!is.finite(wn) ||
        wn <= 1e3 * .Machine$double.eps * max(w_before, .Machine$double.xmin)) {
      break
    }
    q <- w / wn
  }
  if (is.null(Z)) {
    return(0)
  }
  best <- max(best, ratio(Z, Q))
  # Rayleigh-Ritz on A^H A over span(Q): the top eigenvector s of Z^H Z gives
  # y = Q s with A y = Z s (A is linear), so ||Z s|| / ||Q s|| <= ||A||_2.
  G <- cert_inner(Z, Z)
  G <- (G + Conj(t(G))) / 2
  ev <- tryCatch(eigen(G, symmetric = TRUE), error = function(e) NULL)
  if (!is.null(ev)) {
    s <- ev$vectors[, 1L, drop = FALSE]
    best <- max(best, ratio(Z %*% s, Q %*% s))
  }
  best
}

# Resolve the norm bounds behind a certificate.
#
# `sides` is a named list with entries `A` (and optionally `B`), each a list
# with `target` (operator, matrix, or NULL for the identity) and `free`, a list
# of cert_norm_bound() values derived from the solve's own results.
# `backward_fn(norms)` maps a named numeric vector c(A = , B = ) to
# list(backward = , scale = ). Refinement runs only when it could flip a
# failing pair to passing.
#' @keywords internal
resolve_certificate_norms <- function(sides, backward_fn, tol) {
  infos <- lapply(sides, function(side) {
    if (is.null(side$target)) {
      return(cert_norm_bound(1, exact = TRUE, source = "identity"))
    }
    structural <- side$structural %||% two_norm_structural_bound(side$target)
    do.call(cert_norm_max, c(list(structural), side$free %||% list()))
  })
  norms_of <- function(infos) vapply(infos, function(i) i$value, numeric(1))
  current <- backward_fn(norms_of(infos))
  failing <- !is.finite(current$backward) | current$backward > tol
  refinable <- !vapply(infos, function(i) isTRUE(i$exact), logical(1))
  if (any(failing) && any(refinable) && length(tol) == 1L && is.finite(tol)) {
    uppers <- mapply(function(side, info) {
      if (isTRUE(info$exact)) info$value else two_norm_upper_bound(side$target)
    }, sides, infos)
    names(uppers) <- names(sides)
    could_pass <- if (all(is.finite(uppers))) {
      at_upper <- backward_fn(uppers)$backward
      any(failing & is.finite(at_upper) & at_upper <= tol)
    } else {
      TRUE
    }
    if (isTRUE(could_pass)) {
      for (nm in names(sides)[refinable]) {
        infos[[nm]] <- cert_norm_max(infos[[nm]],
                                     two_norm_krylov_bound(sides[[nm]]$target))
      }
      current <- backward_fn(norms_of(infos))
    }
  }
  types <- vapply(names(infos), function(nm) {
    info <- infos[[nm]]
    if (identical(info$source, "identity")) {
      "identity_exact"
    } else if (isTRUE(info$exact)) {
      "two_norm_exact"
    } else {
      "two_norm_lower_bound"
    }
  }, character(1))
  sources <- vapply(infos, function(i) i$source, character(1))
  frob <- frobenius_if_known(sides$A$target)
  list(
    norms = norms_of(infos),
    backward = current$backward,
    scale = current$scale,
    norm_bound_type = paste(types, collapse = "+"),
    norm_source = paste(sources, collapse = "+"),
    frobenius_norm = frob
  )
}

#' @keywords internal
frobenius_if_known <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  target <- cert_norm_target(x)
  if (identical(target$kind, "operator")) {
    frob <- target$op$metadata[["frobenius_norm", exact = TRUE]] %||% NULL
    if (is.numeric(frob) && length(frob) == 1L) {
      return(as.numeric(frob))
    }
  }
  NA_real_
}

# Free lower bounds from the solve's own results -------------------------

#' @keywords internal
bound_from_ratios <- function(numer, denom, source) {
  ok <- is.finite(numer) & is.finite(denom) & denom > 0
  if (!any(ok)) {
    return(NULL)
  }
  cert_norm_bound(max(0, max(numer[ok] / denom[ok])), source = source)
}

# ||A x|| >= |lambda| ||B x|| - ||r||, so (|lambda| ||Bx|| - r) / ||x|| is a
# lower bound on ||A||_2 needing only the residual.
#' @keywords internal
eigen_ritz_bound <- function(values, residuals, vec_norms, bv_norms = vec_norms) {
  if (!length(values)) {
    return(NULL)
  }
  bound_from_ratios(Mod(values) * bv_norms - residuals, vec_norms, "ritz")
}

#' @keywords internal
svd_ritz_bound <- function(d, left, right, u_norms, v_norms) {
  if (!length(d)) {
    return(NULL)
  }
  d <- Mod(d)
  bound_from_ratios(
    c(d * u_norms - left, d * v_norms - right),
    c(v_norms, u_norms),
    "ritz"
  )
}

#' @keywords internal
applied_bound <- function(value) {
  if (is.null(value) || !length(value)) {
    return(NULL)
  }
  value <- suppressWarnings(as.numeric(value[[1L]]))
  if (!is.finite(value) || value <= 0) {
    return(NULL)
  }
  cert_norm_bound(value, source = "applied_vectors")
}

# Backward-error assembly ---------------------------------------------------

#' @keywords internal
eigen_two_norm_backward <- function(A, values, residuals, vec_norms, tol,
                                    B = NULL, has_B = !is.null(B),
                                    free_A = list(), free_B = list(),
                                    structural_A = NULL, structural_B = NULL,
                                    finite = NULL) {
  values <- as.vector(values)
  residuals <- as.numeric(residuals)
  idx <- if (is.null(finite)) seq_along(values) else which(finite)
  sides <- list(A = list(target = A, free = free_A, structural = structural_A))
  if (isTRUE(has_B)) {
    sides$B <- list(target = B, free = free_B, structural = structural_B)
  } else {
    sides$B <- list(target = NULL)
  }
  backward_fn <- function(norms) {
    scale <- rep(Inf, length(values))
    backward <- rep(Inf, length(values))
    if (length(idx)) {
      scale[idx] <- pmax(
        (norms[["A"]] + Mod(values[idx]) * norms[["B"]]) *
          pmax(vec_norms[idx], .Machine$double.eps),
        .Machine$double.eps
      )
      backward[idx] <- residuals[idx] / scale[idx]
    }
    list(backward = backward, scale = scale)
  }
  resolve_certificate_norms(sides, backward_fn, tol)
}

#' @keywords internal
svd_two_norm_backward <- function(A, combined, tol, free = list(),
                                  structural = NULL) {
  combined <- as.numeric(combined)
  backward_fn <- function(norms) {
    scale <- rep(max(norms[["A"]], .Machine$double.eps), length(combined))
    list(backward = combined / scale, scale = scale)
  }
  resolve_certificate_norms(
    list(A = list(target = A, free = free, structural = structural)),
    backward_fn, tol
  )
}

# A structural column-norm bound computed by a native kernel (NULL if absent).
#' @keywords internal
native_column_bound <- function(value) {
  if (is.null(value) || !length(value)) {
    return(NULL)
  }
  cert_norm_bound(value, source = "column_norms")
}

# Build an SVD certificate from residual norms (computed from A in original
# coordinates by the caller) with two-norm scaling.
#' @keywords internal
svd_certificate_from_residuals <- function(A, d, left, right, orthogonality,
                                           tol, u = NULL, v = NULL,
                                           applied = NULL, notes = character(),
                                           certificate_type = "residual_backward_error",
                                           structural = NULL) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  combined <- sqrt(left^2 + right^2)
  free <- list(applied_bound(applied))
  if (!is.null(u) && !is.null(v)) {
    free <- c(free, list(svd_ritz_bound(d, left, right, col_norms(u), col_norms(v))))
  }
  norm <- svd_two_norm_backward(A, combined, tol, free = free,
                                structural = structural)
  new_certificate(
    tol = tol,
    residuals = list(left = left, right = right, combined = combined),
    backward_error = norm$backward,
    orthogonality = orthogonality,
    converged = is.finite(norm$backward) & norm$backward <= tol,
    scale = norm$scale,
    notes = notes,
    certificate_type = certificate_type,
    norm_bound_type = norm$norm_bound_type,
    norm_source = norm$norm_source,
    norm_values = norm$norms,
    frobenius_norm = norm$frobenius_norm
  )
}

# Same for native SVD diagnostics lists (left/right/orthogonality and,
# when the kernel reports it, norm_A_applied_bound). Native scale fields are
# ignored: some kernels still compute internal scales differently.
#' @keywords internal
svd_certificate_from_native_diagnostics <- function(A, d, diag, tol,
                                                    u = NULL, v = NULL,
                                                    swap_sides = FALSE,
                                                    notes = character()) {
  if (is.null(diag)) {
    return(NULL)
  }
  left <- diag$left %||% diag$residuals$left
  right <- diag$right %||% diag$residuals$right
  orthogonality <- diag$orthogonality
  if (isTRUE(swap_sides)) {
    tmp <- left
    left <- right
    right <- tmp
    if (length(orthogonality) >= 2L) {
      orthogonality <- c(
        U = unname(orthogonality[[2L]]),
        V = unname(orthogonality[[1L]])
      )
    }
  }
  svd_certificate_from_residuals(
    A, d, left, right, orthogonality, tol,
    u = u, v = v, applied = diag$norm_A_applied_bound, notes = notes,
    structural = native_column_bound(diag$norm_A_column_bound)
  )
}
