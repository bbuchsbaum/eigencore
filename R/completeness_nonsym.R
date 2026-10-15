# Target-completeness probe for nonsymmetric (general real) eigensolvers.
#
# A right-residual certificate proves that each returned pair (lambda_j, x_j)
# is an eigenpair of A + E_j with small ||E_j||. It does not prove that the
# returned SET is the requested one: a restarted Arnoldi (Krylov-Schur) run
# can converge k wanted Ritz values while a more-preferred eigenvalue never
# emerged in its subspace (C60: random sparse n = 20000, largest_magnitude
# k = 6, the largest-modulus conjugate pair is missed by eigencore and by
# RSpectra alike, and every returned pair certifies).
#
# Theory. Let Q (n x d) be an orthonormal real basis of the span of the
# returned vectors (real and imaginary parts, so conjugate partners are
# included), T = Q' A Q and R = A Q - Q T. In the basis [Q, Q_perp]
#
#     A  =  [ T   G ]  +  [ 0  0 ]        E = Q_perp' A Q,  ||E|| <= ||R||,
#           [ 0   C ]     [ E  0 ]        C = Q_perp' A Q_perp,
#
# so when R = 0 (an exactly invariant subspace) spec(A) = spec(T) U spec(C):
# the eigenvalues missing from the returned set are exactly the eigenvalues
# of the compression C = (I - QQ') A (I - QQ') on the complement. With
# R != 0 the spectrum of A is the spectrum of the block-triangular matrix
# perturbed by ||E|| <= ||R||. For a normal A (more generally, for
# well-conditioned eigenvalues) the perturbation moves eigenvalues by
# O(||R||); for a non-normal A it is bounded only by the ||R||-pseudospectrum
# of the block-triangular matrix, which can be far larger than ||R|| (an
# eigenvalue with condition number kappa moves by up to about kappa * ||R||,
# a Jordan block of size m by ||R||^(1/m)). The probe's margin covers the
# normal/well-conditioned case only (see nonsym_completeness_margin()); it
# does not claim to resolve pseudospectral effects.
#
# The probe runs the native Krylov-Schur kernel (src/arnoldi.cpp) on C, i.e.
# Arnoldi with every operator product projected against Q, from a fixed-seed
# start that never touches the caller's RNG stream, and ranks its Ritz values
# with the request's own target. A complement Ritz value theta with Ritz
# residual r (for C) that is more preferred than the least-preferred returned
# value by more than margin + r indicates an intruder:
#   * converged theta (r below the solve tolerance): for normal A, C has an
#     eigenvalue within r of theta, so A has an eigenvalue more preferred than
#     the returned edge -- proof up to the normal-case perturbation bound;
#   * unconverged theta: a suspect only (Ritz values of a non-normal C can lie
#     anywhere in its field of values, outside its spectrum).
# Either way the result is repaired by a Schur-Rayleigh-Ritz merge: the
# Rayleigh-Ritz problem of A on span(Q, V_c), with V_c the probe's
# (deflated Krylov-Schur) complement basis, selected by the target, then
# re-certified from scratch and probed again, up to a small round budget. A
# merge that does not certify, or does not improve the returned set, leaves
# the result unchanged and reports "failed" (converged intruder) or
# "inconclusive" (suspect only).
#
# A clean probe -- the deflated Krylov-Schur run converged its wanted values
# (or exhausted the complement) and none is beyond the edge -- reports
# "probed". It is evidence, not proof: a short Krylov run on C can miss an
# eigenvalue exactly as the original solve did, and for non-normal A the
# margin does not cover pseudospectral shifts. When the complement is small
# enough the probe spans all of it ("exhausted"): its Ritz values are then the
# complete spectrum of C (computed by a backward-stable Schur decomposition),
# and only the R != 0 perturbation above stands between the probe and a
# proof. A probe that neither converges nor finds a suspect reports
# "inconclusive"; neither "inconclusive" nor "failed" passes.
#
# Shift-invert routes probe the transformed operator M = (A - sigma I)^-1 with
# the largest-magnitude target (nearest sigma in A's spectrum). The sparse
# general-pencil route probes B^-1 A, the operator its Arnoldi run used.

#' @keywords internal
nonsym_completeness_controls <- function() {
  as_count <- function(value, default, lo) {
    value <- suppressWarnings(as.integer(value %||% default))
    if (length(value) != 1L || is.na(value) || value < lo) default else value
  }
  list(
    # Wanted Ritz values of the complement probe and its Krylov-Schur
    # subspace; complements no larger than `exhaust` are spanned completely.
    wanted = as_count(getOption("eigencore.nonsym_probe_wanted"), 6L, 1L),
    subspace = as_count(getOption("eigencore.nonsym_probe_subspace"), 60L, 3L),
    exhaust = as_count(getOption("eigencore.nonsym_probe_exhaust"), 64L, 0L),
    restarts = as_count(getOption("eigencore.nonsym_probe_restarts"), 1000L, 1L),
    # Screening tolerance of the first (clean-verdict) probe pass.
    screen_tol = {
      value <- suppressWarnings(as.numeric(
        getOption("eigencore.nonsym_probe_screen_tol", 1e-4)))
      if (length(value) != 1L || !is.finite(value) || value <= 0) 1e-4 else value
    },
    max_rounds = as_count(getOption("eigencore.completeness_max_rounds"), 3L, 0L)
  )
}

# Targets the probe can rank (the native Krylov-Schur target codes).
#' @keywords internal
nonsym_completeness_target_kind <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else NA_character_
  switch(
    kind,
    largest = ,
    largest_real = "largest_real",
    smallest = ,
    smallest_real = "smallest_real",
    largest_magnitude = "largest_magnitude",
    smallest_magnitude = "smallest_magnitude",
    largest_imaginary = "largest_imaginary",
    smallest_imaginary = "smallest_imaginary",
    NULL
  )
}

# Preference key (larger = more preferred). Each key is 1-Lipschitz in the
# eigenvalue, so a margin in eigenvalue units is a margin in key units.
#' @keywords internal
nonsym_completeness_key <- function(z, kind) {
  switch(
    kind,
    largest_real = Re(z),
    smallest_real = -Re(z),
    largest_magnitude = Mod(z),
    smallest_magnitude = -Mod(z),
    largest_imaginary = Im(z),
    smallest_imaginary = -Im(z)
  )
}

#' @keywords internal
nonsym_completeness_target <- function(kind) {
  switch(
    kind,
    largest_real = largest_real(),
    smallest_real = smallest_real(),
    largest_magnitude = largest_magnitude(),
    smallest_magnitude = smallest_magnitude(),
    largest_imaginary = largest_imaginary(),
    smallest_imaginary = smallest_imaginary()
  )
}

# Margin in eigenvalue units: the solve tolerance relative to the spectrum
# scale, twice the invariance defect ||A Q - Q T||_F (the normal-case
# eigenvalue perturbation of the deflation), and a rounding floor.
#' @keywords internal
nonsym_completeness_margin <- function(values, theta, tol, coupling,
                                       norm_scale = NULL) {
  finite_abs <- function(x) {
    if (!is.numeric(x) && !is.complex(x)) {
      return(numeric())
    }
    x <- Mod(as.vector(x))
    x[is.finite(x)]
  }
  scale <- max(c(finite_abs(values), finite_abs(theta), finite_abs(norm_scale), 0))
  max(tol * scale, 2 * coupling, 64 * .Machine$double.eps * scale)
}

# Two-pass Gram-Schmidt orthonormalisation with a relative drop floor.
#' @keywords internal
nonsym_orthonormalize <- function(X, B = NULL) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  x0 <- sqrt(colSums(X * X))
  if (!is.null(B) && ncol(B)) {
    for (pass in 1:2) {
      X <- X - B %*% crossprod(B, X)
    }
  }
  kept <- 0L
  for (j in seq_len(ncol(X))) {
    if (!is.finite(x0[[j]]) || x0[[j]] == 0) {
      next
    }
    x <- X[, j]
    if (kept) {
      P <- X[, seq_len(kept), drop = FALSE]
      for (pass in 1:2) {
        x <- x - drop(P %*% crossprod(P, x))
      }
    }
    nx <- sqrt(sum(x * x))
    if (is.finite(nx) && nx > 1e-10 * x0[[j]]) {
      kept <- kept + 1L
      X[, kept] <- x / nx
    }
  }
  X[, seq_len(kept), drop = FALSE]
}

# Expected real dimension of span(Re X, Im X): one per real value, two per
# conjugate pair (a pair returned with both members, or one member alone).
#' @keywords internal
nonsym_expected_dimension <- function(values, scale) {
  values <- as.complex(values)
  real <- abs(Im(values)) <= 100 * .Machine$double.eps * max(scale, Mod(values))
  cplx <- values[!real]
  used <- rep(FALSE, length(cplx))
  dim <- sum(real)
  for (i in seq_along(cplx)) {
    if (used[[i]]) next
    used[[i]] <- TRUE
    partner <- which(!used & Mod(cplx - Conj(cplx[[i]])) <=
                       1e-8 * Mod(cplx[[i]]))
    if (length(partner)) used[[partner[[1L]]]] <- TRUE
    dim <- dim + 2L
  }
  as.integer(dim)
}

#' @keywords internal
nonsym_apply <- function(op, X) {
  out <- as.matrix(apply_operator(op, X))
  if (is.complex(out)) {
    if (max(abs(Im(out))) > 0) {
      stop("nonsymmetric completeness probe requires a real operator.", call. = FALSE)
    }
    out <- Re(out)
  }
  storage.mode(out) <- "double"
  out
}

# Orthonormal real basis Q of the returned (approximately invariant)
# subspace, T = Q' A Q and the invariance defect ||A Q - Q T||_F. `ok` is
# FALSE when the vectors do not span a subspace of the expected dimension
# whose Rayleigh quotient reproduces the returned values.
#' @keywords internal
nonsym_returned_basis <- function(op, values, vectors, tol) {
  X <- as.matrix(vectors)
  Xr <- if (is.complex(X)) cbind(Re(X), Im(X)) else X
  scale <- max(c(Mod(values), 0), na.rm = TRUE)
  Q <- nonsym_orthonormalize(Xr)
  expected <- nonsym_expected_dimension(values, scale)
  out <- list(Q = Q, ok = FALSE, columns = 0L, reason = NULL)
  if (ncol(Q) != expected) {
    out$reason <- sprintf(
      "returned vectors span %d dimensions, %d expected (linearly dependent or ill-conditioned eigenvector basis)",
      ncol(Q), expected)
    return(out)
  }
  AQ <- nonsym_apply(op, Q)
  out$columns <- ncol(Q)
  T <- crossprod(Q, AQ)
  R <- AQ - Q %*% T
  coupling <- sqrt(sum(R * R))
  ev <- eigen(T, only.values = TRUE)$values
  match_tol <- max(sqrt(max(tol, .Machine$double.eps)) * max(scale, 1e-300),
                   10 * coupling)
  # Multiset match: each returned value claims a distinct eigenvalue of T.
  free <- rep(TRUE, length(ev))
  dist <- numeric(length(values))
  for (j in seq_along(values)) {
    gap <- Mod(ev - as.complex(values[[j]]))
    gap[!free] <- Inf
    best <- which.min(gap)
    dist[[j]] <- gap[[best]]
    free[[best]] <- FALSE
  }
  if (any(!is.finite(dist)) || any(dist > match_tol)) {
    out$reason <- "the Rayleigh quotient of the returned basis does not reproduce the returned values"
    return(out)
  }
  out$ok <- TRUE
  out$coupling <- coupling
  out$unmatched <- ev[free]
  out
}

# Native kernel spec for the deflated Krylov-Schur entry point.
#' @keywords internal
nonsym_probe_kernel <- function(op) {
  kernel <- tryCatch(native_arnoldi_kernel(op), error = function(e) NULL)
  if (identical(kernel$kind, "dense") && is.double(kernel$matrix)) {
    return(list("dense", kernel$matrix))
  }
  if (identical(kernel$kind, "csc")) {
    source <- tryCatch(native_arnoldi_kernel_transpose(kernel),
                       error = function(e) NULL)
    if (inherits(source, "dgCMatrix")) {
      return(list("csc", methods::slot(source, "i"), methods::slot(source, "p"),
                  methods::slot(source, "x"), methods::slot(source, "Dim"), TRUE))
    }
  }
  closure <- function(X, alpha = 1, beta = 0, Y = NULL) {
    out <- alpha * nonsym_apply(op, X)
    if (is.null(Y) || beta == 0) out else out + beta * Y
  }
  list("r", as.integer(op$dim), closure)
}

# Deflated Krylov-Schur run on the compression of A to the complement of Q.
# Returns the complement Ritz values (most preferred first), their Ritz
# residuals for C, the Krylov-Schur basis of the kept complement subspace and
# whether the run converged / spanned the whole complement.
#' @keywords internal
nonsym_deflated_probe <- function(op, Q, kind, tol, controls, stream = 0L) {
  n <- as.integer(op$dim[[1L]])
  nc <- n - ncol(Q)
  exhausted <- nc <= controls$exhaust
  m <- if (exhausted) nc else min(nc, max(controls$subspace, 2L * controls$wanted + 1L))
  exhausted <- exhausted || m >= nc
  wanted <- min(controls$wanted, m)
  state <- saved_random_seed()
  on.exit(restore_random_seed(state), add = TRUE)
  set.seed(20261010L + as.integer(stream), kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  start <- stats::rnorm(n)
  if (ncol(Q)) {
    for (pass in 1:2) start <- start - drop(Q %*% crossprod(Q, start))
  }
  cycle <- .Call(
    "eigencore_arnoldi_ks_deflated", nonsym_probe_kernel(op), as.numeric(start),
    as.integer(wanted), as.integer(m),
    arnoldi_target_code(nonsym_completeness_target(kind)),
    as.numeric(max(tol, 100 * .Machine$double.eps)),
    as.integer(if (exhausted) 1L else controls$restarts),
    Q, PACKAGE = "eigencore"
  )
  p <- cycle$iterations
  Tp <- cycle$H[seq_len(p), seq_len(p), drop = FALSE]
  b <- cycle$H[p + 1L, seq_len(p)]
  eig <- eigen(Tp)
  yn <- sqrt(colSums(Mod(eig$vectors)^2))
  res <- Mod(colSums(eig$vectors * b)) / pmax(yn, .Machine$double.eps)
  if (exhausted) {
    res[] <- 0
  }
  key <- nonsym_completeness_key(eig$values, kind)
  ord <- order(key, decreasing = TRUE)
  list(
    theta = eig$values[ord],
    residuals = res[ord],
    basis = cycle$V[, seq_len(p), drop = FALSE],
    converged = isTRUE(cycle$converged) || exhausted,
    exhausted = exhausted,
    matvecs = as.integer(cycle$matvecs),
    restarts = as.integer(cycle$krylov_schur_restarts),
    subspace = as.integer(m)
  )
}

# Schur-Rayleigh-Ritz merge: Rayleigh-Ritz of A on span(Q, V_c), the k most
# preferred Ritz pairs by `target`, re-certified by `certify_fn`.
#' @keywords internal
nonsym_completeness_merge <- function(op, Q, Vc, k, target, tol, certify_fn) {
  Z <- cbind(Q, nonsym_orthonormalize(Vc, Q))
  AZ <- nonsym_apply(op, Z)
  H <- crossprod(Z, AZ)
  eig <- eigen(H)
  idx <- utils::head(order_indices(eig$values, target), k)
  values <- eig$values[idx]
  vectors <- Z %*% eig$vectors[, idx, drop = FALSE]
  norms <- sqrt(colSums(Mod(vectors)^2))
  vectors <- sweep(vectors, 2L, pmax(norms, .Machine$double.eps), `/`)
  if (all(arnoldi_real_values(values, arnoldi_projected_scale(H)))) {
    values <- Re(values)
    vectors <- arnoldi_realify_vectors(vectors)
  }
  list(values = values, vectors = vectors,
       certificate = certify_fn(values, vectors),
       columns = ncol(Z))
}

# Probe (and if needed repair) a nonsymmetric partial eigensolution of `op`.
# `values`/`vectors` are eigenpairs of `op` (the probed operator); `target`
# ranks them; `certify_fn(values, vectors)` re-certifies a merged set in the
# caller's coordinates. Returns the verdict, the record, and the possibly
# repaired values, vectors and certificate.
#' @keywords internal
nonsym_completeness_check <- function(op, values, vectors, target, tol,
                                      certify_fn, norm_scale = NULL,
                                      controls = nonsym_completeness_controls()) {
  kind <- nonsym_completeness_target_kind(target)
  started <- proc.time()[["elapsed"]]
  k <- length(values)
  record <- list(
    method = "deflated_krylov_schur_probe",
    wanted = controls$wanted,
    subspace = NA_integer_,
    restarts = 0L,
    operator_columns = 0L,
    operator_block_calls = 0L,
    rounds = 0L,
    edge = NA_real_,
    most_preferred_complement = NA_complex_,
    margin = NA_real_,
    invariance_defect = NA_real_,
    intruder_found = FALSE,
    probe_converged = NA,
    exhausted = FALSE,
    evidence = "probe",
    non_normal_caveat = paste(
      "margin covers eigenvalue shifts of the invariance defect for normal or",
      "well-conditioned eigenvalues only; non-normal pseudospectral shifts are not bounded"
    )
  )
  certificate <- NULL
  repaired <- FALSE
  status <- "inconclusive"
  round <- 0L
  repeat {
    basis <- nonsym_returned_basis(op, values, vectors, tol)
    record$operator_columns <- record$operator_columns + basis$columns
    record$operator_block_calls <- record$operator_block_calls + (basis$columns > 0L)
    if (!isTRUE(basis$ok)) {
      record$reason <- basis$reason
      status <- if (repaired) "failed" else "inconclusive"
      break
    }
    Q <- basis$Q
    record$invariance_defect <- basis$coupling
    edge <- min(nonsym_completeness_key(values, kind))
    record$edge <- edge
    # Candidates: the complement Ritz values, and the eigenvalues of T that no
    # returned value claimed (conjugate partners of returned complex values,
    # deflated with them; for imaginary-part targets a partner can be more
    # preferred than the returned edge). A suspect must lie beyond the edge by
    # the margin plus its own Ritz residual.
    assess <- function(probe) {
      theta <- c(as.complex(basis$unmatched), as.complex(probe$theta))
      theta_res <- c(rep(0, length(basis$unmatched)), probe$residuals)
      margin <- nonsym_completeness_margin(values, theta, tol, basis$coupling,
                                           norm_scale)
      key <- nonsym_completeness_key(theta, kind)
      suspect <- key > edge + margin + theta_res
      theta_scale <- pmax(Mod(theta), .Machine$double.eps^(1 / 3) *
                            max(c(Mod(theta), 1e-300)))
      converged <- theta_res <= max(tol, 100 * .Machine$double.eps) * theta_scale
      list(theta = theta, margin = margin, key = key, suspect = suspect,
           confirmed = suspect & converged)
    }
    run_probe <- function(probe_tol, restarts = controls$restarts) {
      probe_controls <- controls
      probe_controls$restarts <- restarts
      probe <- if (ncol(Q) >= op$dim[[1L]]) {
        # The returned subspace is the whole space: the complement is empty.
        list(theta = complex(), residuals = numeric(),
             basis = matrix(0, op$dim[[1L]], 0L), converged = TRUE,
             exhausted = TRUE, matvecs = 0L, restarts = 0L, subspace = 0L)
      } else {
        nonsym_deflated_probe(op, Q, kind, probe_tol, probe_controls, stream = round)
      }
      record$subspace <<- probe$subspace
      record$restarts <<- record$restarts + probe$restarts
      record$operator_columns <<- record$operator_columns + probe$matvecs
      record$operator_block_calls <<- record$operator_block_calls + probe$matvecs
      record$probe_converged <<- probe$converged
      record$exhausted <<- probe$exhausted
      probe
    }
    # Screen at a loose tolerance with a short restart budget (a clean
    # verdict only needs the complement values resolved to within the gap);
    # re-run at the solve tolerance when the screen finds a suspect (so a
    # repair merges accurate complement vectors) or does not converge (the
    # Krylov-Schur restart selection can stagnate at a loose tolerance on
    # heavily tied keys, e.g. imaginary-part targets over a real complement).
    loose_tol <- max(tol, min(sqrt(tol), controls$screen_tol))
    screened <- loose_tol > tol
    probe <- run_probe(loose_tol, if (screened) min(controls$restarts, 100L) else controls$restarts)
    verdict <- assess(probe)
    # P23: when a suspect of the screen already converged at the solve
    # tolerance (Krylov-Schur resolves the dominant complement values far
    # below the screening tolerance), the screen's basis carries an accurate
    # intruder vector: try the merge with it first, and re-run the probe at
    # the solve tolerance only if that merge does not certify and improve.
    # The merged set is re-certified from scratch and probed again by the
    # next round, which re-examines any suspect the merge left out.
    screen_merge <- screened && !isTRUE(probe$exhausted) &&
      isTRUE(probe$converged) && any(verdict$suspect) &&
      any(verdict$confirmed)
    if (screened && !screen_merge && !isTRUE(probe$exhausted) &&
        (any(verdict$suspect) || !isTRUE(probe$converged))) {
      probe <- run_probe(tol)
      verdict <- assess(probe)
      screened <- FALSE
    }
    theta <- verdict$theta
    margin <- verdict$margin
    suspect <- verdict$suspect
    confirmed <- verdict$confirmed
    record$margin <- margin
    record$probe_tolerance <- if (screened) loose_tol else tol
    record$most_preferred_complement <- if (length(theta)) {
      theta[[which.max(verdict$key)]]
    } else {
      NA_complex_
    }
    if (!any(suspect)) {
      status <- if (!probe$converged) "inconclusive" else if (repaired) "repaired" else "probed"
      if (!probe$converged) {
        record$reason <- "complement probe did not converge within its restart budget"
      }
      break
    }
    record$intruder_found <- TRUE
    record$intruder_confirmed <- any(confirmed)
    if (round >= controls$max_rounds) {
      status <- if (any(confirmed)) "failed" else "inconclusive"
      record$reason <- "intruder still present after the repair round budget"
      break
    }
    round <- round + 1L
    # Accept a merge when its sorted preference keys dominate the old ones
    # and at least one improves by more than the margin (the edge itself may
    # stay, e.g. when k cuts a conjugate pair).
    try_merge <- function(probe, margin) {
      merged <- nonsym_completeness_merge(op, Q, probe$basis, k,
                                          nonsym_completeness_target(kind), tol,
                                          certify_fn)
      record$operator_columns <<- record$operator_columns + merged$columns + k
      record$operator_block_calls <<- record$operator_block_calls + 2L
      merged$improved <- length(merged$values) == k && {
        old_keys <- sort(nonsym_completeness_key(values, kind), decreasing = TRUE)
        new_keys <- sort(nonsym_completeness_key(merged$values, kind), decreasing = TRUE)
        all(new_keys >= old_keys - margin) && any(new_keys > old_keys + margin)
      }
      merged
    }
    merged <- try_merge(probe, margin)
    if ((!isTRUE(merged$certificate$passed) || !isTRUE(merged$improved)) &&
        screen_merge) {
      # The screen-basis merge fell short: fall back to the solve-tolerance
      # probe and its merge, as without the shortcut.
      probe <- run_probe(tol)
      verdict <- assess(probe)
      margin <- verdict$margin
      confirmed <- verdict$confirmed
      record$margin <- margin
      record$probe_tolerance <- tol
      if (length(verdict$theta)) {
        record$most_preferred_complement <- verdict$theta[[which.max(verdict$key)]]
      }
      if (!any(verdict$suspect)) {
        round <- round - 1L
        record$intruder_found <- FALSE
        record$intruder_confirmed <- NULL
        status <- if (!probe$converged) "inconclusive" else if (repaired) "repaired" else "probed"
        if (!probe$converged) {
          record$reason <- "complement probe did not converge within its restart budget"
        }
        break
      }
      record$intruder_confirmed <- any(confirmed)
      merged <- try_merge(probe, margin)
    }
    improved <- isTRUE(merged$improved)
    if (!isTRUE(merged$certificate$passed) || !improved) {
      status <- if (any(confirmed)) "failed" else "inconclusive"
      record$reason <- if (!isTRUE(merged$certificate$passed)) {
        "the Schur-Rayleigh-Ritz merge with the complement basis did not certify"
      } else {
        "the Schur-Rayleigh-Ritz merge did not improve the returned set"
      }
      break
    }
    values <- merged$values
    vectors <- merged$vectors
    certificate <- merged$certificate
    repaired <- TRUE
  }
  record$rounds <- round
  record$seconds <- proc.time()[["elapsed"]] - started
  list(status = status, record = record, repaired = repaired,
       values = values, vectors = vectors, certificate = certificate)
}

# Attach the verdict to a certificate (residual_passed keeps the residual
# meaning; passed is withheld for failed / inconclusive).
#' @keywords internal
nonsym_certificate_with_completeness <- function(cert, status, record) {
  if (is.null(cert)) {
    return(cert)
  }
  cert$residual_passed <- cert$residual_passed %||% isTRUE(cert$passed)
  cert$target_completeness <- status
  cert$target_passed <- switch(status, probed = TRUE, repaired = TRUE,
                               failed = FALSE, NA)
  cert$completeness <- record
  note <- switch(
    status,
    failed = "nonsymmetric completeness probe found a more-preferred eigenvalue outside the returned set that could not be resolved",
    inconclusive = paste0("nonsymmetric completeness probe inconclusive: ",
                          record$reason %||% "no verdict"),
    repaired = "nonsymmetric completeness probe found a more-preferred eigenvalue outside the returned set; the set was repaired by a Schur-Rayleigh-Ritz merge with the deflated complement",
    NULL
  )
  if (status %in% c("failed", "inconclusive")) {
    cert$passed <- FALSE
  }
  if (!is.null(note)) {
    cert$notes <- unique(c(cert$notes, note))
  }
  cert
}

# Is the nonsymmetric probe applicable to this problem / result?
#' @keywords internal
nonsym_completeness_eligible <- function(op, values, vectors, k, target, certify,
                                         certificate) {
  mode <- tryCatch(target_completeness_mode(NULL), error = function(e) "auto")
  isTRUE(certify) && !identical(mode, "none") &&
    !is.null(nonsym_completeness_target_kind(target)) &&
    isTRUE(certificate$passed) && is.null(certificate$target_completeness) &&
    !is.null(vectors) && length(values) == k &&
    identical(op$dtype %||% "double", "double")
}

#' @keywords internal
nonsym_completeness_mode_allows <- function(plan) {
  mode <- tryCatch(target_completeness_mode(plan$method_descriptor),
                   error = function(e) "auto")
  !identical(mode, "none")
}

# Probe a nonsymmetric result of the probed operator `op` (eigenvalues
# `values` of `op`, ranked by `target`) and return the possibly repaired
# values, vectors and certificate with the verdict attached. `certify_fn`
# re-certifies a merged set in the caller's coordinates. Results that are
# not eligible (certification off, completeness mode "none", residual
# certificate not passed, unsupported target, complex operator) are returned
# unchanged and stay "not_checked".
#' @keywords internal
nonsym_completeness_apply <- function(op, values, vectors, certificate, k,
                                      target, tol, certify, plan, certify_fn,
                                      norm_scale = certificate$scale) {
  out <- list(values = values, vectors = vectors, certificate = certificate,
              repaired = FALSE, columns = 0L)
  if (!nonsym_completeness_mode_allows(plan) ||
      !nonsym_completeness_eligible(op, values, vectors, k, target, certify,
                                    certificate)) {
    return(out)
  }
  check <- tryCatch(
    nonsym_completeness_check(op, values, vectors, target, tol,
                              certify_fn = certify_fn,
                              norm_scale = norm_scale),
    error = function(e) list(
      status = "inconclusive", repaired = FALSE,
      record = list(method = "deflated_krylov_schur_probe",
                    reason = paste0("probe error: ", conditionMessage(e)))
    )
  )
  if (isTRUE(check$repaired) && !is.null(check$certificate)) {
    out$values <- check$values
    out$vectors <- check$vectors
    out$certificate <- check$certificate
    out$repaired <- TRUE
  }
  out$certificate <- nonsym_certificate_with_completeness(
    out$certificate, check$status, check$record)
  out$columns <- as.integer(check$record$operator_columns %||% 0L)
  out
}

# Hook for solve_eigen_arnoldi(): probe the native / reference Arnoldi result
# on the operator it solved. `iter` must carry vectors.
#' @keywords internal
nonsym_completeness_after_arnoldi <- function(op, iter, k, target, tol, certify,
                                              plan) {
  cert_op <- arnoldi_certificate_operator(op)
  probed <- nonsym_completeness_apply(
    op, iter$values, iter$vectors, iter$certificate, k, target, tol, certify,
    plan,
    certify_fn = function(values, vectors) {
      certify_general_eigen_operator(cert_op, values, vectors, tol = tol)
    }
  )
  iter$values <- probed$values
  iter$vectors <- probed$vectors
  iter$certificate <- probed$certificate
  iter$matvecs <- as.integer((iter$matvecs %||% 0L) + probed$columns)
  iter
}
