# Target-completeness probe for Hermitian Krylov eigensolvers (C50).
#
# A residual certificate proves that every returned pair is an eigenpair to
# the stated backward error. It does not prove that the returned SET is the
# requested one: a single-vector Krylov method (and a block method whose block
# is smaller than an eigenvalue's multiplicity) only sees, in exact
# arithmetic, the component of its start vector in each eigenspace, so a third
# copy of an exactly repeated eigenvalue can be missing while the k returned
# pairs (9, 9, 7, 7, 5 instead of 9, 9, 7, 7, 7) all certify.
#
# The probe runs a short block Lanczos process, from a fixed-seed start that
# never touches the global RNG stream, on the operator compressed to the
# orthogonal complement of the returned vectors V. Every Ritz value of that
# compression is a Rayleigh quotient of a unit vector orthogonal to V, so it
# lies inside the spectrum of the compression; when V spans an (approximately)
# invariant subspace the compression's spectrum is the rest of A's spectrum up
# to the residual norm. A complement Ritz value beyond the least-preferred
# returned value (by more than the residual and tolerance margin) is therefore
# proof that a more-preferred eigenvalue is missing. The converse does not
# hold: a short Krylov run can fail to resolve an intruder, so a clean probe is
# evidence, not proof. The deterministic answer is the inertia (LDL')
# eigenvalue-counting certificate in R/completeness_inertia.R, used instead
# of the probe whenever the operator has an explicit matrix source and the
# factorisation is affordable (completeness mode "auto"); the probe remains
# the check for matrix-free operators.
#
# When the probe finds an intruder the solve is repaired: the deflated
# operator P A P (P = I - V V') is solved with the native block Lanczos kernel
# started from the probe's intruding Ritz vectors, the complement pairs are
# merged with V by a Rayleigh-Ritz step on span(V, X_c), the merged set is
# re-certified from scratch and probed again, up to a small round budget. A
# result whose probe still finds an unresolved intruder reports
# target_completeness = "failed" and passed = FALSE.

#' @keywords internal
target_completeness_states <- function() {
  c("probed", "repaired", "failed", "not_checked", "exact",
    "inertia_verified", "inertia_failed", "inertia_inconclusive")
}

#' @keywords internal
validate_completeness_mode <- function(mode, arg = "completeness") {
  if (is.null(mode)) {
    return(NULL)
  }
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% completeness_modes()) {
    stop(arg, " must be one of \"auto\", \"inertia\", \"probe\" or \"none\".",
         call. = FALSE)
  }
  mode
}

# Resolve the completeness mode for a solve: an explicit method descriptor
# setting wins, then the eigencore.target_completeness option, then "auto".
#' @keywords internal
target_completeness_mode <- function(method = NULL) {
  mode <- if (inherits(method, "eigencore_method")) method$completeness else NULL
  mode <- mode %||% getOption("eigencore.target_completeness", "auto")
  validate_completeness_mode(mode, "option eigencore.target_completeness")
}

#' @keywords internal
target_completeness_controls <- function() {
  as_count <- function(value, default, lo) {
    value <- suppressWarnings(as.integer(value %||% default))
    if (length(value) != 1L || is.na(value) || value < lo) default else value
  }
  list(
    block = as_count(getOption("eigencore.completeness_probe_block"), 2L, 1L),
    steps = as_count(getOption("eigencore.completeness_probe_steps"), 8L, 1L),
    max_rounds = as_count(getOption("eigencore.completeness_max_rounds"), 3L, 0L)
  )
}

# Targets whose completeness a Krylov probe can witness: the algebraic and
# magnitude edges of a real spectrum.
#' @keywords internal
completeness_target_kind <- function(target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else NA_character_
  switch(
    kind,
    largest = ,
    largest_real = "largest",
    smallest = ,
    smallest_real = "smallest",
    largest_magnitude = "largest_magnitude",
    NULL
  )
}

#' @keywords internal
completeness_edge <- function(values, kind) {
  switch(
    kind,
    largest = min(values),
    smallest = max(values),
    largest_magnitude = min(abs(values))
  )
}

# Is a complement value more preferred than the returned edge by > margin?
#' @keywords internal
completeness_beyond <- function(theta, edge, margin, kind) {
  switch(
    kind,
    largest = theta > edge + margin,
    smallest = theta < edge - margin,
    largest_magnitude = abs(theta) > edge + margin
  )
}

#' @keywords internal
completeness_most_preferred <- function(theta, kind) {
  if (!length(theta)) {
    return(NA_real_)
  }
  switch(
    kind,
    largest = max(theta),
    smallest = min(theta),
    largest_magnitude = theta[[which.max(abs(theta))]]
  )
}

# Deterministic start columns from a fixed Mersenne-Twister stream. The
# caller's .Random.seed (and with it the RNG kind) is restored on exit, so the
# probe never consumes or perturbs the global random stream.
#' @keywords internal
completeness_probe_start <- function(n, cols, stream = 0L) {
  state <- saved_random_seed()
  on.exit(restore_random_seed(state), add = TRUE)
  set.seed(20261009L + as.integer(stream), kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
  matrix(stats::rnorm(n * cols), nrow = n, ncol = cols)
}

# Two-pass classical Gram-Schmidt of the columns of X against basis Bas, then
# column-by-column normalisation with a relative breakdown floor. Columns that
# fall below the floor are dropped.
#' @keywords internal
completeness_orthonormalize <- function(X, Bas = NULL, floor_abs = 0,
                                        Bas2 = NULL) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  x0 <- sqrt(colSums(X * X))
  # Native two-pass block classical Gram-Schmidt against each basis (zero
  # columns in a preallocated basis are harmless).
  for (basis in list(Bas, Bas2)) {
    if (!is.null(basis) && ncol(basis)) {
      X <- reorthogonalize_against(X, basis, passes = 2L)
    }
  }
  kept <- 0L
  for (j in seq_len(ncol(X))) {
    x <- X[, j]
    if (!is.finite(x0[[j]]) || x0[[j]] == 0) {
      next
    }
    if (kept) {
      for (pass in 1:2) {
        for (i in seq_len(kept)) {
          q <- X[, i]
          x <- x - sum(q * x) * q
        }
      }
    }
    nx <- sqrt(sum(x * x))
    if (is.finite(nx) && nx > max(1e-10 * x0[[j]], floor_abs)) {
      kept <- kept + 1L
      X[, kept] <- x / nx
    }
  }
  X[, seq_len(kept), drop = FALSE]
}

# Short block Lanczos (full reorthogonalisation; Rayleigh-Ritz on the
# explored basis) on A compressed to the orthogonal complement of V, run by
# the native kernel in src/completeness_probe.cpp: built-in dense and CSC
# operators apply natively, every other operator through its R apply closure.
# Stops early once an intruder is visible. Returns the complement Ritz values
# and vectors, the operator work spent and whether the complement was
# exhausted.
#' @keywords internal
completeness_probe <- function(op, values, V, kind, margin, block, steps,
                               stream = 0L) {
  n <- as.integer(op$dim[[1L]])
  V <- as.matrix(V)
  storage.mode(V) <- "double"
  finite_values <- abs(values[is.finite(values)])
  scale <- if (length(finite_values) && max(finite_values) > 0) max(finite_values) else 1
  params <- c(
    as.numeric(block), as.numeric(steps),
    switch(kind, largest = 1, smallest = 2, largest_magnitude = 3),
    completeness_edge(values, kind), as.numeric(margin), as.numeric(stream),
    scale
  )
  storage <- op$metadata$storage %||% NULL
  source <- source_or_null(op)
  if (is.matrix(source) && is.double(source) && !is.complex(source)) {
    .Call("eigencore_completeness_probe_dense", source, V, params,
          PACKAGE = "eigencore")
  } else if (identical(storage, "dgCMatrix")) {
    A <- op$metadata$matrix
    .Call("eigencore_completeness_probe_csc",
          methods::slot(A, "i"), methods::slot(A, "p"), methods::slot(A, "x"),
          methods::slot(A, "Dim"), V, params, PACKAGE = "eigencore")
  } else {
    apply_closure <- function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- as.matrix(apply_operator(op, X))
      storage.mode(out) <- "double"
      out <- alpha * out
      if (is.null(Y) || beta == 0) out else out + beta * Y
    }
    .Call("eigencore_completeness_probe_r_operator", as.integer(op$dim),
          apply_closure, V, params, PACKAGE = "eigencore")
  }
}

#' @keywords internal
completeness_margin <- function(values, residuals, theta, tol, scale = NULL) {
  finite_abs <- function(x) {
    x <- abs(as.numeric(x))
    x[is.finite(x)]
  }
  scale <- max(c(finite_abs(values), finite_abs(theta), finite_abs(scale), 0))
  residual_mass <- sqrt(sum(finite_abs(residuals)^2))
  max(tol * scale, 2 * residual_mass, 64 * .Machine$double.eps * scale)
}

# A Hermitian operator equal to P A P + shift * V V' (P = I - V V'): the
# complement problem with the returned directions moved to `shift`, on the
# non-preferred side of the target, so a complement solve does not return them.
#' @keywords internal
completeness_deflated_operator <- function(op, V, shift) {
  V <- as.matrix(V)
  linear_operator(
    dim = op$dim,
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      X <- as.matrix(X)
      C <- crossprod(V, X)
      PX <- X - V %*% C
      APX <- as.matrix(apply_operator(op, PX))
      Z <- APX - V %*% crossprod(V, APX)
      if (shift != 0) {
        Z <- Z + shift * (V %*% C)
      }
      Z <- alpha * Z
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    structure = hermitian(),
    name = "eigencore deflated complement operator"
  )
}

#' @keywords internal
completeness_shift <- function(values, kind, scale) {
  switch(
    kind,
    largest = min(values) - max(scale, 1e-300),
    smallest = max(values) + max(scale, 1e-300),
    largest_magnitude = 0
  )
}

# One repair round: solve the deflated complement problem from the probe's
# intruding Ritz vectors, then Rayleigh-Ritz on span(V, X_c).
#' @keywords internal
completeness_repair_round <- function(op, values, V, target, kind, tol, probe,
                                      margin) {
  n <- as.integer(op$dim[[1L]])
  k <- ncol(V)
  nc <- n - k
  edge <- completeness_edge(values, kind)
  intruding <- which(completeness_beyond(probe$theta, edge, margin, kind))
  scale <- max(abs(c(values, probe$theta)))
  if (isTRUE(probe$exhausted)) {
    # The probe spanned the entire complement: its Ritz pairs are the
    # complement's exact eigenpairs.
    Xc <- probe$Y
    columns <- 0L
    block_calls <- 0L
  } else {
    block <- min(max(2L, length(intruding) + 1L), 8L, nc)
    kc <- min(k, nc)
    m_max <- min(n, default_block_lanczos_max_subspace(kc, block))
    if (m_max < kc + block) {
      kc <- max(1L, m_max - block)
    }
    order_theta <- order_indices(probe$theta, target)
    start_cols <- probe$Y[, utils::head(order_theta, block), drop = FALSE]
    if (ncol(start_cols) < block) {
      extra <- completeness_probe_start(n, block - ncol(start_cols), 977L)
      start_cols <- cbind(start_cols, extra)
    }
    Gop <- completeness_deflated_operator(
      op, V, completeness_shift(values, kind, scale)
    )
    comp <- native_block_lanczos_hermitian(
      Gop, k = kc, target = target, tol = tol, maxit = m_max, block = block,
      max_restarts = 100L, vectors = TRUE, full_subspace = FALSE,
      certificate_fallback = FALSE, start = start_cols
    )
    Xc <- comp$vectors
    columns <- as.integer(comp$operator_columns %||% comp$matvecs %||% 0L)
    block_calls <- as.integer(comp$operator_block_calls %||% comp$matvecs %||% 0L)
  }
  Qc <- completeness_orthonormalize(Xc, V)
  Z <- cbind(V, Qc)
  AZ <- as.matrix(apply_operator(op, Z))
  H <- crossprod(Z, AZ)
  H <- (H + t(H)) / 2
  eig <- eigen(H, symmetric = TRUE)
  sel <- utils::head(order_indices(eig$values, target), k)
  list(
    values = eig$values[sel],
    vectors = Z %*% eig$vectors[, sel, drop = FALSE],
    columns = columns + ncol(Z),
    block_calls = block_calls + 1L
  )
}

# Probe (and if needed repair) a standard Hermitian partial eigensolution.
# Returns the possibly updated values/vectors, the completeness record and the
# operator work spent. `residuals` are the certificate's absolute residuals.
#' @keywords internal
target_completeness_check <- function(op, values, vectors, target, tol,
                                      residuals = NULL, norm_scale = NULL,
                                      controls = target_completeness_controls()) {
  kind <- completeness_target_kind(target)
  V <- as.matrix(vectors)
  record <- list(
    method = "deflated_complement_probe",
    block = controls$block,
    max_steps = controls$steps,
    steps = 0L,
    operator_columns = 0L,
    operator_block_calls = 0L,
    rounds = 0L,
    edge = completeness_edge(values, kind),
    most_preferred_complement = NA_real_,
    margin = NA_real_,
    intruder_found = FALSE,
    exhausted = FALSE
  )
  status <- "probed"
  repaired <- FALSE
  round <- 0L
  repeat {
    margin <- completeness_margin(values, residuals, numeric(), tol, norm_scale)
    probe <- completeness_probe(
      op, values, V, kind, margin,
      block = controls$block, steps = controls$steps, stream = round
    )
    margin <- completeness_margin(values, residuals, probe$theta, tol, norm_scale)
    record$steps <- record$steps + probe$steps
    record$operator_columns <- record$operator_columns + probe$columns
    record$operator_block_calls <- record$operator_block_calls + probe$block_calls
    record$edge <- completeness_edge(values, kind)
    record$margin <- margin
    record$most_preferred_complement <- completeness_most_preferred(probe$theta, kind)
    record$exhausted <- isTRUE(probe$exhausted)
    intruder <- length(probe$theta) &&
      any(completeness_beyond(probe$theta, record$edge, margin, kind))
    if (!intruder) {
      status <- if (repaired) "repaired" else "probed"
      break
    }
    record$intruder_found <- TRUE
    if (round >= controls$max_rounds) {
      status <- "failed"
      break
    }
    round <- round + 1L
    fixed <- completeness_repair_round(op, values, V, target, kind, tol, probe,
                                       margin)
    record$operator_columns <- record$operator_columns + fixed$columns
    record$operator_block_calls <- record$operator_block_calls + fixed$block_calls
    values <- fixed$values
    V <- fixed$vectors
    AV <- as.matrix(apply_operator(op, V))
    residuals <- sqrt(colSums((AV - sweep(V, 2L, values, `*`))^2))
    record$operator_columns <- record$operator_columns + ncol(V)
    record$operator_block_calls <- record$operator_block_calls + 1L
    repaired <- TRUE
  }
  record$rounds <- round
  list(
    status = status,
    record = record,
    repaired = repaired,
    values = values,
    vectors = V
  )
}

# Attach a completeness verdict to a certificate. `passed` requires the
# residual certificate and, when a probe ran, no unresolved intruder.
#' @keywords internal
certificate_with_completeness <- function(certificate, status, record = NULL) {
  if (is.null(certificate)) {
    return(certificate)
  }
  certificate$residual_passed <- certificate$residual_passed %||%
    isTRUE(certificate$passed)
  certificate$target_completeness <- status
  certificate$target_passed <- switch(
    status,
    probed = TRUE,
    repaired = TRUE,
    exact = TRUE,
    inertia_verified = TRUE,
    failed = FALSE,
    inertia_failed = FALSE,
    NA
  )
  certificate$completeness <- record
  if (identical(status, "failed")) {
    certificate$passed <- FALSE
    certificate$notes <- unique(c(
      certificate$notes,
      "target completeness probe found a more-preferred eigenvalue outside the returned set that could not be resolved"
    ))
  } else if (identical(status, "inertia_failed")) {
    certificate$passed <- FALSE
    certificate$notes <- unique(c(
      certificate$notes,
      "inertia count proves a more-preferred eigenvalue is missing from the returned set"
    ))
  } else if (identical(status, "repaired") ||
             (identical(status, "inertia_verified") && isTRUE(record$repaired))) {
    certificate$notes <- unique(c(
      certificate$notes,
      "target completeness check found a missing eigenvalue copy; the returned set was repaired by a deflated complement solve"
    ))
  } else if (identical(status, "inertia_inconclusive")) {
    certificate$notes <- unique(c(
      certificate$notes,
      paste0("inertia completeness check inconclusive: ",
             record$reason %||% "no reliable separation")
    ))
  }
  certificate
}

# Route classes for the completeness verdict: full-spectrum routes are
# complete by construction; Krylov / LOBPCG / shift-invert routes are probed.
#' @keywords internal
target_completeness_route_class <- function(plan, problem) {
  method <- plan$method %||% ""
  if (is_transform_method(problem$transform)) {
    if (plan$method %in% shift_invert_arnoldi_labels()) {
      return("unsupported")
    }
    return("krylov")
  }
  if (plan_dispatches_lobpcg(plan)) {
    # Constrained LOBPCG targets the spectrum on the constraint complement;
    # the probe would (correctly but unhelpfully) find the deflated values.
    if (!is.null(plan$method_descriptor$constraints)) {
      return("unsupported")
    }
    return("krylov")
  }
  if (plan_dispatches_structured_grid_laplacian_2d(plan) ||
      plan_dispatches_native_tridiagonal_eigen(plan)) {
    return("exact")
  }
  if (plan_dispatches_lanczos(plan)) {
    return("krylov")
  }
  if (plan_dispatches_sparse_general_pencil_arnoldi(plan) ||
      plan_dispatches_arnoldi(plan)) {
    return("unsupported")
  }
  "exact"
}

# Is a standard (B = NULL) real Hermitian result eligible for the probe?
#' @keywords internal
target_completeness_eligible <- function(problem, k) {
  op <- problem$A
  is.null(problem$metric) &&
    inherits(op, "eigencore_operator") &&
    identical(op$structure$kind, "hermitian") &&
    identical(op$dtype %||% "double", "double") &&
    !is.null(completeness_target_kind(problem$target)) &&
    k < op$dim[[1L]]
}

# Post-dispatch hook for execute_eigen_plan(): probe eligible Krylov results,
# label full-spectrum results, and leave others "not_checked".
# Is a real Hermitian result (standard, or generalized with SPD B) eligible
# for the inertia certificate? The matrix-source and cost checks are in
# inertia_completeness_gate().
#' @keywords internal
inertia_completeness_eligible <- function(problem, k) {
  op <- problem$A
  inherits(op, "eigencore_operator") &&
    identical(op$structure$kind, "hermitian") &&
    identical(op$dtype %||% "double", "double") &&
    !is.null(inertia_completeness_kind(problem$target)) &&
    k < op$dim[[1L]] &&
    (is.null(problem$metric) || isTRUE(tryCatch(
      generalized_spd_metric_known(problem$metric), error = function(e) FALSE)))
}

#' @keywords internal
apply_target_completeness <- function(result, plan, problem, k, mode,
                                      vectors_requested, solve_seconds = NA_real_,
                                      seed = NULL) {
  cert <- result$certificate
  if (is.null(cert) || !is.null(cert$target_completeness)) {
    if (!isTRUE(vectors_requested)) {
      result["vectors"] <- list(NULL)
    }
    return(result)
  }
  route <- target_completeness_route_class(plan, problem)
  certify <- isTRUE(plan$execution$certify)
  krylov_ok <- identical(route, "krylov") && certify &&
    isTRUE(cert$passed) &&
    is.numeric(result$values) && !is.complex(result$values) &&
    length(result$values) == k
  gate <- NULL
  if (!identical(route, "exact") && mode %in% c("auto", "inertia") && krylov_ok &&
      inertia_completeness_eligible(problem, k)) {
    gate <- inertia_completeness_gate(problem, mode, k, solve_seconds = solve_seconds,
                                      seed = seed)
    if (isTRUE(gate$use)) {
      check <- inertia_completeness_run(
        problem, result$values, result$vectors, cert,
        tol = cert$tolerance %||% plan$execution$tol,
        ctx = gate$ctx, gate = gate
      )
      result <- result_with_completeness(result, check, problem, plan)
      if (!isTRUE(vectors_requested)) {
        result["vectors"] <- list(NULL)
      }
      return(result)
    }
  }
  if (identical(route, "exact")) {
    result$certificate <- certificate_with_completeness(cert, "exact")
  } else if (mode %in% c("auto", "inertia", "probe") && krylov_ok &&
             !is.null(result$vectors) &&
             target_completeness_eligible(problem, k)) {
    started <- proc.time()[["elapsed"]]
    check <- target_completeness_check(
      problem$A, result$values, result$vectors, problem$target,
      tol = cert$tolerance %||% plan$execution$tol,
      residuals = cert$residuals,
      norm_scale = cert$scale
    )
    check$record$seconds <- proc.time()[["elapsed"]] - started
    if (!is.null(gate)) {
      check$record$inertia_gate <- gate$reason
      check$record$inertia_predicted_seconds <- gate$predicted_seconds
    }
    result <- result_with_completeness(result, check, problem, plan)
  } else {
    result$certificate <- certificate_with_completeness(
      cert, "not_checked",
      if (is.null(gate)) NULL else list(inertia_gate = gate$reason,
                                        inertia_predicted_seconds = gate$predicted_seconds)
    )
  }
  if (!isTRUE(vectors_requested)) {
    result["vectors"] <- list(NULL)
  }
  result
}

#' @keywords internal
result_with_completeness <- function(result, check, problem, plan) {
  record <- check$record
  cert <- result$certificate
  extra_cert_columns <- 0L
  if (isTRUE(check$repaired)) {
    tol <- cert$tolerance %||% plan$execution$tol
    new_cert <- certify_eigen_operator(problem$A, check$values, check$vectors,
                                       tol = tol)
    new_cert$notes <- unique(c(cert$notes, new_cert$notes))
    cert <- new_cert
    extra_cert_columns <- length(check$values)
    result$values <- check$values
    result$vectors <- check$vectors
    result$residuals <- cert$residuals
    result$backward_error <- cert$backward_error
    result$orthogonality <- cert$orthogonality
    result$nconv <- sum(cert$converged)
    result$locked <- which(cert$converged)
  }
  cert <- certificate_with_completeness(cert, check$status, record)
  result$certificate <- cert
  if (identical(check$status, "failed")) {
    result$warnings <- c(
      result$warnings,
      paste0(
        "target completeness probe found an eigenvalue (",
        format(record$most_preferred_complement, digits = 8),
        ") more preferred than the returned edge (",
        format(record$edge, digits = 8),
        ") outside the returned set; certificate withheld"
      )
    )
  } else if (identical(check$status, "inertia_failed")) {
    result$warnings <- c(
      result$warnings,
      paste0(
        "inertia count proves an eigenvalue more preferred than the returned edge is ",
        "missing (", record$reason %||% "count mismatch", "); certificate withheld"
      )
    )
  } else if (identical(check$status, "repaired") ||
             (identical(check$status, "inertia_verified") && isTRUE(record$repaired))) {
    result$warnings <- c(
      result$warnings,
      "target completeness probe found a missing eigenvalue copy; result repaired by a deflated complement solve"
    )
  }
  record$operator_block_calls <- record$operator_block_calls %||% 0L
  record$operator_columns <- record$operator_columns %||% 0L
  result$operator_block_calls <- as.integer(
    (result$operator_block_calls %||% 0L) + record$operator_block_calls +
      (extra_cert_columns > 0L)
  )
  result$operator_columns <- as.integer(
    (result$operator_columns %||% 0L) + record$operator_columns + extra_cert_columns
  )
  if (extra_cert_columns > 0L) {
    result$certification_operator_columns <- as.integer(
      (result$certification_operator_columns %||% 0L) + extra_cert_columns
    )
  }
  if (inherits(result$work, "eigencore_work")) {
    work <- unclass(result$work)
    add <- function(field, amount) {
      if (!is.na(work[[field]] %||% NA)) {
        work[[field]] <<- as.integer(work[[field]] + amount)
      }
    }
    add("operator_block_calls", record$operator_block_calls)
    add("operator_columns", record$operator_columns)
    if (extra_cert_columns > 0L) {
      add("certification_operator_block_calls", 1L)
      add("certification_operator_columns", extra_cert_columns)
    }
    result$work <- new_typed_work_record(work)
  }
  result
}
