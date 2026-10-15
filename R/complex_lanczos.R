# Native complex Hermitian thick-restart Lanczos for matrix-free complex
# operators.
#
# The realified route (R/complex_hermitian.R) solves H = A + iB as the real
# symmetric 2n embedding for 2k pairs. Every complex eigendirection is then
# an invariant real plane, found twice, and the real block Lanczos needs
# about twice (in practice 2-4x, block 2 and the restarts included) the
# operator columns of a scalar Lanczos run on H itself. For explicit sources
# that is cheap enough (BLAS-3 kernels on the embedding, CHOLMOD LDL'); for
# matrix-free complex operators every real column is one application of the
# caller's complex function, which then dominates. Those problems run this
# scalar Lanczos in complex arithmetic instead.
#
# Method: Krylov-Schur form of thick-restart Lanczos (Wu & Simon 2000;
# Stewart 2001) for a Hermitian operator. A cycle extends
#   H V_j = V_j T_j + beta_j v_{j+1} e_j^H
# one vector at a time with full two-pass classical Gram-Schmidt against the
# basis (the projections fill T, so the arrowhead coupling after a restart
# needs no special case); the Rayleigh-Ritz pairs (theta_i, V_j y_i) have
# residual norms beta_j |e_j^H y_i|. At the end of a cycle the p preferred
# Ritz vectors are kept, T becomes diag(theta_keep), and v_{j+1} continues
# the basis. Convergence is declared when the k wanted pairs have estimated
# backward error residual / (||H||_est + |theta|) below tol / 2, with
# ||H||_est = max |theta| (a lower bound of ||H||_2, so the certificate's
# backward error, computed from an upper norm bound, is smaller still).
#
# Scope: standard problems with a matrix-free complex operator (no explicit
# matrix or real/imaginary split), targets largest(), smallest() and
# largest_magnitude(), method auto() or lanczos(block = 1). The result is
# certified on the complex operator and its completeness checked like every
# other complex result (realified probe / repair,
# complex_hermitian_target_completeness()). A scalar Krylov method sees one
# copy of a repeated eigenvalue per start vector (exact arithmetic); a missed
# copy is found by that completeness check, as for the real scalar routes.
# getOption("eigencore.complex_native_lanczos", TRUE) = FALSE restores the
# realified route.

#' @keywords internal
complex_native_lanczos_label <- function() {
  "native complex Hermitian thick-restart Lanczos (matrix-free)"
}

#' @keywords internal
complex_native_lanczos_enabled <- function() {
  !isFALSE(getOption("eigencore.complex_native_lanczos", TRUE))
}

#' @keywords internal
complex_native_lanczos_targets <- function() {
  c("largest", "smallest", "largest_magnitude")
}

# Does a planned complex problem run the native complex Lanczos?
#' @keywords internal
complex_native_lanczos_eligible <- function(problem, method, initial_subspace = NULL) {
  if (!complex_native_lanczos_enabled() || !is.null(problem$metric) ||
      !is.null(problem$transform)) {
    return(FALSE)
  }
  op <- problem$A
  if (!identical(op$dtype %||% "double", "complex") ||
      !is.null(source_or_null(op)) || !is.null(op$metadata$matrix %||% NULL) ||
      !is.null(op$metadata$complex_split %||% NULL)) {
    return(FALSE)
  }
  target_kind <- if (inherits(problem$target, "eigencore_target")) problem$target$kind else ""
  if (!target_kind %in% complex_native_lanczos_targets()) {
    return(FALSE)
  }
  kind <- if (inherits(method, "eigencore_method")) method$kind else "auto"
  if (identical(kind, "auto")) {
    return(TRUE)
  }
  identical(kind, "lanczos") && identical(as.integer(method$block %||% 1L), 1L)
}

# Preference order of Ritz values for a target kind.
#' @keywords internal
complex_lanczos_order <- function(theta, kind) {
  switch(
    kind,
    largest = order(theta, decreasing = TRUE),
    smallest = order(theta),
    largest_magnitude = order(abs(theta), decreasing = TRUE)
  )
}

# Thick-restart Lanczos core. `apply_fn(z)` returns H z for an n x 1
# complex matrix z. Returns values (target order), vectors, residual
# estimates, convergence flag and work counters.
#' @keywords internal
complex_hermitian_trlan <- function(apply_fn, n, k, kind, tol, m = NULL,
                                    max_restarts = 100L, start = NULL) {
  n <- as.integer(n)
  k <- as.integer(k)
  if (is.null(m)) {
    m <- max(k + 1L, 3L * k + 20L, 5L * k)
  }
  m <- as.integer(min(n, m))
  if (m < k + 1L && m < n) {
    stop_lanczos_subspace_exhausted(k, n, 1L, m)
  }
  keep_count <- function() {
    as.integer(min(m - 1L, max(k, k + (m - k) %/% 2L)))
  }
  rand_vec <- function() {
    complex(real = stats::rnorm(n), imaginary = stats::rnorm(n))
  }
  V <- matrix(0i, n, m + 1L)
  T <- matrix(0i, m, m)
  v <- if (is.null(start)) rand_vec() else as.complex(as.matrix(start)[, 1L])
  if (length(v) != n || !all(is.finite(Re(v))) || !all(is.finite(Im(v)))) {
    stop("start vector must be finite of length ", n, ".", call. = FALSE)
  }
  nv <- sqrt(sum(Mod(v)^2))
  if (!(nv > 0)) {
    v <- rand_vec()
    nv <- sqrt(sum(Mod(v)^2))
  }
  V[, 1L] <- v / nv
  p <- 0L
  matvecs <- 0L
  restarts <- 0L
  norm_est <- 0
  eta <- 1 / sqrt(2)
  result <- NULL
  check_every <- 5L
  repeat {
    j <- p
    beta <- 0
    while (j < m) {
      j <- j + 1L
      w <- as.vector(apply_fn(V[, j, drop = FALSE]))
      matvecs <- matvecs + 1L
      if (!all(is.finite(Re(w))) || !all(is.finite(Im(w)))) {
        stop("complex Lanczos: the operator returned non-finite values.", call. = FALSE)
      }
      Vj <- V[, seq_len(j), drop = FALSE]
      pre <- sqrt(sum(Mod(w)^2))
      h <- as.vector(crossprod(Conj(Vj), w))
      w <- w - as.vector(Vj %*% h)
      if (sqrt(sum(Mod(w)^2)) < eta * pre) {
        h2 <- as.vector(crossprod(Conj(Vj), w))
        w <- w - as.vector(Vj %*% h2)
        h <- h + h2
      }
      # Hermitian projection: real diagonal; the column fills the upper
      # triangle (kept vectors' couplings included), the row its conjugate.
      h[[j]] <- Re(h[[j]])
      T[seq_len(j), j] <- h
      T[j, seq_len(j)] <- Conj(h)
      norm_est <- max(norm_est, abs(Re(h[[j]])))
      beta <- sqrt(sum(Mod(w)^2))
      if (!(beta > 100 * .Machine$double.eps * max(norm_est, .Machine$double.xmin))) {
        # Invariant subspace: continue with a random direction orthogonal
        # to the basis (zero coupling), or stop when the space is full.
        beta <- 0
        if (j >= n) {
          V[, j + 1L] <- 0
        } else {
          w <- rand_vec()
          for (pass in 1:2) {
            w <- w - as.vector(Vj %*% crossprod(Conj(Vj), w))
          }
          V[, j + 1L] <- w / sqrt(sum(Mod(w)^2))
        }
      } else {
        V[, j + 1L] <- w / beta
      }
      if (j >= k && (j == m || j >= n || beta == 0 || (j - p) %% check_every == 0L)) {
        e <- eigen(T[seq_len(j), seq_len(j), drop = FALSE], symmetric = TRUE)
        theta <- e$values
        norm_est <- max(norm_est, abs(theta))
        res <- beta * Mod(e$vectors[j, ])
        ord <- complex_lanczos_order(theta, kind)
        want <- ord[seq_len(min(k, j))]
        threshold <- 0.5 * tol * (norm_est + abs(theta[want]))
        converged <- res[want] <= threshold
        if ((length(want) == k && all(converged)) || j >= n) {
          result <- list(theta = theta, Y = e$vectors, res = res, j = j, want = want,
                         converged = all(converged))
          break
        }
        if (j == m) {
          result <- list(theta = theta, Y = e$vectors, res = res, j = j, want = want,
                         converged = FALSE, ord = ord)
        }
      }
      if (j >= n) {
        break
      }
    }
    if (isTRUE(result$converged) || result$j >= n || restarts >= max_restarts) {
      break
    }
    # Thick restart: keep the p preferred Ritz vectors.
    restarts <- restarts + 1L
    p <- keep_count()
    keep <- result$ord[seq_len(p)]
    Yk <- result$Y[, keep, drop = FALSE]
    Vnew <- V[, seq_len(m), drop = FALSE] %*% Yk
    vnext <- V[, m + 1L]
    V[] <- 0i
    V[, seq_len(p)] <- Vnew
    V[, p + 1L] <- vnext
    T[] <- 0i
    # The kept Ritz vectors' residuals all lie along v_{p+1}; their
    # couplings to it are filled by the projection of the next column.
    diag(T)[seq_len(p)] <- result$theta[keep]
    result <- NULL
  }
  j <- result$j
  sel <- result$want
  list(
    values = result$theta[sel],
    vectors = V[, seq_len(j), drop = FALSE] %*% result$Y[, sel, drop = FALSE],
    residual_estimates = result$res[sel],
    converged = isTRUE(result$converged),
    matvecs = matvecs,
    restarts = restarts,
    subspace = m
  )
}

# Planner hook (from complex_hermitian_plan()).
#' @keywords internal
complex_native_lanczos_plan <- function(problem, k, method, method_descriptor,
                                        execution, planner_policy, maxit) {
  n <- as.integer(problem$A$dim[[1L]])
  k <- as.integer(k)
  subspace <- if (inherits(method, "eigencore_method")) method$max_subspace else NULL
  subspace <- subspace %||% maxit %||% max(k + 1L, 3L * k + 20L, 5L * k)
  subspace <- as.integer(min(n, subspace))
  max_restarts <- if (inherits(method, "eigencore_method")) method$max_restarts else NULL
  controls <- list(
    embedding = "none",
    subspace = subspace,
    max_restarts = as.integer(max_restarts %||% 100L),
    iteration_limit = subspace,
    iteration_limit_kind = "subspace"
  )
  new_plan(
    problem,
    k = k,
    method = complex_native_lanczos_label(),
    method_descriptor = method_descriptor,
    reasons = c(
      paste0("structure: ", problem$structure$kind, " (complex)"),
      paste0("target: ", target_label(problem$target)),
      "standard eigenproblem",
      paste0("matrix-free complex Hermitian operator: scalar thick-restart Lanczos ",
             "in complex arithmetic (subspace ", subspace, "), not the real 2n ",
             "embedding, which needs about twice the operator columns"),
      operator_kernel_reason(problem$A)
    ),
    fallback = "realified complex Hermitian route (options(eigencore.complex_native_lanczos = FALSE))",
    controls = controls,
    execution = execution,
    planner_policy = planner_policy
  )
}

#' @keywords internal
solve_complex_hermitian_native <- function(plan) {
  started <- proc.time()[["elapsed"]]
  problem <- plan$problem
  execution <- plan$execution
  tol <- execution$tol
  k <- as.integer(plan$requested)
  op <- problem$A
  n <- as.integer(op$dim[[1L]])
  start <- execution$initial_subspace %||% NULL
  if (!is.null(start)) {
    start <- as.matrix(start)
    if (nrow(start) != n) {
      stop("initial_subspace must have ", n, " rows.", call. = FALSE)
    }
    start <- start[, 1L, drop = FALSE] + 0i
  }
  run <- complex_hermitian_trlan(
    function(z) as.matrix(apply_operator(op, z)),
    n = n, k = k, kind = problem$target$kind, tol = tol,
    m = plan$controls$subspace, max_restarts = plan$controls$max_restarts,
    start = start
  )
  values <- run$values
  vectors <- run$vectors
  warnings <- character()
  if (!isTRUE(run$converged)) {
    warnings <- c(warnings, sprintf(paste0(
      "complex Lanczos stopped after %d restarts before all %d pairs met the ",
      "tolerance (max residual estimate %.3g)"),
      run$restarts, k, max(run$residual_estimates)))
  }
  cert <- if (isTRUE(execution$certify)) {
    complex_hermitian_certify(problem, values, vectors, tol)
  } else {
    NULL
  }
  polished <- complex_hermitian_polish(problem, values, vectors, cert, tol)
  if (!is.null(polished)) {
    values <- polished$values
    vectors <- polished$vectors
    cert <- polished$certificate
  }
  if (!is.null(cert)) {
    cert$notes <- unique(c(cert$notes, paste0(
      "complex Hermitian matrix-free operator solved by scalar thick-restart ",
      "Lanczos in complex arithmetic; pairs certified on the complex operator")))
  }
  iter <- list(
    iterations = run$restarts + 1L,
    matvecs = run$matvecs,
    operator_block_calls = as.integer(run$matvecs + 1L),
    operator_columns = as.integer(run$matvecs + length(values)),
    certification_operator_columns = as.integer(length(values))
  )
  result <- make_eigen_result(
    values = values,
    vectors = vectors,
    certificate = cert %||% empty_certificate(tol, "certification disabled"),
    iter = iter,
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(problem$target),
    plan = plan,
    warnings = warnings,
    extras = list(
      operator_block_calls = iter$operator_block_calls,
      operator_columns = iter$operator_columns,
      certification_operator_columns = iter$certification_operator_columns,
      restart = list(
        kind = "thick_restart_complex",
        subspace = run$subspace,
        restarts = run$restarts,
        converged = run$converged,
        polished = if (is.null(polished)) NA_character_ else polished$method,
        seconds = proc.time()[["elapsed"]] - started
      ),
      locked = if (is.null(cert)) integer() else which(cert$converged)
    )
  )
  if (is.null(cert)) {
    result$certificate <- NULL
  }
  if (isTRUE(cert$passed) && length(values) == k &&
      !identical(target_completeness_mode(plan$method_descriptor), "none")) {
    checked <- complex_hermitian_target_completeness(
      result, plan, problem, k, target_completeness_mode(plan$method_descriptor), TRUE,
      solve_seconds = proc.time()[["elapsed"]] - started
    )
    if (!is.null(checked)) {
      result <- checked
    }
  }
  if (!is.null(result$certificate) && is.null(result$certificate$target_completeness)) {
    result$certificate <- certificate_with_completeness(
      result$certificate, "not_checked",
      list(method = "none", reason = "certificate not passed or completeness disabled"))
  }
  if (!isTRUE(execution$vectors)) {
    result["vectors"] <- list(NULL)
  }
  result
}
