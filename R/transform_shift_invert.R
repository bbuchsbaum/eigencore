#' @keywords internal
native_dense_shift_invert_label <- function() {
  "native dense Hermitian shift-invert (factorized Lanczos)"
}

#' @keywords internal
native_dense_generalized_shift_invert_label <- function() {
  "native dense generalized SPD shift-invert (factorized Lanczos)"
}

#' @keywords internal
native_tridiagonal_shift_invert_label <- function() {
  "native tridiagonal Hermitian shift-invert (factorized Lanczos)"
}

#' @keywords internal
native_tridiagonal_generalized_shift_invert_label <- function() {
  "native tridiagonal generalized SPD shift-invert (factorized Lanczos)"
}

#' @keywords internal
#' Nonsymmetric shift-invert: native Krylov-Schur Arnoldi on the factorised
#' (A - sigma I)^{-1} through the matrix-free callback entry.
shift_invert_arnoldi_label <- function(kind = c("dense_qr", "sparse_lu", "user_solve")) {
  kind <- match.arg(kind)
  suffix <- switch(
    kind,
    dense_qr = "dense QR solve callback",
    sparse_lu = "sparse LU solve callback",
    user_solve = "user solve callback"
  )
  paste0("native Krylov-Schur Arnoldi shift-invert (", suffix, ")")
}

#' @keywords internal
shift_invert_arnoldi_labels <- function() {
  vapply(c("dense_qr", "sparse_lu", "user_solve"), shift_invert_arnoldi_label,
         character(1L), USE.NAMES = FALSE)
}

#' @keywords internal
shift_invert_arnoldi_plan_controls <- function(problem, k) {
  n <- as.integer(problem$A$dim[1L])
  k <- as.integer(k)
  requested_subspace <- method_requested_max_subspace(problem$transform, problem)
  list(
    max_subspace = if (is.null(requested_subspace)) {
      native_krylov_schur_default_ncv(n, k)
    } else {
      min(n, requested_subspace)
    },
    max_restarts = 2L,
    krylov_schur_max_iterations = native_krylov_schur_default_maxit(),
    transform = "shift_invert",
    transformed_operator_target = "largest_magnitude",
    eigenvalue_recovery = "lambda = sigma + 1 / theta",
    certified_in_original_coordinates = TRUE,
    certification_policy = "right (and, when computed, left) residual certificate on the original nonsymmetric eigenproblem"
  )
}

#' @keywords internal
shift_invert_tridiagonal_parts <- function(A, shift = 0) {
  if (!inherits(A, "Matrix")) {
    return(NULL)
  }
  n <- nrow(A)
  if (is.null(n) || ncol(A) != n || n < 1L) {
    return(NULL)
  }
  shift <- as.numeric(shift)
  if (length(shift) != 1L || !is.finite(shift)) {
    return(NULL)
  }

  diag <- numeric(n)
  lower <- numeric(max(n - 1L, 0L))
  upper <- numeric(max(n - 1L, 0L))
  if (inherits(A, "diagonalMatrix")) {
    diag <- as.numeric(Matrix::diag(A))
  } else {
    if (!inherits(A, "CsparseMatrix")) {
      A <- tryCatch(methods::as(A, "CsparseMatrix"), error = function(e) NULL)
      if (is.null(A)) {
        return(NULL)
      }
    }
    row_of <- methods::slot(A, "i") + 1L
    p_slot <- methods::slot(A, "p")
    x_slot <- methods::slot(A, "x")
    if (length(x_slot) && any(!is.finite(x_slot))) {
      return(NULL)
    }
    col_of <- rep.int(seq_len(n), diff(p_slot))
    band <- row_of - col_of
    if (any(band > 1L | band < -1L)) {
      return(NULL)
    }
    # Valid CsparseMatrix objects carry unique sorted (row, column) entries,
    # so plain indexed assignment is the loop's accumulation.
    on_diag <- band == 0L
    on_lower <- band == 1L
    on_upper <- band == -1L
    diag[col_of[on_diag]] <- x_slot[on_diag]
    lower[col_of[on_lower]] <- x_slot[on_lower]
    upper[row_of[on_upper]] <- x_slot[on_upper]
  }
  if (any(!is.finite(diag)) || any(!is.finite(lower)) || any(!is.finite(upper))) {
    return(NULL)
  }
  if (n > 1L && !isTRUE(all.equal(lower, upper, tolerance = 1e-12))) {
    return(NULL)
  }
  list(lower = lower, diag = diag + shift, upper = upper)
}

#' @keywords internal
shift_invert_problem_tridiagonal_parts <- function(problem) {
  cached <- problem$A$metadata$shift_invert_tridiagonal_parts %||% NULL
  if (!is.null(cached)) {
    return(cached)
  }
  A <- problem$A$metadata$matrix %||% source_or_null(problem$A)
  if (!(inherits(A, "CsparseMatrix") || inherits(A, "diagonalMatrix"))) {
    return(NULL)
  }
  shift_invert_tridiagonal_parts(A, shift = 0)
}

#' Parse the sparse tridiagonal source once for one planned solve. The cached
#' vectors are immutable solver metadata: shifted diagonals are derived from
#' them, never written back into the operator.
#' @keywords internal
shift_invert_prepare_tridiagonal <- function(problem) {
  if (!is.null(problem$A$metadata$shift_invert_tridiagonal_parts)) {
    return(problem)
  }
  parts <- shift_invert_problem_tridiagonal_parts(problem)
  if (!is.null(parts)) {
    problem$A$metadata$shift_invert_tridiagonal_parts <- parts
  }
  problem
}

#' @keywords internal
shift_invert_is_native_tridiagonal <- function(problem) {
  A <- problem$A$metadata$matrix %||% source_or_null(problem$A)
  (inherits(A, "CsparseMatrix") || inherits(A, "diagonalMatrix")) &&
    !is.null(shift_invert_problem_tridiagonal_parts(problem))
}

#' @keywords internal
shift_invert_diagonal_metric_values <- function(Bop) {
  if (is.null(Bop)) {
    return(NULL)
  }
  Bstorage <- Bop$metadata$storage %||% NULL
  if (!identical(Bstorage, "ddiMatrix")) {
    return(NULL)
  }
  B <- Bop$metadata$matrix
  values <- if (identical(methods::slot(B, "diag"), "U")) {
    rep(1, Bop$dim[1L])
  } else {
    methods::slot(B, "x")
  }
  values <- as.numeric(values)
  if (length(values) != Bop$dim[1L] ||
      any(!is.finite(values)) ||
      any(values <= 0)) {
    return(NULL)
  }
  values
}

#' @keywords internal
shift_invert_plan_label <- function(problem, has_metric, is_hermitian,
                                    is_dense_source, is_native_csc) {
  user_solve <- problem$transform$solve
  if (!is_hermitian) {
    if (has_metric) {
      return("shift-invert requested (nonsymmetric generalized shift-invert is not implemented)")
    }
    if (!is.null(user_solve)) {
      return(shift_invert_arnoldi_label("user_solve"))
    }
    if (is_dense_source) {
      return(shift_invert_arnoldi_label("dense_qr"))
    }
    if (is_native_csc || inherits(problem$A$metadata$matrix, "CsparseMatrix")) {
      return(shift_invert_arnoldi_label("sparse_lu"))
    }
    return("shift-invert requested (provide method$solve for matrix-free A)")
  }
  # Every non-native-kernel route below runs the native thick-restart Lanczos
  # kernel on the factorised (A - sigma B)^{-1} through the matrix-free
  # callback ABI (shift_invert_transformed_lanczos), so the labels say so (C35).
  if (has_metric) {
    Bstorage <- problem$metric$metadata$storage %||% NULL
    Bsource <- source_or_null(problem$metric)
    dense_metric <- is.matrix(Bsource) && is.double(Bsource)
    diagonal_metric <- identical(Bstorage, "ddiMatrix")
    if (!is.null(user_solve)) {
      return("native thick-restart generalized SPD Lanczos shift-invert (user solve callback)")
    }
    if (is_dense_source && dense_metric) {
      return(native_dense_generalized_shift_invert_label())
    }
    if (is_dense_source && diagonal_metric) {
      return("native thick-restart generalized SPD Lanczos shift-invert (dense QR solve callback)")
    }
    if (diagonal_metric && shift_invert_is_native_tridiagonal(problem)) {
      return(native_tridiagonal_generalized_shift_invert_label())
    }
    csc_available <- inherits(problem$A$metadata$matrix, "CsparseMatrix")
    if ((is_native_csc || csc_available) && diagonal_metric) {
      return("native thick-restart generalized SPD Lanczos shift-invert (sparse LU solve callback)")
    }
    return("shift-invert requested (generalized SPD shift-invert requires dense A/B or sparse A with diagonal B)")
  }
  if (!is.null(user_solve)) {
    return("native thick-restart Hermitian Lanczos shift-invert (user solve callback)")
  }
  csc_available <- inherits(problem$A$metadata$matrix, "CsparseMatrix")
  diagonal_available <- inherits(problem$A$metadata$matrix, "diagonalMatrix")
  if (is_native_csc || csc_available || diagonal_available) {
    if (shift_invert_is_native_tridiagonal(problem)) {
      return(native_tridiagonal_shift_invert_label())
    }
    return("native thick-restart Hermitian Lanczos shift-invert (sparse LU solve callback)")
  }
  if (is_dense_source) {
    return(native_dense_shift_invert_label())
  }
  "shift-invert requested (provide method$solve for matrix-free A)"
}

#' @keywords internal
shift_invert_apply_factory <- function(solve_fn) {
  function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- solve_fn(X)
    if (is.null(Y) || beta == 0) {
      if (alpha == 1) Z else alpha * Z
    } else {
      alpha * Z + beta * Y
    }
  }
}

#' @keywords internal
#' Small non-cryptographic digest of a numeric vector. Captures length,
#' running sums, range, and a deterministic head/tail sample so two
#' mathematically identical matrices map to identical fingerprints. Two
#' distinct matrices with all summary statistics matching would collide,
#' which is acceptable for shift-invert cache keys (collisions only matter
#' if a user tries to reuse a factorization across truly distinct A).
shift_invert_double_digest <- function(values) {
  values <- as.numeric(values)
  n <- length(values)
  if (n == 0L) {
    return("empty|0|0|0|0|0|")
  }
  finite <- is.finite(values)
  fin_values <- values[finite]
  head_n <- min(8L, n)
  tail_n <- min(8L, n)
  paste(
    n,
    format(sum(fin_values), digits = 17),
    format(sum(fin_values * fin_values), digits = 17),
    format(if (length(fin_values)) min(fin_values) else NA_real_, digits = 17),
    format(if (length(fin_values)) max(fin_values) else NA_real_, digits = 17),
    paste(format(values[seq_len(head_n)], digits = 17), collapse = "_"),
    paste(format(values[seq.int(n - tail_n + 1L, n)], digits = 17), collapse = "_"),
    sum(!finite),
    sep = "|"
  )
}

#' @keywords internal
shift_invert_operator_fingerprint <- function(op) {
  op <- as_operator(op)
  source <- source_or_null(op)
  storage <- op$metadata$storage %||% NULL
  matrix <- op$metadata$matrix %||% NULL
  if (is.matrix(source)) {
    return(list(
      kind = "dense",
      dim = dim(source),
      storage_mode = storage.mode(source),
      digest = shift_invert_double_digest(as.numeric(source))
    ))
  }
  if (inherits(matrix, "sparseMatrix")) {
    if (!inherits(matrix, "CsparseMatrix")) {
      matrix <- methods::as(matrix, "CsparseMatrix")
    }
    i_slot <- methods::slot(matrix, "i")
    p_slot <- methods::slot(matrix, "p")
    x_slot <- if ("x" %in% methods::slotNames(matrix)) {
      methods::slot(matrix, "x")
    } else {
      numeric(0)
    }
    return(list(
      kind = storage %||% class(matrix)[[1L]],
      class = class(matrix),
      dim = methods::slot(matrix, "Dim"),
      nnz = length(x_slot),
      i_digest = shift_invert_double_digest(as.numeric(i_slot)),
      p_digest = shift_invert_double_digest(as.numeric(p_slot)),
      x_digest = shift_invert_double_digest(x_slot)
    ))
  }
  list(
    kind = "operator",
    dim = op$dim,
    name = op$name %||% NA_character_,
    storage = storage,
    source_available = !is.null(source)
  )
}

#' @keywords internal
shift_invert_factorization_cache_key <- function(Aop, sigma, Bop = NULL) {
  Aop <- as_operator(Aop)
  sigma <- as.numeric(sigma)
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("shift-invert cache key requires a single finite sigma.", call. = FALSE)
  }
  key <- list(
    transform = "shift_invert",
    sigma = sigma,
    structure = Aop$structure$kind %||% NA_character_,
    A = shift_invert_operator_fingerprint(Aop),
    B = if (is.null(Bop)) NULL else shift_invert_operator_fingerprint(Bop),
    standard_problem = is.null(Bop)
  )
  class(key) <- "eigencore_shift_invert_cache_key"
  key
}

#' @keywords internal
shift_invert_factorization_cache_info <- function(Aop, sigma, Bop = NULL,
                                                  label_kind = NA_character_) {
  list(
    key = shift_invert_factorization_cache_key(Aop, sigma, Bop = Bop),
    label_kind = label_kind,
    native = FALSE,
    reusable_within_operator = TRUE,
    external_cache = FALSE
  )
}

#' @keywords internal
shift_invert_factorization_contract <- function(cache) {
  label_kind <- cache$label_kind %||% NA_character_
  native <- isTRUE(cache$native)
  external <- isTRUE(cache$external_cache)
  provider <- if (native) {
    "eigencore_native_factorization"
  } else if (external || identical(label_kind, "user_solve")) {
    "user_supplied_solve"
  } else if (isTRUE(grepl("sparse_lu", label_kind, fixed = TRUE))) {
    "Matrix::lu_reference_factorization"
  } else {
    "eigencore_reference_factorization"
  }
  memory_policy <- if (native) {
    "native_factorized_apply_no_dense_fallback"
  } else if (external || identical(label_kind, "user_solve")) {
    "external_cache_user_owned_no_dense_fallback"
  } else if (isTRUE(grepl("sparse_lu", label_kind, fixed = TRUE))) {
    "sparse_factorization_no_dense_rcond"
  } else {
    "reference_factorization_no_silent_densification"
  }
  list(
    contract_version = "shift_invert_factorization_contract_v1",
    label_kind = label_kind,
    provider = provider,
    promotion_status = if (native) "promoted_native" else "reference_boundary",
    owned_by_eigencore = native,
    external_cache = external,
    generalized = isTRUE(cache$generalized),
    cache_key_scope = "A_fingerprint+B_fingerprint+sigma+structure",
    cache_invalidation = "cache key changes when A, B, sigma, or structure changes",
    memory_policy = memory_policy,
    certificate_policy = "original_coordinate_residual_required",
    native_label_requires_owned_factorized_apply = TRUE
  )
}

#' @keywords internal
shift_invert_factorization_cache_merge <- function(cache_info, label_kind,
                                                   diagnostics = list()) {
  cache <- modifyList(
    modifyList(cache_info, list(label_kind = label_kind)),
    diagnostics
  )
  cache$contract <- shift_invert_factorization_contract(cache)
  cache
}

#' @keywords internal
shift_invert_solver_dense <- function(A, sigma, B = NULL) {
  # Form M = A - sigma * B without materialising an n x n identity for the
  # standard problem, and factor it exactly once.
  if (is.null(B)) {
    M <- A
    diag(M) <- diag(M) - sigma
  } else {
    M <- A - sigma * B
  }
  # LAPACK QR (column pivoting) with an explicit diag(R) tolerance gives a
  # stricter rank check than base::qr(LINPACK), whose qr.coef silently returns
  # NA on rank-deficient columns instead of erroring.
  factor <- tryCatch(qr(M, LAPACK = TRUE), error = function(e) NULL)
  R_factor <- if (is.null(factor)) NULL else qr.R(factor)
  # Reciprocal condition estimate from the same factorisation: M P = Q R with
  # Q orthogonal, so cond(M) and cond(R) agree up to the norm-equivalence
  # factor; dtrcon on R costs O(n^2) instead of a second O(n^3) LU.
  cond <- if (is.null(R_factor)) {
    NA_real_
  } else {
    tryCatch(base::rcond(R_factor, triangular = TRUE), error = function(e) NA_real_)
  }
  if (!is.finite(cond) || cond <= sqrt(.Machine$double.eps)) {
    stop(
      "shift_invert(sigma = ", sigma, ") produced a singular or near-singular ",
      "dense shifted operator; perturb sigma or supply a stable solve function.",
      call. = FALSE
    )
  }
  R_diag <- abs(diag(R_factor))
  rank_tol <- max(dim(M)) * .Machine$double.eps * max(R_diag, 1)
  if (any(R_diag <= rank_tol)) {
    stop(
      "shift_invert(sigma = ", sigma, ") produced a rank-deficient dense ",
      "shifted operator (LAPACK QR): smallest |R[i,i]| = ", min(R_diag),
      ". Perturb sigma or supply a stable solve function.",
      call. = FALSE
    )
  }
  list(
    solve_fn = function(X) {
      Z <- solve(factor, X)
      if (anyNA(Z)) {
        stop(
          "shift_invert(sigma = ", sigma, ") solve returned NA values for the ",
          "dense shifted operator; the factorization is silently rank-deficient. ",
          "Perturb sigma or supply a stable solve function.",
          call. = FALSE
        )
      }
      Z
    },
    label = "dense_qr",
    factor = factor,
    M = M,
    cache = list(
      factorization = "base::qr(LAPACK=TRUE)",
      factorization_cached = TRUE,
      condition_estimate = cond,
      condition_estimate_type = "dense_qr_triangular_rcond",
      near_singular = FALSE
    )
  )
}

#' @keywords internal
shift_invert_solver_csc <- function(A, sigma, B = NULL) {
  n <- nrow(A)
  if (is.null(B)) {
    B <- Matrix::Diagonal(n)
  }
  M <- methods::as(A - sigma * B, "CsparseMatrix")
  # Shift-invert is Hermitian-only and B is diagonal here, so M has a
  # symmetric pattern: order = 1 (AMD on A + A') is the matching fill-reducing
  # ordering. Matrix's default (AMD on A'A, order = 2/NA) produced 2.2x the fill
  # and 3x the factor time on random sparse symmetric matrices. Partial
  # pivoting (tol = 1) is unchanged. Older Matrix versions without an integer
  # `order` fall back to the default.
  factor <- tryCatch(
    tryCatch(Matrix::lu(M, order = 1L), error = function(e) Matrix::lu(M)),
    error = function(e) {
      stop(
        "shift_invert(sigma = ", sigma, ") could not factor the sparse shifted ",
        "operator; perturb sigma or supply a stable solve function. ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  # Cheap non-densifying near-singular diagnostic: min/max of |diag(U)| from
  # the sparseLU factorization is O(n) work and stays sparse. A small ratio
  # indicates a near-singular shifted operator without forming a dense rcond.
  # This is an LU-pivot estimate, not a true condition number — it bounds
  # the smallest singular value from above relative to the largest pivot.
  u_diag_estimate <- tryCatch({
    U_factor <- methods::slot(factor, "U")
    u_diag <- abs(as.numeric(Matrix::diag(U_factor)))
    if (length(u_diag) == 0L) {
      list(min = NA_real_, max = NA_real_, ratio = NA_real_)
    } else {
      max_u <- max(u_diag, na.rm = TRUE)
      min_u <- min(u_diag, na.rm = TRUE)
      ratio <- if (is.finite(max_u) && max_u > 0) min_u / max_u else NA_real_
      list(min = min_u, max = max_u, ratio = ratio)
    }
  }, error = function(e) {
    list(min = NA_real_, max = NA_real_, ratio = NA_real_)
  })
  near_singular <- if (is.finite(u_diag_estimate$ratio)) {
    u_diag_estimate$ratio <= sqrt(.Machine$double.eps)
  } else {
    NA
  }
  # Matrix::solve() refuses a sparse LU whose pivot ratio is below eps
  # ("computationally singular"); detect that here, at factor time, so the
  # caller's singular-shift handling (the implicit smallest_magnitude route
  # perturbs sigma) sees it instead of a callback failure deep inside the
  # Krylov-Schur kernel ("operator apply failed with status=-8", oracle sweep).
  if (isTRUE(u_diag_estimate$ratio < .Machine$double.eps) ||
      isTRUE(u_diag_estimate$min == 0)) {
    stop("shift_invert(sigma = ", sigma, "): the shifted sparse operator is ",
         "numerically singular (LU pivot ratio ",
         format(u_diag_estimate$ratio, digits = 3), "); perturb sigma.",
         call. = FALSE)
  }
  list(
    solve_fn = function(X) {
      Z <- Matrix::solve(factor, X)
      if (inherits(Z, "Matrix")) as.matrix(Z) else Z
    },
    label = "sparse_lu",
    factor = factor,
    M = M,
    cache = list(
      factorization = "Matrix::lu",
      factorization_cached = TRUE,
      condition_estimate = u_diag_estimate$ratio,
      condition_estimate_type = if (is.finite(u_diag_estimate$ratio)) {
        "sparse_lu_pivot_ratio"
      } else {
        "uncomputed_sparse_no_dense_rcond"
      },
      condition_estimate_min_pivot = u_diag_estimate$min,
      condition_estimate_max_pivot = u_diag_estimate$max,
      near_singular = near_singular
    )
  )
}

#' @keywords internal
shift_invert_metric_factor <- function(Bop) {
  Bop <- as_operator(Bop)
  if (!identical(Bop$structure$kind, "hermitian")) {
    stop("generalized shift_invert() requires a Hermitian metric B.", call. = FALSE)
  }
  if (!generalized_spd_metric_known(Bop)) {
    stop("generalized shift_invert() requires known positive definite B.", call. = FALSE)
  }
  Bsource <- source_or_null(Bop)
  Bstorage <- Bop$metadata$storage %||% NULL
  if (is.matrix(Bsource) && is.double(Bsource)) {
    B <- (Bsource + t(Bsource)) / 2
    R <- chol(B)
    return(list(
      kind = "dense",
      matrix_dense = B,
      to_rhs = function(X) crossprod(R, X),
      from_solution = function(X) R %*% X,
      to_original = function(Y) backsolve(R, Y),
      factorization = "base::chol(B)"
    ))
  }
  if (identical(Bstorage, "ddiMatrix")) {
    B <- Bop$metadata$matrix
    values <- if (identical(methods::slot(B, "diag"), "U")) {
      rep(1, Bop$dim[1L])
    } else {
      methods::slot(B, "x")
    }
    if (length(values) != Bop$dim[1L] ||
        any(!is.finite(values)) ||
        any(values <= 0)) {
      stop("generalized shift_invert() requires positive finite diagonal B.", call. = FALSE)
    }
    sqrt_values <- sqrt(values)
    return(list(
      kind = "diagonal",
      matrix_dense = diag(values),
      matrix_sparse = Matrix::Diagonal(x = values),
      to_rhs = function(X) sqrt_values * X,
      from_solution = function(X) sqrt_values * X,
      to_original = function(Y) Y / sqrt_values,
      factorization = "diagonal sqrt(B)"
    ))
  }
  stop(
    "generalized shift_invert() supports dense B or diagonal B only; ",
    "unsupported metric operators are rejected to avoid silent densification.",
    call. = FALSE
  )
}

#' @keywords internal
prepare_shift_invert_operator <- function(problem, sigma, user_solve = NULL) {
  Aop <- problem$A
  Bop <- problem$metric
  n <- Aop$dim[1L]
  cache_info <- shift_invert_factorization_cache_info(Aop, sigma, Bop = Bop)

  if (!identical(problem$structure$kind, "hermitian")) {
    stop("shift_invert() currently requires a Hermitian eigenproblem.", call. = FALSE)
  }

  metric_factor <- if (is.null(Bop)) NULL else shift_invert_metric_factor(Bop)

  # User-supplied solve operator for matrix-free A
  if (!is.null(user_solve)) {
    if (!is.function(user_solve)) {
      stop("shift_invert(solve = ...) must be a function.", call. = FALSE)
    }
    solve_fn <- if (is.null(metric_factor)) {
      user_solve
    } else {
      function(X) metric_factor$from_solution(user_solve(metric_factor$to_rhs(X)))
    }
    op <- linear_operator(
      dim = c(n, n),
      apply = shift_invert_apply_factory(solve_fn),
      apply_adjoint = NULL,
      structure = hermitian(),
      name = "shift_invert_user_solve",
      metadata = list(
        native = FALSE,
        factorization_cache = shift_invert_factorization_cache_merge(
          cache_info,
          "user_solve",
          list(
            factorization = "user_solve",
            factorization_cached = NA,
            condition_estimate = NA_real_,
            condition_estimate_type = "user_supplied",
            near_singular = NA,
            external_cache = TRUE,
            generalized = !is.null(Bop),
            metric_factorization = metric_factor$factorization %||% NA_character_
          )
        )
      )
    )
    cache <- op$metadata$factorization_cache
    return(list(
      operator = op,
      label_kind = "user_solve",
      factorization_cache = cache,
      recover_vectors = if (is.null(metric_factor)) {
        function(Y) Y
      } else {
        metric_factor$to_original
      }
    ))
  }

  # Build factorization from source matrix; CSC operators carry the matrix in
  # metadata$matrix rather than metadata$source.
  source_A <- source_or_null(Aop)
  csc_A    <- if (identical(Aop$metadata$storage %||% NULL, "dgCMatrix")) {
    Aop$metadata$matrix
  } else if (inherits(Aop$metadata$matrix, "CsparseMatrix")) {
    Aop$metadata$matrix
  } else {
    NULL
  }

  prep <- if (is.matrix(source_A) && is.double(source_A)) {
    shift_invert_solver_dense(
      source_A,
      sigma,
      B = metric_factor$matrix_dense %||% NULL
    )
  } else if (!is.null(csc_A) && (is.null(metric_factor) || identical(metric_factor$kind, "diagonal"))) {
    shift_invert_solver_csc(
      methods::as(csc_A, "generalMatrix"),
      sigma,
      B = metric_factor$matrix_sparse %||% NULL
    )
  } else if ((inherits(source_A, "CsparseMatrix") || inherits(source_A, "dgCMatrix")) &&
      (is.null(metric_factor) || identical(metric_factor$kind, "diagonal"))) {
    shift_invert_solver_csc(
      source_A,
      sigma,
      B = metric_factor$matrix_sparse %||% NULL
    )
  } else {
    stop(
      if (is.null(Bop)) {
        "shift_invert() supports dense double matrices and dgCMatrix/dsCMatrix sources, or a user-supplied solve operator."
      } else {
        "generalized shift_invert() supports dense A/B or sparse A with diagonal B; unsupported combinations are rejected to avoid silent densification."
      },
      call. = FALSE
    )
  }

  solve_fn <- if (is.null(metric_factor)) {
    prep$solve_fn
  } else {
    function(X) metric_factor$from_solution(prep$solve_fn(metric_factor$to_rhs(X)))
  }
  label_kind <- if (is.null(metric_factor)) prep$label else paste0(prep$label, "_generalized")

  op <- linear_operator(
    dim = c(n, n),
    apply = shift_invert_apply_factory(solve_fn),
    apply_adjoint = NULL,
    structure = hermitian(),
    name = paste0("shift_invert_", label_kind),
    metadata = list(
      native = FALSE,
      factorization_cache = shift_invert_factorization_cache_merge(
        cache_info,
        label_kind,
        modifyList(
          prep$cache,
          list(
            generalized = !is.null(Bop),
            metric_factorization = metric_factor$factorization %||% NA_character_
          )
        )
      )
    )
  )
  cache <- op$metadata$factorization_cache
  list(
    operator = op,
    label_kind = label_kind,
    factorization_cache = cache,
    recover_vectors = if (is.null(metric_factor)) {
      function(Y) Y
    } else {
      metric_factor$to_original
    }
  )
}

#' @keywords internal
native_dense_shift_invert_lanczos <- function(problem, k, sigma, tol, maxit,
                                              vectors, certify, plan) {
  Aop <- problem$A
  source_A <- source_or_null(Aop)
  if (!(is.matrix(source_A) && is.double(source_A))) {
    stop("native dense shift-invert requires a dense double source.", call. = FALSE)
  }
  if (!is.null(problem$metric)) {
    stop("native dense shift-invert currently supports standard eigenproblems only.", call. = FALSE)
  }

  n <- Aop$dim[1L]
  effective_maxit <- maxit %||% min(n, max(20L, 4L * as.integer(k) + 20L))
  start <- stats::rnorm(n)
  native <- .Call(
    "eigencore_shift_invert_lanczos_dense",
    source_A,
    as.numeric(sigma),
    as.integer(effective_maxit),
    as.numeric(start),
    as.integer(k),
    as.integer(lanczos_target_kind(largest_magnitude())),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )

  iterations <- as.integer(native$iterations)
  alpha <- native$alpha[seq_len(iterations)]
  beta <- native$beta[seq_len(iterations)]
  eig <- native_tridiagonal_eigen(alpha, beta)
  idx <- order_indices(eig$values, largest_magnitude())
  idx <- idx[seq_len(min(as.integer(k), length(idx)))]
  mu <- eig$values[idx]
  if (any(abs(mu) < .Machine$double.eps)) {
    stop(
      "native shift_invert(sigma = ", sigma, ") produced a zero-magnitude ",
      "eigenvalue of the inverted operator; sigma is too close to a true ",
      "eigenvalue. Perturb sigma or use a tighter tolerance.",
      call. = FALSE
    )
  }

  vec <- native$Q[, seq_len(iterations), drop = FALSE] %*%
    eig$vectors[, idx, drop = FALSE]
  lambda <- sigma + 1 / mu
  ord <- order_indices(lambda, problem$target)
  if (length(ord) > k) ord <- ord[seq_len(k)]
  lambda <- lambda[ord]
  vec <- vec[, ord, drop = FALSE]

  cert <- if (isTRUE(certify) && ncol(vec) > 0L) {
    certify_eigen_operator(Aop, lambda, vec, tol = tol)
  } else {
    empty_certificate(
      tol,
      note = if (!isTRUE(certify)) {
        "native shift-invert: certification disabled by caller"
      } else {
        "native shift-invert: no eigenpairs returned; residual certificate not computed"
      }
    )
  }

  cache_info <- shift_invert_factorization_cache_info(
    Aop,
    sigma,
    label_kind = "dense_lu_native"
  )
  cache <- shift_invert_factorization_cache_merge(
    cache_info,
    "dense_lu_native",
    modifyList(
      native$factorization_cache,
      list(
        native = TRUE,
        condition_estimate_type = "dense_lu_pivot_ratio",
        near_singular = FALSE,
        external_cache = FALSE,
        generalized = FALSE,
        metric_factorization = NA_character_
      )
    )
  )

  result <- list(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = iterations,
    matvecs = as.integer(native$matvecs),
    method = plan$method,
    target = target_label(problem$target),
    plan = plan,
    certificate = cert,
    sigma = sigma,
    transform = list(
      kind = "shift_invert",
      sigma = sigma,
      label_kind = "dense_lu_native",
      factorization_cache = cache,
      certification = list(
        problem = "original",
        residual_formula = "A * x - lambda * x",
        transformed_residuals_used = FALSE
      )
    ),
    warnings = character(),
    restart = list(
      kind = "native_dense_shift_invert_lanczos",
      native = TRUE,
      factorization_native = TRUE,
      factorization = cache$factorization,
      max_subspace = effective_maxit,
      transformed_operator_target = "largest_magnitude",
      eigenvalue_recovery = "lambda = sigma + 1 / mu",
      history_nconv = native$history_nconv,
      history_max_residual = native$history_max_residual
    )
  )
  result <- finalize_workflow_result(result, plan)
  class(result) <- "eigencore_eigen_result"
  result
}

#' @keywords internal
native_tridiagonal_shift_invert_lanczos <- function(problem, k, sigma, tol,
                                                    maxit, vectors, certify,
                                                    plan) {
  Aop <- problem$A
  if (!is.null(problem$metric)) {
    stop("native tridiagonal shift-invert supports standard eigenproblems only.", call. = FALSE)
  }
  A <- Aop$metadata$matrix %||% source_or_null(Aop)
  if (!(inherits(A, "CsparseMatrix") || inherits(A, "diagonalMatrix"))) {
    stop("native tridiagonal shift-invert requires a CSC sparse or diagonal source.", call. = FALSE)
  }
  parts <- shift_invert_problem_tridiagonal_parts(problem)
  if (is.null(parts)) {
    stop("native tridiagonal shift-invert requires a symmetric tridiagonal CSC source.", call. = FALSE)
  }
  parts$diag <- parts$diag - as.numeric(sigma)

  n <- Aop$dim[1L]
  effective_maxit <- maxit %||% min(n, max(20L, 4L * as.integer(k) + 20L))
  start <- stats::rnorm(n)
  native <- .Call(
    "eigencore_shift_invert_lanczos_tridiagonal",
    as.numeric(parts$lower),
    as.numeric(parts$diag),
    as.numeric(parts$upper),
    as.integer(effective_maxit),
    as.numeric(start),
    as.integer(k),
    as.integer(lanczos_target_kind(largest_magnitude())),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )

  iterations <- as.integer(native$iterations)
  mu <- native$ritz_values
  if (any(abs(mu) < .Machine$double.eps)) {
    stop(
      "native tridiagonal shift_invert(sigma = ", sigma, ") produced a zero-magnitude ",
      "eigenvalue of the inverted operator; sigma is too close to a true ",
      "eigenvalue. Perturb sigma or use a tighter tolerance.",
      call. = FALSE
    )
  }

  vec <- native$ritz_vectors
  lambda <- sigma + 1 / mu
  ord <- order_indices(lambda, problem$target)
  if (length(ord) > k) ord <- ord[seq_len(k)]
  lambda <- lambda[ord]
  vec <- vec[, ord, drop = FALSE]

  cert <- if (isTRUE(certify) && ncol(vec) > 0L) {
    certify_eigen_operator(Aop, lambda, vec, tol = tol)
  } else {
    empty_certificate(
      tol,
      note = if (!isTRUE(certify)) {
        "native tridiagonal shift-invert: certification disabled by caller"
      } else {
        "native tridiagonal shift-invert: no eigenpairs returned; residual certificate not computed"
      }
    )
  }

  cache_info <- shift_invert_factorization_cache_info(
    Aop,
    sigma,
    label_kind = "tridiagonal_lu_native"
  )
  cache <- shift_invert_factorization_cache_merge(
    cache_info,
    "tridiagonal_lu_native",
    modifyList(
      native$factorization_cache,
      list(
        native = TRUE,
        condition_estimate_type = "tridiagonal_lu_pivot_ratio",
        near_singular = FALSE,
        external_cache = FALSE,
        generalized = FALSE,
        metric_factorization = NA_character_
      )
    )
  )

  result <- list(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = iterations,
    matvecs = as.integer(native$matvecs),
    method = plan$method,
    target = target_label(problem$target),
    plan = plan,
    certificate = cert,
    sigma = sigma,
    transform = list(
      kind = "shift_invert",
      sigma = sigma,
      label_kind = "tridiagonal_lu_native",
      factorization_cache = cache,
      certification = list(
        problem = "original",
        residual_formula = "A * x - lambda * x",
        transformed_residuals_used = FALSE
      )
    ),
    warnings = character(),
    restart = list(
      kind = "native_tridiagonal_shift_invert_lanczos",
      native = TRUE,
      factorization_native = TRUE,
      factorization = cache$factorization,
      max_subspace = effective_maxit,
      transformed_operator_target = "largest_magnitude",
      eigenvalue_recovery = "lambda = sigma + 1 / mu",
      history_nconv = native$history_nconv,
      history_max_residual = native$history_max_residual
    )
  )
  result <- finalize_workflow_result(result, plan)
  class(result) <- "eigencore_eigen_result"
  result
}

#' @keywords internal
native_tridiagonal_shift_invert_retryable_error <- function(error) {
  message <- conditionMessage(error)
  grepl("zero .*pivot|near-singular|zero-magnitude", message)
}

#' @keywords internal
native_tridiagonal_shift_invert_candidate_sigmas <- function(parts, sigma) {
  sigma <- as.numeric(sigma)
  if (length(sigma) != 1L || !is.finite(sigma)) {
    return(numeric())
  }
  bounds <- tridiagonal_gershgorin_bounds(parts)
  diag <- as.numeric(parts$diag)
  offdiag <- abs(as.numeric(parts$upper))
  scale <- max(1, abs(sigma), abs(bounds$lower), abs(bounds$upper), abs(diag), offdiag)
  span <- max(bounds$upper - bounds$lower, scale)
  margin <- max(1e-8 * span, 10 * sqrt(.Machine$double.eps) * scale)
  offsets <- margin * c(1, -1, 10, -10, 100, -100, 1000, -1000, 10000, -10000)
  candidates <- unique(c(sigma, sigma + offsets))
  candidates[is.finite(candidates)]
}

#' @keywords internal
native_tridiagonal_shift_invert_lanczos_with_perturbation <- function(problem, k,
                                                                      sigma, tol,
                                                                      maxit,
                                                                      vectors,
                                                                      certify,
                                                                      plan) {
  parts <- shift_invert_problem_tridiagonal_parts(problem)
  if (is.null(parts)) {
    return(native_tridiagonal_shift_invert_lanczos(
      problem, k = k, sigma = sigma, tol = tol, maxit = maxit,
      vectors = vectors, certify = certify, plan = plan
    ))
  }

  candidates <- native_tridiagonal_shift_invert_candidate_sigmas(parts, sigma)
  last_error <- NULL
  for (candidate in candidates) {
    result <- tryCatch(
      native_tridiagonal_shift_invert_lanczos(
        problem, k = k, sigma = candidate, tol = tol, maxit = maxit,
        vectors = vectors, certify = certify, plan = plan
      ),
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      if (!isTRUE(all.equal(candidate, as.numeric(sigma)))) {
        delta <- candidate - as.numeric(sigma)
        note <- paste0(
          "native tridiagonal shift-invert perturbed requested sigma from ",
          format(as.numeric(sigma), digits = 17),
          " to ",
          format(candidate, digits = 17),
          " after a singular or near-singular pivoted tridiagonal LU (dgttrf) factorization"
        )
        result$warnings <- c(result$warnings, note)
        result$transform$requested_sigma <- as.numeric(sigma)
        result$transform$sigma_perturbed <- TRUE
        result$transform$sigma_perturbation <- delta
        result$restart$requested_sigma <- as.numeric(sigma)
        result$restart$sigma_perturbed <- TRUE
        result$restart$sigma_perturbation <- delta
        result$restart$perturbation_reason <- conditionMessage(last_error)
      }
      return(result)
    }
    if (!native_tridiagonal_shift_invert_retryable_error(result)) {
      stop(result)
    }
    last_error <- result
  }

  stop(
    "native tridiagonal shift_invert(sigma = ", sigma,
    ") failed at the requested shift and all perturbation retries. Last error: ",
    conditionMessage(last_error),
    call. = FALSE
  )
}

#' @keywords internal
native_tridiagonal_generalized_shift_invert_lanczos <- function(problem, k,
                                                               sigma, tol,
                                                               maxit, vectors,
                                                               certify, plan) {
  Aop <- problem$A
  Bop <- problem$metric
  if (is.null(Bop)) {
    stop("native generalized tridiagonal shift-invert requires a metric B.", call. = FALSE)
  }
  A <- Aop$metadata$matrix %||% source_or_null(Aop)
  if (!(inherits(A, "CsparseMatrix") || inherits(A, "diagonalMatrix"))) {
    stop("native generalized tridiagonal shift-invert requires a CSC sparse or diagonal A source.", call. = FALSE)
  }
  metric_values <- shift_invert_diagonal_metric_values(Bop)
  if (is.null(metric_values)) {
    stop("native generalized tridiagonal shift-invert requires positive diagonal B.", call. = FALSE)
  }
  if (length(metric_values) != Aop$dim[1L]) {
    stop("native generalized tridiagonal shift-invert requires conformable diagonal B.", call. = FALSE)
  }
  parts <- shift_invert_problem_tridiagonal_parts(problem)
  if (is.null(parts)) {
    stop("native generalized tridiagonal shift-invert requires a symmetric tridiagonal CSC source.", call. = FALSE)
  }
  shifted_diag <- parts$diag - as.numeric(sigma) * metric_values
  sqrt_metric <- sqrt(metric_values)

  n <- Aop$dim[1L]
  effective_maxit <- maxit %||% min(n, max(20L, 4L * as.integer(k) + 20L))
  start <- stats::rnorm(n)
  native <- .Call(
    "eigencore_shift_invert_lanczos_tridiagonal_generalized",
    as.numeric(parts$lower),
    as.numeric(shifted_diag),
    as.numeric(parts$upper),
    as.numeric(sqrt_metric),
    as.integer(effective_maxit),
    as.numeric(start),
    as.integer(k),
    as.integer(lanczos_target_kind(largest_magnitude())),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )

  iterations <- as.integer(native$iterations)
  alpha <- native$alpha[seq_len(iterations)]
  beta <- native$beta[seq_len(iterations)]
  eig <- native_tridiagonal_eigen(alpha, beta)
  idx <- order_indices(eig$values, largest_magnitude())
  idx <- idx[seq_len(min(as.integer(k), length(idx)))]
  mu <- eig$values[idx]
  if (any(abs(mu) < .Machine$double.eps)) {
    stop(
      "native generalized tridiagonal shift_invert(sigma = ", sigma, ") produced a zero-magnitude ",
      "eigenvalue of the inverted operator; sigma is too close to a true ",
      "eigenvalue. Perturb sigma or use a tighter tolerance.",
      call. = FALSE
    )
  }

  vec_transformed <- native$Q[, seq_len(iterations), drop = FALSE] %*%
    eig$vectors[, idx, drop = FALSE]
  vec <- vec_transformed / sqrt_metric
  lambda <- sigma + 1 / mu
  ord <- order_indices(lambda, problem$target)
  if (length(ord) > k) ord <- ord[seq_len(k)]
  lambda <- lambda[ord]
  vec <- vec[, ord, drop = FALSE]

  cert <- if (isTRUE(certify) && ncol(vec) > 0L) {
    certify_eigen_operator(Aop, lambda, vec, Bop = Bop, tol = tol)
  } else {
    empty_certificate(
      tol,
      note = if (!isTRUE(certify)) {
        "native generalized tridiagonal shift-invert: certification disabled by caller"
      } else {
        "native generalized tridiagonal shift-invert: no eigenpairs returned; residual certificate not computed"
      }
    )
  }

  cache_info <- shift_invert_factorization_cache_info(
    Aop,
    sigma,
    Bop = Bop,
    label_kind = "tridiagonal_lu_generalized_native"
  )
  cache <- shift_invert_factorization_cache_merge(
    cache_info,
    "tridiagonal_lu_generalized_native",
    modifyList(
      native$factorization_cache,
      list(
        native = TRUE,
        condition_estimate_type = "tridiagonal_lu_pivot_ratio",
        near_singular = FALSE,
        external_cache = FALSE,
        generalized = TRUE,
        metric_factorization = "diagonal sqrt(B)"
      )
    )
  )

  result <- list(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = iterations,
    matvecs = as.integer(native$matvecs),
    method = plan$method,
    target = target_label(problem$target),
    plan = plan,
    certificate = cert,
    sigma = sigma,
    transform = list(
      kind = "shift_invert",
      sigma = sigma,
      label_kind = "tridiagonal_lu_generalized_native",
      factorization_cache = cache,
      certification = list(
        problem = "original",
        residual_formula = "A * x - lambda * B * x",
        transformed_residuals_used = FALSE
      )
    ),
    warnings = character(),
    restart = list(
      kind = "native_tridiagonal_generalized_shift_invert_lanczos",
      native = TRUE,
      generalized = TRUE,
      factorization_native = TRUE,
      factorization = cache$factorization,
      max_subspace = effective_maxit,
      transformed_operator_target = "largest_magnitude",
      eigenvalue_recovery = "lambda = sigma + 1 / mu",
      history_nconv = native$history_nconv,
      history_max_residual = native$history_max_residual
    )
  )
  result <- finalize_workflow_result(result, plan)
  class(result) <- "eigencore_eigen_result"
  result
}

#' @keywords internal
native_dense_generalized_shift_invert_lanczos <- function(problem, k, sigma,
                                                         tol, maxit, vectors,
                                                         certify, plan) {
  Aop <- problem$A
  Bop <- problem$metric
  source_A <- source_or_null(Aop)
  source_B <- source_or_null(Bop)
  if (!(is.matrix(source_A) && is.double(source_A))) {
    stop("native dense generalized shift-invert requires a dense double A source.", call. = FALSE)
  }
  if (!(is.matrix(source_B) && is.double(source_B))) {
    stop("native dense generalized shift-invert requires a dense double B source.", call. = FALSE)
  }

  n <- Aop$dim[1L]
  effective_maxit <- maxit %||% min(n, max(20L, 4L * as.integer(k) + 20L))
  start <- stats::rnorm(n)
  native <- .Call(
    "eigencore_shift_invert_lanczos_dense_generalized",
    source_A,
    source_B,
    as.numeric(sigma),
    as.integer(effective_maxit),
    as.numeric(start),
    as.integer(k),
    as.integer(lanczos_target_kind(largest_magnitude())),
    as.numeric(tol),
    PACKAGE = "eigencore"
  )

  iterations <- as.integer(native$iterations)
  alpha <- native$alpha[seq_len(iterations)]
  beta <- native$beta[seq_len(iterations)]
  eig <- native_tridiagonal_eigen(alpha, beta)
  idx <- order_indices(eig$values, largest_magnitude())
  idx <- idx[seq_len(min(as.integer(k), length(idx)))]
  mu <- eig$values[idx]
  if (any(abs(mu) < .Machine$double.eps)) {
    stop(
      "native generalized shift_invert(sigma = ", sigma, ") produced a zero-magnitude ",
      "eigenvalue of the inverted operator; sigma is too close to a true ",
      "eigenvalue. Perturb sigma or use a tighter tolerance.",
      call. = FALSE
    )
  }

  vec_transformed <- native$Q[, seq_len(iterations), drop = FALSE] %*%
    eig$vectors[, idx, drop = FALSE]
  vec <- backsolve(native$chol_factor, vec_transformed)
  lambda <- sigma + 1 / mu
  ord <- order_indices(lambda, problem$target)
  if (length(ord) > k) ord <- ord[seq_len(k)]
  lambda <- lambda[ord]
  vec <- vec[, ord, drop = FALSE]

  cert <- if (isTRUE(certify) && ncol(vec) > 0L) {
    certify_eigen_operator(Aop, lambda, vec, Bop = Bop, tol = tol)
  } else {
    empty_certificate(
      tol,
      note = if (!isTRUE(certify)) {
        "native generalized shift-invert: certification disabled by caller"
      } else {
        "native generalized shift-invert: no eigenpairs returned; residual certificate not computed"
      }
    )
  }

  cache_info <- shift_invert_factorization_cache_info(
    Aop,
    sigma,
    Bop = Bop,
    label_kind = "dense_lu_generalized_native"
  )
  cache <- shift_invert_factorization_cache_merge(
    cache_info,
    "dense_lu_generalized_native",
    modifyList(
      native$factorization_cache,
      list(
        native = TRUE,
        condition_estimate_type = "dense_lu_pivot_ratio",
        near_singular = FALSE,
        external_cache = FALSE,
        generalized = TRUE,
        metric_factorization = native$factorization_cache$metric_factorization %||%
          "LAPACK dpotrf(B)"
      )
    )
  )

  result <- list(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = iterations,
    matvecs = as.integer(native$matvecs),
    method = plan$method,
    target = target_label(problem$target),
    plan = plan,
    certificate = cert,
    sigma = sigma,
    transform = list(
      kind = "shift_invert",
      sigma = sigma,
      label_kind = "dense_lu_generalized_native",
      factorization_cache = cache,
      certification = list(
        problem = "original",
        residual_formula = "A * x - lambda * B * x",
        transformed_residuals_used = FALSE
      )
    ),
    warnings = character(),
    restart = list(
      kind = "native_dense_generalized_shift_invert_lanczos",
      native = TRUE,
      generalized = TRUE,
      factorization_native = TRUE,
      factorization = cache$factorization,
      max_subspace = effective_maxit,
      transformed_operator_target = "largest_magnitude",
      eigenvalue_recovery = "lambda = sigma + 1 / mu",
      history_nconv = native$history_nconv,
      history_max_residual = native$history_max_residual
    )
  )
  result <- finalize_workflow_result(result, plan)
  class(result) <- "eigencore_eigen_result"
  result
}

#' @keywords internal
#' Largest-magnitude eigenpairs of the factorised shift-invert operator M.
#' M is a matrix-free Hermitian callback (its apply is a factorised solve), so
#' it drives the native thick-restart Lanczos kernel through the callback ABI
#' (restarts, locking, active subspace bounded by `maxit`, which here is the
#' internal subspace size). The reference scalar Lanczos remains for subspaces
#' too small for a thick restart (k + 1 > maxit).
shift_invert_transformed_lanczos <- function(M, k, tol, maxit,
                                             max_restarts = 100L) {
  n <- M$dim[1L]
  k <- as.integer(k)
  m_max <- min(n, as.integer(maxit))
  if (!is.na(m_max) && m_max >= k + 1L &&
      native_matrix_free_block_lanczos_available(M)) {
    iter <- native_block_lanczos_hermitian(
      M,
      k = k,
      target = largest_magnitude(),
      tol = tol,
      maxit = m_max,
      block = 1L,
      max_restarts = as.integer(max_restarts),
      vectors = TRUE,
      full_subspace = FALSE,
      certificate_fallback = FALSE
    )
    iter$restart$kind <- "native_thick_restart_shift_invert_callback"
    iter$restart$native_shift_invert_callback <- TRUE
    return(iter)
  }
  reference_lanczos_hermitian(
    M,
    k = k,
    target = largest_magnitude(),
    tol = tol,
    maxit = maxit,
    vectors = TRUE,
    reorthogonalize = TRUE
  )
}

#' @keywords internal
solve_shift_invert_hermitian <- function(problem, k, method, tol, maxit,
                                          vectors, certify, plan) {
  sigma <- method$sigma
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("shift_invert(sigma) requires a single finite numeric shift.", call. = FALSE)
  }
  if (!is.null(method$factorization)) {
    stop(
      "shift_invert(factorization = ...) is not implemented yet; ",
      "supply shift_invert(solve = ...) for a user-managed factorization cache.",
      call. = FALSE
    )
  }
  # `maxit` is an iteration limit (C15): the planner has already resolved it
  # into controls$max_restarts (thick-restart callback routes) or capped
  # controls$max_subspace (unrestarted native Lanczos routes, where an
  # iteration is one Lanczos step). The subspace size comes only from
  # controls$max_subspace.
  controls <- plan$controls %||% list()
  subspace0 <- as.integer(
    controls$max_subspace %||%
      default_shift_invert_max_subspace(problem$A$dim[1L], k)
  )
  restart_limit <- as.integer(controls$max_restarts %||% maxit %||% 100L)

  if (identical(plan$method, native_dense_shift_invert_label()) &&
      is.null(problem$metric) &&
      is.null(method$solve)) {
    return(native_dense_shift_invert_lanczos(
      problem, k = k, sigma = sigma, tol = tol, maxit = subspace0,
      vectors = vectors, certify = certify, plan = plan
    ))
  }
  if (identical(plan$method, native_tridiagonal_shift_invert_label()) &&
      is.null(problem$metric) &&
      is.null(method$solve)) {
    return(native_tridiagonal_shift_invert_lanczos_with_perturbation(
      problem, k = k, sigma = sigma, tol = tol, maxit = subspace0,
      vectors = vectors, certify = certify, plan = plan
    ))
  }
  if (identical(plan$method, native_tridiagonal_generalized_shift_invert_label()) &&
      !is.null(problem$metric) &&
      is.null(method$solve)) {
    return(native_tridiagonal_generalized_shift_invert_lanczos(
      problem, k = k, sigma = sigma, tol = tol, maxit = subspace0,
      vectors = vectors, certify = certify, plan = plan
    ))
  }
  if (identical(plan$method, native_dense_generalized_shift_invert_label()) &&
      !is.null(problem$metric) &&
      is.null(method$solve)) {
    return(native_dense_generalized_shift_invert_lanczos(
      problem, k = k, sigma = sigma, tol = tol, maxit = subspace0,
      vectors = vectors, certify = certify, plan = plan
    ))
  }

  prep <- prepare_shift_invert_operator(problem, sigma, user_solve = method$solve)
  M <- prep$operator
  n <- M$dim[1L]

  effective_maxit <- subspace0
  Aop <- problem$A
  Bop <- problem$metric

  # Transformed-operator solve followed by the original-coordinate
  # certificate. Inner convergence on M = (A - sigma B)^{-1} does not imply
  # original-coordinate convergence (the residuals scale by ||A - sigma B|| /
  # |mu|), so an unconverged original certificate triggers a retry with a
  # tighter inner tolerance and a larger restart subspace instead of
  # returning the unconverged pairs silently.
  attempt <- function(inner_tol, subspace) {
    iter <- shift_invert_transformed_lanczos(M, k = k, tol = inner_tol,
                                             maxit = subspace,
                                             max_restarts = restart_limit)
    mu <- iter$values
    vec <- iter$vectors
    if (any(abs(mu) < .Machine$double.eps)) {
      stop(
        "shift_invert(sigma = ", sigma, ") produced a zero-magnitude eigenvalue ",
        "of the inverted operator; sigma is too close to a true eigenvalue. ",
        "Perturb sigma or use a tighter tolerance.",
        call. = FALSE
      )
    }
    lambda <- sigma + 1 / mu
    ord <- order_indices(lambda, problem$target)
    if (length(ord) > k) ord <- ord[seq_len(k)]
    lambda <- lambda[ord]
    vec <- prep$recover_vectors(vec[, ord, drop = FALSE])
    cert <- if (isTRUE(certify) && !is.null(vec) && ncol(vec) > 0L) {
      certify_eigen_operator(Aop, lambda, vec, Bop = Bop, tol = tol)
    } else {
      empty_certificate(
        tol,
        note = if (!isTRUE(certify)) {
          "shift-invert: certification disabled by caller"
        } else {
          "shift-invert: no eigenpairs returned; residual certificate not computed"
        }
      )
    }
    list(iter = iter, lambda = lambda, vec = vec, cert = cert)
  }

  inner_tol <- tol
  subspace <- effective_maxit
  current <- attempt(inner_tol, subspace)
  total_iterations <- as.integer(current$iter$iterations %||% 0L)
  total_matvecs <- as.integer(current$iter$matvecs %||% 0L)
  retries <- 0L
  native_inner <- isTRUE(current$iter$restart$native_shift_invert_callback)
  while (native_inner && isTRUE(certify) && retries < 2L &&
         length(current$cert$converged) &&
         !all(current$cert$converged)) {
    retries <- retries + 1L
    inner_tol <- max(inner_tol * 1e-3, 10 * .Machine$double.eps)
    subspace <- min(n, max(subspace, 2L * subspace))
    candidate <- attempt(inner_tol, subspace)
    total_iterations <- total_iterations +
      as.integer(candidate$iter$iterations %||% 0L)
    total_matvecs <- total_matvecs + as.integer(candidate$iter$matvecs %||% 0L)
    if (sum(candidate$cert$converged) >= sum(current$cert$converged)) {
      current <- candidate
    }
  }
  iter <- current$iter
  lambda <- current$lambda
  vec <- current$vec
  cert <- current$cert
  iter$iterations <- total_iterations
  iter$matvecs <- total_matvecs

  result <- list(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    residuals = cert$residuals,
    backward_error = cert$backward_error,
    orthogonality = cert$orthogonality,
    nconv = sum(cert$converged),
    requested = k,
    iterations = iter$iterations,
    matvecs = iter$matvecs,
    method = plan$method,
    target = target_label(problem$target),
    plan = plan,
    certificate = cert,
    sigma = sigma,
    transform = list(
      kind = "shift_invert",
      sigma = sigma,
      label_kind = prep$label_kind,
      factorization_cache = prep$factorization_cache,
      certification = list(
        problem = "original",
        residual_formula = if (is.null(Bop)) {
          "A * x - lambda * x"
        } else {
          "A * x - lambda * B * x"
        },
        transformed_residuals_used = FALSE
      )
    ),
    restart = list(
      kind = iter$restart$kind %||% "reference_hermitian_lanczos_shift_invert",
      native = isTRUE(iter$restart$native_shift_invert_callback),
      factorization_native = FALSE,
      max_subspace = subspace,
      restarts_used = as.integer(iter$restarts %||% 0L),
      inner_tolerance = inner_tol,
      certificate_retries = retries,
      transformed_operator_target = "largest_magnitude",
      eigenvalue_recovery = "lambda = sigma + 1 / mu"
    ),
    warnings = if (isTRUE(iter$restart$native_shift_invert_callback)) {
      if (!isTRUE(cert$passed) && isTRUE(certify)) {
        paste0(
          plan$method, " did not certify all ", k, " requested pairs after ",
          retries, " tightened restart(s); subspace ", subspace
        )
      } else {
        character()
      }
    } else {
      paste0(
        "using reference Hermitian Lanczos shift-invert (",
        prep$label_kind, "); native shift-invert hot loop not yet implemented"
      )
    }
  )
  result <- finalize_workflow_result(result, plan)
  class(result) <- "eigencore_eigen_result"
  result
}

#' @keywords internal
#' Factorised solve for the nonsymmetric shifted operator A - sigma I (dense
#' LAPACK QR or sparse LU with the AMD ordering used by the Hermitian path).
#' The transposed factorisation needed by the left (adjoint) solve is built
#' lazily, only when the adjoint operator is applied.
shift_invert_general_solver <- function(Aop, sigma, user_solve = NULL) {
  n <- Aop$dim[1L]
  if (!is.null(user_solve)) {
    if (!is.function(user_solve)) {
      stop("shift_invert(solve = ...) must be a function.", call. = FALSE)
    }
    return(list(
      solve_fn = user_solve,
      adjoint_solve_fn = NULL,
      label_kind = "user_solve",
      cache = list(
        factorization = "user_solve",
        factorization_cached = NA,
        condition_estimate = NA_real_,
        condition_estimate_type = "user_supplied",
        near_singular = NA,
        external_cache = TRUE,
        generalized = FALSE
      )
    ))
  }
  source_A <- source_or_null(Aop)
  csc_A <- if (inherits(Aop$metadata$matrix, "CsparseMatrix")) {
    methods::as(Aop$metadata$matrix, "generalMatrix")
  } else {
    NULL
  }
  dense <- is.matrix(source_A) && is.double(source_A)
  if (!dense && is.null(csc_A)) {
    stop(
      "nonsymmetric shift_invert() supports dense double matrices and sparse ",
      "CSC sources, or a user-supplied solve operator.",
      call. = FALSE
    )
  }
  forward <- if (dense) {
    shift_invert_solver_dense(source_A, sigma)
  } else {
    shift_invert_solver_csc(csc_A, sigma)
  }
  factor <- forward$factor
  # Transposed solves reuse the forward factorisation (no second factor):
  #   dense  M P = Q R          =>  M^T x = b  <=>  x = Q R^{-T} b[pivot]
  #   sparse M[p, q] = L U      =>  M^T x = b  <=>  x[p] = L^{-T} U^{-T} b[q]
  adjoint_solve_fn <- if (dense) {
    R_factor <- NULL
    function(X) {
      X <- as.matrix(X)
      if (is.null(R_factor)) {
        R_factor <<- qr.R(factor)
      }
      Y <- backsolve(R_factor, X[factor$pivot, , drop = FALSE], transpose = TRUE)
      as.matrix(qr.qy(factor, Y))
    }
  } else if (all(c("L", "U", "p", "q") %in% methods::slotNames(factor))) {
    Lt <- NULL
    Ut <- NULL
    function(X) {
      X <- as.matrix(X)
      if (is.null(Lt)) {
        Lt <<- Matrix::t(methods::slot(factor, "L"))
        Ut <<- Matrix::t(methods::slot(factor, "U"))
      }
      p1 <- methods::slot(factor, "p") + 1L
      q1 <- methods::slot(factor, "q") + 1L
      Z <- Matrix::solve(Ut, X[q1, , drop = FALSE])
      Y <- as.matrix(Matrix::solve(Lt, Z))
      out <- matrix(0, nrow(X), ncol(X))
      out[p1, ] <- Y
      out
    }
  } else {
    transposed <- NULL
    function(X) {
      if (is.null(transposed)) {
        transposed <<- shift_invert_solver_csc(
          methods::as(Matrix::t(csc_A), "CsparseMatrix"), sigma
        )
      }
      transposed$solve_fn(X)
    }
  }
  list(
    solve_fn = forward$solve_fn,
    adjoint_solve_fn = adjoint_solve_fn,
    label_kind = paste0(forward$label, "_general"),
    cache = c(forward$cache, list(external_cache = FALSE, generalized = FALSE))
  )
}

#' @keywords internal
#' Nonsymmetric shift-invert (C41): Krylov-Schur Arnoldi on the factorised
#' M = (A - sigma I)^{-1} through the native matrix-free callback, largest
#' magnitude theta of M back-transformed to lambda = sigma + 1 / theta, and the
#' right (optionally left) residual certificate computed on the original A.
solve_shift_invert_general <- function(problem, k, method, tol, vectors,
                                       certify, plan,
                                       left_vectors = "auto") {
  sigma <- method$sigma
  if (length(sigma) != 1L || !is.finite(sigma)) {
    stop("shift_invert(sigma) requires a single finite shift.", call. = FALSE)
  }
  if (is.complex(sigma)) {
    if (Im(sigma) != 0) {
      stop(
        "nonsymmetric shift_invert() currently requires a real sigma; the ",
        "native Krylov-Schur kernel runs in real arithmetic.",
        call. = FALSE
      )
    }
    sigma <- Re(sigma)
  }
  sigma <- as.numeric(sigma)
  if (!is.null(method$factorization)) {
    stop(
      "shift_invert(factorization = ...) is not implemented yet; ",
      "supply shift_invert(solve = ...) for a user-managed factorization cache.",
      call. = FALSE
    )
  }
  if (!is.null(problem$metric)) {
    stop("nonsymmetric generalized shift_invert() is not implemented.", call. = FALSE)
  }
  Aop <- problem$A
  n <- Aop$dim[1L]
  controls <- plan$controls %||% list()
  requested_sigma <- sigma
  solver <- tryCatch(
    shift_invert_general_solver(Aop, sigma, user_solve = method$solve),
    error = function(e) e
  )
  perturbation_reason <- NULL
  if (inherits(solver, "error")) {
    # Only the planner's implicit smallest-magnitude route (sigma = 0 on a
    # possibly singular A) perturbs the shift; an explicit sigma is the
    # caller's choice and its singularity is reported.
    if (!isTRUE(method$perturb_on_singular) ||
        !grepl("singular|rank-deficient|could not factor", conditionMessage(solver))) {
      stop(solver)
    }
    perturbation_reason <- conditionMessage(solver)
    scale <- tryCatch(operator_norm_for_certificate_info(Aop)$value,
                      error = function(e) NA_real_)
    if (!is.finite(scale) || scale <= 0) scale <- 1
    for (offset in scale * c(1e-6, -1e-6, 1e-4, -1e-4, 1e-3, -1e-3)) {
      solver <- tryCatch(
        shift_invert_general_solver(Aop, requested_sigma + offset),
        error = function(e) e
      )
      if (!inherits(solver, "error")) {
        sigma <- requested_sigma + offset
        break
      }
    }
    if (inherits(solver, "error")) {
      stop("shift_invert(sigma = ", requested_sigma, ") failed at the requested ",
           "shift and all perturbation retries: ", conditionMessage(solver),
           call. = FALSE)
    }
  }
  # Ritz vectors of a real nonsymmetric operator can be complex (conjugate
  # pairs); the real factorisations solve the real and imaginary parts as one
  # real block.
  split_complex_solve <- function(f) {
    force(f)
    function(X) {
      if (!is.complex(X)) {
        return(f(X))
      }
      X <- as.matrix(X)
      p <- ncol(X)
      out <- as.matrix(f(cbind(Re(X), Im(X))))
      matrix(complex(real = out[, seq_len(p), drop = FALSE],
                     imaginary = out[, p + seq_len(p), drop = FALSE]),
             nrow(out), p)
    }
  }
  M <- linear_operator(
    dim = c(n, n),
    apply = shift_invert_apply_factory(split_complex_solve(solver$solve_fn)),
    apply_adjoint = if (is.null(solver$adjoint_solve_fn)) {
      NULL
    } else {
      shift_invert_apply_factory(split_complex_solve(solver$adjoint_solve_fn))
    },
    structure = general(),
    name = paste0("shift_invert_", solver$label_kind)
  )
  cert_op <- arnoldi_certificate_operator(Aop)
  subspace <- as.integer(controls$max_subspace %||% native_krylov_schur_default_ncv(n, k))
  ks_maxit <- as.integer(controls$krylov_schur_max_iterations %||%
                           native_krylov_schur_default_maxit())
  inner_attempts <- as.integer(controls$max_restarts %||% 2L)

  attempt <- function(inner_tol, subspace) {
    iter <- native_arnoldi_general(
      M,
      k = k,
      target = largest_magnitude(),
      tol = inner_tol,
      maxit = subspace,
      max_restarts = inner_attempts,
      vectors = TRUE,
      extraction = "projected_ritz",
      krylov_schur_maxit = ks_maxit
    )
    theta <- iter$values
    if (any(Mod(theta) < .Machine$double.eps)) {
      stop(
        "shift_invert(sigma = ", sigma, ") produced a zero-magnitude eigenvalue ",
        "of the inverted operator; sigma is too close to a true eigenvalue. ",
        "Perturb sigma or use a tighter tolerance.",
        call. = FALSE
      )
    }
    lambda <- sigma + 1 / theta
    ord <- order_indices(lambda, problem$target)
    if (length(ord) > k) ord <- ord[seq_len(k)]
    theta <- theta[ord]
    lambda <- lambda[ord]
    vec <- iter$vectors[, ord, drop = FALSE]
    cert <- if (isTRUE(certify) && ncol(vec) > 0L) {
      certify_general_eigen_operator(cert_op, lambda, vec, tol = tol)
    } else {
      empty_certificate(
        tol,
        note = if (!isTRUE(certify)) {
          "shift-invert Arnoldi: certification disabled by caller"
        } else {
          "shift-invert Arnoldi: no eigenpairs returned; residual certificate not computed"
        }
      )
    }
    list(iter = iter, theta = theta, lambda = lambda, vec = vec, cert = cert)
  }

  inner_tol <- tol
  current <- attempt(inner_tol, subspace)
  total_iterations <- as.integer(current$iter$iterations %||% 0L)
  total_matvecs <- as.integer(current$iter$matvecs %||% 0L)
  retries <- 0L
  # Inner convergence on M does not imply original-coordinate convergence
  # (residuals scale by ||A - sigma I|| / |theta|): tighten and enlarge.
  while (isTRUE(certify) && retries < 2L && length(current$cert$converged) &&
         !all(current$cert$converged)) {
    retries <- retries + 1L
    inner_tol <- max(inner_tol * 1e-3, 10 * .Machine$double.eps)
    subspace <- min(n, max(subspace + k, 2L * subspace))
    candidate <- attempt(inner_tol, subspace)
    total_iterations <- total_iterations +
      as.integer(candidate$iter$iterations %||% 0L)
    total_matvecs <- total_matvecs + as.integer(candidate$iter$matvecs %||% 0L)
    if (sum(candidate$cert$converged) >= sum(current$cert$converged)) {
      current <- candidate
    }
  }
  iter <- current$iter
  lambda <- current$lambda
  vec <- current$vec
  cert <- current$cert

  left_contract <- if (identical(left_vectors, "none")) {
    list(supported = FALSE, reason = "not requested (left_vectors = \"none\")")
  } else {
    contract <- arnoldi_left_eigen_contract(
      M, current$theta, vec,
      target = largest_magnitude(),
      tol = inner_tol,
      maxit = subspace,
      max_restarts = inner_attempts,
      extraction = "projected_ritz"
    )
    if (isTRUE(contract$supported)) {
      # Left eigenvectors of M are left eigenvectors of A; certify them against
      # the original operator and the back-transformed eigenvalues.
      contract$certificate <- certify_left_eigen_operator(
        cert_op, lambda, contract$vectors, right_vectors = vec, tol = tol
      )
    }
    contract
  }
  if (identical(left_vectors, "compute") && !isTRUE(left_contract$supported)) {
    stop(
      "left_vectors = \"compute\" but left eigenvectors are unavailable on ",
      plan$method, ": ", left_contract$reason,
      call. = FALSE
    )
  }
  warning_msg <- if (!isTRUE(cert$passed) && isTRUE(certify)) {
    paste0(plan$method, " did not certify all ", k, " requested pairs after ",
           retries, " tightened restart(s); subspace ", subspace)
  } else {
    character()
  }
  if (!is.null(perturbation_reason)) {
    warning_msg <- c(warning_msg, paste0(
      "shift-invert perturbed sigma from ", format(requested_sigma, digits = 17),
      " to ", format(sigma, digits = 17), " after a singular factorization"
    ))
  }
  if (!identical(left_vectors, "none") && !isTRUE(left_contract$supported)) {
    warning_msg <- c(warning_msg,
                     paste0("left eigenvectors unavailable: ", left_contract$reason))
  }
  iter$iterations <- total_iterations
  iter$matvecs <- total_matvecs
  iter$adjoint_block_calls <- left_contract$adjoint_block_calls %||% 0L
  iter$adjoint_columns <- left_contract$adjoint_columns %||% 0L
  restart <- iter$restart
  restart$kind <- "native_krylov_schur_shift_invert_callback"
  restart$max_subspace <- subspace
  restart$inner_tolerance <- inner_tol
  restart$certificate_retries <- retries
  restart$transformed_operator_target <- "largest_magnitude"
  restart$eigenvalue_recovery <- "lambda = sigma + 1 / theta"
  iter$restart <- restart
  make_eigen_result(
    values = lambda,
    vectors = if (isTRUE(vectors)) vec else NULL,
    certificate = cert,
    iter = iter,
    requested = k,
    method_label = plan$method,
    target_label_value = target_label(problem$target),
    plan = plan,
    warnings = warning_msg,
    extras = list(
      sigma = sigma,
      transform = list(
        kind = "shift_invert",
        sigma = sigma,
        requested_sigma = requested_sigma,
        sigma_perturbed = !is.null(perturbation_reason),
        perturbation_reason = perturbation_reason,
        label_kind = solver$label_kind,
        factorization_cache = solver$cache,
        certification = list(
          problem = "original",
          residual_formula = "A * x - lambda * x",
          transformed_residuals_used = FALSE
        )
      ),
      restart = restart,
      right_vectors = if (isTRUE(vectors)) vec else NULL,
      left_vectors = if (isTRUE(left_contract$supported)) left_contract$vectors else NULL,
      left_certificate = if (isTRUE(left_contract$supported)) left_contract$certificate else NULL,
      biorthogonality = if (isTRUE(left_contract$supported)) left_contract$biorthogonality else NULL
    )
  )
}
