# Certified real-double positive-semidefinite factors and image actions.

.psd_factor_fields <- c(
  "schema_version", "dim", "dtype", "representation", "method", "policy",
  "scale", "thresholds", "spectrum", "classification", "rank", "nullity",
  "algebraic_nullity", "rank_bounds", "evidence", "capabilities",
  "operator_identity", "source_identity", "source_semantics",
  "materialization", "certificate", "work", "memory", "serialization",
  "warnings"
)

.psd_capability_names <- c(
  "form", "sqrt", "inverse_sqrt", "pseudoinverse", "image_projector",
  "null_projector", "reduction", "lift", "strict_solve", "gram",
  "orthonormalize", "reduced_operator", "numerical_rank",
  "numerical_nullity", "algebraic_nullity", "serialization", "cache_reuse"
)

.psd_condition_fields <- c(
  "message", "call", "code", "field", "expected", "actual",
  "source_identity", "factor_identity", "representation", "capability",
  "evidence", "scale", "threshold", "defect", "indices", "details"
)

#' Create a scale-aware PSD tolerance
#'
#' @param abs Non-negative finite absolute tolerance.
#' @param rel Non-negative finite relative tolerance.
#' @return An immutable `eigencore_psd_tolerance` record.
#' @export
psd_tolerance <- function(abs = 0, rel = sqrt(.Machine$double.eps)) {
  validate_psd_tolerance_value(abs, "abs")
  validate_psd_tolerance_value(rel, "rel")
  structure(
    list(schema_version = 1L, abs = abs, rel = rel),
    class = "eigencore_psd_tolerance"
  )
}

#' Create a certified PSD classification policy
#'
#' @param symmetry Scale-aware symmetry tolerance.
#' @param positivity Scale-aware negative-eigenvalue acceptance tolerance.
#' @param rank Scale-aware numerical-rank tolerance.
#' @param rhs Scale-aware strict-solve compatibility tolerance.
#' @param scale Matrix scale. Version 1.3 supports only `"frobenius"`.
#' @param symmetry_repair Average an admitted asymmetric source, or reject it.
#' @param negative_repair Clip admitted negative eigenvalues, or reject them.
#' @param structure_repair Canonicalize an admitted structural defect, or
#'   reject it.
#' @return An immutable `eigencore_psd_policy` record.
#' @export
psd_policy <- function(
    symmetry = psd_tolerance(),
    positivity = psd_tolerance(rel = 64 * .Machine$double.eps),
    rank = psd_tolerance(),
    rhs = psd_tolerance(),
    scale = "frobenius",
    symmetry_repair = c("average", "reject"),
    negative_repair = c("clip", "reject"),
    structure_repair = c("canonicalize", "reject")) {
  symmetry <- validate_psd_tolerance(symmetry, "symmetry")
  positivity <- validate_psd_tolerance(positivity, "positivity")
  rank <- validate_psd_tolerance(rank, "rank")
  rhs <- validate_psd_tolerance(rhs, "rhs")
  if (!is.character(scale) || length(scale) != 1L || is.na(scale) ||
      !identical(scale, "frobenius")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "scale",
      "frobenius", scale,
      message = "PSD policy scale must be exactly \"frobenius\" in eigencore 1.3."
    )
  }
  symmetry_repair <- psd_match_enum(
    symmetry_repair, c("average", "reject"), "symmetry_repair"
  )
  negative_repair <- psd_match_enum(
    negative_repair, c("clip", "reject"), "negative_repair"
  )
  structure_repair <- psd_match_enum(
    structure_repair, c("canonicalize", "reject"), "structure_repair"
  )
  structure(
    list(
      schema_version = 1L,
      scale = scale,
      symmetry = symmetry,
      positivity = positivity,
      rank = rank,
      rhs = rhs,
      symmetry_repair = symmetry_repair,
      negative_repair = negative_repair,
      structure_repair = structure_repair
    ),
    class = "eigencore_psd_policy"
  )
}

#' Construct a certified identity PSD factor
#'
#' @param dim Non-negative whole dimension.
#' @param policy A PSD policy from [psd_policy()].
#' @return An `eigencore_psd_factor`.
#' @export
psd_identity <- function(dim, policy = psd_policy()) {
  started <- proc.time()[["elapsed"]]
  policy <- validate_psd_policy(policy)
  n <- psd_dimension(dim, "dim")
  values <- rep(1, n)
  source_identity <- psd_builtin_identity(
    list(kind = "identity", dim = n), c(n, n), "hermitian"
  )
  psd_construct_complete_factor(
    original_values = values,
    vectors = NULL,
    source_index = seq_len(n),
    representation = "identity",
    method = "analytic identity PSD factor",
    source_identity = source_identity,
    policy = policy,
    state_kind = "identity",
    source_dimnames = list(NULL, NULL),
    symmetry_defect = 0,
    original_source = NULL,
    algebraic_nullity = 0L,
    materialization = list(
      source = "identity", factor = "analytic_identity",
      dense_n_by_n = FALSE, notes = character()
    ),
    started = started,
    scale_override = sqrt(as.double(n))
  )
}

#' Construct a certified real-double PSD factor
#'
#' A double vector is interpreted as a diagonal. Base double matrices and the
#' admitted dense/diagonal Matrix classes are validated as complete sources.
#'
#' @param x A supported real-double diagonal or square matrix source.
#' @param policy A PSD policy from [psd_policy()].
#' @param source Source ownership, currently immutable snapshots for complete
#'   identity/diagonal/dense paths.
#' @return An `eigencore_psd_factor`.
#' @export
psd_factor <- function(x, policy = psd_policy(), source = c("snapshot", "live")) {
  started <- proc.time()[["elapsed"]]
  policy <- validate_psd_policy(policy)
  source <- psd_match_enum(source, c("snapshot", "live"), "source")

  if (inherits(x, "eigencore_psd_factor")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_structure", "x",
      "an uncertified source", class(x),
      message = "psd_factor() does not silently reclassify an existing certified factor."
    )
  }
  if (inherits(x, "eigencore_operator")) {
    source_semantics <- if (identical(source, "live")) {
      "live_versioned_source"
    } else {
      "uncacheable_opaque_source"
    }
    psd_abort(
      "eigencore_psd_incomplete_evidence", "incomplete_evidence", "x",
      "an admitted complete or structural PSD provider", class(x),
      source_identity = operator_identity(x),
      representation = "opaque_operator",
      capability = "construction",
      evidence = psd_unavailable_evidence(source_semantics),
      message = paste(
        "A versioned operator identity is not PSD or complete-spectrum evidence;",
        "opaque callbacks are never probed into a certified PSD factor."
      )
    )
  }
  if (identical(source, "live")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "source",
      "snapshot for matrix-backed sources", source,
      message = "source = \"live\" requires an admitted versioned operator provider."
    )
  }

  if (is.double(x) && is.atomic(x) && is.null(dim(x))) {
    if (any(!is.finite(x))) {
      psd_nonfinite_input(x, "x")
    }
    return(psd_construct_diagonal(
      values = x,
      names = names(x),
      policy = policy,
      started = started,
      source_class = class(x)
    ))
  }

  if (inherits(x, "ddiMatrix")) {
    values <- if (identical(methods::slot(x, "diag"), "U")) {
      rep(1, nrow(x))
    } else {
      as.numeric(methods::slot(x, "x"))
    }
    if (any(!is.finite(values))) {
      psd_nonfinite_input(values, "x")
    }
    dn <- dimnames(x)
    axis_names <- dn[[1L]] %||% dn[[2L]]
    return(psd_construct_diagonal(
      values = values,
      names = axis_names,
      policy = policy,
      started = started,
      source_class = class(x)
    ))
  }

  if (inherits(x, "Matrix")) {
    admitted <- vapply(
      c("dgeMatrix", "dsyMatrix", "dpoMatrix"),
      function(cl) methods::is(x, cl),
      logical(1L)
    )
    if (!any(admitted)) {
      psd_abort(
        "eigencore_psd_incomplete_evidence", "incomplete_evidence", "x",
        "ddiMatrix, dgeMatrix, dsyMatrix, or dpoMatrix", class(x),
        representation = "generic_sparse",
        evidence = psd_unavailable_evidence("immutable_snapshot"),
        message = paste(
          "Generic sparse Matrix storage is not complete PSD evidence.",
          "Use an explicitly admitted structural constructor."
        )
      )
    }
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dtype", "x",
      "a real-double vector or admitted square matrix", class(x),
      message = "psd_factor() requires a supported real-double source."
    )
  }
  if (!is.double(x) || is.complex(x)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dtype", "x",
      "a real-double matrix", typeof(x),
      message = "PSD factors in eigencore 1.3 require a real-double matrix without implicit coercion."
    )
  }
  if (length(dim(x)) != 2L || nrow(x) != ncol(x)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", "x",
      "a square matrix", dim(x),
      message = "A PSD source must be square."
    )
  }
  if (any(!is.finite(x))) {
    psd_nonfinite_input(x, "x")
  }
  psd_validate_square_dimnames(x)
  psd_construct_dense(x, policy = policy, started = started)
}

#' Construct a supplied Gram PSD factor
#'
#' Dense base and `dgeMatrix` inputs receive a complete compact-SVD
#' certificate and canonical spectral actions. A `dgCMatrix` remains sparse
#' and exposes only the form and Gram actions justified by `K = L L^T` or
#' `K = L^T L`; sparse storage alone does not certify rank or a principal
#' square root.
#'
#' @param x A finite real-double base matrix, `dgeMatrix`, or `dgCMatrix`.
#' @param orientation Whether columns define `K = x x^T` or rows define
#'   `K = x^T x`.
#' @param policy A PSD policy.
#' @param source Source ownership.
#' @return A certified PSD factor when the representation is admitted.
#' @export
psd_gram_factor <- function(x, orientation = c("columns", "rows"),
                            policy = psd_policy(),
                            source = c("snapshot", "live")) {
  started <- proc.time()[["elapsed"]]
  orientation <- psd_match_enum(orientation, c("columns", "rows"), "orientation")
  source <- psd_match_enum(source, c("snapshot", "live"), "source")
  policy <- validate_psd_policy(policy)
  if (identical(source, "live")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "source",
      "snapshot for the admitted matrix Gram providers", source,
      representation = if (inherits(x, "Matrix")) "gram_sparse" else "gram_dense",
      capability = "construction",
      message = "source = \"live\" requires a separately admitted versioned Gram provider."
    )
  }
  if (is.matrix(x)) {
    return(psd_construct_dense_gram(
      x, orientation = orientation, policy = policy, started = started,
      source_class = class(x)
    ))
  }
  if (methods::is(x, "dgeMatrix")) {
    return(psd_construct_dense_gram(
      as.matrix(x), orientation = orientation, policy = policy,
      started = started, source_class = class(x)
    ))
  }
  if (methods::is(x, "dgCMatrix")) {
    return(psd_construct_sparse_gram(
      x, orientation = orientation, policy = policy, started = started
    ))
  }
  psd_abort(
    "eigencore_psd_invalid_input", "invalid_structure", "x",
    "a real-double base matrix, dgeMatrix, or dgCMatrix", class(x),
    representation = if (inherits(x, "Matrix")) "gram_sparse" else "gram_dense",
    capability = "construction",
    details = list(orientation = orientation, source = source),
    message = "psd_gram_factor() accepts only the explicitly admitted 1.3 Gram classes."
  )
}

#' Construct a structural sparse graph-Laplacian PSD factor
#'
#' The constructor admits only CSC `dgCMatrix` and `dsCMatrix` inputs. It
#' validates symmetry, non-positive off-diagonals, and zero row sums without
#' dense conversion. Connected components certify algebraic nullity; numerical
#' rank and spectral actions remain unavailable without separate spectral
#' evidence.
#'
#' @param x A finite real-double `dgCMatrix` or `dsCMatrix` graph Laplacian.
#' @param policy A PSD policy.
#' @param source Source ownership.
#' @return A certified PSD factor when the representation is admitted.
#' @export
psd_laplacian <- function(x, policy = psd_policy(),
                          source = c("snapshot", "live")) {
  started <- proc.time()[["elapsed"]]
  source <- psd_match_enum(source, c("snapshot", "live"), "source")
  policy <- validate_psd_policy(policy)
  if (identical(source, "live")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "source",
      "snapshot for the admitted sparse Laplacian providers", source,
      representation = "laplacian_sparse",
      capability = "construction",
      message = "source = \"live\" requires a separately admitted versioned Laplacian provider."
    )
  }
  if (!methods::is(x, "dgCMatrix") && !methods::is(x, "dsCMatrix")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_structure", "x",
      "a real-double dgCMatrix or dsCMatrix", class(x),
      representation = "laplacian_sparse",
      capability = "construction",
      message = "psd_laplacian() accepts only the explicitly admitted CSC classes."
    )
  }
  psd_construct_sparse_laplacian(x, policy = policy, started = started)
}

#' Inspect certified PSD capabilities
#'
#' @param x A certified PSD factor.
#' @return An `eigencore_psd_capabilities` record.
#' @export
psd_capabilities <- function(x) {
  x <- validate_psd_factor(x)
  deep_copy_record(x$capabilities)
}

#' Return a complete certified PSD spectrum
#'
#' @param x A certified PSD factor.
#' @param repaired Return the repaired/action spectrum instead of the admitted
#'   source spectrum.
#' @return A nonincreasing numeric spectrum.
#' @export
psd_spectrum <- function(x, repaired = FALSE) {
  x <- validate_psd_factor(x)
  if (!is.logical(repaired) || length(repaired) != 1L || is.na(repaired)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "repaired",
      "one TRUE or FALSE value", repaired,
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation
    )
  }
  if (!identical(x$spectrum$coverage, "complete") ||
      is.null(x$spectrum$original) || is.null(x$spectrum$repaired)) {
    psd_incomplete_evidence(x, "spectrum")
  }
  value <- if (isTRUE(repaired)) x$spectrum$repaired else x$spectrum$original
  as.numeric(value)
}

#' Return certified tolerance-relative PSD rank
#'
#' @param x A certified PSD factor.
#' @return One integer numerical rank.
#' @export
psd_rank <- function(x) {
  x <- validate_psd_factor(x)
  if (is.na(x$rank) || !isTRUE(x$capabilities$numerical_rank$available)) {
    psd_incomplete_evidence(x, "numerical_rank")
  }
  x$rank
}

#' Return certified numerical or algebraic nullity
#'
#' @param x A certified PSD factor.
#' @param type Tolerance-relative numerical nullity or structurally proven
#'   algebraic nullity.
#' @return One integer nullity.
#' @export
psd_nullity <- function(x, type = c("numerical", "algebraic")) {
  x <- validate_psd_factor(x)
  type <- psd_match_enum(type, c("numerical", "algebraic"), "type")
  field <- if (identical(type, "numerical")) "nullity" else "algebraic_nullity"
  capability <- paste0(type, "_nullity")
  if (is.na(x[[field]]) || !isTRUE(x$capabilities[[capability]]$available)) {
    psd_incomplete_evidence(x, capability)
  }
  x[[field]]
}

#' Apply a certified PSD action to a vector or block
#'
#' @param x A certified PSD factor.
#' @param X A finite real-double vector or matrix.
#' @param action Certified action to apply.
#' @return A vector or matrix matching the input block shape.
#' @export
psd_apply <- function(
    x, X,
    action = c("form", "sqrt", "inverse_sqrt", "pseudoinverse",
               "image_projector", "null_projector")) {
  x <- validate_psd_factor(x)
  action <- psd_match_enum(
    action,
    c("form", "sqrt", "inverse_sqrt", "pseudoinverse",
      "image_projector", "null_projector"),
    "action"
  )
  psd_require_capability(x, action)
  block <- psd_prepare_block(X, x$dim[[2L]], "X", x)
  out <- psd_apply_matrix(x, block$value, action)
  psd_restore_block_shape(out, block, x)
}

#' Expose a certified PSD action as an eigencore operator
#'
#' @param x A certified PSD factor.
#' @param action Certified square action.
#' @return An `eigencore_operator`.
#' @export
psd_operator <- function(
    x,
    action = c("form", "sqrt", "inverse_sqrt", "pseudoinverse",
               "image_projector", "null_projector")) {
  x <- validate_psd_factor(x)
  action <- psd_match_enum(
    action,
    c("form", "sqrt", "inverse_sqrt", "pseudoinverse",
      "image_projector", "null_projector"),
    "action"
  )
  psd_require_capability(x, action)
  frozen <- deep_copy_record(x)
  action_identity <- psd_action_identity(frozen, action)
  apply <- function(X, alpha = 1, beta = 0, Y = NULL) {
    value <- psd_apply(frozen, X, action = action)
    out <- alpha * as.matrix(value)
    if (!is.null(Y) && beta != 0) {
      out <- out + beta * Y
    }
    out
  }
  op <- linear_operator(
    dim = frozen$dim,
    apply = apply,
    apply_adjoint = apply,
    dtype = "double",
    structure = hermitian(),
    name = paste("certified PSD", action),
    metadata = list(
      parent = frozen$operator_identity,
      psd_action = action,
      frobenius_norm = psd_action_frobenius_norm(frozen, action)
    ),
    operator_id = action_identity$operator_id,
    revision = action_identity$revision,
    portable = action_identity$portable
  )
  op$identity <- action_identity
  op
}

#' Reduce an original-coordinate block to canonical PSD image coordinates
#'
#' @param x A certified PSD factor.
#' @param X An original-coordinate vector or matrix.
#' @return Canonical reduced coordinates.
#' @export
psd_reduce <- function(x, X) {
  x <- validate_psd_factor(x)
  psd_require_capability(x, "reduction")
  block <- psd_prepare_block(X, x$dim[[2L]], "X", x)
  out <- psd_reduce_matrix(x, block$value)
  psd_restore_reduced_shape(out, block)
}

#' Lift canonical PSD image coordinates
#'
#' @param x A certified PSD factor.
#' @param Z A reduced-coordinate vector or matrix.
#' @return Original-coordinate minimum-norm representatives.
#' @export
psd_lift <- function(x, Z) {
  x <- validate_psd_factor(x)
  psd_require_capability(x, "lift")
  block <- psd_prepare_reduced_block(Z, x$rank, "Z", x)
  out <- psd_lift_matrix(x, block$value)
  psd_restore_block_shape(out, block, x)
}

#' Strictly solve a compatible singular PSD system
#'
#' @param x A certified PSD factor.
#' @param B A right-hand-side vector or matrix.
#' @param tolerance Optional typed RHS compatibility override.
#' @return An `eigencore_psd_solve_result`.
#' @export
psd_solve <- function(x, B, tolerance = NULL) {
  started <- proc.time()[["elapsed"]]
  x <- validate_psd_factor(x)
  psd_require_capability(x, "strict_solve")
  tolerance <- if (is.null(tolerance)) {
    x$policy$rhs
  } else {
    validate_psd_tolerance(tolerance, "tolerance")
  }
  block <- psd_prepare_block(B, x$dim[[2L]], "B", x)
  null_part <- psd_apply_matrix(x, block$value, "null_projector")
  defect <- psd_column_norms(null_part)
  rhs_norm <- psd_column_norms(block$value)
  threshold <- tolerance$abs + tolerance$rel * rhs_norm
  compatible <- defect <= threshold
  if (any(!compatible)) {
    bad <- which(!compatible)
    psd_abort(
      "eigencore_psd_incompatible_rhs", "incompatible_rhs", "B",
      "every column in image(K)", defect,
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation,
      capability = "strict_solve",
      evidence = x$evidence,
      scale = rhs_norm,
      threshold = threshold,
      defect = defect,
      indices = bad,
      details = list(
        offending_columns = bad,
        compatibility_defect = defect,
        compatibility_threshold = threshold
      ),
      message = paste0(
        "Strict PSD solve requires every RHS column to lie in image(K); ",
        "incompatible column(s): ", paste(bad, collapse = ", "), "."
      )
    )
  }
  solution_matrix <- psd_apply_matrix(x, block$value, "pseudoinverse")
  residual <- psd_apply_matrix(x, solution_matrix, "form") - block$value
  equation_defect <- psd_column_norms(residual)
  solution_norm <- psd_column_norms(solution_matrix)
  equation_threshold <- threshold + 64 * .Machine$double.eps *
    (x$scale * solution_norm + rhs_norm)
  elapsed <- proc.time()[["elapsed"]] - started
  work <- psd_work_record(solve_seconds = elapsed, total_seconds = elapsed)
  solution <- psd_restore_block_shape(solution_matrix, block, x)
  certificate <- psd_operation_certificate(
    type = "psd_strict_solve",
    x = x,
    thresholds = list(
      compatibility = threshold,
      equation = equation_threshold
    ),
    residuals = list(compatibility = defect, equation = equation_defect),
    passed = compatible & equation_defect <= equation_threshold,
    action_bounds = list(
      compatibility = threshold,
      equation = equation_threshold
    ),
    notes = "strict solve certified per RHS column"
  )
  structure(
    list(
      schema_version = 1L,
      solution = solution,
      compatibility_defect = defect,
      compatibility_threshold = threshold,
      compatible = compatible,
      factor_identity = deep_copy_record(x$operator_identity),
      certificate = certificate,
      work = work,
      warnings = character()
    ),
    class = "eigencore_psd_solve_result"
  )
}

#' Compute a PSD Gram or cross-Gram block
#'
#' @param x A certified PSD factor.
#' @param X A finite real-double vector or matrix.
#' @param Y Optional second block.
#' @return `crossprod(X, K Y)` as a dense matrix.
#' @export
psd_gram <- function(x, X, Y = NULL) {
  x <- validate_psd_factor(x)
  psd_require_capability(x, "gram")
  left <- psd_prepare_block(X, x$dim[[2L]], "X", x)
  right <- if (is.null(Y)) left else psd_prepare_block(Y, x$dim[[2L]], "Y", x)
  crossprod(left$value, psd_apply_matrix(x, right$value, "form"))
}

#' PSD-orthonormalize a block modulo the certified null space
#'
#' @param x A certified PSD factor.
#' @param X A finite real-double vector or matrix.
#' @param required_rank Optional exact output width.
#' @param tolerance Optional typed Gram-rank tolerance.
#' @return An `eigencore_psd_block_result`.
#' @export
psd_orthonormalize <- function(x, X, required_rank = NULL, tolerance = NULL) {
  started <- proc.time()[["elapsed"]]
  x <- validate_psd_factor(x)
  psd_require_capability(x, "orthonormalize")
  block <- psd_prepare_block(X, x$dim[[2L]], "X", x)
  tolerance <- if (is.null(tolerance)) {
    x$policy$rank
  } else {
    validate_psd_tolerance(tolerance, "tolerance")
  }
  if (!is.null(required_rank)) {
    required_rank <- psd_dimension(required_rank, "required_rank")
    if (required_rank > ncol(block$value)) {
      psd_abort(
        "eigencore_psd_infeasible_block", "infeasible_block_rank",
        "required_rank", paste0("at most ", ncol(block$value)), required_rank,
        source_identity = x$source_identity,
        factor_identity = x$operator_identity,
        representation = x$representation,
        capability = "orthonormalize",
        evidence = x$evidence
      )
    }
  }
  gram <- psd_gram(x, block$value)
  gram <- (gram + t(gram)) / 2
  m <- ncol(gram)
  eig <- if (m == 0L) {
    list(values = numeric(), vectors = matrix(numeric(), 0L, 0L))
  } else {
    native_dense_symmetric_eigen(gram)
  }
  ord <- order(eig$values, decreasing = TRUE)
  vals <- eig$values[ord]
  vecs <- eig$vectors[, ord, drop = FALSE]
  scale <- psd_frobenius_norm(gram)
  threshold <- tolerance$abs + tolerance$rel * scale
  negative_threshold <- x$policy$positivity$abs +
    x$policy$positivity$rel * scale
  materially_negative <- which(vals < -negative_threshold)
  if (length(materially_negative)) {
    psd_abort(
      "eigencore_psd_infeasible_block", "infeasible_block_rank", "gram",
      "a PSD block Gram", vals[materially_negative],
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation,
      capability = "orthonormalize",
      evidence = x$evidence,
      scale = scale,
      threshold = negative_threshold,
      defect = min(vals),
      indices = materially_negative,
      message = "The PSD block Gram is materially indefinite relative to its certified scale."
    )
  }
  retained <- which(vals > threshold)
  discovered_rank <- length(retained)
  requested <- required_rank
  if (!is.null(requested) && discovered_rank < requested) {
    psd_abort(
      "eigencore_psd_infeasible_block", "infeasible_block_rank",
      "required_rank", paste0("at most ", discovered_rank), requested,
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation,
      capability = "orthonormalize",
      evidence = x$evidence,
      scale = scale,
      threshold = threshold,
      defect = vals,
      details = list(discovered_rank = discovered_rank, eigenvalues = vals)
    )
  }
  use_rank <- if (is.null(requested)) discovered_rank else requested
  use <- if (use_rank == 0L) integer() else retained[seq_len(use_rank)]
  basis <- if (use_rank == 0L) {
    matrix(numeric(), nrow(block$value), 0L)
  } else {
    block$value %*% vecs[, use, drop = FALSE] %*%
      diag(1 / sqrt(vals[use]), nrow = use_rank)
  }
  rownames(basis) <- rownames(block$value)
  post <- if (use_rank == 0L) {
    0
  } else {
    psd_frobenius_norm(psd_gram(x, basis) - diag(use_rank))
  }
  post_threshold <- sqrt(.Machine$double.eps) * max(1, sqrt(use_rank))
  condition <- if (use_rank == 0L) NA_real_ else max(vals[use]) / min(vals[use])
  dropped <- setdiff(seq_len(m), use)
  elapsed <- proc.time()[["elapsed"]] - started
  work <- psd_work_record(solve_seconds = elapsed, total_seconds = elapsed)
  cert <- psd_operation_certificate(
    type = "psd_block_orthonormalization",
    x = x,
    thresholds = list(
      gram_rank = threshold,
      postcondition = post_threshold
    ),
    residuals = list(postcondition = post),
    passed = post <= post_threshold,
    action_bounds = list(postcondition = post_threshold),
    notes = "PSD block orthonormalization postcondition"
  )
  structure(
    list(
      schema_version = 1L,
      basis = basis,
      rank = as.integer(use_rank),
      required_rank = if (is.null(requested)) NULL else as.integer(requested),
      condition = condition,
      dropped = dropped,
      gram = gram,
      postcondition_error = post,
      factor_identity = deep_copy_record(x$operator_identity),
      certificate = cert,
      work = work,
      warnings = character()
    ),
    class = "eigencore_psd_block_result"
  )
}

#' Construct the Euclidean operator induced on a certified PSD image
#'
#' @param x A certified PSD factor.
#' @param A A compatible real-double eigencore operator or matrix.
#' @return An `eigencore_operator` for `R+^T A R+`.
#' @export
psd_reduced_operator <- function(x, A) {
  x <- validate_psd_factor(x)
  psd_require_capability(x, "reduced_operator")
  Aop <- as_operator(A)
  if (!identical(Aop$dtype, "double") ||
      !identical(as.integer(Aop$dim), as.integer(x$dim))) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", "A",
      list(dim = x$dim, dtype = "double"),
      list(dim = Aop$dim, dtype = Aop$dtype),
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation,
      capability = "reduced_operator"
    )
  }
  frozen <- deep_copy_record(x)
  parent_identity <- operator_identity(Aop)
  payload <- list(
    kind = "psd_reduced_operator", factor = frozen$operator_identity,
    parent = parent_identity
  )
  operator_id <- paste0("psd-reduced-", stable_raw_hash(list(
    frozen$operator_identity$operator_id, parent_identity$operator_id
  )))
  revision <- stable_raw_hash(payload)
  portable <- isTRUE(frozen$operator_identity$portable) && isTRUE(parent_identity$portable)
  apply <- function(Z, alpha = 1, beta = 0, Y = NULL) {
    Z <- as.matrix(Z)
    lifted <- psd_lift_matrix(frozen, Z)
    applied <- Aop$apply(lifted)
    out <- alpha * psd_inverse_reduce_matrix(frozen, applied)
    if (!is.null(Y) && beta != 0) out <- out + beta * Y
    out
  }
  apply_adjoint <- if (is.null(Aop$apply_adjoint)) {
    NULL
  } else {
    function(Z, alpha = 1, beta = 0, Y = NULL) {
      Z <- as.matrix(Z)
      lifted <- psd_lift_matrix(frozen, Z)
      applied <- Aop$apply_adjoint(lifted)
      out <- alpha * psd_inverse_reduce_matrix(frozen, applied)
      if (!is.null(Y) && beta != 0) out <- out + beta * Y
      out
    }
  }
  linear_operator(
    dim = c(frozen$rank, frozen$rank),
    apply = apply,
    apply_adjoint = apply_adjoint,
    dtype = "double",
    structure = Aop$structure,
    name = "PSD image-reduced operator",
    metadata = list(parent = parent_identity, psd_factor = frozen$operator_identity),
    operator_id = operator_id,
    revision = revision,
    portable = portable
  )
}

#' @export
as_operator.eigencore_psd_factor <- function(x, ...) {
  psd_operator(x, action = "form")
}

#' @export
print.eigencore_psd_policy <- function(x, ...) {
  cat("eigencore PSD policy\n")
  cat("  scale:", x$scale, "\n")
  cat("  symmetry: abs", format(x$symmetry$abs), "rel", format(x$symmetry$rel), "\n")
  cat("  positivity: abs", format(x$positivity$abs), "rel", format(x$positivity$rel), "\n")
  cat("  rank: abs", format(x$rank$abs), "rel", format(x$rank$rel), "\n")
  cat("  rhs: abs", format(x$rhs$abs), "rel", format(x$rhs$rel), "\n")
  invisible(x)
}

#' @export
print.eigencore_psd_factor <- function(x, ...) {
  x <- validate_psd_factor(x)
  cat("eigencore certified PSD factor\n")
  cat("  representation:", x$representation, "\n")
  cat("  dimension:", x$dim[[1L]], "\n")
  cat("  numerical rank:", x$rank, "\n")
  cat("  numerical nullity:", x$nullity, "\n")
  cat("  fidelity:", x$evidence$action_fidelity, "\n")
  cat("  certificate passed:", x$certificate$passed, "\n")
  invisible(x)
}

#' @export
print.eigencore_psd_capabilities <- function(x, ...) {
  cat("eigencore PSD capabilities (", x$representation, ")\n", sep = "")
  available <- vapply(x[.psd_capability_names], function(y) isTRUE(y$available), logical(1L))
  cat("  available:", paste(names(available)[available], collapse = ", "), "\n")
  unavailable <- names(available)[!available]
  if (length(unavailable)) cat("  unavailable:", paste(unavailable, collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.eigencore_psd_certificate <- function(x, ...) {
  cat("eigencore PSD certificate\n")
  cat("  scope:", x$scope, "\n")
  cat("  passed:", x$passed, "\n")
  cat("  representation:", x$representation, "\n")
  cat("  fidelity:", x$evidence$action_fidelity, "\n")
  cat("  symmetry defect:", format(x$symmetry_defect), "\n")
  cat("  repair defect:", format(x$repair_defect), "\n")
  cat("  source/action defect:", format(x$source_action_defect), "\n")
  invisible(x)
}

#' @export
print.eigencore_psd_solve_result <- function(x, ...) {
  cat("eigencore strict PSD solve\n")
  cat("  columns:", length(x$compatible), "\n")
  cat("  compatible:", all(x$compatible), "\n")
  cat("  certificate passed:", x$certificate$passed, "\n")
  invisible(x)
}

#' @export
print.eigencore_psd_block_result <- function(x, ...) {
  cat("eigencore PSD-orthonormal block\n")
  cat("  rank:", x$rank, "\n")
  cat("  condition:", format(x$condition), "\n")
  cat("  postcondition error:", format(x$postcondition_error), "\n")
  invisible(x)
}

# Internal validation and construction -------------------------------------

validate_psd_tolerance_value <- function(value, field) {
  if (!is.double(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", field,
      "one finite non-negative double", value,
      message = paste0(field, " must be one finite non-negative double.")
    )
  }
  invisible(value)
}

validate_psd_tolerance <- function(x, field = "tolerance") {
  if (!inherits(x, "eigencore_psd_tolerance") ||
      !identical(names(x), c("schema_version", "abs", "rel")) ||
      !identical(x$schema_version, 1L)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", field,
      "eigencore_psd_tolerance schema version 1", class(x)
    )
  }
  validate_psd_tolerance_value(x$abs, paste0(field, "$abs"))
  validate_psd_tolerance_value(x$rel, paste0(field, "$rel"))
  deep_copy_record(x)
}

validate_psd_policy <- function(x) {
  required <- c(
    "schema_version", "scale", "symmetry", "positivity", "rank", "rhs",
    "symmetry_repair", "negative_repair", "structure_repair"
  )
  if (!inherits(x, "eigencore_psd_policy") || !identical(names(x), required) ||
      !identical(x$schema_version, 1L)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "policy",
      "eigencore_psd_policy schema version 1", class(x)
    )
  }
  if (!identical(x$scale, "frobenius")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", "policy$scale",
      "frobenius", x$scale
    )
  }
  for (field in c("symmetry", "positivity", "rank", "rhs")) {
    validate_psd_tolerance(x[[field]], paste0("policy$", field))
  }
  psd_match_enum(x$symmetry_repair, c("average", "reject"), "policy$symmetry_repair")
  psd_match_enum(x$negative_repair, c("clip", "reject"), "policy$negative_repair")
  psd_match_enum(x$structure_repair, c("canonicalize", "reject"), "policy$structure_repair")
  deep_copy_record(x)
}

psd_match_enum <- function(value, choices, field) {
  if (!is.character(value) || !length(value) || anyNA(value)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", field,
      choices, value
    )
  }
  choice <- tryCatch(match.arg(value, choices), error = identity)
  if (inherits(choice, "error")) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_policy", field,
      choices, value,
      message = conditionMessage(choice)
    )
  }
  choice
}

psd_dimension <- function(x, field) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 0 || x != floor(x) || x > .Machine$integer.max) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", field,
      "one non-negative whole dimension", x
    )
  }
  as.integer(x)
}

psd_abort <- function(subclass, code, field = NULL, expected = NULL,
                      actual = NULL, source_identity = NULL,
                      factor_identity = NULL, representation = NULL,
                      capability = NULL, evidence = NULL, scale = NULL,
                      threshold = NULL, defect = NULL, indices = integer(),
                      details = list(), message = NULL) {
  if (is.null(message)) {
    message <- paste0(
      "Eigencore PSD failure (", code, ")",
      if (is.null(field)) "." else paste0(": ", field, ".")
    )
  }
  values <- list(
    message = message,
    call = NULL,
    code = code,
    field = field,
    expected = expected,
    actual = actual,
    source_identity = source_identity,
    factor_identity = factor_identity,
    representation = representation,
    capability = capability,
    evidence = evidence,
    scale = scale,
    threshold = threshold,
    defect = defect,
    indices = as.integer(indices),
    details = details
  )
  values <- values[.psd_condition_fields]
  stop(structure(
    values,
    class = c(subclass, "eigencore_psd_error", "error", "condition")
  ))
}

psd_nonfinite_input <- function(x, field, factor = NULL) {
  bad <- which(!is.finite(x))
  psd_abort(
    "eigencore_psd_invalid_input", "nonfinite_input", field,
    "finite real-double values", x[bad], indices = bad,
    source_identity = if (is.null(factor)) NULL else factor$source_identity,
    factor_identity = if (is.null(factor)) NULL else factor$operator_identity,
    representation = if (is.null(factor)) NULL else factor$representation,
    evidence = if (is.null(factor)) NULL else factor$evidence,
    message = paste0("PSD input ", field, " contains non-finite values.")
  )
}

psd_validate_square_dimnames <- function(x) {
  dn <- dimnames(x)
  if (!is.null(dn[[1L]]) && !is.null(dn[[2L]]) &&
      !identical(dn[[1L]], dn[[2L]])) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", "dimnames",
      "matching row and column names", dn,
      message = "PSD source row and column names must agree when both are present."
    )
  }
  invisible(TRUE)
}

psd_frobenius_norm <- function(x) {
  if (!length(x)) return(0)
  magnitude <- max(abs(x))
  if (!is.finite(magnitude)) return(Inf)
  if (magnitude == 0) return(0)
  value <- magnitude * sqrt(sum((x / magnitude)^2))
  if (is.finite(value)) value else Inf
}

psd_symmetry_defect <- function(A) {
  magnitude <- max(abs(A))
  if (!length(A) || magnitude == 0) return(0)
  scaled <- A / magnitude
  magnitude * psd_frobenius_norm(scaled - t(scaled))
}

psd_sparse_frobenius <- function(x) {
  if (!length(x) || any(dim(x) == 0L)) return(0)
  magnitude <- as.numeric(max(abs(x)))
  if (!is.finite(magnitude)) return(Inf)
  if (magnitude == 0) return(0)
  scaled <- x / magnitude
  value <- magnitude * sqrt(as.numeric(sum(scaled * scaled)))
  if (is.finite(value)) value else Inf
}

psd_sparse_gram_frobenius <- function(x, orientation) {
  small_gram <- if (identical(orientation, "columns")) {
    Matrix::crossprod(x)
  } else {
    Matrix::tcrossprod(x)
  }
  psd_sparse_frobenius(small_gram)
}

psd_sparse_components <- function(laplacian) {
  n <- nrow(laplacian)
  entries <- Matrix::summary(laplacian)
  edges <- which(entries$i != entries$j & entries$x < 0)
  if (n == 0L) {
    return(list(count = 0L, membership = integer(), edge_count = 0L))
  }
  parent <- seq_len(n)
  find_root <- function(index) {
    while (parent[[index]] != index) index <- parent[[index]]
    index
  }
  if (length(edges)) {
    for (edge in edges) {
      left <- find_root(entries$i[[edge]])
      right <- find_root(entries$j[[edge]])
      if (left != right) parent[[right]] <- left
    }
  }
  roots <- vapply(seq_len(n), find_root, integer(1L))
  unique_roots <- unique(roots)
  membership <- match(roots, unique_roots)
  list(
    count = as.integer(length(unique_roots)),
    membership = as.integer(membership),
    edge_count = as.integer(length(edges))
  )
}

psd_builtin_identity <- function(payload, dim, structure) {
  digest <- stable_raw_hash(list(schema_version = 1L, payload = payload))
  new_operator_identity(
    operator_id = paste0("builtin-", digest),
    revision = digest,
    origin = "builtin",
    dim = dim,
    dtype = "double",
    structure = structure,
    portable = TRUE
  )
}

psd_factor_identity <- function(source_identity, policy, representation,
                                state_token, dim) {
  lineage <- stable_raw_hash(list(
    schema_version = 1L,
    kind = "certified_psd_factor",
    source_operator_id = source_identity$operator_id,
    representation = representation,
    policy = unclass(policy)
  ))
  revision <- stable_raw_hash(list(
    schema_version = 1L,
    source_revision = source_identity$revision,
    state_token = state_token
  ))
  new_operator_identity(
    operator_id = paste0("psd-", lineage),
    revision = revision,
    origin = "composite",
    dim = dim,
    dtype = "double",
    structure = "hermitian",
    portable = isTRUE(source_identity$portable)
  )
}

psd_construct_diagonal <- function(values, names, policy, started,
                                   source_class) {
  n <- length(values)
  source_identity <- psd_builtin_identity(
    list(kind = "diagonal", values = values, names = names,
         source_class = source_class),
    c(n, n), "hermitian"
  )
  original <- values
  ord <- order(original, decreasing = TRUE)
  psd_construct_complete_factor(
    original_values = original[ord],
    vectors = NULL,
    source_index = ord,
    representation = "diagonal",
    method = "analytic diagonal PSD factor",
    source_identity = source_identity,
    policy = policy,
    state_kind = "diagonal",
    source_dimnames = list(names, names),
    symmetry_defect = 0,
    original_source = original,
    algebraic_nullity = as.integer(sum(original == 0)),
    materialization = list(
      source = "diagonal_vector", factor = "analytic_diagonal",
      dense_n_by_n = FALSE, notes = character()
    ),
    started = started,
    scale_override = psd_frobenius_norm(original)
  )
}

psd_construct_dense <- function(A, policy, started) {
  n <- nrow(A)
  defect <- psd_symmetry_defect(A)
  S <- 0.5 * A + 0.5 * t(A)
  scale <- psd_frobenius_norm(S)
  if (!is.finite(scale) || !is.finite(defect)) {
    psd_abort(
      "eigencore_psd_invalid_input", "nonfinite_input", "scale",
      "a finite Frobenius scale and symmetry defect",
      list(scale = scale, symmetry_defect = defect),
      message = "The finite PSD input overflows its required Frobenius diagnostics."
    )
  }
  tau_sym <- policy$symmetry$abs + policy$symmetry$rel * scale
  source_structure <- if (defect == 0) "hermitian" else "general"
  source_identity <- psd_builtin_identity(
    list(kind = "dense", source = A, source_class = class(A)),
    c(n, n), source_structure
  )
  if (defect > tau_sym ||
      (defect > 0 && identical(policy$symmetry_repair, "reject"))) {
    psd_abort(
      "eigencore_psd_asymmetry_error", "asymmetric_input", "x",
      if (identical(policy$symmetry_repair, "reject")) {
        "exact symmetry"
      } else {
        paste0("symmetry defect <= ", format(tau_sym))
      },
      defect,
      source_identity = source_identity,
      representation = "dense_spectral",
      scale = scale,
      threshold = tau_sym,
      defect = defect,
      message = "Dense PSD source is not admitted by the declared symmetry policy."
    )
  }
  eig <- if (n == 0L) {
    list(values = numeric(), vectors = matrix(numeric(), 0L, 0L))
  } else {
    native_dense_symmetric_eigen(S)
  }
  ord <- order(eig$values, decreasing = TRUE)
  values <- eig$values[ord]
  vectors <- eig$vectors[, ord, drop = FALSE]
  dn <- dimnames(A)
  axis_names <- dn[[1L]] %||% dn[[2L]]
  rownames(vectors) <- axis_names
  psd_construct_complete_factor(
    original_values = values,
    vectors = vectors,
    source_index = seq_len(n),
    representation = "dense_spectral",
    method = "complete dense symmetric eigen PSD factor",
    source_identity = source_identity,
    policy = policy,
    state_kind = "dense",
    source_dimnames = list(axis_names, axis_names),
    symmetry_defect = defect,
    original_source = A,
    symmetric_source = S,
    algebraic_nullity = NA_integer_,
    materialization = list(
      source = "dense_matrix", factor = "dense_eigen",
      dense_n_by_n = TRUE, notes = character()
    ),
    started = started,
    scale_override = scale
  )
}

psd_construct_dense_gram <- function(x, orientation, policy, started,
                                     source_class) {
  if (!is.double(x) || is.complex(x) || length(dim(x)) != 2L) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dtype", "x",
      "a finite real-double matrix", typeof(x),
      representation = "gram_dense",
      capability = "construction"
    )
  }
  if (any(!is.finite(x))) psd_nonfinite_input(x, "x")
  defining <- if (identical(orientation, "columns")) x else t(x)
  n <- nrow(defining)
  q <- min(dim(defining))
  decomp <- if (q == 0L) {
    list(
      d = numeric(),
      u = matrix(numeric(), n, 0L),
      v = matrix(numeric(), ncol(defining), 0L)
    )
  } else {
    native_dense_svd(defining)
  }
  values <- c(as.numeric(decomp$d)^2, rep(0, n - length(decomp$d)))
  vectors <- decomp$u
  axis_names <- rownames(defining)
  rownames(vectors) <- axis_names
  gram <- tcrossprod(defining)
  source_identity <- psd_builtin_identity(
    list(
      kind = "dense_gram_factor", factor = x, orientation = orientation,
      source_class = source_class
    ),
    c(n, n), "hermitian"
  )
  psd_construct_complete_factor(
    original_values = values,
    vectors = vectors,
    source_index = seq_len(n),
    representation = "gram_dense",
    method = "complete dense Gram compact-SVD PSD factor",
    source_identity = source_identity,
    policy = policy,
    state_kind = "dense",
    source_dimnames = list(axis_names, axis_names),
    symmetry_defect = 0,
    original_source = gram,
    symmetric_source = gram,
    algebraic_nullity = NA_integer_,
    materialization = list(
      source = "dense_matrix", factor = "compact_svd",
      dense_n_by_n = TRUE,
      notes = "supplied Gram array validated by complete compact SVD"
    ),
    started = started,
    scale_override = psd_frobenius_norm(values)
  )
}

psd_construct_sparse_gram <- function(x, orientation, policy, started) {
  if (any(!is.finite(methods::slot(x, "x")))) psd_nonfinite_input(methods::slot(x, "x"), "x")
  defining_dim <- dim(x)
  n <- if (identical(orientation, "columns")) defining_dim[[1L]] else defining_dim[[2L]]
  axis_names <- if (identical(orientation, "columns")) rownames(x) else colnames(x)
  scale <- psd_sparse_gram_frobenius(x, orientation)
  if (!is.finite(scale)) {
    psd_abort(
      "eigencore_psd_invalid_input", "nonfinite_input", "scale",
      "a finite Gram Frobenius scale", scale,
      representation = "gram_sparse",
      capability = "construction"
    )
  }
  snapshot <- deep_copy_record(x)
  state <- list(
    schema_version = 1L,
    kind = "gram_sparse",
    rank = NA_integer_,
    dim = as.integer(c(n, n)),
    dimnames = list(axis_names, axis_names),
    factor = snapshot,
    orientation = orientation,
    action_frobenius_norm = scale
  )
  source_identity <- psd_builtin_identity(
    list(kind = "sparse_gram_factor", factor = snapshot, orientation = orientation),
    c(n, n), "hermitian"
  )
  evidence <- psd_structural_evidence(
    repaired = FALSE,
    theorem = if (identical(orientation, "columns")) {
      "K = L L^T is positive semidefinite"
    } else {
      "K = L^T L is positive semidefinite"
    },
    details = list(
      orientation = orientation,
      factor_dim = as.integer(defining_dim),
      algebraic_rank = "not inferred from sparse storage"
    )
  )
  psd_construct_structural_factor(
    dim = c(n, n),
    representation = "gram_sparse",
    method = "structural sparse Gram PSD factor",
    source_identity = source_identity,
    policy = policy,
    state = state,
    source_dimnames = list(axis_names, axis_names),
    scale = scale,
    symmetry_defect = 0,
    repair_defect = 0,
    source_action_defect = 0,
    algebraic_nullity = NA_integer_,
    rank_bounds = list(lower = 0L, upper = as.integer(min(defining_dim))),
    evidence = evidence,
    materialization = list(
      source = "sparse_csc", factor = "sparse_gram",
      dense_n_by_n = FALSE,
      notes = "retains the supplied CSC factor and no square dense action"
    ),
    started = started,
    spectrum_upper_bound = scale
  )
}

psd_construct_sparse_laplacian <- function(x, policy, started) {
  n <- nrow(x)
  if (n != ncol(x)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", "x",
      "a square sparse Laplacian", dim(x),
      representation = "laplacian_sparse",
      capability = "construction"
    )
  }
  if (any(!is.finite(methods::slot(x, "x")))) psd_nonfinite_input(methods::slot(x, "x"), "x")
  psd_validate_square_dimnames(x)
  snapshot <- deep_copy_record(x)
  transposed <- Matrix::t(snapshot)
  symmetry_defect <- psd_sparse_frobenius(snapshot - transposed)
  symmetric <- Matrix::drop0(0.5 * snapshot + 0.5 * transposed)
  scale <- psd_sparse_frobenius(symmetric)
  if (!is.finite(scale) || !is.finite(symmetry_defect)) {
    psd_abort(
      "eigencore_psd_invalid_input", "nonfinite_input", "scale",
      "finite sparse Frobenius diagnostics",
      list(scale = scale, symmetry_defect = symmetry_defect),
      representation = "laplacian_sparse",
      capability = "construction"
    )
  }
  tau_sym <- policy$symmetry$abs + policy$symmetry$rel * scale
  tau_psd <- policy$positivity$abs + policy$positivity$rel * scale
  source_identity <- psd_builtin_identity(
    list(kind = "sparse_laplacian", source = snapshot, source_class = class(x)),
    c(n, n), if (symmetry_defect == 0) "hermitian" else "general"
  )
  if (symmetry_defect > tau_sym ||
      (symmetry_defect > 0 && identical(policy$symmetry_repair, "reject"))) {
    psd_abort(
      "eigencore_psd_asymmetry_error", "asymmetric_input", "x",
      if (identical(policy$symmetry_repair, "reject")) {
        "exact symmetry"
      } else {
        paste0("symmetry defect <= ", format(tau_sym))
      },
      symmetry_defect,
      source_identity = source_identity,
      representation = "laplacian_sparse",
      capability = "construction",
      scale = scale,
      threshold = tau_sym,
      defect = symmetry_defect
    )
  }
  entries <- Matrix::summary(symmetric)
  positive_off_diagonal <- which(entries$i != entries$j & entries$x > 0)
  if (length(positive_off_diagonal)) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_structure", "off_diagonal",
      "exactly non-positive off-diagonal entries",
      entries$x[positive_off_diagonal],
      source_identity = source_identity,
      representation = "laplacian_sparse",
      capability = "construction",
      scale = scale,
      threshold = 0,
      defect = max(entries$x[positive_off_diagonal]),
      indices = positive_off_diagonal,
      message = "A structural graph Laplacian cannot have a positive off-diagonal weight."
    )
  }
  row_sums <- as.numeric(Matrix::rowSums(symmetric))
  row_sum_defect <- if (length(row_sums)) max(abs(row_sums)) else 0
  if (row_sum_defect > tau_psd ||
      (row_sum_defect > 0 && identical(policy$structure_repair, "reject"))) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_structure", "row_sums",
      if (identical(policy$structure_repair, "reject")) {
        "exact zero row sums"
      } else {
        paste0("maximum absolute row sum <= ", format(tau_psd))
      },
      row_sums,
      source_identity = source_identity,
      representation = "laplacian_sparse",
      capability = "construction",
      scale = scale,
      threshold = tau_psd,
      defect = row_sum_defect,
      indices = which(abs(row_sums) == row_sum_defect),
      message = "Sparse source does not satisfy the admitted graph-Laplacian row-sum structure."
    )
  }
  repaired <- symmetric
  if (row_sum_defect > 0) {
    off_diagonal <- symmetric
    diag(off_diagonal) <- 0
    repaired <- off_diagonal
    diag(repaired) <- -as.numeric(Matrix::rowSums(off_diagonal))
    repaired <- Matrix::drop0(repaired)
  }
  repaired <- Matrix::forceSymmetric(repaired, uplo = "U")
  repair_defect <- psd_sparse_frobenius(symmetric - repaired)
  source_action_defect <- psd_sparse_frobenius(snapshot - repaired)
  components <- psd_sparse_components(repaired)
  axis_names <- rownames(repaired) %||% colnames(repaired)
  action_norm <- psd_sparse_frobenius(repaired)
  state <- list(
    schema_version = 1L,
    kind = "laplacian_sparse",
    rank = NA_integer_,
    dim = as.integer(c(n, n)),
    dimnames = list(axis_names, axis_names),
    matrix = repaired,
    components = components$membership,
    action_frobenius_norm = action_norm
  )
  evidence <- psd_structural_evidence(
    repaired = symmetry_defect > 0 || row_sum_defect > 0,
    theorem = "symmetric non-positive-edge graph Laplacian component theorem",
    details = list(
      component_count = as.integer(components$count),
      edge_count = as.integer(components$edge_count),
      admitted_row_sum_defect = row_sum_defect,
      repaired_row_sum_defect = if (n) {
        max(abs(as.numeric(Matrix::rowSums(repaired))))
      } else {
        0
      }
    )
  )
  psd_construct_structural_factor(
    dim = c(n, n),
    representation = "laplacian_sparse",
    method = "structural sparse graph-Laplacian PSD factor",
    source_identity = source_identity,
    policy = policy,
    state = state,
    source_dimnames = list(axis_names, axis_names),
    scale = scale,
    symmetry_defect = symmetry_defect,
    repair_defect = repair_defect,
    source_action_defect = source_action_defect,
    algebraic_nullity = as.integer(components$count),
    rank_bounds = list(
      lower = 0L,
      upper = as.integer(n - components$count)
    ),
    evidence = evidence,
    materialization = list(
      source = "sparse_csc", factor = "sparse_laplacian",
      dense_n_by_n = FALSE,
      notes = "retains a canonical sparse graph Laplacian"
    ),
    started = started,
    spectrum_upper_bound = action_norm
  )
}

psd_construct_complete_factor <- function(
    original_values, vectors, source_index, representation, method,
    source_identity, policy, state_kind, source_dimnames, symmetry_defect,
    original_source, algebraic_nullity, materialization, started,
    symmetric_source = NULL, scale_override = NULL) {
  n <- length(original_values)
  scale <- scale_override %||% psd_frobenius_norm(original_values)
  thresholds <- list(
    scale_definition = "frobenius",
    symmetry = policy$symmetry$abs + policy$symmetry$rel * scale,
    positivity = policy$positivity$abs + policy$positivity$rel * scale,
    rank = policy$rank$abs + policy$rank$rel * scale
  )
  category <- rep("retained_positive", n)
  category[original_values < -thresholds$positivity] <- "materially_negative"
  category[original_values >= -thresholds$positivity & original_values < 0] <-
    "accepted_negative"
  category[original_values == 0] <- "exact_zero"
  category[original_values > 0 & original_values <= thresholds$rank] <-
    "numerical_null"

  materially_negative <- which(category == "materially_negative")
  accepted_negative <- which(category == "accepted_negative")
  if (length(materially_negative) ||
      (length(accepted_negative) && identical(policy$negative_repair, "reject"))) {
    offending <- if (length(materially_negative)) materially_negative else accepted_negative
    failure_code <- "indefinite_input"
    evidence <- psd_complete_evidence(
      repaired = FALSE,
      source_semantics = "immutable_snapshot"
    )
    failed_cert <- psd_new_certificate(
      passed = FALSE,
      certificate_type = "psd_factor_validation",
      scope = "source_validation",
      thresholds = thresholds,
      source_identity = source_identity,
      factor_identity = NULL,
      representation = representation,
      evidence = evidence,
      repair_applied = FALSE,
      symmetry_defect = symmetry_defect,
      repair_defect = NA_real_,
      source_action_defect = NA_real_,
      original_spectrum = original_values,
      repaired_spectrum = NULL,
      classification = psd_classification_counts(category),
      rank = NA_integer_,
      nullity = NA_integer_,
      algebraic_nullity = algebraic_nullity,
      rank_bounds = list(lower = 0L, upper = as.integer(n)),
      capabilities = NULL,
      action_bounds = list(),
      residuals = list(
        symmetry = symmetry_defect,
        most_negative = if (length(original_values)) min(original_values) else 0
      ),
      orthogonality = numeric(),
      notes = if (length(materially_negative)) {
        "materially negative spectrum"
      } else {
        "negative repair rejected by policy"
      }
    )
    psd_abort(
      "eigencore_psd_indefinite_error", failure_code, "spectrum",
      if (length(materially_negative)) {
        paste0("lambda >= -", format(thresholds$positivity))
      } else {
        "no negative eigenvalue when negative_repair = reject"
      },
      original_values[offending],
      source_identity = source_identity,
      representation = representation,
      evidence = evidence,
      scale = scale,
      threshold = thresholds$positivity,
      defect = min(original_values[offending]),
      indices = offending,
      details = list(certificate = failed_cert, category = category),
      message = if (length(materially_negative)) {
        "PSD source has a materially negative eigenvalue."
      } else {
        "PSD source has a negative eigenvalue and the declared policy rejects clipping."
      }
    )
  }

  repaired <- original_values
  repaired[category %in% c("accepted_negative", "exact_zero", "numerical_null")] <- 0
  retained <- which(category == "retained_positive")
  rank <- as.integer(length(retained))
  nullity <- as.integer(n - rank)
  classification <- psd_classification_counts(category)
  repair_applied <- symmetry_defect > 0 ||
    length(accepted_negative) > 0L || any(category == "numerical_null")
  fidelity <- if (repair_applied) {
    "repaired_with_defect"
  } else {
    "exact_for_certified_factor"
  }
  evidence <- psd_complete_evidence(
    repaired = repair_applied,
    source_semantics = "immutable_snapshot"
  )

  if (identical(state_kind, "identity")) {
    state <- list(
      schema_version = 1L,
      kind = "identity",
      rank = rank,
      dim = c(n, n),
      dimnames = source_dimnames,
      active = identical(rank, as.integer(n))
    )
    repair_defect <- psd_frobenius_norm(original_values - repaired)
    source_action_defect <- repair_defect
    orthogonality <- 0
  } else if (identical(state_kind, "diagonal")) {
    repaired_source <- numeric(n)
    repaired_source[source_index] <- repaired
    state <- list(
      schema_version = 1L,
      kind = "diagonal",
      rank = rank,
      dim = c(n, n),
      dimnames = source_dimnames,
      repaired = repaired_source,
      retained_source_index = source_index[retained],
      retained_values = repaired[retained]
    )
    repair_defect <- psd_frobenius_norm(original_values - repaired)
    source_action_defect <- repair_defect
    orthogonality <- 0
  } else {
    Qr <- vectors[, retained, drop = FALSE]
    retained_values <- repaired[retained]
    state <- list(
      schema_version = 1L,
      kind = "dense",
      rank = rank,
      dim = c(n, n),
      dimnames = source_dimnames,
      basis = Qr,
      retained_values = retained_values
    )
    Kr <- if (rank == 0L) {
      matrix(0, n, n)
    } else {
      tcrossprod(sweep(Qr, 2L, retained_values, `*`), Qr)
    }
    repair_defect <- psd_frobenius_norm(symmetric_source - Kr)
    source_action_defect <- psd_frobenius_norm(original_source - Kr)
    orthogonality <- if (rank == 0L) 0 else {
      max(abs(crossprod(Qr) - diag(rank)))
    }
  }

  state_token <- stable_raw_hash(state)
  factor_identity <- psd_factor_identity(
    source_identity, policy, representation, state_token, c(n, n)
  )
  capabilities <- psd_complete_capabilities(
    representation = representation,
    evidence = evidence,
    materialization = materialization$factor,
    algebraic_available = !is.na(algebraic_nullity)
  )
  spectrum <- list(
    coverage = "complete",
    original = as.numeric(original_values),
    repaired = as.numeric(repaired),
    category = category,
    source_index = as.integer(source_index),
    retained_indices = as.integer(retained),
    exact_zero_indices = as.integer(which(category == "exact_zero")),
    accepted_negative_indices = as.integer(accepted_negative),
    numerical_null_indices = as.integer(which(category == "numerical_null")),
    signed_zero_count = as.integer(psd_signed_zero_count(original_values)),
    lower_bound = if (n) min(original_values) else 0,
    upper_bound = if (n) max(original_values) else 0
  )
  rank_bounds <- list(lower = rank, upper = rank)
  certificate <- psd_new_certificate(
    passed = TRUE,
    certificate_type = "psd_factor_validation",
    scope = "source_validation_and_factor_actions",
    thresholds = thresholds,
    source_identity = source_identity,
    factor_identity = factor_identity,
    representation = representation,
    evidence = evidence,
    repair_applied = repair_applied,
    symmetry_defect = symmetry_defect,
    repair_defect = repair_defect,
    source_action_defect = source_action_defect,
    original_spectrum = original_values,
    repaired_spectrum = repaired,
    classification = classification,
    rank = rank,
    nullity = nullity,
    algebraic_nullity = algebraic_nullity,
    rank_bounds = rank_bounds,
    capabilities = capabilities,
    action_bounds = list(
      source_action_frobenius = source_action_defect,
      projector_orthogonality = orthogonality
    ),
    residuals = list(
      symmetry = symmetry_defect,
      repair = repair_defect,
      source_action = source_action_defect
    ),
    orthogonality = c(factor_basis = orthogonality),
    notes = if (repair_applied) "factor action includes recorded tolerance repair" else character()
  )
  if (!isTRUE(certificate$passed)) {
    psd_abort(
      "eigencore_psd_incomplete_evidence", "incomplete_evidence", "certificate",
      "a passing factor-action certificate", certificate,
      source_identity = source_identity,
      factor_identity = factor_identity,
      representation = representation,
      capability = "factor_actions",
      evidence = evidence,
      scale = scale,
      threshold = certificate$orthogonality_tolerance,
      defect = certificate$max_orthogonality_loss,
      indices = certificate$failed_indices,
      details = list(certificate = certificate),
      message = "The candidate PSD factor failed its construction certificate."
    )
  }
  elapsed <- proc.time()[["elapsed"]] - started
  work <- psd_work_record(setup_seconds = elapsed, total_seconds = elapsed)
  memory <- new_memory_record(
    list(
      source_snapshot = list(),
      factor_state = state,
      cached_actions = list(),
      metadata = list(
        spectrum = spectrum,
        classification = classification,
        policy = policy,
        evidence = evidence,
        certificate = certificate
      )
    ),
    native_bytes = 0
  )
  serialization <- structure(
    list(
      schema_version = 1L,
      portable = isTRUE(factor_identity$portable),
      originating_session = eigencore_session_id(),
      incompatibility_reason = NULL,
      integrity_token = NULL,
      source_identity_token = stable_raw_hash(source_identity),
      operator_identity_token = stable_raw_hash(factor_identity),
      policy_token = stable_raw_hash(policy),
      state_token = state_token,
      reconstruction = representation
    ),
    class = "eigencore_psd_serialization"
  )
  factor <- structure(
    list(
      schema_version = 1L,
      dim = as.integer(c(n, n)),
      dtype = "double",
      representation = representation,
      method = method,
      policy = policy,
      scale = scale,
      thresholds = thresholds,
      spectrum = spectrum,
      classification = classification,
      rank = rank,
      nullity = nullity,
      algebraic_nullity = algebraic_nullity,
      rank_bounds = rank_bounds,
      evidence = evidence,
      capabilities = capabilities,
      operator_identity = factor_identity,
      source_identity = source_identity,
      source_semantics = "immutable_snapshot",
      materialization = materialization,
      certificate = certificate,
      work = work,
      memory = memory,
      serialization = serialization,
      warnings = character()
    ),
    class = "eigencore_psd_factor"
  )
  attr(factor, "eigencore_psd_state") <- state
  factor$serialization$integrity_token <- stable_raw_hash(psd_integrity_payload(factor))
  factor
}

psd_construct_structural_factor <- function(
    dim, representation, method, source_identity, policy, state,
    source_dimnames, scale, symmetry_defect, repair_defect,
    source_action_defect, algebraic_nullity, rank_bounds, evidence,
    materialization, started, spectrum_upper_bound) {
  dim <- as.integer(dim)
  state$dimnames <- source_dimnames
  thresholds <- list(
    scale_definition = "frobenius",
    symmetry = policy$symmetry$abs + policy$symmetry$rel * scale,
    positivity = policy$positivity$abs + policy$positivity$rel * scale,
    rank = policy$rank$abs + policy$rank$rel * scale
  )
  classification <- psd_classification_counts(character())
  spectrum <- list(
    coverage = "structural",
    original = NULL,
    repaired = NULL,
    category = NULL,
    source_index = NULL,
    retained_indices = NULL,
    exact_zero_indices = NULL,
    accepted_negative_indices = NULL,
    numerical_null_indices = NULL,
    signed_zero_count = 0L,
    lower_bound = 0,
    upper_bound = as.numeric(spectrum_upper_bound)
  )
  state_token <- stable_raw_hash(state)
  factor_identity <- psd_factor_identity(
    source_identity, policy, representation, state_token, dim
  )
  capabilities <- psd_structural_capabilities(
    representation = representation,
    evidence = evidence,
    materialization = materialization$factor,
    algebraic_available = !is.na(algebraic_nullity)
  )
  repair_applied <- symmetry_defect > 0 || repair_defect > 0 ||
    source_action_defect > 0
  certificate <- psd_new_certificate(
    passed = TRUE,
    certificate_type = "psd_factor_validation",
    scope = "structural_source_validation_and_form_action",
    thresholds = thresholds,
    source_identity = source_identity,
    factor_identity = factor_identity,
    representation = representation,
    evidence = evidence,
    repair_applied = repair_applied,
    symmetry_defect = symmetry_defect,
    repair_defect = repair_defect,
    source_action_defect = source_action_defect,
    original_spectrum = NULL,
    repaired_spectrum = NULL,
    classification = classification,
    rank = NA_integer_,
    nullity = NA_integer_,
    algebraic_nullity = algebraic_nullity,
    rank_bounds = rank_bounds,
    capabilities = capabilities,
    action_bounds = list(form_frobenius = spectrum_upper_bound),
    residuals = list(
      symmetry = symmetry_defect,
      repair = repair_defect,
      source_action = source_action_defect
    ),
    orthogonality = numeric(),
    notes = evidence$theorem,
    certificate_scale = scale
  )
  elapsed <- proc.time()[["elapsed"]] - started
  work <- psd_work_record(setup_seconds = elapsed, total_seconds = elapsed)
  memory <- new_memory_record(
    list(
      source_snapshot = list(),
      factor_state = state,
      cached_actions = list(),
      metadata = list(
        spectrum = spectrum,
        classification = classification,
        policy = policy,
        evidence = evidence,
        certificate = certificate
      )
    ),
    native_bytes = 0
  )
  serialization <- structure(
    list(
      schema_version = 1L,
      portable = TRUE,
      originating_session = eigencore_session_id(),
      incompatibility_reason = NULL,
      integrity_token = NULL,
      source_identity_token = stable_raw_hash(source_identity),
      operator_identity_token = stable_raw_hash(factor_identity),
      policy_token = stable_raw_hash(policy),
      state_token = state_token,
      reconstruction = representation
    ),
    class = "eigencore_psd_serialization"
  )
  factor <- structure(
    list(
      schema_version = 1L,
      dim = dim,
      dtype = "double",
      representation = representation,
      method = method,
      policy = policy,
      scale = scale,
      thresholds = thresholds,
      spectrum = spectrum,
      classification = classification,
      rank = NA_integer_,
      nullity = NA_integer_,
      algebraic_nullity = algebraic_nullity,
      rank_bounds = rank_bounds,
      evidence = evidence,
      capabilities = capabilities,
      operator_identity = factor_identity,
      source_identity = source_identity,
      source_semantics = "immutable_snapshot",
      materialization = materialization,
      certificate = certificate,
      work = work,
      memory = memory,
      serialization = serialization,
      warnings = character()
    ),
    class = "eigencore_psd_factor"
  )
  attr(factor, "eigencore_psd_state") <- state
  factor$serialization$integrity_token <- stable_raw_hash(psd_integrity_payload(factor))
  factor
}

psd_classification_counts <- function(category) {
  keys <- c(
    "retained_positive", "exact_zero", "accepted_negative", "numerical_null",
    "materially_negative", "user_truncated"
  )
  out <- stats::setNames(as.list(integer(length(keys))), keys)
  counts <- table(factor(category, levels = keys))
  for (key in keys) out[[key]] <- as.integer(counts[[key]])
  out
}

psd_complete_evidence <- function(repaired, source_semantics) {
  structure(list(
    schema_version = 1L,
    spectrum_coverage = "complete",
    validation = "computed",
    action_fidelity = if (repaired) "repaired_with_defect" else "exact_for_certified_factor",
    source_semantics = source_semantics,
    theorem = "complete real symmetric spectrum",
    bound_type = "exact_frobenius_diagnostics",
    details = list()
  ), class = "eigencore_psd_evidence")
}

psd_structural_evidence <- function(repaired, theorem, details) {
  structure(list(
    schema_version = 1L,
    spectrum_coverage = "structural",
    validation = "computed",
    action_fidelity = if (repaired) "repaired_with_defect" else "exact_for_certified_factor",
    source_semantics = "immutable_snapshot",
    theorem = theorem,
    bound_type = "exact_structural_identity",
    details = details
  ), class = "eigencore_psd_evidence")
}

psd_unavailable_evidence <- function(source_semantics) {
  structure(list(
    schema_version = 1L,
    spectrum_coverage = "partial",
    validation = "unavailable",
    action_fidelity = "approximate_with_bound",
    source_semantics = source_semantics,
    theorem = NULL,
    bound_type = "none",
    details = list()
  ), class = "eigencore_psd_evidence")
}

psd_capability_entry <- function(available, fidelity, evidence_required,
                                 materialization, reason = NULL) {
  list(
    available = isTRUE(available),
    fidelity = fidelity,
    evidence_required = evidence_required,
    materialization = materialization,
    reason = if (isTRUE(available)) NULL else reason
  )
}

psd_complete_capabilities <- function(representation, evidence,
                                      materialization, algebraic_available) {
  out <- list(
    schema_version = 1L,
    representation = representation,
    evidence = evidence
  )
  for (name in .psd_capability_names) {
    available <- !identical(name, "algebraic_nullity") || algebraic_available
    out[[name]] <- psd_capability_entry(
      available = available,
      fidelity = evidence$action_fidelity,
      evidence_required = if (identical(name, "algebraic_nullity")) {
        "structural exact-zero evidence"
      } else {
        "complete classified spectrum"
      },
      materialization = materialization,
      reason = if (available) NULL else "incomplete_evidence"
    )
  }
  structure(out, class = "eigencore_psd_capabilities")
}

psd_structural_capabilities <- function(representation, evidence,
                                        materialization, algebraic_available) {
  out <- list(
    schema_version = 1L,
    representation = representation,
    evidence = evidence
  )
  available_names <- c("form", "gram", "serialization", "cache_reuse")
  if (algebraic_available) {
    available_names <- c(available_names, "algebraic_nullity")
  }
  for (name in .psd_capability_names) {
    available <- name %in% available_names
    evidence_required <- if (identical(name, "algebraic_nullity")) {
      "a structural null-space theorem"
    } else if (name %in% c("form", "gram")) {
      "a supplied Gram identity or admitted graph-Laplacian theorem"
    } else if (name %in% c("serialization", "cache_reuse")) {
      "an immutable portable sparse snapshot"
    } else {
      "complete or threshold-separating spectral-factor evidence"
    }
    out[[name]] <- psd_capability_entry(
      available = available,
      fidelity = evidence$action_fidelity,
      evidence_required = evidence_required,
      materialization = materialization,
      reason = if (available) NULL else "incomplete_evidence"
    )
  }
  structure(out, class = "eigencore_psd_capabilities")
}

psd_new_certificate <- function(
    passed, certificate_type, scope, thresholds, source_identity,
    factor_identity, representation, evidence, repair_applied,
    symmetry_defect, repair_defect, source_action_defect, original_spectrum,
    repaired_spectrum, classification, rank, nullity, algebraic_nullity,
    rank_bounds, capabilities, action_bounds, residuals, orthogonality,
    notes = character(), certificate_scale = NULL) {
  residual_values <- unlist(residuals, recursive = TRUE, use.names = FALSE)
  residual_values <- residual_values[is.finite(residual_values)]
  max_residual <- if (length(residual_values)) max(abs(residual_values)) else NA_real_
  orthogonality <- as.numeric(orthogonality)
  max_orthogonality <- if (length(orthogonality)) max(abs(orthogonality)) else NA_real_
  orth_tol <- sqrt(.Machine$double.eps)
  orth_passed <- is.na(max_orthogonality) || max_orthogonality <= orth_tol
  scale <- if (!is.null(certificate_scale)) {
    as.numeric(certificate_scale)
  } else if (!is.null(original_spectrum)) {
    psd_frobenius_norm(original_spectrum)
  } else {
    NA_real_
  }
  denom <- if (is.finite(scale) && scale > 0) scale else 1
  backward_error <- lapply(residuals, function(value) as.numeric(value) / denom)
  backward_values <- unlist(backward_error, recursive = TRUE, use.names = FALSE)
  max_backward <- if (length(backward_values)) max(abs(backward_values)) else NA_real_
  threshold_values <- psd_numeric_leaves(thresholds)
  threshold_values <- threshold_values[is.finite(threshold_values)]
  tolerance <- if (length(threshold_values)) max(threshold_values) else NA_real_
  pass_flags <- as.logical(passed)
  scope_passed <- all(!is.na(pass_flags) & pass_flags)
  failed_indices <- which(is.na(pass_flags) | !pass_flags)
  if (!orth_passed && !length(failed_indices)) failed_indices <- 1L
  cert <- list(
    passed = scope_passed && orth_passed,
    tolerance = as.numeric(tolerance),
    orthogonality_tolerance = orth_tol,
    orthogonality_required = TRUE,
    certificate_type = certificate_type,
    norm_bound_type = "exact_frobenius",
    scale_is_estimate = FALSE,
    max_backward_error = max_backward,
    max_residual = max_residual,
    max_orthogonality_loss = max_orthogonality,
    orthogonality_passed = orth_passed,
    failed_indices = as.integer(failed_indices),
    scale = scale,
    notes = as.character(notes),
    residuals = residuals,
    backward_error = backward_error,
    orthogonality = orthogonality,
    converged = scope_passed && orth_passed,
    schema_version = 1L,
    scope = scope,
    thresholds = thresholds,
    source_identity = source_identity,
    factor_identity = factor_identity,
    representation = representation,
    evidence = evidence,
    repair_applied = isTRUE(repair_applied),
    symmetry_defect = symmetry_defect,
    repair_defect = repair_defect,
    source_action_defect = source_action_defect,
    original_spectrum = original_spectrum,
    repaired_spectrum = repaired_spectrum,
    classification = classification,
    rank = rank,
    nullity = nullity,
    algebraic_nullity = algebraic_nullity,
    rank_bounds = rank_bounds,
    capabilities = capabilities,
    action_bounds = action_bounds
  )
  class(cert) <- c("eigencore_psd_certificate", "eigencore_certificate")
  cert
}

psd_operation_certificate <- function(type, x, thresholds, residuals, passed,
                                      notes, action_bounds = list()) {
  psd_new_certificate(
    passed = passed,
    certificate_type = type,
    scope = if (identical(type, "psd_strict_solve")) "right_hand_side" else "block_action",
    thresholds = thresholds,
    source_identity = x$source_identity,
    factor_identity = x$operator_identity,
    representation = x$representation,
    evidence = x$evidence,
    repair_applied = x$certificate$repair_applied,
    symmetry_defect = x$certificate$symmetry_defect,
    repair_defect = x$certificate$repair_defect,
    source_action_defect = x$certificate$source_action_defect,
    original_spectrum = x$spectrum$original,
    repaired_spectrum = x$spectrum$repaired,
    classification = x$classification,
    rank = x$rank,
    nullity = x$nullity,
    algebraic_nullity = x$algebraic_nullity,
    rank_bounds = x$rank_bounds,
    capabilities = x$capabilities,
    action_bounds = action_bounds,
    residuals = residuals,
    orthogonality = numeric(),
    notes = notes,
    certificate_scale = x$scale
  )
}

psd_numeric_leaves <- function(x) {
  if (is.list(x)) {
    return(unlist(lapply(x, psd_numeric_leaves), recursive = TRUE, use.names = FALSE))
  }
  if (is.numeric(x)) as.numeric(x) else numeric()
}

psd_work_record <- function(setup_seconds = 0, solve_seconds = 0,
                            certification_seconds = 0, total_seconds = NULL) {
  counters <- stats::setNames(
    as.list(integer(length(work_counter_fields()))),
    work_counter_fields()
  )
  if (is.null(total_seconds)) {
    total_seconds <- setup_seconds + solve_seconds + certification_seconds
  }
  new_typed_work_record(c(
    counters,
    list(
      iterations = 0L,
      restarts = 0L,
      setup_seconds = setup_seconds,
      solve_seconds = solve_seconds,
      certification_seconds = certification_seconds,
      total_seconds = total_seconds,
      legacy_matvecs = 0L
    )
  ))
}

psd_signed_zero_count <- function(x) {
  zeros <- x == 0
  if (!any(zeros)) return(0L)
  reciprocals <- 1 / x[zeros]
  as.integer(sum(is.infinite(reciprocals) & reciprocals < 0))
}

psd_integrity_payload <- function(x) {
  public <- unclass(x)
  public$serialization$integrity_token <- NULL
  if (inherits(public$work, "eigencore_work")) {
    for (field in c(
      "setup_seconds", "solve_seconds", "certification_seconds", "total_seconds"
    )) {
      public$work[[field]] <- 0
    }
  }
  list(public = public, state = attr(x, "eigencore_psd_state", exact = TRUE))
}

validate_psd_factor <- function(x) {
  if (!inherits(x, "eigencore_psd_factor") ||
      !identical(names(x), .psd_factor_fields)) {
    psd_abort(
      "eigencore_psd_corrupt_state", "corrupt_factor", "class/fields",
      list(class = "eigencore_psd_factor", fields = .psd_factor_fields),
      list(class = class(x), fields = names(x))
    )
  }
  if (!identical(x$schema_version, 1L) ||
      !inherits(x$serialization, "eigencore_psd_serialization") ||
      !identical(x$serialization$schema_version, 1L)) {
    psd_abort(
      "eigencore_psd_corrupt_state", "unsupported_schema", "schema_version",
      1L, x$schema_version %||% NULL,
      source_identity = x$source_identity %||% NULL,
      factor_identity = x$operator_identity %||% NULL,
      representation = x$representation %||% NULL
    )
  }
  state <- attr(x, "eigencore_psd_state", exact = TRUE)
  if (!is.list(state) ||
      !identical(stable_raw_hash(state), x$serialization$state_token)) {
    psd_abort(
      "eigencore_psd_corrupt_state", "corrupt_factor", "factor_state",
      x$serialization$state_token, if (is.list(state)) stable_raw_hash(state) else class(state),
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation
    )
  }
  expected <- stable_raw_hash(psd_integrity_payload(x))
  if (!identical(expected, x$serialization$integrity_token)) {
    psd_abort(
      "eigencore_psd_corrupt_state", "corrupt_factor", "integrity_token",
      expected, x$serialization$integrity_token,
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation
    )
  }
  if (!isTRUE(x$serialization$portable) &&
      !identical(x$serialization$originating_session, eigencore_session_id())) {
    psd_abort(
      "eigencore_psd_nonportable", "session_incompatible", "session",
      x$serialization$originating_session, eigencore_session_id(),
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation
    )
  }
  x
}

psd_incomplete_evidence <- function(x, capability) {
  psd_abort(
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", capability,
    "complete or threshold-separating structural evidence", x$evidence,
    source_identity = x$source_identity,
    factor_identity = x$operator_identity,
    representation = x$representation,
    capability = capability,
    evidence = x$evidence,
    scale = x$scale,
    threshold = x$thresholds$rank
  )
}

psd_require_capability <- function(x, capability) {
  entry <- x$capabilities[[capability]] %||% NULL
  if (is.null(entry)) {
    psd_abort(
      "eigencore_psd_unsupported_action", "unsupported_action", capability,
      .psd_capability_names, capability,
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation,
      capability = capability,
      evidence = x$evidence
    )
  }
  if (!isTRUE(entry$available)) {
    if (identical(entry$reason, "incomplete_evidence")) {
      psd_incomplete_evidence(x, capability)
    }
    psd_abort(
      "eigencore_psd_unsupported_action", "unsupported_action", capability,
      "an available certified action", entry,
      source_identity = x$source_identity,
      factor_identity = x$operator_identity,
      representation = x$representation,
      capability = capability,
      evidence = x$evidence
    )
  }
  invisible(entry)
}

psd_prepare_block <- function(X, n, field, factor) {
  vector <- is.atomic(X) && is.null(dim(X))
  if (!is.double(X) || is.complex(X) ||
      (!vector && (!is.matrix(X) || length(dim(X)) != 2L))) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dtype", field,
      "a finite real-double vector or matrix", class(X),
      source_identity = factor$source_identity,
      factor_identity = factor$operator_identity,
      representation = factor$representation
    )
  }
  value <- if (vector) matrix(X, ncol = 1L) else X
  if (vector && !is.null(names(X))) rownames(value) <- names(X)
  if (nrow(value) != n) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", field,
      paste0(n, " rows"), dim(value),
      source_identity = factor$source_identity,
      factor_identity = factor$operator_identity,
      representation = factor$representation
    )
  }
  if (any(!is.finite(value))) psd_nonfinite_input(value, field, factor)
  expected_names <- attr(factor, "eigencore_psd_state", exact = TRUE)$dimnames[[1L]]
  if (!is.null(expected_names) && !is.null(rownames(value)) &&
      !identical(expected_names, rownames(value))) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", paste0(field, "$rownames"),
      expected_names, rownames(value),
      source_identity = factor$source_identity,
      factor_identity = factor$operator_identity,
      representation = factor$representation
    )
  }
  list(value = value, vector = vector, names = names(X), colnames = colnames(value))
}

psd_prepare_reduced_block <- function(Z, rank, field, factor) {
  vector <- is.atomic(Z) && is.null(dim(Z))
  if (!is.double(Z) || is.complex(Z) ||
      (!vector && (!is.matrix(Z) || length(dim(Z)) != 2L))) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dtype", field,
      "a finite real-double reduced vector or matrix", class(Z),
      source_identity = factor$source_identity,
      factor_identity = factor$operator_identity,
      representation = factor$representation
    )
  }
  value <- if (vector) matrix(Z, ncol = 1L) else Z
  if (nrow(value) != rank) {
    psd_abort(
      "eigencore_psd_invalid_input", "invalid_dimension", field,
      paste0(rank, " rows"), dim(value),
      source_identity = factor$source_identity,
      factor_identity = factor$operator_identity,
      representation = factor$representation
    )
  }
  if (any(!is.finite(value))) psd_nonfinite_input(value, field, factor)
  list(value = value, vector = vector, names = names(Z), colnames = colnames(value))
}

psd_restore_block_shape <- function(out, block, factor) {
  state <- attr(factor, "eigencore_psd_state", exact = TRUE)
  rownames(out) <- state$dimnames[[1L]]
  colnames(out) <- block$colnames
  if (isTRUE(block$vector)) {
    value <- as.numeric(out[, 1L])
    names(value) <- state$dimnames[[1L]] %||% block$names
    value
  } else {
    out
  }
}

psd_restore_reduced_shape <- function(out, block) {
  rownames(out) <- NULL
  colnames(out) <- block$colnames
  if (isTRUE(block$vector)) as.numeric(out[, 1L]) else out
}

psd_apply_matrix <- function(x, X, action) {
  state <- attr(x, "eigencore_psd_state", exact = TRUE)
  if (identical(state$kind, "identity")) {
    if (isTRUE(state$active)) {
      if (identical(action, "null_projector")) matrix(0, nrow(X), ncol(X)) else X
    } else {
      if (identical(action, "null_projector")) X else matrix(0, nrow(X), ncol(X))
    }
  } else if (identical(state$kind, "diagonal")) {
    values <- state$repaired
    weights <- switch(
      action,
      form = values,
      sqrt = sqrt(values),
      inverse_sqrt = ifelse(values > 0, 1 / sqrt(values), 0),
      pseudoinverse = ifelse(values > 0, 1 / values, 0),
      image_projector = as.numeric(values > 0),
      null_projector = as.numeric(values == 0)
    )
    weights * X
  } else if (identical(state$kind, "gram_sparse")) {
    if (!identical(action, "form")) {
      stop("Internal error: unavailable sparse Gram action reached execution.", call. = FALSE)
    }
    out <- if (identical(state$orientation, "columns")) {
      state$factor %*% Matrix::crossprod(state$factor, X)
    } else {
      Matrix::crossprod(state$factor, state$factor %*% X)
    }
    as.matrix(out)
  } else if (identical(state$kind, "laplacian_sparse")) {
    if (!identical(action, "form")) {
      stop("Internal error: unavailable sparse Laplacian action reached execution.", call. = FALSE)
    }
    as.matrix(state$matrix %*% X)
  } else {
    Q <- state$basis
    if (identical(action, "null_projector")) {
      if (state$rank == 0L) return(X)
      return(X - Q %*% crossprod(Q, X))
    }
    if (identical(action, "image_projector")) {
      if (state$rank == 0L) return(matrix(0, nrow(X), ncol(X)))
      return(Q %*% crossprod(Q, X))
    }
    if (state$rank == 0L) return(matrix(0, nrow(X), ncol(X)))
    weights <- switch(
      action,
      form = state$retained_values,
      sqrt = sqrt(state$retained_values),
      inverse_sqrt = 1 / sqrt(state$retained_values),
      pseudoinverse = 1 / state$retained_values
    )
    Q %*% (weights * crossprod(Q, X))
  }
}

psd_reduce_matrix <- function(x, X) {
  state <- attr(x, "eigencore_psd_state", exact = TRUE)
  if (identical(state$kind, "identity")) {
    if (isTRUE(state$active)) return(X)
    return(matrix(numeric(), 0L, ncol(X)))
  }
  if (identical(state$kind, "diagonal")) {
    if (state$rank == 0L) return(matrix(numeric(), 0L, ncol(X)))
    return(sqrt(state$retained_values) *
      X[state$retained_source_index, , drop = FALSE])
  }
  if (state$rank == 0L) return(matrix(numeric(), 0L, ncol(X)))
  sqrt(state$retained_values) * crossprod(state$basis, X)
}

psd_lift_matrix <- function(x, Z) {
  state <- attr(x, "eigencore_psd_state", exact = TRUE)
  if (identical(state$kind, "identity")) {
    if (isTRUE(state$active)) return(Z)
    return(matrix(0, state$dim[[1L]], ncol(Z)))
  }
  if (identical(state$kind, "diagonal")) {
    out <- matrix(0, state$dim[[1L]], ncol(Z))
    if (state$rank) {
      out[state$retained_source_index, ] <-
        (1 / sqrt(state$retained_values)) * Z
    }
    return(out)
  }
  if (state$rank == 0L) return(matrix(0, state$dim[[1L]], ncol(Z)))
  state$basis %*% ((1 / sqrt(state$retained_values)) * Z)
}

psd_inverse_reduce_matrix <- function(x, X) {
  state <- attr(x, "eigencore_psd_state", exact = TRUE)
  if (identical(state$kind, "identity")) {
    if (isTRUE(state$active)) return(X)
    return(matrix(numeric(), 0L, ncol(X)))
  }
  if (identical(state$kind, "diagonal")) {
    if (state$rank == 0L) return(matrix(numeric(), 0L, ncol(X)))
    return((1 / sqrt(state$retained_values)) *
      X[state$retained_source_index, , drop = FALSE])
  }
  if (state$rank == 0L) return(matrix(numeric(), 0L, ncol(X)))
  (1 / sqrt(state$retained_values)) * crossprod(state$basis, X)
}

psd_action_identity <- function(x, action) {
  if (identical(action, "form")) return(x$operator_identity)
  token <- stable_raw_hash(list(
    schema_version = 1L,
    parent = unclass(x$operator_identity),
    action = action
  ))
  new_operator_identity(
    operator_id = paste0("psd-action-", stable_raw_hash(list(
      x$operator_identity$operator_id, action
    ))),
    revision = token,
    origin = "composite",
    dim = x$dim,
    dtype = "double",
    structure = "hermitian",
    portable = x$operator_identity$portable
  )
}

psd_action_frobenius_norm <- function(x, action) {
  state <- attr(x, "eigencore_psd_state", exact = TRUE)
  if (state$kind %in% c("gram_sparse", "laplacian_sparse")) {
    return(state$action_frobenius_norm)
  }
  values <- x$spectrum$repaired
  transformed <- switch(
    action,
    form = values,
    sqrt = sqrt(values),
    inverse_sqrt = ifelse(values > 0, 1 / sqrt(values), 0),
    pseudoinverse = ifelse(values > 0, 1 / values, 0),
    image_projector = as.numeric(values > 0),
    null_projector = as.numeric(values == 0)
  )
  psd_frobenius_norm(transformed)
}

psd_column_norms <- function(X) {
  if (!ncol(X)) return(numeric())
  vapply(seq_len(ncol(X)), function(j) psd_frobenius_norm(X[, j]), numeric(1L))
}
