# Reusable restart-state records and the version-1 Hermitian Lanczos adapter.

#' @keywords internal
restart_state_error <- function(code, field = NULL, expected = NULL,
                                actual = NULL, message = NULL) {
  if (is.null(message)) {
    message <- paste0(
      "Invalid eigencore restart state (", code, ")",
      if (is.null(field)) "." else paste0(": ", field, ".")
    )
  }
  stop(structure(
    list(
      message = message,
      call = NULL,
      code = code,
      field = field,
      expected = expected,
      actual = actual
    ),
    class = c("eigencore_restart_state_error", "error", "condition")
  ))
}

#' @keywords internal
required_restart_state_fields <- function() {
  c(
    "schema_version", "basis", "operator_identity", "problem_signature",
    "requested", "method", "method_state", "provenance", "serialization",
    "memory"
  )
}

#' @keywords internal
restart_target_family <- function(problem) {
  kind <- problem$target$kind %||% target_label(problem$target)
  if (identical(problem$type, "svd")) {
    if (kind %in% c("nearest", "interval")) "svd_interior" else "svd_extremal"
  } else if (!identical(problem$structure$kind, "hermitian")) {
    "general_eigen"
  } else if (kind %in% c("nearest", "interval")) {
    "hermitian_interior"
  } else {
    "hermitian_extremal"
  }
}

#' @keywords internal
restart_coordinate_id <- function(problem, side) {
  op <- if (identical(side, "B") && !is.null(problem$metric)) {
    problem$metric
  } else {
    problem$A
  }
  stable_raw_hash(list(
    schema_version = 1L,
    problem_type = problem$type,
    dim = as.integer(op$dim),
    dtype = op$dtype,
    structure = op$structure$kind,
    metric_present = !is.null(problem$metric),
    side = side
  ))
}

#' @keywords internal
new_problem_signature <- function(plan, basis_sides = NULL) {
  problem <- plan$problem
  if (is.null(basis_sides)) {
    basis_sides <- if (identical(problem$type, "eigen")) "eigen" else character()
  }
  structure(list(
    schema_version = 1L,
    problem_type = problem$type,
    dim = as.integer(problem$A$dim),
    dtype = as.character(problem$A$dtype),
    structure = as.character(problem$structure$kind),
    target_family = restart_target_family(problem),
    metric_present = !is.null(problem$metric),
    coordinate_id_A = restart_coordinate_id(problem, "A"),
    coordinate_id_B = if (is.null(problem$metric)) {
      NULL
    } else {
      restart_coordinate_id(problem, "B")
    },
    basis_sides = as.character(basis_sides),
    basis_metric = if (is.null(problem$metric)) "euclidean" else "B"
  ), class = "eigencore_problem_signature")
}

#' @keywords internal
restart_basis_q <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (is.null(x)) {
    restart_state_error(
      "corrupt_state", "basis", "a non-empty finite numeric matrix", NULL
    )
  }
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || !all(is.finite(x))) {
    restart_state_error(
      "corrupt_state", "basis",
      "a non-empty finite numeric matrix", dim(x)
    )
  }
  q <- matrix(x[0], nrow = nrow(x), ncol = 0L)
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    scale <- max(Mod(v))
    if (!is.finite(scale) || scale == 0) {
      next
    }
    v <- v / scale
    if (ncol(q)) {
      v <- v - q %*% crossprod(q, v)
      v <- v - q %*% crossprod(q, v)
    }
    norm_v <- sqrt(Re(drop(crossprod(v))))
    if (is.finite(norm_v) && norm_v > tol) {
      q <- cbind(q, v / norm_v)
    }
  }
  if (!ncol(q)) {
    restart_state_error(
      "corrupt_state", "basis", "positive numerical rank", 0L
    )
  }
  q
}

#' @keywords internal
restart_basis_error <- function(q) {
  gram <- crossprod(as.matrix(q))
  max(Mod(gram - diag(ncol(q))))
}

#' @keywords internal
restart_original_metric_basis <- function(x) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || !all(is.finite(x))) {
    restart_state_error(
      "corrupt_state", "basis", "a non-empty finite numeric matrix", dim(x)
    )
  }
  rank <- ncol(restart_basis_q(x))
  if (rank != ncol(x)) {
    restart_state_error(
      "corrupt_state", "basis", "full numerical column rank", rank
    )
  }
  x
}

#' @keywords internal
deterministic_restart_matrix <- function(n, p, token, offset = 0L) {
  if (p < 1L) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  token_value <- sum(utf8ToInt(substr(token %||% "restart", 1L, 32L))) + offset
  i <- seq_len(n) + (token_value %% 101L)
  j <- seq_len(p) + (token_value %% 37L)
  outer(i, j, function(ii, jj) {
    sin(ii * jj * sqrt(2)) + cos((ii + jj) * sqrt(3))
  })
}

#' @keywords internal
lanczos_restart_adapter_supported <- function(plan) {
  identical(plan$problem_type, "eigen") &&
    identical(plan$problem$A$dtype, "double") &&
    identical(plan$problem$structure$kind, "hermitian") &&
    is.null(plan$problem$metric) &&
    !is_transform_method(plan$problem$transform) &&
    warm_start_plan_consumes_start(plan$problem, plan)
}

#' @keywords internal
lanczos_restart_start_width <- function(plan) {
  if (plan_dispatches_native_lanczos(plan)) {
    as.integer(plan$controls$block %||% 1L)
  } else {
    1L
  }
}

#' @keywords internal
lanczos_restart_controls_token <- function(plan) {
  stable_raw_hash(list(
    schema_version = 1L,
    requested = plan$requested,
    planned_method = plan$planned_method,
    width = lanczos_restart_start_width(plan),
    block = plan$controls$block %||% 1L,
    reorthogonalize = plan$controls$reorthogonalize %||% TRUE
  ))
}

#' @keywords internal
lanczos_restart_target_token <- function(plan) {
  stable_raw_hash(list(schema_version = 1L, target = plan$problem$target))
}

#' @keywords internal
fit_restart_basis <- function(basis, plan, require_exploration = FALSE) {
  q <- restart_basis_q(basis)
  n <- nrow(q)
  width <- lanczos_restart_start_width(plan)
  if (width > n) {
    restart_state_error(
      "method_incompatible", "method_state$payload$start_block",
      paste0("at most ", n, " columns"), width
    )
  }
  token <- stable_raw_hash(list(
    q = q,
    controls = lanczos_restart_controls_token(plan),
    target = lanczos_restart_target_token(plan)
  ))
  coefficient <- deterministic_restart_matrix(ncol(q), width, token, 11L)
  start <- restart_basis_q(q %*% coefficient)
  if (ncol(start) < width) {
    filler <- deterministic_restart_matrix(n, width, token, 29L)
    for (j in seq_len(ncol(filler))) {
      candidate <- filler[, j, drop = FALSE]
      against <- restart_basis_q_or_empty(cbind(q, start))
      if (ncol(against)) {
        candidate <- candidate - against %*% crossprod(against, candidate)
        candidate <- candidate - against %*% crossprod(against, candidate)
      }
      norm_candidate <- sqrt(Re(drop(crossprod(candidate))))
      if (is.finite(norm_candidate) && norm_candidate > sqrt(.Machine$double.eps)) {
        start <- cbind(start, candidate / norm_candidate)
      }
      if (ncol(start) == width) break
    }
  }
  if (ncol(start) < width) {
    restart_state_error(
      "method_incompatible", "basis",
      paste0("enough rank to fit a ", width, "-column Lanczos start"),
      ncol(start)
    )
  }
  start <- start[, seq_len(width), drop = FALSE]

  exploration <- deterministic_restart_matrix(n, width, token, 53L)
  exploration <- exploration - q %*% crossprod(q, exploration)
  exploration <- restart_basis_q_or_empty(exploration)
  if (!ncol(exploration) && isTRUE(require_exploration)) {
    return(NULL)
  }
  if (ncol(exploration)) {
    mixing <- deterministic_restart_matrix(
      ncol(exploration), width, token, 71L
    )
    mixed <- restart_basis_q(start + 0.05 * exploration %*% mixing)
    if (ncol(mixed) == width) {
      start <- mixed
    }
  }
  list(
    start = start,
    public_basis_columns = seq_len(ncol(q)),
    compression = list(
      source_rank = ncol(q),
      start_width = width,
      compressed = ncol(q) > width,
      augmented = max(0L, width - ncol(q)),
      exploration_columns = ncol(exploration),
      deterministic_token = token
    )
  )
}

#' @keywords internal
restart_basis_q_or_empty <- function(x, tol = sqrt(.Machine$double.eps)) {
  x <- as.matrix(x)
  if (!ncol(x)) return(x)
  tryCatch(
    restart_basis_q(x, tol = tol),
    eigencore_restart_state_error = function(e) {
      matrix(x[0], nrow = nrow(x), ncol = 0L)
    }
  )
}

#' @keywords internal
method_state_token_payload <- function(method_state) {
  payload <- deep_copy_record(method_state)
  payload$payload$integrity_token <- NULL
  payload
}

#' @keywords internal
new_lanczos_method_state <- function(basis, plan) {
  if (!lanczos_restart_adapter_supported(plan)) {
    return(NULL)
  }
  fitted <- fit_restart_basis(basis, plan, require_exploration = TRUE)
  if (is.null(fitted)) {
    return(NULL)
  }
  identities_portable <- all(vapply(
    plan$operator_identity, function(x) isTRUE(x$portable), logical(1L)
  ))
  out <- structure(list(
    kind = "hermitian_lanczos_start_block",
    schema_version = 1L,
    adapter_version = 1L,
    problem_type = "eigen",
    method_family = "standard_real_hermitian_lanczos",
    controls_token = lanczos_restart_controls_token(plan),
    target_token = lanczos_restart_target_token(plan),
    operator_identity_token = stable_raw_hash(plan$operator_identity),
    portable = identities_portable,
    payload = list(
      start_block = fitted$start,
      public_basis_columns = fitted$public_basis_columns,
      compression = fitted$compression,
      integrity_token = NULL
    )
  ), class = "eigencore_method_state")
  out$payload$integrity_token <- stable_raw_hash(method_state_token_payload(out))
  out
}

#' @keywords internal
restart_serialization_without_integrity <- function(serialization) {
  out <- serialization
  out$integrity_token <- NULL
  out
}

#' @keywords internal
restart_state_integrity_payload <- function(state) {
  known <- state[required_restart_state_fields()]
  known$memory <- NULL
  known$serialization <- restart_serialization_without_integrity(
    known$serialization
  )
  known
}

#' @keywords internal
new_restart_memory <- function(state) {
  metadata <- list(
    operator_identity = state$operator_identity,
    problem_signature = state$problem_signature,
    requested = state$requested,
    method = state$method,
    provenance = state$provenance,
    serialization = state$serialization
  )
  sizes <- c(
    basis = as.numeric(utils::object.size(state$basis)),
    method_state = if (is.null(state$method_state)) 0 else {
      as.numeric(utils::object.size(state$method_state))
    },
    cached_operator_actions = 0,
    metadata = as.numeric(utils::object.size(metadata))
  )
  structure(list(
    schema_version = 1L,
    r_bytes = sum(sizes),
    native_bytes = 0,
    total_bytes = sum(sizes),
    complete = TRUE,
    components = as.list(sizes)
  ), class = "eigencore_memory")
}

#' @keywords internal
new_restart_serialization <- function(plan, method_state) {
  identity_portable <- all(vapply(
    plan$operator_identity, function(x) isTRUE(x$portable), logical(1L)
  ))
  structure(list(
    schema_version = 1L,
    portable = identity_portable,
    method_state_portable = if (is.null(method_state)) {
      TRUE
    } else {
      isTRUE(method_state$portable)
    },
    originating_session = eigencore_session_id(),
    integrity_token = NULL
  ), class = "eigencore_serialization")
}

#' @keywords internal
construct_restart_state <- function(result, retention = c("basis", "same_operator")) {
  retention <- match.arg(retention)
  if (!inherits(result, c("eigencore_eigen_result", "eigencore_svd_result"))) {
    restart_state_error(
      "corrupt_state", "class", "certified eigencore eigen or SVD result",
      class(result)
    )
  }
  plan <- result$plan %||% NULL
  validate_eigencore_plan(plan)
  certificate <- result$certificate %||% NULL
  if (is.null(certificate) || !isTRUE(certificate$passed)) {
    restart_state_error(
      "corrupt_state", "certificate$passed", TRUE,
      certificate$passed %||% NULL,
      "A reusable restart state requires a result with a passed current-operator certificate."
    )
  }

  if (inherits(result, "eigencore_eigen_result")) {
    basis <- if (is.null(plan$problem$metric)) {
      restart_basis_q(result$vectors %||% NULL)
    } else {
      # A passed generalized certificate already established B-orthonormality.
      # Preserve those original-coordinate vectors: Euclidean QR would destroy
      # the public basis_metric = "B" invariant.
      restart_original_metric_basis(result$vectors %||% NULL)
    }
    if (nrow(basis) != plan$problem$A$dim[[1L]]) {
      restart_state_error(
        "corrupt_state", "basis", plan$problem$A$dim[[1L]], nrow(basis)
      )
    }
    basis_sides <- "eigen"
  } else {
    left <- result$u %||% NULL
    right <- result$v %||% NULL
    if (!is.null(left)) left <- restart_basis_q(left)
    if (!is.null(right)) right <- restart_basis_q(right)
    if (is.null(left) && is.null(right)) {
      restart_state_error(
        "corrupt_state", "basis", "at least one computed SVD vector side", NULL
      )
    }
    if (!is.null(left) && nrow(left) != plan$problem$A$dim[[1L]]) {
      restart_state_error(
        "corrupt_state", "basis$left", plan$problem$A$dim[[1L]], nrow(left)
      )
    }
    if (!is.null(right) && nrow(right) != plan$problem$A$dim[[2L]]) {
      restart_state_error(
        "corrupt_state", "basis$right", plan$problem$A$dim[[2L]], nrow(right)
      )
    }
    basis <- structure(
      list(schema_version = 1L, left = left, right = right),
      class = "eigencore_svd_basis"
    )
    basis_sides <- c(if (!is.null(left)) "left", if (!is.null(right)) "right")
  }

  method_state <- if (
    identical(retention, "same_operator") &&
      inherits(result, "eigencore_eigen_result") &&
      identical(result$actual_method %||% result$method, plan$planned_method)
  ) {
    new_lanczos_method_state(basis, plan)
  } else {
    NULL
  }
  provenance <- structure(list(
    schema_version = 1L,
    retention = retention,
    source_method = result$actual_method %||% result$method,
    source_relation = result$state_transition$relation %||% "cold_start",
    certificate_passed = TRUE,
    method_state_available = !is.null(method_state),
    basis_rank = if (inherits(basis, "eigencore_svd_basis")) {
      c(
        left = if (is.null(basis$left)) 0L else ncol(basis$left),
        right = if (is.null(basis$right)) 0L else ncol(basis$right)
      )
    } else {
      ncol(basis)
    }
  ), class = "eigencore_restart_provenance")
  state <- structure(list(
    schema_version = 1L,
    basis = basis,
    operator_identity = deep_copy_record(plan$operator_identity),
    problem_signature = new_problem_signature(plan, basis_sides),
    requested = as.integer(plan$requested),
    method = as.character(result$actual_method %||% result$method),
    method_state = method_state,
    provenance = provenance,
    serialization = new_restart_serialization(plan, method_state),
    memory = NULL
  ), class = "eigencore_restart_state")
  state$serialization$integrity_token <- stable_raw_hash(
    restart_state_integrity_payload(state)
  )
  state$memory <- new_restart_memory(state)
  deep_copy_record(state)
}

#' @keywords internal
validate_problem_signature <- function(signature) {
  required <- c(
    "schema_version", "problem_type", "dim", "dtype", "structure",
    "target_family", "metric_present", "coordinate_id_A", "coordinate_id_B",
    "basis_sides", "basis_metric"
  )
  if (!inherits(signature, "eigencore_problem_signature") ||
      !identical(signature$schema_version, 1L)) {
    restart_state_error(
      "unsupported_schema", "problem_signature$schema_version", 1L,
      signature$schema_version %||% NULL
    )
  }
  missing <- setdiff(required, names(signature))
  if (length(missing)) {
    restart_state_error("corrupt_state", paste0("problem_signature$", missing[[1L]]),
                        "present", NULL)
  }
  invisible(signature)
}

#' @keywords internal
validate_restart_identity <- function(identity, field) {
  if (!inherits(identity, "eigencore_operator_identity") ||
      !identical(identity$schema_version, 1L)) {
    restart_state_error(
      "unsupported_schema", paste0(field, "$schema_version"), 1L,
      identity$schema_version %||% NULL
    )
  }
  required <- c(
    "operator_id", "revision", "origin", "dim", "dtype", "structure",
    "portable", "session_id"
  )
  missing <- setdiff(required, names(identity))
  if (length(missing)) {
    restart_state_error(
      "corrupt_state", paste0(field, "$", missing[[1L]]), "present", NULL
    )
  }
  invisible(identity)
}

#' @keywords internal
validate_restart_basis <- function(basis, signature) {
  if (identical(signature$problem_type, "eigen")) {
    if (!identical(signature$basis_sides, "eigen")) {
      restart_state_error(
        "corrupt_state", "problem_signature$basis_sides", "eigen",
        signature$basis_sides
      )
    }
    if (!is.matrix(basis) || !is.numeric(basis) || !all(is.finite(basis)) ||
        ncol(basis) < 1L || nrow(basis) != signature$dim[[1L]]) {
      restart_state_error(
        "corrupt_state", "basis",
        paste0(signature$dim[[1L]], "-row finite eigen basis"), dim(basis)
      )
    }
    if (identical(signature$basis_metric, "euclidean") &&
        restart_basis_error(basis) > 1e-7) {
      restart_state_error(
        "corrupt_state", "basis", "orthonormal columns",
        restart_basis_error(basis)
      )
    }
    if (identical(signature$basis_metric, "B") &&
        ncol(restart_basis_q(basis)) != ncol(basis)) {
      restart_state_error(
        "corrupt_state", "basis", "full numerical column rank", ncol(basis)
      )
    }
  } else {
    if (!inherits(basis, "eigencore_svd_basis") ||
        !identical(basis$schema_version, 1L)) {
      restart_state_error(
        "unsupported_schema", "basis$schema_version", 1L,
        basis$schema_version %||% NULL
      )
    }
    if (is.null(basis$left) && is.null(basis$right)) {
      restart_state_error(
        "corrupt_state", "basis", "at least one SVD side", NULL
      )
    }
    actual_sides <- c(
      if (!is.null(basis$left)) "left",
      if (!is.null(basis$right)) "right"
    )
    if (!identical(signature$basis_sides, actual_sides)) {
      restart_state_error(
        "corrupt_state", "problem_signature$basis_sides", actual_sides,
        signature$basis_sides
      )
    }
    for (side in c("left", "right")) {
      value <- basis[[side]]
      if (is.null(value)) next
      expected_rows <- signature$dim[[if (identical(side, "left")) 1L else 2L]]
      if (!is.matrix(value) || !is.numeric(value) || !all(is.finite(value)) ||
          nrow(value) != expected_rows || ncol(value) < 1L ||
          restart_basis_error(value) > 1e-7) {
        restart_state_error(
          "corrupt_state", paste0("basis$", side),
          paste0(expected_rows, "-row finite orthonormal basis"), dim(value)
        )
      }
    }
  }
  invisible(basis)
}

#' @keywords internal
validate_restart_state <- function(state) {
  if (!inherits(state, "eigencore_restart_state")) {
    restart_state_error(
      "corrupt_state", "class", "eigencore_restart_state", class(state)
    )
  }
  if (!identical(state$schema_version, 1L)) {
    restart_state_error(
      "unsupported_schema", "schema_version", 1L, state$schema_version
    )
  }
  missing <- setdiff(required_restart_state_fields(), names(state))
  if (length(missing)) {
    restart_state_error("corrupt_state", missing[[1L]], "present", NULL)
  }
  validate_problem_signature(state$problem_signature)
  validate_restart_basis(state$basis, state$problem_signature)
  if (!is.list(state$operator_identity) ||
      !"A" %in% names(state$operator_identity)) {
    restart_state_error(
      "corrupt_state", "operator_identity", "named identity list containing A",
      names(state$operator_identity)
    )
  }
  for (name in names(state$operator_identity)) {
    validate_restart_identity(
      state$operator_identity[[name]], paste0("operator_identity$", name)
    )
  }
  if (!is.numeric(state$requested) || length(state$requested) != 1L ||
      is.na(state$requested) || state$requested < 1L) {
    restart_state_error(
      "corrupt_state", "requested", "one positive integer", state$requested
    )
  }
  if (!is.character(state$method) || length(state$method) != 1L ||
      is.na(state$method) || !nzchar(state$method)) {
    restart_state_error(
      "corrupt_state", "method", "one non-empty method label", state$method
    )
  }
  if (!is.null(state$method_state) && !is.list(state$method_state)) {
    restart_state_error(
      "corrupt_state", "method_state", "NULL or a versioned list",
      class(state$method_state)
    )
  }
  if (!inherits(state$serialization, "eigencore_serialization") ||
      !identical(state$serialization$schema_version, 1L)) {
    restart_state_error(
      "unsupported_schema", "serialization$schema_version", 1L,
      state$serialization$schema_version %||% NULL
    )
  }
  expected_token <- stable_raw_hash(restart_state_integrity_payload(state))
  if (!identical(state$serialization$integrity_token, expected_token)) {
    restart_state_error(
      "corrupt_state", "serialization$integrity_token", expected_token,
      state$serialization$integrity_token,
      "The restart state was mutated or corrupted after construction."
    )
  }
  if (!inherits(state$memory, "eigencore_memory") ||
      !identical(state$memory$schema_version, 1L)) {
    restart_state_error(
      "corrupt_state", "memory", "eigencore_memory schema version 1",
      class(state$memory)
    )
  }
  expected_memory <- new_restart_memory(state)
  if (!identical(unclass(state$memory), unclass(expected_memory))) {
    restart_state_error(
      "corrupt_state", "memory", unclass(expected_memory), unclass(state$memory)
    )
  }
  deep_copy_record(state)
}

#' Extract a reusable spectral restart state.
#'
#' A restart state is an immutable acceleration hint built only from a result
#' with a passed current-operator certificate. Its public basis is expressed in
#' original problem coordinates. `retention = "same_operator"` may additionally
#' retain a versioned method-fitted start block when the producing route has an
#' admitted adapter; it never retains convergence, locks, cached operator
#' actions, or a certificate.
#'
#' @param x A certified eigencore eigen/SVD result or an existing restart state.
#' @param retention Retain only the public basis, or also request eligible
#'   same-operator method state.
#' @return An `eigencore_restart_state` schema-version-1 record.
#' @export
restart_state <- function(x, retention = c("basis", "same_operator")) {
  retention <- match.arg(retention)
  if (inherits(x, "eigencore_restart_state")) {
    state <- validate_restart_state(x)
    if (identical(retention, "basis") && !is.null(state$method_state)) {
      state$method_state <- NULL
      state$provenance$retention <- "basis"
      state$provenance$method_state_available <- FALSE
      state$serialization$method_state_portable <- TRUE
      state$serialization$integrity_token <- stable_raw_hash(
        restart_state_integrity_payload(state)
      )
      state$memory <- new_restart_memory(state)
    }
    return(deep_copy_record(state))
  }
  construct_restart_state(x, retention = retention)
}

#' Return retained-memory accounting.
#'
#' @param x An eigencore plan, restart state, result, or `eigencore_memory`
#'   record.
#' @return The known retained byte count. The full `eigencore_memory` record is
#'   attached as the `memory` attribute; when `complete` is false the value is a
#'   documented lower bound.
#' @export
retained_bytes <- function(x) {
  memory <- if (inherits(x, "eigencore_memory")) x else x$memory %||% NULL
  if (!inherits(memory, "eigencore_memory")) {
    stop("retained_bytes() requires an object with eigencore memory metadata.",
         call. = FALSE)
  }
  value <- as.numeric(memory$total_bytes)
  attr(value, "memory") <- deep_copy_record(memory)
  value
}

#' @keywords internal
coordinate_signature_fields <- function() {
  c(
    "problem_type", "dim", "dtype", "structure", "metric_present",
    "coordinate_id_A", "coordinate_id_B", "basis_metric"
  )
}

#' @keywords internal
validate_restart_coordinate_compatibility <- function(state, plan) {
  destination <- new_problem_signature(
    plan, basis_sides = state$problem_signature$basis_sides
  )
  for (field in coordinate_signature_fields()) {
    if (!identical(state$problem_signature[[field]], destination[[field]])) {
      restart_state_error(
        "coordinate_incompatible", paste0("problem_signature$", field),
        destination[[field]], state$problem_signature[[field]]
      )
    }
  }
  destination
}

#' @keywords internal
restart_operator_relation <- function(source, destination) {
  if (!identical(sort(names(source)), sort(names(destination)))) {
    restart_state_error(
      "coordinate_incompatible", "operator_identity",
      sort(names(destination)), sort(names(source))
    )
  }
  source <- source[sort(names(source))]
  destination <- destination[sort(names(destination))]
  same_id <- vapply(seq_along(source), function(i) {
    identical(source[[i]]$operator_id, destination[[i]]$operator_id)
  }, logical(1L))
  same_revision <- vapply(seq_along(source), function(i) {
    identical(source[[i]]$revision, destination[[i]]$revision) &&
      identical(source[[i]]$dim, destination[[i]]$dim) &&
      identical(source[[i]]$dtype, destination[[i]]$dtype) &&
      identical(source[[i]]$structure, destination[[i]]$structure)
  }, logical(1L))
  if (all(same_id) && all(same_revision)) {
    "same_operator"
  } else if (all(same_id)) {
    "changed_revision"
  } else {
    "changed_operator"
  }
}

#' @keywords internal
validate_restart_identity_session <- function(state) {
  nonportable <- Filter(
    function(x) !isTRUE(x$portable), state$operator_identity
  )
  if (length(nonportable) && any(vapply(
    nonportable,
    function(x) !identical(x$session_id, eigencore_session_id()),
    logical(1L)
  ))) {
    restart_state_error(
      "session_incompatible", "operator_identity", eigencore_session_id(),
      vapply(nonportable, `[[`, character(1L), "session_id")
    )
  }
  invisible(state)
}

#' @keywords internal
validate_lanczos_method_state <- function(method_state, state, plan) {
  if (is.null(method_state)) {
    restart_state_error(
      "method_incompatible", "method_state",
      "a retained Hermitian Lanczos start-block payload", NULL
    )
  }
  required <- c(
    "kind", "schema_version", "adapter_version", "problem_type",
    "method_family", "controls_token", "target_token",
    "operator_identity_token", "portable", "payload"
  )
  missing <- setdiff(required, names(method_state))
  if (length(missing)) {
    restart_state_error(
      "stale_method_state", paste0("method_state$", missing[[1L]]),
      "present", NULL
    )
  }
  if (!is.list(method_state$payload) ||
      restart_method_payload_has_forbidden_object(method_state$payload)) {
    restart_state_error(
      "stale_method_state", "method_state$payload",
      "ordinary serializable R data without functions, environments, or external pointers",
      typeof(method_state$payload)
    )
  }
  if (!identical(method_state$schema_version, 1L) ||
      !identical(method_state$adapter_version, 1L) ||
      !identical(method_state$kind, "hermitian_lanczos_start_block")) {
    restart_state_error(
      "stale_method_state", "method_state$schema_version",
      "Hermitian Lanczos start-block adapter version 1",
      list(
        kind = method_state$kind %||% NULL,
        schema_version = method_state$schema_version %||% NULL,
        adapter_version = method_state$adapter_version %||% NULL
      )
    )
  }
  expected <- list(
    problem_type = "eigen",
    method_family = "standard_real_hermitian_lanczos",
    controls_token = lanczos_restart_controls_token(plan),
    target_token = lanczos_restart_target_token(plan),
    operator_identity_token = stable_raw_hash(plan$operator_identity)
  )
  for (field in names(expected)) {
    if (!identical(method_state[[field]], expected[[field]])) {
      restart_state_error(
        if (field %in% c("controls_token", "target_token")) {
          "method_incompatible"
        } else {
          "stale_method_state"
        },
        paste0("method_state$", field), expected[[field]], method_state[[field]]
      )
    }
  }
  if (!isTRUE(method_state$portable)) {
    restart_state_error(
      "session_incompatible", "method_state$portable", TRUE, FALSE
    )
  }
  expected_integrity <- stable_raw_hash(method_state_token_payload(method_state))
  if (!identical(method_state$payload$integrity_token, expected_integrity)) {
    restart_state_error(
      "stale_method_state", "method_state$payload$integrity_token",
      expected_integrity, method_state$payload$integrity_token
    )
  }
  start <- method_state$payload$start_block %||% NULL
  basis_columns <- method_state$payload$public_basis_columns %||% NULL
  if (!is.integer(basis_columns) ||
      !identical(basis_columns, seq_len(ncol(state$basis)))) {
    restart_state_error(
      "stale_method_state", "method_state$payload$public_basis_columns",
      seq_len(ncol(state$basis)), basis_columns
    )
  }
  expected_dim <- c(plan$problem$A$dim[[1L]], lanczos_restart_start_width(plan))
  if (!is.matrix(start) || !is.numeric(start) || !all(is.finite(start)) ||
      !identical(dim(start), as.integer(expected_dim)) ||
      restart_basis_error(start) > 1e-7) {
    restart_state_error(
      "stale_method_state", "method_state$payload$start_block",
      expected_dim, dim(start)
    )
  }
  deep_copy_record(start)
}

#' @keywords internal
restart_method_payload_has_forbidden_object <- function(x) {
  if (is.function(x) || is.environment(x) ||
      typeof(x) %in% c("externalptr", "weakref")) {
    return(TRUE)
  }
  if (!is.list(x)) {
    return(FALSE)
  }
  any(vapply(x, restart_method_payload_has_forbidden_object, logical(1L)))
}

#' @keywords internal
restart_prepared_provenance <- function(state, fitted, method_state_used) {
  rank <- ncol(state$basis)
  list(
    start_source = if (method_state_used) {
      "restart_state_method_state"
    } else {
      "restart_state_basis"
    },
    supplied = rank,
    accepted = rank,
    rejected = 0L,
    augmented = fitted$compression$augmented %||% 0L,
    rank = rank,
    compressed = isTRUE(fitted$compression$compressed),
    invariant_guard_used = FALSE,
    invariant_relative_residual = NA_real_,
    guard_operator_block_calls = 0L,
    guard_operator_columns = 0L,
    restart_state = TRUE
  )
}

#' @keywords internal
prepare_restart_state_for_plan <- function(state, plan, reuse) {
  if (is.null(state)) {
    relation <- if (is.null(plan$execution$initial_subspace)) {
      "cold_start"
    } else {
      "initial_subspace"
    }
    return(list(
      state_supplied = FALSE,
      prepared_start = NULL,
      transition = new_state_transition(
        relation = relation,
        reason = new_state_reason(
          relation,
          if (identical(relation, "cold_start")) {
            "no restart state supplied"
          } else {
            "plan contains an initial_subspace basis hint"
          }
        ),
        reuse = reuse
      )
    ))
  }
  if (!is.null(plan$execution$initial_subspace)) {
    restart_state_error(
      "method_incompatible", "initial_subspace", NULL,
      "plan already contains initial_subspace",
      "A plan cannot consume both initial_subspace and restart_state hints."
    )
  }
  state <- validate_restart_state(state)
  validate_restart_coordinate_compatibility(state, plan)
  validate_restart_identity_session(state)
  relation <- restart_operator_relation(
    state$operator_identity, plan$operator_identity
  )
  if (!lanczos_restart_adapter_supported(plan) ||
      !identical(state$problem_signature$problem_type, "eigen")) {
    restart_state_error(
      "method_incompatible", "planned_method",
      "standard real Hermitian Lanczos restart-basis adapter",
      plan$planned_method
    )
  }
  if (identical(reuse, "same_operator") &&
      !identical(relation, "same_operator")) {
    restart_state_error(
      "operator_incompatible", "operator_identity",
      state$operator_identity, plan$operator_identity
    )
  }

  invalidated <- c(
    "recurrence", "locked", "cached_operator_actions", "projection",
    "residuals", "convergence", "ordering", "certificate"
  )
  method_state_used <- FALSE
  adapter <- NULL
  downgrade_reason <- NULL
  fitted <- NULL
  start <- NULL
  if (identical(reuse, "same_operator")) {
    start <- validate_lanczos_method_state(state$method_state, state, plan)
    fitted <- list(
      start = start,
      compression = state$method_state$payload$compression
    )
    method_state_used <- TRUE
  } else if (identical(reuse, "auto") &&
             identical(relation, "same_operator") &&
             !is.null(state$method_state)) {
    start <- tryCatch(
      validate_lanczos_method_state(state$method_state, state, plan),
      eigencore_restart_state_error = function(e) {
        downgrade_reason <<- new_state_reason(
          e$code,
          paste0("same-operator method state was invalidated: ", e$message),
          list(field = e$field)
        )
        NULL
      }
    )
    if (!is.null(start)) {
      fitted <- list(
        start = start,
        compression = state$method_state$payload$compression
      )
      method_state_used <- TRUE
    }
  }
  if (is.null(start)) {
    fitted <- fit_restart_basis(state$basis, plan, require_exploration = FALSE)
    start <- fitted$start
    invalidated <- unique(c("method_state", invalidated))
  }
  adapter <- list(
    kind = "hermitian_lanczos_start_block",
    schema_version = 1L,
    adapter_version = 1L
  )
  reason <- downgrade_reason %||% new_state_reason(
    if (method_state_used) "method_state_accepted" else paste0(relation, "_basis"),
    if (method_state_used) {
      "validated same-operator Lanczos start block accepted"
    } else {
      paste0("public restart basis accepted for ", relation, " reuse")
    }
  )
  transition <- new_state_transition(
    relation = relation,
    basis_used = TRUE,
    method_state_used = method_state_used,
    invalidated = invalidated,
    reason = reason,
    source_operator_identity = state$operator_identity,
    destination_operator_identity = plan$operator_identity,
    reuse = reuse,
    adapter = adapter
  )
  list(
    state_supplied = TRUE,
    state = state,
    prepared_start = list(
      start = start,
      guard_basis = start,
      provenance = restart_prepared_provenance(
        state, fitted, method_state_used
      ),
      strict_same_operator = identical(reuse, "same_operator"),
      method_state_candidate = method_state_used
    ),
    transition = transition
  )
}

#' @keywords internal
finalize_restart_transition <- function(result, preparation) {
  transition <- preparation$transition
  if (!isTRUE(preparation$state_supplied)) {
    return(transition)
  }
  source <- result$start_source %||% ""
  consumed <- startsWith(source, "restart_state_") &&
    !grepl("discarded", source, fixed = TRUE)
  transition$basis_used <- consumed
  transition$method_state_used <- consumed &&
    isTRUE(preparation$prepared_start$method_state_candidate)
  if (!consumed) {
    transition$adapter <- NULL
    transition$reason <- new_state_reason(
      "restart_start_discarded",
      paste0("restart state was validated but the fitted start was not consumed: ",
             source %||% "unknown reason")
    )
  }
  transition
}

#' @keywords internal
result_memory_record <- function(result) {
  new_memory_record(
    list(
      metadata = list(
        planned_method = result$planned_method,
        actual_method = result$actual_method,
        state_transition = result$state_transition
      ),
      restart_state = result$restart_state
    )
  )
}
