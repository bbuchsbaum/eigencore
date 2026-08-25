# Executable-plan and reusable-workflow records introduced in eigencore 1.2.

#' @keywords internal
eigencore_session_id <- local({
  id <- NULL
  function() {
    if (is.null(id)) {
      payload <- serialize(
        list(pid = Sys.getpid(), started = Sys.time(), tempdir = tempdir()),
        NULL,
        version = 3L
      )
      id <<- paste0("r-session-", stable_raw_hash(payload))
    }
    id
  }
})

#' @keywords internal
stable_raw_hash <- function(x) {
  if (!is.raw(x)) {
    x <- serialize(x, NULL, version = 3L)
  }
  .Call("eigencore_stable_raw_hash", x, PACKAGE = "eigencore")
}

#' @keywords internal
deep_copy_record <- function(x) {
  unserialize(serialize(x, NULL, version = 3L))
}

#' @keywords internal
new_operator_identity <- function(operator_id, revision, origin, dim, dtype,
                                  structure, portable, session_id = NULL) {
  structure(
    list(
      schema_version = 1L,
      operator_id = operator_id,
      revision = revision,
      origin = origin,
      dim = as.integer(dim),
      dtype = as.character(dtype),
      structure = as.character(structure),
      portable = isTRUE(portable),
      session_id = if (isTRUE(portable)) NULL else session_id %||% eigencore_session_id()
    ),
    class = "eigencore_operator_identity"
  )
}

#' @keywords internal
canonical_identity_value <- function(x) {
  if (inherits(x, "eigencore_operator")) {
    return(list(identity = unclass(refresh_operator_identity(x))))
  }
  if (is.function(x) || is.environment(x) || typeof(x) == "externalptr") {
    return(NULL)
  }
  if (is.list(x)) {
    out <- lapply(x, canonical_identity_value)
    names(out) <- names(x)
    return(out)
  }
  x
}

#' @keywords internal
builtin_operator_identity_payload <- function(op) {
  metadata <- op$metadata %||% list()
  source <- metadata$source %||% metadata$matrix %||% NULL
  structural_metadata <- metadata
  structural_metadata$source <- NULL
  structural_metadata$matrix <- NULL
  structural_metadata$native <- NULL
  structural_metadata$frobenius_norm <- NULL
  list(
    dim = as.integer(op$dim),
    dtype = op$dtype,
    structure = op$structure$kind,
    source_class = class(source),
    source = source,
    metadata = canonical_identity_value(structural_metadata)
  )
}

#' @keywords internal
operator_has_builtin_provenance <- function(metadata) {
  keys <- c(
    "source", "matrix", "parent", "left", "right", "terms",
    "structured_grid_laplacian_2d", "tridiagonal"
  )
  any(keys %in% names(metadata %||% list()))
}

#' @keywords internal
next_opaque_identity_token <- local({
  counter <- 0L
  function(dim, dtype, structure) {
    counter <<- counter + 1L
    stable_raw_hash(list(
      session = eigencore_session_id(),
      counter = counter,
      dim = as.integer(dim),
      dtype = dtype,
      structure = structure
    ))
  }
})

#' @keywords internal
make_operator_identity <- function(dim, dtype, structure, metadata,
                                   operator_id = NULL, revision = NULL,
                                   portable = FALSE) {
  supplied <- !is.null(operator_id) || !is.null(revision)
  if (xor(is.null(operator_id), is.null(revision))) {
    stop("operator_id and revision must be supplied together.", call. = FALSE)
  }
  if (supplied) {
    if (!is.character(operator_id) || length(operator_id) != 1L ||
        is.na(operator_id) || !nzchar(operator_id)) {
      stop("operator_id must be one non-empty character value.", call. = FALSE)
    }
    if (!is.character(revision) || length(revision) != 1L ||
        is.na(revision) || !nzchar(revision)) {
      stop("revision must be one non-empty character value.", call. = FALSE)
    }
    if (!is.logical(portable) || length(portable) != 1L || is.na(portable)) {
      stop("portable must be TRUE or FALSE.", call. = FALSE)
    }
    return(new_operator_identity(
      operator_id = operator_id,
      revision = revision,
      origin = "callback_explicit",
      dim = dim,
      dtype = dtype,
      structure = structure$kind,
      portable = portable
    ))
  }
  if (isTRUE(portable)) {
    stop(
      "portable = TRUE requires explicit operator_id and revision provenance.",
      call. = FALSE
    )
  }
  if (operator_has_builtin_provenance(metadata)) {
    proto <- list(
      dim = as.integer(dim),
      apply = NULL,
      apply_adjoint = NULL,
      dtype = dtype,
      structure = structure,
      metadata = metadata
    )
    digest <- stable_raw_hash(builtin_operator_identity_payload(proto))
    return(new_operator_identity(
      operator_id = paste0("builtin-", digest),
      revision = digest,
      origin = "builtin",
      dim = dim,
      dtype = dtype,
      structure = structure$kind,
      portable = TRUE
    ))
  }
  token <- next_opaque_identity_token(dim, dtype, structure$kind)
  new_operator_identity(
    operator_id = paste0("opaque-", token),
    revision = "opaque-revision-1",
    origin = "callback_opaque",
    dim = dim,
    dtype = dtype,
    structure = structure$kind,
    portable = FALSE
  )
}

#' @keywords internal
refresh_operator_identity <- function(op) {
  identity <- op$identity %||% NULL
  if (!inherits(identity, "eigencore_operator_identity")) {
    stop("Operator is missing eigencore 1.2 identity metadata.", call. = FALSE)
  }
  if (!identical(identity$origin, "builtin")) {
    return(identity)
  }
  digest <- stable_raw_hash(builtin_operator_identity_payload(op))
  new_operator_identity(
    operator_id = paste0("builtin-", digest),
    revision = digest,
    origin = "builtin",
    dim = op$dim,
    dtype = op$dtype,
    structure = op$structure$kind,
    portable = TRUE
  )
}

#' Return operator identity and revision provenance.
#'
#' @param x An eigencore operator, problem, executable plan, restart state, or
#'   result.
#' @return An `eigencore_operator_identity` for an operator, or a named list
#'   with entries `A` and, when present, `B` for larger workflow objects.
#' @export
operator_identity <- function(x) {
  if (inherits(x, "eigencore_operator")) {
    return(refresh_operator_identity(x))
  }
  if (inherits(x, "eigencore_eigen_problem")) {
    out <- list(A = refresh_operator_identity(x$A))
    if (!is.null(x$metric)) {
      out$B <- refresh_operator_identity(x$metric)
    }
    return(out)
  }
  if (inherits(x, "eigencore_svd_problem")) {
    return(list(A = refresh_operator_identity(x$A)))
  }
  if (inherits(x, "eigencore_plan")) {
    return(x$operator_identity %||% NULL)
  }
  if (inherits(x, c("eigencore_eigen_result", "eigencore_svd_result"))) {
    return(operator_identity(x$plan))
  }
  if (inherits(x, "eigencore_restart_state")) {
    return(x$operator_identity %||% NULL)
  }
  stop("operator_identity() does not support class ",
       paste(class(x), collapse = "/"), ".", call. = FALSE)
}

#' @keywords internal
normalize_scalar_option <- function(name, default, mode = c("numeric", "integer", "logical"),
                                    min = -Inf, max = Inf, allow_infinite = TRUE,
                                    invalid = c("default", "error")) {
  mode <- match.arg(mode)
  invalid <- match.arg(invalid)
  value <- getOption(name, default)
  valid <- length(value) == 1L && !is.na(value)
  if (mode == "logical") {
    valid <- valid && is.logical(value)
  } else {
    valid <- valid && is.numeric(value) && value >= min && value <= max &&
      (allow_infinite || is.finite(value))
  }
  if (!valid) {
    if (invalid == "error") {
      stop("Option ", name, " has an invalid value.", call. = FALSE)
    }
    value <- default
  }
  switch(
    mode,
    integer = if (is.finite(value)) as.integer(value) else Inf,
    numeric = as.numeric(value),
    logical = isTRUE(value)
  )
}

#' @keywords internal
planner_policy_keys <- function() {
  c(
    "eigencore.arnoldi_max_restarts",
    "eigencore.block_dense_full_subspace_max_n",
    "eigencore.dense_partial_lanczos_max_fraction",
    "eigencore.dense_partial_lanczos_min_n",
    "eigencore.promote_sparse_block_lanczos",
    "eigencore.lobpcg_maxit",
    "eigencore.dense_fallback_mb",
    "eigencore.max_dense_fallback_bytes",
    "eigencore.gram_svd_max_dimension",
    "eigencore.gram_svd_max_dimension_wide",
    "eigencore.gram_svd_memory_mb",
    "eigencore.gram_svd_rank_fraction_limit",
    "eigencore.gram_svd_min_aspect_ratio",
    "eigencore.gram_svd_work_budget",
    "eigencore.implicit_gram_svd_min_dimension",
    "eigencore.promote_retained_golub_kahan",
    "eigencore.golub_kahan_projected_stop",
    "eigencore.golub_kahan_prefix_diagnostics",
    "eigencore.randomized_stage_timing",
    "eigencore.randomized_adaptive_stop"
  )
}

#' @keywords internal
planner_policy_snapshot <- function() {
  dense_mb <- getOption("eigencore.dense_fallback_mb", NULL)
  legacy_bytes <- getOption("eigencore.max_dense_fallback_bytes", NULL)
  if (is.null(dense_mb)) {
    if (is.null(legacy_bytes)) {
      dense_mb <- 256
      legacy_bytes <- 256e6
    } else {
      if (!is.numeric(legacy_bytes) || length(legacy_bytes) != 1L ||
          is.na(legacy_bytes) || legacy_bytes < 0) {
        stop("Option eigencore.max_dense_fallback_bytes must be one non-negative numeric value.",
             call. = FALSE)
      }
      dense_mb <- as.numeric(legacy_bytes) / 1e6
    }
  } else {
    if (!is.numeric(dense_mb) || length(dense_mb) != 1L ||
        is.na(dense_mb) || dense_mb < 0) {
      stop("Option eigencore.dense_fallback_mb must be one non-negative numeric value.",
           call. = FALSE)
    }
    dense_mb <- as.numeric(dense_mb)
    legacy_bytes <- dense_mb * 1e6
  }
  values <- list(
    "eigencore.arnoldi_max_restarts" = normalize_scalar_option(
      "eigencore.arnoldi_max_restarts", 5L, "integer", min = 0
    ),
    "eigencore.block_dense_full_subspace_max_n" = normalize_scalar_option(
      "eigencore.block_dense_full_subspace_max_n", 256L, "integer", min = 1
    ),
    "eigencore.dense_partial_lanczos_max_fraction" = normalize_scalar_option(
      "eigencore.dense_partial_lanczos_max_fraction", 0.25, "numeric",
      min = .Machine$double.eps, max = 1
    ),
    "eigencore.dense_partial_lanczos_min_n" = normalize_scalar_option(
      "eigencore.dense_partial_lanczos_min_n", 128L, "integer", min = 1
    ),
    "eigencore.promote_sparse_block_lanczos" = isTRUE(
      getOption("eigencore.promote_sparse_block_lanczos", FALSE)
    ),
    "eigencore.lobpcg_maxit" = normalize_scalar_option(
      "eigencore.lobpcg_maxit", 200L, "integer", min = 1
    ),
    "eigencore.dense_fallback_mb" = dense_mb,
    "eigencore.max_dense_fallback_bytes" = as.numeric(legacy_bytes),
    "eigencore.gram_svd_max_dimension" = normalize_scalar_option(
      "eigencore.gram_svd_max_dimension", 512, "integer", min = 1,
      invalid = "error"
    ),
    "eigencore.gram_svd_max_dimension_wide" = normalize_scalar_option(
      "eigencore.gram_svd_max_dimension_wide", 1024, "integer", min = 1,
      invalid = "error"
    ),
    "eigencore.gram_svd_memory_mb" = normalize_scalar_option(
      "eigencore.gram_svd_memory_mb", 64, "numeric", min = 0,
      invalid = "error"
    ),
    "eigencore.gram_svd_rank_fraction_limit" = normalize_scalar_option(
      "eigencore.gram_svd_rank_fraction_limit", 0.5, "numeric",
      min = .Machine$double.eps, max = 1, allow_infinite = FALSE,
      invalid = "error"
    ),
    "eigencore.gram_svd_min_aspect_ratio" = normalize_scalar_option(
      "eigencore.gram_svd_min_aspect_ratio", 2, "numeric", min = 1,
      allow_infinite = FALSE, invalid = "error"
    ),
    "eigencore.gram_svd_work_budget" = normalize_scalar_option(
      "eigencore.gram_svd_work_budget", Inf, "numeric", min = 0,
      invalid = "error"
    ),
    "eigencore.implicit_gram_svd_min_dimension" = normalize_scalar_option(
      "eigencore.implicit_gram_svd_min_dimension", 64L, "integer", min = 1
    ),
    "eigencore.promote_retained_golub_kahan" = isTRUE(
      getOption("eigencore.promote_retained_golub_kahan", FALSE)
    ),
    "eigencore.golub_kahan_projected_stop" = isTRUE(
      getOption("eigencore.golub_kahan_projected_stop", FALSE)
    ),
    "eigencore.golub_kahan_prefix_diagnostics" = isTRUE(
      getOption("eigencore.golub_kahan_prefix_diagnostics", FALSE)
    ),
    "eigencore.randomized_stage_timing" = isTRUE(
      getOption("eigencore.randomized_stage_timing", FALSE)
    ),
    "eigencore.randomized_adaptive_stop" = isTRUE(
      getOption("eigencore.randomized_adaptive_stop", TRUE)
    )
  )
  stopifnot(identical(names(values), planner_policy_keys()))
  structure(values, schema_version = 1L, class = "eigencore_planner_policy")
}

#' @keywords internal
planner_policy_value <- function(plan, name) {
  policy <- plan$planner_policy
  if (!inherits(policy, "eigencore_planner_policy") || is.null(policy[[name]])) {
    plan_error(
      "corrupt_plan", "planner_policy", name, names(policy),
      paste0("Plan is missing frozen policy key ", name, ".")
    )
  }
  policy[[name]]
}

# Execution policy is dynamically scoped without mutating base R options.
# Nested solves restore the prior snapshot on exit.
.eigencore_execution_policy <- new.env(parent = emptyenv())
.eigencore_execution_policy$current <- NULL

#' @keywords internal
with_execution_policy <- function(policy, code) {
  previous <- .eigencore_execution_policy$current
  .eigencore_execution_policy$current <- policy
  on.exit({
    .eigencore_execution_policy$current <- previous
  }, add = TRUE)
  force(code)
}

#' @keywords internal
execution_policy_option <- function(name, default = NULL) {
  policy <- .eigencore_execution_policy$current
  if (inherits(policy, "eigencore_planner_policy") && name %in% names(policy)) {
    return(policy[[name]])
  }
  getOption(name, default)
}

#' @keywords internal
new_plan_execution <- function(problem_type, tol = 1e-8, maxit = NULL,
                               vectors = TRUE, certify = TRUE,
                               allow_dense_fallback = "auto",
                               initial_subspace = NULL) {
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive numeric value.", call. = FALSE)
  }
  if (!is.null(maxit)) {
    maxit <- as.integer(maxit)
    if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
      stop("maxit must be NULL or one positive integer.", call. = FALSE)
    }
  }
  allow_dense_fallback <- match.arg(
    allow_dense_fallback,
    c("auto", "never", "always")
  )
  if (identical(problem_type, "eigen")) {
    if (!is.logical(vectors) || length(vectors) != 1L || is.na(vectors)) {
      stop("vectors must be TRUE or FALSE for an eigen plan.", call. = FALSE)
    }
  } else {
    vectors <- match.arg(vectors, c("both", "left", "right", "none"))
  }
  if (!is.logical(certify) || length(certify) != 1L || is.na(certify)) {
    stop("certify must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(initial_subspace)) {
    initial_subspace <- as.matrix(initial_subspace)
    storage.mode(initial_subspace) <- "double"
  }
  list(
    tol = as.numeric(tol),
    maxit = maxit,
    vectors = vectors,
    certify = certify,
    allow_dense_fallback = allow_dense_fallback,
    dense_fallback_budget_bytes = as.numeric(NA),
    initial_subspace = initial_subspace
  )
}

#' @keywords internal
new_memory_record <- function(components, native_bytes = NA_real_) {
  components <- lapply(components, function(x) as.numeric(utils::object.size(x)))
  r_bytes <- sum(unlist(components, use.names = FALSE))
  native_known <- length(native_bytes) == 1L && is.finite(native_bytes) && native_bytes >= 0
  structure(
    list(
      schema_version = 1L,
      r_bytes = r_bytes,
      native_bytes = if (native_known) as.numeric(native_bytes) else NA_real_,
      total_bytes = r_bytes + if (native_known) as.numeric(native_bytes) else 0,
      complete = native_known,
      components = components
    ),
    class = "eigencore_memory"
  )
}

#' @keywords internal
new_plan_serialization <- function(identities, planned_method,
                                   method_descriptor = NULL,
                                   controls = list(), execution = list(),
                                   planner_policy = list()) {
  portable <- all(vapply(identities, function(x) isTRUE(x$portable), logical(1L)))
  structure(list(
    schema_version = 1L,
    portable = portable,
    originating_session = eigencore_session_id(),
    incompatibility_reason = if (portable) NULL else "opaque callback identity is session-local",
    dispatch_token = stable_raw_hash(list(1L, planned_method)),
    method_descriptor_token = stable_raw_hash(method_descriptor),
    controls_token = stable_raw_hash(controls),
    execution_token = stable_raw_hash(execution),
    planner_policy_token = stable_raw_hash(planner_policy),
    operator_identity_token = stable_raw_hash(identities)
  ), class = "eigencore_serialization")
}

#' @keywords internal
plan_error <- function(code, field = NULL, expected = NULL, actual = NULL,
                       message = NULL) {
  if (is.null(message)) {
    message <- paste0("Invalid eigencore plan (", code, ")", if (!is.null(field)) {
      paste0(": ", field)
    } else {
      ""
    }, ".")
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
    class = c("eigencore_plan_error", "error", "condition")
  ))
}

#' @keywords internal
required_plan_fields <- function() {
  c(
    "schema_version", "problem", "problem_type", "requested", "target",
    "method_descriptor", "planned_method", "method", "reasons", "fallback",
    "controls", "execution", "planner_policy", "operator_identity",
    "serialization", "memory"
  )
}

#' @keywords internal
validate_eigencore_plan <- function(plan) {
  if (!inherits(plan, "eigencore_plan")) {
    plan_error("corrupt_plan", "class", "eigencore_plan", class(plan))
  }
  if (is.null(plan$schema_version) || !identical(plan$schema_version, 1L)) {
    plan_error(
      "unsupported_schema", "schema_version", 1L, plan$schema_version,
      "This plan is not an executable eigencore 1.2 schema. Replan from the original problem."
    )
  }
  missing <- setdiff(required_plan_fields(), names(plan))
  if (length(missing)) {
    plan_error("missing_field", missing[[1L]], "present", NULL)
  }
  if (!plan$problem_type %in% c("eigen", "svd") ||
      !identical(plan$problem$type, plan$problem_type)) {
    plan_error("corrupt_plan", "problem_type", plan$problem$type, plan$problem_type)
  }
  if (!is.numeric(plan$requested) || length(plan$requested) != 1L ||
      is.na(plan$requested) || plan$requested < 1) {
    plan_error("corrupt_plan", "requested", "one positive integer", plan$requested)
  }
  if (!is.character(plan$planned_method) || length(plan$planned_method) != 1L ||
      is.na(plan$planned_method) || !nzchar(plan$planned_method) ||
      !identical(plan$method, plan$planned_method)) {
    plan_error("corrupt_plan", "planned_method", plan$method, plan$planned_method)
  }
  expected_dispatch <- stable_raw_hash(list(1L, plan$planned_method))
  if (!identical(plan$serialization$dispatch_token, expected_dispatch)) {
    plan_error(
      "dispatch_unavailable", "planned_method", "registered frozen route",
      plan$planned_method
    )
  }
  if (!isTRUE(plan_dispatch_available(plan))) {
    plan_error(
      "dispatch_unavailable", "planned_method", "registered 1.2 dispatch",
      plan$planned_method
    )
  }
  token_fields <- list(
    method_descriptor_token = plan$method_descriptor,
    controls_token = plan$controls,
    execution_token = plan$execution,
    planner_policy_token = plan$planner_policy,
    operator_identity_token = plan$operator_identity
  )
  for (token_name in names(token_fields)) {
    actual_token <- stable_raw_hash(token_fields[[token_name]])
    if (!identical(plan$serialization[[token_name]], actual_token)) {
      plan_error(
        "corrupt_plan", sub("_token$", "", token_name),
        plan$serialization[[token_name]], actual_token
      )
    }
  }
  if (!inherits(plan$method_descriptor, "eigencore_method")) {
    plan_error("corrupt_plan", "method_descriptor", "eigencore_method",
               class(plan$method_descriptor))
  }
  expected_target <- target_label(plan$problem$target)
  if (!identical(plan$target, expected_target)) {
    plan_error("corrupt_plan", "target", expected_target, plan$target)
  }
  if (!inherits(plan$planner_policy, "eigencore_planner_policy") ||
      !identical(attr(plan$planner_policy, "schema_version"), 1L)) {
    plan_error("unsupported_schema", "planner_policy", 1L,
               attr(plan$planner_policy, "schema_version"))
  }
  if (!identical(names(plan$planner_policy), planner_policy_keys())) {
    plan_error(
      "corrupt_plan", "planner_policy", planner_policy_keys(),
      names(plan$planner_policy)
    )
  }
  current_identity <- operator_identity(plan$problem)
  if (!identical(current_identity, plan$operator_identity)) {
    plan_error(
      "operator_incompatible", "operator_identity", plan$operator_identity,
      current_identity,
      "The embedded operator no longer matches the identity frozen in this plan. Replan from the current problem."
    )
  }
  nonportable <- Filter(function(x) !isTRUE(x$portable), current_identity)
  if (length(nonportable) && any(vapply(
    nonportable,
    function(x) !identical(x$session_id, eigencore_session_id()),
    logical(1L)
  ))) {
    plan_error(
      "session_incompatible", "operator_identity", eigencore_session_id(),
      vapply(nonportable, `[[`, character(1L), "session_id"),
      "This plan contains an opaque callback operator from another R session. Replan with the restored operator or supply explicit portable identity provenance."
    )
  }
  execution <- plan$execution
  if (!is.list(execution) || !is.numeric(execution$tol) ||
      length(execution$tol) != 1L || !is.finite(execution$tol) ||
      execution$tol <= 0) {
    plan_error("corrupt_plan", "execution$tol", "one finite positive value",
               execution$tol)
  }
  if (!execution$allow_dense_fallback %in% c("auto", "never", "always")) {
    plan_error("corrupt_plan", "execution$allow_dense_fallback",
               c("auto", "never", "always"), execution$allow_dense_fallback)
  }
  invisible(plan)
}

#' @keywords internal
new_fallback_reason <- function(code, message, planned_method, attempted_method,
                                details = list()) {
  structure(list(
    schema_version = 1L,
    code = code,
    message = message,
    planned_method = planned_method,
    attempted_method = attempted_method,
    details = details
  ), class = "eigencore_fallback_reason")
}

#' @keywords internal
new_state_transition <- function(relation = "cold_start", reason = "no restart state supplied") {
  structure(list(
    schema_version = 1L,
    relation = relation,
    basis_used = FALSE,
    method_state_used = FALSE,
    invalidated = character(),
    reason = reason
  ), class = "eigencore_state_transition")
}

#' @keywords internal
runtime_actual_method <- function(plan, method_label, iter) {
  planned <- plan$planned_method %||% plan$method %||% method_label
  if (!identical(method_label, planned)) {
    return(method_label)
  }
  restart <- iter$restart %||% list()
  if (isTRUE(restart$fallback_used)) {
    route <- restart$fallback_method %||% restart$kind %||% "certification fallback"
    if (length(route) != 1L || is.na(route) || !nzchar(route)) {
      route <- restart$kind %||% "certification fallback"
    }
    return(paste0(route, " fallback from ", planned))
  }
  if (isTRUE(restart$refinement_passed)) {
    route <- restart$refinement_kind %||% "certified refinement"
    if (length(route) != 1L || is.na(route) || !nzchar(route)) {
      route <- "certified refinement"
    }
    return(paste0(route, " refinement from ", planned))
  }
  method_label
}

#' @keywords internal
result_workflow_fields <- function(plan, actual_method, iter, result = NULL) {
  planned_method <- plan$planned_method %||% plan$method %||% actual_method
  fallback_used <- !identical(actual_method, planned_method)
  fallback_reason <- if (fallback_used) {
    new_fallback_reason(
      "certificate_failed",
      paste0("The planned route did not certify; returned result uses ", actual_method, "."),
      planned_method,
      actual_method
    )
  } else {
    NULL
  }
  list(
    planned_method = planned_method,
    actual_method = actual_method,
    fallback_used = fallback_used,
    fallback_reason = fallback_reason,
    work = typed_work_record(iter, result, plan),
    memory = new_memory_record(list(metadata = list(plan = planned_method))),
    state_transition = new_state_transition(),
    restart_state = NULL
  )
}

#' @keywords internal
finalize_workflow_result <- function(result, plan) {
  actual_method <- runtime_actual_method(plan, result$method, result)
  result$method <- actual_method
  fields <- result_workflow_fields(plan, actual_method, result, result)
  for (name in names(fields)) {
    result[[name]] <- fields[[name]]
  }
  result
}
