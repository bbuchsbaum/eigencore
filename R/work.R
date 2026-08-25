# Typed logical work accounting for eigencore 1.2 workflows.

.eigencore_work_context <- new.env(parent = emptyenv())
.eigencore_work_context$current <- NULL

#' @keywords internal
work_counter_fields <- function() {
  c(
    "operator_block_calls", "operator_columns",
    "adjoint_block_calls", "adjoint_columns",
    "metric_block_calls", "metric_columns",
    "preconditioner_calls", "preconditioner_columns",
    "certification_operator_block_calls",
    "certification_operator_columns",
    "certification_adjoint_block_calls",
    "certification_adjoint_columns"
  )
}

#' @keywords internal
work_record_fields <- function() {
  c(
    "schema_version", "complete", work_counter_fields(),
    "iterations", "restarts", "setup_seconds", "solve_seconds",
    "certification_seconds", "total_seconds", "legacy_matvecs"
  )
}

#' @keywords internal
new_work_context <- function(plan) {
  context <- new.env(parent = emptyenv())
  context$phase <- "solve"
  context$operator_depth <- 0L
  context$phase_seconds <- c(solve = 0, certification = 0)
  context$counters <- stats::setNames(
    as.list(integer(length(work_counter_fields()))),
    work_counter_fields()
  )
  context$identity_A <- plan$operator_identity$A %||% NULL
  context$identity_B <- plan$operator_identity$B %||% NULL
  context$operator_observations_authoritative <-
    !identical(context$identity_A$origin %||% NULL, "builtin")
  context$metric_observations_authoritative <-
    !is.null(context$identity_B) &&
    !identical(context$identity_B$origin %||% NULL, "builtin")
  context
}

#' @keywords internal
with_work_context <- function(context, code) {
  previous <- .eigencore_work_context$current
  .eigencore_work_context$current <- context
  on.exit({
    .eigencore_work_context$current <- previous
  }, add = TRUE)
  force(code)
}

#' @keywords internal
work_identity_equal <- function(left, right) {
  if (is.null(left) || is.null(right)) {
    return(FALSE)
  }
  identical(left$operator_id, right$operator_id) &&
    identical(left$revision, right$revision) &&
    identical(left$dim, right$dim) &&
    identical(left$dtype, right$dtype) &&
    identical(left$structure, right$structure)
}

#' @keywords internal
work_block_columns <- function(X) {
  dims <- dim(X)
  if (length(dims) == 2L) as.integer(dims[[2L]]) else 1L
}

#' @keywords internal
work_operator_enter <- function(identity, kind = c("operator", "adjoint"), X) {
  kind <- match.arg(kind)
  context <- .eigencore_work_context$current
  if (is.null(context)) {
    return(NULL)
  }
  outermost <- identical(context$operator_depth, 0L)
  context$operator_depth <- context$operator_depth + 1L
  if (outermost) {
    columns <- work_block_columns(X)
    metric <- work_identity_equal(identity, context$identity_B)
    phase <- context$phase
    if (metric) {
      block_field <- "metric_block_calls"
      column_field <- "metric_columns"
    } else if (identical(phase, "certification")) {
      block_field <- if (identical(kind, "adjoint")) {
        "certification_adjoint_block_calls"
      } else {
        "certification_operator_block_calls"
      }
      column_field <- if (identical(kind, "adjoint")) {
        "certification_adjoint_columns"
      } else {
        "certification_operator_columns"
      }
    } else {
      block_field <- paste0(kind, "_block_calls")
      column_field <- paste0(kind, "_columns")
    }
    context$counters[[block_field]] <- context$counters[[block_field]] + 1L
    context$counters[[column_field]] <-
      context$counters[[column_field]] + columns
  }
  list(context = context, outermost = outermost)
}

#' @keywords internal
work_operator_exit <- function(token) {
  if (is.null(token)) {
    return(invisible(NULL))
  }
  context <- token$context
  context$operator_depth <- max(0L, context$operator_depth - 1L)
  invisible(NULL)
}

#' @keywords internal
work_phase_enter <- function(phase = c("solve", "certification")) {
  phase <- match.arg(phase)
  context <- .eigencore_work_context$current
  if (is.null(context) || identical(context$phase, phase)) {
    return(NULL)
  }
  token <- list(
    context = context,
    previous = context$phase,
    phase = phase,
    started = proc.time()[["elapsed"]]
  )
  context$phase <- phase
  token
}

#' @keywords internal
work_phase_exit <- function(token) {
  if (is.null(token)) {
    return(invisible(NULL))
  }
  elapsed <- proc.time()[["elapsed"]] - token$started
  token$context$phase_seconds[[token$phase]] <-
    token$context$phase_seconds[[token$phase]] + elapsed
  token$context$phase <- token$previous
  invisible(NULL)
}

#' @keywords internal
work_record_preconditioner_call <- function(X) {
  context <- .eigencore_work_context$current
  if (is.null(context)) {
    return(invisible(NULL))
  }
  context$counters$preconditioner_calls <-
    context$counters$preconditioner_calls + 1L
  context$counters$preconditioner_columns <-
    context$counters$preconditioner_columns + work_block_columns(X)
  invisible(NULL)
}

#' @keywords internal
whole_count_or_na <- function(x) {
  if (is.null(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    return(NA_integer_)
  }
  as.integer(round(x))
}

#' @keywords internal
first_work_value <- function(...) {
  candidates <- list(...)
  for (value in candidates) {
    if (!is.null(value) && length(value) == 1L && !is.na(value)) {
      return(value)
    }
  }
  NULL
}

#' @keywords internal
new_typed_work_record <- function(values) {
  missing <- setdiff(work_record_fields(), names(values))
  for (field in missing) {
    values[[field]] <- NA
  }
  values$schema_version <- 1L
  count_fields <- c(work_counter_fields(), "iterations", "restarts",
                    "legacy_matvecs")
  for (field in count_fields) {
    values[[field]] <- whole_count_or_na(values[[field]])
  }
  timing_fields <- c(
    "setup_seconds", "solve_seconds", "certification_seconds", "total_seconds"
  )
  for (field in timing_fields) {
    value <- values[[field]]
    values[[field]] <- if (is.null(value) || length(value) != 1L ||
      is.na(value) || !is.finite(value) || value < 0) {
      NA_real_
    } else {
      as.numeric(value)
    }
  }
  required_counters <- work_counter_fields()
  values$complete <- all(!vapply(values[required_counters], is.na, logical(1L)))
  values <- values[work_record_fields()]
  structure(values, class = "eigencore_work")
}

#' @keywords internal
typed_work_record <- function(iter, result = NULL, plan = NULL) {
  result <- result %||% list()
  plan <- plan %||% result$plan %||% NULL
  problem_type <- plan$problem_type %||% NULL
  if (is.null(problem_type)) {
    problem_type <- if (inherits(result, "eigencore_eigen_result")) {
      "eigen"
    } else if (inherits(result, "eigencore_svd_result")) {
      "svd"
    } else {
      NULL
    }
  }
  execution <- plan$execution %||% list()
  certificate <- result$certificate %||% iter$certificate %||% NULL
  vectors_available <- if (identical(problem_type, "eigen")) {
    !is.null(result$vectors %||% iter$vectors)
  } else {
    !is.null(result$u %||% iter$u) && !is.null(result$v %||% iter$v)
  }
  matvecs <- whole_count_or_na(first_work_value(iter$matvecs, result$matvecs))
  iterations <- whole_count_or_na(first_work_value(
    iter$iterations, result$iterations, 1L
  ))
  restarts <- whole_count_or_na(first_work_value(
    iter$restarts, iter$restart$restarts_used,
    result$restarts, result$restart$restarts_used, 0L
  ))
  if (!inherits(plan, "eigencore_plan") || !identical(plan$schema_version, 1L)) {
    return(new_typed_work_record(list(
      schema_version = 1L,
      complete = FALSE,
      preconditioner_calls = whole_count_or_na(first_work_value(
        iter$preconditioner_calls, result$preconditioner_calls
      )),
      iterations = iterations,
      restarts = restarts,
      legacy_matvecs = matvecs
    )))
  }
  cert_columns <- whole_count_or_na(first_work_value(
    iter$certification_operator_columns,
    result$certification_operator_columns
  ))
  cert_block_calls <- whole_count_or_na(first_work_value(
    iter$certification_operator_block_calls,
    result$certification_operator_block_calls
  ))
  if (is.na(cert_columns) && isTRUE(execution$certify) && vectors_available) {
    cert_columns <- if (identical(problem_type, "eigen")) {
      ncol(as.matrix(result$vectors %||% iter$vectors))
    } else {
      ncol(as.matrix(result$v %||% iter$v))
    }
  }
  if (is.na(cert_block_calls) && !is.na(cert_columns) && cert_columns > 0L) {
    certificate_width <- if (identical(problem_type, "eigen")) {
      ncol(as.matrix(result$vectors %||% iter$vectors))
    } else {
      ncol(as.matrix(result$v %||% iter$v))
    }
    if (is.null(certificate_width) || !length(certificate_width) ||
        is.na(certificate_width) || certificate_width < 1L) {
      certificate_width <- as.integer(plan$requested %||% cert_columns)
    }
    cert_block_calls <- as.integer(ceiling(cert_columns / certificate_width))
  }
  if (identical(cert_columns, 0L)) {
    cert_block_calls <- 0L
  }
  if (!isTRUE(execution$certify) || !vectors_available || is.null(certificate)) {
    cert_columns <- 0L
    cert_block_calls <- 0L
  }

  operator_block_calls <- adjoint_block_calls <- NA_integer_
  operator_columns <- adjoint_columns <- NA_integer_
  certification_adjoint_block_calls <- whole_count_or_na(first_work_value(
    iter$certification_adjoint_block_calls,
    result$certification_adjoint_block_calls,
    0L
  ))
  certification_adjoint_columns <- whole_count_or_na(first_work_value(
    iter$certification_adjoint_columns,
    result$certification_adjoint_columns,
    0L
  ))
  metric_block_calls <- metric_columns <- if (is.null(plan$operator_identity$B)) {
    0L
  } else {
    NA_integer_
  }

  if (identical(problem_type, "eigen")) {
    total_blocks <- whole_count_or_na(first_work_value(
      iter$operator_block_calls, result$operator_block_calls
    ))
    total_columns <- whole_count_or_na(first_work_value(
      iter$operator_columns, result$operator_columns
    ))
    operator_block_calls <- if (!is.na(total_blocks) && !is.na(cert_block_calls)) {
      max(0L, total_blocks - cert_block_calls)
    } else {
      matvecs
    }
    operator_columns <- if (!is.na(total_columns) && !is.na(cert_columns)) {
      max(0L, total_columns - cert_columns)
    } else if (!is.na(matvecs)) {
      method <- plan$planned_method %||% result$planned_method %||% result$method %||% ""
      if (grepl("LOBPCG", method, fixed = TRUE)) {
        matvecs * as.integer(plan$requested %||% 1L)
      } else {
        matvecs
      }
    } else {
      NA_integer_
    }
    adjoint_block_calls <- whole_count_or_na(first_work_value(
      iter$adjoint_block_calls, result$adjoint_block_calls, 0L
    ))
    adjoint_columns <- whole_count_or_na(first_work_value(
      iter$adjoint_columns, result$adjoint_columns, 0L
    ))
  } else if (identical(problem_type, "svd")) {
    adjoint_matvecs <- whole_count_or_na(first_work_value(
      iter$adjoint_matvecs, result$adjoint_matvecs,
      iter$restart$adjoint_matvecs, result$restart$adjoint_matvecs
    ))
    method <- plan$planned_method %||% result$planned_method %||% result$method %||% ""
    if (!is.na(adjoint_matvecs) && !is.na(matvecs)) {
      adjoint_block_calls <- adjoint_matvecs
      operator_block_calls <- max(0L, matvecs - adjoint_matvecs)
    } else if (!is.na(matvecs) && grepl("Golub-Kahan", method, fixed = TRUE)) {
      adjoint_block_calls <- as.integer(floor(matvecs / 2))
      operator_block_calls <- matvecs - adjoint_block_calls
    } else if (identical(matvecs, 0L)) {
      operator_block_calls <- 0L
      adjoint_block_calls <- 0L
    }
    if (!is.na(operator_block_calls)) {
      operator_columns <- operator_block_calls
    }
    if (!is.na(adjoint_block_calls)) {
      adjoint_columns <- adjoint_block_calls
    }
    if (isTRUE(execution$certify) && vectors_available) {
      certification_adjoint_block_calls <- 1L
      certification_adjoint_columns <- ncol(as.matrix(result$u %||% iter$u))
    }
  }

  preconditioner_calls <- whole_count_or_na(first_work_value(
    iter$preconditioner_calls, result$preconditioner_calls,
    iter$restart$preconditioner_calls, result$restart$preconditioner_calls
  ))
  has_preconditioner <- !is.null(plan$method_descriptor$preconditioner)
  if (is.na(preconditioner_calls) && !has_preconditioner) {
    preconditioner_calls <- 0L
  }
  preconditioner_columns <- if (!is.na(preconditioner_calls)) {
    if (preconditioner_calls == 0L) 0L else {
      preconditioner_calls * as.integer(plan$requested %||% 1L)
    }
  } else {
    NA_integer_
  }

  stage_seconds <- iter$stage_seconds %||% result$stage_seconds %||%
    iter$restart$stage_seconds %||% result$restart$stage_seconds %||% numeric()
  cert_seconds <- if ("certificate" %in% names(stage_seconds)) {
    as.numeric(stage_seconds[["certificate"]])
  } else {
    NA_real_
  }

  new_typed_work_record(list(
    schema_version = 1L,
    complete = FALSE,
    operator_block_calls = operator_block_calls,
    operator_columns = operator_columns,
    adjoint_block_calls = adjoint_block_calls,
    adjoint_columns = adjoint_columns,
    metric_block_calls = metric_block_calls,
    metric_columns = metric_columns,
    preconditioner_calls = preconditioner_calls,
    preconditioner_columns = preconditioner_columns,
    certification_operator_block_calls = cert_block_calls,
    certification_operator_columns = cert_columns,
    certification_adjoint_block_calls = certification_adjoint_block_calls,
    certification_adjoint_columns = certification_adjoint_columns,
    iterations = iterations,
    restarts = restarts,
    setup_seconds = NA_real_,
    solve_seconds = NA_real_,
    certification_seconds = cert_seconds,
    total_seconds = NA_real_,
    legacy_matvecs = matvecs
  ))
}

#' @keywords internal
finalize_work_record <- function(result, context, setup_seconds,
                                 execution_seconds, total_seconds) {
  current <- result$work
  if (!inherits(current, "eigencore_work")) {
    current <- typed_work_record(result, result, result$plan)
  }
  values <- unclass(current)
  counters <- context$counters
  for (kind in c("operator", "adjoint")) {
    block_field <- paste0(kind, "_block_calls")
    column_field <- paste0(kind, "_columns")
    cert_block_field <- paste0("certification_", block_field)
    cert_column_field <- paste0("certification_", column_field)
    observed_blocks <- counters[[block_field]]
    observed_columns <- counters[[column_field]]
    observed_cert_blocks <- counters[[cert_block_field]]
    observed_cert_columns <- counters[[cert_column_field]]
    derived_cert_blocks <- values[[cert_block_field]] %||% 0L
    derived_cert_columns <- values[[cert_column_field]] %||% 0L
    total_cert_blocks <- max(derived_cert_blocks, observed_cert_blocks)
    total_cert_columns <- max(derived_cert_columns, observed_cert_columns)
    values[[cert_block_field]] <- total_cert_blocks
    values[[cert_column_field]] <- total_cert_columns
    if (isTRUE(context$operator_observations_authoritative)) {
      # Native callback kernels may perform current-certificate applications
      # inside their C/C++ execution boundary, so the R callback wrapper sees
      # them while the dynamic phase is still "solve". Route diagnostics give
      # the certified subset; remove it once so typed solve + certificate work
      # equals the independently observed callback total.
      hidden_cert_blocks <- max(0L, total_cert_blocks - observed_cert_blocks)
      hidden_cert_columns <- max(0L, total_cert_columns - observed_cert_columns)
      values[[block_field]] <- max(0L, observed_blocks - hidden_cert_blocks)
      values[[column_field]] <- max(0L, observed_columns - hidden_cert_columns)
    }
  }
  if (isTRUE(context$metric_observations_authoritative)) {
    values$metric_block_calls <- counters$metric_block_calls
    values$metric_columns <- counters$metric_columns
  }
  for (field in c("preconditioner_calls", "preconditioner_columns")) {
    if (counters[[field]] > 0L) {
      values[[field]] <- counters[[field]]
    }
  }
  if (is.null(result$plan$operator_identity$B)) {
    values$metric_block_calls <- 0L
    values$metric_columns <- 0L
  }
  if (is.null(result$plan$method_descriptor$preconditioner) &&
      counters$preconditioner_calls == 0L) {
    values$preconditioner_calls <- 0L
    values$preconditioner_columns <- 0L
  }
  certification_seconds <- context$phase_seconds[["certification"]]
  if (!is.finite(certification_seconds) || certification_seconds <= 0) {
    certification_seconds <- values$certification_seconds
  }
  values$setup_seconds <- setup_seconds
  values$certification_seconds <- certification_seconds
  values$solve_seconds <- if (is.finite(certification_seconds)) {
    max(0, execution_seconds - certification_seconds)
  } else {
    NA_real_
  }
  values$total_seconds <- total_seconds
  new_typed_work_record(values)
}

#' Extract typed logical work diagnostics.
#'
#' Returns shared work counters and phase timings without redefining the
#' route-specific legacy `matvecs` field. Operator, adjoint, metric,
#' preconditioner, and certification work are reported separately. Unknown
#' legacy counters are `NA` and make `complete` false.
#'
#' @details `schema_version` identifies the frozen record layout and `complete`
#' is true only when every logical counter is known. Each `*_block_calls` field
#' counts logical applications, while its paired `*_columns` field counts the
#' vector columns processed by those applications. Certification has distinct
#' forward-operator and adjoint ledgers. `iterations` and `restarts` retain
#' route-level progress, the four `*_seconds` fields report measured phase or
#' total elapsed time when available, and `legacy_matvecs` preserves the result's
#' historical route-specific `matvecs` value without treating it as a common
#' unit across solver families.
#'
#' @param x An eigencore result or an existing `eigencore_work` record.
#' @param ... Reserved for future methods.
#' @return An `eigencore_work` schema-version-1 record.
#' @examples
#' fit <- eig_partial(diag(c(5, 3, 1)), k = 2)
#' work(fit)$operator_columns
#' work(fit)$legacy_matvecs
#' @export
work <- function(x, ...) {
  if (inherits(x, "eigencore_work")) {
    return(x)
  }
  if (inherits(x$work, "eigencore_work")) {
    return(x$work)
  }
  typed_work_record(x, x, x$plan %||% NULL)
}
