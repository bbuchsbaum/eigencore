#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

raw_args <- commandArgs(trailingOnly = TRUE)
args <- benchmark_args(raw_args)

gram_max_arg <- benchmark_arg_value(raw_args, "--gram-max=")
gram_max <- if (is.null(gram_max_arg)) {
  1024L
} else {
  as.integer(gram_max_arg)
}
if (is.na(gram_max) || gram_max < 1L) {
  stop("--gram-max must be a positive integer.", call. = FALSE)
}

iterations <- if (is.na(args$iterations)) {
  if (args$quick) 1L else 3L
} else {
  args$iterations
}

methods <- args$methods %||% c(
  "eigencore",
  "RSpectra",
  "PRIMME",
  "irlba",
  if (args$quick) "base"
)
gate_subject <- args$subject %||% "eigencore"
if (!gate_subject %in% methods) {
  stop(
    "Gram cutoff gate subject `", gate_subject, "` is not in the selected methods.",
    call. = FALSE
  )
}

old_options <- options(
  eigencore.gram_svd_max_dimension = gram_max,
  eigencore.gram_svd_max_dimension_wide = gram_max
)
on.exit(options(old_options), add = TRUE)

gram_cutoff_case <- function(case, m, n, density, rank, seed,
                             base_meaningful = FALSE) {
  A <- tall_skinny_sparse(m, n, density = density, seed = seed)
  list(
    case = case,
    id = paste0(case, ":", m, "x", n),
    A = A,
    rank = rank,
    small_side = min(m, n),
    density = density,
    base_meaningful = base_meaningful
  )
}

wide_gram_cutoff_case <- function(case, m, n, density, rank, seed,
                                  base_meaningful = FALSE) {
  tall <- gram_cutoff_case(
    case = case,
    m = n,
    n = m,
    density = density,
    rank = rank,
    seed = seed,
    base_meaningful = base_meaningful
  )
  tall$id <- paste0(case, ":", m, "x", n)
  tall$A <- Matrix::t(tall$A)
  tall
}

cases <- if (args$quick) {
  list(
    gram_cutoff_case("tall_sparse_600", 1200L, 600L, 0.004, 4L, 9201L, TRUE),
    wide_gram_cutoff_case("wide_sparse_600", 600L, 1200L, 0.004, 4L, 9202L, TRUE),
    gram_cutoff_case("tall_sparse_768", 1536L, 768L, 0.003, 4L, 9203L),
    wide_gram_cutoff_case("wide_sparse_768", 768L, 1536L, 0.003, 4L, 9204L),
    gram_cutoff_case("tall_sparse_1024", 2048L, 1024L, 0.002, 4L, 9205L),
    wide_gram_cutoff_case("wide_sparse_1024", 1024L, 2048L, 0.002, 4L, 9206L)
  )
} else {
  list(
    gram_cutoff_case("tall_sparse_600", 50000L, 600L, 0.0010, 10L, 9301L),
    wide_gram_cutoff_case("wide_sparse_600", 600L, 50000L, 0.0010, 10L, 9302L),
    gram_cutoff_case("tall_sparse_768", 60000L, 768L, 0.0008, 12L, 9303L),
    wide_gram_cutoff_case("wide_sparse_768", 768L, 60000L, 0.0008, 12L, 9304L),
    gram_cutoff_case("tall_sparse_1024", 80000L, 1024L, 0.0006, 12L, 9305L),
    wide_gram_cutoff_case("wide_sparse_1024", 1024L, 80000L, 0.0006, 12L, 9306L)
  )
}
cases <- filter_benchmark_cases(cases, args$cases)

quick_reference_contract_gate <- function(gate, quick = FALSE) {
  if (!isTRUE(quick) || !nrow(gate)) {
    return(gate)
  }
  gate$speed_gate <- TRUE
  gate$memory_gate <- TRUE
  gate$passed <- gate$subject_certified
  note <- "quick smoke contract-only; speed and memory ratios are diagnostics"
  gate$note <- ifelse(nzchar(gate$note), paste(gate$note, note, sep = "; "), note)
  gate
}

gram_cutoff_contract <- function(case_rows, case, subject, gram_max) {
  subject_row <- case_rows[case_rows$method == subject, , drop = FALSE]
  if (nrow(subject_row) != 1L) {
    stop("Gram cutoff contract requires exactly one subject row.", call. = FALSE)
  }
  gram_path <- identical(
    subject_row$solver_label,
    "native certified Gram SVD special case"
  )
  certified <- isTRUE(subject_row$certificate_passed) &&
    isTRUE(subject_row$nconv >= case$rank)
  original_coordinate_certificate <- isTRUE(
    subject_row$certified_in_original_coordinates
  )
  original_residuals_finite <- all(is.finite(c(
    subject_row$max_left_residual,
    subject_row$max_right_residual,
    subject_row$max_cyclic_residual,
    subject_row$max_backward_error
  )))
  smaller_gram_only <- isTRUE(subject_row$materialized_gram) &&
    !isTRUE(subject_row$normal_operator_implicit) &&
    isTRUE(subject_row$gram_dimension == case$small_side)
  no_original_sparse_densification <- inherits(case$A, "sparseMatrix") &&
    gram_path &&
    smaller_gram_only
  within_cutoff <- !is.na(subject_row$gram_dimension) &&
    subject_row$gram_dimension <= gram_max
  passed <- certified &&
    gram_path &&
    original_coordinate_certificate &&
    original_residuals_finite &&
    smaller_gram_only &&
    no_original_sparse_densification &&
    within_cutoff &&
    !isTRUE(subject_row$fallback_used)

  data.frame(
    case = case$case,
    m = nrow(case$A),
    n = ncol(case$A),
    rank = case$rank,
    small_side = case$small_side,
    subject = subject,
    solver_label = subject_row$solver_label,
    gram_max_dimension = gram_max,
    gram_dimension = subject_row$gram_dimension,
    gram_path = gram_path,
    within_cutoff = within_cutoff,
    materialized_gram = subject_row$materialized_gram,
    normal_operator_implicit = subject_row$normal_operator_implicit,
    native_gram_kernel = subject_row$native_gram_kernel,
    native_gram_eigensolver = subject_row$native_gram_eigensolver,
    gram_certificate_passed = subject_row$gram_certificate_passed,
    certified = certified,
    certified_in_original_coordinates = original_coordinate_certificate,
    original_residuals_finite = original_residuals_finite,
    no_original_sparse_densification = no_original_sparse_densification,
    fallback_used = subject_row$fallback_used,
    max_backward_error = subject_row$max_backward_error,
    passed = passed,
    stringsAsFactors = FALSE
  )
}

rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  message_benchmark_case("bench-svd-gram-cutoff", case)
  active_methods <- methods
  if (!isTRUE(args$include_dense) && !isTRUE(case$base_meaningful)) {
    active_methods <- setdiff(active_methods, "base")
  }
  out <- benchmark_svd_case(
    case$A,
    rank = case$rank,
    methods = active_methods,
    iterations = iterations,
    tol = 1e-8,
    seed = 9400L + i
  )
  out$case <- case$case
  out$m <- nrow(case$A)
  out$n <- ncol(case$A)
  out$rank <- case$rank
  out$small_side <- case$small_side
  out$gram_max_dimension <- gram_max
  out$input_density <- case$density
  out
})

result <- do.call(rbind, rows)
row.names(result) <- NULL
print(result)

external_references <- setdiff(unique(result$method), svd_internal_methods())
can_evaluate_gates <- gate_subject %in% result$method &&
  length(setdiff(external_references, gate_subject)) > 0L

gates <- if (isTRUE(can_evaluate_gates)) {
  lapply(seq_along(cases), function(i) {
    case <- cases[[i]]
    case_rows <- result[result$case == case$case, , drop = FALSE]
    gate_rows <- case_rows[
      case_rows$method == gate_subject |
        !case_rows$method %in% svd_internal_methods(),
      ,
      drop = FALSE
    ]
    gate <- evaluate_reference_gate(
      gate_rows,
      subject = gate_subject,
      references = setdiff(unique(gate_rows$method), gate_subject),
      requested = case$rank,
      speed_ratio_required = release_speed_gate("svd"),
      memory_ratio_required = release_memory_gate("svd")
    )
    gate <- quick_reference_contract_gate(gate, quick = args$quick)
    gate$case <- case$case
    gate$m <- nrow(case$A)
    gate$n <- ncol(case$A)
    gate$rank <- case$rank
    gate$small_side <- case$small_side
    gate$gram_max_dimension <- gram_max
    gate
  })
} else {
  list()
}
gates <- if (length(gates)) do.call(rbind, gates) else data.frame()
row.names(gates) <- NULL
print(gates)

memory_diagnostics <- if (isTRUE(can_evaluate_gates)) {
  lapply(seq_along(cases), function(i) {
    case <- cases[[i]]
    case_rows <- result[result$case == case$case, , drop = FALSE]
    gate_rows <- case_rows[
      case_rows$method == gate_subject |
        !case_rows$method %in% svd_internal_methods(),
      ,
      drop = FALSE
    ]
    diagnostics <- evaluate_memory_diagnostics(
      gate_rows,
      subject = gate_subject,
      references = setdiff(unique(gate_rows$method), gate_subject),
      requested = case$rank
    )
    diagnostics$case <- case$case
    diagnostics$m <- nrow(case$A)
    diagnostics$n <- ncol(case$A)
    diagnostics$rank <- case$rank
    diagnostics$small_side <- case$small_side
    diagnostics$gram_max_dimension <- gram_max
    diagnostics
  })
} else {
  list()
}
memory_diagnostics <- if (length(memory_diagnostics)) {
  do.call(rbind, memory_diagnostics)
} else {
  data.frame()
}
row.names(memory_diagnostics) <- NULL
print(memory_diagnostics)

contracts <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  case_rows <- result[result$case == case$case, , drop = FALSE]
  gram_cutoff_contract(case_rows, case, gate_subject, gram_max)
})
contracts <- do.call(rbind, contracts)
row.names(contracts) <- NULL
print(contracts)

policy <- data.frame(
  subject = gate_subject,
  gram_max_dimension = gram_max,
  cases = paste(contracts$case, collapse = ","),
  contracts_passed = all(contracts$passed),
  reference_gates_passed = if (nrow(gates)) all(gates$passed) else FALSE,
  recommendation = if (all(contracts$passed) && nrow(gates) && all(gates$passed)) {
    "raise_cutoff_candidate"
  } else {
    "keep_current_cutoff_or_collect_more_evidence"
  },
  stringsAsFactors = FALSE
)
print(policy)

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "svd-gram-cutoff-rows"))
  message("saved gates: ", save_benchmark_result(gates, "svd-gram-cutoff-gates"))
  message(
    "saved memory diagnostics: ",
    save_benchmark_result(memory_diagnostics, "svd-gram-cutoff-memory")
  )
  message(
    "saved contracts: ",
    save_benchmark_result(contracts, "svd-gram-cutoff-contracts")
  )
  message("saved policy: ", save_benchmark_result(policy, "svd-gram-cutoff-policy"))
}

if (args$strict) {
  if (!nrow(gates)) {
    stop("Gram cutoff strict mode requires at least one external reference row.", call. = FALSE)
  }
  if (!all(contracts$passed)) {
    stop("Gram cutoff strict mode failed certification/provenance contract.", call. = FALSE)
  }
}
