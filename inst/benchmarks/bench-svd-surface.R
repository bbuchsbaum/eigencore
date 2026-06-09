#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (is.na(args$iterations)) {
  if (args$quick) 1L else 5L
} else {
  args$iterations
}
methods <- svd_surface_default_methods(args)
gate_subject <- svd_surface_gate_subject(args, methods)

rank_deficient_sparse <- function(m, n, intrinsic_rank, density = 0.1, seed = 1L) {
  set.seed(seed)
  L <- Matrix::rsparsematrix(m, intrinsic_rank, density = density)
  R <- Matrix::rsparsematrix(intrinsic_rank, n, density = density)
  L %*% R
}

clustered_dense_svd <- function(m, n, rank, seed = 1L) {
  set.seed(seed)
  r <- min(m, n)
  U <- qr.Q(qr(matrix(stats::rnorm(m * r), nrow = m, ncol = r)))
  V <- qr.Q(qr(matrix(stats::rnorm(n * r), nrow = n, ncol = r)))
  d <- c(seq(10, 9.99, length.out = min(rank + 2L, r)),
         if (r > rank + 2L) seq(1, 0.1, length.out = r - rank - 2L) else numeric())
  U %*% (d * t(V))
}

complex_dense_svd <- function(m, n, seed = 1L) {
  set.seed(seed)
  real <- matrix(stats::rnorm(m * n), nrow = m, ncol = n)
  imag <- matrix(stats::rnorm(m * n), nrow = m, ncol = n)
  real + 1i * imag
}

smallest_sparse_svd <- function(m, n, seed = 1L) {
  set.seed(seed)
  values <- sort(stats::runif(n, min = 0.1, max = 10), decreasing = TRUE)
  Matrix::sparseMatrix(
    i = seq_len(n),
    j = seq_len(n),
    x = values,
    dims = c(m, n)
  )
}

interior_sparse_svd <- function(m, n, seed = 1L) {
  set.seed(seed)
  values <- sort(c(10, 5, 1, 0.2, 0.1, stats::runif(max(0L, n - 5L), 0.02, 8)),
                 decreasing = TRUE)
  Matrix::sparseMatrix(
    i = seq_len(n),
    j = seq_len(n),
    x = values,
    dims = c(m, n)
  )
}

projected_stop_comparison <- function(result) {
  rows <- lapply(split(result, result$case), function(case_rows) {
    plain <- case_rows[case_rows$method == "eigencore_golub_kahan", , drop = FALSE]
    projected <- case_rows[
      case_rows$method == "eigencore_golub_kahan_projected",
      ,
      drop = FALSE
    ]
    if (nrow(plain) != 1L || nrow(projected) != 1L) {
      return(NULL)
    }
    data.frame(
      case = unique(case_rows$case),
      rank = unique(case_rows$rank),
      plain_certified = plain$certificate_passed,
      projected_certified = projected$certificate_passed,
      plain_median = plain$median,
      projected_median = projected$median,
      projected_speed_ratio = plain$median / projected$median,
      plain_mem_alloc = plain$mem_alloc,
      projected_mem_alloc = projected$mem_alloc,
      projected_memory_ratio = plain$mem_alloc / projected$mem_alloc,
      plain_iterations = plain$final_iterations,
      projected_iterations = projected$final_iterations,
      iteration_savings = plain$final_iterations - projected$final_iterations,
      iteration_savings_fraction = benchmark_safe_ratio(
        plain$final_iterations - projected$final_iterations,
        plain$final_iterations
      ),
      plain_matvecs = plain$final_matvecs,
      projected_matvecs = projected$final_matvecs,
      matvec_savings = plain$final_matvecs - projected$final_matvecs,
      matvec_savings_fraction = benchmark_safe_ratio(
        plain$final_matvecs - projected$final_matvecs,
        plain$final_matvecs
      ),
      plain_reorthogonalization_seconds = plain$stage_reorthogonalization_seconds,
      projected_reorthogonalization_seconds = projected$stage_reorthogonalization_seconds,
      reorthogonalization_speed_ratio = benchmark_safe_ratio(
        plain$stage_reorthogonalization_seconds,
        projected$stage_reorthogonalization_seconds
      ),
      plain_reorthogonalization_passes = plain$reorthogonalization_passes,
      projected_reorthogonalization_passes = projected$reorthogonalization_passes,
      reorthogonalization_pass_savings =
        plain$reorthogonalization_passes - projected$reorthogonalization_passes,
      reorthogonalization_pass_savings_fraction = benchmark_safe_ratio(
        plain$reorthogonalization_passes - projected$reorthogonalization_passes,
        plain$reorthogonalization_passes
      ),
      projected_stop_requested = projected$projected_stop_requested,
      projected_stop_enabled = projected$projected_stop_enabled,
      projected_stop_disable_reason = projected$projected_stop_disable_reason,
      projected_stop = projected$projected_stop,
      projected_nconv = projected$projected_nconv,
      projected_max_residual = projected$projected_max_residual,
      plain_projected_checks = plain$projected_checks,
      projected_checks = projected$projected_checks,
      projected_seconds = projected$projected_seconds,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

svd_target_contract <- function(rows, quick = FALSE) {
  internal <- rows[
    rows$method %in% c("eigencore_smallest", "eigencore_interior"),
    ,
    drop = FALSE
  ]
  if (!nrow(internal)) {
    return(data.frame())
  }
  out <- lapply(seq_len(nrow(internal)), function(i) {
    row <- internal[i, , drop = FALSE]
    ref_method <- if (identical(row$method, "eigencore_smallest")) {
      "base_smallest"
    } else {
      "base_interior"
    }
    ref <- rows[
      rows$case == row$case & rows$method == ref_method,
      ,
      drop = FALSE
    ]
    speed_ratio <- if (nrow(ref) == 1L && isTRUE(ref$certificate_passed)) {
      ref$median / row$median
    } else {
      NA_real_
    }
    memory_ratio <- if (nrow(ref) == 1L && isTRUE(ref$certificate_passed)) {
      ref$mem_alloc / row$mem_alloc
    } else {
      NA_real_
    }
    certificate_gate <- isTRUE(row$certificate_passed) &&
      isTRUE(row$nconv >= row$rank) &&
      !isTRUE(row$scale_is_estimate)
    provenance_gate <- if (identical(row$method, "eigencore_smallest")) {
      identical(row$solver_label, "native certified Gram SVD special case") &&
        isTRUE(row$materialized_gram) &&
        identical(row$native_gram_eigensolver, "native_dense_symmetric_eigen")
    } else {
      identical(row$solver_label, eigencore:::native_interior_golub_kahan_label()) &&
        isTRUE(row$final_iterations == min(row$m, row$n))
    }
    performance_gate <- if (isTRUE(quick)) {
      TRUE
    } else {
      isTRUE(speed_ratio >= release_speed_gate("svd")) &&
        isTRUE(memory_ratio >= release_memory_gate("svd"))
    }
    data.frame(
      case = row$case,
      method = row$method,
      rank = row$rank,
      solver_label = row$solver_label,
      speed_ratio_vs_base_target = speed_ratio,
      memory_ratio_vs_base_target = memory_ratio,
      certificate_gate = certificate_gate,
      provenance_gate = provenance_gate,
      performance_gate = performance_gate,
      passed = certificate_gate && provenance_gate && performance_gate,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

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

cases <- if (args$quick) {
  list(
    list(case = "tall_sparse", id = "tall_sparse:600x90", A = tall_skinny_sparse(600L, 90L, density = 0.03, seed = 701), rank = 5L),
    list(case = "wide_sparse", id = "wide_sparse:90x600", A = Matrix::t(tall_skinny_sparse(600L, 90L, density = 0.03, seed = 702)), rank = 5L),
    list(case = "rank_deficient_sparse", id = "rank_deficient_sparse:160x50", A = rank_deficient_sparse(160L, 50L, 4L, seed = 703), rank = 6L),
    list(case = "smallest_sparse", id = "smallest_sparse:160x50", A = smallest_sparse_svd(160L, 50L, seed = 708), rank = 3L, methods = c("eigencore_smallest", "base_smallest"), gate = FALSE),
    list(case = "interior_sparse", id = "interior_sparse:160x50", A = interior_sparse_svd(160L, 50L, seed = 709), rank = 3L, methods = c("eigencore_interior", "base_interior"), gate = FALSE),
    list(case = "clustered_dense", id = "clustered_dense:120x60", A = clustered_dense_svd(120L, 60L, 6L, seed = 704), rank = 6L),
    list(case = "complex_dense", id = "complex_dense:80x40", A = complex_dense_svd(80L, 40L, seed = 707), rank = 5L, methods = c("eigencore", "base"), gate = FALSE),
    list(case = "slow_decay_dense", id = "slow_decay_dense:120x60", A = slow_decay_svd_matrix(120L, 60L, decay = 0.35, seed = 705), rank = 6L),
    list(case = "low_rank_sparse", id = "low_rank_sparse:300x80", A = rank_deficient_sparse(300L, 80L, 8L, density = 0.08, seed = 706), rank = 5L)
  )
} else {
  list(
    list(case = "tall_sparse", id = "tall_sparse:100000x500", A = tall_skinny_sparse(100000L, 500L, density = 0.002, seed = 701), rank = 20L),
    list(case = "wide_sparse", id = "wide_sparse:500x100000", A = Matrix::t(tall_skinny_sparse(100000L, 500L, density = 0.002, seed = 702)), rank = 20L),
    list(case = "rank_deficient_sparse", id = "rank_deficient_sparse:5000x500", A = rank_deficient_sparse(5000L, 500L, 20L, density = 0.01, seed = 703), rank = 30L),
    list(case = "smallest_sparse", id = "smallest_sparse:5000x500", A = smallest_sparse_svd(5000L, 500L, seed = 708), rank = 10L, methods = c("eigencore_smallest", "base_smallest"), gate = FALSE),
    list(case = "interior_sparse", id = "interior_sparse:5000x500", A = interior_sparse_svd(5000L, 500L, seed = 709), rank = 10L, methods = c("eigencore_interior", "base_interior"), gate = FALSE),
    list(case = "clustered_dense", id = "clustered_dense:2000x500", A = clustered_dense_svd(2000L, 500L, 20L, seed = 704), rank = 20L),
    list(case = "complex_dense", id = "complex_dense:240x120", A = complex_dense_svd(240L, 120L, seed = 707), rank = 10L, methods = c("eigencore", "base"), gate = FALSE),
    list(case = "slow_decay_dense", id = "slow_decay_dense:2000x500", A = slow_decay_svd_matrix(2000L, 500L, decay = 0.35, seed = 705), rank = 20L),
    list(case = "low_rank_sparse", id = "low_rank_sparse:10000x500", A = rank_deficient_sparse(10000L, 500L, 40L, density = 0.01, seed = 706), rank = 20L)
  )
}
cases <- filter_benchmark_cases(cases, args$cases)

rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  message_benchmark_case("bench-svd-surface", case)
  active_methods <- case$methods %||% if (args$quick && inherits(case$A, "matrix")) {
    methods
  } else {
    setdiff(methods, "base")
  }
  out <- benchmark_svd_case(
    case$A,
    rank = case$rank,
    methods = active_methods,
    iterations = iterations,
    tol = 1e-8,
    seed = 700L + i
  )
  out$case <- case$case
  out$m <- nrow(case$A)
  out$n <- ncol(case$A)
  out$rank <- case$rank
  out$gate <- isTRUE(case$gate %||% TRUE)
  out
})

result <- do.call(rbind, rows)
row.names(result) <- NULL
print(result)

gated_result <- result[result$gate %in% TRUE, , drop = FALSE]
can_evaluate_gates <- nrow(gated_result) > 0L &&
  gate_subject %in% gated_result$method &&
  any(!gated_result$method %in% svd_internal_methods())
gates <- if (isTRUE(can_evaluate_gates)) lapply(split(gated_result, gated_result$case), function(case_rows) {
  internal_methods <- svd_internal_methods()
  gate_rows <- case_rows[
    case_rows$method == gate_subject | !case_rows$method %in% internal_methods,
    ,
    drop = FALSE
  ]
  gate <- evaluate_reference_gate(
    gate_rows,
    subject = gate_subject,
    references = setdiff(unique(gate_rows$method), gate_subject),
    requested = unique(case_rows$rank)[[1L]],
    speed_ratio_required = release_speed_gate("svd"),
    memory_ratio_required = release_memory_gate("svd")
  )
  gate <- quick_reference_contract_gate(gate, quick = args$quick)
  gate$case <- unique(case_rows$case)
  gate$m <- unique(case_rows$m)
  gate$n <- unique(case_rows$n)
  gate$rank <- unique(case_rows$rank)
  gate
}) else {
  list()
}
gates <- if (length(gates)) do.call(rbind, gates) else data.frame()
row.names(gates) <- NULL
print(gates)

target_contract <- svd_target_contract(result, quick = args$quick)
cat("\nSVD target contracts\n")
print(target_contract)

memory_diagnostics <- if (isTRUE(can_evaluate_gates)) lapply(split(gated_result, gated_result$case), function(case_rows) {
  internal_methods <- svd_internal_methods()
  gate_rows <- case_rows[
    case_rows$method == gate_subject | !case_rows$method %in% internal_methods,
    ,
    drop = FALSE
  ]
  diagnostics <- evaluate_memory_diagnostics(
    gate_rows,
    subject = gate_subject,
    references = setdiff(unique(gate_rows$method), gate_subject),
    requested = unique(case_rows$rank)[[1L]]
  )
  diagnostics$case <- unique(case_rows$case)
  diagnostics$m <- unique(case_rows$m)
  diagnostics$n <- unique(case_rows$n)
  diagnostics$rank <- unique(case_rows$rank)
  diagnostics
}) else {
  list()
}
memory_diagnostics <- if (length(memory_diagnostics)) {
  do.call(rbind, memory_diagnostics)
} else {
  data.frame()
}
row.names(memory_diagnostics) <- NULL
print(memory_diagnostics)

projected_comparison <- NULL
has_projected_pair <- all(c(
  "eigencore_golub_kahan",
  "eigencore_golub_kahan_projected"
) %in% result$method)
if (isTRUE(args$svd_projected_stop) || isTRUE(args$h_candidate) || has_projected_pair) {
  projected_comparison <- projected_stop_comparison(result)
  print(projected_comparison)
}

if (args$save) {
  message("saved rows: ", save_benchmark_result(result, "svd-surface-rows"))
  message("saved gates: ", save_benchmark_result(gates, "svd-surface-gates"))
  message("saved target contracts: ", save_benchmark_result(target_contract, "svd-surface-target-contracts"))
  message("saved memory diagnostics: ", save_benchmark_result(memory_diagnostics, "svd-surface-memory"))
  if (!is.null(projected_comparison)) {
    message(
      "saved projected stop comparison: ",
      save_benchmark_result(projected_comparison, "svd-projected-stop-comparison")
    )
  }
}

if (args$strict) {
  if (!nrow(gates) && !nrow(target_contract)) {
    stop(
      "SVD surface strict mode requires the gate subject, a target contract, or at least one external reference.",
      call. = FALSE
    )
  }
  if (nrow(gates) && !all(gates$passed)) {
    stop("SVD surface benchmark failed PRD release gate.", call. = FALSE)
  }
  if (nrow(target_contract) && !all(target_contract$passed)) {
    stop("SVD target contract failed.", call. = FALSE)
  }
}
