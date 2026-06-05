#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (!is.na(args$iterations)) args$iterations else if (args$quick) 1L else 5L
tol <- 1e-8

failure_rows <- function(case, target, methods, requested, message) {
  data.frame(
    case = case,
    kind = "generalized_eigen",
    target = target,
    method = methods,
    median = NA_real_,
    min = NA_real_,
    mem_alloc = NA_real_,
    max_residual = NA_real_,
    max_backward_error = NA_real_,
    orthogonality_loss = NA_real_,
    certificate_passed = FALSE,
    certificate_type = NA_character_,
    norm_bound_type = NA_character_,
    scale_is_estimate = NA,
    nconv = 0L,
    requested = requested,
    iterations = NA_integer_,
    matvecs = NA_integer_,
    restarts = NA_integer_,
    ortho_passes = NA_integer_,
    locking_events = NA_integer_,
    block_size = NA_integer_,
    native = NA,
    native_kernels = NA,
    generalized = NA,
    orthogonalization_native = NA,
    orthogonalization_methods = NA_character_,
    q_rank_final = NA_integer_,
    constrained = NA,
    constraints_rank = NA_integer_,
    preconditioner_kind = NA_character_,
    preconditioner_native = NA,
    preconditioner_calls = NA_integer_,
    metric_solve_kind = NA_character_,
    metric_solve_label = NA_character_,
    metric_solves = NA_integer_,
    seed = NA_integer_,
    pkg_version = as.character(utils::packageVersion("eigencore")),
    solver_label = NA_character_,
    error = message,
    stringsAsFactors = FALSE
  )
}

with_case <- function(case, target, methods, requested, expr) {
  tryCatch({
    rows <- expr
    rows$case <- case
    rows$kind <- "generalized_eigen"
    rows$target <- target
    rows$requested <- requested
    rows$error <- ""
    rows
  }, error = function(e) {
    failure_rows(case, target, methods, requested, conditionMessage(e))
  })
}

generalized_lobpcg_tridiagonal_shift <- function(A, target) {
  kind <- if (inherits(target, "eigencore_target")) target$kind else "largest"
  if (!identical(kind, "largest")) {
    return(1e-3)
  }
  n <- nrow(A)
  if (!is.finite(n) || n < 1L) {
    return(1e-3)
  }
  max(1e-3, as.numeric(n))
}

benchmark_generalized_lobpcg_case <- function(A, B, k, target = smallest(),
                                              methods = c("eigencore", "base"),
                                              iterations = 3L, tol = 1e-8,
                                              seed = 1L, maxit = 200L,
                                              constraints = NULL,
                                              A_dense = NULL, B_dense = NULL) {
  methods <- intersect(
    methods,
    c(
      "eigencore_auto",
      "eigencore",
      "eigencore_lanczos_reference",
      "eigencore_shifted_diagonal",
      "eigencore_shifted_tridiagonal",
      "eigencore_constrained",
      "base"
    )
  )
  A_cert <- A_dense %||% as.matrix(A)
  B_cert <- B_dense %||% as.matrix(B)
  rows <- lapply(methods, function(method) {
    timed <- run_timed({
      if (identical(method, "eigencore_auto")) {
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          tol = tol,
          allow_dense_fallback = "auto"
        )
      } else if (identical(method, "eigencore")) {
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          method = lobpcg(maxit = maxit),
          tol = tol,
          allow_dense_fallback = "never"
        )
      } else if (identical(method, "eigencore_lanczos_reference")) {
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          method = lanczos(max_subspace = max(maxit, k)),
          tol = tol,
          allow_dense_fallback = "never"
        )
      } else if (identical(method, "eigencore_shifted_diagonal")) {
        preconditioner <- shifted_diagonal_preconditioner(Matrix::diag(A), shift = 1e-3)
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          method = lobpcg(maxit = maxit, preconditioner = preconditioner),
          tol = tol,
          allow_dense_fallback = "never"
        )
      } else if (identical(method, "eigencore_shifted_tridiagonal")) {
        preconditioner <- shifted_tridiagonal_preconditioner(
          A,
          shift = generalized_lobpcg_tridiagonal_shift(A, target)
        )
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          method = lobpcg(maxit = maxit, preconditioner = preconditioner),
          tol = tol,
          allow_dense_fallback = "never"
        )
      } else if (identical(method, "eigencore_constrained")) {
        if (is.null(constraints)) {
          constraints <- matrix(c(1, rep(0, nrow(A) - 1L)), ncol = 1L)
        }
        eig_partial(
          A,
          B = B,
          k = k,
          target = target,
          method = lobpcg(maxit = maxit, constraints = constraints),
          tol = tol,
          allow_dense_fallback = "never"
        )
      } else {
        eig <- eigencore:::dense_generalized_spd_eigen(A_cert, B_cert)
        idx <- eigencore:::order_indices(eig$values, target)
        idx <- idx[seq_len(k)]
        list(values = eig$values[idx], vectors = eig$vectors[, idx, drop = FALSE])
      }
    }, iterations = iterations, seed = seed)
    cert <- eigencore:::certify_eigen(
      A_cert,
      eigencore:::method_values(timed$value, kind = "eigen"),
      timed$value$vectors,
      B = B_cert,
      tol = tol
    )
    data.frame(
      method = method,
      median = timed$median,
      min = timed$min,
      mem_alloc = timed$mem_alloc,
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      orthogonality_loss = cert$max_orthogonality_loss,
      certificate_passed = cert$passed,
      certificate_type = cert$certificate_type,
      norm_bound_type = cert$norm_bound_type,
      scale_is_estimate = cert$scale_is_estimate,
      nconv = sum(cert$converged),
      iterations = result_iterations(timed$value),
      matvecs = result_matvecs(timed$value),
      restarts = result_restarts(timed$value),
      ortho_passes = result_ortho_passes(timed$value),
      locking_events = result_locking_events(timed$value),
      block_size = result_block_size(timed$value),
      native = result_restart_logical(timed$value, "native"),
      native_kernels = result_restart_logical(timed$value, "native_kernels"),
      generalized = result_restart_logical(timed$value, "generalized"),
      orthogonalization_native = result_restart_logical(timed$value, "orthogonalization_native"),
      orthogonalization_methods = result_restart_character(timed$value, "orthogonalization_methods"),
      q_rank_final = result_restart_integer(timed$value, "q_rank_final"),
      constrained = result_restart_logical(timed$value, "constrained"),
      constraints_rank = result_restart_integer(timed$value, "constraints_rank"),
      preconditioner_kind = result_preconditioner_field(timed$value, "kind"),
      preconditioner_native = result_preconditioner_field(timed$value, "native"),
      preconditioner_calls = result_preconditioner_calls(timed$value),
      metric_solve_kind = {
        restart <- timed$value$restart %||% list()
        metric_solve <- restart$metric_solve %||% NA_character_
        if (is.list(metric_solve)) {
          metric_solve$kind %||% NA_character_
        } else {
          NA_character_
        }
      },
      metric_solve_label = {
        restart <- timed$value$restart %||% list()
        metric_solve <- restart$metric_solve %||% NA_character_
        if (is.list(metric_solve)) {
          metric_solve$label %||% NA_character_
        } else {
          metric_solve
        }
      },
      metric_solves = result_restart_integer(timed$value, "metric_solves"),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      solver_label = timed$value$method %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

generalized_lobpcg_native_contract <- function(rows) {
  internal <- rows[rows$method %in% c(
    "eigencore",
    "eigencore_lanczos_reference",
    "eigencore_shifted_diagonal",
    "eigencore_shifted_tridiagonal",
    "eigencore_constrained"
  ), , drop = FALSE]
  if ("solver_label" %in% names(internal)) {
    native_label <- eigencore:::native_generalized_lobpcg_label()
    internal <- internal[
      internal$solver_label == native_label,
      ,
      drop = FALSE
    ]
  }
  if (!nrow(internal)) {
    return(data.frame())
  }
  out <- lapply(seq_len(nrow(internal)), function(i) {
    row <- internal[i, , drop = FALSE]
    certified <- isTRUE(row$certificate_passed) && row$nconv >= row$requested
    native_gate <- certified &&
      isTRUE(row$native) &&
      isTRUE(row$native_kernels) &&
      isTRUE(row$generalized) &&
      isTRUE(row$orthogonalization_native)
    preconditioner_gate <- if (identical(row$method, "eigencore_shifted_diagonal")) {
      identical(row$preconditioner_kind, "shifted_diagonal") &&
        isTRUE(row$preconditioner_native) &&
        isTRUE(row$preconditioner_calls > 0L)
    } else if (identical(row$method, "eigencore_shifted_tridiagonal")) {
      identical(row$preconditioner_kind, "shifted_tridiagonal") &&
        isTRUE(row$preconditioner_native) &&
        isTRUE(row$preconditioner_calls > 0L)
    } else {
      TRUE
    }
    constraint_gate <- if (identical(row$method, "eigencore_constrained")) {
      isTRUE(row$constrained) && isTRUE(row$constraints_rank > 0L)
    } else {
      TRUE
    }
    data.frame(
      case = row$case,
      method = row$method,
      requested = row$requested,
      nconv = row$nconv,
      native_gate = native_gate,
      preconditioner_gate = preconditioner_gate,
      constraint_gate = constraint_gate,
      passed = native_gate && preconditioner_gate && constraint_gate,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

generalized_lanczos_reference_contract <- function(rows) {
  internal <- rows[rows$method == "eigencore_lanczos_reference", , drop = FALSE]
  if (!nrow(internal)) {
    return(data.frame())
  }
  out <- lapply(seq_len(nrow(internal)), function(i) {
    row <- internal[i, , drop = FALSE]
    certificate_gate <- isTRUE(row$certificate_passed) && row$nconv >= row$requested
    label_gate <- identical(row$solver_label, eigencore:::generalized_lanczos_label())
    generalized_gate <- isTRUE(row$generalized)
    reference_gate <- !isTRUE(row$native) &&
      !isTRUE(row$native_kernels) &&
      isTRUE(row$metric_solves > 0L) &&
      isTRUE(nzchar(row$metric_solve_label))
    orthogonality_gate <- certificate_gate &&
      isTRUE(row$orthogonality_loss <= sqrt(.Machine$double.eps) * 10)
    data.frame(
      case = row$case,
      method = row$method,
      requested = row$requested,
      nconv = row$nconv,
      label_gate = label_gate,
      generalized_gate = generalized_gate,
      reference_gate = reference_gate,
      orthogonality_gate = orthogonality_gate,
      certificate_gate = certificate_gate,
      passed = label_gate && generalized_gate && reference_gate &&
        orthogonality_gate && certificate_gate,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

generalized_lobpcg_adversarial_b_specs <- function() {
  csc_n <- 10L
  csc_A <- Matrix::bandSparse(
    csc_n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, csc_n - 1L), rep(2.5, csc_n), rep(-1, csc_n - 1L))
  )
  csc_B <- methods::as(Matrix::bandSparse(
    csc_n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-0.05, csc_n - 1L),
      seq(1.2, 2, length.out = csc_n),
      rep(-0.05, csc_n - 1L)
    )
  ), "dgCMatrix")

  mf_n <- 9L
  mf_bdiag <- seq(0.75, 2.5, length.out = mf_n)
  mf_B_dense <- diag(mf_bdiag)
  mf_B <- linear_operator(
    dim = c(mf_n, mf_n),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (mf_bdiag * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (mf_bdiag * X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    metadata = list(
      frobenius_norm = sqrt(sum(mf_bdiag^2)),
      positive_definite = TRUE
    )
  )

  cases <- list(
    list(
      name = "ill_conditioned_diagonal_b",
      A = Matrix::Diagonal(x = c(1, 3, 7, 12, 20, 33, 50, 80)),
      B = Matrix::Diagonal(x = 10^seq(-3, 3, length.out = 8L)),
      B_dense = diag(10^seq(-3, 3, length.out = 8L)),
      expected_orthogonalization = "native_diagonal_b_mgs2"
    ),
    list(
      name = "sparse_csc_b",
      A = csc_A,
      B = csc_B,
      B_dense = as.matrix(csc_B),
      expected_orthogonalization = "native_csc_b_mgs2"
    ),
    list(
      name = "explicit_spd_matrix_free_b",
      A = Matrix::Diagonal(x = seq(2, 18, by = 2)),
      B = mf_B,
      B_dense = mf_B_dense,
      expected_orthogonalization = "native_matrix_free_b_mgs2"
    )
  )
  targets <- list(smallest = smallest(), largest = largest())
  specs <- list()
  for (case in cases) {
    for (target_name in names(targets)) {
      specs[[length(specs) + 1L]] <- list(
        case = paste("adversarial", case$name, target_name, sep = "_"),
        n = nrow(case$A),
        k = 2L,
        target = targets[[target_name]],
        methods = c("eigencore", "base"),
        subject = "eigencore",
        A = case$A,
        B = case$B,
        B_dense = case$B_dense,
        maxit = 220L,
        adversarial_b = TRUE,
        performance_gate = FALSE,
        expected_orthogonalization = case$expected_orthogonalization
      )
    }
  }
  specs
}

generalized_lobpcg_adversarial_b_contract <- function(rows, case_specs) {
  specs <- Filter(function(spec) isTRUE(spec$adversarial_b), case_specs)
  if (!length(specs)) {
    return(data.frame())
  }
  out <- lapply(specs, function(spec) {
    row <- rows[
      rows$case == spec$case & rows$method == "eigencore",
      ,
      drop = FALSE
    ]
    if (!nrow(row)) {
      return(data.frame(
        case = spec$case,
        target = eigencore:::target_label(spec$target),
        requested = spec$k,
        nconv = 0L,
        native_gate = FALSE,
        orthogonalization_gate = FALSE,
        certificate_gate = FALSE,
        passed = FALSE,
        stringsAsFactors = FALSE
      ))
    }
    row <- row[1L, , drop = FALSE]
    certificate_gate <- isTRUE(row$certificate_passed) && row$nconv >= spec$k
    native_gate <- certificate_gate &&
      isTRUE(row$native) &&
      isTRUE(row$native_kernels) &&
      isTRUE(row$generalized) &&
      isTRUE(row$orthogonalization_native)
    orthogonalization_gate <- identical(
      row$orthogonalization_methods,
      spec$expected_orthogonalization
    )
    data.frame(
      case = spec$case,
      target = row$target,
      requested = spec$k,
      nconv = row$nconv,
      native_gate = native_gate,
      orthogonalization_gate = orthogonalization_gate,
      certificate_gate = certificate_gate,
      passed = native_gate && orthogonalization_gate && certificate_gate,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

case_specs <- if (args$quick) {
  list(
    list(
      case = "sparse_generalized_path_smallest",
      n = 80L, k = 3L, target = smallest(), sparse = TRUE,
      methods = c(
        "eigencore",
        "eigencore_shifted_diagonal",
        "eigencore_shifted_tridiagonal",
        "base"
      ),
      subject = "eigencore"
    ),
    list(
      case = "sparse_generalized_path_largest",
      n = 80L, k = 3L, target = largest(), sparse = TRUE,
      methods = c(
        "eigencore",
        "eigencore_shifted_diagonal",
        "eigencore_shifted_tridiagonal",
        "base"
      ),
      subject = "eigencore"
    ),
    list(
      case = "diagonal_generalized_constrained_largest",
      n = 8L, k = 2L, target = largest(), sparse = FALSE,
      methods = c("eigencore_constrained", "base"),
      subject = "eigencore_constrained",
      A = Matrix::Diagonal(x = c(0, 1, 4, 9, 16, 25, 36, 49)),
      B = Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6, 7, 8)),
      constraints = matrix(c(1, rep(0, 7L)), ncol = 1L)
    ),
    list(
      case = "diagonal_generalized_lanczos_ref_smallest",
      n = 8L, k = 2L, target = smallest(), sparse = FALSE,
      methods = c("eigencore_lanczos_reference", "eigencore", "base"),
      subject = "eigencore_lanczos_reference",
      A = Matrix::Diagonal(x = c(1, 4, 9, 16, 25, 36, 49, 64)),
      B = Matrix::Diagonal(x = c(1, 2, 3, 4, 5, 6, 7, 8)),
      performance_gate = FALSE
    ),
    list(
      case = "sparse_csc_generalized_lanczos_ref_smallest",
      n = 10L, k = 2L, target = smallest(), sparse = TRUE,
      methods = c("eigencore_lanczos_reference", "eigencore", "base"),
      subject = "eigencore_lanczos_reference",
      A = Matrix::bandSparse(
        10L,
        k = c(-1, 0, 1),
        diagonals = list(rep(-1, 9L), rep(2.5, 10L), rep(-1, 9L))
      ),
      B = methods::as(Matrix::bandSparse(
        10L,
        k = c(-1, 0, 1),
        diagonals = list(
          rep(-0.05, 9L),
          seq(1.2, 2, length.out = 10L),
          rep(-0.05, 9L)
        )
      ), "dgCMatrix"),
      performance_gate = FALSE
    )
  )
} else {
  list(
    list(
      case = "sparse_generalized_path_smallest",
      n = 1000L, k = 10L, target = smallest(), sparse = TRUE,
      methods = c(
        "eigencore_shifted_tridiagonal",
        "base"
      ),
      subject = "eigencore_shifted_tridiagonal"
    ),
    list(
      case = "sparse_generalized_path_largest",
      n = 1000L, k = 10L, target = largest(), sparse = TRUE,
      methods = c(
        "eigencore_shifted_tridiagonal",
        "base"
      ),
      subject = "eigencore_shifted_tridiagonal",
      maxit = 300L
    ),
    list(
      case = "dense_generalized_partial_smallest",
      n = 180L, k = 8L, target = smallest(), sparse = FALSE,
      methods = c("eigencore_auto", "base"),
      subject = "eigencore_auto",
      performance_gate = FALSE
    ),
    list(
      case = "dense_generalized_partial_largest",
      n = 180L, k = 8L, target = largest(), sparse = FALSE,
      methods = c("eigencore_auto", "base"),
      subject = "eigencore_auto",
      performance_gate = FALSE
    ),
    list(
      case = "diagonal_generalized_constrained_largest",
      n = 80L, k = 8L, target = largest(), sparse = FALSE,
      methods = c("eigencore_constrained", "base"),
      subject = "eigencore_constrained",
      A = Matrix::Diagonal(x = c(0, seq(1, 79)^2)),
      B = Matrix::Diagonal(x = seq(1, 80)),
      constraints = matrix(c(1, rep(0, 79L)), ncol = 1L),
      performance_gate = FALSE
    ),
    list(
      case = "diagonal_generalized_lanczos_ref_smallest",
      n = 80L, k = 8L, target = smallest(), sparse = FALSE,
      methods = c("eigencore_lanczos_reference", "eigencore", "base"),
      subject = "eigencore_lanczos_reference",
      A = Matrix::Diagonal(x = seq_len(80L)^2),
      B = Matrix::Diagonal(x = seq_len(80L)),
      performance_gate = FALSE
    ),
    list(
      case = "sparse_csc_generalized_lanczos_ref_smallest",
      n = 80L, k = 8L, target = smallest(), sparse = TRUE,
      methods = c("eigencore_lanczos_reference", "eigencore", "base"),
      subject = "eigencore_lanczos_reference",
      A = Matrix::bandSparse(
        80L,
        k = c(-1, 0, 1),
        diagonals = list(rep(-1, 79L), rep(2.5, 80L), rep(-1, 79L))
      ),
      B = methods::as(Matrix::bandSparse(
        80L,
        k = c(-1, 0, 1),
        diagonals = list(
          rep(-0.05, 79L),
          seq(1.2, 2, length.out = 80L),
          rep(-0.05, 79L)
        )
      ), "dgCMatrix"),
      performance_gate = FALSE
    )
  )
}
case_specs <- c(case_specs, generalized_lobpcg_adversarial_b_specs())
case_specs <- filter_benchmark_cases(case_specs, args$cases)

rows <- lapply(seq_along(case_specs), function(i) {
  spec <- case_specs[[i]]
  message_benchmark_case("bench-generalized-lobpcg", spec)
  pair <- if (!is.null(spec$A) && !is.null(spec$B)) {
    list(A = spec$A, B = spec$B)
  } else {
    generalized_spd_pair(
      spec$n,
      rank = min(12L, spec$n),
      sparse = spec$sparse,
      seed = 12000L + spec$n + i
    )
  }
  methods <- spec$methods %||% c("eigencore", "base")
  if (!is.null(args$methods)) {
    methods <- intersect(methods, args$methods)
    if (!length(methods)) {
      stop(
        "--methods did not match any methods for case ",
        benchmark_case_id(spec),
        call. = FALSE
      )
    }
  }
  case_rows <- with_case(
    spec$case,
    eigencore:::target_label(spec$target),
    methods,
    spec$k,
    benchmark_generalized_lobpcg_case(
      pair$A,
      pair$B,
      k = spec$k,
      target = spec$target,
      methods = methods,
      iterations = iterations,
      tol = tol,
      seed = 12100L + spec$n + spec$k + i,
      constraints = spec$constraints %||% NULL,
      maxit = spec$maxit %||% 200L,
      A_dense = spec$A_dense %||% NULL,
      B_dense = spec$B_dense %||% NULL
    )
  )
  case_rows$adversarial_b <- isTRUE(spec$adversarial_b)
  case_rows$expected_orthogonalization <- spec$expected_orthogonalization %||% NA_character_
  case_rows
})
rows <- do.call(rbind, rows)
row.names(rows) <- NULL
case_by_name <- stats::setNames(case_specs, vapply(case_specs, `[[`, character(1), "case"))

gate_rows <- lapply(split(rows, rows$case), function(case_rows) {
  requested <- unique(case_rows$requested)
  spec <- case_by_name[[unique(case_rows$case)]]
  subject <- args$subject %||% spec$subject %||% "eigencore"
  if (!subject %in% unique(case_rows$method)) {
    stop(
      "--subject ", subject, " is not present for case ",
      unique(case_rows$case),
      call. = FALSE
    )
  }
  internal_methods <- c(
    "eigencore_auto",
    "eigencore",
    "eigencore_lanczos_reference",
    "eigencore_shifted_diagonal",
    "eigencore_shifted_tridiagonal",
    "eigencore_constrained"
  )
  gate <- evaluate_reference_gate(
    case_rows,
    subject = subject,
    references = setdiff(unique(case_rows$method), internal_methods),
    requested = requested[[1L]],
    speed_ratio_required = if (args$quick || isFALSE(spec$performance_gate %||% TRUE)) {
      0
    } else {
      release_speed_gate("generalized_eigen")
    },
    memory_ratio_required = if (args$quick || isFALSE(spec$performance_gate %||% TRUE)) {
      0
    } else {
      release_memory_gate("generalized_eigen")
    }
  )
  gate$case <- unique(case_rows$case)
  gate$target <- unique(case_rows$target)
  gate$kind <- "generalized_eigen"
  gate
})
gates <- do.call(rbind, gate_rows)
row.names(gates) <- NULL
native_contract <- generalized_lobpcg_native_contract(rows)
lanczos_contract <- generalized_lanczos_reference_contract(rows)
adversarial_contract <- generalized_lobpcg_adversarial_b_contract(rows, case_specs)

cat("Generalized SPD LOBPCG benchmark rows\n")
print(rows)
cat("\nGeneralized SPD LOBPCG gates\n")
print(gates)
cat("\nGeneralized SPD LOBPCG native contracts\n")
print(native_contract)
cat("\nGeneralized SPD B-orthogonal Lanczos reference contracts\n")
print(lanczos_contract)
cat("\nGeneralized SPD LOBPCG adversarial B contracts\n")
print(adversarial_contract)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "generalized-lobpcg-rows"))
  message("saved gates: ", save_benchmark_result(gates, "generalized-lobpcg-gates"))
  message("saved native contracts: ", save_benchmark_result(native_contract, "generalized-lobpcg-native-contracts"))
  message("saved generalized Lanczos contracts: ", save_benchmark_result(lanczos_contract, "generalized-lanczos-reference-contracts"))
  message("saved adversarial B contracts: ", save_benchmark_result(adversarial_contract, "generalized-lobpcg-adversarial-b-contracts"))
}

if (args$strict && (
  !all(gates$passed) ||
    !all(native_contract$passed) ||
    !all(lanczos_contract$passed) ||
    !all(adversarial_contract$passed)
)) {
  stop("Generalized SPD LOBPCG benchmark failed release gate.", call. = FALSE)
}
