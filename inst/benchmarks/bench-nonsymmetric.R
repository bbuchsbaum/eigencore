#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (!is.na(args$iterations)) args$iterations else if (args$quick) 1L else 5L
tol <- 1e-10

nonsymmetric_oracle_label <- "dense LAPACK general eigen oracle (prototype fallback)"
nonsymmetric_reference_arnoldi_label <- "reference Arnoldi (prototype/oracle fallback)"
nonsymmetric_native_arnoldi_label <- "native Arnoldi cycle + native Ritz extraction (compatibility)"

nonsymmetric_real_nonnormal <- function(n) {
  values <- seq(n, 1)
  A <- diag(values, n)
  A[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- 25
  A
}

nonsymmetric_complex_blocks <- function(n) {
  A <- matrix(0, n, n)
  A[1:2, 1:2] <- matrix(c(0, -2, 2, 0), 2, byrow = TRUE)
  if (n >= 4L) {
    A[3:4, 3:4] <- matrix(c(0, -1, 1, 0), 2, byrow = TRUE)
  }
  if (n > 4L) {
    diag(A)[5:n] <- seq(0.2, 0.2 + 0.1 * (n - 5L), length.out = n - 4L)
  }
  A
}

nonsymmetric_sparse_real <- function(n) {
  A <- Matrix::sparseMatrix(
    i = c(seq_len(n), seq_len(n - 1L)),
    j = c(seq_len(n), seq_len(n - 1L) + 1L),
    x = c(seq(n, 1), rep(2, n - 1L)),
    dims = c(n, n)
  )
  methods::as(A, "dgCMatrix")
}

nonsymmetric_sparse_complex_blocks <- function(n) {
  methods::as(Matrix::Matrix(nonsymmetric_complex_blocks(n), sparse = TRUE), "dgCMatrix")
}

nonsymmetric_cases <- function(quick = FALSE) {
  if (quick) {
    return(list(
      list(
        case = "dense_native_arnoldi_lm",
        n = 8L, k = 3L, api = "eig_partial",
        target = largest_magnitude(),
        which = NA_character_,
        build = nonsymmetric_real_nonnormal
      ),
      list(
        case = "dense_native_arnoldi_li",
        n = 6L, k = 2L, api = "eig_partial",
        target = largest_imaginary(),
        which = NA_character_,
        build = nonsymmetric_complex_blocks
      ),
      list(
        case = "dense_eigs_native_arnoldi_li",
        n = 6L, k = 2L, api = "eigs",
        target = largest_imaginary(),
        which = "LI",
        build = nonsymmetric_complex_blocks
      ),
      list(
        case = "sparse_native_arnoldi_lr",
        n = 8L, k = 3L, api = "eig_partial",
        target = largest_real(),
        which = NA_character_,
        build = nonsymmetric_sparse_real
      ),
      list(
        case = "sparse_native_arnoldi_li",
        n = 8L, k = 2L, api = "eigs",
        target = largest_imaginary(),
        which = "LI",
        build = nonsymmetric_sparse_complex_blocks
      )
    ))
  }

  list(
    list(
      case = "dense_native_arnoldi_lm",
      n = 80L, k = 8L, api = "eig_partial",
      target = largest_magnitude(),
      which = NA_character_,
      build = nonsymmetric_real_nonnormal
    ),
    list(
      case = "dense_native_arnoldi_li",
      n = 40L, k = 4L, api = "eig_partial",
      target = largest_imaginary(),
      which = NA_character_,
      build = nonsymmetric_complex_blocks
    ),
    list(
      case = "dense_eigs_native_arnoldi_li",
      n = 40L, k = 4L, api = "eigs",
      target = largest_imaginary(),
      which = "LI",
      build = nonsymmetric_complex_blocks
    ),
    list(
      case = "sparse_native_arnoldi_lr",
      n = 80L, k = 8L, api = "eig_partial",
      target = largest_real(),
      which = NA_character_,
      build = nonsymmetric_sparse_real
    ),
    list(
      case = "sparse_native_arnoldi_li",
      n = 40L, k = 2L, api = "eigs",
      target = largest_imaginary(),
      which = "LI",
      build = nonsymmetric_sparse_complex_blocks
    )
  )
}

run_nonsymmetric_case <- function(case) {
  A <- case$build(case$n)
  if (identical(case$api, "eigs")) {
    eigs(A, k = case$k, which = case$which, tol = tol)
  } else {
    eig_partial(A, k = case$k, target = case$target, tol = tol)
  }
}

nonsymmetric_fit_diagnostics <- function(fit) {
  if (!is.null(fit$diagnostics) && is.list(fit$diagnostics)) {
    fit$diagnostics
  } else {
    diagnostics(fit)
  }
}

nonsymmetric_fit_values <- function(fit) {
  fit$values
}

nonsymmetric_fit_vectors <- function(fit) {
  fit$vectors
}

nonsymmetric_stage_second <- function(stage_seconds, name) {
  if (!length(stage_seconds) || is.null(names(stage_seconds)) ||
      !name %in% names(stage_seconds)) {
    return(NA_real_)
  }
  unname(stage_seconds[[name]])
}

benchmark_nonsymmetric_case <- function(case, iterations = 3L, seed = 1L) {
  timed <- run_timed(run_nonsymmetric_case(case), iterations = iterations, seed = seed)
  fit <- timed$value
  cert <- fit$certificate
  diag <- nonsymmetric_fit_diagnostics(fit)
  restart <- diag$restart %||% list()
  stage_seconds <- diag$stage_seconds %||% numeric()
  vals <- nonsymmetric_fit_values(fit)
  vecs <- nonsymmetric_fit_vectors(fit)
  data.frame(
    case = case$case,
    kind = "nonsymmetric_eigen",
    api = case$api,
    target = eigencore:::target_label(case$target),
    which = case$which,
    method = diag$method,
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
    orthogonality_required = cert$orthogonality_required,
    nconv = sum(cert$converged),
    requested = case$k,
    iterations = diag$iterations %||% NA_integer_,
    matvecs = diag$matvecs %||% NA_integer_,
    restart_count = restart$restart_count %||% NA_integer_,
    max_restarts = restart$max_restarts %||% NA_integer_,
    certified_attempt = restart$certified_attempt %||% NA_integer_,
    selected_attempt = restart$selected_attempt %||% NA_integer_,
    stage_arnoldi_cycle_seconds =
      nonsymmetric_stage_second(stage_seconds, "cycle"),
    stage_ritz_extraction_seconds =
      nonsymmetric_stage_second(stage_seconds, "ritz_extraction"),
    ritz_extraction_native = restart$ritz_extraction_native %||% NA,
    value_real_1 = Re(vals[[1L]]),
    value_imag_1 = Im(vals[[1L]]),
    complex_values = is.complex(vals),
    complex_vectors = is.complex(vecs),
    right_residual_certified = identical(cert$certificate_type, "right_residual_backward_error"),
    dense_oracle_label = identical(diag$method, nonsymmetric_oracle_label),
    reference_arnoldi_label = identical(diag$method, nonsymmetric_reference_arnoldi_label),
    native_arnoldi_label = identical(diag$method, nonsymmetric_native_arnoldi_label),
    arnoldi_native = grepl("native.*arnoldi|arnoldi.*native", diag$method, ignore.case = TRUE),
    warnings = paste(diag$warnings %||% character(), collapse = "; "),
    seed = seed,
    pkg_version = as.character(utils::packageVersion("eigencore")),
    stringsAsFactors = FALSE
  )
}

nonsymmetric_contract <- function(rows) {
  out <- lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, , drop = FALSE]
    certificate_gate <- isTRUE(row$certificate_passed) && row$nconv >= row$requested
    right_residual_gate <- isTRUE(row$right_residual_certified) &&
      !isTRUE(row$orthogonality_required)
    label_gate <- isTRUE(row$dense_oracle_label) ||
      isTRUE(row$reference_arnoldi_label) ||
      isTRUE(row$native_arnoldi_label)
    warning_gate <- (
      isTRUE(row$dense_oracle_label) &&
        grepl("dense general eigen oracle", row$warnings, fixed = TRUE) &&
        grepl("right residuals certified", row$warnings, fixed = TRUE)
    ) || (
      isTRUE(row$reference_arnoldi_label) &&
        grepl("reference Arnoldi prototype", row$warnings, fixed = TRUE) &&
        grepl("native nonsymmetric Arnoldi not yet implemented", row$warnings, fixed = TRUE)
    ) || (
      isTRUE(row$native_arnoldi_label) &&
        grepl("native Arnoldi cycle", row$warnings, fixed = TRUE) &&
        grepl("right residuals certified", row$warnings, fixed = TRUE)
    )
    data.frame(
      case = row$case,
      api = row$api,
      requested = row$requested,
      nconv = row$nconv,
      certificate_gate = certificate_gate,
      right_residual_gate = right_residual_gate,
      label_gate = label_gate,
      reference_arnoldi_label = isTRUE(row$reference_arnoldi_label),
      native_arnoldi_label = isTRUE(row$native_arnoldi_label),
      arnoldi_native = isTRUE(row$arnoldi_native),
      restart_gate = if (isTRUE(row$native_arnoldi_label)) {
        !is.na(row$max_restarts) && row$max_restarts >= 1L &&
          !is.na(row$restart_count) &&
          isTRUE(row$ritz_extraction_native)
      } else {
        TRUE
      },
      warning_gate = warning_gate,
      passed = certificate_gate && right_residual_gate && label_gate &&
        warning_gate &&
        (if (isTRUE(row$native_arnoldi_label)) {
          !is.na(row$max_restarts) && row$max_restarts >= 1L &&
            !is.na(row$restart_count) &&
            isTRUE(row$ritz_extraction_native)
        } else {
          TRUE
        }),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

cases <- filter_benchmark_cases(nonsymmetric_cases(args$quick), args$cases)
rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  message_benchmark_case("bench-nonsymmetric", case)
  benchmark_nonsymmetric_case(
    case,
    iterations = iterations,
    seed = 15000L + case$n + case$k + i
  )
})
rows <- do.call(rbind, rows)
row.names(rows) <- NULL
contracts <- nonsymmetric_contract(rows)

cat("Nonsymmetric benchmark rows\n")
print(rows)
cat("\nNonsymmetric contracts\n")
print(contracts)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "nonsymmetric-rows"))
  message("saved contracts: ", save_benchmark_result(contracts, "nonsymmetric-contracts"))
}

if (args$strict && !all(contracts$passed)) {
  stop("Nonsymmetric benchmark failed compatibility contract gate.", call. = FALSE)
}
