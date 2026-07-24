#!/usr/bin/env Rscript

# Cold-versus-warm spectral continuation benchmark.
#
# Exercises dense and dgCMatrix Hermitian continuation across multiple sizes,
# requested ranks, spectral-gap scales, and tolerances. Every row reports wall
# time, exact operator block calls/columns, certification columns, restarts,
# native operator workspace bytes, and agreement with a common certified answer.
# Bounded cases are also checked against an independent dense eigen() oracle.
#
# Usage:
#   Rscript inst/benchmarks/bench-warm-start-continuation.R
#   Rscript inst/benchmarks/bench-warm-start-continuation.R --quick
#   Rscript inst/benchmarks/bench-warm-start-continuation.R --save --strict

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
base_seed <- 2026L
rho_seq <- c(0, 0.5, 1, 1.5, 2, 60, 60.5, 61)
overlap_loss_step <- 6L

as_dgc <- function(M) {
  if (methods::is(M, "dgCMatrix")) return(M)
  methods::as(methods::as(M, "generalMatrix"), "CsparseMatrix")
}

laplacian_1d <- function(n, gap_scale) {
  i <- c(seq_len(n), seq_len(n - 1L), 2:n)
  j <- c(seq_len(n), 2:n, seq_len(n - 1L))
  x <- gap_scale * c(rep(2, n), rep(-1, n - 1L), rep(-1, n - 1L))
  as_dgc(Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n)))
}

case_specs <- if (args$quick) {
  list(
    list(id = "csc_quick", storage = "dgCMatrix", n = 200L, k = 6L,
         gap_scale = 1, tol = 1e-8, reps = 2L),
    list(id = "dense_quick", storage = "dense", n = 160L, k = 4L,
         gap_scale = 4, tol = 1e-6, reps = 2L)
  )
} else {
  list(
    list(id = "csc_n1000_k10_tight", storage = "dgCMatrix",
         n = 1000L, k = 10L, gap_scale = 1, tol = 1e-8, reps = 5L),
    list(id = "csc_n2000_k20_loose", storage = "dgCMatrix",
         n = 2000L, k = 20L, gap_scale = 4, tol = 1e-6, reps = 5L),
    list(id = "dense_n600_k8_tight", storage = "dense",
         n = 600L, k = 8L, gap_scale = 1, tol = 1e-8, reps = 5L),
    list(id = "dense_n1000_k12_loose", storage = "dense",
         n = 1000L, k = 12L, gap_scale = 4, tol = 1e-6, reps = 5L)
  )
}

median_time <- function(thunk, reps) {
  thunk()
  stats::median(vapply(seq_len(reps), function(i) {
    as.numeric(system.time(thunk())[["elapsed"]])
  }, numeric(1)))
}

subspace_overlap <- function(X, Y) {
  qx <- qr.Q(qr(as.matrix(X)))
  qy <- qr.Q(qr(as.matrix(Y)))
  min(svd(crossprod(qx, qy), nu = 0, nv = 0)$d)
}

certified <- function(fit) {
  cert <- certificate(fit)
  isTRUE(cert$passed) ||
    (all(cert$converged) &&
       cert$max_backward_error <= cert$tolerance &&
       isTRUE(cert$orthogonality_passed))
}

all_rows <- list()
row_index <- 0L
trace_fixture <- NULL

for (case_index in seq_along(case_specs)) {
  spec <- case_specs[[case_index]]
  n <- spec$n
  k <- spec$k
  set.seed(base_seed + case_index)

  L_sparse <- laplacian_1d(n, spec$gap_scale)
  grid <- seq_len(n) / n
  potential <- 1 + grid + 0.5 * sin(4 * pi * grid)
  V_sparse <- Matrix::sparseMatrix(
    i = seq_len(n), j = seq_len(n), x = potential, dims = c(n, n)
  )
  shifted_operator <- if (identical(spec$storage, "dense")) {
    L <- as.matrix(L_sparse)
    V <- as.matrix(V_sparse)
    function(rho) L - rho * V
  } else {
    function(rho) as_dgc(L_sparse - rho * V_sparse)
  }
  if (is.null(trace_fixture) && identical(spec$storage, "dgCMatrix")) {
    trace_fixture <- list(spec = spec, shifted_operator = shifted_operator)
  }

  method <- lanczos(
    block = k,
    max_subspace = min(n - 1L, 8L * k),
    max_restarts = 400L
  )
  solve_once <- function(op, start = NULL) {
    eig_partial(
      op, k = k, target = smallest(), method = method, tol = spec$tol,
      seed = base_seed, initial_subspace = start
    )
  }
  oracle_values <- function(op) {
    if (n > 1200L) return(NULL)
    sort(eigen(as.matrix(op), symmetric = TRUE)$values)[seq_len(k)]
  }

  prev_vectors <- NULL
  for (step in seq_along(rho_seq)) {
    rho <- rho_seq[[step]]
    op <- shifted_operator(rho)
    message_benchmark_case(
      "bench-warm-start-continuation",
      list(case = spec$id, n = n, k = k)
    )

    cold <- solve_once(op)
    warm <- if (is.null(prev_vectors)) NULL else solve_once(op, prev_vectors)
    cold_sec <- median_time(function() solve_once(op), spec$reps)
    warm_sec <- if (is.null(prev_vectors)) NA_real_ else {
      start <- prev_vectors
      median_time(function() solve_once(op, start), spec$reps)
    }

    cold_values <- sort(values(cold))
    warm_values <- if (is.null(warm)) NULL else sort(values(warm))
    oracle <- oracle_values(op)
    overlap <- if (is.null(prev_vectors)) NA_real_ else {
      subspace_overlap(prev_vectors, vectors(cold))
    }
    oracle_gap <- if (is.null(oracle) || length(oracle) < k) {
      NA_real_
    } else {
      full <- sort(eigen(as.matrix(op), symmetric = TRUE,
                         only.values = TRUE)$values)
      full[[k + 1L]] - full[[k]]
    }

    row_index <- row_index + 1L
    all_rows[[row_index]] <- data.frame(
      case = spec$id,
      storage = spec$storage,
      n = n,
      k = k,
      gap_scale = spec$gap_scale,
      observed_gap = oracle_gap,
      tol = spec$tol,
      step = step,
      rho = rho,
      overlap_loss = identical(step, overlap_loss_step),
      start_overlap = overlap,
      cold_certified = certified(cold),
      warm_certified = if (is.null(warm)) NA else certified(warm),
      warm_start_source = if (is.null(warm)) "cold" else warm$start_source,
      cold_operator_block_calls = cold$operator_block_calls,
      warm_operator_block_calls =
        if (is.null(warm)) NA_integer_ else warm$operator_block_calls,
      cold_operator_columns = cold$operator_columns,
      warm_operator_columns =
        if (is.null(warm)) NA_integer_ else warm$operator_columns,
      operator_column_ratio = if (is.null(warm)) NA_real_ else {
        warm$operator_columns / max(cold$operator_columns, 1L)
      },
      cold_certification_columns = cold$certification_operator_columns,
      warm_certification_columns = if (is.null(warm)) {
        NA_integer_
      } else {
        warm$certification_operator_columns
      },
      cold_restarts = cold$restart$restarts_used %||% 0L,
      warm_restarts = if (is.null(warm)) {
        NA_integer_
      } else {
        warm$restart$restarts_used %||% 0L
      },
      cold_operator_workspace_bytes =
        cold$restart$operator_bytes_allocated %||% NA_real_,
      warm_operator_workspace_bytes = if (is.null(warm)) {
        NA_real_
      } else {
        warm$restart$operator_bytes_allocated %||% NA_real_
      },
      cold_sec = cold_sec,
      warm_sec = warm_sec,
      speedup = if (is.na(warm_sec)) NA_real_ else cold_sec / warm_sec,
      agree_cold_warm = if (is.null(warm_values)) {
        NA_real_
      } else {
        max(abs(cold_values - warm_values))
      },
      agree_oracle = if (is.null(oracle)) {
        NA_real_
      } else {
        max(abs(cold_values - oracle))
      },
      stringsAsFactors = FALSE
    )
    prev_vectors <- vectors(if (is.null(warm)) cold else warm)
  }
}

rows <- do.call(rbind, all_rows)
print(rows, row.names = FALSE, digits = 4)

# Exported-API traceratio-shaped continuation row; no internal access.
trace_fun <- function(fixture) {
  spec <- fixture$spec
  start <- NULL
  total_columns <- 0L
  total_certification_columns <- 0L
  certified_rows <- logical(0)
  for (rho in c(0, 0.5, 1)) {
    fit <- eigencore::eig_partial(
      fixture$shifted_operator(rho),
      k = spec$k,
      target = eigencore::smallest(),
      method = eigencore::lanczos(
        block = spec$k,
        max_subspace = min(spec$n - 1L, 8L * spec$k),
        max_restarts = 400L
      ),
      tol = spec$tol,
      seed = base_seed,
      initial_subspace = start
    )
    start <- eigencore::vectors(fit)
    total_columns <- total_columns + fit$operator_columns
    total_certification_columns <- total_certification_columns +
      fit$certification_operator_columns
    certified_rows <- c(certified_rows, certified(fit))
  }
  list(
    operator_columns = total_columns,
    certification_operator_columns = total_certification_columns,
    all_certified = all(certified_rows)
  )
}
trace <- tryCatch(
  trace_fun(trace_fixture),
  error = function(e) list(
    operator_columns = NA_integer_,
    certification_operator_columns = NA_integer_,
    all_certified = FALSE,
    error = conditionMessage(e)
  )
)
cat(
  "\ntraceratio integration: operator_columns=", trace$operator_columns,
  " certification_operator_columns=", trace$certification_operator_columns,
  " all_certified=", trace$all_certified, "\n", sep = ""
)

if (args$save) {
  message("saved rows: ",
          save_benchmark_result(rows, "warm-start-continuation-rows"))
}

warm_rows <- rows[!is.na(rows$agree_cold_warm), ]
loss_rows <- rows[rows$overlap_loss, ]
correctness_ok <-
  all(rows$cold_certified) &&
  all(rows$warm_certified[!is.na(rows$warm_certified)]) &&
  all(warm_rows$agree_cold_warm < 1e-6) &&
  all(rows$agree_oracle[!is.na(rows$agree_oracle)] < 1e-6) &&
  isTRUE(trace$all_certified)
overlap_ok <-
  nrow(loss_rows) == length(case_specs) &&
  all(loss_rows$warm_certified) &&
  all(loss_rows$agree_cold_warm < 1e-6) &&
  all(loss_rows$start_overlap < 0.5)
accounting_ok <-
  all(rows$cold_operator_columns >= rows$cold_operator_block_calls) &&
  all(warm_rows$warm_operator_columns >=
        warm_rows$warm_operator_block_calls) &&
  all(rows$cold_certification_columns <= rows$cold_operator_columns) &&
  all(warm_rows$warm_certification_columns <=
        warm_rows$warm_operator_columns)

cat(
  "\ncorrectness_ok=", correctness_ok,
  " overlap_loss_real_and_honest=", overlap_ok,
  " accounting_ok=", accounting_ok, "\n", sep = ""
)
if (args$strict && !(correctness_ok && overlap_ok && accounting_ok)) {
  stop("Warm-start continuation benchmark failed strict invariants.",
       call. = FALSE)
}
