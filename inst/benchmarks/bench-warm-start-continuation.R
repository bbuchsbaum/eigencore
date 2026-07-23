#!/usr/bin/env Rscript

# Cold-versus-warm spectral continuation benchmark (V2.1 warm-start seam).
#
# Models the traceratio access pattern: repeatedly solve the k smallest
# eigenpairs of a slowly drifting standard Hermitian operator H(rho) = L - rho*V,
# a 1-D discrete Laplacian L plus a smooth diagonal potential V. Because L has
# delocalized (sinusoidal) eigenvectors, increasing rho genuinely rotates and
# localizes the leading eigenspace, so consecutive small-rho steps overlap
# strongly (warm start pays off) while a large-rho jump destroys the overlap
# (the adversary: warm must stay correct and not claim a false saving).
#
# Each step is solved twice:
#   * cold: a fresh random start;
#   * warm: initial_subspace = the previous step's computed eigenvectors.
# We report wall time AND operator columns (matvecs) to a common certified
# answer, cross-checked against an independent dense eigen() oracle. Per
# planning/prd.json no performance claim is made from a single quick run;
# --strict gates the correctness invariants only.
#
# Usage:
#   Rscript inst/benchmarks/bench-warm-start-continuation.R            # non-quick
#   Rscript inst/benchmarks/bench-warm-start-continuation.R --quick
#   Rscript inst/benchmarks/bench-warm-start-continuation.R --save --strict

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
n    <- if (args$quick) 200L else 2000L
k    <- if (args$quick) 6L else 20L
reps <- if (args$quick) 3L else 5L
seed <- 2026L

set.seed(seed)

as_dgc <- function(M) {
  if (methods::is(M, "dgCMatrix")) return(M)
  methods::as(methods::as(M, "generalMatrix"), "CsparseMatrix")
}

# 1-D discrete Laplacian: symmetric, delocalized (sinusoidal) eigenvectors.
laplacian_1d <- function(n) {
  i <- c(seq_len(n), seq_len(n - 1L), 2:n)
  j <- c(seq_len(n), 2:n, seq_len(n - 1L))
  x <- c(rep(2, n), rep(-1, n - 1L), rep(-1, n - 1L))
  as_dgc(Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n)))
}

L <- laplacian_1d(n)
# Strictly positive smooth ramp-plus-bump potential; a large rho localizes the
# lowest eigenstates near its maximum, far from L's smooth modes.
grid <- seq_len(n) / n
Vdiag <- 1 + grid + 0.5 * sin(4 * pi * grid)          # in (0.5, 2.5), > 0
V <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = Vdiag, dims = c(n, n))

shifted_operator <- function(rho) as_dgc(L - rho * V)

# rho schedule: a slow drift, a deliberate overlap-loss jump, then resume.
rho_seq <- c(0, 0.5, 1.0, 1.5, 2.0, 60, 60.5, 61)
overlap_loss_step <- 6L

target <- smallest()
method <- lanczos(block = k, max_subspace = min(n - 1L, 8L * k), max_restarts = 400L)

median_time <- function(thunk, reps) {
  thunk()                                       # warmup
  ts <- vapply(seq_len(reps), function(i) {
    as.numeric(system.time(thunk())[["elapsed"]])
  }, numeric(1))
  stats::median(ts)
}

solve_once <- function(op, start = NULL) {
  eig_partial(op, k = k, target = target, method = method, seed = seed,
              initial_subspace = start)
}

sorted_vals <- function(fit) sort(values(fit))     # ascending: k smallest

oracle_values <- function(op) {
  if (n > 2500L) return(NULL)
  sort(eigen(as.matrix(op), symmetric = TRUE)$values)[seq_len(k)]
}

# Subspace overlap between two orthonormal-ish bases (largest principal angle
# cosine); ~1 means high overlap, ~0 means the warm start is useless.
subspace_overlap <- function(X, Y) {
  qx <- qr.Q(qr(as.matrix(X)))
  qy <- qr.Q(qr(as.matrix(Y)))
  min(svd(crossprod(qx, qy), nu = 0, nv = 0)$d)
}

rows <- vector("list", length(rho_seq))
prev_vectors <- NULL

for (t in seq_along(rho_seq)) {
  rho <- rho_seq[[t]]
  op <- shifted_operator(rho)
  message_benchmark_case("bench-warm-start-continuation",
                         list(case = sprintf("step %d (rho=%.2f)", t, rho)))

  cold <- solve_once(op, start = NULL)
  warm <- if (is.null(prev_vectors)) NULL else solve_once(op, start = prev_vectors)

  cold_time <- median_time(function() solve_once(op, start = NULL), reps)
  warm_time <- if (is.null(prev_vectors)) NA_real_ else {
    pv <- prev_vectors
    median_time(function() solve_once(op, start = pv), reps)
  }

  oracle <- oracle_values(op)
  cold_vals <- sorted_vals(cold)
  warm_vals <- if (is.null(warm)) NULL else sorted_vals(warm)

  overlap <- if (is.null(prev_vectors)) NA_real_ else {
    subspace_overlap(prev_vectors, vectors(cold))
  }

  rows[[t]] <- data.frame(
    step = t,
    rho = rho,
    overlap_loss = identical(t, overlap_loss_step),
    start_overlap = overlap,
    cold_certified = certificate(cold)$passed,
    warm_certified = if (is.null(warm)) NA else certificate(warm)$passed,
    warm_start_source = if (is.null(warm)) "cold" else warm$start_source,
    cold_matvecs = cold$matvecs,          # operator columns to certified answer
    warm_matvecs = if (is.null(warm)) NA_integer_ else warm$matvecs,
    matvec_ratio = if (is.null(warm)) NA_real_ else warm$matvecs / max(cold$matvecs, 1L),
    cold_sec = cold_time,
    warm_sec = warm_time,
    speedup = if (is.na(warm_time)) NA_real_ else cold_time / warm_time,
    agree_cold_warm = if (is.null(warm_vals)) NA_real_ else max(abs(cold_vals - warm_vals)),
    agree_oracle = if (is.null(oracle)) NA_real_ else max(abs(cold_vals - oracle)),
    stringsAsFactors = FALSE
  )

  prev_vectors <- vectors(if (is.null(warm)) cold else warm)
}

rows <- do.call(rbind, rows)

cat("\nCold vs warm spectral continuation (L - rho*V), n=", n, " k=", k, "\n", sep = "")
print(rows, row.names = FALSE, digits = 4)

# --- traceratio integration row: exported API only, no ::: --------------------
trace_fun <- function() {
  Vs <- NULL
  total_matvecs <- 0L
  certified <- logical(0)
  for (rho in c(0, 0.5, 1.0)) {
    op <- as_dgc(L - rho * V)
    fit <- eigencore::eig_partial(
      op, k = k, target = eigencore::smallest(),
      method = eigencore::lanczos(block = k, max_subspace = min(n - 1L, 8L * k)),
      seed = seed, initial_subspace = Vs
    )
    Vs <- eigencore::vectors(fit)                 # exported accessor, no :::
    total_matvecs <- total_matvecs + fit$matvecs
    certified <- c(certified, eigencore::certificate(fit)$passed)
  }
  list(matvecs = total_matvecs, all_certified = all(certified))
}
trace <- tryCatch(trace_fun(),
                  error = function(e) list(matvecs = NA, all_certified = FALSE,
                                           error = conditionMessage(e)))
cat("\ntraceratio integration (exported API, warm continuation): ",
    "matvecs=", trace$matvecs, " all_certified=", trace$all_certified, "\n", sep = "")

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "warm-start-continuation-rows"))
}

# --- Correctness gate (the only thing --strict enforces) ----------------------
warm_rows <- rows[!is.na(rows$agree_cold_warm), ]
correctness_ok <-
  all(rows$cold_certified) &&
  all(rows$warm_certified[!is.na(rows$warm_certified)]) &&
  all(warm_rows$agree_cold_warm < 1e-6) &&
  all(rows$agree_oracle[!is.na(rows$agree_oracle)] < 1e-6) &&
  isTRUE(trace$all_certified)

# Overlap-loss honesty: at the jump, overlap must actually be low, yet warm must
# still certify and agree with cold on the answer.
loss_row <- rows[rows$overlap_loss, ]
overlap_ok <- nrow(loss_row) == 1L &&
  isTRUE(loss_row$warm_certified) &&
  loss_row$agree_cold_warm < 1e-6 &&
  isTRUE(loss_row$start_overlap < 0.5)

cat("\ncorrectness_ok=", correctness_ok, " overlap_loss_real_and_honest=", overlap_ok,
    "\n", sep = "")

if (args$strict && !(correctness_ok && overlap_ok)) {
  stop("Warm-start continuation benchmark failed its correctness invariants.",
       call. = FALSE)
}
