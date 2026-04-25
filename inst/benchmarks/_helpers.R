library(eigencore)

`%||%` <- function(x, y) if (is.null(x)) y else x

benchmark_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  list(
    quick = "--quick" %in% args,
    save = "--save" %in% args,
    strict = "--strict" %in% args
  )
}

release_speed_gate <- function(kind) {
  switch(
    kind,
    hermitian = 1.25,
    svd = 1.5,
    randomized_svd = 2.0,
    generalized_eigen = 1.0,
    lobpcg = 1.0,
    1.0
  )
}

release_memory_gate <- function(kind) {
  switch(
    kind,
    generalized_eigen = 0.25,
    lobpcg = 0.25,
    1.0
  )
}

path_laplacian <- function(n) {
  Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
}

tall_skinny_sparse <- function(m, n, density = 0.01, seed = 1L) {
  set.seed(seed)
  Matrix::rsparsematrix(m, n, density = density)
}

dense_low_rank_spd <- function(n, rank = 20L, diagonal_shift = 1, seed = 1L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * rank), nrow = n, ncol = rank)
  tcrossprod(X) + diag(diagonal_shift, n)
}

slow_decay_svd_matrix <- function(m, n, decay = 0.35, seed = 1L) {
  set.seed(seed)
  r <- min(m, n)
  U <- qr.Q(qr(matrix(stats::rnorm(m * r), nrow = m, ncol = r)))
  V <- qr.Q(qr(matrix(stats::rnorm(n * r), nrow = n, ncol = r)))
  d <- seq_len(r)^(-decay)
  U %*% (d * t(V))
}

generalized_spd_pair <- function(n, rank = 12L, sparse = FALSE, seed = 1L) {
  set.seed(seed)
  if (isTRUE(sparse)) {
    A <- path_laplacian(n) + Matrix::Diagonal(n, 0.25)
    weights <- stats::runif(n, min = 0.5, max = 2)
    B <- Matrix::Diagonal(n, weights)
  } else {
    X <- matrix(stats::rnorm(n * rank), nrow = n, ncol = rank)
    Y <- matrix(stats::rnorm(n * rank), nrow = n, ncol = rank)
    A <- tcrossprod(X) + diag(0.5, n)
    B <- tcrossprod(Y) + diag(1, n)
  }
  list(A = A, B = B)
}

bench_methods <- function(kind, requested = NULL) {
  methods <- if (identical(kind, "eigen")) {
    eigencore:::available_eigen_methods()
  } else {
    eigencore:::available_svd_methods()
  }
  if (!is.null(requested)) intersect(requested, methods) else methods
}

run_timed <- function(expr, iterations = 3L) {
  expr <- substitute(expr)
  env <- parent.frame()
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("bench is required for release benchmarks.", call. = FALSE)
  }
  value <- NULL
  mark <- bench::mark(
    value <- eval(expr, envir = env),
    iterations = iterations,
    check = FALSE,
    time_unit = "s",
    memory = TRUE,
    filter_gc = FALSE
  )
  list(
    value = value,
    median = stats::median(as.numeric(mark$time[[1L]])),
    min = min(as.numeric(mark$time[[1L]])),
    mem_alloc = as.numeric(mark$mem_alloc[[1L]]),
    bench = mark
  )
}

certify_eigen_result <- function(A, fit, tol = 1e-8) {
  vals <- eigencore:::method_values(fit, kind = "eigen")
  vecs <- fit$vectors
  if (is.null(vals) || is.null(vecs)) return(eigencore:::empty_certificate(tol, "vectors unavailable"))
  eigencore:::certify_eigen_operator(as_operator(A), vals, vecs, tol = tol)
}

certify_svd_result <- function(A, fit, tol = 1e-8) {
  d <- eigencore:::method_values(fit, kind = "svd")
  u <- fit$u
  v <- fit$v
  if (is.null(d) || is.null(u) || is.null(v)) {
    return(eigencore:::empty_certificate(tol, "both singular-vector sides unavailable"))
  }
  eigencore:::certify_svd_operator(as_operator(A), d, u, v, tol = tol)
}

result_preconditioner_field <- function(x, field) {
  info <- x$preconditioner %||% x$restart$preconditioner %||% NULL
  if (is.null(info) || is.null(info[[field]])) {
    return(NA)
  }
  info[[field]]
}

result_preconditioner_calls <- function(x) {
  x$preconditioner_calls %||% x$restart$preconditioner_calls %||% NA_integer_
}

benchmark_eigen_case <- function(A, k, target = largest(), methods = NULL,
                                 iterations = 3L, tol = 1e-8, seed = 1L) {
  methods <- bench_methods("eigen", methods)
  rows <- lapply(methods, function(method) {
    set.seed(seed)
    timed <- run_timed(eigencore:::run_eigen_method(method, A, k = k, target = target, tol = tol),
                       iterations = iterations)
    cert <- certify_eigen_result(A, timed$value, tol = tol)
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
      preconditioner_kind = result_preconditioner_field(timed$value, "kind"),
      preconditioner_native = result_preconditioner_field(timed$value, "native"),
      preconditioner_calls = result_preconditioner_calls(timed$value),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

benchmark_generalized_eigen_case <- function(A, B, k, target = smallest(),
                                             methods = "eigencore",
                                             iterations = 3L, tol = 1e-8,
                                             seed = 1L) {
  methods <- intersect(methods, c("eigencore", "base"))
  rows <- lapply(methods, function(method) {
    set.seed(seed)
    timed <- run_timed({
      if (identical(method, "eigencore")) {
        eig_partial(A, B = B, k = k, target = target, tol = tol)
      } else {
        eig <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
        idx <- eigencore:::order_indices(eig$values, target)
        idx <- idx[seq_len(k)]
        list(values = eig$values[idx], vectors = eig$vectors[, idx, drop = FALSE])
      }
    }, iterations = iterations)
    cert <- eigencore:::certify_eigen(as.matrix(A), eigencore:::method_values(timed$value, kind = "eigen"),
                                      timed$value$vectors, B = as.matrix(B), tol = tol)
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
      preconditioner_kind = result_preconditioner_field(timed$value, "kind"),
      preconditioner_native = result_preconditioner_field(timed$value, "native"),
      preconditioner_calls = result_preconditioner_calls(timed$value),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

evaluate_native_hermitian_gate <- function(rows, k,
                                           speed_ratio_required = release_speed_gate("hermitian"),
                                           memory_ratio_required = 1.0) {
  eig <- rows[rows$method == "eigencore", , drop = FALSE]
  refs <- rows[rows$method != "eigencore" & rows$certificate_passed, , drop = FALSE]
  if (nrow(eig) != 1L) {
    stop("gate requires exactly one eigencore row", call. = FALSE)
  }
  if (!nrow(refs)) {
    stop("gate requires at least one certified reference row", call. = FALSE)
  }

  best_time <- min(refs$median)
  best_mem <- min(refs$mem_alloc)
  speed_ratio <- best_time / eig$median
  memory_ratio <- best_mem / eig$mem_alloc
  eig_certified <- isTRUE(eig$certificate_passed) && eig$nconv >= k
  speed_gate <- isTRUE(speed_ratio >= speed_ratio_required)
  memory_gate <- isTRUE(memory_ratio >= memory_ratio_required)

  data.frame(
    eigencore_certified = eig_certified,
    eigencore_nconv = eig$nconv,
    requested = k,
    speed_ratio_vs_best_reference = speed_ratio,
    memory_ratio_vs_best_reference = memory_ratio,
    speed_gate = speed_gate,
    memory_gate = memory_gate,
    passed = eig_certified && speed_gate && memory_gate,
    stringsAsFactors = FALSE
  )
}

evaluate_preconditioned_lobpcg_gate <- function(rows, k,
                                                speed_ratio_required = release_speed_gate("lobpcg"),
                                                memory_ratio_required = release_memory_gate("lobpcg")) {
  lob <- rows[rows$method %in% c(
    "eigencore_lobpcg_preconditioned",
    "eigencore_lobpcg_tridiagonal"
  ), , drop = FALSE]
  scalar <- rows[rows$method == "eigencore", , drop = FALSE]
  refs <- rows[
    rows$method %in% c("RSpectra", "PRIMME") & rows$certificate_passed,
    ,
    drop = FALSE
  ]
  if (nrow(lob) < 1L || nrow(scalar) != 1L) {
    stop("gate requires exactly one scalar and at least one preconditioned LOBPCG row", call. = FALSE)
  }
  if (!nrow(refs)) {
    stop("gate requires at least one certified reference row", call. = FALSE)
  }

  lob <- lob[order(lob$median), , drop = FALSE][1L, , drop = FALSE]
  best_time <- min(refs$median)
  best_mem <- min(refs$mem_alloc)
  speed_ratio_vs_scalar <- scalar$median / lob$median
  speed_ratio_vs_best_reference <- best_time / lob$median
  memory_ratio_vs_best_reference <- best_mem / lob$mem_alloc
  lob_certified <- isTRUE(lob$certificate_passed) && lob$nconv >= k
  scalar_speed_gate <- isTRUE(speed_ratio_vs_scalar >= speed_ratio_required)
  reference_speed_gate <- isTRUE(speed_ratio_vs_best_reference >= speed_ratio_required)
  memory_gate <- isTRUE(memory_ratio_vs_best_reference >= memory_ratio_required)

  data.frame(
    lobpcg_certified = lob_certified,
    lobpcg_nconv = lob$nconv,
    requested = k,
    selected_method = lob$method,
    preconditioner_kind = lob$preconditioner_kind,
    preconditioner_native = lob$preconditioner_native,
    preconditioner_calls = lob$preconditioner_calls,
    speed_ratio_vs_scalar = speed_ratio_vs_scalar,
    speed_ratio_vs_best_reference = speed_ratio_vs_best_reference,
    memory_ratio_vs_best_reference = memory_ratio_vs_best_reference,
    scalar_speed_gate = scalar_speed_gate,
    reference_speed_gate = reference_speed_gate,
    memory_gate = memory_gate,
    passed = lob_certified && scalar_speed_gate && reference_speed_gate && memory_gate,
    stringsAsFactors = FALSE
  )
}

benchmark_native_hermitian_gate <- function(A, k, target = smallest(),
                                            iterations = 3L, tol = 1e-8,
                                            seed = 1L,
                                            methods = c("eigencore", "RSpectra", "PRIMME")) {
  rows <- benchmark_eigen_case(
    A,
    k = k,
    target = target,
    methods = methods,
    iterations = iterations,
    tol = tol,
    seed = seed
  )
  gate <- evaluate_native_hermitian_gate(rows, k = k)
  list(rows = rows, gate = gate)
}

evaluate_reference_gate <- function(rows, subject = "eigencore", references = setdiff(unique(rows$method), subject),
                                    requested, speed_ratio_required = release_speed_gate("hermitian"),
                                    memory_ratio_required = 1.0) {
  eig <- rows[rows$method == subject, , drop = FALSE]
  refs <- rows[rows$method %in% references & rows$certificate_passed, , drop = FALSE]
  if (nrow(eig) != 1L) {
    stop("gate requires exactly one subject row", call. = FALSE)
  }
  if (!nrow(refs)) {
    return(data.frame(
      subject_certified = isTRUE(eig$certificate_passed) && eig$nconv >= requested,
      subject_nconv = eig$nconv,
      requested = requested,
      speed_ratio_vs_best_reference = NA_real_,
      memory_ratio_vs_best_reference = NA_real_,
      speed_gate = FALSE,
      memory_gate = FALSE,
      passed = FALSE,
      note = "no certified reference rows",
      stringsAsFactors = FALSE
    ))
  }
  best_time <- min(refs$median)
  best_mem <- min(refs$mem_alloc)
  speed_ratio <- best_time / eig$median
  memory_ratio <- best_mem / eig$mem_alloc
  subject_certified <- isTRUE(eig$certificate_passed) && eig$nconv >= requested
  speed_gate <- isTRUE(speed_ratio >= speed_ratio_required)
  memory_gate <- isTRUE(memory_ratio >= memory_ratio_required)
  data.frame(
    subject_certified = subject_certified,
    subject_nconv = eig$nconv,
    requested = requested,
    speed_ratio_vs_best_reference = speed_ratio,
    memory_ratio_vs_best_reference = memory_ratio,
    speed_gate = speed_gate,
    memory_gate = memory_gate,
    passed = subject_certified && speed_gate && memory_gate,
    note = "",
    stringsAsFactors = FALSE
  )
}

benchmark_svd_case <- function(A, rank, methods = NULL, iterations = 3L,
                               tol = 1e-8, seed = 1L) {
  methods <- bench_methods("svd", methods)
  rows <- lapply(methods, function(method) {
    set.seed(seed)
    timed <- run_timed(eigencore:::run_svd_method(method, A, rank = rank, tol = tol, seed = seed),
                       iterations = iterations)
    cert <- certify_svd_result(A, timed$value, tol = tol)
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
      preconditioner_kind = result_preconditioner_field(timed$value, "kind"),
      preconditioner_native = result_preconditioner_field(timed$value, "native"),
      preconditioner_calls = result_preconditioner_calls(timed$value),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

save_benchmark_result <- function(result, name) {
  dir.create("inst/benchmarks/results", recursive = TRUE, showWarnings = FALSE)
  path <- file.path(
    "inst/benchmarks/results",
    paste0(format(Sys.Date(), "%Y%m%d"), "-", name, ".rds")
  )
  saveRDS(result, path)
  path
}
