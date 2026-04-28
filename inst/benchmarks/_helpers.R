library(eigencore)

`%||%` <- function(x, y) if (is.null(x)) y else x

benchmark_arg_value <- function(args, prefix) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) {
    return(NULL)
  }
  sub(paste0("^", prefix), "", hit[[1L]])
}

benchmark_arg_csv <- function(args, prefix) {
  value <- benchmark_arg_value(args, prefix)
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

benchmark_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  iterations <- benchmark_arg_value(args, "--iterations=")
  iterations <- if (is.null(iterations)) {
    NA_integer_
  } else {
    as.integer(iterations)
  }
  if (!is.na(iterations) && iterations < 1L) {
    stop("--iterations must be a positive integer.", call. = FALSE)
  }
  list(
    quick = "--quick" %in% args,
    save = "--save" %in% args,
    strict = "--strict" %in% args,
    block_candidate = "--block-candidate" %in% args,
    include_dense = "--include-dense" %in% args,
    svd_projected_stop = "--projected-stop" %in% args,
    h_candidate = "--h-candidate" %in% args,
    iterations = iterations,
    subject = benchmark_arg_value(args, "--subject="),
    methods = benchmark_arg_csv(args, "--methods="),
    cases = benchmark_arg_csv(args, "--cases=")
  )
}

svd_surface_default_methods <- function(args) {
  methods <- c(
    "eigencore",
    "eigencore_golub_kahan",
    "eigencore_randomized",
    "RSpectra",
    "PRIMME",
    "irlba",
    "rsvd",
    "base"
  )
  if (isTRUE(args$svd_projected_stop)) {
    methods <- append(methods, "eigencore_golub_kahan_projected", after = 2L)
  }
  if (isTRUE(args$h_candidate)) {
    methods <- c(
      "eigencore_golub_kahan",
      "eigencore_golub_kahan_projected",
      "RSpectra",
      "PRIMME",
      "irlba",
      "rsvd",
      "base"
    )
  }
  if (!is.null(args$methods)) {
    methods <- args$methods
  }
  methods
}

svd_surface_gate_subject <- function(args, methods) {
  subject <- args$subject %||% if (isTRUE(args$h_candidate)) {
    "eigencore_golub_kahan_projected"
  } else {
    "eigencore"
  }
  if (!subject %in% methods) {
    stop(
      "SVD gate subject `", subject, "` is not in the selected methods. ",
      "Use --methods=... to include it, --projected-stop, or --h-candidate.",
      call. = FALSE
    )
  }
  subject
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

dense_hermitian_with_spectrum <- function(values, seed = 1L) {
  n <- length(values)
  set.seed(seed)
  Q <- qr.Q(qr(matrix(stats::rnorm(n * n), nrow = n, ncol = n)))
  Q %*% (values * t(Q))
}

g1_candidate_baseline_cases <- function(quick = FALSE) {
  path_n <- if (isTRUE(quick)) 80L else 120L
  list(
    list(
      case = "path_laplacian",
      A = path_laplacian(path_n),
      k = 4L,
      target = smallest(),
      seed = 1200L + path_n
    ),
    list(
      case = "dense_hermitian",
      A = dense_hermitian_with_spectrum(c(9, 6, 4, 2, 1, 0.5, 0.25, 0.1), seed = 1208L),
      k = 3L,
      target = largest(),
      seed = 1208L
    ),
    list(
      case = "clustered",
      A = dense_hermitian_with_spectrum(c(9, 8.999, 8.998, 2, 1, 0.5, 0.25, 0.1), seed = 1210L),
      k = 3L,
      target = largest(),
      seed = 1210L
    ),
    list(
      case = "ill_conditioned_diag",
      A = diag(c(1e8, 1e4, 1, 1e-4, 1e-8, 1e-12)),
      k = 3L,
      target = largest(),
      seed = 1206L
    )
  )
}

benchmark_g1_candidate_baseline <- function(quick = FALSE,
                                            iterations = if (isTRUE(quick)) 1L else 5L,
                                            methods = c("eigencore", "eigencore_block_candidate", "RSpectra", "PRIMME"),
                                            tol = 1e-8) {
  rows <- lapply(g1_candidate_baseline_cases(quick = quick), function(case) {
    out <- benchmark_eigen_case(
      case$A,
      k = case$k,
      target = case$target,
      methods = methods,
      iterations = iterations,
      tol = tol,
      seed = case$seed
    )
    out$case <- case$case
    out$n <- nrow(case$A)
    out$k <- case$k
    out
  })
  rows <- do.call(rbind, rows)
  row.names(rows) <- NULL
  rows
}

bench_methods <- function(kind, requested = NULL) {
  methods <- if (identical(kind, "eigen")) {
    eigencore:::available_eigen_methods()
  } else {
    eigencore:::available_svd_methods()
  }
  if (!is.null(requested)) intersect(requested, methods) else methods
}

run_timed <- function(expr, iterations = 3L, seed = NULL) {
  expr <- substitute(expr)
  env <- parent.frame()
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("bench is required for release benchmarks.", call. = FALSE)
  }
  value <- NULL
  eval_once <- function() {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    eval(expr, envir = env)
  }
  mark <- bench::mark(
    value <- eval_once(),
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
  if (!is.null(fit$certificate)) {
    return(fit$certificate)
  }
  vals <- eigencore:::method_values(fit, kind = "eigen")
  vecs <- fit$vectors
  if (is.null(vals) || is.null(vecs)) return(eigencore:::empty_certificate(tol, "vectors unavailable"))
  eigencore:::certify_eigen_operator(as_operator(A), vals, vecs, tol = tol)
}

certify_svd_result <- function(A, fit, tol = 1e-8) {
  if (!is.null(fit$certificate)) {
    return(fit$certificate)
  }
  d <- eigencore:::method_values(fit, kind = "svd")
  u <- fit$u
  v <- fit$v
  if (is.null(d) || is.null(u) || is.null(v)) {
    return(eigencore:::empty_certificate(tol, "both singular-vector sides unavailable"))
  }
  eigencore:::certify_svd_operator(as_operator(A), d, u, v, tol = tol)
}

time_certified_eigen_method <- function(method, A, k, target, tol, seed,
                                        iterations) {
  set.seed(seed)
  warm_fit <- eigencore:::run_eigen_method(method, A, k = k, target = target, tol = tol)
  warm_cert <- certify_eigen_result(A, warm_fit, tol = tol)
  rm(warm_fit, warm_cert)
  gc()

  set.seed(seed)
  solver <- run_timed(
    eigencore:::run_eigen_method(method, A, k = k, target = target, tol = tol),
    iterations = iterations,
    seed = seed
  )
  fit_for_cert <- solver$value
  certificate <- run_timed(
    certify_eigen_result(A, fit_for_cert, tol = tol),
    iterations = iterations
  )
  set.seed(seed)
  total <- run_timed({
    fit <- eigencore:::run_eigen_method(method, A, k = k, target = target, tol = tol)
    cert <- certify_eigen_result(A, fit, tol = tol)
    list(fit = fit, cert = cert)
  }, iterations = iterations, seed = seed)
  list(
    fit = total$value$fit,
    cert = total$value$cert,
    solver = solver,
    certificate = certificate,
    total = total
  )
}

time_certified_svd_method <- function(method, A, rank, tol, seed, iterations) {
  set.seed(seed)
  warm_fit <- eigencore:::run_svd_method(method, A, rank = rank, tol = tol, seed = seed)
  warm_cert <- certify_svd_result(A, warm_fit, tol = tol)
  rm(warm_fit, warm_cert)
  gc()

  set.seed(seed)
  solver <- run_timed(
    eigencore:::run_svd_method(method, A, rank = rank, tol = tol, seed = seed),
    iterations = iterations,
    seed = seed
  )
  fit_for_cert <- solver$value
  certificate <- run_timed(
    certify_svd_result(A, fit_for_cert, tol = tol),
    iterations = iterations
  )
  set.seed(seed)
  total <- run_timed({
    fit <- eigencore:::run_svd_method(method, A, rank = rank, tol = tol, seed = seed)
    cert <- certify_svd_result(A, fit, tol = tol)
    list(fit = fit, cert = cert)
  }, iterations = iterations, seed = seed)
  list(
    fit = total$value$fit,
    cert = total$value$cert,
    solver = solver,
    certificate = certificate,
    total = total
  )
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

result_restart_field <- function(x, field) {
  info <- x$restart %||% NULL
  if (is.null(info) || is.null(info[[field]])) {
    return(NA)
  }
  info[[field]]
}

result_restart_logical <- function(x, field) {
  value <- result_restart_field(x, field)
  if (length(value) != 1L || is.na(value)) {
    return(NA)
  }
  isTRUE(value)
}

result_restart_numeric <- function(x, field) {
  value <- result_restart_field(x, field)
  if (length(value) != 1L || is.na(value)) {
    return(NA_real_)
  }
  as.numeric(value)
}

result_iterations <- function(x) {
  x$iterations %||% x$niter %||% NA_integer_
}

result_matvecs <- function(x) {
  x$matvecs %||% x$nops %||% NA_integer_
}

result_restarts <- function(x) {
  x$restarts %||% result_restart_field(x, "restarts_used")
}

result_ortho_passes <- function(x) {
  x$ortho_passes %||% result_restart_field(x, "ortho_passes")
}

result_locking_events <- function(x) {
  x$locking_events %||% result_restart_field(x, "locking_events")
}

result_block_size <- function(x) {
  x$block %||% result_restart_field(x, "block")
}

result_stage_seconds <- function(x, field) {
  stages <- x$stage_seconds %||% result_restart_field(x, "stage_seconds") %||% NULL
  if (is.null(stages) || !field %in% names(stages)) {
    return(NA_real_)
  }
  as.numeric(stages[[field]])
}

benchmark_failed_eigen_row <- function(method, seed, error) {
  data.frame(
    method = method,
    median = Inf,
    min = Inf,
    mem_alloc = NA_real_,
    solver_median = Inf,
    solver_min = Inf,
    solver_mem_alloc = NA_real_,
    certificate_median = Inf,
    certificate_min = Inf,
    certificate_mem_alloc = NA_real_,
    total_median = Inf,
    total_min = Inf,
    total_mem_alloc = NA_real_,
    max_residual = Inf,
    max_backward_error = Inf,
    orthogonality_loss = Inf,
    certificate_passed = FALSE,
    certificate_type = "method_error",
    norm_bound_type = NA_character_,
    scale_is_estimate = NA,
    nconv = 0L,
    iterations = NA_integer_,
    matvecs = NA_integer_,
    restarts = NA_integer_,
    ortho_passes = NA_integer_,
    locking_events = NA_integer_,
    block_size = NA_integer_,
    stage_apply_seconds = NA_real_,
    stage_recurrence_seconds = NA_real_,
    stage_reorthogonalization_seconds = NA_real_,
    stage_projected_solve_seconds = NA_real_,
    stage_projection_update_seconds = NA_real_,
    stage_projection_copy_seconds = NA_real_,
    stage_projected_eigensolve_seconds = NA_real_,
    stage_selected_vector_copy_seconds = NA_real_,
    stage_ritz_residual_seconds = NA_real_,
    stage_ritz_vector_form_seconds = NA_real_,
    stage_ritz_operator_apply_seconds = NA_real_,
    stage_ritz_norm_seconds = NA_real_,
    stage_ritz_final_polish_seconds = NA_real_,
    stage_locking_seconds = NA_real_,
    stage_restart_seconds = NA_real_,
    preconditioner_kind = NA_character_,
    preconditioner_native = NA,
    preconditioner_calls = NA_integer_,
    seed = seed,
    pkg_version = as.character(utils::packageVersion("eigencore")),
    error = conditionMessage(error),
    stringsAsFactors = FALSE
  )
}

benchmark_eigen_case <- function(A, k, target = largest(), methods = NULL,
                                 iterations = 3L, tol = 1e-8, seed = 1L) {
  methods <- bench_methods("eigen", methods)
  rows <- lapply(methods, function(method) {
    timed <- tryCatch(
      time_certified_eigen_method(method, A, k, target, tol, seed, iterations),
      error = function(e) e
    )
    if (inherits(timed, "error")) {
      return(benchmark_failed_eigen_row(method, seed, timed))
    }
    fit <- timed$fit
    cert <- timed$cert
    data.frame(
      method = method,
      median = timed$total$median,
      min = timed$total$min,
      mem_alloc = timed$total$mem_alloc,
      solver_median = timed$solver$median,
      solver_min = timed$solver$min,
      solver_mem_alloc = timed$solver$mem_alloc,
      certificate_median = timed$certificate$median,
      certificate_min = timed$certificate$min,
      certificate_mem_alloc = timed$certificate$mem_alloc,
      total_median = timed$total$median,
      total_min = timed$total$min,
      total_mem_alloc = timed$total$mem_alloc,
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      orthogonality_loss = cert$max_orthogonality_loss,
      certificate_passed = cert$passed,
      certificate_type = cert$certificate_type,
      norm_bound_type = cert$norm_bound_type,
      scale_is_estimate = cert$scale_is_estimate,
      nconv = sum(cert$converged),
      iterations = result_iterations(fit),
      matvecs = result_matvecs(fit),
      restarts = result_restarts(fit),
      ortho_passes = result_ortho_passes(fit),
      locking_events = result_locking_events(fit),
      block_size = result_block_size(fit),
      stage_apply_seconds = result_stage_seconds(fit, "apply"),
      stage_recurrence_seconds = result_stage_seconds(fit, "recurrence"),
      stage_reorthogonalization_seconds = result_stage_seconds(fit, "reorthogonalization"),
      stage_projected_solve_seconds = result_stage_seconds(fit, "projected_solve"),
      stage_projection_update_seconds = result_stage_seconds(fit, "projection_update"),
      stage_projection_copy_seconds = result_stage_seconds(fit, "projection_copy"),
      stage_projected_eigensolve_seconds = result_stage_seconds(fit, "projected_eigensolve"),
      stage_selected_vector_copy_seconds = result_stage_seconds(fit, "selected_vector_copy"),
      stage_ritz_residual_seconds = result_stage_seconds(fit, "ritz_residual"),
      stage_ritz_vector_form_seconds = result_stage_seconds(fit, "ritz_vector_form"),
      stage_ritz_operator_apply_seconds = result_stage_seconds(fit, "ritz_operator_apply"),
      stage_ritz_norm_seconds = result_stage_seconds(fit, "ritz_norm"),
      stage_ritz_final_polish_seconds = result_stage_seconds(fit, "ritz_final_polish"),
      stage_locking_seconds = result_stage_seconds(fit, "locking"),
      stage_restart_seconds = result_stage_seconds(fit, "restart"),
      preconditioner_kind = result_preconditioner_field(fit, "kind"),
      preconditioner_native = result_preconditioner_field(fit, "native"),
      preconditioner_calls = result_preconditioner_calls(fit),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      error = "",
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
    timed <- run_timed({
      if (identical(method, "eigencore")) {
        eig_partial(A, B = B, k = k, target = target, tol = tol)
      } else {
        eig <- eigencore:::dense_generalized_spd_eigen(as.matrix(A), as.matrix(B))
        idx <- eigencore:::order_indices(eig$values, target)
        idx <- idx[seq_len(k)]
        list(values = eig$values[idx], vectors = eig$vectors[, idx, drop = FALSE])
      }
    }, iterations = iterations, seed = seed)
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
      iterations = result_iterations(timed$value),
      matvecs = result_matvecs(timed$value),
      restarts = result_restarts(timed$value),
      ortho_passes = result_ortho_passes(timed$value),
      locking_events = result_locking_events(timed$value),
      block_size = result_block_size(timed$value),
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
                                           subject = "eigencore",
                                           speed_ratio_required = release_speed_gate("hermitian"),
                                           memory_ratio_required = 1.0,
                                           reference_methods = "RSpectra",
                                           parity_methods = "PRIMME",
                                           parity_ratio_required = 1.0) {
  eig <- rows[rows$method == subject, , drop = FALSE]
  speed_refs <- rows[
    rows$method %in% reference_methods & rows$certificate_passed,
    ,
    drop = FALSE
  ]
  parity_present <- any(rows$method %in% parity_methods)
  parity_refs <- rows[
    rows$method %in% parity_methods & rows$certificate_passed,
    ,
    drop = FALSE
  ]
  memory_ref_methods <- unique(c(reference_methods, parity_methods))
  memory_refs <- rows[
    rows$method %in% memory_ref_methods & rows$certificate_passed,
    ,
    drop = FALSE
  ]
  if (nrow(eig) != 1L) {
    stop("gate requires exactly one subject row", call. = FALSE)
  }

  speed_ratio <- if (nrow(speed_refs)) {
    min(speed_refs$median) / eig$median
  } else {
    NA_real_
  }
  memory_ratio <- if (nrow(memory_refs)) {
    min(memory_refs$mem_alloc) / eig$mem_alloc
  } else {
    NA_real_
  }
  parity_required <- parity_present && length(parity_methods) > 0L
  parity_ratio <- if (nrow(parity_refs)) {
    min(parity_refs$median) / eig$median
  } else {
    NA_real_
  }
  eig_certified <- isTRUE(eig$certificate_passed) && eig$nconv >= k
  speed_gate <- nrow(speed_refs) > 0L && isTRUE(speed_ratio >= speed_ratio_required)
  memory_gate <- nrow(memory_refs) > 0L && isTRUE(memory_ratio >= memory_ratio_required)
  parity_gate <- if (isTRUE(parity_required)) {
    isTRUE(parity_ratio >= parity_ratio_required)
  } else {
    TRUE
  }
  notes <- character()
  if (!nrow(speed_refs)) {
    notes <- c(notes, "no certified speed reference rows")
  }
  if (!nrow(memory_refs)) {
    notes <- c(notes, "no certified memory reference rows")
  }
  if (isTRUE(parity_required) && !nrow(parity_refs)) {
    notes <- c(notes, "parity reference present but uncertified")
  } else if (!isTRUE(parity_required)) {
    notes <- c(notes, "parity reference not present")
  }
  note <- paste(notes, collapse = "; ")

  data.frame(
    eigencore_certified = eig_certified,
    eigencore_nconv = eig$nconv,
    subject = subject,
    requested = k,
    speed_ratio_vs_best_reference = speed_ratio,
    memory_ratio_vs_best_reference = memory_ratio,
    parity_ratio_vs_best_reference = parity_ratio,
    speed_reference_methods = paste(reference_methods, collapse = ","),
    memory_reference_methods = paste(memory_ref_methods, collapse = ","),
    parity_reference_methods = paste(parity_methods, collapse = ","),
    parity_required = parity_required,
    speed_gate = speed_gate,
    memory_gate = memory_gate,
    parity_gate = parity_gate,
    passed = eig_certified && speed_gate && memory_gate && parity_gate,
    note = note,
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
                                            methods = c("eigencore", "RSpectra", "PRIMME"),
                                            subject = "eigencore",
                                            reference_methods = "RSpectra",
                                            parity_methods = "PRIMME") {
  rows <- benchmark_eigen_case(
    A,
    k = k,
    target = target,
    methods = methods,
    iterations = iterations,
    tol = tol,
    seed = seed
  )
  gate <- evaluate_native_hermitian_gate(
    rows,
    k = k,
    subject = subject,
    reference_methods = reference_methods,
    parity_methods = parity_methods
  )
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
      subject = subject,
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
    subject = subject,
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
    timed <- time_certified_svd_method(method, A, rank, tol, seed, iterations)
    fit <- timed$fit
    cert <- timed$cert
    data.frame(
      method = method,
      median = timed$total$median,
      min = timed$total$min,
      mem_alloc = timed$total$mem_alloc,
      solver_median = timed$solver$median,
      solver_min = timed$solver$min,
      solver_mem_alloc = timed$solver$mem_alloc,
      certificate_median = timed$certificate$median,
      certificate_min = timed$certificate$min,
      certificate_mem_alloc = timed$certificate$mem_alloc,
      total_median = timed$total$median,
      total_min = timed$total$min,
      total_mem_alloc = timed$total$mem_alloc,
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      orthogonality_loss = cert$max_orthogonality_loss,
      certificate_passed = cert$passed,
      certificate_type = cert$certificate_type,
      norm_bound_type = cert$norm_bound_type,
      scale_is_estimate = cert$scale_is_estimate,
      nconv = sum(cert$converged),
      iterations = result_iterations(fit),
      matvecs = result_matvecs(fit),
      restarts = result_restarts(fit),
      ortho_passes = result_ortho_passes(fit),
      locking_events = result_locking_events(fit),
      block_size = result_block_size(fit),
      stage_apply_seconds = result_stage_seconds(fit, "apply"),
      stage_recurrence_seconds = result_stage_seconds(fit, "recurrence"),
      stage_reorthogonalization_seconds = result_stage_seconds(fit, "reorthogonalization"),
      stage_projected_solve_seconds = result_stage_seconds(fit, "projected_solve"),
      stage_projection_update_seconds = result_stage_seconds(fit, "projection_update"),
      stage_projection_copy_seconds = result_stage_seconds(fit, "projection_copy"),
      stage_projected_eigensolve_seconds = result_stage_seconds(fit, "projected_eigensolve"),
      stage_selected_vector_copy_seconds = result_stage_seconds(fit, "selected_vector_copy"),
      stage_ritz_residual_seconds = result_stage_seconds(fit, "ritz_residual"),
      stage_ritz_vector_form_seconds = result_stage_seconds(fit, "ritz_vector_form"),
      stage_ritz_operator_apply_seconds = result_stage_seconds(fit, "ritz_operator_apply"),
      stage_ritz_norm_seconds = result_stage_seconds(fit, "ritz_norm"),
      stage_ritz_final_polish_seconds = result_stage_seconds(fit, "ritz_final_polish"),
      stage_locking_seconds = result_stage_seconds(fit, "locking"),
      stage_restart_seconds = result_stage_seconds(fit, "restart"),
      restart_attempts = result_restart_field(fit, "attempts"),
      final_iterations = result_restart_field(fit, "final_iterations"),
      final_matvecs = result_restart_field(fit, "final_matvecs"),
      projected_stop_requested = result_restart_logical(fit, "projected_stop_requested"),
      projected_stop_enabled = result_restart_logical(fit, "projected_stop_enabled"),
      projected_stop_disable_reason = result_restart_field(fit, "projected_stop_disable_reason"),
      projected_stop = result_restart_logical(fit, "projected_stop"),
      projected_nconv = result_restart_field(fit, "projected_nconv"),
      projected_max_residual = result_restart_numeric(fit, "projected_max_residual"),
      projected_checks = result_restart_field(fit, "projected_checks"),
      projected_seconds = result_restart_numeric(fit, "projected_seconds"),
      first_certified_prefix = result_restart_field(fit, "first_certified_prefix"),
      final_prefix_iteration_overshoot = result_restart_field(fit, "final_prefix_iteration_overshoot"),
      final_prefix_matvec_overshoot = result_restart_field(fit, "final_prefix_matvec_overshoot"),
      stage_native_iteration_seconds = result_stage_seconds(fit, "native_iteration"),
      stage_golub_kahan_ritz_seconds = result_stage_seconds(fit, "ritz"),
      stage_retry_overhead_seconds = result_stage_seconds(fit, "retry_overhead"),
      fallback_attempted = result_restart_logical(fit, "fallback_attempted"),
      fallback_used = result_restart_logical(fit, "fallback_used"),
      fallback_method = result_restart_field(fit, "fallback_method"),
      gram_max_backward_error = result_restart_numeric(fit, "gram_max_backward_error"),
      fallback_max_backward_error = result_restart_numeric(fit, "fallback_max_backward_error"),
      preconditioner_kind = result_preconditioner_field(fit, "kind"),
      preconditioner_native = result_preconditioner_field(fit, "native"),
      preconditioner_calls = result_preconditioner_calls(fit),
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
