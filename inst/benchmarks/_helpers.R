library(eigencore)

`%||%` <- function(x, y) if (is.null(x)) y else x

benchmark_safe_ratio <- function(numerator, denominator) {
  out <- rep(NA_real_, length(numerator))
  ok <- !is.na(numerator) & !is.na(denominator) & denominator != 0
  out[ok] <- numerator[ok] / denominator[ok]
  out
}

benchmark_row_sum <- function(...) {
  values <- cbind(...)
  out <- rep(NA_real_, nrow(values))
  ok <- rowSums(!is.na(values)) > 0L
  out[ok] <- rowSums(values[ok, , drop = FALSE], na.rm = TRUE)
  out
}

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

benchmark_case_field <- function(case, name) {
  case[[name, exact = TRUE]]
}

benchmark_case_id <- function(case) {
  id <- benchmark_case_field(case, "id")
  if (!is.null(id)) {
    return(id)
  }
  case_name <- benchmark_case_field(case, "case")
  m <- benchmark_case_field(case, "m")
  n <- benchmark_case_field(case, "n")
  if (!is.null(m) && !is.null(n)) {
    return(paste0(case_name, ":", m, "x", n))
  }
  if (is.null(n)) {
    case_name
  } else {
    paste0(case_name, ":", n)
  }
}

filter_benchmark_cases <- function(cases, selected = NULL) {
  if (is.null(selected) || !length(selected)) {
    return(cases)
  }
  ids <- vapply(cases, benchmark_case_id, character(1))
  names <- vapply(cases, benchmark_case_field, character(1), name = "case")
  keep <- ids %in% selected | names %in% selected
  if (!any(keep)) {
    stop(
      "--cases did not match any benchmark cases. Available cases: ",
      paste(ids, collapse = ", "),
      call. = FALSE
    )
  }
  cases[keep]
}

message_benchmark_case <- function(label, case) {
  requested <- case$k %||% case$rank %||% NA_integer_
  requested_name <- if (!is.na(requested)) {
    if (!is.null(case$rank)) " rank=" else " k="
  } else {
    " requested="
  }
  message(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    label, ": ", benchmark_case_id(case), requested_name, requested
  )
  flush.console()
}

svd_surface_default_methods <- function(args) {
  methods <- c(
    "eigencore",
    "eigencore_golub_kahan",
    "eigencore_randomized",
    "RSpectra",
    "irlba",
    "rsvd",
    "base"
  )
  if (isTRUE(args$svd_projected_stop)) {
    methods <- append(methods, "eigencore_golub_kahan_projected", after = 2L)
  }
  if (isTRUE(args$h_candidate)) {
    methods <- c(
      "eigencore",
      "eigencore_golub_kahan",
      "eigencore_golub_kahan_one_sided",
      "eigencore_irlba_lbd_one_sided",
      "eigencore_irlba_lbd_retained_native",
      "eigencore_irlba_lbd_retained_bpro",
      "eigencore_irlba_lbd_retained_bpro_one_sided_guarded",
      "eigencore_irlba_lbd_retained_bpro_block_guarded",
      "eigencore_irlba_lbd_normal_scout",
      "eigencore_golub_kahan_projected",
      "eigencore_implicit_normal_lanczos",
      "eigencore_gram_dsyevx",
      "eigencore_block_golub_kahan_cycle",
      "eigencore_block_golub_kahan_cycle_cached",
      "eigencore_block_golub_kahan_cycle_cached_random",
      "eigencore_block_golub_kahan_cycle_residual",
      "eigencore_block_golub_kahan_cycle_lean",
      "eigencore_block_golub_kahan_retained",
      "eigencore_block_golub_kahan_retained_cached",
      "eigencore_block_golub_kahan_retained_deflated",
      "RSpectra",
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
    "eigencore"
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

svd_internal_methods <- function() {
  c(
    "eigencore",
    "eigencore_smallest",
    "eigencore_interior",
    "eigencore_golub_kahan",
    "eigencore_golub_kahan_smallest",
    "eigencore_golub_kahan_interior",
    "eigencore_golub_kahan_one_sided",
    "eigencore_irlba_lbd_one_sided",
    "eigencore_irlba_lbd_retained_native",
    "eigencore_irlba_lbd_retained_bpro",
    "eigencore_irlba_lbd_retained_bpro_one_sided_guarded",
    "eigencore_irlba_lbd_retained_bpro_block_guarded",
    "eigencore_irlba_lbd_normal_scout",
    "eigencore_golub_kahan_projected",
    "eigencore_implicit_normal_lanczos",
    "eigencore_gram_dsyevx",
    "eigencore_block_golub_kahan_cycle",
    "eigencore_block_golub_kahan_cycle_cached",
    "eigencore_block_golub_kahan_cycle_cached_random",
    "eigencore_block_golub_kahan_cycle_residual",
    "eigencore_block_golub_kahan_cycle_lean",
    "eigencore_block_golub_kahan_retained",
    "eigencore_block_golub_kahan_retained_cached",
    "eigencore_block_golub_kahan_retained_deflated",
    "eigencore_randomized"
  )
}

svd_crossover_evidence_contract <- function(rows, summary,
                                            subject = "eigencore",
                                            orientations = c("tall", "wide")) {
  required_row_fields <- c(
    "case", "method", "rank", "certificate_passed", "nconv"
  )
  required_summary_fields <- c(
    "case", "orientation", "long_side", "best_certified_reference",
    "raw_speed_ratio", "certified_total_speed_ratio"
  )
  missing_fields <- c(
    setdiff(required_row_fields, names(rows)),
    setdiff(required_summary_fields, names(summary))
  )
  if (length(missing_fields)) {
    stop(
      "crossover evidence is missing required fields: ",
      paste(unique(missing_fields), collapse = ", "),
      call. = FALSE
    )
  }

  subject_cases <- unique(rows$case[rows$method == subject])
  external <- rows[
    rows$method != subject & !rows$method %in% svd_internal_methods(),
    ,
    drop = FALSE
  ]
  external_certified <- nrow(external) > 0L && all(
    !is.na(external$certificate_passed) & external$certificate_passed &
      !is.na(external$nconv) & external$nconv >= external$rank
  )
  summary_complete <- length(subject_cases) > 0L &&
    nrow(summary) == length(subject_cases) &&
    identical(sort(as.character(summary$case)), sort(as.character(subject_cases))) &&
    all(!is.na(summary$best_certified_reference)) &&
    all(nzchar(summary$best_certified_reference)) &&
    all(is.finite(summary$raw_speed_ratio) & summary$raw_speed_ratio > 0) &&
    all(
      is.finite(summary$certified_total_speed_ratio) &
        summary$certified_total_speed_ratio > 0
    )
  orientations_complete <- all(orientations %in% unique(summary$orientation))

  sampled_crossover <- function(field) {
    orientations_complete && all(vapply(orientations, function(orientation) {
      slice <- summary[summary$orientation == orientation, , drop = FALSE]
      slice <- slice[order(slice$long_side), , drop = FALSE]
      if (nrow(slice) < 2L) {
        return(FALSE)
      }
      ratios <- slice[[field]]
      faster <- which(ratios > 1)
      any(vapply(faster, function(index) {
        index > 1L && any(ratios[seq_len(index - 1L)] <= 1)
      }, logical(1)))
    }, logical(1)))
  }

  raw_crossover <- summary_complete && sampled_crossover("raw_speed_ratio")
  certified_crossover <- summary_complete &&
    sampled_crossover("certified_total_speed_ratio")
  data.frame(
    subject = subject,
    external_reference_rows = nrow(external),
    external_references_certified = external_certified,
    summary_complete = summary_complete,
    orientations_complete = orientations_complete,
    raw_sampled_crossover = raw_crossover,
    certified_total_sampled_crossover = certified_crossover,
    passed = external_certified && summary_complete &&
      orientations_complete && certified_crossover,
    stringsAsFactors = FALSE
  )
}

release_speed_gate <- function(kind) {
  switch(
    kind,
    hermitian = 1.25,
    svd = 1.1,
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

nearly_low_rank_svd_matrix <- function(m, n, rank, noise = 1e-3, seed = 1L) {
  set.seed(seed)
  U <- qr.Q(qr(matrix(stats::rnorm(m * rank), nrow = m, ncol = rank)))
  V <- qr.Q(qr(matrix(stats::rnorm(n * rank), nrow = n, ncol = rank)))
  signal <- U %*% (seq(rank, 1) * t(V))
  signal + noise * matrix(stats::rnorm(m * n), nrow = m, ncol = n)
}

randomized_rsvd_benchmark_cases <- function(quick = FALSE) {
  if (isTRUE(quick)) {
    return(list(
      list(
        case = "exact_low_rank_dense",
        id = "exact_low_rank_dense:120x80",
        A = nearly_low_rank_svd_matrix(120L, 80L, rank = 8L, noise = 0, seed = 1601L),
        rank = 6L,
        seed = 1601L,
        release_gate_required = FALSE,
        release_gate_note = "small dense diagnostic; fixed public-result overhead dominates this size"
      ),
      list(
        case = "slow_decay_dense",
        id = "slow_decay_dense:140x90",
        A = slow_decay_svd_matrix(140L, 90L, decay = 0.25, seed = 1602L),
        rank = 8L,
        seed = 1602L,
        release_gate_required = FALSE,
        release_gate_note = "diagnostic row; rsvd baseline fails eigencore certification"
      ),
      list(
        case = "low_rank_sparse",
        id = "low_rank_sparse:140x90",
        A = {
          set.seed(1603L)
          Matrix::rsparsematrix(140L, 12L, density = 0.12) %*%
            Matrix::rsparsematrix(12L, 90L, density = 0.12)
        },
        rank = 8L,
        seed = 1603L,
        release_gate_required = FALSE,
        release_gate_note = "sparse native-path diagnostic; quick timing is fixed-overhead sensitive"
      )
    ))
  }
  list(
    list(
      case = "exact_low_rank_dense",
      id = "exact_low_rank_dense:2000x500",
      A = nearly_low_rank_svd_matrix(2000L, 500L, rank = 60L, noise = 0, seed = 1701L),
      rank = 50L,
      seed = 1701L,
      release_gate_required = TRUE,
      release_gate_note = "large exact-low-rank randomized planner candidate"
    ),
    list(
      case = "nearly_low_rank_dense",
      id = "nearly_low_rank_dense:2000x500",
      A = nearly_low_rank_svd_matrix(2000L, 500L, rank = 80L, noise = 1e-3, seed = 1702L),
      rank = 50L,
      seed = 1702L,
      release_gate_required = FALSE,
      release_gate_note = "diagnostic row; rsvd baseline fails eigencore certification"
    ),
    list(
      case = "slow_decay_dense",
      id = "slow_decay_dense:2000x500",
      A = slow_decay_svd_matrix(2000L, 500L, decay = 0.25, seed = 1703L),
      rank = 50L,
      seed = 1703L,
      release_gate_required = FALSE,
      release_gate_note = "diagnostic row; rsvd baseline fails eigencore certification"
    ),
    list(
      case = "low_rank_sparse",
      id = "low_rank_sparse:2000x500",
      A = {
        set.seed(1705L)
        Matrix::rsparsematrix(2000L, 40L, density = 0.05) %*%
          Matrix::rsparsematrix(40L, 500L, density = 0.05)
      },
      rank = 30L,
      seed = 1705L,
      release_gate_required = FALSE,
      release_gate_note = "native sparse CSC controller contract; sparse baselines are diagnostic rather than promoted release gates"
    )
  )
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
    benchmark_available_eigen_methods()
  } else {
    benchmark_available_svd_methods()
  }
  if (is.null(requested)) {
    return(methods)
  }
  missing <- setdiff(requested, methods)
  if (length(missing)) {
    label <- if (identical(kind, "eigen")) "eigen" else "SVD"
    stop(
      "Requested ", label, " benchmark method(s) unavailable in the loaded ",
      "eigencore namespace: ", paste(missing, collapse = ", "), ". ",
      "Available methods: ", paste(methods, collapse = ", "), ". ",
      "If this is a source checkout benchmark, reinstall eigencore or run the ",
      "script under pkgload::load_all().",
      call. = FALSE
    )
  }
  requested
}

benchmark_available_eigen_methods <- function() {
  c(
    eigencore:::available_eigen_methods(),
    if (requireNamespace("PRIMME", quietly = TRUE)) "PRIMME"
  )
}

benchmark_available_svd_methods <- function() {
  c(
    eigencore:::available_svd_methods(),
    if (requireNamespace("PRIMME", quietly = TRUE)) "PRIMME"
  )
}

run_benchmark_eigen_method <- function(method, A, k, target, tol) {
  if (identical(method, "PRIMME")) {
    which <- eigencore:::target_to_rspectra_which(target, symmetric = TRUE)
    return(PRIMME::eigs_sym(A, NEig = k, which = which, tol = tol))
  }
  eigencore:::run_eigen_method(method, A, k = k, target = target, tol = tol)
}

run_benchmark_svd_method <- function(method, A, rank, tol, seed = NULL) {
  if (identical(method, "PRIMME")) {
    return(PRIMME::svds(A, NSvals = rank, tol = tol))
  }
  eigencore:::run_svd_method(method, A, rank = rank, tol = tol, seed = seed)
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
    # bench::mark(memory = TRUE) calls utils::Rprofmem(), which errors on R
    # builds without --enable-memory-profiling (e.g. r-devel fedora).
    memory = isTRUE(capabilities("profmem")),
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

svd_certificate_source <- function(fit) {
  if (is.null(fit$certificate)) {
    "eigencore_recomputed"
  } else {
    "method_certificate"
  }
}

certificate_residual_max <- function(cert, side) {
  residuals <- cert$residuals %||% list()
  value <- residuals[[side]]
  if (is.null(value) || !length(value)) {
    return(NA_real_)
  }
  max(value, na.rm = TRUE)
}

certificate_orthogonality_component <- function(cert, side) {
  orth <- cert$orthogonality %||% numeric()
  if (!length(orth)) {
    return(NA_real_)
  }
  if (!is.null(names(orth)) && side %in% names(orth)) {
    return(as.numeric(orth[[side]]))
  }
  idx <- match(side, c("U", "V"))
  if (!is.na(idx) && length(orth) >= idx) {
    return(as.numeric(orth[[idx]]))
  }
  NA_real_
}

svd_result_sorted_descending <- function(fit) {
  d <- eigencore:::method_values(fit, kind = "svd")
  if (is.null(d) || length(d) < 2L) {
    return(TRUE)
  }
  all(diff(as.numeric(d)) <= sqrt(.Machine$double.eps) * max(1, abs(d), na.rm = TRUE))
}

svd_subspace_error <- function(observed, oracle) {
  if (is.null(observed) || is.null(oracle) || !ncol(observed) || !ncol(oracle)) {
    return(NA_real_)
  }
  observed <- qr.Q(qr(as.matrix(observed)))
  oracle <- qr.Q(qr(as.matrix(oracle)))
  sqrt(sum((tcrossprod(observed) - tcrossprod(oracle))^2))
}

svd_oracle_accuracy <- function(A, fit, rank, oracle = NULL) {
  if (is.null(oracle)) {
    oracle <- svd(as.matrix(A), nu = rank, nv = rank)
  }
  idx <- seq_len(min(rank, length(oracle$d)))
  observed <- eigencore:::method_values(fit, kind = "svd")
  observed <- observed[seq_len(min(length(observed), length(idx)))]
  expected <- oracle$d[idx][seq_along(observed)]
  singular_value_relative_error <- if (length(observed)) {
    max(abs(observed - expected) / pmax(abs(expected), .Machine$double.eps))
  } else {
    Inf
  }
  u <- fit$u
  v <- fit$v
  list(
    singular_value_relative_error = singular_value_relative_error,
    left_subspace_error = svd_subspace_error(u, oracle$u[, idx, drop = FALSE]),
    right_subspace_error = svd_subspace_error(v, oracle$v[, idx, drop = FALSE])
  )
}

time_certified_eigen_method <- function(method, A, k, target, tol, seed,
                                        iterations) {
  set.seed(seed)
  warm_fit <- run_benchmark_eigen_method(method, A, k = k, target = target, tol = tol)
  warm_cert <- certify_eigen_result(A, warm_fit, tol = tol)
  rm(warm_fit, warm_cert)
  gc()

  set.seed(seed)
  solver <- run_timed(
    run_benchmark_eigen_method(method, A, k = k, target = target, tol = tol),
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
    fit <- run_benchmark_eigen_method(method, A, k = k, target = target, tol = tol)
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
  warm_fit <- run_benchmark_svd_method(method, A, rank = rank, tol = tol, seed = seed)
  warm_cert <- certify_svd_result(A, warm_fit, tol = tol)
  rm(warm_fit, warm_cert)
  gc()

  set.seed(seed)
  solver <- run_timed(
    run_benchmark_svd_method(method, A, rank = rank, tol = tol, seed = seed),
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
    fit <- run_benchmark_svd_method(method, A, rank = rank, tol = tol, seed = seed)
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

result_restart_character <- function(x, field) {
  value <- result_restart_field(x, field)
  if (!length(value) || all(is.na(value))) {
    return(NA_character_)
  }
  paste(as.character(value), collapse = ",")
}

result_restart_integer <- function(x, field) {
  value <- result_restart_field(x, field)
  if (length(value) != 1L || is.na(value)) {
    return(NA_integer_)
  }
  as.integer(value)
}

result_attempt_history <- function(x) {
  history <- result_restart_field(x, "attempt_history")
  if (is.data.frame(history)) {
    history
  } else {
    NULL
  }
}

result_attempt_history_max <- function(x, field) {
  history <- result_attempt_history(x)
  if (is.null(history) || !field %in% names(history) || !nrow(history)) {
    return(NA_integer_)
  }
  value <- max(history[[field]], na.rm = TRUE)
  if (!is.finite(value)) NA_integer_ else as.integer(value)
}

result_attempt_history_count_true <- function(x, field) {
  history <- result_attempt_history(x)
  if (is.null(history) || !field %in% names(history) || !nrow(history)) {
    return(NA_integer_)
  }
  value <- history[[field]]
  if (!is.logical(value)) {
    return(NA_integer_)
  }
  as.integer(sum(value, na.rm = TRUE))
}

result_certified_attempt <- function(x) {
  history <- result_attempt_history(x)
  if (is.null(history) || !"certificate_passed" %in% names(history) || !nrow(history)) {
    return(NA_integer_)
  }
  if (!is.logical(history$certificate_passed)) {
    return(NA_integer_)
  }
  passed <- which(history$certificate_passed)
  if (length(passed)) as.integer(history$attempt[[passed[[1L]]]]) else NA_integer_
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
    native = NA,
    native_kernels = NA,
    generalized = NA,
    orthogonalization_native = NA,
    orthogonalization_methods = NA_character_,
    q_rank_final = NA_integer_,
    constrained = NA,
    constraints_rank = NA_integer_,
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
      native = result_restart_logical(fit, "native"),
      native_kernels = result_restart_logical(fit, "native_kernels"),
      generalized = result_restart_logical(fit, "generalized"),
      orthogonalization_native = result_restart_logical(fit, "orthogonalization_native"),
      orthogonalization_methods = result_restart_character(fit, "orthogonalization_methods"),
      q_rank_final = result_restart_integer(fit, "q_rank_final"),
      constrained = result_restart_logical(fit, "constrained"),
      constraints_rank = result_restart_integer(fit, "constraints_rank"),
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

evaluate_memory_diagnostics <- function(rows, subject = "eigencore",
                                        references = setdiff(unique(rows$method), subject),
                                        requested) {
  eig <- rows[rows$method == subject, , drop = FALSE]
  refs <- rows[rows$method %in% references & rows$certificate_passed, , drop = FALSE]
  if (nrow(eig) != 1L) {
    stop("memory diagnostics require exactly one subject row", call. = FALSE)
  }
  subject_certified <- isTRUE(eig$certificate_passed) && eig$nconv >= requested
  if (!nrow(refs)) {
    return(data.frame(
      subject = subject,
      subject_certified = subject_certified,
      subject_nconv = eig$nconv,
      requested = requested,
      subject_total_mem_alloc = eig$mem_alloc,
      subject_solver_mem_alloc = eig$solver_mem_alloc,
      subject_certificate_mem_alloc = eig$certificate_mem_alloc,
      subject_solver_memory_fraction = eig$solver_mem_alloc / eig$mem_alloc,
      subject_certificate_memory_fraction = eig$certificate_mem_alloc / eig$mem_alloc,
      best_reference = NA_character_,
      best_reference_total_mem_alloc = NA_real_,
      total_memory_gap_bytes = NA_real_,
      solver_memory_gap_bytes = NA_real_,
      certificate_memory_gap_bytes = NA_real_,
      total_memory_ratio_vs_best_reference = NA_real_,
      solver_memory_ratio_vs_best_reference = NA_real_,
      certificate_memory_ratio_vs_best_reference = NA_real_,
      note = "no certified reference rows",
      stringsAsFactors = FALSE
    ))
  }

  ref <- refs[order(refs$mem_alloc), , drop = FALSE][1L, , drop = FALSE]
  data.frame(
    subject = subject,
    subject_certified = subject_certified,
    subject_nconv = eig$nconv,
    requested = requested,
    subject_total_mem_alloc = eig$mem_alloc,
    subject_solver_mem_alloc = eig$solver_mem_alloc,
    subject_certificate_mem_alloc = eig$certificate_mem_alloc,
    subject_solver_memory_fraction = eig$solver_mem_alloc / eig$mem_alloc,
    subject_certificate_memory_fraction = eig$certificate_mem_alloc / eig$mem_alloc,
    best_reference = ref$method,
    best_reference_total_mem_alloc = ref$mem_alloc,
    total_memory_gap_bytes = eig$mem_alloc - ref$mem_alloc,
    solver_memory_gap_bytes = eig$solver_mem_alloc - ref$solver_mem_alloc,
    certificate_memory_gap_bytes = eig$certificate_mem_alloc - ref$certificate_mem_alloc,
    total_memory_ratio_vs_best_reference = ref$mem_alloc / eig$mem_alloc,
    solver_memory_ratio_vs_best_reference = ref$solver_mem_alloc / eig$solver_mem_alloc,
    certificate_memory_ratio_vs_best_reference = if (eig$certificate_mem_alloc > 0) {
      ref$certificate_mem_alloc / eig$certificate_mem_alloc
    } else {
      NA_real_
    },
    note = "",
    stringsAsFactors = FALSE
  )
}

benchmark_randomized_rsvd_case <- function(A, rank, methods = c("eigencore_randomized", "rsvd"),
                                           iterations = 3L, tol = 1e-8,
                                           seed = 1L) {
  methods <- bench_methods("svd", methods)
  oracle <- svd(as.matrix(A), nu = rank, nv = rank)
  rows <- lapply(methods, function(method) {
    timed <- time_certified_svd_method(method, A, rank, tol, seed, iterations)
    fit <- timed$fit
    cert <- timed$cert
    accuracy <- svd_oracle_accuracy(A, fit, rank = rank, oracle = oracle)
    data.frame(
      method = method,
      median = timed$total$median,
      min = timed$total$min,
      mem_alloc = timed$total$mem_alloc,
      solver_median = timed$solver$median,
      solver_mem_alloc = timed$solver$mem_alloc,
      certificate_median = timed$certificate$median,
      certificate_mem_alloc = timed$certificate$mem_alloc,
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      orthogonality_loss = cert$max_orthogonality_loss,
      certificate_passed = cert$passed,
      certificate_type = cert$certificate_type,
      norm_bound_type = cert$norm_bound_type,
      scale_is_estimate = cert$scale_is_estimate,
      nconv = sum(cert$converged),
      stage_random_seconds = result_stage_seconds(fit, "random"),
      stage_apply_seconds = result_stage_seconds(fit, "apply"),
      stage_normalize_seconds = result_stage_seconds(fit, "normalize"),
      stage_small_svd_seconds = result_stage_seconds(fit, "small_svd"),
      stage_vector_form_seconds = result_stage_seconds(fit, "vector_form"),
      stage_internal_certificate_seconds = result_stage_seconds(fit, "certificate"),
      stage_refinement_seconds = result_stage_seconds(fit, "refinement"),
      randomized_native_sketch = result_restart_logical(fit, "native_sketch"),
      randomized_sketch_kind = result_restart_character(fit, "sketch_kind"),
      randomized_controller_native = result_restart_logical(fit, "controller_native"),
      randomized_controller_kind = result_restart_character(fit, "controller_kind"),
      randomized_dense_native_controller = result_restart_logical(fit, "dense_native_controller"),
      randomized_sparse_native_controller = result_restart_logical(fit, "sparse_native_controller"),
      randomized_native_certificate_diagnostics =
        result_restart_logical(fit, "native_certificate_diagnostics"),
      randomized_adaptive_stop_used = result_restart_logical(fit, "adaptive_stop_used"),
      randomized_core_solver = result_restart_character(fit, "core_solver"),
      randomized_projection_kind = result_restart_character(fit, "projection_kind"),
      randomized_projection_transposed = result_restart_logical(fit, "projection_transposed"),
      singular_value_relative_error = accuracy$singular_value_relative_error,
      left_subspace_error = accuracy$left_subspace_error,
      right_subspace_error = accuracy$right_subspace_error,
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

evaluate_randomized_rsvd_gate <- function(rows, subject = "eigencore_randomized",
                                          baseline = "rsvd", requested,
                                          speed_ratio_required = release_speed_gate("randomized_svd"),
                                          accuracy_multiplier = 1.05,
                                          accuracy_floor = 1e-12) {
  eig <- rows[rows$method == subject, , drop = FALSE]
  ref <- rows[rows$method == baseline, , drop = FALSE]
  if (nrow(eig) != 1L) {
    stop("randomized rsvd gate requires exactly one subject row", call. = FALSE)
  }
  if (nrow(ref) != 1L) {
    return(data.frame(
      subject = subject,
      baseline = baseline,
      subject_certified = isTRUE(eig$certificate_passed) && eig$nconv >= requested,
      baseline_certified = FALSE,
      subject_nconv = eig$nconv,
      baseline_nconv = NA_integer_,
      requested = requested,
      speed_ratio_vs_rsvd = NA_real_,
      singular_value_error_ratio_vs_rsvd = NA_real_,
      left_subspace_error_ratio_vs_rsvd = NA_real_,
      right_subspace_error_ratio_vs_rsvd = NA_real_,
      speed_gate = FALSE,
      accuracy_gate = FALSE,
      passed = FALSE,
      note = "rsvd baseline row unavailable",
      stringsAsFactors = FALSE
    ))
  }

  sv_ratio <- eig$singular_value_relative_error /
    max(ref$singular_value_relative_error, accuracy_floor)
  left_ratio <- eig$left_subspace_error / max(ref$left_subspace_error, accuracy_floor)
  right_ratio <- eig$right_subspace_error / max(ref$right_subspace_error, accuracy_floor)
  speed_ratio <- ref$median / eig$median
  subject_certified <- isTRUE(eig$certificate_passed) && eig$nconv >= requested
  baseline_certified <- isTRUE(ref$certificate_passed) && ref$nconv >= requested
  accuracy_gate <- isTRUE(sv_ratio <= accuracy_multiplier) &&
    (is.na(left_ratio) || isTRUE(left_ratio <= accuracy_multiplier)) &&
    (is.na(right_ratio) || isTRUE(right_ratio <= accuracy_multiplier))
  speed_gate <- baseline_certified && isTRUE(speed_ratio >= speed_ratio_required)
  data.frame(
    subject = subject,
    baseline = baseline,
    subject_certified = subject_certified,
    baseline_certified = baseline_certified,
    subject_nconv = eig$nconv,
    baseline_nconv = ref$nconv,
    requested = requested,
    speed_ratio_vs_rsvd = speed_ratio,
    singular_value_error_ratio_vs_rsvd = sv_ratio,
    left_subspace_error_ratio_vs_rsvd = left_ratio,
    right_subspace_error_ratio_vs_rsvd = right_ratio,
    speed_gate = speed_gate,
    accuracy_gate = accuracy_gate,
    passed = subject_certified && baseline_certified && speed_gate && accuracy_gate,
    note = if (!baseline_certified) {
      "rsvd baseline did not satisfy eigencore certificate"
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
}

randomized_controller_contract <- function(rows) {
  eig <- rows[rows$method == "eigencore_randomized", , drop = FALSE]
  if (!nrow(eig)) {
    return(data.frame())
  }
  out <- lapply(seq_len(nrow(eig)), function(idx) {
    row <- eig[idx, , drop = FALSE]
    case_name <- if ("case" %in% names(row)) as.character(row$case) else ""
    is_sparse_case <- grepl("sparse", case_name, fixed = TRUE)
    expected_kind <- if (is_sparse_case) {
      "native_csc_randomized_controller"
    } else {
      "native_dense_randomized_controller"
    }
    provenance_gate <- isTRUE(row$randomized_controller_native) &&
      isTRUE(row$randomized_native_certificate_diagnostics) &&
      identical(as.character(row$randomized_controller_kind), expected_kind) &&
      if (is_sparse_case) {
        isTRUE(row$randomized_sparse_native_controller) &&
          !isTRUE(row$randomized_dense_native_controller)
      } else {
        isTRUE(row$randomized_dense_native_controller) &&
          !isTRUE(row$randomized_sparse_native_controller)
      }
    certificate_gate <- isTRUE(row$certificate_passed) &&
      isTRUE(row$nconv >= row$rank) &&
      !isTRUE(row$scale_is_estimate)
    data.frame(
      case = case_name,
      rank = row$rank,
      controller_kind = as.character(row$randomized_controller_kind),
      native_certificate_diagnostics =
        isTRUE(row$randomized_native_certificate_diagnostics),
      certificate_gate = certificate_gate,
      provenance_gate = provenance_gate,
      passed = certificate_gate && provenance_gate,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

benchmark_svd_case <- function(A, rank, methods = NULL, iterations = 3L,
                               tol = 1e-8, seed = 1L) {
  methods <- bench_methods("svd", methods)
  rows <- lapply(methods, function(method) {
    timed <- time_certified_svd_method(method, A, rank, tol, seed, iterations)
    fit <- timed$fit
    cert <- timed$cert
    cert_source <- svd_certificate_source(fit)
    data.frame(
      method = method,
      solver_label = fit$method %||% NA_character_,
      median = timed$total$median,
      min = timed$total$min,
      mem_alloc = timed$total$mem_alloc,
      raw_solver_median = timed$solver$median,
      raw_solver_min = timed$solver$min,
      raw_solver_mem_alloc = timed$solver$mem_alloc,
      eigencore_certificate_median = timed$certificate$median,
      eigencore_certificate_min = timed$certificate$min,
      eigencore_certificate_mem_alloc = timed$certificate$mem_alloc,
      eigencore_certified_total_median = timed$total$median,
      eigencore_certified_total_min = timed$total$min,
      eigencore_certified_total_mem_alloc = timed$total$mem_alloc,
      certificate_source = cert_source,
      certificate_recomputed_by_eigencore =
        identical(cert_source, "eigencore_recomputed"),
      solver_median = timed$solver$median,
      solver_min = timed$solver$min,
      solver_mem_alloc = timed$solver$mem_alloc,
      certificate_median = timed$certificate$median,
      certificate_min = timed$certificate$min,
      certificate_mem_alloc = timed$certificate$mem_alloc,
      total_median = timed$total$median,
      total_min = timed$total$min,
      total_mem_alloc = timed$total$mem_alloc,
      max_left_residual = certificate_residual_max(cert, "left"),
      max_right_residual = certificate_residual_max(cert, "right"),
      max_cyclic_residual = certificate_residual_max(cert, "combined"),
      max_residual = cert$max_residual,
      max_backward_error = cert$max_backward_error,
      orthogonality_loss = cert$max_orthogonality_loss,
      orthogonality_U = certificate_orthogonality_component(cert, "U"),
      orthogonality_V = certificate_orthogonality_component(cert, "V"),
      singular_values_sorted = svd_result_sorted_descending(fit),
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
      stage_gram_seconds = result_stage_seconds(fit, "gram"),
      stage_gram_eigensolve_seconds = result_stage_seconds(fit, "eigensolve"),
      stage_gram_vector_form_seconds = result_stage_seconds(fit, "vector_form"),
      stage_gram_diagnostics_seconds = result_stage_seconds(fit, "diagnostics"),
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
      native_workspace_bytes = result_restart_numeric(fit, "native_workspace_bytes"),
      native_workspace_allocator = result_restart_character(fit, "native_workspace_allocator"),
      basis_returned = result_restart_logical(fit, "basis_returned"),
      retained_restart = result_restart_logical(fit, "retained_restart"),
      retained_restart_native = result_restart_logical(fit, "retained_restart_native"),
      retained_av_cache = result_restart_logical(fit, "retained_av_cache"),
      retained_deflation = result_restart_logical(fit, "retained_deflation"),
      retained_locked_count = result_restart_integer(fit, "retained_locked_count"),
      native_attempt_certification = result_restart_logical(fit, "native_attempt_certification"),
      native_early_stop = result_restart_logical(fit, "native_early_stop"),
      retained_converged_count = result_attempt_history_max(fit, "converged_count"),
      retained_leading_converged_count = result_attempt_history_max(fit, "leading_converged_count"),
      reorthogonalization_passes = result_restart_field(fit, "reorthogonalization_passes"),
      first_certified_prefix = result_restart_field(fit, "first_certified_prefix"),
      final_prefix_iteration_overshoot = result_restart_field(fit, "final_prefix_iteration_overshoot"),
      final_prefix_matvec_overshoot = result_restart_field(fit, "final_prefix_matvec_overshoot"),
      stage_native_iteration_seconds = result_stage_seconds(fit, "native_iteration"),
      stage_golub_kahan_ritz_seconds = result_stage_seconds(fit, "ritz"),
      stage_retry_overhead_seconds = result_stage_seconds(fit, "retry_overhead"),
      attempted_subspaces = result_restart_character(fit, "attempted_subspaces"),
      max_attempted_subspace = result_attempt_history_max(fit, "max_subspace"),
      max_start_cols = result_attempt_history_max(fit, "start_cols"),
      warm_started_attempts = result_attempt_history_count_true(fit, "warm_started"),
      cached_start_attempts = result_attempt_history_count_true(fit, "cached_start_used"),
      irlba_lbd_small_work_matvecs =
        result_restart_integer(fit, "irlba_lbd_small_work_matvecs"),
      irlba_lbd_fallback_matvecs =
        result_restart_integer(fit, "irlba_lbd_fallback_matvecs"),
      irlba_lbd_scout_matvec_overhead_fraction =
        result_restart_numeric(fit, "irlba_lbd_scout_matvec_overhead_fraction"),
      irlba_lbd_small_work_accounted_seconds =
        result_restart_numeric(fit, "irlba_lbd_small_work_accounted_seconds"),
      irlba_lbd_fallback_accounted_seconds =
        result_restart_numeric(fit, "irlba_lbd_fallback_accounted_seconds"),
      irlba_lbd_retained_native_attempted =
        result_restart_logical(fit, "irlba_lbd_retained_native_attempted"),
      irlba_lbd_retained_matvecs =
        result_restart_integer(fit, "irlba_lbd_retained_matvecs"),
      irlba_lbd_total_matvecs =
        result_restart_integer(fit, "irlba_lbd_total_matvecs"),
      irlba_lbd_retained_native_fallback_reason =
        result_restart_character(fit, "irlba_lbd_retained_native_fallback_reason"),
      irlba_lbd_locking_policy =
        result_restart_character(fit, "irlba_lbd_locking_policy"),
      irlba_lbd_lock_source =
        result_restart_character(fit, "irlba_lbd_lock_source"),
      irlba_lbd_soft_locked_count =
        result_restart_integer(fit, "irlba_lbd_soft_locked_count"),
      irlba_lbd_hard_locked_count =
        result_restart_integer(fit, "irlba_lbd_hard_locked_count"),
      irlba_lbd_locked_triplets_certified =
        result_restart_logical(fit, "irlba_lbd_locked_triplets_certified"),
      irlba_lbd_locked_orthogonality_loss =
        result_restart_numeric(fit, "irlba_lbd_locked_orthogonality_loss"),
      irlba_lbd_future_vectors_orthogonal_to_locks =
        result_restart_logical(fit, "irlba_lbd_future_vectors_orthogonal_to_locks"),
      irlba_lbd_lock_fallback_reason =
        result_restart_character(fit, "irlba_lbd_lock_fallback_reason"),
      irlba_lbd_restart_state_kind =
        result_restart_character(fit, "irlba_lbd_restart_state_kind"),
      irlba_lbd_recurrence_available =
        result_restart_logical(fit, "irlba_lbd_recurrence_available"),
      irlba_lbd_augmented_recurrence =
        result_restart_logical(fit, "irlba_lbd_augmented_recurrence"),
      irlba_lbd_residual_augmented_cols =
        result_restart_integer(fit, "irlba_lbd_residual_augmented_cols"),
      irlba_lbd_augmented_tail_steps =
        result_restart_integer(fit, "irlba_lbd_augmented_tail_steps"),
      irlba_lbd_augmented_basis_cols =
        result_restart_integer(fit, "irlba_lbd_augmented_basis_cols"),
      irlba_lbd_augmented_restart_cycles =
        result_restart_integer(fit, "irlba_lbd_augmented_restart_cycles"),
      irlba_lbd_augmented_kept_vectors =
        result_restart_integer(fit, "irlba_lbd_augmented_kept_vectors"),
      irlba_lbd_augmented_small_svds =
        result_restart_integer(fit, "irlba_lbd_augmented_small_svds"),
      irlba_lbd_augmented_cached_aq_cols =
        result_restart_integer(fit, "irlba_lbd_augmented_cached_aq_cols"),
      irlba_lbd_augmented_from_scratch_matvecs =
        result_restart_integer(fit, "irlba_lbd_augmented_from_scratch_matvecs"),
      irlba_lbd_augmented_matvec_savings =
        result_restart_integer(fit, "irlba_lbd_augmented_matvec_savings"),
      irlba_lbd_augmented_min_cheap_residual =
        result_restart_numeric(fit, "irlba_lbd_augmented_min_cheap_residual"),
      irlba_lbd_augmented_final_cheap_residual =
        result_restart_numeric(fit, "irlba_lbd_augmented_final_cheap_residual"),
      irlba_lbd_augmented_reduces_from_scratch_work =
        result_restart_logical(fit, "irlba_lbd_augmented_reduces_from_scratch_work"),
      irlba_lbd_native_certificate_diagnostics_reused =
        result_restart_logical(fit, "irlba_lbd_native_certificate_diagnostics_reused"),
      irlba_lbd_native_certificate_diagnostics_swapped =
        result_restart_logical(fit, "irlba_lbd_native_certificate_diagnostics_swapped"),
      irlba_lbd_bpro_policy =
        result_restart_logical(fit, "irlba_lbd_bpro_policy"),
      irlba_lbd_bpro_passes_per_append =
        result_restart_integer(fit, "irlba_lbd_bpro_passes_per_append"),
      irlba_lbd_bpro_monitoring_threshold =
        result_restart_numeric(fit, "irlba_lbd_bpro_monitoring_threshold"),
      irlba_lbd_bpro_monitored_appends =
        result_restart_integer(fit, "irlba_lbd_bpro_monitored_appends"),
      irlba_lbd_bpro_threshold_reorthogonalizations =
        result_restart_integer(fit, "irlba_lbd_bpro_threshold_reorthogonalizations"),
      irlba_lbd_bpro_max_estimated_orthogonality_loss =
        result_restart_numeric(fit, "irlba_lbd_bpro_max_estimated_orthogonality_loss"),
      irlba_lbd_bpro_max_post_append_orthogonality_loss =
        result_restart_numeric(fit, "irlba_lbd_bpro_max_post_append_orthogonality_loss"),
      irlba_lbd_bpro_basis_orthogonality_loss =
        result_restart_numeric(fit, "irlba_lbd_bpro_basis_orthogonality_loss"),
      irlba_lbd_bpro_escalation_recommended =
        result_restart_logical(fit, "irlba_lbd_bpro_escalation_recommended"),
      irlba_lbd_reorth_mode =
        result_restart_character(fit, "irlba_lbd_reorth_mode"),
      irlba_lbd_one_sided_reorth_used =
        result_restart_logical(fit, "irlba_lbd_one_sided_reorth_used"),
      irlba_lbd_bpro_block_size =
        result_restart_integer(fit, "irlba_lbd_bpro_block_size"),
      irlba_lbd_bpro_exact_orthogonality_loss =
        result_restart_numeric(fit, "irlba_lbd_bpro_exact_orthogonality_loss"),
      irlba_lbd_bpro_exact_orthogonality_passed =
        result_restart_logical(fit, "irlba_lbd_bpro_exact_orthogonality_passed"),
      irlba_lbd_bpro_guard_fallback_reason =
        result_restart_character(fit, "irlba_lbd_bpro_guard_fallback_reason"),
      irlba_lbd_retained_seed_strategy =
        result_restart_character(fit, "irlba_lbd_retained_seed_strategy"),
      irlba_lbd_retained_from_scout =
        result_restart_integer(fit, "irlba_lbd_retained_from_scout"),
      irlba_lbd_retained_padding =
        result_restart_integer(fit, "irlba_lbd_retained_padding"),
      irlba_lbd_retained_fixed_work_attempts =
        result_restart_integer(fit, "irlba_lbd_retained_fixed_work_attempts"),
      irlba_lbd_scout_matvecs =
        result_restart_integer(fit, "irlba_lbd_scout_matvecs"),
      irlba_lbd_scout_certificate_passed =
        result_restart_logical(fit, "irlba_lbd_scout_certificate_passed"),
      irlba_lbd_normal_scout_attempted =
        result_restart_logical(fit, "irlba_lbd_normal_scout_attempted"),
      irlba_lbd_normal_scout_steps =
        result_restart_character(fit, "irlba_lbd_normal_scout_steps"),
      irlba_lbd_normal_scout_chosen_steps =
        result_restart_integer(fit, "irlba_lbd_normal_scout_chosen_steps"),
      irlba_lbd_normal_scout_count =
        result_restart_integer(fit, "irlba_lbd_normal_scout_count"),
      irlba_lbd_normal_scout_side =
        result_restart_character(fit, "irlba_lbd_normal_scout_side"),
      irlba_lbd_normal_scout_materialized =
        result_restart_logical(fit, "irlba_lbd_normal_scout_materialized"),
      irlba_lbd_normal_scout_certificate_trusted =
        result_restart_logical(fit, "irlba_lbd_normal_scout_certificate_trusted"),
      irlba_lbd_normal_scout_matvecs =
        result_restart_integer(fit, "irlba_lbd_normal_scout_matvecs"),
      irlba_lbd_normal_scout_operator_matvecs =
        result_restart_integer(fit, "irlba_lbd_normal_scout_operator_matvecs"),
      irlba_lbd_normal_scout_iterations =
        result_restart_integer(fit, "irlba_lbd_normal_scout_iterations"),
      irlba_lbd_normal_scout_accounted_seconds =
        result_restart_numeric(fit, "irlba_lbd_normal_scout_accounted_seconds"),
      irlba_lbd_normal_scout_polish_matvecs =
        result_restart_integer(fit, "irlba_lbd_normal_scout_polish_matvecs"),
      irlba_lbd_normal_scout_polish_iterations =
        result_restart_integer(fit, "irlba_lbd_normal_scout_polish_iterations"),
      irlba_lbd_normal_scout_polish_accounted_seconds =
        result_restart_numeric(fit, "irlba_lbd_normal_scout_polish_accounted_seconds"),
      certified_attempt = result_certified_attempt(fit),
      final_attempt_matvecs = result_restart_integer(fit, "final_attempt_matvecs"),
      final_attempt_ortho_passes = result_restart_integer(fit, "final_attempt_ortho_passes"),
      total_ortho_passes = result_restart_integer(fit, "total_ortho_passes"),
      fallback_attempted = result_restart_logical(fit, "fallback_attempted"),
      fallback_used = result_restart_logical(fit, "fallback_used"),
      fallback_method = result_restart_field(fit, "fallback_method"),
      gram_max_backward_error = result_restart_numeric(fit, "gram_max_backward_error"),
      gram_certificate_passed = result_restart_logical(fit, "gram_certificate_passed"),
      gram_dimension = result_restart_integer(fit, "gram_dimension"),
      native_gram_kernel = result_restart_character(fit, "native_gram_kernel"),
      native_gram_eigensolver = result_restart_character(fit, "native_gram_eigensolver"),
      native_gram_subspace_max_backward_error =
        result_restart_numeric(fit, "native_gram_subspace_max_backward_error"),
      normal_operator_implicit = result_restart_logical(fit, "normal_operator_implicit"),
      materialized_gram = result_restart_logical(fit, "materialized_gram"),
      certified_in_original_coordinates =
        result_restart_logical(fit, "certified_in_original_coordinates"),
      native_implicit_normal_lanczos_max_backward_error =
        result_restart_numeric(fit, "native_implicit_normal_lanczos_max_backward_error"),
      native_implicit_normal_lanczos_iterations =
        result_restart_integer(fit, "native_implicit_normal_lanczos_iterations"),
      explicit_gram_retry_used =
        result_restart_logical(fit, "explicit_gram_retry_used"),
      fallback_max_backward_error = result_restart_numeric(fit, "fallback_max_backward_error"),
      preconditioner_kind = result_preconditioner_field(fit, "kind"),
      preconditioner_native = result_preconditioner_field(fit, "native"),
      preconditioner_calls = result_preconditioner_calls(fit),
      seed = seed,
      pkg_version = as.character(utils::packageVersion("eigencore")),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$native_stage_accounted_seconds <- benchmark_row_sum(
    out$stage_apply_seconds,
    out$stage_recurrence_seconds,
    out$stage_reorthogonalization_seconds,
    out$stage_projected_solve_seconds
  )
  out$stage_reorthogonalization_fraction <- benchmark_safe_ratio(
    out$stage_reorthogonalization_seconds,
    out$native_stage_accounted_seconds
  )
  out$reorthogonalization_seconds_per_pass <- benchmark_safe_ratio(
    out$stage_reorthogonalization_seconds,
    out$reorthogonalization_passes
  )
  out$reorthogonalization_passes_per_iteration <- benchmark_safe_ratio(
    out$reorthogonalization_passes,
    out$final_iterations
  )
  out$native_seconds_per_matvec <- benchmark_safe_ratio(
    out$native_stage_accounted_seconds,
    out$final_matvecs
  )
  out$projected_seconds_per_check <- benchmark_safe_ratio(
    out$projected_seconds,
    out$projected_checks
  )
  out
}

tiny_gram_fixture <- function(n, width = max(4L * n, n + 16L),
                              density = 0.03, seed = 1L) {
  set.seed(seed)
  A <- Matrix::rsparsematrix(n, width, density = density)
  gram <- as.matrix(Matrix::tcrossprod(A))
  0.5 * (gram + t(gram))
}

tiny_gram_backend_fit <- function(gram, k, backend) {
  backend <- match.arg(backend, c(
    "lapack_dsyevr_selected",
    "lapack_dsyevx_selected",
    "lapack_dsyev_full",
    "lapack_dsyevd_full"
  ))
  if (identical(backend, "lapack_dsyevr_selected")) {
    eigencore:::native_dense_symmetric_eigen_selected(gram, k, largest())
  } else if (identical(backend, "lapack_dsyevx_selected")) {
    .Call(
      "eigencore_dense_symmetric_eigen_dsyevx_selected",
      gram,
      as.integer(k),
      as.integer(1L),
      PACKAGE = "eigencore"
    )
  } else {
    fit <- if (identical(backend, "lapack_dsyevd_full")) {
      eigencore:::native_dense_symmetric_eigen_dsyevd(gram)
    } else {
      eigencore:::native_dense_symmetric_eigen(gram)
    }
    idx <- order(fit$values, decreasing = TRUE)
    idx <- idx[seq_len(min(k, length(idx)))]
    list(
      values = fit$values[idx],
      vectors = fit$vectors[, idx, drop = FALSE]
    )
  }
}

benchmark_tiny_gram_eigensolvers <- function(dimensions = c(32L, 64L, 90L, 128L),
                                             ranks = c(5L, 8L, 16L),
                                             backends = c(
                                               "lapack_dsyevr_selected",
                                               "lapack_dsyevx_selected",
                                               "lapack_dsyev_full",
                                               "lapack_dsyevd_full"
                                             ),
                                             iterations = 20L,
                                             seed = 810L) {
  if (!requireNamespace("bench", quietly = TRUE)) {
    stop("bench is required for tiny Gram eigensolver benchmarks.", call. = FALSE)
  }
  rows <- list()
  out_idx <- 1L
  for (n in as.integer(dimensions)) {
    gram <- tiny_gram_fixture(n, seed = seed + n)
    oracle_rank <- min(max(as.integer(ranks)), n)
    oracle <- eigencore:::native_dense_symmetric_eigen_selected(
      gram,
      oracle_rank,
      largest()
    )
    for (k in as.integer(ranks)) {
      k <- min(k, n)
      oracle_values <- oracle$values[seq_len(k)]
      for (backend in backends) {
        fit <- tiny_gram_backend_fit(gram, k, backend)
        mark <- bench::mark(
          tiny_gram_backend_fit(gram, k, backend),
          iterations = iterations,
          check = FALSE,
          time_unit = "s",
          # bench::mark(memory = TRUE) calls utils::Rprofmem(), which errors on
          # R builds without --enable-memory-profiling (e.g. r-devel fedora).
          memory = isTRUE(capabilities("profmem")),
          filter_gc = FALSE
        )
        rows[[out_idx]] <- data.frame(
          dimension = n,
          rank = k,
          backend = backend,
          median = stats::median(as.numeric(mark$time[[1L]])),
          min = min(as.numeric(mark$time[[1L]])),
          mem_alloc = as.numeric(mark$mem_alloc[[1L]]),
          max_value_error = max(abs(fit$values - oracle_values)),
          values_sorted = all(diff(fit$values) <=
                                sqrt(.Machine$double.eps) *
                                  max(1, abs(fit$values), na.rm = TRUE)),
          stringsAsFactors = FALSE
        )
        out_idx <- out_idx + 1L
      }
    }
  }
  out <- do.call(rbind, rows)
  out$winner <- as.logical(ave(
    out$median,
    out$dimension,
    out$rank,
    FUN = function(x) as.integer(x == min(x, na.rm = TRUE))
  ))
  row.names(out) <- NULL
  out
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
