#!/usr/bin/env Rscript

source("inst/benchmarks/_helpers.R")

args <- benchmark_args()
iterations <- if (!is.na(args$iterations)) args$iterations else if (args$quick) 1L else 5L
tol <- 1e-8

shift_invert_path_laplacian <- function(n) {
  Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, n - 1L), c(1, rep(2, n - 2L), 1), rep(-1, n - 1L))
  )
}

shift_invert_symmetric_with_spectrum <- function(values, seed = 1L) {
  set.seed(seed)
  q <- qr.Q(qr(matrix(rnorm(length(values)^2), length(values))))
  q %*% diag(values, length(values)) %*% t(q)
}

shift_invert_general_sparse <- function(n, seed = 1L) {
  methods::as(
    Matrix::Matrix(
      shift_invert_symmetric_with_spectrum(seq_len(n), seed = seed),
      sparse = TRUE
    ),
    "CsparseMatrix"
  )
}

shift_invert_matrix_free_user_solve <- function(n, sigma, seed = 1L) {
  A <- shift_invert_symmetric_with_spectrum(seq(0.1, n - 0.5, length.out = n), seed = seed)
  shifted <- A - sigma * diag(n)
  Aop <- linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (A %*% X)
      if (!is.null(Y) && beta != 0) {
        out <- out + beta * Y
      }
      out
    },
    structure = hermitian(),
    metadata = list(frobenius_norm = sqrt(sum(A^2)))
  )
  list(
    A = Aop,
    B = NULL,
    solve = function(X) base::solve(shifted, X)
  )
}

shift_invert_cases <- function(quick = FALSE) {
  if (quick) {
    return(list(
      list(
        case = "dense_standard_mid",
        n = 40L, k = 3L, sigma = 0.35,
        expected_native = TRUE,
        expected_label_kind = "dense_lu_native",
        build = function(n) list(A = as.matrix(shift_invert_path_laplacian(n)), B = NULL)
      ),
      list(
        case = "dense_generalized_mid",
        n = 24L, k = 3L, sigma = 2.75,
        expected_native = TRUE,
        expected_label_kind = "dense_lu_generalized_native",
        build = function(n) {
          list(
            A = shift_invert_symmetric_with_spectrum(seq(1, n), seed = 1701L + n),
            B = diag(seq(1.5, 2.5, length.out = n), n)
          )
        }
      ),
      list(
        case = "sparse_tridiagonal_native",
        n = 50L, k = 3L, sigma = 0.01,
        expected_native = TRUE,
        expected_label_kind = "tridiagonal_thomas_native",
        build = function(n) list(A = shift_invert_path_laplacian(n), B = NULL)
      ),
      list(
        case = "diagonal_standard_native",
        n = 60L, k = 3L, sigma = 21.5,
        expected_native = TRUE,
        expected_label_kind = "tridiagonal_thomas_native",
        build = function(n) list(A = Matrix::Diagonal(n, x = seq_len(n)), B = NULL)
      ),
      list(
        case = "sparse_tridiagonal_generalized_native",
        n = 40L, k = 3L, sigma = 0.01,
        expected_native = TRUE,
        expected_label_kind = "tridiagonal_thomas_generalized_native",
        build = function(n) {
          list(
            A = shift_invert_path_laplacian(n),
            B = Matrix::Diagonal(n, x = seq(1.1, 2.1, length.out = n))
          )
        }
      ),
      list(
        case = "sparse_general_reference",
        n = 18L, k = 3L, sigma = 8.4,
        expected_native = FALSE,
        expected_label_kind = "sparse_lu",
        expected_certificate = "estimated_converged",
        build = function(n) {
          list(A = shift_invert_general_sparse(n, seed = 1801L + n), B = NULL)
        }
      ),
      list(
        case = "sparse_general_diagonal_b_reference",
        n = 16L, k = 3L, sigma = 4.2,
        expected_native = FALSE,
        expected_label_kind = "sparse_lu_generalized",
        expected_certificate = "estimated_converged",
        build = function(n) {
          list(
            A = shift_invert_general_sparse(n, seed = 1901L + n),
            B = Matrix::Diagonal(n, x = seq(1.1, 2.0, length.out = n))
          )
        }
      ),
      list(
        case = "matrix_free_user_solve_reference",
        n = 18L, k = 3L, sigma = 8.4,
        expected_native = FALSE,
        expected_label_kind = "user_solve",
        expected_certificate = "passed",
        expected_external_cache = TRUE,
        build = function(n) {
          shift_invert_matrix_free_user_solve(n, sigma = 8.4, seed = 2001L + n)
        }
      )
    ))
  }

  list(
    list(
      case = "dense_standard_mid",
      n = 160L, k = 6L, sigma = 0.25,
      expected_native = TRUE,
      expected_label_kind = "dense_lu_native",
      build = function(n) list(A = as.matrix(shift_invert_path_laplacian(n)), B = NULL)
    ),
    list(
      case = "dense_generalized_mid",
      n = 80L, k = 6L, sigma = 2.75,
      expected_native = TRUE,
      expected_label_kind = "dense_lu_generalized_native",
      build = function(n) {
        list(
          A = shift_invert_symmetric_with_spectrum(seq(1, n), seed = 1701L + n),
          B = diag(seq(1.5, 2.5, length.out = n), n)
        )
      }
    ),
    list(
      case = "sparse_tridiagonal_native",
      n = 300L, k = 6L, sigma = 0.01,
      expected_native = TRUE,
      expected_label_kind = "tridiagonal_thomas_native",
      build = function(n) list(A = shift_invert_path_laplacian(n), B = NULL)
    ),
    list(
      case = "diagonal_standard_native",
      n = 300L, k = 6L, sigma = 121.5,
      expected_native = TRUE,
      expected_label_kind = "tridiagonal_thomas_native",
      build = function(n) list(A = Matrix::Diagonal(n, x = seq_len(n)), B = NULL)
    ),
    list(
      case = "sparse_tridiagonal_generalized_native",
      n = 160L, k = 6L, sigma = 0.01,
      expected_native = TRUE,
      expected_label_kind = "tridiagonal_thomas_generalized_native",
      build = function(n) {
        list(
          A = shift_invert_path_laplacian(n),
          B = Matrix::Diagonal(n, x = seq(1.1, 2.1, length.out = n))
        )
      }
    ),
    list(
      case = "sparse_general_reference",
      n = 60L, k = 6L, sigma = 25.4,
      expected_native = FALSE,
      expected_label_kind = "sparse_lu",
      expected_certificate = "estimated_converged",
      build = function(n) {
        list(A = shift_invert_general_sparse(n, seed = 1801L + n), B = NULL)
      }
    ),
    list(
      case = "sparse_general_diagonal_b_reference",
      n = 50L, k = 6L, sigma = 9.5,
      expected_native = FALSE,
      expected_label_kind = "sparse_lu_generalized",
      expected_certificate = "estimated_converged",
      build = function(n) {
        list(
          A = shift_invert_general_sparse(n, seed = 1901L + n),
          B = Matrix::Diagonal(n, x = seq(1.1, 2.0, length.out = n))
        )
      }
    ),
    list(
      case = "matrix_free_user_solve_reference",
      n = 60L, k = 6L, sigma = 25.4,
      expected_native = FALSE,
      expected_label_kind = "user_solve",
      expected_certificate = "passed",
      expected_external_cache = TRUE,
      build = function(n) {
        shift_invert_matrix_free_user_solve(n, sigma = 25.4, seed = 2001L + n)
      }
    )
  )
}

benchmark_shift_invert_case <- function(case, iterations = 3L, tol = 1e-8,
                                        seed = 1L) {
  built <- case$build(case$n)
  timed <- run_timed({
    eig_partial(
      built$A,
      B = built$B,
      k = case$k,
      target = nearest(case$sigma),
      method = shift_invert(case$sigma, solve = built$solve %||% NULL),
      tol = tol,
      allow_dense_fallback = "never"
    )
  }, iterations = iterations, seed = seed)
  fit <- timed$value
  cert <- fit$certificate
  cache <- fit$transform$factorization_cache %||% list()
  contract <- cache$contract %||% list()
  data.frame(
    case = case$case,
    kind = "shift_invert_eigen",
    target = "nearest",
    sigma = case$sigma,
    method = fit$method,
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
    requested = case$k,
    iterations = result_iterations(fit),
    matvecs = result_matvecs(fit),
    native = result_restart_logical(fit, "native"),
    factorization_native = result_restart_logical(fit, "factorization_native"),
    cache_native = cache$native %||% NA,
    cache_label_kind = cache$label_kind %||% NA_character_,
    factorization = cache$factorization %||% NA_character_,
    factorization_cached = cache$factorization_cached %||% NA,
    condition_estimate_type = cache$condition_estimate_type %||% NA_character_,
    external_cache = cache$external_cache %||% NA,
    generalized = cache$generalized %||% !is.null(built$B),
    metric_factorization = cache$metric_factorization %||% NA_character_,
    factorization_contract_version = contract$contract_version %||% NA_character_,
    factorization_contract_provider = contract$provider %||% NA_character_,
    factorization_contract_promotion_status = contract$promotion_status %||% NA_character_,
    factorization_contract_owned_by_eigencore = contract$owned_by_eigencore %||% NA,
    factorization_contract_external_cache = contract$external_cache %||% NA,
    factorization_contract_memory_policy = contract$memory_policy %||% NA_character_,
    factorization_contract_certificate_policy = contract$certificate_policy %||% NA_character_,
    certification_problem = fit$transform$certification$problem %||% NA_character_,
    seed = seed,
    pkg_version = as.character(utils::packageVersion("eigencore")),
    stringsAsFactors = FALSE
  )
}

shift_invert_contract <- function(rows, cases) {
  out <- lapply(cases, function(case) {
    row <- rows[rows$case == case$case, , drop = FALSE]
    if (!nrow(row)) {
      return(data.frame())
    }
    row <- row[1L, , drop = FALSE]
    expected_certificate <- case$expected_certificate %||% "passed"
    certification_gate <- if (identical(expected_certificate, "estimated_converged")) {
      !isTRUE(row$certificate_passed) &&
        isTRUE(row$scale_is_estimate) &&
        row$nconv >= case$k &&
        isTRUE(row$max_backward_error <= tol)
    } else {
      isTRUE(row$certificate_passed) && row$nconv >= case$k
    }
    native_gate <- if (isTRUE(case$expected_native)) {
      isTRUE(row$native) && isTRUE(row$factorization_native) && isTRUE(row$cache_native)
    } else if (identical(case$expected_label_kind, "user_solve")) {
      !isTRUE(row$native) &&
        !isTRUE(row$cache_native) &&
        identical(row$factorization, "user_solve") &&
        isTRUE(row$external_cache)
    } else {
      !isTRUE(row$native) &&
        !isTRUE(row$cache_native) &&
        identical(row$factorization, "Matrix::lu") &&
        isTRUE(row$factorization_cached) &&
        !isTRUE(row$external_cache)
    }
    label_gate <- identical(row$cache_label_kind, case$expected_label_kind)
    original_gate <- identical(row$certification_problem, "original")
    contract_gate <- identical(
      row$factorization_contract_version,
      "shift_invert_factorization_contract_v1"
    ) &&
      identical(row$factorization_contract_certificate_policy,
                "original_coordinate_residual_required") &&
      if (isTRUE(case$expected_native)) {
        identical(row$factorization_contract_provider,
                  "eigencore_native_factorization") &&
          identical(row$factorization_contract_promotion_status,
                    "promoted_native") &&
          isTRUE(row$factorization_contract_owned_by_eigencore)
      } else if (identical(case$expected_label_kind, "user_solve")) {
        identical(row$factorization_contract_provider, "user_supplied_solve") &&
          identical(row$factorization_contract_promotion_status,
                    "reference_boundary") &&
          !isTRUE(row$factorization_contract_owned_by_eigencore) &&
          isTRUE(row$factorization_contract_external_cache)
      } else {
        identical(row$factorization_contract_provider,
                  "Matrix::lu_reference_factorization") &&
          identical(row$factorization_contract_promotion_status,
                    "reference_boundary") &&
          !isTRUE(row$factorization_contract_owned_by_eigencore) &&
          !isTRUE(row$factorization_contract_external_cache)
      }
    data.frame(
      case = case$case,
      requested = case$k,
      nconv = row$nconv,
      expected_native = isTRUE(case$expected_native),
      cache_label_kind = row$cache_label_kind,
      certification_gate = certification_gate,
      native_gate = native_gate,
      label_gate = label_gate,
      original_coordinate_gate = original_gate,
      factorization_contract_gate = contract_gate,
      passed = certification_gate && native_gate && label_gate && original_gate &&
        contract_gate,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  row.names(out) <- NULL
  out
}

cases <- filter_benchmark_cases(shift_invert_cases(args$quick), args$cases)

rows <- lapply(seq_along(cases), function(i) {
  case <- cases[[i]]
  message_benchmark_case("bench-shift-invert", case)
  benchmark_shift_invert_case(
    case,
    iterations = iterations,
    tol = tol,
    seed = 14000L + case$n + case$k + i
  )
})
rows <- do.call(rbind, rows)
row.names(rows) <- NULL
contracts <- shift_invert_contract(rows, cases)

cat("Shift-invert benchmark rows\n")
print(rows)
cat("\nShift-invert contracts\n")
print(contracts)

if (args$save) {
  message("saved rows: ", save_benchmark_result(rows, "shift-invert-rows"))
  message("saved contracts: ", save_benchmark_result(contracts, "shift-invert-contracts"))
}

if (args$strict && !all(contracts$passed)) {
  stop("Shift-invert benchmark failed native/reference contract gate.", call. = FALSE)
}
