# exported function signatures are frozen

    Code
      writeLines(vapply(exports, api_signature, character(1L)))
    Output
      adjoint(x, ...)
      alpha_beta(x, ...)
      as_operator(x, ...)
      auto()
      backward_error(x, ...)
      both_ends(k_low, k_high)
      center(A, rows = FALSE, columns = TRUE, row_means = NULL, col_means = NULL, name = NULL)
      certificate(x, ...)
      check_adjoint(A, trials = 20, tol = 1e-12, seed = NULL)
      compose(A, B, name = NULL)
      crossprod_operator(A, name = NULL)
      diagnostics(x, ...)
      eig_full(A, B = NULL, structure = NULL, vectors = TRUE, tol = 1e-08, allow_dense_fallback = c("auto", "never", "always"), ...)
      eig_partial(A, k, target = largest(), B = NULL, method = auto(), tol = 1e-08, maxit = NULL, vectors = TRUE, seed = NULL, certify = TRUE, allow_dense_fallback = c("auto", "never", "always"), initial_subspace = NULL)
      eigen_problem(A, metric = NULL, structure = NULL, target = largest(), transform = NULL)
      eigs(A, k, which = "LM", opts = list(), ...)
      eigs_sym(A, k, which = "LA", opts = list(), ...)
      euclidean(dim, dtype = "double")
      general()
      generalized_schur(A, B, sort = NULL, vectors = TRUE, ...)
      generalized_svd(A, B, tol = 1e-08, ...)
      golub_kahan(max_subspace = NULL, reorthogonalize = TRUE)
      hermitian()
      lanczos(max_subspace = NULL, max_restarts = NULL, block = 1L, check_stride = 0L, reorthogonalize = TRUE)
      largest()
      largest_imaginary()
      largest_magnitude()
      largest_real()
      left_vectors(x, ...)
      linear_operator(dim, apply, apply_adjoint = NULL, dtype = "double", structure = general(), name = NULL, metadata = list())
      lobpcg(maxit = 200L, preconditioner = NULL, constraints = NULL)
      nearest(sigma)
      plan_solver(problem, ...)
      randomized(oversample = 10, n_iter = 2, block = NULL, normalizer = c("qr", "lu", "none"), refine = TRUE)
      right_vectors(x, ...)
      scale_cols(A, weights, name = NULL)
      scale_rows(A, weights, name = NULL)
      shift_invert(sigma, solve = NULL, factorization = NULL)
      shifted_cholesky_preconditioner(A, shift = 0)
      shifted_diagonal_preconditioner(A, shift = 0)
      shifted_tridiagonal_preconditioner(A, shift = 0)
      smallest()
      smallest_imaginary()
      smallest_magnitude()
      smallest_real()
      svd_partial(A, rank, target = largest(), method = auto(), tol = 1e-08, vectors = c("both", "left", "right", "none"), seed = NULL, certify = TRUE, allow_dense_fallback = c("auto", "never", "always"))
      svd_problem(A, domain = NULL, codomain = NULL, target = largest())
      svds(A, k, nu = k, nv = k, opts = list(), ...)
      symmetric_operator(A, validate = TRUE, tol = 1e-10)
      values(x, ...)
      vectors(x, ...)

# S3 method registrations are frozen

    Code
      writeLines(regs)
    Output
      adjoint.eigencore_operator
      as_operator.default
      as_operator.eigencore_operator
      as_operator.matrix
      plan_solver.eigencore_eigen_problem
      plan_solver.eigencore_svd_problem
      print.eigencore_benchmark
      print.eigencore_certificate
      print.eigencore_eigen_result
      print.eigencore_gsvd_result
      print.eigencore_operator
      print.eigencore_plan
      print.eigencore_svd_result
      print.eigencore_validation
      residuals.eigencore_certificate
      residuals.eigencore_eigen_result
      residuals.eigencore_svd_result
      solve.eigencore_eigen_problem
      solve.eigencore_svd_problem

# result and certificate field names are frozen

    Code
      cat("eigen result:\n")
    Output
      eigen result:
    Code
      writeLines(sort(names(efit)))
    Output
      backward_error
      certificate
      iterations
      matvecs
      method
      nconv
      orthogonality
      plan
      requested
      residuals
      restart
      target
      values
      vectors
      warnings
    Code
      cat("svd result:\n")
    Output
      svd result:
    Code
      writeLines(sort(names(sfit)))
    Output
      backward_error
      certificate
      d
      iterations
      matvecs
      method
      nconv
      orthogonality
      plan
      requested
      residuals
      stage_seconds
      target
      u
      v
      values
      warnings
    Code
      cat("certificate:\n")
    Output
      certificate:
    Code
      writeLines(sort(names(efit$certificate)))
    Output
      backward_error
      certificate_type
      converged
      failed_indices
      max_backward_error
      max_orthogonality_loss
      max_residual
      norm_bound_type
      notes
      orthogonality
      orthogonality_passed
      orthogonality_required
      orthogonality_tolerance
      passed
      residuals
      scale
      scale_is_estimate
      tolerance

