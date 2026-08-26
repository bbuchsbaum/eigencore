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
      linear_operator(dim, apply, apply_adjoint = NULL, dtype = "double", structure = general(), name = NULL, metadata = list(), operator_id = NULL, revision = NULL, portable = FALSE)
      lobpcg(maxit = 200L, preconditioner = NULL, constraints = NULL)
      nearest(sigma)
      operator_identity(x)
      plan_solver(problem, ...)
      psd_apply(x, X, action = c("form", "sqrt", "inverse_sqrt", "pseudoinverse", "image_projector", "null_projector"))
      psd_capabilities(x)
      psd_factor(x, policy = psd_policy(), source = c("snapshot", "live"))
      psd_gram(x, X, Y = NULL)
      psd_gram_factor(x, orientation = c("columns", "rows"), policy = psd_policy(), source = c("snapshot", "live"))
      psd_identity(dim, policy = psd_policy())
      psd_laplacian(x, policy = psd_policy(), source = c("snapshot", "live"))
      psd_lift(x, Z)
      psd_nullity(x, type = c("numerical", "algebraic"))
      psd_operator(x, action = c("form", "sqrt", "inverse_sqrt", "pseudoinverse", "image_projector", "null_projector"))
      psd_orthonormalize(x, X, required_rank = NULL, tolerance = NULL)
      psd_policy(symmetry = psd_tolerance(), positivity = psd_tolerance(rel = 64 * .Machine$double.eps), rank = psd_tolerance(), rhs = psd_tolerance(), scale = "frobenius", symmetry_repair = c("average", "reject"), negative_repair = c("clip", "reject"), structure_repair = c("canonicalize", "reject"))
      psd_rank(x)
      psd_reduce(x, X)
      psd_reduced_operator(x, A)
      psd_solve(x, B, tolerance = NULL)
      psd_spectrum(x, repaired = FALSE)
      psd_tolerance(abs = 0, rel = sqrt(.Machine$double.eps))
      randomized(oversample = 10, n_iter = 2, block = NULL, normalizer = c("qr", "lu", "none"), refine = TRUE)
      restart_state(x, retention = c("basis", "same_operator"))
      retained_bytes(x)
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
      work(x, ...)

# S3 method registrations are frozen

    Code
      writeLines(regs)
    Output
      adjoint.eigencore_operator
      as_operator.default
      as_operator.eigencore_operator
      as_operator.eigencore_psd_factor
      as_operator.matrix
      plan_solver.eigencore_eigen_problem
      plan_solver.eigencore_svd_problem
      print.eigencore_benchmark
      print.eigencore_certificate
      print.eigencore_eigen_result
      print.eigencore_gsvd_result
      print.eigencore_operator
      print.eigencore_plan
      print.eigencore_psd_block_result
      print.eigencore_psd_capabilities
      print.eigencore_psd_certificate
      print.eigencore_psd_factor
      print.eigencore_psd_policy
      print.eigencore_psd_solve_result
      print.eigencore_svd_result
      print.eigencore_validation
      residuals.eigencore_certificate
      residuals.eigencore_eigen_result
      residuals.eigencore_svd_result
      solve.eigencore_eigen_problem
      solve.eigencore_plan
      solve.eigencore_svd_problem

# result and certificate field names are frozen

    Code
      cat("eigen result:\n")
    Output
      eigen result:
    Code
      writeLines(sort(names(efit)))
    Output
      actual_method
      backward_error
      certificate
      fallback_reason
      fallback_used
      iterations
      matvecs
      memory
      method
      nconv
      orthogonality
      plan
      planned_method
      requested
      residuals
      restart
      restart_state
      state_transition
      target
      values
      vectors
      warnings
      work
    Code
      cat("svd result:\n")
    Output
      svd result:
    Code
      writeLines(sort(names(sfit)))
    Output
      actual_method
      backward_error
      certificate
      d
      fallback_reason
      fallback_used
      iterations
      matvecs
      memory
      method
      nconv
      orthogonality
      plan
      planned_method
      requested
      residuals
      restart_state
      stage_seconds
      state_transition
      target
      u
      v
      values
      warnings
      work
    Code
      cat("dense general-pencil result:\n")
    Output
      dense general-pencil result:
    Code
      writeLines(sort(names(pencil_fit)))
    Output
      actual_method
      alpha
      backward_error
      beta
      biorthogonality
      certificate
      classification
      classification_policy
      conditioning
      fallback_reason
      fallback_used
      finite
      generalized
      infinite
      iterations
      left_certificate
      left_vectors
      matvecs
      memory
      method
      nconv
      orthogonality
      plan
      planned_method
      requested
      residuals
      restart_state
      state_transition
      target
      undefined
      values
      vectors
      warnings
      work
    Code
      cat("generalized Schur result:\n")
    Output
      generalized Schur result:
    Code
      writeLines(sort(names(qz_fit)))
    Output
      Q
      S
      T
      Z
      alpha
      beta
      certificate
      classification
      classification_policy
      finite
      infinite
      method
      plan
      sdim
      sort
      undefined
      values
      warnings
    Code
      cat("sparse general-pencil result:\n")
    Output
      sparse general-pencil result:
    Code
      writeLines(sort(names(sparse_pencil_fit)))
    Output
      actual_method
      alpha
      backward_error
      beta
      certificate
      classification
      classification_policy
      fallback_reason
      fallback_used
      finite
      generalized
      infinite
      iterations
      locked
      matvecs
      memory
      method
      nconv
      orthogonality
      plan
      planned_method
      requested
      residuals
      restart
      restart_state
      restarts
      right_hand_pencil
      state_transition
      target
      transform
      undefined
      values
      vectors
      warnings
      work
    Code
      cat("GSVD result:\n")
    Output
      GSVD result:
    Code
      writeLines(sort(names(gsvd_fit)))
    Output
      A_factor
      B_factor
      D1
      D2
      Q
      R
      U
      V
      alpha
      backward_error
      beta
      certificate
      classification
      classification_policy
      dimensions
      finite
      infinite
      iterations
      k
      l
      matvecs
      method
      nconv
      orthogonality
      plan
      rank
      requested
      residuals
      target
      undefined
      values
      warnings
      zero_R
    Code
      cat("nonsymmetric left-vector result:\n")
    Output
      nonsymmetric left-vector result:
    Code
      writeLines(sort(names(nonsymmetric_fit)))
    Output
      actual_method
      backward_error
      biorthogonality
      certificate
      fallback_reason
      fallback_used
      iterations
      left_certificate
      left_vectors
      locked
      matvecs
      memory
      method
      nconv
      orthogonality
      plan
      planned_method
      requested
      residuals
      restart
      restart_state
      restarts
      right_vectors
      state_transition
      target
      values
      vectors
      warnings
      work
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

