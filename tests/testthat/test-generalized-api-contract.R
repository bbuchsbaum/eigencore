generalized_api_contract_text <- function() {
  path <- testthat::test_path("..", "..", "docs", "generalized-eigen-api.md")
  if (!file.exists(path)) {
    skip("repository docs/generalized-eigen-api.md is excluded from package tarballs")
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

expect_contract_contains <- function(text, pattern) {
  expect_true(
    grepl(pattern, text, fixed = TRUE),
    info = paste("Missing generalized API contract text:", pattern)
  )
}

test_that("generalized-eigen API contract freezes eigencore naming", {
  contract <- generalized_api_contract_text()

  expect_contract_contains(contract, "`eig_partial(A, B = NULL, ...)`")
  expect_contract_contains(
    contract,
    "`eig_full(A, B = NULL, structure = NULL, vectors = TRUE, ...)`"
  )
  expect_contract_contains(
    contract,
    "`generalized_schur(A, B, sort = NULL, vectors = TRUE, ...)`"
  )
  expect_contract_contains(contract, "`generalized_svd(A, B, ...)`")
  expect_contract_contains(contract, "There is no primary `geigen()` export.")
  expect_contract_contains(
    contract,
    "`eigen_problem(metric = )` remains SPD/Hermitian-metric-only when the problem"
  )
  expect_contract_contains(
    contract,
    eigencore:::sparse_general_pencil_diagonal_arnoldi_label()
  )
  expect_contract_contains(
    contract,
    eigencore:::sparse_general_pencil_unsupported_label()
  )
  # eig_full full-decomposition native labels must stay documented and must
  # match the exact strings returned by the *_label() helpers in R/solve.R.
  expect_contract_contains(
    contract,
    eigencore:::native_dense_generalized_spd_full_label()
  )
  expect_contract_contains(
    contract,
    eigencore:::native_dense_generalized_pencil_full_label()
  )
})

test_that("metric= rejects nonsymmetric and non-Hermitian B before dispatch", {
  A <- diag(c(1, 4, 9))

  # Real nonsymmetric B whose UPPER triangle is PD (2*I): the native SPD kernel
  # runs dpotrf(uplo="U") and references only the upper triangle, so without a
  # guard this nonsymmetric B (lower-triangle entry 9) is silently solved as 2I.
  B_nonsym <- matrix(c(2, 0, 0, 0, 2, 9, 0, 0, 2), 3, 3)
  expect_false(isSymmetric(B_nonsym))
  expect_error(eigen_problem(A, metric = B_nonsym), "symmetric")
  expect_error(
    eig_partial(A, B = B_nonsym, k = 2L, target = smallest()),
    "symmetric"
  )

  # Sparse nonsymmetric dgCMatrix B is rejected too.
  B_sparse <- methods::as(B_nonsym, "CsparseMatrix")
  expect_error(
    eig_partial(A, B = B_sparse, k = 2L, target = smallest()),
    "symmetric"
  )

  # Complex non-Hermitian B (B[2,3]=1i, B[3,2]=0) is rejected.
  Ac <- diag(as.complex(c(1, 4, 9)))
  Bc <- matrix(as.complex(c(2, 0, 0, 0, 2, 1i, 0, 0, 2)), 3, 3)
  expect_false(isSymmetric(Bc))
  expect_error(eigen_problem(Ac, metric = Bc), "symmetric|Hermitian")

  # Control: a symmetric (even if not yet definiteness-checked) B is accepted by
  # the symmetry guard; definiteness is enforced downstream, not here.
  expect_no_error(eigen_problem(A, metric = diag(c(1, 2, 3))))
  expect_no_error(eigen_problem(A, metric = methods::as(diag(c(1, 2, 3)), "CsparseMatrix")))
})

test_that("generalized replacement does not export rejected aliases", {
  exports <- getNamespaceExports("eigencore")

  expect_false("geigen" %in% exports)
  expect_false("gqz" %in% exports)
  expect_false("gsvd" %in% exports)
  expect_false("qz" %in% exports)
})

test_that("current generalized planner labels match the API contract", {
  A <- diag(c(1, 4, 9, 16))
  B <- diag(c(1, 2, 3, 4))
  dense_problem <- eigen_problem(A, metric = B, target = smallest())

  expect_equal(
    plan_solver(dense_problem, k = 2L)$method,
    "native dense generalized SPD LAPACK fallback"
  )
  expect_equal(
    plan_solver(dense_problem, k = 2L, method = lobpcg(maxit = 50L))$method,
    eigencore:::native_generalized_lobpcg_label()
  )
  expect_equal(
    plan_solver(
      eigen_problem(A, metric = B, target = largest()),
      k = 2L,
      method = lanczos(max_subspace = 4L)
    )$method,
    eigencore:::native_generalized_lanczos_label()
  )
  expect_equal(
    plan_solver(
      eigen_problem(A, metric = B, target = nearest(2), transform = shift_invert(2)),
      k = 2L
    )$method,
    eigencore:::native_dense_generalized_shift_invert_label()
  )

  S <- Matrix::Diagonal(x = c(1, 4, 9, 16))
  D <- Matrix::Diagonal(x = c(1, 2, 3, 4))
  expect_equal(
    plan_solver(
      eigen_problem(S, metric = D, target = nearest(2), transform = shift_invert(2)),
      k = 2L
    )$method,
    eigencore:::native_tridiagonal_generalized_shift_invert_label()
  )

  G <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 3),
    j = c(1, 2, 2, 3),
    x = c(4, 1, 2, -1),
    dims = c(3, 3)
  )
  expect_equal(
    plan_solver(eigen_problem(G, metric = Matrix::Diagonal(3),
                              structure = general()), k = 2L)$method,
    eigencore:::sparse_general_pencil_diagonal_arnoldi_label()
  )
})
