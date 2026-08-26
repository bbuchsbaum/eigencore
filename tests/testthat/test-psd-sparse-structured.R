psd_path_laplacian <- function(n) {
  if (n == 1L) {
    return(Matrix::sparseMatrix(
      i = integer(), j = integer(), x = numeric(), dims = c(1, 1)
    ))
  }
  incidence <- Matrix::sparseMatrix(
    i = c(seq_len(n - 1L), seq_len(n - 1L) + 1L),
    j = rep(seq_len(n - 1L), 2L),
    x = c(rep(-1, n - 1L), rep(1, n - 1L)),
    dims = c(n, n - 1L)
  )
  Matrix::tcrossprod(incidence)
}

test_that("dense Gram factors use complete compact-SVD geometry in both orientations", {
  columns <- matrix(c(
    1, 0,
    0, 1,
    1, 1
  ), 3, 2, byrow = TRUE)
  colnames(columns) <- c("u", "v")
  rownames(columns) <- c("a", "b", "c")
  factor <- psd_gram_factor(columns, orientation = "columns")
  K <- tcrossprod(columns)
  direct <- psd_factor(K)

  expect_identical(factor$representation, "gram_dense")
  expect_identical(factor$method, "complete dense Gram compact-SVD PSD factor")
  expect_identical(factor$materialization$factor, "compact_svd")
  expect_true(factor$materialization$dense_n_by_n)
  expect_identical(psd_rank(factor), 2L)
  expect_identical(psd_nullity(factor), 1L)
  expect_psd_condition(
    psd_nullity(factor, "algebraic"),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence",
    "algebraic_nullity"
  )
  expect_equal(unname(psd_apply(factor, diag(3), "form")), unname(K), tolerance = 1e-12)
  for (action in c(
    "sqrt", "inverse_sqrt", "pseudoinverse",
    "image_projector", "null_projector"
  )) {
    expect_equal(
      psd_apply(factor, diag(3), action),
      psd_apply(direct, diag(3), action),
      tolerance = 1e-11
    )
  }
  expect_identical(rownames(psd_apply(factor, diag(3))), rownames(columns))

  row_factor <- psd_gram_factor(t(columns), orientation = "rows")
  expect_equal(unname(psd_apply(row_factor, diag(3), "form")), unname(K), tolerance = 1e-12)
  expect_equal(
    psd_apply(row_factor, diag(3), "image_projector"),
    psd_apply(factor, diag(3), "image_projector"),
    tolerance = 1e-12
  )

  matrix_factor <- psd_gram_factor(Matrix::Matrix(columns, sparse = FALSE))
  expect_identical(matrix_factor$representation, "gram_dense")
  expect_equal(psd_spectrum(matrix_factor), psd_spectrum(factor), tolerance = 1e-13)
})

test_that("dense Gram validation is typed and snapshot ownership is immutable", {
  expect_psd_condition(
    psd_gram_factor(matrix(1L, 2, 2)),
    "eigencore_psd_invalid_input", "invalid_dtype", "x"
  )
  expect_psd_condition(
    psd_gram_factor(matrix(complex(real = 1:4), 2, 2)),
    "eigencore_psd_invalid_input", "invalid_dtype", "x"
  )
  bad <- matrix(c(1, 2, Inf, 4), 2, 2)
  expect_psd_condition(
    psd_gram_factor(bad),
    "eigencore_psd_invalid_input", "nonfinite_input", "x"
  )
  expect_psd_condition(
    psd_gram_factor(Matrix::Diagonal(2)),
    "eigencore_psd_invalid_input", "invalid_structure", "x"
  )
  expect_psd_condition(
    psd_gram_factor(diag(2), source = "live"),
    "eigencore_psd_invalid_input", "invalid_policy", "source"
  )

  source <- matrix(c(1, 0, 1, 1), 2, 2)
  factor <- psd_gram_factor(source)
  before <- psd_apply(factor, diag(2), "form")
  source[,] <- 100
  expect_identical(psd_apply(factor, diag(2), "form"), before)
})

test_that("Gram rank classification and structural scale are rescaling equivariant", {
  base <- diag(c(1, 1e-5))
  dense_factors <- lapply(c(1e-6, 1, 1e6), function(multiplier) {
    psd_gram_factor(multiplier * base)
  })
  expect_identical(
    vapply(dense_factors, psd_rank, integer(1L)),
    rep(1L, 3L)
  )
  expect_identical(
    lapply(dense_factors, function(x) x$spectrum$category),
    rep(list(c("retained_positive", "numerical_null")), 3L)
  )

  sparse <- Matrix::sparseMatrix(
    i = 1:2, j = 1:2, x = diag(base), dims = c(2, 2)
  )
  sparse_factors <- lapply(c(1e-6, 1, 1e6), function(multiplier) {
    psd_gram_factor(multiplier * sparse)
  })
  normalized_scale <- vapply(
    seq_along(sparse_factors),
    function(i) sparse_factors[[i]]$scale / c(1e-6, 1, 1e6)[[i]]^2,
    numeric(1L)
  )
  expect_equal(normalized_scale, rep(sqrt(1 + 1e-20), 3L), tolerance = 1e-13)
  expect_true(all(vapply(
    sparse_factors,
    function(x) !x$capabilities$numerical_rank$available,
    logical(1L)
  )))
})

test_that("sparse Gram factors expose only theorem-supported capabilities", {
  dense_factor <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    0, 2
  ), 4, 2, byrow = TRUE)
  sparse_factor <- methods::as(Matrix::Matrix(dense_factor, sparse = TRUE), "dgCMatrix")
  factor <- psd_gram_factor(sparse_factor)
  K <- tcrossprod(dense_factor)
  X <- matrix(c(1, 2, 3, 4, -1, 0, 2, 1), 4, 2)

  expect_identical(factor$representation, "gram_sparse")
  expect_identical(factor$method, "structural sparse Gram PSD factor")
  expect_identical(factor$spectrum$coverage, "structural")
  expect_identical(factor$evidence$validation, "computed")
  expect_identical(factor$materialization$factor, "sparse_gram")
  expect_false(factor$materialization$dense_n_by_n)
  expect_true(psd_capabilities(factor)$form$available)
  expect_true(psd_capabilities(factor)$gram$available)
  expect_false(psd_capabilities(factor)$sqrt$available)
  expect_equal(factor$scale, norm(K, "F"), tolerance = 1e-13)
  expect_equal(psd_apply(factor, X), K %*% X, tolerance = 1e-13)
  expect_equal(psd_gram(factor, X), crossprod(X, K %*% X), tolerance = 1e-12)
  expect_equal(as_operator(factor)$apply(X), K %*% X, tolerance = 1e-13)

  for (request in list(
    function() psd_spectrum(factor),
    function() psd_rank(factor),
    function() psd_nullity(factor),
    function() psd_nullity(factor, "algebraic"),
    function() psd_apply(factor, X, "sqrt"),
    function() psd_reduce(factor, X),
    function() psd_solve(factor, X)
  )) {
    expect_s3_class(tryCatch(request(), error = identity), "eigencore_psd_incomplete_evidence")
  }

  row_factor <- psd_gram_factor(Matrix::t(sparse_factor), orientation = "rows")
  expect_equal(psd_apply(row_factor, X), K %*% X, tolerance = 1e-13)

  expected <- psd_apply(factor, diag(4))
  sparse_factor@x[] <- 50
  expect_identical(psd_apply(factor, diag(4)), expected)
})

test_that("advertised sparse Gram construction does not densify or retain square dense state", {
  n <- 800L
  sparse_factor <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = rep(seq_len(20L), length.out = n),
    x = rep(c(1, -0.5), length.out = n),
    dims = c(n, 20L)
  )
  sparse_factor <- methods::as(sparse_factor, "dgCMatrix")
  sparse_as_matrix_calls <- 0L
  invisible(trace(
    "as.matrix",
    tracer = quote({
      if (inherits(x, "sparseMatrix")) {
        sparse_as_matrix_calls <<- sparse_as_matrix_calls + 1L
      }
    }),
    print = FALSE,
    where = asNamespace("Matrix")
  ))
  on.exit(untrace("as.matrix", where = asNamespace("Matrix")), add = TRUE)

  factor <- psd_gram_factor(sparse_factor)
  laplacian_factor <- psd_laplacian(psd_path_laplacian(n))
  expect_identical(sparse_as_matrix_calls, 0L)
  state <- attr(factor, "eigencore_psd_state")
  expect_s4_class(state$factor, "dgCMatrix")
  expect_false(any(vapply(
    state,
    function(value) is.matrix(value) && identical(dim(value), c(n, n)),
    logical(1L)
  )))
  expect_lt(retained_bytes(factor), 2e6)
  expect_lt(retained_bytes(laplacian_factor), 2e6)
  expect_false(laplacian_factor$materialization$dense_n_by_n)
})

test_that("path, disconnected, and isolated Laplacians certify structural nullity", {
  path <- psd_path_laplacian(5)
  factor <- psd_laplacian(path)
  X <- matrix(seq_len(10), 5, 2)
  storage.mode(X) <- "double"
  expect_identical(factor$representation, "laplacian_sparse")
  expect_identical(factor$method, "structural sparse graph-Laplacian PSD factor")
  expect_identical(psd_nullity(factor, "algebraic"), 1L)
  expect_equal(psd_apply(factor, X), as.matrix(path %*% X), tolerance = 1e-14)
  expect_false(factor$materialization$dense_n_by_n)
  expect_true(certificate(factor)$passed)
  expect_psd_condition(
    psd_rank(factor),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", "numerical_rank"
  )
  expect_psd_condition(
    psd_apply(factor, X, "pseudoinverse"),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", "pseudoinverse"
  )

  disconnected <- Matrix::bdiag(psd_path_laplacian(3), psd_path_laplacian(2))
  disconnected <- methods::as(disconnected, "generalMatrix")
  disconnected_factor <- psd_laplacian(disconnected)
  expect_identical(psd_nullity(disconnected_factor, "algebraic"), 2L)
  expect_identical(disconnected_factor$rank_bounds, list(lower = 0L, upper = 3L))
  block <- matrix(seq_len(10), 5, 2)
  storage.mode(block) <- "double"
  component_constants <- cbind(c(rep(3, 3), rep(-2, 2)), c(rep(1, 3), rep(4, 2)))
  expect_equal(
    psd_apply(disconnected_factor, block + component_constants),
    psd_apply(disconnected_factor, block),
    tolerance = 1e-14
  )

  isolated <- Matrix::sparseMatrix(
    i = integer(), j = integer(), x = numeric(), dims = c(4, 4)
  )
  isolated_factor <- psd_laplacian(isolated)
  expect_identical(psd_nullity(isolated_factor, "algebraic"), 4L)
  expect_equal(psd_apply(isolated_factor, diag(4)), matrix(0, 4, 4))
})

test_that("Laplacian component evidence and actions are permutation equivariant", {
  source <- methods::as(
    Matrix::bdiag(psd_path_laplacian(4), psd_path_laplacian(3)),
    "generalMatrix"
  )
  permutation <- c(5, 1, 7, 3, 2, 6, 4)
  permuted <- methods::as(source[permutation, permutation], "dgCMatrix")
  factor <- psd_laplacian(source)
  permuted_factor <- psd_laplacian(permuted)
  X <- matrix(seq_len(14), 7, 2)
  storage.mode(X) <- "double"

  expect_identical(psd_nullity(factor, "algebraic"), 2L)
  expect_identical(psd_nullity(permuted_factor, "algebraic"), 2L)
  expect_equal(
    psd_apply(permuted_factor, X[permutation, , drop = FALSE]),
    psd_apply(factor, X)[permutation, , drop = FALSE],
    tolerance = 1e-14
  )
  expect_equal(permuted_factor$scale, factor$scale, tolerance = 1e-14)
})

test_that("Laplacian structural claims are invariant across weight scales", {
  base <- psd_path_laplacian(5)
  multipliers <- c(1e-10, 1, 1e10)
  factors <- lapply(multipliers, function(multiplier) {
    psd_laplacian(multiplier * base)
  })
  expect_identical(
    vapply(factors, function(x) psd_nullity(x, "algebraic"), integer(1L)),
    rep(1L, length(multipliers))
  )
  normalized <- vapply(
    seq_along(factors),
    function(i) factors[[i]]$scale / multipliers[[i]],
    numeric(1L)
  )
  expect_equal(normalized, rep(norm(as.matrix(base), "F"), 3L), tolerance = 1e-13)
})

test_that("Laplacian canonicalization records each source-to-action defect", {
  source <- methods::as(psd_path_laplacian(3), "generalMatrix")
  source[1, 1] <- source[1, 1] + 5e-7
  policy <- psd_policy(
    positivity = psd_tolerance(abs = 1e-6, rel = 0),
    structure_repair = "canonicalize"
  )
  factor <- psd_laplacian(source, policy = policy)
  repaired <- psd_path_laplacian(3)
  source_defect <- norm(as.matrix(source - repaired), "F")
  row_sum_defect <- max(abs(Matrix::rowSums(source)))

  expect_true(certificate(factor)$repair_applied)
  expect_equal(certificate(factor)$repair_defect, source_defect, tolerance = 1e-14)
  expect_equal(certificate(factor)$source_action_defect, source_defect, tolerance = 1e-14)
  expect_equal(psd_apply(factor, diag(3)), as.matrix(repaired), tolerance = 1e-14)
  expect_equal(
    factor$evidence$details$admitted_row_sum_defect,
    row_sum_defect,
    tolerance = 1e-14
  )

  rejecting <- psd_policy(
    positivity = psd_tolerance(abs = 1e-6, rel = 0),
    structure_repair = "reject"
  )
  expect_psd_condition(
    psd_laplacian(source, policy = rejecting),
    "eigencore_psd_invalid_input", "invalid_structure", "row_sums"
  )
})

test_that("non-Laplacian sparse structures fail with quantitative evidence", {
  positive_edge <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2), j = c(1, 2, 1, 2),
    x = c(-1, 1, 1, -1), dims = c(2, 2)
  )
  positive_error <- expect_psd_condition(
    psd_laplacian(methods::as(positive_edge, "dgCMatrix")),
    "eigencore_psd_invalid_input", "invalid_structure", "off_diagonal"
  )
  expect_gt(positive_error$defect, 0)

  bad_rows <- methods::as(psd_path_laplacian(3), "generalMatrix")
  bad_rows[1, 1] <- bad_rows[1, 1] + 1e-3
  row_error <- expect_psd_condition(
    psd_laplacian(bad_rows),
    "eigencore_psd_invalid_input", "invalid_structure", "row_sums"
  )
  expect_gt(row_error$defect, row_error$threshold)

  asymmetric <- methods::as(psd_path_laplacian(3), "generalMatrix")
  asymmetric[1, 2] <- asymmetric[1, 2] - 1e-3
  asymmetry_error <- expect_psd_condition(
    psd_laplacian(asymmetric),
    "eigencore_psd_asymmetry_error", "asymmetric_input", "x"
  )
  expect_gt(asymmetry_error$defect, asymmetry_error$threshold)

  expect_psd_condition(
    psd_laplacian(Matrix::Diagonal(3)),
    "eigencore_psd_invalid_input", "invalid_structure", "x"
  )
  expect_psd_condition(
    psd_laplacian(methods::as(psd_path_laplacian(3), "generalMatrix"), source = "live"),
    "eigencore_psd_invalid_input", "invalid_policy", "source"
  )
})

test_that("sparse snapshots survive caller mutation and portable serialization", {
  source <- methods::as(psd_path_laplacian(4), "generalMatrix")
  factor <- psd_laplacian(source)
  expected <- psd_apply(factor, diag(4))
  source@x[] <- 100
  expect_identical(psd_apply(factor, diag(4)), expected)

  path <- tempfile(fileext = ".rds")
  saveRDS(factor, path)
  restored <- readRDS(path)
  expect_identical(psd_apply(restored, diag(4)), expected)
  expect_identical(operator_identity(restored), operator_identity(factor))
})

test_that("opaque matrix-free identity is not promoted to PSD evidence", {
  calls <- 0L
  operator <- linear_operator(
    dim = c(3, 3),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      alpha * X
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      calls <<- calls + 1L
      alpha * X
    },
    structure = hermitian(),
    operator_id = "opaque-psd-claim",
    revision = "r1"
  )
  error <- expect_psd_condition(
    psd_factor(operator),
    "eigencore_psd_incomplete_evidence", "incomplete_evidence", "x"
  )
  expect_identical(calls, 0L)
  expect_identical(error$evidence$validation, "unavailable")
  expect_identical(error$representation, "opaque_operator")
})
