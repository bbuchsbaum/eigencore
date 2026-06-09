test_that("explicit 2D grid Laplacian metadata routes to separable diagnostic prototype", {
  op <- eigencore:::grid_laplacian_2d_operator(5, 4)
  A <- op$metadata$matrix

  plan <- plan_solver(eigen_problem(op, target = smallest()), k = 4)
  expect_identical(plan$method, eigencore:::structured_grid_laplacian_2d_label())
  expect_identical(plan$controls$structured_operator, "separable_2d_grid_laplacian")
  expect_identical(plan$controls$recognition_policy, "explicit_internal_metadata_only")
  expect_false(plan$controls$arbitrary_sparse_claim)
  expect_false(plan$controls$materializes_dense_operator)

  fit <- eig_partial(op, k = 4, target = smallest(), tol = 1e-10)
  dense <- eigen(as.matrix(A), symmetric = TRUE)
  idx <- order(dense$values)[seq_len(4)]
  dense_values <- dense$values[idx]
  dense_vectors <- dense$vectors[, idx, drop = FALSE]

  expect_identical(fit$restart$kind, "separable_2d_grid_laplacian")
  expect_false(fit$restart$materialized_dense_operator)
  expect_true(fit$restart$certificate_in_original_coordinates)
  expect_true(fit$certificate$passed)
  expect_equal(fit$values, dense_values, tolerance = 1e-10)
  expect_equal(
    fit$vectors %*% t(fit$vectors),
    dense_vectors %*% t(dense_vectors),
    tolerance = 1e-8
  )
})

test_that("plain sparse 2D grid-looking matrices are not auto-claimed as structured grids", {
  op <- eigencore:::grid_laplacian_2d_operator(5, 4)
  A <- op$metadata$matrix

  plan <- plan_solver(eigen_problem(A, structure = hermitian(), target = smallest()), k = 4)
  expect_false(identical(plan$method, eigencore:::structured_grid_laplacian_2d_label()))
  expect_match(plan$method, "Lanczos")
})

test_that("structured 2D grid Laplacian values match RSpectra on small grids", {
  skip_if_not_installed("RSpectra")

  op <- eigencore:::grid_laplacian_2d_operator(6, 5)
  A <- op$metadata$matrix
  fit <- eig_partial(op, k = 5, target = smallest(), tol = 1e-9)
  rs <- RSpectra::eigs_sym(A, k = 5, which = "SA")

  expect_true(fit$certificate$passed)
  expect_equal(fit$values, sort(rs$values), tolerance = 1e-8)
})

test_that("path graph Laplacian primitives are well-formed and validated", {
  expect_error(
    eigencore:::path_graph_laplacian_matrix(1L),
    "integer >= 2",
    fixed = TRUE
  )
  expect_error(
    eigencore:::path_graph_laplacian_eigenbasis(NA_integer_),
    "integer >= 2",
    fixed = TRUE
  )

  L <- as.matrix(eigencore:::path_graph_laplacian_matrix(4L))
  expect_equal(diag(L), c(1, 2, 2, 1))
  expect_equal(L[cbind(2:4, 1:3)], rep(-1, 3))
  expect_equal(L[cbind(1:3, 2:4)], rep(-1, 3))
  expect_equal(rowSums(L), rep(0, 4))

  eb <- eigencore:::path_graph_laplacian_eigenbasis(4L)
  expect_equal(eb$values[1], 0, tolerance = 1e-12)
  expect_equal(diag(crossprod(eb$vectors)), rep(1, 4), tolerance = 1e-10)
  residual <- L %*% eb$vectors - sweep(eb$vectors, 2L, eb$values, `*`)
  expect_equal(sqrt(colSums(residual^2)), rep(0, 4), tolerance = 1e-10)
})

test_that("structured grid metadata rejects forged or inconsistent tags", {
  op <- eigencore:::grid_laplacian_2d_operator(4L, 5L)
  expect_equal(
    eigencore:::structured_grid_laplacian_2d_metadata(op),
    list(nx = 4L, ny = 5L, n = 20L)
  )

  op_no_tag <- op
  op_no_tag$metadata$structured_grid_laplacian_2d <- FALSE
  expect_null(eigencore:::structured_grid_laplacian_2d_metadata(op_no_tag))

  op_bad_nx <- op
  op_bad_nx$metadata$grid_nx <- 99L
  expect_null(eigencore:::structured_grid_laplacian_2d_metadata(op_bad_nx))

  op_bad_ny <- op
  op_bad_ny$metadata$grid_ny <- 2L
  expect_null(eigencore:::structured_grid_laplacian_2d_metadata(op_bad_ny))
})

test_that("structured grid support gate rejects invalid planner contexts", {
  op <- eigencore:::grid_laplacian_2d_operator(3L, 3L)
  problem <- eigen_problem(op, target = smallest())
  method <- auto()

  expect_true(eigencore:::structured_grid_laplacian_2d_supported(problem, method, k = 9L))
  expect_false(eigencore:::structured_grid_laplacian_2d_supported(problem, method, k = 10L))
  expect_false(eigencore:::structured_grid_laplacian_2d_supported(problem, method, k = 0L))
  expect_false(eigencore:::structured_grid_laplacian_2d_supported(problem, method, k = NA_integer_))

  expect_false(eigencore:::structured_grid_laplacian_2d_supported(
    eigen_problem(op, target = largest()),
    method,
    k = 3L
  ))
  expect_false(eigencore:::structured_grid_laplacian_2d_supported(
    eigen_problem(op, target = nearest(0.5)),
    method,
    k = 3L
  ))
  expect_false(eigencore:::structured_grid_laplacian_2d_supported(
    eigen_problem(op, metric = Matrix::Diagonal(9), target = smallest()),
    method,
    k = 3L
  ))
  expect_false(eigencore:::structured_grid_laplacian_2d_supported(
    problem,
    lanczos(max_subspace = 6L),
    k = 3L
  ))
})

test_that("structured grid planner negatives route away from separable prototype", {
  op <- eigencore:::grid_laplacian_2d_operator(3L, 3L)
  grid_label <- eigencore:::structured_grid_laplacian_2d_label()

  expect_identical(
    plan_solver(eigen_problem(op, target = smallest_magnitude()), k = 4L)$method,
    grid_label
  )

  expect_false(identical(
    plan_solver(eigen_problem(op, target = largest()), k = 4L)$method,
    grid_label
  ))
  expect_false(identical(
    plan_solver(eigen_problem(op, target = nearest(0.5)), k = 4L)$method,
    grid_label
  ))
  expect_false(identical(
    plan_solver(
      eigen_problem(op, metric = Matrix::Diagonal(9), target = smallest()),
      k = 4L
    )$method,
    grid_label
  ))
  expect_false(identical(
    plan_solver(eigen_problem(op, target = smallest()), k = 10L)$method,
    grid_label
  ))
  expect_false(identical(
    plan_solver(
      eigen_problem(op, target = smallest()),
      k = 4L,
      method = lanczos(max_subspace = 6L)
    )$method,
    grid_label
  ))
})

test_that("structured_grid_laplacian_2d_eigen truncates k and returns zero mode first", {
  nx <- 4L
  ny <- 5L
  n <- nx * ny
  eig <- eigencore:::structured_grid_laplacian_2d_eigen(nx, ny, k = n + 5L)

  expect_length(eig$values, n)
  expect_equal(nrow(eig$mode_pairs), n)
  expect_equal(ncol(eig$vectors), n)
  expect_equal(eig$values[1], 0, tolerance = 1e-12)
  expect_equal(eig$mode_pairs[1, ], c(1L, 1L))

  A <- as.matrix(eigencore:::grid_laplacian_2d_matrix(nx, ny))
  residual <- A %*% eig$vectors - sweep(eig$vectors, 2L, eig$values, `*`)
  expect_equal(sqrt(colSums(residual^2)), rep(0, n), tolerance = 1e-10)
})
