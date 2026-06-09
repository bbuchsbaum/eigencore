# Internal diagnostic support for explicitly tagged separable 2D path-grid
# Laplacians. This is intentionally unexported: arbitrary sparse Hermitian
# matrices must not be auto-claimed as grid Laplacians.

#' @keywords internal
path_graph_laplacian_matrix <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("path graph Laplacian dimension must be an integer >= 2.", call. = FALSE)
  }
  Matrix::bandSparse(
    n,
    k = c(-1, 0, 1),
    diagonals = list(
      rep(-1, n - 1L),
      c(1, rep(2, n - 2L), 1),
      rep(-1, n - 1L)
    )
  )
}

#' @keywords internal
path_graph_laplacian_eigenbasis <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("path graph Laplacian dimension must be an integer >= 2.", call. = FALSE)
  }
  modes <- 0:(n - 1L)
  values <- 2 - 2 * cos(pi * modes / n)
  i <- seq_len(n) - 0.5
  vectors <- outer(i, modes, function(x, mode) cos(pi * x * mode / n))
  vectors[, 1L] <- 1 / sqrt(n)
  if (n > 1L) {
    vectors[, -1L] <- sqrt(2 / n) * vectors[, -1L, drop = FALSE]
  }
  list(values = values, vectors = vectors)
}

#' @keywords internal
grid_laplacian_2d_matrix <- function(nx, ny) {
  nx <- as.integer(nx)
  ny <- as.integer(ny)
  Lx <- path_graph_laplacian_matrix(nx)
  Ly <- path_graph_laplacian_matrix(ny)
  Ix <- Matrix::Diagonal(nx)
  Iy <- Matrix::Diagonal(ny)
  Matrix::kronecker(Iy, Lx) + Matrix::kronecker(Ly, Ix)
}

#' @keywords internal
grid_laplacian_2d_operator <- function(nx, ny) {
  nx <- as.integer(nx)
  ny <- as.integer(ny)
  A <- methods::as(grid_laplacian_2d_matrix(nx, ny), "dgCMatrix")
  linear_operator(
    dim = dim(A),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
      Z <- alpha * (A %*% X)
      if (is.null(Y) || beta == 0) Z else Z + beta * Y
    },
    dtype = "double",
    structure = hermitian(),
    name = "separable_2d_grid_laplacian",
    metadata = list(
      matrix = A,
      storage = "dgCMatrix",
      native = FALSE,
      structured_grid_laplacian_2d = TRUE,
      grid_laplacian_kind = "path_kronecker_sum",
      grid_nx = nx,
      grid_ny = ny,
      frobenius_norm = sqrt(sum(methods::slot(A, "x")^2))
    )
  )
}

#' @keywords internal
structured_grid_laplacian_2d_metadata <- function(op) {
  op <- as_operator(op)
  if (!isTRUE(op$metadata$structured_grid_laplacian_2d)) {
    return(NULL)
  }
  nx <- as.integer(op$metadata$grid_nx)
  ny <- as.integer(op$metadata$grid_ny)
  if (length(nx) != 1L || is.na(nx) || nx < 2L ||
      length(ny) != 1L || is.na(ny) || ny < 2L ||
      prod(c(nx, ny)) != op$dim[1L] || op$dim[1L] != op$dim[2L]) {
    return(NULL)
  }
  list(nx = nx, ny = ny, n = nx * ny)
}

#' @keywords internal
structured_grid_laplacian_2d_supported <- function(problem, method, k = NULL) {
  if (!inherits(method, "eigencore_method") || !identical(method$kind, "auto")) {
    return(FALSE)
  }
  if (!identical(problem$structure$kind, "hermitian") || !is.null(problem$metric)) {
    return(FALSE)
  }
  target_kind <- if (inherits(problem$target, "eigencore_target")) {
    problem$target$kind
  } else {
    "largest"
  }
  if (!target_kind %in% c("smallest", "smallest_magnitude")) {
    return(FALSE)
  }
  meta <- structured_grid_laplacian_2d_metadata(problem$A)
  if (is.null(meta)) {
    return(FALSE)
  }
  k <- as.integer(k %||% 1L)
  length(k) == 1L && !is.na(k) && k >= 1L && k <= meta$n
}

#' @keywords internal
structured_grid_laplacian_2d_controls <- function(problem, k) {
  meta <- structured_grid_laplacian_2d_metadata(problem$A)
  list(
    structured_operator = "separable_2d_grid_laplacian",
    grid_laplacian_kind = problem$A$metadata$grid_laplacian_kind,
    grid_nx = meta$nx,
    grid_ny = meta$ny,
    grid_vertices = meta$n,
    requested = as.integer(k),
    prototype = TRUE,
    promotion_status = "diagnostic_v2_future_only",
    recognition_policy = "explicit_internal_metadata_only",
    arbitrary_sparse_claim = FALSE,
    materializes_dense_operator = FALSE,
    certificate_policy = "original-coordinate residual certificate against sparse operator"
  )
}

#' @keywords internal
structured_grid_laplacian_2d_eigen <- function(nx, ny, k) {
  ex <- path_graph_laplacian_eigenbasis(nx)
  ey <- path_graph_laplacian_eigenbasis(ny)
  values_2d <- outer(ex$values, ey$values, "+")
  ord <- order(as.vector(values_2d), decreasing = FALSE)
  ord <- ord[seq_len(min(as.integer(k), length(ord)))]
  pairs <- arrayInd(ord, dim(values_2d))

  values <- values_2d[ord]
  vectors <- matrix(0, nrow = nx * ny, ncol = length(ord))
  for (j in seq_along(ord)) {
    ix <- pairs[j, 1L]
    iy <- pairs[j, 2L]
    vectors[, j] <- kronecker(ey$vectors[, iy], ex$vectors[, ix])
  }
  list(values = values, vectors = vectors, mode_pairs = pairs)
}
