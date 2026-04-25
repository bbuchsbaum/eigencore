#' Modified Gram-Schmidt with two passes.
mgs2 <- function(X, against = NULL, tol = sqrt(.Machine$double.eps)) {
  X <- as.matrix(X)
  if (!is.null(against)) {
    against <- as.matrix(against)
    X <- reorthogonalize_against(X, against, passes = 2L)
  }

  decomp <- native_mgs2(X, tol = tol)
  Q <- decomp$Q
  R <- decomp$R
  loss <- orthogonality_loss(Q)
  if (loss > max(10 * tol, 100 * .Machine$double.eps)) {
    warning("MGS2 orthogonality loss is ", format(loss), ".", call. = FALSE)
  }

  structure(
    list(Q = Q, R = R, rank = ncol(Q), orthogonality = loss),
    class = "eigencore_orthogonalization"
  )
}

#' @keywords internal
native_mgs2 <- function(X, tol = sqrt(.Machine$double.eps)) {
  .Call("eigencore_mgs2", as.matrix(X), as.numeric(tol), PACKAGE = "eigencore")
}

#' Cholesky QR with a corrective second pass.
cholqr2 <- function(X, tol = sqrt(.Machine$double.eps)) {
  X <- as.matrix(X)
  decomp <- native_cholqr2(X)
  Q <- decomp$Q
  R <- decomp$R
  loss <- orthogonality_loss(Q)
  if (loss > max(10 * tol, 100 * .Machine$double.eps)) {
    warning("Cholesky QR2 orthogonality loss is ", format(loss), ".", call. = FALSE)
  }
  structure(
    list(Q = Q, R = R, rank = decomp$rank, orthogonality = loss),
    class = "eigencore_orthogonalization"
  )
}

#' @keywords internal
native_cholqr2 <- function(X) {
  .Call("eigencore_cholqr2", as.matrix(X), PACKAGE = "eigencore")
}

#' Orthogonalize columns in the B-inner product.
b_orthogonalize <- function(X, B, against = NULL, tol = sqrt(.Machine$double.eps)) {
  X <- as.matrix(X)
  B <- as.matrix(B)
  if (nrow(B) != nrow(X) || ncol(B) != nrow(X)) {
    stop("B must be square with dimension matching nrow(X).", call. = FALSE)
  }

  if (!is.null(against)) {
    against <- as.matrix(against)
    X <- b_reorthogonalize_against(X, against, B, passes = 2L)
  }

  decomp <- native_b_cholqr2(X, B)
  Q <- decomp$Q
  R <- decomp$R

  loss <- orthogonality_loss(Q, B = B)
  if (loss > max(10 * tol, 100 * .Machine$double.eps)) {
    warning("B-orthogonalization loss is ", format(loss), ".", call. = FALSE)
  }
  structure(
    list(Q = Q, R = R, rank = decomp$rank, orthogonality = loss, B = B),
    class = "eigencore_orthogonalization"
  )
}

#' @keywords internal
native_b_cholqr2 <- function(X, B) {
  .Call("eigencore_b_cholqr2", as.matrix(X), as.matrix(B), PACKAGE = "eigencore")
}

#' Measure orthogonality loss.
orthogonality_loss <- function(Q, B = NULL) {
  Q <- as.matrix(Q)
  if (is.null(B)) {
    .Call("eigencore_orthogonality_loss", Q, NULL, PACKAGE = "eigencore")
  } else {
    .Call("eigencore_orthogonality_loss", Q, as.matrix(B), PACKAGE = "eigencore")
  }
}

#' Rayleigh-Ritz projection in a trial basis.
rayleigh_ritz <- function(A, Q, B = NULL, target = largest(), symmetric = TRUE) {
  A <- as.matrix(A)
  Q <- as.matrix(Q)

  if (is.null(B) && isTRUE(symmetric)) {
    native <- native_rayleigh_ritz_symmetric(A, Q)
    projected_A <- native$projected
    values <- native$values
    vectors <- native$vectors
  } else {
    AQ <- A %*% Q
    projected_A <- crossprod(Q, AQ)
    projected_A <- if (isTRUE(symmetric)) (projected_A + t(projected_A)) / 2 else projected_A

    if (is.null(B)) {
      eig <- eigen(projected_A, symmetric = isTRUE(symmetric))
      values <- eig$values
      vectors <- Q %*% eig$vectors
    } else {
      B <- as.matrix(B)
      BQ <- B %*% Q
      projected_B <- crossprod(Q, BQ)
      projected_B <- (projected_B + t(projected_B)) / 2
      eig <- dense_generalized_spd_eigen(projected_A, projected_B, vectors = TRUE)
      values <- eig$values
      vectors <- Q %*% eig$vectors
    }
  }

  idx <- order_indices(values, target)
  values <- values[idx]
  vectors <- vectors[, idx, drop = FALSE]
  residuals <- if (is.null(B)) {
    dense_eigen_residuals(A, values, vectors)
  } else {
    dense_eigen_residuals(A, values, vectors, B = B)
  }

  structure(
    list(
      values = values,
      vectors = vectors,
      residuals = residuals,
      projected = projected_A,
      target = target_label(target)
    ),
    class = "eigencore_ritz"
  )
}

#' @keywords internal
native_rayleigh_ritz_symmetric <- function(A, Q) {
  .Call(
    "eigencore_rayleigh_ritz_symmetric",
    as.matrix(A),
    as.matrix(Q),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
cholqr_once <- function(X) {
  gram <- crossprod(X)
  gram <- (gram + t(gram)) / 2
  R <- chol(gram)
  Q <- X %*% backsolve(R, diag(ncol(R)))
  list(Q = Q, R = R)
}

#' @keywords internal
reorthogonalize_against <- function(X, Q, passes = 2L, workspace = NULL) {
  if (is.null(workspace)) {
    .Call(
      "eigencore_reorthogonalize_against",
      as.matrix(X),
      as.matrix(Q),
      as.integer(passes),
      PACKAGE = "eigencore"
    )
  } else {
    .Call(
      "eigencore_reorthogonalize_against_workspace",
      as.matrix(X),
      as.matrix(Q),
      as.integer(passes),
      workspace,
      PACKAGE = "eigencore"
    )
  }
}

#' @keywords internal
basis_workspace <- function(rows, basis_cols, block_cols) {
  .Call(
    "eigencore_basis_workspace_create",
    as.numeric(rows),
    as.numeric(basis_cols),
    as.numeric(block_cols),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
basis_workspace_info <- function(workspace) {
  .Call("eigencore_basis_workspace_info", workspace, PACKAGE = "eigencore")
}

#' @keywords internal
b_reorthogonalize_against <- function(X, Q, B, passes = 2L) {
  out <- X
  BQ <- B %*% Q
  for (i in seq_len(passes)) {
    out <- out - Q %*% crossprod(BQ, out)
  }
  out
}

#' @keywords internal
mgs_once <- function(X, tol = sqrt(.Machine$double.eps)) {
  n <- nrow(X)
  p <- ncol(X)
  Q <- matrix(0, n, p)
  R <- matrix(0, p, p)
  rank <- 0L

  for (j in seq_len(p)) {
    v <- X[, j]
    if (rank > 0L) {
      for (i in seq_len(rank)) {
        r <- sum(Q[, i] * v)
        R[i, j] <- R[i, j] + r
        v <- v - r * Q[, i]
      }
      for (i in seq_len(rank)) {
        r <- sum(Q[, i] * v)
        R[i, j] <- R[i, j] + r
        v <- v - r * Q[, i]
      }
    }
    rjj <- sqrt(sum(v^2))
    if (rjj > tol) {
      rank <- rank + 1L
      Q[, rank] <- v / rjj
      R[rank, j] <- rjj
    }
  }

  if (rank == 0L) {
    list(Q = matrix(0, n, 0L), R = matrix(0, 0L, p))
  } else {
    list(Q = Q[, seq_len(rank), drop = FALSE], R = R[seq_len(rank), , drop = FALSE])
  }
}
