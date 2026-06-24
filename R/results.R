#' @export
print.eigencore_eigen_result <- function(x, ...) {
  cat("Partial eigen decomposition\n")
  cat("  requested:", x$requested, "\n")
  cat("  converged:", x$nconv, "\n")
  cat("  method:", x$method, "\n")
  cat("  target:", x$target, "\n")
  if (!is.null(x$restart)) {
    if (!is.null(x$restart$locking)) {
      cat("  restart:", x$restart$kind, "(", x$restart$locking, ")\n", sep = "")
    } else {
      cat("  restart:", x$restart$kind, "\n")
    }
    cat("  locked:", length(x$locked %||% integer(0)), "\n")
  }
  cat("  max residual:", format(x$certificate$max_residual), "\n")
  cat("  max backward error:", format(x$certificate$max_backward_error), "\n")
  cat("  max orthogonality loss:", format(x$certificate$max_orthogonality_loss), "\n")
  cat("  norm bound:", x$certificate$norm_bound_type, "\n")
  cat("  scale estimated:", x$certificate$scale_is_estimate, "\n")
  cat("  certificate:", if (isTRUE(x$certificate$passed)) "passed" else "failed", "\n")
  invisible(x)
}

#' @export
print.eigencore_svd_result <- function(x, ...) {
  cat("Partial SVD\n")
  cat("  requested rank:", x$requested, "\n")
  cat("  converged rank:", x$nconv, "\n")
  cat("  method:", x$method, "\n")
  cat("  target:", x$target, "\n")
  cat("  max residual:", format(x$certificate$max_residual), "\n")
  cat("  max backward error:", format(x$certificate$max_backward_error), "\n")
  cat("  max orthogonality loss:", format(x$certificate$max_orthogonality_loss), "\n")
  cat("  norm bound:", x$certificate$norm_bound_type, "\n")
  cat("  scale estimated:", x$certificate$scale_is_estimate, "\n")
  cat("  certificate:", if (isTRUE(x$certificate$passed)) "passed" else "failed", "\n")
  invisible(x)
}

#' @export
print.eigencore_gsvd_result <- function(x, ...) {
  cat("Generalized SVD\n")
  cat("  dimensions:", paste(x$dimensions, collapse = " x "), "\n")
  cat("  rank:", x$rank, "(k =", x$k, ", l =", x$l, ")\n")
  cat("  method:", x$method, "\n")
  cat("  finite values:", sum(x$finite), "\n")
  cat("  infinite values:", sum(x$infinite), "\n")
  cat("  undefined values:", sum(x$undefined), "\n")
  cat("  max residual:", format(x$certificate$max_residual), "\n")
  cat("  max backward error:", format(x$certificate$max_backward_error), "\n")
  cat("  max orthogonality loss:", format(x$certificate$max_orthogonality_loss), "\n")
  cat("  certificate:", if (isTRUE(x$certificate$passed)) "passed" else "failed", "\n")
  invisible(x)
}
