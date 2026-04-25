#' @export
print.eigencore_eigen_result <- function(x, ...) {
  cat("Partial eigen decomposition\n")
  cat("  requested:", x$requested, "\n")
  cat("  converged:", x$nconv, "\n")
  cat("  method:", x$method, "\n")
  cat("  target:", x$target, "\n")
  if (!is.null(x$restart)) {
    cat("  restart:", x$restart$kind, "(", x$restart$locking, ")\n", sep = "")
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
