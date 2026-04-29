# Result builders shared by per-method solve_*() helpers.
# These take the iterate produced by an iterative kernel and the certificate
# computed against it, and produce a classed eigencore_eigen_result or
# eigencore_svd_result. Each per-method helper still owns its own
# method-specific restart/diagnostics fields; the builders supply the common
# skeleton (values/vectors/residuals/backward_error/orthogonality/nconv,
# requested/method/target/plan/certificate/warnings).

#' @keywords internal
make_eigen_result <- function(values,
                              vectors,
                              certificate,
                              iter,
                              requested,
                              method_label,
                              target_label_value,
                              plan,
                              warnings,
                              extras = list()) {
  base <- list(
    values = values,
    vectors = vectors,
    residuals = certificate$residuals,
    backward_error = certificate$backward_error,
    orthogonality = certificate$orthogonality,
    nconv = sum(certificate$converged),
    requested = requested,
    iterations = iter$iterations %||% 1L,
    matvecs = iter$matvecs %||% 0L,
    method = method_label,
    target = target_label_value,
    plan = plan,
    certificate = certificate,
    warnings = warnings
  )
  result <- modifyList(base, extras %||% list())
  class(result) <- "eigencore_eigen_result"
  result
}

#' @keywords internal
make_svd_result <- function(d,
                            u,
                            v,
                            certificate,
                            iter,
                            requested,
                            method_label,
                            target_label_value,
                            plan,
                            warnings,
                            extras = list()) {
  base <- list(
    d = d,
    u = u,
    v = v,
    values = d,
    residuals = certificate$residuals,
    backward_error = certificate$backward_error,
    orthogonality = certificate$orthogonality,
    nconv = sum(certificate$converged),
    requested = requested,
    iterations = iter$iterations %||% 1L,
    matvecs = iter$matvecs %||% 0L,
    stage_seconds = iter$stage_seconds %||% iter$restart$stage_seconds %||% numeric(),
    method = method_label,
    target = target_label_value,
    plan = plan,
    certificate = certificate,
    restart = iter$restart,
    warnings = warnings
  )
  result <- modifyList(base, extras %||% list())
  class(result) <- "eigencore_svd_result"
  result
}
