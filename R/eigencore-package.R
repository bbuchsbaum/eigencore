#' eigencore
#'
#' Computes the top-k singular triplets or eigenpairs of large sparse and
#' structured matrices, with a numerical certificate attached to every result.
#' See [svd_partial()], [eig_partial()], and `vignette("eigencore")`.
#'
#' @keywords internal
#' @useDynLib eigencore, .registration = TRUE
#' @importFrom methods is
#' @importFrom Matrix crossprod tcrossprod
"_PACKAGE"
