#' eigencore
#'
#' Computes the top-k singular triplets or eigenpairs of large sparse and
#' structured matrices, with a numerical certificate attached to every result.
#' It also validates and factors finite real symmetric positive-semidefinite
#' forms, including singular image-space geometry and structural sparse Gram
#' and graph-Laplacian paths. See [svd_partial()], [eig_partial()],
#' [psd_factor()], `vignette("eigencore")`, and
#' `vignette("psd-geometry")`.
#'
#' @keywords internal
#' @useDynLib eigencore, .registration = TRUE
#' @importFrom methods is
#' @importFrom Matrix crossprod tcrossprod
"_PACKAGE"
