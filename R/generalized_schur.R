#' Compute a dense generalized Schur decomposition
#'
#' `generalized_schur()` is eigencore's dense QZ surface for general matrix
#' pencils `A x = lambda B x`. It computes the generalized Schur pair `S`, `T`
#' and, when requested, left/right Schur vectors `Q`, `Z` such that
#' `A = Q S Z*` and `B = Q T Z*` for complex inputs, with transpose replacing
#' conjugate-transpose for real inputs. Sparse and operator inputs are not
#' silently densified.
#'
#' @param A Base dense square matrix.
#' @param B Base dense square matrix with the same dimension as `A`.
#' @param sort Optional LAPACK sorting class. Use `NULL` or `"none"` for no
#'   sorting, `"finite"` to move finite generalized eigenvalues first,
#'   `"infinite"` to move beta-zero nonzero-alpha eigenvalues first. Custom
#'   predicates and undefined alpha-zero/beta-zero sorting are not part of the
#'   public contract.
#' @param vectors Whether to compute Schur vectors `Q` and `Z`.
#' @param ... Reserved for future options.
#' @return A classed generalized Schur result with fields `S`, `T`, `Q`, `Z`,
#'   `alpha`, `beta`, `values`, `classification`, `sdim`, `method`, `plan`, and
#'   `certificate`.
#' @examples
#' A <- matrix(c(0, -1, 1, 0), 2, 2)
#' B <- diag(2)
#' qz <- generalized_schur(A, B)
#' values(qz)
#'
#' pencil <- generalized_schur(diag(c(2, 3, 0)), diag(c(1, 0, 0)),
#'                             sort = "infinite")
#' pencil$classification
#' @export
generalized_schur <- function(A, B, sort = NULL, vectors = TRUE, ...) {
  dots <- list(...)
  if (length(dots)) {
    stop(
      "unused generalized_schur() options: ",
      paste(names(dots), collapse = ", "),
      call. = FALSE
    )
  }
  A <- generalized_schur_dense_input(A, "A")
  B <- generalized_schur_dense_input(B, "B")
  eig_full_validate_square(A, B)
  sort_key <- generalized_schur_sort_key(sort)
  vectors <- generalized_schur_vectors_flag(vectors)

  raw <- if (is.complex(A) || is.complex(B)) {
    native_dense_complex_generalized_schur(
      eig_full_as_complex_matrix(A),
      eig_full_as_complex_matrix(B),
      vectors,
      sort_key$code
    )
  } else {
    native_dense_generalized_schur(A, B, vectors, sort_key$code)
  }
  generalized_schur_result(raw, vectors = vectors, sort_key = sort_key)
}

#' @keywords internal
generalized_schur_result <- function(raw, vectors, sort_key) {
  pencil <- generalized_pencil_values(raw$alpha, raw$beta)
  method <- native_dense_generalized_schur_label()
  n <- nrow(raw$S)
  plan <- list(
    problem_type = "generalized_schur",
    requested = n,
    method = method,
    target = "all",
    reasons = c("full dense generalized Schur decomposition requested"),
    fallback = "none",
    controls = list(
      full = TRUE,
      dense = TRUE,
      qz = TRUE,
      sort = sort_key$name,
      schur_vectors = isTRUE(vectors)
    )
  )
  class(plan) <- "eigencore_plan"
  certificate <- empty_certificate(
    sqrt(.Machine$double.eps),
    note = "generalized Schur decomposition returned; per-eigenvector residual certificate not computed"
  )
  out <- list(
    S = raw$S,
    T = raw$T,
    Q = if (isTRUE(vectors)) raw$Q else NULL,
    Z = if (isTRUE(vectors)) raw$Z else NULL,
    alpha = raw$alpha,
    beta = raw$beta,
    values = pencil$values,
    classification = pencil$classification,
    finite = pencil$finite,
    infinite = pencil$infinite,
    undefined = pencil$undefined,
    sdim = raw$sdim,
    sort = sort_key$name,
    method = method,
    plan = plan,
    certificate = certificate,
    warnings = "using native dense generalized Schur QZ LAPACK full decomposition"
  )
  class(out) <- "eigencore_generalized_schur_result"
  out
}

#' @keywords internal
generalized_schur_dense_input <- function(x, name) {
  if (is.matrix(x)) {
    return(x)
  }
  stop(
    name, " must be a base dense matrix for generalized_schur(); ",
    "sparse/operator full decompositions are not silently densified",
    call. = FALSE
  )
}

#' @keywords internal
generalized_schur_vectors_flag <- function(vectors) {
  if (!is.logical(vectors) || length(vectors) != 1L || is.na(vectors)) {
    stop("vectors must be TRUE or FALSE.", call. = FALSE)
  }
  vectors
}

#' @keywords internal
generalized_schur_sort_key <- function(sort) {
  if (is.null(sort)) {
    return(list(code = 0L, name = "none"))
  }
  if (!is.character(sort) || length(sort) != 1L || is.na(sort)) {
    stop(
      "unsupported generalized_schur() sort; use NULL, 'finite', ",
      "or 'infinite'",
      call. = FALSE
    )
  }
  key <- match.arg(sort, c("none", "finite", "infinite"))
  code <- switch(
    key,
    none = 0L,
    finite = 1L,
    infinite = 2L
  )
  list(code = code, name = key)
}

#' @keywords internal
native_dense_generalized_schur <- function(A, B, vectors, sort_code) {
  .Call(
    "eigencore_dense_generalized_schur",
    as.matrix(A),
    as.matrix(B),
    as.logical(vectors),
    as.integer(sort_code),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_complex_generalized_schur <- function(A, B, vectors, sort_code) {
  .Call(
    "eigencore_dense_complex_generalized_schur",
    as.matrix(A),
    as.matrix(B),
    as.logical(vectors),
    as.integer(sort_code),
    PACKAGE = "eigencore"
  )
}

#' @keywords internal
native_dense_generalized_schur_label <- function() {
  "native dense generalized Schur QZ LAPACK full"
}
