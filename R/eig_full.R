#' Compute a full dense eigendecomposition
#'
#' `eig_full()` is the dense/full eigencore surface for standard and
#' generalized eigenproblems. Sparse and operator inputs are not silently
#' densified.
#'
#' @param A Base dense square matrix.
#' @param B Optional base dense square matrix for generalized problems.
#' @param structure Optional structure descriptor. Use [general()] to force the
#'   general-pencil path for dense generalized inputs.
#' @param vectors Whether to return right eigenvectors.
#' @param tol Certification tolerance.
#' @param allow_dense_fallback Reserved dense fallback policy. Sparse/operator
#'   inputs still fail unless a future issue explicitly opens an opt-in dense
#'   fallback contract for `eig_full()`.
#' @param ... Reserved for future options.
#' @return An `eigencore_eigen_result`. Dense general-pencil results
#'   additionally carry `alpha`, `beta`, `classification`,
#'   `classification_policy`, `left_vectors` (left generalized eigenvectors
#'   satisfying `w^H A = lambda w^H B`), and `conditioning`. For real pencils
#'   the decomposition runs through the expert LAPACK driver `DGGEVX` with
#'   balancing, so `conditioning` contains reciprocal condition numbers
#'   `rconde`/`rcondv` and the balanced pencil norms `abnrm`/`bbnrm`. Complex
#'   pencils use `ZGGEV` (R's bundled LAPACK subset has no `ZGGEVX`), so they
#'   return left vectors but `conditioning$available` is `FALSE`.
#' @examples
#' A <- diag(c(1, 4, 9))
#' B <- diag(c(1, 2, 3))
#' fit <- eig_full(A, B = B)
#' values(fit)
#' certificate(fit)$passed
#'
#' # Force the dense general-pencil path when alpha/beta diagnostics matter.
#' pencil <- eig_full(A, B = B, structure = general())
#' pencil$alpha
#' pencil$beta
#' pencil$classification
#' @export
eig_full <- function(A, B = NULL, structure = NULL, vectors = TRUE,
                     tol = 1e-8,
                     allow_dense_fallback = c("auto", "never", "always"),
                     ...) {
  allow_dense_fallback <- match.arg(allow_dense_fallback)
  dots <- list(...)
  if (length(dots)) {
    stop("unused eig_full() options: ", paste(names(dots), collapse = ", "),
         call. = FALSE)
  }
  A <- eig_full_dense_input(A, "A", allow_dense_fallback)
  B <- if (is.null(B)) NULL else eig_full_dense_input(B, "B", allow_dense_fallback)
  eig_full_validate_square(A, B)

  if (!is.null(B) && (is.complex(A) || is.complex(B))) {
    A <- eig_full_as_complex_matrix(A)
    B <- eig_full_as_complex_matrix(B)
  }

  structure <- eig_full_resolve_structure(A, B, structure, tol = sqrt(.Machine$double.eps))
  if (is.null(B)) {
    return(eig_full_standard(A, structure = structure, vectors = vectors, tol = tol))
  }
  eig_full_generalized(A, B, structure = structure, vectors = vectors, tol = tol)
}

#' @keywords internal
eig_full_standard <- function(A, structure, vectors, tol) {
  n <- nrow(A)
  hermitian_path <- identical(structure$kind, "hermitian")
  if (hermitian_path) {
    eig <- if (is.complex(A)) {
      native_dense_complex_hermitian_eigen(A)
    } else {
      native_dense_symmetric_eigen(A)
    }
    method <- if (is.complex(A)) {
      native_dense_complex_hermitian_label()
    } else {
      "native dense Hermitian LAPACK fallback"
    }
    cert <- eig_full_certificate(A, NULL, eig$values, eig$vectors, vectors, tol,
                                 general = FALSE)
    return(eig_full_result(
      values = eig$values,
      vectors = if (isTRUE(vectors)) eig$vectors else NULL,
      certificate = cert,
      method = method,
      n = n,
      warnings = "using native dense full Hermitian LAPACK decomposition",
      extras = list(generalized = FALSE)
    ))
  }

  eig <- if (is.complex(A)) {
    native_dense_complex_general_eigen(A)
  } else {
    eigen(A, symmetric = FALSE)
  }
  cert <- eig_full_certificate(A, NULL, eig$values, eig$vectors, vectors, tol,
                               general = TRUE)
  eig_full_result(
    values = eig$values,
    vectors = if (isTRUE(vectors)) eig$vectors else NULL,
    certificate = cert,
    method = if (is.complex(A)) {
      native_dense_complex_general_label()
    } else {
      "dense LAPACK general eigen oracle (base fallback)"
    },
    n = n,
    warnings = if (is.complex(A)) {
      "using native dense complex general LAPACK full decomposition"
    } else {
      "using base dense general eigen fallback; right residuals certified"
    },
    extras = list(generalized = FALSE)
  )
}

#' @keywords internal
eig_full_generalized <- function(A, B, structure, vectors, tol) {
  n <- nrow(A)
  hermitian_path <- identical(structure$kind, "hermitian")
  if (hermitian_path) {
    eig <- if (is.complex(A) || is.complex(B)) {
      native_dense_complex_generalized_hpd_eigen(A, B)
    } else {
      dense_generalized_spd_eigen(A, B)
    }
    cert <- eig_full_certificate(A, B, eig$values, eig$vectors, vectors, tol,
                                 general = FALSE)
    return(eig_full_result(
      values = eig$values,
      vectors = if (isTRUE(vectors)) eig$vectors else NULL,
      certificate = cert,
      method = native_dense_generalized_spd_full_label(),
      n = n,
      warnings = "using native dense generalized SPD/Hermitian LAPACK full decomposition",
      extras = list(
        generalized = TRUE,
        classification = rep("finite", n),
        finite = rep(TRUE, n),
        infinite = rep(FALSE, n),
        undefined = rep(FALSE, n)
      ),
      controls_extra = list(
        generalized = TRUE,
        target_family = "dense_full_generalized",
        promotion_gate = "dense_full_generalized:spd",
        dense_fallback_policy = "none; dense full surface is the native path",
        alpha_beta_semantics = "all finite; SPD/Hermitian-definite pencil"
      )
    ))
  }

  complex_pencil <- is.complex(A) || is.complex(B)
  eig <- if (complex_pencil) {
    native_dense_complex_generalized_pencil_eigen(A, B)
  } else {
    native_dense_generalized_pencil_eigen(A, B)
  }
  # Real pencils go through DGGEVX, which returns one-norms of the balanced
  # pencil (abnrm/bbnrm) for scale-aware classification plus rconde/rcondv
  # conditioning diagnostics. R's bundled LAPACK subset has no ZGGEVX, so
  # complex pencils use ZGGEV with input one-norms for the same
  # classification policy and no conditioning diagnostics.
  norm_A <- eig$abnrm %||% norm(A, type = "1")
  norm_B <- eig$bbnrm %||% norm(B, type = "1")
  pencil <- generalized_pencil_values(eig$alpha, eig$beta,
                                      norm_A = norm_A, norm_B = norm_B)
  vecs <- eig$vectors
  cert <- if (isTRUE(vectors)) {
    certify_dense_generalized_pencil(A, B, eig$alpha, eig$beta, vecs, tol = tol,
                                     norm_A = norm_A, norm_B = norm_B)
  } else {
    empty_certificate(tol, note = "vectors not returned; residual certificate not computed")
  }
  conditioning <- eig_full_pencil_conditioning(eig, complex_pencil)
  eig_full_result(
    values = pencil$values,
    vectors = if (isTRUE(vectors)) vecs else NULL,
    certificate = cert,
    method = native_dense_generalized_pencil_full_label(),
    n = n,
    warnings = "using native dense general pencil LAPACK full decomposition",
    extras = list(
      generalized = TRUE,
      alpha = eig$alpha,
      beta = eig$beta,
      classification = pencil$classification,
      finite = pencil$finite,
      infinite = pencil$infinite,
      undefined = pencil$undefined,
      classification_policy = list(
        policy = pencil$tolerance_policy,
        tolerance = pencil$tolerance,
        alpha_threshold = pencil$alpha_threshold,
        beta_threshold = pencil$beta_threshold,
        norm_A = pencil$norm_A,
        norm_B = pencil$norm_B
      ),
      left_vectors = if (isTRUE(vectors)) eig$left_vectors else NULL,
      conditioning = conditioning
    ),
    controls_extra = list(
      generalized = TRUE,
      target_family = "dense_full_generalized",
      promotion_gate = if (complex_pencil) {
        "dense_full_generalized:general_pencil_complex"
      } else {
        "dense_full_generalized:general_pencil_real"
      },
      dense_fallback_policy = "none; dense full surface is the native path",
      alpha_beta_semantics = paste(
        "homogeneous alpha/beta with", pencil$tolerance_policy,
        "classification"
      ),
      conditioning_diagnostics = isTRUE(conditioning$available)
    )
  )
}

#' @keywords internal
eig_full_pencil_conditioning <- function(eig, complex_pencil) {
  if (complex_pencil || is.null(eig$rconde)) {
    return(list(
      available = FALSE,
      source = "none",
      note = paste(
        "conditioning diagnostics require the expert LAPACK driver;",
        "R's bundled LAPACK subset provides DGGEVX for real pencils only"
      )
    ))
  }
  list(
    available = TRUE,
    source = "LAPACK dggevx (balanced, sense = 'B')",
    rconde = eig$rconde,
    rcondv = eig$rcondv,
    abnrm = eig$abnrm,
    bbnrm = eig$bbnrm
  )
}

#' @keywords internal
eig_full_result <- function(values, vectors, certificate, method, n, warnings,
                            extras = list(), controls_extra = list()) {
  problem <- list(type = "eigen", target = largest())
  plan <- new_plan(
    problem = problem,
    k = n,
    method = method,
    reasons = c("full dense eigendecomposition requested through eig_full()"),
    fallback = "none",
    controls = c(list(full = TRUE, dense = TRUE), controls_extra)
  )
  plan$target <- "all"
  make_eigen_result(
    values = values,
    vectors = vectors,
    certificate = certificate,
    iter = list(iterations = 1L, matvecs = 0L),
    requested = n,
    method_label = method,
    target_label_value = "all",
    plan = plan,
    warnings = warnings,
    extras = extras
  )
}

#' @keywords internal
eig_full_certificate <- function(A, B, values, vecs, vectors, tol, general) {
  if (!isTRUE(vectors)) {
    return(empty_certificate(tol, note = "vectors not returned; residual certificate not computed"))
  }
  if (isTRUE(general)) {
    certify_dense_general_eigen(A, values, vecs, tol = tol)
  } else {
    certify_eigen(A, values, vecs, B = B, tol = tol)
  }
}

#' @keywords internal
eig_full_dense_input <- function(x, name, allow_dense_fallback) {
  if (is.matrix(x)) {
    return(x)
  }
  stop(
    name, " must be a base dense matrix for eig_full(); sparse/operator ",
    "full decompositions are not silently densified",
    call. = FALSE
  )
}

#' @keywords internal
eig_full_validate_square <- function(A, B = NULL) {
  if (is.null(dim(A)) || nrow(A) != ncol(A)) {
    stop("A must be a square matrix.", call. = FALSE)
  }
  if (!is.null(B) && (is.null(dim(B)) || nrow(B) != ncol(B) ||
      nrow(B) != nrow(A))) {
    stop("B must be a square matrix with the same dimension as A.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
eig_full_resolve_structure <- function(A, B, structure, tol) {
  if (!is.null(structure)) {
    return(structure)
  }
  if (eig_full_is_hermitian(A, tol = tol) &&
      (is.null(B) || eig_full_is_hermitian(B, tol = tol))) {
    hermitian()
  } else {
    general()
  }
}

#' @keywords internal
eig_full_is_hermitian <- function(A, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    return(FALSE)
  }
  if (is.complex(A)) {
    scale <- max(1, Mod(A))
    return(max(Mod(A - Conj(t(A)))) <= tol * scale)
  }
  isTRUE(.Call("eigencore_dense_is_symmetric", as.matrix(A), as.numeric(tol),
               PACKAGE = "eigencore"))
}

#' @keywords internal
eig_full_as_complex_matrix <- function(A) {
  if (is.complex(A)) {
    return(A)
  }
  matrix(as.complex(A), nrow = nrow(A), ncol = ncol(A), dimnames = dimnames(A))
}
