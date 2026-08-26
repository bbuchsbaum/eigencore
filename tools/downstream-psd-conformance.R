#!/usr/bin/env Rscript

# Installed-only downstream conformance gate for the eigencore 1.3 PSD API.
#
# This script deliberately uses exported package APIs only. Run it from an
# arbitrary directory with an isolated library containing exact source builds
# of eigencore, gprocrustes, and DKGE. Provenance values are required so a
# successful run is an artifact receipt rather than an untraceable smoke test.

options(warn = 2)

required_provenance <- c(
  EIGENCORE_SOURCE_SHA = "eigencore source commit",
  EIGENCORE_TARBALL_SHA256 = "eigencore source tarball",
  GPROCRUSTES_SOURCE_SHA = "gprocrustes source commit",
  GPROCRUSTES_TARBALL_SHA256 = "gprocrustes source tarball",
  DKGE_SOURCE_SHA = "DKGE source commit",
  DKGE_TARBALL_SHA256 = "DKGE source tarball"
)

provenance <- Sys.getenv(names(required_provenance), unset = NA_character_)
names(provenance) <- names(required_provenance)
missing_provenance <- names(provenance)[is.na(provenance) | !nzchar(provenance)]
if (length(missing_provenance)) {
  stop(
    "Missing exact-artifact provenance: ",
    paste(missing_provenance, collapse = ", "),
    call. = FALSE
  )
}

sha256_pattern <- "^[[:xdigit:]]{64}$"
sha256_names <- grep("SHA256$", names(provenance), value = TRUE)
if (any(!grepl(sha256_pattern, provenance[sha256_names]))) {
  stop("Every tarball SHA-256 must contain exactly 64 hexadecimal digits.", call. = FALSE)
}

expected_library <- Sys.getenv("EIGENCORE_CONFORMANCE_LIB", unset = NA_character_)
if (is.na(expected_library) || !nzchar(expected_library)) {
  stop("EIGENCORE_CONFORMANCE_LIB must name the isolated installed library.", call. = FALSE)
}
expected_library <- normalizePath(expected_library, mustWork = TRUE)

library(eigencore)
library(gprocrustes)
library(dkge)
library(Matrix)

package_names <- c("eigencore", "gprocrustes", "dkge")
package_paths <- vapply(package_names, find.package, character(1L))
package_paths <- normalizePath(package_paths, mustWork = TRUE)
if (any(dirname(package_paths) != expected_library)) {
  stop(
    "Every tested package must be loaded from EIGENCORE_CONFORMANCE_LIB; got: ",
    paste(paste(names(package_paths), package_paths, sep = "="), collapse = "; "),
    call. = FALSE
  )
}

checks <- 0L
check_true <- function(value, label) {
  checks <<- checks + 1L
  if (!isTRUE(value)) stop("Conformance failure: ", label, call. = FALSE)
  invisible(TRUE)
}

check_close <- function(actual, expected, tolerance, label) {
  checks <<- checks + 1L
  actual <- as.matrix(actual)
  expected <- as.matrix(expected)
  if (!identical(dim(actual), dim(expected))) {
    stop(
      "Conformance failure: ", label, " has dimensions ",
      paste(dim(actual), collapse = "x"), " rather than ",
      paste(dim(expected), collapse = "x"),
      call. = FALSE
    )
  }
  scale <- max(1, sqrt(sum(expected * expected)))
  defect <- sqrt(sum((actual - expected)^2))
  if (!is.finite(defect) || defect > tolerance * scale) {
    stop(
      "Conformance failure: ", label, " relative Frobenius defect ",
      format(defect / scale, digits = 6L), " exceeds ", tolerance,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

check_scalar_close <- function(actual, expected, tolerance, label) {
  checks <<- checks + 1L
  scale <- max(1, abs(expected))
  defect <- abs(actual - expected)
  if (length(actual) != 1L || !is.finite(defect) || defect > tolerance * scale) {
    stop(
      "Conformance failure: ", label, " relative defect ",
      format(defect / scale, digits = 6L), " exceeds ", tolerance,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

dependency_names <- function(package) {
  description <- utils::packageDescription(package)
  fields <- unlist(description[c("Depends", "Imports", "Suggests")], use.names = FALSE)
  if (!length(fields)) return(character())
  fields <- paste(fields, collapse = ",")
  names <- trimws(strsplit(fields, ",", fixed = TRUE)[[1L]])
  unique(sub("[[:space:]]*\\(.*$", "", names))
}

downstream_names <- c("gprocrustes", "dkge", "rfugw")
check_true(
  !any(downstream_names %in% dependency_names("eigencore")),
  "eigencore must not reverse the downstream dependency direction"
)
check_true(
  "eigencore" %in% dependency_names("gprocrustes"),
  "gprocrustes must declare its existing eigencore dependency"
)
check_true(
  !"eigencore" %in% dependency_names("dkge"),
  "DKGE remains an independent native differential oracle for this gate"
)

hadamard4 <- matrix(c(
  1, 1, 1, 1,
  1, -1, 1, -1,
  1, 1, -1, -1,
  1, -1, -1, 1
), 4L, 4L, byrow = TRUE) / 2

seed <- 20260825L
set.seed(seed)
eigenvalues <- c(9, 1, 0, 0)
K <- hadamard4 %*% diag(eigenvalues) %*% t(hadamard4)
factor <- psd_factor(K)
X <- matrix(c(1, 2, -1, 3, 0, 2, 4, -2), 4L, 2L)
Y <- matrix(c(2, -1, 3, 1, -2, 5, 1, 0), 4L, 2L)
P <- hadamard4[, 1:2, drop = FALSE] %*% t(hadamard4[, 1:2, drop = FALSE])
N <- diag(4L) - P
null_X <- N %*% matrix(rnorm(8L), 4L, 2L)
null_Y <- N %*% matrix(rnorm(8L), 4L, 2L)

# Generic exported API fixture: complete singular geometry, action identities,
# quotient reduction, block primitive, serialization, and operator exposure.
check_true(identical(psd_rank(factor), 2L), "generic numerical rank")
check_true(identical(psd_nullity(factor), 2L), "generic numerical nullity")
check_close(psd_apply(factor, X, "form"), K %*% X, 2e-12, "generic form action")
check_close(
  psd_apply(factor, X, "image_projector"), P %*% X, 2e-12,
  "generic image projector"
)
check_close(
  psd_apply(factor, X, "null_projector"), N %*% X, 2e-12,
  "generic null projector"
)
check_close(
  crossprod(psd_reduce(factor, X), psd_reduce(factor, Y)),
  crossprod(X, K %*% Y), 3e-12, "generic quotient Gram identity"
)
check_close(
  psd_reduce(factor, X + null_X), psd_reduce(factor, X), 3e-12,
  "generic null-space invariance"
)
block <- psd_orthonormalize(factor, X, required_rank = 2L)
check_true(isTRUE(certificate(block)$passed), "generic block certificate")
check_close(psd_gram(factor, block$basis), diag(2L), 3e-11, "generic block postcondition")
check_close(
  psd_operator(factor, "form")$apply(X), K %*% X, 2e-12,
  "generic exported operator action"
)
roundtrip_file <- tempfile(fileext = ".rds")
saveRDS(factor, roundtrip_file)
roundtrip <- readRDS(roundtrip_file)
unlink(roundtrip_file)
check_close(
  psd_apply(roundtrip, X, "form"), K %*% X, 2e-12,
  "generic installed serialization round trip"
)

# gprocrustes owns the polar solver and transform semantics. Eigencore supplies
# certified K validation and optional image reduction; the direct moment below
# remains the independent application oracle.
Xr <- psd_reduce(factor, X)
Yr <- psd_reduce(factor, Y)
C <- crossprod(X, K %*% Y)
sv <- base::svd(C)
R_o <- sv$u %*% t(sv$v)
gproc_o <- gprocrustes::procrustes(
  Xr, Yr, transform = gprocrustes::proc_orthogonal("O"), center = FALSE
)
check_close(gproc_o$moments$C, C, 3e-12, "gprocrustes direct K moment")
check_close(gproc_o$transform$R, R_o, 3e-12, "gprocrustes O action oracle")
check_scalar_close(
  gproc_o$objective, sum((Xr %*% R_o - Yr)^2), 3e-12,
  "gprocrustes O objective"
)

gproc_null <- gprocrustes::procrustes(
  psd_reduce(factor, X + null_X),
  psd_reduce(factor, Y + null_Y),
  transform = gprocrustes::proc_orthogonal("O"), center = FALSE
)
check_close(gproc_null$transform$R, gproc_o$transform$R, 3e-12, "gprocrustes null invariance")
check_scalar_close(gproc_null$objective, gproc_o$objective, 3e-12, "gprocrustes null objective")

permutation <- c(3L, 1L, 4L, 2L)
Kp <- K[permutation, permutation]
factor_p <- psd_factor(Kp)
gproc_p <- gprocrustes::procrustes(
  psd_reduce(factor_p, X[permutation, , drop = FALSE]),
  psd_reduce(factor_p, Y[permutation, , drop = FALSE]),
  transform = gprocrustes::proc_orthogonal("O"), center = FALSE
)
check_close(gproc_p$moments$C, C, 4e-12, "gprocrustes permutation-equivariant moment")
check_close(gproc_p$transform$R, gproc_o$transform$R, 4e-12, "gprocrustes permutation action")

diagonal_factor <- psd_factor(c(4, 1, 0))
diagonal_Xr <- diag(2L)
diagonal_Yr <- diag(c(2, -1), 2L)
diagonal_X <- psd_lift(diagonal_factor, diagonal_Xr)
diagonal_Y <- psd_lift(diagonal_factor, diagonal_Yr)
diagonal_C <- crossprod(diagonal_X, diag(c(4, 1, 0)) %*% diagonal_Y)
check_close(diagonal_C, diag(c(2, -1), 2L), 1e-14, "gprocrustes diagonal compatibility moment")

gproc_diag_o <- gprocrustes::procrustes(
  psd_reduce(diagonal_factor, diagonal_X),
  psd_reduce(diagonal_factor, diagonal_Y),
  transform = gprocrustes::proc_orthogonal("O"), center = FALSE
)
gproc_diag_so <- gprocrustes::procrustes(
  psd_reduce(diagonal_factor, diagonal_X),
  psd_reduce(diagonal_factor, diagonal_Y),
  transform = gprocrustes::proc_orthogonal("SO"), center = FALSE
)
manual_o <- diag(c(1, -1), 2L)
manual_so <- diag(2L)
check_close(gproc_diag_o$transform$R, manual_o, 1e-14, "gprocrustes diagonal O action")
check_close(gproc_diag_so$transform$R, manual_so, 1e-14, "gprocrustes proper SO action")
check_true(det(gproc_diag_o$transform$R) < 0, "gprocrustes O fixture distinguishes reflection")
check_scalar_close(det(gproc_diag_so$transform$R), 1, 1e-14, "gprocrustes SO determinant")
check_scalar_close(
  sum(gproc_diag_so$transform$R * diagonal_C), 1, 1e-14,
  "gprocrustes SO constrained optimum"
)

# DKGE remains a public native differential oracle. Compare mathematical
# values, projectors, and quotient coordinates—not arbitrary eigenvector signs.
dkge_roots <- dkge::dkge_kernel_roots(K, jitter = 0, tol = sqrt(.Machine$double.eps))
check_true(identical(as.integer(dkge_roots$rank), psd_rank(factor)), "DKGE root rank")
check_close(
  dkge_roots$Khalf, psd_apply(factor, diag(4L), "sqrt"), 3e-12,
  "DKGE square-root oracle"
)
check_close(
  dkge_roots$Kihalf, psd_apply(factor, diag(4L), "inverse_sqrt"), 3e-12,
  "DKGE inverse-root oracle"
)
check_close(
  dkge_roots$support_projector, psd_apply(factor, diag(4L), "image_projector"),
  3e-12, "DKGE support projector oracle"
)

W <- hadamard4[, 1:2, drop = FALSE] %*% matrix(c(1, 2, -1, 1), 2L, 2L)
eigencore_U <- psd_orthonormalize(factor, W, required_rank = 2L)$basis
dkge_U <- dkge::dkge_k_orthonormalize(W, K)
check_close(crossprod(dkge_U, K %*% dkge_U), diag(2L), 3e-11, "DKGE K-orthonormality")
check_close(
  tcrossprod(qr.Q(qr(dkge_U))), tcrossprod(qr.Q(qr(eigencore_U))), 3e-11,
  "DKGE/eigencore orthonormalized span"
)

theta <- 0.43
R_fixture <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2L, 2L)
dkge_Uref <- dkge_U
dkge_U2 <- dkge_U %*% R_fixture
dkge_alignment <- dkge::dkge_procrustes_K(dkge_Uref, dkge_U2, K)
dkge_C <- crossprod(dkge_Uref, K %*% dkge_U2)
dkge_sv <- base::svd(dkge_C)
dkge_R_oracle <- dkge_sv$v %*% t(dkge_sv$u)
check_close(dkge_alignment$R, dkge_R_oracle, 3e-12, "DKGE K-Procrustes action")
check_close(dkge_alignment$U_aligned, dkge_Uref, 3e-12, "DKGE aligned basis")
check_scalar_close(
  dkge_alignment$d, sum(diag(dkge_C %*% dkge_R_oracle)), 3e-12,
  "DKGE K-Procrustes objective"
)

null_basis <- hadamard4[, 3:4, drop = FALSE]
dkge_Uref_null <- dkge_Uref + null_basis %*% matrix(c(2, -1, 3, 4), 2L, 2L)
dkge_U2_null <- dkge_U2 + null_basis %*% matrix(c(-2, 5, 1, -3), 2L, 2L)
dkge_null_alignment <- dkge::dkge_procrustes_K(dkge_Uref_null, dkge_U2_null, K)
check_close(
  dkge_null_alignment$R, dkge_alignment$R, 3e-12,
  "DKGE K-Procrustes null invariance"
)
check_close(
  psd_reduce(factor, dkge_Uref_null), psd_reduce(factor, dkge_Uref), 3e-12,
  "DKGE reference quotient invariance"
)

projector_native <- dkge::dkge_projector_K(W, K)
gram_W <- crossprod(W, K %*% W)
projector_K_oracle <- W %*% solve(gram_W) %*% crossprod(W, K)
root_K <- psd_apply(factor, diag(4L), "sqrt")
projector_hat_oracle <- root_K %*% W %*% solve(gram_W) %*% t(W) %*% root_K
check_close(projector_native$P_K, projector_K_oracle, 3e-11, "DKGE effect-space projector")
check_close(projector_native$P_hat, projector_hat_oracle, 3e-11, "DKGE root-space projector")

dkge_scale_divergence <- NULL
for (scale in c(1e-12, 1, 1e12)) {
  scaled_K <- scale * K
  scaled_factor <- psd_factor(scaled_K)
  scaled_eigencore <- psd_orthonormalize(scaled_factor, W, required_rank = 2L)$basis
  check_true(
    identical(psd_rank(scaled_factor), 2L),
    paste("scale-invariant eigencore rank at", format(scale, scientific = TRUE))
  )
  scaled_native <- tryCatch(
    dkge::dkge_k_orthonormalize(W, scaled_K),
    error = identity
  )
  if (identical(scale, 1e-12)) {
    check_true(
      inherits(scaled_native, "error") &&
        grepl("full column rank in the K metric", conditionMessage(scaled_native), fixed = TRUE),
      "DKGE tiny-scale differential remains explicit and diagnostically stable"
    )
    dkge_scale_divergence <- conditionMessage(scaled_native)
    next
  }
  check_true(!inherits(scaled_native, "error"), paste("DKGE admits scale", scale))
  check_close(
    crossprod(scaled_native, scaled_K %*% scaled_native), diag(2L), 4e-10,
    paste("DKGE scaled K-orthonormality at", format(scale, scientific = TRUE))
  )
  check_close(
    tcrossprod(qr.Q(qr(scaled_native))), tcrossprod(qr.Q(qr(scaled_eigencore))),
    4e-10, paste("DKGE/eigencore scaled span at", format(scale, scientific = TRUE))
  )
}

cat("eigencore 1.3 installed downstream PSD conformance: PASS\n")
cat("checks:", checks, "\n")
cat("seed:", seed, "\n")
cat("R:", as.character(getRversion()), "\n")
cat(
  "packages:",
  paste(
    paste0(
      package_names, "=",
      vapply(package_names, function(package) as.character(packageVersion(package)), character(1L))
    ),
    collapse = "; "
  ),
  "\n"
)
cat("library:", expected_library, "\n")
cat("DKGE_TINY_SCALE_DIFFERENTIAL=", dkge_scale_divergence, "\n", sep = "")
for (name in names(provenance)) cat(name, "=", provenance[[name]], "\n", sep = "")
