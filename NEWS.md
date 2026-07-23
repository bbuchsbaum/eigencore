# eigencore (development version)

## Performance

* New production `auto` route for largest-target partial SVD: a native
  implicit normal-equations (Gram) thick-restart Lanczos that runs on
  \eqn{A^T A} or \eqn{A A^T} as an operator, without materializing the Gram
  matrix. It covers dense matrices and sparse operators whose smaller side
  exceeds the explicit-Gram caps — regimes that previously fell back to a
  full LAPACK SVD (dense) or a single unrestarted Golub-Kahan sweep
  (sparse). Representative speedups at `tol = 1e-8` with certificates still
  passing: dense 4000x1000 `k = 10` ~18x, dense 2000x2000 `k = 10` ~80x,
  sparse 20000x5000 `k = 10` ~2x, sparse 20000x5000 `k = 50` ~2x. Results
  keep the exact two-sided residual certificate in original coordinates,
  and uncertified results still fall back to the native Golub-Kahan path.
* The scalar Golub-Kahan kernel now uses BLAS (dgemv) classical
  Gram-Schmidt reorthogonalization on the sparse CSC, matrix-free, and
  retained-restart paths, matching the dense path; previously these used
  scalar loops.
* The projected Golub-Kahan convergence check now runs `dbdsqr` directly on
  the projected bidiagonal, tracking only the last row of the left singular
  vectors, instead of a dense `dgesvd` with full vectors; scratch buffers
  are reused across checks. The check drops from O(iter^3) plus per-check
  allocations to O(iter^2).
* The scalar Lanczos convergence estimate (used by the shift-invert paths)
  computes eigenvalues with `dsterf` and recovers only the selected
  eigenvectors with `dstevr` instead of a full `dstev` decomposition each
  iteration, with scratch reuse across iterations.
* The CSC sparse matrix-multiply kernel now processes wide blocks in
  cache-friendly chunks of 10 columns instead of a strided generic loop.
* `as_operator()` no longer forces a full copy of an already-double dense
  matrix (`storage.mode<-` is now conditional), removing an O(m*n) copy
  plus GC churn from every problem construction — operator construction on
  a 2000x2000 input drops from ~52ms to ~6ms.
* The scalar thick-restart Hermitian Lanczos default subspace grows from
  `3k+20` to `max(3k+20, 5k)`: unchanged for small `k`, and at `k = 30`
  on a general sparse operator it cuts operator applications by ~25%.

## Portability

* The package again compiles on R < 4.4: guarded compatibility typedefs for
  `La_INT`/`La_LGL` and declarations for the complex QZ drivers
  (`zggev`, `zgges`) were added for headers that predate them, and
  `crossprod`/`tcrossprod` are now imported from Matrix so sparse-matrix
  dispatch does not rely on the base generics added in R 4.4.

# eigencore 1.0.0

First CRAN release.

* `svd_partial()` and `eig_partial()` compute the top-*k* singular triplets or
  eigenpairs of large dense, sparse (CSC), diagonal, banded/tridiagonal, and
  matrix-free operators through native C++ kernels.
* Every result carries a numerical certificate: residuals for both singular
  relations, a backward-error bound, orthogonality loss, a labeled norm bound,
  and a single `passed` flag. Bounds that can only be estimated (for example
  stochastic norm estimates on centered sparse operators) are reported as
  estimates and never produce an unqualified `passed`.
* Operator algebra — `center()`, `scale_cols()`, `compose()`,
  `crossprod_operator()`, `linear_operator()` — solves centered, scaled, and
  composed problems without forming dense matrices.
* Transparent method selection: `plan_solver()` reports the chosen kernel
  before a solve, and `fit$method` names the path that actually ran. Problem
  classes without a production kernel carry explicit `reference` labels.
* RSpectra-compatible wrappers `eigs()`, `eigs_sym()`, and `svds()` accept the
  same `which` codes and additionally return certificates.
* Benchmarked against 'RSpectra', 'irlba', and 'PRIMME'; reproduce with
  `Rscript inst/benchmarks/bench-readme.R`.
* Generalized eigen support: `eig_full()` for dense SPD and general pencils,
  `generalized_schur()` and `generalized_svd()` for dense QZ/GSVD, partial
  sparse general pencils with nonsingular diagonal `B` via transformed native
  Arnoldi, left eigenvectors and conditioning diagnostics on supported dense
  paths, and `pencil_norm_scaled` alpha/beta classification. Sparse SPD partial
  paths remain under `eig_partial()` / LOBPCG / B-orthogonal Lanczos; general
  sparse QZ and non-diagonal sparse `B` are explicit unsupported boundaries.
  The real dense GSVD path currently requires a linked LAPACK that provides the
  deprecated `dggsvd` routine.
* The exported API is stable as of 1.0.0; breaking changes from here follow
  semantic versioning.
