# Changelog

## eigencore 1.0.0

First CRAN release.

- [`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md)
  and
  [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
  compute the top-*k* singular triplets or eigenpairs of large dense,
  sparse (CSC), diagonal, banded/tridiagonal, and matrix-free operators
  through native C++ kernels.
- Every result carries a numerical certificate: residuals for both
  singular relations, a backward-error bound, orthogonality loss, a
  labeled norm bound, and a single `passed` flag. Bounds that can only
  be estimated (for example stochastic norm estimates on centered sparse
  operators) are reported as estimates and never produce an unqualified
  `passed`.
- Operator algebra —
  [`center()`](https://bbuchsbaum.github.io/eigencore/reference/center.md),
  [`scale_cols()`](https://bbuchsbaum.github.io/eigencore/reference/scale_cols.md),
  [`compose()`](https://bbuchsbaum.github.io/eigencore/reference/compose.md),
  [`crossprod_operator()`](https://bbuchsbaum.github.io/eigencore/reference/crossprod_operator.md),
  [`linear_operator()`](https://bbuchsbaum.github.io/eigencore/reference/linear_operator.md)
  — solves centered, scaled, and composed problems without forming dense
  matrices.
- Transparent method selection:
  [`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md)
  reports the chosen kernel before a solve, and `fit$method` names the
  path that actually ran. Problem classes without a production kernel
  carry explicit `reference` labels.
- RSpectra-compatible wrappers
  [`eigs()`](https://bbuchsbaum.github.io/eigencore/reference/eigs.md),
  [`eigs_sym()`](https://bbuchsbaum.github.io/eigencore/reference/eigs_sym.md),
  and
  [`svds()`](https://bbuchsbaum.github.io/eigencore/reference/svds.md)
  accept the same `which` codes and additionally return certificates.
- Benchmarked against ‘RSpectra’, ‘irlba’, and ‘PRIMME’; reproduce with
  `Rscript inst/benchmarks/bench-readme.R`.
- Generalized eigen support:
  [`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md)
  for dense SPD and general pencils,
  [`generalized_schur()`](https://bbuchsbaum.github.io/eigencore/reference/generalized_schur.md)
  and
  [`generalized_svd()`](https://bbuchsbaum.github.io/eigencore/reference/generalized_svd.md)
  for dense QZ/GSVD, partial sparse general pencils with nonsingular
  diagonal `B` via transformed native Arnoldi, left eigenvectors and
  conditioning diagnostics on supported dense paths, and
  `pencil_norm_scaled` alpha/beta classification. Sparse SPD partial
  paths remain under
  [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
  / LOBPCG / B-orthogonal Lanczos; general sparse QZ and non-diagonal
  sparse `B` are explicit unsupported boundaries. The real dense GSVD
  path currently requires a linked LAPACK that provides the deprecated
  `dggsvd` routine.
- The exported API is stable as of 1.0.0; breaking changes from here
  follow semantic versioning.
