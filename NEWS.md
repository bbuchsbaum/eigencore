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
* The exported API is stable as of 1.0.0; breaking changes from here follow
  semantic versioning.
