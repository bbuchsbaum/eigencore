# eigencore 1.0.2

* Complete the no-memory-profiling fix from 1.0.1. The CRAN benchmark
  vignette and the installed README benchmark now also request
  `bench::mark()` memory measurements only when `capabilities("profmem")` is
  true. The vignette's summary tables also report unavailable results instead
  of failing when a benchmark regime has no successful methods.
* Add a regression test that audits every shipped benchmark entry point for
  unconditional `memory = TRUE` and exercises the package timing helper on the
  current R build.

# eigencore 1.0.1

* Fix an R CMD check ERROR on R builds without memory profiling (for example
  the r-devel Linux fedora flavors, which are configured without
  `--enable-memory-profiling`). The internal benchmark/validation timing
  helpers passed `memory = TRUE` to `bench::mark()`, which calls
  `utils::Rprofmem()` and aborts with "memory profiling is not available on
  this system" on those platforms. Memory measurement is now requested only
  when `capabilities("profmem")` is `TRUE`; elsewhere timing still runs and
  `mem_alloc` is reported as `NA`.
* Benchmark smoke tests now use `testthat::skip_on_cran()` so they are skipped
  on CRAN as intended (the previous `Sys.getenv("CRAN")` guard never fired on
  the CRAN check farm).

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
