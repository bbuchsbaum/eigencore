# cran-comments

## Release summary

This is a major update from 1.3.0 to 2.0.0.

* `certificate$passed` is now strict: it is `TRUE` only when residuals are
  certified and the returned set was verified to be the requested one.
* Backward errors use the standard 2-norm definition.
* Persisted plans from 1.x must be re-planned; they fail with a typed error,
  not a crash.

These behaviour changes are why the major version is bumped. New features:
eigenvalue counting, interval targets, complex Hermitian iterative routes,
and faster single- and multi-threaded kernels.

## Test environments

* local: Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran --no-manual`
* GitHub Actions: macOS R-release, Windows R-release, and Ubuntu R-devel,
  R-release and R-oldrel

## R CMD check results

0 errors | 0 warnings | 1 note

* installed package size: the compiled solver kernels (libs) are about
  7 MB on Linux.

## Downstream dependencies

None: no reverse dependencies are affected by this release.
