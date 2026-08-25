# cran-comments

## Release summary

This is a minor update from 1.0.2 to 1.1.0. It adds certified warm starts for
the supported real Hermitian Lanczos paths, an optional mid-sweep convergence
check, and a native implicit-Gram route for eligible largest-target partial
SVD problems. Every new iterative route retains original-coordinate residual
and backward-error certification. The release also includes the
no-memory-profiling portability fixes from 1.0.1 and 1.0.2.

The public changes are additive. Existing calls with the new arguments omitted
retain their previous behavior, and exported symbols are checked against the
release snapshot.

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1,
  `R CMD check --as-cran --no-manual`
* GitHub Actions release matrix: macOS R-release, Windows R-release, and Ubuntu
  R-devel, R-release, and R-oldrel

## R CMD check results

0 errors | 0 warnings | 0 notes locally.

The exact source artifact, SHA-256 digest, and hosted-matrix receipts are
recorded with the 1.1.0 release rather than embedded into the source tarball.

## Downstream dependencies

None: no reverse dependencies are affected by this release.
