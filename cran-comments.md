# cran-comments

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1
* Fresh GitHub Actions, R-hub, and win-builder checks will be run on this exact
  release candidate before submission. Earlier cross-platform results do not
  sign off this candidate because pre-CRAN solver changes altered the package
  payload.

## R CMD check results

0 errors | 0 warnings | 1 note

* "New submission" — this is the first CRAN release of eigencore.

The local result is from a fresh `eigencore_1.0.0.tar.gz` built and checked
with `R CMD check --as-cran --no-manual` on 2026-07-15.

(A local-only "unable to verify current time" NOTE appears on the development
machine and is an environment artifact; it does not occur on CRAN
infrastructure.)

## Downstream dependencies

None: this is a new package with no reverse dependencies.
