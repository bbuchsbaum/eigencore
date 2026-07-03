# cran-comments

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1
* GitHub Actions, commit 088d3bd, 2026-07-03: R-CMD-check passed on
  macOS release, Windows release, Ubuntu devel, Ubuntu release, and Ubuntu
  oldrel-1; lint, coverage, and pkgdown workflows also passed.
* R-hub, run 28686224852, 2026-07-03: Ubuntu release and Windows completed
  successfully.
* Pending before submission (tracked in bd-01KTV9C30CZ0XTZHEHVNVHPB1H):
  win-builder R-release and R-devel were submitted from commit 088d3bd with
  `--no-manual`; maintainer email results must be confirmed before CRAN
  submission. These checks are important because the package compiles ~17k
  lines of C++17 and x86 runs exercise 80-bit `long double` accumulation paths
  that do not exist on the arm64 development machine.

## R CMD check results

0 errors | 0 warnings | 1 note

* "New submission" — this is the first CRAN release of eigencore.

(A local-only "unable to verify current time" NOTE appears on the development
machine and is an environment artifact; it does not occur on CRAN
infrastructure.)

## Downstream dependencies

None: this is a new package with no reverse dependencies.
