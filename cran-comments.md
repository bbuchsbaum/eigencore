# cran-comments

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1
* TODO before submission (tracked in bd-01KTV9C30CZ0XTZHEHVNVHPB1H): win-builder
  (release + devel), R-hub Linux x86_64 and Windows. The package compiles
  ~17k lines of C++17; x86 runs exercise 80-bit `long double` accumulation
  paths that do not exist on the arm64 development machine.

## R CMD check results

0 errors | 0 warnings | 1 note

* "New submission" — this is the first CRAN release of eigencore.

(A local-only "unable to verify current time" NOTE appears on the development
machine and is an environment artifact; it does not occur on CRAN
infrastructure.)

## Downstream dependencies

None: this is a new package with no reverse dependencies.
