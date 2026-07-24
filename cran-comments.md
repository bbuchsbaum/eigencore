# cran-comments

## Resubmission

This is a patch release (1.0.1) that fixes the R CMD check ERROR reported for
1.0.0 on the r-devel Linux fedora flavors
(<https://cran.r-project.org/web/checks/check_results_eigencore.html>).

Those check machines are configured without `--enable-memory-profiling`. The
package's internal benchmark/validation timing helpers requested memory
measurement from `bench::mark(memory = TRUE)`, which calls `utils::Rprofmem()`
and errors with "memory profiling is not available on this system" on such
builds. Memory measurement is now requested only when
`capabilities("profmem")` is `TRUE`; timing still runs otherwise and
`mem_alloc` is reported as `NA`. The affected benchmark smoke tests also now
use `testthat::skip_on_cran()`.

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1
* win-builder and R-hub (including a configuration without
  `--enable-memory-profiling`) will be re-run on this candidate before
  submission.

## R CMD check results

0 errors | 0 warnings | 2 notes

* "Days since last update: 1" — this 1.0.1 patch is submitted promptly to fix
  the R CMD check ERROR that CRAN reported for 1.0.0 (see Resubmission above).
* "unable to verify current time" — a local environment artifact on the
  development machine (no network time source); it does not occur on CRAN
  infrastructure.

## Downstream dependencies

None: no reverse dependencies are affected by this change.
