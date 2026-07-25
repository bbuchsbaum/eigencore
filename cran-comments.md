# cran-comments

## Resubmission

This is a follow-up patch release (1.0.2) for the R CMD check ERROR on the
r-devel Linux fedora flavors
(<https://cran.r-project.org/web/checks/check_results_eigencore.html>).

The 1.0.1 fix was incomplete. It guarded the package's validation and installed
benchmark helpers, but overlooked a separate `bench::mark(memory = TRUE)` call
in `vignettes/benchmarks.Rmd`. On R builds without memory profiling, that call
recorded every benchmark method as an error; the timing summary then failed
while rebuilding the vignette.

Version 1.0.2 audits every shipped `bench::mark()` entry point. Memory
measurement is requested only when `capabilities("profmem")` is true, while
timing and numerical checks continue when it is false. The vignette summary
also handles a regime with no successful benchmark method without failing.
A regression test rejects unconditional `memory = TRUE` in the shipped
benchmark entry points.

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1,
  `R CMD check --as-cran --no-manual`
* Fedora/no-memory-profiling check: pending on this exact 1.0.2 tarball

## R CMD check results

0 errors | 0 warnings | 0 notes locally. Fedora is pending on the exact source
tarball.

## Downstream dependencies

None: no reverse dependencies are affected by this patch.
