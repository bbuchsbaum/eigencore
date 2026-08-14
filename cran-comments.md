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
* R-hub `gcc16`: Fedora Linux 44, R-devel (2026-07-24 r90297), GCC 16
* On the same Fedora runner, an explicit no-memory-profiling probe forced
  `capabilities("profmem")` to false, exercised both timing helpers and the
  regression test, and rebuilt `vignettes/benchmarks.Rmd`.

## R CMD check results

0 errors | 0 warnings | 0 notes locally.

The Fedora R-hub job completed successfully (3,806 test expectations passed;
47 environment/CRAN skips). The explicit no-memory-profiling probe also
completed successfully, including the benchmark-vignette render.

## Reproducible verification evidence

* Release payload: commit
  `cd32f89a577ef68ceb673da4a70953b74de11d2b`; built tarball SHA-256
  `0264d972b494775a6e1eb7c357b120f69e48016b7b68b21f74b3acbaaf51a2cb`.
* Public Fedora R-hub run:
  <https://github.com/bbuchsbaum/eigencore/actions/runs/30157914407>.
  Its `gcc16` job completed successfully:
  <https://github.com/bbuchsbaum/eigencore/actions/runs/30157914407/job/89678675741>.
  Both the `Exercise no-memory-profiling paths` step and the full R-hub check
  are green.
* The no-memory-profiling step forced `capabilities("profmem")` to false,
  asserted that the package and installed benchmark helpers returned finite
  timings with unavailable allocation recorded as `NA`, ran the portability
  regression test, and rebuilt `vignettes/benchmarks.Rmd`.

The R-hub image itself supports memory profiling. The explicit probe therefore
tests the exact false-capability branch used on CRAN builds without profiling;
it is not presented as an R binary compiled without that optional feature.
The checked CI commit `b50277deb51e357938a8c0f48af50958688d081c` differs
from the release payload only by `.github/scripts/check-no-profmem.R` and
`.github/workflows/rhub.yaml`; `.github` is excluded by `.Rbuildignore`.

## Downstream dependencies

None: no reverse dependencies are affected by this patch.
