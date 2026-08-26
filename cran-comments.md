# cran-comments

## Resubmission

This is a focused patch release (1.0.3) responding to the `clang23`
"Additional issue" reported at
<https://cran.r-project.org/web/checks/check_results_eigencore.html>.

LLVM 23's libc++ removed transitive standard-library includes. The package's
`src/native_operators.cpp` used `std::fill()` but did not include the declaring
`<algorithm>` header directly, so installation failed under LLVM 23.1.0.

Version 1.0.3 adds the direct `<algorithm>` include. It does not define the
temporary `_LIBCPP_KEEP_TRANSITIVE_INCLUDES_LLVM23` compatibility macro. A
regression test now audits every shipped C++ source and header that uses a
standard algorithm and requires it to include `<algorithm>` directly.

## Test environments

* local: macOS 14 (aarch64-apple-darwin20), R 4.5.1, Homebrew clang 20.1.8,
  exact-tarball `R CMD check --as-cran --no-manual`
* GitHub Actions: macOS release, Windows release, Ubuntu release,
  Ubuntu devel, and Ubuntu oldrel-1
* dedicated GitHub Actions check: Ubuntu 24.04 (x86_64), R-devel,
  LLVM 23.1.0, libc++ 23.1.0, and libc++abi 23.1.0

## R CMD check results

0 errors | 0 warnings | 0 notes locally.

The five-platform GitHub Actions matrix passed. The dedicated LLVM 23.1.0/
libc++ 23.1.0 job also passed a full package check with no errors or warnings.
Its compiler probe compiled the repaired `native_operators.cpp`. As a negative
control, the same source with only the new `<algorithm>` line removed
reproduced the reported error: `no member named 'fill' in namespace 'std'`.

## Reproducible verification evidence

* Exact source tarball SHA-256:
  `5e21bf5dd6a76b9cb5599cbbe784f229c91bbe5295280fa34a5267a0ec758fd8`.
* Tested package-source commit:
  `f940c55eb4e160ad780e830b6a28edcb757e61fb`.
* The full local exact-tarball check passed installation, incoming
  feasibility, examples, tests, and vignette rebuilding with `Status: OK`.
* LLVM 23/libc++ 23 check:
  <https://github.com/bbuchsbaum/eigencore/actions/runs/32973329914>.
* Standard five-platform check matrix:
  <https://github.com/bbuchsbaum/eigencore/actions/runs/32973329915>.

## Downstream dependencies

CRAN reports no reverse dependencies for eigencore as of 2026-08-26.
