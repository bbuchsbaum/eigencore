
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# eigencore

**eigencore** is a native engine for certified partial eigenvalue and
singular-value computation in R. The goal is a partial spectrum that is
**fast enough to use in production and trustworthy enough to build on**
— every returned result carries residuals, backward-error estimates,
orthogonality loss, and an inspectable solver plan that records exactly
which kernel ran. The package ships RSpectra-compatible shims (`eigs()`,
`eigs_sym()`, `svds()`) so existing code can migrate gradually.

## Installation

eigencore is in early development. You can install the development
version from GitHub:

``` r
# install.packages("pak")
pak::pak("bbuchsbaum/eigencore")
```

## Quick start

``` r
library(eigencore)

set.seed(1)
n <- 200
A <- crossprod(matrix(rnorm(n * n), n, n)) / n + diag(n)

fit <- eig_partial(A, k = 5, target = largest())
fit
#> Partial eigen decomposition
#>   requested: 5
#>   converged: 5
#>   method: native scalar thick-restart Hermitian Lanczos
#>   target: largest
#>   restart:thick_restart(in_native_loop)
#>   locked: 5
#>   max residual: 1.67046e-07
#>   max backward error: 4.546939e-09
#>   max orthogonality loss: 1.776357e-15
#>   norm bound: frobenius_exact+identity_exact
#>   scale estimated: FALSE
#>   certificate: passed
```

The result prints a one-screen ledger: which solver actually ran, how
many pairs converged, the worst residual and backward error across the
returned basis, and whether the certificate passed.

``` r
fit$certificate
#> eigencore certificate
#>   passed: TRUE
#>   tolerance: 1e-08
#>   type: residual_backward_error
#>   norm bound: frobenius_exact+identity_exact
#>   scale estimated: FALSE
#>   max residual: 1.67046e-07
#>   max backward error: 4.546939e-09
#>   max orthogonality loss: 1.776357e-15
#>   orthogonality tolerance: 1.490116e-08
#>   orthogonality required: TRUE
```

The certificate is the differentiator — eigencore tells you whether the
numbers are trustworthy, not just what they are.

A partial decomposition is a certified *slice* of the spectrum: you
compute the part you need and leave the rest untouched.

<img src="man/figures/README-spectrum-1.png" alt="The five largest eigenvalues highlighted in blue against the full 200-point spectrum of A in grey." width="100%" />

## RSpectra-compatible drop-in

If you have existing code calling `RSpectra::eigs_sym()`, swapping in
`eigencore::eigs_sym()` is a one-line change and gives you a certificate
for free:

``` r
res <- eigs_sym(A, k = 5, which = "LA")
res$values
#> [1] 5.010576 4.769169 4.700866 4.566055 4.504502
res$certificate
#> eigencore certificate
#>   passed: TRUE
#>   tolerance: 1e-08
#>   type: residual_backward_error
#>   norm bound: frobenius_exact+identity_exact
#>   scale estimated: FALSE
#>   max residual: 4.542481e-08
#>   max backward error: 1.236449e-09
#>   max orthogonality loss: 1.290817e-15
#>   orthogonality tolerance: 1.490116e-08
#>   orthogonality required: TRUE
```

`eigs()`, `eigs_sym()`, and `svds()` accept the same `which` codes as
`RSpectra` (`"LM"`, `"SM"`, `"LA"`, `"SA"`, `"LR"`, `"SR"`, `"LI"`,
`"SI"`, `"BE"`).

## Partial SVD

``` r
M <- matrix(rnorm(400 * 50), 400, 50)
svd_fit <- svd_partial(M, rank = 5, target = largest())
svd_fit
#> Partial SVD
#>   requested rank: 5
#>   converged rank: 5
#>   method: native certified Gram SVD special case
#>   target: largest
#>   max residual: 1.22125e-15
#>   max backward error: 8.568377e-18
#>   max orthogonality loss: 5.551115e-16
#>   norm bound: frobenius_exact
#>   scale estimated: FALSE
#>   certificate: passed
```

The `method` field reports which path ran. eigencore takes a fast
certified Gram special case when it can prove the result to tolerance,
and falls back to a Golub–Kahan bidiagonalization when the Gram route
can’t be certified. Either way the certificate covers *both* singular
relations, `||A v - sigma u||` and `||A^T u - sigma v||`, so the
triplets come back with honest backward-error bounds — eigencore never
takes a shortcut it can’t certify.

## What makes eigencore different

- **Certificates by default.** Every result carries `max_residual`,
  `max_backward_error`, `max_orthogonality_loss`, and a
  `passed`/`failed` verdict. The returned norm bound is labeled
  (`frobenius_exact`, `frobenius_hutchinson_estimate`, …) so you know
  how the certificate was computed.
- **Planner honesty.** The `method` field on every result identifies the
  path that actually ran —
  `native scalar thick-restart Hermitian Lanczos`,
  `native certified Gram SVD special case`,
  `native dense LAPACK SVD fallback`, or, when a problem class has no
  production kernel yet,
  `reference Hermitian Lanczos (prototype/oracle fallback)`. Sparse
  inputs do not silently densify; reference oracles are never dressed as
  production paths. Call `plan_solver()` to see the chosen path *before*
  you solve.
- **Block-native engine.** Operators apply to dense blocks through a
  frozen C++17 ABI. Dense, CSC, diagonal, scaled, centered, and composed
  operators go through native kernels without R-level iteration in the
  hot loop.
- **First-class operator algebra.** Build operators with
  `linear_operator()`, combine them with `compose()`,
  `crossprod_operator()`, `scale_cols()`, `center()`, etc. — the planner
  inspects the resulting structure and routes accordingly.

Run `vignette("eigencore", package = "eigencore")` for a guided tour of
operators, problems, plans, and certificates, and
`vignette("certificates", package = "eigencore")` for the deep dive on
reading the numerical evidence — including what to do when a check
fails.

## Status

**eigencore is experimental.** The public surface is settling but now
has a scoped V1 release surface with fresh local gate evidence. The V1
native engine ships scalar Hermitian Lanczos plus explicit and
diagnostic block Hermitian Lanczos, native structured-tridiagonal
Hermitian defaults, scoped sparse generalized LOBPCG, certified
tall/wide sparse SVD special cases, scoped randomized SVD acceleration,
dense LAPACK fallbacks, dense/structured shift-invert, and
dense/sparse-CSC nonsymmetric Arnoldi compatibility. Broader general
sparse and matrix-free SVD, fully native randomized solver control,
general sparse native LU, matrix-free native Arnoldi restart, and
native/block generalized Lanczos remain future scope. See
[`plan_v1.md`](plan_v1.md) and the docs under `docs/` for the
engineering roadmap and current limitations. Start with
[`docs/method-selection-and-workflows.md`](docs/method-selection-and-workflows.md)
for the current API workflow map and
[`docs/v1-benchmark-manifest.md`](docs/v1-benchmark-manifest.md) for the
benchmark gate inventory.

## License

MIT © Bradley Buchsbaum.
