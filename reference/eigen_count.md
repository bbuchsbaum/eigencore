# Count eigenvalues below, at and above a shift (Sylvester inertia)

`eigen_count()` counts the eigenvalues of a real symmetric or complex
Hermitian matrix `A` (or of the symmetric-definite pencil `(A, B)`) that
lie below, at, and above `sigma`, without computing any eigenvalue. It
factors `A - sigma * B` as `P L D L' P'` and reads the counts off the
signs of the pivots of `D` (Sylvester's law of inertia).

## Usage

``` r
eigen_count(A, sigma, B = NULL, perturb = TRUE, pivot_tol = NULL)
```

## Arguments

- A:

  A real symmetric or complex Hermitian matrix: base `matrix`, a
  `Matrix` object (dense, sparse or diagonal), or an
  `eigencore_operator` with an explicit matrix source. The input must be
  symmetric (checked); only one triangle enters the factorisation.

- sigma:

  A single finite shift.

- B:

  Optional symmetric positive definite matrix for the pencil
  `A x = lambda B x` (checked). Counts are then of the generalized
  eigenvalues.

- perturb:

  Logical; retry an unreliable shift at nearby shifts.

- pivot_tol:

  Optional extra relative pivot floor: a count is unreliable when the
  smallest pivot is below `pivot_tol * scale`. Default
  `getOption("eigencore.inertia_pivot_tol", 0)` (the backward-error
  bound alone decides).

## Value

An object of class `eigencore_inertia`, a list with

- `sigma`: the requested shift; `n`: the dimension;

- `below`, `zero`, `above`: eigenvalue counts below `sigma` (below
  `sigma - delta` when perturbed), inside the zero band, and above;

- `zero_band`: `c(sigma - delta, sigma + delta)` (`delta = 0` unless
  perturbed); `perturbed`, `perturbation` (`delta`);

- `reliable`: whether the counts are trustworthy (see Details);

- `method`: `"diagonal"`, `"tridiagonal_sturm"`,
  `"dense_bunch_kaufman"`, `"dense_bunch_kaufman_hermitian_embedding"`
  or `"sparse_cholmod_ldl"`; `generalized`: whether `B` was given;

- `diagnostics`: a list with `scale` (an upper bound of
  `||A - sigma B||_2`), `min_pivot` and `max_pivot` (absolute),
  `min_pivot_relative` (`min_pivot / scale`), `pivot_growth`,
  `backward_bound`, `factor_nnz` (sparse), `attempts` (one row per
  factorisation: shift, counts, min pivot, reliable, error) and
  `seconds`;

- `notes`: character vector of explanations.

## Details

Factorisations by storage:

- diagonal `A` (and diagonal `B`): exact comparison;

- tridiagonal `A` with identity or diagonal `B`: Sturm sequence
  (backward stable for counting);

- dense real symmetric: LAPACK Bunch-Kaufman `dsytrf`; complex
  Hermitian: the same on the real symmetric embedding
  `[Re(A) -Im(A); Im(A) Re(A)]`;

- sparse (`dsCMatrix`, symmetric `dgCMatrix`, other sparse classes):
  CHOLMOD simplicial \\LDL^T\\ with an AMD fill-reducing ordering via
  [`Matrix::Cholesky()`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html).

**Reliability.** The computed inertia is exact for a nearby matrix
`A - sigma B + E`. `eigen_count()` bounds `||E||` by
`c * eps * || |L| |D| |L'| ||` (`backward_bound`) and declares a count
`reliable` only when every pivot is larger than that bound, the pivot
growth `|| |L| |D| |L'| || / ||A - sigma B||` is below `1 / sqrt(eps)`
and the factorisation completed. CHOLMOD's \\LDL^T\\ does not pivot, so
a shift near an eigenvalue (or an unlucky ordering) can produce a tiny
pivot. With `perturb = TRUE` an unreliable shift is retried at
`sigma - delta` and `sigma + delta` for increasing `delta` (relative to
the matrix scale); the counts then refer to the half-lines below
`sigma - delta` and above `sigma + delta`, and `zero` counts the
eigenvalues inside the band `[sigma - delta, sigma + delta]`
(`zero_band`). `reliable = FALSE` is reported whenever no attempt was
trustworthy; such counts must not be used as a certificate.
`reliable = TRUE` means no warning sign was found; it is a numerical
diagnosis, not an interval-arithmetic proof.

## See also

[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
uses these counts to certify that a computed set of extreme eigenvalues
is complete
(`certificate(fit)$target_completeness == "inertia_verified"`); see the
"Certificates" vignette.

## Examples

``` r
A <- diag(c(1, 2, 2, 5, 7))
eigen_count(A, 3)
#> <eigencore_inertia> dense_bunch_kaufman  
#>   sigma: 3 
#>   below: 3  zero: 0  above: 2  (n = 5 )
#>   reliable: TRUE  min pivot / scale: 0.25  pivot growth: 1 
S <- Matrix::sparseMatrix(i = c(1:6, 1:5), j = c(1:6, 2:6),
                          x = c(rep(2, 6), rep(-1, 5)), symmetric = TRUE)
eigen_count(S, 1)
#> <eigencore_inertia> tridiagonal_sturm  
#>   sigma: 1  zero band: sigma -/+ 5e-12 
#>   below: 2  zero: 0  above: 4  (n = 6 )
#>   reliable: TRUE  min pivot / scale: 2e-12  pivot growth: 1 
#>   note: sigma was within the numerical zero band of the factorisation; counts use sigma -/+ 5e-12 and 'zero' counts the eigenvalues inside that band 
sum(eigen(as.matrix(S), only.values = TRUE)$values < 1)
#> [1] 2
```
