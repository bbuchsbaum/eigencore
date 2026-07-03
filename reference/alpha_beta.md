# Extract homogeneous generalized coordinates.

Generalized dense-pencil and generalized Schur results store eigenvalues
in homogeneous form as `alpha / beta`; generalized SVD results use the
same alpha/beta convention for generalized singular values.
`alpha_beta()` exposes those coordinates along with the
finite/infinite/undefined classification computed by eigencore.

## Usage

``` r
alpha_beta(x, ...)
```

## Arguments

- x:

  An eigencore result with homogeneous `alpha` and `beta` fields.

- ...:

  Reserved for future methods.

## Value

A list containing `alpha`, `beta`, and any available `values`,
`classification`, `finite`, `infinite`, and `undefined` fields. Results
that record how the finite/infinite/undefined labels were decided also
include a `classification_policy` list with the policy name, the
tolerance, the per-coordinate zero thresholds, and the pencil norms used
for norm-scaled classification.

## Examples

``` r
A <- diag(c(2, 3, 0))
B <- diag(c(1, 0, 0))
fit <- eig_full(A, B = B, structure = general())
alpha_beta(fit)$classification
#> [1] "finite"    "infinite"  "undefined"
```
