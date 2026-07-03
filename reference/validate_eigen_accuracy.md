# Validate eigencore eigen results against a dense oracle.

Validate eigencore eigen results against a dense oracle.

## Usage

``` r
validate_eigen_accuracy(
  A,
  k,
  target = largest(),
  B = NULL,
  fit = NULL,
  tol = 1e-08
)
```

## Arguments

- A:

  Matrix or eigencore operator to validate.

- k:

  Number of eigenpairs to validate.

- target:

  Eigencore eigenvalue target descriptor.

- B:

  Optional metric matrix or operator for generalized problems.

- fit:

  Optional precomputed eigencore eigen result.

- tol:

  Validation tolerance.
