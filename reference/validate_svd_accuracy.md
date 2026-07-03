# Validate eigencore SVD results against base::svd().

Validate eigencore SVD results against base::svd().

## Usage

``` r
validate_svd_accuracy(A, rank, target = largest(), fit = NULL, tol = 1e-08)
```

## Arguments

- A:

  Matrix or eigencore operator to validate.

- rank:

  Number of singular values to validate.

- target:

  Eigencore singular-value target descriptor.

- fit:

  Optional precomputed eigencore SVD result.

- tol:

  Validation tolerance.
