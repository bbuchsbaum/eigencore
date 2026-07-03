# Orthogonalize columns in the B-inner product.

Orthogonalize columns in the B-inner product.

## Usage

``` r
b_orthogonalize(X, B, against = NULL, tol = sqrt(.Machine$double.eps))
```

## Arguments

- X:

  Numeric matrix whose columns are orthogonalized.

- B:

  Symmetric positive-definite metric matrix defining the inner product.

- against:

  Optional matrix whose columns are projected out in the B-inner
  product.

- tol:

  Orthogonality warning tolerance.
