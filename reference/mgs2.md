# Modified Gram-Schmidt with two passes.

Modified Gram-Schmidt with two passes.

## Usage

``` r
mgs2(X, against = NULL, tol = sqrt(.Machine$double.eps))
```

## Arguments

- X:

  Numeric matrix whose columns are orthogonalized.

- against:

  Optional matrix whose columns are projected out before
  orthogonalization.

- tol:

  Orthogonality warning tolerance.
