# LOBPCG method descriptor.

LOBPCG method descriptor.

## Usage

``` r
lobpcg(maxit = 200L, preconditioner = NULL, constraints = NULL)
```

## Arguments

- maxit:

  Maximum LOBPCG iterations.

- preconditioner:

  Optional function taking a residual block and returning a
  preconditioned block with the same dimensions.

- constraints:

  Optional matrix whose columns span a subspace to deflate. Iterates are
  kept orthogonal to this subspace in the Euclidean or generalized `B`
  inner product. Native constrained LOBPCG is not promoted yet;
  constrained problems use the honest reference path.

## Value

An `eigencore_method` descriptor selecting LOBPCG. Built-in standard
Hermitian dense/CSC operators may use a native prototype; unsupported
cases route to the reference prototype.
