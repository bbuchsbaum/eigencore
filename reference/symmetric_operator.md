# Mark an operator as symmetric/Hermitian.

Mark an operator as symmetric/Hermitian.

## Usage

``` r
symmetric_operator(A, validate = TRUE, tol = 1e-10)
```

## Arguments

- A:

  Operator-like object.

- validate:

  Whether to check the adjoint identity before marking the operator
  symmetric.

- tol:

  Relative tolerance for the adjoint check.

## Value

An `eigencore_operator` with Hermitian structure metadata.
