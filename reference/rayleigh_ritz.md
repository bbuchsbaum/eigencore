# Rayleigh-Ritz projection in a trial basis.

Rayleigh-Ritz projection in a trial basis.

## Usage

``` r
rayleigh_ritz(A, Q, B = NULL, target = largest(), symmetric = TRUE)
```

## Arguments

- A:

  Dense matrix defining the projected eigenproblem.

- Q:

  Trial-basis matrix with basis vectors in columns.

- B:

  Optional symmetric positive-definite metric matrix.

- target:

  Eigencore target descriptor.

- symmetric:

  Whether the projected problem should be treated as
  symmetric/Hermitian.
