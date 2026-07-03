# Check an operator adjoint identity.

Check an operator adjoint identity.

## Usage

``` r
check_adjoint(A, trials = 20, tol = 1e-12, seed = NULL)
```

## Arguments

- A:

  Operator-like object.

- trials:

  Number of random block trials.

- tol:

  Relative tolerance for each adjoint identity check.

- seed:

  Optional random seed for reproducible trials.

## Value

An `eigencore_adjoint_check` list with pass/fail status, tolerance,
maximum relative error, per-trial errors, and trial count.
