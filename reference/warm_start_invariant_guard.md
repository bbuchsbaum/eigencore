# Detect a supplied invariant subspace that cannot establish target identity.

Residual certification establishes that returned pairs are eigenpairs;
it does not prove that an exact invariant start contains the requested
extremal pairs. If a fully supplied basis is already invariant at the
requested tolerance, discard it and use the solver's cold start rather
than risk certifying a non-target invariant block.

## Usage

``` r
warm_start_invariant_guard(op, basis, tol)
```

## Arguments

- op:

  Hermitian operator.

- basis:

  Orthonormal accepted user basis.

- tol:

  Solver tolerance.

## Value

Guard decision, relative invariance residual, and exact operator work.
