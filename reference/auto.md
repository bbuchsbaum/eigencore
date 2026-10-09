# Automatic solver choice.

Automatic solver choice.

## Usage

``` r
auto(max_subspace = NULL)
```

## Arguments

- max_subspace:

  Optional maximum Krylov subspace size (the ARPACK `ncv`). When the
  planner selects a restarted Krylov route (thick-restart or block
  Lanczos, Krylov-Schur Arnoldi, shift-invert Lanczos or Arnoldi,
  Golub-Kahan, or implicit-Gram Lanczos for SVD) this caps the active
  basis; it is capped at the problem dimension and must leave room for
  the wanted pairs (at least `k + 1`, `k + block` for block Lanczos).
  Routes without a Krylov basis (dense LAPACK, explicit Gram SVD,
  LOBPCG) ignore it. `NULL` (default) uses the route's own default. The
  solve's `maxit` argument is the iteration (restart) limit and never
  changes the subspace size.

## Value

An `eigencore_method` descriptor that lets the planner choose a solver
based on problem structure.
