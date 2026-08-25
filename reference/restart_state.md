# Extract a reusable spectral restart state.

A restart state is an immutable acceleration hint built only from a
result with a passed current-operator certificate. Its public basis is
expressed in original problem coordinates. `retention = "same_operator"`
may additionally retain a versioned method-fitted start block when the
producing route has an admitted adapter; it never retains convergence,
locks, cached operator actions, or a certificate.

## Usage

``` r
restart_state(x, retention = c("basis", "same_operator"))
```

## Arguments

- x:

  A certified eigencore eigen/SVD result or an existing restart state.

- retention:

  Retain only the public basis, or also request eligible same-operator
  method state.

## Value

An `eigencore_restart_state` schema-version-1 record.
