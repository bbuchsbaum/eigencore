# Create a certified PSD classification policy

Create a certified PSD classification policy

## Usage

``` r
psd_policy(
  symmetry = psd_tolerance(),
  positivity = psd_tolerance(rel = 64 * .Machine$double.eps),
  rank = psd_tolerance(),
  rhs = psd_tolerance(),
  scale = "frobenius",
  symmetry_repair = c("average", "reject"),
  negative_repair = c("clip", "reject"),
  structure_repair = c("canonicalize", "reject")
)
```

## Arguments

- symmetry:

  Scale-aware symmetry tolerance.

- positivity:

  Scale-aware negative-eigenvalue acceptance tolerance.

- rank:

  Scale-aware numerical-rank tolerance.

- rhs:

  Scale-aware strict-solve compatibility tolerance.

- scale:

  Matrix scale. Version 1.3 supports only `"frobenius"`.

- symmetry_repair:

  Average an admitted asymmetric source, or reject it.

- negative_repair:

  Clip admitted negative eigenvalues, or reject them.

- structure_repair:

  Canonicalize an admitted structural defect, or reject it.

## Value

An immutable `eigencore_psd_policy` record.
