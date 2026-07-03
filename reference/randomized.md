# Randomized SVD method descriptor.

Randomized SVD method descriptor.

## Usage

``` r
randomized(
  oversample = 10,
  n_iter = 2,
  block = NULL,
  normalizer = c("qr", "lu", "none"),
  refine = TRUE
)
```

## Arguments

- oversample:

  Number of extra samples beyond the requested rank.

- n_iter:

  Number of subspace-iteration refinement passes.

- block:

  Optional block size.

- normalizer:

  Basis normalizer to use (`"qr"`, `"lu"`, or `"none"`).

- refine:

  Whether to refine with a certified Lanczos pass.

## Value

An `eigencore_method` descriptor selecting randomized SVD.
