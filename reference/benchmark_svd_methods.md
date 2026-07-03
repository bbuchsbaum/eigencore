# Benchmark SVD methods against base and optional references.

Benchmark SVD methods against base and optional references.

## Usage

``` r
benchmark_svd_methods(
  A,
  rank,
  repeats = 3L,
  include = c("eigencore", "base", "RSpectra", "irlba", "rsvd"),
  tol = 1e-08
)
```

## Arguments

- A:

  Matrix or eigencore operator to benchmark.

- rank:

  Number of singular values to compute.

- repeats:

  Number of timing repetitions.

- include:

  Character vector of method labels to include when available.

- tol:

  Solver tolerance.
