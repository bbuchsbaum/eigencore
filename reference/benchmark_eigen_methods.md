# Benchmark eigen methods against base and optional references.

Benchmark eigen methods against base and optional references.

## Usage

``` r
benchmark_eigen_methods(
  A,
  k,
  target = largest(),
  repeats = 3L,
  include = c("eigencore", "base", "RSpectra"),
  tol = 1e-08
)
```

## Arguments

- A:

  Matrix or eigencore operator to benchmark.

- k:

  Number of eigenpairs to compute.

- target:

  Eigencore eigenvalue target descriptor.

- repeats:

  Number of timing repetitions.

- include:

  Character vector of method labels to include when available.

- tol:

  Solver tolerance.
