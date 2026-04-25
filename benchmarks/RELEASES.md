# eigencore Benchmark Releases

This file records curated release benchmark summaries. Raw benchmark runs are
written under `inst/benchmarks/results/` and are intentionally ignored because
they are machine-dependent.

## Methodology

- Metric: time to certified answer.
- Timer: `bench::mark`, with certificate recomputation after each method.
- Eigen reference category: sparse path graph Laplacian, smallest eigenpairs.
- SVD reference category: tall-skinny sparse matrix, largest singular triplets.
- Reported columns include median time, allocated memory, maximum residual,
  maximum backward error, orthogonality loss, converged count, seed, and
  eigencore package version.

## Unreleased

- Benchmark plumbing added.
- No release gate numbers claimed yet.
