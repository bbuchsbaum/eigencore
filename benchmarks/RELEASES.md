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
- G1 block Hermitian candidate baseline captured in
  `inst/benchmarks/baselines/g1_candidate_pre.csv`; regenerate it with
  `inst/benchmarks/bench-g1-candidate-baseline.R --save`. The current baseline
  certifies the sparse Laplacian, dense Hermitian, clustered, and
  ill-conditioned diagonal block-candidate rows. Strict mode now fails if any
  block-candidate row is missing, errors, or returns without a passing
  certificate.
- G1 quick tuning grid captured in
  `inst/benchmarks/baselines/g1_tuning_quick.csv`.
- H SVD surface benchmark now supports `--subject=<method>`, so projected
  Golub-Kahan candidate rows can be gated directly against external references
  without changing the default `eigencore` SVD release gate. A quick
  `eigencore_golub_kahan_projected` subject check certifies sampled
  wide-sparse and clustered-dense rows, but remains well below the `1.5x` SVD
  speed gate and memory parity; H is still a performance milestone.
- G1 native block Hermitian Lanczos is promoted for benchmark-proven regimes.
  The default `eigencore` path now passes both strict release scripts with dense
  diagnostics included:
  `Rscript inst/benchmarks/bench-native-hermitian-gate.R --include-dense --strict`
  and
  `Rscript inst/benchmarks/bench-hermitian-sparse.R --include-dense --strict`.
  Representative installed-package gate rows:
  path Laplacian `n = 1000, k = 20` certifies with block size `2`, about
  `1.42x` to `1.49x` faster than certified RSpectra, memory ratio about
  `1.79x`, and PRIMME parity above `5x`;
  path Laplacian `n = 10000, k = 20` certifies with block size `4` and
  `max_subspace = 320`, about `2.77x` to `2.99x` faster than certified
  RSpectra, memory ratio about `1.65x`, and PRIMME parity above `3.8x`;
  dense Hermitian `n = 200, k = 20` certifies through the small-dense
  full-subspace native path, about `1.42x` faster than certified RSpectra,
  memory ratio about `1.08x`, and PRIMME parity above `9x`.
  The native block path includes structured projected recurrence, capped Ritz
  padding, sparse CSC Ritz-vector residual application, in-solver certificate
  reuse in the benchmark harness, selected-Ritz workspace shrink, selected
  dense LAPACK extraction for small-dense full-subspace fallback, tiny-block
  native projection updates, one-pass reorthogonalization for `n >= 64`,
  upper-triangle projected-matrix copies, read-only restart Ritz residual
  norming that preserves `A * V` blocks for retained Ritz vectors, guarded
  native final subspace polishing with residual reuse, residual-refresh and
  Rayleigh-Ritz fallback paths, best-so-far restart snapshot recovery, and
  native malloc-based workspaces where zero-fill is unnecessary.
- Native standard LOBPCG with the shifted tridiagonal preconditioner passes
  `bench-lobpcg-preconditioned.R --strict` on path Laplacians
  `n = 200, 1000, 2000`, `k = 5`, with failed references recorded as
  uncertified rows.
- Native generalized SPD LOBPCG now has built-in dense, diagonal, and CSC A/B
  prototype slices for explicit `lobpcg()` plans. Production promotion remains
  blocked on generalized preconditioning and broader adversarial gates.
