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
- Current strict block-candidate gates on path Laplacian `n = 200, k = 5`
  certify all requested pairs but are not releasable. By default,
  `bench-native-hermitian-gate.R --quick --block-candidate --strict` and
  `bench-hermitian-sparse.R --quick --block-candidate --strict` now gate only
  the sparse Hermitian release case from the PRD. Dense Hermitian diagnostics
  are opt-in with `--include-dense` so they remain visible without blocking the
  sparse G1 gate. With the current installed package, the default quick native
  gate reports sparse path-Laplacian speed ratio about `1.0`, memory ratio
  about `1.01`, and PRIMME parity about `1.9`; it still misses the required
  `1.25x` RSpectra speed bar. The opt-in dense diagnostic certifies through
  selected native `dsyevr` extraction and no-allocation dense symmetry
  detection, but still misses the `1.25x` RSpectra speed bar and remains just
  below memory parity at about `0.986x`. The current candidate includes
  structured projected recurrence, capped Ritz padding, sparse CSC Ritz-vector
  residual application, in-solver certificate reuse in the benchmark harness,
  selected-Ritz workspace shrink, selected dense LAPACK extraction for
  small-dense full-subspace fallback, tiny-block native projection updates,
  one-pass reorthogonalization for `n >= 64`, upper-triangle projected-matrix copies,
  read-only restart Ritz residual norming that preserves `A * V` blocks for
  retained Ritz vectors, guarded native final subspace polishing with residual
  reuse, a residual-refresh fast path, and
  Rayleigh-Ritz fallback, and native malloc-based workspaces where zero-fill is
  unnecessary. Certificate pass/fail now includes an
  orthogonality gate, so fast runs with duplicate locked vectors are excluded
  from certified gates; the native lock path now rejects numerically duplicate
  locked Ritz vectors before they enter the returned basis. Default small-dense
  block runs with `n <= getOption("eigencore.block_dense_full_subspace_max_n",
  256)` now use an honestly labeled native full-subspace LAPACK/Rayleigh-Ritz
  path; explicit bounded `max_subspace` still exercises restart. A
  `k = 20, n = 1000` staging row with
  `max_subspace = 160` now certifies and passes memory/PRIMME parity, but
  reaches only about `1.09x` of RSpectra's median time. A larger
  `n = 10000, k = 20` path-Laplacian probe currently exhausts the block
  candidate's restart budget without certification, while the PRIMME-inclusive
  full gate can run too long for tight-loop tuning. The production
  block label remains withheld because
  the promotion gate still fails the required `1.25x` RSpectra speed ratio
  and the large sparse row is not yet certifiable.
- Native standard LOBPCG with the shifted tridiagonal preconditioner passes
  `bench-lobpcg-preconditioned.R --strict` on path Laplacians
  `n = 200, 1000, 2000`, `k = 5`, with failed references recorded as
  uncertified rows.
- Native generalized SPD LOBPCG now has built-in dense, diagonal, and CSC A/B
  prototype slices for explicit `lobpcg()` plans. Production promotion remains
  blocked on generalized preconditioning and broader adversarial gates.
