# Native Block Hermitian Lanczos Design

## Goal

Add a native Hermitian block Krylov path that exercises the V1 block-operator
ABI without weakening planner honesty. The current implementation is a
thick-restart candidate with native locking metadata. It is not yet the
promoted G1 block solver.

## Scope

- Standard Hermitian eigenproblems only.
- Dense double matrices and `Matrix::dgCMatrix` operators.
- Targets supported by the scalar native path: largest, smallest, largest
  magnitude, smallest magnitude.
- Block operator application uses `apply(..., block_cols = b)` in native C++.
- Basis orthogonalization is native and allocation-free after setup. Candidate
  blocks are reorthogonalized with BLAS-3 projections where safe, with
  one-column fallback for residual-tail recovery.

## Algorithm

1. Start from a dense random block `S` with `b` columns.
2. Orthogonalize accepted columns against the existing basis using two-pass
   modified Gram-Schmidt, batched as BLAS-3 projections for normal block
   extension.
3. Apply `A` to accepted columns as one block and cache `AV`.
4. Generate the next candidate block from the block recurrence residual
   `A V_j - V_j alpha_j - V_{j-1} beta_{j-1}^T`, then reorthogonalize against
   locked and active basis vectors as a numerical correction.
5. Maintain the projected Hermitian problem incrementally from block diagonal
   and nearest-block coupling updates, with diagonal-plus-spike reconstruction
   after thick restart. Solve the compact projected problem and extract
   target-ordered Ritz pairs as `V S_k`.
6. Reuse the `AV S - V S theta` residual workspace for residual norms and
   restart tails, then return explicit residuals for certification in the
   existing certificate path.
7. Run a native final subspace polish so the typed certificate sees
   orthonormal locked vectors without adding R-tracked allocations to the
   solver path. The fast path first reuses existing locked residuals when the
   returned basis already satisfies the certificate orthogonality tolerance;
   otherwise it Cholesky-orthonormalizes the returned vectors and refreshes
   residuals directly. If that still does not certify all requested pairs, it
   falls back to a final Rayleigh-Ritz solve in the returned subspace.

## Acceptance

- `lanczos(block > 1)` is accepted by the R API.
- `plan_solver()` labels the path as
  `native block Hermitian Lanczos thick-restart candidate`.
- The block result matches dense oracle values on small Hermitian matrices.
- The block result certifies on the quick sparse Laplacian gate with bounded
  restart subspaces.
- The benchmark script reports an explicit promotion gate. On the current quick
  path-Laplacian gate (`n = 200`, `k = 5`) the installed candidate certifies
  and reaches about `0.98x` to `1.01x` of the RSpectra speed reference on the
  current quick gate scripts while
  passing memory and PRIMME parity. The dense Hermitian quick row (`n = 80`,
  `k = 5`) now certifies through selected native `dsyevr` full-subspace
  extraction and no-allocation dense symmetry detection, but still misses the
  `1.25x` RSpectra speed bar and remains just below memory parity. Current
  tuning uses `block = 2` for smaller rows and `block = 4` with at least
  `16*k` restart space for large sparse `k >= 16` candidate rows, with a
  four-vector capped Ritz pad, sparse
  CSC Ritz-vector residual application, in-solver certificate reuse in the
  benchmark harness, selected-Ritz workspaces, upper-triangle projected-matrix
  copies, selective one-pass reorthogonalization for larger restart spaces, and
  native workspaces that avoid unnecessary zero-fill. The bounded sparse
  restart path now also uses native loops instead of BLAS calls for tiny block
  projection updates, a small-column native combiner for selected Ritz-vector
  rotations, and one reorthogonalization pass for `n >= 64` when the certificate
  remains tight.
  Certificates now require both residual/backward-error convergence and bounded
  orthogonality, so a fast run with nearly duplicate locked vectors is not
  counted as certified. The native lock path now rejects a converged Ritz vector
  if it is numerically dependent on already locked vectors, which fixes the
  duplicate-lock orthogonality failure exposed by ill-conditioned diagonal
  stress cases. Default small-dense block runs with
  `n <= getOption("eigencore.block_dense_full_subspace_max_n", 256)` now take
  an honestly labeled native full-subspace LAPACK/Rayleigh-Ritz path instead of
  driving the restart loop; explicit bounded `max_subspace` still exercises
  restart. Native final subspace polishing now fixes the larger
  `k = 20, n = 1000` orthogonality certificate without R-tracked allocations;
  it also preserves any genuinely locked prefix when the restart budget is
  exhausted, so a noncertified tail cannot rotate away vectors that already met
  the certificate scale. After retuning the larger-k restart space to `8*k`,
  adding guarded final-polish fast paths, splitting Ritz residual timing, and
  adding the small-column selected-vector combiner, that certified
  `n = 1000, k = 20` staging row reaches about `1.45x` to `1.50x` of RSpectra
  in local benchmark runs. The larger `n = 10000, k = 20` path row now
  certifies with the adaptive large-row candidate (`block = 4`,
  `max_subspace = 320`), with local one-iteration timing around `6.0s`.
  With robust RSpectra reference controls (`ncv = 120`, `maxitr = 20000`),
  the same harness certifies RSpectra around `18.2s` and PRIMME around `25.3s`,
  so the large sparse row passes the strict RSpectra speed gate. The candidate
  remains deliberately unpromoted because the dense regression memory gate is
  still just below parity.

## Structured residual sidecar

The boundary-only structured Ritz residual shortcut was tested and rejected:
it failed the direct certificate checks after thick restart. The reason is
that a restarted active basis is no longer a plain unrestarted block Lanczos
basis. Its leading columns are kept Ritz vectors, and those vectors carry their
own small but nonzero residuals. Therefore the valid post-restart invariant is:

```text
A V_active = V_active T_proj + R_sidecar
```

`R_sidecar` must include both kept-Ritz residual columns and the current
block-boundary residual columns. A residual for projected coefficient vector
`s` is `R_sidecar s`, not just `F_boundary s_tail`.

The next structured residual implementation should keep a thin sidecar rather
than a full `n x max_subspace` residual matrix:

- store residual columns for kept Ritz vectors at restart while `A * B_v` is
  still available;
- keep a map from active-basis column to sidecar column for those kept vectors;
- add the current boundary residual block only for the active tail columns;
- form selected residuals as one or two dense products over the sidecar
  columns, then run the same certificate checks as the sparse-apply path.

Acceptance for this optimization is stricter than ordinary speed tuning: it
must pass the thick-restart residual/direct-certificate test, the adversarial
bank, and the benchmark smoke tests before any timing result counts. It should
also reduce `stage_ritz_residual_seconds` on the `k = 20, n = 1000` path case
without increasing total median time or memory enough to hurt the G1 gate.

## Projected solve tuning notes

The projected Rayleigh-Ritz solve is the largest remaining timed stage on the
bounded sparse path, but simple selected-eigenvector substitutions are not a
current win on the tested LAPACK build. Guarded `dsyevr` and `dsyevx` selected
eigenvector paths for algebraic largest/smallest targets both passed the
thick-restart correctness tests, but they increased the `k = 20, n = 1000,
max_subspace = 160` path-Laplacian projected-solve time from about `0.0196s`
to about `0.023s`. Raising the full-solve `dsyevd` threshold so the same row
uses `dsyev` was worse, increasing projected-solve time to about `0.036s`.

Restart-space sweeps now split the small and large sparse regimes. For the
`n = 1000, k = 20` staging row, `block = 2` with `max_subspace = 160` remains
the best tested point; wider block-4 runs certify but are slower. For the
`n = 10000, k = 20` row, `block = 2` only certifies at much wider and slower
subspaces, while `block = 4` certifies at `max_subspace = 320` and is faster
again around `640`. Blocks 5 and 8 certify in some large-row settings but are
slower on the local harness. Small sparse full-subspace runs can remove
restarts at `n = 200`, but their median time is roughly tied with the bounded
default rather than clearly faster.

Native diagnostics now split the legacy `projected_solve` aggregate into
`projection_update`, `projection_copy`, and `projected_eigensolve`. The
selected Ritz-vector gather is reported separately as `selected_vector_copy`;
it is intentionally not included in the `projected_solve` aggregate because
older timing runs counted that work with Ritz residual formation.
Ritz residual diagnostics are also split into `ritz_vector_form`,
`ritz_operator_apply`, `ritz_norm`, and `ritz_final_polish`; the legacy
`ritz_residual` aggregate remains the sum of those four fields.

The accepted projection-update optimization is a fused native small-block
update for the common inner recurrence case. For `block = 2`, the old path
computed the new self block and the previous/current cross block in separate
dot-product scans over the same accepted `A V` columns. The fused path computes
both in one pass for `left_cols <= 4 && right_cols <= 4`, and falls back to the
existing BLAS/general helpers otherwise. On the `k = 20, n = 1000,
max_subspace = 160` path-Laplacian staging row, this reduced
`stage_projection_update_seconds` from about `0.0074s` to about `0.0035s` and
improved the observed RSpectra speed ratio to about `1.18x`, still short of
the `1.25x` G1 gate.

The accepted Ritz-side micro-optimization is a guarded native combiner for
small selected-vector rotations (`selected_count <= 32`) plus reuse of the
already computed final-polish Gram matrix. On the `n = 1000, k = 20` staging
row, this moved `stage_ritz_vector_form_seconds` from roughly `0.011s` to
roughly `0.0042s`. Sparse Ritz application and residual norming are minor.
Restart diagnostics now expose a per-Rayleigh-Ritz history with active subspace
size, selected count, lock counts before/after locking, wanted-pair convergence
count, max residual, and max backward error. This caught the final-polish
partial-lock regression and exposed the large path-Laplacian best-snapshot
recovery issue.
A block-sized Ritz pad was tested and rejected: it reduced selected columns but
increased restarts and total time on both the quick path row and the
`k = 20, n = 1000` staging row. A larger large-`n` pad (`pad = k`) was also
rejected on `n = 10000, k = 20`; it increased Ritz-vector formation time and
did not certify. Retesting `max_subspace` values around `8*k` did not justify a
default change for the `n = 1000` staging row; `160` remains the best current
point there, with `192` roughly tied or worse depending on run noise. The
large row now needs the adaptive block-4 controls instead of a larger block-2
budget.

## G1 blocker

The candidate now forms block-recurrence residuals, maintains a structured
projected problem, uses allocation-free CholQR2 block acceptance when rank
permits, caps the thick-restart Ritz pad, uses CSC block apply for sparse
Ritz-vector residuals, exposes native stage timings, and uses fair certified
benchmark accounting without double-certifying eigencore rows. It is still not
the promoted G1 solver: installed strict quick gates certify and now pass
memory plus PRIMME parity on the path row, but they still fail the RSpectra
`1.25x` speed bar. Dense quick rows are now opt-in benchmark diagnostics via
`--include-dense`; the full dense regression row certifies and passes the
RSpectra speed bar locally, but remains just below memory parity against the
best reference, so it remains a separate release blocker.
The saved G1 pre-promotion baseline is generated by
`inst/benchmarks/bench-g1-candidate-baseline.R`; its strict mode now certifies
the sparse Laplacian, dense Hermitian, clustered, and ill-conditioned diagonal
block-candidate rows. The `n = 1000, k = 20` sparse row passes the RSpectra
speed gate with `block = 2`. The `n = 10000, k = 20` path row now certifies
with `block = 4` and `max_subspace = 320`, and the same benchmark harness shows
RSpectra speed ratio around `3.0x` and PRIMME parity above `4.0x` when RSpectra
uses the robust large-sparse reference controls. Until dense memory parity is
green, the planner must keep the candidate label.
