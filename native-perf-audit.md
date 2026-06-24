# Native Performance Audit: `src/`

Date: 2026-06-09
Scope: full survey of `src/` (~17k lines C++) for provable efficiency improvements.
Status: findings and recommendations only — no changes applied beyond the
already-landed `native_operators.cpp` pass (commit `7e4189d`).

## Background: what was already proven and landed

A first pass on `src/native_operators.cpp` landed five changes, all verified
bit-identical against the previous build under fixed RNG seeds, with
1,570 test assertions passing:

1. Removed a dead `work_B = B` copy before `dgesvd` in both randomized SVD
   candidates (`dgesvd` destroys its input; `B` was never read again).
2. `dgemm` -> `dsyrk` for the Gram matrix in
   `dense_randomized_max_orthogonality` (half the flops; upper-triangle scan).
3. Pre-transposed `Q` panels in both CSC `Q^T A` projection kernels.
   Same-binary kernel benchmark (`clang++ -O2`, identical inputs,
   max|diff| = 0): **1.27x-1.62x** across m = 20k..2M, nnz = 400k..10M.
4. `memcpy` packing of `d`/`U`/`V` into result SEXPs in both controllers.
5. `int64_t` offset in `eigencore_col_norms` (`col * rows` overflowed `int`
   past 2^31 elements — correctness hardening).

Evidence recorded in mote bead `bd-01KTQ165VFNR8K3GB4DE43DBZM` (closed).

---

## New findings, ranked by expected payoff

### 1. `scalar_krylov.cpp` — scalar reorthogonalization loops (largest candidate)

**Where:** the scalar Lanczos cycle, lines ~315-327.

The full reorthogonalization runs every iteration as two passes of
vector-at-a-time modified Gram-Schmidt with hand-rolled scalar loops and
`long double` accumulation:

```cpp
for (int pass = 0; pass < 2; ++pass) {
  for (int prev = 0; prev <= j; ++prev) {
    const double* qprev_basis = Q + prev * n;
    long double dot = 0.0L;
    for (int row = 0; row < n; ++row) dot += ...;   // scalar dot
    for (int row = 0; row < n; ++row) z[row] -= ...; // scalar axpy
  }
}
```

This is O(m^2 * n) scalar work per restart cycle — typically the dominant
non-matvec cost of the solver. Two `dgemv` calls per pass
(`tmp = Q^T z; z -= Q tmp`) perform the same projection at BLAS-2 speed, and
**the same file already uses exactly that pattern** for the Golub-Kahan path
(lines ~498-501 and ~584-587), so block-CGS2 semantics are already accepted
in this codebase.

- **Expected win:** several-fold reduction of the reorthogonalization stage
  for large n and m; this stage often dominates once matvecs are cheap
  (CSC/diagonal operators).
- **Risk:** numerics-bearing. MGS2 -> CGS2 changes rounding behavior.
  CGS with two passes carries the same orthogonality guarantee
  ("twice is enough", Giraud et al.), but this must be arbitrated by the
  adversarial test bank plus a before/after orthogonality-loss comparison
  on clustered/ill-conditioned spectra.
- **Proof plan:** same-binary kernel benchmark (both variants compiled
  together), orthogonality-loss distributions before/after on the
  adversarial bank, full test suite.

### 2. `arnoldi.cpp` — per-eigenpair waste in refined extraction

**Where:** refined Ritz extraction, lines ~160-234, inside `for (col < k)`.

Per eigenpair, the loop currently:

1. Re-allocates `z` ((m+1) x m complex), `s`, `vt`, `rwork`, `work` — k times
   instead of once (sizes are loop-invariant).
2. Re-runs the `zgesvd` workspace query — k times instead of once.
3. Forms the Ritz vector with a scalar O(n * m) loop over complex
   coefficients.

Fixes: hoist allocations and the workspace query out of the loop (free;
bit-identical), and split `coeff` into real/imaginary vectors so vector
formation becomes two `dgemv` calls per eigenpair (last-ulp shifts only).

- **Expected win:** moderate; refined extraction runs once per solve but the
  vector-formation loop is O(k * n * m) scalar work that becomes BLAS-2.
- **Risk:** minimal. The hoisting is mechanically bit-identical.
- **Proof plan:** fixed-seed output comparison plus microbenchmark at
  n ~ 1e5-1e6, m ~ 60-200, k ~ 10-50.

### 3. `certificates.cpp` — eight `dgemm` Gram products -> `dsyrk`

**Where:** lines 76, 252, 422/425, 534, 637/640, 755/758, 899.

Each computes an orthogonality Gram `X^T X` via full `dgemm` (2 * n * k^2
flops) where `dsyrk` computes the upper triangle in n * k^2, and the consumer
`max_orthogonality_loss_cert` only needs that triangle. Identical pattern to
the change already proven in `native_operators.cpp`.

- **Expected win:** halves certificate Gram cost. Runs once per solve, but
  certificates are a visible slice of total time for large n (the rSVD stage
  timers showed this).
- **Risk:** minimal; same-triangle results. The B-inner-product variant at
  line 93 (`Q^T (BQ)`) is *not* a `dsyrk` candidate (two distinct operands) —
  leave it.
- **Proof plan:** same as the landed dsyrk change: bit-level comparison of the
  computed triangle, test suite.

### 4. `gram_svd.cpp` — scalar long-double crossprod (needs a decision)

**Where:** `small_column_crossprod_gram` (line 27), called at lines 1027/1028
and 1372/1373 on full-size `U` (m x rank) and `V` (n x rank).

Hand-rolled O(rows * cols^2 / 2) crossprod with `long double` accumulation.
`dsyrk` would be substantially faster, **but** the extended-precision
accumulation may be a deliberate accuracy choice for certificate-grade
orthogonality checks. Note: on Apple Silicon `long double == double`, so the
extra precision only materializes on x86 builds.

- **Recommendation:** decide policy first (see cross-cutting note below);
  if certificate accuracy on x86 matters, keep as is and accept the cost.

---

## Cross-cutting observation: `long double` accumulation loops

There are ~100 `long double` accumulation loops across `src/`
(gram_svd 30, block_lanczos 16, scalar_krylov 14, lobpcg 10, retained_svd 10,
certificates 8, orthogonalization 6, others few). Two facts:

- On AArch64 (Apple Silicon — the primary dev platform), `long double` is
  64-bit: these loops gain **zero** precision while inhibiting SIMD
  vectorization.
- On x86 Linux/Windows, `long double` is 80-bit extended: converting such a
  loop to BLAS (`ddot`/`dnrm2`/`dsyrk`) genuinely trades accuracy away there.

Any conversion should be a per-call-site decision weighing whether the
quantity feeds a certificate (keep precision) or an internal iterate
(speed is fine). A blanket sweep is not recommended.

---

## What already looks tight (do not touch without a profile)

- `block_lanczos.cpp`: preallocated buffer struct, hoisted `dsyev`/`dsyevd`
  workspace, custom small-block kernels with a `dgemm` crossover at 32
  columns, staged timers throughout.
- `native_operators.cpp` randomized controllers: addressed in the landed pass.
- `small_dense.cpp`: one-shot LAPACK drivers with standard workspace queries;
  nothing per-iteration.

## Suggested execution order

| Step | Item | Risk | Expected effect |
|------|------|------|-----------------|
| 1 | certificates.cpp dsyrk (#3) | minimal | -50% certificate Gram flops |
| 2 | arnoldi.cpp hoist + dgemv (#2) | minimal | refined extraction O(knm) scalar -> BLAS-2 |
| 3 | scalar_krylov.cpp CGS2 dgemv (#1) | numerics-bearing | multi-x on reorthogonalization stage |
| 4 | gram_svd crossprod policy (#4) | decision needed | conditional |

Steps 1-2 are mechanical and provable by bit-identical (or
triangle-identical) outputs. Step 3 requires adversarial-bank arbitration and
should ship with before/after orthogonality-loss evidence in the PR per the
repo guardrails (planner labels and certificate provenance are unaffected;
no public API changes).
