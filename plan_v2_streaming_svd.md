# plan_v2_streaming_svd.md — Streaming Bidiagonal SVD Update (V2)

> Status: design note (V2 scope). Post-V1. Does not modify any V1
> solver, certificate, or planner code. Subordinate to `vision.md`
> (goals), `prd.json` (scope), and `plan_v1.md` (execution discipline).
> Read `AGENTS.md` before any implementation phase.

Source: Brust, J. J. & Saunders, M. A. (2025), *Fast and Accurate
SVD-Type Updating in Streaming Data*, arXiv:2509.02840v1.

## 1. Goal

Add a **certified online low-rank SVD tracker** to eigencore. Given an
existing rank-`r` bidiagonal factorization `A ≈ Q B Pᵀ`, fold in a
rank-1 / rank-`s` change `A⁺ = A + bcᵀ` in `O(n²)` per update (BGU)
instead of `O(mn·t)` recomputation, and return certified singular
triplets with an honest streaming error certificate.

This makes concrete the existing `prd.json` `v2_scope` entries **“warm
starts and continuation”** and **“dynamic rank selection”**. It is not
an unplanned surface; it is the operational form of those entries.

## 2. Why this paper fits eigencore

1.  **It updates eigencore’s native intermediate.** eigencore’s SVD core
    is committed to Golub–Kahan **bidiagonalization, not normal
    equations** (`vision.md`; `prd.json:222/468/617`; `golub_kahan.hpp`;
    `R/*_golub_kahan.R`). Most streaming-SVD literature updates `UΣVᵀ`;
    this paper updates `Q B Pᵀ` — exactly the factorization eigencore
    produces.
2.  **Matrix-free aligns with the BlockOperator ABI.** BHU needs only
    matvecs; eigencore is operator-native by design.
3.  **It has computable, honest error accounting.** Exact BD error (eqs
    1.4 / 1.6) and two-sided SVD↔︎BD bounds (eq 1.7) map directly onto
    eigencore’s “certification by default” promise.

## 3. Scope / non-goals

In scope:

- **BGU** (Givens, paper Algorithm 2) — production streaming engine,
  `O(n²)` per update.
- **BHU** (Householder / compact-WY, paper Algorithm 1) — secondary
  matrix-free / low-memory + partial-factorization path (~`(m+n+2)k`
  memory vs standard `~mn`).
- §4 subspace tracker (eqs 4.1–4.5): `(r+1)×(r+1)` middle-matrix
  update + lowest-column-norm truncation.
- **Bidiagonal→SVD finish** on the small `(r+1)×(r+1)` core to return
  certified singular triplets. The paper explicitly defers SVD-updating
  to future work (§4.4); eigencore closes that gap on the tiny core
  only.
- Streaming certificate from exact BD error (1.6) + SVD↔︎BD bounds
  (1.7) + standard residual/backward-error/orthogonality diagnostics.

Out of scope (defer to V2.x / V3):

- **RBD** (randomized bidiagonal). Low marginal value — eigencore
  already has a randomized SVD path. Document, do not build.
- Distributed / GPU streaming, file-backed streams.
- Generalized / weighted streaming SVD.

## 4. Product surface (new public API)

Follows the existing problem → plan → solve → result → certificate
shape.

- `svd_stream(A0 = NULL, rank, target = ec_largest(), tol = 1e-8, ...)`
  → an `eigencore_svd_stream` object holding the truncated bidiagonal
  state (`Q, B, P`), rank, target, certificate-of-record, and a planner
  label.
- `update(stream, b, c)` (rank-1) and `update(stream, B_lr, C_lr)`
  (rank-`s`, looped per the paper) → updated stream; S3
  [`update()`](https://rdrr.io/r/stats/update.html) method.
- Accessors reuse existing generics:
  [`values()`](https://bbuchsbaum.github.io/eigencore/reference/values.md),
  [`left_vectors()`](https://bbuchsbaum.github.io/eigencore/reference/left_vectors.md),
  [`right_vectors()`](https://bbuchsbaum.github.io/eigencore/reference/right_vectors.md),
  [`certificate()`](https://bbuchsbaum.github.io/eigencore/reference/certificate.md),
  [`diagnostics()`](https://bbuchsbaum.github.io/eigencore/reference/diagnostics.md),
  `plan()`.
- `as_svd_result(stream)` → materialize a standard
  `eigencore_svd_result` for interop with the rest of the package.

Non-negotiable (`AGENTS.md`): every new exported symbol carries real
public semantics, supporting tests, a planner label, and certificate
provenance. The stream object carries an honest planner label
distinguishing `reference` from native, exactly like current solvers.

## 5. Layering (reference → native promotion, per plan_v1 discipline)

1.  `R/reference_bidiagonal_update.R` — unexported reference prototypes:
    `reference_bgu_rank1()`, `reference_bhu_rank1()`,
    `reference_svd_stream_track()`. Honest `reference`/oracle planner
    label. Permanent adversarial-test baseline; never presented as
    production.
2.  `src/bidiagonal_update.{hpp,cpp}` — native BGU bulge-chasing kernel
    (Givens rotations stored as `(c, s, i, j)` 4-vectors per the paper’s
    BGU memory layout) and native BHU compact-WY kernel; matrix-free
    `b/c` projection (eq 4.2) through the existing `BlockOperator` ABI.
3.  `R/solve_svd_stream.R` — dispatcher + stream state management,
    structured like the `R/solve_svd.R` helpers.
4.  Planner promotion only after the native path passes its gate. Until
    then the public path runs reference code with an honest label
    (`vision.md` planner honesty). Sparse inputs must not silently
    densify.

## 6. Algorithm work breakdown

**BGU core (highest implementation risk):**

- Phase 1: eliminate `b⁺`/`c⁺` spikes; transform `B` into banded `B_pq`
  (subdiagonal `p=1`, two superdiagonals `q=2`) via Givens rotations.
  Handle the zero-permutation shortcut: if `b⁺ = γₙeₙ + γₙ₋₁eₙ₋₁` then
  `B + b⁺c⁺ᵀ` is already bidiagonal.
- Bulge chasing: each rotation is `O(1)` (~10 flops, ≤ 2×5 element
  touch). ~`2n²` to clear `b⁺/c⁺`; `~4n²` total including phase-2
  quaddiagonal→bidiagonal reduction (`2n−3` extra bulges).
- Permutation bookkeeping (paper Fig 3, labels `rs`/`cs`/`bl`) is the
  main source of subtle bugs — exercise heavily against the reference
  oracle.
- Trailing spike reflector `G_{n+1n}···G_{mm−1} b⁺(n+1:m)` for the
  `m > n` block (does not touch `B`).

**Subspace tracker (§4):**

- Project `b, c` onto the orthogonal complement (eq 4.2): `b⁺ = Q⁽¹⁾ᵀb`,
  `b⊥ = b − Q⁽¹⁾b⁺`, `δ = ‖b⊥‖`; symmetric for `c`. Handle `δ = 0` /
  `γ = 0` degenerate cases (no augmentation; stays rank-`r`,
  automatically exact).
- Run BGU on the `(r+1)×(r+1)` middle matrix (eq 4.4); accumulate
  `Q̄⁽¹⁾Q⁽²⁾` and `P̄⁽¹⁾P⁽²⁾`.
- Truncate the lowest-column-norm direction of `B̄⁺` back to rank `r`.
  This is **dynamic rank selection**: expose `max_rank` and a
  column-norm threshold as `svd_until`-style controls.
- Periodic reorthogonalization of accumulated `Q̄⁽¹⁾Q⁽²⁾` / `P̄⁽¹⁾P⁽²⁾`
  (paper monitors orthogonality loss). Reuse `R/orthogonalization.R`.

**SVD finish (closes the §4.4 gap):**

- Each update, run the existing dense bidiagonal→SVD on the tiny `r×r` /
  `(r+1)×(r+1)` core only (`O(r³)` on the small matrix, negligible
  against the stream). Converts the BD tracker into certified singular
  triplets without re-touching `A`.

## 7. Certificate design (the differentiator)

The streaming certificate uses the paper’s own exact accounting:

- BD reconstruction error: `‖A − A_r^BD‖²_F = Σ_{i>r}(αᵢ² + βᵢ₋₁²)` (eq
  1.6) — exact, computed from the maintained `B`.
- SVD↔︎BD gap bounds (eq 1.7) reported as a certificate field. A flat or
  slowly-decaying spectrum yields an honest **wide-bound** certificate
  rather than a false “exact” claim — the method’s main accuracy caveat
  becomes a transparency feature, directly on-vision.
- Standard eigencore residual / backward-error on the finished triplets
  (reuse the `certify` path / `new_certificate()`), plus orthogonality
  drift diagnostics on accumulated `Q, P`.
- One shared scale definition across reference and native paths
  (`AGENTS.md` certificate non-negotiable). Tests must compare both
  paths.

## 8. Milestones (phased, plan_v1 style)

| Phase | Deliverable | Gate |
|----|----|----|
| **S0** | This design note; `prd.json`/`vision.md` deltas; mote issue tree | Reviewed; scope agreed |
| **S1** | `reference_bgu_rank1` + `reference_svd_stream_track` (R prototype, honest label) | Matches LAPACK recompute to 1e-12 on the SuiteSparse rank-1 bank (paper Table 1 set) |
| **S2** | Streaming certificate (1.6/1.7 + residual/orthogonality), reference path | Scale parity reference↔︎dense; adversarial bank incl. flat-spectrum wide-bound case |
| **S3** | `svd_stream()` / [`update()`](https://rdrr.io/r/stats/update.html) public API + planner labels + result interop | API tests; planner honesty tests; `as_svd_result()` round-trip |
| **S4** | Native BGU C++ kernel (Givens bulge-chasing) via BlockOperator ABI | Bit-comparable to reference within tol; sparse never densifies |
| **S5** | Native BHU kernel (matrix-free / low-memory / partial) | Memory ≈ `(m+n+2)k` vs standard `~mn`; matrix-free correctness |
| **S6** | Promotion + benchmarks vs Brand iSVD / RPI (Deng et al.) / `zgebrd` recompute | BGU `O(n²)` scaling confirmed; faster than recompute at machine-precision residuals |

S1–S3 are the minimum shippable V2 increment (certified streaming
tracker on the honest reference engine). S4–S6 are the native
performance promotion.

## 9. Risks & mitigations

- **BGU bulge-chasing / permutation bookkeeping is intricate** (Fig 3).
  Mitigate: reference prototype first as the correctness oracle; native
  kernel must be bit-comparable before promotion.
- **BD ≠ truncated SVD on flat spectra.** Mitigate: SVD finish on the
  small core + eq 1.7 bounds surfaced in the certificate. The honesty
  contract is a feature here, not a workaround.
- **Orthogonality drift over thousands of updates** (paper’s known
  issue). Mitigate: reuse existing reorthogonalization + drift
  diagnostic in the certificate; expose reorth cadence as a control.
- **No liftable code** (paper ships Fortran 90 + C/Matlab/Python
  wrappers). Accept: clean-room C++ in eigencore’s kernel layer; the
  paper fully specifies the algorithms.

## 10. Doc / spec deltas required before build (S0)

- `prd.json`: add a `v2_scope.product_additions` entry, e.g.
  `svd_stream(A0, rank)` + `update(stream, b, c)`, explicitly tied to
  the existing “warm starts and continuation” / “dynamic rank selection”
  algorithmic additions.
- `vision.md`: one sentence acknowledging streaming/online certified SVD
  as a V2 product surface (currently absent).
- `plan_v1.md`: no change (post-V1). This document is referenced from V2
  planning, not from the V1 execution order.

## 11. References

- Brust & Saunders 2025, arXiv:2509.02840v1 — Algorithms 1 (BHU) and 2
  (BGU), Theorem 3.1, §3.4 (Givens Low Rank), §4 (subspace tracking),
  eqs 1.4/1.6/1.7 (error accounting), Table 1 / Figs 4–6 (benchmarks).
- Brand 2006 (incremental SVD) — comparison baseline.
- Deng et al. 2024 (RPI) — comparison baseline.
- eigencore: `vision.md`, `prd.json` (`v2_scope`), `AGENTS.md`
  (non-negotiables), `plan_v1.md` (execution discipline),
  `golub_kahan.hpp`, `R/reference_golub_kahan.R`, `R/solve_svd.R`,
  `R/certification.R`, `R/orthogonalization.R`.
