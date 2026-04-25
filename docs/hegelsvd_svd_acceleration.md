# HegelSVD Lessons for eigencore SVD Acceleration

`~/code/HegelSVD` is a useful research prototype for eigencore's randomized
SVD roadmap. It should not be copied into eigencore as-is, but its algorithms,
planner heuristics, and benchmark regimes should inform the native randomized
SVD milestone.

## What HegelSVD Contributes

- **Adaptive randomized QB / PCS strategy.** HegelSVD builds a range basis
  `Q`, then solves a compressed SVD core. This maps naturally to eigencore's
  `BlockOperator` model as `Y = A Omega`, optional power iterations,
  orthogonalization, and a small core SVD.
- **Regime-aware planning.** Its benchmark notes show randomized PCS can beat
  `irlba` and RSpectra in moderate/high-rank, slow-decay, nearly low-rank, and
  ill-conditioned regimes, while `irlba` remains better for low-rank problems
  with sharp singular-value gaps.
- **Core solver heuristics.** HegelSVD switches among QR, Gram, full small-core
  SVD, and `irlba` on the compressed core depending on core size and requested
  rank fraction.
- **Benchmark regimes.** The exact-low-rank, nearly-low-rank, slow-decay, and
  ill-conditioned benchmark families should be ported into eigencore's
  reproducible SVD benchmark suite.

## What Not To Copy Directly

- HegelSVD currently assumes dense numeric matrices, while eigencore needs
  dense, sparse, diagonal, and matrix-free operators.
- HegelSVD's C++ implementation uses RcppArmadillo; eigencore's native engine
  should stay aligned with its direct BLAS/LAPACK and operator ABI.
- HegelSVD's accuracy checks are approximation-oriented. eigencore must keep
  residual, backward-error, orthogonality, and certificate provenance as the
  public contract.
- HegelSVD's own notes show some C++ PCS core variants were slower than
  R/LAPACK. Reuse the strategy and benchmark evidence, not every implementation
  detail.

## eigencore Integration Plan

1. Add a native `randomized()` SVD path over `BlockOperator`.
2. Generate sketch matrices deterministically from eigencore's seed plumbing.
3. Compute `Q` using native block applies and native orthogonalization.
4. Solve the compressed core with a planner-selected QR, Gram, or small-SVD
   path.
5. Recompute true SVD residuals:
   `||A v - sigma u||` and `||A* u - sigma v||`.
6. Mark certificates as estimated unless the residual certification/refinement
   path satisfies the requested tolerance with non-stochastic scale provenance.
7. Use HegelSVD-style planner rules:
   - prefer deterministic Golub-Kahan for low `k` and sharp low-rank cases;
   - prefer randomized PCS for larger `k`, slow spectral decay, nearly
     low-rank workloads without sharp gaps, ill-conditioned spectra, and
     tall/wide dense workloads.

## Benchmark Gate

The randomized SVD milestone should include HegelSVD-derived benchmark cases
and report:

- elapsed time;
- memory allocation;
- singular-value error against a dense oracle where feasible;
- true SVD residuals;
- backward error;
- `U` and `V` orthogonality;
- certificate type, norm provenance, and whether scale is estimated;
- speedup against RSpectra, PRIMME, `irlba`, and deterministic eigencore
  Golub-Kahan where available.

The target is not merely to match HegelSVD's approximation behavior. The target
is to convert the regimes where HegelSVD is fast into eigencore paths that are
operator-native, certified, and faster by time-to-certified-answer.
