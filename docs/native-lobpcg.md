# Native LOBPCG Design

## Goal

Keep standard Hermitian LOBPCG as a diagnostic/prototype path unless a future
design beats the current eigencore scalar/tridiagonal default. Generalized SPD
LOBPCG has its own promoted scoped surface and is documented separately.

## Current Evidence

On the 2026-06-06 path-Laplacian probe (`n = 200, 1000, 2000`, `k = 5`,
`tol = 1e-8`):

- built-in preconditioners are typed functions with inspectable `kind`,
  `native`, `factorization`, `shift`, and call-count diagnostics, so the
  planner and benchmark gates can distinguish opaque R callbacks from
  native-backed solves.
- a native standard Hermitian LOBPCG prototype now runs the iteration loop,
  trial-basis orthogonalization, Rayleigh-Ritz projection, and shifted
  tridiagonal preconditioner application in C++ for dense double and
  `dgCMatrix` operators.
- the native shifted-tridiagonal LOBPCG row certifies `5/5` and beats certified
  external references by about `2.44x`, `8.11x`, and `6.44x`, with memory gates
  green.
- the same rows fail the scalar-speed gate against the current eigencore
  scalar/tridiagonal default, with speed ratios of about `0.65`, `0.32`, and
  `0.14`.
- standard Hermitian LOBPCG/preconditioner promotion is therefore closed as a
  documented no-promotion decision under `bd-01KTEH4G1QPR4RT14B4G78PF1M`; keep
  the native label diagnostic/prototype-only.
- a direct large path-Laplacian probe (`n = 10000`, `k = 20`, `tol = 1e-8`,
  shifted tridiagonal preconditioner) certifies 20/20 requested pairs in about
  `4.8s` locally. This is evidence for the LOBPCG/M milestone, not a substitute
  for G1 block-Lanczos promotion.

## Native ABI

The native solver should consume:

- `A` as an `EigencoreApplyFn`;
- optional `B` as an `EigencoreApplyFn` for generalized SPD;
- optional preconditioner `M` as a block solve/apply callback;
- a typed preconditioner descriptor used to select native-backed solves and to
  report setup/application diagnostics;
- dense block workspaces `X`, `AX`, `R`, `W`, `P`, `Q`, `AQ`;
- projected dense matrices and LAPACK workspaces allocated once at setup.

The preconditioner callback contract mirrors operator application:

```c
int apply_preconditioner(void* impl,
                         int64_t block_cols,
                         const double* R,
                         int64_t ldr,
                         double* W,
                         int64_t ldw,
                         EigencoreWorkspace* workspace);
```

## Algorithm

1. Build a random block `X` and orthonormalize it.
2. Compute `AX`, Rayleigh quotients, and residual block `R = AX - X Lambda`.
3. Stop when certified residual/backward-error tolerance is met.
4. Apply preconditioner `W = M(R)` or use `R` if none is supplied.
5. Build trial subspace `[X, W, P]`, orthogonalize natively, and apply `A`.
6. Solve the small projected Hermitian eigenproblem.
7. Ritz-rotate `X`, update conjugate directions `P`, and repeat.
8. Return final residuals, backward errors, orthogonality, matvec/preconditioner
   counts, and convergence history.

## Acceptance

- The diagnostic path must keep certifying supported path-Laplacian probes.
- Planner labels must stay explicit: `native standard Hermitian LOBPCG
  prototype` is not a production-promotion label.
- Any future promotion attempt must beat the current eigencore
  scalar/tridiagonal default, not only external references, while preserving
  no-densification and certificate checks.
