# Native LOBPCG Design

## Goal

Turn the current `reference LOBPCG prototype` into a production native path for
standard Hermitian and generalized SPD problems. The first performance target is
smallest graph-Laplacian eigenpairs with a solve-style preconditioner.

## Current Evidence

On the quick path-Laplacian gate (`n = 200`, `k = 5`, `tol = 1e-8`):

- scalar native thick-restart Lanczos certifies but is slower than RSpectra;
- native block Lanczos prototype certifies only with a full subspace and is
  slower than the scalar path;
- reference preconditioned LOBPCG certifies in about 10 iterations and is
  faster than scalar eigencore and PRIMME;
- reference preconditioned LOBPCG still loses to RSpectra and allocates much
  more memory because factorization/preconditioner application currently runs
  through R/Matrix.
- a native shifted tridiagonal preconditioner preserves the useful iteration
  count on path-Laplacian-style gates and reduces preconditioner memory versus
  the Cholesky-backed setup, but the R-level LOBPCG loop still keeps the memory
  gate open.
- built-in preconditioners are typed functions with inspectable `kind`,
  `native`, `factorization`, `shift`, and call-count diagnostics, so the
  planner and benchmark gates can distinguish opaque R callbacks from
  native-backed solves.
- a native standard Hermitian LOBPCG prototype now runs the iteration loop,
  trial-basis orthogonalization, Rayleigh-Ritz projection, and shifted
  tridiagonal preconditioner application in C++ for dense double and
  `dgCMatrix` operators. On the path-Laplacian staging gate
  (`n = 200, 1000, 2000`, `k = 5`) it certifies, beats scalar eigencore and the
  best certified external reference, and passes the restored `0.25` memory
  ratio gate. Reference failures are recorded as uncertified benchmark rows
  rather than aborting the gate.
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

- No R allocations inside the native iteration loop for the native standard
  Hermitian path with built-in operators and the shifted tridiagonal
  preconditioner.
- Preconditioned path certifies the Laplacian staging gate faster than scalar
  eigencore and the best certified external reference.
- Memory ratio versus the best sparse reference improves materially over the
  reference R/Matrix path and stays inside the preconditioned LOBPCG staging
  ceiling.
- Planner label changes from `reference LOBPCG prototype` to
  `native standard Hermitian LOBPCG prototype` only for supported built-in
  operator/preconditioner combinations.
