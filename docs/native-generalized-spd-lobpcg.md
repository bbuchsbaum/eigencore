# Native Generalized SPD LOBPCG Plan

## Goal

Extend the native standard Hermitian LOBPCG prototype to first-class
generalized SPD problems:

```text
A x = lambda B x,    B > 0
```

The production path must keep the standard-problem speed work while adding
B-inner-product basis management and original-coordinate certification.

## Native ABI

The native generalized solver should consume:

- `A` as an `EigencoreApplyFn`;
- `B` as an `EigencoreApplyFn`;
- optional preconditioner `M` as a typed block apply/solve callback;
- dense block workspaces for `X`, `AX`, `BX`, `R`, `W`, `P`, `Q`, `AQ`, `BQ`;
- projected matrices `H = Q' A Q` and `G = Q' B Q`;
- LAPACK workspaces allocated once at setup.

No native generalized path may materialize sparse `A` or sparse `B`.

## Algorithm Contract

1. Build a random block `X`.
2. B-orthonormalize `X` so `X' B X = I`.
3. Compute `AX`, `BX`, Rayleigh quotients, and residuals:

```text
R = A X - B X Lambda
```

4. Stop only when original-coordinate residual/backward-error certification
   passes.
5. Apply a typed preconditioner `W = M(R)` when available.
6. Build trial subspace `[X, W, P]`.
7. B-orthonormalize the trial subspace, caching both `Q` and `BQ`.
8. Solve the small projected Hermitian eigenproblem in the B-orthonormal basis.
9. Ritz-rotate `X`, update directions `P`, and repeat.

## Certificate Contract

All convergence decisions and returned certificates use original coordinates:

```text
r_i = A v_i - lambda_i B v_i
eta_i = ||r_i|| / ((||A|| + |lambda_i| ||B||) ||v_i||)
```

Returned vectors must satisfy:

```text
V' B V = I
```

The result must report `orthogonality` as B-orthogonality loss, not Euclidean
orthogonality.

## Planner Gate

The planner may use a native generalized LOBPCG label only when all are true:

- `A` and `B` are built-in native operators or native apply callbacks;
- `B` is known SPD by factorization, typed metadata, or a successful native
  setup check;
- the target is `largest`, `smallest`, `largest_magnitude`, or
  `smallest_magnitude`;
- the preconditioner is `NULL` or a typed native preconditioner compatible with
  the generalized residual;
- certification can compute residuals in original coordinates.

The current implementation may use a native generalized prototype label for
built-in dense, diagonal, and CSC operator slices that satisfy those
constraints. Matrix-free `B`, generalized preconditioners, and production
promotion remain outside the gate.

## First Implementation Slices

1. Factor the standard native LOBPCG loop so orthogonalization accepts an
   optional `B` apply callback and returns both `Q` and `BQ`.
2. Add native dense-B B-CholQR2 inside the loop, reusing the current R-facing
   `native_b_cholqr2()` numerical contract without R allocations.
3. Add diagonal and `dgCMatrix` B apply paths without densification. Native
   diagonal and CSC `B` paths are now wired into the generalized LOBPCG
   prototype for supported built-in `A` operators.
4. Add original-coordinate generalized residual/certificate kernels. The
   R-facing residual-backed certificate helper and dense native generalized
   certificate path are now present; the built-in A/B C++ loop now produces
   residuals against `A X - B X Lambda`.
5. Promote planner labels only after dense and sparse generalized adversarial
   tests pass.
