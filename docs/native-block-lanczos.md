# Native Block Hermitian Lanczos Design

## Goal

Add a native Hermitian block Krylov path that exercises the V1 block-operator
ABI without weakening planner honesty. The first implementation is a prototype
block Arnoldi/Lanczos basis builder with Rayleigh-Ritz extraction. It is not the
final thick-restart, locking block solver.

## Scope

- Standard Hermitian eigenproblems only.
- Dense double matrices and `Matrix::dgCMatrix` operators.
- Targets supported by the scalar native path: largest, smallest, largest
  magnitude, smallest magnitude.
- Block operator application uses `apply(..., block_cols = b)` in native C++.
- Basis orthogonalization is native and allocation-free after setup, but may
  orthogonalize candidate columns one at a time inside the native loop.

## Algorithm

1. Start from a dense random block `S` with `b` columns.
2. Orthogonalize accepted columns against the existing basis using two-pass
   modified Gram-Schmidt.
3. Apply `A` to accepted columns as one block and cache `AV`.
4. Generate the next candidate block from `A * V_last`, orthogonalize against
   the full active basis, and repeat until `max_subspace` columns are built or
   the block breaks down.
5. Form `H = V' A V`, symmetrize, solve the projected dense Hermitian problem,
   and extract target-ordered Ritz pairs as `V S_k`.
6. Return explicit residuals `||A v_i - lambda_i v_i||` for certification in
   the existing certificate path.

## Acceptance

- `lanczos(block > 1)` is accepted by the R API.
- `plan_solver()` labels the path as `native block Hermitian Lanczos prototype`.
- The block result matches dense oracle values on small Hermitian matrices.
- The block result certifies on the quick sparse Laplacian gate when enough
  subspace is requested.
- If benchmarks underperform the scalar path or references, this remains a
  prototype path and does not count as Milestone G1.

