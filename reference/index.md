# Package index

## User entry points

High-level functions for partial eigenvalue and SVD computation.

- [`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
  : Compute a partial eigendecomposition.
- [`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md)
  : Compute a full dense eigendecomposition
- [`svd_partial()`](https://bbuchsbaum.github.io/eigencore/reference/svd_partial.md)
  : Compute a partial singular-value decomposition.
- [`generalized_schur()`](https://bbuchsbaum.github.io/eigencore/reference/generalized_schur.md)
  : Compute a dense generalized Schur decomposition
- [`generalized_svd()`](https://bbuchsbaum.github.io/eigencore/reference/generalized_svd.md)
  : Compute a dense generalized singular value decomposition

## RSpectra-compatible shims

Drop-in replacements for
[`RSpectra::eigs()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
[`eigs_sym()`](https://bbuchsbaum.github.io/eigencore/reference/eigs_sym.md),
and
[`svds()`](https://bbuchsbaum.github.io/eigencore/reference/svds.md).

- [`eigs()`](https://bbuchsbaum.github.io/eigencore/reference/eigs.md) :
  RSpectra-compatible eigen shim.
- [`eigs_sym()`](https://bbuchsbaum.github.io/eigencore/reference/eigs_sym.md)
  : RSpectra-compatible symmetric eigen shim.
- [`svds()`](https://bbuchsbaum.github.io/eigencore/reference/svds.md) :
  RSpectra-compatible SVD shim.

## Problems and planning

Build problem descriptors and inspect the solver plan.

- [`eigen_problem()`](https://bbuchsbaum.github.io/eigencore/reference/eigen_problem.md)
  : Define an eigenproblem.
- [`svd_problem()`](https://bbuchsbaum.github.io/eigencore/reference/svd_problem.md)
  : Define an SVD problem.
- [`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md)
  : Plan a solver for a problem.
- [`solve(`*`<eigencore_eigen_problem>`*`)`](https://bbuchsbaum.github.io/eigencore/reference/solve.eigencore_eigen_problem.md)
  : Solve a planned eigenproblem.
- [`solve(`*`<eigencore_svd_problem>`*`)`](https://bbuchsbaum.github.io/eigencore/reference/solve.eigencore_svd_problem.md)
  : Solve a planned SVD problem.

## Operators

Block-native linear operators and operator algebra.

- [`linear_operator()`](https://bbuchsbaum.github.io/eigencore/reference/linear_operator.md)
  : Create a block-native linear operator.
- [`as_operator()`](https://bbuchsbaum.github.io/eigencore/reference/as_operator.md)
  : Convert an object to an eigencore operator.
- [`adjoint()`](https://bbuchsbaum.github.io/eigencore/reference/adjoint.md)
  : Return the adjoint operator.
- [`check_adjoint()`](https://bbuchsbaum.github.io/eigencore/reference/check_adjoint.md)
  : Check an operator adjoint identity.
- [`compose()`](https://bbuchsbaum.github.io/eigencore/reference/compose.md)
  : Compose two operators.
- [`crossprod_operator()`](https://bbuchsbaum.github.io/eigencore/reference/crossprod_operator.md)
  : Create A^\* A as an operator.
- [`symmetric_operator()`](https://bbuchsbaum.github.io/eigencore/reference/symmetric_operator.md)
  : Mark an operator as symmetric/Hermitian.
- [`scale_cols()`](https://bbuchsbaum.github.io/eigencore/reference/scale_cols.md)
  : Scale operator columns.
- [`scale_rows()`](https://bbuchsbaum.github.io/eigencore/reference/scale_rows.md)
  : Scale operator rows.
- [`center()`](https://bbuchsbaum.github.io/eigencore/reference/center.md)
  : Center an operator by rows or columns.

## Targets

Selectors describing which part of the spectrum you want.

- [`largest()`](https://bbuchsbaum.github.io/eigencore/reference/largest.md)
  : Target the largest algebraic values.
- [`smallest()`](https://bbuchsbaum.github.io/eigencore/reference/smallest.md)
  : Target the smallest algebraic values.
- [`largest_magnitude()`](https://bbuchsbaum.github.io/eigencore/reference/largest_magnitude.md)
  : Target the largest values by magnitude.
- [`smallest_magnitude()`](https://bbuchsbaum.github.io/eigencore/reference/smallest_magnitude.md)
  : Target the smallest values by magnitude.
- [`largest_real()`](https://bbuchsbaum.github.io/eigencore/reference/largest_real.md)
  : Target the largest real part.
- [`smallest_real()`](https://bbuchsbaum.github.io/eigencore/reference/smallest_real.md)
  : Target the smallest real part.
- [`largest_imaginary()`](https://bbuchsbaum.github.io/eigencore/reference/largest_imaginary.md)
  : Target the largest imaginary part.
- [`smallest_imaginary()`](https://bbuchsbaum.github.io/eigencore/reference/smallest_imaginary.md)
  : Target the smallest imaginary part.
- [`both_ends()`](https://bbuchsbaum.github.io/eigencore/reference/both_ends.md)
  : Target both algebraic ends.
- [`nearest()`](https://bbuchsbaum.github.io/eigencore/reference/nearest.md)
  : Target values nearest a shift.

## Methods

Solver families and method descriptors.

- [`auto()`](https://bbuchsbaum.github.io/eigencore/reference/auto.md) :
  Automatic solver choice.
- [`lanczos()`](https://bbuchsbaum.github.io/eigencore/reference/lanczos.md)
  : Hermitian Lanczos method descriptor.
- [`lobpcg()`](https://bbuchsbaum.github.io/eigencore/reference/lobpcg.md)
  : LOBPCG method descriptor.
- [`golub_kahan()`](https://bbuchsbaum.github.io/eigencore/reference/golub_kahan.md)
  : Golub-Kahan bidiagonalization method descriptor.
- [`randomized()`](https://bbuchsbaum.github.io/eigencore/reference/randomized.md)
  : Randomized SVD method descriptor.
- [`shift_invert()`](https://bbuchsbaum.github.io/eigencore/reference/shift_invert.md)
  : Shift-invert method descriptor.

## Preconditioners

Preconditioner factories for iterative solvers.

- [`shifted_cholesky_preconditioner()`](https://bbuchsbaum.github.io/eigencore/reference/shifted_cholesky_preconditioner.md)
  : Shifted Cholesky preconditioner.
- [`shifted_diagonal_preconditioner()`](https://bbuchsbaum.github.io/eigencore/reference/shifted_diagonal_preconditioner.md)
  : Shifted diagonal preconditioner.
- [`shifted_tridiagonal_preconditioner()`](https://bbuchsbaum.github.io/eigencore/reference/shifted_tridiagonal_preconditioner.md)
  : Shifted tridiagonal preconditioner.

## Spaces and structure tags

- [`euclidean()`](https://bbuchsbaum.github.io/eigencore/reference/euclidean.md)
  : Euclidean vector space descriptor.
- [`hermitian()`](https://bbuchsbaum.github.io/eigencore/reference/hermitian.md)
  : Hermitian/symmetric operator structure descriptor.
- [`general()`](https://bbuchsbaum.github.io/eigencore/reference/general.md)
  : General operator structure descriptor.

## Results, certificates, diagnostics

Inspecting solver output and the numerical evidence.

- [`certificate()`](https://bbuchsbaum.github.io/eigencore/reference/certificate.md)
  : Extract a result certificate.
- [`residuals(`*`<eigencore_eigen_result>`*`)`](https://bbuchsbaum.github.io/eigencore/reference/residuals.md)
  [`residuals(`*`<eigencore_svd_result>`*`)`](https://bbuchsbaum.github.io/eigencore/reference/residuals.md)
  [`residuals(`*`<eigencore_certificate>`*`)`](https://bbuchsbaum.github.io/eigencore/reference/residuals.md)
  : Extract residual diagnostics.
- [`backward_error()`](https://bbuchsbaum.github.io/eigencore/reference/backward_error.md)
  : Extract backward-error diagnostics.
- [`values()`](https://bbuchsbaum.github.io/eigencore/reference/values.md)
  : Extract computed values.
- [`alpha_beta()`](https://bbuchsbaum.github.io/eigencore/reference/alpha_beta.md)
  : Extract homogeneous generalized coordinates.
- [`vectors()`](https://bbuchsbaum.github.io/eigencore/reference/vectors.md)
  : Extract eigenvectors.
- [`left_vectors()`](https://bbuchsbaum.github.io/eigencore/reference/left_vectors.md)
  : Extract left singular vectors or left eigenvectors.
- [`right_vectors()`](https://bbuchsbaum.github.io/eigencore/reference/right_vectors.md)
  : Extract right singular vectors.
- [`diagnostics()`](https://bbuchsbaum.github.io/eigencore/reference/diagnostics.md)
  : Extract diagnostics.
