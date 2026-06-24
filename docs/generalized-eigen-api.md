# Generalized Eigen API Contract

This note freezes the eigencore-native generalized eigen replacement surface
before new public solver exports are added. It is a contract for the staged
implementation under the geigen-replacement epic, not a claim that every
function listed here is already exported.

## Public Names

The primary public names are eigencore names, not geigen-compatible names:

- `eig_partial(A, B = NULL, ...)` is the existing partial-spectrum surface.
  `B` denotes a symmetric/Hermitian positive-definite metric for
  generalized SPD/Hermitian problems.
- `eig_full(A, B = NULL, structure = NULL, vectors = TRUE, ...)` is the full
  dense decomposition surface. It covers standard dense eigenproblems, dense
  generalized SPD/Hermitian problems, and dense general pencils through native
  LAPACK-backed paths.
- `generalized_schur(A, B, sort = NULL, vectors = TRUE, ...)` is the planned
  dense QZ/generalized Schur surface. A short `qz()` alias is not part of the
  current public contract; it can be added only after a separate export audit.
- `generalized_svd(A, B, ...)` is the planned generalized SVD compatibility
  surface. It is deferred until the dense generalized eigen and QZ contracts
  are stable.

There is no primary `geigen()` export. Migration documentation may map common
`geigen::geigen()`, `geigen::gqz()`, and `geigen::gsvd()` calls to the names
above, but exact geigen names are compatibility aliases only if a future audit
explicitly accepts them.

## Meaning Of B

`eigen_problem(metric = )` remains SPD/Hermitian-metric-only. The same is true
for `eig_partial(A, B = ...)`: sparse and operator partial paths interpret `B`
as a metric defining the generalized problem `A x = lambda B x`, with
B-orthogonality and original-coordinate residual certification.

General right-hand pencils with indefinite, singular, or nonsymmetric `B` are
not accepted through `metric = `. Those belong to the future dense/full
surfaces:

- `eig_full(A, B = ...)` for dense finite/infinite generalized eigenvalues.
- `generalized_schur(A, B, ...)` for dense QZ/Schur factors.

Sparse general-pencil partial support is a separate boundary. It can be
claimed only for explicit transformed/factorized/user-solve configurations
whose planner labels and certificates state the solve provenance. Arbitrary
sparse indefinite pencils are not a native production claim.

## Result Contract

Current generalized SPD partial results returned by `eig_partial(A, B = ...)`
use the normal `eigencore_eigen_result` shape:

- `values`
- `vectors`
- `method`
- `plan`
- `certificate`
- `residuals`
- `backward_error`
- `warnings`
- `restart`
- `transform`, when a spectral transform such as shift-invert is used

General dense-pencil results will extend that shape rather than replace it:

- `values` is `alpha / beta` where finite.
- `alpha` and `beta` are always present for general pencils.
- `classification` records `finite`, `infinite`, or `undefined` per value.
- `vectors` contains right generalized eigenvectors when requested.
- `left_vectors` is present only for paths that compute left generalized
  eigenvectors.
- `certificate` states whether a finite-value right residual
  `A v - lambda B v` was computed, skipped, or undefined.

Internal result helpers use `generalized_pencil_values(alpha, beta)` to classify
homogeneous pairs before any dense full solver result is exposed. Finite pairs
are certified with `certify_dense_generalized_pencil()` or
`certify_generalized_pencil_operator()`, both of which use the shared
`eigen_backward_scale()` definition. Infinite or undefined pairs are explicit
failed certificate entries for now; they require QZ/Schur-specific semantics
before they can be presented as certified.

QZ results will expose the Schur form separately from eigenvalue accessors:

- `S` and `T` for the generalized Schur pair.
- `Q` and `Z` when Schur vectors are requested.
- `alpha`, `beta`, `values`, and `classification` accessors matching the
  dense-pencil result contract.
- `method`, `plan`, `certificate`, and `warnings`.

GSVD results will define their own field names in the GSVD child issue before
export.

## Planner Labels

Planner labels are part of the public diagnostic contract. Current
generalized paths must retain honest labels until a child issue explicitly
changes the implementation boundary:

- `native dense generalized SPD LAPACK fallback`
- `native generalized SPD LOBPCG (B-orthogonal, residual certified)`
- `reference generalized SPD LOBPCG (matrix-free fallback)`
- `native transformed generalized SPD B-orthogonal Lanczos`
- `reference generalized SPD B-orthogonal Lanczos refinement`
- `native dense generalized SPD shift-invert (factorized Lanczos)`
- `native tridiagonal generalized SPD shift-invert (factorized Lanczos)`
- `reference generalized SPD Lanczos shift-invert (user solve)`
- `reference generalized SPD Lanczos shift-invert (dense QR)`
- `reference generalized SPD Lanczos shift-invert (sparse LU)`

The dense full decomposition surface `eig_full()` carries its own
full-decomposition labels, distinct from the partial-spectrum labels above:

- `native dense Hermitian LAPACK fallback`
- `native dense complex Hermitian LAPACK fallback`
- `native dense complex general LAPACK fallback`
- `native dense generalized SPD/Hermitian LAPACK full`
- `native dense general pencil LAPACK full`
- `dense LAPACK general eigen oracle (base fallback)`

The last label is intentionally not a native label: the real dense general
standard path routes through base `eigen()` and must say so. The two
`... LAPACK full` labels are the native real/complex SPD/Hermitian-definite and
general-pencil `eig_full(A, B = ...)` paths backed by eigencore-owned LAPACK
kernels.

Native labels require eigencore-owned native kernels or LAPACK calls routed
through the native layer. Dense, reference, user-solve, Matrix-factorization,
and oracle fallbacks must say so in the label and in result provenance.

## Certificate Scale And Provenance

Generalized certificates use the shared backward-error scale
`(||A|| + |lambda| ||B||) ||v||`. Dense native and dense LAPACK-backed paths
use exact Frobenius scale metadata. Operator and sparse paths may combine exact
and estimated norm bounds; if any part of the scale is estimated,
`scale_is_estimate = TRUE` withholds `passed` even when residual convergence is
numerically small. Factorized, transformed, Matrix-backed, and user-solve paths
must preserve their provenance in `plan`, `method`, `restart`, `transform`, or
certificate notes.

## Export Audit

Accepted current exports for generalized eigen work:

- `eig_partial`
- `eigen_problem`
- `plan_solver`
- method descriptors such as `auto()`, `lanczos()`, `lobpcg()`, and
  `shift_invert()`
- target descriptors such as `largest()`, `smallest()`, and `nearest()`
- result/certificate accessors

Exported full dense surface:

- `eig_full`

Planned but not yet exported:

- `generalized_schur`
- `generalized_svd`

Rejected unless a future audit reopens them:

- `geigen`
- `gqz`
- `gsvd`
- `qz`
- helper names that expose implementation details rather than public
  semantics

## Dense And Sparse Fallback Policy

Dense full paths may use LAPACK through eigencore-owned native registration
and must report the exact LAPACK-backed planner label. R prototypes remain
reference code only.

Sparse and operator paths must not silently densify. If a sparse generalized
problem can run only by densifying, the planner must reject it or carry a
reference/oracle label that makes the dense boundary explicit. Test or
migration oracles may densify small fixtures for comparison, but solver paths
may not hide that conversion.
