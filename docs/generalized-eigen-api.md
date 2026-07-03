# Generalized Eigen API Contract

This note freezes the eigencore-native generalized eigen replacement surface
before new public solver exports are added. It is a contract for the staged
implementation under the geigen-replacement epic, not a claim that every
function listed here is already exported.

## Public Names

The primary public names are eigencore names, not geigen-compatible names:

- `eig_partial(A, B = NULL, ...)` is the existing partial-spectrum surface.
  For Hermitian problems, `B` denotes a symmetric/Hermitian positive-definite
  metric. For general sparse problems, `B` is accepted only on explicitly
  labelled right-hand-pencil partial boundaries.
- `eig_full(A, B = NULL, structure = NULL, vectors = TRUE, ...)` is the full
  dense decomposition surface. It covers standard dense eigenproblems, dense
  generalized SPD/Hermitian problems, and dense general pencils through native
  LAPACK-backed paths.
- `generalized_schur(A, B, sort = NULL, vectors = TRUE, ...)` is the dense
  QZ/generalized Schur surface. A short `qz()` alias is not part of the
  current public contract; it can be added only after a separate export audit.
- `alpha_beta(x)` extracts homogeneous generalized coordinates from dense
  general-pencil, generalized Schur, and generalized SVD results.
- `generalized_svd(A, B, ...)` is the dense GSVD compatibility surface. The
  current promoted path is real dense LAPACK `dggsvd`, so it requires a linked
  LAPACK that provides that deprecated routine. Complex and sparse GSVD remain
  explicit future scope.

There is no primary `geigen()` export. Migration documentation may map common
`geigen::geigen()`, `geigen::gqz()`, and `geigen::gsvd()` calls to the names
above, but exact geigen names are compatibility aliases only if a future audit
explicitly accepts them.

## Meaning Of B

`eigen_problem(metric = )` remains SPD/Hermitian-metric-only when the problem
structure is Hermitian. The same is true for Hermitian
`eig_partial(A, B = ...)`: sparse and operator partial paths interpret `B` as a
metric defining the generalized problem `A x = lambda B x`, with
B-orthogonality and original-coordinate residual certification.

For general-structure partial problems, `B` is a right-hand pencil only when an
explicit sparse partial solver boundary accepts it. The current supported
general sparse boundary is:

- sparse `dgCMatrix` `A` with nonsingular diagonal `B`: native Arnoldi on the
  sparse transformed operator `B^{-1} A`, certified in original coordinates as
  `A v - lambda B v`.

General dense right-hand pencils remain full-decomposition surfaces:

- `eig_full(A, B = ...)` for dense finite/infinite generalized eigenvalues.
- `generalized_schur(A, B, ...)` for dense QZ/Schur factors.

Unsupported sparse general-pencil cases fail before dense fallback. That
includes singular diagonal `B`, non-diagonal sparse `B`, full sparse QZ, and
factorized/user-solve `B` configurations that do not yet have explicit
provenance and certificate semantics.

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

General dense-pencil results extend that shape rather than replace it:

- `values` is `alpha / beta` where finite.
- `alpha` and `beta` are always present for general pencils.
- `classification` records `finite`, `infinite`, or `undefined` per value.
- `classification_policy` records how those labels were decided (see
  "Alpha/Beta Classification Tolerance" below).
- `vectors` contains right generalized eigenvectors when requested.
- `left_vectors` contains left generalized eigenvectors satisfying
  `w^H A = lambda w^H B` when the native backend computes them. Both the
  real (`DGGEVX`) and complex (`ZGGEV`) dense pencil paths compute left
  vectors; they are dropped when `vectors = FALSE`.
- `conditioning` carries eigenvalue/eigenvector reciprocal condition
  numbers when the selected LAPACK routine provides them. The real dense
  pencil path runs `DGGEVX` with balancing and `sense = 'B'`, so
  `conditioning` contains `rconde`, `rcondv`, and the balanced pencil
  one-norms `abnrm`/`bbnrm` with `available = TRUE`. R's bundled LAPACK
  subset does not ship `ZGGEVX`, so complex pencils report
  `available = FALSE` with an explanatory note; a conditioning routine for
  complex pencils requires bundling or linking an expert driver and is
  explicitly out of scope for the current surface.
- `certificate` states whether a finite-value right residual
  `A v - lambda B v` was computed, skipped, or undefined. Right-residual
  certification and conditioning diagnostics are deliberately separate:
  a passing residual certificate asserts a small backward error, while
  `conditioning$rconde` near zero warns that the eigenvalue itself is
  sensitive to perturbations even when the backward error is tiny.

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
- `alpha`, `beta`, `values`, and `classification` fields matching the
  dense-pencil result contract.
- `method`, `plan`, `certificate`, and `warnings`.

GSVD results use their own class and field names rather than the partial-SVD
result shape:

- `alpha`, `beta`, `values`, and `classification` describe generalized
  singular values in homogeneous form.
- `U`, `V`, and `Q` are the orthogonal factors.
- `D1`, `D2`, `R`, and `zero_R` expose the reconstructable LAPACK `dggsvd`
  factors: `A = U D1 zero_R t(Q)` and `B = V D2 zero_R t(Q)` for the real
  dense path.
- `A_factor` and `B_factor` retain the overwritten LAPACK factor workspaces.
- `k`, `l`, and `rank = k + l` expose the effective-rank partition.
- `certificate` records exact Frobenius reconstruction residuals and
  orthogonality loss for the returned real dense factors.

Complex GSVD fails clearly until a complex GSVD driver is bundled or otherwise
made available through the eigencore native layer.

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
- `native transformed sparse general-pencil Arnoldi (diagonal B)`
- `unsupported sparse general-pencil partial solver`

The dense full decomposition surface `eig_full()` carries its own
full-decomposition labels, distinct from the partial-spectrum labels above:

- `native dense Hermitian LAPACK fallback`
- `native dense complex Hermitian LAPACK fallback`
- `native dense complex general LAPACK fallback`
- `native dense generalized SPD/Hermitian LAPACK full`
- `native dense general pencil LAPACK full`
- `native dense generalized Schur QZ LAPACK full`
- `native dense real LAPACK dggsvd GSVD full`
- `dense LAPACK general eigen oracle (base fallback)`

The last label is intentionally not a native label: the real dense general
standard path routes through base `eigen()` and must say so. The two
`... LAPACK full` labels are the native real/complex SPD/Hermitian-definite and
general-pencil `eig_full(A, B = ...)` paths backed by eigencore-owned LAPACK
kernels.

Native labels require eigencore-owned native kernels or LAPACK calls routed
through the native layer. Dense, reference, user-solve, Matrix-factorization,
and oracle fallbacks must say so in the label and in result provenance.
The GSVD label intentionally names `dggsvd` because that real dense path is
not portable to LAPACK builds that omit the deprecated routine.

## Alpha/Beta Classification Tolerance

`generalized_pencil_values(alpha, beta)` labels each homogeneous pair as
`finite`, `infinite`, or `undefined`. Two tolerance policies exist, and every
classified result records which one ran in `classification_policy`:

- `pencil_norm_scaled` (dense pencil and QZ surfaces): `|alpha|` is treated
  as zero when `|alpha| <= tol * norm_A` and `|beta|` when
  `|beta| <= tol * norm_B`, with `tol = sqrt(.Machine$double.eps)`. For real
  pencils `norm_A`/`norm_B` are the one-norms of the balanced pencil
  reported by `DGGEVX` (`abnrm`/`bbnrm`); complex pencils and
  `generalized_schur()` use one-norms of the input matrices. Because LAPACK
  returns `alpha` on the scale of `norm(A)` and `beta` on the scale of
  `norm(B)`, this policy is invariant under joint rescaling
  `(A, B) -> (c A, c B)` and does not misclassify uniformly small or large
  pencils. Near-singular `B` has an explicit boundary: an eigenvalue is
  reported as `infinite` exactly when its `beta` falls below
  `tol * norm_B`, rather than returned as an untrustworthy enormous finite
  number.
- `per_pair_magnitude` (fallback): `tol * max(1, |alpha|, |beta|)` per pair.
  This is used when no valid pencil norms are available, for example when
  classifying homogeneous pairs supplied directly by a caller.

`classification_policy` exposes `policy`, `tolerance`, the per-coordinate
zero thresholds, and the norms used, and `alpha_beta()` passes it through so
users can audit any classification decision. The QZ `sort = "finite"` /
`sort = "infinite"` selectors run inside LAPACK's ordering callback where
pencil norms are not available; they use the per-pair rule for ordering
only, and the final labels reported on the result always come from the
norm-scaled policy.

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
- `generalized_schur`
- `generalized_svd`
- `alpha_beta`

Not yet promoted behind that export:

- complex `generalized_svd`

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

Sparse general-pencil partial support is intentionally narrower than dense QZ:
the current production boundary is nonsingular diagonal `B`, transformed
Arnoldi on `B^{-1} A`, and a generalized right-residual certificate in original
coordinates. Sparse QZ, singular/near-singular `B`, and non-diagonal
factorized/user-solve right-hand pencils remain separate issues.
