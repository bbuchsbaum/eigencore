# Eigencore 1.2 reusable restart-space contract

- Status: approved implementation detail for the frozen 1.2 workflow contract
- Detail contract version: 1
- Package target: 1.2.0
- Date: 2026-08-25
- Mote: `bd-01KY8K55H75Z0J40CTRPZ068A7`
- Normative parent: `docs/plans/1.2-reusable-workflow-contract.md`

## Decision summary

Eigencore restart state is an immutable, opt-in acceleration hint. Its stable
public meaning is an orthonormal basis in the original problem coordinates.
It is not a saved solver, a convergence claim, or a certificate. A later solve
may reuse that basis across a coordinate-compatible operator revision, but it
must invalidate every operator-dependent claim and produce a fresh result and
certificate for the current operator.

Version 1 makes five concrete decisions.

1. `restart_state()` accepts certified eigen or SVD results, but a receiving
   solver consumes a state only when its adapter explicitly advertises basis
   reuse. Unsupported receiving routes fail with `method_incompatible`; they
   never ignore a supplied state.
2. The first production basis adapter is the existing standard real Hermitian
   Lanczos start-subspace family, across dense double, `dgCMatrix`, and callback
   operators. This is the route already covered by `initial_subspace`.
3. The first eligible same-operator method payload is a versioned Lanczos start
   block. It contains no recurrence, locked status, cached operator action,
   convergence flag, or certificate. It avoids repeating method-specific basis
   compression while keeping the public basis representation independent of
   the Lanczos implementation.
4. SVD state represents both original coordinate spaces explicitly. Version 1
   permits extraction and serialization, but no SVD route claims consumption
   until a Golub--Kahan adapter is implemented and independently certified.
5. `reuse = "auto"` may downgrade same-operator payload reuse to basis-only
   reuse. It may not downgrade an incompatible public basis to a cold start.
   `reuse = "same_operator"` never downgrades.

These choices are the smallest implementation that fulfills the common 1.2
contract without exposing current native workspace layouts. Future additive
adapters can support retained Golub--Kahan or richer exact-revision payloads
without changing the public state schema.

## Non-negotiable invariants

- A state is immutable under extraction, validation, execution, and
  serialization.
- The public basis is expressed in original coordinates and is always checked
  and re-orthonormalized at the receiving boundary.
- Operator identity determines whether method state may be considered; a
  coordinate signature determines whether the public basis may be considered.
- A changed revision invalidates recurrences, locked claims, cached operator
  actions, projected quantities, residuals, convergence, target ordering, and
  certificates.
- No historical certificate or convergence flag is an input to a new
  certificate decision.
- Every accepted solve runs the current operator and computes current
  projections, Ritz values or singular values, residuals, orthogonality,
  ordering, convergence, and certification.
- State validation finishes before the first application of the current
  operator. Rejected state therefore has zero operator work.
- Retention is off by default and every retained byte is reported as retained
  memory, never as peak allocation or peak RSS.

## Public entry points

```r
restart_state(x, retention = c("basis", "same_operator"))

solve(
  plan,
  restart_state = NULL,
  reuse = c("auto", "basis_only", "same_operator"),
  retain_state = c("none", "basis", "same_operator"),
  replan = FALSE
)
```

`restart_state()` is both an accessor and a constructor from a result. If `x`
is already an `eigencore_restart_state`, it validates the object and returns an
immutable copy with no increase in retention. If `x` is a result, extraction
requires computed vectors in the result's private execution record or public
result fields. It never reruns a solve to manufacture missing vectors.

`retain_state` controls the state attached to the new result. It does not
change `vectors`, which continues to control public result vectors. An adapter
may retain an internal basis needed by `retain_state` even when the caller asks
not to return vectors, but that work and memory must be visible in the result.
Version 1 does not add that hidden-vector behavior: retention errors clearly
when no valid basis is available.

The default `retain_state = "none"` always returns
`result$restart_state = NULL`. Supplying `restart_state` does not imply that a
new state should be retained.

## Schema-version-1 state

An `eigencore_restart_state` is a classed list with exactly these required
known fields. Additive unknown fields are preserved under schema version 1.

| Field | Required type and meaning |
|---|---|
| `schema_version` | Integer scalar `1L`. |
| `basis` | An eigen basis matrix or an SVD basis record, described below. |
| `operator_identity` | Named identity record for `A` and, if present, `B`. |
| `problem_signature` | Classed `eigencore_problem_signature` record. |
| `requested` | Positive integer `k` or rank of the producing solve. |
| `method` | Non-empty actual-method label of the producing result. |
| `method_state` | `NULL` or one classed, versioned adapter payload. |
| `provenance` | Classed audit record; never executable evidence. |
| `serialization` | Portability and originating-session record. |
| `memory` | Complete or lower-bound `eigencore_memory` record. |

The state constructor deep-copies ordinary R containers. Callback environments
retain the same limitations as executable plans: user-supplied identity and
revision are a provenance covenant, not a magical freeze of a mutable closure.

### Eigen basis

For a standard eigenproblem, `basis` is a finite numeric matrix with one row
per operator-domain coordinate and at least one column. It is Euclidean-
orthonormal for a standard problem. For a generalized Hermitian problem it is
stored in original coordinates with `basis_metric = "B"` in the problem
signature and is B-orthonormal with respect to the producing metric revision.

A receiving generalized solver never assumes stored B-orthogonality. It
recomputes `B %*% Q` and B-orthonormalizes against the current metric before
using the basis. A changed metric revision therefore remains eligible for
basis-only reuse when all coordinate fields agree, but it invalidates every
stored metric-dependent quantity.

### SVD basis

For an SVD problem, `basis` is a classed `eigencore_svd_basis` list:

```text
schema_version = 1L
left             # m by p finite numeric matrix, or NULL
right            # n by q finite numeric matrix, or NULL
```

`left` is in the original codomain of `A`; `right` is in its original domain.
Internal transposition of a wide operator does not swap these public meanings.
Each non-`NULL` side is independently orthonormal. At least one side must be
present, and the problem signature records which sides are available.

Version 1 has no receiving SVD adapter, so an SVD state is useful for stable
inspection, memory accounting, and forward-compatible serialization but is
rejected by `solve(plan, restart_state = state)` with
`method_incompatible`. A future Golub--Kahan adapter must declare which side it
consumes and must not reinterpret orientation implicitly.

## Problem and coordinate signature

`problem_signature` has class `eigencore_problem_signature`, schema version 1,
and these fields:

```text
schema_version
problem_type        # "eigen" or "svd"
dim                 # c(n, n) or c(m, n)
dtype
structure
target_family
metric_present
coordinate_id_A
coordinate_id_B     # NULL without a metric
basis_sides         # "eigen", or a subset of c("left", "right")
basis_metric        # "euclidean" or "B"
```

The version-1 coordinate identifiers are deterministic strings derived from
the problem type, dimensions, dtype, structure, metric presence, and public
side. They do not include operator identity or revision. This separation is
intentional:

- identity and revision answer whether operator-dependent method state is
  reusable;
- the coordinate signature answers whether numerical basis columns have the
  same meaning and shape.

Eigencore does not infer coordinate equivalence across permutations,
regridding, basis changes, or unrelated domains. A future explicit transform
would require a new reviewed contract. Matching dimensions alone are not an
authorization to transform coordinates, but they are the version-1 default
coordinate identity for built-in R matrices, whose row/column order is their
coordinate system. Callback authors use their stable `operator_id` lineage to
assert continuity across revisions; a different callback lineage with the same
shape is eligible only for `reuse = "basis_only"` or `"auto"` and is recorded
as `changed_operator`.

Targets are reduced to adapter-level families rather than exact labels:

- Hermitian extremal: largest, smallest, and their magnitude aliases;
- Hermitian interior: nearest or interval targets;
- general eigenvalue targets;
- SVD extremal and SVD interior.

Changing target within a family keeps the basis coordinate-compatible but
invalidates method state unless an adapter explicitly declares that transition.
Version 1 Lanczos method state requires the exact target descriptor token.

## Operator relation

Validation compares the state's named identity record with the plan's identity
record after validating both schemas.

| Relation | Definition |
|---|---|
| `same_operator` | Every identity component has identical `operator_id`, `revision`, dimension, dtype, and structure. |
| `changed_revision` | Each component has the same `operator_id`, at least one revision differs, and coordinate signatures agree. |
| `changed_operator` | At least one `operator_id` differs while coordinate signatures agree. |
| incompatible | Coordinate signatures differ, an identity component is missing, or an opaque session-local identity is no longer valid. |

For generalized problems, `A` and `B` are compared independently. A changed
revision of either component yields `changed_revision`; a changed lineage of
either yields `changed_operator`. A missing or newly introduced metric is
coordinate-incompatible.

## Method-state adapter protocol

`method_state` is never interpreted generically. A registered adapter owns:

```text
kind
schema_version
adapter_version
problem_type
method_family
controls_token
target_token
operator_identity_token
portable
payload
```

The registry key is `(problem_type, method_family, kind, adapter_version)`.
Validation checks the key and all tokens before inspecting `payload`. Unknown
keys or future versions are stale method state, not an invitation to guess.

An adapter must implement four operations:

1. `extract(result)` creates a payload from a completed solve;
2. `validate(payload, plan, identity)` performs no operator application;
3. `prepare(payload, public_basis, plan)` returns an execution start object;
4. `memory(payload)` reports retained R/native bytes and portability.

The execution path records `method_state_used = TRUE` only when `prepare()`
returns a start object that the solver actually consumes. Merely validating or
discarding a payload does not count as use.

### Version-1 Lanczos start-block adapter

The adapter kind is `hermitian_lanczos_start_block`, adapter version 1. It is
eligible only when all of the following hold:

- producing and receiving routes are standard real Hermitian Lanczos routes;
- there is no generalized metric or shift-invert transform;
- dense double, `dgCMatrix`, and callback paths resolve to the existing native
  block/scalar or reference scalar Lanczos start seam;
- the exact target descriptor, block width, reorthogonalization policy, and
  start-consuming route family match;
- operator identity and revision match exactly; and
- the payload is finite, has the expected row count and width, and passes its
  stored integrity token.

Its payload is the already validated and method-width-fitted start block plus
the public-basis column indices and deterministic compression/exploration
metadata needed to reproduce it. Fitting mixes every retained block with a
finite deterministic direction from its orthogonal complement when one exists.
This target-safety tail ensures that a previously invariant Ritz basis remains
a hint rather than a self-confirming target claim. The current invariant-start
guard remains authoritative: if it still rejects the fitted block,
`reuse = "auto"` records the rejection and continues through the safe
basis/cold path, while `reuse = "same_operator"` fails with
`stale_method_state` instead of pretending the payload was consumed.

The fitted block retains a nonzero contribution from the public basis; a pure
random replacement is not reuse. If the basis spans the full coordinate space
and no exploratory complement exists, the adapter reports method state
unavailable. Its payload intentionally contains none of the following:

- a Lanczos recurrence or tridiagonal projection;
- locked eigenvalue or convergence claims;
- cached `A %*% Q`;
- residuals or a certificate;
- an external pointer or native workspace.

The payload is portable whenever the plan and public basis are portable. The
same-operator mode skips only repeated adapter compression; the solver still
runs its ordinary iteration and fresh certificate. Rich recurrence retention
is outside 1.2 and would require a new adapter kind and numerical court.

All other current routes report `method_state_available = FALSE`. In
particular, existing internal retained Golub--Kahan structures are not exposed
as public method state in version 1: their orientation, locked claims, cached
actions, and native ABI remain implementation-specific and are not yet a
stable cross-solve adapter.

## Construction and retention

State construction follows this order:

1. validate the result class and its executable plan;
2. require a current certificate that passed; an uncertified result cannot
   produce reusable state even when its solve used `certify = FALSE`;
3. extract result vectors in original coordinates;
4. reject non-finite, empty, or dimensionally invalid bases;
5. orthonormalize the public basis and record accepted/rejected rank;
6. capture operator identities and the coordinate signature;
7. if `retention = "same_operator"`, ask the producing route adapter for a
   payload; absence is recorded, not fabricated;
8. build serialization and memory records; and
9. deep-copy and return the immutable state.

For eigen results, the public basis is the returned Ritz-vector matrix. For SVD
results, available public `u` and `v` sides populate the SVD basis record. A
result with no reusable vectors errors with `corrupt_state` and names the
missing basis field.

`solve(plan, retain_state = ...)` performs the same construction after final
current-operator certification and attaches the state to the result. It never
constructs state from a failed fallback attempt or a pre-fallback route.

## Validation before operator application

The receiving boundary performs these checks in order:

1. state class and integer schema version;
2. required fields and types;
3. finite basis values, dimensions, nonzero rank, and orthogonality tolerance;
4. problem and coordinate signature compatibility;
5. identity schemas, portability, and current-session compatibility;
6. requested-size and target-family compatibility;
7. receiving route's basis-adapter capability;
8. method payload schema, adapter key, control token, target token, identity
   token, integrity token, and portability; and
9. requested reuse policy.

Only then may basis preparation apply a metric or the solver apply `A`.
Validation failures carry the state error class and machine-readable fields
`code`, `field`, `expected`, and `actual`.

## Reuse decision algorithm

For `reuse = "basis_only"`:

1. require coordinate compatibility and a receiving basis adapter;
2. validate, re-orthonormalize, and target-safely fit the public basis to the
   receiving start seam;
3. ignore method state without validating its internal payload;
4. invalidate all operator-dependent categories; and
5. execute with `basis_used = TRUE`, `method_state_used = FALSE` only when the
   fitted start is actually consumed; a safety-guard discard is recorded.

For `reuse = "same_operator"`:

1. require relation `same_operator`;
2. require a compatible, portable method payload and receiving adapter;
3. fail on every mismatch rather than downgrade; and
4. execute the prepared start with both `basis_used` and
   `method_state_used` true.

For `reuse = "auto"`:

1. when relation is `same_operator`, try compatible method state;
2. if it is absent, unknown, stale, nonportable, target-incompatible, or
   method-incompatible, record the reason and try the public basis;
3. for `changed_revision` or `changed_operator`, skip method-state inspection
   and use only the coordinate-compatible public basis; and
4. if the basis or receiving adapter is incompatible, error rather than run a
   silent cold start.

Every accepted state transition invalidates at least `projection`, `residuals`,
`convergence`, and `certificate`. Basis-only transitions additionally
invalidate `recurrence`, `locked`, `cached_operator_actions`, and
`method_state`. The version-1 start-block payload does not claim any of those
categories, so its same-operator transition still invalidates them and records
that only the method-fitted start block was reused.

## State-transition record

Every result has one `eigencore_state_transition`, schema version 1:

```text
schema_version
relation
basis_used
method_state_used
invalidated
reason
source_operator_identity
destination_operator_identity
reuse
adapter
```

`adapter` is `NULL` for cold starts and contains adapter kind/version when a
state is consumed. `reason` is a classed code/message record, not prose parsed
by downstream code.

Without a restart state, relation is `cold_start` or `initial_subspace`.
With a restart state, relation is the operator relation even if `auto` discarded
method state and used only the basis.

## Adversarial transition table

| State and destination | `auto` | `basis_only` | `same_operator` | Required invalidation or failure |
|---|---|---|---|---|
| Same operator/revision, matching Lanczos payload | Use fitted start block | Refit public basis | Use fitted start block | Fresh projection, residuals, convergence, ordering, certificate. |
| Same operator/revision, no payload | Use public basis | Use public basis | Error `method_incompatible` | All method state categories invalidated. |
| Same operator/revision, unknown payload version | Use public basis with reason | Use public basis without payload inspection | Error `stale_method_state` | Payload never deserialized as executable state. |
| Same operator/revision, changed block/control token | Use public basis with reason | Use public basis | Error `method_incompatible` | Old fitted start invalidated. |
| Same lineage, changed `A` revision | Use public basis | Use public basis | Error `operator_incompatible` | Recurrence, locks, cached actions, projection, residuals, convergence, certificate. |
| Same generalized `A`, changed `B` revision | Re-B-orthonormalize basis | Re-B-orthonormalize basis | Error `operator_incompatible` | All B-dependent quantities plus all certificate claims. |
| Different lineage, same coordinate signature | Use public basis | Use public basis | Error `operator_incompatible` | All operator-dependent state. |
| Same shape, changed dtype or structure | Error `coordinate_incompatible` | Same error | Same error | No operator application. |
| Dimension or metric-presence change | Error `coordinate_incompatible` | Same error | Same error | No operator application. |
| Exact target changes within family | Use public basis with reason | Use public basis | Error `method_incompatible` | Method payload and target ordering. |
| Target-family or method changes to another basis-capable route | Use public basis with reason | Use public basis | Error `method_incompatible` | Method payload, recurrence, locks, ordering. |
| Receiving route has no basis adapter | Error `method_incompatible` | Same error | Same error | State is never silently ignored. |
| Opaque callback state restored in another session | Error `session_incompatible` | Same error | Same error | Failure before callback invocation. |
| Portable basis plus nonportable payload after restore | Use public basis with reason | Use public basis | Error `session_incompatible` | Payload invalidated; basis remains eligible. |
| Non-finite or rank-zero basis | Error `corrupt_state` | Same error | Same error | Failure before identity or operator work. |
| Unknown state schema | Error `unsupported_schema` | Same error | Same error | No upgrade or operator work. |
| Caller mutates a state field after construction | Error `corrupt_state` | Same error | Same error | Integrity-token mismatch before apply. |

## Worked same-operator example

```r
A <- diag(seq(20, 1))
plan <- plan_solver(
  eigen_problem(A, target = largest()),
  k = 4,
  method = lanczos(block = 4),
  tol = 1e-10
)

first <- solve(plan, retain_state = "same_operator")
state <- restart_state(first, retention = "same_operator")
second <- solve(plan, restart_state = state, reuse = "same_operator")

stopifnot(second$state_transition$relation == "same_operator")
stopifnot(second$state_transition$basis_used)
stopifnot(second$state_transition$method_state_used)
stopifnot("certificate" %in% second$state_transition$invalidated)
stopifnot(work(second)$certification_operator_columns > 0)
stopifnot(isTRUE(certificate(second)$passed))
```

The method payload is only the fitted Lanczos start block. The second solve
does not inherit locks, convergence, Ritz values, residuals, or certification.

## Worked changed-operator `A - rho B` example

```r
A0 <- diag(seq(12, 1))
B0 <- diag(seq(0.5, 1.3, length.out = 12))

make_shifted <- function(rho) {
  M <- A0 - rho * B0
  linear_operator(
    dim = dim(M),
    apply = function(X, alpha = 1, beta = 0, Y = NULL) {
      out <- alpha * (M %*% X)
      if (is.null(Y) || beta == 0) out else out + beta * Y
    },
    structure = hermitian(),
    operator_id = "A-rho-B-continuation",
    revision = sprintf("rho=%.17g", rho),
    portable = TRUE,
    metadata = list(frobenius_norm = sqrt(sum(M^2)))
  )
}

p1 <- plan_solver(eigen_problem(make_shifted(0.10)), k = 3,
                  method = lanczos(block = 3), tol = 1e-10)
f1 <- solve(p1, retain_state = "same_operator")

p2 <- plan_solver(eigen_problem(make_shifted(0.11)), k = 3,
                  method = lanczos(block = 3), tol = 1e-10)
f2 <- solve(p2, restart_state = restart_state(f1), reuse = "auto")

stopifnot(f2$state_transition$relation == "changed_revision")
stopifnot(f2$state_transition$basis_used)
stopifnot(!f2$state_transition$method_state_used)
stopifnot(all(c(
  "method_state", "recurrence", "locked", "cached_operator_actions",
  "projection", "residuals", "convergence", "certificate"
) %in% f2$state_transition$invalidated))
stopifnot(isTRUE(certificate(f2)$passed))
```

The basis is a starting hint for `rho = 0.11`. No value, lock, cached action,
or certificate from `rho = 0.10` survives the transition.

## Serialization and portability

- Base RDS is the only persistence format.
- Basis matrices and the version-1 Lanczos start block are ordinary R numeric
  objects and can be portable.
- State portability requires every operator identity component and the public
  basis to be portable. Payload portability is recorded separately.
- A session-local callback identity always fails after restore, including under
  `basis_only`; users must provide explicit portable lineage and revision if
  cross-session continuation is intended.
- External pointers and native workspace addresses are invalid payloads.
- Integrity tokens cover all known state fields except the token itself and
  the derived memory record. Mutation therefore fails closed.
- Unknown schema versions never auto-upgrade. Unknown additive fields within
  schema version 1 round-trip unchanged.

## Retained memory

The state memory record has components:

```text
basis
method_state
cached_operator_actions
metadata
```

For version 1, `cached_operator_actions` is zero because no public adapter
retains them. `basis` includes both SVD sides when present. `method_state`
includes the fitted start block and its adapter metadata. `metadata` includes
identity, signature, provenance, serialization, and integrity records.

`r_bytes` is computed from retained R objects without double-counting shared
references. `native_bytes = 0` and `complete = TRUE` for version-1 basis and
Lanczos payloads because they own no native allocation. If a future adapter
owns unmeasured native state, `complete` becomes false and `total_bytes` is a
documented lower bound.

## Error contract

All state-validation errors inherit from `eigencore_restart_state_error`,
`error`, and `condition`, and carry `code`, `field`, `expected`, and `actual`.

| Code | Use |
|---|---|
| `unsupported_schema` | Unknown or non-integer state, basis, signature, or payload schema. |
| `coordinate_incompatible` | Dimension, dtype, structure, metric presence, coordinate ID, or basis-side mismatch. |
| `operator_incompatible` | Exact operator required but lineage or revision differs. |
| `method_incompatible` | Receiving route cannot consume the basis or the requested method payload. |
| `stale_method_state` | Payload adapter/version/token is unknown, corrupt, or stale. |
| `session_incompatible` | Opaque identity or nonportable required payload belongs to another session. |
| `corrupt_state` | Missing field, invalid type/value, non-finite basis, rank failure, or integrity mismatch. |

Warnings are not used for state invalidation. An accepted `auto` downgrade is
represented in `state_transition$reason`, where callers can inspect it without
parsing console output.

## Reconciliation with the 1.1 and early-1.2 source

This design was checked against the implementation at executable-plan commit
`3497ac9`.

- `R/warm_start.R` already defines original-coordinate basis validation,
  rank reduction, method-width fitting, and the invariant-subspace safety
  guard. Restart-state basis preparation extends that boundary; it does not
  bypass it.
- `R/solve_eigen.R` currently passes a start only to standard real Hermitian
  Lanczos routes. Generalized Lanczos, shift-invert, Arnoldi, LOBPCG, and dense
  fallback routes therefore remain unsupported until they gain explicit
  adapters.
- `R/solve.R` already freezes `restart_state`, `reuse`, and `retain_state` in
  the public `solve(plan)` signature while rejecting nondefault use. The
  implementation Mote replaces that temporary rejection without changing the
  signature.
- `R/reference_golub_kahan.R` and the native block Golub--Kahan code contain
  retained subspaces, locked vectors, projected recurrences, and cached
  actions. Their present ABI is richer than the version-1 public payload and
  is intentionally not serialized or exposed as stable method state.
- Result builders already reserve `state_transition`, `restart_state`, and
  memory fields. The implementation must replace provisional records with the
  validated types here while preserving legacy `restart` diagnostics.

No current source path provides a safe reason to broaden the version-1 adapter
table. Capability claims follow actual consumption and validation, not the
presence of an internal array with a suggestive name.

## Implementation and review gates

The implementation Mote must establish all of the following before this API is
described as supported:

1. constructors, validators, accessors, immutable copies, integrity tokens,
   memory records, and RDS behavior;
2. standard real Hermitian Lanczos basis reuse over dense double,
   `dgCMatrix`, and explicit/opaque callback operators;
3. exact same-operator start-block payload use and strict-mode failures;
4. changed-revision and changed-lineage basis-only transitions with fresh
   projections, residuals, ordering, convergence, and certificates;
5. generalized-metric rejection or correct B-reorthonormalization before any
   claim of generalized support;
6. zero-callback-work assertions for every rejected state;
7. independent callback work counts showing that current operators, not saved
   certificates, produced the new result;
8. adversarial mutation, stale payload, target/method/shape/dtype/structure,
   cross-session, serialization, and non-finite-basis tests;
9. retained-byte accounting with no native-memory or peak-RSS overclaim; and
10. continuation benchmarks comparing cold, public-basis, and exact-revision
    start-block runs at a common certified answer.

The integrated 1.2 gate must additionally verify that no unsupported SVD,
generalized, shift-invert, Arnoldi, LOBPCG, or dense-fallback route silently
accepts state. Adding one of those adapters changes the capability table and
requires route-specific numerical and work-accounting tests, but does not
require a new public state schema when it obeys this contract.
