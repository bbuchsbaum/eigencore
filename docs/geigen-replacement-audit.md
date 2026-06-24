# geigen Replacement Audit

Date: 2026-06-24

This artifact records the evidence behind the geigen migration test bank. The
goal is use-case replacement, not namespace compatibility: eigencore keeps its
own public names and planner/certificate contract.

## Sources Checked

- `tools::CRAN_package_db()` for current CRAN reverse dependencies importing
  or depending on `geigen`.
- Current CRAN source tarballs for those packages, scanned for `geigen`,
  `gqz`, `gevalues`, and `gsvd` calls.
- Installed `geigen` 2.3 manual pages for `geigen()`, `gqz()`, `gevalues()`,
  and `gsvd()`.

The current CRAN reverse dependencies found in the package database were:
`AIDA`, `conicfit`, `CovTools`, `decp`, `HDTSA`, `influenceAUC`, `itdr`,
`iTensor`, `loadings`, `multiCCA`, `multivarious`, `PEIP`, `quadmatrix`,
`rMultiNet`, and `wideRhino`.

## geigen Manual Surface

| geigen surface | Manual pattern | eigencore migration |
| --- | --- | --- |
| `geigen(A, B, symmetric = TRUE)` | Dense symmetric/Hermitian definite generalized eigenproblem | `eig_full(A, B = ...)` for full dense decompositions, or `eig_partial(A, B = ...)` for partial SPD/Hermitian metrics |
| `geigen(A, B, symmetric = FALSE)` | Dense general matrix pencil, values plus homogeneous `alpha`/`beta` | `eig_full(A, B = ..., structure = general())` |
| `gqz(A, B, sort = "N")` | Dense QZ/generalized Schur decomposition | `generalized_schur(A, B, sort = NULL)` |
| `gevalues(gqz_result)` | Values extracted from QZ homogeneous coordinates | `values(qz)` for `alpha / beta`, `alpha_beta(qz)` for homogeneous coordinates and finite/infinite classification |
| `gsvd(A, B)` | Generalized singular value decomposition | Deferred to `generalized_svd(A, B, ...)` in `bd-01KVWKVRK08D6ZERE81ANGN3QG`; use `geigen::gsvd()` until that surface exists |

`geigen::gqz()` also documents sort modes such as `"+"`, `"-"`, `"S"`,
`"B"`, and `"R"`. eigencore currently exposes only no sorting, finite-first,
and infinite-first generalized Schur sorting. Unsupported sort predicates fail
explicitly instead of being silently interpreted as geigen-compatible modes.

## Reverse Dependency Usage

| Package | Observed geigen usage | Use class | eigencore replacement status |
| --- | --- | --- | --- |
| AIDA | `geigen::geigen(est_cov, diag(...), only.values = TRUE)` | Dense SPD/identity metric, values only | Covered by dense SPD full tests |
| CovTools | `geigen(A, B, symmetric = TRUE, only.values = TRUE)` | Dense SPD covariance distance | Covered by dense SPD full tests |
| decp | `geigen::geigen(sigma2, sigma1, symmetric = TRUE)$values` | Dense SPD covariance distance | Covered by dense SPD full tests |
| HDTSA | `geigen::geigen(K2, K1)` | Dense covariance/tensor time-series pencil, no explicit symmetric flag in scanned call | Covered by dense SPD and dense general-pencil tests |
| influenceAUC | `geigen(capE, capD, symmetric = TRUE)` | Dense SPD | Covered by dense SPD full tests |
| iTensor | `geigen(U_m, U_n, symmetric = TRUE)` | Dense SPD/ICA | Covered by dense SPD full tests |
| itdr | `geigen(KernelX, sigma)` | Dense kernel/sigma generalized eigenproblem | Covered by dense SPD and dense general-pencil tests |
| loadings | `geigen::geigen(..., symmetric = FALSE)` and `geigen::geigen(..., symmetric = TRUE)` | Dense general kernel PCA and dense SPD PLS-ROG | Covered by dense general-pencil and dense SPD tests |
| multiCCA | `geigen::geigen(A, B, symmetric = TRUE)` | Dense SPD/kernel CCA | Covered by dense SPD full tests |
| multivarious | `geigen::geigen(A, B, symmetric = geigen_symmetric)` | Dense SPD or dense general depending on caller flag | Covered by dense SPD and dense general-pencil tests |
| conicfit | `geigen(P, Q)` | Small dense general/indefinite pencil | Covered by dense general-pencil tests |
| quadmatrix | `geigen::geigen(E1, E2)` | Dense companion pencil for quadratic matrix equations | Covered by dense general-pencil tests |
| rMultiNet | `geigen(L, D, symmetric = TRUE)` | Graph Laplacian with diagonal degree matrix | Covered by sparse diagonal-B SPD tests and dense SPD full tests |
| PEIP | `geigen::gsvd(A, B)` plus GSVD helper accessors | GSVD | Deferred to the GSVD child issue |
| wideRhino | `geigen::gsvd(A, B)` plus GSVD helper accessors | GSVD | Deferred to the GSVD child issue |

No current CRAN reverse dependency scanned here directly called
`geigen::gqz()` or `geigen::gevalues()`. They remain in the migration bank
because they are documented geigen entry points and are covered by the QZ
surface tests.

## Test Bank Coverage

`tests/testthat/test-geigen-migration-bank.R` is package-check-safe and does not
depend on `geigen`. It covers:

- Dense real SPD `geigen(..., symmetric = TRUE)` migration to certified
  `eig_full(A, B = ...)`.
- Dense complex Hermitian SPD migration.
- Dense real general pencils with nontrivial complex eigenvalues and
  `alpha_beta()` evidence.
- Singular `B` / beta-zero classification into finite, infinite, and undefined
  values.
- Sparse diagonal-B SPD partial solves with `allow_dense_fallback = "never"`.
- Sparse transformed shift-invert generalized SPD solves with explicit native
  planner labels.
- `gqz()` / `gevalues()` migration to `generalized_schur()`, `values()`, and
  `alpha_beta()`.
- Deferred GSVD migration with no accidental `gsvd` or `generalized_svd`
  export.

`tests/testthat/test-geigen-parity.R` is intentionally excluded from package
tarballs because `geigen` is archived/optional. When `geigen` is locally
installed and `NOT_CRAN=true`, it compares eigencore against live
`geigen::geigen()` and `geigen::gqz()` oracles for dense SPD, dense general,
dense complex general, sparse SPD oracle comparisons, and QZ/gevalues parity.

## Sparse Migration Boundary

Sparse users should not migrate by wrapping sparse matrices in a full dense
solver. eigencore's full dense surfaces (`eig_full()` and
`generalized_schur()`) require base dense matrices and reject sparse/operator
inputs instead of silently densifying.

Sparse SPD/Hermitian generalized problems should use `eig_partial(A, B = ...)`
with an explicit method and, when densification would be unacceptable,
`allow_dense_fallback = "never"`. Result `method`, `plan`, `restart`, and
certificate fields carry the transform, factorization, or metric-solve
provenance.

Sparse indefinite or nonsymmetric general pencils are not covered by this
migration bank. They belong to the sparse general-pencil boundary issue, where
factorized or user-solve `B` semantics must be explicit in planner labels and
certificates.
