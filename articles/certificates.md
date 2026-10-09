# Certificates: reading the numerical evidence

Every result returned by eigencore carries a **certificate**: a compact
record of the residual, backward-error scale, orthogonality check, and
tolerance used for that result. The certificate lets you inspect those
checks without recomputing the residual. It does not certify spectral
ordering or any property that its `certificate_type` does not name.

This vignette walks through the certificate’s fields, the three things
they measure, and what to do when one of them fails.

``` r

library(eigencore)
```

## Anatomy of a certificate

``` r

set.seed(1)
n <- 200
A <- crossprod(matrix(rnorm(n * n), n, n)) / n + diag(n)

fit  <- eig_partial(A, k = 5, target = largest())
cert <- fit$certificate
cert
#> eigencore certificate
#>   passed: TRUE 
#>   tolerance: 1e-08 
#>   type: residual_backward_error 
#>   norm bound: two_norm_lower_bound+identity_exact 
#>   norm source: ritz+identity 
#>   scale estimated: FALSE 
#>   max residual: 6.259907e-11 
#>   max backward error: 6.578934e-12 
#>   max orthogonality loss: 8.881784e-16 
#>   orthogonality tolerance: 1.490116e-08 
#>   orthogonality required: TRUE 
#>   target completeness: probed
```

The fields you act on, in order:

- **`passed`** — overall verdict. `TRUE` means every returned pair
  satisfies the checks required by this certificate type. The
  backward-error scale is exact or a lower bound on the matrix 2-norm,
  never an estimate, so a pass is never optimistic (see “Norm bounds”
  below).
- **`tolerance`** — the user-requested tolerance. Defaults to `1e-8`.
- **`max_residual`** — the worst absolute residual `||A v - lambda v||`
  (or `||A v - lambda B v||` for generalized problems) across the
  returned basis.
- **`max_backward_error`** — the worst residual divided by the labelled
  scale. This is the number that has to be below `tolerance` for a pair
  to be considered converged.
- **`max_orthogonality_loss`** — the worst entry of `|V* V - I|` (or
  `|V* B V - I|` in B-inner-product problems).
- **`failed_indices`** — which pairs failed. Empty when `passed = TRUE`.

Three more fields explain *how* the verdict was reached rather than what
it is: `norm_bound_type` says whether the `||A||_2` (and `||B||_2`)
value in the backward-error scale is exact or a lower bound,
`norm_source` says where it came from, and `norm_values` holds the
values used. `scale_is_estimate` is always `FALSE` for built-in
certificates. `certificate_type` and `notes` carry provenance text for
unusual states.
[`?certificate`](https://bbuchsbaum.github.io/eigencore/reference/certificate.md)
documents every field.

The `max_*` fields are summaries. The certificate also keeps the
underlying *per-pair* vectors (`cert$backward_error`, `cert$converged`),
so you can see at a glance whether the whole basis cleared the bar or
just barely missed:

![Stem plot of backward error for five eigenpairs on a log scale, all
falling well below the dashed tolerance
line.](certificates_files/figure-html/anatomy-bars-1.png)

Per-pair backward error for the five returned eigenpairs. The
`max_backward_error` field is simply the height of the tallest stem;
`passed` is TRUE because every stem clears the tolerance line.

## Three things a certificate measures

### 1. Residual

For a Hermitian eigenproblem `A v = lambda v`, the *absolute* residual
is `||A v_i - lambda_i v_i||`. For an SVD, both the left and right
relations matter, so the combined residual is

> sqrt( \|\|A v - sigma u\|\|^2 + \|\|A^T u - sigma v\|\|^2 ).

For a generalized SPD problem `A v = lambda B v`, the residual is
computed in the *original* coordinates: `||A v_i - lambda_i B v_i||`.
Eigencore never reports a residual computed in a transformed problem
space without saying so in `certificate_type`.

### 2. Backward error

Absolute residuals can be misleading when the operator’s norm is large
or small. Backward error is the residual divided by a scale that
captures “how big could a perturbation of `A` (and `B`) be that
*exactly* makes `(lambda, v)` an eigenpair?”. eigencore uses the
standard normwise definition in the matrix 2-norm:

> eta_i = \|\|r_i\|\| / ( (\|\|A\|\|\_2 + \|lambda_i\| \|\|B\|\|\_2)
> \|\|v_i\|\| ),

and for a singular triplet
`eta_i = sqrt(||A v - sigma u||^2 + ||A^T u - sigma v||^2) / ||A||_2`
(unit vectors). A pair “converges” when `eta_i <= tol`. The tolerance
you pass to `eig_partial(tol = ...)` is this backward-error tolerance,
not a residual tolerance.

`||A||_2` is rarely known exactly, so the denominator uses a value that
never exceeds it. A smaller denominator can only make the reported
`eta_i` larger than the true backward error, so a certificate that
passes would also pass with the exact norm. (Before eigencore 1.4 the
scale was the Frobenius norm, which can be up to `sqrt(rank)` times
larger than `||A||_2` and made the check correspondingly more lenient.)

### 3. Orthogonality

Iterative methods drift. After enough restarts the returned basis can
lose orthogonality even when each per-pair residual is small. The
certificate records `max_abs(V* V - I)` (or `max_abs(V* B V - I)`); its
acceptance tolerance is `max(tol, sqrt(eps))`. A passing residual should
be interpreted alongside the orthogonality field. If orthogonality loss
exceeds its tolerance, clustered or repeated eigenvalues may be
represented by nearly duplicate vectors even when individual residuals
are small.

A certificate evaluates the pairs the solver returned. It does not prove
that those pairs are the requested largest, smallest, or nearest
eigenvalues. When target identity is consequential, inspect the method
and target in the plan and compare against an independent oracle on a
tractable instance.

## Reading common certificate states

The single most useful habit is to picture the per-pair backward error
against the tolerance line. A passing certificate is “all stems below
the line”; a failing one is “at least one stem above it.” Here are the
same ten eigenpairs of `A`, computed two ways:

``` r

# Largest eigenvalues are well separated -> easy, converges fast.
fit_pass <- eig_partial(A, k = 10, target = largest())

# Smallest eigenvalues are densely clustered near 1 -> a tight maxit (restart)
# budget leaves them short of tolerance.
fit_fail <- eig_partial(A, k = 10, target = smallest(), maxit = 4)

c(largest_passed  = fit_pass$certificate$passed,
  smallest_passed = fit_fail$certificate$passed)
#>  largest_passed smallest_passed 
#>            TRUE           FALSE
```

![Two side-by-side stem plots of per-pair backward error. Left panel
shows ten blue points below the tolerance line; right panel shows ten
red points above the tolerance
line.](certificates_files/figure-html/pass-vs-fail-1.png)

Same matrix, same k, two targets. Left: the ten largest eigenpairs all
clear the tolerance (blue, passed). Right: the ten smallest stall above
it under a tight iteration budget (red, failed).

### Clean: every box ticked

``` r

fit_pass$certificate$passed
#> [1] TRUE
fit_pass$certificate$norm_bound_type
#> [1] "two_norm_lower_bound+identity_exact"
fit_pass$certificate$scale_is_estimate
#> [1] FALSE
```

This is the easy case. Every residual is below tolerance and
orthogonality is near machine precision. The scale is a lower bound on
`||A||_2`, here the largest computed Ritz value (`norm_source` is
`"ritz+identity"`), which for a largest-eigenvalue target is essentially
the exact norm.

### Failed: residual too large

The right-hand panel above is a genuine failure. The smallest
eigenvalues of `A` sit in a dense cluster just above 1, so they need
many more iterations than the well-separated largest ones. With
`maxit = 4` the solver runs out of restart budget before all pairs
converge, and the verdict flips. The `failed_indices` slot tells you
which Ritz pairs missed:

``` r

fit_fail$certificate$passed
#> [1] FALSE
fit_fail$certificate$failed_indices
#>  [1]  1  2  3  4  5  6  7  8  9 10
fit_fail$certificate$max_backward_error
#> [1] 0.004174378
```

You do not have to guess how far off it was, or whether it was inching
toward convergence. The solver records a `convergence_history`; plotting
the worst backward error per restart shows the failed run plateauing
above the line while a generous budget drives it underneath.

``` r

fit_ok <- eig_partial(A, k = 10, target = smallest(), maxit = 40)
fit_ok$certificate$passed
#> [1] TRUE
```

![Line plot on a log scale of backward error versus restart number. The
red maxit-4 curve plateaus above the tolerance line; the blue maxit-40
curve descends below
it.](certificates_files/figure-html/convergence-1.png)

Worst backward error per restart for the ten smallest eigenpairs. A
tight budget (red) stalls above the tolerance; a generous one (blue)
drives the error under the line and the certificate passes.

What to do: increase `maxit`, raise `tol`, or — when you suspect a
clustered spectrum — request a larger `k` (so the cluster is fully
covered by the returned basis) and slice afterwards. Here, lifting
`maxit` (the restart limit) from 4 to 40 is enough.

### Norm bounds: exact or a lower bound

`norm_bound_type` says what kind of `||A||_2` value scaled the backward
error, as `A` or `A+B` parts (`B = I` appears as `identity_exact`):

- `two_norm_exact` — the exact spectral norm: diagonal matrices
  (`max |d_i|`), a full computed Hermitian spectrum
  ([`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md)
  and dense solves that compute every eigenvalue), or a user-asserted
  `metadata$two_norm`.
- `two_norm_lower_bound` — a value that never exceeds `||A||_2`.

`norm_source` names where the value came from:

- `diagonal`, `full_spectrum`, `metadata`, `identity` — exact values.
- `column_norms` — the largest column norm `||A e_j||` (dense, sparse,
  tridiagonal storage; no operator applies).
- `applied_vectors` — the largest `||A x|| / ||x||` over the certified
  vectors, computed with the same operator applies as the residuals.
- `ritz` — the residual-corrected value
  `|lambda| ||B x|| / ||x|| - ||r|| / ||x||` when only residual norms
  are available.
- `frobenius_rank_bound` — `||A||_F / sqrt(min(m, n))` when an operator
  only carries Frobenius metadata.
- `lanczos` — a 24-step Lanczos estimate on `A` (Hermitian) or `A^T A`,
  verified as `||A y|| / ||y||` for its Ritz vector. It runs only when
  some pair fails with the cheaper bounds but could pass with a larger
  one, uses a fixed start vector (it never touches the R random-number
  stream), and is memoised per operator.

Because every source is a lower bound, matrix-free operators certify
like any other operator. A
[`linear_operator()`](https://bbuchsbaum.github.io/eigencore/reference/linear_operator.md)
with no norm metadata used to fall back to a stochastic Hutchinson
estimate and withhold `passed`; it now certifies with the bound from its
own vectors:

``` r

set.seed(2)
op <- linear_operator(
  dim = c(n, n),
  apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (A %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  },
  apply_adjoint = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (A %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  },
  structure = hermitian(),
  name = "matrix-free Hermitian wrapper"
)
# The matrix-free path runs the reference Lanczos solver, which needs a
# larger subspace than the default to resolve this clustered spectrum.
fit_mf <- eig_partial(op, k = 5, target = largest(),
                      method = lanczos(max_subspace = 100))
fit_mf$certificate$norm_bound_type
#> [1] "two_norm_lower_bound+identity_exact"
fit_mf$certificate$norm_source
#> [1] "applied_vectors+identity"
fit_mf$certificate$scale_is_estimate
#> [1] FALSE
fit_mf$certificate$passed
#> [1] TRUE
```

A lower bound can be weak — for example `||A x|| / ||x||` for the
*smallest* eigenvalues of an operator with no other information. Then
the reported backward error over-states the true one and a pair may fail
although it is accurate; the Lanczos refinement exists to close that gap
before a pair is reported as failed. What to do if a pass matters and
the bound is still weak: supply `metadata = list(two_norm = ...)` when
you know `||A||_2`, or wrap a built-in dense / `dgCMatrix` / `ddiMatrix`
matrix.

### Generalized SPD: B-orthogonality matters

For `A v = lambda B v`, the certificate’s orthogonality field is in the
**B-inner product**: `max_abs(V* B V - I)`. For this certificate type,
`passed = TRUE` requires both the residual threshold and B-orthogonality
of the returned vectors.

``` r

set.seed(4)
B <- diag(seq(1, 5, length.out = n))
fit_gen <- eig_partial(A, k = 5, target = largest(), B = B,
                       method = lobpcg(maxit = 200))
fit_gen$certificate$norm_bound_type
#> [1] "two_norm_lower_bound+two_norm_lower_bound"
fit_gen$certificate$max_orthogonality_loss
#> [1] 2.442491e-15
fit_gen$certificate$passed
#> [1] TRUE
```

If a B-orthogonality value comes back near machine precision, the
B-inner-product Cholesky-QR refinement inside the solver did its job. If
it comes back loose (say, `1e-4`), increase `maxit` or lower `tol` —
orthogonality loss is usually the first thing to surface in
ill-conditioned-B problems.

### Target completeness: the right set, not just right pairs

Residuals prove that each returned pair is an eigenpair. They cannot
prove that the returned *set* is the one you asked for. A single-vector
Krylov method sees, in exact arithmetic, only one direction of each
eigenspace, so it can miss a copy of a repeated eigenvalue: for the
spectrum `9, 9, 7, 7, 7, 5, ...` it may return `9, 9, 7, 7, 5`, and
every one of those five pairs certifies.

After a certified Hermitian Krylov solve with a
[`largest()`](https://bbuchsbaum.github.io/eigencore/reference/largest.md),
[`smallest()`](https://bbuchsbaum.github.io/eigencore/reference/smallest.md)
or
[`largest_magnitude()`](https://bbuchsbaum.github.io/eigencore/reference/largest_magnitude.md)
target, eigencore therefore runs a short **deflated complement probe**:
a few block Lanczos steps on the operator restricted to the orthogonal
complement of the returned vectors (B-orthogonal for generalized
problems), from a fixed-seed start that does not touch R’s random-number
stream. Any Ritz value of that restricted operator that lies beyond the
least-preferred returned value (by more than the residual and tolerance
margin) proves a more-preferred eigenvalue is missing. The solver then
repairs the result with a deflated complement solve, merges the two sets
by a Rayleigh-Ritz step, re-certifies and probes again.

``` r

d <- c(9, 9, 7, 7, 7, seq(5, 0.1, length.out = 55))
set.seed(1001)
p <- sample(60)
S <- Matrix::sparseMatrix(i = p, j = p, x = d)
bare <- eig_partial(S, k = 5, method = lanczos(completeness = "none"), seed = 1)
sort(values(bare), decreasing = TRUE)
#> [1] 9 9 7 7 5
bare$certificate$passed
#> [1] TRUE
fit_c <- eig_partial(S, k = 5, method = lanczos(), seed = 1)
sort(values(fit_c), decreasing = TRUE)
#> [1] 9 9 7 7 7
fit_c$certificate$target_completeness
#> [1] "repaired"
fit_c$certificate$completeness[c("steps", "operator_columns", "rounds")]
#> $steps
#> [1] 10
#> 
#> $operator_columns
#> [1] 124
#> 
#> $rounds
#> [1] 1
```

`certificate$target_completeness` is one of:

- `"probed"`: the probe found nothing beyond the returned edge;
- `"repaired"`: the probe found a missing value and the repaired set
  then probed clean;
- `"failed"`: an intruder remained after the repair budget
  (`options(eigencore.completeness_max_rounds = 3)`); `passed` is then
  `FALSE` even though every residual passes (`residual_passed` keeps
  that verdict, `target_passed` is `FALSE`);
- `"exact"`: a full-spectrum route (dense LAPACK, tridiagonal, analytic
  grid Laplacian) selected from every eigenvalue, so the set is complete
  by construction;
- `"not_checked"`: the probe did not apply (interior or both-ends
  targets, nonsymmetric problems, a residual certificate that already
  failed, `certify = FALSE`, or `completeness = "none"`).

The probe is cheap (by default at most 8 steps of block size 2, set by
the options `eigencore.completeness_probe_steps` and
`eigencore.completeness_probe_block`) and is switched off with
`lanczos(completeness = "none")` or
`options(eigencore.target_completeness = "none")`. It is a probabilistic
check: it can prove a set incomplete, not complete — an intruder whose
Ritz value a short run does not resolve goes unnoticed. The
deterministic answer is an eigenvalue-counting (inertia, `LDL'`)
certificate, planned for a later release.

## Comparing eigencore certificates to RSpectra diagnostics

[`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html)
returns `nconv` and `niter` but does not return residuals, backward
errors, or an orthogonality measure. The eigencore shim
([`eigencore::eigs_sym()`](https://bbuchsbaum.github.io/eigencore/reference/eigs_sym.md))
returns the same RSpectra-shaped list with two added fields:

``` r

res <- eigs_sym(A, k = 5, which = "LA")
names(res)
#> [1] "values"      "vectors"     "nconv"       "niter"       "nops"       
#> [6] "certificate" "diagnostics"
res$certificate
#> eigencore certificate
#>   passed: TRUE 
#>   tolerance: 1e-08 
#>   type: residual_backward_error 
#>   norm bound: two_norm_lower_bound+identity_exact 
#>   norm source: ritz+identity 
#>   scale estimated: FALSE 
#>   max residual: 6.019888e-09 
#>   max backward error: 6.007181e-10 
#>   max orthogonality loss: 2.677869e-15 
#>   orthogonality tolerance: 1.490116e-08 
#>   orthogonality required: TRUE 
#>   target completeness: probed
```

Code already written against
[`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html)
ignores `certificate` and `diagnostics` silently; new code can opt in to
certified results without changing call sites.

## Cheat sheet

| You see | What it means | What to do |
|----|----|----|
| `passed = TRUE` | All checks required by this certificate type passed; the 2-norm scale was exact or a lower bound, so the backward error is not under-stated. | Use the values/vectors within the stated certificate scope. |
| `passed = FALSE`, `failed_indices` non-empty | Some pairs hit `maxit` before converging. | Increase `maxit`, raise `tol`, or request a wider `k`. |
| `passed = FALSE`, `norm_bound_type` is `two_norm_lower_bound`, residuals tiny | The norm lower bound may be weak (e.g. smallest targets of a matrix-free operator). | Supply `metadata$two_norm` if `||A||_2` is known, or use a built-in matrix class. |
| `max_orthogonality_loss` near `sqrt(eps)` but residuals tiny | Iterative drift; clustered eigenvalues at risk of duplicates. | Increase `maxit`; check whether there are repeated eigenvalues. |
| `failed_indices` includes the first returned pairs | Some leading returned pairs exceed the tolerance. | Inspect `convergence_history`; increase `maxit` if the errors are still declining. |
| `passed = FALSE`, `residual_passed = TRUE`, `target_completeness = "failed"` | Every pair is accurate, but a more-preferred eigenvalue (e.g. a missed copy of a repeated value) lies outside the returned set. | Use a block method with a block at least the multiplicity (`lanczos(block = )`), raise `eigencore.completeness_max_rounds`, or a dense solve. |
| `failed_indices` includes the last returned pairs | Some trailing returned pairs exceed the tolerance. | Inspect the spectral gap; a wider `k` may help when the target boundary cuts through a cluster. |

The certificate records which numerical checks were run, the scale used
for each backward error, and whether those checks passed. Use those
fields together with the solver plan and the requirements of your
downstream analysis.

## Where to go next

- [`vignette("sparse-pca")`](https://bbuchsbaum.github.io/eigencore/articles/sparse-pca.md)
  shows certificates for centered and scaled sparse operators in
  context, including operators that carry no norm metadata.
- [`vignette("eigencore")`](https://bbuchsbaum.github.io/eigencore/articles/eigencore.md)
  and
  [`vignette("generalized-eigenproblems")`](https://bbuchsbaum.github.io/eigencore/articles/generalized-eigenproblems.md)
  cover the workflows that produce the certificates read here.
