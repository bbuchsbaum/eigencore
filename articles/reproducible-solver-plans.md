# Freeze and replay a solver decision

Sometimes you need more than the answer to an eigenproblem: you need a
record of the solver decision that produced it. A problem describes what
to compute, whereas an executable plan freezes how eigencore will
compute it. Saving the plan lets you inspect, run, and replay the same
route without asking the planner to decide again.

This vignette follows one small Hermitian problem from construction
through execution and RDS persistence. It also shows how operator
identity prevents a saved workflow from quietly drifting to different
data.

``` r

library(eigencore)
```

## What does a plan add to a problem?

Start with a matrix and the numerical question: its two largest
eigenpairs. The problem records that question.
[`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md)
then chooses a route and freezes the controls, tolerance, operator
identity, and planner policy used in that choice.

``` r

A <- diag(seq(200, 1))
problem <- eigen_problem(A, target = largest())
plan <- plan_solver(problem, k = 2, tol = 1e-10)
plan
#> eigencore solver plan
#>   problem: eigen 
#>   requested: 2 
#>   target: largest 
#>   method: native scalar thick-restart Hermitian Lanczos 
#>   reasons:
#>    - structure: hermitian 
#>    - target: largest 
#>    - standard eigenproblem 
#>    - built-in dense operator has native block apply 
#>   controls:
#>    - block : 1 
#>    - max_subspace : 26 
#>    - max_restarts : 100 
#>    - check_stride : 0 
#>    - reorthogonalize : TRUE 
#>   fallback: dense oracle prototype if unsupported
```

Keep the distinction simple:

| Object | What it preserves | What happens when you solve it |
|:---|:---|:---|
| problem | operator, structure, and target | eigencore plans again under current policy |
| plan | problem plus a frozen route and execution contract | eigencore validates and runs that route |
| result | values, vectors, certificate, work, and governing plan | nothing is recomputed until you ask |

## How do you run and interpret the frozen decision?

`solve(plan)` validates the plan before applying the operator. The
result keeps both the planned route and the route that actually returned
the answer. They differ only when an allowed runtime fallback is needed.

``` r

fit <- solve(plan)
data.frame(
  planned = fit$planned_method,
  actual = fit$actual_method,
  fallback = fit$fallback_used,
  certified = certificate(fit)$passed
)
#>                                         planned
#> 1 native scalar thick-restart Hermitian Lanczos
#>                                          actual fallback certified
#> 1 native scalar thick-restart Hermitian Lanczos    FALSE      TRUE
```

If `fallback` is `TRUE`, inspect `fit$fallback_reason`. It names the
failed attempt, the route that returned the result, and a typed reason
such as a failed certificate. The certificate always describes the
returned result, not the abandoned attempt.

The work record uses logical units that can be compared without
pretending that every solver’s historical `matvecs` field means the same
thing:

``` r

w <- work(fit)
data.frame(
  complete = w$complete,
  solve_columns = w$operator_columns,
  certificate_columns = w$certification_operator_columns,
  iterations = w$iterations,
  restarts = w$restarts
)
#>   complete solve_columns certificate_columns iterations restarts
#> 1     TRUE           130                   2        110        4
```

`complete = TRUE` means every logical counter is known. When it is
`FALSE`, an `NA` counter means unknown, not zero; compare only the
phases that both routes report.

## How does the plan identify its operator?

Matrix-backed operators receive deterministic identities derived from
their dimensions, structure, and values. The revision therefore changes
when the matrix changes.

``` r

identity_record <- operator_identity(plan)$A
data.frame(
  operator_id = identity_record$operator_id,
  revision = identity_record$revision,
  origin = identity_record$origin,
  portable = identity_record$portable
)
#>                operator_id         revision  origin portable
#> 1 builtin-e695432add8b2285 e695432add8b2285 builtin     TRUE
```

For a callback, eigencore cannot infer whether two functions represent
the same external data. Supply provenance only when your application can
maintain that promise:

``` r

callback <- linear_operator(
  dim = dim(A),
  apply = function(X, alpha = 1, beta = 0, Y = NULL) {
    Z <- alpha * (A %*% X)
    if (is.null(Y) || beta == 0) Z else Z + beta * Y
  },
  structure = hermitian(),
  metadata = list(frobenius_norm = sqrt(sum(A^2))),
  operator_id = "diagonal-example",
  revision = "values-v1",
  portable = TRUE
)
operator_identity(callback)[c("operator_id", "revision", "portable")]
#> $operator_id
#> [1] "diagonal-example"
#> 
#> $revision
#> [1] "values-v1"
#> 
#> $portable
#> [1] TRUE
```

`operator_id` names the logical lineage; `revision` names its current
values and behavior. Marking a callback portable is an assertion by its
owner, not a numerical test performed by eigencore. Without explicit
provenance, a callback receives an opaque session-local identity and its
restored plan cannot run in a different R session.

## How do you save and replay a plan?

Plans use ordinary RDS persistence. A restored portable plan validates
its schema, route, controls, execution arguments, policy snapshot, and
operator identity before the operator is called.

``` r

plan_path <- tempfile(fileext = ".rds")
saveRDS(plan, plan_path)
restored_plan <- readRDS(plan_path)
restored_fit <- solve(restored_plan)

data.frame(
  portable = restored_plan$serialization$portable,
  same_route = identical(restored_fit$planned_method, fit$planned_method),
  same_values = isTRUE(all.equal(values(restored_fit), values(fit)))
)
#>   portable same_route same_values
#> 1     TRUE       TRUE        TRUE
```

Replay means the same validated route and execution contract, not
necessarily bit-for-bit equality from an iterative method. Compare
numerical answers under the requested tolerance and inspect the fresh
certificate on every result.

## What changes require a new plan?

Execution arguments are part of the frozen contract. An attempted
override fails with a typed error instead of producing a result under a
route that no longer matches the inspected plan.

``` r

data.frame(
  class = class(override_error)[1],
  code = override_error$code,
  field = override_error$field
)
#>                  class               code field
#> 1 eigencore_plan_error execution_override   tol
```

Use `solve(plan, replan = TRUE)` only when you deliberately want
eigencore to make a new decision from the embedded problem under current
policy. The original plan remains unchanged. If instead you want to
reuse a certified basis across related operators, use a restart state;
that is a different contract because the destination operator must
receive a fresh certificate.

## Where should you go next?

- [`vignette("warm-start-continuation")`](https://bbuchsbaum.github.io/eigencore/articles/warm-start-continuation.md)
  shows basis reuse across related operators and the exact-revision
  versus changed-revision transition.
- [`vignette("certificates")`](https://bbuchsbaum.github.io/eigencore/articles/certificates.md)
  explains the evidence carried by each result.
- [`vignette("sparse-pca")`](https://bbuchsbaum.github.io/eigencore/articles/sparse-pca.md)
  shows how planning exposes a non-densifying sparse route before
  computation begins.
- [`?plan_solver`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md),
  [`?operator_identity`](https://bbuchsbaum.github.io/eigencore/reference/operator_identity.md),
  and
  [`?work`](https://bbuchsbaum.github.io/eigencore/reference/work.md)
  document the complete records.
