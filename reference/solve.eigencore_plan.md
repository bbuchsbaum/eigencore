# Execute a frozen eigencore solver plan.

`solve(plan)` validates and runs the problem, method, controls, policy,
and execution arguments captured by
[`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md).
It never invokes the planner unless `replan = TRUE`. Execution arguments
cannot be overridden through `...`.

## Usage

``` r
# S3 method for class 'eigencore_plan'
solve(
  a,
  b,
  restart_state = NULL,
  reuse = c("auto", "basis_only", "same_operator"),
  retain_state = c("none", "basis", "same_operator"),
  replan = FALSE,
  ...
)
```

## Arguments

- a:

  An executable `eigencore_plan`.

- b:

  Unused second argument reserved by the base
  [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html)
  generic.

- restart_state:

  Optional
  [`restart_state()`](https://bbuchsbaum.github.io/eigencore/reference/restart_state.md)
  built from a certified result. Version 1 accepts the public basis only
  on admitted standard real Hermitian Lanczos routes and validates it
  before applying the destination operator.

- reuse:

  Reuse policy for an optional restart state. `"basis_only"` always fits
  the public basis to the receiving method, `"same_operator"` requires a
  matching operator identity and compatible retained method state, and
  `"auto"` uses eligible same-operator method state or downgrades to the
  public basis. No policy reuses convergence, locked vectors, operator
  actions, residuals, or a certificate.

- retain_state:

  State-retention request for the returned result. `"basis"` retains
  only the certified public basis; `"same_operator"` also requests an
  eligible versioned method payload. Unsupported producer routes still
  return a basis-only state.

- replan:

  Whether to create and execute a new plan from the embedded problem
  under current policy. The original plan is not modified.

- ...:

  Must be empty; execution controls are frozen in the plan.

## Value

An `eigencore_eigen_result` or `eigencore_svd_result` containing the
exact plan that governed execution.
