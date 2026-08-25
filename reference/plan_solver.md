# Plan a solver for a problem.

Plan a solver for a problem.

## Usage

``` r
plan_solver(problem, ...)
```

## Arguments

- problem:

  Eigencore eigen or SVD problem object.

- ...:

  Additional planning arguments passed to methods.

## Value

A schema-versioned executable `eigencore_plan` containing the problem,
requested size, original method descriptor, frozen route, canonical
controls and planner policy, execution arguments, operator identity,
serialization capability, and retained-memory metadata.

## Details

The eigen and SVD methods accept their usual request and method
arguments plus execution controls through `...`. Eigen plans freeze
`tol`, `maxit`, `vectors`, `certify`, `allow_dense_fallback`, and an
optional `initial_subspace`; SVD plans freeze `tol`, `vectors`,
`certify`, and `allow_dense_fallback`. Use `solve(plan)` to execute
those values or `solve(plan, replan = TRUE)` to make a fresh decision
under current policy.

## Examples

``` r
A <- diag(c(4, 3, 2, 1))
plan <- plan_solver(eigen_problem(A), k = 2)
plan$method
#> [1] "native dense Hermitian LAPACK fallback"
plan$reasons
#> [1] "structure: hermitian"                          
#> [2] "target: largest"                               
#> [3] "standard eigenproblem"                         
#> [4] "built-in dense operator has native block apply"
```
