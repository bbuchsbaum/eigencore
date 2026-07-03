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

An `eigencore_plan` list describing the requested problem, chosen method
label, target, planner reasons, fallback label, and control metadata
used by solver dispatch.

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
