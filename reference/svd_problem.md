# Define an SVD problem.

Define an SVD problem.

## Usage

``` r
svd_problem(A, domain = NULL, codomain = NULL, target = largest())
```

## Arguments

- A:

  Matrix or operator defining the rectangular linear map.

- domain:

  Optional domain-space descriptor.

- codomain:

  Optional codomain-space descriptor.

- target:

  Eigencore singular-value target descriptor.

## Value

An `eigencore_svd_problem` object containing the operator, domain and
codomain descriptors, and singular-value target consumed by
[`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md)
and [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html).

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(40), 8, 5)
S <- svd_problem(X, target = largest())
fit <- solve(S, rank = 2)
values(fit)
#> [1] 3.309631 3.044806
```
