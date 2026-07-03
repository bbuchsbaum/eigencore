# Define an eigenproblem.

Define an eigenproblem.

## Usage

``` r
eigen_problem(
  A,
  metric = NULL,
  structure = NULL,
  target = largest(),
  transform = NULL
)
```

## Arguments

- A:

  Matrix or operator defining the linear map.

- metric:

  Optional metric operator for generalized eigenproblems.

- structure:

  Optional structure descriptor; defaults to the operator structure.

- target:

  Eigencore target descriptor.

- transform:

  Optional transform method such as
  [`shift_invert()`](https://bbuchsbaum.github.io/eigencore/reference/shift_invert.md).

## Value

An `eigencore_eigen_problem` object containing the operator, optional
metric, structure, target, and transform metadata consumed by
[`plan_solver()`](https://bbuchsbaum.github.io/eigencore/reference/plan_solver.md)
and [`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html).

## Examples

``` r
A <- diag(c(4, 3, 2, 1))
P <- eigen_problem(A, target = largest())
fit <- solve(P, k = 2)
values(fit)
#> [1] 4 3
```
