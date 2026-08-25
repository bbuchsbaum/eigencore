# Whether the resolved plan consumes a user-supplied starting subspace.

The boundary is standard (non-generalized, non-transformed) real
Hermitian Lanczos: native dense double / dgCMatrix CSC paths, the native
matrix-free callback path selected by `lanczos(block > 1)`, and the
scalar matrix-free reference path selected by `lanczos(block = 1)`.
Every other dispatch is out of scope and must reject `initial_subspace`
rather than ignore it, densify, or borrow a production label.

## Usage

``` r
warm_start_plan_consumes_start(problem, plan)
```

## Value

`TRUE` if the plan consumes a user-supplied start block, else `FALSE`.
