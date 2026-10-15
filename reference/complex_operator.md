# Complex operator from its real and imaginary parts.

Builds an eigencore operator for the complex matrix `re + 1i * im` from
two real matrices of the same dimensions, each dense or sparse (any
`Matrix` class). This is how complex *sparse* matrices enter eigencore:
the `Matrix` package has no complex sparse classes (`zgCMatrix` and
friends are not implemented), so a complex sparse Hermitian matrix is
passed as a symmetric real part and a skew-symmetric imaginary part.

## Usage

``` r
complex_operator(re, im = NULL, structure = NULL)
```

## Arguments

- re:

  Real part: a numeric matrix or a real `Matrix` object.

- im:

  Imaginary part with the dimensions of `re`; `NULL` means zero.

- structure:

  Optional structure descriptor. By default the operator is
  [`hermitian()`](https://bbuchsbaum.github.io/eigencore/reference/hermitian.md)
  when `re` is symmetric and `im` skew-symmetric (to a relative
  tolerance of `sqrt(.Machine$double.eps)`), else
  [`general()`](https://bbuchsbaum.github.io/eigencore/reference/general.md).

## Value

An `eigencore_operator` of `dtype = "complex"` whose apply takes real or
complex blocks.

## Details

A Hermitian operator (`re` symmetric, `im` skew-symmetric) is solved by
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
on the real symmetric `2n x 2n` embedding `[re -im; im re]`, kept
sparse: native block Lanczos, CHOLMOD LDL' shift-invert and inertia
counts all run on the real sparse embedding, and the eigenpairs are
mapped back to complex vectors and certified on the complex operator.

## Examples

``` r
re <- Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 2, 3, 2),
                           x = c(2, 3, 4, 1), dims = c(3, 3), symmetric = TRUE)
im <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = c(0.5, -0.5),
                           dims = c(3, 3))
op <- complex_operator(re, im)
fit <- eig_partial(op, k = 1, target = largest())
values(fit)
#> [1] 4
```
