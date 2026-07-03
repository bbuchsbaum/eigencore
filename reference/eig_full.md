# Compute a full dense eigendecomposition

`eig_full()` is the dense/full eigencore surface for standard and
generalized eigenproblems. Sparse and operator inputs are not silently
densified.

## Usage

``` r
eig_full(
  A,
  B = NULL,
  structure = NULL,
  vectors = TRUE,
  tol = 1e-08,
  allow_dense_fallback = c("auto", "never", "always"),
  ...
)
```

## Arguments

- A:

  Base dense square matrix.

- B:

  Optional base dense square matrix for generalized problems.

- structure:

  Optional structure descriptor. Use
  [`general()`](https://bbuchsbaum.github.io/eigencore/reference/general.md)
  to force the general-pencil path for dense generalized inputs.

- vectors:

  Whether to return right eigenvectors.

- tol:

  Certification tolerance.

- allow_dense_fallback:

  Reserved dense fallback policy. Sparse/operator inputs still fail
  unless a future issue explicitly opens an opt-in dense fallback
  contract for `eig_full()`.

- ...:

  Reserved for future options.

## Value

An `eigencore_eigen_result`. Dense general-pencil results additionally
carry `alpha`, `beta`, `classification`, `classification_policy`,
`left_vectors` (left generalized eigenvectors satisfying
`w^H A = lambda w^H B`), and `conditioning`. For real pencils the
decomposition runs through the expert LAPACK driver `DGGEVX` with
balancing, so `conditioning` contains reciprocal condition numbers
`rconde`/`rcondv` and the balanced pencil norms `abnrm`/`bbnrm`. Complex
pencils use `ZGGEV` (R's bundled LAPACK subset has no `ZGGEVX`), so they
return left vectors but `conditioning$available` is `FALSE`.

## Examples

``` r
A <- diag(c(1, 4, 9))
B <- diag(c(1, 2, 3))
fit <- eig_full(A, B = B)
values(fit)
#> [1] 1 2 3
certificate(fit)$passed
#> [1] TRUE

# Force the dense general-pencil path when alpha/beta diagnostics matter.
pencil <- eig_full(A, B = B, structure = general())
pencil$alpha
#> [1] 1+0i 4+0i 9+0i
pencil$beta
#> [1] 1+0i 2+0i 3+0i
pencil$classification
#> [1] "finite" "finite" "finite"
```
