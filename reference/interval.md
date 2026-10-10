# Target every eigenvalue in an interval.

`interval(a, b)` selects all eigenvalues of a Hermitian problem (or of a
symmetric-definite pencil `A x = lambda B x`) that lie in the closed
interval `[a, b]`. The number of such eigenvalues is not supplied: the
solver counts them first with Sylvester's law of inertia (see
[`eigen_count()`](https://bbuchsbaum.github.io/eigencore/reference/eigen_count.md)),
so `eig_partial(A, target = interval(a, b))` needs no `k`. A supplied
`k` is an upper bound and must be at least the count.

## Usage

``` r
interval(a, b)
```

## Arguments

- a, b:

  Interval end points, `a < b`. One of them may be infinite
  (`interval(-Inf, b)` selects every eigenvalue up to `b`); both
  infinite is the full spectrum, use
  [`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md).

## Value

An `eigencore_target` descriptor of kind `"interval"`.

## Details

**End points.** The interval is closed. Counting factors `A - t B` at
`t = a` and `t = b`; when a factorisation there is not reliable (an
eigenvalue numerically at the end point) the end point is moved outward
by a small recorded perturbation, so an eigenvalue inside that numerical
zero band is counted as inside. Eigenvalues whose residual bound
straddles an end point are reported in
`certificate(fit)$completeness$boundary_ambiguous` rather than dropped.

**Routes.** Dense matrices use LAPACK `dsyevr` with `RANGE = "V"`
(exact: the tridiagonal Sturm bisection selects the eigenvalues);
generalized dense pencils reduce with the Cholesky factor of `B` first.
Sparse matrices are counted and solved with \\LDL^T\\ shift-invert
Lanczos at the interval centre; wide intervals are split into slices by
inertia counts (spectrum slicing), each slice is solved at its own
centre (reusing the symbolic factorisation) and the slices are merged.
Every result is certified by residuals and by the inertia count:
`certificate(fit)$target_completeness` is `"inertia_verified"` when the
returned values provably are all eigenvalues in the interval (`"exact"`
for the dense LAPACK route).

Options: `eigencore.interval_dense_limit` (default 800) is the largest
sparse dimension solved by the dense route;
`eigencore.interval_slice_size` fixes the number of eigenvalues per
slice (default: chosen from the predicted factorisation cost, between 60
and 200). Per-slice diagnostics (bounds, counts, the free count at each
slice centre, operator applies, block size used for repeated
eigenvalues) are in `fit$interval`.

## Examples

``` r
A <- diag(c(1, 2, 2, 5, 7, 9))
fit <- eig_partial(A, target = interval(1.5, 6))
values(fit)
#> [1] 2 2 5
certificate(fit)$target_completeness
#> [1] "exact"
```
