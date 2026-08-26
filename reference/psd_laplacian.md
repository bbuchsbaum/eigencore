# Construct a structural sparse graph-Laplacian PSD factor

The constructor admits only CSC `dgCMatrix` and `dsCMatrix` inputs. It
validates symmetry, non-positive off-diagonals, and zero row sums
without dense conversion. Connected components certify algebraic
nullity; numerical rank and spectral actions remain unavailable without
separate spectral evidence.

## Usage

``` r
psd_laplacian(x, policy = psd_policy(), source = c("snapshot", "live"))
```

## Arguments

- x:

  A finite real-double `dgCMatrix` or `dsCMatrix` graph Laplacian.

- policy:

  A PSD policy.

- source:

  Source ownership.

## Value

A certified PSD factor when the representation is admitted.
