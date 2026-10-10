# Parse the sparse tridiagonal source once for one planned solve. The cached vectors are immutable solver metadata: shifted diagonals are derived from them, never written back into the operator.

Parse the sparse tridiagonal source once for one planned solve. The
cached vectors are immutable solver metadata: shifted diagonals are
derived from them, never written back into the operator.

## Usage

``` r
shift_invert_prepare_tridiagonal(problem)
```
