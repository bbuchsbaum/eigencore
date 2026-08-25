# Implicit normal-equations (Gram) partial SVD.

Runs the production block thick-restart Lanczos on the smaller-side
normal operator (\\A^T A\\ or \\A A^T\\) without materializing the Gram
matrix, then recovers the opposite singular factor and certifies the
triplets with the exact two-sided residual in original coordinates. This
removes the explicit-Gram memory/dimension caps: cost per operator
application is one forward and one adjoint apply of \\A\\.

## Usage

``` r
native_implicit_gram_svd(
  op,
  rank,
  target = largest(),
  tol = 1e-08,
  vectors = c("both", "left", "right", "none"),
  block = NULL,
  max_subspace = NULL,
  max_restarts = 100L
)
```
