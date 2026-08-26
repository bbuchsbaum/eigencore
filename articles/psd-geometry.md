# Work with singular PSD geometry

``` r

library(eigencore)
library(Matrix)
```

A positive-semidefinite matrix `K` can measure energy even when it is
singular. The singular directions have zero energy, so `K` is a seminorm
on the original coordinates and a genuine metric only after restricting
to `image(K)` or identifying points that differ by `null(K)`.

Eigencore 1.3 turns that geometry into a certified factor. You begin
with four facts:

1.  the source shape and representation;
2.  the real symmetric PSD form `K`;
3.  the action you need, such as a form, root, projector, or
    pseudoinverse; and
4.  the evidence the representation can supply.

The factor validates those facts once, records any repair, and fails
before a fallback when the requested action is not certified.

## Certify one singular dense form

This four-coordinate example has two positive directions and a
two-dimensional null space. The orthogonal rotation makes the null space
nontrivial in the original coordinate system.

``` r

Q <- matrix(c(
  1,  1,  1,  1,
  1, -1,  1, -1,
  1,  1, -1, -1,
  1, -1, -1,  1
), 4, 4, byrow = TRUE) / 2

K <- Q %*% diag(c(9, 1, 0, 0)) %*% t(Q)
K_factor <- psd_factor(K)
K_factor
#> eigencore certified PSD factor
#>   representation: dense_spectral 
#>   dimension: 4 
#>   numerical rank: 2 
#>   numerical nullity: 2 
#>   fidelity: repaired_with_defect 
#>   certificate passed: TRUE
c(rank = psd_rank(K_factor), nullity = psd_nullity(K_factor))
#>    rank nullity 
#>       2       2
```

[`psd_rank()`](https://bbuchsbaum.github.io/eigencore/reference/psd_rank.md)
and the default
[`psd_nullity()`](https://bbuchsbaum.github.io/eigencore/reference/psd_nullity.md)
are numerical claims under the factor’s recorded tolerance. A generic
floating-point eigendecomposition does not prove an exact-arithmetic
nullity theorem, so `psd_nullity(K_factor, "algebraic")` is deliberately
unavailable. Structural constructors can make that stronger claim when
their proof supports it.

Inspect the manifest before asking reusable code for an optional action:

``` r

psd_capabilities(K_factor)
#> eigencore PSD capabilities (dense_spectral)
#>   available: form, sqrt, inverse_sqrt, pseudoinverse, image_projector, null_projector, reduction, lift, strict_solve, gram, orthonormalize, reduced_operator, numerical_rank, numerical_nullity, serialization, cache_reuse 
#>   unavailable: algebraic_nullity
certificate(K_factor)
#> eigencore PSD certificate
#>   scope: source_validation_and_factor_actions 
#>   passed: TRUE 
#>   representation: dense_spectral 
#>   fidelity: repaired_with_defect 
#>   symmetry defect: 0 
#>   repair defect: 2.945755e-15 
#>   source/action defect: 2.945755e-15
```

The certificate distinguishes the source from the action that will run.
Tiny accepted negative eigenvalues, numerical-null positive eigenvalues,
and symmetry averaging are repairs, not erased history: their categories
and defects remain in `K_factor$spectrum` and `certificate(K_factor)`.

## Apply roots, projectors, and the pseudoinverse

All complete factors expose the same six square actions.
[`psd_apply()`](https://bbuchsbaum.github.io/eigencore/reference/psd_apply.md)
accepts a vector or a block with one row per original coordinate.

``` r

X <- matrix(c(
  1, 2, 3, 4,
  4, 3, 2, 1
), 4, 2)

KX <- psd_apply(K_factor, X, "form")
K_half_X <- psd_apply(K_factor, X, "sqrt")
image_X <- psd_apply(K_factor, X, "image_projector")
null_X <- psd_apply(K_factor, X, "null_projector")

stopifnot(
  isTRUE(all.equal(KX, K %*% X, tolerance = 1e-12)),
  isTRUE(all.equal(image_X + null_X, X, tolerance = 1e-12))
)
```

The inverse square root and pseudoinverse invert only retained positive
modes. They are not ordinary inverses of singular `K`:

``` r

K_half <- psd_apply(K_factor, diag(4), "sqrt")
K_ihalf <- psd_apply(K_factor, diag(4), "inverse_sqrt")
K_plus <- psd_apply(K_factor, diag(4), "pseudoinverse")
P <- psd_apply(K_factor, diag(4), "image_projector")

stopifnot(
  isTRUE(all.equal(K_half %*% K_half, K, tolerance = 1e-12)),
  isTRUE(all.equal(K_half %*% K_ihalf, P, tolerance = 1e-12)),
  isTRUE(all.equal(K %*% K_plus, P, tolerance = 1e-12))
)
```

Use
[`psd_operator()`](https://bbuchsbaum.github.io/eigencore/reference/psd_operator.md)
when an eigencore solver or another operator consumer needs the action
without a dense materialization. The operator and direct interfaces
represent the same certified map:

``` r

sqrt_operator <- psd_operator(K_factor, "sqrt")
sqrt_operator
#> <eigencore operator>
#>   name: certified PSD sqrt 
#>   dim: 4 x 4 
#>   dtype: double 
#>   structure: hermitian
max(abs(
  sqrt_operator$apply(X) - psd_apply(K_factor, X, "sqrt")
))
#> [1] 0
```

## Move explicitly to the metric space

[`psd_reduce()`](https://bbuchsbaum.github.io/eigencore/reference/psd_reduce.md)
maps original-coordinate columns into canonical image coordinates.
[`psd_lift()`](https://bbuchsbaum.github.io/eigencore/reference/psd_lift.md)
returns their minimum-Euclidean-norm representatives. Null-space
contamination disappears under reduction, while every K-Gram is
preserved.

``` r

reduced_X <- psd_reduce(K_factor, X)
contaminated_X <- X + 3 * null_X

stopifnot(
  isTRUE(all.equal(
    psd_reduce(K_factor, contaminated_X),
    reduced_X,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    crossprod(reduced_X),
    psd_gram(K_factor, X),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    psd_lift(K_factor, reduced_X),
    image_X,
    tolerance = 1e-12
  ))
)
```

This is the right boundary for algorithms that require an honest inner
product.
[`psd_orthonormalize()`](https://bbuchsbaum.github.io/eigencore/reference/psd_orthonormalize.md)
performs the corresponding block operation and reports discovered rank,
conditioning, dropped directions, work, and a fresh postcondition
certificate.

``` r

K_basis <- psd_orthonormalize(K_factor, X, required_rank = 2)
K_basis
#> eigencore PSD-orthonormal block
#>   rank: 2 
#>   condition: 225 
#>   postcondition error: 1.566058e-14
psd_gram(K_factor, K_basis$basis)
#>               [,1]          [,2]
#> [1,]  1.000000e+00 -4.348374e-16
#> [2,] -1.110223e-16  1.000000e+00
certificate(K_basis)$passed
#> [1] TRUE
```

## Distinguish pseudoinverse application from solving an equation

The Moore–Penrose action `K^+ b` is defined for every finite `b`; it
silently drops `b`’s null component. A strict solution of `K x = b`
exists only when `b` lies in `image(K)`.
[`psd_solve()`](https://bbuchsbaum.github.io/eigencore/reference/psd_solve.md)
checks that compatibility per right-hand side and certifies the original
equation.

``` r

b <- c(1, 0, 0, 0)
psd_apply(K_factor, b, "pseudoinverse")
#> [1]  0.2777778 -0.2222222  0.2777778 -0.2222222
```

``` r

c(
  class = class(incompatible)[1],
  code = incompatible$code,
  field = incompatible$field
)
#>                            class                             code 
#> "eigencore_psd_incompatible_rhs"               "incompatible_rhs" 
#>                            field 
#>                              "B"
```

Projecting the right-hand side onto `image(K)` makes the equation
compatible:

``` r

b_image <- psd_apply(K_factor, b, "image_projector")
strict <- psd_solve(K_factor, b_image)
strict
#> eigencore strict PSD solve
#>   columns: 1 
#>   compatible: TRUE 
#>   certificate passed: TRUE
```

Relax the RHS tolerance only when that change is part of your model; the
factor records the default scale-relative policy, while an
operation-specific override is recorded in the solve certificate.

## Reduce a singular generalized eigenproblem

The existing generalized-eigen `B`/`metric=` surface remains SPD-only.
Do not pass a singular PSD form to it. Instead, reduce the problem
explicitly to the finite image space, solve the ordinary Euclidean
problem there, then lift the vectors.

``` r

A <- Q %*% diag(c(8, 3, 5, 2)) %*% t(Q)
A_reduced <- psd_reduced_operator(K_factor, A)
A_image <- A_reduced$apply(diag(psd_rank(K_factor)))
fit <- eig_full(A_image)

largest_index <- which.max(values(fit))
value <- values(fit)[largest_index]
vector_original <- psd_lift(K_factor, vectors(fit)[, largest_index])
residual <- A %*% vector_original -
  value * psd_apply(K_factor, vector_original, "form")

c(value = value, residual_norm = sqrt(sum(residual^2)))
#>         value residual_norm 
#>  3.000000e+00  3.607798e-15
stopifnot(certificate(fit)$passed, sqrt(sum(residual^2)) < 1e-9)
```

Here the certified image has dimension two, so materializing its
two-by-two operator gives
[`eig_full()`](https://bbuchsbaum.github.io/eigencore/reference/eig_full.md)
an exact norm for its certificate. For a large image, pass `A_reduced`
to
[`eig_partial()`](https://bbuchsbaum.github.io/eigencore/reference/eig_partial.md)
instead; an operator without exact norm metadata can return accurate
residuals while withholding an unqualified `passed` flag because its
certificate scale is estimated.

The caller still owns the scientific question: whether discarding the
null space preserves the requested estimand, and whether the retained
rank is adequate. Eigencore owns the numerical reduction, provenance,
and certificate.

## Choose a constructor whose evidence matches the action

Identity and diagonal sources use analytic complete paths. Dense square
sources use a complete symmetric eigendecomposition. A supplied dense
Gram factor uses a complete compact SVD. These paths expose roots,
projectors, rank, reduction, strict solve, block primitives, and reduced
operators.

``` r

identity_factor <- psd_identity(4)
diagonal_factor <- psd_factor(c(9, 1, 0, 0))

L_dense <- matrix(c(
  1, 0,
  0, 1,
  1, 1
), 3, 2, byrow = TRUE)
gram_dense <- psd_gram_factor(L_dense, orientation = "columns")

c(
  identity = psd_rank(identity_factor),
  diagonal = psd_rank(diagonal_factor),
  dense_gram = psd_rank(gram_dense)
)
#>   identity   diagonal dense_gram 
#>          4          2          2
```

Sparse Gram and graph-Laplacian constructors preserve sparse state and
make only structural claims. A sparse Gram factor proves `K = L L^T` or
`K = L^T L`, so form and Gram actions are available without an `n` by
`n` dense matrix. It does not label `L` as the principal square root.

``` r

L_sparse <- Matrix(L_dense, sparse = TRUE)
gram_sparse <- psd_gram_factor(L_sparse, orientation = "columns")
psd_capabilities(gram_sparse)
#> eigencore PSD capabilities (gram_sparse)
#>   available: form, gram, serialization, cache_reuse 
#>   unavailable: sqrt, inverse_sqrt, pseudoinverse, image_projector, null_projector, reduction, lift, strict_solve, orthonormalize, reduced_operator, numerical_rank, numerical_nullity, algebraic_nullity
psd_apply(gram_sparse, diag(3), "form")
#>      [,1] [,2] [,3]
#> [1,]    1    0    1
#> [2,]    0    1    1
#> [3,]    1    1    2
```

A sparse graph Laplacian proves positive semidefiniteness through its
graph structure and proves algebraic nullity by counting connected
components. Numerical rank, roots, projectors, and reduction remain
unavailable because the constructor never computes a complete spectrum.

``` r

L_path <- bandSparse(
  5, 5,
  k = c(-1, 0, 1),
  diagonals = list(rep(-1, 4), c(1, 2, 2, 2, 1), rep(-1, 4))
)
laplacian <- psd_laplacian(L_path)

psd_nullity(laplacian, type = "algebraic")
#> [1] 1
psd_apply(laplacian, rep(1, 5), "form")
#> [1] 0 0 0 0 0
```

``` r

c(class = class(rank_error)[1], code = rank_error$code)
#>                               class                                code 
#> "eigencore_psd_incomplete_evidence"               "incomplete_evidence"
```

The admitted release manifest is:

| Representation | Certified form/Gram | Root, pseudoinverse, projectors | Reduction, solve, block | Rank/nullity evidence | Dense `n` by `n` state |
|:---|:--:|:--:|:--:|:---|:--:|
| identity | yes | yes | yes | numerical and algebraic complete | no |
| diagonal | yes | yes | yes | numerical and algebraic complete | no |
| dense spectral | yes | yes | yes | numerical complete; algebraic unavailable | yes |
| dense Gram | yes | yes | yes | numerical complete; algebraic unavailable | yes |
| sparse Gram | yes | no | no | unavailable | no |
| sparse Laplacian | yes | no | no | algebraic component nullity only | no |
| generic sparse symmetric | rejected | no | no | unavailable | no |
| opaque callback | rejected without probing | no | no | unavailable | no reusable factor |

[`psd_capabilities()`](https://bbuchsbaum.github.io/eigencore/reference/psd_capabilities.md)
is the runtime authority. An unavailable cell raises a typed
`eigencore_psd_error` with a stable `code` and `capability`; it does not
densify or run an undocumented approximation.

## Understand tolerances, ownership, and persistence

All default thresholds use the admitted matrix’s exact Frobenius scale:

``` text
threshold = absolute tolerance + relative tolerance * ||K||_F.
```

No `max(1, ||K||_F)` floor is inserted, so zero-absolute-tolerance
classification is invariant under finite positive rescaling. Symmetry,
positivity, rank, and RHS compatibility have distinct typed tolerances
because they answer different questions. Customize them with
[`psd_policy()`](https://bbuchsbaum.github.io/eigencore/reference/psd_policy.md)
and
[`psd_tolerance()`](https://bbuchsbaum.github.io/eigencore/reference/psd_tolerance.md),
not unnamed scalar cutoffs.

### How does a rank tolerance change the factor?

Consider a diagonal source with one large mode, one small positive mode,
and a tiny accepted negative value. The default rank policy treats the
small positive mode as numerical null. A finer, explicitly recorded rank
tolerance retains it:

``` r

classification_source <- c(4, 1e-9, -1e-15)
default_factor <- psd_factor(classification_source)
fine_policy <- psd_policy(
  rank = psd_tolerance(rel = 1e-12)
)
fine_factor <- psd_factor(classification_source, policy = fine_policy)

data.frame(
  policy = c("default", "finer rank tolerance"),
  rank_threshold = c(
    certificate(default_factor)$thresholds$rank,
    certificate(fine_factor)$thresholds$rank
  ),
  retained_rank = c(psd_rank(default_factor), psd_rank(fine_factor))
)
#>                 policy rank_threshold retained_rank
#> 1              default   5.960464e-08             1
#> 2 finer rank tolerance   4.000000e-12             2
```

[`psd_spectrum()`](https://bbuchsbaum.github.io/eigencore/reference/psd_spectrum.md)
keeps the admitted source separate from the repaired spectrum that
powers factor actions:

``` r

data.frame(
  source = psd_spectrum(default_factor),
  action = psd_spectrum(default_factor, repaired = TRUE)
)
#>   source action
#> 1  4e+00      4
#> 2  1e-09      0
#> 3 -1e-15      0
```

The negative entry remains visible in `source` and is zero in `action`.
The small positive entry is also zero under the default rank policy.
Changing that policy changes the numerical rank claim; it does not erase
the original spectrum or convert numerical rank into an algebraic
theorem.

### What survives an RDS round trip?

Every 1.3 factor is an immutable snapshot.
[`operator_identity()`](https://bbuchsbaum.github.io/eigencore/reference/operator_identity.md),
[`work()`](https://bbuchsbaum.github.io/eigencore/reference/work.md),
and
[`retained_bytes()`](https://bbuchsbaum.github.io/eigencore/reference/retained_bytes.md)
expose lineage, logical work, and retained memory. Base
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
is the persistence format; integrity and schema checks run before a
restored factor acts. Live and opaque factors are not admitted in 1.3,
even when a callback carries an operator ID and revision.

``` r

factor_path <- tempfile(fileext = ".rds")
saveRDS(K_factor, factor_path)
restored_factor <- readRDS(factor_path)

data.frame(
  portable = operator_identity(restored_factor)$portable,
  same_identity = identical(
    operator_identity(restored_factor),
    operator_identity(K_factor)
  ),
  same_action = isTRUE(all.equal(
    psd_apply(restored_factor, X, "form"),
    psd_apply(K_factor, X, "form")
  )),
  work_complete = work(restored_factor)$complete,
  retained_bytes = as.numeric(retained_bytes(restored_factor))
)
#>   portable same_identity same_action work_complete retained_bytes
#> 1     TRUE          TRUE        TRUE          TRUE          51504
```

The byte count describes the R object retained by this factor, not peak
memory during construction. A restored factor acts only after its
recorded identity, schema, and integrity token agree with its state.

``` r

callback_calls <- 0L
opaque <- linear_operator(
  dim = c(4, 4),
  apply = function(X, ...) {
    callback_calls <<- callback_calls + 1L
    X
  },
  dtype = "double",
  structure = hermitian(),
  name = "declared Hermitian callback"
)
opaque
#> <eigencore operator>
#>   name: declared Hermitian callback 
#>   dim: 4 x 4 
#>   dtype: double 
#>   structure: hermitian
```

``` r

c(
  class = class(opaque_error)[1],
  code = opaque_error$code,
  callback_calls = callback_calls
)
#>                               class                                code 
#> "eigencore_psd_incomplete_evidence"               "incomplete_evidence" 
#>                      callback_calls 
#>                                 "0"
```

This boundary keeps evidence separate from assertion: Hermitian metadata
says how an operator is intended to behave, but it is not a PSD proof, a
complete spectrum, or a certified square root.
