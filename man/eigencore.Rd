\name{eigencore}
\alias{eigencore}
\alias{adjoint}
\alias{as_operator}
\alias{auto}
\alias{backward_error}
\alias{certificate}
\alias{center}
\alias{check_adjoint}
\alias{compose}
\alias{crossprod_operator}
\alias{diagnostics}
\alias{eig_partial}
\alias{eigen_problem}
\alias{eigs}
\alias{eigs_sym}
\alias{euclidean}
\alias{general}
\alias{golub_kahan}
\alias{hermitian}
\alias{lanczos}
\alias{largest}
\alias{largest_real}
\alias{largest_imaginary}
\alias{largest_magnitude}
\alias{left_vectors}
\alias{linear_operator}
\alias{nearest}
\alias{plan_solver}
\alias{randomized}
\alias{residuals}
\alias{right_vectors}
\alias{scale_cols}
\alias{scale_rows}
\alias{shift_invert}
\alias{smallest}
\alias{smallest_real}
\alias{smallest_imaginary}
\alias{smallest_magnitude}
\alias{symmetric_operator}
\alias{svd_partial}
\alias{svd_problem}
\alias{svds}
\alias{values}
\alias{vectors}
\title{Certified partial spectral computation}
\description{
Foundations for eigencore's certified partial eigenvalue and singular-value
computation API. This initial implementation provides operator, problem,
planner, result, certificate, and compatibility surfaces, with dense oracle
prototype solvers used until native iterative solvers are added.
}
\usage{
eig_partial(A, k, target = largest(), B = NULL, method = auto(),
  tol = 1e-8, maxit = NULL, vectors = TRUE, seed = NULL,
  certify = TRUE)

svd_partial(A, rank, target = largest(), method = auto(), tol = 1e-8,
  vectors = c("both", "left", "right", "none"), seed = NULL,
  certify = TRUE)

eigen_problem(A, metric = NULL, structure = NULL, target = largest(),
  transform = NULL)

svd_problem(A, domain = NULL, codomain = NULL, target = largest())

linear_operator(dim, apply, apply_adjoint = NULL, dtype = "double",
  structure = general(), name = NULL, metadata = list())

compose(A, B, name = NULL)
scale_rows(A, weights, name = NULL)
scale_cols(A, weights, name = NULL)
center(A, rows = FALSE, columns = TRUE, row_means = NULL,
  col_means = NULL, name = NULL)
symmetric_operator(A, validate = TRUE, tol = 1e-10)
crossprod_operator(A, name = NULL)
check_adjoint(A, trials = 20, tol = 1e-12, seed = NULL)

as_operator(x, ...)
adjoint(x, ...)
plan_solver(problem, ...)
certificate(x, ...)
diagnostics(x, ...)
values(x, ...)
vectors(x, ...)
left_vectors(x, ...)
right_vectors(x, ...)
residuals(x, ...)
backward_error(x, ...)

largest()
smallest()
largest_real()
smallest_real()
largest_imaginary()
smallest_imaginary()
largest_magnitude()
smallest_magnitude()
nearest(sigma)
auto()
lanczos(max_subspace = NULL, reorthogonalize = TRUE)
golub_kahan(max_subspace = NULL, reorthogonalize = TRUE)
randomized(oversample = 10, n_iter = 2, block = NULL,
  normalizer = c("qr", "lu", "none"), refine = TRUE)
shift_invert(sigma, solve = NULL, factorization = NULL)
general()
hermitian()
euclidean(dim, dtype = "double")

eigs(A, k, which = "LM", opts = list(), ...)
eigs_sym(A, k, which = "LA", opts = list(), ...)
svds(A, k, nu = k, nv = k, opts = list(), ...)
}
\arguments{
\item{A}{A matrix, Matrix object, or eigencore operator.}
\item{k, rank}{Number of eigenpairs or singular triplets requested.}
\item{target}{Target descriptor such as \code{largest()} or \code{smallest()}.}
\item{B, metric}{Optional SPD metric for generalized eigenproblems.}
\item{max_subspace}{Maximum subspace size for prototype Lanczos or Golub-Kahan solvers.}
\item{reorthogonalize}{Whether to use full reorthogonalization in prototype iterative solvers.}
\item{domain, codomain}{Optional vector-space descriptors for SVD problems.}
\item{method, transform}{Solver or spectral transform descriptor.}
\item{tol}{Requested certification tolerance.}
\item{maxit}{Optional iteration limit, reserved for native iterative solvers.}
\item{vectors}{Whether and which vectors to return.}
\item{seed}{Optional random seed.}
\item{certify}{Whether to compute certificates when enough vectors are available.}
\item{dim}{Operator dimensions.}
\item{apply, apply_adjoint}{Block apply functions.}
\item{dtype}{Operator scalar type.}
\item{structure}{Operator structure descriptor.}
\item{name}{Optional operator name.}
\item{metadata}{Optional implementation metadata.}
\item{weights}{Row or column scaling weights.}
\item{rows, columns}{Whether to center rows or columns.}
\item{row_means, col_means}{Explicit means for matrix-free centering.}
\item{validate}{Whether to validate adjoint consistency before marking a symmetric operator.}
\item{trials}{Number of randomized adjoint-check trials.}
\item{x, problem}{Objects consumed by generic accessors and planners.}
\item{sigma}{Shift or target value.}
\item{solve, factorization}{Optional user-supplied shift-invert solve state.}
\item{oversample, n_iter, block, normalizer, refine}{Randomized SVD controls.}
\item{which, opts, nu, nv, ...}{Compatibility-shim arguments.}
}
\value{
Problem constructors return problem objects. Solver functions return structured
result objects with values, vectors, diagnostics, a solver plan, and a
certificate. Compatibility shims return RSpectra-like lists with additional
certificate and diagnostics fields.
}
\examples{
A <- diag(c(4, 1, 3, 2))
fit <- eig_partial(A, k = 2)
values(fit)
certificate(fit)

sf <- svd_partial(A, rank = 2)
values(sf)
certificate(sf)
}
