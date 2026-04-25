# eigencore Mission

eigencore's mission is to give R users fast, native, elegant, and certified tools for partial eigenvalue and singular-value computation.

We serve package authors, applied researchers, and numerical power users who need reliable spectra from dense, sparse, generalized, and matrix-free problems. They should be able to compute the eigenpairs or singular triplets they need, understand whether the computation converged, and inspect the numerical evidence without leaving R.

The standard is non-negotiable: eigencore must be faster and unambiguously better than the Python equivalents on its core workloads. Better means not only lower time to answer, but stronger mathematical design, clearer diagnostics, cleaner extensibility, and more trustworthy certification.

To do this, eigencore will:

- provide native block-operator solvers for eigenproblems, generalized SPD/Hermitian eigenproblems, SVD, randomized SVD, LOBPCG, and shift-invert workflows;
- return residuals, backward-error estimates, orthogonality diagnostics, convergence metadata, and certificate status by default;
- avoid default normal-equation SVD when a true SVD method is the right numerical tool;
- make automatic solver choices inspectable through clear solver plans;
- support RSpectra-compatible entry points so existing code can migrate gradually;
- measure success by time to certified answer, not by raw speed alone;
- preserve clean extension points for V2 and V3 solvers, transforms, preconditioners, operators, and backends.

Our engineering standard is simple: every abstraction must either express the mathematics more clearly, make the computation faster, make the system more flexible, or make the result more trustworthy. Convenience matters, but not at the expense of numerical honesty, mathematical beauty, or efficiency.

eigencore exists so R users can build production workflows, research pipelines, and package APIs on spectral results that are fast enough to use and transparent enough to trust.
