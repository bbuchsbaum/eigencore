# eigencore Vision

eigencore will be the fastest and most trustworthy path to a partial spectrum in R, and it must be unambiguously better than the Python equivalents users would otherwise reach for.

The package exists to make large-scale eigenvalue and singular-value computations fast, accountable, elegant, and extensible. R users should not have to choose between convenient interfaces, native numerical performance, mathematical clarity, and clear evidence that the result actually converged. eigencore will provide all of these through a native spectral engine built around operators, spaces, metrics, transforms, solvers, and certificates.

The first version should feel familiar enough for RSpectra users to migrate quickly, but its foundation is more ambitious. It should also set a higher bar than SciPy-style sparse linear algebra by being faster for the intended workloads, clearer in its mathematical model, and more explicit about correctness. Instead of treating spectral computation as a matrix-in, vectors-out routine, eigencore will model the mathematical problem directly: standard and generalized eigenproblems, rectangular SVD maps, weighted spaces, shift-invert transforms, preconditioners, and certified residuals.

The core product promise is certification by default. Every returned eigenpair or singular triplet should carry residuals, backward-error estimates, orthogonality diagnostics, convergence metadata, and an inspectable solver plan. A result should not merely contain numbers; it should explain whether those numbers are trustworthy.

eigencore will also be block-native from the start. The fundamental operation is applying an operator to dense blocks, not repeatedly crossing the R/native boundary for scalar callbacks. Built-in dense, sparse, diagonal, symmetric, centered, scaled, and composed operators should run through native kernels and avoid R-level iteration overhead in solver hot loops.

Elegance is not cosmetic. The design should expose the mathematical beauty of the problem while remaining brutally efficient in implementation. Good abstractions must make the solver more expressive, more optimizable, or easier to certify. If an abstraction obscures performance, weakens diagnostics, or blocks future solver families, it does not belong in the core.

V1 should establish eigencore as a disruptive replacement for common production uses of RSpectra by delivering:

- certified Hermitian and symmetric partial eigensolvers;
- compatibility-grade nonsymmetric eigensolvers;
- first-class generalized SPD/Hermitian eigenproblems;
- shift-invert with factorization-aware solving;
- true partial SVD based on Golub-Kahan bidiagonalization, not default normal equations;
- randomized SVD with honest residual certification and optional refinement;
- LOBPCG for symmetric and generalized problems;
- an inspectable automatic solver planner;
- RSpectra-compatible `eigs()`, `eigs_sym()`, and `svds()` shims.

The benchmark that matters is time to certified answer. eigencore should be faster than RSpectra and Python equivalents on its core workloads, but never by hiding weaker convergence or poorer diagnostics. If a result is slower because it is more rigorously certified, the certificate should make that tradeoff visible. If it is slower without stronger correctness, that is a performance failure.

The long-term vision is a mathematically expressive spectral computation system for R: native at its core, matrix-free by design, extensible to harder solvers and larger backends, and disciplined about numerical evidence. V1 must be designed with V2 and V3 extension points in mind: new solvers, spectral transforms, preconditioners, matrix ecosystems, external backends, GPU kernels, distributed memory, and PRIMME or SLEPc plugins should plug into the same problem, operator, solver, and certificate model rather than replace it.

eigencore should become the R package users reach for when they need not just some eigenvalues or singular values, but a reliable partial spectrum they can inspect, trust, and build on.
