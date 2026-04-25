# eigencore Vision

eigencore will be the fastest and most trustworthy path to a partial spectrum in R, and it must be unambiguously better than the Python equivalents users would otherwise reach for.

The package exists to make large-scale eigenvalue and singular-value computations fast, accountable, elegant, and extensible. R users should not have to choose between convenient interfaces, native numerical performance, mathematical clarity, and clear evidence that the result actually converged. eigencore will provide all of these through a native spectral engine built around operators, spaces, metrics, transforms, solvers, and certificates.

The first version should feel familiar enough for RSpectra users to migrate quickly, but its foundation is more ambitious. It should also set a higher bar than SciPy-style sparse linear algebra by being faster for the intended workloads, clearer in its mathematical model, and more explicit about correctness. Instead of treating spectral computation as a matrix-in, vectors-out routine, eigencore will model the mathematical problem directly: standard and generalized eigenproblems, rectangular SVD maps, weighted spaces, shift-invert transforms, preconditioners, and certified residuals.

eigencore is not trying to be merely another binding to a fast eigensolver.
PRIMME and RSpectra are important performance reference points, and eigencore's
end state is to beat both on its core workloads by time-to-certified-answer.
The claim to being better, however, must come from the whole product contract:
safer defaults, inspectable planning, certified numerical evidence,
operator-native workflows, and an R interface that makes trustworthy spectral
computation easier to build into packages and analyses. If eigencore is faster
but less explainable, it has missed the point. If it is more explainable but
permanently slower than PRIMME or RSpectra on the core workloads, it has also
missed the point.

The core product promise is certification by default. Every returned eigenpair or singular triplet should carry residuals, backward-error estimates, orthogonality diagnostics, convergence metadata, and an inspectable solver plan. A result should not merely contain numbers; it should explain whether those numbers are trustworthy.

Planner honesty is part of that promise. Every result should say what actually
ran: native production code, scalar staging code, dense fallback, reference
oracle, external plugin, or unsupported path. Sparse inputs must not silently
become dense; unsupported targets must not be silently remapped; approximate
answers must not be dressed as certified production results. The solver plan is
not decoration. It is part of the numerical evidence.

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

The benchmark that matters is time to certified answer. eigencore should be
faster than RSpectra, PRIMME, and Python equivalents on its core workloads, but
never by hiding weaker convergence or poorer diagnostics. If a result is slower
because it is more rigorously certified, the certificate should make that
tradeoff visible during development. In the finished system, stronger
certification and faster runtime should converge: if eigencore is slower without
stronger correctness, that is a performance failure; if it remains slower even
with stronger correctness on the workloads it claims to own, that is also a
product failure.

This makes validation and benchmarking first-class product features. Comparisons
against RSpectra, PRIMME, irlba, SciPy, and randomized methods should report not
only elapsed time, but residuals, backward error, orthogonality loss, convergence
status, memory use, target semantics, and certificate provenance. The goal is
not fastest numbers; it is fastest trustworthy numbers.

The long-term vision is a mathematically expressive spectral computation system for R: native at its core, matrix-free by design, extensible to harder solvers and larger backends, and disciplined about numerical evidence. V1 must be designed with V2 and V3 extension points in mind: new solvers, spectral transforms, preconditioners, matrix ecosystems, external backends, GPU kernels, distributed memory, and PRIMME or SLEPc plugins should plug into the same problem, operator, solver, and certificate model rather than replace it.

eigencore should become the R package users reach for when they need not just some eigenvalues or singular values, but a reliable partial spectrum they can inspect, trust, and build on.
