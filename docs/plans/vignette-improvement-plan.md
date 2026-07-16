# eigencore vignette improvement plan

Date: 2026-07-16. Reviewed against the r-vignette-craft scorecard, on branch
`release/1.0.0` (4 vignettes: eigencore, generalized-eigenproblems,
certificates, benchmarks).

## Diagnosis

### The headline feature has no vignette

DESCRIPTION and README both lead with: *"the computation behind PCA on big
sparse data … Centered, scaled, and composed operators are solved through
native C++ kernels without forming dense matrices."* That is the package's
reason to exist over RSpectra/irlba — and **none** of `center()`,
`scale_cols()`, `scale_rows()`, `compose()`, `crossprod_operator()` appears in
any vignette. A reader sold by the README pitch clicks "Articles" and finds
nothing that shows it. This is the single biggest gap.

### Getting started teaches architecture before it gives a win

`eigencore.Rmd` opens with the five-step workflow list (operator → problem →
plan → solve → certificate) and walks `eigen_problem()`/`plan_solver()`/
`solve()` first; `eig_partial()` then appears in §2 with no introduction. The
reader is taught the plumbing before the tap. The high-level verbs
(`eig_partial()`, `svd_partial()`) should carry the first win; problem/plan/
solve is the "when you want control" layer.

### Maintainer language leaks into user docs

"For the V2 CRAN release, planner labels are part of the contract: promoted
paths are native and benchmark-backed…" (eigencore.Rmd) and "…keeps promoted
native paths, reference fallbacks, and V3 deferrals honest" (certificates.Rmd).
V2/V3, "promoted", "deferrals" are release-engineering vocabulary; a user has
no referent for any of it. Cut or translate.

### Every example is an abstract random matrix

All four vignettes compute on `crossprod(matrix(rnorm(...)))` or `diag()`.
Nothing is *about* anything. One realistic, consistent dataset (a simulated
sparse document–term or spectral-embedding matrix) would anchor the whole
suite and make the "why partial + certified" story land.

### Per-vignette notes

- **certificates.Rmd** — the strongest of the four: motivating opening,
  pass-vs-fail paired panels, convergence-history plot, cheat-sheet table.
  Needs only the jargon cut and a link from every other vignette.
- **generalized-eigenproblems.Rmd** — good task-first headers and honest
  routing table, but zero figures, no motivating context (where do
  generalized problems arise — LDA, CCA, normalized-Laplacian embeddings?),
  and 3×3 diagonal toys with no scale context.
- **benchmarks.Rmd** — 260 of 567 lines are hidden harness code; it runs
  `bench::mark()` at build time and then spends most of its prose explaining
  why the numbers don't mean anything ("not evidence of package rankings or
  scaling behavior"). A vignette that mostly disclaims itself is a pkgdown
  article, not a CRAN vignette. It also rebuilds on CRAN's machines, which is
  both slow and noisy. Title "Benchmark smoke tests" is internal QA language.
- **YAML/theme**: `toc_depth: 2.0` should be `2`; shipping `albers.css` (37K),
  `albers.js`, and a `fonts/` directory inside `vignettes/` risks an
  installed-size NOTE and duplicates what the pkgdown template already
  provides — verify what actually lands in `inst/doc`.

## Plan (ranked by reader benefit)

### 1. New flagship workflow vignette: `sparse-pca.Rmd` — "PCA on sparse data without densifying"

The missing centerpiece. One realistic dataset built once (e.g., a simulated
100k × 2k document–term `dgCMatrix` with planted low-rank structure), then:

1. The problem: centering a sparse matrix densifies it — show the memory
   arithmetic (one sentence + inline R, not a wall).
2. First win: `svd_partial(center(A), rank = 10)` — certified scree plot.
3. Scaling too: `scale_cols()`, composition with `compose()`.
4. What the planner did: `plan_solver()` on the composed operator — this is
   where the operator/plan layer gets introduced *in context*, doing real work.
5. `linear_operator()` for the fully matrix-free case, and the R-callback
   performance warning story.
6. Scale context sentence + pointer to certificates vignette.

Covers: `center`, `scale_cols`/`scale_rows`, `compose`, `crossprod_operator`,
`linear_operator`, `svd_partial`, `plan_solver`. This becomes the vignette the
README links to.

### 2. Restructure `eigencore.Rmd` (Get started)

- Open with the problem ("you need the top k of a matrix too large to
  eigendecompose — and you need to know the answer is right"), then a first
  win inside the first two chunks: `eig_partial(A, k = 5)` → print → spectrum
  plot. Keep the existing plot helpers; they're good.
- Demote operator/problem/plan/solve to one later section: "Under the hood:
  problems and plans" (or push most of it into the sparse-pca vignette per #1
  and link).
- Delete the "V2 CRAN release" paragraph.
- Keep the RSpectra shim section — it's a genuine on-ramp — and the SVD
  section.
- End with a routing table by reader intent: sparse PCA → `sparse-pca`,
  trust/validation → `certificates`, `B`-metric problems →
  `generalized-eigenproblems`, performance → benchmarks article.

### 3. Move `benchmarks.Rmd` out of CRAN vignettes

Make it a pkgdown-only article (`vignettes/articles/benchmarks.Rmd`, listed
under Articles but not built by CRAN). Keeps the honest content, drops the
CRAN build-time risk and the `bench`/`irlba` soft-dependency dance. Retitle
"Benchmarks: what we measure and what it means". If a CRAN-visible pointer is
wanted, one paragraph in the README suffices.

### 4. Lift `generalized-eigenproblems.Rmd` to the standard of the others

- Add a 2–3 sentence motivating opening: where `A v = λ B v` shows up
  (whitened PCA/CCA, LDA, normalized graph Laplacians).
- Add one figure using the shared spectrum-plot helper (extract the helpers
  to a single sourced file under `vignettes/` instead of copy-pasting them
  into three vignettes).
- One scale-context sentence after the 3×3 toys.
- Keep the decision table and the `include=FALSE` validation chunks (they're
  invisible and keep the prose honest), but mirror the key assertions in
  testthat so the contract lives in tests, not only in doc builds.

### 5. Polish pass on `certificates.Rmd`

- Cut the "V2 CRAN release" paragraph.
- Compress the anatomy field list: keep the six fields a user acts on
  (`passed`, `tolerance`, the three `max_*`, `failed_indices`); fold
  `certificate_type`/`norm_bound_type` details into the "Withheld" section
  where they matter (details live in `?certificate`).
- Otherwise leave it alone — it's the model the others should match.

### 6. Cross-cutting cleanups

- Shared plot helpers → one file, `source()`d from an `include=FALSE` chunk.
- `toc_depth: 2.0` → `2` in all four YAML headers.
- Audit vignette asset weight (css/js/fonts) against what albersdown injects;
  confirm no size NOTE in `R CMD check`.
- Every vignette ends with "where to go next" links; every mention of
  certificates links `vignette("certificates")`.

### Deferred (post-1.0)

- `advanced.Rmd`: interior eigenvalues (`nearest()`, `shift_invert()`),
  `both_ends()`, preconditioners, `randomized()` — currently zero vignette
  coverage but too much for the 1.0 window. The functions have help pages;
  a vignette can follow evidence of user demand.

## Suggested order of execution

1 (new flagship) → 2 (get started restructure) → 3 (benchmarks move, small) →
4, 5 (parallel, small) → 6 (sweep) → rebuild site, `R CMD check`, read the
four articles top-to-bottom in the rendered site before calling it done.
