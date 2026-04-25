# eigencore Engineering Guardrails

This repository follows `vision.md` for goals, `prd.json` for product scope,
and `plan_v1.md` for engineering method. When they appear to conflict:

1. `vision.md` wins on goals.
2. `prd.json` wins on scope.
3. `plan_v1.md` wins on execution discipline.

## Non-negotiables

- Public solver paths must run native engine code or carry an honest planner
  label identifying the path as a reference/oracle fallback.
- R solver prototypes are reference code only. Keep them named
  `reference_*`, keep them unexported, and do not present them as production
  algorithms.
- Sparse inputs must not silently densify in solver paths.
- Certificates must carry one shared scale definition across dense and
  operator paths; changes to certificate semantics require tests that compare
  both paths.
- Every exported symbol is a compatibility commitment. Do not export a helper
  unless it has real public semantics and supporting tests.
- Any solver PR must update planner labels, certificate provenance, and the
  adversarial/reference test bank.

## Current Execution Order

Work through `plan_v1.md` milestones A-M. A-E are prerequisites to production
native solvers:

- A: layer split and public-surface contraction
- B: frozen `BlockOperator` ABI plus dense/CSC/diagonal native operators
- C: adversarial test bank
- D: typed certificates and full target taxonomy
- E: dense-fallback memory policy

Do not add new solver families before these are in place.
