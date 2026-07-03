# CLAUDE.md - eigencore

Use **mote** (`mote`) for issue tracking and agent coordination in this
repo, not beads. The repo-level `~/code/CLAUDE.md` references `bd`/beads
as the default, but eigencore uses mote — see `AGENTS.md` for the
workflow.

Quick pointers:

- `mote ready` — find actionable work
- `mote show <bd-id>` — inspect an issue
- `mote begin <bd-id> --paths <path>` — reserve paths before editing
- `mote done <bd-id>` / `mote handoff <bd-id> --to <actor>` — finish or
  hand off

Engineering guardrails (non-negotiables, plan order, layer split) live
in `AGENTS.md`. Read it before editing solver, certificate, planner, or
operator code.
