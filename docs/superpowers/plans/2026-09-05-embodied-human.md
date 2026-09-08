# Embodied Human Implementation Plan

> Independent domain work uses the dispatching-parallel-agents skill. Prior user authorization covers source-constrained synthesis and autonomous implementation.

**Goal:** Extend the canonical body with visible physical breathing, functioning peripheral sensory/motor interfaces and fine anatomy, and make the workbench a fixed 3D workspace.
**Architecture:** Keep one body registry and native physiological owner. Add reproducible domain artifacts and tested runtime components; bind their outputs into canonical shared state and renderer.
**Spec:** docs/superpowers/specs/2026-09-05-embodied-human.md

- [ ] Workspace: app/src/style.css/main.js and tests. Fixed viewport, toggles and independent scrolling, responsive controls, no document scroll, existing evidence/scenario workflows accessible.
- [ ] Respiration: assembly respiratory mechanics + builder/tests/docs. Native volume/pressure input, chest/diaphragm/rib displacement, explicit reduced mechanics/strain/energy assumptions. Report render contract and interaction scope.
- [ ] Peripheral system: assembly peripheral module + builder/tests/docs and targeted brain adapter input. Canonical skin/muscle/spinal/brain bindings, receptor stimuli and motor outputs with conduction delays. Preserve IBM sources, document completeness and reject unsupported anatomical certainty.
- [ ] Fine anatomy: assembly detail synthesis + builder/tests/docs. Source-constrained skin hair and representative microvascular topology with physical flow parameters; separate underlying detail resolution from display.
- [ ] Root integration: canonical body/manifest/API, scenarios/stimulus controls with shared clock, source/executor provenance, replay artifacts and coverage audit.
- [ ] Verify scientific invariants and frontend behavior, independent review, complete required checks, preserve research outputs and leave local app running.
