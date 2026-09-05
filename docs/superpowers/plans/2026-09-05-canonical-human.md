# Canonical Human Implementation Plan

> For agentic workers: use subagent-driven-development and independent domain reviews. User authorizes continuous implementation and source-constrained synthesis.

**Goal:** One executable generic human with shared coordinates, connected anatomy, mechanics, physiology and IBM-derived brain, exposed as the application's primary body.

**Architecture:** Canonical anatomical entities and connections form the shared substrate. Acquired sources constrain this substrate; source geometry, registered additions and synthesized representations retain distinct evidence labels. Native and reduced physical executors bind through explicit body-state ports.

1. Anatomy assembly: implement `ihm/assembly/anatomy.py`, canonical builder and verification. BP3D scaffold; fit additional lymphatic/soft anatomy using recorded correspondences; avoid duplicate competing organs; output anatomy registry, assumption ledger and single-family display fragment.
2. Mechanics: implement independently tested rigid/soft tissue and muscle/tendon/ligament mechanisms, derive source-informed constitutive/attachment parameters, and output canonical bindings plus executable mechanical state. Coordinate anatomy schema before integration.
3. Brain: inspect IBM-1 actual assets/runtime, reuse a source-backed brain representation and executable neural model; register to canonical brain support and record input/output coupling and inherited limitations.
4. Root integration: assemble shared canonical state, connect native physiological drivers and domain runners with an explicit clock, expose `/api/body` and materialization, make IHM body default and move source selector to evidence inspection.
5. Validate source/assumption lineage, attachment and registration consistency, constitutive behavior and coupled updates; run API/browser/body scenario tests, review findings, update documentation and leave app running.

Ownership: anatomy, mechanics and brain workers own separate assembly modules/scripts. Root owns shared API/CLI/app/build orchestration. Work stays in the existing feature branch and retains all raw assets.
