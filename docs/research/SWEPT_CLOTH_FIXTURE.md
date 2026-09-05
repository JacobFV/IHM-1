# Transactional elastic cloth first impact

The opt-in `ihm/assembly/swept_cloth.py` module advances existing finite-mass elastic body owners on one clock. `step_first_impact(bodies, dt_s, contacts=[SweptEdgePair(...)])` computes an elastic force kick, searches the specified edge pairs, drifts every owner to the common event time, applies the swept finite-mass Coulomb impulse, and drifts the remainder. It validates all proposed states before changing any owner. Invalid final elasticity, unsupported initial/coplanar contact, and multiple candidate impacts reject the entire step.

The independent fixture uses two actual triangular `Cloth` patches with area-derived nodal masses and central elastic edge springs. An isolated pair crosses at 99 µs during a 100 µs interval. Contact changes both patches' velocities and subsequently stores spring energy; disabling friction changes tangential motion. The accepted run in `data/derived/swept-cloth-yz772thn` records a 77.99999847 nJ impact loss and 20.34397 fJ final spring energy. Every trajectory includes actual calculated positions, velocities and energies. Exact source files, source hashes and copied-input hashes accompany the fresh report.

100, 50 and 25 µs integration steps agree on this first impact and final positions. The incident motion is rigid translation and the impact occurs in the final substep, so this comparison checks event scheduling parity, **not convergence of subsequent elastic relaxation**. The tiny spring energy reflects only the 1 µs post-impact remainder. The report separates the integration energy defect from contact dissipation and records linear/angular momentum residuals. Kernel tests independently check fixed and moving support work.

Only the specified pair is scheduled; the fixture does not claim collision coverage of every edge or triangle. There is no sustained-contact state, post-impact event rescan, general multi-event schedule, or self-contact. The ballistic remainder is explicitly reported. This path remains separate from the default body viewer and existing face-interior panel artifacts. No garment containment result follows from it.

All fixture density, spring stiffness and friction values are explicit engineering priors. Next whole-garment work targets source-pinned shirt/shorts connectivity, seams and support over the torso, shoulders, hips and thighs, with physical contact coverage still a required gate.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_swept_edge_contact.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_swept_cloth.py
```

The second command creates a fresh immutable run directory. Its negative tests demonstrate that simultaneous unscheduled impacts and invalid final elasticity leave every body's arrays and clock unchanged.
