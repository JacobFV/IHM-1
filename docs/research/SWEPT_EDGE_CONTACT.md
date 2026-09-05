# Opt-in generic swept edge impact

`ihm/assembly/swept_edge_contact.py` supplies a pure generic collision-event and finite-mass impact kernel. It has no anatomical specialization and does not change the existing face-interior garment panel solver or its retained trajectories. Subsequent garment integration targets the torso, shoulders, hips and thighs as whole-body contact interfaces.

`swept_edge_event(xa, va, xb, vb, duration_s)` searches two linear endpoint trajectories. Coplanarity is a cubic polynomial in normalized interval time. Candidate roots must place both segment parameters within `[0,1]` and produce coincident contact points. The returned record is `impact`, `no_impact`, or `unresolved`; an impact includes its time, endpoint barycentric weights, common point, normal and numerical point-separation residual.

`resolve_edge_impact(xa_at_event, va, ma, xb_at_event, vb, mb, event, friction_static=..., friction_kinetic=..., mobile_a=..., mobile_b=...)` is a pure zero-restitution operator. It returns both endpoint velocity arrays and impulses, normal/tangential impulse, static/sliding regime, dissipation, support impulse/work, momentum residual, angular impulse residual and energy residual. Fixed endpoints retain their prescribed velocities, including moving supports; the corresponding support work is explicit. All lengths, times, masses and impulses use SI units.

The impulse occurs at a shared zero-thickness contact point. There is no artificial attraction inside a proximity band and no unmodeled finite-clearance friction couple. The implementation uses floating-point polynomial roots and checked coincidence, not a certified exact collision predicate. Initial/persistent contact, coplanar/parallel motion and grazing degeneracies are reported unresolved. This kernel alone supplies neither a multi-event schedule, self-contact, contact thickness nor garment containment.

The analytic midpoint fixture has four 2 g endpoints, 0.1 m/s relative normal approach, a 10 ms event time, 0.0002 N·s impulse and 10 µJ kinetic loss. Tests additionally cover endpoint hits, both crossing directions and endpoint orders, high-speed tunneling missed by endpoint sampling, nonuniform endpoint velocities, translated geometry, event-time agreement across subdivided search windows, frictional sticking/sliding, fixed and moving supports, no-impact cases, invalid input and explicit geometric degeneracy. Input arrays remain unchanged.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/verify_swept_edge_contact.py
```

The numerical receipt is retained under `artifacts/verification/swept-edge-contact/kernel-report.json`, with implementation/verifier hashes. The next slice is transactional common-clock cloth-patch impact and source-pinned full-garment seam/support work; no display animation is substituted for that mechanical calculation.
