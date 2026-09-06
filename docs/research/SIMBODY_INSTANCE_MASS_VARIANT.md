# Opt-in Simbody instance-mass source variant

This is a separately built correction to the concrete native limitation measured
in `LOCAL_INTAKE_MASS_INSTANCE.md`. It does not replace the installed library or
automatically enable intake mass changes in embodied runtime.

`patch_simbody_instance_mass.py` requires exact SHA-256 identities for seven
source files from retained Simbody revision
`59c6e7b89b3bdf266a2f3d54c611599d964205f7`. It patches a copied source directory:

- Public mass properties now read the State's existing instance variables.
- Generic position kinematics receives those mass properties; spatial inertia
  and COM derive from them. Gyroscopic force uses the resulting spatial inertia.
- Weld and LoneParticle specialized paths use instance mass/inertia too.
- An additive exported C ABI, declared by `native_simbody_instance_mass.h`,
  validates and writes one body instance and invalidates Instance/higher caches.

No object layout changes, const-casts, topology mutation, model initialization
or per-node mutable mass cache are introduced. The first exported setter requires
nondegenerate physical rigid-body inertia for both original and proposed body.
It rejects zero/negative/nonfinite mass and invalid central inertia. It does not
change speed, apply a mass impulse or invent transport velocity; those belong
to the caller's explicit material flux transaction.

The setter changes Simbody State properties. OpenSim Model::getTotalMass and
Model::calcMassCenterPosition delegate to Simbody and therefore observe them.
Static OpenSim Body model/XML properties remain the baseline; continuing mass
must be retained with the State and an explicit material-port receipt.

## Resumable build and resource scope

Prepare only (no compiler):

```sh
.venv/bin/python scripts/build_simbody_instance_mass_variant.py \
  --output data/runtime/opensim/variants/instance_mass_v1
```

The dependency inventory identifies 19 replacement translation units from the
retained shared-library build. Other donor objects are copied and hashed. Each
native invocation requires the coordinated heavy slot and compiles at most three
objects using O1, single threads, nice 10 and a 4 GiB address-space ceiling.
The Python caller should run with a 1 GiB soft / 4 GiB hard address-space limit.
No intermediate library is linked or published. Resume a chunk with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
prlimit --as=1073741824:4294967296 -- nice -n 10 \
.venv/bin/python scripts/build_simbody_instance_mass_variant.py \
  --output data/runtime/opensim/variants/instance_mass_v1 \
  --build-native --max-objects 3
```

Use `--only-object MobilizedBody.cpp.o --max-objects 1` for one pending object.
`--link` requires all replacement objects to exist; it never changes the installed
library, symlinks or a latest pointer. Source/header/library/donor hashes and
completed object hashes are checked on every resume and before linking.
Per-object command, output, wall time and maximum child RSS remain beside the
manifest. A failure retains the object log without marking that object complete.

`verify_simbody_instance_mass_recipe.py` is a three-test source-only check of
exact source patching, mismatch rejection and frozen donor/overlay guards.
`verify_simbody_instance_mass_variant.py --variant ... --run-native` separately
compiles/runs the positive native fixture only after the library is linked,
checks loaded-library selection, and retains its result. That fixture checks
mass/COM/M/gravity/kinetic consistency, localized capture and passive outflow,
momentum and energy flux, concurrent States, restore/replay and invalid mass.
It is deliberately limited to independent free bodies; full constrained OpenSim
muscle and contact acceptance is still required before a whole-body material port.

## Material ownership remains a separate integration task

A corrected mass matrix does not locate ingested contents or close physiology
mass exits. Native stomach/chyme/vascular/bladder transfers need nonoverlapping
baseline inventories and registered local mass distributions. Ingested mass
must be consumed once, internal transfer must conserve system mass, and urine,
sweat/gas exits need actual boundary receipts. Never uniformly scale segments,
double-count GI contents, or erase performed mechanical work to make a transfer
appear balanced. Capture dissipation and carried enthalpy need explicitly named
thermal owners distinct from the native active-muscle energy ledger.

## Retained positive acceptance, 2026-09-05

All 19 replacement translation units compiled in serialized chunks (38.511 s
summed compiler wall time, 522,384 KiB maximum compiler RSS). The 59-object
inventory includes 40 unchanged donor objects. The opt-in library at
`data/runtime/opensim/variants/instance_mass_v1/libSimTKsimbody.so.3.9` has SHA-256
`a0cc25c0543aba0a16ca80e7a2b111db78f7ed0a0ec83439266ffd6e90625da7`.

`data/derived/instance-mass-positive-9hi1mldg/report.json` passed the actual native
positive fixture. Unlike the original library's zero change, the new mass matrix
changed by norm 0.9485778829 and the co-moving payload added exactly 0.2938 J
kinetic energy. Arbitrary-velocity capture recorded 0.3194725347 J dissipation
with momentum residual norm 1.47e-15; passive outflow had dissipation 8.88e-16 J
and momentum residual norm 1.01e-15. Local mass, COM, gravitational potential,
COM acceleration under gravity, concurrent independent States, restore/replay
and invalid-mass rejection passed. Fixture compile/run took 5.590 s with
455,200 KiB maximum child RSS. The receipt binds fixture source, selected variant
and every library reported by the loader.

This establishes the isolated corrected engine behavior. No default library,
mechanical stream or physiology mass coupling was changed. The next acceptance
must retain actual OpenSim muscle activation/fiber states and signed energy
ledgers while updating an explicitly registered local GI payload, then verify
constraint/contact impulse and per-owner material incidence.
