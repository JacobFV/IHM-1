# Cellular voltage and measured skin fields

The workbench's **BETSE · cellular bioelectricity** model shows 212 actual solver
cells and 1,221 membrane edges from the retained native BETSE simulation. The
source is planar; the app embeds those real polygons in a 3D viewport without
inventing tissue depth or registering them to the reference human's skin.

Its 34 saved states contain native cell-averaged membrane voltage, sodium,
potassium, protein and anion concentrations, and a mean incident-membrane gap
junction open fraction. Picking a cell exposes its values at the selected time.
Each field keeps one color range across the entire recording. Constant fields
remain constant; color normalization does not invent variation.

Membrane voltage is in **V**. The separately fitted human wound measurement is
an extracellular lateral field in **V/m**. The observation operator is the
extracellular potential gradient, not cellular membrane voltage. These are
separate materializations with separate source evidence.

```bash
.venv-physiology/bin/python scripts/export_betse_tissue.py
.venv/bin/python scripts/verify_bioelectric.py
```

The exporter will unpickle only the SHA-256-pinned artifact generated locally
by the reviewed BETSE solver. It writes portable JSON; neither the browser nor
the HTTP API unpickles data. The exact generated random geometry and configuration
are retained, although the initial random seed was not fixed. A fresh random
simulation requires a new reviewed artifact receipt; it cannot impersonate the
saved source run.

The 0.035-second generic tissue simulation checks execution and exposes cellular
heterogeneity. It is not a human wound-healing replication, does not identify
human skin channel/pump kinetics, and supplies no healing-time prediction.
