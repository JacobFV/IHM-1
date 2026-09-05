# Muscle bulk intersection screening

The two source-preserving dimensional partitions expose424 candidate bulk
surfaces. Closed edge incidence alone is insufficient for volume construction:
the next native geometric check finds self-intersections in180 of them.
**244 pass both intersection and vertex-link screening.** One of425 original
muscles remains blocked by the earlier degeneracy gate and is not submitted to
the native intersection routine. None of these results changes source geometry.

The retained run is `data/derived/muscle-intersections-v1/`. It takes18.99 seconds
on the local machine. Every face-pair result maps through the retained bulk
partition to original source-face indices. Exact numeric coordinate welding
joins authored seams without moving coordinates; inputs, mappings, native
module bytes and the executing script are copied and hashed. The separate
verifier checks every recorded coordinate and face mapping without rerunning
the native computation.

The backend is the pinned libigl2.6.2 CGAL `remesh_self_intersections` binding,
called with `detect_only=True` and `first_only=False`. Despite the routine's
name, this audit performs no remeshing. Its
[upstream API](https://libigl.github.io/libigl-python-bindings/api/igl_copyleft_cgal/)
describes the detection mode. Native fixtures test a proper tetrahedron,
duplicated raw seam vertices, crossing triangles, coplanar overlap, pinched
vertex links, degenerate rejection and nested disjoint shells. Nesting passes
intersection detection by design; an empty intersection list does not resolve
cavity conventions, overlap of enclosed volumes, or physical tissue occupancy.

No retained candidate is yet certified as a whole-body finite-element domain.
Cross-entity collisions, source registration, shell orientation/nesting, thin
structures, material assignment and exclusive mass ownership remain distinct
questions. In particular, deleting intersecting faces would discard evidence;
conditional reconstruction needs a specified volume/field interpretation and
an uncertainty record.

The ARM host uses an isolated Python3.10.20 geometry-tool environment because
the stable2.6.2 ARM wheel is built for CPython3.10. The main IHM environment is
unchanged. Retained PyPI metadata and wheel receipts are under
`data/raw/geometry/libigl/2.6.2/`; versions are libigl2.6.2, NumPy2.2.6 and
SciPy1.15.3. An initially acquired x86 wheel was rejected as incompatible by the
installer; those bytes remain separate from the actually used ARM wheel.

```sh
# Small native fixtures, single thread and lower scheduling priority:
nice -n 10 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  data/runtime/geometry/libigl-2.6.2/venv/bin/python \
  scripts/audit_muscle_intersections.py --self-test
# Retained artifact/mapping check; does not launch a browser or native solver:
nice -n 10 env OPENBLAS_NUM_THREADS=1 .venv/bin/python \
  scripts/verify_muscle_intersections.py data/derived/muscle-intersections-v1
```

A fresh full audit uses `--output` with a new directory. Original partitions,
failed attempts, source meshes and all canonical material ownership remain
unchanged.
