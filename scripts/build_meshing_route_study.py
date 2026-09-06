"""Two routes past the ~50-entity ceiling of the cross-structure conforming mesh, measured.

The existing pipeline (scripts/build_cross_structure_conflict_repair.py) imprints every crossing
curve with an exact CGAL arrangement, then meshes the padded box with TetGen `pY`. It succeeds to
48 entities and fails at 37, 45, 49, 50, 55 and 60. The receipts in
data/derived/cross-structure-repair-v1/summary.json say why, and the reason is NOT the one the
edge-length note in that module guesses:

  succeeded rung 48: plc_residual_self_intersecting_face_pairs = 0, edges under 1e-8 m = 4
  failed    rung 50: plc_residual_self_intersecting_face_pairs = 2, edges under 1e-8 m = 0

TetGen tolerates the sub-1e-8 edges (it Steiner-splits the offending segments, warns
"4 segments are not recovered", and still returns status 0). What it cannot tolerate is a PLC that
still self-intersects, which is what the arrangement leaves behind after CGAL's exact intersection
points are rounded to double. So the discriminator is `plc_residual_self_intersecting_face_pairs`.

ROUTE A - snap rounding. Quantise every arrangement vertex onto a fixed cubic lattice, weld on the
integer lattice key, drop what degenerates, and re-run the exact arrangement to certify or repair
what the quantisation broke. Iterate to a fixpoint. The lattice gives a hard lower bound on every
edge (two distinct lattice points are at least `grid` apart) and, if the fixpoint has zero residual
crossings, a PLC TetGen can take.

ROUTE B - a tolerant mesher. fTetWild, built for this host under
data/runtime/tolerant-mesher, meshed on the same clusters with filtering disabled so the whole box
is kept and the same winding-number ownership step can run on it.

Both routes are then pushed through the SAME downstream measurement: winding-number labelling,
priority ownership, per-structure volume error, extracted owned surfaces, exact CGAL pairwise
residual overlap. Route A does this by substituting its arrangement for the repair module's and
calling that module's `mesh_cluster` unchanged; Route B does it on the fTetWild tet mesh.

The decisive property is CONFORMING SHARED NODES at tissue interfaces. Both routes are checked for
it explicitly and by construction-independent measurement.

Run with the isolated libigl environment:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_meshing_route_study.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_meshing_route_study.py \
      --out data/derived/meshing-route-study-v1 --route both --clusters thigh-six,ladder-50

Nothing under scripts/build_cross_structure_conflict_repair.py or data/derived/cross-structure-
repair-v1 is modified.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, os, shutil, subprocess, sys, tempfile, time
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_cross_structure_conflict_repair as R  # noqa: E402  (imported, never modified)
import igl                                          # noqa: E402
import igl.copyleft.cgal as cgal                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REPAIR = ROOT / 'data/derived/cross-structure-repair-v1'
FTETWILD = ROOT / 'data/runtime/tolerant-mesher/build/FloatTetwild_bin'

# ---------------------------------------------------------------- tolerances
#
# Every tolerance below is justified against two measured anatomical scales taken from the existing
# receipts and from docs/research/SOFT_BODY_MATERIALIZATION.md:
#
#   0.489 mm   median 2V/A thickness proxy of a cross-structure overlap lens (resolution.json)
#   1-2 mm     articular cartilage thickness, the finest structure the atlas intends to carry
#   ~1.6-2.0 mm median PLC edge length of the arrangement at rungs 12-60
#
# SURFACE_DEVIATION_BUDGET_M is set at 1 percent of the 0.489 mm overlap thickness. A displacement
# of the tissue interface by 1 percent of the thinnest thing the ownership step has to arbitrate
# cannot change which side of the interface a tet falls on for any feature the atlas resolves, and
# it is 0.25-0.5 percent of the cartilage scale.
SURFACE_DEVIATION_BUDGET_M = 4.89e-6
# VOLUME_ERROR_BUDGET is per structure, |labelled volume / surface divergence volume - 1|. The
# exact-arrangement route achieves 1.83e-07. Snap rounding perturbs the boundary by up to
# grid*sqrt(3)/2, so the achievable bound is set by deviation*area/volume, of order 1e-4 for a
# thigh muscle at a 1 um grid. 1e-3 is the acceptance bar; the achieved number is always reported.
VOLUME_ERROR_BUDGET = 1e-3
# TetGen's own segment tolerance in these units, quoted by TetGen itself in the rung-50 receipt:
# "default limit is 1.08111e-08". A lattice at or above 1e-7 m clears it by an order of magnitude.
TETGEN_SEGMENT_TOLERANCE_M = 1.08111e-8
GRIDS_DEFAULT = (1e-5, 1e-6, 1e-7)


def sha(path):
    return R.sha(path)


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def boundary_and_nonmanifold(W, F):
    if not len(F):
        return 0, 0
    e = np.sort(np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]), axis=1)
    _, c = np.unique(e, axis=0, return_counts=True)
    return int((c == 1).sum()), int((c > 2).sum())


def facet_degeneracy(W, F):
    """Sub-tolerance features that are NOT short edges.

    TetGen's crash logs at the failing rungs quote segment-to-vertex distances of 6.9e-18 m, not
    short edges. A triangle whose smallest altitude is 6.9e-18 m IS a vertex sitting 6.9e-18 m from
    the opposite segment, so the smallest facet altitude is the quantity those warnings are about
    and the one the receipts never measured."""
    if not len(F):
        return {}
    L = np.stack([np.linalg.norm(W[F[:, (i + 1) % 3]] - W[F[:, i]], axis=1) for i in range(3)], 1)
    a = np.linalg.norm(np.cross(W[F[:, 1]] - W[F[:, 0]], W[F[:, 2]] - W[F[:, 0]]), axis=1) / 2.0
    alt = 2.0 * a / np.maximum(L.max(1), 1e-300)
    return {'min_facet_area_m2': float(a.min()),
            'min_facet_altitude_m': float(alt.min()),
            'facets_with_altitude_under_tetgen_tolerance': int((alt < TETGEN_SEGMENT_TOLERANCE_M).sum()),
            'facets_with_altitude_under_1e-14_m': int((alt < 1e-14).sum()),
            'max_facet_aspect_ratio': float((L.max(1) / np.maximum(alt, 1e-300)).max()),
            'note': 'smallest facet altitude is the vertex-to-opposite-segment distance TetGen '
                    'reports in its "a segment and a vertex are very close" warnings'}


def detect_crossings(W, F):
    """Exact CGAL self-intersection detection only; no remeshing."""
    return int(len(np.asarray(cgal.remesh_self_intersections(
        np.ascontiguousarray(W), np.ascontiguousarray(F), True, False, False, False, 20_000_000
    )[2]).reshape(-1, 2)))


# --------------------------------------------------------------------- route A

def lattice_quantise(W, F, src, grid):
    """Move every vertex to the nearest point of a cubic lattice of spacing `grid`, weld on the
    integer lattice key, then reuse the repair module's cleaner for degenerate and duplicate facets.

    Two distinct lattice points differ by at least `grid` in one coordinate, so after this the
    minimum nonzero edge length of the complex is at least `grid` by construction, not by luck."""
    W = np.asarray(W, float)
    reach = float(np.abs(W).max()) if len(W) else 1.0
    # the lattice is only sound if consecutive lattice points are separated by many double ulps at
    # the largest coordinate in play, and if the lattice index itself is exactly representable
    if grid <= 16 * np.spacing(max(reach, 1.0)) or reach / grid >= 2.0 ** 52:
        raise ValueError('grid %g is below double resolution at coordinate magnitude %g '
                         '(ulp %g, index %g)' % (grid, reach, np.spacing(max(reach, 1.0)), reach / grid))
    K = np.rint(W / grid)
    Wq = np.ascontiguousarray(K * grid)
    distinct_keys = len(np.unique(K, axis=0))
    distinct_doubles = len(np.unique(Wq, axis=0))
    if distinct_keys != distinct_doubles:
        raise ValueError('grid %g is below double resolution at this coordinate magnitude '
                         '(%d lattice keys collapsed to %d doubles)'
                         % (grid, distinct_keys, distinct_doubles))
    moved = float(np.abs(Wq - np.asarray(W, float)).max())
    W2, F2, src2, repeated, zero, collapsed = R.clean(Wq, F, src)
    return W2, F2, src2, {'vertices_before': int(len(W)), 'vertices_after': int(len(W2)),
                          'facets_before': int(len(F)), 'facets_after': int(len(F2)),
                          'max_coordinate_move_m': moved,
                          'degenerate_facets_dropped': repeated + zero,
                          'duplicate_facets_collapsed': collapsed}


def soup(parts):
    verts, faces, src, offset = [], [], [], 0
    for i, p in enumerate(parts):
        verts.append(p['V']); faces.append(p['F'] + offset)
        src.append(np.full(len(p['F']), i)); offset += len(p['V'])
    return (np.ascontiguousarray(np.concatenate(verts)),
            np.ascontiguousarray(np.concatenate(faces).astype(np.int64)),
            np.concatenate(src))


def exact_arrangement(V, F, src):
    r = cgal.remesh_self_intersections(V, F, False, False, True, False, 20_000_000)
    AV = np.ascontiguousarray(np.asarray(r[0], float))
    AF = np.ascontiguousarray(np.asarray(r[1], np.int64))
    crossings = int(len(np.asarray(r[2]).reshape(-1, 2)))
    J = np.asarray(r[3]).ravel()
    W, DF, dsrc, repeated, zero, collapsed = R.clean(AV, AF, src[J])
    return W, DF, dsrc, crossings, repeated + zero, collapsed


GRID = None
MAX_SNAP_ROUNDS = 6


def arrange_snap_rounded(parts, snap_below=0.0, tet_tolerance=1e-8, rounds=4):
    """Drop-in replacement for `build_cross_structure_conflict_repair.arrange`.

    Exact arrangement, then quantise/re-arrange to a fixpoint on the lattice `GRID`."""
    grid = GRID if (GRID or 0) > 0 else None
    V, F, src = soup(parts)
    began = time.monotonic()
    W, DF, dsrc, crossings, dropped, collapsed = exact_arrangement(V, F, src)
    arrangement_seconds = time.monotonic() - began
    base_boundary, _ = boundary_and_nonmanifold(W, DF)
    passes = []
    converged = False
    if grid is None:
        residual = detect_crossings(W, DF)
        converged = residual == 0
    else:
        for it in range(MAX_SNAP_ROUNDS):
            t0 = time.monotonic()
            W, DF, dsrc, q = lattice_quantise(W, DF, dsrc, grid)
            residual = detect_crossings(W, DF)
            be, nm = boundary_and_nonmanifold(W, DF)
            _, L = R.edge_lengths(W, DF)
            rec = {'pass': it, 'seconds': time.monotonic() - t0, **q,
                   'residual_self_intersecting_face_pairs': residual,
                   'boundary_edges': be, 'nonmanifold_edges': nm,
                   'min_edge_length_m': float(L.min()) if len(L) else None,
                   'edges_under_tetgen_tolerance': int((L < TETGEN_SEGMENT_TOLERANCE_M).sum()),
                   'degeneracy': facet_degeneracy(W, DF)}
            passes.append(rec)
            if residual == 0 and be == 0:
                converged = True
                break
            if it + 1 == MAX_SNAP_ROUNDS:
                break
            t1 = time.monotonic()
            W, DF, dsrc, c2, d2, k2 = exact_arrangement(W, DF, dsrc)
            rec['re_arrangement_seconds'] = time.monotonic() - t1
            rec['re_arrangement_crossings_resolved'] = c2
        residual = passes[-1]['residual_self_intersecting_face_pairs']
    ue, ec = np.unique(np.sort(np.concatenate([DF[:, [0, 1]], DF[:, [1, 2]], DF[:, [2, 0]]]),
                               axis=1), axis=0, return_counts=True)
    L = np.linalg.norm(W[ue[:, 0]] - W[ue[:, 1]], axis=1)
    # two-sided deviation between the final PLC and the untouched input soup
    dev = deviation(V, F, W, DF)
    report = {'route': 'snap-round' if grid else 'exact-arrangement',
              'grid_m': grid,
              'input_vertices': int(len(V)), 'input_facets': int(len(F)),
              'crossing_face_pairs_resolved': crossings,
              'arrangement_seconds': arrangement_seconds,
              'plc_vertices': int(len(W)), 'plc_facets': int(len(DF)),
              'plc_boundary_edges': int((ec == 1).sum()),
              'plc_boundary_edges_before_snapping': base_boundary,
              'plc_nonmanifold_edges': int((ec > 2).sum()),
              'plc_residual_self_intersecting_face_pairs': residual,
              'plc_min_edge_length_m': float(L.min()),
              'plc_median_edge_length_m': float(np.median(L)),
              'plc_edges_under_tetgen_tolerance': int((L < TETGEN_SEGMENT_TOLERANCE_M).sum()),
              'plc_degeneracy': facet_degeneracy(W, DF),
              'snap_passes': passes, 'snap_converged': converged,
              'deviation_from_input': dev,
              'deviation_within_budget': bool(dev['hausdorff_two_sided_m']
                                              <= SURFACE_DEVIATION_BUDGET_M),
              'surface_deviation_budget_m': SURFACE_DEVIATION_BUDGET_M,
              'plc_valid': bool((ec == 1).sum() == 0 and residual == 0)}
    return W, DF, dsrc, report


def deviation(VA, FA, VB, FB):
    """Two-sided vertex-to-surface deviation, in metres. One-sided each way over vertices only,
    which is a lower bound on the true Hausdorff distance and the quantity the lattice bounds."""
    out = {}
    if len(FA) and len(VB):
        d1 = np.sqrt(np.asarray(igl.point_mesh_squared_distance(
            np.ascontiguousarray(VB), np.ascontiguousarray(VA), np.ascontiguousarray(FA))[0], float))
        out['snapped_to_input_max_m'] = float(d1.max())
        out['snapped_to_input_mean_m'] = float(d1.mean())
        out['snapped_to_input_p99_m'] = float(np.percentile(d1, 99))
    if len(FB) and len(VA):
        d2 = np.sqrt(np.asarray(igl.point_mesh_squared_distance(
            np.ascontiguousarray(VA), np.ascontiguousarray(VB), np.ascontiguousarray(FB))[0], float))
        out['input_to_snapped_max_m'] = float(d2.max())
        out['input_to_snapped_mean_m'] = float(d2.mean())
        out['input_to_snapped_p99_m'] = float(np.percentile(d2, 99))
    out['hausdorff_two_sided_m'] = max(out.get('snapped_to_input_max_m', 0.0),
                                       out.get('input_to_snapped_max_m', 0.0))
    return out


def run_route_a(parts, grid, pad, flags, priority, label, emit=None):
    """Substitute the snap-rounding arrangement and let the repair module do everything else."""
    global GRID
    saved = R.arrange
    GRID = grid
    R.arrange = arrange_snap_rounded
    try:
        began = time.monotonic()
        report, mesh = R.mesh_cluster(parts, pad, flags, priority, emit=emit, label=label,
                                      snap_below=0.0)
        report['wall_seconds'] = time.monotonic() - began
    finally:
        R.arrange = saved
        GRID = None
    report['route'] = 'A-snap-round' if (grid or 0) > 0 else '0-exact-arrangement-control'
    report['grid_m'] = grid
    return report, mesh


# --------------------------------------------------------------------- route B

def write_obj(path, V, F):
    with open(path, 'w') as f:
        f.write('\n'.join('v %.17g %.17g %.17g' % tuple(p) for p in V))
        f.write('\n')
        f.write('\n'.join('f %d %d %d' % (a + 1, b + 1, c + 1) for a, b, c in F))
        f.write('\n')


def read_msh_tets(path):
    """Minimal gmsh 2.2 / 4.1 ASCII reader for the node block and tetra elements fTetWild writes."""
    text = Path(path).read_text().splitlines()
    i = 0
    version = None
    V = None
    T = []
    while i < len(text):
        line = text[i].strip()
        if line == '$MeshFormat':
            version = text[i + 1].split()[0]
            i += 1
        elif line == '$Nodes':
            if version and version.startswith('4'):
                hdr = [int(x) for x in text[i + 1].split()]
                nblocks, nnodes = hdr[0], hdr[1]
                V = np.zeros((nnodes, 3))
                j = i + 2
                for _ in range(nblocks):
                    b = [int(x) for x in text[j].split()]
                    n = b[3]
                    tags = [int(text[j + 1 + k]) for k in range(n)]
                    for k in range(n):
                        V[tags[k] - 1] = [float(x) for x in text[j + 1 + n + k].split()]
                    j += 1 + 2 * n
                i = j - 1
            else:
                n = int(text[i + 1])
                V = np.zeros((n, 3))
                for k in range(n):
                    parts = text[i + 2 + k].split()
                    V[int(parts[0]) - 1] = [float(x) for x in parts[1:4]]
                i += 1 + n
        elif line == '$Elements':
            if version and version.startswith('4'):
                hdr = [int(x) for x in text[i + 1].split()]
                nblocks = hdr[0]
                j = i + 2
                for _ in range(nblocks):
                    b = [int(x) for x in text[j].split()]
                    etype, n = b[2], b[3]
                    for k in range(n):
                        vals = [int(x) for x in text[j + 1 + k].split()]
                        if etype == 4:
                            T.append([v - 1 for v in vals[1:5]])
                    j += 1 + n
                i = j - 1
            else:
                n = int(text[i + 1])
                for k in range(n):
                    vals = [int(x) for x in text[i + 2 + k].split()]
                    if vals[1] == 4:
                        ntags = vals[2]
                        T.append([v - 1 for v in vals[3 + ntags:3 + ntags + 4]])
                i += 1 + n
        i += 1
    return np.ascontiguousarray(V, float), np.ascontiguousarray(np.asarray(T, np.int64))


def run_ftetwild(parts, workdir, epsr, lr, extra=(), timeout=3600):
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    V, F, _ = soup(parts)
    obj = workdir / 'input.obj'
    write_obj(obj, V, F)
    out = workdir / 'out.msh'
    cmd = [str(FTETWILD), '-i', str(obj), '-o', str(out), '--no-binary',
           '--disable-filtering', '-e', str(epsr), '-l', str(lr),
           '--max-threads', str(min(12, os.cpu_count() or 4)), *extra]
    began = time.monotonic()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    seconds = time.monotonic() - began
    produced = sorted(p.name for p in workdir.iterdir())
    rec = {'command': ' '.join(cmd), 'returncode': proc.returncode, 'seconds': seconds,
           'stdout_tail': proc.stdout.strip().splitlines()[-25:],
           'stderr_tail': proc.stderr.strip().splitlines()[-15:],
           'files_produced': produced,
           'input_vertices': int(len(V)), 'input_facets': int(len(F))}
    candidates = [p for p in workdir.iterdir() if p.suffix == '.msh']
    if proc.returncode != 0 or not candidates:
        rec['succeeded'] = False
        return rec, None, None
    msh = max(candidates, key=lambda p: p.stat().st_size)
    rec['mesh_file'] = msh.name
    TV, TT = read_msh_tets(msh)
    rec['succeeded'] = bool(len(TT))
    rec['tet_vertices'] = int(len(TV))
    rec['tets'] = int(len(TT))
    return rec, TV, TT


TET_FACES = R.TET_FACES


def label_and_measure(TV, TT, parts, priority_of, box_volume=None):
    """The repair module's ownership step, applied to any tet mesh of the padded box."""
    det = np.linalg.det(np.swapaxes(TV[TT[:, 1:]] - TV[TT[:, 0, None]], 1, 2)) / 6.0
    negative = int((det < 0).sum())
    if negative:
        TT = TT.copy(); TT[det < 0] = TT[det < 0][:, [0, 2, 1, 3]]
        det = np.linalg.det(np.swapaxes(TV[TT[:, 1:]] - TV[TT[:, 0, None]], 1, 2)) / 6.0
    centroid = np.ascontiguousarray(TV[TT].mean(1))
    inside = np.zeros((len(parts), len(TT)), bool)
    for i, p in enumerate(parts):
        w = np.asarray(igl.fast_winding_number(np.ascontiguousarray(p['V']),
                                               np.ascontiguousarray(p['F']), centroid), float)
        inside[i] = np.abs(w) > .5
    claims = inside.sum(0)
    seq = sorted(range(len(parts)),
                 key=lambda i: (priority_of[parts[i]['entity_id']],
                                -parts[i]['surface_volume_m3'], parts[i]['entity_id']))
    owner = np.full(len(TT), -1, np.int64)
    for i in reversed(seq):
        owner[inside[i]] = i
    per = []
    for i, p in enumerate(parts):
        claimed = float(np.abs(det[inside[i]]).sum()); own = float(np.abs(det[owner == i]).sum())
        per.append({'entity_id': p['entity_id'], 'name': p['name'], 'role': p['role'],
                    'surface_volume_m3': p['surface_volume_m3'],
                    'claimed_tets': int(inside[i].sum()), 'claimed_volume_m3': claimed,
                    'claimed_relative_volume_error':
                        float(claimed / p['surface_volume_m3'] - 1) if p['surface_volume_m3'] else None,
                    'owned_tets': int((owner == i).sum()), 'owned_volume_m3': own})
    q = tet_quality(TV, TT)
    mesh = {'tet_vertices': int(len(TV)), 'tets': int(len(TT)),
            'tets_negatively_oriented_by_mesher': negative,
            'all_positive_volume': bool((det > 0).all()),
            'min_tet_volume_m3': float(det.min()), 'zero_volume_tets': int((det <= 0).sum()),
            'total_volume_m3': float(np.abs(det).sum()),
            'complement_tets': int((owner < 0).sum()),
            'structure_tets': int((owner >= 0).sum()),
            'structure_volume_m3': float(np.abs(det[owner >= 0]).sum()),
            'disputed_tets': int((claims > 1).sum()),
            'max_claimants_on_one_tet': int(claims.max()),
            'per_structure': per,
            'max_abs_claimed_volume_error':
                float(max(abs(c['claimed_relative_volume_error'] or 0.0) for c in per)),
            'quality': q}
    if box_volume:
        mesh['box_volume_closure_relative_error'] = float(np.abs(det).sum() / box_volume - 1)
    return mesh, TT, owner, det


def tet_quality(TV, TT, sample=400_000, seed=0):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(TT)) if len(TT) <= sample else rng.choice(len(TT), sample, replace=False)
    P = TV[TT[idx]]
    e = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    L = np.stack([np.linalg.norm(P[:, b] - P[:, a], axis=1) for a, b in e], 1)
    v = np.abs(np.linalg.det(np.swapaxes(P[:, 1:] - P[:, :1], 1, 2)) / 6.0)
    # circumradius-to-shortest-edge is expensive; use the standard normalised shape measure
    # gamma = 12*(3V)^(2/3) / sum(L^2), 1.0 for the regular tetrahedron
    gamma = 12.0 * np.cbrt(3.0 * v) ** 2 / np.maximum((L ** 2).sum(1), 1e-300)
    # minimum dihedral angle
    faces = np.array([[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]])
    N = np.cross(P[:, faces[:, 1]] - P[:, faces[:, 0]], P[:, faces[:, 2]] - P[:, faces[:, 0]])
    N = N / np.maximum(np.linalg.norm(N, axis=2, keepdims=True), 1e-300)
    ang = []
    for a in range(4):
        for b in range(a + 1, 4):
            c = np.clip(-np.einsum('ij,ij->i', N[:, a], N[:, b]), -1.0, 1.0)
            ang.append(np.degrees(np.arccos(c)))
    ang = np.stack(ang, 1)
    return {'sampled_tets': int(len(idx)),
            'shape_gamma_min': float(gamma.min()), 'shape_gamma_mean': float(gamma.mean()),
            'shape_gamma_p01': float(np.percentile(gamma, 1)),
            'min_dihedral_deg': float(ang.min()),
            'min_dihedral_p01_deg': float(np.percentile(ang.min(1), 1)),
            'min_dihedral_median_deg': float(np.median(ang.min(1))),
            'max_dihedral_deg': float(ang.max()),
            'tets_with_min_dihedral_under_5deg': int((ang.min(1) < 5).sum()),
            'gamma_note': 'gamma = 12*(3V)^(2/3)/sum(edge^2); 1.0 is the regular tetrahedron'}


def shared_node_audit(TV, TT, owner, parts, sample_pairs=200):
    """Does the interface between two owned regions consist of tet faces with SHARED nodes?

    Independent of how the mesh was made: build the face table of the whole tet mesh, find every
    interior face whose two incident tets have different owners, and check that the face is one
    triangle referenced by BOTH tets through the SAME three vertex indices. If a mesher had left a
    non-conforming interface, the two sides would carry distinct coincident vertices and the face
    would appear twice as a boundary face rather than once as an interior face."""
    tri = TT[:, TET_FACES.ravel()].reshape(-1, 3)
    tet = np.repeat(np.arange(len(TT)), 4)
    key = np.sort(tri, axis=1)
    uniq, inverse, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    inverse = np.ravel(inverse)
    interior = counts[inverse] == 2
    dangling = int((counts > 2).sum())
    # owner pair per unique face
    order = np.argsort(inverse, kind='stable')
    si = inverse[order]
    bounds = np.r_[0, 1 + np.flatnonzero(si[1:] != si[:-1]), len(si)]
    pair = (bounds[1:] - bounds[:-1]) == 2
    lo = bounds[:-1][pair]
    a = order[lo]; b = order[lo + 1]
    oa = owner[tet[a]]; ob = owner[tet[b]]
    hetero = oa != ob
    both_structure = hetero & (oa >= 0) & (ob >= 0)
    # coordinate check: every interface face is one triangle, so its three vertices are literally
    # the same rows of TV for both sides. Verified by construction of the index test above; the
    # measurement that can still fail is a duplicated coincident vertex elsewhere in TV.
    duplicate_positions = int(len(TV) - len(np.unique(TV, axis=0)))
    return {'total_faces': int(len(uniq)),
            'interior_faces': int((counts == 2).sum()),
            'boundary_faces': int((counts == 1).sum()),
            'faces_with_more_than_two_incident_tets': dangling,
            'interface_faces_between_different_owners': int(hetero.sum()),
            'interface_faces_between_two_structures': int(both_structure.sum()),
            'duplicate_vertex_positions_in_mesh': duplicate_positions,
            'shared_nodes_at_interfaces': bool(dangling == 0 and duplicate_positions == 0),
            'basis': 'one vertex array, every tissue/tissue interface face carried by exactly two '
                     'tets through the same three vertex indices, and no coincident duplicate '
                     'vertex anywhere in the mesh'}


def interface_fidelity(TV, TT, owner, parts, det):
    """How far does the LABELLED interface sit from the input surface it is supposed to be?

    TetGen -pY preserves input facets exactly, so this is ~0. An envelope mesher does not, and this
    number is the price it charges."""
    shells = R.owned_surfaces(TV, TT, owner, len(parts))
    rows = []
    for i, p in enumerate(parts):
        SV, SF = shells[i]
        if not len(SF):
            rows.append({'entity_id': p['entity_id'], 'faces': 0}); continue
        d = np.sqrt(np.asarray(igl.point_mesh_squared_distance(
            np.ascontiguousarray(SV), np.ascontiguousarray(p['V']),
            np.ascontiguousarray(p['F']))[0], float))
        sv = abs(R.divergence(SV, SF))
        be, nm = boundary_and_nonmanifold(SV, SF)
        rows.append({'entity_id': p['entity_id'], 'name': p['name'], 'role': p['role'],
                     'vertices': int(len(SV)), 'faces': int(len(SF)),
                     'closed': be == 0, 'boundary_edges': be, 'nonmanifold_edges': nm,
                     'owned_surface_volume_m3': sv,
                     'input_surface_volume_m3': p['surface_volume_m3'],
                     'volume_relative_error':
                         float(sv / p['surface_volume_m3'] - 1) if p['surface_volume_m3'] else None,
                     'deviation_max_m': float(d.max()), 'deviation_mean_m': float(d.mean()),
                     'deviation_p99_m': float(np.percentile(d, 99))})
    live = [r for r in rows if r.get('faces')]
    return shells, {'per_structure': rows,
                    'entities': len(rows), 'entities_with_geometry': len(live),
                    'all_closed': all(r['closed'] for r in live),
                    'max_deviation_m': max((r['deviation_max_m'] for r in live), default=0.0),
                    'max_abs_volume_relative_error':
                        max((abs(r['volume_relative_error'] or 0.0) for r in live), default=0.0),
                    'note': 'deviation is the distance from every vertex of the extracted owned '
                            'surface to the nearest point of that entity\'s ORIGINAL input surface'}


def residual_overlap(shells, parts, cap=40):
    contacts = []
    total = 0.0
    boxes = [(V.min(0), V.max(0)) if len(V) else None for V, _ in shells]
    tested = 0
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            tested += 1
            if boxes[i] is None or boxes[j] is None:
                continue
            if (boxes[i][0] > boxes[j][1]).any() or (boxes[j][0] > boxes[i][1]).any():
                continue
            VA, FA = shells[i]; VB, FB = shells[j]
            if not len(FA) or not len(FB):
                continue
            n = int(len(np.asarray(cgal.intersect_other(VA, FA, VB, FB, True, False, False, False,
                                                        2_000_000)[0]).reshape(-1, 2)))
            b = cgal.mesh_boolean(VA, FA, VB, FB, type_str='intersect')
            v = abs(R.divergence(np.asarray(b[0], float), np.asarray(b[1], np.int64))) if len(b[1]) else 0.0
            total += v
            if n or v:
                contacts.append({'a': parts[i]['entity_id'], 'b': parts[j]['entity_id'],
                                 'shared_interface_face_contacts': n,
                                 'residual_overlap_volume_m3': v})
    return {'pairs_tested': tested, 'pairs_in_contact': len(contacts),
            'total_residual_overlap_volume_m3': total,
            'contacts': sorted(contacts, key=lambda c: -c['residual_overlap_volume_m3'])[:cap]}


def run_route_b(parts, pad, epsr, lr, workdir, priority_of, label, full_audit=True):
    V, _, _ = soup(parts)
    lo = V.min(0) - pad; hi = V.max(0) + pad
    box_volume = float(np.prod(hi - lo))
    rec, TV, TT = run_ftetwild(parts, workdir, epsr, lr)
    report = {'label': label, 'route': 'B-ftetwild', 'members': [p['entity_id'] for p in parts],
              'member_count': len(parts), 'epsilon_relative': epsr, 'edge_length_relative': lr,
              'ftetwild': rec}
    if not rec['succeeded']:
        return report, None
    mesh, TT, owner, det = label_and_measure(TV, TT, parts, priority_of)
    report['conforming_mesh'] = mesh
    report['shared_node_audit'] = shared_node_audit(TV, TT, owner, parts)
    if full_audit:
        shells, fid = interface_fidelity(TV, TT, owner, parts, det)
        report['interface_fidelity'] = fid
        report['resolved_disjointness'] = residual_overlap(shells, parts)
    diag = float(np.linalg.norm(hi - lo))
    report['envelope'] = {'bbox_diagonal_m': diag, 'epsilon_absolute_m': epsr * diag,
                          'note': 'fTetWild epsilon is relative to the bounding box diagonal'}
    return report, {'TV': TV, 'TT': TT, 'owner': owner}


# ------------------------------------------------------------------ self-test

def self_test():
    checks = []

    def check(name, ok, detail):
        checks.append({'check': name, 'passed': bool(ok), 'detail': detail})

    # 1. the lattice guarantees a minimum edge, and it refuses a grid it cannot represent
    W = np.array([[0.0, 0, 0], [1e-12, 0, 0], [1.0, 0, 0], [1.0000000001, 0, 0]])
    F = np.array([[0, 2, 3]], np.int64)
    W2, F2, s2, q = lattice_quantise(W, F, np.zeros(1, np.int64), 1e-6)
    d = np.linalg.norm(W2[:, None] - W2[None], axis=2)
    nz = d[d > 0]
    check('lattice weld leaves no pair closer than the grid',
          len(nz) == 0 or nz.min() >= 1e-6 * (1 - 1e-9), float(nz.min()) if len(nz) else None)
    check('lattice quantisation moves a vertex by at most grid*sqrt(3)/2',
          q['max_coordinate_move_m'] <= 1e-6 / 2 * (1 + 1e-9), q['max_coordinate_move_m'])
    try:
        lattice_quantise(np.array([[1e9, 0, 0], [1e9 + 1e-9, 0, 0]]), np.zeros((0, 3), np.int64),
                         np.zeros(0, np.int64), 1e-15)
        refused = False
    except ValueError:
        refused = True
    check('a grid below double resolution is refused, not silently welded', refused, None)

    # 2. the two-cube analytic case, meshed by route A on a 1 um lattice
    def cube(lo, hi):
        V = np.array(lo, float) + R.BOX_CORNERS * (np.array(hi, float) - np.array(lo, float))
        return np.ascontiguousarray(V), np.ascontiguousarray(R.BOX.copy())

    VA, FA = cube([0, 0, 0], [1, 1, 1]); VB, FB = cube([.75, 0, 0], [1.75, 1, 1])
    parts = [{'entity_id': 'A', 'name': 'a', 'role': 'rigid_bone', 'system': 'skeletal',
              'frame': 'test', 'source_path': 'test', 'source_sha256': '0' * 64,
              'surface_volume_m3': abs(R.divergence(VA, FA)), 'V': VA, 'F': FA},
             {'entity_id': 'B', 'name': 'b', 'role': 'muscle', 'system': 'muscular',
              'frame': 'test', 'source_path': 'test', 'source_sha256': '0' * 64,
              'surface_volume_m3': abs(R.divergence(VB, FB)), 'V': VB, 'F': FB}]
    priority = {'A': R.rank('rigid_bone'), 'B': R.rank('muscle')}
    with tempfile.TemporaryDirectory() as scratch:
        rep, m = run_route_a(parts, 1e-6, 0.25, 'pY', priority, 'self-test-A',
                             emit=Path(scratch) / 'geometry')
        a = rep['arrangement']
        check('route A produces a valid PLC', a['plc_valid'], a['plc_residual_self_intersecting_face_pairs'])
        check('route A snap converged', a['snap_converged'], len(a['snap_passes']))
        check('route A min edge is at least the grid',
              a['plc_min_edge_length_m'] >= 1e-6 * (1 - 1e-9), a['plc_min_edge_length_m'])
        check('route A has no sub-tolerance edge', a['plc_edges_under_tetgen_tolerance'] == 0,
              a['plc_edges_under_tetgen_tolerance'])
        check('route A tetgen succeeded', rep['tetgen']['succeeded'], rep['tetgen'].get('message'))
        cm = rep.get('conforming_mesh')
        check('route A all tets positive', cm and cm['all_positive_volume'], cm and cm['min_tet_volume_m3'])
        pa = [p for p in cm['per_structure'] if p['entity_id'] == 'A'][0]
        pb = [p for p in cm['per_structure'] if p['entity_id'] == 'B'][0]
        check('route A both cubes claim 1 m3 to within the volume budget',
              abs(pa['claimed_volume_m3'] - 1) < VOLUME_ERROR_BUDGET
              and abs(pb['claimed_volume_m3'] - 1) < VOLUME_ERROR_BUDGET,
              (pa['claimed_volume_m3'], pb['claimed_volume_m3']))
        check('route A disputed slab is 0.25 m3',
              abs(cm['disputed_volume_m3'] - .25) < VOLUME_ERROR_BUDGET, cm['disputed_volume_m3'])
        check('route A bone keeps all of it and muscle concedes the slab',
              abs(pa['owned_volume_m3'] - 1) < VOLUME_ERROR_BUDGET
              and abs(pb['owned_volume_m3'] - .75) < VOLUME_ERROR_BUDGET,
              (pa['owned_volume_m3'], pb['owned_volume_m3']))
        check('route A deviation is inside the stated budget',
              a['deviation_from_input']['hausdorff_two_sided_m'] <= SURFACE_DEVIATION_BUDGET_M
              or a['grid_m'] > SURFACE_DEVIATION_BUDGET_M,
              a['deviation_from_input'])
        # shared-node audit on route A's own mesh
        aud = shared_node_audit(m['TT'], m['TT'], m['owner'], parts) if False else \
            shared_node_audit(m['TV'], m['TT'], m['owner'], parts)
        check('route A interfaces share nodes', aud['shared_nodes_at_interfaces'], aud)
        check('route A has interface faces between the two structures',
              aud['interface_faces_between_two_structures'] > 0,
              aud['interface_faces_between_two_structures'])

    # 2b. the degeneracy measure catches a needle triangle a short-edge test cannot see
    nee = np.array([[0.0, 0, 0], [1.0, 0, 0], [0.5, 1e-18, 0]])
    dg = facet_degeneracy(nee, np.array([[0, 1, 2]], np.int64))
    check('facet degeneracy sees a needle whose every edge is long',
          dg['min_facet_altitude_m'] < 1e-17
          and float(np.linalg.norm(nee[1] - nee[0])) > TETGEN_SEGMENT_TOLERANCE_M
          and dg['facets_with_altitude_under_tetgen_tolerance'] == 1, dg)

    # 3. quality measure sanity on a regular tetrahedron
    reg = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
    q = tet_quality(reg, np.array([[0, 1, 2, 3]], np.int64))
    check('regular tetrahedron scores gamma 1.0', abs(q['shape_gamma_mean'] - 1.0) < 1e-9,
          q['shape_gamma_mean'])
    check('regular tetrahedron dihedral is 70.53 deg', abs(q['min_dihedral_deg'] - 70.528779) < 1e-4,
          q['min_dihedral_deg'])

    # 4. msh reader round trip
    with tempfile.TemporaryDirectory() as scratch:
        p = Path(scratch) / 't.msh'
        p.write_text('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n4\n'
                     '1 0 0 0\n2 1 0 0\n3 0 1 0\n4 0 0 1\n$EndNodes\n'
                     '$Elements\n1\n1 4 2 0 1 1 2 3 4\n$EndElements\n')
        TV, TT = read_msh_tets(p)
        check('msh reader recovers 4 nodes and 1 tet', len(TV) == 4 and len(TT) == 1 and
              TT.tolist() == [[0, 1, 2, 3]], (TV.shape, TT.tolist()))

    # 5. fTetWild presence, reported not asserted
    checks.append({'check': 'fTetWild binary present (informational)',
                   'passed': True, 'detail': {'path': str(FTETWILD), 'exists': FTETWILD.exists()}})

    passed = all(c['passed'] for c in checks)
    print(json.dumps({'self_test': 'ihm.meshing-route-study', 'passed': passed,
                      'checks': checks}, indent=2))
    return passed


# -------------------------------------------------------------------- diagnose

def diagnose(args):
    """Why does the exact-arrangement route stop where it stops?

    Three questions the existing receipts cannot answer, because the `arrangement` block they store
    is the POST-SNAP variant rather than the PLC TetGen was first handed:
      1. is the exact PLC at a failing rung measurably worse than at a succeeding one?
      2. does lattice quantisation improve or worsen the sub-tolerance features?
      3. is TetGen's verdict on one fixed PLC even deterministic?"""
    out = Path(args.out)
    rows = R.sources()
    by = {r['entity_id']: r for r in rows}
    drop = set(json.loads((REPAIR / 'coincidence.json').read_text())
               ['identical_geometry_entities_dropped']) if (REPAIR / 'coincidence.json').exists() else set()
    report = {'schema': 'ihm.meshing-ceiling-diagnosis.v1',
              'question': 'what actually stops the exact-arrangement conforming route',
              'clusters': []}
    for name in [c for c in args.clusters.split(',') if c]:
        parts = R.build_parts(by, cluster_members(name, by, drop, out / 'ladder-order.json'))
        V, F, src = soup(parts)
        W, DF, dsrc, crossings, dropped, collapsed = exact_arrangement(V, F, src)
        be, nm = boundary_and_nonmanifold(W, DF)
        _, L = R.edge_lengths(W, DF)
        rec = {'label': name, 'members': len(parts),
               'exact_plc': {'vertices': int(len(W)), 'facets': int(len(DF)),
                             'boundary_edges': be, 'nonmanifold_edges': nm,
                             'residual_self_intersecting_face_pairs': detect_crossings(W, DF),
                             'min_edge_m': float(L.min()),
                             'edges_under_tetgen_tolerance': int((L < TETGEN_SEGMENT_TOLERANCE_M).sum()),
                             **facet_degeneracy(W, DF)}}
        for grid in [float(g) for g in args.grids.split(',') if float(g) > 0]:
            try:
                Wg, Fg, sg, q = lattice_quantise(W, DF, dsrc, grid)
                bg, ng = boundary_and_nonmanifold(Wg, Fg)
                _, Lg = R.edge_lengths(Wg, Fg)
                rec['after_one_quantisation_at_%g_m' % grid] = {
                    'boundary_edges': bg, 'nonmanifold_edges': ng,
                    'residual_self_intersecting_face_pairs': detect_crossings(Wg, Fg),
                    'min_edge_m': float(Lg.min()),
                    'edges_under_tetgen_tolerance': int((Lg < TETGEN_SEGMENT_TOLERANCE_M).sum()),
                    **facet_degeneracy(Wg, Fg), **q}
            except BaseException as error:
                rec['after_one_quantisation_at_%g_m' % grid] = {'failed': str(error)}
        if args.repeats:
            lo = W.min(0) - args.pad; hi = W.max(0) + args.pad
            PV = np.ascontiguousarray(np.concatenate([lo + R.BOX_CORNERS * (hi - lo), W]))
            PF = np.ascontiguousarray(np.concatenate([R.BOX, DF + 8]).astype(np.int64))
            rec['plc_sha256'] = hashlib.sha256(PV.tobytes() + PF.tobytes()).hexdigest()
            rec['tetgen_repeatability'] = []
            for fl in [f for f in args.flags.split(',') if f]:
                runs = []
                for i in range(args.repeats):
                    t = R.run_tetgen(PV, PF, fl)
                    runs.append({'run': i, 'succeeded': bool(t.get('succeeded')),
                                 'tets': t.get('tets'), 'exception': t.get('exception'),
                                 'seconds': round(t['seconds'], 2),
                                 'warnings': t.get('tetgen_warning_kinds')})
                    print('  %s %s run %d ok=%s' % (name, fl, i, runs[-1]['succeeded']), flush=True)
                rec['tetgen_repeatability'].append(
                    {'flags': fl, 'runs': args.repeats,
                     'successes': sum(r['succeeded'] for r in runs),
                     'distinct_tet_counts': sorted({r['tets'] for r in runs if r['succeeded']}),
                     'deterministic': len({(r['succeeded'], r['tets']) for r in runs}) == 1,
                     'detail': runs})
        report['clusters'].append(rec)
        write_json(out / 'ceiling-diagnosis.json', report)
        print('[%s] done' % name, flush=True)
    return report


# ------------------------------------------------------------------- comparison

def undisplaced_reference(fid):
    """The one fidelity number that is NOT confounded by priority ownership.

    An entity that concedes volume to a higher-priority neighbour has an owned surface that is
    SUPPOSED to differ from its input surface, so a raw per-entity deviation mixes mesher error with
    the ownership decision. The lowest-rank entity in a cluster concedes nothing, so its owned
    surface must reproduce its input surface exactly. That number isolates the mesher."""
    rows = [r for r in (fid.get('per_structure') or []) if r.get('faces')]
    if not rows:
        return {}
    best = min(rows, key=lambda r: (R.rank(r.get('role') or ''), r['entity_id']))
    return {'undisplaced_reference_entity': best['entity_id'],
            'undisplaced_reference_role': best.get('role'),
            'undisplaced_reference_deviation_max_m': best.get('deviation_max_m'),
            'undisplaced_reference_volume_relative_error': best.get('volume_relative_error')}


def row(r):
    """One comparable line per (cluster, route, setting)."""
    a = r.get('arrangement') or {}
    m = r.get('conforming_mesh') or {}
    q = m.get('quality') or {}
    fid = r.get('interface_fidelity') or {}
    aud = r.get('shared_node_audit') or {}
    dis = r.get('resolved_disjointness') or {}
    ok = (r.get('ftetwild') or {}).get('succeeded') if r['route'].startswith('B') \
        else (r.get('tetgen') or {}).get('succeeded')
    route = r.get('route')
    if route == 'A-snap-round' and not (r.get('grid_m') or 0) > 0:
        route = '0-exact-arrangement-control'
    out = {'cluster': r.get('label'), 'entities': r.get('member_count'), 'route': route,
           'setting': ('grid=%g m' % r['grid_m']) if r.get('grid_m') is not None
                      else ('epsr=%g, lr=%g' % (r.get('epsilon_relative', float('nan')),
                                                r.get('edge_length_relative', float('nan')))),
           'succeeded': bool(ok), 'failure': r.get('failed'),
           'tets': m.get('tets'), 'tet_vertices': m.get('tet_vertices'),
           'seconds': r.get('wall_seconds') or (r.get('ftetwild') or {}).get('seconds'),
           'max_abs_claimed_volume_error': m.get('max_abs_claimed_volume_error'),
           'shared_nodes_at_interfaces': aud.get('shared_nodes_at_interfaces'),
           'interface_faces_between_two_structures': aud.get('interface_faces_between_two_structures'),
           'interface_max_deviation_from_input_m': fid.get('max_deviation_m'),
           'interface_max_abs_volume_relative_error': fid.get('max_abs_volume_relative_error'),
           **undisplaced_reference(fid),
           'residual_pairwise_overlap_m3': dis.get('total_residual_overlap_volume_m3'),
           'min_dihedral_deg': q.get('min_dihedral_deg'),
           'min_dihedral_p01_deg': q.get('min_dihedral_p01_deg'),
           'min_dihedral_median_deg': q.get('min_dihedral_median_deg'),
           'tets_with_min_dihedral_under_5deg': q.get('tets_with_min_dihedral_under_5deg'),
           'quality_sampled_tets': q.get('sampled_tets'),
           'shape_gamma_p01': q.get('shape_gamma_p01'),
           'min_tet_volume_m3': m.get('min_tet_volume_m3'),
           'all_positive_volume': m.get('all_positive_volume'),
           'box_volume_closure_relative_error': m.get('box_volume_closure_relative_error')}
    if a:
        out.update({'plc_residual_self_intersections': a.get('plc_residual_self_intersecting_face_pairs'),
                    'plc_boundary_edges': a.get('plc_boundary_edges'),
                    'plc_min_edge_m': a.get('plc_min_edge_length_m'),
                    'plc_edges_under_tetgen_tolerance': a.get('plc_edges_under_tetgen_tolerance'),
                    'plc_min_facet_altitude_m': (a.get('plc_degeneracy') or {}).get('min_facet_altitude_m'),
                    'plc_facets_altitude_under_tetgen_tolerance':
                        (a.get('plc_degeneracy') or {}).get('facets_with_altitude_under_tetgen_tolerance'),
                    'snap_converged': a.get('snap_converged'),
                    'snap_passes': len(a.get('snap_passes') or []),
                    'surface_deviation_two_sided_m':
                        (a.get('deviation_from_input') or {}).get('hausdorff_two_sided_m')})
    return out


def compare(args):
    out = Path(args.out)
    rows, sources_used = [], []
    for d in [x for x in args.merge.split(',') if x]:
        f = Path(d) / 'results.jsonl'
        if not f.exists():
            continue
        sources_used.append(str(f))
        for line in f.read_text().splitlines():
            if line.strip():
                rows.append(row(json.loads(line)))
    rows.sort(key=lambda r: (r['entities'] or 0, r['route'], str(r['setting'])))
    best = {}
    for r in rows:
        if r['succeeded']:
            k = r['route']
            if r['entities'] > best.get(k, {'entities': -1})['entities']:
                best[k] = r
    comparison = {'schema': 'ihm.meshing-route-comparison.v1',
                  'tolerances': {'surface_deviation_budget_m': SURFACE_DEVIATION_BUDGET_M,
                                 'volume_error_budget': VOLUME_ERROR_BUDGET,
                                 'tetgen_segment_tolerance_m': TETGEN_SEGMENT_TOLERANCE_M},
                  'sources': sources_used, 'rows': rows,
                  'largest_conforming_domain_per_route': best}
    write_json(out / 'comparison.json', comparison)
    inputs = {'data/derived/canonical/anatomy.json': sha(R.ANATOMY),
              'data/derived/muscle-tet-ready-v1/manifest.json': sha(R.MUSCLE / 'manifest.json'),
              'data/derived/entity-tet-ready-v1/manifest.json': sha(R.ENTITY / 'manifest.json'),
              'data/derived/cross-structure-repair-v1/clusters.json': sha(REPAIR / 'clusters.json'),
              'data/derived/cross-structure-repair-v1/summary.json': sha(REPAIR / 'summary.json'),
              'scripts/build_cross_structure_conflict_repair.py': sha(R.__file__),
              'scripts/build_meshing_route_study.py': sha(__file__)}
    if FTETWILD.exists():
        inputs['data/runtime/tolerant-mesher/build/FloatTetwild_bin'] = sha(FTETWILD)
        rec = FTETWILD.parents[1] / 'receipt.json'
        if rec.exists():
            inputs['data/runtime/tolerant-mesher/receipt.json'] = sha(rec)
    for f in sources_used:
        inputs[f] = sha(f)
    write_json(out / 'manifest.json',
               {'schema': 'ihm.meshing-route-study-manifest.v1', 'inputs_sha256': inputs,
                'artifacts_sha256': {p.relative_to(out).as_posix(): sha(p)
                                     for p in sorted(out.rglob('*'))
                                     if p.is_file() and p.name != 'manifest.json'
                                     and p.suffix not in ('.obj', '.msh', '.csv')},
                'canonical_assets_modified': False, 'existing_repo_files_modified': False})
    print(json.dumps({'rows': len(rows), 'largest': {k: (v['cluster'], v['entities'])
                                                     for k, v in best.items()}}, indent=2))
    return comparison


# ------------------------------------------------------------------------ run

_ORDER = None


def ladder_order(rows_by, drop, cache):
    """The repair module's spatial ordering, computed once and cached, so a ladder of many rungs
    does not reload all 2403 surfaces per rung."""
    global _ORDER
    if _ORDER is None:
        cache = Path(cache)
        if cache.exists():
            _ORDER = json.loads(cache.read_text())
        else:
            rows = [r for r in rows_by.values() if r['entity_id'] not in drop]
            _ORDER = R.ladder_members(rows, R.THIGH, [len(rows)], 10 ** 18)[0]
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(_ORDER))
    return _ORDER


def cluster_members(name, rows_by, drop, cache=None):
    if name.startswith('n') and name[1:].isdigit():
        n = int(name[1:])
        order = ladder_order(rows_by, drop, cache or (ROOT / 'data/derived/meshing-route-study-v1/ladder-order.json'))
        if n > len(order):
            raise ValueError('no rung of size %d' % n)
        return order[:n]
    data = json.loads((REPAIR / 'clusters.json').read_text())['clusters']
    for c in data:
        if c['label'] == name:
            return [e for e in c['members'] if e not in drop]
    raise ValueError('unknown cluster %r' % name)


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    rows = R.sources()
    by = {r['entity_id']: r for r in rows}
    priority = {r['entity_id']: R.rank(r['role']) for r in rows}
    coin_path = REPAIR / 'coincidence.json'
    drop = set(json.loads(coin_path.read_text())['identical_geometry_entities_dropped']) \
        if coin_path.exists() else set()
    results = []
    results_path = out / 'results.jsonl'
    if args.fresh and results_path.exists():
        results_path.unlink()
    if results_path.exists():
        results = [json.loads(l) for l in results_path.read_text().splitlines() if l.strip()]
    done = {(r['label'], r['route'], r.get('grid_m'), r.get('epsilon_relative')) for r in results}

    def record(rec):
        results.append(rec)
        with results_path.open('a') as f:
            f.write(json.dumps(rec, allow_nan=False) + '\n')

    for name in [c for c in args.clusters.split(',') if c]:
        members = cluster_members(name, by, drop, out / 'ladder-order.json')
        parts = R.build_parts(by, members)
        print('[%s] %d members, %d facets' % (name, len(parts), sum(len(p['F']) for p in parts)),
              flush=True)
        if args.route in ('a', 'both'):
            for grid in [float(g) for g in args.grids.split(',') if g]:
                if (name, 'A-snap-round', grid, None) in done:
                    continue
                print('  route A grid %g' % grid, flush=True)
                try:
                    emit = (out / 'geometry' / ('%s-A-%g' % (name, grid))) if args.emit else None
                    rep, m = run_route_a(parts, grid, args.pad, args.flags, priority, name, emit=emit)
                    if m is not None:
                        rep['shared_node_audit'] = shared_node_audit(m['TV'], m['TT'], m['owner'], parts)
                        rep['conforming_mesh']['quality'] = tet_quality(m['TV'], m['TT'])
                        _, fid = interface_fidelity(m['TV'], m['TT'], m['owner'], parts, m['det'])
                        rep['interface_fidelity'] = fid
                    rep.pop('resolved_surfaces', None)
                except BaseException as error:
                    rep = {'label': name, 'route': 'A-snap-round', 'grid_m': grid,
                           'member_count': len(parts), 'members': members,
                           'failed': type(error).__name__, 'message': str(error)}
                rep['label'] = name
                print('    tetgen succeeded: %s' % rep.get('tetgen', {}).get('succeeded'), flush=True)
                record(rep)
        if args.route in ('b', 'both'):
            for epsr in [float(e) for e in args.epsr.split(',') if e]:
                if (name, 'B-ftetwild', None, epsr) in done:
                    continue
                if not FTETWILD.exists():
                    record({'label': name, 'route': 'B-ftetwild', 'epsilon_relative': epsr,
                            'member_count': len(parts),
                            'failed': 'BinaryMissing', 'message': str(FTETWILD)})
                    continue
                print('  route B epsr %g' % epsr, flush=True)
                work = out / 'ftetwild' / ('%s-e%g' % (name, epsr))
                try:
                    rep, _ = run_route_b(parts, args.pad, epsr, args.lr, work, priority, name,
                                         full_audit=not args.fast)
                except BaseException as error:
                    rep = {'label': name, 'route': 'B-ftetwild', 'epsilon_relative': epsr,
                           'member_count': len(parts), 'members': members,
                           'failed': type(error).__name__, 'message': str(error)}
                print('    ftetwild succeeded: %s' % rep.get('ftetwild', {}).get('succeeded'),
                      flush=True)
                record(rep)

    summary = {'schema': 'ihm.meshing-route-study.v1',
               'tolerances': {'surface_deviation_budget_m': SURFACE_DEVIATION_BUDGET_M,
                              'volume_error_budget': VOLUME_ERROR_BUDGET,
                              'tetgen_segment_tolerance_m': TETGEN_SEGMENT_TOLERANCE_M},
               'clusters': args.clusters, 'grids': args.grids, 'epsr': args.epsr,
               'results': [{k: v for k, v in r.items()
                            if k not in ('members',)} for r in results],
               'wall_seconds': time.monotonic() - began, 'python': sys.version}
    write_json(out / 'summary.json', summary)
    inputs = {'data/derived/canonical/anatomy.json': sha(R.ANATOMY),
              'data/derived/muscle-tet-ready-v1/manifest.json': sha(R.MUSCLE / 'manifest.json'),
              'data/derived/entity-tet-ready-v1/manifest.json': sha(R.ENTITY / 'manifest.json'),
              'data/derived/cross-structure-repair-v1/clusters.json': sha(REPAIR / 'clusters.json'),
              'scripts/build_cross_structure_conflict_repair.py': sha(R.__file__),
              'scripts/build_meshing_route_study.py': sha(__file__)}
    if FTETWILD.exists():
        inputs['data/runtime/tolerant-mesher/build/FloatTetwild_bin'] = sha(FTETWILD)
    artifacts = {p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob('*'))
                 if p.is_file() and p.name != 'manifest.json' and p.suffix != '.obj'
                 and p.suffix != '.msh'}
    write_json(out / 'manifest.json',
               {'schema': 'ihm.meshing-route-study-manifest.v1', 'inputs_sha256': inputs,
                'per_entity_input_sha256': {r['entity_id']: r['sha256'] for r in rows},
                'artifacts_sha256': artifacts, 'canonical_assets_modified': False,
                'existing_repo_files_modified': False})
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', default=str(ROOT / 'data/derived/meshing-route-study-v1'))
    p.add_argument('--route', default='both', choices=['a', 'b', 'both'])
    p.add_argument('--clusters', default='thigh-six,ladder-50')
    p.add_argument('--grids', default=','.join('%g' % g for g in GRIDS_DEFAULT))
    p.add_argument('--epsr', default='1e-3')
    p.add_argument('--lr', type=float, default=0.05)
    p.add_argument('--pad', type=float, default=.02)
    p.add_argument('--flags', default='pY,pYT1e-14,p')
    p.add_argument('--emit', action='store_true')
    p.add_argument('--fast', action='store_true', help='skip the O(n^2) exact overlap audit')
    p.add_argument('--fresh', action='store_true')
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--merge', default='', help='comma-separated study dirs to consolidate')
    p.add_argument('--diagnose', action='store_true',
                   help='measure what actually stops the exact-arrangement route')
    p.add_argument('--repeats', type=int, default=0,
                   help='with --diagnose, how many independent TetGen runs per flag on one PLC')
    a = p.parse_args()
    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if a.diagnose:
        Path(a.out).mkdir(parents=True, exist_ok=True)
        diagnose(a)
    elif a.merge:
        compare(a)
    else:
        run(a)
