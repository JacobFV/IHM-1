"""per-segment atlas -> scaffold registration.

WHY.  every anatomical structure enters the scaffold through ONE similarity
(`binding.json`, 24.7 mm RMS at bone-group centroids).  measured this session:
  * the atlas's articular centres land 10.7-26.3 mm from the scaffold's own
    (same closest-approach rule on both skeletons; scripts/measure_joint_centre_registration.py)
  * translating each failing ligament by its joint's displacement alone takes
    0/37 admissible to 19/37 and the median peak strain 0.579 -> 0.170
    (scripts/measure_ligament_registration_probe.py)
  * the whole skin encloses the registered OpenSim bones at 0.888, failing at
    the hands, feet, humerus and torso (scripts/measure_skin_enclosure_whole.py)
all three are the same defect: one rigid map for a different specimen.  so fit
one similarity PER SEGMENT, atlas bone group -> the OpenSim bone mesh of the same
body at the binding's own reference pose.

HOW.  area-weighted surface samples on both meshes (the OpenSim meshes are sparse
-- femur_r has 132 vertices -- so vertex-to-vertex is not a surface distance),
initialised from the global similarity, refined by SYMMETRIC trimmed ICP: nearest
neighbours in both directions, the farthest 10% of each rejected, one Umeyama
similarity on the union.  one-way point-to-point with a free scale is degenerate
-- shrinking the source onto one spot of the target zeroes its residual -- and
the reverse correspondences are what forbid that.

GATES, each printed with the answer it must give:
  a  recovery: a cloud fitted onto a copy of itself under a known similarity must
     recover it (transform error < 1e-6, scale error < 1e-9).
  b  joint centres: the GLOBAL map must reproduce the probe's displacements
     before the per-segment ones are read.
  c  ligaments: both ends through the GLOBAL map must reproduce the stored peak
     strains (< 1e-6) before either end is moved.
  d  skin: linear blend skinning with the GLOBAL transform on every segment must
     reproduce the binding-map enclosure column before the per-segment one is read.
"""
import gzip, hashlib, importlib.util, json, sys, time
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/anatomy-segment-registration"
FREE_REGISTRATION = OUT / "registration.json"
# ROTATION MODE.  Point-to-point ICP cannot determine an elongated bone's spin about
# its own long axis, and the free fit shows it: femur 19.8 deg of rotation relative to
# the global map with only 1.3 deg off-axis, radius 32-36 deg almost all spin.  Those
# spins move a ligament's attachment round the bone, so any ligament result that rides
# on them is noise.  Two controls:
#   global    rotation held at the global map's; only translation and scale are fitted
#   no-twist  each free fit's rotation with its twist about the bone's long axis removed
#             (swing kept), then translation and scale refitted with that rotation held
import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--rotation", choices=("free", "global", "no-twist"), default="free")
_ap.add_argument("--out", type=Path, default=None)
ARGS = _ap.parse_args()
if ARGS.out is not None: OUT = ARGS.out.resolve()   # provenance paths are taken relative to ROOT
OUT.mkdir(parents=True, exist_ok=True)
N_FIT, N_FIT_TORSO, N_EVAL, TRIM, ITERS, TOL = 8000, 20000, 8000, 0.10, 400, 1e-12

def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
B = load("btfe", "scripts/build_tissue_force_elements.py")
bscm = load("bscm", "scripts/build_skin_contact_meshes.py")
bind = B.load_module("bind_anatomy", ROOT / "scripts/bind_anatomy_to_segments.py")
render = B.load_module("render_body_3d", ROOT / "scripts/render_body_3d.py")
crawl = B.load_module("crawl", ROOT / "scripts/crawl.py")
t0 = time.time()
def say(*a): print(*a, flush=True)
def jdump(o): return json.dumps(o, indent=2, default=lambda x: x.item() if hasattr(x, "item") else str(x)) + "\n"

# --------------------------------------------------------------- inputs
entities = json.loads(B.ANATOMY.read_text())["entities"]
gpath = {e["id"]: ROOT / e["reference_geometry"]["path"] for e in entities if e.get("reference_geometry")}
bind.GEOMETRY_PATH.update(gpath)
groups = {}
for e in entities:
    if e["role"] == "rigid_bone": groups.setdefault(bind.named_segment(e["name"]), []).append(e["id"])
segments = sorted(groups)

def atlas_mesh(seg):
    V, F, off = [], [], 0
    for i in groups[seg]:
        g = json.loads(gzip.decompress(gpath[i].read_bytes()))
        v = np.asarray(g["positions"], float).reshape(-1, 3); f = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        V.append(v); F.append(f + off); off += len(v)
    return np.concatenate(V), np.concatenate(F)

from ihm.spatial.vtk import surface as read_surface
GEOM = ROOT / "data/raw/anatomy/opensim-models/source/Geometry"
def osim_meshes():
    root = ET.parse(B.MODEL).getroot().find("Model"); out = {}
    for body in root.iter("Body"):
        V, F, off = [], [], 0
        for mesh in body.iter("Mesh"):
            fac = np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            p, f = read_surface(GEOM / mesh.findtext("mesh_file")); V.append(p * fac); F.append(f + off); off += len(p)
        if V: out[body.get("name")] = (np.concatenate(V), np.concatenate(F))
    return out

binding = json.loads(B.BINDING.read_text())
G = np.linalg.inv(np.asarray(binding["similarity_atlas_from_opensim_ground"], float))   # atlas -> ground, global
RG = G[:3, :3] / np.cbrt(np.linalg.det(G[:3, :3]))   # the global map's rotation
ref = binding["reference_pose_rad"]; model = render.OsimModel(B.MODEL); rest = model.forward(ref)
parent = {j["child"]: j["parent"] for j in model.joints}

def apply(M, p): return p @ M[:3, :3].T + M[:3, 3]
def sim(s, R, t): M = np.eye(4); M[:3, :3] = s * R; M[:3, 3] = t; return M
def scale_of(M): return float(np.cbrt(np.linalg.det(M[:3, :3])))

def sample(V, F, n, rng):
    tri = V[F]; a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    k = rng.choice(len(F), n, p=a / a.sum()); r1 = np.sqrt(rng.random(n)); r2 = rng.random(n); t = tri[k]
    return (1 - r1)[:, None] * t[:, 0] + (r1 * (1 - r2))[:, None] * t[:, 1] + (r1 * r2)[:, None] * t[:, 2]

def icp(src, tgt, M0, R_fixed=None):
    """symmetric trimmed ICP; returns the atlas->target similarity and iterations used.

    with R_fixed, the rotation is held and only scale and translation are solved:
    s = sum((R a_c) . b_c) / sum(|a_c|^2), t = b_mean - s R a_mean.
    """
    tt = cKDTree(tgt); M = M0.copy()
    for k in range(1, ITERS + 1):
        S = apply(M, src)
        d1, j1 = tt.query(S); d2, j2 = cKDTree(S).query(tgt)
        k1 = d1 <= np.quantile(d1, 1 - TRIM); k2 = d2 <= np.quantile(d2, 1 - TRIM)
        A = np.vstack([src[k1], src[j2[k2]]]); Bt = np.vstack([tgt[j1[k1]], tgt[k2]])
        if R_fixed is None:
            s, R, t = bind.umeyama(A, Bt)
        else:
            R = R_fixed; ac, bc = A - A.mean(0), Bt - Bt.mean(0)
            s = float(((ac @ R.T) * bc).sum() / (ac ** 2).sum()); t = Bt.mean(0) - s * R @ A.mean(0)
        Mn = sim(s, R, t)
        if np.abs(Mn - M).max() < TOL: return Mn, k
        M = Mn
    return M, ITERS

def residual(M, src, tgt):
    S = apply(M, src); d1 = cKDTree(tgt).query(S)[0]; d2 = cKDTree(S).query(tgt)[0]; d = np.concatenate([d1, d2])
    return float(np.sqrt((d ** 2).mean())), float(np.sqrt((np.sort(d)[: int(len(d) * (1 - TRIM))] ** 2).mean()))

# --------------------------------------------------------------- gate a
say("== gate a: recovery of a known similarity ==")
rng = np.random.default_rng(0)
Va, Fa = atlas_mesh("femur_l"); P = sample(Va, Fa, N_FIT, rng)
# in a fixed-rotation mode the gate exercises the fixed-rotation solver: recover a known
# scale and translation with the rotation held at the global map's
if ARGS.rotation == "free":
    M_true = sim(1.04, Rotation.from_rotvec(np.deg2rad(6) * np.array([1, 2, 3]) / np.sqrt(14)).as_matrix(), np.array([0.015, -0.010, 0.008]))
    M_fit, it = icp(P, apply(M_true, P), np.eye(4))
else:
    M_true = sim(1.04, RG, np.array([0.015, -0.010, 0.008]))
    M_fit, it = icp(P, apply(M_true, P), sim(1.0, RG, np.zeros(3)), R_fixed=RG)
terr, serr = float(np.abs(M_fit - M_true).max()), abs(scale_of(M_fit) - 1.04)
say(f"GATE a: transform error {terr:.2e} (< 1e-6), scale error {serr:.2e} (< 1e-9), {it} iterations -> "
    + ("PASS" if terr < 1e-6 and serr < 1e-9 else "FAIL"))
if not (terr < 1e-6 and serr < 1e-9):
    say("gate a failed; stopping before any registration is fitted"); sys.exit(1)

# --------------------------------------------------------------- fit
say("\n== per-segment fit (atlas bone group -> OpenSim bone mesh at the binding reference pose) ==")
om = osim_meshes()
fits = {}
Rfix = {}
if ARGS.rotation == "global":
    Rfix = {seg: RG for seg in segments}
elif ARGS.rotation == "no-twist":
    free = json.loads(FREE_REGISTRATION.read_text())["segments"]
    say(f"{'segment':10s} {'removed twist':>14s} {'kept swing':>11s}")
    for seg in segments:
        Va_, _ = atlas_mesh(seg); g_ = apply(G, Va_)
        axis = np.linalg.eigh(np.cov((g_ - g_.mean(0)).T))[1][:, -1]
        Mf = np.asarray(free[seg]["atlas_to_ground"], float)
        D = Rotation.from_matrix((Mf[:3, :3] / scale_of(Mf)) @ RG.T); q = D.as_quat()
        twist = Rotation.from_quat([*((q[:3] @ axis) * axis), q[3]]); swing = D * twist.inv()
        Rfix[seg] = swing.as_matrix() @ RG
        say(f"{seg:10s} {np.degrees(twist.magnitude()):11.1f} deg {np.degrees(swing.magnitude()):8.1f} deg")
say(f"{'segment':10s} {'scale':>7s} {'iters':>5s} {'RMS global':>11s} {'RMS per-seg':>12s} {'trimmed90 g':>12s} {'trimmed90 s':>12s}   mm")
for seg in segments:
    r = np.random.default_rng(segments.index(seg) + 1)
    Va, Fa = atlas_mesh(seg); Vo, Fo = om[seg]; Vo_g = apply(rest[seg], Vo)
    n = N_FIT_TORSO if seg == "torso" else N_FIT
    src, tgt = sample(Va, Fa, n, r), sample(Vo_g, Fo, n, r)
    if seg in Rfix:
        Rf = Rfix[seg]; s0 = scale_of(G); c = src.mean(0)
        M, it = icp(src, tgt, sim(s0, Rf, apply(G, c[None])[0] - s0 * Rf @ c), R_fixed=Rf)
    else:
        M, it = icp(src, tgt, G)
    es, et = sample(Va, Fa, N_EVAL, r), sample(Vo_g, Fo, N_EVAL, r)
    rg, rgt = residual(G, es, et); rs, rst = residual(M, es, et)
    fits[seg] = dict(M=M, iters=it, scale=scale_of(M), rms_global_m=rg, rms_segment_m=rs, trimmed_global_m=rgt, trimmed_segment_m=rst)
    say(f"{seg:10s} {scale_of(M):7.4f} {it:5d} {1000*rg:11.2f} {1000*rs:12.2f} {1000*rgt:12.2f} {1000*rst:12.2f}"
        + ("   (hit iteration cap)" if it == ITERS else ""))
say(f"global map scale {scale_of(G):.4f}")

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
inputs = {str(Path(p).relative_to(ROOT)): sha(p) for p in
          (B.ANATOMY, B.BINDING, B.MODEL, ROOT / bscm.BINDING, ROOT / bscm.EVIDENCE)}
(OUT / "registration.json").write_text(json.dumps(dict(
    schema="ihm.anatomy-segment-registration.v1", rotation_mode=ARGS.rotation, frame_from="bodyparts3d-display-m (atlas canonical)",
    frame_to="OpenSim ground at binding.json reference_pose_rad", reference_pose_rad=ref,
    method=("per-segment similarity, symmetric trimmed point-to-point ICP on area-weighted surface samples "
            "(atlas bone group -> OpenSim bone mesh of the same body), initialised from the global binding similarity"),
    parameters=dict(samples_fit=N_FIT, samples_fit_torso=N_FIT_TORSO, samples_eval=N_EVAL, trim=TRIM, max_iterations=ITERS, tolerance=TOL),
    global_atlas_to_ground=G.tolist(),
    segments={s: dict(atlas_to_ground=f["M"].tolist(), scale=f["scale"], iterations=f["iters"],
                      rms_nearest_surface_global_m=f["rms_global_m"], rms_nearest_surface_segment_m=f["rms_segment_m"],
                      trimmed90_rms_global_m=f["trimmed_global_m"], trimmed90_rms_segment_m=f["trimmed_segment_m"])
              for s, f in fits.items()},
    recovery_gate=dict(transform_error=terr, scale_error=serr, iterations=it),
    provenance=dict(git_sha=B.git_sha(), script=str(Path(__file__).resolve().relative_to(ROOT)),
                    script_sha256=sha(__file__), inputs=inputs)), indent=2) + "\n")
say(f"wrote {(OUT/'registration.json').relative_to(ROOT)}   [{time.time()-t0:.0f}s]")
Mseg = {s: f["M"] for s, f in fits.items()}
report = dict(gate_a=dict(transform_error=terr, scale_error=serr))

# --------------------------------------------------------------- gate b: joint centres + gaps
say("\n== gate b: joint centres, same closest-approach rule on both skeletons ==")
atlas_v = {s: np.concatenate([bind.read_geometry(i) for i in groups[s]]) for s in segments}
osim_g = {s: apply(rest[s], p) for s, p in bscm.bone_clouds().items()}
def centre(pa, pc, margin=0.008):
    d, j = cKDTree(pc).query(pa); ok = d <= d.min() + margin
    return ((pa[ok] + pc[j[ok]]) / 2).mean(0)
PROBE = {"tibia_l": 26.0, "tibia_r": 26.3, "talus_l": 21.4, "talus_r": 21.5, "ulna_l": 22.7, "ulna_r": 10.7, "femur_l": 19.8, "femur_r": 23.5}
JOINTS = [("femur_l", "tibia_l"), ("femur_r", "tibia_r"), ("tibia_l", "talus_l"), ("tibia_r", "talus_r"),
          ("humerus_l", "ulna_l"), ("humerus_r", "ulna_r"), ("pelvis", "femur_l"), ("pelvis", "femur_r")]
say(f"{'joint':20s} {'global':>8s} {'probe':>6s} {'per-seg':>8s}   {'gap glob':>8s} {'gap seg':>8s}   {'child-in-parent glob':>21s} {'seg':>5s}")
rows_b, gate_b = [], []
for par, ch in JOINTS:
    so = centre(osim_g[par], osim_g[ch])
    Pg, Cg = apply(G, atlas_v[par]), apply(G, atlas_v[ch])
    Ps, Cs = apply(Mseg[par], atlas_v[par]), apply(Mseg[ch], atlas_v[ch])
    dg = 1000 * np.linalg.norm(so - centre(Pg, Cg)); ds = 1000 * np.linalg.norm(so - centre(Ps, Cs))
    gg = 1000 * cKDTree(Cg).query(Pg)[0].min(); gs = 1000 * cKDTree(Cs).query(Ps)[0].min()
    Vp, Fp = atlas_mesh(par)
    ing = bscm.enclosure(apply(G, Vp), Fp, Cg); ins = bscm.enclosure(apply(Mseg[par], Vp), Fp, Cs)
    gate_b.append(abs(round(dg, 1) - PROBE[ch]))
    rows_b.append(dict(joint=f"{par}-{ch}", global_mm=dg, per_segment_mm=ds, min_gap_global_mm=gg, min_gap_segment_mm=gs,
                       child_vertices_inside_parent_global=ing, child_vertices_inside_parent_segment=ins))
    say(f"{par+'-'+ch:20s} {dg:8.1f} {PROBE[ch]:6.1f} {ds:8.1f}   {gg:8.2f} {gs:8.2f}   {ing:21.3f} {ins:5.3f}")
say(f"GATE b: max |global - probe| = {max(gate_b):.2f} mm (must be <= 0.1, the probe's print precision) -> "
    + ("PASS" if max(gate_b) <= 0.1 + 1e-9 else "FAIL"))
report["gate_b"] = dict(max_abs_global_minus_probe_mm=max(gate_b), rows=rows_b)
if max(gate_b) > 0.1 + 1e-9: say("gate b failed; stopping"); (OUT / "report.json").write_text(jdump(report)); sys.exit(1)

# --------------------------------------------------------------- gate c: ligaments
say("\n== gate c: the 51 inadmissible ligaments, each end through its own bone's map ==")
ranges = crawl.declared_ranges(B.MODEL); jc = {}
for j in model.joints:
    for c in j["coords"]:
        if c in ranges: jc.setdefault(j["child"], []).append(c)
def chain(x):
    out = [x]
    while x in parent: x = parent[x]; out.append(x)
    return out
def spanning(a, b):
    ca, cb = chain(a), chain(b); common = next(x for x in ca if x in cb)
    return [c for body in ca[:ca.index(common)] + cb[:cb.index(common)] for c in jc.get(body, [])]
cache = {}
def T(c, v):
    k = (c, round(float(v), 12))
    if k not in cache: p = dict(ref); p[c] = float(v); cache[k] = model.forward(p)
    return cache[k]
def peak(ga, gb, b1, b2):
    la = np.linalg.solve(rest[b1][:3, :3], ga - rest[b1][:3, 3]); lb = np.linalg.solve(rest[b2][:3, :3], gb - rest[b2][:3, 3])
    def L(tr): return np.linalg.norm((tr[b1][:3, :3] @ la + tr[b1][:3, 3]) - (tr[b2][:3, :3] @ lb + tr[b2][:3, 3]))
    s = L(rest)
    return max((L(T(c, v)) - s) / s for c in spanning(b1, b2) for v in np.linspace(*ranges[c], B.ADMISSIBILITY_SAMPLES))
trees = [cKDTree(atlas_v[s]) for s in segments]
elements = json.loads((B.OUT / "ligaments.json").read_text())["elements"]
bad = [r for r in elements if r.get("status") == "two_segment" and not r["kinematically_admissible"]]
probe = {r["name"]: r for r in json.loads((B.OUT / "registration_probe.json").read_text())["rows"]} \
    if (B.OUT / "registration_probe.json").exists() else {}
u = B.ULTIMATE_STRAIN; rows_c, gate_c = [], []
for r in bad:
    g = json.loads(gzip.decompress(gpath[r["id"]].read_bytes())); V = np.asarray(g["positions"], float).reshape(-1, 3)
    nearest = np.argmin(np.stack([tr.query(V)[0] for tr in trees], 1), 1)
    va, vb = V[nearest == segments.index(r["body1"])], V[nearest == segments.index(r["body2"])]
    def tip(side, other):
        rad = np.linalg.norm(side - other, axis=1); ch = side[rad >= np.percentile(rad, 100 - B.TIP_QUANTILE)]
        return (ch if len(ch) else side).mean(0)
    ta, tb = tip(va, vb.mean(0)), tip(vb, va.mean(0))
    glob = peak(apply(G, ta), apply(G, tb), r["body1"], r["body2"])
    segp = peak(apply(Mseg[r["body1"]], ta), apply(Mseg[r["body2"]], tb), r["body1"], r["body2"])
    gate_c.append(abs(glob - r["peak_strain_over_declared_range"]))
    pr = probe.get(r["name"], {}).get("corrected")
    rows_c.append(dict(name=r["name"], body1=r["body1"], body2=r["body2"], stored=r["peak_strain_over_declared_range"],
                       global_recomputed=glob, per_segment=segp, joint_translation_probe=pr))
say(f"GATE c: max |global recomputed - stored| = {max(gate_c):.2e} (< 1e-6) -> " + ("PASS" if max(gate_c) < 1e-6 else "FAIL"))
report["gate_c_max_abs"] = max(gate_c)
if max(gate_c) >= 1e-6: say("gate c failed; stopping"); (OUT / "report.json").write_text(jdump(report)); sys.exit(1)
adm = int(sum(x["per_segment"] <= u for x in rows_c))
say(f"admissible at {u}: per-segment {adm}/{len(rows_c)}   (global 0/{len(rows_c)}; one translation per joint 19/37)")
say(f"peak strain fell for {sum(x['per_segment'] < x['stored'] for x in rows_c)}/{len(rows_c)}; median "
    f"{np.median([x['stored'] for x in rows_c]):.3f} -> {np.median([x['per_segment'] for x in rows_c]):.3f}")
GROUPS = [("cruciates", "cruciate"), ("collaterals", "collateral"), ("talofibulars", "talofibular"),
          ("ligament of head of femur", "ligament of head of femur"), ("transverse acetabular", "transverse acetabular")]
say(f"\n{'element':46s} {'stored':>7s} {'joint-shift':>11s} {'per-seg':>8s}")
for label, key in GROUPS:
    say(f"-- {label}")
    for x in rows_c:
        if key in x["name"]:
            pr = x["joint_translation_probe"]
            say(f"   {x['name'][:43]:43s} {x['stored']:7.3f} {('%11.3f' % pr) if pr is not None else '          -'} {x['per_segment']:8.3f}"
                + ("  ok" if x["per_segment"] <= u else ""))
say("-- everything else")
named = [k for _, k in GROUPS]
for x in rows_c:
    if not any(k in x["name"] for k in named):
        say(f"   {x['name'][:43]:43s} {x['stored']:7.3f} {'':11s} {x['per_segment']:8.3f}" + ("  ok" if x["per_segment"] <= u else ""))
report["gate_c"] = dict(admissible_per_segment=adm, total=len(rows_c), rows=rows_c)
(OUT / "report.json").write_text(jdump(report))

# --------------------------------------------------------------- gate d: skin by linear blend skinning
say(f"\n== gate d: skin carried by linear blend skinning   [{time.time()-t0:.0f}s] ==")
mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
skin = next(e for e in mech["entities"] if e["role"] == "skin")
sg = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
SV = np.asarray(sg["positions"], float).reshape(-1, 3); SF = np.asarray(sg["indices"], np.int64).reshape(-1, 3)
ext = np.asarray(json.loads((ROOT / bscm.EVIDENCE).read_text())["contact_eligible_triangle_ids"], np.int64)
used, inv = np.unique(SF[ext], return_inverse=True); sv, sf = SV[used], inv.reshape(-1, 3)
if bscm.simtk_precondition(sv, sf) is not None:
    cv, cf = bscm.repair(sv, sf); sv, sf, _ = bscm.cap_boundaries(cv, cf)
sb = json.loads(gzip.decompress((ROOT / bscm.BINDING).read_bytes()))
skin_segs = [s["id"] for s in sb["segments"]]
W = np.asarray(sb["weights"], np.float64)[used]
W = W[cKDTree(SV[used]).query(sv)[1]]            # cap vertices take their nearest skin vertex's weights
W /= W.sum(1, keepdims=True)                      # exactly 1, so the global gate is exact
def lbs(Ms):
    out = np.zeros_like(sv)
    for k, s in enumerate(skin_segs): out += W[:, k:k + 1] * apply(Ms[s], sv)
    return out
osim_v = bscm.bone_clouds()
BINDING_MAP = dict(calcn_l=0.754, calcn_r=0.756, femur_l=1.000, femur_r=1.000, hand_l=0.634, hand_r=0.572, humerus_l=0.821,
                   humerus_r=0.674, patella_l=1.000, patella_r=1.000, pelvis=0.948, radius_l=1.000, radius_r=1.000,
                   talus_l=0.990, talus_r=0.990, tibia_l=0.994, tibia_r=0.987, toes_l=0.784, toes_r=0.796, torso=0.844,
                   ulna_l=1.000, ulna_r=1.000)
skin_g, skin_s = lbs({s: G for s in skin_segs}), lbs(Mseg)
say(f"{'segment':10s} {'binding map':>12s} {'LBS global':>11s} {'LBS per-seg':>12s}")
rows_d = []
for s in sorted(BINDING_MAP):
    bones = apply(rest[s], osim_v[s])
    eg = bscm.enclosure(skin_g, sf, bones, samples=500); es = bscm.enclosure(skin_s, sf, bones, samples=500)
    rows_d.append((s, BINDING_MAP[s], eg, es)); say(f"{s:10s} {BINDING_MAP[s]:12.3f} {eg:11.3f} {es:12.3f}")
R = np.array([x[1:] for x in rows_d])
say(f"{'mean':10s} {R[:,0].mean():12.3f} {R[:,1].mean():11.3f} {R[:,2].mean():12.3f}")
say(f"{'>= 0.99':10s} {int((R[:,0]>=.99).sum()):>9d}/22 {int((R[:,1]>=.99).sum()):>8d}/22 {int((R[:,2]>=.99).sum()):>9d}/22   (ceiling 0.997, 20/22)")
gd = float(np.abs(R[:, 1] - R[:, 0]).max())
say(f"GATE d: max |LBS global - binding map| = {gd:.4f} (<= 0.0005, the binding column's print precision) -> "
    + ("PASS" if gd <= 0.0005 + 1e-9 else "FAIL -- the LBS path changes the answer by itself; the per-segment column is not interpretable"))
report["gate_d"] = dict(max_abs_lbs_global_minus_binding=gd,
                        rows=[dict(segment=s, binding_map=b, lbs_global=g_, lbs_segment=e_) for s, b, g_, e_ in rows_d])
(OUT / "report.json").write_text(jdump(report))
say(f"wrote {(OUT/'report.json').relative_to(ROOT)}   [{time.time()-t0:.0f}s]")
