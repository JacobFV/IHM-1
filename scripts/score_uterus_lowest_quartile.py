"""Lowest-quartile containment for every registered uterus -- a POST-HOC proposal, not a gate.

docs/BODY_PARAMETERS.md records that the containment gate as written ("uterus surface inside the
pelvic ring >= 0.99") fails a correctly seated ENLARGED uterus: D1-017 (726.5 mL) reads 0.688
with its lowest 25% entirely inside the ring. It proposed a lowest-quartile test "to be applied
to every subject and labelled post-hoc if adopted". This script applies it to every subject
and changes nothing: the gate, the manifests and the batch are untouched.

Same construction as register_pelvic_organs.py, imported rather than copied: this body's hip
bones and sacrum, the Delaunay ring over their vertices, 20,000 area samples of the registered
uterus with seed 0. KNOWN ANSWER: the full-surface share recomputed here must equal each
manifest's uterus_in_ring exactly, or nothing else is printed.

"Lowest" is along this body's own vertical: the axis and sign that put the calcanei below the
skull, asserted rather than assumed. The lowest q of the uterus is the q-quantile of its own
surface samples' heights.
"""
import gzip, importlib.util, json, sys
from pathlib import Path
import numpy as np
from scipy.spatial import Delaunay

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("rpo", ROOT / "scripts/register_pelvic_organs.py")
R = importlib.util.module_from_spec(spec); spec.loader.exec_module(R)
MANIFESTS = ["data/derived/ut-endomri-registered-v1/manifest.json",
             "data/derived/ut-endomri-registered-v2-world/manifest.json"]
QUANTILES = (0.25, 0.10)

def centroid(name):
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    c = [e for e in ents if e["id"].startswith("body-bp3d-") and e.get("reference_geometry") and e["name"] == name]
    if len(c) != 1: raise SystemExit(f"expected one entity named {name!r}, found {len(c)}")
    g = json.loads(gzip.decompress((ROOT / c[0]["reference_geometry"]["path"]).read_bytes()))
    return np.asarray(g["positions"], float).reshape(-1, 3).mean(0)

def vertical():
    feet = (centroid("left calcaneus") + centroid("right calcaneus")) / 2
    head = centroid("frontal bone"); d = head - feet
    axis = int(np.argmax(np.abs(d))); sign = float(np.sign(d[axis]))
    assert abs(d[axis]) > 1.2, f"calcaneus-to-skull distance along the vertical is {abs(d[axis]):.3f} m"
    assert abs(d[axis]) > 5 * np.delete(np.abs(d), axis).max(), "the body is not upright along one axis"
    return axis, sign

def main():
    body = R.body_meshes()
    ring = Delaunay(np.vstack([body["left hip bone"][0], body["right hip bone"][0], body["sacrum"][0]]))
    axis, sign = vertical()
    print(f"vertical: axis {axis}, sign {sign:+.0f} (calcanei below the frontal bone)")
    rows, mismatch = [], []
    for rel in MANIFESTS:
        man = json.loads((ROOT / rel).read_text()); base = (ROOT / rel).parent
        for sid, rec in sorted(man["subjects"].items()):
            obj = base / sid / "uterus.obj"
            if not obj.exists(): continue
            V, F = R.read_obj(obj)
            P = R.area_samples(V, F, 20000)
            inside = ring.find_simplex(P) >= 0
            full = float(inside.mean())
            if full != rec["uterus_in_ring"]: mismatch.append((rel, sid, full, rec["uterus_in_ring"]))
            h = sign * P[:, axis]
            low = {q: float(inside[h <= np.quantile(h, q)].mean()) for q in QUANTILES}
            vol = float(abs(np.einsum("ij,ij->i", V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum()) / 6e-6)
            rows.append((Path(rel).parent.name, sid, vol, full, low, rec["gate_containment"], max(rec["uterus_in_bone"].values())))
    if mismatch:
        print("KNOWN ANSWER FAILED: recomputed ring share differs from the manifest:")
        for m in mismatch: print("  ", m)
        raise SystemExit(1)
    print(f"known answer: full-surface ring share reproduced exactly for all {len(rows)} registered uteri\n")
    print(f"{'run':30s} {'subject':8s} {'volume':>9s} {'as written':>11s} {'lowest 25%':>11s} {'lowest 10%':>11s}  gate  proposal")
    for run, sid, vol, full, low, gate, bone in rows:
        prop = low[0.25] >= R.RING_MIN and bone <= R.BONE_MAX
        flag = "  <- differs" if prop != gate else ""
        print(f"{run:30s} {sid:8s} {vol:7.1f} mL {full:11.4f} {low[0.25]:11.4f} {low[0.10]:11.4f}  {'PASS' if gate else 'FAIL'}  {'PASS' if prop else 'FAIL'}{flag}")
    print("\nThe proposal (lowest 25% >= 0.99, bone <= 0.01) is POST-HOC: defined after D1-017 failed the gate "
          "as written. The gate as written remains the registered verdict for every subject.")

if __name__ == "__main__": main()
