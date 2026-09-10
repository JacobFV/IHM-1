"""which ligaments DRIVE hip flexion out of its declared range in tissue-force-elements-v2?

measure_ligament_moments.py showed v2 fixes the knee (v1 drove it out of range at up to
101 N.m; v2 restores on every out-of-range sample) and introduces a hip-flexion DRIVE:
0 restoring / 18 driving samples on the left, 0 / 9 on the right, to 37 N.m, where v1
mostly restored. This breaks that moment down per element, v1 against v2.

A line element under tension F makes a generalized force -F dl/dq on q. Evaluated per
element at 10 and 30 deg past each end of hip_flexion's declared range, everything else
held at the binding's reference pose -- the same arithmetic, the same spring law and
the same held pose as measure_ligament_moments.py.

GATE: at each evaluated angle the per-element moments must sum to what the moment
check itself computes there, or this breakdown is not of the same quantity.
"""
import importlib.util, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
LM = load("measure_ligament_moments", "scripts/measure_ligament_moments.py")
render = load("render_body_3d", "scripts/render_body_3d.py")
crawl = load("crawl", "scripts/crawl.py")
MODEL = ROOT / "data/models/engineering_stance_v1/model.osim"
binding = json.loads((ROOT / "data/derived/anatomy-segment-binding/binding.json").read_text())
pose0 = binding["reference_pose_rad"]
model = render.OsimModel(MODEL)
ranges = crawl.declared_ranges(MODEL)
H = 1e-5

def elements(version):
    rows = [r for r in json.loads((ROOT / f"data/derived/tissue-force-elements-{version}/ligaments.json").read_text())["elements"]
            if r.get("status") == "two_segment"]
    for r in rows: r["_p1"], r["_p2"] = np.asarray(r["point1_m"]), np.asarray(r["point2_m"])
    return rows

def per_element(rows, coord, q):
    p = dict(pose0); p[coord] = q
    lp = LM.lengths(model, rows, {**p, coord: q + H}); lm = LM.lengths(model, rows, {**p, coord: q - H}); l0 = LM.lengths(model, rows, p)
    dl = (lp - lm) / (2 * H)
    F = np.array([LM.spring_force(l0[i], r["slack_length_m"], r["linear_stiffness_n"], r["transition_strain"]) for i, r in enumerate(rows)])
    return -F * dl, F, (l0 - np.array([r["slack_length_m"] for r in rows])) / np.array([r["slack_length_m"] for r in rows])

def main():
    total_check = {}
    for version in ("v1", "v2"):
        rows = elements(version)
        print(f"\n################ {version} ################")
        for coord in ("hip_flexion_l", "hip_flexion_r"):
            lo, hi = ranges[coord]
            for side, sign, edge in (("above high", +1, hi), ("below low", -1, lo)):
                for past in (10, 30):
                    q = edge + sign * np.radians(past)
                    m, F, strain = per_element(rows, coord, q)
                    want = -sign                                   # restoring moment points back toward the range
                    tot = float(m.sum()); total_check[(version, coord, round(q, 6))] = tot
                    verdict = "RESTORES" if tot * want > 0 else ("DRIVES" if abs(tot) > 1e-6 else "nothing")
                    print(f"{coord} {past:2d} deg {side:10s}: total {tot:+8.2f} N.m  {verdict}")
                    order = np.argsort(-np.abs(m))[:4]
                    for i in order:
                        if abs(m[i]) < 0.05: continue
                        r = rows[i]; role = "restores" if m[i] * want > 0 else "DRIVES"
                        pel = "  [pelvis]" if "pelvis" in (r["body1"], r["body2"]) else ""
                        print(f"     {r['name'][:44]:44s} {m[i]:+8.2f} N.m  {role:8s} F {F[i]:7.1f} N  strain {strain[i]:+.3f}{pel}")
    # GATE against the moment check's own curve at matching angles
    worst = 0.0; compared = 0
    for version in ("v1", "v2"):
        rep = LM.__dict__.get("_last_report")
        path = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if path is not None:
        for version, f in (("v1", path / "moments_v1.json"), ("v2", path / "moments_v2.json")):
            d = json.loads(f.read_text())
            coords = next(v for v in d.values() if isinstance(v, list) and v and isinstance(v[0], dict) and "coordinate" in v[0])
            for c in coords:
                if c["coordinate"] not in ("hip_flexion_l", "hip_flexion_r"): continue
                rows = elements(version)
                for smp in c["samples"][::10]:
                    m, _, _ = per_element(rows, c["coordinate"], smp["value_rad"])
                    worst = max(worst, abs(float(m.sum()) - smp["ligament_moment_nm"])); compared += 1
        print(f"\nGATE: per-element sum vs the moment check's own total, {compared} samples: max |diff| = {worst:.2e} N.m -> "
              + ("PASS" if worst < 1e-3 else "FAIL -- this is not the same quantity"))

if __name__ == "__main__": main()
