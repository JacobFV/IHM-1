"""How much of this body has a volume at all? A watertightness census over every mesh.

`scripts/audit_organ_volumes.py` found that seven of ten solid organs are not closed surfaces, so
no volume exists for them. Ten organs cannot say whether that is an organ problem or a BODY
problem, and the difference matters: the mass ledger in `mechanics.json` allocates mass from tissue
VOLUMES, and a volume computed by the divergence theorem on an open surface returns a number rather
than an error. If open meshes are widespread, every mass in that ledger inherits the defect.

This script asks the question of all 4,000 entities at once, and reports by tissue class so an
answer like "the viscera are open and the bones are closed" is visible rather than averaged away.

KNOWN ANSWERS, run before any census figure is printed:
  1. a closed cube reports 0 boundary edges;
  2. the SAME cube with one triangle deleted reports exactly 3 -- the detector must find a hole of
     known size, not merely agree that a good mesh is good. A census tool that only ever sees closed
     meshes cannot distinguish "all closed" from "detector broken", which is the failure this
     control exists to exclude;
  3. a mesh with one duplicated triangle reports exactly 3 non-manifold edges and 0 open ones.

THE THIRD KNOWN ANSWER WAS WRITTEN WRONG THE FIRST TIME, and the detector was right. It expected
"2 open, 1 non-manifold" from a duplicated face. Duplicating triangle (0,2,1) raises each of its
three edges from two uses to three, so the correct answer is 3 non-manifold and 0 open -- which is
what the detector printed. The EXPECTATION was corrected, not the detector, and the census had not
been run, so nothing was tuned to a result. Recorded because a known answer that fails is worth
exactly as much as the care taken to find out which side of it was wrong.

NON-MANIFOLD IS REPORTED SEPARATELY FROM OPEN. An edge used by three or more triangles is a
different defect from an edge used by one, and lumping them would hide which repair is needed.
"""
import gzip, json, sys
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANATOMY = ROOT / "data/derived/canonical/anatomy.json"


def edge_use(F):
    """Counts of how many triangles use each undirected edge, computed without a Python loop."""
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]], axis=0)
    e.sort(axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return counts


def classify(F):
    c = edge_use(F)
    return int((c == 1).sum()), int((c >= 3).sum())


def cube():
    V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    F = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                  [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], np.int64)
    return V, F


def gate():
    ok = True
    _, F = cube()
    for name, Ft, want_open, want_nm in (
        ("closed cube", F, 0, 0),
        ("cube less one triangle", F[1:], 3, 0),
        ("cube with a duplicated face", np.vstack([F, F[:1]]), 0, 3),
    ):
        o, n = classify(Ft)
        good = (o == want_open and n == want_nm)
        ok &= good
        print(f"  {name:28s} open edges {o:3d} (want {want_open})  non-manifold {n:3d}"
              f" (want {want_nm})  {'PASS' if good else 'FAIL'}")
    return ok


def main():
    print("KNOWN ANSWERS: the detector must find a hole of known size, not just bless a good mesh")
    if not gate():
        sys.exit("known answer FAILED -- no census figure is printed")

    ents = json.loads(ANATOMY.read_text())["entities"]
    geo = [e for e in ents if e.get("reference_geometry")]
    print(f"\nentities carrying geometry: {len(geo)}")

    tally = Counter(); by_class = {}
    worst = []
    for e in geo:
        try:
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        except Exception:
            tally["unreadable"] += 1
            continue
        if len(F) == 0:
            tally["empty"] += 1
            continue
        o, n = classify(F)
        state = "closed" if (o == 0 and n == 0) else ("non-manifold" if n else "open")
        tally[state] += 1
        cls = (e.get("tissue_class") or e.get("system") or "unlabelled")
        d = by_class.setdefault(cls, Counter()); d[state] += 1
        if o: worst.append((o, e.get("name", "?")))

    total = sum(tally.values())
    print(f"\n{'state':16s} {'count':>7s} {'share':>8s}")
    for k in ("closed", "open", "non-manifold", "empty", "unreadable"):
        if tally[k]: print(f"{k:16s} {tally[k]:7d} {tally[k]/total:8.1%}")

    print(f"\nby tissue class / system (classes with >= 20 entities):")
    print(f"{'class':34s} {'n':>6s} {'closed':>8s}")
    for cls, d in sorted(by_class.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(d.values())
        if n >= 20: print(f"{str(cls)[:34]:34s} {n:6d} {d['closed']/n:8.1%}")

    worst.sort(reverse=True)
    print(f"\nlargest holes, by count of boundary edges:")
    for o, name in worst[:12]: print(f"  {o:6d}  {name}")

    closed = tally["closed"]
    print(f"\n{closed} of {total} meshes ({closed/total:.1%}) are closed surfaces and have a volume.")
    print("Any mass or volume property computed over the rest is reading a number the divergence")
    print("theorem returns for an open surface, which is not that surface's volume.")


if __name__ == "__main__":
    main()
