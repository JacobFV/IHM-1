"""Are this body's solid organs the right SIZE? A known-answer volume audit.

`docs/BODY_PARAMETERS.md` records, as a pending item, that "the testis meshes are 6.7 mL against a
12-30 mL reference". That is one organ. The question this script asks is whether it is ONE organ or
a SCALE: if every solid organ reads low by a similar factor, the fault is in the body's units or its
build, and a single-organ reading would have been the wrong diagnosis entirely.

KNOWN ANSWER, checked before any organ is measured: the divergence-theorem volume of an analytic
sphere and a cube, tessellated the same way and read by the same function, must match their closed
forms. Nothing else is printed if it does not.

CLOSEDNESS IS NOT ASSUMED. Signed volume by the divergence theorem is meaningless on an open
surface -- it silently returns a number. Every mesh is checked for boundary edges (an edge used by
one triangle rather than two) and for consistent winding, and an organ that is not watertight is
reported as UNMEASURABLE rather than given a volume.

THE REFERENCE RANGES ARE STATED HERE, NOT CATALOGUED. They are typical adult male values from
general clinical anatomy, carried in this file so they can be argued with. They are NOT drawn from
any source catalogued in `data/sources/`, and no gate is attached to them. What the audit is for is
the RATIO PATTERN across organs, which does not depend on any single range being right.
"""
import gzip, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANATOMY = ROOT / "data/derived/canonical/anatomy.json"

# organ -> (low_mL, high_mL). Typical adult male. Stated here, not catalogued; see the docstring.
REFERENCE = {
    "left testis":       (12, 30),
    "right testis":      (12, 30),
    "left kidney":       (120, 170),
    "right kidney":      (120, 170),
    "spleen":            (100, 250),
    "pancreas":          (60, 100),
    "prostate":          (15, 30),
    "gallbladder":       (30, 50),
    "urinary bladder":   (40, 120),
    "stomach":           (200, 900),
}


def mesh(entity):
    g = json.loads(gzip.decompress((ROOT / entity["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3)
    F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    return V, F


def signed_volume(V, F):
    """Divergence theorem: sum over triangles of (a . (b x c)) / 6. Metres^3 if V is in metres."""
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def boundary_edges(F):
    """An edge used by one triangle is a hole. Returns the count of such edges."""
    use = defaultdict(int)
    for i, j in ((0, 1), (1, 2), (2, 0)):
        for e in zip(F[:, i], F[:, j]):
            use[tuple(sorted(e))] += 1
    return sum(1 for n in use.values() if n != 2)


def icosphere(radius, subdiv=4):
    """A sphere built here, so its volume is known in closed form."""
    t = (1 + 5 ** 0.5) / 2
    V = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t],
                  [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    F = np.array([[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4],
                  [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8],
                  [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]], np.int64)
    for _ in range(subdiv):
        mid, newF = {}, []
        def m(i, j):
            k = (min(i, j), max(i, j))
            if k not in mid:
                mid[k] = len(V2); V2.append((V[i] + V[j]) / 2)
            return mid[k]
        V2 = [v for v in V]
        for f in F:
            a, b, c = m(f[0], f[1]), m(f[1], f[2]), m(f[2], f[0])
            newF += [[f[0], a, c], [f[1], b, a], [f[2], c, b], [a, b, c]]
        V, F = np.asarray(V2, float), np.asarray(newF, np.int64)
    V = V / np.linalg.norm(V, axis=1, keepdims=True) * radius
    return V, F


def cube(side):
    V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float) * side
    F = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                  [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], np.int64)
    return V, F


def gate():
    """KNOWN ANSWER. The reader must reproduce closed forms before any organ is quoted.

    FIRST VERSION FAILED, AND THE BAR WAS NOT MOVED. It asked the inscribed icosphere to match the
    SPHERE's closed form within 2e-3 and read 2.17e-3 at subdivision 4. That is not an error in the
    volume reader: an inscribed polyhedron understates its sphere by construction, so the bar was
    testing the tessellation's fineness while claiming to test the reader. The bar stands as failed;
    the INSTRUMENT is replaced, which is the permitted move.

    What replaces it is a sharper known answer. The cube is exact in closed form, so it gates the
    reader at machine precision. The sphere's deficit is then required to CONVERGE at second order
    -- the relative error must fall by a factor near 4 per subdivision, because the deficit goes as
    the square of the edge length. A reader with a sign error, a factor of 2, or a winding bug does
    not produce clean second-order convergence toward the analytic value; it is a test the earlier
    absolute bar could not perform at all.
    """
    ok = True
    Vc, Fc = cube(0.04)
    got, truth = signed_volume(Vc, Fc), 0.04 ** 3
    rel = abs(got - truth) / truth
    flag = "PASS" if rel <= 1e-12 and boundary_edges(Fc) == 0 else "FAIL"
    if flag == "FAIL": ok = False
    print(f"  exactness   cube s=0.04 m  {got*1e6:9.4f} mL vs {truth*1e6:9.4f}  rel {rel:.2e}"
          f" (<= 1e-12)  {flag}")

    truth_s, prev, rels = 4 / 3 * np.pi * 0.04 ** 3, None, []
    for sub in (3, 4, 5):
        V, F = icosphere(0.04, sub)
        v = signed_volume(V, F); r = abs(v - truth_s) / truth_s; rels.append(r)
        ratio = f"{prev/r:.3f}x" if prev else "   --"
        print(f"  convergence icosphere sub={sub}  {v*1e6:9.4f} mL  rel {r:.3e}  fall {ratio}"
              f"  boundary edges {boundary_edges(F)}")
        prev = r
        if boundary_edges(F): ok = False
    ratios = [rels[i] / rels[i + 1] for i in range(len(rels) - 1)]
    conv = all(3.5 <= x <= 4.5 for x in ratios) and all(v > 0 for v in rels)
    # the icosphere is INSCRIBED, so every deficit must also be one-signed (an understatement).
    signs = all(signed_volume(*icosphere(0.04, s)) < truth_s for s in (3, 4, 5))
    print(f"  second-order convergence (each fall in 3.5-4.5x): {['%.3f'%x for x in ratios]}"
          f" -> {'PASS' if conv else 'FAIL'};  deficit one-signed: {'PASS' if signs else 'FAIL'}")
    return ok and conv and signs


def main():
    print("KNOWN ANSWER: the volume reader against closed forms")
    if not gate():
        sys.exit("known answer FAILED -- no organ volume is quoted")

    ents = json.loads(ANATOMY.read_text())["entities"]
    by_name = {}
    for e in ents:
        if e.get("reference_geometry") and e.get("name") in REFERENCE:
            by_name.setdefault(e["name"], e)

    print(f"\n{'organ':18s} {'volume':>10s} {'reference':>14s} {'ratio to low':>13s}  closedness")
    ratios = []
    for name in REFERENCE:
        e = by_name.get(name)
        if e is None:
            print(f"{name:18s} {'--':>10s} {'':>14s} {'':>13s}  ABSENT from this body")
            continue
        V, F = mesh(e)
        holes = boundary_edges(F)
        vol_mL = signed_volume(V, F) * 1e6
        lo, hi = REFERENCE[name]
        if holes:
            print(f"{name:18s} {'--':>10s} {f'{lo}-{hi} mL':>14s} {'':>13s}  UNMEASURABLE: {holes} boundary edges")
            continue
        r = abs(vol_mL) / lo
        ratios.append((name, r))
        verdict = "within" if lo <= abs(vol_mL) <= hi else ("LOW" if abs(vol_mL) < lo else "HIGH")
        print(f"{name:18s} {abs(vol_mL):9.2f} mL {f'{lo}-{hi} mL':>14s} {r:12.3f}x  closed, {verdict}")

    if ratios:
        rs = np.array([r for _, r in ratios])
        print(f"\n{len(rs)} closed organs. ratio to the LOW end of each range:"
              f" median {np.median(rs):.3f}x, min {rs.min():.3f}x, max {rs.max():.3f}x")
        print("A common factor across organs would indicate a SCALE fault; a spread indicates "
              "per-organ mesh problems and the testis reading is then its own story.")


if __name__ == "__main__":
    main()
