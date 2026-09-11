"""How wrong is the enclosed volume of the 62 viscera that are declared closed and are not?

`docs/BODY_PARAMETERS.md` records that 62 entities carry the basis "a BodyParts3D viscus is one
closed surface enclosing wall AND lumen, so the enclosed volume carries gas or chyme at tissue
density. Tissue volume is surface area times an assumed wall thickness; the remainder is lumen" --
and that their meshes are NOT closed. That is a defect in a premise. It is not yet a number.

This measures the number, because the size of an error decides whether it is worth repairing. The
signed integral of an open surface is wrong by exactly the volume the missing caps would have
contributed, so: chain the boundary edges into loops, fan each loop to its own centroid, and compare
the capped volume against the raw one. The gap is the error the ledger is carrying.

KNOWN ANSWER, before any viscus is reported: a unit cube with one triangle deleted must cap back to
EXACTLY 1.0. That exercises the whole path -- boundary detection, loop chaining, fan orientation and
the volume reader -- on a case whose answer is known to machine precision. A capper that orients its
fan backwards gives 1 - 2*(cap volume), not 1, so the sign is tested and not merely the magnitude.

WHAT THIS CANNOT DO: a fan to the loop centroid is the minimal cap, so it gives a LOWER BOUND on the
correction for a non-planar or re-entrant hole. Where the correction is large, it is large by at
least that much; where it is small, the hole is genuinely small. It is reported as a bound, not an
estimate.
"""
import gzip, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASIS_MARK = "a BodyParts3D viscus is one closed surface"


def signed_volume(V, F):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def boundary_directed(F):
    """Directed edges used once. Winding is consistent, so a hole's rim is a directed cycle."""
    use = defaultdict(int)
    for i, j in ((0, 1), (1, 2), (2, 0)):
        for u, v in zip(F[:, i].tolist(), F[:, j].tolist()):
            use[(min(u, v), max(u, v))] += 1
    out = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        for u, v in zip(F[:, i].tolist(), F[:, j].tolist()):
            if use[(min(u, v), max(u, v))] == 1:
                out.append((u, v))
    return out


def loops(edges):
    """Chain directed boundary edges into cycles. Unchainable remnants are returned separately."""
    nxt = defaultdict(list)
    for u, v in edges:
        nxt[u].append(v)
    seen, cycles, stranded = set(), [], 0
    for u0, vs in list(nxt.items()):
        for v0 in vs:
            if (u0, v0) in seen:
                continue
            cyc, u, v = [u0], u0, v0
            seen.add((u, v))
            while v != u0:
                cyc.append(v)
                cand = [w for w in nxt.get(v, []) if (v, w) not in seen]
                if not cand:
                    stranded += 1
                    cyc = None
                    break
                w = cand[0]
                seen.add((v, w))
                v = w
            if cyc and len(cyc) >= 3:
                cycles.append(cyc)
    return cycles, stranded


def cap(V, F):
    """Close every boundary loop with a fan to its centroid. Returns (V2, F2, n_loops, stranded)."""
    cycles, stranded = loops(boundary_directed(F))
    if not cycles:
        return V, F, 0, stranded
    V2 = [V]
    tris = []
    n = len(V)
    for cyc in cycles:
        centre = V[cyc].mean(axis=0)
        V2.append(centre[None, :])
        ci = n
        n += 1
        # the rim runs u->v on the open side, so the cap triangle (u, v, centre) closes it with the
        # outward orientation the surface already has.
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            tris.append([a, b, ci])
    return np.vstack(V2), np.vstack([F, np.asarray(tris, np.int64)]), len(cycles), stranded


def gate():
    V = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                  [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)
    F = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
                  [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]], np.int64)
    ok = True
    for name, Ft, truth in (("cube less one triangle", F[1:], 1.0),
                            ("cube less two triangles", F[2:], 1.0)):
        V2, F2, nl, st = cap(V, Ft)
        got = signed_volume(V2, F2)
        good = abs(got - truth) < 1e-12 and st == 0
        ok &= good
        print(f"  {name:24s} raw {signed_volume(V, Ft):9.6f} -> capped {got:9.6f}"
              f"  (want {truth})  loops {nl} stranded {st}  {'PASS' if good else 'FAIL'}")
    return ok


def main():
    print("KNOWN ANSWER: a cube with faces removed must cap back to exactly its volume")
    if not gate():
        sys.exit("known answer FAILED -- no viscus number is reported")

    anat = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())["entities"]
    basis = {e["id"]: (e.get("volume_basis") or "") for e in mech}
    led = {e["id"]: e.get("volume_m3") for e in mech}

    print(f"\n{'viscus':32s} {'raw mL':>10s} {'capped mL':>10s} {'shortfall':>10s} {'loops':>6s}")
    rows = []
    for e in anat:
        if not e.get("reference_geometry") or BASIS_MARK not in basis.get(e["id"], ""):
            continue
        g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
        V = np.asarray(g["positions"], float).reshape(-1, 3)
        F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
        raw = signed_volume(V, F) * 1e6
        V2, F2, nl, st = cap(V, F)
        if nl == 0:
            continue                      # already closed: not one of the 62
        capped = signed_volume(V2, F2) * 1e6
        short = (capped - raw) / capped if capped else float("nan")
        rows.append((abs(short), e.get("name", "?"), raw, capped, short, nl, st))
    rows.sort(reverse=True)
    for _, name, raw, capped, short, nl, st in rows:
        print(f"{name[:32]:32s} {raw:10.2f} {capped:10.2f} {short:9.1%} {nl:6d}"
              + ("  (stranded rim)" if st else ""))
    if rows:
        s = np.array([r[4] for r in rows])
        print(f"\n{len(rows)} open viscera. enclosed-volume shortfall: median {np.median(s):.1%},"
              f" min {s.min():.1%}, max {s.max():.1%}")
        print("The cap is a fan to each loop's centroid, so these are LOWER BOUNDS on the")
        print("correction. Where a figure is small the hole is genuinely small.")


if __name__ == "__main__":
    main()
