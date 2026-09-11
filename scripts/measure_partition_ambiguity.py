#!/usr/bin/env python3
"""How much does the hard partition throw away, and where?

`build_skin_contact_meshes.py` cuts the canonical skin into one rigid piece per segment by taking
the ARGMAX of `continuous_surface_binding`'s graph-diffused skinning weights -- and says so in its
own docstring: "a HARD partition of a surface that is really continuous, so every segment boundary
is a seam that in the real body does not exist", forced because "Simbody is a rigid multibody engine
and there is no deformable continuum anywhere in it".

The weights are continuous; the partition is not. A triangle whose winning weight is 0.98 is
unambiguously one segment's. One at 0.36 is being assigned by a coin-toss and then carried rigidly
as though it belonged. This measures the size of that population and which segments carry it, which
is a property of the partition alone -- no warp, no solver, no gate.

KNOWN ANSWERS, both checked before any figure is printed:
  1. the winning weight is the row maximum of a partition of unity, so it can never be below
     1/n_segments and never above 1. Violations mean the weights are not what they claim;
  2. a synthetic hard one-hot binding must report ambiguity exactly zero -- a detector that cannot
     distinguish a genuinely hard binding from a soft one measures nothing.
"""
import gzip, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BINDING = ROOT / "data/derived/canonical/continuous_surface_binding.json.gz"
BARS = (0.5, 0.6, 0.8)


def winning(weights, faces):
    """Per triangle: the mean vertex weight row, its winner and the winner's share."""
    tri = (weights[faces[:, 0]] + weights[faces[:, 1]] + weights[faces[:, 2]]) / 3.0
    owner = tri.argmax(axis=1)
    return owner, tri[np.arange(len(tri)), owner]


def gate(n_seg):
    """FIRST VERSION OF THIS GATE WAS WRONG AND THE DETECTOR WAS RIGHT.

    It made each VERTEX one-hot on an independently drawn segment, so a triangle's three vertices
    were hot on three different segments and the mean row's maximum was 1/3 -- correctly. A hard
    binding means the three vertices of a triangle agree, not that each vertex is confident. The
    corrected case assigns one-hot PER TRIANGLE.

    That failure is itself the finding the metric needed: a low winning share conflates two things --
    weights that are genuinely soft, and a triangle that straddles a boundary between two confident
    vertices. The second is reported separately below as the spanning count, which needs no
    threshold.
    """
    ok = True
    rng = np.random.default_rng(0)
    n_tri = 100
    f = np.arange(3 * n_tri).reshape(-1, 3)
    hot = np.zeros((3 * n_tri, n_seg), np.float32)
    per_tri = rng.integers(0, n_seg, n_tri)
    hot[np.arange(3 * n_tri), np.repeat(per_tri, 3)] = 1.0
    _, share = winning(hot, f)
    hard_ok = float(share.min()) == 1.0
    ok &= hard_ok
    print(f"  one-hot PER TRIANGLE -> minimum winning share {share.min():.4f} (want 1.0)"
          f"  {'PASS' if hard_ok else 'FAIL'}")
    split = np.zeros((3 * n_tri, n_seg), np.float32)          # every triangle straddles a boundary
    split[np.arange(3 * n_tri), np.repeat(per_tri, 3)] = 1.0
    third = np.arange(2, 3 * n_tri, 3)                        # explicit ROW indices, not a slice:
    split[third] = 0.0                                        # `split[2::3, idx]` mixes a slice with
    split[third, (per_tri + 1) % n_seg] = 1.0                 # an array and BROADCASTS, setting
                                                              # len(idx) columns in every selected
                                                              # row instead of pairing them one to
                                                              # one. That made every third vertex
                                                              # hot on all 22 segments and the gate
                                                              # read 1.0 where 2/3 is right.
    _, s3 = winning(split, f)
    span_ok = abs(float(s3.max()) - 2 / 3) < 1e-6
    ok &= span_ok
    print(f"  confident vertices straddling a seam -> winning share {s3.max():.4f} (want 0.6667)"
          f"  {'PASS' if span_ok else 'FAIL'}  <- why a low share is not by itself soft weights")
    flat = np.full((300, n_seg), 1.0 / n_seg, np.float32)
    _, s2 = winning(flat, f)
    flat_ok = abs(float(s2.min()) - 1.0 / n_seg) < 1e-6
    ok &= flat_ok
    print(f"  uniform binding -> winning share {s2.min():.4f} (want {1/n_seg:.4f})"
          f"  {'PASS' if flat_ok else 'FAIL'}")
    return ok


def main():
    b = json.loads(gzip.decompress(BINDING.read_bytes()))
    segments = b["segments"]
    weights = np.asarray(b["weights"], np.float32)
    # the binding carries weights and reference positions but no topology; the faces come from the
    # canonical skin entity it was built on. KNOWN ANSWER: that entity's vertices must reproduce the
    # binding's own reference_positions_m, or the two are not the same mesh and the weights do not
    # index these triangles.
    ents = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
    want = set(b["surface_entity_ids"])
    V, faces, off = [], [], 0
    for e in ents:
        if e["id"] in want and e.get("reference_geometry"):
            g = json.loads(gzip.decompress((ROOT / e["reference_geometry"]["path"]).read_bytes()))
            v = np.asarray(g["positions"], float).reshape(-1, 3)
            f = np.asarray(g["indices"], np.int64).reshape(-1, 3)
            V.append(v); faces.append(f + off); off += len(v)
    V = np.vstack(V); faces = np.vstack(faces)
    ref = np.asarray(b["reference_positions_m"], float).reshape(-1, 3)
    if V.shape != ref.shape or not np.allclose(V, ref, atol=1e-9):
        d = np.abs(V - ref).max() if V.shape == ref.shape else float("nan")
        sys.exit(f"skin entity does not reproduce the binding's reference positions "
                 f"(shapes {V.shape} vs {ref.shape}, max |dv| {d}) -- not the same mesh")
    print(f"known answer: skin entity reproduces the binding's reference positions exactly "
          f"({len(V)} vertices)")
    print(f"binding: {weights.shape[0]} skin vertices x {len(segments)} segments, {len(faces)} triangles")

    rows = weights.sum(1)
    print(f"partition of unity: row sums min {rows.min():.6f} max {rows.max():.6f}")
    print("KNOWN ANSWERS")
    if not gate(len(segments)):
        sys.exit("known answer FAILED -- no ambiguity figure is printed")

    owner, share = winning(weights, faces)
    if share.min() < 1.0 / len(segments) - 1e-6 or share.max() > 1.0 + 1e-6:
        sys.exit(f"winning share out of range [{1/len(segments):.4f}, 1]: "
                 f"{share.min():.4f}..{share.max():.4f} -- the weights are not a partition of unity")

    print(f"\nwinning share over all {len(faces)} triangles: "
          f"median {np.median(share):.3f}, p10 {np.percentile(share, 10):.3f}, min {share.min():.3f}")
    for bar in BARS:
        n = int((share < bar).sum())
        print(f"  triangles assigned on a winning share below {bar}: {n} ({n/len(faces):.2%})")

    # the threshold-free half: a triangle whose three vertices do not agree on their argmax is
    # straddling a seam, whatever the weights' confidence. No bar to choose.
    vote = weights.argmax(1)
    spanning = (vote[faces[:, 0]] != vote[faces[:, 1]]) | (vote[faces[:, 1]] != vote[faces[:, 2]])
    print(f"\ntriangles whose three vertices disagree on their own argmax (a seam, no threshold): "
          f"{int(spanning.sum())} ({spanning.mean():.2%})")
    vmax = weights.max(1)
    print(f"per-VERTEX confidence, independent of topology: median {np.median(vmax):.3f}, "
          f"p10 {np.percentile(vmax, 10):.3f}, below 0.6 {(vmax < 0.6).mean():.2%}")

    print(f"\n{'segment':12s} {'triangles':>10s} {'median share':>13s} {'< 0.6':>8s} {'< 0.6 %':>9s}")
    worst = []
    for i, seg in enumerate(segments):
        name = seg["id"] if isinstance(seg, dict) else str(seg)
        m = owner == i
        if not m.any():
            continue
        sh = share[m]
        amb = int((sh < 0.6).sum())
        worst.append((amb / len(sh), name, len(sh), float(np.median(sh)), amb))
    for frac, name, n, med, amb in sorted(worst, key=lambda r: -r[0]):
        print(f"{name:12s} {n:10d} {med:13.3f} {amb:8d} {frac:8.1%}")

    print("\nA triangle with a low winning share is assigned to one segment and carried rigidly by")
    print("it, while the binding says it belongs partly to another. That is the seam the real body")
    print("does not have, measured rather than described.")


if __name__ == "__main__":
    main()
