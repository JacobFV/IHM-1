"""Bake a gimbal-sized low-poly copy of this body's own outer envelope.

The orientation widget shows THIS body, not a generic figure, so it is
decimated from data/derived/outer-envelope/outer-envelope.json.gz rather than
modelled. Only the display copy is written; the canonical asset is untouched.

    .venv/bin/python scripts/bake_gimbal_body.py [--faces 480]
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "data/derived/outer-envelope/outer-envelope.json.gz"
TARGET = ROOT / "app/src/gimbal-body.json"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--faces", type=int, default=480)
    options = parser.parse_args()

    raw = SOURCE.read_bytes()
    envelope = json.loads(gzip.decompress(raw))
    if envelope["representation"] != "triangular_surface":
        raise SystemExit("Outer envelope is not a triangular surface")
    points = np.asarray(envelope["positions"], dtype=np.float64).reshape(-1, 3)
    faces = np.asarray(envelope["indices"], dtype=np.int64).reshape(-1, 3)

    import fast_simplification

    keep = options.faces / len(faces)
    if not 0 < keep < 1:
        raise SystemExit("Target face count must be below the source face count")
    out_points, out_faces = fast_simplification.simplify(
        points.astype(np.float32), faces.astype(np.int32), target_reduction=1 - keep
    )
    out_points = np.asarray(out_points, dtype=np.float64)
    out_faces = np.asarray(out_faces, dtype=np.int64)

    # Centre on the bounding box and scale so the tallest axis spans 1 unit; the
    # widget only needs orientation, never absolute size.
    low, high = out_points.min(axis=0), out_points.max(axis=0)
    centre = (low + high) / 2
    scale = float((high - low).max())
    normalized = (out_points - centre) / scale

    payload = {
        "schema": "ihm.gimbal-body.v1",
        "builder": "scripts/bake_gimbal_body.py",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "frame": envelope["frame"],
            "units": envelope["units"],
            "source_faces": int(len(faces)),
        },
        "display_faces": int(len(out_faces)),
        "normalization": {
            "centre_m": [round(float(v), 6) for v in centre],
            "scale_m": round(scale, 6),
            "note": "Orientation only; the widget carries no measurement.",
        },
        "positions": [round(float(v), 4) for v in normalized.reshape(-1)],
        "indices": [int(v) for v in out_faces.reshape(-1)],
    }
    TARGET.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    print(
        f"{len(faces)} -> {len(out_faces)} faces, {len(out_points)} vertices, "
        f"{TARGET.stat().st_size / 1024:.1f} kB at {TARGET.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    main()
