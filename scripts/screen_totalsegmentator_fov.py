"""screen coverage candidates by their CT's actual field of view, read from the NIfTI
header inside the zip -- one small HTTP Range request per subject, no CT download.

why: of four subjects chosen by label coverage, three failed on field of view. an
in-plane extent of 330-380 mm cuts the breasts off laterally and a short scan cuts
the ribs. both are in the NIfTI header (dim, pixdim), which sits in the first 348
bytes of the decompressed member, so the first few KB of the gzip stream are enough.

writes fov_selection.json: candidates with in-plane >= --min-inplane-mm and
superior-inferior >= --min-si-mm, ages closest to 40, in the same 'candidates'
format fetch_totalsegmentator_subjects.py --subjects-from takes.
"""
import argparse, json, struct, zlib, importlib.util
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fetch", ROOT / "scripts/fetch_totalsegmentator_subjects.py")
F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)

def header(cd, sid):
    m = cd[f"{sid}/ct.nii.gz"]
    head, _ = F.rng(m["offset"], m["offset"] + 29); nl, xl = struct.unpack("<HH", head[26:30])
    start = m["offset"] + 30 + nl + xl
    # a PARTIAL deflate stream: one-shot zlib.decompress refuses a stream cut mid-way
    # ("incomplete or truncated stream"), which failed every subject on the first run.
    # a decompressobj yields whatever the chunk holds, and 64 KB is ample for 348 bytes.
    chunk, _ = F.rng(start, start + min(m["compressed"], 65536) - 1)
    raw = zlib.decompressobj(-15).decompress(chunk) if m["method"] == 8 else chunk   # the zip layer
    nii = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw, 400)                # the .gz layer
    if len(nii) < 348: raise ValueError(f"only {len(nii)} header bytes decoded")
    endian = "<" if struct.unpack("<i", nii[:4])[0] == 348 else ">"
    dim = struct.unpack(endian + "8h", nii[40:56]); pixdim = struct.unpack(endian + "8f", nii[76:108])
    return [dim[i + 1] for i in range(3)], [abs(pixdim[i + 1]) for i in range(3)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-inplane-mm", type=float, default=420.0)
    ap.add_argument("--min-si-mm", type=float, default=450.0)
    ap.add_argument("--exclude", nargs="*", default=["s0897", "s0790", "s1157", "s1067"])
    ap.add_argument("--want", type=int, default=None,
                    help="stop once this many candidates pass, taking them in the selection's order")
    a = ap.parse_args()
    cd = F.central_directory()["members"]
    cands = json.loads((F.OUT / "coverage_selection.json").read_text())["candidates"]
    kept, rows, unreadable = [], [], []
    def screen(sid):
        dims, pix = header(cd, sid)
        ext = [d * p for d, p in zip(dims, pix)]
        ok = min(ext[0], ext[1]) >= a.min_inplane_mm and ext[2] >= a.min_si_mm
        rows.append(dict(subject=sid, dims=dims, spacing_mm=pix, extent_mm=ext, passes=ok))
        print(f"  {sid}: {ext[0]:.0f} x {ext[1]:.0f} x {ext[2]:.0f} mm  {'KEEP' if ok else '-'}", flush=True)
        if ok: kept.append(sid)
    for sid in cands:
        if sid in a.exclude: continue
        if a.want is not None and len(kept) >= a.want: break
        try: screen(sid)
        except Exception as e: print(f"  {sid}: header unreadable ({e}) -- retried in a second pass", flush=True); unreadable.append(sid)
    # a server timeout is NOT a field-of-view failure: the first version dropped these silently
    still = []
    for sid in unreadable:
        if a.want is not None and len(kept) >= a.want: still.append(sid); continue
        try: screen(sid)
        except Exception as e: print(f"  {sid}: still unreadable ({e})", flush=True); still.append(sid)
    unreadable = still
    out = F.OUT / "fov_selection.json"
    out.write_text(json.dumps(dict(rule=f"coverage candidates with in-plane >= {a.min_inplane_mm} mm and SI >= {a.min_si_mm} mm",
                                   screened=rows, unreadable=unreadable, candidates=kept), indent=2) + "\n")
    print(f"{len(kept)} of {len(rows)} screened pass; wrote {out.relative_to(ROOT)}")

if __name__ == "__main__": main()
