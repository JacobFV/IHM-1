"""fetch SELECTED subjects out of the TotalSegmentator v2.0.1 archive by HTTP Range.

the archive is one 23.58 GB zip of 147,361 members.  a female breast needs a handful
of female chest CTs, not the corpus, so read the central directory once, choose
subjects by a stated rule from the archive's OWN meta.csv, and range-fetch only
their members.  every member is checked against the CRC32 the zip stores for it --
the integrity gate, with an answer this script never computes itself.

licence: CC-BY-4.0 (Zenodo record 10047292).  cite Wasserthal et al., Radiology: AI
2023, doi 10.1148/ryai.230024.  data lands under data/raw/ (gitignored).
"""
import argparse, csv, hashlib, io, json, struct, zlib, urllib.request
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
URL = "https://zenodo.org/records/10047292/files/Totalsegmentator_dataset_v201.zip?download=1"
OUT = ROOT / "data/raw/anatomy/totalsegmentator"
# the members registration needs: the SAME 39 labels build_female_torso_from_totalsegmentator.py
# registers by.  a subject carries 118 members; the breast and skin are produced by running
# subtask models on the CT, so the other 78 (sex-neutral organs) are not needed for this step.
REGISTRATION_LABELS = ([f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)]
                       + ["sternum", "clavicula_left", "clavicula_right"]
                       + [f"vertebrae_T{i}" for i in range(1, 13)])

def rng(a, b):
    req = urllib.request.Request(URL, headers={"Range": f"bytes={a}-{b}"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                if r.status != 206: raise SystemExit(f"server ignored Range (HTTP {r.status})")
                return r.read(), int(r.headers["Content-Range"].split("/")[1])
        except (OSError, TimeoutError) as e:
            if attempt == 4: raise
            print(f"  retry {attempt+1}: {e}", flush=True)

def central_directory():
    cache = OUT / "central_directory.json"
    if cache.exists(): return json.loads(cache.read_text())
    _, size = rng(0, 0); tail, _ = rng(size - 131072, size - 1)
    e = tail.rfind(b"PK\x05\x06"); n = struct.unpack("<H", tail[e+10:e+12])[0]
    cd_size, cd_off = struct.unpack("<II", tail[e+12:e+20]); z = tail.rfind(b"PK\x06\x06")
    if z >= 0: n, cd_size, cd_off = struct.unpack("<QQQ", tail[z+32:z+56])
    cd, _ = rng(cd_off, cd_off + cd_size - 1); members, p = {}, 0
    while p < len(cd) and cd[p:p+4] == b"PK\x01\x02":
        meth, = struct.unpack("<H", cd[p+10:p+12]); crc, csz, usz = struct.unpack("<III", cd[p+16:p+28])
        nl, xl, cl = struct.unpack("<HHH", cd[p+28:p+34]); off, = struct.unpack("<I", cd[p+42:p+46])
        name = cd[p+46:p+46+nl].decode(); extra = cd[p+46+nl:p+46+nl+xl]; q = 0
        while q < len(extra):
            hid, hl = struct.unpack("<HH", extra[q:q+4])
            if hid == 1:
                it = iter(struct.unpack("<" + "Q"*(hl//8), extra[q+4:q+4+hl]))
                if usz == 0xFFFFFFFF: usz = next(it)
                if csz == 0xFFFFFFFF: csz = next(it)
                if off == 0xFFFFFFFF: off = next(it)
            q += 4 + hl
        members[name] = dict(method=meth, crc32=crc, compressed=csz, size=usz, offset=off); p += 46 + nl + xl + cl
    if len(members) != n: raise SystemExit(f"central directory parsed {len(members)} of {n} members")
    OUT.mkdir(parents=True, exist_ok=True); cache.write_text(json.dumps(dict(archive_bytes=size, members=members)))
    return dict(archive_bytes=size, members=members)

def fetch(name, m, dest):
    # already on disk and the zip's own CRC32 agrees: nothing to fetch
    if dest.exists() and dest.stat().st_size == m["size"]:
        raw = dest.read_bytes()
        if (zlib.crc32(raw) & 0xFFFFFFFF) == m["crc32"]:
            return hashlib.sha256(raw).hexdigest()
    head, _ = rng(m["offset"], m["offset"] + 29); nl, xl = struct.unpack("<HH", head[26:30])
    start = m["offset"] + 30 + nl + xl
    data, _ = rng(start, start + m["compressed"] - 1) if m["compressed"] else (b"", 0)
    raw = zlib.decompress(data, -15) if m["method"] == 8 else data
    if len(raw) != m["size"] or (zlib.crc32(raw) & 0xFFFFFFFF) != m["crc32"]:
        raise SystemExit(f"CRC/size mismatch on {name}: the zip's own checksum refuses this member")
    dest.parent.mkdir(parents=True, exist_ok=True); dest.write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--subjects", type=int, default=4)
    ap.add_argument("--target-age", type=float, default=40.0)
    ap.add_argument("--members", choices=("all", "registration"), default="all",
                    help="'registration': only ct.nii.gz and the 39 registration labels")
    a = ap.parse_args()
    cd = central_directory()["members"]
    if not (OUT / "meta.csv").exists(): fetch("meta.csv", cd["meta.csv"], OUT / "meta.csv")
    rows = list(csv.DictReader(io.open(OUT / "meta.csv", encoding="utf-8-sig"), delimiter=";"))
    # THE RULE, stated so it can be argued with: female; no pathology recorded; a study
    # that spans neck to pelvis (so the thorax is whole AND the ribs, sternum,
    # clavicles and vertebrae this body shares are all labelled for registration);
    # then the ages closest to --target-age, because the median female age here is 65
    # and a breast is not age-invariant.
    pool = [r for r in rows if r["gender"] == "f" and r["pathology"] == "no_pathology"
            and "neck" in r["study_type"] and "thorax" in r["study_type"] and "pelvis" in r["study_type"]
            and r["age"] not in ("", "nan")]
    pool.sort(key=lambda r: abs(float(r["age"]) - a.target_age))
    chosen = pool[:a.subjects]
    print(f"{len(pool)} subjects satisfy the rule; taking {len(chosen)}:")
    for r in chosen: print(f"  {r['image_id']}  age {float(r['age']):.0f}  {r['study_type']}  split {r['split']}")
    manifest = dict(source="totalsegmentator", zenodo_record="10047292", version="2.0.1", licence="CC-BY-4.0",
                    citation="Wasserthal et al., Radiology: Artificial Intelligence 2023, doi:10.1148/ryai.230024",
                    selection_rule="female; pathology == no_pathology; study_type spans neck, thorax and pelvis; ages closest to %g" % a.target_age,
                    pool_size=len(pool), members=a.members, subjects={})
    for r in chosen:
        sid = r["image_id"]; names = sorted(n for n in cd if n.startswith(sid + "/") and not n.endswith("/"))
        if a.members == "registration":
            keep = {f"{sid}/ct.nii.gz"} | {f"{sid}/segmentations/{l}.nii.gz" for l in REGISTRATION_LABELS}
            if not keep <= set(names): raise SystemExit(f"{sid}: archive lacks {sorted(keep - set(names))}")
            names = sorted(keep)
        total = sum(cd[n]["size"] for n in names); print(f"{sid}: {len(names)} members, {total/1e6:.0f} MB", flush=True)
        files = {n: fetch(n, cd[n], OUT / n) for n in names}
        manifest["subjects"][sid] = dict(meta=r, members=len(names), bytes=total, sha256=files)
        print(f"  {sid}: all {len(names)} members pass the zip's CRC32", flush=True)
        (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("done:", OUT / "manifest.json")

if __name__ == "__main__": main()
