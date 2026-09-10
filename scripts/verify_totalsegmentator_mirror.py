"""is a HuggingFace copy of the TotalSegmentator zip byte-identical to Zenodo's?

zenodo's file backend is answering every range request with HTTP 504 at 30.5 s (six of
six, on both the /records and /api endpoints). a mirror is only usable if its bytes ARE
zenodo's, and that can be checked without zenodo: the central directory cached from
zenodo this morning records every member's CRC32, compressed size, size and offset. read
the mirror's central directory and compare all 147,361 entries. provenance stays zenodo;
the mirror is recorded as transport only.
"""
import importlib.util, json, struct, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fetch", ROOT / "scripts/fetch_totalsegmentator_subjects.py")
F = importlib.util.module_from_spec(spec); spec.loader.exec_module(F)
MIRROR = "https://huggingface.co/datasets/HajihajihaJimmy/TotalSegmentator_v2/resolve/main/Totalsegmentator_dataset_v201.zip"

def main():
    cached = json.loads((F.OUT / "central_directory.json").read_text())
    F.URL = MIRROR                                   # this process only; the cache is never written
    _, size = F.rng(0, 0)
    print(f"archive bytes: zenodo {cached['archive_bytes']:,}  mirror {size:,}  -> "
          + ("same" if size == cached["archive_bytes"] else "DIFFERENT"))
    if size != cached["archive_bytes"]: sys.exit("sizes differ; not the same archive")
    tail, _ = F.rng(size - 131072, size - 1)
    e = tail.rfind(b"PK\x05\x06"); n = struct.unpack("<H", tail[e+10:e+12])[0]
    cd_size, cd_off = struct.unpack("<II", tail[e+12:e+20]); z = tail.rfind(b"PK\x06\x06")
    if z >= 0: n, cd_size, cd_off = struct.unpack("<QQQ", tail[z+32:z+56])
    cd, _ = F.rng(cd_off, cd_off + cd_size - 1)
    got, p = {}, 0
    while p < len(cd) and cd[p:p+4] == b"PK\x01\x02":
        crc, csz, usz = struct.unpack("<III", cd[p+16:p+28]); nl, xl, cl = struct.unpack("<HHH", cd[p+28:p+34])
        off, = struct.unpack("<I", cd[p+42:p+46]); name = cd[p+46:p+46+nl].decode(); extra = cd[p+46+nl:p+46+nl+xl]; q = 0
        while q < len(extra):
            hid, hl = struct.unpack("<HH", extra[q:q+4])
            if hid == 1:
                it = iter(struct.unpack("<" + "Q" * (hl // 8), extra[q+4:q+4+hl]))
                if usz == 0xFFFFFFFF: usz = next(it)
                if csz == 0xFFFFFFFF: csz = next(it)
                if off == 0xFFFFFFFF: off = next(it)
            q += 4 + hl
        got[name] = (crc, csz, usz, off); p += 46 + nl + xl + cl
    want = {k: (v["crc32"], v["compressed"], v["size"], v["offset"]) for k, v in cached["members"].items()}
    diff = [k for k in want if got.get(k) != want[k]]
    extra_names = set(got) - set(want)
    print(f"members: zenodo {len(want):,}  mirror {len(got):,}  differing {len(diff)}  mirror-only {len(extra_names)}")
    ok = not diff and not extra_names and len(got) == len(want)
    print("VERDICT:", "BYTE-IDENTICAL central directory -- every member's CRC32, sizes and offset match" if ok
          else f"NOT identical (first differences: {diff[:5]})")
    sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
