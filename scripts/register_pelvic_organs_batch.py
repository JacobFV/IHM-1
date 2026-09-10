"""register every extracted UT-EndoMRI subject into this body's pelvis, one subject at a time.

register_pelvic_organs.py takes one --case per subject: SUBJECT T2 SEGDIR ORGANDIR. this
drives it over the whole extraction. for each subject with a uterus mesh:
  1. pull the image the organ labels were drawn on -- the grid extract_ut_endomri.py recorded
     for the uterus (T2 where it exists; 92 of 124 subjects have one) -- from the VERIFIED zip,
     so bones and organs share one frame;
  2. run TotalSegmentator total_mr on it (skipped if already done: the batch resumes);
  3. register that subject alone, so one failing case never stops the others.
it ends with a count per gate. the caveat travels with every artefact: an endometriosis
cohort, pathology-selected, NOT a typical-anatomy reference.
"""
import json, re, subprocess, sys, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / "data/raw/anatomy/ut-endomri/UT-EndoMRI.zip"
ORGANS = ROOT / "data/derived/ut-endomri-organs-v1"
MR = ROOT / "data/derived/ut-endomri-mr"
OUT = ROOT / "data/derived/ut-endomri-registered-v1"
PY = ROOT / ".venv-totalseg/bin/python"; TS = ROOT / ".venv-totalseg/bin/TotalSegmentator"

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--organs", type=Path, default=ORGANS, help="an extract_ut_endomri.py output directory")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--subjects", nargs="*", default=None, help="only these subjects")
    a = ap.parse_args()
    organs_dir, out_dir = a.organs.resolve(), a.out.resolve()
    subjects = json.loads((organs_dir / "manifest.json").read_text())["subjects"]
    if a.subjects: subjects = {k: v for k, v in subjects.items() if k in set(a.subjects)}
    zf = zipfile.ZipFile(ZIP)
    members = {n.replace(" ", ""): n for n in zf.namelist() if "__MACOSX" not in n}
    out_dir.mkdir(parents=True, exist_ok=True); MR.mkdir(parents=True, exist_ok=True)
    tally = {"registered": 0, "no uterus mesh": 0, "no image": 0, "total_mr failed": 0, "registration failed": 0}
    for sid, rec in sorted(subjects.items()):
        ut = rec.get("organs", {}).get("ut", {})
        if not (organs_dir / sid / "uterus.obj").exists() or "grid" not in ut:
            tally["no uterus mesh"] += 1; continue
        grid = ut["grid"]
        key = next((k for k in members if re.search(rf"/{sid}_{grid}\.nii(\.gz)?$", k)), None)
        if key is None: print(f"{sid}: no {grid} image in the zip", flush=True); tally["no image"] += 1; continue
        d = MR / sid; d.mkdir(exist_ok=True)
        img = d / f"{grid}.nii.gz" if key.endswith(".gz") else d / f"{grid}.nii"
        if not img.exists(): img.write_bytes(zf.read(members[key]))
        seg = d / "seg_mr"
        if not (seg / "hip_left.nii.gz").exists():
            r = subprocess.run([str(TS), "-i", str(img), "-o", str(seg), "-ta", "total_mr", "--device", "gpu"],
                               capture_output=True, text=True)
            if r.returncode != 0 or not (seg / "hip_left.nii.gz").exists():
                print(f"{sid}: total_mr failed ({r.returncode})", flush=True); tally["total_mr failed"] += 1; continue
        r = subprocess.run([str(PY), "-u", str(ROOT / "scripts/register_pelvic_organs.py"),
                            "--case", sid, str(img), str(seg), str(organs_dir / sid), "--out", str(out_dir)],
                           capture_output=True, text=True)
        (d / "register.log").write_text(r.stdout + r.stderr)
        gates = [l for l in r.stdout.splitlines() if "GATE" in l or l.startswith(sid)]
        ok = r.returncode == 0
        tally["registered" if ok else "registration failed"] += 1
        print(f"{sid} [{grid}]: {'ok' if ok else 'FAILED (exit %d)' % r.returncode} | " + " | ".join(gates[-3:]), flush=True)
    print("\nTALLY:", json.dumps(tally))

if __name__ == "__main__": main()
