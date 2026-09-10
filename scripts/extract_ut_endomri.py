"""uterus and ovary geometry from UT-EndoMRI, with the caveat that it is a pathology cohort.

SOURCE. UTHealth Endometriosis MRI Dataset, Zenodo record 13749613, 7,975,894,471 bytes, md5
7ace6e1b08efa10d0a1967073b0ba41c. Terms: free use exclusively in non-commercial scientific
research; cite Liang et al., "A Multi-Modal Pelvic MRI Dataset for Deep Learning-Based Pelvic
Organ Segmentation in Endometriosis". The programme owner confirmed on 2026-09-10 that this
programme qualifies. Raw and derived data stay out of git.

THE CAVEAT, which travels with every artefact. This is an ENDOMETRIOSIS cohort. The disease
distorts pelvic anatomy -- adhesions, endometriomas, adenomyosis -- and the record itself says
ovaries may be deformed or surgically absent. These are not typical-anatomy references.

LAYOUT, from the archive's own central directory: D1_MHS 51 subjects, up to three raters
(ut_r1/r2/r3, ov_r*, em_r*), rater coverage uneven; D2_TCPW 73 subjects, one rater (ut, ov,
em, cy). One label is misnamed ut_re3 (read as rater 3, flagged); three are named pat
(undocumented, excluded, flagged).

GATES, each with an answer this script does not compute:
  md5          the zip matches the md5 Zenodo records for it
  grid         every label shares EXACTLY one image's grid (shape and affine); otherwise it
               is excluded, never resampled onto a guess
  consensus    D1 uterus/ovary: a voxel is in when >= 2 of 3 raters mark it; with two raters,
               both; with one, that one. Pairwise inter-rater Dice is reported per subject.
  ovary count  connected pieces above OVARY_MIN_ML; 0-2 expected (surgical absence is real)
  volume       COARSE sanity bounds, typo catchers not anatomy: uterus UTERUS_ML, ovary OVARY_ML
"""
import argparse, hashlib, io, json, re, zipfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / "data/raw/anatomy/ut-endomri/UT-EndoMRI.zip"
OUT = ROOT / "data/derived/ut-endomri-organs-v1"
MD5 = "7ace6e1b08efa10d0a1967073b0ba41c"
CAVEAT = ("UT-EndoMRI: endometriosis cohort, pathology-selected -- NOT a typical-anatomy reference. "
          "Non-commercial research use only; cite Liang et al.")
UTERUS_ML, OVARY_ML, OVARY_MIN_ML = (10.0, 1000.0), (0.5, 400.0), 0.5
IMAGE_KINDS = ("T2", "T2FS", "T1FS", "T1")          # preference order for the reference image

def md5(p, block=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""): h.update(chunk)
    return h.hexdigest()

def load(zf, name):
    import nibabel as nib
    img = nib.Nifti1Image.from_bytes(zf.read(name)) if not name.endswith(".gz") else \
          nib.Nifti1Image.from_bytes(__import__("gzip").decompress(zf.read(name)))
    return img

def dice(a, b):
    s = a.sum() + b.sum(); return float(2 * (a & b).sum() / s) if s else float("nan")

def mesh(mask, affine):
    from skimage import measure
    v, f, _, _ = measure.marching_cubes(np.pad(mask.astype(np.uint8), 1), 0.5)
    v = v - 1
    return (v @ affine[:3, :3].T + affine[:3, 3]) / 1000.0, f

def write_obj(path, v, f, header):
    with open(path, "w") as h:
        h.write(f"# {header}\n# {CAVEAT}\n")
        for x in v: h.write("v %.7g %.7g %.7g\n" % tuple(x))
        for t in f: h.write("f %d %d %d\n" % (t[0] + 1, t[1] + 1, t[2] + 1))

def main():
    from scipy import ndimage
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-md5", action="store_true", help="only when the download step already verified it")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    if not a.skip_md5:
        got = md5(ZIP); print(f"GATE md5: {got} {'PASS' if got == MD5 else 'FAIL'}")
        if got != MD5: raise SystemExit("the zip is not the archive Zenodo records; nothing extracted")
    zf = zipfile.ZipFile(ZIP)
    pat = re.compile(r"UT-EndoMRI/(D[12]_\w+)/(D\d-\d+)/\s*(D\d-\d+)_\s*(\w+)\.nii(\.gz)?$")
    subjects = {}
    for n in zf.namelist():
        if "__MACOSX" in n or n.endswith(".DS_Store"): continue
        m = pat.match(n.replace(" ", ""))
        if m: subjects.setdefault(m.group(2), {})[m.group(4)] = n
    OUT.mkdir(parents=True, exist_ok=True)
    report = dict(schema="ihm.ut-endomri-organs.v1", source="Zenodo 13749613", md5=MD5, caveat=CAVEAT,
                  consensus_rule=">=2 of 3 raters; both of 2; the one of 1", subjects={})
    for i, (sid, files) in enumerate(sorted(subjects.items())):
        if a.limit is not None and i >= a.limit: break
        flags = []
        if "pat" in files: flags.append("'pat' label present: undocumented, excluded"); files.pop("pat")
        if "ut_re3" in files: flags.append("'ut_re3' read as rater 3"); files.setdefault("ut_r3", files.pop("ut_re3"))
        images = {k: load(zf, files[k]) for k in IMAGE_KINDS if k in files}
        def grid_of(img):
            for k, im in images.items():
                if im.shape[:3] == img.shape[:3] and np.allclose(im.affine, img.affine, atol=1e-3): return k
            return None
        rec = dict(flags=flags, images=sorted(images), organs={})
        for organ in ("ut", "ov"):
            raters = sorted(k for k in files if re.fullmatch(organ + r"(_r\d)?", k))
            masks, grids = {}, set()
            for k in raters:
                img = load(zf, files[k]); g = grid_of(img)
                if g is None: flags.append(f"{k}: matches no image grid -- excluded"); continue
                masks[k] = np.asanyarray(img.dataobj) > 0; grids.add(g); ref = img
            if not masks: continue
            if len(grids) != 1: flags.append(f"{organ}: raters on different grids {sorted(grids)} -- excluded"); continue
            stack = np.stack(list(masks.values()))
            need = 2 if len(masks) >= 2 else 1
            cons = stack.sum(0) >= (need if len(masks) != 2 else 2)
            pair = {f"{x}~{y}": dice(masks[x], masks[y]) for ix, x in enumerate(masks) for y in list(masks)[ix + 1:]}
            vox_ml = abs(np.linalg.det(ref.affine[:3, :3])) / 1000.0
            o = dict(grid=next(iter(grids)), raters=sorted(masks), inter_rater_dice=pair)
            d = OUT / sid; d.mkdir(parents=True, exist_ok=True)
            if organ == "ut":
                vol = float(cons.sum() * vox_ml); o["volume_ml"] = vol
                o["volume_in_bound"] = UTERUS_ML[0] <= vol <= UTERUS_ML[1]
                if cons.any(): v, f = mesh(cons, ref.affine); write_obj(d / "uterus.obj", v, f, f"{sid} uterus, MRI world frame, metres")
            else:
                lab, k = ndimage.label(cons, structure=np.ones((3, 3, 3), bool))
                sizes = np.bincount(lab.ravel())[1:] * vox_ml
                keep = [j + 1 for j in np.argsort(-sizes) if sizes[j] >= OVARY_MIN_ML]
                o["pieces_ml"] = [round(float(sizes[j - 1]), 2) for j in keep]
                o["count"] = len(keep); o["count_expected"] = len(keep) <= 2
                o["volume_in_bound"] = all(OVARY_ML[0] <= sizes[j - 1] <= OVARY_ML[1] for j in keep)
                for n_, j in enumerate(keep[:2], 1):
                    v, f = mesh(lab == j, ref.affine); write_obj(d / f"ovary_{n_}.obj", v, f, f"{sid} ovary piece {n_}, MRI world frame, metres")
            rec["organs"][organ] = o
        report["subjects"][sid] = rec
        u = rec["organs"].get("ut", {}); ov = rec["organs"].get("ov", {})
        print(f"{sid}: uterus {u.get('volume_ml', float('nan')):7.1f} mL (raters {u.get('raters', [])}) | "
              f"ovaries {ov.get('count', '-')} {ov.get('pieces_ml', [])} | flags {flags or '-'}", flush=True)
        (OUT / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{len(report['subjects'])} subjects written to {OUT.relative_to(ROOT)}")

if __name__ == "__main__": main()
