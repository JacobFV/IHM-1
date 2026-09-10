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

WORLD ROUTE (--route world), added AFTER the strict route's result and recorded as such. Across
D1, 35 of 83 uterus rater labels were drawn on MRI series not in the release, and some subjects
split raters across two released grids. Rater pairs compared in scanner world space
(world_frame_agreement.json) agree across released grids (uterus 0.798, ovary 0.626) as well as
on one grid (0.812, 0.602); off-release uterus labels mostly agree (0.708) with gross failures;
off-release ovary labels do not (0.162). So, in this route:
  * raters on different RELEASED grids are merged onto one reference grid (T2 where present) by
    their affines, nearest neighbour;
  * an off-release UTERUS label is accepted only when an on-grid uterus label of the same subject
    confirms it at world Dice >= PER_PAIR (0.406 = the committed criterion, 0.5 x the same-grid
    median 0.812, applied per pair);
  * an off-release OVARY label is never recovered.
The strict route is the default and is unchanged.
"""
import argparse, hashlib, io, json, re, zipfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / "data/raw/anatomy/ut-endomri/UT-EndoMRI.zip"
OUT_DEFAULT = ROOT / "data/derived/ut-endomri-organs-v1"
MD5 = "7ace6e1b08efa10d0a1967073b0ba41c"
CAVEAT = ("UT-EndoMRI: endometriosis cohort, pathology-selected -- NOT a typical-anatomy reference. "
          "Non-commercial research use only; cite Liang et al.")
UTERUS_ML, OVARY_ML, OVARY_MIN_ML = (10.0, 1000.0), (0.5, 400.0), 0.5
IMAGE_KINDS = ("T2", "T2FS", "T1FS", "T1")          # preference order for the reference image
PER_PAIR = 0.406                                     # world route: 0.5 x the same-grid uterus median 0.812

def into(src_img, src_mask, ref_img):
    """nearest neighbour: each voxel of ref_img's grid takes the src value at the same world point"""
    M = np.linalg.inv(src_img.affine) @ ref_img.affine
    idx = np.indices(ref_img.shape[:3]).reshape(3, -1).T
    s_ = np.rint(idx @ M[:3, :3].T + M[:3, 3]).astype(int)
    ok = np.all((s_ >= 0) & (s_ < np.array(src_mask.shape)), axis=1)
    out = np.zeros(len(idx), bool); out[ok] = src_mask[s_[ok, 0], s_[ok, 1], s_[ok, 2]]
    return out.reshape(ref_img.shape[:3])

def world_gather(organ, raters, files, images, grid_of, zf, flags):
    """the world route's masks for one organ, all on one reference grid, or None"""
    on, off = {}, {}
    for k in raters:
        img = load(zf, files[k]); g = grid_of(img); m = np.asanyarray(img.dataobj) > 0
        if m.any(): (on if g else off)[k] = (img, m, g)
    if not on:
        if off: flags.append(f"{organ}: only off-release labels; nothing on the release can confirm them -- excluded")
        return None
    grids_on = sorted({v[2] for v in on.values()})
    ref_grid = "T2" if "T2" in grids_on else grids_on[0]
    ref = images[ref_grid]; out = {}
    for k, (img, m, g) in on.items():
        out[k] = m if g == ref_grid else into(img, m, ref)
        if g != ref_grid: flags.append(f"{k}: on released grid {g}, merged onto {ref_grid} in world space")
    for k, (img, m, g) in off.items():
        if organ != "ut":
            flags.append(f"{k}: off the release -- ovary labels are not recovered"); continue
        best, by = -1.0, None
        for kk, (img2, m2, g2) in on.items():
            d = dice(m2, into(img, m, img2))
            if d > best: best, by = d, kk
        if best >= PER_PAIR:
            out[k] = into(img, m, ref); flags.append(f"{k}: off the release, confirmed by {by} at Dice {best:.3f}, merged onto {ref_grid}")
        else:
            flags.append(f"{k}: off the release, NOT confirmed (best Dice {best:.3f} < {PER_PAIR}) -- excluded")
    return out, ref, ref_grid

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
    ap.add_argument("--zip", type=Path, default=ZIP, help="the archive (default: the verified download)")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--route", choices=("strict", "world"), default="strict")
    a = ap.parse_args()
    out = a.out
    if not a.skip_md5:
        got = md5(a.zip); print(f"GATE md5: {got} {'PASS' if got == MD5 else 'FAIL'}")
        if got != MD5: raise SystemExit("the zip is not the archive Zenodo records; nothing extracted")
    zf = zipfile.ZipFile(a.zip)
    pat = re.compile(r"UT-EndoMRI/(D[12]_\w+)/(D\d-\d+)/\s*(D\d-\d+)_\s*(\w+)\.nii(\.gz)?$")
    subjects = {}
    for n in zf.namelist():
        if "__MACOSX" in n or n.endswith(".DS_Store"): continue
        m = pat.match(n.replace(" ", ""))
        if m: subjects.setdefault(m.group(2), {})[m.group(4)] = n
    out.mkdir(parents=True, exist_ok=True)
    report = dict(schema="ihm.ut-endomri-organs.v1", route=a.route, source="Zenodo 13749613", md5=MD5, caveat=CAVEAT,
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
            if a.route == "world":
                got = world_gather(organ, raters, files, images, grid_of, zf, flags)
                if got is None or not got[0]: continue
                masks, ref, gname = got; grids = {gname}
            else:
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
            d = out / sid; d.mkdir(parents=True, exist_ok=True)
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
        (out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{len(report['subjects'])} subjects written to {out}")

if __name__ == "__main__": main()
