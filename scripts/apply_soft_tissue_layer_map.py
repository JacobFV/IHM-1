"""Give every skin contact patch this body's own soft-tissue depth and an in vivo modulus.

The rule was pre-registered in docs/SEGMENT_CONTACT_SURFACES.md ("Depth alone would make it
worse", commit 97eba07) before the engine could carry per-segment stiffness:

  h = the segment's median skin-to-nearest-bone/muscle depth (soft-tissue-depth-v1), a vertex
      belonging to its argmax binding weight;
  E = the in vivo APPARENT compressive modulus -- the median of three heel-pad values for calcn,
      the Linder-Ganz buttock-fat secant for pelvis (transferred) and for every other segment
      (unsourced for that site);
  k = E/h, with no confined-layer factor, because the heel moduli are pressure over thickness
      strain (checked below against Yang's own pressure-strain pair).

The meshes are copied byte for byte; only the manifest gains a `layer` per record and a
`layer_map` block.  Nothing here is tuned: every number is read from the depth map or the card.
"""
import argparse, gzip, json, shutil, statistics
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
CARD = ROOT / "data/sources/in-vivo-soft-tissue-compression.json"
DEPTH = ROOT / "data/derived/soft-tissue-depth-v1/depth.npz"
BINDING = ROOT / "data/derived/canonical/continuous_surface_binding.json.gz"
HEEL_SEGMENTS = ("calcn_l", "calcn_r")

def moduli():
    card = json.loads(CARD.read_text()); teng, yang, linder = card["entries"]
    assert "Teng" in teng["citation"] and "Yang" in yang["citation"] and "Linder-Ganz" in linder["citation"]
    # Gefen et al. 2001 is carried only as Teng's citation ("as high as 175 kPa"); it is the third
    # in vivo heel value the pre-registration named.
    assert "175 kPa" in teng["also_cites"]
    heel = [175.0, teng["values"]["elastic_modulus_kpa"]["median"], yang["values"]["elastic_modulus_kpa"]["median"]]
    # The reading that licenses k = E/h: Yang's pressure over strain must reproduce Yang's modulus.
    v = yang["values"]; secant = v["peak_heel_pressure_kpa"]["mean"] / v["peak_gait_strain"]["mean"]
    agree = abs(secant - v["elastic_modulus_kpa"]["median"]) / v["elastic_modulus_kpa"]["median"]
    assert agree <= 0.05, f"heel modulus is not a layer modulus ({secant:.1f} vs {v['elastic_modulus_kpa']['median']})"
    f = linder["values"]
    fat = f["fat_peak_principal_compressive_stress_kpa"]["mean"] / f["fat_peak_principal_compressive_strain"]["mean"]
    return 1e3 * statistics.median(heel), 1e3 * fat, dict(heel_values_kpa=heel, yang_secant_kpa=secant,
                                                           yang_secant_vs_reported=agree, fat_secant_kpa=fat)

def segment_depths():
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3)
    b = json.loads(gzip.decompress(BINDING.read_bytes()))
    segments = [s["id"] for s in b["segments"]]; owner = np.asarray(b["weights"], np.float32).argmax(1)
    z = np.load(DEPTH); gap, vertex = cKDTree(V).query(z["skin_points_m"])
    assert gap.max() == 0.0, "depth points are not vertices of this skin"
    seg = owner[vertex]
    return {name: (float(np.median(z["depth_m"][seg == i])), int((seg == i).sum())) for i, name in enumerate(segments)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", default="data/derived/segment-contact-meshes/skin")
    ap.add_argument("--out", default="data/derived/segment-contact-meshes/skin-layer-map-v1")
    args = ap.parse_args()
    src, out = ROOT / args.bundle, ROOT / args.out
    if out.exists(): raise SystemExit(f"{out} exists; a layer map is written once")
    manifest = json.loads((src / "manifest.json").read_text())
    assert manifest["layer"] == "skin"
    heel_pa, fat_pa, checks = moduli(); depth = segment_depths()
    for r in manifest["records"]:
        h, n = depth[r["body"]]
        if n < 100: raise SystemExit(f"{r['body']}: {n} depth points, too few for a median")
        heel = r["body"] in HEEL_SEGMENTS
        E = heel_pa if heel else fat_pa
        r["layer"] = dict(thickness_m=h, depth_points=n, apparent_modulus_pa=E, stiffness_pa_per_m=E / h,
            modulus_basis="in vivo heel pad, median of Gefen 2001 / Teng 2022 / Yang 2022" if heel else
            ("Linder-Ganz 2007 buttock-fat secant, transferred" if r["body"] == "pelvis" else
             "Linder-Ganz 2007 buttock-fat secant, UNSOURCED for this site"))
        print(f"  {r['body']:10s} h {1e3*h:5.1f} mm  E {E/1e3:6.1f} kPa  k {E/h/1e6:6.2f} MPa/m  ({r['layer']['modulus_basis']})")
    manifest["layer_map"] = dict(rule="per-segment k = E/h: h from this body's soft-tissue depth map, E the in vivo apparent modulus; pre-registered in docs/SEGMENT_CONTACT_SURFACES.md (97eba07)",
        source_card=str(CARD.relative_to(ROOT)), depth_map=str(DEPTH.relative_to(ROOT)),
        from_bundle=str(src.relative_to(ROOT)), checks=checks,
        replaces="skin_material's uniform layer (E 3 kPa, h 6.6 mm, k 1.72 MPa/m) as the plant's stiffness; skin_material is kept only as the record of what the bundle was built with")
    out.mkdir(parents=True); shutil.copytree(src / "meshes", out / "meshes")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"heel E {heel_pa/1e3:.1f} kPa (from {checks['heel_values_kpa']}), fat E {fat_pa/1e3:.1f} kPa; "
          f"Yang secant {checks['yang_secant_kpa']:.1f} kPa, {100*checks['yang_secant_vs_reported']:.1f}% from reported -> layer modulus")
    print(out)

if __name__ == "__main__": main()
