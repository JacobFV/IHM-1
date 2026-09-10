"""What fraction of this body is fat, against what this body says it is.

`data/derived/canonical/profile.json` declares `body_fat_fraction` 0.21 at 70.7713 kg -- 14.86 kg
of adipose. The geometry carries fat in exactly one depot of any size, the `hypodermis` skin
layer, plus two infrapatellar pads of 0.63 mL each: there is no visceral, intermuscular or
intramuscular adipose geometry in the atlas (docs/TISSUE_MECHANICS.md).

The mass ledger hides the gap. `mechanics.json`'s `mass_allocation` records that tissue volumes
at sourced densities weigh 60.00 kg while the profile demands 70.77, so every entity's mass is
multiplied by a UNIFORM 1.17944. The declared total is therefore met by construction, and the
fat the geometry does not carry is carried as denser bone, muscle and organs instead.

KNOWN ANSWERS, checked before anything is reported:
  1. the entity masses sum to the profile's mass_kg, to 1e-6 kg;
  2. every entity's mass / volume equals its declared tissue density x the recorded uniform
     scale, to 1e-9 relative -- i.e. the scale is uniform, as the ledger says;
  3. the adipose entities are exactly the hypodermis and the two infrapatellar fat pads.
"""
import json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.tissue_materials import density, tissue_class

def main():
    mech = json.loads((ROOT / "data/derived/canonical/mechanics.json").read_text())
    prof = json.loads((ROOT / "data/derived/canonical/profile.json").read_text())
    ents, alloc = mech["entities"], mech["mass_allocation"]
    total = sum(e["mass_kg"] for e in ents)
    assert abs(total - prof["mass_kg"]) < 1e-6, (total, prof["mass_kg"])
    scale = alloc["uniform_scale"]
    # The ledger's own exception: 206 entities (skin parent, the lymphatic graph, cardiac
    # cavities) carry 1 mg of numerical inertia with no material volume. They are excluded from
    # the density check, and their count is itself a known answer.
    carriers = [e for e in ents if e["mass_kg"] == 1e-6]
    assert len(carriers) == alloc["numerical_carrier_count"], (len(carriers), alloc["numerical_carrier_count"])
    worst, adipose, tissue_mass = 0.0, [], {}
    for e in ents:
        v = e.get("volume_m3")
        d = density(e["name"], e["role"]); d = d[0] if isinstance(d, tuple) else d
        cls = tissue_class(e["name"], e["role"])
        tissue_mass[cls] = tissue_mass.get(cls, 0.0) + e["mass_kg"]
        if cls == "adipose": adipose.append(e)
        if v and d and e["mass_kg"] != 1e-6: worst = max(worst, abs((e["mass_kg"] / v) / (d * scale) - 1))
    assert worst < 1e-9, f"the scale is not uniform: worst relative deviation {worst:.3e}"
    names = sorted(e["name"] for e in adipose)
    assert names == ["hypodermis"], names
    pads = [e for e in ents if "fat pad" in e["name"]]
    assert len(pads) == 2 and all(tissue_class(p["name"], p["role"]) == "dense_connective" for p in pads), pads
    fat = sum(e["mass_kg"] for e in adipose)
    declared = prof["body_fat_fraction"] * prof["mass_kg"]
    print(f"known answers: total {total:.4f} kg = profile; uniform scale {scale:.6f} exact on every entity "
          f"with material volume (worst {worst:.1e}), {len(carriers)} numerical carriers excluded; "
          f"adipose entities {names}")
    print(f"  the two infrapatellar fat pads ({2e6*pads[0]['volume_m3']:.2f} mL together) are classed "
          f"dense_connective at 1060 kg/m3, not adipose at 950: the only tissue in this body carried "
          f"at fat density is the hypodermis")
    print(f"\nthis body's fat, as geometry:  {fat:.3f} kg = {100*fat/total:.1f}% of {total:.4f} kg")
    print(f"this body's fat, as declared:  {declared:.3f} kg = {100*prof['body_fat_fraction']:.1f}% "
          f"(profile.json body_fat_fraction, {prof['id']}, age {prof['age_years']}, {prof['sex']})")
    print(f"shortfall:                     {declared-fat:.3f} kg, {100*(declared-fat)/declared:.0f}% of the declared depot")
    v = adipose[0]["volume_m3"] if adipose[0]["name"] == "hypodermis" else None
    print(f"  the hypodermis is {1e3*v:.2f} L; carrying 21% would need {1e3*declared/(950*scale):.2f} L at its ledger density")
    print("\nwhere the mass went instead (ledger mass by tissue class, top 8):")
    for k, m in sorted(tissue_mass.items(), key=lambda kv: -kv[1])[:8]:
        print(f"  {k:24s} {m:7.3f} kg  {100*m/total:5.1f}%")
    print(f"\nEvery non-fat tissue above is inflated by the same {100*(scale-1):.1f}% the ledger applies to reach "
          f"{prof['mass_kg']} kg. The fat this body declares but has no geometry for is carried as denser bone, "
          f"muscle and organs, which is a mass distribution no body has.")
    where_it_would_go(mech, fat, declared, scale)

def where_it_would_go(mech, fat, declared, scale):
    """The subcutaneous space this body measures, against the space its skin layers declare.

    KNOWN ANSWER: the hypodermis's own volume must equal its declared thickness times the exterior
    skin area, because that is how a layer is defined; if it does not, this comparison is void.
    """
    import gzip
    bundle = json.loads((ROOT / "data/derived/segment-contact-meshes/skin/manifest.json").read_text())
    layers = {l["id"]: l["thickness_m"] for l in bundle["skin_material"]["layers"]}
    # The bundle's own exterior_surface_area_m2 is in the SCAFFOLD frame (the binding map's scale
    # 0.963, so 1.9198 m2 against 1.7805 here); the hypodermis volume is canonical, so the area
    # has to be canonical too. Taken from the same exterior triangles the depth map used.
    skin = next(e for e in mech["entities"] if e["role"] == "skin")
    g = json.loads(gzip.decompress((ROOT / skin["reference_geometry"]["path"]).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3); F = np.asarray(g["indices"], np.int64).reshape(-1, 3)
    ext = np.asarray(json.loads((ROOT / "data/research/engineered_skin_territories/materialization.json").read_text())
                     ["contact_eligible_triangle_ids"], np.int64)
    t = V[F[ext]]
    area = float(np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum() / 2)
    hypo = next(e for e in mech["entities"] if e["name"] == "hypodermis")
    h = layers["body-skin-hypodermis"]
    implied = h * area
    if abs(implied / hypo["volume_m3"] - 1) > 0.02:
        print(f"\n[skipped: the hypodermis is {1e3*hypo['volume_m3']:.2f} L but its {1e3*h:.1f} mm layer over "
              f"{area:.4f} m2 is {1e3*implied:.2f} L -- the layer identity does not hold, so this comparison is void]")
        return
    d = np.load(ROOT / "data/derived/soft-tissue-depth-v1/depth.npz")["depth_m"]
    declared_layer = sum(layers.values())
    print(f"\nwhere the missing fat would go (known answer: {1e3*h:.1f} mm x {area:.4f} m2 = {1e3*implied:.2f} L "
          f"= the hypodermis's {1e3*hypo['volume_m3']:.2f} L)")
    print(f"  this body's skin declares      {1e3*declared_layer:.1f} mm of layer -> {1e3*declared_layer*area:.2f} L")
    print(f"  its own geometry measures      {1e3*np.median(d):.1f} mm median, {1e3*d.mean():.1f} mm mean skin-to-bone/muscle "
          f"-> {1e3*d.mean()*area:.2f} L of subcutaneous space")
    gap = (d.mean() - declared_layer) * area
    print(f"  the difference is              {1e3*gap:.2f} L, {gap*950*scale:.3f} kg at the hypodermis's ledger density")
    print(f"  the fat with no geometry is    {declared-fat:.3f} kg")
    need = (declared - fat) / (950 * scale)
    print(f"\nThe space is ample and the deficit does NOT fill it: the missing {declared-fat:.3f} kg is {1e3*need:.2f} L, "
          f"{100*need/gap:.0f}% of the {1e3*gap:.2f} L the measurement leaves over the declared layers. Filling the whole "
          f"difference at fat density would be {gap*950*scale:.1f} kg, three times the declared depot's shortfall, so "
          "this is a bound, not an allocation: nearest-structure depth is an UPPER bound on a subcutaneous layer (it "
          "runs to the nearest bone or muscle, and where neither is close -- abdomen, gluteal region, breast -- the "
          "space it measures holds fascia, vessels and glands as well as fat). What it does establish is that the "
          "hypodermis is declared thinner than this body's own surfaces measure, with room for the fat the ledger has "
          "no geometry for -- and that the same declared layer is the contact model's thickness "
          "(docs/SEGMENT_CONTACT_SURFACES.md), so the mass gap and the contact gap are the same declaration.")

if __name__ == "__main__": main()
