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

if __name__ == "__main__": main()
