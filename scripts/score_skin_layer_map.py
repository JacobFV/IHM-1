"""Score the per-segment skin layer against the gates pre-registered in 97eba07.

Reads a measure_segment_contact_meshes.py report.  Gates (docs/SEGMENT_CONTACT_SURFACES.md,
"Depth alone would make it worse"), fixed before the engine carried per-segment stiffness:

  1. momentum balance: worst residual over the run <= 1e-5 N;
  2. never bone: every contacting patch's worst compression < that patch's h;
  3. heel strain: calcn worst compression / h <= 0.73 (top of Teng 2022's in vivo gait range).

The declared-layer `skin` arm, if present in the same report, is scored against its own uniform
h beside it -- reported, not gated.
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESIDUAL_N, HEEL_STRAIN = 1e-5, 0.73

def thicknesses(arm):
    manifest = json.loads((ROOT / arm["bundle"] / "manifest.json").read_text())
    if "layer_map" in manifest:
        return {r["element"]: r["layer"]["thickness_m"] for r in manifest["records"]}
    h = manifest["skin_material"]["layer_thickness_m"]
    return {r["element"]: h for r in manifest["records"]}

def score(arm):
    if "failed" in arm: return dict(arm=arm["arm"], failed=arm["failed"])
    h = thicknesses(arm); comp = arm["max_compression_m"]
    touching = {e: c for e, c in comp.items() if c > 0}
    strain = {e: c / h[e] for e, c in touching.items()}
    heel = {e: s for e, s in strain.items() if e.startswith("skin_calcn")}
    g1 = arm["max_momentum_balance_residual_norm_n"] <= RESIDUAL_N
    g2 = all(s < 1 for s in strain.values())
    g3 = all(s <= HEEL_STRAIN for s in heel.values())
    return dict(arm=arm["arm"], residual_n=arm["max_momentum_balance_residual_norm_n"], gate_momentum=g1,
                strain_by_patch={e: round(s, 4) for e, s in sorted(strain.items(), key=lambda kv: -kv[1])},
                gate_never_bone=g2, heel_strain=heel, gate_heel=g3,
                vertical_force_final_n=arm["final"]["vertical_contact_force_n"], weight_n=arm["weight_n"],
                pelvis_ty_final_m=arm["final"]["pelvis_ty_m"])

def main(path):
    report = json.loads(Path(path).read_text()); out = {}
    for arm in report["arms"]:
        if arm["arm"] not in ("skin_layer_map", "skin"): continue
        s = score(arm); out[arm["arm"]] = s
        tag = "GATED" if arm["arm"] == "skin_layer_map" else "reported (declared uniform layer)"
        print(f"== {arm['arm']}  [{tag}]")
        if "failed" in s: print("   FAILED TO RUN:", s["failed"]); continue
        print(f"   1 momentum: worst residual {s['residual_n']:.3e} N (<= {RESIDUAL_N}) -> {'PASS' if s['gate_momentum'] else 'FAIL'}")
        print(f"   2 never bone: {len(s['strain_by_patch'])} patches touched; worst compression/h "
              + ", ".join(f"{e} {v:.2f}" for e, v in list(s['strain_by_patch'].items())[:6])
              + f" -> {'PASS' if s['gate_never_bone'] else 'FAIL'}")
        print(f"   3 heel strain: {s['heel_strain']} (<= {HEEL_STRAIN}) -> {'PASS' if s['gate_heel'] else 'FAIL'}")
        print(f"   reported: final vertical contact {s['vertical_force_final_n']:.1f} N vs weight {s['weight_n']:.1f} N; pelvis_ty {s['pelvis_ty_final_m']:.4f} m")
    lm = out.get("skin_layer_map")
    verdict = None if lm is None or "failed" in lm else all((lm["gate_momentum"], lm["gate_never_bone"], lm["gate_heel"]))
    print(f"\nVERDICT (skin_layer_map, all three gates): {'PASS' if verdict else 'FAIL' if verdict is False else 'NOT SCORED'}")
    Path(path).with_suffix(".verdict.json").write_text(json.dumps(dict(report=str(path), verdict=verdict, arms=out), indent=2) + "\n")

if __name__ == "__main__": main(sys.argv[1])
