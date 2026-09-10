"""does the global atlas->scaffold map put the atlas's JOINTS where the scaffold's are?

a straight ligament across a joint sweeps about the joint centre.  if the global
similarity (24.7 mm RMS at bone centroids) lands the atlas knee ~2 cm from the
scaffold knee, a cruciate attached to atlas geometry sweeps about the wrong centre
and tears inside the joint's range no matter where in its footprint it attaches
or what it wraps around.  that would make per-segment registration -- the same
fix the skin needs -- the ligament fix too.

atlas joint centre: midpoint of the closest-approach pairs between the two bone
groups that meet there (within 12 mm), mapped to ground by the binding's own
inverse similarity at its own reference pose.  scaffold joint centre: the child
body's origin at that pose, which in this model sits on the joint.
"""
import importlib.util, json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
ROOT = Path(__file__).resolve().parents[1]
def load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
bind = load("bind_anatomy", "scripts/bind_anatomy_to_segments.py")
render = load("render_body_3d", "scripts/render_body_3d.py")
entities = json.loads((ROOT / "data/derived/canonical/anatomy.json").read_text())["entities"]
bind.GEOMETRY_PATH.update({e["id"]: ROOT / e["reference_geometry"]["path"] for e in entities if e.get("reference_geometry")})
groups = {}
for e in entities:
    if e["role"] == "rigid_bone": groups.setdefault(bind.named_segment(e["name"]), []).append(e["id"])
binding = json.loads((ROOT / "data/derived/anatomy-segment-binding/binding.json").read_text())
inv = np.linalg.inv(np.asarray(binding["similarity_atlas_from_opensim_ground"], float))
model = render.OsimModel(ROOT / "data/models/engineering_stance_v1/model.osim")
rest = model.forward(binding["reference_pose_rad"])
cloud = {s: np.concatenate([bind.read_geometry(i) for i in ids]) @ inv[:3, :3].T + inv[:3, 3] for s, ids in groups.items()}
res = binding["registration"]["segment_centroid_residual_m"]

# THE CONTROL.  "scaffold joint centre = child body origin" is an assumption, and
# the atlas centre is taken at the articular surface, so part of any offset could
# be the definition rather than the registration.  apply the IDENTICAL rule to the
# scaffold's OWN bone meshes: if its own articular centre lands near its own child
# origin, the definition is sound and the atlas offset is registration error.
spec = importlib.util.spec_from_file_location("bscm", ROOT / "scripts/build_skin_contact_meshes.py")
bscm = importlib.util.module_from_spec(spec); spec.loader.exec_module(bscm)
osim = {s: p @ rest[s][:3, :3].T + rest[s][:3, 3] for s, p in bscm.bone_clouds().items() if s in rest}

def centre(pa, pc, margin=0.008):
    """midpoint of every closest-approach pair within `margin` of the minimum gap."""
    d, j = cKDTree(pc).query(pa)
    ok = d <= d.min() + margin
    return ((pa[ok] + pc[j[ok]]) / 2).mean(0), float(d.min()), int(ok.sum())

print(f"{'joint':20s} {'scaffold own bones':>19s} {'atlas bones':>12s}    atlas gap  centroid residuals")
print(f"{'':20s} {'-> child origin':>19s} {'-> origin':>12s}")
for parent, child in [("femur_l", "tibia_l"), ("femur_r", "tibia_r"), ("tibia_l", "talus_l"), ("tibia_r", "talus_r"),
                      ("humerus_l", "ulna_l"), ("humerus_r", "ulna_r"), ("pelvis", "femur_l"), ("pelvis", "femur_r")]:
    origin = rest[child][:3, 3]
    so, sgap, _ = centre(osim[parent], osim[child])
    ao, agap, _ = centre(cloud[parent], cloud[child])
    print(f"{parent+'-'+child:20s} {1000*np.linalg.norm(so-origin):16.1f} mm {1000*np.linalg.norm(ao-origin):9.1f} mm"
          f"   {1000*agap:6.1f} mm   {1000*res.get(parent, float('nan')):.1f}, {1000*res.get(child, float('nan')):.1f} mm")
print(f"\nfor scale: cruciate slack lengths are 33-37 mm; the similarity's own RMS residual is "
      f"{1000*binding['registration']['segment_centroid_residual_rms_m']:.1f} mm")
