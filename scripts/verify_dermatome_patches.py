"""Coverage, identity and honesty checks on the dermatomal skin patches.

The point of these tests is not that the file parses.  It is that the numbers a
reader would quote out of it -- how much of the skin is covered, how many patches
are innervated, which nerve carries which patch -- cannot drift away from the
geometry and the peripheral model they claim to be about.
"""
from pathlib import Path
import json
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "data/derived/canonical"


class DermatomeVerification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = OUT / "dermatomes.json"
        if not path.exists():
            raise unittest.SkipTest("run scripts/build_dermatome_patches.py first")
        cls.d = json.loads(path.read_text())
        cls.per = json.loads((OUT / "peripheral.json").read_text())
        cls.anat = json.loads((OUT / "anatomy.json").read_text())

    def test_area_is_the_measured_exterior_and_not_the_raw_mesh(self):
        """The denominator is stated and is the one the skin entity declares.

        The raw mesh is 3.50 m2 and includes interior and orifice surfaces; the
        exterior component is 1.78 m2.  Quoting coverage against the raw area
        would understate it by nearly half, and quoting the raw area as body
        surface area would overstate the skin by the same factor.
        """
        c = self.d["coverage"]
        skin = next(e for e in self.anat["entities"] if e["id"] == c["skin_entity_id"])
        self.assertAlmostEqual(c["exterior_component_area_m2"],
                               skin["physical_surface_support"]["area_m2"], places=9)
        self.assertAlmostEqual(c["raw_source_area_m2"], skin["surface_area_m2"], places=9)
        self.assertGreater(c["raw_source_area_m2"], c["exterior_component_area_m2"])

    def test_patches_tile_the_exterior_without_gap_or_overlap(self):
        c = self.d["coverage"]
        total = sum(p["area_m2"] for p in self.d["patches"])
        self.assertAlmostEqual(total, c["exterior_component_area_m2"], places=9)
        self.assertEqual(c["unassigned_triangle_count"], 0)
        self.assertEqual(sum(p["triangle_count"] for p in self.d["patches"]),
                         c["exterior_triangle_count"])

    def test_every_patch_names_a_nerve_that_the_body_actually_declares(self):
        """A patch routed to a nerve id that does not exist is not innervated.

        This is the check the programme keeps needing: the trunk table and the
        body were separate objects, so a name could look wired and reach nothing.
        """
        nerves = {n["id"] for n in self.per["nerves"]}
        relays = {r["id"] for r in self.per["relays"]}
        for p in self.d["patches"]:
            self.assertIn(p["nerve_id"], nerves, p["id"])
            self.assertIn(p["relay_id"], relays, p["id"])
        self.assertEqual(self.d["coverage"]["trunks_absent_from_peripheral_json"], [])

    def test_root_levels_are_spinal_or_explicitly_absent(self):
        """A trigeminal patch must never carry a spinal root level."""
        levels = {f"{a}{i}" for a, n in (("c", 8), ("t", 12), ("l", 5), ("s", 5))
                  for i in range(1, n + 1)}
        cranial = 0
        for p in self.d["patches"]:
            if p["root_level"] is None:
                self.assertTrue(p["root_gap_reason"], p["id"])
                self.assertEqual(p["nerve_name"], "trigeminal", p["id"])
                cranial += 1
            else:
                self.assertIn(p["root_level"], levels, p["id"])
                self.assertEqual(p["root_level"], p["dermatome"], p["id"])
        self.assertEqual(cranial, self.d["coverage"]["patches_without_spinal_root"])
        self.assertGreater(cranial, 0, "the face cannot have vanished")

    def test_patch_positions_lie_on_the_real_skin_surface(self):
        import gzip
        skin = next(e for e in self.anat["entities"]
                    if e["id"] == self.d["skin_geometry"]["entity_id"])
        g = json.loads(gzip.decompress(
            (ROOT / skin["reference_geometry"]["path"]).read_bytes()))
        V = np.asarray(g["positions"], float).reshape(-1, 3)
        P = np.asarray([p["position_m"] for p in self.d["patches"]], float)
        from scipy.spatial import cKDTree
        self.assertLess(float(cKDTree(V).query(P)[0].max()), 1e-9)

    def test_path_lengths_are_positive_and_reach_their_own_relay(self):
        relays = {r["id"]: np.asarray(r["position_m"], float)
                  for r in self.per["relays"]}
        for p in self.d["patches"]:
            straight = float(np.linalg.norm(
                np.asarray(p["position_m"]) - relays[p["relay_id"]]))
            self.assertGreater(p["path_length_m"], 0.0)
            # a routed path through a proximal waypoint is never shorter than
            # the straight line it replaces
            self.assertGreaterEqual(p["path_length_m"] + 1e-12, straight, p["id"])

    def test_the_measured_levels_really_are_measured(self):
        """T2-T12 on the trunk wall must cite the ribs, not a coordinate band."""
        wall = [p for p in self.d["patches"] if p["region"] == "thorax_wall"]
        self.assertTrue(wall)
        for p in wall:
            self.assertIn("nearest rib", p["assignment_rule"], p["id"])
        levels = sorted({p["dermatome"] for p in wall}, key=lambda s: int(s[1:]))
        self.assertEqual(levels, [f"t{i}" for i in range(2, 13)],
                         "the rib cage should produce a complete T2-T12 ladder")

    def test_the_evidence_flags_stay_honest(self):
        self.assertFalse(self.d["measured_dermatome_atlas"])
        self.assertFalse(self.d["biological_validation"])
        for p in self.d["patches"]:
            self.assertFalse(p["measured_dermatome_atlas"])
            self.assertTrue(p["assignment_rule"])

    def test_sides_are_balanced_within_a_few_percent(self):
        """A gross left/right asymmetry means the side rule or a frame is wrong."""
        a = sum(p["area_m2"] for p in self.d["patches"] if p["side"] == "left")
        b = sum(p["area_m2"] for p in self.d["patches"] if p["side"] == "right")
        self.assertLess(abs(a - b) / (a + b), 0.02, f"left {a:.4f} right {b:.4f}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
