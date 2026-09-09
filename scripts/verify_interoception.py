"""Visceral afferent transduction fixtures and recorded-trajectory checks; no native process.

Four kinds of failure, all silent, all with a case here.

*A source key that does not resolve.* IHM keys proprioceptors by the bare muscle
id, and reading 'proprio:' + id returned 0.0 for every muscle every step -- a
stretch reflex that never fired looked like one with nothing to do. A visceral
channel has the same shape: an empty stomach and a broken key are both 0.0 Hz.
So every read raises, and `assert_sources` reports the whole missing set at once.

*A receptor law fitted to the corpus.* Thresholds and saturations are declared
from physiology, so a high-threshold channel is allowed to sit silent through a
protocol that never drives it. The test that a channel CAN be silent is as
important as the test that it can fire.

*A rate outside its declared range.* A saturating law that returns 46 Hz on a
45 Hz channel is a law that is not saturating.

*A trunk or fibre class that has drifted from IBM-1's declaration.* The two
tables are in two repositories.
"""
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.interoception import (  # noqa: E402
    CHANNELS, FIBRE_AVERSIVE_WEIGHT, SCALAR_NAMES, assert_sources,
    channel_table, derived_scalars, discomfort_index, endurance_hours,
    source_keys, visceral_afferent_rates)

CORPUS = ROOT / "data/derived/intero-corpus"
RECORDED = ROOT / "data/derived/systemic/exertion_v3"


def frame(path):
    with path.open() as fh:
        return json.loads(fh.readline())["values"]


class TransductionTests(unittest.TestCase):

    def test_a_missing_or_null_source_raises_and_never_returns_zero(self):
        v = {k: 1.0 for k in source_keys()}
        v.update({"liver_glycogen_g": 100.0, "muscle_glycogen_g": 400.0,
                  "metabolic_rate_w": 90.0, "stored_fat_g": 16000.0,
                  "fatigue_fraction": 0.0, "core_temperature_c": 37.0})
        visceral_afferent_rates(v)                    # the intact case works
        for k in list(source_keys())[:4]:
            broken = dict(v); broken.pop(k)
            with self.assertRaises(KeyError):
                visceral_afferent_rates(broken)
            nulled = dict(v); nulled[k] = None
            with self.assertRaises((KeyError, ValueError)):
                visceral_afferent_rates(nulled)
            nan = dict(v); nan[k] = float("nan")
            with self.assertRaises((KeyError, ValueError)):
                visceral_afferent_rates(nan)

    def test_assert_sources_names_every_missing_key_at_once(self):
        v = {k: 1.0 for k in list(source_keys())[3:]}
        with self.assertRaises(KeyError) as cm:
            assert_sources(v)
        self.assertIn("3 key(s) absent", str(cm.exception))

    def test_every_rate_stays_inside_its_declared_range(self):
        # sweep each channel's own quantity well past both ends of its declared
        # operating range.  a law that leaves the range is not saturating.
        for c in CHANNELS:
            for q in (-1e6, 0.0, c.threshold, c.saturation, 1e6):
                v = {k: q for k in source_keys()}
                r = c.rate_hz(v)
                self.assertTrue(math.isfinite(r))
                self.assertGreaterEqual(r, 0.0)
                self.assertLessEqual(r, c.max_hz + 1e-9)

    def test_a_high_threshold_channel_is_exactly_silent_below_threshold(self):
        # not "nearly zero": a floor that is nearly zero is indistinguishable
        # from a bug that returns nearly zero.  the splanchnic foregut
        # mechanoreceptor needs a stomach fuller than 700 mL and must read 0.0
        # below it, which is what makes its silence over a one-hour protocol a
        # correct reading rather than a dead wire.
        c = next(x for x in CHANNELS if x.name == "foregut_mechano")
        v = {k: 0.0 for k in source_keys()}
        v["stomach_water_ml"] = 600.0
        self.assertEqual(c.rate_hz(v), 0.0)
        v["stomach_water_ml"] = 2500.0
        self.assertEqual(c.rate_hz(v), c.max_hz)

    def test_the_inverted_channel_is_inverted(self):
        c = next(x for x in CHANNELS if x.name == "aortic_chemo_hypoxia")
        v = {k: 0.0 for k in source_keys()}
        v["arterial_o2_mmhg"] = 100.0
        self.assertEqual(c.rate_hz(v), 0.0)           # normal oxygen, silent
        v["arterial_o2_mmhg"] = 30.0
        self.assertEqual(c.rate_hz(v), c.max_hz)      # hypoxic, maximal

    def test_the_same_organ_reports_on_two_trunks_at_two_thresholds(self):
        # the pair the trunk/fibre-class split exists to express.  one gastric
        # volume, one vagal low-threshold channel, one splanchnic high-threshold
        # channel, and a range in which exactly one of them fires.
        vag = next(x for x in CHANNELS if x.name == "gastric_distension")
        spl = next(x for x in CHANNELS if x.name == "foregut_mechano")
        self.assertEqual(vag.trunk, "vagus")
        self.assertEqual(spl.trunk, "greater_splanchnic")
        self.assertEqual(vag.sources, spl.sources)
        v = {k: 0.0 for k in source_keys()}
        v["stomach_water_ml"] = 400.0
        self.assertGreater(vag.rate_hz(v), 1.0)
        self.assertEqual(spl.rate_hz(v), 0.0)


class ScalarTests(unittest.TestCase):

    def test_endurance_is_substrate_over_demand_and_refuses_a_zero_demand(self):
        v = {"liver_glycogen_g": 100.0, "muscle_glycogen_g": 500.0,
             "metabolic_rate_w": 100.0}
        # 600 g x 17 kJ/g = 10,200 kJ at 100 W -> 102,000 s -> 28.33 h
        self.assertAlmostEqual(endurance_hours(v), 10200e3 / 100.0 / 3600.0, 6)
        v["metabolic_rate_w"] = 200.0
        self.assertAlmostEqual(endurance_hours(v), 10200e3 / 200.0 / 3600.0, 6)
        v["metabolic_rate_w"] = 0.0
        with self.assertRaises(ValueError):
            endurance_hours(v)

    def test_discomfort_is_a_fraction_of_its_own_declared_ceiling(self):
        zero = {f"{c.trunk}/{c.name}": 0.0 for c in CHANNELS}
        full = {f"{c.trunk}/{c.name}": c.max_hz for c in CHANNELS}
        self.assertEqual(discomfort_index(zero), 0.0)
        self.assertAlmostEqual(discomfort_index(full), 1.0, 12)
        with self.assertRaises(KeyError):
            discomfort_index({k: 0.0 for k in list(zero)[1:]})

    def test_discomfort_weights_the_unmyelinated_classes_above_the_rest(self):
        # the fibre class is load-bearing: C fibres carry the great majority of
        # visceral nociceptive traffic and A-beta essentially none, so the same
        # firing rate on a C channel must contribute more than on an A-beta one.
        self.assertGreater(FIBRE_AVERSIVE_WEIGHT["c"],
                           FIBRE_AVERSIVE_WEIGHT["adelta"])
        self.assertGreater(FIBRE_AVERSIVE_WEIGHT["adelta"],
                           FIBRE_AVERSIVE_WEIGHT["abeta"])
        base = {f"{c.trunk}/{c.name}": 0.0 for c in CHANNELS}
        c_only = dict(base); c_only["greater_splanchnic/foregut_ischaemia"] = 10.0
        a_only = dict(base); a_only["vagus/aortic_baroreceptor"] = 10.0
        self.assertGreater(discomfort_index(c_only), discomfort_index(a_only))


class RecordedTrajectoryTests(unittest.TestCase):

    def test_every_recorded_protocol_transduces_without_a_default(self):
        if not RECORDED.exists():
            self.skipTest("no recorded systemic trajectories")
        for d in sorted(RECORDED.iterdir()):
            f = d / "frames.jsonl"
            if not f.exists():
                continue
            v = frame(f)
            assert_sources(v)
            r = visceral_afferent_rates(v)
            self.assertEqual(len(r), len(CHANNELS), d.name)
            s = derived_scalars(v, r)
            for k in SCALAR_NAMES:
                self.assertTrue(math.isfinite(s[k]), f"{d.name}: {k}")
            self.assertGreater(s["endurance_h"], 0.0)

    def test_the_corpus_channel_table_matches_this_module(self):
        meta_path = CORPUS / "meta.json"
        if not meta_path.exists():
            self.skipTest("no intero corpus built")
        meta = json.loads(meta_path.read_text())
        self.assertEqual(sorted(meta["channels"]),
                         sorted(f"{c.trunk}/{c.name}" for c in CHANNELS))
        self.assertEqual(meta["channel_table"], channel_table())
        # the corpus must record which split its sensation basis was fitted on,
        # because a basis fitted on data it is then evaluated against is a leak
        # that no downstream metric can see.
        self.assertIn("sensation_basis_fitted_on", meta)
        self.assertIn("held_out_protocols", meta)


if __name__ == "__main__":
    unittest.main()
