"""Visceral afferent firing rates from the native systemic state.

The body already simulates digestion, absorption, substrate stores and exertion.
Nothing read that state as *afference*: `ihm/assembly/systemic.py` records 180
native quantities per frame and hands them to an audit, and the brain on the
other side of the wire has vision, hearing and a somatic port and is blind to its
own gut.  This module is the transducer between the two.

**Every rate here comes from a native quantity that BioGears integrated.**  There
is no oscillator, no noise process and no hand-drawn trajectory: the input is
`NativeSession.snapshot()['values']` live, or a recorded `frames.jsonl` row
offline, and those are the same dict.  What this module adds is the receptor law
-- a declared operating range and a saturating transfer -- which is the one thing
the native engine does not compute for the gut.

**What is already native and is NOT re-derived here.**  `native_afferents.py`
carries four channels the engine's own nervous model computes as Hz --
`chemoreceptor`, `baroreceptor_carotid`, `baroreceptor_aortic`,
`pulmonary_stretch`.  Those are cardiorespiratory.  There is no native gastric,
intestinal, hepatoportal, renal or bladder afferent, and no splanchnic route of
any kind, which is exactly the gap this fills.  A caller that has the afferent
adapter built should prefer the native rate for the two channels that overlap;
the overlap is declared per channel in `NATIVE_EQUIVALENT`.

**The trunk and the fibre class are load-bearing, not labels.**  The same organ
reports on two routes at two speeds: the stomach's low-threshold mechanoreceptors
are vagal and myelinated, and the high-threshold ones that signal painful
distension are splanchnic and unmyelinated.  IBM-1's `ibm/topologies/nerve.py`
resolves conduction delay per fibre class, and over IHM's own authored route
lengths that is 9 ms for the vagal A-beta channel against 508 ms for the vagal C
channel -- a half-second separation between the same meal's mechanical and
chemical report.  Declaring the trunk and the class per channel is what lets the
receiving model reproduce that instead of lumping it.

**The key trap this module is written against.**  IHM keys proprioceptor rates by
the BARE muscle id, and reading `"proprio:" + id` returned 0.0 for every muscle
every step -- a silent zero, not an exception, so a stretch reflex that never
fired looked like a stretch reflex with nothing to do.  Every source key here is
therefore checked before it is read: `visceral_afferent_rates` raises on a
missing key, on a null and on a non-finite value, and `assert_sources` will tell
a caller up front which of the 15 channels a given snapshot can support.  A
channel that cannot be sourced is an error, never a zero.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# receptor laws
# ---------------------------------------------------------------------------

def _saturating(x: float, threshold: float, saturation: float, max_hz: float,
                invert: bool = False) -> float:
    """A receptor with a threshold, a saturation and a smooth knee between them.

    The form is the one IBM-1's ontology already declares for the baroreceptor
    (`ibm/processes/transduction.py:baroreceptor_arterial`): silent below
    threshold, saturated above, steepest in between, with the operating range
    placed where the quantity normally lives.  A raised cosine rather than a
    logistic so that "silent" and "saturated" are exact rather than asymptotic --
    a channel that is supposed to be quiet on an empty stomach should read 0.0
    and not 0.7 Hz, because a floor that is nearly zero is indistinguishable from
    a bug that returns nearly zero.
    """
    if saturation <= threshold:
        raise ValueError("saturation must exceed threshold")
    u = (x - threshold) / (saturation - threshold)
    u = 0.0 if u < 0.0 else (1.0 if u > 1.0 else u)
    if invert:
        u = 1.0 - u
    return max_hz * 0.5 * (1.0 - math.cos(math.pi * u))


@dataclass(frozen=True)
class VisceralChannel:
    """One afferent channel: where it originates, what it reads, how it fires.

    `sources` are literal keys of the native systemic snapshot.  `combine` maps
    them to the physical quantity the receptor sits in, in `units`; the law then
    maps that to Hz.  Keeping the two apart is deliberate -- the quantity is
    physiology and is checkable against the engine, the law is a receptor model
    and is not.
    """
    name: str
    trunk: str                       # an IBM-1 declared trunk, and an IHM route
    fibre: str                       # an IBM-1 declared fibre class on that trunk
    units: str
    sources: tuple[str, ...]
    threshold: float
    saturation: float
    max_hz: float
    doc: str
    invert: bool = False
    combine: str = "sum"             # sum | first
    weights: tuple[float, ...] = ()
    native_equivalent: str = ""      # a channel `native_afferents.py` also computes

    def quantity(self, values: dict) -> float:
        w = self.weights or (1.0,) * len(self.sources)
        if len(w) != len(self.sources):
            raise ValueError(f"{self.name}: {len(w)} weights for "
                             f"{len(self.sources)} sources")
        xs = [_read(values, k, self.name) for k in self.sources]
        if self.combine == "first":
            return xs[0]
        return math.fsum(a * b for a, b in zip(xs, w))

    def rate_hz(self, values: dict) -> float:
        return _saturating(self.quantity(values), self.threshold,
                           self.saturation, self.max_hz, self.invert)


def _read(values: dict, key: str, channel: str) -> float:
    """Read one native key, or raise.  NEVER a default.

    The whole reason this function exists rather than `values.get(key, 0.0)`:
    a wrong key that returns 0.0 produces a channel that is silent, and a silent
    interoceptive channel is exactly what a body with an empty stomach looks
    like.  The failure and the correct answer are the same number.
    """
    if key not in values:
        raise KeyError(f"channel {channel!r}: native snapshot has no {key!r}")
    v = values[key]
    if v is None:
        raise ValueError(f"channel {channel!r}: native {key!r} is null "
                         f"(absent or nonfinite native quantity)")
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(f"channel {channel!r}: native {key!r} is not a finite "
                         f"number ({v!r})")
    return float(v)


# ---------------------------------------------------------------------------
# the channels
# ---------------------------------------------------------------------------
# Operating ranges are declared from physiology and NOT fitted to the recorded
# corpus, so that a channel is allowed to sit silent through a whole run if the
# body never entered its range.  Two of them do: the splanchnic foregut
# mechanoreceptor needs a stomach fuller than any protocol here produces, and
# that silence is a correct reading of a high-threshold nociceptor, not a dead
# wire.  Fitting the ranges to the data would have hidden that.

CHANNELS: tuple[VisceralChannel, ...] = (
    # -- vagus: ~80% afferent, the main gut-to-brain route ------------------
    VisceralChannel(
        "gastric_distension", "vagus", "abeta", "mL",
        ("stomach_water_ml", "stomach_carbohydrate_g", "stomach_protein_g",
         "stomach_fat_g"),
        threshold=50.0, saturation=1200.0, max_hz=45.0,
        weights=(1.0, 1.0, 1.0, 1.11),
        doc="low-threshold vagal mechanoreceptors in the gastric wall.  the "
            "solid weights convert grams to millilitres at the macronutrient "
            "densities (carbohydrate and protein ~1.0 g/mL, fat ~0.9), so the "
            "quantity is a volume and not a mixed sum.  threshold at 50 mL "
            "because a stomach at resting tone is not silent-empty"),
    VisceralChannel(
        "gastric_nutrient", "vagus", "c", "g",
        ("stomach_carbohydrate_g", "stomach_protein_g", "stomach_fat_g"),
        threshold=0.5, saturation=90.0, max_hz=20.0,
        doc="vagal chemosensitive endings reporting macronutrient load in the "
            "stomach.  unmyelinated, so this is the half-second-late chemical "
            "report of the same meal the A-beta channel reported in 9 ms"),
    VisceralChannel(
        "intestinal_distension", "vagus", "abeta", "mL",
        ("SmallIntestineChyme.volume_ml",),
        threshold=20.0, saturation=600.0, max_hz=40.0,
        doc="vagal mechanoreceptors on the small intestine, reading the chyme "
            "volume the native lumped digestion model actually carries"),
    VisceralChannel(
        "intestinal_nutrient", "vagus", "c", "mg/dL",
        ("SmallIntestineChyme.Glucose.concentration_mg_per_dl",),
        threshold=10.0, saturation=5000.0, max_hz=25.0,
        doc="vagal glucose sensing in the intestinal wall.  the corpus reaches "
            "14,700 mg/dL post-meal, three times the saturation, so this "
            "channel is pinned for much of a digesting hour -- which is the "
            "correct behaviour of a saturating chemoreceptor and the reason "
            "saturation is declared rather than fitted"),
    VisceralChannel(
        "gi_absorption", "vagus", "adelta", "mL/min",
        ("chyme_water_absorption_ml_per_min",),
        threshold=0.3, saturation=14.0, max_hz=18.0,
        doc="absorptive activity at the intestinal wall.  A-delta puts it "
            "between the fast mechanical and the slow chemical report, which is "
            "the only reason the trunk needs three classes rather than two"),
    VisceralChannel(
        "hepatoportal_glucose", "vagus", "c", "mg/dL",
        ("LiverVasculature.Glucose.concentration_mg_per_dl",),
        threshold=40.0, saturation=160.0, max_hz=20.0,
        doc="the hepatoportal glucose sensor: vagal afferents in the portal "
            "wall that report absorbed glucose before it reaches the systemic "
            "circulation.  it is the one channel here that sees a meal as "
            "energy rather than as bulk"),
    VisceralChannel(
        "pulmonary_stretch", "vagus", "abeta", "mL",
        ("lung_volume_ml",),
        threshold=1800.0, saturation=4200.0, max_hz=50.0,
        doc="slowly adapting pulmonary stretch receptors, myelinated and vagal.",
        native_equivalent="pulmonary_stretch"),
    VisceralChannel(
        "aortic_baroreceptor", "vagus", "abeta", "mmHg",
        ("mean_arterial_pressure_mmhg",),
        threshold=60.0, saturation=180.0, max_hz=50.0,
        doc="aortic arch baroreceptors.  threshold and saturation are taken "
            "from IBM-1's declared `baroreceptor_arterial` parameters (60 and "
            "180 mmHg) rather than chosen here, so the two models state the "
            "same operating range",
        native_equivalent="baroreceptor_aortic"),
    VisceralChannel(
        "aortic_chemo_hypoxia", "vagus", "c", "mmHg",
        ("arterial_o2_mmhg",),
        threshold=40.0, saturation=95.0, max_hz=30.0, invert=True,
        doc="the hypoxic limb of aortic body chemoreception: silent when "
            "arterial oxygen is normal, maximal when it is not.  inverted, "
            "which is why `invert` is a declared field and not a minus sign "
            "buried in a lambda"),
    VisceralChannel(
        "aortic_chemo_hypercapnia", "vagus", "c", "mmHg",
        ("arterial_co2_mmhg",),
        threshold=35.0, saturation=60.0, max_hz=30.0,
        doc="the hypercapnic limb of the same body"),

    # -- splanchnic: the high-threshold arm, unmyelinated ------------------
    # AUTONOMIC trunks declare ("b_preganglionic", "c_postganglionic", "c") and
    # only "c" of those is afferent, so every splanchnic channel is a C fibre.
    # That is not a simplification -- it is what the trunk composition says.
    VisceralChannel(
        "foregut_mechano", "greater_splanchnic", "c", "mL",
        ("stomach_water_ml", "stomach_carbohydrate_g", "stomach_protein_g",
         "stomach_fat_g"),
        threshold=700.0, saturation=2000.0, max_hz=20.0,
        weights=(1.0, 1.0, 1.0, 1.11),
        doc="the SAME gastric volume as `gastric_distension`, on a different "
            "trunk with a different threshold.  this is the pair the whole "
            "trunk/class split exists to express: the vagal channel is a "
            "low-threshold volume report arriving in 9 ms, this one is a "
            "high-threshold nociceptive report arriving in 168 ms and staying "
            "silent unless the stomach is genuinely overfull"),
    VisceralChannel(
        "foregut_ischaemia", "greater_splanchnic", "c", "mg/dL",
        ("Aorta.Lactate.concentration_mg_per_dl",),
        threshold=20.0, saturation=150.0, max_hz=25.0,
        doc="metabolite accumulation as visceral nociceptive drive.  arterial "
            "lactate is the mismatch indicator -- glycolytic flux exceeding "
            "what oxidation can consume -- which is the quantity that makes an "
            "ischaemic or heavily exercising viscus hurt"),
    VisceralChannel(
        "midgut_distension", "lesser_splanchnic", "c", "mL",
        ("SmallIntestineChyme.volume_ml",),
        threshold=500.0, saturation=1200.0, max_hz=20.0,
        doc="the high-threshold twin of `intestinal_distension`, on the trunk "
            "IHM routes to the upper abdominal visceral endpoint"),
    VisceralChannel(
        "renal_afferent", "least_splanchnic", "c", "mL/min",
        ("glomerular_filtration_ml_per_min",),
        threshold=60.0, saturation=220.0, max_hz=15.0,
        doc="renal afferents on the trunk IHM routes to the renal visceral "
            "endpoint, reading filtration rate"),
    VisceralChannel(
        "bladder_distension", "pelvic_splanchnic", "c", "mL",
        ("Bladder.volume_ml",),
        threshold=150.0, saturation=500.0, max_hz=25.0,
        doc="bladder filling on the pelvic splanchnic route, which is the "
            "correct afferent nerve for it.  threshold at first sensation "
            "(~150 mL) and saturation near urgency (~500 mL)"),
)

#: channels whose signal the native afferent adapter also computes, when it is
#: built.  a caller with `native_biogears_afferent_signed` available should
#: prefer the native rate for these two and keep the other thirteen from here.
NATIVE_EQUIVALENT = {c.name: c.native_equivalent for c in CHANNELS
                     if c.native_equivalent}

CHANNEL_NAMES = tuple(c.name for c in CHANNELS)
TRUNKS = tuple(dict.fromkeys(c.trunk for c in CHANNELS))


def source_keys() -> tuple[str, ...]:
    """Every native key the channel set reads, deduplicated."""
    seen: dict[str, None] = {}
    for c in CHANNELS:
        for k in c.sources:
            seen[k] = None
    return tuple(seen)


def assert_sources(values: dict) -> None:
    """Raise unless every channel can be sourced from this snapshot.

    Call this ONCE before a run rather than discovering a bad key on frame
    40,000, and call it at all rather than trusting that a `.get` returned a
    real zero.  The message names every missing key at once, because fixing
    them one exception at a time is how a morning goes.
    """
    missing = [k for k in source_keys() if k not in values]
    null = [k for k in source_keys() if k in values and values[k] is None]
    if missing or null:
        raise KeyError(
            f"visceral afference cannot be sourced: "
            f"{len(missing)} key(s) absent {sorted(missing)}; "
            f"{len(null)} key(s) null {sorted(null)}. "
            f"the snapshot has {len(values)} keys.")


def visceral_afferent_rates(values: dict) -> dict[str, float]:
    """Afferent firing rate in Hz per channel, keyed `<trunk>/<name>`.

    The key carries the trunk because the receiving model routes by trunk and
    delays by fibre class, and a bare channel name would make it guess.
    """
    assert_sources(values)
    out = {}
    for c in CHANNELS:
        r = c.rate_hz(values)
        if not math.isfinite(r) or r < 0.0 or r > c.max_hz + 1e-9:
            raise ValueError(f"channel {c.name!r} produced {r!r}, outside "
                             f"[0, {c.max_hz}]")
        out[f"{c.trunk}/{c.name}"] = r
    return out


def channel_table() -> list[dict]:
    """The declaration, as rows -- for a receipt, a doc or the other repo."""
    return [dict(name=c.name, trunk=c.trunk, fibre=c.fibre, units=c.units,
                 sources=list(c.sources), threshold=c.threshold,
                 saturation=c.saturation, max_hz=c.max_hz, invert=c.invert,
                 native_equivalent=c.native_equivalent, doc=c.doc)
            for c in CHANNELS]


# ---------------------------------------------------------------------------
# derived scalars
# ---------------------------------------------------------------------------
#
# READ THIS BEFORE USING THE NAMES.
#
# `endurance`, `discomfort` and `sensation` below are NAMED PROJECTIONS OF
# MEASURED PHYSIOLOGICAL STATE.  They are arithmetic on quantities BioGears
# integrated, and nothing about that arithmetic is evidence that the simulation
# feels tired, uncomfortable, or anything at all.  The names are chosen because
# they say which physiological question each projection answers, and they are a
# liability precisely because they are readable: a reader who sees "discomfort
# 0.4" will supply a meaning the number does not have.
#
# What each one IS, stated so the name cannot smuggle more:
#   endurance   hours of carbohydrate substrate remaining at the current
#               metabolic rate.  A ratio of two engine outputs.
#   discomfort  a weighted sum of the afferent rates above, weighted toward the
#               unmyelinated classes.  A summary of traffic on a wire.
#   sensation   a low-dimensional linear projection of the same afferent vector,
#               with the basis fitted on a training split.  A compression.
#
# None of the three is a claim about experience, and none of them should be
# reported as one.  They are here because a model that receives fifteen visceral
# rates and is asked to predict nothing has not been asked a question.

#: kJ per gram of the two substrate pools.  glycogen is stored hydrated, so the
#: usable figure per gram of stored glycogen is the standard 17 kJ/g of the
#: carbohydrate itself; fat is 37 kJ/g.
KJ_PER_G_GLYCOGEN = 17.0
KJ_PER_G_FAT = 37.0

#: fibre-class weights for `discomfort`.  unmyelinated C fibres carry the great
#: majority of visceral nociceptive traffic, A-delta a minority, A-beta
#: essentially none -- so this is a statement about which classes report
#: aversive state, and it is why the fibre class has to be a declared field.
FIBRE_AVERSIVE_WEIGHT = {"c": 1.0, "adelta": 0.4, "abeta": 0.05}


def endurance_hours(values: dict) -> float:
    """Hours of carbohydrate substrate remaining at the current metabolic rate.

    A named projection of measured state, not a claim about felt endurance.
    Carbohydrate only: fat oxidation cannot be recruited fast enough to be the
    binding constraint on the timescale this matters, and mixing the two pools
    into one number would report a body that can run for four days.
    `stored_fat_g` is returned separately by `derived_scalars` for anyone who
    wants the other bound.
    """
    liver = _read(values, "liver_glycogen_g", "endurance")
    muscle = _read(values, "muscle_glycogen_g", "endurance")
    rate_w = _read(values, "metabolic_rate_w", "endurance")
    if rate_w <= 0.0:
        raise ValueError("metabolic_rate_w must be positive to divide by")
    kj = KJ_PER_G_GLYCOGEN * (liver + muscle)
    return kj * 1000.0 / rate_w / 3600.0


def discomfort_index(rates_hz: dict[str, float]) -> float:
    """Visceral afferent load, weighted toward the unmyelinated classes, in [0,1].

    A named projection of measured afferent traffic, not a claim about felt
    discomfort.  Normalised by the load the same weighting would produce if
    every channel were saturated, so the number is a fraction of a declared
    ceiling rather than an unbounded sum whose scale depends on how many
    channels happen to be declared.
    """
    num = den = 0.0
    for c in CHANNELS:
        w = FIBRE_AVERSIVE_WEIGHT[c.fibre]
        key = f"{c.trunk}/{c.name}"
        if key not in rates_hz:
            raise KeyError(f"discomfort: no rate for {key!r}")
        num += w * rates_hz[key]
        den += w * c.max_hz
    return num / den


def derived_scalars(values: dict, rates_hz: dict[str, float] | None = None) -> dict:
    """The scalars, with their units, from one snapshot.

    `sensation` is not here: it is a projection whose basis has to be fitted on
    a training split, so it belongs to the corpus builder and not to a
    single-frame function.  Putting it here would have meant fitting a basis on
    whatever frame happened to be passed.
    """
    r = visceral_afferent_rates(values) if rates_hz is None else rates_hz
    return {
        "endurance_h": endurance_hours(values),
        "discomfort": discomfort_index(r),
        # context, reported alongside rather than folded in
        "fat_reserve_kj": KJ_PER_G_FAT * _read(values, "stored_fat_g", "context"),
        "metabolic_rate_w": _read(values, "metabolic_rate_w", "context"),
        "fatigue_fraction": _read(values, "fatigue_fraction", "context"),
        "core_temperature_c": _read(values, "core_temperature_c", "context"),
    }


SCALAR_NAMES = ("endurance_h", "discomfort")


def describe() -> str:
    lines = [f"{len(CHANNELS)} visceral afferent channels on {len(TRUNKS)} trunks",
             ""]
    lines.append(f"  {'channel':24s} {'trunk':20s} {'fibre':7s} {'units':8s} "
                 f"{'range':>22s} {'max Hz':>7s}")
    for c in CHANNELS:
        rng = f"{c.threshold:g}-{c.saturation:g}" + (" (inv)" if c.invert else "")
        lines.append(f"  {c.name:24s} {c.trunk:20s} {c.fibre:7s} {c.units:8s} "
                     f"{rng:>22s} {c.max_hz:7.1f}")
    lines += ["", f"reads {len(source_keys())} native systemic keys; every one is "
                  f"checked before it is read.", "",
              "endurance / discomfort are NAMED PROJECTIONS OF MEASURED STATE.",
              "they are not claims about experience.  see the module comment."]
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
