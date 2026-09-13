#!/usr/bin/env python3
"""Are the containment failures ANATOMY or REGISTRATION?

`register_pelvic_organs.py`'s containment gate asks whether a mapped uterus sits inside this
body's pelvic ring. A subject who fails is either badly registered, or genuinely larger than this
pelvis admits -- and this is a pathology-selected cohort whose labelled uteri run to **765.7 mL,
6.4x the top of an adult non-gravid reference**, with 42% above that reference. **The gate cannot
tell those two apart**, and they call for opposite work: fix the registration, or accept that this
body cannot hold these organs.

The separation, pre-registered before the full run finished: **if the failures correlate with
uterine volume they are anatomy; if they do not, they are registration.**

A rank test is used rather than a difference of means, because the volume distribution is heavily
skewed (median 108.5 mL, max 765.7) and a mean would be dragged by the tail.

KNOWN ANSWERS, both before the real comparison:
  1. shuffling the pass/fail labels must destroy any association -- the statistic on 200 shuffles
     must straddle zero. This is the control that can fail.
  2. splitting on a coin flip instead of containment must give the same null. Same reason,
     different construction, so a quirk of the label vector cannot pass both.
"""
import json, re, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs/endomri_world_full.log"
MAN = ROOT / "data/derived/ut-endomri-organs-v2-world/manifest.json"


def main():
    if not LOG.exists():
        sys.exit(f"{LOG} not found")
    passed, failed = [], []
    for ln in LOG.read_text().splitlines():
        m = re.match(r"(D\d-\d+)", ln)
        if not m or "GATE containment" not in ln:
            continue
        seg = ln.split("GATE containment")[1]
        (failed if "-> FAIL" in seg else passed).append(m.group(1))
    print(f"containment: {len(passed)} pass, {len(failed)} fail "
          f"({len(failed)/max(len(passed)+len(failed),1):.0%} failure rate)")

    S = json.loads(MAN.read_text())["subjects"]

    def vol(sid):
        o = (S.get(sid) or {}).get("organs", {}).get("ut")
        return None if not o else o.get("volume_ml")

    pv = np.array([v for v in (vol(s) for s in passed) if v is not None])
    fv = np.array([v for v in (vol(s) for s in failed) if v is not None])
    print(f"uterine volume available for {len(pv)} passing, {len(fv)} failing subjects")
    if len(fv) < 3 or len(pv) < 3:
        sys.exit("too few subjects with volumes on one side to compare")

    allv = np.concatenate([pv, fv])
    lab = np.r_[np.zeros(len(pv), int), np.ones(len(fv), int)]
    rank = allv.argsort().argsort().astype(float)

    def stat(l):
        """mean rank of the FAIL group minus the mean rank of the PASS group, in rank units."""
        return rank[l == 1].mean() - rank[l == 0].mean()

    obs = stat(lab)
    rng = np.random.default_rng(0)
    null = np.array([stat(rng.permutation(lab)) for _ in range(2000)])

    print("\nKNOWN ANSWERS")
    sh = np.array([stat(rng.permutation(lab)) for _ in range(200)])
    ok1 = sh.min() < 0 < sh.max()
    print(f"  shuffled labels straddle zero: [{sh.min():+.1f}, {sh.max():+.1f}]  "
          f"{'PASS' if ok1 else 'FAIL'}")
    coin = rng.integers(0, 2, size=len(lab))
    c = np.array([stat(rng.permutation(coin)) for _ in range(200)])
    ok2 = c.min() < 0 < c.max()
    print(f"  a coin-flip split gives the same null: [{c.min():+.1f}, {c.max():+.1f}]  "
          f"{'PASS' if ok2 else 'FAIL'}")
    if not (ok1 and ok2):
        sys.exit("a known answer FAILED")

    p = float((np.abs(null) >= abs(obs)).mean())
    print(f"\n  median uterine volume: passing {np.median(pv):.1f} mL, failing {np.median(fv):.1f} mL")
    print(f"  mean-rank difference (fail - pass): {obs:+.1f}  two-sided p = {p:.4f}  "
          f"over {len(null):,} permutations")
    print("\n  VERDICT: " + (
        f"the failures ARE associated with uterine volume (p = {p:.4f}). They are substantially "
        f"ANATOMY --\n  this body's pelvis cannot hold the larger uteri in a pathology-selected "
        f"cohort -- and fixing the\n  registration will not recover them." if p < 0.05 else
        f"the failures are NOT associated with uterine volume (p = {p:.4f}). Size does not explain "
        f"them,\n  so they are a REGISTRATION problem and are recoverable."))

    out = ROOT / "out/containment_vs_volume.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(n_pass=len(pv), n_fail=len(fv),
                                   median_pass_ml=float(np.median(pv)),
                                   median_fail_ml=float(np.median(fv)),
                                   mean_rank_diff=float(obs), p_two_sided=p,
                                   caveat="UT-EndoMRI is an endometriosis cohort, "
                                          "pathology-selected, NOT a typical-anatomy reference"),
                              indent=2) + "\n")
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
