#!/usr/bin/env python3
"""Known answers for the convergence summary `seat_on_bed` returns.

Added 2026-09-11 with the fields themselves. They are only exercised when a drive runs — which
takes hours on the anatomical bed — so without this they would first be trusted at the moment a
result depends on them.

WHAT THEY ARE FOR. `ihm/assembly/sliding_contact.py` used to return
`converged=record[-1]['converged']`: the LAST load step's flag alone. A drive that pushed through
several `Newton 300 UNCONVERGED` steps and happened to converge on its final one reported
`converged: true`, and nothing downstream could tell — while the per-breast gates (a)–(d) never
ask about convergence at all, so a state assembled from non-equilibria could pass all four.

The second case below is the whole point: the old field alone says `converged=1` on a run with an
unconverged step in the middle of it.

The fourth is the one that would otherwise raise rather than mislead — `max()` over an empty
sequence — and it is here because a drive where NOTHING converges is exactly the case this
summary exists to describe.
"""
import sys


def summarise(record):
    """The expressions `seat_on_bed` returns, isolated so they can be checked.

    Kept textually identical to the ones in `sliding_contact.py`; if that changes, change this
    and re-run, because a test of a copy that has drifted is worse than no test.
    """
    return dict(
        converged=record[-1]["converged"],
        all_steps_converged=all(bool(s["converged"]) for s in record),
        n_unconverged_steps=sum(1 for s in record if not s["converged"]),
        last_converged_fraction=max([s["fraction"] for s in record if s["converged"]],
                                    default=None),
        max_fraction_reached=max([s["fraction"] for s in record], default=None))


def rec(*pairs):
    return [dict(fraction=f, converged=c) for f, c in pairs]


CASES = [
    ("every step converged",
     rec((0.1, 1), (0.2, 1), (0.3, 1)),
     dict(all_steps_converged=True, n_unconverged_steps=0,
          last_converged_fraction=0.3, max_fraction_reached=0.3)),

    ("THE BUG: the last step converged, an earlier one did not",
     rec((0.1, 1), (0.2, 0), (0.3, 1)),
     dict(all_steps_converged=False, n_unconverged_steps=1,
          last_converged_fraction=0.3, max_fraction_reached=0.3)),

    ("the 2026-09-11 s1159-left shape: converged to 0.1094, then not",
     rec((0.0625, 1), (0.1094, 1), (0.1797, 0), (0.2324, 0), (0.3115, 0)),
     dict(all_steps_converged=False, n_unconverged_steps=3,
          last_converged_fraction=0.1094, max_fraction_reached=0.3115)),

    ("EDGE: nothing converged -- None, not a raised ValueError",
     rec((0.1, 0), (0.2, 0)),
     dict(all_steps_converged=False, n_unconverged_steps=2,
          last_converged_fraction=None, max_fraction_reached=0.2)),
]


def main():
    ok = True
    for name, record, want in CASES:
        got = summarise(record)
        good = all(got[k] == v for k, v in want.items())
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {name}")
        if not good:
            for k, v in want.items():
                if got[k] != v:
                    print(f"        {k}: want {v!r}, got {got[k]!r}")
    print("\n  The second case is the one that matters: the old `converged` field alone reports")
    print("  True on a run with an unconverged step in it. `all_steps_converged` is what catches it.")
    print(f"\n{'ALL KNOWN ANSWERS PASS' if ok else 'A KNOWN ANSWER FAILED'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
