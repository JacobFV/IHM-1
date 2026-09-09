#!/usr/bin/env python3
"""Render the full-stack embodied run against its severed control, side by side.

Draws the recorded joint kinematics from `unified-world-ofb2dp7z` -- the run in
which the trained cortex, the 1 ms segmental cord, the native OpenSim body,
BioGears physiology and the bedroom completed five seconds and 250 exchanges
with synchronised clocks, under a 5 N pelvis push.

Two arms, same push, same everything else:

  full     cortical dynamics in the loop
  severed  cortical dynamics removed

The severed arm is the control that makes the left panel mean something. Without
it a standing figure is just a standing figure.

**What this video does and does not show.** It shows that the cortical dynamics
are load-bearing for balance: severed falls, full holds. It does NOT show that
the kernel's learned content matters -- a matched control with the kernel's site
rows permuted recovers the same push at 1.7452 mm against 1.7344 mm, so what is
doing the work is the E/I network being a well-conditioned filter, not what it
learned from vision and audio. The caption says so.
"""
import argparse, json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# a planar stick figure in the sagittal plane, driven by the joints the trace has
SEGMENTS = [
    ("pelvis", "torso", 0.60), ("torso", "head", 0.24),
    ("pelvis", "hip_r", 0.09), ("hip_r", "knee_r", 0.43), ("knee_r", "ankle_r", 0.43),
    ("ankle_r", "toe_r", 0.16),
    ("pelvis", "hip_l", 0.09), ("hip_l", "knee_l", 0.43), ("knee_l", "ankle_l", 0.43),
    ("ankle_l", "toe_l", 0.16),
    ("torso", "shoulder_r", 0.16), ("shoulder_r", "elbow_r", 0.30),
    ("elbow_r", "hand_r", 0.26),
    ("torso", "shoulder_l", 0.16), ("shoulder_l", "elbow_l", 0.30),
    ("elbow_l", "hand_l", 0.26),
]


def pose(frame):
    """joint angles -> 2D points in the sagittal plane, metres."""
    j = frame["joints"]
    g = lambda k, d=0.0: float(j[k]["value"]) if k in j else d
    px, py = g("pelvis_tx"), g("pelvis_ty")
    tilt = g("pelvis_tilt")
    P = {}
    P["pelvis"] = np.array([px, py])
    P["torso"] = P["pelvis"] + 0.60 * np.array([math.sin(-tilt + g("lumbar_extension")),
                                                math.cos(-tilt + g("lumbar_extension"))])
    P["head"] = P["torso"] + 0.24 * np.array([math.sin(-tilt), math.cos(-tilt)])
    # arms, hung from the shoulders, for readability as a human figure
    for side, sgn in (("r", 1), ("l", -1)):
        sh = P["torso"] + np.array([0.16 * sgn, -0.04])
        P[f"shoulder_{side}"] = sh
        a_s = -tilt + 0.12 * sgn
        el = sh + 0.30 * np.array([math.sin(a_s), -math.cos(a_s)])
        P[f"elbow_{side}"] = el
        P[f"hand_{side}"] = el + 0.26 * np.array([math.sin(a_s), -math.cos(a_s)])
    for side, sgn in (("r", 1), ("l", -1)):
        # a small lateral offset per leg.  the sagittal projection puts both legs
        # on the same line in stance, which renders a standing body as a single
        # stroke and makes it look like nothing is there.  the offset is drawing,
        # not kinematics -- the joint angles are untouched.
        hip = P["pelvis"] + np.array([0.055 * sgn, -0.05])
        P[f"hip_{side}"] = hip
        a_h = -tilt - g(f"hip_flexion_{side}")
        knee = hip + 0.43 * np.array([math.sin(a_h), -math.cos(a_h)])
        P[f"knee_{side}"] = knee
        a_k = a_h + g(f"knee_angle_{side}")
        ank = knee + 0.43 * np.array([math.sin(a_k), -math.cos(a_k)])
        P[f"ankle_{side}"] = ank
        a_a = a_k + g(f"ankle_angle_{side}") + math.pi / 2
        P[f"toe_{side}"] = ank + 0.16 * np.array([math.sin(a_a), -math.cos(a_a)])
    return P


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default="data/derived/unified-world-ofb2dp7z")
    ap.add_argument("--out", default="artifacts/embodied_full_vs_severed.gif")
    ap.add_argument("--fps", type=int, default=25)
    ap.add_argument("--stride", type=int, default=2)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    run = ROOT / a.run
    full = json.load(open(run / "full/trace.json"))
    sev = json.load(open(run / "sever/trace.json"))
    n = max(len(full), len(sev))
    idx = list(range(0, n, a.stride))

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 5.4), facecolor="#0f151b")
    titles = ("cortex in the loop", "cortex severed")
    cols = ("#4fc3ae", "#f2879f")
    for ax, t, c in zip(axes, titles, cols):
        ax.set_facecolor("#0f151b")
        ax.set_xlim(-0.9, 0.9); ax.set_ylim(-0.15, 2.0)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(t, color=c, fontsize=12, family="monospace", pad=8)
        ax.axhline(0, color="#2a3742", lw=1.5)

    lines, dots, labels = [], [], []
    for ax, c in zip(axes, cols):
        lines.append([ax.plot([], [], color=c, lw=3, solid_capstyle="round")[0]
                      for _ in SEGMENTS])
        dots.append(ax.plot([], [], "o", color=c, ms=6)[0])
        labels.append(ax.text(-0.85, 1.85, "", color="#a9b8c5", fontsize=9,
                              family="monospace", va="top"))

    def draw(k):
        for panel, trace in enumerate((full, sev)):
            i = min(k, len(trace) - 1)
            f = trace[i]
            P = pose(f)
            for ln, (a_, b_, _) in zip(lines[panel], SEGMENTS):
                if a_ in P and b_ in P:
                    ln.set_data([P[a_][0], P[b_][0]], [P[a_][1], P[b_][1]])
            dots[panel].set_data([P["head"][0]], [P["head"][1]])
            d = f["com_horizontal_displacement_m"] * 1000
            state = "FALLEN" if f.get("fallen") else "upright"
            labels[panel].set_text(
                f"t {f['time_s']:5.2f} s\nCOM {d:8.2f} mm\npelvis {f['pelvis_height_m']:.3f} m\n{state}")
        return [l for p in lines for l in p] + dots + labels

    fig.suptitle("IBM-1 cortex + segmental cord + native body + physiology + world",
                 color="#e6edf3", fontsize=12, family="monospace", y=0.97)
    fig.text(0.5, 0.03,
             "5 N pelvis push, 250 exchanges, synchronised clocks.  the severed arm is the control.\n"
             "the cortical DYNAMICS are load-bearing; the kernel's learned CONTENT is not "
             "(a permuted kernel recovers the same push).",
             ha="center", color="#7c8b99", fontsize=7.5, family="monospace")
    anim = FuncAnimation(fig, draw, frames=idx, interval=1000 // a.fps, blit=False)
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out, writer=PillowWriter(fps=a.fps))
    print(f"wrote {out}  ({len(idx)} frames)")


if __name__ == "__main__":
    main()
