#!/usr/bin/env python3
"""Collect a motor corpus from the engineered LQR, for IBM-1's supervision soup.

The IBM curriculum trains every materialization the same way: a fixed paired
corpus, a task head, and a loss against an explicit baseline. Vision is
(image, evoked EEG). Hearing is (cochleagram, MEG). This makes the body the same
shape: (muscle state, motor command), so it can be one more objective in the
soup rather than a separate pipeline with its own physics loop.

**The teacher has to be the LQR, and that is the whole design constraint.**
Three candidates were measured:

  bare postural servo, no baselines     falls at 1.19 s
  postural servo + equilibrium          falls at 1.97 s
  engineered LQR                        holds

The first two produce a corpus of a body falling over, which is what the earlier
fine-tuning attempt cloned -- and every ablation arm then lost to predicting the
training mean, because there is nothing to learn from a controller that is
failing. The LQR is also the only candidate INDEPENDENT of the IBM kernel: the
cortical stance controller holds five seconds, but it is built from that kernel
by an offline decoder fit, so cloning it would be circular.

Physics runs once, here. The curriculum then trains on the arrays like any other
corpus, which is what keeps a 1.45 s/step training loop from becoming a 300 s one.
"""
import argparse, json, shutil, sys, tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.native.mechanical_stream import NativeMechanicalStream

BUNDLE = ROOT / "data/models/engineering_stance_v1"


def load_lqr():
    src = ROOT / "ihm/native/stance_lqr.py"
    ns = {"__name__": "_frozen_stance_lqr", "__file__": str(src)}
    exec(compile(src.read_bytes(), str(src), "exec"), ns)
    reg = json.loads((BUNDLE / "registration.json").read_bytes())
    with tempfile.TemporaryDirectory(prefix="lqr-") as tmp:
        p = Path(tmp) / "linearization.npz"
        p.write_bytes((BUNDLE / "linearization.npz").read_bytes())
        return ns["NativeStanceLQR"](p, model_sha256=reg["model_sha256"], dt_s=.01)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seconds", type=float, default=20.0)
    ap.add_argument("--push-every", type=int, default=50,
                    help="steps between pelvis perturbations.  a corpus recorded "
                         "at equilibrium has almost no command variance and "
                         "cloning it is cloning a constant")
    ap.add_argument("--push-n", type=float, default=6.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="data/derived/body-corpus")
    a = ap.parse_args()

    out = ROOT / a.out
    shutil.rmtree(out, ignore_errors=True); out.mkdir(parents=True)
    pose = json.loads((BUNDLE / "initial_pose.json").read_bytes())
    lqr = load_lqr()
    rng = np.random.default_rng(a.seed)

    # the LQR artifact was identified against the stance bundle's augmented
    # model, so the stream has to be built on the SAME registration or its muscle
    # catalog differs and the controller refuses -- correctly, since a gain
    # matrix indexed by a different muscle set is meaningless.
    stream = NativeMechanicalStream(
        ROOT, out / "stream", environment="upright", target_mass_kg=70,
        initial_pose=pose,
        augmented_registration="data/models/engineering_stance_v1/registration.json")
    try:
        names = list(lqr.muscle_names)
        X, Y, T = [], [], []
        force = ()
        while stream.state["time_s"] < a.seconds - 1e-9:
            st = stream.state
            cmd, _ = lqr.commands(st)
            feat = []
            for n in names:
                m = st["muscles"][n]
                o = float(m["optimal_fiber_length_m"])
                feat += [m["fiber_length_m"] / o, m["fiber_velocity_m_s"] / o]
            X.append(feat); Y.append([cmd[n] for n in names]); T.append(st["time_s"])
            k = len(X)
            if k % a.push_every == 1:
                force = ({"body": "pelvis", "point_m": [0., 0., 0.],
                          "force_n": [float(rng.uniform(-a.push_n, a.push_n)), 0.,
                                      float(rng.uniform(-a.push_n, a.push_n))]},)
            elif k % a.push_every == 12:
                force = ()
            stream.advance(.01, actuation=cmd, forces=force)
            if st["coordinates"]["pelvis_ty"]["value"] < .6:
                print(f"  FELL at step {k}", flush=True); break
        X = np.array(X, np.float32); Y = np.array(Y, np.float32)
        np.save(out / "state.npy", X); np.save(out / "command.npy", Y)
        (out / "meta.json").write_text(json.dumps({
            "muscles": names, "n_steps": len(X), "dt_s": .01,
            "teacher": "engineered LQR (engineering_stance_v1), independent of the IBM kernel",
            "command_variance": float(Y.var(0).mean()),
            "state_variance": float(X.var(0).mean()),
            "pushed": True, "push_n": a.push_n, "push_every": a.push_every,
        }, indent=2))
        print(f"{len(X)} steps, {len(names)} muscles, {X.shape[1]} features")
        print(f"  command variance {Y.var(0).mean():.8f}   state variance {X.var(0).mean():.6f}")
        print(f"  final pelvis_ty {stream.state['coordinates']['pelvis_ty']['value']:.4f}")
        print(f"wrote {out}")
    finally:
        stream.close()


if __name__ == "__main__":
    main()
