#!/usr/bin/env python3
"""Fine-tune the IBM cortical kernel on real body states, then ablate it.

The permuted-kernel control already showed that the IBM kernel's learned content
does not carry the motor task. That result is correct and it is also expected:
the kernel was trained on vision, audio and EEG retrieval, and nothing in that
corpus is motor. A kernel that has never seen a muscle has no reason to encode
one, and a permutation of it has no reason to be worse.

So the question the control leaves open is whether the kernel CAN carry motor
content, and that requires training it on one. This does that:

  1. Run the native upright body under the working postural servo, recording
     (muscle state -> servo command) at every step. Real mechanics, real fiber
     lengths and velocities, not synthetic coverage.
  2. Train the IBM cortex to reproduce those commands from that state, with the
     kernel in the gradient path. Afference enters postcentral sites, the
     dynamics run, the command is read from precentral sites -- so the signal
     must cross between two disjoint populations through the association kernel.
  3. Re-run the permuted control on the fine-tuned kernel.

The prediction that makes this falsifiable: BEFORE fine-tuning, permuted and
trained kernels should score the same, reproducing the existing negative result.
AFTER, the trained kernel should beat its own permutation -- and if it does not,
the kernel cannot carry motor content and that is a real limit rather than a
missing corpus.

Behavioural cloning is deliberately the objective rather than reward learning:
the teacher is a servo that demonstrably balances, so the target is known-good,
and the question here is representational capacity, not credit assignment.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.postural_control import PosturalController, PosturalConfig

IBM = Path.home() / "Documents/IBM-1"


def collect(root, out, seconds, dt, pushes=True, seed=0):
    """Real body, real servo, recorded state->command pairs.

    UNDER PERTURBATION by default.  the first version of this recorded a servo
    holding a body that was already balanced, and a stabilising controller at
    equilibrium barely moves: 119 steps of near-identical states with almost no
    variance in the command.  cloning that is cloning a constant, and every
    ablation arm lost to predicting the training mean.

    random pelvis pushes give the servo something to correct, which is the only
    way the recorded commands carry information about the state.
    """
    stream = NativeMechanicalStream(root, out, environment="upright", target_mass_kg=70)
    rng = np.random.default_rng(seed)
    try:
        initial = stream.snapshot()
        policy = PosturalController(initial, PosturalConfig())
        names = sorted(initial["muscles"])
        X, Y = [], []
        force = None
        while stream.state["time_s"] < seconds - 1e-9:
            st = stream.state
            cmd = policy.commands(st)
            feat = []
            for n in names:
                m = st["muscles"][n]
                o = float(m["optimal_fiber_length_m"])
                feat += [(m["fiber_length_m"] - policy.reference[n]) / o,
                         m["fiber_velocity_m_s"] / o]
            X.append(feat); Y.append([cmd[n] for n in names])
            # a fresh push every ~0.3 s, alternating direction, small enough that
            # the servo recovers rather than falling
            if pushes and len(X) % 30 == 1:
                mag = float(rng.uniform(-6.0, 6.0))
                force = [{"body": "pelvis", "point_m": [0.0, 0.0, 0.0],
                          "force_n": [mag, 0.0, float(rng.uniform(-3.5, 3.5))]}]
            elif pushes and len(X) % 30 == 8:
                force = None
            stream.advance(min(dt, seconds - st["time_s"]), actuation=cmd,
                           forces=force or ())
            if st["coordinates"]["pelvis_ty"]["value"] < .6:
                break
        return names, np.array(X, np.float32), np.array(Y, np.float32)
    finally:
        stream.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--checkpoint",
                    default=str(IBM / "ckpt/ibm1_embodied_s2k_final.pt"))
    ap.add_argument("--sites", type=int, default=2048)
    ap.add_argument("--seconds", type=float, default=3.0)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--dyn-steps", type=int, default=32,
                    help="integration steps before the motor readout.  must be "
                         "large enough that the drive crosses from postcentral "
                         "to precentral; at 2 it does not arrive at all")
    ap.add_argument("--no-pushes", action="store_true",
                    help="record the undisturbed servo.  the target variance then "
                         "collapses and every arm loses to the mean -- kept only to "
                         "reproduce that")
    ap.add_argument("--out", default="data/derived/ibm-body-finetune")
    a = ap.parse_args()

    root = Path(a.root)
    work = root / a.out
    if work.exists():
        import shutil; shutil.rmtree(work)
    work.mkdir(parents=True)

    print("collecting real body states under the postural servo...", flush=True)
    names, X, Y = collect(root, work / "stream", a.seconds, a.dt,
                          pushes=not a.no_pushes)
    print(f"  {len(X)} steps, {len(names)} muscles, {X.shape[1]} features", flush=True)
    print(f"  command variance across steps: {Y.var(0).mean():.8f} "
          f"(the thing cloning has to have)", flush=True)
    if len(X) < 32:
        raise SystemExit("too few steps collected; the body fell immediately")

    sys.path.insert(0, str(IBM))
    import importlib.util
    sp = importlib.util.spec_from_file_location(
        "ptrain", str(IBM / "scripts/pretrain_video_loop.py"))
    P = importlib.util.module_from_spec(sp); sp.loader.exec_module(P)

    ck = torch.load(a.checkpoint, map_location="cpu", weights_only=False)
    kernel = ck["dyn.embed"]
    torch.manual_seed(0)

    xt = torch.from_numpy(X); yt = torch.from_numpy(Y)
    n_tr = int(len(xt) * 0.8)
    # contiguous split with a guard band: consecutive body states are nearly
    # identical, so a random split leaks the test set into training.
    gap = max(1, len(xt) // 50)
    tr = slice(0, n_tr); te = slice(n_tr + gap, len(xt))

    def build(embed_init):
        dyn = P.CorticalDynamics(a.sites, kernel.shape[1], 48, "cpu")
        dyn.embed.data.copy_(embed_init)
        return P.SensorimotorLoop(dyn, muscles=names,
                                  afferent_channels=X.shape[1])

    # the drive enters postcentral and the command is read from precentral, and
    # those populations are DISJOINT, so the signal has to traverse the sheet.
    # measured: at n_steps=2 the precentral across-batch sd is exactly 0.0 --
    # the drive has not arrived, every kernel produces the same constant, and
    # all three ablation arms then train an identical constant-input decoder and
    # report identical numbers to six decimal places.  that is not an ablation
    # result, it is an ablation that never ran.
    n_steps = a.dyn_steps

    def train(model, epochs):
        opt = torch.optim.Adam(model.parameters(), lr=3e-3)
        for e in range(epochs):
            pred, _ = model(xt[tr], n_steps=n_steps, substeps=4)
            loss = (pred - yt[tr]).square().mean()
            opt.zero_grad(); loss.backward(); opt.step()
        with torch.no_grad():
            p, _ = model(xt[te], n_steps=n_steps, substeps=4)
            held = float((p - yt[te]).square().mean())
        return held

    base = float(yt[te].var())
    zero = float(yt[te].square().mean())
    mean_b = float((yt[te] - yt[tr].mean(0)).square().mean())
    print(f"\nbaselines on held-out: predict-mean {mean_b:.6f}  "
          f"predict-zero {zero:.6f}", flush=True)

    g = torch.Generator().manual_seed(0)
    arms = {
        "trained kernel": kernel.clone(),
        "permuted sites": kernel[torch.randperm(kernel.shape[0], generator=g)],
        "random kernel": torch.randn(kernel.shape, generator=g) * 0.02,
    }
    res = {}
    print(f"\n{'arm':18s} {'held-out MSE':>13s} {'skill vs mean':>14s}")
    for name, emb in arms.items():
        torch.manual_seed(0)
        e = emb if emb.shape[0] == a.sites else emb[:a.sites]
        held = train(build(e), a.epochs)
        res[name] = {"held_mse": held, "skill_vs_mean": 1 - held / mean_b}
        print(f"{name:18s} {held:13.6f} {1-held/mean_b:+14.4f}", flush=True)

    t, p_ = res["trained kernel"]["held_mse"], res["permuted sites"]["held_mse"]
    print(f"\n  trained vs permuted: {100*(p_-t)/p_:+.2f}% MSE difference")
    verdict = ("the fine-tuned kernel's learned content CARRIES the motor task"
               if p_ > t * 1.05 else
               "trained and permuted are equivalent -- the kernel's content still "
               "does not carry the motor task even after fine-tuning on it")
    print(f"  -> {verdict}")
    res["verdict"] = verdict
    res["n_steps"] = len(X); res["muscles"] = len(names)
    (work / "receipt.json").write_text(json.dumps(res, indent=2))
    print(f"\nwrote {work/'receipt.json'}")


if __name__ == "__main__":
    main()
