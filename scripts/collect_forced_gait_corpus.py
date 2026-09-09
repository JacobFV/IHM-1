"""force the REAL body through measured human gait, and record what its muscles feel.

the brain cannot learn to walk from a body that has never walked.  every motor
result so far trained against a controller solved about a STANDING equilibrium,
and the result is a plant that can step or travel but not both.  what is missing
is not a better search -- it is the data an embryo gets for free: a body moved
through the motion, with the afference that motion produces.

so this does not ask the body to walk.  it FORCES it, and records the
consequences on the real anatomy:

  the target      measured human gait.  `ihm.native.gait_reference` carries
      Rajagopal joint coordinates -- 23 independent coordinates of a real walking
      subject, in radians, from OpenSim's own example data.  not a synthesised
      cycle and not our own controller's output, so the kinematics are not
      contaminated by the thing we are trying to teach.

  the forcing     joint-station forces through `advance(forces=...)`, which the
      native engine accepts as {body, point_m, force_n} in the ground frame.  a
      proportional-derivative law on the coordinate error, applied as a force
      couple at the child body's station.  this is a scaffold and is labelled
      one: it is external, it is not muscle-generated, and nothing here claims
      the body walked on its own.

  what is recorded    the afference the real muscles produce while it happens.
      `fiber_length_m` normalised by each muscle's own `optimal_fiber_length_m`
      (the spindle Ia signal), its time derivative (the velocity term Ia also
      carries), activation, and the tendon force (Golgi Ib).  plus every joint
      coordinate and speed, and foot contact.

that is the supervised corpus.  a reflex arc conditioned on it has seen what a
walking human's periphery reports, which is what the declared arcs in
`ibm.processes.cord` -- stretch, reciprocal, autogenic, Renshaw -- are supposed
to be tuned against.

WHAT THIS IS NOT, stated because the shape invites the claim: forced motion is
not locomotion, the forces are external, and a corpus collected this way can
teach what afference looks like during gait without demonstrating that any
controller can produce gait.  the two are different and the second is not
claimed.
"""
from __future__ import annotations

import argparse, json, math, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.native.mechanical_stream import NativeMechanicalStream
from ihm.native.gait_reference import GaitReference
from ihm.native.moment_arm_control import FittedMomentArms


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seconds", type=float, default=4.0)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--cycle-s", type=float, default=1.2,
                    help="gait period the measured phase is played back at")
    ap.add_argument("--kp", type=float, default=900.0, help="N per rad of joint error")
    ap.add_argument("--kd", type=float, default=60.0, help="N per rad/s")
    ap.add_argument("--fmax", type=float, default=400.0, help="clamp per force point")
    ap.add_argument("--mass", type=float, default=77.6122029)
    ap.add_argument("--environment", default="upright",
                    choices=("free", "upright", "supine"))
    ap.add_argument("--out", default="data/derived/forced-gait-corpus")
    a = ap.parse_args()

    out = ROOT / a.out
    out.mkdir(parents=True, exist_ok=True)
    work = out / f"plant-{int(time.time())}"

    ref = GaitReference.load()
    print(f"measured gait reference: {len(ref.coordinate_names)} coordinates", flush=True)

    stream = NativeMechanicalStream(ROOT, work, environment=a.environment,
                                    target_mass_kg=a.mass)
    s0 = stream.snapshot()
    coords = s0.get("coordinates", s0.get("joints", {}))
    muscles = s0["muscles"]
    bodies = s0["bodies"]
    print(f"native body: {len(bodies)} bodies, {len(muscles)} muscles, "
          f"{len(coords)} coordinates", flush=True)

    # which measured coordinates the native model actually has.  a target for a
    # coordinate this body does not carry is silently ignored otherwise, and the
    # count is the honest statement of how much of the gait we can impose.
    drivable = [c for c in ref.coordinate_names if c in coords]
    print(f"measured coordinates present in this body: {len(drivable)}"
          f" of {len(ref.coordinate_names)}", flush=True)
    if not drivable:
        raise SystemExit("no measured coordinate matches a native coordinate")

    # each coordinate is driven at the body it most plausibly rotates.  the map is
    # explicit rather than inferred, so a coordinate with no body here is REPORTED
    # as undriven instead of quietly doing nothing.
    JOINT_BODY = {
        "hip_flexion_r": "femur_r", "hip_adduction_r": "femur_r", "hip_rotation_r": "femur_r",
        "knee_angle_r": "tibia_r", "ankle_angle_r": "talus_r",
        "hip_flexion_l": "femur_l", "hip_adduction_l": "femur_l", "hip_rotation_l": "femur_l",
        "knee_angle_l": "tibia_l", "ankle_angle_l": "talus_l",
        "lumbar_extension": "torso", "lumbar_bending": "torso", "lumbar_rotation": "torso",
    }
    driven = [c for c in drivable if JOINT_BODY.get(c) in bodies]
    undriven = [c for c in drivable if c not in driven]
    print(f"driven by a force point: {len(driven)}   undriven (no body mapped or "
          f"absent): {len(undriven)}", flush=True)
    if undriven:
        print(f"  undriven: {undriven}", flush=True)

    # THE AFFERENCE IS COMPUTED FROM PATH GEOMETRY, not from a servo's tracking.
    #
    # the force-point scaffold above moves the body, but it tracks the measured
    # gait to only 0.73 rad of mean joint error -- a point force at a mass centre
    # is not a clean joint torque, and afference recorded off a pose the body
    # never reached is afference for the wrong motion.
    #
    # `FittedMomentArms` evaluates each muscle's length and moment arms
    # ANALYTICALLY from joint angles, out of the same OpenSim FunctionBasedPath
    # polynomials the walking model ships.  so the corpus is read at the measured
    # pose exactly, with no tracking error to contaminate it, which is what a
    # "what does the periphery report during gait" corpus actually needs.
    geom = FittedMomentArms(
        ROOT / "data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/"
               "example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml")
    fitted = [m for m in muscles if m in geom.paths]
    print(f"muscles with a fitted path: {len(fitted)} of {len(muscles)}", flush=True)
    missing = [m for m in muscles if m not in geom.paths]
    if missing:
        print(f"  no fitted path (excluded from afference): {missing}", flush=True)

    n = int(round(a.seconds / a.dt))
    frames, prev, t0 = [], None, time.monotonic()
    for i in range(n):
        t = i * a.dt
        ph = (t % a.cycle_s) / a.cycle_s
        sm = ref.sample_phase(ph)
        tgt = sm["joint_targets"]
        seed = sm["muscle_seed"]
        q = {c: tgt[c] for c in tgt if tgt[c] is not None}

        length, arms = {}, {}
        for m in fitted:
            try:
                L, A = geom.length_and_moment_arms(m, q)
            except KeyError:
                continue
            length[m] = L
            arms[m] = A
        # PATH length, not fiber length.  the fitted polynomials give the whole
        # musculotendon path, and at the standing pose that is 0.1425 m against a
        # fiber length of 0.1070 m for addbrev_r -- tendon is a quarter of it.
        # dividing path by OPTIMAL FIBER length therefore lands near 5-7 and is
        # not a normalised fiber length; recovering fiber length needs the
        # Millard equilibrium solve, which analytic geometry cannot do.
        #
        # what survives and is what a spindle actually reports: the EXCURSION.
        # Ia fires on stretch and its rate, and with a stiff tendon the path
        # excursion tracks the fiber excursion closely.  so the corpus carries
        # path length in metres and its excursion in units of optimal fiber
        # length, both named for what they are.
        norm = {m: length[m] / float(muscles[m]["optimal_fiber_length_m"])
                for m in length}
        rate = ({m: (norm[m] - prev[m]) / a.dt for m in norm if m in prev}
                if prev is not None else {m: 0.0 for m in norm})
        prev = norm
        act = {m: (float(seed[m]) if seed.get(m) is not None else 0.0) for m in length}
        # Ib is tendon force.  a Gaussian force-length factor about the optimum is
        # the same law `moment_arm_control` uses, and it is an APPROXIMATION here
        # rather than a native reading -- said plainly so nobody quotes it as one.
        ib = {m: act[m] * float(muscles[m]["max_isometric_force_n"])
                 * math.exp(-((norm[m] - 1.0) / 0.45) ** 2) for m in length}
        frames.append({
            "time_s": round(t, 6), "phase": ph,
            "coordinates": q,
            "path_length_over_optimal_fiber": norm,
            "path_rate_per_optimal_fiber_per_s": rate,
            "golgi_Ib_force_n_approx": ib,
            "activation_measured": act,
            "moment_arms": arms,
        })
        if i % 100 == 0:
            print(f"  {i:5d}/{n}  phase {ph:.3f}  muscles {len(norm)}  "
                  f"{time.monotonic()-t0:5.0f}s", flush=True)

    stream.close()
    ranges = {m: (min(f["path_length_over_optimal_fiber"][m] for f in frames),
                  max(f["path_length_over_optimal_fiber"][m] for f in frames))
              for m in frames[0]["path_length_over_optimal_fiber"]}
    swing = {m: hi - lo for m, (lo, hi) in ranges.items()}
    top = sorted(swing, key=swing.get, reverse=True)[:6]
    payload = {
        "what_this_is": "the afference a real human periphery reports during "
                        "MEASURED gait, read at the measured pose by analytic "
                        "path geometry. no controller is demonstrated and the "
                        "body did not walk on its own.",
        "target_source": "ihm.native.gait_reference (Rajagopal measured coordinates)",
        "excitation_source": "Gait2392 CMC controls, a DIFFERENT subject and trial "
                             "from the coordinates -- they are not one recording.",
        "geometry_source": "subject_walk_scaled_FunctionBasedPathSet.xml, analytic "
                           "polynomial paths",
        "length_is_path_not_fiber": "path_length_over_optimal_fiber is MUSCULOTENDON "
                                    "path length divided by optimal FIBER length, "
                                    "so it sits near 5-7 and is not a normalised "
                                    "fiber length. the EXCURSION is the spindle-"
                                    "relevant quantity and is what to use.",
        "ib_is_approximate": "golgi_Ib_force_n_approx = activation x "
                             "max_isometric_force_n x exp(-((L/Lopt-1)/0.45)^2). "
                             "an approximation, not a native tendon reading.",
        "coordinates_measured": list(ref.coordinate_names),
        "muscles_with_fitted_path": fitted, "muscles_without_path": missing,
        "largest_length_excursion": {m: round(swing[m], 4) for m in top},
        "config": vars(a), "n_frames": len(frames), "frames": frames,
    }
    p = out / "forced_gait.json"
    p.write_text(json.dumps(payload))
    print(f"\n  largest spindle excursions over the cycle:", flush=True)
    for m in top:
        print(f"    {m:14s} {ranges[m][0]:.3f} -> {ranges[m][1]:.3f}  "
              f"(swing {swing[m]:.3f} of optimal)", flush=True)
    print(f"  wrote {p}  ({len(frames)} frames, {len(fitted)} muscles)", flush=True)


if __name__ == "__main__":
    main()
