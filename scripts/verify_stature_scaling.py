"""Gate the stature scaling, including the case where it is deliberately broken.

    python scripts/verify_stature_scaling.py

A gate that has never been seen to fail is not a gate.  The load-bearing arm
here is ``sabotage``: it scales the model and leaves the fitted path set at the
source subject's, which is exactly the mistake this whole mechanism exists to
prevent, and it reports which gates notice and which do not.  The answer is the
point.  Stature still passes; the muscles are wrong.

Also run:

* an identity arm at scale 1.0, whose answer is known in advance -- every
  measured quantity must come back bit-identical to the base model's;
* a monotonicity arm across five scales, so a sign error cannot hide;
* the refusal arm, which asks ``NativeMechanicalStream`` to accept a scaled
  model with no scaled path set and requires it to say no.

No native process is started.  Everything below is measured on the XML.
"""
from pathlib import Path
import json, shutil, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.body_parameters import MECHANICAL_STATURE_M  # noqa: E402
from ihm.native.model_scaling import head_marker_height_m  # noqa: E402
from ihm.native.moment_arm_control import FittedMomentArms  # noqa: E402
from ihm.body_parameters import MECHANICAL_TARGET_MASS_KG
import scripts.materialize_stature_variant as msv  # noqa: E402

SCALES = (0.80, 0.90, 1.0, 1.05, 1.13)


def build(scale, output):
    shutil.rmtree(output, ignore_errors=True)
    return msv.materialize(ROOT, output, scale=scale)


def sabotage(scale, output):
    """Scale the model; leave the fitted path set at the source subject's."""
    record = build(scale, output)
    shutil.copy(ROOT / msv.RAW / msv.PATHSET, output / msv.PATHSET)
    base = json.loads((ROOT / msv.BASE / 'registration.json').read_text())
    catalog = json.loads((ROOT / base['catalog_path']).read_text())
    return msv.gate(ROOT / base['model_path'], ROOT / msv.RAW / msv.PATHSET, catalog,
                    output / 'subject_scaled_stature.osim', output / msv.PATHSET,
                    json.loads((output / 'catalog.json').read_text()), scale)


def refusal(output):
    """The stream must refuse a scaled model that carries no scaled path set."""
    from ihm.native.mechanical_stream import NativeMechanicalStream
    record = json.loads((output / 'registration.json').read_text())
    record['source_overrides'] = {}
    broken = output / 'registration_without_pathset.json'
    broken.write_text(json.dumps(record, indent=2) + '\n')
    fresh = output / 'refused-native-run'
    shutil.rmtree(fresh, ignore_errors=True)
    try:
        NativeMechanicalStream(ROOT, fresh, environment='supine',
                               target_mass_kg=MECHANICAL_TARGET_MASS_KG,
                               augmented_registration=str(broken.relative_to(ROOT)))
    except ValueError as error:
        return str(error)
    finally:
        shutil.rmtree(fresh, ignore_errors=True)
    return None


def main():
    work = Path(tempfile.mkdtemp(dir=ROOT / 'data/derived', prefix='stature-verify-'))
    summary = {'schema': 'ihm.stature-scaling-verification.v1', 'arms': {}}
    try:
        # --- identity ---------------------------------------------------
        record = build(1.0, work / 'identity')
        base = json.loads((ROOT / msv.BASE / 'registration.json').read_text())
        base_stature = head_marker_height_m(ROOT / base['model_path'])
        identity_stature = head_marker_height_m(work / 'identity/subject_scaled_stature.osim')
        base_paths = FittedMomentArms(ROOT / msv.RAW / msv.PATHSET)
        identity_paths = FittedMomentArms(work / 'identity' / msv.PATHSET)
        pose = msv.reference_pose(ROOT / base['model_path'])
        drift = max(abs(identity_paths.length_and_moment_arms(name, pose)[0]
                        / base_paths.length_and_moment_arms(name, pose)[0] - 1)
                    for name in base_paths.paths)
        summary['arms']['identity'] = {
            'note': 'scale 1.0. The answer is known in advance: nothing moves.',
            'stature_relative_error': abs(identity_stature / base_stature - 1),
            'worst_path_length_relative_error': drift,
            'passed': abs(identity_stature / base_stature - 1) < 1e-12 and drift < 1e-12}

        # --- monotonicity -----------------------------------------------
        rows = []
        for scale in SCALES:
            out = work / ('scale_%s' % ('%.2f' % scale).replace('.', 'p'))
            record = build(scale, out)
            paths = FittedMomentArms(out / msv.PATHSET)
            rows.append({
                'scale': scale,
                'stature_m': head_marker_height_m(out / 'subject_scaled_stature.osim'),
                'soleus_r_path_m': paths.length_and_moment_arms('soleus_r', pose)[0],
                'gates_passed': all(g['passed'] for g in record['gates'])})
        statures = [row['stature_m'] for row in rows]
        soleus = [row['soleus_r_path_m'] for row in rows]
        summary['arms']['monotonic'] = {
            'note': 'five scales; stature and one muscle path must both rise with it',
            'rows': rows,
            'passed': (all(b > a for a, b in zip(statures, statures[1:]))
                       and all(b > a for a, b in zip(soleus, soleus[1:]))
                       and all(row['gates_passed'] for row in rows))}

        # --- sabotage ---------------------------------------------------
        broken = sabotage(1.13, work / 'sabotage')
        caught = [g['gate'] for g in broken if not g['passed']]
        missed = [g['gate'] for g in broken if g['passed']]
        summary['arms']['sabotage'] = {
            'note': ('model scaled to 1.13x, fitted path set left at the source '
                     "subject's. This is the failure the mechanism exists to "
                     'prevent, injected on purpose.'),
            'scale': 1.13,
            'gates_that_caught_it': caught,
            'gates_that_did_not': missed,
            'detail': broken,
            'passed': len(caught) >= 1 and 'stature scales' in missed}

        # --- refusal ----------------------------------------------------
        message = refusal(work / 'scale_1p13')
        summary['arms']['refusal'] = {
            'note': ('NativeMechanicalStream given a registration with '
                     'geometric_scale 1.13 and no scaled path set'),
            'raised': message, 'passed': message is not None}

        # The finding the sabotage arm exists to produce, stated rather than
        # left for a reader to notice in the JSON.
        band = next(g for g in broken
                    if g['gate'].startswith('path length over optimal fibre'))
        summary['findings'] = [
            'The 0.2-20 path-over-optimal-fibre band is NOT sufficient to catch '
            'a mis-scaled body. On the sabotage arm the ratio drifted %.1f%% and '
            'stayed inside the band at [%.3f, %.3f]. The band is a corpus '
            'quality check -- it catches a motion file holding accelerations '
            'under coordinate names -- and it was never a scaling check. Treat '
            'it as necessary and nowhere near sufficient.'
            % (100 * band['drift'], band['measured_min'], band['measured_max']),
            'Stature is also not sufficient: it passed exactly on the sabotage '
            'arm. A body can be the requested height and have every lower-limb '
            'muscle at the wrong length.',
            'What catches it is a quantity that reads BOTH files: the fitted '
            'path length ratio, and the ratio of that length to the model\'s own '
            'path-point polyline. The second is the one that would also survive '
            'the two files being scaled consistently and both being wrong.']

        summary['passed'] = all(arm['passed'] for arm in summary['arms'].values())
        print(json.dumps(summary, indent=2))
        return 0 if summary['passed'] else 1
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    sys.exit(main())
