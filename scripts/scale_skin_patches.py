"""The 1,326 skin patches: area scales, count does not, and why that is a choice.

    .venv/bin/python scripts/scale_skin_patches.py --stature-m 2.03
    .venv/bin/python scripts/scale_skin_patches.py --scale 1.0        # identity

Area is not the interesting question.  A patch is a piece of the skin surface,
the surface takes the area exponent, and every patch takes it too; at 2.03 m the
exterior goes from 1.7805 m2 to 2.2718 m2.  There is nothing to decide there.

**The decision is the patch COUNT, and it is a real one.**

`scripts/build_dermatome_patches.py` makes patches by recursive bisection until
each is under an area budget -- 0.0008 m2 on hand, foot and head, 0.003 m2
elsewhere.  Re-running that rule on a larger body with the budget unchanged
would give **more patches**: about `s**2` times as many, 1,692 at 2.03 m instead
of 1,326.  Re-running it with the budget scaled as `s**2` gives the same 1,326,
larger.  Both are defensible readings of the same script and they are different
bodies, so the choice has to be made explicitly and defended.

**Chosen: count fixed at 1,326, area scales as s**2, receptor density therefore
falls as s**-2.**  Three reasons, in the order they actually decide it:

1. *A patch is an afferent channel, not a square centimetre.*  Every patch is
   innervated, carries a nerve and a relay and a cortical target, and IBM-1
   consumes them as a fixed set of inputs.  Letting the count follow body size
   would change the dimension of the brain's sensory input with stature -- a
   1.4 m body and a 2.0 m body could not run the same trained model, and the
   parameter would stop being a property of the body and start being a property
   of the interface.  The anatomy agrees with the engineering here: a taller
   person has the same 31 spinal nerve pairs and the same 32 dermatomes, and
   `docs/MILESTONES.md` counts 29 of 30 dermatome levels reached.

2. *Real receptor counts do not follow body size.*  Tactile spatial acuity is
   better on smaller fingers, and the reason is that Merkel-cell density is
   higher on them -- the innervation is laid down as a roughly fixed count and
   spread over whatever surface the person grows.  Holding the count and
   letting density fall is the reading that matches that; scaling the count is
   the reading that contradicts it.

3. *The alternative is available and is not lost.*  `--count-follows-area`
   reports what the other choice would give, in the same run, so the decision
   is visible as a decision rather than as the only thing the code could do.

**What this choice is NOT a claim about.**  It does not make the patch layout
right.  `dermatomes.json` already records that patch density here is two-tier by
AREA and not by measured receptor density, and that real fingertip innervation
exceeds trunk innervation by more than an order of magnitude.  Holding the count
fixed preserves that limitation exactly; it neither improves nor worsens it.
The honest summary is: **count fixed, area scales, density falls as s**-2, and
the density that falls was never a measured density in the first place.**
"""
from pathlib import Path
import argparse, gzip, hashlib, json, sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm import body_scaling as bs                                  # noqa: E402
from ihm.body_parameters import MECHANICAL_STATURE_M, resolve       # noqa: E402

DERMATOMES = 'data/derived/canonical/dermatomes.json'
SKIN_GEOMETRY = 'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz'
EXTERIOR = 'data/research/engineered_skin_territories/materialization.json'

#: Patch fields, by quantity.  `triangle_count` is a count of mesh triangles and
#: is invariant because the mesh is not retriangulated; `surface_normal` is a
#: unit direction and a uniform scale preserves directions.
PATCH_FIELDS = {'position_m': 'centroid', 'area_m2': 'area',
                'path_length_m': 'route_length'}
PATCH_INVARIANT = {
    'triangle_count': 'a count of mesh triangles; the mesh is not retriangulated',
    'surface_normal': 'a unit direction, preserved by a uniform scale',
}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exterior_area(root, scale):
    """The exterior skin area, integrated over the SCALED triangles themselves.

    The independent arm for this stage.  It opens the gzipped skin mesh that
    `dermatomes.json` merely names, selects the same exterior triangles by the
    same recorded id list, scales every vertex, and sums cross products.  The
    patch records know nothing about this computation and it knows nothing about
    them; they have to agree anyway, because they are two accounts of one
    surface.
    """
    with gzip.open(root / SKIN_GEOMETRY) as handle:
        mesh = json.load(handle)
    P = np.asarray(mesh['positions'], float).reshape(-1, 3) * float(scale)
    F = np.asarray(mesh['indices'], np.int64).reshape(-1, 3)
    ids = json.loads((root / EXTERIOR).read_text())['contact_eligible_triangle_ids']
    F = F[np.asarray(ids, np.int64)]
    a, b, c = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()), len(F)


def scale_patch(patch, s, origin=(0.0, 0.0, 0.0)):
    out = dict(patch)
    for field, quantity in PATCH_FIELDS.items():
        if out.get(field) is None:
            continue
        out[field] = (bs.scale_point(out[field], s, origin) if quantity == 'centroid'
                      else bs.scale(quantity, out[field], s))
    return out


def build(root, scale):
    root = Path(root)
    dermatomes = json.loads((root / DERMATOMES).read_text())
    patches = dermatomes['patches']
    scaled = [scale_patch(p, scale) for p in patches]

    coverage = dict(dermatomes['coverage'])
    for field in ('raw_source_area_m2', 'exterior_component_area_m2',
                  'patch_area_m2', 'patches_with_spinal_root_area_m2',
                  'patches_without_spinal_root_area_m2', 'unassigned_area_m2'):
        if coverage.get(field) is not None:
            coverage[field] = bs.scale('area', coverage[field], scale)
    coverage['by_dermatome'] = [dict(row, area_m2=bs.scale('area', row['area_m2'], scale))
                                for row in coverage['by_dermatome']]

    method = dict(dermatomes['method'])
    budget = {k: bs.scale('area', v, scale) for k, v in method['target_area_m2'].items()}

    base_area = dermatomes['coverage']['patch_area_m2']
    scaled_area = coverage['patch_area_m2']
    base_count = len(patches)
    alternative_count = base_count * scale ** bs.exponent('area')

    decision = dict(
        chosen='count_fixed_area_scales',
        patch_count=base_count, patch_count_at_every_stature=base_count,
        alternative='count_follows_area',
        alternative_patch_count=alternative_count,
        alternative_patch_count_rounded=int(round(alternative_count)),
        area_exponent=bs.exponent('area'),
        count_exponent=bs.exponent('patch_count'),
        receptor_density_exponent=bs.exponent('receptor_density'),
        base_density_per_m2=base_count / base_area,
        scaled_density_per_m2=base_count / scaled_area,
        density_change_percent=100 * (scale ** bs.exponent('receptor_density') - 1.0),
        area_budget_treated_as='scaled as s**%g, so re-running the bisection rule '
                               'on the scaled body reproduces the same 1,326 patches '
                               'rather than more of them' % bs.exponent('area'),
        scaled_area_budget_m2=budget,
        reasons=[
            'A patch is an afferent channel, not a square centimetre. Every one is '
            'innervated and carries a nerve, a relay and a cortical target, and '
            'IBM-1 consumes them as a fixed input set. A count that followed body '
            'size would change the dimension of the brain\'s sensory input with '
            'stature, so a 1.4 m and a 2.0 m body could not run the same trained '
            'model and the parameter would become a property of the interface.',
            'The anatomy agrees: a taller person has the same 31 spinal nerve pairs '
            'and the same 32 dermatomes. Nothing about being tall adds a root.',
            'Real receptor counts do not follow body size. Tactile spatial acuity '
            'is better on smaller fingers because Merkel-cell density is higher on '
            'them -- innervation is laid down as a roughly fixed count and spread '
            'over whatever surface the person grows.'],
        not_a_claim=(
            'This does not make the layout right. dermatomes.json already records '
            'that patch density here is two-tier by AREA and not by measured '
            'receptor density, and that real fingertip innervation exceeds trunk '
            'innervation by more than an order of magnitude. Holding the count '
            'preserves that limitation exactly. The density that falls as s**-2 '
            'was never a measured density.'))

    gates = run_gates(root, dermatomes, patches, scaled, coverage, scale, decision)
    report = dict(
        schema='ihm.skin-patch-scaling.v1', stature_scale=scale,
        patches=base_count,
        base_exterior_area_m2=base_area, scaled_exterior_area_m2=scaled_area,
        decision=decision,
        fields_scaled={f: dict(quantity=q, exponent=bs.exponent(q))
                       for f, q in PATCH_FIELDS.items()},
        fields_invariant=dict(PATCH_INVARIANT),
        gates=gates)
    out = dict(dermatomes, patches=scaled, coverage=coverage,
               method=dict(method, target_area_m2=budget),
               stature_scale=scale, scaling=decision)
    return out, report


def run_gates(root, dermatomes, patches, scaled, coverage, scale, decision):
    gates = []

    def record(name, error, tolerance, note, **extra):
        gates.append(dict(gate=name, relative_error=error, tolerance=tolerance,
                          passed=(error <= tolerance), note=note, **extra))

    if scale == 1.0:
        record('identity at s=1.0 reproduces the base exactly',
               float(sum(1 for a, b in zip(patches, scaled) if a != b)), 0.0,
               'whole-record equality over all 1,326 patches')

    # 1. Every patch area takes the area exponent.
    worst = max(abs(b['area_m2'] / (a['area_m2'] * scale ** bs.exponent('area')) - 1.0)
                for a, b in zip(patches, scaled))
    record('every patch area scales by exactly s**%g' % bs.exponent('area'),
           worst, 1e-12, 'area, not length', patches=len(patches))

    # 2. The count is the decision, and it must be visibly unchanged.
    record('patch count is unchanged', float(len(scaled) - len(patches)), 0.0,
           'the choice this script exists to make explicit: count fixed, area '
           'scales, density falls as s**%g. The alternative would give %.1f '
           'patches at this scale.' % (bs.exponent('receptor_density'),
                                       decision['alternative_patch_count']),
           chosen=decision['chosen'], alternative=decision['alternative'])

    # 3. THE INDEPENDENT ARM. Sum of 1,326 patch areas against the area of the
    #    same surface integrated over its own SCALED triangles -- a different
    #    file, a different code path, and 109,183 triangles that have never
    #    heard of a patch. Coverage is exactly 1.0 at the base, so any drift is
    #    the scaling and not the partition.
    mesh_area, triangles = exterior_area(root, scale)
    patch_sum = sum(p['area_m2'] for p in scaled)
    record('1,326 scaled patch areas sum to the SCALED mesh exterior',
           abs(patch_sum / mesh_area - 1.0), 1e-9,
           'INDEPENDENT ARM. The right-hand side is a sum of cross products over '
           '%d scaled triangles of the gzipped skin mesh, selected by the same '
           'recorded exterior id list; the left is 1,326 records that know nothing '
           'about it. Give patch area the length exponent and these differ by '
           'exactly s.' % triangles,
           patch_sum_m2=patch_sum, mesh_area_m2=mesh_area, triangles=triangles)

    # 4. Coverage is a fraction and must be exactly invariant. It is the one
    #    number here that would move under an inconsistent scaling of the two
    #    sides and cannot move under a consistent one.
    record('coverage fraction is invariant',
           abs(patch_sum / mesh_area
               - dermatomes['coverage']['patch_area_fraction_of_exterior']), 1e-9,
           'a fraction takes exponent %g; 100.00%% before, 100.00%% after'
           % bs.exponent('fraction'))

    # 5. Re-running the bisection rule with the SCALED budget must reproduce the
    #    same partition. Checked as the thing the rule actually tests: every
    #    patch under its own budget, before and after, with the same margin.
    acral = set(dermatomes['method']['acral_regions'])
    base_budget = dermatomes['method']['target_area_m2']
    scaled_budget = {k: bs.scale('area', v, scale) for k, v in base_budget.items()}
    worst, over = 0.0, 0
    for a, b in zip(patches, scaled):
        key = 'acral' if a['region'] in acral else 'default'
        r0 = a['area_m2'] / base_budget[key]
        r1 = b['area_m2'] / scaled_budget[key]
        worst = max(worst, abs(r1 / r0 - 1.0))
        over += int(r1 > 1.0)
    record('every patch sits at the same fraction of its own scaled budget',
           worst, 1e-12,
           'the bisection rule stops when a patch is under its area budget. '
           'Scale the budget as s**%g and every patch sits at exactly the '
           'fraction it did, so the rule re-run on the scaled body returns this '
           'partition and not a finer one. %d patches over budget, before and '
           'after alike.' % (bs.exponent('area'), over), patches=len(patches))

    # 6. Patch route lengths are conduction routes and take the length exponent,
    #    so this stage and the nerve stage agree by construction rather than by
    #    coincidence -- the same quantity, from the same table.
    worst = max(abs(b['path_length_m'] / (a['path_length_m'] * scale) - 1.0)
                for a, b in zip(patches, scaled) if a.get('path_length_m'))
    record('patch afferent routes take the same exponent as the nerve routes',
           worst, 1e-12,
           'route_length, s**%g, the same table entry scale_nerve_conduction.py '
           'uses; a patch afferent is a conduction route like any other'
           % bs.exponent('route_length'))

    # 7. Completeness over every numeric leaf, as in stage 1.
    numeric = set()
    for patch in patches:
        for key, value in patch.items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                numeric.add(key)
            elif isinstance(value, list) and value and all(
                    isinstance(x, (int, float)) and not isinstance(x, bool) for x in value):
                numeric.add(key)
    unhandled = numeric - set(PATCH_FIELDS) - set(PATCH_INVARIANT)
    record('every numeric patch field is scaled or declared invariant',
           float(len(unhandled)), 0.0, 'unhandled: %r' % (sorted(unhandled),),
           numeric_fields=sorted(numeric))
    return gates


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stature-m', type=float)
    group.add_argument('--scale', type=float)
    parser.add_argument('--output', default=None)
    parser.add_argument('--sabotage', choices=('area-as-length',), default=None)
    args = parser.parse_args()

    stature = args.stature_m if args.stature_m is not None else args.scale * MECHANICAL_STATURE_M
    scale = resolve({'stature_m': stature})['derived']['stature_scale']
    out = ROOT / (args.output or ('data/derived/body-variants/stature_%s'
                                  % ('%.6f' % scale).replace('.', 'p')))
    out.mkdir(parents=True, exist_ok=True)

    if args.sabotage == 'area-as-length':
        PATCH_FIELDS['area_m2'] = 'length'

    scaled, report = build(ROOT, scale)
    report['stature_m'] = stature
    report['sources'] = {p: _sha(ROOT / p) for p in (DERMATOMES, SKIN_GEOMETRY, EXTERIOR)}
    report['script_sha256'] = _sha(Path(__file__))
    (out / 'dermatomes_scaled.json').write_text(json.dumps(scaled, indent=1) + '\n')
    (out / 'skin_patch_report.json').write_text(json.dumps(report, indent=2) + '\n')

    for gate in report['gates']:
        print('%-5s %-62s %.3g (tol %.3g)'
              % ('PASS' if gate['passed'] else 'FAIL', gate['gate'],
                 gate['relative_error'], gate['tolerance']))
    d = report['decision']
    print('\npatches           %d, at every stature' % d['patch_count'])
    print('exterior          %.4f m2 -> %.4f m2  (s**%g)'
          % (report['base_exterior_area_m2'], report['scaled_exterior_area_m2'],
             d['area_exponent']))
    print('density           %.1f -> %.1f patches/m2  (%+.1f%%, s**%g)'
          % (d['base_density_per_m2'], d['scaled_density_per_m2'],
             d['density_change_percent'], d['receptor_density_exponent']))
    print('not chosen        count_follows_area would give %.1f patches'
          % d['alternative_patch_count'])
    print('\nwrote', out)
    if not all(g['passed'] for g in report['gates']):
        raise SystemExit('scaled skin patches failed their own gates')


if __name__ == '__main__':
    main()
