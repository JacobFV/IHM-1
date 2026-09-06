"""Source-only verification/feasibility of separate mattress and skin compression."""
from pathlib import Path
import hashlib,json,sys,tempfile
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.bed_compression import load_bed,series_response,maximum_approach


def main():
    manifest_path=ROOT/'data/derived/supine-surface-contact-exmzq9pq/manifest.json'
    manifest=json.loads(manifest_path.read_text());data=np.load(ROOT/manifest['arrays_path']);skin=manifest['material']
    points=data['reference_points_source_m'];areas=data['area_m2'];minimum=float(points[:,0].min());weight=77.6122029*9.81
    report=dict(passed=True,native_run=False,accepted_equilibrium=False,
        scope='Source constitutive/domain and rigid-translation feasibility; no native support or moment equilibrium acceptance',materials={})
    for material in ('SM','MM','HM'):
        bed=load_bed(ROOT,material);limit=maximum_approach(skin,bed)
        samples=np.linspace(0,limit*.99,31);result=series_response(samples,skin,bed)
        assert np.max(np.abs(result['pressure_balance_residual_pa']))<1e-7
        assert np.allclose(result['skin_indentation_m']+result['bed_indentation_m'],samples,atol=1e-12)
        assert np.max(result['skin_indentation_m'])<=.0033+1e-12
        eps=1e-7;probe=limit*.35
        derivative=(series_response(probe+eps,skin,bed)['total_energy_j_m2']-series_response(probe-eps,skin,bed)['total_energy_j_m2'])/(2*eps)
        assert np.isclose(derivative,series_response(probe,skin,bed)['pressure_pa'],rtol=1e-6)
        try:series_response(limit+1e-5,skin,bed)
        except ValueError:pass
        else:raise AssertionError('Bed/skin extrapolation accepted')
        def support(translation):
            delta=np.maximum(0,manifest['plane_source_x_m']-points[:,0]+translation)
            values=series_response(delta,skin,bed)
            return float(np.dot(areas,values['pressure_pa'])),values
        max_translation=minimum-manifest['plane_source_x_m']+limit*(1-1e-10)
        max_force,_=support(max_translation)
        result_row=dict(curve=bed,maximum_admissible_total_approach_m=limit,
            maximum_rigid_translation_m=max_translation,maximum_support_n=max_force,
            weight_n=weight,force_balance_reachable=max_force>=weight)
        if max_force>=weight:
            lo,hi=0.,max_translation
            for _ in range(50):
                middle=(lo+hi)/2
                if support(middle)[0]>weight:hi=middle
                else:lo=middle
            translation=(lo+hi)/2;force,values=support(translation);point_forces=areas*values['pressure_pa']
            cop=(points*point_forces[:,None]).sum(0)/force
            result_row.update(force_balanced_translation_m=translation,support_n=force,center_of_pressure_source_m=cop.tolist(),
                maximum_skin_indentation_m=float(values['skin_indentation_m'].max()),maximum_bed_indentation_m=float(values['bed_indentation_m'].max()),
                warning='Only net normal force solved; contact moment and all generalized residuals remain unverified')
        report['materials'][material]=result_row
    output=Path(tempfile.mkdtemp(prefix='bed-skin-series-',dir=ROOT/'data/derived'))
    report['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (manifest_path,ROOT/'data/research/bed_material/manifest.json',Path(__file__).resolve(),ROOT/'ihm/assembly/bed_compression.py')}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({**report,'output':str(output)},indent=2))


if __name__=='__main__':main()
