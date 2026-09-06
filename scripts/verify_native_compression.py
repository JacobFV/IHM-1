"""Actual optional external-pressure topology: zero-boundary parity and drainage."""
from pathlib import Path
import json
import tempfile
from ihm.native.coupled_session import CoupledNativeSession
from ihm.native.session import SessionConfig


def main():
    root=Path(__file__).resolve().parents[1]
    state=root/'data/derived/canonical/native_baseline_v1/states/native_stabilized.xml'
    config=SessionConfig(state_path=state,engine_variant='whole_body_integrity_depletion',horizon_s=.4)
    output=Path(tempfile.mkdtemp(prefix='native-skin-compression-',dir=root/'data/derived/audits'))
    results={}
    for label,pressure in [('unchanged',None),('zero',0.),('loaded',133.322387415)]:
        with CoupledNativeSession(config,output/label) as session:
            if pressure is not None:session.skin_compression(pressure)
            results[label]=[session.step(.02) for _ in range(10)]
            if pressure is not None:session.skin_compression(0)
            results[label].extend(session.step(.02) for _ in range(10))
    residuals={}
    for original,zero in zip(results['unchanged'],results['zero']):
        for key,value in original['values'].items():
            other=zero['values'].get(key)
            if value is None or other is None or key.startswith('coupling.'):continue
            error=abs(other-value);residuals[key]=max(residuals.get(key,0),error)
            assert error<=1e-8+1e-9*max(abs(value),abs(other)),(key,value,other)
    pressure_key='tissue.node.SkinE3.pressure_mmhg'
    flow_key='tissue.path.SkinE3ToSkinL1.flow_ml_per_s'
    pressure_change=max(abs(a['values'][pressure_key]-b['values'][pressure_key]) for a,b in zip(results['loaded'],results['zero']))
    flow_change=max(abs(a['values'][flow_key]-b['values'][flow_key]) for a,b in zip(results['loaded'],results['zero']))
    assert pressure_change>1e-6 and flow_change>1e-10
    assert all(f['values']['tissue.compression.Skin.requested_pa']==0 for f in results['loaded'][10:])
    record={'passed':True,'directory':str(output.relative_to(root)),'duration_s':.4,
        'maximum_zero_boundary_absolute_residual':max(residuals.values()),
        'skin_interstitial_pressure_change_mmhg':pressure_change,'skin_lymph_flow_change_ml_per_s':flow_change,
        'scope':'Uniform native whole-skin extracellular pressure boundary; intracellular reference unchanged. Not local contact or long-horizon validation.'}
    (output/'verification.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))

if __name__=='__main__':main()
