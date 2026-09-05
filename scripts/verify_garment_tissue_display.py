"""Boundary display must preserve computed positions, contact and source ownership."""
from pathlib import Path
import sys,json,gzip,tempfile,argparse
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

def verify(experiment):
    try:from scripts.export_garment_tissue_display import export_display
    except ImportError:raise AssertionError('Verified mechanics display export missing') from None
    out=Path(tempfile.mkdtemp(prefix='garment-display-check-',dir='data/derived'))/'display'
    receipt=export_display(experiment,out)
    d=json.loads(gzip.decompress((out/'display.json.gz').read_bytes()))
    raw=np.load(Path(experiment)/'trajectory.npz')
    ids=np.array(d['tissue']['source_volume_node_indices'])
    assert np.max(np.abs(np.array(d['tissue']['positions_m'])-raw['tissue_positions_m'][:,ids]))<1e-7
    assert len(d['tissue']['triangle_material_index'])==len(d['tissue']['triangles'])
    assert set(d['replaced_body_entity_ids'])=={'body-bp3d-FJ3132','body-bp3d-FJ3133','body-bp3d-FJ3134'}
    assert np.max(np.abs(np.array(d['panel']['positions_m'])-raw['panel_positions_m']))<1e-7
    force_a=np.array(d['tissue']['contact_force_n']).sum(axis=1)
    force_b=np.array(d['panel']['contact_force_n']).sum(axis=1)
    assert np.max(np.linalg.norm(force_a+force_b,axis=1))<1e-7
    assert len(d['time_s'])==len(d['frames'])+1 and d['time_s'][0]==0
    assert d['report']['whole_garment_containment_validated'] is False
    assert receipt['position_quantization_max_error_m']<1e-7
    try:export_display(experiment,out)
    except ValueError:pass
    else:raise AssertionError('Display export overwrote a retained artifact')
    report={'status':'passed','output_dir':str(out),'receipt':receipt}
    (out.parent/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('experiment',type=Path);verify(p.parse_args().experiment)
