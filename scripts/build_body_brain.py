"""Assemble IBM's concrete brain law and cached DK substrate into the BP3D body.

Runtime needs only NumPy; asset ingestion uses the explicitly selected IBM Python
with nibabel. The source workspace/cache are read-only; every used byte is copied.
"""
from __future__ import annotations
import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / 'data/derived/canonical'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def build(ibm_root, asset_root=None, asset_python=None):
    ibm_root = Path(ibm_root).expanduser().resolve()
    out = OUT
    out.mkdir(parents=True, exist_ok=True)
    preserved = out / 'brain-sources'
    preserved.mkdir(exist_ok=True)
    sources = []

    def preserve(path, name, role):
        path = Path(path)
        target = preserved / name
        shutil.copyfile(path, target)
        sources.append({'source_path': str(path), 'preserved_path': str(target.relative_to(ROOT)),
                        'sha256': digest(path), 'bytes': path.stat().st_size, 'role': role})
        return target

    for rel, name, role in [
        ('ibm/processes/neural.py', 'ibm-neural.py', 'executed_neural_law'),
        ('ibm/anatomy/systems.py', 'ibm-anatomy-systems.py', 'enumerated_DK_parcellation'),
        ('ibm/anatomy/sources.py', 'ibm-anatomy-sources.py', 'anatomical_source_scope'),
        ('ibm/topologies/association.py', 'ibm-association.py', 'distance_prior_definition'),
        ('ibm/materialize/geometry.py', 'ibm-geometry.py', 'source_geometry_loading_reference'),
        ('data/sources/mne-sample/raw/.location.yaml', 'mne-sample-location.yaml', 'asset_location'),
        ('data/sources/mne-sample/card.yaml', 'mne-sample-card.yaml', 'asset_source_card'),
    ]:
        preserve(ibm_root / rel, name, role)
    if asset_root is None:
        location = (preserved / 'mne-sample-location.yaml').read_text()
        local = next(line.split(':',1)[1].strip().strip('"\'') for line in location.splitlines() if line.startswith('local_root:'))
        asset_root = Path(local) / 'processed-v6/MNE-sample-data/subjects/fsaverage'
    asset_root = Path(asset_root).expanduser().resolve()
    for hemi in ['lh','rh']:
        preserve(asset_root / 'surf' / f'{hemi}.pial', f'{hemi}.pial', 'original_full_cortical_surface')
        preserve(asset_root / 'label' / f'{hemi}.aparc.annot', f'{hemi}.aparc.annot', 'original_DK_annotation')
    # Invoke the declared asset environment without sys.path/env changes to IBM.
    interpreter = Path(asset_python or ibm_root / '.venv/bin/python').absolute()
    extraction = '''import json,sys,numpy as np,nibabel as nib
from pathlib import Path
base=Path(sys.argv[1]); out=Path(sys.argv[2]); arrays={}; info={"nibabel_version":nib.__version__,"numpy_version":np.__version__}
for h in ["lh","rh"]:
 v,f=nib.freesurfer.read_geometry(base/(h+".pial"))
 labels,table,names=nib.freesurfer.read_annot(base/(h+".aparc.annot"))
 arrays[h+"_vertices"]=v;arrays[h+"_faces"]=f;arrays[h+"_labels"]=labels
 info[h+"_names"]=[n.decode() for n in names]
np.savez_compressed(out,**arrays)
print(json.dumps(info))
'''
    with tempfile.TemporaryDirectory(prefix='ihm-brain-ingest-') as directory:
        archive = Path(directory) / 'source.npz'
        result = subprocess.run([str(interpreter), '-c', extraction, str(preserved), str(archive)],
                                check=True, capture_output=True, text=True)
        asset_info = json.loads(result.stdout)
        loaded = np.load(archive)
        arrays = {key: loaded[key].copy() for key in loaded.files}
    asset_info['python'] = str(interpreter)
    asset_info['python_sha256'] = digest(interpreter)
    systems = ast.parse((preserved / 'ibm-anatomy-systems.py').read_text())
    dk = next(ast.literal_eval(n.value) for n in systems.body if isinstance(n, ast.Assign)
              and any(isinstance(t, ast.Name) and t.id == 'DK_GYRI' for t in n.targets))
    manifest = json.loads((ROOT / 'data/derived/app/manifest.json').read_text())
    bp_model = next(m for m in manifest['models'] if m['id'] == 'bodyparts3d')
    bp = [s for s in manifest['structures'] if s['model_id']=='bodyparts3d']
    cortical = [s for s in bp if s['system']=='nervous' and any(
        term in s['name'] for term in ('gyrus', 'lobule', 'occipital lobe', 'insula'))]
    tissue = {}
    for s in cortical + [s for s in bp if s['name'] in
            ('medulla oblongata','pons','midbrain','hypothalamus','cerebellum','left thalamus','right thalamus')]:
        path = ROOT / 'data/derived/app/geometry' / (s['id'] + '.json.gz')
        mesh = json.loads(gzip.decompress(path.read_bytes()))
        positions = np.asarray(mesh['positions']).reshape(-1,3)
        tissue[s['id']] = {'structure': s, 'center': positions.mean(axis=0),
            'min': positions.min(axis=0), 'max': positions.max(axis=0), 'sha256': digest(path)}
    # Homologous named gyri provide observable fit residuals, not an invented accuracy.
    matches = {'fusiform':'fusiform gyrus','inferiortemporal':'inferior temporal gyrus',
               'middletemporal':'middle temporal gyrus','parahippocampal':'parahippocampal gyrus',
               'postcentral':'postcentral gyrus','precentral':'precentral gyrus',
               'superiorfrontal':'superior frontal gyrus','superiorparietal':'superior parietal lobule',
               'supramarginal':'supramarginal gyrus','insula':'insula'}
    # FreeSurfer surface tkRAS millimeters: right, anterior, superior.
    rotation = np.array([[-1.,0.,0.],[0.,0.,1.],[0.,1.,0.]])
    source_nodes = []
    landmark_source, landmark_target, landmark_ids = [], [], []
    for hemi, hemisphere in [('lh','left'),('rh','right')]:
        names = asset_info[hemi+'_names']
        for label in dk:
            index = names.index(label)
            selected = arrays[hemi+'_labels']==index
            xyz = arrays[hemi+'_vertices'][selected].mean(axis=0)
            node = {'id':f'brain-{hemi}-{label}', 'kind':'cortical_population', 'name':f'{hemisphere} {label}',
                    'hemisphere':hemisphere, 'ibm_partition':{'system':'cortical_areas','label':label},
                    'source_position_tkras_mm':xyz.tolist(), 'source_vertex_count':int(selected.sum())}
            source_nodes.append(node)
            if label in matches:
                target_name = hemisphere + ' ' + matches[label]
                target = next(t for t in tissue.values() if t['structure']['name']==target_name)
                landmark_source.append(xyz @ rotation.T * .001)
                landmark_target.append(target['center'])
                landmark_ids.append({'brain_node_id':node['id'],'body_entity_id':'body-'+target['structure']['id']})
    x, y = np.array(landmark_source), np.array(landmark_target)
    scale = np.sum((x-x.mean(0))*(y-y.mean(0)),axis=0)/np.sum((x-x.mean(0))**2,axis=0)
    if np.any(scale <= 0):
        raise ValueError('Registration would reverse a declared anatomical axis')
    translation = y.mean(0)-x.mean(0)*scale
    linear = np.diag(scale) @ rotation * .001
    residuals = np.linalg.norm(x*scale+translation-y,axis=1)
    matrix = np.eye(4);matrix[:3,:3]=linear;matrix[:3,3]=translation
    for node in source_nodes:
        position = np.array(node['source_position_tkras_mm']) @ linear.T + translation
        candidates = [t for t in tissue.values() if t['structure'] in cortical and
                      t['structure']['name'].startswith(node['hemisphere']+' ')]
        preferred = matches.get(node['ibm_partition']['label'])
        if preferred:
            target = next(t for t in candidates if t['structure']['name']==node['hemisphere']+' '+preferred)
            correspondence = 'homologous_gyral_name'
        else:
            target = min(candidates,key=lambda t:np.linalg.norm(t['center']-position))
            correspondence = 'nearest_coarse_cortical_support'
        node.update(position_m=position.tolist(), body_entity_id='body-'+target['structure']['id'],
                    tissue_support_name=target['structure']['name'], support_transfer=correspondence,
                    support_centroid_distance_m=float(np.linalg.norm(position-target['center'])),
                    state_scope='IBM cortical adaptive population; atlas label transferred between templates')
    nodes = list(source_nodes)
    for t in tissue.values():
        s=t['structure']
        if s in cortical:
            continue
        role = 'medulla' if s['name']=='medulla oblongata' else s['name']
        nodes.append({'id':'brain-'+s['id'], 'kind':'subcortical_population_proxy', 'name':s['name'],
                      'position_m':t['center'].tolist(), 'body_entity_id':'body-'+s['id'],
                      'autonomic_role':role, 'state_scope':'Reduced cortical-law proxy on actual subcortical tissue; nucleus-specific dynamics uncalibrated'})
    points = np.array([n['position_m'] for n in nodes])
    distances = np.linalg.norm(points[:,None]-points[None,:],axis=-1)
    edges=[]
    for i,node in enumerate(nodes):
        for j in np.argsort(distances[i])[1:9]:
            edges.append({'source':node['id'],'target':nodes[j]['id'],
                          'distance_m':float(distances[i,j]),'weight':float(np.exp(-distances[i,j]/.04)),
                          'evidence_kind':'geometric_prior','conduction':'instantaneous regional coupling; no tract-delay model'})
    surface_arrays = {}
    for h in ['lh','rh']:
        surface_arrays[h+'_positions_m'] = arrays[h+'_vertices'] @ linear.T + translation
        surface_arrays[h+'_faces'] = arrays[h+'_faces']
        surface_arrays[h+'_DK_label_indices'] = arrays[h+'_labels']
    np.savez_compressed(out/'brain-registered-surfaces.npz', **surface_arrays)
    peripheral = [s for s in bp if s['system']=='nervous' and 'nerve' in s['name']]
    medulla_ids = [n['id'] for n in nodes if n.get('autonomic_role')=='medulla']
    nerve_connections = []
    for s in peripheral:
        side = 'left' if 'left' in s['name'] else 'right' if 'right' in s['name'] else None
        if 'optic nerve' in s['name']:
            target_name = f'{side} thalamus' if side else 'left thalamus'
            route_scope = 'visual relay represented by coarse thalamic proxy'
        elif any(term in s['name'] for term in ('oculomotor', 'trochlear')):
            target_name = 'midbrain'
            route_scope = 'ocular motor midbrain interface prior'
        else:
            target_name = 'pons'
            route_scope = 'ophthalmic/trigeminal branch represented by coarse pontine interface'
        origins = [n for n in nodes if n['name'] == target_name]
        if side and len(origins) > 1:
            origin = max(origins, key=lambda n: n['position_m'][0]) if side == 'left' else min(origins, key=lambda n: n['position_m'][0])
        else:
            origin = origins[0]
        nerve_connections.append({'brain_node_id':origin['id'], 'body_entity_id':'body-'+s['id'],
            'name':s['name'], 'relation':'cranial_neural_interface', 'scope':route_scope,
            'evidence_kind':'anatomical_interface_prior', 'routed_axons':False,
            'dynamic_feedback_applied':False})
    parameters = {'ibm_wilson_cowan':{'tau_membrane_s':.015,'v_rest_mv':-65.,'v_half_mv':-55.,
                  'slope_mv':4.,'r_max_hz':100.,'tau_adaptation_s':.5,'adaptation_gain':.05,'drive_gain':1.},
                  'body_transfer_priors':{'tonic_drive_nS':4.,'recurrent_drive_nS_per_Hz':.08,
                  'inhibitory_conductance_per_Hz':.001,'map_reference_mmHg':90.,'temperature_Q10':2.}}
    data = {'schema_version':1, 'id':'canonical-ibm-brain', 'frame':'bodyparts3d-display-m',
            'axes':{'x':'left','y':'superior','z':'anterior'}, 'units':'m', 'nodes':nodes,'edges':edges,
            'nerve_connections':nerve_connections,'parameters':parameters,'sources':sources,
            'source_model_identity':{'module':'ibm.processes.neural','implementation':'wilson_cowan_adaptive',
                'executed_symbols':['_sigmoid','wilson_cowan_excitatory','shunting_inhibition_rate'],
                'neural_source_sha256':digest(preserved/'ibm-neural.py')},
            'parameter_provenance':{'ibm_wilson_cowan':'exact defaults of preserved IBM law; literature/weak priors, not fitted human measurements',
                                    'body_transfer_priors':'IHM synthesis priors; uncalibrated body coupling'},
            'registration':{'source_frame':'fsaverage-surface-tkRAS-mm','target_frame':'bodyparts3d-display-m',
                'matrix_4x4':matrix.tolist(),'method':'axis-preserving affine scale/translation fitted to 20 named gyral centroids',
                'landmarks':[dict(pair,residual_m=float(r)) for pair,r in zip(landmark_ids,residuals)],
                'landmark_rms_m':float(np.sqrt(np.mean(residuals**2))), 'maximum_landmark_residual_m':float(residuals.max()),
                'residual_scope':'in-sample inter-template centroid residual; not held-out anatomical accuracy',
                'bp3d_source_display_transform':bp_model['display_transform'],'same_individual':False,
                'full_source_surface_archive':'data/derived/canonical/brain-registered-surfaces.npz',
                'surface_archive_sha256':digest(out/'brain-registered-surfaces.npz'),
                'surface_role':'registered source evidence; canonical body retains BP cortex geometry without duplicate organs'},
            'body_tissue_sources':{key:{'sha256':value['sha256'],'centroid_m':value['center'].tolist()} for key,value in tissue.items()},
            'asset_ingestion':asset_info,
            'peripheral_coverage':'Available BP3D eye/orbital nerves only; spinal cord and systemic peripheral/autonomic nerve trajectories are not present in this source inventory.',
            'ports':{'afferent':[{'id':key,'unit':unit} for key,unit in [('mean_arterial_pressure_mmHg','mmHg'),('oxygen_saturation','fraction'),('core_temperature_C','C')]],
                     'efferent':[{'id':'autonomic-sympathetic','origin_brain_nodes':medulla_ids,'signal':'sympathetic_fraction','applied_to_body':False},
                                  {'id':'autonomic-parasympathetic','origin_brain_nodes':medulla_ids,'signal':'parasympathetic_fraction','applied_to_body':False}]},
            'assumptions':[
                {'id':'brain-registration','kind':'cross_source_registration','detail':'fsaverage and BP3D are different reference anatomies. Gyral centroid fit transfers cortical partitions; residuals quantify only fitted centroid mismatch.'},
                {'id':'brain-topology','kind':'source_prior','detail':'8 nearest regional neighbors with IBM 40 mm exponential distance prior. Edges are not measured tracts; instantaneous rate coupling omits axonal delays.'},
                {'id':'brain-state','kind':'source_model','detail':'Exact IBM Wilson-Cowan adaptive and shunting law ASTs execute from preserved, verified original source. IHM RK4 integrates coupled regional state; IBM spectral full-model runtime is not claimed.'},
                {'id':'brain-physiology','kind':'synthesis_prior','detail':'MAP/oxygen availability scales synaptic drive; temperature Q10 scales kinetics. These transfers are generic priors and not calibrated brain perfusion/metabolic models.'},
                {'id':'brain-autonomic','kind':'synthesis_prior','detail':'Medulla rate and physiology stress produce bounded autonomic commands. Commands are returned, never silently applied to native engines; no validated closed loop is claimed.'},
                {'id':'brain-subcortex','kind':'synthesis_prior','detail':'Subcortical source geometry supports coarse activity proxies using the cortical law; nucleus-specific neural dynamics are not supplied by this integration.'},
                {'id':'brain-nerves','kind':'anatomical_interface_prior','detail':'Cranial nerve tissue interfaces identify canonical endpoints; axonal routing, nucleus specificity, and peripheral conduction are not measured or simulated here.'}],
            'calibration':{'source_fidelity':'hash-preserved IBM functions, source annotations and full surfaces',
                           'anatomical_fit':'generic inter-template label transfer', 'physiological_validation':'not empirically calibrated',
                           'closed_loop_validation':False}}
    dump(out/'brain.json',data)
    from ihm.assembly.brain import BodyBrain
    brain=BodyBrain.from_dict(data)
    baseline=brain.step(.5, {})
    low=brain.step(.5, {'mean_arterial_pressure_mmHg':45.,'oxygen_saturation':.65,'core_temperature_C':37.})
    dump(out/'brain-state.json',{'baseline':baseline,'hypoperfusion_hypoxia':low,
        'scope':'executed synthetic perturbation after baseline; not native engine feedback or empirical validation'})
    print(json.dumps({'nodes':len(nodes),'cortical_regions':len(source_nodes),'edges':len(edges),
                      'landmark_rms_m':data['registration']['landmark_rms_m'],'preserved_sources':len(sources)}))
    return data


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ibm-root',default=str(Path.home()/'Documents/IBM-1'))
    parser.add_argument('--asset-root',help='Explicit fsaverage subject directory; default follows IBM mne-sample location')
    parser.add_argument('--asset-python',help='Python with nibabel for build-only asset ingestion; default IBM .venv')
    args=parser.parse_args()
    build(args.ibm_root,args.asset_root,args.asset_python)
