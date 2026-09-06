"""Exact retained 80-muscle source catalog; cortical assignments are priors."""
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
SOURCE_MODEL='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'

def native_muscle_catalog(root):
    path=Path(root)/SOURCE_MODEL
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    rows=[]
    for muscle in ET.parse(path).getroot().iter('Millard2012EquilibriumMuscle'):
        name=muscle.attrib['name'];side=name[-1]
        if side not in ('r','l') or name[-2]!='_':raise ValueError('Unresolved source laterality')
        bodies=sorted({p.text.rsplit('/',1)[-1] for p in muscle.findall('.//PathPointSet/objects/*/socket_parent_frame')})
        if not bodies:raise ValueError('Missing source attachment bodies')
        distal=any(b.startswith(('calcn','talus','toes')) for b in bodies)
        knee=any(b.startswith(('tibia','patella')) for b in bodies)
        group='ankle_foot' if distal else 'knee' if knee else 'hip_pelvis'
        hemi='lh' if side=='r' else 'rh'
        rows.append({'id':name,'side':side,'body_group':group,'attachment_bodies':bodies,
            'sensory_region':f'brain-{hemi}-postcentral','motor_region':f'brain-{hemi}-precentral',
            'max_isometric_force_n':float(muscle.findtext('max_isometric_force')),
            'optimal_fiber_length_m':float(muscle.findtext('optimal_fiber_length')),
            'source_path':SOURCE_MODEL,'source_sha256':digest,
            'assignment_basis':'contralateral regional cortical assignment engineering prior; body group from source attachment names; no measured somatotopic recruitment'})
    if not rows or len({m['id'] for m in rows})!=len(rows):raise ValueError('Invalid source muscle catalog')
    return rows
