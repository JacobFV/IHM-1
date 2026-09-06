"""Write a fresh isolated native source with an explicit pre-init hair partition."""
from pathlib import Path
import argparse,json,sys,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256
from ihm.assembly.hair_source_ownership import partition_represented_guides
from ihm.native.mechanical_stream import SOURCE_FILES


def materialize(destination,root=ROOT):
    root=Path(root).resolve();destination=Path(destination).resolve();dependencies={}
    if destination.exists() or not destination.is_relative_to(root):raise ValueError('Fresh isolated owned destination required')
    def read(path,expected=None):
        path=root/path;raw=path.read_bytes();digest=file_sha256(raw)
        if expected is not None and digest!=expected:raise ValueError('Input hash mismatch: '+str(path))
        dependencies[str(path.relative_to(root))]=digest;return raw
    partition=json.loads(read('data/research/hair_native_partition.json'));base=json.loads(read('data/research/hair_source_factory.json'))
    for path,digest in partition['source_dependencies'].items():read(path,digest)
    native=base['actual_frozen_native_inertia'];prepared={k:v['prepared'] for k,v in partition['groups'].items()}
    recomputed=partition_represented_guides(prepared,native['bodies'])
    if recomputed!=partition['partition'] or recomputed['represented_guide_count']!=94:raise ValueError('Partition receipt mismatch')
    selected={r['owner']:r['native_after'] for r in recomputed['partitions']};bodies={b['owner']:b for b in native['bodies']}
    raw=read(base['native_source_inertia']['model_path'],base['native_source_inertia']['model_sha256']);tree=ET.fromstring(raw);xml_bodies=tree.findall('.//BodySet/objects/Body')
    if [b.attrib['name'] for b in xml_bodies]!=list(bodies):raise ValueError('Source native body topology mismatch')
    before={};after={}
    for b in xml_bodies:
        name=b.attrib['name'];p=selected.get(name,bodies[name]);mass=float(p['mass_kg']);center=p['centroid_m'];tensor=np.asarray(p['inertia_com_kg_m2'])
        if mass<=0 or np.linalg.eigvalsh(.5*np.trace(tensor)*np.eye(3)-tensor).min()<=0:raise ValueError('Invalid residual native inertia')
        before[name]={k:b.findtext(k) for k in ('mass','mass_center','inertia')}
        b.find('mass').text=format(mass,'.17g');b.find('mass_center').text=' '.join(format(x,'.17g') for x in center)
        b.find('inertia').text=' '.join(format(tensor[i,j],'.17g') for i,j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2)))
        after[name]={'mass_kg':mass,'centroid_m':center,'inertia_com_kg_m2':tensor.tolist()}
    # Nonzero engineering initialization tests rigid velocity matching without an
    # integrated interval or pose change. These are not held human observations.
    speeds={'pelvis_tx':.013,'pelvis_rotation':.007}
    for name,speed in speeds.items():
        coordinate=tree.find(f'.//Coordinate[@name="{name}"]')
        if coordinate is None or coordinate.find('default_speed_value') is None:raise ValueError('Expected free pelvis coordinate missing')
        coordinate.find('default_speed_value').text=str(speed)
    files={};original=Path('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking')
    for name in SOURCE_FILES:files[name]=ET.tostring(tree,encoding='utf-8',xml_declaration=True) if name=='subject_walk_scaled.osim' else read(original/name)
    read(Path(__file__).relative_to(root));destination.mkdir(parents=True)
    for name,data in files.items():(destination/name).write_bytes(data)
    target=sum(b['mass_kg'] for b in after.values())
    result={'schema':'ihm.hair-residual-native-input.v1','native_started':False,'live_enabled':False,'destination':str(destination.relative_to(root)),
        'dependencies':dependencies,'files':{name:file_sha256(raw) for name,raw in files.items()},'target_native_mass_kg':target,'original_actual_mass_kg':native['effective_mass_kg'],
        'represented_hair_mass_kg':recomputed['represented_mass_kg'],'mass_scale_required':1.,'mass_transfer_port_enabled':False,'expected_native_bodies':after,
        'source_xml_before':before,'engineering_default_speeds':speeds,'partition_conservation':recomputed['conservation'],
        'installation_order':'Explicit actual residual mass/COM/full tensor XML properties before model load/initSystem/muscle reference; target equals residual XML total; no uniform restoration to full subject mass',
        'limitations':partition['limits']}
    (destination/'hair_residual_identity.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('destination',type=Path);args=parser.parse_args();report=materialize(args.destination)
    print(json.dumps({k:report[k] for k in ('destination','target_native_mass_kg','represented_hair_mass_kg','native_started')}))
