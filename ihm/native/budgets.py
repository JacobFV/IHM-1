"""Source-state fluid inventory audits without double counting circuit views."""
from pathlib import Path
import json
from .circuits import native_circuit_graph


def fluid_inventory(state_path):
    graph=native_circuit_graph(state_path)
    circuit=next(c for c in graph['circuits'] if c['name']=='FullCardiovascular')
    selected=set(circuit['nodes']);volumes={};boundaries=[]
    for node in graph['nodes']:
        if node['id'] not in selected:continue
        item=node['properties'].get('Volume')
        if item is None:continue
        if item['status']=='infinite_boundary':boundaries.append(node['id']);continue
        if item['si_value'] is None or item['si_unit']!='m^3':raise ValueError('volume lacks valid native units')
        volumes[node['id']]=item['si_value']
    system=next(s for s in graph['systems'] if s['type']=='BioGearsGastrointestinalSystemData')
    water=system['properties']['StomachContents/Water']
    if water['si_unit']!='m^3' or water['si_value'] is None:raise ValueError('stomach water lacks valid native units')
    return dict(source=graph['source'],simulation_time=graph['simulation_time'],circuit_volumes_m3=volumes,
        stomach_water_m3=water['si_value'],infinite_boundaries=boundaries,
        combined_finite_volume_m3=sum(volumes.values())+water['si_value'])


def audit_fluid_budget(run_dir):
    run_dir=Path(run_dir);a=fluid_inventory(run_dir/'states/native_stabilized.xml');b=fluid_inventory(run_dir/'states/native_final.xml')
    if set(a['circuit_volumes_m3'])!=set(b['circuit_volumes_m3']):raise ValueError('circuit inventory changed; explicit remapping required')
    changes={name:b['circuit_volumes_m3'][name]-volume for name,volume in a['circuit_volumes_m3'].items()}
    stomach_change=b['stomach_water_m3']-a['stomach_water_m3']
    residual=sum(changes.values())+stomach_change
    return dict(schema_version=1,kind='native_finite_volume_inventory',initial=a,final=b,
        node_volume_changes_m3=changes,circuit_volume_change_m3=sum(changes.values()),stomach_water_change_m3=stomach_change,
        combined_volume_change_m3=residual,relative_combined_change=residual/a['combined_finite_volume_m3'],
        limitations=['FullCardiovascular selected once; globally shared circuit nodes are counted once. Alternative circuit views and parent compartments are not added.',
            'Stomach water is an explicit reservoir outside these circuit nodes; its transfer can change blood and interstitial volumes without adding external water.',
            'An inventory difference is not automatically a numerical error. External fluids, hemorrhage, urination and other sources/sinks require their own cumulative flux accounting.',
            'Finite volume accounting is not chemical mass, energy or independent patient validation.'])
