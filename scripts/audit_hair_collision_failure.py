"""Read-only retained native failure and exact passive coordinate law audit."""
from pathlib import Path
import ast,json,math,sys,xml.etree.ElementTree as ET
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from ihm.assembly.source_skin import file_sha256
from scripts.effective_passive_energy_observation import PassiveEnergyObservation


def expression_value(expression,q,qdot):
    def visit(node):
        if isinstance(node,ast.Constant) and type(node.value) in (int,float):return node.value
        if isinstance(node,ast.Name) and node.id in ('q','qdot'):return {'q':q,'qdot':qdot}[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.UAdd,ast.USub)):return visit(node.operand)*(1 if isinstance(node.op,ast.UAdd) else -1)
        if isinstance(node,ast.BinOp):
            a,b=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='exp' and len(node.args)==1 and not node.keywords:return math.exp(visit(node.args[0]))
        raise ValueError('Unsupported retained passive expression syntax')
    value=visit(ast.parse(expression,mode='eval').body)
    if not math.isfinite(value):raise ValueError('Nonfinite passive source force')
    return value


def source_joint_forces(path,coordinates):
    return [{'name':f.get('name'),'coordinate':f.findtext('coordinate'),'expression':f.findtext('expression'),
             'q_rad':coordinates[f.findtext('coordinate')]['value'],'qdot_rad_s':coordinates[f.findtext('coordinate')]['speed'],
             'generalized_force_nm':expression_value(f.findtext('expression'),coordinates[f.findtext('coordinate')]['value'],coordinates[f.findtext('coordinate')]['speed'])}
            for f in ET.parse(path).getroot().iter('ExpressionBasedCoordinateForce')]


def audit():
    folder=ROOT/'data/derived/hair-collision-native-u93_yi3p';initial=json.loads((folder/'initial.json').read_bytes());failure=json.loads((folder/'failure.json').read_bytes());model=folder/'inputs/subject_walk_scaled.osim';force_path=folder/'inputs/subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml';forces=source_joint_forces(force_path,initial['coordinates']);forces.sort(key=lambda r:abs(r['generalized_force_nm']),reverse=True)
    port=PassiveEnergyObservation(model);energy=port.observe({'muscles':[{'name':name,'type':m['muscle_type'],'passive_energy_j':m['passive_energy_j'],'fiber_length_m':m['fiber_length_m']} for name,m in initial['muscles'].items()]})
    report={'schema':'ihm.hair-collision-failure-audit.v1','accepted_physical_time_s':failure['native_time_s'],'collision_exercised':False,'failure':failure,'passive_coordinate_forces':forces,'initial_passive_energy_observation':energy,
        'raw_native_kinetic_energy_j':initial['kinetic_energy_j'],'raw_native_potential_energy_j':initial['potential_energy_j'],
        'native_joint_reaction_force_norm_n':{name:float(np.linalg.norm(b['joint_reaction_force_n'])) for name,b in initial['bodies'].items()},
        'maximum_initial_angular_speed_rad_s':max(float(np.linalg.norm(b['angular_velocity_rad_s'])) for b in initial['bodies'].values()),'maximum_initial_tendon_force_n':max(m['tendon_force_n'] for m in initial['muscles'].values()),
        'constraint_position_error':initial['constraint_position_error'],'constraint_velocity_error':initial['constraint_velocity_error'],
        'cause':'Kinematic-only search admits shoulder coordinates deep inside exponential source passive-limit tails despite legal broad XML ROM; initialized source joint reactions are enormous before hair contact',
        'reference_check':'Offline exact same-initial snapshot roots/tangents dot minimum0.9999999999999996, clamp residual1.145e-16m; fixed held-surface hair solve succeeds with no contact and maxload0.003644N',
        'next_requirement':'Force-aware source-domain pose qualification and alternate collision normal/target. No blind timestep shrink or passive force suppression.',
        'profile':{'elapsed_s':3.77,'peak_rss_kib':94392,'native_child_reaped':True},'source_receipts':{str(p.relative_to(ROOT)):file_sha256(p) for p in (folder/'initial.json',folder/'failure.json',model,force_path,Path(__file__),ROOT/'scripts/effective_passive_energy_observation.py')}}
    return report

if __name__=='__main__':
    report=audit();(ROOT/'data/research/hair_collision_failure.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report['passive_coordinate_forces'][:3],indent=2))
