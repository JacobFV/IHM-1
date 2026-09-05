#!/usr/bin/env python3
"""Run native FEBio compression and check an independent analytical solution."""
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/derived/mechanics-reference"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model():
    nodes = [(x, y, z) for z in (0, .5, 1) for x, y in ((0, 0), (1, 0), (1, 1), (0, 1))]
    def node_set(name, ids):
        return f'<NodeSet name="{name}">' + ''.join(f'<node id="{i}"/>' for i in ids) + '</NodeSet>'
    return '''<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_spec version="2.5">
<Module type="solid"/>
<Control><analysis type="static"/><time_steps>10</time_steps><step_size>0.1</step_size></Control>
<Material><material id="1" name="solid" type="neo-Hookean"><E>1000</E><v>0</v></material></Material>
<Geometry><Nodes>''' + ''.join(f'<node id="{i}">{x},{y},{z}</node>' for i, (x,y,z) in enumerate(nodes, 1)) + '''</Nodes>
<Elements type="hex8" mat="1" name="block"><elem id="1">1,2,3,4,5,6,7,8</elem><elem id="2">5,6,7,8,9,10,11,12</elem></Elements>
''' + node_set('bottom', range(1,5)) + node_set('top', range(9,13)) + node_set('xzero', [i for i,(x,y,z) in enumerate(nodes,1) if x==0]) + node_set('yzero', [i for i,(x,y,z) in enumerate(nodes,1) if y==0]) + '''</Geometry>
<Boundary><fix bc="z" node_set="bottom"/><fix bc="x" node_set="xzero"/><fix bc="y" node_set="yzero"/>
<prescribe bc="z" node_set="top"><scale lc="1">-0.1</scale><relative>0</relative></prescribe></Boundary>
<LoadData><loadcurve id="1" type="linear"><point>0,0</point><point>1,1</point></loadcurve></LoadData>
<Output><logfile>
<node_data data="ux;uy;uz;Rz" file="nodes.txt" delim=",">1,2,3,4,5,6,7,8,9,10,11,12</node_data>
<element_data data="sx;sy;sz" file="stress.txt" delim=",">1,2</element_data>
</logfile><plotfile type="febio"><var type="displacement"/><var type="stress"/></plotfile></Output>
</febio_spec>
'''


def records(path):
    blocks = []
    for line in path.read_text().splitlines():
        if line.startswith('*Time'):
            blocks.append({"time": float(line.split('=')[1]), "rows": []})
        elif line and not line.startswith('*') and blocks:
            blocks[-1]['rows'].append([float(v) for v in line.split(',')])
    return blocks


def main():
    manifest = json.loads((OUT / 'build-manifest.json').read_text())
    executable = ROOT / manifest['executable']
    assert digest(executable) == manifest['executable_sha256'], 'Executable hash mismatch'
    for name, expected_hash in manifest['shared_library_sha256'].items():
        assert digest(ROOT / name) == expected_hash, f'Library hash mismatch: {name}'
    run = OUT / 'compression'
    run.mkdir(exist_ok=True)
    (run / 'compression.feb').write_text(model())
    command = [str(executable), '-i', 'compression.feb', '-nosplash']
    with (run / 'stdout.log').open('w') as log:
        result = subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT, timeout=120,
                                env={**os.environ, 'OMP_NUM_THREADS': '1'})
    report = {'backend': 'FEBio', 'version': manifest['version'], 'commit': manifest['commit'],
              'executable_sha256': manifest['executable_sha256'],
              'build_manifest_sha256': digest(OUT / 'build-manifest.json'),
              'command': command, 'returncode': result.returncode, 'status': 'failed',
              'model': 'two hex8 elements; unit cube; E=1000 Pa; nu=0; 10% axial compression; traction-free lateral faces',
              'expected_solution': 'F=diag(1,1,lambda); sigma_zz=500*(lambda^2-1)/lambda; u_z=(lambda-1)*Z',
              'tolerance': {'stress_absolute_Pa': .002, 'displacement_absolute_m': 2e-6},
              'scope': 'homogeneous material and equilibrium reference only; no contact, anatomy or local reduced-model equivalence'}
    try:
        assert result.returncode == 0, 'FEBio execution failed; inspect stdout.log'
        assert 'N O R M A L   T E R M I N A T I O N' in (run / 'compression.log').read_text(), 'No normal termination'
        stresses, nodes = records(run/'stress.txt'), records(run/'nodes.txt')
        for blocks in (stresses, nodes):
            assert len(blocks) == 11, 'Expected initial state and ten increments'
            assert all(abs(b['time']-i*.1)<1e-12 for i,b in enumerate(blocks)), 'Unexpected time sequence'
        assert stresses and abs(stresses[-1]['time']-1)<1e-12, 'Final stress time missing'
        assert nodes and abs(nodes[-1]['time']-1)<1e-12, 'Final node time missing'
        max_stress = max_displacement = 0.
        for block in stresses:
            stretch = 1-.1*block['time']
            expected = 500*(stretch**2-1)/stretch
            assert len(block['rows']) == 2
            assert [row[0] for row in block['rows']] == [1,2], 'Element IDs mismatch'
            for _,sx,sy,sz in block['rows']:
                assert all(math.isfinite(v) for v in (sx,sy,sz))
                max_stress = max(max_stress, abs(sx), abs(sy), abs(sz-expected))
        for block in nodes:
            assert len(block['rows']) == 12
            assert [row[0] for row in block['rows']] == list(range(1,13)), 'Node IDs mismatch'
            for i,ux,uy,uz,rz in block['rows']:
                z = ((int(i)-1)//4)*.5
                assert all(math.isfinite(v) for v in (ux,uy,uz,rz))
                max_displacement = max(max_displacement,abs(ux),abs(uy),abs(uz+.1*block['time']*z))
        report.update(max_stress_error_Pa=max_stress, max_displacement_error_m=max_displacement,
                      final_stress_Pa=stresses[-1]['rows'], final_nodes=nodes[-1]['rows'])
        assert max_stress < .002 and max_displacement < 2e-6, 'Analytical comparison failed'
        report['status'] = 'passed'
    finally:
        report['files_sha256'] = {p.name:digest(p) for p in sorted(run.iterdir()) if p.is_file() and p.name!='benchmark.json'}
        (run/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
