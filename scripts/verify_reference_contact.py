#!/usr/bin/env python3
"""Native FEBio rigid-wall benchmark against confined neo-Hookean compression.

SI units. Full-face, laterally confined analytical limit; not regional indentation.
No third-party Python dependencies. Run with --self-test to exercise rejection checks.
"""
import argparse
import copy
import json
import math
import os
from pathlib import Path
import subprocess
from verify_reference_mechanics import digest, records

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/derived/mechanics-reference'
OUT = BASE / 'contact'
MU, LAM, H, A = 2800., 25200., .003, .0001
GAP, DEPTH = .00005, .0001


def mesh(n):
    nodes = [(i*.01/n, j*.01/n, k*H/n) for k in range(n+1) for j in range(n+1) for i in range(n+1)]
    def idx(i,j,k): return 1+i+(n+1)*j+(n+1)**2*k
    elements = [[idx(i,j,k),idx(i+1,j,k),idx(i+1,j+1,k),idx(i,j+1,k),
                 idx(i,j,k+1),idx(i+1,j,k+1),idx(i+1,j+1,k+1),idx(i,j+1,k+1)]
                for k in range(n) for j in range(n) for i in range(n)]
    faces = [e[4:] for e in elements if nodes[e[4]-1][2] == H]
    return nodes, elements, faces


def model(n, penalty):
    nodes, elems, faces = mesh(n)
    def tags(name, rows):
        return ''.join(f'<{name} id="{i}">' + ','.join(map(str,row)) + f'</{name}>' for i,row in enumerate(rows,1))
    def ns(name, ids): return f'<NodeSet name="{name}">' + ''.join(f'<node id="{i}"/>' for i in ids) + '</NodeSet>'
    return f'''<?xml version="1.0"?>
<febio_spec version="2.5"><Module type="solid"/>
<Control><analysis type="static"/><time_steps>20</time_steps><step_size>0.1</step_size><dtol>1e-8</dtol><etol>1e-10</etol><rtol>1e-8</rtol></Control>
<Material><material id="1" name="block" type="neo-Hookean"><E>8120</E><v>0.45</v></material></Material>
<Geometry><Nodes>{tags('node',nodes)}</Nodes><Elements type="hex8" mat="1" name="block">{tags('elem',elems)}</Elements>
{ns('all',range(1,len(nodes)+1))}{ns('bottom',[i for i,p in enumerate(nodes,1) if p[2]==0])}
<Surface name="top">{tags('quad4',faces)}</Surface></Geometry>
<Boundary><fix bc="x,y" node_set="all"/><fix bc="z" node_set="bottom"/></Boundary>
<Contact><contact type="rigid_wall" surface="top"><laugon>0</laugon><penalty>{penalty}</penalty><plane>0,0,-1,{-H-GAP}</plane><offset lc="1">{DEPTH+GAP}</offset></contact></Contact>
<LoadData><loadcurve id="1" type="linear"><point>0,0</point><point>1,1</point><point>2,0</point></loadcurve></LoadData>
<Output><logfile><node_data data="ux;uy;uz;Rz" file="nodes.txt" delim=",">{','.join(map(str,range(1,len(nodes)+1)))}</node_data><element_data data="J;sz;sed" file="elements.txt" delim=",">{','.join(map(str,range(1,len(elems)+1)))}</element_data></logfile><plotfile type="febio"><var type="displacement"/><var type="stress"/><var type="strain energy density"/></plotfile></Output></febio_spec>
'''


def force(s):
    return -A*(MU*(s-1/s)+LAM*math.log(s)/s)


def energy(s):
    return A*H*(MU/2*(s*s-1)-MU*math.log(s)+LAM/2*math.log(s)**2)


def check(n, penalty, nb, eb):
    nodes, elems, _ = mesh(n)
    assert len(nb)==len(eb)==21, 'Missing time states'
    top = [i for i,p in enumerate(nodes) if p[2]==H]
    steps=[]
    for j,(nblock,eblock) in enumerate(zip(nb,eb)):
        t = j*.1
        assert abs(nblock['time']-t)<1e-10 and abs(eblock['time']-t)<1e-10, 'Wrong time grid'
        nr, er = nblock['rows'], eblock['rows']
        for rows,count,width in ((nr,len(nodes),5),(er,len(elems),4)):
            assert len(rows)==count and [r[0] for r in rows]==list(range(1,count+1)), 'Wrong IDs/counts'
            assert all(len(r)==width and all(math.isfinite(v) for v in r) for r in rows), 'Invalid numeric data'
        u = sum(nr[i][3] for i in top)/len(top)
        s = 1+u/H
        assert 0<s<=1.000001, 'Invalid stretch'
        # Affine displacement verifies every element, not only averaged Jacobians.
        assert max(abs(r[1])+abs(r[2])+abs(r[3]-u*nodes[i][2]/H) for i,r in enumerate(nr))<2e-9, 'Non-affine deformation'
        assert min(r[1] for r in er)>0, 'Nonpositive Jacobian'
        assert max(abs(r[1]-s) for r in er)<2e-6, 'Jacobian inconsistent with displacement'
        # FEBio Rz is node.get_load(z), so sum loaded contact nodes, not fixed nodes.
        f = sum(nr[i][4] for i in top)
        exact = force(s)
        assert abs(f-exact)<2e-6, 'Reaction fails constitutive solution'
        assert max(abs(r[2]+exact/A) for r in er)<.03, 'Stress fails analytical limit'
        native_energy = sum(r[3] for r in er)*A*H/len(elems)
        assert abs(native_energy-energy(s))<1e-10, 'Native strain energy fails analytical limit'
        plane = H+GAP-(DEPTH+GAP)*(t if t<=1 else 2-t)
        penetration=max(0.,H+u-plane)
        assert abs(f-penalty*A*penetration)<3e-6, 'Contact force/penetration imbalance'
        assert f>=-1e-8 and penetration<2e-6, 'Invalid unilateral contact'
        if plane>=H:
            assert abs(f)<1e-8 and abs(u)<1e-10, 'Open gap has residual force/displacement'
        steps.append({'time':t,'reaction_N':f,'top_displacement_m':u,'penetration_m':penetration,
                      'minimum_J':min(r[1] for r in er),'native_strain_energy_J':native_energy,'analytical_energy_from_measured_stretch_J':energy(s)})
    assert steps[10]['reaction_N']>0 and steps[10]['top_displacement_m']<0, 'No loaded contact'
    assert abs(steps[-1]['reaction_N'])<1e-8 and abs(steps[-1]['top_displacement_m'])<1e-10, 'Release residual'
    ideal=force(1-DEPTH/H)
    return {'mesh_divisions':n,'elements':len(elems),'penalty_Pa_per_m':penalty,'steps':steps,
            'peak_force_N':steps[10]['reaction_N'],'ideal_hard_contact_force_N':ideal,
            'relative_force_error':abs(steps[10]['reaction_N']-ideal)/ideal,
            'peak_penetration_m':steps[10]['penetration_m']}


def check_refinement(results):
    assert len(results)==5, 'Expected three meshes and three penalties (shared baseline)'
    assert [(r['mesh_divisions'],r['penalty_Pa_per_m']) for r in results]==[(1,1e10),(2,1e10),(4,1e10),(2,1e9),(2,1e11)], 'Unexpected refinement cases'
    forces=[r['peak_force_N'] for r in results[:3]]
    assert max(forces)-min(forces)<2e-6, 'Mesh force invariance failed'
    errors=[]
    for i in (3,1,4):
        r=results[i]
        actual=abs(r['peak_force_N']-force(1-DEPTH/H))/force(1-DEPTH/H)
        assert math.isfinite(actual) and abs(actual-r['relative_force_error'])<1e-12, 'Reported force error inconsistent'
        errors.append(actual)
    assert 0<errors[2]<errors[1]/5<errors[0]/25, 'Penalty convergence failed'
    assert errors[2]<.0002, 'Hard-contact force limit failed'
    return {'mesh_force_spread_N':max(forces)-min(forces), 'penalty_relative_force_errors':errors,
            'mesh_interpretation':'Affine solution is represented exactly by all meshes; invariance, not an observed spatial order.'}


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--self-test', action='store_true'); args=parser.parse_args()
    if args.self_test:
        # Reuse native raw outputs, mutate independently, and ensure fail-closed validation.
        path=OUT/'n2-p1e10'; nb=records(path/'nodes.txt'); eb=records(path/'elements.txt')
        check(2,1e10,nb,eb)
        mutations=[]
        bad=copy.deepcopy(eb); bad[10]['rows'][0][1]=-1; mutations.append((nb,bad))
        bad=copy.deepcopy(nb); bad[10]['rows'][0][4]=float('nan'); mutations.append((bad,eb))
        mutations.append((nb[:-1],eb))
        bad=copy.deepcopy(eb); bad[10]['rows'][0][3]+=1.; mutations.append((nb,bad))
        bad=copy.deepcopy(nb); bad[-1]['rows'][-1][4]=.01; mutations.append((bad,eb))
        for a,b in mutations:
            try: check(2,1e10,a,b)
            except AssertionError: pass
            else: raise AssertionError('Invalid native record accepted')
        rr=json.loads((OUT/'benchmark.json').read_text())['cases']; check_refinement(rr)
        bad=copy.deepcopy(rr); bad[4]['relative_force_error']=bad[1]['relative_force_error']
        try: check_refinement(bad)
        except AssertionError: pass
        else: raise AssertionError('False convergence accepted')
        print('PASS: negative Jacobian, NaN, missing state, release residual, bad native energy, false convergence rejected'); return
    manifest=json.loads((BASE/'build-manifest.json').read_text()); exe=ROOT/manifest['executable']
    assert digest(exe)==manifest['executable_sha256'], 'Executable hash mismatch'
    for p,h in manifest['shared_library_sha256'].items(): assert digest(ROOT/p)==h, f'Library hash mismatch: {p}'
    OUT.mkdir(exist_ok=True)
    report={'status':'failed','backend':'FEBio','version':manifest['version'],'commit':manifest['commit'],
            'executable_sha256':digest(exe),'build_manifest_sha256':digest(BASE/'build-manifest.json'),
            'scope':'Full-face frictionless rigid plane, confined neo-Hookean block; not regional indentation or anatomical collision.', 'cases':[]}
    try:
        for n,p in [(1,1e10),(2,1e10),(4,1e10),(2,1e9),(2,1e11)]:
            run=OUT/f'n{n}-p{p:.0e}'.replace('+',''); run.mkdir(exist_ok=True)
            for name in ('nodes.txt','elements.txt','contact.log','contact.xplt','stdout.log'):
                (run/name).unlink(missing_ok=True)
            (run/'contact.feb').write_text(model(n,p))
            command=[str(exe),'-i','contact.feb','-nosplash']
            with (run/'stdout.log').open('w') as log:
                result=subprocess.run(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=120,env={**os.environ,'OMP_NUM_THREADS':'1'})
            assert result.returncode==0, f'Native execution failed: {run}'
            assert 'N O R M A L   T E R M I N A T I O N' in (run/'contact.log').read_text(), 'No normal termination'
            case=check(n,p,records(run/'nodes.txt'),records(run/'elements.txt'))
            case.update(command=command,returncode=result.returncode,files_sha256={f.name:digest(f) for f in sorted(run.iterdir()) if f.is_file()})
            report['cases'].append(case)
        report['refinement']=check_refinement(report['cases']); report['status']='passed'
    finally:
        (OUT/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'refinement':report['refinement']},indent=2))


if __name__=='__main__': main()
