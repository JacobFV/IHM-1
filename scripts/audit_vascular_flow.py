"""Audit original cerebral CFD boundary flow, clock and P1 divergence identity."""
from pathlib import Path
import argparse
import hashlib
import json
import mmap
import re
import time
import xml.etree.ElementTree as ET
import numpy as np
from ihm.spatial.vascular import TetrahedralAuditMesh, surface_integrals
from ihm.spatial.vtk import read_arrays, surface

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'data/raw/vascular/vmr/extracted/0050_H_CERE_H'
SIM=CASE/'Simulations/0078_0001'
VTU=ROOT/'data/raw/vascular/vmr/extracted/0050_H_CERE_H_3D_RIGID.vtu'


def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def metadata(path):
    with Path(path).open('rb') as f,mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as mm:
        end=mm.find(b'<AppendedData')
        return ET.fromstring(mm[:end]+b'</VTKFile>' if end>=0 else mm[:])


def error_metrics(actual,expected):
    actual=np.asarray(actual);expected=np.asarray(expected);difference=actual-expected
    rms=float(np.sqrt(np.mean(difference**2)));scale=float(np.sqrt(np.mean(expected**2)))
    return {'rmse':rms,'maximum_absolute_error':float(np.max(np.abs(difference))),
            'mean_error':float(np.mean(difference)),
            'relative_rmse_to_expected_rms':rms/scale if scale>0 else None}


def read_bct(path):
    """Read the actual solver boundary table, retaining its saved time precision."""
    with path.open() as f:
        node_count,time_count=map(int,f.readline().split())
        coordinates=[];ids=[];vectors=[];times=None
        for _ in range(node_count):
            header=f.readline().split()
            if len(header)!=5 or int(header[3])!=time_count:raise ValueError('Unexpected saved BCT node header')
            coordinates.append([float(x) for x in header[:3]]);ids.append(int(header[4]))
            values=np.fromstring(''.join(f.readline() for _ in range(time_count)),sep=' ').reshape(time_count,4)
            if times is None:times=values[:,3]
            elif not np.array_equal(times,values[:,3]):raise ValueError('BCT nodes do not share the same clock')
            vectors.append(values[:,:3])
        if f.read().strip():raise ValueError('Unexpected trailing BCT data')
    vectors=np.stack(vectors,axis=1)
    if not np.isfinite(vectors).all() or not np.all(np.diff(times)>0):raise ValueError('Invalid saved BCT data')
    return np.array(coordinates),np.array(ids),times,vectors


def plot_audit(out,flow):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=np.array(flow['time_s']);caps={c['id']:c for c in flow['caps']}
    fig,axes=plt.subplots(2,2,figsize=(12,7.8),layout='constrained')
    ax=axes[0,0]
    for name in ['inflow','inflow_2']:
        if caps.get(name,{}).get('matching_complete'):ax.plot(t,caps[name]['outward_flow'],label=name)
    other=[c['outward_flow'] for c in flow['caps'] if c['id'].startswith('outlet_') and c['matching_complete']]
    if other:ax.plot(t,np.sum(other,axis=0),label='12 named outlets combined')
    ax.axhline(0,color='.6',lw=.8);ax.set_title('Signed cap flow: positive = outward');ax.set_ylabel('Source length³ / s');ax.legend(fontsize=8)
    ax=axes[0,1];ax.semilogy(t,np.maximum(flow['global']['relative_global_mass_imbalance'],1e-16));ax.set_title('Global flux imbalance / throughput');ax.set_ylabel('Dimensionless (log scale)')
    ax=axes[1,0]
    if 'inflow' in caps and 'bct_expected_outward_flow' in caps['inflow']:
        c=caps['inflow'];ax.plot(t,c['outward_flow'],label='Archived 3D integral');ax.plot(t,c['bct_expected_outward_flow'],'--',label='Saved BCT, periodic source clock');ax.legend(fontsize=8)
    ax.set_title('Prescribed boundary consistency');ax.set_ylabel('Source length³ / s')
    ax=axes[1,1];ax.plot(t,flow['global']['divergence_volume_rms']);ax.set_title('Local P1 divergence is not identically zero');ax.set_ylabel('P1 divergence RMS / source time⁻¹')
    for ax in axes.flat:ax.set_xlabel('Archived solver time / s');ax.grid(alpha=.2)
    fig.suptitle('Cerebral CFD numerical audit · archived source simulation\nUnits not independently confirmed; job metadata says “Simulation failed”',fontsize=12)
    fig.savefig(out/'audit.svg');fig.savefig(out/'audit.png',dpi=150);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stride',type=int,default=1,help='Audit every nth archived frame (default all 201)');parser.add_argument('--output',type=Path,default=ROOT/'data/derived/vascular');args=parser.parse_args()
    if args.stride<1:parser.error('stride must be positive')
    out=args.output;out.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    fields=metadata(VTU)
    steps=sorted(int(n.get('Name').split('_')[-1]) for n in fields.findall('.//PointData/DataArray') if re.fullmatch(r'velocity_\d+',n.get('Name','')))
    selected=steps[::args.stride]
    text=(SIM/'solver.inp').read_text();settings=dict(line.split(':',1) for line in text.splitlines() if ':' in line)
    dt=float(settings['Time Step Size']);times=np.array(selected)*dt
    source=read_arrays(VTU,{'Points/Points','Cells/connectivity','Cells/offsets','Cells/types'})
    if not np.all(source['Cells/types']==10) or not np.array_equal(source['Cells/offsets'],np.arange(1,len(source['Cells/types'])+1)*4):raise ValueError('Original cells are not uniformly linear tetrahedra')
    mesh=TetrahedralAuditMesh(source['Points/Points'],source['Cells/connectivity'].reshape(-1,4))
    simulation_mesh=read_arrays(SIM/'mesh-complete/mesh-complete.mesh.vtu',{'Points/Points','Cells/connectivity','PointData/GlobalNodeID'})
    coherent=np.array_equal(simulation_mesh['Points/Points'],mesh.points) and np.array_equal(np.sort(simulation_mesh['Cells/connectivity'].reshape(-1,4),axis=1),np.sort(mesh.tets,axis=1)) and np.array_equal(simulation_mesh['PointData/GlobalNodeID'],np.arange(1,len(mesh.points)+1))
    if not coherent:raise ValueError('Saved simulation mesh does not match archived CFD point ordering/connectivity')
    print(f"Verified original mesh: {len(mesh.points)} nodes, {len(mesh.tets)} tets, {len(mesh.boundary_faces)} boundary faces",flush=True)
    svpre=(SIM/'0078_0001.svpre').read_text()
    assigned={Path(path).stem:int(id) for path,id in re.findall(r'^set_surface_id_vtp\s+(\S+)\s+(\d+)',svpre,re.M)}
    prescribed={Path(p).stem for p in re.findall(r'^prescribed_velocities_vtp\s+(\S+)',svpre,re.M)}
    resistance={}
    for line in (SIM/'resistance.dat').read_text().splitlines():
        name,r,p=line.split();resistance[name]={'resistance':float(r),'reference_pressure':float(p)}
    matching={};caps=[];owners=np.zeros(len(mesh.boundary_faces),dtype=int)
    for file in sorted((SIM/'mesh-complete/mesh-surfaces').glob('*.vtp'))+[SIM/'mesh-complete/walls_combined.vtp']:
        points,faces=surface(file);ids=read_arrays(file,{'PointData/GlobalNodeID'}).get('PointData/GlobalNodeID')
        matched=mesh.match_surface(points,faces,ids)
        name=file.stem;matching[name]={k:v for k,v in matched.items() if k!='boundary_face_indices'}
        matching[name].update(source=str(file.relative_to(ROOT)),sha256=sha(file))
        indices=np.array(matched['boundary_face_indices'],dtype=int)
        if matched['complete']:owners[indices]+=1
        cap={'id':name,'matching_complete':matched['complete'],'_indices':indices,
             'area_source_length2':float(mesh.face_areas[indices].sum()) if matched['complete'] else None,
             'boundary_condition':{'type':'no slip' if name=='walls_combined' else 'prescribed velocity' if name in prescribed else 'resistance' if name in resistance else 'not declared in saved svpre',
                                   'surface_id':assigned.get(name),**resistance.get(name,{})},
             'outward_flow':[],'mean_pressure':[],'normal_velocity_area_rms':[],'maximum_nodal_speed':[]}
        caps.append(cap)
    # The BCT fixes signed flow without assuming that a cap named 'inflow' points inward.
    bp,bids,bt,bv=read_bct(SIM/'bct.dat');bct_reference={}
    original_bct_nodes=bids-1
    if np.max(np.linalg.norm(mesh.points[original_bct_nodes]-bp,axis=1))>np.linalg.norm(np.ptp(mesh.points,axis=0))*5e-7:raise ValueError('BCT point IDs do not match the original coordinates')
    for cap in caps:
        if cap['id'] not in prescribed or not cap['matching_complete']:continue
        local=np.full(len(mesh.points),-1,dtype=int);local[original_bct_nodes]=np.arange(len(bp))
        faces=local[mesh.boundary_faces[cap['_indices']]]
        if (faces<0).any():continue
        vectors=mesh.area_vectors[cap['_indices']]
        bq=np.array([surface_integrals(faces,vectors,u)['outward_flux'] for u in bv])
        period=float(re.search(r'^bct_period\s+(\S+)',svpre,re.M)[1]);phase=(times-bt[0])%period+bt[0]
        cap['bct_expected_outward_flow']=np.interp(phase,bt,bq).tolist()
        cap['bct_clock_note']='Original bct.dat node times, periodically extended by saved bct_period; no fitted phase shift.'
        bct_reference[cap['id']]={'points':len(bp),'saved_samples':len(bt),'period_s':period,'flow_sign':'integrated with outward tetrahedral boundary normals','source':'bct.dat'}
        waveform=np.loadtxt(SIM/'inflow.flow')
        cap['saved_waveform_value']=np.interp(phase,waveform[:,0],waveform[:,1]).tolist()
    global_series={key:[] for key in ['outward_flux','volume_divergence_integral','absolute_volume_divergence_integral','divergence_volume_rms','divergence_theorem_residual','divergence_theorem_relative_residual','relative_global_mass_imbalance','total_inward_flow','total_outward_flow','wall_outward_flow']}
    for index,(step,t) in enumerate(zip(selected,times)):
        names={f'PointData/velocity_{step:05d}',f'PointData/pressure_{step:05d}'}
        arrays=read_arrays(VTU,names);velocity=arrays[f'PointData/velocity_{step:05d}'];pressure=arrays[f'PointData/pressure_{step:05d}']
        if index==0:first_velocity=velocity.copy();first_pressure=pressure.copy()
        audit=mesh.audit(velocity,pressure);cap_flows=[];wall_flow=0.
        for cap in caps:
            if not cap['matching_complete']:continue
            result=mesh.integrate(velocity,pressure,cap['_indices'])
            for key,source_key in [('outward_flow','outward_flux'),('mean_pressure','mean_pressure'),('normal_velocity_area_rms','normal_velocity_area_rms'),('maximum_nodal_speed','maximum_nodal_speed')]:cap[key].append(result[source_key])
            if cap['id']=='walls_combined':wall_flow=result['outward_flux']
            else:cap_flows.append(result['outward_flux'])
        outward=sum(q for q in cap_flows if q>0);inward=-sum(q for q in cap_flows if q<0)
        audit.update(relative_global_mass_imbalance=abs(audit['outward_flux'])/max(outward,inward,np.finfo(float).tiny),total_inward_flow=inward,total_outward_flow=outward,wall_outward_flow=wall_flow)
        for key in global_series:global_series[key].append(audit[key])
        if index%20==0 or index==len(selected)-1:print(f'Frame {index+1}/{len(selected)} at {t:.3f}s · relative net flow {audit["relative_global_mass_imbalance"]:.3g}',flush=True)
    period=float(re.search(r'^bct_period\s+(\S+)',svpre,re.M)[1])
    endpoint_periodicity={'interval_s':float(times[-1]-times[0]),'prescribed_period_s':period,'same_prescribed_phase':bool(np.isclose((times[-1]-times[0])/period,round((times[-1]-times[0])/period),rtol=0,atol=1e-12) and len(times)>1),'velocity_nodal_difference':error_metrics(velocity,first_velocity),'pressure_nodal_difference':error_metrics(pressure,first_pressure),'interpretation':'Unweighted nodal endpoint differences, without fitted phase or pressure-gauge adjustment; matching inlet phase does not imply whole-state periodic convergence.'}
    condition_summary={}
    for cap in caps:
        cap.pop('_indices')
        if not cap['matching_complete']:continue
        bc=cap['boundary_condition'];summary={'condition':bc,'flow_minimum':min(cap['outward_flow']),'flow_maximum':max(cap['outward_flow'])}
        if 'bct_expected_outward_flow' in cap:
            summary['bct_comparison']=error_metrics(cap['outward_flow'],cap['bct_expected_outward_flow'])
            summary['raw_waveform_comparison']=error_metrics(cap['outward_flow'],cap['saved_waveform_value'])
        if bc['type']=='resistance':
            expected=bc['reference_pressure']+bc['resistance']*np.array(cap['outward_flow'])
            cap['resistance_expected_pressure']=expected.tolist()
            summary['pressure_minus_resistance_flow']=error_metrics(cap['mean_pressure'],expected)
            summary['interpretation']='Area-mean nodal pressure versus saved p_ref + R*Q_out; traction/viscous discretization may prevent exact pressure equality.'
        condition_summary[cap['id']]=summary
    units={'absolute_si_confirmed':False,'length':'original source length unit (unconfirmed)','velocity':'original source length per source time','pressure':'original pressure gauge and unit (unconfirmed absolute SI)','flow':'source length^3 / source time','time':'s inferred directly from solver.inp dt and original array step names',
           'evidence':['solver.inp density 1.06 and viscosity 0.04 match documented CGS examples','Original image filename OSMSC0078-cm.vti suggests centimeters','SimVascular documentation states the solver does not enforce a unit system'],
           'status':'CGS is supported as an inference, not a verified case-wide dimensional declaration; no SI conversion applied',
           'reference':'https://github.com/SimVascular/simvascular.github.io-archive/blob/master/docsFlowSolver.html'}
    limitations=['The archived project says Simulation failed. Conservation identities do not override that provenance or certify solver convergence.',
                 'The divergence theorem checks geometry, winding and integration of the same P1 field; it is not an independent validation of incompressibility.',
                 'Global net flow can cancel while local piecewise-linear divergence remains nonzero. Both are reported.',
                 'Cap names are source labels, not inferred flow directions; positive is outward from the tetrahedral volume.',
                 'inflow_2 exists and matches the archive boundary but has no assigned BC in the saved svpre; saved setup files may not fully describe the archived result.',
                 'Pressure values and resistance comparisons retain the archived gauge. No absolute physiologic pressure is inferred.',
                 'There is no mesh-refinement, time-refinement, full momentum-residual or independent in-vivo validation in this audit.']
    files=[VTU,SIM/'solver.inp',SIM/'0078_0001.svpre',SIM/'resistance.dat',SIM/'inflow.flow',SIM/'bct.dat',SIM.parent/'0078_0001.sjb',SIM/'mesh-complete/mesh-complete.mesh.vtu']
    report={'schema_version':1,'case_id':'0050_H_CERE_H','source_kind':'archived_source_cfd','archived_job_status':re.search(r'<mitk_job[^>]+status="([^"]+)"',(SIM.parent/'0078_0001.sjb').read_text())[1],
            'mesh':{**mesh.metadata,'original_simulation_nodes_and_tetrahedron_sets_equal':bool(coherent),'original_cell_winding_identical':bool(np.array_equal(simulation_mesh['Cells/connectivity'],mesh.tets.ravel()))},'units':units,
            'clock':{'solver_dt_s':dt,'archive_frames':len(steps),'audited_frames':len(selected),'stride':args.stride,'start_s':float(times[0]),'end_s':float(times[-1]),'array_step_names':selected},
            'matching':matching,'boundary_partition':{'unassigned_triangles':int((owners==0).sum()),'multiply_assigned_triangles':int((owners>1).sum()),'complete_disjoint_partition':bool(np.all(owners==1))},
            'conservation':{'maximum_divergence_theorem_relative_residual':max(global_series['divergence_theorem_relative_residual']),
                            'maximum_absolute_divergence_theorem_residual':max(abs(x) for x in global_series['divergence_theorem_residual']),
                            'maximum_relative_global_mass_imbalance':max(global_series['relative_global_mass_imbalance']),
                            'mean_relative_global_mass_imbalance':float(np.mean(global_series['relative_global_mass_imbalance'])),
                            'maximum_volume_divergence_rms':max(global_series['divergence_volume_rms']),
                            'maximum_wall_nodal_speed':max(next(c for c in caps if c['id']=='walls_combined')['maximum_nodal_speed'],default=None)},
            'periodic_endpoint_check':endpoint_periodicity,'definitions':{'outward_flow':'sum over oriented boundary triangles of mean nodal velocity dotted with outward area vector','relative_global_mass_imbalance':'absolute net full-boundary flow divided by max(total positive matched-cap flow, absolute total negative matched-cap flow); uniform density cancels','volume_divergence':'independently computed P1 tetrahedral shape gradients and exact tetrahedron volumes','normal_velocity_area_rms':'exact triangular integral of squared P1 normal velocity, divided by area, then square root'},'boundary_conditions':condition_summary,'bct_reference':bct_reference,'limitations':limitations,
            'provenance':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size} for p in files],
            'wall_seconds':time.monotonic()-started}
    flow={'schema_version':1,'case_id':'0050_H_CERE_H','source_kind':'archived_source_cfd','time_s':times.tolist(),'step':selected,'units':units,'caps':caps,'global':global_series,'limitations':limitations}
    (out/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');(out/'cap-flow.json').write_text(json.dumps(flow,separators=(',',':'),allow_nan=False)+'\n')
    plot_audit(out,flow)
    print(json.dumps({'audited_frames':len(selected),'conservation':report['conservation'],'boundary_partition':report['boundary_partition'],'output':str(out)},indent=2),flush=True)


if __name__=='__main__':main()
