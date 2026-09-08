"""A single generic body with explicit source ownership and numerical projections.

Native physiology owns fluid/thermal balances. Neural and mechanical reductions
consume its state; these transfers are assumptions, not fitted whole-body laws.
"""
from pathlib import Path
import csv, hashlib, json, math, io
import numpy as np
from .certainty import entity_certainty, summarize


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def write_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.partial.json');temp.write_text(json.dumps(data,allow_nan=False,separators=(',',':'))+'\n');temp.replace(path)


def build_bindings(anatomy):
    bindings=[]
    for e in anatomy['entities']:
        column=None;assumption=None
        if e['id'] in ('body-bp3d-FJ2422','body-bp3d-FJ2423'):
            side='Left' if e['id'].endswith('2422') else 'Right'
            column=side+'HeartVolume(mL)'
            assumption='Lumped native heart volume ratio projected onto ventricular cavity surface; atrial/ventricular partition and local shape are unresolved. This is not myocardium volume change.'
        elif 'lobe' in e['name'].lower() and 'lung' in e['name'].lower():
            side='Left' if 'left' in e['name'].lower() else 'Right' if 'right' in e['name'].lower() else None
            if side:
                column=side+'LungPulmonaryGasVolume(mL)'
                assumption='Uniform isotropic lobe expansion following ipsilateral gas volume ratio. Gas-to-total-tissue volume transfer, pleural sliding and lobar ventilation partition are unresolved.'
        if column:bindings.append(dict(entity_id=e['id'],column=column,evidence_kind='synthesized_transfer',assumption=assumption,calibrated_probability=None,reference='native stabilized t=0 volume equals canonical reference shape',owner='BioGears',direction='physiology_to_mechanics'))
    return bindings


def project_volumes(bindings,row,reference):
    result={}
    for b in bindings:
        key=b['column'];base=reference[key];value=row[key]
        if not math.isfinite(base) or base<=0 or not math.isfinite(value) or value<=0:raise ValueError('Positive finite native volume required')
        result[b['entity_id']]=value/base
    return result


def read_native(directory):
    directory=Path(directory)
    snapshots={name:(directory/name).read_bytes() for name in ['native_multisystem.csv','body_compartments.csv','summary.json']}
    input_hashes={name:hashlib.sha256(data).hexdigest() for name,data in snapshots.items()}
    summary=json.loads(snapshots['summary.json'])
    receipt=summary.get('compartment_telemetry',{})
    if input_hashes['native_multisystem.csv']!=summary.get('csv_sha256') or input_hashes['body_compartments.csv']!=receipt.get('csv_sha256'):
        raise ValueError('Native telemetry does not match its execution receipt; rerun or explicitly audit missing receipts')
    def read(name):
        with io.StringIO(snapshots[name].decode('utf-8')) as f:
            reader=csv.reader(f);header=next(reader,[])
            values=[[float(v) for v in row] for row in reader]
        if not header or len(set(header))!=len(header) or not values or any(len(row)!=len(header) for row in values):raise ValueError('Native CSV has missing, duplicate or ragged columns')
        rows=[dict(zip(header,row)) for row in values]
        if not all(math.isfinite(v) for row in rows for v in row.values()):raise ValueError('Native telemetry must be finite')
        return rows
    physiology=read('native_multisystem.csv');compartments=read('body_compartments.csv')
    times=np.array([r['#Time(s)'] for r in physiology]);ct=np.array([r['Time(s)'] for r in compartments])
    if len(ct)!=len(times)+1 or abs(ct[0])>1e-10 or not np.allclose(ct[1:],times,atol=1e-8,rtol=0):raise ValueError('Native physiology and compartments have different clocks')
    if times[0]<=0 or not (np.diff(times)>0).all() or np.max(np.diff(ct))>1:raise ValueError('Invalid native integration clock')
    return dict(input_hashes=input_hashes,time_s=times.tolist(),physiology=physiology,compartments=compartments[1:],initial_compartments=compartments[0],summary=summary)


SURFACE_OWNERSHIP_LANES=(
 ('cross_structure_repair','data/derived/cross-structure-repair-v1',
  ('manifest.json','summary.json','ownership-ledger.jsonl','coincidence.json','resolution.json',
   'conflict-graph.json','conflict-components.json','entity-operations.jsonl','clusters.json'),
  'every pairwise surface conflict in the acquired atlas, and which structure owns the disputed volume'),
 ('conflict_free_atlas','data/derived/conflict-free-atlas-v1',
  ('manifest.json','summary.json','conformity.json','residual-conflicts.json','annihilation.json',
   'duplicates.json','hole-fill.json','entities.jsonl','conflict-census.json'),
  'the repaired, mutually non-overlapping surfaces those ownership decisions imply'),
)


def surface_ownership(root,anatomy):
    """Declare the promoted surface-ownership lanes as canonical sources.

    Neither lane rewrites a reference surface, so neither can be declared through
    ``sources``; they are a second, conforming representation of the same anatomy plus the
    ledger that says who owns each disputed millilitre. Declaring them here with their
    digests is what makes them part of the model rather than an unpromoted candidate, and
    what lets any reader follow a structure back to the decision that shaped it.

    The coverage numbers are computed, not asserted. The ledger resolves the surfaces its
    two tet-ready parents produced; every canonical entity outside that set has no
    ownership decision, and the count of those is reported rather than left implied.
    """
    root=Path(root)
    lanes={}
    for key,relative,artifacts,role in SURFACE_OWNERSHIP_LANES:
        base=root/relative
        lanes[key]={'path':relative,'role':role,
                    'artifacts':{name:{'path':relative+'/'+name,'sha256':digest(base/name)}
                                 for name in artifacts if (base/name).exists()}}
    ledger=root/'data/derived/cross-structure-repair-v1/ownership-ledger.jsonl'
    decided=set()
    rows=0
    for line in ledger.read_text().splitlines():
        if not line.strip():continue
        record=json.loads(line);rows+=1
        decided.add(record['a']);decided.add(record['b'])
    present={e['id'] for e in anatomy['entities']}
    covered=json.loads((root/'data/derived/cross-structure-repair-v1/manifest.json').read_text())['per_entity_input_sha256']
    return {'lanes':lanes,'decisions':rows,
            'entities_with_a_decision':len(decided&present),
            'entities_with_a_resolved_surface':len(set(covered)&present),
            'canonical_entities':len(present),
            'entities_without_a_resolved_surface':len(present-set(covered)),
            'decisions_naming_an_entity_the_model_no_longer_has':len(decided-present),
            'not_regenerated':'the lanes are promoted as shipped. Their inputs are the tet-ready '
                              'surfaces, which are byte-identical to what they recorded, so a rerun '
                              'is a chance to introduce drift rather than a source of new evidence. '
                              'scripts/verify_cross_structure_ownership_reproduces.py replays every '
                              'decision through the shipped rule at the roles canonical carries now.',
            'coverage_limit':'entities without a resolved surface have no ownership decision at all. '
                             'They are not certified conflict-free; they are unresolved. The '
                             'tet-ready lanes have not been re-run over the promoted structures.'}


def provenance_index(root,anatomy):
    """Declare the promoted structure-provenance index as a canonical source.

    Every structure in this model came from a dataset, a script or an explicit prior, and
    anyone should be able to ask which. That answer lives in one index with one record
    shape, and pointing at it from the body manifest is what makes traceability part of
    the model rather than a side artifact. The completeness numbers are recomputed here
    from the records themselves, so this cannot claim coverage the records do not have.
    """
    root=Path(root)
    base=root/'data/derived/structure-provenance-v1'
    if not (base/'index.json').is_file():
        return {'present':False,
                'reason':'run scripts/build_structure_provenance.py --output '
                         'data/derived/structure-provenance-v1 --promoted'}
    index=json.loads((base/'index.json').read_text())
    complete=partial=absent=0
    missing_fields={}
    for e in anatomy['entities']:
        path=base/'records'/(e['id']+'.json')
        if not path.is_file():
            absent+=1;continue
        record=json.loads(path.read_text())
        gaps=record.get('completeness',{}).get('missing',[])
        if gaps:
            partial+=1
            for field in gaps:missing_fields[field]=missing_fields.get(field,0)+1
        else:
            complete+=1
    return {'present':True,
            'path':'data/derived/structure-provenance-v1',
            'schema':index['schema'],
            'index':{'path':'data/derived/structure-provenance-v1/index.json','sha256':digest(base/'index.json')},
            'manifest':{'path':'data/derived/structure-provenance-v1/manifest.json','sha256':digest(base/'manifest.json')},
            'schema_record':{'path':'data/derived/structure-provenance-v1/schema.json','sha256':digest(base/'schema.json')},
            'audit':{'path':'data/derived/structure-provenance-v1/audit.json','sha256':digest(base/'audit.json')},
            'records':index['record_count'],
            'canonical_entity_coverage':{'entities':len(anatomy['entities']),'complete':complete,
                                         'partial':partial,'no_record':absent,
                                         'missing_field_counts':missing_fields,
                                         'required_fields':'listed in the index schema; a record is '
                                                           'complete when none of them is empty'},
            'candidate_of_record':'data/derived/structure-provenance-candidate-v1'}


def build(root):
    root=Path(root);directory=root/'data/derived/canonical'
    anatomy=json.loads((directory/'anatomy.json').read_text());profile=json.loads((directory/'profile.json').read_text())
    sources={}
    for name in ['anatomy','profile','mechanics','brain','respiration','peripheral']:
        path=directory/(name+'.json');sources[name]={'path':str(path.relative_to(root)),'sha256':digest(path)}
    runtime_sources={name:{'path':'ihm/assembly/'+name+'.py','sha256':digest(root/'ihm/assembly'/f'{name}.py')} for name in ['body','body_protocol','body_runtime','body_states','contracts','interfaces','cosimulation','evidence','mechanics','brain','respiration','peripheral','certainty','temporal']}
    payload=dict(runtime_sources=runtime_sources,schema_version=1,model_id='ihm-body',name='IHM · one generic human',status='executable generic research assembly',validated_digital_twin=False,
        entity_count=len(anatomy['entities']),profile=profile,frame=anatomy['frame'],sources=sources,surface_ownership=surface_ownership(root,anatomy),provenance_index=provenance_index(root,anatomy),certainty=summarize(anatomy),volume_bindings=build_bindings(anatomy),
        coupling_contract=[
            {'owner':'BioGears','state':'cardiorespiratory, blood, renal, endocrine, tissue fluid, lymph aggregate and thermal balances','direction':'native integrated physiology','spatial_resolution':'lumped compartments; vascular surfaces are not solved 3D flow lumens'},
            {'owner':'canonical mechanics','state':'rigid translation, affine tissue deformation and muscle attachment forces','direction':'native volumes → explicit synthesized shape transfer','feedback_applied':False},
            {'owner':'thoracic mechanics','state':'three chest/diaphragm modes and regional skin displacement','direction':'native gas volume → constrained thorax → prescribed mechanics boundaries','feedback_applied':False,'mode':'native_replay'},
            {'owner':'somatic peripheral','state':'receptor adaptation, conduction queues and motor activation','direction':'stimuli/mechanics → delayed sensory inputs; explicit motor commands → delayed muscle activation','autonomic_control':False,'motor_policy':'explicit commands only; cortical activity is not an inferred recruitment policy'},
            {'owner':'IBM-derived brain','state':'80 neural populations','direction':'native perfusion/oxygen/temperature → neural rates','feedback_applied':False,'reason':'Native nervous system already owns autonomic control; applying a second uncalibrated controller would double count it.'},
            {'owner':'BETSE / human wound evidence','state':'local membrane and ion field; wound parameters','direction':'retained regional evidence and separately executable materialization','whole_body_voltage_field_calibrated':False},
            {'owner':'published lymphatic graph','state':'registered topology and estimated source lengths','direction':'anatomical prior','flow_solution_available':False},
            {'owner':'JOS3 / CSF / reproductive / vascular CFD','state':'retained domain materializations and parameter evidence','direction':'explicit alternative or regional models; not silently added to native conserved storage'}],
        limitations=[
            'Generic source-constrained body, not a validated individual digital twin or complete account of every microscopic component.',
            'Registration residuals, source dependence, unmeasured coefficients and synthesized anatomy remain explicit. Source participation does not establish calibrated interactions.',
            'Reduced mechanics does not resolve volumetric contact, joint articulation, vessel-wall fluid interaction or moving-body physiology.',
            'One-way partitioned coupling: mechanics and IBM neural outputs do not yet feed back into native physiology.',
            'Supine is a reference modeling condition; native engine posture limitations remain. The body display retains anatomical coordinates.',
            'Temporal spectra describe selected finite trajectories; they are not independently validated predictor transfer functions.',
            'Surface ownership is resolved for the acquired atlas only. Promoted display structures carry no ownership decision yet, so overlap between a promoted surface and an acquired one is unresolved rather than adjudicated.'])
    write_json(directory/'body.json',payload)
    return payload


class CanonicalBody:
    def __init__(self,root,payload,assets):
        self.root=Path(root);self.payload=payload;self.assets=assets;self.entities={e['id']:e for e in assets['anatomy']['entities']}
        if len(self.entities)!=payload['entity_count']:raise ValueError('Canonical entity count mismatch')
        if set(e['id'] for e in assets['mechanics']['entities'])!=set(self.entities):raise ValueError('Anatomy and mechanics identity mismatch; rebuild mechanics')
        if any(b['entity_id'] not in self.entities for b in payload['volume_bindings']):raise ValueError('Unknown volume binding')
    @classmethod
    def from_workspace(cls,root):
        root=Path(root).resolve();payload=json.loads((root/'data/derived/canonical/body.json').read_text());assets={}
        for name,record in payload['sources'].items():
            path=(root/record['path']).resolve()
            if not path.is_relative_to(root) or digest(path)!=record['sha256']:raise ValueError('Canonical source hash mismatch: '+name)
            assets[name]=json.loads(path.read_text())
        for name,record in payload['runtime_sources'].items():
            path=(root/record['path']).resolve()
            if not path.is_relative_to(root) or digest(path)!=record['sha256']:raise ValueError('Canonical executor changed: '+name)
        profile=assets['profile'];patient=(root/profile['native_patient_path']).resolve()
        if not patient.is_relative_to(root) or digest(patient)!=profile['native_patient_sha256']:raise ValueError('Canonical patient changed')
        return cls(root,payload,assets)
    def describe(self):return json.loads(json.dumps(self.payload))
    def certainty(self,entity_id):
        if entity_id not in self.entities:raise ValueError('Unknown canonical entity')
        return entity_certainty(self.entities[entity_id])
    def simulate(self,native_directory,output,output_hz=10,body_interventions=()):
        from .body_runtime import BodyRuntime
        from .body_protocol import validate_interventions,inputs_at
        if isinstance(output_hz,bool) or not isinstance(output_hz,(int,float)) or not math.isfinite(output_hz) or not 1<=output_hz<=50:raise ValueError('Output frequency must be 1..50 Hz')
        native_directory=Path(native_directory);series=read_native(native_directory)
        summary=series['summary'];profile=self.assets['profile']
        if summary['configuration']['patient']!=profile['native_patient'] or summary.get('patient_sha256')!=profile['native_patient_sha256']:raise ValueError('Native run does not match canonical patient identity')
        if summary.get('input_state_sha256') is not None:raise ValueError('Canonical profile matching requires fresh patient initialization, not an unverified saved state')
        protocol=validate_interventions(body_interventions,self.assets['peripheral'],series['time_s'][-1])
        runtime=BodyRuntime(self.root,self.assets,self.payload['volume_bindings'],series['initial_compartments'])
        frames=[];previous=0.;next_output=0.;max_force_residual=0.;max_torque_residual=0.
        # Sparse export drops only transforms indistinguishable at this explicit
        # numerical threshold; source geometry and full final state are retained.
        tolerance=1e-9;identity=np.eye(3)
        for time,row,comp in zip(series['time_s'],series['physiology'],series['compartments']):
            # Split at intervention boundaries; interpolate retained native
            # samples only for these subintervals and record that approximation.
            boundaries=sorted({event[key] for event in protocol for key in ['start_s','end_s'] if previous<event[key]<time})+[time]
            native_start=series['initial_compartments'] if previous==0 else prior_comp
            interval_start=previous
            for boundary in boundaries:
                alpha=(boundary-interval_start)/(time-interval_start)
                boundary_comp={key:native_start[key]+alpha*(comp[key]-native_start[key]) for key in comp}
                state=runtime.step(boundary-previous,row,boundary_comp,inputs_at(protocol,previous))
                previous=boundary
            prior_comp=comp
            mechanical,neural=state['mechanics'],state['brain']
            drivers={'volume_ratios':state['volume_ratios']}
            if abs(mechanical['time_s']-time)>1e-7 or abs(neural['time_s']-time)>1e-7:raise RuntimeError('Canonical component clocks diverged')
            max_force_residual=max(max_force_residual,mechanical['audit']['internal_force_residual_n']);max_torque_residual=max(max_torque_residual,mechanical['audit']['internal_torque_residual_nm'])
            if time+1e-9>=next_output or time==series['time_s'][-1]:
                explicit={}
                for id,entity_state in mechanical['entities'].items():
                    if max(np.max(np.abs(entity_state['translation_m'])),np.max(np.abs(np.array(entity_state['rotation_matrix'])-identity)),np.max(np.abs(np.array(entity_state['deformation_gradient'])-identity)))>tolerance:
                        explicit[id]={k:entity_state[k] for k in ['translation_m','rotation_matrix','deformation_gradient']}
                frames.append(dict(time_s=time,entities=explicit,physiology=row,compartments=comp,brain=neural,respiration=state['respiration'],peripheral=state['peripheral'],muscle_forces_n=mechanical['muscle_forces_n'],mechanical_audit=mechanical['audit'],volume_ratios=drivers['volume_ratios']))
                next_output=(math.floor(time*output_hz+1e-9)+1)/output_hz
                if len(frames)%50==0:print(f'Canonical integration {time:.2f}/{series["time_s"][-1]:.2f}s',flush=True)
        result=dict(runtime_sources=self.payload['runtime_sources'],schema_version=2,model_id='ihm-body',body_interventions=protocol,frames=frames,centroids_m={id:e['centroid_m'] for id,e in self.entities.items()},
            clock={'start_s':series['time_s'][0],'end_s':series['time_s'][-1],'native_samples':len(series['time_s']),'requested_output_hz':output_hz,'integration':'transactional reduced body intervals; motor/sensory outputs applied on subsequent interval','coupling_mode':'native_replay','interpolation':'native compartments linear only at intervention boundaries; physiology held from native interval endpoint; output on first native sample at/after display interval','playback':'recorded finite trajectory; no extrapolation'},
            sources={**self.payload['sources'],'native_run':{'directory':str(native_directory.resolve()),'files':series['input_hashes']}},
            executed_neural_source=dict(runtime.brain.source_identity),
            executed_neural_source_basis='The brain asset declares a retained baseline path; this records the law that actually integrated these frames.',
            volume_bindings=self.payload['volume_bindings'],coupling_contract=self.payload['coupling_contract'],
            display_numerics={'sparse_transform_component_tolerance':tolerance,'missing_entity_transform':'identity at reference centroid, never hold last value','geometry_reduction_applied':False,'biological_accuracy':None},
            audit={'maximum_internal_force_residual_n':max_force_residual,'maximum_internal_torque_residual_nm':max_torque_residual,'validated_digital_twin':False},limitations=self.payload['limitations'])
        output=Path(output);write_json(output,result);write_json(output.with_name(output.stem+'-final-state.json'),{'time_s':previous,**state,'physiology':series['physiology'][-1],'compartments':series['compartments'][-1]})
        return result
