"""Persistent IBM implicit cortex and 1 ms segmental cord on native plant ports.

Encoder/decoder initialization is an engineering materialization, not learned
motor control. Physiology modulates excitatory availability and neural timescales through explicit IHM engineering priors.
"""
from copy import deepcopy
import hashlib
import ast
import io
import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from .sensorimotor import SensorimotorController, SensorimotorParameters
from .reflexes import finite


def _module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class IBMImplicitController(SensorimotorController):
    @classmethod
    def from_root(cls, root, *, muscle_catalog, ibm_root=None, checkpoint_path=None,
                  sites=1024, seed=0, sever=False, no_cord=False, kind='implicit', **unused):
        import torch
        if not isinstance(kind, str) or not kind:
            raise ValueError('Controller kind must be a nonempty name')
        root = Path(root)
        ibm_root = Path(ibm_root or Path.home() / 'Documents/IBM-1').resolve()
        path = Path(checkpoint_path or ibm_root / 'ckpt/ibm1_implicit.pt').resolve()
        checkpoint_bytes = path.read_bytes()
        checkpoint = torch.load(io.BytesIO(checkpoint_bytes), map_location='cpu', weights_only=True)
        embed = checkpoint['dyn.embed']
        if not torch.isfinite(embed).all() or embed.ndim != 2:
            raise ValueError('Invalid implicit embedding')
        if type(sites) is not int or sites < 128 or sites > embed.shape[0]:
            raise ValueError('sites must be between128 and checkpoint site count')
        source = ibm_root / 'scripts/pretrain_video_loop.py'
        source_bytes = source.read_bytes()
        namespace = {'__name__':'_ihm_ibm_sensorimotor'}
        exec(compile(source_bytes, str(source), 'exec'), namespace)
        module = SimpleNamespace(**namespace)
        # Resolve cord's anatomy dependency without importing or replacing the
        # pinned IBM package used by the regional controller elsewhere.
        anatomy_path = ibm_root / 'ibm/anatomy/muscles.py'
        anatomy_bytes = anatomy_path.read_bytes()
        tree = ast.parse(anatomy_bytes)
        innervation = next(ast.literal_eval(node.value) for node in tree.body
            if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
            and node.target.id == 'INNERVATION')
        cord_path = ibm_root / 'ibm/processes/cord.py'
        cord_bytes = cord_path.read_bytes()
        cord_source = cord_bytes.decode().replace('from ibm.anatomy.muscles import INNERVATION', '')
        cord_namespace = {'INNERVATION': innervation, '__name__': '_ihm_ibm_cord'}
        exec(compile(cord_source, str(cord_path), 'exec'), cord_namespace)
        catalog = deepcopy(muscle_catalog)
        muscles = tuple(row['id'] for row in catalog)
        brain_bytes = (root / 'data/derived/canonical/brain.json').read_bytes()
        data = json.loads(brain_bytes)
        channels = tuple(sorted(row['id'] for row in data['nodes']))
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(seed)
            dyn = module.CorticalDynamics(sites, embed.shape[1], min(48, sites-1), 'cpu')
            sampled = torch.nn.functional.interpolate(embed.T[None], size=sites,
                        mode='linear', align_corners=False)[0].T
            dyn.embed.data.copy_(sampled)
            loop = module.SensorimotorLoop(dyn, list(muscles), afferent_channels=len(channels)).eval()
        if not len(loop.sense_idx) or not len(loop.motor_idx):
            raise ValueError('Materialization lacks sensory or motor sites')
        identity = {'checkpoint_sha256': hashlib.sha256(checkpoint_bytes).hexdigest(),
            'implementation_sha256': hashlib.sha256(source_bytes).hexdigest(),
            'innervation_sha256': hashlib.sha256(anatomy_bytes).hexdigest(),
            'cord_sha256': hashlib.sha256(cord_bytes).hexdigest(),
            'adapter_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'source_sites': int(embed.shape[0]), 'sites': sites, 'seed': seed,
            'resampling': 'linear embedding interpolation, align_corners=False',
            'muscle_catalog':catalog, 'brain_sha256':hashlib.sha256(brain_bytes).hexdigest(),
            'physiology_basis':'IHM body-transfer MAP/O2 availability scales excitatory drive; Q10 scales neural derivatives; uncalibrated engineering extension',
            'muscles': list(muscles), 'channels': list(channels), 'sever': sever, 'no_cord': no_cord}
        self = cls.__new__(cls)
        self.loop = loop
        self.source_bytes = {path.name:checkpoint_bytes,source.name:source_bytes,cord_path.name:cord_bytes,
            anatomy_path.name:anatomy_bytes,'brain.json':brain_bytes,Path(__file__).name:Path(__file__).read_bytes()}
        self.channels = channels
        self.catalog = catalog
        self.muscles = muscles
        self.bindings = {row['id']: row for row in catalog}
        if len(set(muscles)) != len(muscles) or not muscles:
            raise ValueError('Unique native muscle identities required')
        self.parameters = SensorimotorParameters()
        self.identity = identity
        self.model_sha256 = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        self.brain = SimpleNamespace(data=data, ids=list(channels), source_identity=identity,
                                     time_s=0., state=np.zeros((len(channels), 3)))
        self.time_s = 0.
        self.sever, self.no_cord = bool(sever), bool(no_cord)
        self.controller_metadata = {'kind':kind,'sever':self.sever,'no_cord':self.no_cord,
            'sites':sites,'trained_motor_policy':False,'model_sha256':self.model_sha256,
            # Which association kernel this cortex was materialized from. The digest
            # is the identity; the name is only there to be readable in a receipt.
            'kernel_name':path.name}
        self.cord = cord_namespace['SegmentalCord'](list(muscles), dt=.001, muscle_bindings=[dict(row, muscle_id=row['id']) for row in catalog])
        self.controller_metadata.update({
            'checkpoint_sha256':identity['checkpoint_sha256'],
            'muscle_count':len(muscles),
            'mapping_count':len(muscles)-len(self.cord.unmapped),
            'mapped_muscles':[m for m in muscles if m not in self.cord.unmapped],
            'unmapped_muscles':list(self.cord.unmapped)})
        self.antagonist = np.full(len(muscles), -1, dtype=np.int64)
        self.reciprocal_pairs = {}
        for side in ('r','l'):
            for target, origin in (('tibant','soleus'),('soleus','tibant'),('gasmed','tibant'),('gaslat','tibant')):
                target, origin = target+'_'+side, origin+'_'+side
                if target in muscles and origin in muscles:
                    self.antagonist[muscles.index(target)] = muscles.index(origin)
                    self.reciprocal_pairs[target] = origin
        self.controller_metadata['arc_availability'] = {
            name: not no_cord for name in ('stretch','autogenic','renshaw')}
        self.controller_metadata['arc_availability']['reciprocal'] = bool(self.reciprocal_pairs) and not no_cord
        self.controller_metadata['reciprocal_pairs'] = dict(self.reciprocal_pairs)
        self.controller_metadata['reciprocal_basis'] = 'Explicit ankle antagonist prior; tibant samples soleus, plantarflexors sample tibant; other channels unpaired'
        self.cortical_state = tuple(s.detach().clone() for s in dyn.init_state(1, 'cpu'))
        with torch.no_grad():
            self.weights = dyn.edge_weights().detach()
            if sever: self.weights = torch.zeros_like(self.weights)
        self.excitations = {m: 0. for m in muscles}
        return self

    def retain_sources(self, output):
        target = Path(output) / 'implicit-controller'
        target.mkdir(parents=True, exist_ok=True)
        for name, raw in self.source_bytes.items(): (target / name).write_bytes(raw)
        (target / 'manifest.json').write_text(json.dumps(self.identity, indent=2))
        return str(target)

    def step(self, dt_s, mechanical_observation, *, descending=None, sensory_blocks=(),
             motor_blocks=(), physiology=None, additional_sensory_inputs_hz=None):
        import torch
        dt = finite(dt_s, .001, .1, 'dt_s')
        ticks = round(dt / .001)
        if abs(dt - ticks * .001) > 1e-10:
            raise ValueError('IBM controller exchanges require integral1ms ticks')
        sensors = self._validate_observation(mechanical_observation)
        sb = self._blocks(sensory_blocks, 'sensory block')
        mb = self._blocks(motor_blocks, 'motor block')
        descending = {} if descending is None else descending
        if not isinstance(descending, dict) or set(descending)-set(self.muscles):
            raise ValueError('Unknown descending effector')
        descending = {k: finite(v, 0, 1, 'descending') for k,v in descending.items()}
        additional = {} if additional_sensory_inputs_hz is None else additional_sensory_inputs_hz
        if not isinstance(additional, dict) or set(additional)-set(self.channels):
            raise ValueError('Unknown additional sensory channel')
        rates = {k: 0. for k in self.channels}
        for k,v in additional.items(): rates[k] = finite(v, 0, 1000, 'sensory Hz')
        stretch = np.array([0. if m in sb else min(1., max(0., sensors[m]['length']-1)*2)
                            for m in self.muscles], np.float32)
        force = np.array([0. if m in sb else min(1., sensors[m]['force']) for m in self.muscles], np.float32)
        for i,m in enumerate(self.muscles):
            region = self.bindings[m]['sensory_region']
            rates[region] = min(1000., rates[region] + float(100*stretch[i]+40*force[i]))
        for m,v in descending.items():
            region = self.bindings[m]['motor_region']
            rates[region] = min(1000., rates[region] + 100*v)
        from .brain import INPUT_BASELINES
        supplied = {} if physiology is None else physiology
        if not isinstance(supplied,dict) or set(supplied)-set(INPUT_BASELINES):
            raise ValueError('Unknown physiology input')
        physiological = dict(INPUT_BASELINES, **supplied)
        pressure = finite(physiological['mean_arterial_pressure_mmHg'],0,300,'MAP')
        oxygen = finite(physiological['oxygen_saturation'],0,1,'oxygen saturation')
        temperature = finite(physiological['core_temperature_C'],20,45,'temperature')
        priors = self.brain.data['parameters']['body_transfer_priors']
        availability = min(1., pressure/priors['map_reference_mmHg'])*min(1.,oxygen/.98)
        temperature_factor = priors['temperature_Q10']**((temperature-37.)/10.)
        saved = self.checkpoint()
        try:
            with torch.no_grad():
                aff = torch.tensor([[rates[k] for k in self.channels]], dtype=torch.float32)
                if aff.shape[1] != self.loop.enc[0].in_features:
                    raise ValueError('Afferent width differs from encoder')
                drive = torch.zeros(1, self.loop.dyn.n)
                drive[:, self.loop.sense_idx] = self.loop.to_cortex(self.loop.enc(aff))
                # A block invalidates in-flight afference as well as the sample.
                for name in ('stretch', 'autogenic', 'reciprocal'):
                    for sample in self.cord._delay_buf.get(name, []):
                        for i,m in enumerate(self.muscles):
                            if (name != 'reciprocal' and m in sb or
                                name == 'reciprocal' and self.reciprocal_pairs.get(m) in sb): sample[i] = 0.
                arcs = {name: 0. for name in ('stretch','reciprocal','autogenic','renshaw')}
                for _ in range(ticks):
                    old = self.cortical_state
                    h = .001*temperature_factor
                    updated = self.loop.dyn.step(old, drive*availability, h, self.weights*availability)
                    missing_local = h*20*self.loop.dyn.w_ee*old[1]/self.loop.dyn.r_max*(1-availability)/self.loop.dyn.tau_m
                    self.cortical_state = (updated[0]-missing_local, *updated[1:])
                    if not all(torch.isfinite(x).all() for x in self.cortical_state):
                        raise ValueError('Nonfinite cortical state')
                    read = self.cortical_state[1][:, self.loop.motor_idx][:, self.loop.motor_read]
                    cortical = torch.sigmoid(self.loop.dec(read))[0].numpy().copy()
                    if self.no_cord: alpha = cortical
                    else:
                        result = self.cord.step(cortical, stretch=stretch, force=force, antagonist=self.antagonist)
                        alpha = result['alpha'].copy()
                        for name in arcs: arcs[name] = max(arcs[name], float(np.abs(result[name]).max()))
                    for i,m in enumerate(self.muscles):
                        if m in mb: alpha[i] = 0.
                self.excitations = dict(zip(self.muscles, map(float, alpha)))
                self.time_s += dt
                self.brain.time_s = self.time_s
                # Render only the actually materialized pooled strips. This
                # model has no hemisphere-resolved or medullary population map.
                indices = (self.loop.sense_idx, self.loop.motor_idx)
                region_ids = ['implicit-postcentral', 'implicit-precentral']
                region = {'node_ids': region_ids}
                for name, state in (('potential_mV',0),('activity_hz',1),('adaptation_mV',2)):
                    region[name] = [float(self.cortical_state[state][0,ix].mean()) for ix in indices]
                return {'schema':'ihm.sensorimotor.v1', 'time_s':self.time_s,
                    'motor_excitations':dict(self.excitations), 'requested_excitations':dict(self.excitations),
                    'cortical_commands':dict(zip(self.muscles,map(float,cortical))),
                    'sensors':deepcopy(sensors), 'arc_max':arcs,
                    'brain':{'time_s':self.time_s,'model_id':'ibm1-implicit','source_identity':self.identity,
                             'regional_state':region,'physiology_inputs':physiological,
                             'physiology_coupling_applied':True,'oxygen_perfusion_availability':availability,
                             'temperature_factor':temperature_factor,'physiology_basis':self.identity['physiology_basis']},
                    'controller':deepcopy(self.controller_metadata),
                    'brain_sensory_inputs_hz':rates,'model_sha256':self.model_sha256,
                    'neural_delay_s':.030,'exchange_interval_s':dt,
                    'activation_owner':'mechanical_plant','biological_validation':False,
                    'scope':'Persistent implicit cortical kernel; untrained sensory/motor heads; no walking policy',
                    'mapping':{'afferent_channels':list(self.channels),'muscles':list(self.muscles),
                               'unmapped_spinal_muscles':list(self.cord.unmapped)}}
        except Exception:
            self.restore(saved)
            raise

    def checkpoint(self):
        return {'schema':'ihm.ibm-implicit-state.v1','model_sha256':self.model_sha256,
                'time_s':self.time_s,'cortical_state':[x.tolist() for x in self.cortical_state],
                'alpha':self.cord.alpha.tolist(),'gamma':self.cord.gamma.tolist(),
                'delays':{k:[v.tolist() for v in values] for k,values in self.cord._delay_buf.items()},
                'excitations':dict(self.excitations)}

    def restore(self, checkpoint):
        import torch
        c = deepcopy(checkpoint)
        if not isinstance(c,dict) or c.get('schema')!='ihm.ibm-implicit-state.v1' or c.get('model_sha256')!=self.model_sha256:
            raise ValueError('Implicit checkpoint identity mismatch')
        time = finite(c.get('time_s'),0,1e9,'checkpoint time')
        states = [torch.tensor(x,dtype=torch.float32) for x in c['cortical_state']]
        if len(states)!=4 or any(x.shape!=(1,self.loop.dyn.n) or not torch.isfinite(x).all() for x in states):
            raise ValueError('Invalid cortical checkpoint')
        arrays = [np.asarray(c[k],np.float32) for k in ('alpha','gamma')]
        if any(x.shape!=(len(self.muscles),) or not np.isfinite(x).all() or (x<0).any() or (x>1).any() for x in arrays):
            raise ValueError('Invalid cord checkpoint')
        delays = {}
        expected = {'stretch':30,'autogenic':34,'reciprocal':32,'renshaw':4}
        required = set() if time == 0 or self.no_cord else set(expected)
        if not isinstance(c.get('delays'),dict) or set(c['delays']) != required:
            raise ValueError('Missing or unexpected active delay histories')
        for k,values in c['delays'].items():
            if k not in expected or len(values)!=expected[k]: raise ValueError('Invalid delay checkpoint')
            delays[k] = [np.asarray(x,np.float32) for x in values]
            if any(x.shape!=(len(self.muscles),) or not np.isfinite(x).all() for x in delays[k]):
                raise ValueError('Invalid delayed sample')
        exc = c['excitations']
        if set(exc)!=set(self.muscles): raise ValueError('Invalid checkpoint muscle identities')
        for v in exc.values(): finite(v,0,1,'checkpoint excitation')
        self.time_s = self.brain.time_s = time
        self.cortical_state = tuple(states)
        self.cord.alpha,self.cord.gamma = arrays
        self.cord._delay_buf = delays
        self.excitations = exc
