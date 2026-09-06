"""Evidence-backed implicit body and explicit, source-qualified materializations.

Source identities and parameter families remain traceable. The canonical body
registers and assembles them with explicit fits and assumptions; shared names
alone never imply measured covariance or anatomical registration.
"""
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import numpy as np
from ihm.runtime.gaussian import covariance

ASSETS={
 'coverage':'data/derived/system-coverage.json',
 'anatomy':'data/derived/app/manifest.json',
 'anatomy_coverage':'data/derived/anatomy/anatomy_coverage.json',
 'anatomy_fidelity':'data/derived/anatomy/fidelity.json',
 'lymph_network':'data/derived/lymphatic/graph.json',
 'native':'data/derived/native-circuits/graph.json',
 'systemic_backend':'data/runtime/physiology/native_biogears_stream.manifest.json',
 'population':'data/derived/population/nhanes-2017-2018/joint-population-prior.json',
 'skin_field':'data/derived/calibration/skin-fit.json',
 'skin_lymph':'data/derived/coupling/native-skin-circuit.json',
 'temporal':'data/derived/temporal/index.json',
 'opensim':'data/derived/opensim/native_corrected/baseline/summary.json',
 'reproductive':'data/derived/reproductive/index.json',
 'csf':'data/derived/csf/index.json',
 'bioelectric':'data/derived/bioelectric/tissue.json',
 'native_targets':'data/derived/calibration/native-target-audit.json',
 'thermal':'data/derived/thermal/index.json',
 'canonical_anatomy':'data/derived/canonical/anatomy.json',
 'canonical_profile':'data/derived/canonical/profile.json',
 'canonical_brain':'data/derived/canonical/brain.json',
 'canonical_mechanics':'data/derived/canonical/mechanics.json',
 'canonical_body':'data/derived/canonical/body.json',
 'canonical_respiration':'data/derived/canonical/respiration.json',
 'canonical_peripheral':'data/derived/canonical/peripheral.json',
 'hair_strands':'data/derived/hair/elastic_v3/manifest_fragment.json',
 'hair_dynamics':'app/src/hair_dynamics.js',
 'canonical_details':'data/derived/canonical/details.json',
 'canonical_microvascular':'data/derived/canonical/microvascular.json',
 'ibm_source':'data/derived/canonical/ibm-backend/manifest.json',
 'kidney_card':'data/sources/hipct-kidney-arterial.json',
 'kidney_graph':'data/derived/microstructure/kidney/example_graph.npz',
 'kidney_statistics':'data/derived/microstructure/kidney/statistics.json',
 'kidney_slab':'data/derived/microstructure/kidney/slab_validation.json',
 'penile_constitutive':'data/measurements/biomechanics/khorshidi_2024.json',
 'penile_constitutive_source':'data/raw/biomechanics/human-penile-mechanics-2024/paper.pdf',
 'penile_volume':'data/derived/material-domains/pelvis-0.004m/manifest.json',
 'penile_volume_data':'data/derived/material-domains/pelvis-0.004m/pelvic-domain.npz',
}
CANONICAL_ASSETS = ('canonical_anatomy', 'canonical_profile', 'canonical_brain',
                    'canonical_mechanics', 'canonical_body')

def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

@dataclass
class PopulationBelief:
    components: list
    units: list
    mean: np.ndarray
    cov: np.ndarray
    provenance: dict
    def __post_init__(self):
        self.components=list(self.components);self.units=list(self.units)
        self.mean=np.array(self.mean,float,copy=True)
        n=len(self.components)
        if len(set(self.components))!=n or len(self.units)!=n or self.mean.shape!=(n,) or not np.isfinite(self.mean).all():raise ValueError('invalid population components, units or mean')
        self.cov=covariance(self.cov,n).copy()
    def condition(self,observations):
        """Return a new belief from {component:(value, independent variance, unit)}.

        This predicts concurrent population state; it introduces no dynamics or
        calibrated intervention response. Correlated assay error requires a
        separately specified joint observation model.
        """
        if not isinstance(observations,dict) or not observations:raise ValueError('nonempty observation mapping required')
        ids=[];values=[];variances=[]
        for name,item in observations.items():
            if name not in self.components or len(item)!=3:raise ValueError('unknown component or invalid observation')
            value,var,unit=item;i=self.components.index(name)
            if not np.isfinite(value) or not np.isfinite(var) or var<=0 or unit!=self.units[i]:raise ValueError('finite observation, positive variance and exact declared unit required')
            ids.append(i);values.append(value);variances.append(var)
        # Factor on correlation scale so heterogeneous physical units do not
        # determine numerical rank. Whiten observations rather than adding tiny
        # noise to a singular C_obs, where floating point can erase it entirely.
        scale=np.sqrt(np.diag(self.cov));safe=np.where(scale>0,scale,1.)
        correlation=self.cov/safe[:,None]/safe[None,:]
        eigenvalues,eigenvectors=np.linalg.eigh(correlation)
        # Discard correlation-scale roundoff modes, consistent with covariance
        # validation; otherwise tiny assay noise falsely resolves exact aliases.
        rank_tolerance=64*np.finfo(float).eps*len(scale)*max(1.,float(eigenvalues.max()))
        active=eigenvalues>rank_tolerance
        factor=safe[:,None]*eigenvectors[:,active]*np.sqrt(eigenvalues[active])
        if factor.shape[1]:
            noise_sd=np.sqrt(variances)
            design=factor[ids]/noise_sd[:,None]
            residual=(np.array(values)-self.mean[ids])/noise_sd
            if not np.isfinite(design).all() or not np.isfinite(residual).all():raise ValueError('observation scaling overflow')
            u,s,vt=np.linalg.svd(design,full_matrices=True)
            root=np.hypot(1.,s);gain=(s/root)/root
            latent=vt[:len(s)].T@(gain*(u[:,:len(s)].T@residual))
            inverse_root=np.ones(factor.shape[1]);inverse_root[:len(s)]=1/root
            posterior_factor=(factor@vt.T)*inverse_root
            mean=self.mean+factor@latent
            cov=posterior_factor@posterior_factor.T
        else:
            mean=self.mean.copy();cov=self.cov.copy()
        provenance={**self.provenance,'conditioning':[*self.provenance.get('conditioning',[]),dict(observations)]}
        return PopulationBelief(self.components,self.units,mean,(cov+cov.T)/2,provenance)
    def to_dict(self):
        return dict(kind='concurrent_population_gaussian',components=self.components,units=self.units,mean=self.mean.tolist(),covariance=self.cov.tolist(),provenance=self.provenance,temporal_dynamics=None)

@dataclass(frozen=True)
class SkinFieldPredictor:
    evidence: dict
    def predict(self,age_group,sex,site):
        from ihm.calibration.skin import predict_skin
        return predict_skin(self.evidence,age_group,sex,site)

@dataclass(frozen=True)
class NativePredictor:
    config: object
    def run(self,output_dir):
        from ihm.native import run_native
        return run_native(self.config,output_dir)

@dataclass(frozen=True)
class SystemicPredictor:
    root: Path
    config: object
    sources: dict
    def run(self,output_dir):
        for path,sha in self.sources.items():
            if digest(path)!=sha:raise ValueError('Systemic input changed after materialization: '+str(path))
        from ihm.assembly.systemic import run_systemic
        return run_systemic(self.root,output_dir,self.config)

@dataclass(frozen=True)
class OpenSimPredictor:
    config: object
    def run(self,output_dir):
        from ihm.native.opensim_backend import run_opensim
        return run_opensim(self.config,output_dir)

@dataclass(frozen=True)
class ReproductivePredictor:
    root: Path
    options: dict
    def run(self):
        from ihm.native.reproductive import run_reproductive
        return run_reproductive(self.root,**self.options)

@dataclass(frozen=True)
class CSFPredictor:
    root: Path
    options: dict
    def run(self):
        from ihm.native.csf import run_csf
        return run_csf(self.root,**self.options)

@dataclass(frozen=True)
class ThermalPredictor:
    root: Path
    options: dict
    def run(self):
        from ihm.native.thermal import run_thermal
        return run_thermal(self.root,**self.options)

@dataclass(frozen=True)
class RegionalTouchPredictor:
    root: Path
    options: dict
    sources: dict
    def run(self):
        for path,sha in self.sources.items():
            if digest(self.root/path)!=sha:raise ValueError('Body evidence changed after materialization: '+path)
        from ihm.assembly.regional_touch import run_touch
        return run_touch(self.root,**self.options)

@dataclass(frozen=True)
class RegionalElectricPredictor:
    root: Path
    options: dict
    sources: dict
    def run(self):
        for path,sha in self.sources.items():
            if digest(self.root/path)!=sha:raise ValueError('Body evidence changed after materialization: '+path)
        from ihm.assembly.skin_bioelectric import build_skin_electric
        return build_skin_electric(self.root,**self.options)

class ImplicitHuman:
    @classmethod
    def open(cls,root=None):
        self=cls();self.root=Path(root or Path(__file__).resolve().parents[1]).resolve()
        self.assets={k:dict(path=v,sha256=digest(self.root/v)) for k,v in ASSETS.items() if (self.root/v).is_file()}
        self._cache={};return self
    def _read(self,key):
        if key not in self.assets:raise ValueError('Required evidence unavailable: '+key)
        if key not in self._cache:
            asset=self.assets[key];path=(self.root/asset['path']).resolve()
            if not path.is_relative_to(self.root) or digest(path)!=asset['sha256']:raise ValueError('Evidence changed: '+key+'; reopen or rebuild the substrate')
            self._cache[key]=json.loads(path.read_text())
        return deepcopy(self._cache[key])
    def _executed_neural_source(self):
        """Identity of the law the brain materializations execute, not the baseline asset."""
        from ihm.brain.active_source import resolve_source
        try:pin=resolve_source(self.root)
        except ValueError as error:return {'available':False,'reason':str(error)}
        return {'available':True,**pin.to_dict()}

    def describe(self):
        result=dict(schema_version=1,kind='heterogeneous_implicit_human',coverage=self._read('coverage')['summary'] if 'coverage' in self.assets else {},
            materializations=[name for name,key in [('population','population'),('skin-field','skin_field'),('skin-lymph','skin_lymph'),('temporal','temporal'),('native','native'),('opensim','opensim'),('reproductive','reproductive'),('csf','csf'),('thermal','thermal'),('body','canonical_body'),('ibm-causal','ibm_source'),('body-touch','canonical_microvascular'),('body-skin-transport','canonical_microvascular'),('body-skin-electric','canonical_microvascular'),('kidney-arterial-geometry','kidney_graph')] if key in self.assets and all(k in self.assets for k in {'body':CANONICAL_ASSETS,'body-touch':('ibm_source','canonical_anatomy'),'body-skin-transport':('skin_lymph','canonical_anatomy'),'body-skin-electric':('skin_field','canonical_anatomy'),'kidney-arterial-geometry':('kidney_card','kidney_statistics')}.get(name,()))],
            temporal_runs=[r['id'] for r in self._read('temporal')['runs']] if 'temporal' in self.assets else [],
            assets=deepcopy(self.assets),independently_validated_whole_human=False,
            executed_neural_source=self._executed_neural_source(),
            coupling='Native systems are coupled inside their source engine. The canonical body assembles registered anatomy, mechanics and brain state with explicit physiological drivers and source assumptions.',
            limitations=['No universal coefficient calibration or complete patient digital twin.',
                'Population covariance predicts concurrent measured states; it does not identify causal dynamics.',
                'Canonical anatomy uses recorded inter-template fits and synthesis priors; source families retain distinct specimen identities.',
                'Frozen circuit responses and fitted temporal spectra have different meanings and validity domains.'])
        if all(k in self.assets for k in ('hair_strands','hair_dynamics')):result['materializations'].append('hair-strands')
        if 'systemic_backend' in self.assets:result['materializations'].append('body-systemic')
        if all(key in self.assets for key in ('penile_constitutive','penile_constitutive_source')):
            result['materializations'].append('penile-constitutive')
            if all(key in self.assets for key in ('penile_volume','penile_volume_data')):
                result['materializations'].append('penile-volume')
        return result
    def microstructure_evidence(self):
        """Acquired organ evidence; availability does not confer population validity."""
        if 'kidney_card' not in self.assets:return {}
        card=self._read('kidney_card')
        return {'kidney':dict(tier=card['source_tier'],body_registered=False,
            usable_for_population_priors=card['acquired_graph']['usable_for_population_priors'],
            source_card=card,statistics=self._read('kidney_statistics') if 'kidney_statistics' in self.assets else None,
            paired_slab=self._read('kidney_slab') if 'kidney_slab' in self.assets else None,
            evidence_receipts={k:deepcopy(v) for k,v in self.assets.items() if k.startswith('kidney_')})}
    def fields(self):
        """Concrete native scalar state/parameter fields with units and support IDs."""
        graph=self._read('native');fields=[]
        for category in ('nodes','paths','systems'):
            for owner in graph[category]:
                for name,value in owner.get('properties',{}).items():
                    if not isinstance(value,dict):continue
                    fields.append(dict(id=f"native/{category}/{owner.get('id',owner.get('type'))}/{name}",support=owner.get('id',owner.get('type')),quantity=name,**value,
                        evidence_kind='native_initialized_state_or_parameter',source=graph['source'],independently_calibrated=False))
        return fields
    def materialize(self,kind,**options):
        if kind=='hair-strands':
            if options:raise ValueError('Hair materialization returns its explicit retained guide/render resolution')
            for key in ('hair_strands','hair_dynamics'):
                if key not in self.assets or digest(self.root/self.assets[key]['path'])!=self.assets[key]['sha256']:
                    raise ValueError('Hair evidence changed; reopen the implicit body')
            from ihm.app.experiments import read_experiment
            return read_experiment(self.root,'hair-strands')
        if kind=='penile-volume':
            if options:raise ValueError('The source-volume materialization uses its explicit retained CC/CS partition')
            for key in ('penile_volume','penile_volume_data','penile_constitutive','penile_constitutive_source'):
                if key not in self.assets or digest(self.root/self.assets[key]['path'])!=self.assets[key]['sha256']:
                    raise ValueError('Volume or constitutive evidence changed; reopen the implicit body')
            from ihm.assembly.penile_volume import materialize_penile_volume
            return materialize_penile_volume(self.root)
        if kind=='penile-constitutive':
            for key in ('penile_constitutive','penile_constitutive_source'):
                if key not in self.assets or digest(self.root/self.assets[key]['path'])!=self.assets[key]['sha256']:
                    raise ValueError('Constitutive evidence changed or is unavailable; reopen the implicit body')
            from ihm.calibration.penile import penile_material
            return penile_material(evidence_root=self.root,**options)
        if kind=='body-systemic':
            from ihm.assembly.systemic import SystemicConfig
            from ihm.native.session import RUNTIME
            config=SystemicConfig(**options)
            library=(RUNTIME/'biogears-build/outputs/Release/lib' if config.engine_variant=='upstream'
                     else RUNTIME/'variants'/config.engine_variant)/'libbiogears.so.8.0.0'
            paths=[Path(config.state_path).resolve(), RUNTIME/'native_biogears_stream', library,
                   self.root/'ihm/assembly/systemic.py',self.root/'ihm/assembly/systemic_evidence.py',
                   self.root/'ihm/assembly/native_environment_evidence.py',
                   self.root/'ihm/native/session.py',self.root/'ihm/native/__init__.py',
                   self.root/'scripts/native_body_ports.h',self.root/'scripts/native_biogears_stream.cpp']
            return SystemicPredictor(self.root,config,{str(path):digest(path) for path in paths})
        if kind=='kidney-arterial-geometry':
            if options:raise ValueError('Measured donor geometry accepts no synthesis or registration options')
            for key in ('kidney_card','kidney_graph','kidney_statistics'):
                if key not in self.assets:raise ValueError('Required evidence unavailable: '+key)
                if digest(self.root/self.assets[key]['path'])!=self.assets[key]['sha256']:
                    raise ValueError('Evidence changed: '+key+'; reopen the implicit body')
            from ihm.anatomy.kidney_graph import validate_graph
            with np.load(self.root/self.assets['kidney_graph']['path'],allow_pickle=False) as archive:
                graph={k:archive[k] for k in archive.files}
            validate_graph(graph)
            return dict(kind=kind,coordinate_frame='original donor graph local coordinates',
                nodes_m=graph['nodes_mm']*.001,points_m=graph['points_mm']*.001,
                edges=graph['edges'],point_offsets=graph['point_offsets'],
                thickness_native=graph['thickness_native'],radius_m=None,flow_solution=None,
                source_transform=graph.get('source_transform'),source_transform_applied=False,
                body_registered=False,evidence=self.microstructure_evidence()['kidney'])
        if kind in ('ibm-causal','body-touch','body-skin-transport','body-skin-electric'):
            required={'ibm-causal':('ibm_source',),'body-touch':('ibm_source','canonical_anatomy','canonical_microvascular'),
                'body-skin-transport':('skin_lymph','canonical_anatomy','canonical_microvascular'),
                'body-skin-electric':('skin_field','canonical_anatomy','canonical_microvascular')}[kind]
            for key in required:
                if key not in self.assets:raise ValueError('Required evidence unavailable: '+key)
                if digest(self.root/self.assets[key]['path'])!=self.assets[key]['sha256']:
                    raise ValueError('Evidence changed: '+key+'; reopen the implicit body')
            if kind in ('ibm-causal','body-touch'):
                # The ibm_source asset is the retained baseline; these two kinds execute
                # the active source, so gate on the artifact that actually runs as well.
                from ihm.brain.active_source import resolve_source
                resolve_source(self.root)
            if kind=='ibm-causal':
                from ihm.brain.ibm_backend import IBMBackend
                from ihm.brain.causal import CausalIBM
                response_kind=options.pop('response_kind','rapid')
                return CausalIBM(IBMBackend(root=self.root),kind=response_kind,**options)
            sources={self.assets[key]['path']:self.assets[key]['sha256'] for key in required}
            if kind=='body-touch':return RegionalTouchPredictor(self.root,deepcopy(options),sources)
            if kind=='body-skin-electric':return RegionalElectricPredictor(self.root,deepcopy(options),sources)
            from ihm.assembly.skin_transport import RegionalSkinTransport
            return RegionalSkinTransport.from_sources(self.root,**options)
        if kind=='body':
            if options:raise ValueError('body options belong to simulate()')
            # from_workspace reads current files; recheck even previously cached
            # assets so this materialization cannot silently switch evidence.
            for key in CANONICAL_ASSETS:
                if key not in self.assets:raise ValueError('Required evidence unavailable: '+key)
                asset=self.assets[key]
                path=(self.root/asset['path']).resolve()
                if not path.is_relative_to(self.root) or digest(path)!=asset['sha256']:
                    raise ValueError('Evidence changed: '+key+'; reopen or rebuild the substrate')
            from ihm.assembly.body import CanonicalBody
            return CanonicalBody.from_workspace(self.root)
        if kind=='population':
            if options:raise ValueError('population materialization has no implicit dynamics options')
            d=self._read('population')
            if set(d['fit_subject_ids'])&set(d['holdout_subject_ids']):raise ValueError('population holdout leakage')
            return PopulationBelief(d['components'],d['units'],d['mean'],d['covariance'],{**self.assets['population'],'source':d['source'],'basis':d['basis'],'limitations':d['limitations']})
        if kind=='skin-field':
            if options:raise ValueError('skin-field options belong to predict()')
            return SkinFieldPredictor(self._read('skin_field'))
        if kind=='skin-lymph':
            if options:raise ValueError('boundary perturbations belong to step() or response()')
            from ihm.coupling.circuit import FluidCircuit
            return FluidCircuit.from_dict(self._read('skin_lymph')['model'])
        if kind=='temporal':
            if set(options)!={'run_id'}:raise ValueError('temporal requires exactly run_id')
            from ihm.temporal.predictor import ReducedPredictor
            matches=[r for r in self._read('temporal')['runs'] if r['id']==options['run_id']]
            if len(matches)!=1:raise ValueError('Unknown temporal source run')
            return ReducedPredictor.from_dict(matches[0]['predictor']['model'])
        if kind=='opensim':
            from ihm.native.opensim_backend import OpenSimConfig
            options.setdefault('engine_variant','wrap_8_0.0005_cache')
            return OpenSimPredictor(OpenSimConfig(**options))
        if kind=='thermal':
            return ThermalPredictor(self.root,dict(options))
        if kind=='csf':
            return CSFPredictor(self.root,dict(options))
        if kind=='reproductive':
            return ReproductivePredictor(self.root,dict(options))
        if kind=='native':
            from ihm.native import NativeConfig
            return NativePredictor(NativeConfig.from_dict(options))
        raise ValueError('Unknown materialization: '+kind)
    def save(self,path):
        # Recheck even cached assets when freezing a reproducible manifest.
        for key,asset in self.assets.items():
            if digest(self.root/asset['path'])!=asset['sha256']:raise ValueError('Evidence changed: '+key)
        Path(path).parent.mkdir(parents=True,exist_ok=True)
        Path(path).write_text(json.dumps(self.describe(),indent=2,allow_nan=False)+'\n')
    @classmethod
    def load(cls,path,root=None):
        data=json.loads(Path(path).read_text())
        if data.get('schema_version')!=1 or data.get('kind')!='heterogeneous_implicit_human':raise ValueError('unsupported human manifest')
        self=cls.open(root);self.assets=data['assets'];self._cache={}
        for key in self.assets:
            if key not in ASSETS or self.assets[key]['path']!=ASSETS[key]:raise ValueError('unsupported evidence location')
            if key in ('kidney_graph','penile_constitutive_source','penile_volume_data'):
                asset=self.assets[key];path=(self.root/asset['path']).resolve()
                if not path.is_relative_to(self.root) or digest(path)!=asset['sha256']:
                    raise ValueError('Evidence changed: '+key+'; reopen or rebuild the substrate')
            else:self._read(key)
        return self
