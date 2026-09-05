"""Evidence-backed implicit body and explicit, source-qualified materializations.

Distinct specimens and parameter families remain distinct. Missing cross-family
covariance or anatomical registration is never inferred from a shared name.
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
 'native':'data/derived/native-circuits/graph.json',
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
}

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
    def describe(self):
        return dict(schema_version=1,kind='heterogeneous_implicit_human',coverage=self._read('coverage')['summary'] if 'coverage' in self.assets else {},
            materializations=[name for name,key in [('population','population'),('skin-field','skin_field'),('skin-lymph','skin_lymph'),('temporal','temporal'),('native','native'),('opensim','opensim'),('reproductive','reproductive'),('csf','csf'),('thermal','thermal')] if key in self.assets],
            temporal_runs=[r['id'] for r in self._read('temporal')['runs']] if 'temporal' in self.assets else [],
            assets=deepcopy(self.assets),independently_validated_whole_human=False,
            coupling='Native systems are coupled inside their source engine. Cross-source anatomy, human measurements and reduced predictors retain explicit identities and bindings.',
            limitations=['No universal coefficient calibration or complete patient digital twin.',
                'Population covariance predicts concurrent measured states; it does not identify causal dynamics.',
                'Geometric families are separate specimens; names alone do not register anatomy.',
                'Frozen circuit responses and fitted temporal spectra have different meanings and validity domains.'])
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
            self._read(key)
        return self
