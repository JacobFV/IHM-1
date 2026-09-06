"""Shadow partition of existing fluid owners; no native mutation or ATP source.

Every call is one alternative snapshot, not a second simultaneously live pool.
Native commit is deliberately unavailable until native diffusion exclusions,
atomic epoch checks, and electrical/thermal state persistence exist.
"""
import copy
import gzip
import hashlib
import json
import numpy as np
from .epithelial_electrodiffusion import EpithelialPatch, F, Z, simulate
from .skin_layers import physical_skin_support

IONS=('Sodium','Potassium','Chloride')


def exact_patch_support(reference, geometry_bytes, evidence_bytes, triangle_ids):
    """Area from explicit source triangles within reviewed exterior proxy."""
    support=physical_skin_support(reference,geometry_bytes,evidence_bytes)
    evidence=json.loads(evidence_bytes)
    ids=list(triangle_ids)
    if not ids or any(type(i) is not int for i in ids) or ids!=sorted(set(ids)) or not set(ids)<=set(evidence['contact_eligible_triangle_ids']):
        raise ValueError('Unique ordered exterior source triangles required')
    geometry=json.loads(gzip.decompress(geometry_bytes))
    p=np.asarray(geometry['positions'],float).reshape(-1,3)
    f=np.asarray(geometry['indices'],int).reshape(-1,3)[ids]
    area=float(np.linalg.norm(np.cross(p[f[:,1]]-p[f[:,0]],p[f[:,2]]-p[f[:,0]]),axis=1).sum()/2)
    if not area>0:raise ValueError('Positive patch area required')
    return dict(area_m2=area,support_identity=hashlib.sha256(json.dumps([support,ids],sort_keys=True).encode()).hexdigest(),
                source_support=support,triangle_ids=ids)


def bind_inventory(*,owners,owner_ids,epoch,area_m2,support_identity,epithelial_depth_m,
                   water_fraction,cellular_fraction,film_depth_m,molar_mass_g_mol,
                   capacitance_f,initial_potential_v,temperature_k,prior_source,
                   countercharge_owner,initial_field_energy_owner,heat_owner):
    """Partition one snapshot by homogeneous-fluid engineering priors.

    Cell volume=A*depth*water*cellular_fraction; basal volume uses its complement.
    Apical volume=A*film_depth belongs to an explicitly supplied external film.
    All supplied species are partitioned. Missing species remain unobserved;
    a present unknown mass cannot be silently converted to zero. The countercharge
    is an unresolved electrical prior, NOT an inferred mass of another ion.
    """
    for x in (epoch,support_identity,prior_source,countercharge_owner,initial_field_energy_owner,heat_owner):
        if not isinstance(x,str) or not x.strip():raise ValueError('Explicit source, epoch and energy/charge owners required')
    p=np.array([area_m2,epithelial_depth_m,water_fraction,cellular_fraction,film_depth_m],float)
    if not np.isfinite(p).all() or (p<=0).any() or water_fraction>1 or cellular_fraction>=1:
        raise ValueError('Invalid layer geometry prior')
    mm=np.array(molar_mass_g_mol,float)
    if mm.shape!=(3,) or not np.isfinite(mm).all() or (mm<=0).any():raise ValueError('Explicit positive ion molar masses required')
    if len(owner_ids)!=3 or len(set(owner_ids))!=3:raise ValueError('Three distinct owning pools required')
    volumes=area_m2*np.array([film_depth_m,epithelial_depth_m*water_fraction*cellular_fraction,epithelial_depth_m*water_fraction*(1-cellular_fraction)])
    represented={};complement={};concentrations=[]
    for i,key in enumerate(owner_ids):
        owner=owners[key]
        expected=('external_surface_film','native_skin_intracellular','native_skin_extracellular')[i]
        if owner.get('accounting_owner') is not True or owner.get('origin')!=expected:
            raise ValueError('Actual owning leaf with correct compartment origin required')
        volume=owner['volume_m3']
        if not np.isfinite(volume) or not volume>0 or volumes[i]>volume:raise ValueError('Represented volume exceeds owner')
        masses=owner['mass_g']
        if any(v is None or isinstance(v,bool) or not np.isfinite(v) or v<0 for v in masses.values()):raise ValueError('Known finite nonnegative species masses required')
        if any(masses.get(s,0)<=0 for s in IONS):raise ValueError('All three ion inventories required')
        fraction=float(volumes[i]/volume)
        allocated={s:float(m*fraction) for s,m in masses.items()}
        represented[key]=dict(volume_m3=float(volumes[i]),fraction=fraction,mass_g=allocated)
        complement[key]=dict(volume_m3=float(volume-volumes[i]),mass_g={s:float(m-allocated[s]) for s,m in masses.items()})
        concentrations.append([allocated[s]/mm[j]/volumes[i] for j,s in enumerate(IONS)])
    patch=EpithelialPatch(volumes,concentrations,capacitance_f,np.zeros((3,3)),initial_potential_v,temperature_k,
        dict(preparation='native-snapshot-derived offline patch',parameters=prior_source))
    fixed=patch.capacitance_matrix@patch.initial_potential_v-F*(patch.initial_moles@Z)
    return dict(patch=patch,owner_ids=list(owner_ids),epoch=epoch,support_identity=support_identity,
        owners=copy.deepcopy(owners),represented=represented,complement=complement,molar_mass_g_mol=mm,
        fixed_countercharge_c=fixed,countercharge_owner=countercharge_owner,
        initial_field_energy_owner=initial_field_energy_owner,heat_owner=heat_owner,
        initial_field_energy_j=patch.free_energy_change(np.zeros((3,3)))[1],
        native_commit_available=False,ownership='alternative shadow partition, not additional native pools')


def run_passive(binding,conductance_s,duration_s,samples=11):
    """Offline passive evolution and proposed owner reassembly, never a commit.

    No imposed active flux is exposed: ATP/native energetic supply is unresolved.
    Capacitor energy already belongs to the declared initial-field owner.
    All transport is integrated by the existing instantaneous-force operator;
    this does not validate the prior conductances or native runtime coupling.
    """
    p=binding['patch']
    patch=EpithelialPatch(p.volumes_m3,p.concentrations_mol_m3,p.capacitance_f,
        conductance_s,p.initial_potential_v,p.temperature_k,p.provenance)
    r=simulate(patch,duration_s,samples)
    delta=np.asarray(r['moles'][-1])-p.initial_moles
    proposed={k:dict(binding['owners'][k]['mass_g']) for k in binding['owner_ids']}
    for i,k in enumerate(binding['owner_ids']):
        for j,s in enumerate(IONS):
            proposed[k][s]=binding['complement'][k]['mass_g'][s]+float(np.asarray(r['moles'][-1])[i,j]*binding['molar_mass_g_mol'][j])
    return dict(epoch=binding['epoch'],support_identity=binding['support_identity'],
        proposed_owner_mass_g=proposed,species_delta_total_mol=delta.sum(axis=0).tolist(),
        fixed_countercharge_c=binding['fixed_countercharge_c'].tolist(),
        tep_v=r['tep_v'][-1],basal_membrane_v=r['basal_membrane_v'][-1],apical_membrane_v=r['apical_membrane_v'][-1],
        heat_credit_j=r['dissipated_energy_j'][-1],heat_owner=binding['heat_owner'],
        electrical_energy_j=r['electrical_energy_j'][-1],external_active_work_j=0.,
        audit=r['audit'],native_commit_available=False,
        limitation='Alternative terminal snapshot; does not reserve native pools, resolve ATP or advance native time')


def native_skin_owner(observed, *, compartment):
    """Adapt NativeTissueExchange owning readout without importing its runtime.

    Regional EC leaves are accepted only when explicitly supplied as accounting
    owners. Their geometry allocation remains a separate source mapping prior.
    The observer's mass list is partial, not a claim that other species are absent.
    """
    if compartment not in ('intracellular','extracellular') or observed.get('accounting_owner') is not True:
        raise ValueError('Explicit native Skin owning IC or EC observation required')
    name=observed.get('native_owner','')
    if not (name=='Skin.'+compartment or (compartment=='extracellular' and name in ('Skin.region_a.extracellular','Skin.region_b.extracellular','Skin.residual.extracellular'))):
        raise ValueError('Unrecognized native Skin owner identity')
    return dict(origin='native_skin_'+compartment,accounting_owner=True,
        native_owner=name,volume_m3=observed['volume_ml']*1e-6,mass_g=copy.deepcopy(observed['mass_g']),
        inventory_scope='all supplied observed species; unreported species unresolved')
