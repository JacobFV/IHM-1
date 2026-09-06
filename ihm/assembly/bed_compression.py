"""Finite measured mattress compression curves in series with bounded skin.

The bed table is a digitized nominal-stress interpretation of Hong et al. 2022,
not a universal foam law. Skin and bed compression remain separate observables.
"""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np


def load_bed(root, material):
    directory=Path(root)/'data/research/bed_material';manifest_path=directory/'manifest.json'
    manifest=json.loads(manifest_path.read_text());table=directory/manifest['digitization']['path']
    digest=hashlib.sha256(table.read_bytes()).hexdigest()
    if digest!=manifest['derived_files'][table.name]['sha256']:raise ValueError('Bed curve identity mismatch')
    rows=[r for r in csv.DictReader(table.open()) if r['material']==material]
    if not rows:raise ValueError('Unknown retained mattress curve')
    strain=np.array([float(r['compressive_strain']) for r in rows]);pressure=np.array([float(r['stress_pa']) for r in rows])
    if not np.isfinite(strain).all() or not np.isfinite(pressure).all() or strain[0]!=0 or pressure[0]!=0 or (np.diff(strain)<=0).any() or (np.diff(pressure)<0).any():raise ValueError('Invalid monotone compression table')
    return dict(material=material,thickness_m=manifest['engineering_law_candidate']['thickness_m'],
        strain=strain.tolist(),pressure_pa=pressure.tolist(),table_sha256=digest,
        evidence_manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        source_doi=manifest['doi'],stress_basis=manifest['engineering_law_candidate']['nominal_stress_interpretation'],
        normal_rate_law='conservative elastic curve; bed damping/hysteresis uncalibrated')


def skin_pressure(compression,skin):
    stretch=1-np.asarray(compression)/skin['total_layer_thickness_m']
    return -(skin['shear_modulus_pa']*(stretch-1/stretch)+skin['lame_lambda_pa']*np.log(stretch)/stretch)


def bed_energy(strain,bed):
    x=np.asarray(bed['strain']);y=np.asarray(bed['pressure_pa']);strain=np.asarray(strain)
    index=np.clip(np.searchsorted(x,strain,side='right')-1,0,len(x)-2)
    integral=np.r_[0,np.cumsum(.5*(y[:-1]+y[1:])*np.diff(x))]
    offset=strain-x[index];slope=(y[index+1]-y[index])/(x[index+1]-x[index])
    return bed['thickness_m']*(integral[index]+y[index]*offset+.5*slope*offset**2)


def series_response(approach_m,skin,bed):
    approach=np.maximum(0.,np.asarray(approach_m,float))
    if not np.isfinite(approach).all():raise ValueError('Finite contact approach required')
    hs=skin['total_layer_thickness_m'];hb=bed['thickness_m']
    if hs<=0 or hb<=0:raise ValueError('Positive layer thickness required')
    skin_max=hs*(1-skin['minimum_thickness_ratio']);bed_max=hb*bed['strain'][-1]
    low=np.maximum(0.,approach-bed_max);high=np.minimum(skin_max,approach)
    if (low>high+1e-12).any():raise ValueError('Combined skin/bed compression outside retained domains')
    def residual(delta):return skin_pressure(delta,skin)-np.interp((approach-delta)/hb,bed['strain'],bed['pressure_pa'])
    if (residual(low)>1e-7).any() or (residual(high)<-1e-7).any():raise ValueError('Equal-pressure series solution exceeds retained material domain')
    for _ in range(50):
        middle=(low+high)/2;positive=residual(middle)>0
        high=np.where(positive,middle,high);low=np.where(positive,low,middle)
    delta_skin=(low+high)/2;delta_bed=approach-delta_skin
    pressure=skin_pressure(delta_skin,skin);stretch=1-delta_skin/hs;log=np.log(stretch)
    energy_skin=hs*(.5*skin['shear_modulus_pa']*(stretch**2-1)-skin['shear_modulus_pa']*log+.5*skin['lame_lambda_pa']*log**2)
    energy_bed=bed_energy(delta_bed/hb,bed)
    return dict(pressure_pa=pressure,skin_indentation_m=delta_skin,bed_indentation_m=delta_bed,
                skin_energy_j_m2=energy_skin,bed_energy_j_m2=energy_bed,
                total_energy_j_m2=energy_skin+energy_bed,
                pressure_balance_residual_pa=pressure-np.interp(delta_bed/hb,bed['strain'],bed['pressure_pa']))


def maximum_approach(skin,bed):
    # Domain boundary is whichever constituent reaches its retained range first.
    hs=skin['total_layer_thickness_m'];skin_max=hs*(1-skin['minimum_thickness_ratio'])
    pressure=min(float(skin_pressure(skin_max,skin)),bed['pressure_pa'][-1])
    lo,hi=0.,skin_max
    for _ in range(60):
        middle=(lo+hi)/2
        if skin_pressure(middle,skin)>pressure:hi=middle
        else:lo=middle
    # At plateau pressure choose largest admissible bed compression, preserving
    # every measured strain interval; do not delete plateau rows.
    bed_strain=float(np.interp(pressure,bed['pressure_pa'],bed['strain']))
    return float((lo+hi)/2+bed['thickness_m']*bed_strain)
