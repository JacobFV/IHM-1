"""Stateless GI native-boundary transaction planner; never owns native inventory.

Owner order: existing SI chyme, proposed colon lumen, proposed rectum lumen,
existing large-intestine vascular pool, external cumulative fecal ledger.
Ion order Na,K,Cl. Actual native scalar bindings and electrical/chemical paths
are not supplied by this module. A plan is reviewable, not permission to commit.
"""
import hashlib
import json
import numpy as np
from .gi_lumen_transit import transit
from .gi_colon_perfusion import exchange


def _snapshot(snapshot):
    if snapshot.get('phase')!='after_native_SI_absorption':raise ValueError('Bridge must follow existing native absorption once')
    owners=snapshot.get('owner_ids',[])
    if len(owners)!=5 or len(set(owners))!=5 or any(not isinstance(o,str) or not o for o in owners):raise ValueError('Five distinct authoritative owners required')
    if not isinstance(snapshot.get('epoch'),int):raise ValueError('Native timestep epoch required')
    if not snapshot.get('current_ledger_owner'):raise ValueError('Explicit current/charge residual ledger owner required')
    w=np.array(snapshot['water_ml'],float);n=np.array(snapshot['ion_mol'],float)
    if w.shape!=(5,) or n.shape!=(5,3) or not np.isfinite(w).all() or not np.isfinite(n).all() or (w<0).any() or (n<0).any():raise ValueError('Finite nonnegative owner inventories required')
    serial=dict(epoch=snapshot['epoch'],phase=snapshot['phase'],owner_ids=owners,water_ml=w.tolist(),ion_mol=n.tolist(),current_ledger_owner=snapshot['current_ledger_owner'])
    digest=hashlib.sha256(json.dumps(serial,sort_keys=True,allow_nan=False).encode()).hexdigest()
    return w,n,digest


def prepare(snapshot,transit_requested_ml,duration_s,colon_condition):
    """Plan donor-arbitrated signed absorption followed by aqueous transit.

    Priority is an explicit engineered splitting choice: colon wall exchange
    first, transit from remaining lumen inventories second. Received upstream
    fluid cannot traverse a second lumen edge in this operation. SI wall uptake
    is already owned by the existing native GI method and is never repeated.
    """
    w,n,digest=_snapshot(snapshot);oldw=w.copy();oldn=n.copy()
    wall=exchange(w[[1,3]],n[[1,3]],duration_s,colon_condition)
    w[[1,3]]=wall['water_ml'];n[[1,3]]=wall['ion_mol']
    flow=transit(w[:3],n[:3],transit_requested_ml,provenance='caller-supplied native transit request; engineered split priority wall then transit')
    w[:3]=flow['water_ml'];n[:3]=flow['solute_mol']
    w[4]+=flow['fecal_water_ml'];n[4]+=flow['fecal_solute_mol']
    dw=w-oldw;dn=n-oldn
    if not np.allclose(dw.sum(),0,rtol=0,atol=1e-10) or not np.allclose(dn.sum(axis=0),0,rtol=0,atol=1e-14):raise ArithmeticError('Transaction conservation failure')
    return dict(snapshot_sha256=digest,epoch=snapshot['epoch'],owner_ids=list(snapshot['owner_ids']),
                water_delta_ml=dw.tolist(),ion_delta_mol=dn.tolist(),
                wall_signed_payload_ml_mol=wall['signed_lumen_to_serosa'].tolist(),
                transit_water_payload_ml=flow['edge_water_ml'].tolist(),transit_ion_payload_mol=flow['edge_solute_mol'].tolist(),
                colon_uncompensated_charge_mol=wall['uncompensated_charge_mol'],
                current_ledger_owner=snapshot['current_ledger_owner'],
                split_order='native SI uptake already committed to snapshot; colon wall; snapshot lumen transit',
                empirical_transfer_prior=colon_condition.get('transfer_prior'),native_stores_mutated=False,native_commit_ready=False,
                missing_native_interfaces=['Distinct colon/rectum fluid nodes and substance quantities with serialization',
                    'Single volume mutation owner across GI PreProcess and circuit NextVolume',
                    'H/bicarbonate/electrical current path for empirical Na+K-Cl residual',
                    'Native-bound complete species payload including nutrients, drugs and nonionic solutes',
                    'Native atomic scalar writer and consumed-epoch transaction guard'])


def validate_commit(snapshot,transaction):
    """Read-only precondition check, not a native scalar writer.

    A future writer must validate inside the native mutation boundary, atomically
    mark the epoch consumed, and update every volume/species/ledger once. Merely
    invoking this Python check cannot make a subsequent native write atomic.
    """
    _,_,digest=_snapshot(snapshot)
    if transaction.get('snapshot_sha256')!=digest or transaction.get('epoch')!=snapshot['epoch'] or transaction.get('owner_ids')!=snapshot['owner_ids']:
        raise ValueError('Stale snapshot or changed owner bindings; rebuild transaction')
    return True
