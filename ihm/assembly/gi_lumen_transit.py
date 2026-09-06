"""Conservative aqueous SI→colon→rectum→fecal-boundary transit transaction.

No motility rate, epithelial absorption, secretion, water source or native pool
is created. The caller supplies a measured/declared volume request per edge.
All dissolved species follow each donor's well-mixed aqueous fraction. Dry
residue stays put; particulate/solid transit needs a separate physical model.
"""
import numpy as np


def transit(water_ml,solute_mol,requested_ml,*,provenance):
    """Return a new snapshot plus paired edge payloads and external fecal loss.

    Each edge uses the starting donor snapshot, so a single call cannot carry
    a newly received packet through several regions. An eventual native adapter
    must commit these deltas atomically once to existing pools; this function
    never owns persistent duplicates. No wall absorption is implied by transit.
    """
    water=np.array(water_ml,float,copy=True);mass=np.array(solute_mol,float,copy=True);request=np.array(requested_ml,float,copy=True)
    if water.shape!=(3,) or request.shape!=(3,) or mass.ndim!=2 or mass.shape[0]!=3 or mass.shape[1]<1:
        raise ValueError('Require three lumen regions and at least one declared dissolved species')
    if any(not np.isfinite(a).all() or (a<0).any() for a in (water,mass,request)):
        raise ValueError('Finite nonnegative snapshots and requests required')
    if not isinstance(provenance,str) or not provenance.strip():raise ValueError('Request provenance required')
    moved_water=np.minimum(water,request)
    fraction=np.divide(moved_water,water,out=np.zeros(3),where=water>0)
    moved_mass=mass*fraction[:,None]
    remaining_water=water-moved_water;remaining_mass=mass-moved_mass
    remaining_water[1:]+=moved_water[:-1];remaining_mass[1:]+=moved_mass[:-1]
    return dict(water_ml=remaining_water,solute_mol=remaining_mass,
                edge_water_ml=moved_water,edge_solute_mol=moved_mass,
                fecal_water_ml=float(moved_water[-1]),fecal_solute_mol=moved_mass[-1].copy(),
                request_provenance=provenance,native_stores_mutated=False,
                scope='prescribed-volume aqueous transit only; no motility or absorption prediction')
