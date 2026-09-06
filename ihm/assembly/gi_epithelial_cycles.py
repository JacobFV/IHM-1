"""Finite, uncalibrated GI cycle preparation, independent of native physiology.

Species are Na+, K+, Cl-, glucose and an explicitly NEUTRAL AA. Native lumped
AminoAcids is not mapped here. Requests are finite cycle extents, not kinetics.
"""
import numpy as np
from .epithelial_electrodiffusion import R, F, INCIDENCE

SPECIES = ('Na', 'K', 'Cl', 'glucose', 'neutral_AA')
Z = np.array([1., 1., -1., 0., 0.])


def _cycle(transfers):
    result=np.zeros((3,5))
    for source,dest,species,coefficient in transfers:
        result[source,species]-=coefficient
        result[dest,species]+=coefficient
    return result

CYCLES = {
    'sglt1': _cycle([(0,1,0,2),(0,1,3,1)]),
    'b0at1': _cycle([(0,1,0,1),(0,1,4,1)]),
    'nak_pump': _cycle([(1,2,0,3),(2,1,1,2)]),
    'k_recycling': _cycle([(1,2,1,1)]),
    'glucose_efflux': _cycle([(1,2,3,1)]),
    'aa_efflux': _cycle([(1,2,4,1)]),
    'paracellular_na': _cycle([(0,2,0,1)]),
    'paracellular_cl': _cycle([(0,2,2,1)]),
}
for _v in CYCLES.values(): _v.setflags(write=False)


class GIPreparation:
    """Ideal fixed-volume pools with two membrane capacitors and an ATP work port.

    ATP is a finite external cycle allowance with a supplied hydrolysis free
    energy, not an intracellular ATP/ADP/Pi chemical model. No native metabolic
    demand is charged. Pump work minus free-energy increase is dissipated heat
    in this isothermal ledger; it is not exported to a native thermal object.
    """
    def __init__(self,moles,volumes_m3,capacitance_f,initial_potential_v,
                 temperature_k,provenance,*,atp_cycles_mol,atp_free_energy_j_mol):
        def array(value,shape,positive=False):
            a=np.array(value,dtype=float,copy=True)
            if a.shape!=shape or not np.isfinite(a).all() or (positive and (a<=0).any()):
                raise ValueError('Invalid finite preparation array')
            return a
        self.moles=array(moles,(3,5))
        if (self.moles<0).any():raise ValueError('Negative initial inventory')
        self.volumes=array(volumes_m3,(3,),True)
        cap=array(capacitance_f,(2,),True)
        self.phi0=array(initial_potential_v,(3,))
        if self.phi0[2]!=0:raise ValueError('Serosal reference must be zero')
        for val in (temperature_k,atp_free_energy_j_mol):
            if not np.isfinite(val) or val<=0:raise ValueError('Positive temperature and ATP free energy required')
        if not np.isfinite(atp_cycles_mol) or atp_cycles_mol<0:raise ValueError('Invalid ATP allowance')
        if not provenance.get('preparation') or not provenance.get('parameters'):raise ValueError('Explicit provenance required')
        self.provenance=dict(provenance);self.temperature=temperature_k
        self.atp_cycles_mol=atp_cycles_mol;self.atp_free_energy_j_mol=atp_free_energy_j_mol
        b=INCIDENCE[:,:2];self.capacitance=(b*cap)@b.T
        self.initial=self.moles.copy();self.work_j=0.;self.dissipation_j=0.;self.history=[]
        self.initial_energy=self.free_energy(self.initial)

    def electrical_state(self,n=None):
        n=self.moles if n is None else n
        q=self.capacitance@self.phi0+F*((n-self.initial)@Z)
        phi=np.r_[np.linalg.solve(self.capacitance[:2,:2],q[:2]),0.]
        return q,phi

    def free_energy(self,n):
        # lim n log(n/V)=0 admits exact donor exhaustion.
        positive=n>0
        concentrations=n/self.volumes[:,None]
        chemical=R*self.temperature*np.sum(n[positive]*(np.log(concentrations[positive])-1))
        _,phi=self.electrical_state(n)
        return float(chemical+.5*phi@self.capacitance@phi)

    def marginal_cost(self, n, nu):
        """Directional derivative of convex ideal-mixture plus capacitor energy.

        At a zero pool, creation has -infinite and withdrawal +infinite
        chemical marginal cost. Zero-inventory donors are removed by the
        extent budget before this derivative is evaluated.
        """
        used=nu!=0
        if (n[used]<=0).any():
            if ((n<=0)&(nu<0)).any():return float('inf')
            return -float('inf')
        _,phi=self.electrical_state(n)
        mu=R*self.temperature*np.log(n[used]/np.broadcast_to(self.volumes[:,None],n.shape)[used])
        electrical=np.broadcast_to(phi[:,None],n.shape)[used]*np.broadcast_to(Z,n.shape)[used]*F
        return float(np.sum(nu[used]*(mu+electrical)))

    def apply(self,pathway,requested_extent_mol):
        """Transactional donor-limited extent. Uphill unpowered requests fail.

        Convex marginal-affinity stopping prevents passage beyond equilibrium.
        An imposed extent is still not a fitted rate. Reverse passive/cotransport extents are allowed;
        reverse ATP synthesis requires a separate explicit model and is rejected.
        """
        if pathway not in CYCLES or not np.isfinite(requested_extent_mol):raise ValueError('Unknown pathway or invalid extent')
        pump=pathway=='nak_pump'
        if pump and requested_extent_mol<0:raise ValueError('ATP synthesis not modeled')
        sign=1 if requested_extent_mol>=0 else -1
        nu=CYCLES[pathway]*sign
        donor=nu<0
        extent=min(abs(requested_extent_mol),float(np.min(self.moles[donor]/-nu[donor])))
        if pump:extent=min(extent,self.atp_cycles_mol)
        trial=self.moles+extent*nu
        while (trial<0).any():
            extent=np.nextafter(extent,0.);trial=self.moles+extent*nu
        supplied=self.atp_free_energy_j_mol if pump else 0.
        if extent>0:
            initial_cost=self.marginal_cost(self.moles,nu)
            if initial_cost>supplied:
                raise ValueError('Initial marginal cost exceeds supplied driving work')
            if self.marginal_cost(trial,nu)>supplied:
                # Along a fixed stoichiometric column, ideal chemical energy
                # and positive capacitor energy are convex. Thus the marginal
                # cost is monotone and its first crossing is the only stop.
                lower,upper=0.,extent
                for _ in range(100):
                    mid=lower+(upper-lower)/2
                    if mid==lower or mid==upper:break
                    if self.marginal_cost(self.moles+mid*nu,nu)<=supplied:lower=mid
                    else:upper=mid
                extent=lower;trial=self.moles+extent*nu
        dg=self.free_energy(trial)-self.free_energy(self.moles)
        work=extent*self.atp_free_energy_j_mol if pump else 0.
        if dg>work:raise ValueError('Requested finite extent is uphill beyond available work; reduce extent or supply an explicit driver')
        self.moles=trial;self.work_j+=work;self.dissipation_j+=work-dg
        if pump:self.atp_cycles_mol-=extent
        row=dict(pathway=pathway,extent_mol=sign*extent,free_energy_change_j=dg,
                 external_work_j=work,dissipation_j=work-dg,atp_cycles_mol=extent if pump else 0.,
                 pump_net_outward_charge_c=F*extent if pump else 0.,
                 final_marginal_cost_j_mol=self.marginal_cost(trial,nu) if extent else None)
        self.history.append(row)
        return row

    def audit(self):
        return dict(maximum_species_residual_mol=float(np.max(np.abs((self.moles-self.initial).sum(axis=0)))),
                    global_charge_change_c=float(F*np.sum((self.moles-self.initial)@Z)),
                    energy_residual_j=self.free_energy(self.moles)-self.initial_energy+self.dissipation_j-self.work_j,
                    external_work_j=self.work_j,dissipation_j=self.dissipation_j,
                    native_inventory_mutated=False,calibrated_kinetics=False)
