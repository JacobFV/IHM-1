"""Finite three-compartment epithelial ion/charge/energy primitive.

Compartments: apical, cell, basal. Ions: Na+, K+, Cl-. Transport interfaces:
apical->cell, cell->basal, apical->basal (paracellular). All parameters must be
supplied with preparation/prior provenance. No defaults imply human calibration.
No native simulator objects are accepted or mutated.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

R = 8.314462618
F = 96485.33212
Z = np.array([1., 1., -1.])
# Each column removes one mole at its source and adds one at its destination.
INCIDENCE = np.array([[-1., 0., -1.], [1., -1., 0.], [0., 1., 1.]])


def _array(value, shape, name, positive=False, nonnegative=False):
    a = np.array(value, dtype=float, copy=True)
    if a.shape != shape or not np.isfinite(a).all() or (positive and (a <= 0).any()) or (nonnegative and (a < 0).any()):
        raise ValueError('Invalid '+name)
    a.setflags(write=False)
    return a


@dataclass(frozen=True)
class EpithelialPatch:
    volumes_m3: object
    concentrations_mol_m3: object
    capacitance_f: object
    conductance_s: object
    initial_potential_v: object
    temperature_k: float
    provenance: dict
    active_flux_mol_s: object = None

    def __post_init__(self):
        for key, shape, positive, nonnegative in (
            ('volumes_m3',(3,),True,False), ('concentrations_mol_m3',(3,3),True,False),
            ('capacitance_f',(2,),True,False), ('conductance_s',(3,3),False,True),
            ('initial_potential_v',(3,),False,False)):
            object.__setattr__(self,key,_array(getattr(self,key),shape,key,positive,nonnegative))
        object.__setattr__(self,'active_flux_mol_s',_array(np.zeros((3,3)) if self.active_flux_mol_s is None else self.active_flux_mol_s,(3,3),'active_flux_mol_s'))
        if isinstance(self.temperature_k,bool) or not np.isfinite(self.temperature_k) or self.temperature_k <= 0:
            raise ValueError('Positive temperature required')
        if self.initial_potential_v[2] != 0:
            raise ValueError('Basal potential must be grounded at zero')
        if not isinstance(self.provenance,dict) or not self.provenance.get('preparation') or not self.provenance.get('parameters'):
            raise ValueError('Explicit preparation and parameter provenance required')
        object.__setattr__(self,'provenance',dict(self.provenance))

    @property
    def capacitance_matrix(self):
        b=INCIDENCE[:,:2]
        return (b*self.capacitance_f)@b.T

    @property
    def initial_moles(self):
        return self.concentrations_mol_m3*self.volumes_m3[:,None]

    def electrical_state(self, delta_moles):
        """Fixed countercharge sets initial Q; subsequent charge is ionic transfer."""
        delta=_array(delta_moles,(3,3),'delta_moles')
        q=self.capacitance_matrix@self.initial_potential_v+F*(delta@Z)
        phi=np.r_[np.linalg.solve(self.capacitance_matrix[:2,:2],q[:2]),0.]
        return q,phi

    def transport(self, delta_moles):
        """Return species rate, interface fluxes and thermodynamic powers.

        Active flux is an externally powered signed transfer, not ATP kinetics.
        Its work is signed flux times the destination-source electrochemical
        potential. Negative active work means work returned to the actuator.
        """
        n=self.initial_moles+delta_moles
        if (n <= 0).any() or not np.isfinite(n).all():
            raise ValueError('Finite reservoir exhausted; reduce flux/horizon or supply depletion kinetics')
        _,phi=self.electrical_state(delta_moles)
        mu=R*self.temperature_k*np.log(n/self.volumes_m3[:,None])+F*phi[:,None]*Z
        downhill=-INCIDENCE.T@mu
        passive=self.conductance_s/(F*Z)**2*downhill
        flux=passive+self.active_flux_mol_s
        return dict(rate_mol_s=INCIDENCE@flux,passive_flux_mol_s=passive,
                    active_flux_mol_s=self.active_flux_mol_s,
                    dissipation_w=float(np.sum(passive*downhill)),
                    active_power_w=float(-np.sum(self.active_flux_mol_s*downhill)))

    def free_energy_change(self, delta):
        n0=self.initial_moles
        if ((n0+delta)<=0).any():raise ValueError('Nonpositive inventory')
        chemical=R*self.temperature_k*np.sum((n0+delta)*np.log1p(delta/n0)+delta*(np.log(self.concentrations_mol_m3)-1))
        _,phi=self.electrical_state(delta)
        electrical=.5*phi@self.capacitance_matrix@phi
        initial=.5*self.initial_potential_v@self.capacitance_matrix@self.initial_potential_v
        return float(chemical+electrical-initial),float(electrical)


def simulate(patch, duration_s, samples=51):
    """Small bounded independent preparation; integrate transfer extents.

    Extents directly preserve each ionic species. Output includes capacitor
    energy, electrochemical dissipation and signed active source work.
    """
    if isinstance(duration_s,bool) or not np.isfinite(duration_s) or not 0 < duration_s <= 10 or not isinstance(samples,int) or isinstance(samples,bool) or not 2 <= samples <= 2001:
        raise ValueError('Require duration (0,10] seconds and 2..2001 samples')
    times=np.linspace(0,duration_s,samples)
    # States: 9 net transport extents, dissipated energy, active work.
    def rhs(_,state):
        delta=INCIDENCE@state[:9].reshape(3,3)
        t=patch.transport(delta)
        return np.r_[(t['passive_flux_mol_s']+t['active_flux_mol_s']).ravel(),t['dissipation_w'],t['active_power_w']]
    sol=solve_ivp(rhs,(0,duration_s),np.zeros(11),t_eval=times,rtol=2e-9,atol=np.r_[np.full(9,1e-25),1e-22,1e-22],max_step=duration_s/100)
    if not sol.success:raise RuntimeError(sol.message)
    extents=sol.y[:9].T.reshape(-1,3,3)
    delta=np.array([INCIDENCE@e for e in extents])
    n=patch.initial_moles+delta
    if (n<=0).any():raise ValueError('Reservoir depletion')
    qp=[patch.electrical_state(d) for d in delta];q=np.array([x[0] for x in qp]);phi=np.array([x[1] for x in qp])
    energies=np.array([patch.free_energy_change(d) for d in delta]);loss=sol.y[9];work=sol.y[10]
    residual=energies[:,0]+loss-work
    scale=max(float(np.max(np.abs(energies[:,0]))),float(np.max(np.abs(loss))),float(np.max(np.abs(work))),1e-30)
    return dict(time_s=times.tolist(),moles=n.tolist(),transport_extents_mol=extents.tolist(),potential_v=phi.tolist(),
        capacitor_charge_c=q.tolist(),electrical_energy_j=energies[:,1].tolist(),free_energy_change_j=energies[:,0].tolist(),
        dissipated_energy_j=loss.tolist(),active_work_j=work.tolist(),tep_v=(-phi[:,0]).tolist(),
        basal_membrane_v=phi[:,1].tolist(),apical_membrane_v=(phi[:,1]-phi[:,0]).tolist(),
        fixed_countercharge_c=(patch.capacitance_matrix@patch.initial_potential_v-F*(patch.initial_moles@Z)).tolist(),
        audit=dict(max_species_residual_mol=float(np.max(np.abs(delta.sum(axis=1)))),
            max_charge_residual_c=float(np.max(np.abs(q-phi@patch.capacitance_matrix.T))),
            max_energy_residual_j=float(np.max(np.abs(residual))),relative_energy_residual=float(np.max(np.abs(residual))/scale)),
        ownership=dict(inventories='independent epithelial preparation; includes finite apical and basal baths',native_blood_mutated=False,
            external_source='prescribed active transfer exchanges work, not ion mass, with an explicitly unmodelled actuator'),
        provenance=dict(patch.provenance),limitations=['Ideal dilute activities; fixed volumes, temperature and capacitances; no water/osmotic mechanics.',
            'Linear electrochemical transport coefficients and prescribed active flux are unidentified unless separately evidenced.',
            'No ATP reservoir or channel gating; externally powered transport is not an ATP consumption prediction.',
            'One epithelial cell layer, not stratified intact epidermis; no nerve, perfusion, migration or regeneration coupling.'])


def prior_ensemble(patches,duration_s,samples=51):
    """Equal treatment of explicitly supplied scenarios; no posterior weights."""
    patches=list(patches)
    if not 2 <= len(patches) <= 32:raise ValueError('Supply 2..32 explicit prior scenarios')
    if len({p.provenance['preparation'] for p in patches}) != 1:raise ValueError('Do not pool different preparations')
    runs=[simulate(p,duration_s,samples) for p in patches]
    tep=np.array([r['tep_v'] for r in runs])
    return dict(member_count=len(runs),time_s=runs[0]['time_s'],tep_min_v=tep.min(axis=0).tolist(),tep_max_v=tep.max(axis=0).tolist(),
        uncertainty_kind='explicit_prior_scenario_envelope_not_posterior',members=runs)


def fit_voltage_observations(rows):
    """Identify only a same-preparation/quantity voltage location, never RC.

    An observed spread is descriptive; it is not automatically a likelihood
    uncertainty or an independent-donor standard error.
    """
    rows=list(rows)
    if not rows or any(not r.get('source_id') or not r.get('quantity') or not r.get('preparation') for r in rows):raise ValueError('Source/quantity/preparation required')
    if len({(r['quantity'],r['preparation']) for r in rows})!=1:raise ValueError('Mixed observables/preparations cannot be fitted together')
    values=np.array([r['value_v'] for r in rows],float)
    if not np.isfinite(values).all():raise ValueError('Finite voltage observations required')
    return dict(voltage_location_v=float(values.mean()),identified_parameters=['voltage_location_v'],
        observed_range_v=[float(values.min()),float(values.max())],n_rows=len(rows),
        parameter_uncertainty=None,uncertainty_note='No independent-donor or measurement-error model supplied',
        unidentifiable_parameters=['conductances','capacitances','active_fluxes','intracellular_concentrations'],
        quantity=rows[0]['quantity'],preparation=rows[0]['preparation'],sources=sorted({r['source_id'] for r in rows}))


def condition_tep_magnitude(patch, observation, *, basal_positive):
    """Condition only the initial observed magnitude, with explicit polarity.

    A source whose voltage wiring sign is unresolved can constrain magnitude.
    Polarity and cell Vm then remain explicit priors, not fitted measurements.
    Initial fixed countercharge is re-established, not an instantaneous flux.
    """
    from dataclasses import replace
    value=observation.get('value_v')
    if observation.get('quantity')!='tep_magnitude' or observation.get('preparation')!=patch.provenance['preparation'] or not observation.get('source_id'):
        raise ValueError('Matching preparation and sourced TEP magnitude required')
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not np.isfinite(value) or value<0 or not isinstance(basal_positive,bool):
        raise ValueError('Nonnegative finite magnitude and explicit Boolean polarity required')
    phi=patch.initial_potential_v.copy();phi[0]=-value if basal_positive else value
    return replace(patch,initial_potential_v=phi,provenance=dict(patch.provenance,
        initial_tep_magnitude_source=observation['source_id'],polarity_prior='basal_positive' if basal_positive else 'apical_positive',
        conditioning_scope='Initial TEP magnitude only; initial countercharge reset; no current, RC or cell Vm fit'))


def fit_tape_strip_response(rows):
    """Descriptive linear voltage-magnitude response; last strip count held out.

    Repeated interventions on one participant cannot identify population error,
    channel kinetics, conductance or a time constant. Tape count is not time.
    """
    rows=sorted(list(rows),key=lambda r:r['tape_strips'])
    if len(rows)<4:raise ValueError('At least three training states and one held state required')
    x=np.array([r['tape_strips'] for r in rows],float)
    y=np.array([r['value_v'] for r in rows],float)
    if not np.isfinite(np.r_[x,y]).all() or (x<0).any() or (y<0).any() or len(set(x))!=len(x):raise ValueError('Distinct nonnegative strip counts and voltage magnitudes required')
    design=np.c_[np.ones(len(x)-1),x[:-1]]
    beta,_,rank,_=np.linalg.lstsq(design,y[:-1],rcond=None)
    if rank!=2:raise ValueError('Unidentifiable descriptive design')
    prediction=float(beta@[1,x[-1]])
    return dict(intercept_v=float(beta[0]),slope_v_per_strip=float(beta[1]),
        training_strip_counts=x[:-1].tolist(),holdout_strip_count=float(x[-1]),holdout_observed_v=float(y[-1]),
        holdout_prediction_v=prediction,holdout_error_v=float(prediction-y[-1]),
        training_residual_v=(y[:-1]-design@beta).tolist(),identified_parameters=['intercept_v','slope_v_per_strip'],
        rc_parameters_identified=False,population_uncertainty=None,
        limitation='Within-participant descriptive response; no independent donor validation or time constant')
