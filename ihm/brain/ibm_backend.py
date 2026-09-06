"""Execute pinned IBM source through an explicit, bounded materializer/runtime adapter.

The selected runtime equation is dy/dt=H*x-y. A periodic window solves
y_hat=H*x_hat/(i*omega+1). It is not a causal online filter or a calibrated
receptor voltage. No copied IHM transfer formula is used to generate output.
"""
from dataclasses import dataclass
from ihm.brain.active_source import DEFAULT_SOURCE,resolve_source
from pathlib import Path
import hashlib
import inspect
import json
import sys
import numpy as np

DEFAULT_ARTIFACT = Path(__file__).resolve().parents[2] / 'data/derived/canonical/ibm-backend'
PINNED_PACKAGE_SHA256 = '0bcb65606a350d586c1c5fc1634913aee7f90c7de4201eef57770c6af715346f'


@dataclass(frozen=True)
class WindowResult:
    values: np.ndarray
    residual: float
    audit: dict


class IBMBackend:
    def __init__(self, artifact_dir=None, *, source_pin=DEFAULT_SOURCE, root=None):
        # Naming an artifact directory selects that artifact. Only an unaddressed
        # backend takes the active default, so a positional legacy caller keeps
        # reaching the legacy verification instead of a pin-mismatch error.
        if source_pin is DEFAULT_SOURCE and artifact_dir is not None:source_pin=None
        source_pin=resolve_source(root or Path(__file__).resolve().parents[2],source_pin)
        self.source_pin = source_pin
        if source_pin is not None:
            from .candidate import verify_pin
            if artifact_dir is not None and Path(artifact_dir).resolve() != source_pin.artifact_dir:
                raise ValueError('Artifact directory differs from explicit source pin')
            self.identity, snapshots = verify_pin(source_pin)
            self.artifact_dir = source_pin.artifact_dir
            source = self.artifact_dir / 'source'
            expected = source_pin.package_sha256
        else:
            self.artifact_dir = Path(artifact_dir or DEFAULT_ARTIFACT).resolve()
            manifest = self.artifact_dir / 'manifest.json'
            if not manifest.is_file():
                raise FileNotFoundError('IBM artifact missing; run scripts/vendor_ibm_backend.py')
            self.identity = json.loads(manifest.read_text())
            source = self.artifact_dir / 'source'
            actual={p.relative_to(source).as_posix() for p in source.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc'}
            if actual!=set(self.identity['files']):
                raise ValueError('IBM source contains missing or unmanifested files')
            snapshots={}
            for rel, digest in self.identity['files'].items():
                path=(source/rel).resolve()
                if not path.is_relative_to(source):raise ValueError('IBM source path escapes artifact')
                snapshots[rel]=path.read_bytes()
                if hashlib.sha256(snapshots[rel]).hexdigest() != digest:
                    raise ValueError(f'IBM artifact hash mismatch: {rel}')
            expected = hashlib.sha256(json.dumps(self.identity['files'], sort_keys=True).encode()).hexdigest()
            if expected != self.identity['package_sha256']:
                raise ValueError('IBM package manifest identity mismatch')
            if expected != PINNED_PACKAGE_SHA256:
                raise ValueError('IBM artifact differs from the package pinned by this IHM adapter')

        from .source_loader import install_snapshot
        install_snapshot(source,snapshots,expected)
        import ibm
        ibm.load_all()
        self.registry = ibm.REGISTRY

    def identity_summary(self):
        """The law this backend actually executed, independent of any declared path."""
        return {'selection':'legacy independently preserved regional law' if self.source_pin is None
                    else 'active immutable candidate',
                'artifact_dir':str(self.artifact_dir),
                'package_sha256':self.identity['package_sha256'],
                'source_pin':None if self.source_pin is None else self.source_pin.to_dict()}

    @staticmethod
    def basis(n, dt_s):
        from ibm.fields.uncertainty.spectral import TemporalBasis
        if n < 2 or not np.isfinite(dt_s) or dt_s <= 0:
            raise ValueError('window needs n >= 2 and positive finite dt_s')
        return TemporalBasis(n, dt_s)

    def materialize_tactile(self, sites_m, n=512, dt_s=.002, kind='rapid'):
        from ibm.materialize import MaterializationRequest, Window, build
        from ibm.vocabulary import OnSupport, sel
        from .body_supports import body_surface_geometry
        if kind not in ('rapid', 'slow'):
            raise ValueError('tactile materialization supports rapid or slow')
        self.basis(n, dt_s)
        geometry = body_surface_geometry(sites_m)
        request = MaterializationRequest(name='ihm-body-tactile-'+kind,
            targets=(sel('transduction.mechanoreceptor'),),
            regions=(('canonical_body_surface', OnSupport('body_surface')),),
            window=Window(n=n, dt=dt_s), frame='ihm_canonical_body',
            notes='IHM explicit body support subset; omitted dependencies are not executable')
        build_options = {'substrate': False} if self.source_pin is not None and 'substrate' in inspect.signature(build).parameters else {}
        model = build(request, geometry=geometry, strict=False,
                      implementation_override={'transduction': 'mechanoreceptor_'+kind}, **build_options)
        if 'transduction.mechanoreceptor' not in model.components or not model.layout['transduction.mechanoreceptor'].n_sites:
            raise ValueError('IBM materializer produced no receptor sites')
        return model

    def _transfer_spec(self, kind):
        from ibm.processes.effector import twitch_transfer
        names = {'rapid': 'mechanoreceptor_rapid', 'slow': 'mechanoreceptor_slow',
                 'activation': 'twitch_activation_lti'}
        if kind == 'twitch':
            return twitch_transfer, {}, 'ibm.processes.effector.twitch_transfer'
        if kind not in names:
            raise ValueError(f'unsupported IBM transfer kind: {kind}')
        implementation = next(i for i in self.registry.implementations.values() if i.name == names[kind])
        from ibm.forge.priors import median_of
        accepted = inspect.signature(implementation.transfer).parameters
        theta = {k: median_of(v) for k, v in implementation.params.items() if k in accepted}
        return implementation.transfer, theta, names[kind]

    def transfer(self, kind, basis, parameters=None):
        fn, theta, _ = self._transfer_spec(kind)
        supplied = dict(parameters or {})
        invalid = set(supplied) - (set(inspect.signature(fn).parameters) - {'basis'})
        if invalid:
            raise ValueError(f'parameters not accepted by selected IBM transfer: {sorted(invalid)}')
        theta.update(supplied)
        return fn(basis, **theta)

    def run_window(self, inputs, *, kind='rapid', dt_s=.002, parameters=None,
                   delay_s=0., sites_m=None):
        from ibm.runtime import Block, Layout, State, Coupling, Solve, Clamp, solve_window
        from ibm.vocabulary import sel
        from .body_supports import SUPPORT_GAPS
        values = np.asarray(inputs, dtype=float)
        if values.ndim == 1:
            values = values[None, :]
        if values.ndim != 2 or not values.shape[0] or not np.isfinite(values).all():
            raise ValueError('inputs must be finite (sites, samples) trajectories')
        if not np.isfinite(delay_s) or delay_s < 0:
            raise ValueError('delay_s must be finite and nonnegative')
        count, n = values.shape
        basis = self.basis(n, dt_s)
        model = None
        caller_sites = sites_m is not None
        if kind in ('rapid', 'slow'):
            if sites_m is None:
                # Explicit test support. Caller supplies real canonical sites for body integration.
                sites_m = np.column_stack((np.arange(count)*.001, np.zeros((count, 2))))
            model = self.materialize_tactile(sites_m, n, dt_s, kind)
            count_sites = model.layout['transduction.mechanoreceptor'].n_sites
            if count_sites != count:
                raise ValueError('materialized receptor count differs from input count')
        output = 'transduction.mechanoreceptor' if model is not None else 'ihm.ibm.effector_response'
        receptor = model.layout[output] if model is not None else Block(output, 'spectral', count, basis=basis)
        # Normalize selected materialized state into this full-band runtime window.
        layout = Layout((Block('ihm.indentation_um', 'spectral', count, basis=basis),
                         Block(output, 'spectral', receptor.n_sites, basis=basis, support=receptor.support)))
        H = self.transfer(kind, basis, parameters)
        couplings = (Coupling('ihm.selected_ibm_transfer', output, ('ihm.indentation_um',),
                              transfer=H, delay_s=delay_s),
                     Coupling('ihm.explicit_unit_rate_leak', output, (output,), gain=-1.))
        state, report = solve_window(State.zeros(layout), couplings, basis,
            solve=Solve(damping=1., newton=False),
            clamps=(Clamp(sel('ihm.indentation_um'), values),))
        result = basis.synthesize(state[output].mean)
        if not np.isfinite(result).all() or report.residual > 1e-7:
            raise RuntimeError(f'IBM selected runtime solve did not converge: {report.residual}')
        return WindowResult(result, float(report.residual), {
            'backend': 'ibm.runtime.solve_window', 'package_sha256': self.identity['package_sha256'],
            'implementation': self._transfer_spec(kind)[2], 'materialized': model is not None,
            'materialized_components': list(model.components) if model else [],
            'support_origin': 'caller canonical coordinates' if caller_sites else 'synthetic probe',
            'request': model.request.name if model else None,
            'materializer_notes': list(model.notes) if model else [],
            'materializer_missing': list(model.provenance.missing) if model else [],
            'materializer_unimplemented': list(model.provenance.unimplemented) if model else [],
            'full_materialized_graph_executable': False,
            'semantics': 'periodic single-window drift H/(i*omega+1); no carry or causal streaming',
            'input_units': 'um indentation' if model else 'unit drive',
            'output_units': 'donor drift response; not calibrated receptor mV or muscle N',
            'adapter_leak_rate_per_s': 1., 'delay_s': delay_s,
            'delay_origin': 'explicit caller delay; no reconstructed axonal topology',
            'support_gaps': list(SUPPORT_GAPS)})
