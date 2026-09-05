"""inspect, condition, materialize and export views of the body substrate."""
import argparse
import json
from pathlib import Path
from itertools import groupby
import numpy as np
from ihm import body, materialize, Request
from ihm.materialize import LIBRARY
from ihm.materialize.model import Model
from ihm.materialize.cardiopulmonary import cardiopulmonary_model, CardiopulmonaryModel
from ihm.processes.cardiopulmonary import Parameters
from ihm.materialize.skin import SkinPatch, Cell, SurfaceSite, skin_model, electric_field
from ihm.materialize.fluids import FluidNetwork, FluidCompartment, FluidEdge, fluid_model, flows
from ihm.forge.ingest import SourceCard, load


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)


def geometry(path, kind):
    if path is None:
        return SkinPatch.line() if kind == 'skin' else FluidNetwork.example()
    data = json.loads(Path(path).read_text())
    if kind == 'skin':
        data['cells'] = tuple(Cell(**c) for c in data['cells'])
        data['surface'] = tuple(SurfaceSite(**s) for s in data['surface'])
        data['gap_edges'] = tuple(tuple(e) for e in data['gap_edges'])
        data['surface_edges'] = tuple(tuple(e) for e in data['surface_edges'])
        return SkinPatch(**data)
    return FluidNetwork(tuple(FluidCompartment(**c) for c in data['compartments']),
                        tuple(FluidEdge(**e) for e in data['edges']))


def read_model(path):
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data['metadata']))
    return CardiopulmonaryModel.load(path) if metadata.get('model_type') == 'cardiopulmonary' else Model.load(path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('inventory'); sub.add_parser('models')
    data=sub.add_parser('data',help='inspect acquired real source data and coefficients')
    data.add_argument('action',choices=('summary','search'))
    data.add_argument('--term');data.add_argument('--source');data.add_argument('--limit',type=int,default=20)
    data.add_argument('--catalog',default='data/derived/evidence-catalog.sqlite')
    integrated=sub.add_parser('integrated',help='inspect or freeze the acquired evidence-backed human substrate')
    integrated.add_argument('--root',default='.')
    integrated.add_argument('--output')
    integrated.add_argument('--fields',action='store_true')
    native=sub.add_parser('native',help='execute original native multisystem physiology')
    native.add_argument('--config',required=True);native.add_argument('--output',required=True)
    serve=sub.add_parser('serve',help='open the local 3D scientific workbench')
    serve.add_argument('--port',type=int,default=8765);serve.add_argument('--root',default='.')
    build = sub.add_parser('build')
    build.add_argument('model', choices=[*LIBRARY, 'whole-body', 'skin', 'fluids', 'cardiopulmonary'])
    build.add_argument('--subject', required=True); build.add_argument('--output', required=True)
    build.add_argument('--population-prior', help='measured joint population-state prior JSON')
    build.add_argument('--parameters', help='cardiopulmonary Parameters JSON; supine only')
    build.add_argument('--geometry', help='explicit SkinPatch or FluidNetwork JSON')
    build.add_argument('--source', nargs=2, action='append', default=[], metavar=('CARD', 'DATA'))
    build.add_argument('--time', type=float, default=None)
    forecast = sub.add_parser('forecast'); forecast.add_argument('artifact')
    forecast.add_argument('--times', nargs='+', type=float, required=True)
    forecast.add_argument('--clamp', action='append', default=[], metavar='COMPONENT=VALUE')
    inspect = sub.add_parser('inspect'); inspect.add_argument('artifact')
    demo = sub.add_parser('demo'); demo.add_argument('--output', default='artifacts')
    args = parser.parse_args(argv)
    try:
        if args.command == 'integrated':
            from ihm.human import ImplicitHuman
            human=ImplicitHuman.open(args.root)
            if args.output:human.save(args.output)
            print(json.dumps(human.fields() if args.fields else human.describe(),indent=2,allow_nan=False))
        elif args.command == 'native':
            from ihm.native import NativeConfig,run_native
            print(json.dumps(run_native(NativeConfig.from_dict(json.loads(Path(args.config).read_text())),args.output),indent=2,allow_nan=False))
        elif args.command == 'serve':
            from ihm.app import serve
            serve(args.root,args.port)
        elif args.command == 'data':
            from ihm.forge.catalog import EvidenceCatalog
            catalog=EvidenceCatalog(args.catalog)
            if args.action=='search' and not args.term: raise ValueError('--term is required for search')
            result=catalog.summary() if args.action=='summary' else catalog.search(args.term,args.limit,args.source)
            print(json.dumps(result,indent=2))
        elif args.command == 'inventory':
            r = body()
            print(json.dumps({'counts': {k: len(getattr(r, k)) for k in ('components', 'supports', 'anatomy', 'topologies', 'processes')},
                              'components': r.describe()}, indent=2))
        elif args.command == 'models':
            print(json.dumps({**LIBRARY, 'cardiopulmonary': 'beating four-chamber heart and tidal breathing, supine', 'whole-body': 'all registered organ-scale components', 'skin': 'explicit heterogeneous electrical patch',
                              'fluids': 'explicit conservative vascular/interstitial/lymphatic network'}, indent=2))
        elif args.command == 'build':
            if args.geometry and args.model not in ('skin', 'fluids'):
                raise ValueError('--geometry applies only to skin or fluids')
            if args.parameters and args.model != 'cardiopulmonary':
                raise ValueError('--parameters applies only to cardiopulmonary')
            if args.model == 'cardiopulmonary':
                if args.source:
                    raise ValueError('the nonlinear cardiopulmonary view does not yet support evidence assimilation')
                parameters = Parameters(**json.loads(Path(args.parameters).read_text())) if args.parameters else Parameters()
                m = cardiopulmonary_model(args.subject, parameters); r = body(cardiopulmonary=parameters)
            elif args.model == 'skin':
                patch = geometry(args.geometry, 'skin'); r = body(skin=patch); m = skin_model(patch, args.subject)
            elif args.model == 'fluids':
                network = geometry(args.geometry, 'fluids'); r = body(fluids=network); m = fluid_model(network, args.subject)
            else:
                r = body(); targets = tuple(r.components) if args.model == 'whole-body' else LIBRARY[args.model]
                m = materialize(r, Request(targets, args.subject))
            if args.population_prior:
                if isinstance(m,CardiopulmonaryModel):
                    raise ValueError('population Gaussian initialization is not implemented for the nonlinear view')
                from ihm.forge.population import initialize_population
                initialize_population(m,args.population_prior)
            evidence = []
            for card_path, data_path in args.source:
                evidence.extend(load(data_path, SourceCard.read(card_path), r))
            # Validate complete input before partial conditioning of the artifact.
            if len({e.key for e in evidence}) != len(evidence):
                raise ValueError('duplicate evidence across source inputs')
            if any(e.subject != args.subject for e in evidence):
                raise ValueError('data contains other subjects; subset explicitly before building')
            for _, batch in groupby(sorted(evidence, key=lambda e: e.time), key=lambda e: e.time):
                m.assimilate(list(batch))
            if args.time is not None:
                m.advance(args.time)
            m.save(args.output)
            print(json.dumps({'artifact': args.output, 'states': len(m.components), 'time': m.time,
                              'evidence_records': len(getattr(m, 'seen', ())), 'validated_biology': False}))
        elif args.command == 'forecast':
            clamps = {}
            for item in args.clamp:
                key, value = item.split('=', 1)
                if key in clamps:
                    raise ValueError('duplicate clamp')
                clamps[key] = float(value)
            m = read_model(args.artifact)
            if isinstance(m, CardiopulmonaryModel):
                if clamps:
                    raise ValueError('nonlinear cardiac/respiratory clamps are not implemented; use declared parameters')
                result = m.forecast(args.times)
            else:
                result = m.forecast(args.times, clamps)
            print(json.dumps(result, default=jsonable, indent=2, allow_nan=False))
        elif args.command == 'inspect':
            m = read_model(args.artifact)
            print(json.dumps({'subject': m.subject, 'time': m.time, 'components': m.components,
                              'targets': getattr(m, 'targets', m.components), 'provenance': m.provenance}, indent=2))
        else:
            out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
            summary = {}
            r = body()
            for name, targets in LIBRARY.items():
                m = materialize(r, Request(targets, 'synthetic-demo')); m.save(out/f'{name}.npz')
                summary[name] = {'states': len(m.components), 'forecast': m.forecast([0, 1, 60])}
            for wounded in (False, True):
                patch = SkinPatch.line(wound=wounded); m = skin_model(patch); m.advance(5)
                name = 'skin-wound' if wounded else 'skin-intact'
                m.save(out/f'{name}.npz'); summary[name] = electric_field(m, patch)
            network = FluidNetwork.example(); m = fluid_model(network); m.advance(10)
            m.save(out/'fluids.npz'); summary['fluids'] = flows(m, network)
            cardio = cardiopulmonary_model(); cardio.save(out/'cardiopulmonary.npz')
            cycle = cardio.forecast(np.linspace(0, 30, 1501))
            (out/'cardiopulmonary-trace.json').write_text(json.dumps(cycle, default=jsonable, allow_nan=False))
            summary['cardiopulmonary'] = {'trace': 'cardiopulmonary-trace.json', 'posture': 'supine', 'validated_biology': False}
            (out/'demo.json').write_text(json.dumps(summary, default=jsonable, indent=2, allow_nan=False))
            print(f'Wrote {len(summary)} illustrative predictor artifacts and demo.json to {out}')
    except (ValueError, KeyError, OSError, TypeError) as exc:
        parser.exit(2, f'ihm: {exc}\n')


if __name__ == '__main__':
    main()
