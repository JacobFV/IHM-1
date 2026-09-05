"""Materialize a typed graph and compact index from an actual native saved state."""
from pathlib import Path
from collections import defaultdict
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from ihm.native.circuits import native_circuit_graph


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser();parser.add_argument('--state',type=Path,default=root/'data/derived/physiology/native_baseline_v2/states/native_stabilized.xml');parser.add_argument('--output',type=Path,default=root/'data/derived/native-circuits')
    args=parser.parse_args();graph=native_circuit_graph(args.state);args.output.mkdir(parents=True,exist_ok=True)
    graphpath=args.output/'graph.json';graphpath.write_text(json.dumps(graph,separators=(',',':'),allow_nan=False))
    compnames=defaultdict(list)
    for c in graph['compartments']:compnames[c['kind']].append(c['name'])
    circuits=[dict(id=c['id'],family=c['family'],name=c['name'],node_count=len(c['nodes']),path_count=len(c['paths']),
        reference_nodes=c['reference_nodes'],balance_audit={k:v for k,v in c['balance_audit'].items() if k!='nodes'}) for c in graph['circuits']]
    index=dict(schema_version=1,source=graph['source'],simulation_time=graph['simulation_time'],parameter_status=graph['parameter_status'],
        graph_path=str(graphpath.relative_to(root)) if graphpath.is_relative_to(root) else str(graphpath),
        graph_bytes=graphpath.stat().st_size,summary=graph['summary'],circuits=circuits,
        systems=[dict(type=s['type'],scalar_fields=len(s['properties'])) for s in graph['systems']],
        compartment_names_by_kind=dict(compnames),limitations=graph['limitations'])
    (args.output/'index.json').write_text(json.dumps(index,separators=(',',':'),allow_nan=False))
    print(json.dumps(dict(graph= str(graphpath),bytes=graphpath.stat().st_size,summary=graph['summary']),indent=2))

if __name__=='__main__':main()
