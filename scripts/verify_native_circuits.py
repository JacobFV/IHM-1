"""Native circuit import, SI conversion and snapshot storage accounting."""
from pathlib import Path
import tempfile,sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.native.circuits import native_circuit_graph

FIXTURE='''<BioGearsState><CircuitManager>
<FluidNode><Name>A</Name><Pressure unit="mmHg" value="1"/><Volume unit="mL" value="1"/></FluidNode>
<FluidNode><Name>B</Name><Pressure unit="Pa" value="0"/></FluidNode>
<FluidPath><Name>R</Name><SourceNode>A</SourceNode><TargetNode>B</TargetNode><Resistance unit="mmHg s/mL" value="1"/><Flow unit="mL/s" value="1"/></FluidPath>
<FluidPath><Name>C</Name><SourceNode>B</SourceNode><TargetNode>A</TargetNode><Compliance unit="mL/mmHg" value="2"/><Flow unit="L/min" value="0.06"/></FluidPath>
<FluidCircuit><Name>Main</Name><Node>A</Node><Node>B</Node><Path>R</Path><Path>C</Path></FluidCircuit>
<FluidCircuit><Name>Shared</Name><Node>A</Node><Node>B</Node><Path>R</Path><Path>C</Path></FluidCircuit>
<ThermalNode><Name>A</Name><Temperature unit="degC" value="20"/></ThermalNode>
<ThermalNode><Name>B</Name><Temperature unit="K" value="290"/></ThermalNode>
<ThermalPath><Name>C</Name><SourceNode>A</SourceNode><TargetNode>B</TargetNode><Capacitance unit="kcal/degC" value="1"/></ThermalPath>
</CircuitManager></BioGearsState>'''

def main():
    with tempfile.TemporaryDirectory() as d:
        path=Path(d)/'fixture.xml';path.write_text(FIXTURE);graph=native_circuit_graph(path)
        assert len(graph['nodes'])==4 and len(graph['paths'])==3
        byid={p['id']:p for p in graph['paths']}
        assert len(byid['fluid:R']['circuits'])==2
        np.testing.assert_allclose(byid['fluid:R']['properties']['Flow']['si_value'],1e-6)
        np.testing.assert_allclose(byid['fluid:C']['properties']['Compliance']['si_value'],2e-6/133.322387415)
        np.testing.assert_allclose(byid['thermal:C']['properties']['Capacitance']['si_value'],4184)
        np.testing.assert_allclose(graph['circuits'][0]['balance_audit']['max_abs_net_inflow_si'],0,atol=1e-18)
        assert graph['circuits'][0]['balance_audit']['nodes'][0]['transport_net_inflow_si'] != 0
        path.write_text(FIXTURE.replace('unit="mL/s"','unit="mystery"'))
        try:native_circuit_graph(path)
        except ValueError:pass
        else:raise AssertionError('unknown circuit unit accepted')
        path.write_text(FIXTURE.replace('unit="mmHg" value="1"','unit="K" value="1"'))
        try:native_circuit_graph(path)
        except ValueError:pass
        else:raise AssertionError('pressure accepted temperature dimensions')
        path.write_text(FIXTURE.replace('unit="mmHg" value="1"','unit="mmHg" value="NaN"'))
        try:native_circuit_graph(path)
        except ValueError:pass
        else:raise AssertionError('NaN pressure accepted')
        path.write_text(FIXTURE.replace('<Volume unit="mL" value="1"/>','<Volume unit="mL" value="INF"/>'))
        graph=native_circuit_graph(path);v=graph['nodes'][0]['properties']['Volume'];assert v['si_value'] is None and v['status']=='infinite_boundary'
        json.dumps(graph,allow_nan=False)
    root=Path(__file__).resolve().parents[1]
    real=root/'data/derived/physiology/native_baseline_v2/states/native_stabilized.xml'
    if real.exists():
        graph=native_circuit_graph(real)
        assert len(graph['nodes'])==222 and len(graph['paths'])==395 and len(graph['circuits'])==13
        assert len(graph['systems'])==15
        assert len({x['id'] for x in graph['nodes']})==222
        assert all(p['source'] in {n['id'] for n in graph['nodes']} for p in graph['paths'])
        json.dumps(graph,allow_nan=False)
    print('native circuits: mixed units, shared circuit deduplication, storage flow, unknown units, infinite boundary and actual state counts PASS')

if __name__=='__main__':main()
