"""Inventory deduplication, unit-equivalent reservoirs and source transfer balance."""
from pathlib import Path
import sys,tempfile,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ihm.native.budgets import fluid_inventory,audit_fluid_budget


def fixture(blood_ml,water_ml,extra_node=False):
    extra='<FluidNode><Name>Other</Name><Volume value="8" unit="mL"/></FluidNode>' if extra_node else ''
    selected='<Node>Other</Node>' if extra_node else ''
    return f'''<BioGearsState xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <CircuitManager>
    <FluidNode><Name>Blood</Name><Volume value="{blood_ml}" unit="mL"/></FluidNode>
    <FluidNode><Name>Tissue</Name><Volume value="0.03" unit="L"/></FluidNode>
    <FluidNode><Name>Ground</Name><Volume value="INF" unit="mL"/></FluidNode>
    <FluidNode><Name>Unselected</Name><Volume value="999" unit="L"/></FluidNode>{extra}
    <FluidCircuit><Name>FullCardiovascular</Name><Node>Blood</Node><Node>Blood</Node><Node>Tissue</Node><Node>Ground</Node>{selected}</FluidCircuit>
    <FluidCircuit><Name>DuplicateView</Name><Node>Blood</Node><Node>Tissue</Node></FluidCircuit>
    </CircuitManager>
    <System xsi:type="BioGearsGastrointestinalSystemData"><StomachContents><Water value="{water_ml}" unit="mL"/></StomachContents></System>
    </BioGearsState>'''


def main():
    with tempfile.TemporaryDirectory() as directory:
        run=Path(directory);states=run/'states';states.mkdir()
        a=states/'native_stabilized.xml';b=states/'native_final.xml'
        a.write_text(fixture(100,500));b.write_text(fixture(300,300))
        inv=fluid_inventory(a)
        assert len(inv['circuit_volumes_m3'])==2 and inv['infinite_boundaries']==['fluid:Ground']
        np.testing.assert_allclose(inv['combined_finite_volume_m3'],630e-6)
        budget=audit_fluid_budget(run)
        np.testing.assert_allclose(budget['circuit_volume_change_m3'],200e-6)
        np.testing.assert_allclose(budget['stomach_water_change_m3'],-200e-6)
        assert abs(budget['combined_volume_change_m3'])<1e-18
        b.write_text(fixture(301,300))
        np.testing.assert_allclose(audit_fluid_budget(run)['combined_volume_change_m3'],1e-6)
        b.write_text(fixture(300,300,True))
        try:audit_fluid_budget(run)
        except ValueError:pass
        else:raise AssertionError('changed circuit inventory accepted')
        b.write_text(fixture(300,300).replace('unit="L"','unit="Pa"'))
        try:fluid_inventory(b)
        except ValueError:pass
        else:raise AssertionError('wrong volume dimension accepted')
    root=Path(__file__).resolve().parents[1]
    actual=root/'data/derived/physiology/native_hour_rest'
    if actual.exists():
        budget=audit_fluid_budget(actual)
        np.testing.assert_allclose(budget['initial']['stomach_water_m3'],.0005)
        assert budget['stomach_water_change_m3']<0 and budget['circuit_volume_change_m3']>0
        json.dumps(budget,allow_nan=False)
        print('Actual hour combined finite-volume change (m3):',budget['combined_volume_change_m3'])
    print('Native budgets: duplicate views/nodes, excluded external nodes, SI reservoirs, internal transfer, nonzero residual and topology/unit rejection PASS')

if __name__=='__main__':main()
