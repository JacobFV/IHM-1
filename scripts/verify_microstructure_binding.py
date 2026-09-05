"""Unregistered donor microstructure remains evidence, never inferred body truth."""
from pathlib import Path
import tempfile
import numpy as np
from ihm.human import ImplicitHuman

ROOT = Path(__file__).resolve().parents[1]
body = ImplicitHuman.open(ROOT)
inventory = body.microstructure_evidence()
assert inventory['kidney']['tier'] == 'measured_other_donor'
assert inventory['kidney']['body_registered'] is False
assert inventory['kidney']['usable_for_population_priors'] is False
graph = body.materialize('kidney-arterial-geometry')
raw = dict(np.load(ROOT/'data/derived/microstructure/kidney/example_graph.npz'))
np.testing.assert_allclose(graph['nodes_m'], raw['nodes_mm']*.001, rtol=0, atol=0)
np.testing.assert_array_equal(graph['edges'], raw['edges'])
np.testing.assert_array_equal(graph['points_m'], raw['points_mm']*.001)
assert graph['radius_m'] is None and graph['flow_solution'] is None
assert graph['source_transform_applied'] is False
assert graph['body_registered'] is False
assert len(graph['points_m']) == 19924
graph['nodes_m'][0, 0] += 1
np.testing.assert_array_equal(body.materialize('kidney-arterial-geometry')['nodes_m'],raw['nodes_mm']*.001)
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    for key in ('kidney_card','kidney_graph','kidney_statistics'):
        path = root/body.assets[key]['path']; path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes((ROOT/body.assets[key]['path']).read_bytes())
    opened = ImplicitHuman.open(root)
    assert 'kidney-arterial-geometry' in opened.describe()['materializations']
    (root/body.assets['kidney_graph']['path']).write_bytes(b'changed')
    try: opened.materialize('kidney-arterial-geometry')
    except ValueError as error: assert 'Evidence changed' in str(error)
    else: raise AssertionError('Changed derived graph accepted')
    empty = ImplicitHuman.open(root/'empty')
    assert empty.microstructure_evidence() == {}
    assert 'kidney-arterial-geometry' not in empty.describe()['materializations']
print('PASS: measured donor graph, exact SI conversion, unresolved radius/registration, immutable evidence and source changes')
