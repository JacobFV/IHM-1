"""data integration checks; domain collectors separately verify original assets."""
import json
from pathlib import Path
import tempfile
import numpy as np
from ihm.forge.catalog import EvidenceCatalog
from ihm.forge.acquisition import fetch,sha256

catalog=EvidenceCatalog();summary=catalog.summary()
assert summary['raw_bytes_excluding_git']>5_000_000_000
assert summary['whole_body_calibrated'] is False and summary['whole_body_registered_3d'] is False
found=catalog.search('BloodDensity',source='biogears')
assert any(row['basis']=='literature_validation_target' and float(row['value'])==1050 for row in found)
assert catalog.search('e50_W',source='biogears')
skin=[json.loads(line) for line in Path('data/derived/integumentary/human-wound-observations.jsonl').read_text().splitlines()]
assert len(skin)==160 and sum(r['value'] is not None for r in skin)==155
assert sum(r['excluded_by_authors'] for r in skin)==5
assert len({(r['age_group_table'],r['sex'],r['participant_within_group']) for r in skin})==40
assert all(r['si_scale']==1 and r['age_discrepancy'] for r in skin)
target=catalog.search('wound_field/LA-1/18-29/male',source='human-wound-study-2011')
assert len(target)==1 and target[0]['value']==107
literal_percent=catalog.search('%',source='nhanes-2017-2018')
assert len(literal_percent)==3 and all('blood.hematocrit' in row['name'] for row in literal_percent)
semantics=json.loads(Path('data/derived/semantics/index.json').read_text())
assert semantics['asctb_rows']>1_700_000 and semantics['geometry_measurements']==193
geom=[json.loads(line) for line in Path('data/derived/semantics/vascular-geometry-measurements.jsonl').read_text().splitlines()]
assert sum(row['citation_columns_malformed'] for row in geom)==191
assert all(row['citation_tail_raw'] for row in geom)
with np.load('data/derived/vascular/0050_H_CERE_H_first_cfd_frame.npz',allow_pickle=False) as frame:
    points=frame['Points/Points'];velocity=frame['PointData/velocity_01010']
    assert points.shape==velocity.shape==(153082,3)
    assert np.isfinite(points).all() and np.isfinite(velocity).all()
with tempfile.TemporaryDirectory() as root:
    path=Path(root)/'asset';path.write_bytes(b'known')
    receipt={'url':'https://example.invalid/asset','bytes':5,'sha256':sha256(path)}
    path.with_name('asset.receipt.json').write_text(json.dumps(receipt))
    for kwargs in ({'max_bytes':4},{'prefix':b'wrong'}):
        try:fetch(receipt['url'],path,**kwargs)
        except ValueError:pass
        else:raise AssertionError('cached acquisition bypassed new validation constraints')
print('verified: real-asset catalog, literature target values, semantic parsing flags, decoded original flow, acquisition constraints')
native=catalog.search('Skin1ToSkin2',source='biogears-native')
assert native and all(x['basis']=='native_initialized_scalar_state_or_parameter' for x in native)
assert catalog.search('V0_LH',source='schlosser-selgrade-2000')[0]['value']==1400
fit=catalog.search('intercept',source='human-wound-phenotype-fit')
assert fit[0]['basis']=='human_observation_phenotypic_fit'
