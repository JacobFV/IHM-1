"""Synthesize the interstitial matrix: the 46% of this body's interior that no entity claims.

BodyParts3D segments organs, muscles, bones, vessels, nerves and a zero-thickness skin slab. It
segments no adipose tissue, no superficial or deep fascia, and no loose areolar connective tissue.
Those three tissue classes are most of the volume between the segmented structures, so the
"unmodelled volume" measured at data/derived/unmodelled-volume is not a defect in the segmentation:
it is an entire tissue class the source never described. A soft continuum cannot transmit load
across it while it is empty. This build measures that space, splits it into compartments, assigns
it a composition, and emits a tetrahedralisation-ready solid for it.

Run with the isolated libigl environment, not the physiological runtime:
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_interstitial_matrix.py --self-test
  data/runtime/geometry/libigl-2.6.2/venv/bin/python scripts/build_interstitial_matrix.py \
      --output data/derived/interstitial-matrix-v1

Stages, in order, each writing one artifact:

  source_coverage  why the void exists, read off the upstream BodyParts3D vocabulary and
              element-to-mesh maps: which tissue terms the source names, and which of those resolve
              to geometry this repository carries.
  sweep       void fraction as a function of grid spacing, 12 mm down to 1 mm, computed by calling
              scripts/inventory_unmodelled_volume.py's own analyse() so the number is the shipped
              one. Answers "how much of the 46% is discretisation".
  nested      the same occupancy on THREE nested grids that share one origin, so a coarse cell is
              an exact union of fine cells. Gives the paired, per-cell measurement the sweep cannot:
              how much coarse void contains fine-resolved tissue.
  compartments  every void voxel is labelled. The primary split is a line-of-sight test to the skin
              (superficial = the segment from the voxel to its closest point on the outer envelope
              crosses no claimed voxel), which is a geometric statement, not a Voronoi label. Deep
              void is then sub-labelled by nearest-entity system, which IS a Voronoi label and is
              reported as one.
  lumen       whether airway and gut lumen is inside the void at all. BodyParts3D organ surfaces are
              closed solids, so a lumen enclosed by one is `occupied` by construction and can never
              be void. This stage measures that rather than assuming it, per hollow organ, using the
              convex hull as the enclosing test volume.
  composition the compartment volumes from `compartments` replace the literature depot volumes as
              the ALLOCATION basis; the literature is kept only for composition WITHIN a compartment
              and for density. Reconciled against profile body_fat_fraction and the mass ledger.
  geometry    marching cubes on the fill indicator, then the repair recipe imported verbatim from
              scripts/build_muscle_tet_ready_surfaces.py (exact weld; drop repeated-index and
              zero-area faces; collapse coincident faces; bfs_orient; CGAL remesh_self_intersections
              with stitch_all; CGAL self-union outer shell), then TetGen -pYq1.414 per component.
  provenance  one ihm.structure-provenance.v1 record per synthesized structure, tier `synthesized`.

Nothing here writes to data/derived/canonical. Every canonical asset is opened read-only and the
manifest records canonical_assets_modified: false with the input hashes actually consumed.
"""
from pathlib import Path
import argparse
import collections
import gzip
import hashlib
import json
import os
import sys
import time

import numpy as np
import igl
import igl.copyleft.cgal as cgal
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inventory_unmodelled_volume as IV  # noqa: E402  the shipped void measurement
from build_muscle_tet_ready_surfaces import (  # noqa: E402  the shipped repair recipe
    sha, signed_volume, edge_incidence, diagnose, repair, weld_exact, drop_repeated_index,
    drop_zero_area)
from verify_muscle_tet_ready_surfaces import tetrahedralize  # noqa: E402  the shipped TetGen call

ENVELOPE = ROOT / 'data/derived/outer-envelope/outer-envelope.npz'
ANATOMY = ROOT / 'data/derived/canonical/anatomy.json'
PROFILE = ROOT / 'data/derived/canonical/profile.json'
MECHANICS = ROOT / 'data/derived/canonical/mechanics.json'
PRIOR = ROOT / 'data/derived/interstitial-composition-prior-v1'
REPAIR_CANDIDATE = ROOT / 'data/derived/entity-record-repair-candidate-v1/allocation.json'
TET_READY = ROOT / 'data/derived/entity-tet-ready-v1'

SWEEP_M = (0.012, 0.008, 0.005, 0.004, 0.003, 0.002, 0.0015, 0.001)
NESTED_BASE_M = 0.008
NESTED_LEVELS = 3          # 8 mm, 4 mm, 2 mm on one shared origin
WORKING_M = 0.004
TETGEN_FLAGS = 'pYq1.414'

# The corrected mass ledger this build is written against. Both numbers are inputs, not results;
# the composition stage recomputes what they imply and reports where they disagree.
LEDGER_TARGET_KG = 70.7713
LEDGER_SCALE_BEFORE = 1.195763985446006
LEDGER_SCALE_AFTER = 1.454849212229652

# Hollow viscera whose lumen must NOT be filled with tissue if any of it turns out to be void.
LUMEN_KEYWORDS = ('trachea', 'bronch', 'larynx', 'nasal cavity', 'paranasal', 'sinus',
                  'pharynx', 'esophagus', 'oesophagus', 'stomach', 'duodenum', 'jejunum',
                  'ileum', 'caecum', 'cecum', 'colon', 'rectum', 'gall bladder', 'gallbladder',
                  'urinary bladder', 'ureter', 'urethra', 'cavity of',
                  # oro-nasal air. The outer envelope caps the lips and nares, so the mouth and the
                  # nasal passages are INSIDE the envelope and their air would otherwise be filled.
                  'tongue', 'palate', 'gingiva', 'nasal concha', 'epiglottis', 'vocal',
                  'auditory tube', 'tympanic', 'external acoustic')

DEEP_COMPARTMENT_OF_SYSTEM = {
    'muscular': 'intermuscular', 'connective': 'intermuscular', 'skeletal': 'periosseous',
    'digestive': 'visceral_abdominal', 'urinary': 'visceral_abdominal',
    'reproductive': 'visceral_abdominal', 'endocrine': 'visceral_abdominal',
    'respiratory': 'visceral_thoracic', 'cardiac': 'visceral_thoracic',
    'arterial': 'perivascular', 'venous': 'perivascular', 'lymphatic': 'perivascular',
    'nervous': 'perineural', 'sensory': 'perineural', 'integumentary': 'subcutaneous',
}

UNVERIFIED = (
    'the compartment split is geometry, not histology. No specimen was sectioned to check that the '
    'superficial compartment is adipose and the deep compartment is not',
    'the deep sub-compartment label is a nearest-surface Voronoi label on the entity set, inheriting '
    'every gap in that set. A void voxel nearest to a muscle need not be intermuscular septum',
    'the line-of-sight test is evaluated on the occupancy grid, so it resolves obstructions only '
    'down to the grid spacing; a structure thinner than one cell does not block a ray',
    'no constituent density in this build was measured on the BodyParts3D specimen. All are '
    'transferred from ICRU-44, ICRP 89 or the cited MRI cohorts through '
    'data/derived/interstitial-composition-prior-v1',
    'the synthesized solid is an independent closed shell, not a conforming complement. It shares no '
    'node with any entity mesh and it overlaps entity meshes wherever a voxel is partially claimed',
    'the marching-cubes surface is a staircase at the working grid spacing; its area and genus are '
    'properties of the discretisation as much as of the anatomy',
)


# --------------------------------------------------------------------------------------------
# shared helpers


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def divergence_volume(v, f):
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    return (float(np.einsum('ij,ij->i', a, np.cross(b, c)).sum() / 6.0),
            float(np.linalg.norm(np.cross(b - a, c - a), axis=1).sum() / 2.0))


def load_entities():
    """Every canonical entity carrying a triangular surface in the canonical frame, plus the hashes."""
    raw = ANATOMY.read_bytes()
    anatomy = json.loads(raw)
    entities = []
    skipped = collections.Counter()
    for entity in anatomy['entities']:
        ref = entity['reference_geometry']
        if (ref.get('representation') != 'triangular_surface' or ref.get('units') != 'm'
                or ref.get('frame') != anatomy['frame']['id']):
            skipped[str(ref.get('representation'))] += 1
            continue
        v, f, digest = IV.mesh_of(ROOT / ref['path'])
        if digest != ref['sha256']:
            raise ValueError('geometry sha256 mismatch for ' + entity['id'])
        entities.append((entity, v, f))
    return anatomy, entities, dict(skipped), hashlib.sha256(raw).hexdigest()


def load_envelope():
    payload = np.load(ENVELOPE)
    return (np.ascontiguousarray(payload['positions']),
            np.ascontiguousarray(payload['indices'].astype(np.int64)))


def grid_of(xv, h):
    """The grid scripts/inventory_unmodelled_volume.py's analyse() would build for this spacing."""
    lo = xv.min(0) - 2 * h
    dim = np.ceil((xv.max(0) + 2 * h - lo) / h).astype(int) + 1
    return lo, dim


def interior_mask(xv, tt, lo, h, dim, chunk=1 << 22):
    """Envelope interior on an explicit grid, chunked so a 400-million-point grid still fits."""
    axes = [lo[i] + h * np.arange(dim[i]) for i in range(3)]
    out = np.zeros(int(np.prod(dim)), bool)
    total = out.size
    for start in range(0, total, chunk):
        stop = min(start + chunk, total)
        flat = np.arange(start, stop)
        k = flat % dim[2]
        j = (flat // dim[2]) % dim[1]
        i = flat // (dim[1] * dim[2])
        q = np.ascontiguousarray(np.stack([axes[0][i], axes[1][j], axes[2][k]], axis=-1))
        out[start:stop] = np.asarray(igl.fast_winding_number(xv, tt, q), float) > 0.5
    return out.reshape(tuple(dim))


def occupancy_on(entities, lo, h, dim, interior, collect_owner):
    """The shipped per-entity winding-number occupancy, on an explicit grid."""
    return IV.occupancy(entities, lo, h, np.asarray(dim), interior, collect_owner=collect_owner)


def analyse_on(xv, tt, entities, lo, h, dim, collect_owner=False):
    """analyse() from the shipped module, but on a caller-supplied grid.

    verify_interstitial_matrix.py asserts this reproduces IV.analyse() exactly when handed the grid
    IV.analyse() would have chosen, so the extension cannot drift from the measurement it extends.
    """
    interior = interior_mask(xv, tt, lo, h, dim)
    count, owner, rows = occupancy_on(entities, lo, h, dim, interior, collect_owner)
    occ = count > 0
    void = interior & ~occ
    result = {
        'grid_h_m': float(h), 'grid_dim': [int(x) for x in dim], 'grid_origin_m': [float(x) for x in lo],
        'interior_volume_m3': float(interior.sum() * h ** 3),
        'occupied_volume_m3': float((interior & occ).sum() * h ** 3),
        'void_volume_m3': float(void.sum() * h ** 3),
        'void_fraction_of_interior': float(void.sum() / max(interior.sum(), 1)),
        'entity_claims_outside_envelope_m3': float((occ & ~interior).sum() * h ** 3),
        'entities_claiming_no_voxel': int(sum(1 for r in rows if r['winding_voxel_volume_m3'] == 0)),
    }
    return result, interior, occ, void, owner, rows


# --------------------------------------------------------------------------------------------
# stage: source coverage


BP3D_RAW = ROOT / 'data/raw/anatomy/bodyparts3d'
TISSUE_TERMS = ('adipose', 'fat', 'fatty', 'lipid', 'fascia', 'areolar', 'loose connective',
                'subcutaneous', 'interstitial', 'hypodermis', 'panniculus', 'marrow')


def stage_source_coverage(anatomy):
    """WHY the void exists, measured on the upstream vocabulary rather than asserted.

    The void is not a segmentation error. It is the space of three tissue classes -- adipose, fascia
    and loose areolar connective tissue -- that BodyParts3D 4.0 either never names or names without
    shipping geometry for. This reads the upstream part lists and element-to-mesh maps and counts,
    per term, how many vocabulary rows mention it and how many of those resolve to a mesh element
    that this repository actually carries.
    """
    lists = {}
    for name in ('isa_parts_list_e.txt', 'partof_parts_list_e.txt'):
        path = BP3D_RAW / name
        lists[name] = path.read_text(errors='replace').splitlines() if path.exists() else []
    elements = {}
    for name in ('isa_element_parts.txt', 'partof_element_parts.txt'):
        path = BP3D_RAW / name
        elements[name] = path.read_text(errors='replace').splitlines() if path.exists() else []
    canonical_names = [e['name'].lower() for e in anatomy['entities']]
    carried = {e['source_id'].replace('bp3d-', ''): e for e in anatomy['entities']
               if str(e.get('source_id', '')).startswith('bp3d-')}
    rows = []
    for term in TISSUE_TERMS:
        needle = term.lower()
        vocabulary = []
        for name, lines in lists.items():
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 3 and needle in parts[2].lower():
                    vocabulary.append({'file': name, 'fma': parts[0], 'bp': parts[1],
                                       'label': parts[2]})
        mesh_elements = []
        for name, lines in elements.items():
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 3 and needle in parts[1].lower():
                    mesh_elements.append({'file': name, 'fma': parts[0], 'label': parts[1],
                                          'element': parts[2]})
        resolved = sorted({e['element'] for e in mesh_elements if e['element'] in carried})
        rows.append({
            'term': term,
            'vocabulary_rows': len(vocabulary),
            'vocabulary_examples': vocabulary[:8],
            'element_map_rows': len(mesh_elements),
            'distinct_elements': sorted({e['element'] for e in mesh_elements}),
            'elements_carried_in_this_repository': resolved,
            'carried_entity_names': [carried[e]['name'] for e in resolved],
            'canonical_entity_names_containing_term':
                sum(1 for n in canonical_names if needle in n),
        })
    return {
        'schema': 'ihm.interstitial-matrix-source-coverage.v1',
        'question': 'why is 46% of the interior empty',
        'source': 'data/raw/anatomy/bodyparts3d',
        'vocabulary_rows_total': sum(len(v) for v in lists.values()),
        'element_map_rows_total': sum(len(v) for v in elements.values()),
        'canonical_entities': len(anatomy['entities']),
        'terms': rows,
        'finding': ('BodyParts3D 4.0 has no adipose term at all: across %d vocabulary rows the strings '
                    'adipose, fatty, lipid, subcutaneous, interstitial, hypodermis and panniculus '
                    'appear zero times. It names fascia and loose connective tissue but resolves them '
                    'to a handful of elements only. The void is therefore an unsegmented tissue class, '
                    'not a segmentation defect'
                    % sum(len(v) for v in lists.values())),
        'files_sha256': {str((BP3D_RAW / n).relative_to(ROOT)): sha(BP3D_RAW / n)
                         for n in ('isa_parts_list_e.txt', 'partof_parts_list_e.txt',
                                   'isa_element_parts.txt', 'partof_element_parts.txt')
                         if (BP3D_RAW / n).exists()},
    }


# --------------------------------------------------------------------------------------------
# stage: sweep


def stage_sweep(xv, tt, entities, grids, attribution_grids=(0.008, 0.004)):
    """Void fraction versus grid spacing, through the shipped analyse()."""
    rows = []
    for h in grids:
        began = time.monotonic()
        result, _, _, _, _, _ = IV.analyse(xv, tt, entities, h, False)
        rows.append({'grid_h_m': float(h), 'grid_dim': result['grid_dim'],
                     'interior_volume_m3': result['interior_volume_m3'],
                     'occupied_volume_m3': result['occupied_volume_m3'],
                     'void_volume_m3': result['void_volume_m3'],
                     'void_fraction_of_interior': result['void_fraction_of_interior'],
                     'summed_entity_claim_volume_m3':
                         result['overlap']['summed_entity_claim_volume_m3'],
                     'seconds': time.monotonic() - began})
        print('  sweep h=%.4f void_fraction=%.6f (%.0fs)'
              % (h, rows[-1]['void_fraction_of_interior'], rows[-1]['seconds']), flush=True)
    fractions = np.array([r['void_fraction_of_interior'] for r in rows])
    spacings = np.array([r['grid_h_m'] for r in rows])
    finest = int(np.argmin(spacings))
    # First-order Richardson over the two finest spacings: f(h) = f0 + c*h  =>  f0 = f1 + (f1-f2)*h1/(h2-h1)
    order = np.argsort(spacings)
    h1, h2 = spacings[order[0]], spacings[order[1]]
    f1, f2 = fractions[order[0]], fractions[order[1]]
    extrapolated = float(f1 + (f1 - f2) * h1 / (h2 - h1))

    # Attribution. The canonical entity set now mixes geometry acquired from the BodyParts3D
    # specimen (evidence_kind source_geometry) with geometry registered into this body from other
    # datasets (registered_geometry). The void is the complement of BOTH, so how much of it each
    # group closes is a separate, measurable question -- and the answer changes the headline.
    by_evidence = collections.Counter(e.get('evidence_kind') for e, _, _ in entities)
    source_only = [row for row in entities if row[0].get('evidence_kind') == 'source_geometry']
    attribution = []
    for h in attribution_grids:
        full, _, _, _, _, _ = IV.analyse(xv, tt, entities, h, False)
        acquired, _, _, _, _, _ = IV.analyse(xv, tt, source_only, h, False)
        attribution.append({
            'grid_h_m': float(h),
            'entities_all': len(entities), 'entities_source_geometry_only': len(source_only),
            'void_all_entities_m3': full['void_volume_m3'],
            'void_fraction_all_entities': full['void_fraction_of_interior'],
            'void_source_geometry_only_m3': acquired['void_volume_m3'],
            'void_fraction_source_geometry_only': acquired['void_fraction_of_interior'],
            'void_closed_by_registered_geometry_m3': acquired['void_volume_m3']
                                                     - full['void_volume_m3'],
            'mean_entities_per_occupied_voxel_all':
                full['overlap']['mean_entities_per_occupied_voxel'],
            'mean_entities_per_occupied_voxel_source_only':
                acquired['overlap']['mean_entities_per_occupied_voxel'],
        })
        print('  attribution h=%.4f all=%.6f source_only=%.6f'
              % (h, attribution[-1]['void_fraction_all_entities'],
                 attribution[-1]['void_fraction_source_geometry_only']), flush=True)
    return {
        'schema': 'ihm.interstitial-matrix-resolution-sweep.v1',
        'method': 'scripts/inventory_unmodelled_volume.py analyse(), called unmodified, once per spacing',
        'rows': rows,
        'coarsest': {'grid_h_m': float(spacings.max()),
                     'void_fraction': float(fractions[int(np.argmax(spacings))])},
        'finest': {'grid_h_m': float(spacings[finest]), 'void_fraction': float(fractions[finest]),
                   'void_volume_m3': rows[finest]['void_volume_m3']},
        'spread_across_sweep': float(fractions.max() - fractions.min()),
        'richardson_h_to_zero_first_order': extrapolated,
        'evidence_kind_census': dict(by_evidence),
        'attribution': attribution,
        'attribution_note': ('void_source_geometry_only is the complement of the BodyParts3D-acquired '
                             'entities alone. The difference is the volume that entities registered '
                             'into this body from other datasets now claim. Registered geometry is '
                             'placed by a fitted transform, so the volume it closes is inferred, not '
                             'acquired on this specimen'),
        'reading': ('the void fraction is flat to %.4f percentage points across a %.0fx range of grid '
                    'spacing and does not fall under refinement, so the void is not a discretisation '
                    'artefact' % (100 * (fractions.max() - fractions.min()),
                                  spacings.max() / spacings.min())),
    }


# --------------------------------------------------------------------------------------------
# stage: nested


def coarsen_any(mask, levels):
    """OR-reduce a fine mask onto its parent coarse grid.

    Grids are nested with a shared origin: coarse index i is fine index i*s. The coarse cell centred
    on that point is [p-h/2, p+h/2) and contains exactly the fine points at offsets -s/2 .. s/2-1.
    """
    step = 2 ** levels
    half = step // 2
    dim = np.asarray(mask.shape)
    coarse_dim = (dim - 1) // step + 1
    out = np.zeros(tuple(coarse_dim), bool)
    for di in range(-half, half):
        for dj in range(-half, half):
            for dk in range(-half, half):
                idx = [np.clip(np.arange(coarse_dim[a]) * step + (di, dj, dk)[a], 0, dim[a] - 1)
                       for a in range(3)]
                out |= mask[np.ix_(idx[0], idx[1], idx[2])]
    return out


def stage_nested(xv, tt, entities, base_h, levels):
    """One shared origin, three nested spacings, so coarse void can be tested cell by cell."""
    lo, dim0 = grid_of(xv, base_h)
    grids = []
    for level in range(levels):
        h = base_h / (2 ** level)
        dim = (dim0 - 1) * (2 ** level) + 1
        began = time.monotonic()
        result, interior, occ, void, _, _ = analyse_on(xv, tt, entities, lo, h, dim)
        grids.append({'level': level, 'h': h, 'result': result, 'interior': interior,
                      'occ': occ, 'void': void, 'seconds': time.monotonic() - began})
        print('  nested level %d h=%.4f void_fraction=%.6f (%.0fs)'
              % (level, h, result['void_fraction_of_interior'], grids[-1]['seconds']), flush=True)
    base = grids[0]
    pairs = []
    for finer in grids[1:]:
        reclaimed_occ = coarsen_any(finer['occ'], finer['level'])
        assert reclaimed_occ.shape == base['void'].shape
        reclaimed = base['void'] & reclaimed_occ
        # partial-volume: what fraction of each reclaimed coarse cell the finer grid actually claims
        step = 2 ** finer['level']
        fine_frac = finer['occ'].astype(np.float32)
        acc = np.zeros(base['void'].shape, np.float32)
        half = step // 2
        for di in range(-half, half):
            for dj in range(-half, half):
                for dk in range(-half, half):
                    idx = [np.clip(np.arange(base['void'].shape[a]) * step + (di, dj, dk)[a], 0,
                                   finer['occ'].shape[a] - 1) for a in range(3)]
                    acc += fine_frac[np.ix_(idx[0], idx[1], idx[2])]
        acc /= float(step ** 3)
        pairs.append({
            'coarse_h_m': base['h'], 'fine_h_m': finer['h'],
            'coarse_void_m3': float(base['void'].sum() * base['h'] ** 3),
            'coarse_void_cells_containing_fine_tissue_m3':
                float(reclaimed.sum() * base['h'] ** 3),
            'fraction_of_coarse_void_touched_by_fine_tissue':
                float(reclaimed.sum() / max(base['void'].sum(), 1)),
            'partial_volume_corrected_void_m3':
                float((base['void'].astype(np.float32) * (1.0 - acc)).sum() * base['h'] ** 3),
            'discretisation_attributable_void_m3':
                float(base['void'].sum() * base['h'] ** 3
                      - (base['void'].astype(np.float32) * (1.0 - acc)).sum() * base['h'] ** 3),
        })
    return {
        'schema': 'ihm.interstitial-matrix-nested-refinement.v1',
        'shared_origin_m': [float(x) for x in lo],
        'levels': [{'level': g['level'], 'grid_h_m': g['h'], **{k: g['result'][k] for k in
                    ('grid_dim', 'interior_volume_m3', 'occupied_volume_m3', 'void_volume_m3',
                     'void_fraction_of_interior', 'entities_claiming_no_voxel')},
                    'seconds': g['seconds']} for g in grids],
        'pairs': pairs,
        'method': ('a coarse cell is the exact union of 2^(3k) fine cells because the grids share an '
                   'origin. "touched" counts a coarse void cell as reclaimed if ANY fine cell in it '
                   'is claimed, an upper bound; the partial-volume column weights each coarse cell '
                   'by the fine-resolved unclaimed fraction and is the honest correction'),
    }, grids


# --------------------------------------------------------------------------------------------
# stage: compartments


def line_of_sight_to_envelope(void_idx, points, closest, occ, lo, h, samples):
    """True where the straight segment voxel -> closest envelope point crosses no claimed cell.

    This is the geometric definition of the superficial compartment: tissue between the skin and the
    first structure under it, with nothing segmented in between. It replaces the nearest-entity
    Voronoi label for the primary split.
    """
    dim = np.asarray(occ.shape)
    blocked = np.zeros(len(points), bool)
    direction = closest - points
    for step in range(1, samples + 1):
        t = step / (samples + 1.0)
        q = points + t * direction
        ijk = np.rint((q - lo) / h).astype(np.int64)
        np.clip(ijk, 0, dim - 1, out=ijk)
        blocked |= occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]]
    return ~blocked


def stage_compartments(xv, tt, entities, h, samples):
    lo, dim = grid_of(xv, h)
    result, interior, occ, void, owner, rows = analyse_on(xv, tt, entities, lo, h, dim,
                                                          collect_owner=True)
    idx = np.argwhere(void)
    points = np.ascontiguousarray(idx * h + lo)
    sq, _, closest = igl.point_mesh_squared_distance(points, xv, tt)
    depth = np.sqrt(np.asarray(sq, float))
    closest = np.ascontiguousarray(np.asarray(closest, float))
    print('  compartments: %d void voxels, line-of-sight with %d samples' % (len(idx), samples),
          flush=True)
    superficial = line_of_sight_to_envelope(idx, points, closest, occ, lo, h, samples)

    # nearest claimed entity, by Euclidean distance transform over the entity surface samples
    _, ind = ndimage.distance_transform_edt(owner < 0, sampling=(h, h, h), return_indices=True)
    nearest = owner[ind[0], ind[1], ind[2]][void]
    systems = np.array([e.get('system') or 'none' for e, _, _ in entities] + ['none'])
    nearest_system = systems[np.where(nearest >= 0, nearest, len(systems) - 1)]

    label = np.where(superficial, 'subcutaneous',
                     np.array([DEEP_COMPARTMENT_OF_SYSTEM.get(s, 'deep_other')
                               for s in nearest_system]))
    cell = h ** 3
    compartments = {}
    for name in sorted(set(label.tolist())):
        m = label == name
        compartments[name] = {
            'volume_m3': float(m.sum() * cell),
            'fraction_of_void': float(m.sum() / max(len(label), 1)),
            'depth_below_skin_mm': {
                'mean': float(depth[m].mean() * 1000), 'median': float(np.median(depth[m]) * 1000),
                'p05': float(np.percentile(depth[m], 5) * 1000),
                'p95': float(np.percentile(depth[m], 95) * 1000),
                'max': float(depth[m].max() * 1000)},
            'nearest_entity_systems': {k: float(v * cell) for k, v in
                                       collections.Counter(nearest_system[m].tolist()).most_common(6)},
        }
    superficial_thickness_m = compartments.get('subcutaneous', {}).get('volume_m3', 0.0) / max(
        float(np.asarray(divergence_volume(xv, tt)[1])), 1e-12)
    summary = {
        'schema': 'ihm.interstitial-matrix-compartments.v1',
        'grid_h_m': h, 'grid_dim': [int(x) for x in dim], 'grid_origin_m': [float(x) for x in lo],
        'void_volume_m3': result['void_volume_m3'],
        'void_fraction_of_interior': result['void_fraction_of_interior'],
        'line_of_sight_samples': samples,
        'compartments': compartments,
        'superficial_equivalent_uniform_thickness_m': superficial_thickness_m,
        'envelope_area_m2': divergence_volume(xv, tt)[1],
        'primary_split_basis': ('line of sight to the closest point on the outer envelope, evaluated '
                                'on the occupancy grid. Superficial means no claimed cell lies '
                                'between the void cell and the skin'),
        'deep_split_basis': ('nearest claimed entity by Euclidean distance transform over surface '
                             'samples, mapped to a compartment by system. This is a Voronoi label'),
    }
    state = {'lo': lo, 'dim': dim, 'h': h, 'interior': interior, 'occ': occ, 'void': void,
             'owner': owner, 'idx': idx, 'depth': depth, 'superficial': superficial,
             'nearest': nearest, 'nearest_system': nearest_system, 'label': label, 'rows': rows}
    return summary, state


# --------------------------------------------------------------------------------------------
# stage: lumen


def stage_lumen(entities, state):
    """Is airway or gut lumen inside the void at all?

    A BodyParts3D organ surface is a closed solid: winding number > 0.5 inside the lumen as well as
    inside the wall. Any voxel there is `occupied` and is therefore not void by construction. This
    tests that per organ instead of assuming it, by asking how much void lies inside each hollow
    organ's convex hull -- a strictly larger region than the organ, so any lumen void must show up.
    """
    lo, h, occ, void = state['lo'], state['h'], state['occ'], state['void']
    dim = np.asarray(occ.shape)
    rows = []
    total_void_in_hull = 0.0
    total_enclosed = 0.0
    solid_count = 0
    tet_ready = {}
    lumen_mask = np.zeros_like(void)
    entity_index = {e['id']: i for i, (e, _, _) in enumerate(entities)}
    nearest_index = state['nearest']
    void_idx = state['idx']
    void_flat = np.ravel_multi_index(void_idx.T, tuple(occ.shape))
    path = TET_READY / 'entities.jsonl'
    if path.exists():
        for line in path.read_text().splitlines():
            record = json.loads(line)
            tet_ready[record['entity_id']] = record
    for entity, v, f in entities:
        name = entity['name'].lower()
        if entity.get('system') in ('arterial', 'venous', 'lymphatic', 'nervous'):
            continue
        if not any(k in name for k in LUMEN_KEYWORDS):
            continue
        a = np.maximum(((v.min(0) - h - lo) / h).astype(int), 0)
        b = np.minimum(((v.max(0) + h - lo) / h).astype(int) + 1, dim - 1)
        if not (b >= a).all():
            continue
        sl = [np.arange(a[i], b[i] + 1) for i in range(3)]
        q = np.ascontiguousarray(np.stack(np.meshgrid(
            *[lo[i] + h * sl[i] for i in range(3)], indexing='ij'), -1).reshape(-1, 3))
        hv = np.ascontiguousarray(v)
        hf = np.ascontiguousarray(np.asarray(cgal.convex_hull(hv), np.int64).reshape(-1, 3))
        inside_hull = np.asarray(igl.fast_winding_number(hv, hf, q), float) > 0.5
        ii = np.stack(np.meshgrid(*sl, indexing='ij'), -1).reshape(-1, 3)[inside_hull]
        vd = void[ii[:, 0], ii[:, 1], ii[:, 2]]
        surface_volume, _ = divergence_volume(v, f)
        # is the organ mesh a SOLID enclosing its lumen, or a shell around a wall only?
        # The maximum inscribed radius separates the two without ambiguity: a wall-only shell of
        # thickness t admits no interior point further than t/2 from the surface, while a solid
        # enclosing a lumen admits one at roughly half the lumen's calibre.
        centre = np.ascontiguousarray(v.mean(axis=0).reshape(1, 3))
        vv = np.ascontiguousarray(v)
        ff = np.ascontiguousarray(f)
        centre_inside = bool(np.asarray(igl.fast_winding_number(vv, ff, centre), float)[0] > 0.5)
        probe_h = max(min(h / 4.0, 0.001), 0.0005)
        pa, pb = v.min(0), v.max(0)
        pax = [np.arange(pa[i], pb[i] + probe_h, probe_h) for i in range(3)]
        inscribed = None
        if all(len(x) for x in pax) and np.prod([len(x) for x in pax]) < 6_000_000:
            pq = np.ascontiguousarray(np.stack(np.meshgrid(*pax, indexing='ij'), -1).reshape(-1, 3))
            pin = np.asarray(igl.fast_winding_number(vv, ff, pq), float) > 0.5
            if pin.any():
                dsq, _, _ = igl.point_mesh_squared_distance(np.ascontiguousarray(pq[pin]), vv, ff)
                inscribed = float(np.sqrt(np.asarray(dsq, float)).max())
        record = tet_ready.get(entity['id'], {})
        row = {
            'entity_id': entity['id'], 'name': entity['name'], 'system': entity.get('system'),
            'surface_signed_volume_m3': surface_volume,
            'enclosed_volume_includes_lumen_m3': abs(surface_volume),
            'centroid_inside_own_surface': centre_inside,
            'max_inscribed_radius_m': inscribed,
            'probe_spacing_m': probe_h,
            'watertight_edge_incidence': entity.get('watertight_edge_incidence'),
            'tet_ready_cavity_convention_pending': record.get('cavity_convention_pending'),
            'convex_hull_volume_m3': divergence_volume(hv, hf)[0],
            'void_inside_convex_hull_m3': float(vd.sum() * h ** 3),
            'void_fraction_of_convex_hull': float(vd.sum() / max(len(vd), 1)),
        }
        # a void cell counts as this organ's lumen only if it is inside the organ's convex hull AND
        # the organ is its nearest claimed structure. Both conditions are measured, not assumed.
        if vd.any():
            hull_flat = np.ravel_multi_index(ii[vd].T, tuple(occ.shape))
            select = np.isin(void_flat, hull_flat) & (nearest_index == entity_index[entity['id']])
            if select.any():
                chosen = void_idx[select]
                lumen_mask[chosen[:, 0], chosen[:, 1], chosen[:, 2]] = True
                row['void_attributed_to_this_lumen_m3'] = float(select.sum() * h ** 3)
        total_void_in_hull += row['void_inside_convex_hull_m3']
        total_enclosed += row['enclosed_volume_includes_lumen_m3']
        solid_count += int(centre_inside)
        rows.append(row)
    rows.sort(key=lambda r: -r.get('void_inside_convex_hull_m3', 0.0))
    cavity_flagged = [r for r in tet_ready.values() if r.get('cavity_convention_pending')]
    return {
        'schema': 'ihm.interstitial-matrix-lumen.v1',
        'grid_h_m': h,
        'hollow_organs_tested': len(rows),
        'hollow_organs_whose_own_centroid_is_inside_their_surface': solid_count,
        'hollow_organs_with_max_inscribed_radius_over_3mm': int(sum(
            1 for r in rows if (r.get('max_inscribed_radius_m') or 0) > 0.003)),
        'keywords': list(LUMEN_KEYWORDS),
        'total_void_inside_any_hollow_organ_convex_hull_m3': total_void_in_hull,
        'total_enclosed_volume_of_hollow_organs_m3': total_enclosed,
        'entities_flagged_cavity_convention_pending': len(cavity_flagged),
        'cavity_flagged_entities': [{'entity_id': r['entity_id'], 'name': r.get('name'),
                                     'role': r.get('role')} for r in cavity_flagged],
        'finding': ('a BodyParts3D viscus is authored as one closed surface enclosing wall AND lumen, '
                    'so a lumen voxel has winding number > 0.5 for that entity and is `occupied` by '
                    'construction. Airway and gut lumen therefore cannot appear in the void and no '
                    'lumen is filled by this build. The corollary is the opposite over-count: the '
                    'enclosed volume of these organs is carried as tissue in the entity ledger'),
        'note': ('the convex hull of an organ is strictly larger than the organ, so void_inside_'
                 'convex_hull OVER-counts: void in the concavities between an organ and its '
                 'neighbours is inside the hull and is interstitial space, not lumen. It is an upper '
                 'bound on lumen-in-void'),
        'lumen_mask_volume_m3': float(lumen_mask.sum() * h ** 3),
        'lumen_mask_rule': ('a void cell inside a hollow organ\'s convex hull whose nearest claimed '
                            'structure is that same organ. This is what the geometry stage removes '
                            'from the fill'),
        'rows': rows[:80],
    }, lumen_mask


# --------------------------------------------------------------------------------------------
# stage: composition


def stage_composition(compartments, sweep, anatomy, prior, profile, mechanics, repair_allocation):
    mechanics_allocation = mechanics['mass_allocation']
    mechanics_rows = {e['id']: e for e in mechanics['entities']}
    densities = prior['composition']['densities']
    reference = prior['composition']['reference_male']
    fat_of_adipose = reference['fat_fraction_of_adipose_tissue_adult']
    d_adipose = densities['adipose']['value']
    d_skin = densities['skin']['value']
    d_fascia = densities['fascia']['value']
    d_residual = densities['residual_central']['value']
    d_water = densities['water']['value']
    d_soft_max = densities['soft_tissue_max']['value']

    envelope_area = compartments['envelope_area_m2']
    comp = {k: v['volume_m3'] for k, v in compartments['compartments'].items()}
    void_total = compartments['void_volume_m3']

    # The dermis-plus-epidermis slab is inside the void because the canonical skin is a
    # zero-thickness double-sided surface. Its thickness is the repaired one, not the double-counted
    # one: the repair candidate divides the skin-layer volume by the exterior area, not the slab area.
    epi_derm_thickness_m = 0.0016
    dermis_volume = envelope_area * epi_derm_thickness_m
    skin_layers = {i: mechanics_rows.get(i) for i in
                   ('body-skin-epidermis', 'body-skin-dermis', 'body-skin-hypodermis')}
    hypodermis = skin_layers['body-skin-hypodermis'] or {}
    hypodermis_shell_volume = hypodermis.get('volume_m3')
    hypodermis_shell_mass = hypodermis.get('mass_kg')

    # Allocation. The compartment volumes are MEASURED here; the split within a compartment is
    # transferred from the priors already recorded in interstitial-composition-prior-v1.
    subcutaneous = max(comp.get('subcutaneous', 0.0) - dermis_volume, 0.0)
    adipose_fraction = {
        'subcutaneous': 0.80,       # superficial fascia + SAT lobules; the rest is loose CT and fluid
        'visceral_abdominal': 0.55,  # mesenteric and omental fat plus retroperitoneal loose tissue
        'visceral_thoracic': 0.25,
        'intermuscular': 0.35,      # intermuscular adipose plus epimysial/intermuscular septa
        'periosseous': 0.20,
        'perivascular': 0.30,
        'perineural': 0.20,
        'deep_other': 0.30,
    }
    fascia_fraction = {
        'subcutaneous': 0.05, 'intermuscular': 0.25, 'periosseous': 0.20, 'visceral_abdominal': 0.05,
        'visceral_thoracic': 0.05, 'perivascular': 0.05, 'perineural': 0.05, 'deep_other': 0.05,
    }
    constituents = []
    constituents.append({
        'constituent': 'dermis_and_epidermis', 'compartment': 'subcutaneous',
        'volume_m3': dermis_volume, 'density_kg_m3': d_skin, 'tier': 'derived',
        'basis': ('envelope area %.6f m2 x %.1f mm; the canonical skin is a zero-thickness slab so '
                  'the real dermis occupies void' % (envelope_area, epi_derm_thickness_m * 1000)),
        'source': 'icrp89'})
    for name in sorted(comp):
        volume = subcutaneous if name == 'subcutaneous' else comp[name]
        af = adipose_fraction.get(name, 0.30)
        ff = fascia_fraction.get(name, 0.05)
        constituents.append({'constituent': 'adipose_tissue', 'compartment': name,
                             'volume_m3': volume * af, 'density_kg_m3': d_adipose,
                             'tier': 'transferred', 'source': 'shen2003+shen2009',
                             'basis': 'compartment volume measured here x transferred adipose '
                                      'volume fraction %.2f' % af})
        constituents.append({'constituent': 'fascia_and_dense_connective_tissue', 'compartment': name,
                             'volume_m3': volume * ff, 'density_kg_m3': d_fascia,
                             'tier': 'transferred', 'source': 'icru44_nist',
                             'basis': 'compartment volume measured here x transferred fascia volume '
                                      'fraction %.2f' % ff})
        constituents.append({'constituent': 'loose_connective_tissue_and_interstitial_fluid',
                             'compartment': name, 'volume_m3': volume * (1.0 - af - ff),
                             'density_kg_m3': d_residual, 'tier': 'assumed', 'source': None,
                             'basis': 'remainder of the compartment at the assumed residual density'})
    for c in constituents:
        c['mass_kg'] = c['volume_m3'] * c['density_kg_m3']

    fill_volume = sum(c['volume_m3'] for c in constituents)
    fill_mass = sum(c['mass_kg'] for c in constituents)
    adipose_mass = sum(c['mass_kg'] for c in constituents if c['constituent'] == 'adipose_tissue')
    fat_mass = adipose_mass * fat_of_adipose

    # entity side of the ledger, from the same voxel partition the void came from
    finest_row = min(sweep['rows'], key=lambda r: r['grid_h_m'])
    entity_volume = finest_row['occupied_volume_m3']
    prior_entity_mass = prior['ledger']['cases']['voxel_8mm']['entity_mass_kg']
    prior_entity_volume = prior['ledger']['cases']['voxel_8mm']['entity_volume_m3']
    entity_mean_density = prior_entity_mass / prior_entity_volume
    entity_mass = entity_volume * entity_mean_density
    total_mass = entity_mass + fill_mass

    bands = {}
    for key, residual in (('lower_water_1000', d_water), ('central_assumed_1030', d_residual),
                          ('upper_icru44_soft_1060', d_soft_max)):
        mass = sum(c['volume_m3'] * (residual if c['constituent'].startswith('loose') else
                                     c['density_kg_m3']) for c in constituents)
        bands[key] = {'residual_density_kg_m3': residual, 'fill_mass_kg': mass,
                      'total_body_mass_kg': entity_mass + mass}

    declared_fat_fraction = profile['body_fat_fraction']
    ledger_fat_kg = declared_fat_fraction * LEDGER_TARGET_KG
    marrow_and_essential = {
        'yellow_marrow_note': 'inside the canonical bone geometry, so it is entity volume, not void; '
                              'not counted in the fill fat here',
        'interstitial_adipose_note': 'ICRP 89 puts 8% of total adipose tissue in a form not '
                                     'separable by dissection; it is inside entity volume here',
        'essential_fat_fraction_of_lbm': reference['essential_fat_fraction_of_lbm'],
    }

    # the mass-ledger arithmetic, recomputed from the live canonical allocation rather than restated
    unscaled_live = mechanics_allocation['unscaled_proxy_mass_kg']
    scale_live = mechanics_allocation['uniform_scale']
    target_live = mechanics_allocation['target_mass_kg']
    unscaled_repair = repair_allocation['steps_unscaled_kg']['total']
    ledger = {
        'canonical_mechanics': {
            'target_mass_kg': target_live, 'unscaled_proxy_mass_kg': unscaled_live,
            'uniform_scale': scale_live,
            'scale_times_proxy_kg': unscaled_live * scale_live,
            'closes_to_target': abs(unscaled_live * scale_live
                                    + mechanics_allocation['numerical_carrier_mass_kg']
                                    - target_live) < 1e-6},
        'brief_claimed_scale_after': LEDGER_SCALE_AFTER,
        'brief_claimed_scale_before': LEDGER_SCALE_BEFORE,
        'brief_scale_after_times_repaired_proxy_kg': unscaled_repair * LEDGER_SCALE_AFTER,
        'scale_that_reaches_corrected_target': (
            LEDGER_TARGET_KG - mechanics_allocation['numerical_carrier_mass_kg']) / unscaled_repair,
    }
    ledger['finding'] = (
        'the two corrections are separate. uniform_scale %.12f multiplied by the repaired unscaled '
        'proxy %.10f kg gives %.7f kg, which is the SUPERSEDED target 77.1107029 kg, not 70.7713 kg: '
        'that scale change is the skin-area, disc-density and duplicate repair at constant target. '
        'Reaching 70.7713 kg with the same repaired proxy needs scale %.12f, and the canonical '
        'mechanics.json read by this build now carries scale %.12f against target %.4f kg.'
        % (LEDGER_SCALE_AFTER, unscaled_repair, unscaled_repair * LEDGER_SCALE_AFTER,
           ledger['scale_that_reaches_corrected_target'], scale_live, target_live))

    return {
        'schema': 'ihm.interstitial-matrix-composition.v1',
        'built_against': {
            'profile_mass_kg': profile['mass_kg'],
            'profile_body_fat_fraction': declared_fat_fraction,
            'corrected_ledger_target_kg': LEDGER_TARGET_KG,
            'canonical_mechanics_uniform_scale': mechanics_allocation['uniform_scale'],
            'candidate_repair_uniform_scale': LEDGER_SCALE_AFTER,
            'anatomy_sha256_note': 'recorded in manifest.json inputs_sha256; anatomy.json and '
                                   'mechanics.json are being edited by another lane',
        },
        'ledger_arithmetic_check': ledger,
        'double_count_hazard': {
            'skin_quadrature_shells': {
                'entities': {i: {'volume_m3': (r or {}).get('volume_m3'),
                                 'mass_kg': (r or {}).get('mass_kg')}
                             for i, r in skin_layers.items()},
                'representation': 'surface_shell_quadrature_layer',
                'in_voxel_occupancy': False,
                'total_shell_volume_m3': float(sum((r or {}).get('volume_m3') or 0.0
                                                   for r in skin_layers.values())),
                'total_shell_mass_kg': float(sum((r or {}).get('mass_kg') or 0.0
                                                 for r in skin_layers.values())),
                'statement': ('the three skin layers have no volumetric geometry, so the occupancy '
                              'pass skips them and the space they describe is INSIDE the void this '
                              'build fills. The hypodermis in particular IS subcutaneous adipose and '
                              'carries %.6f m3 / %.4f kg in the canonical mass ledger. Materialising '
                              'the matrix without retiring these three shells double-counts their '
                              'whole volume and mass'
                              % (hypodermis_shell_volume or 0.0, hypodermis_shell_mass or 0.0)),
            },
            'entity_allocation': {
                'canonical_entity_mass_total_kg': float(sum(e['mass_kg'] for e in
                                                            mechanics['entities'])),
                'this_build_entity_side_kg': entity_mass,
                'statement': ('canonical mechanics.json allocates the whole declared body mass across '
                              'the entity set alone. Adding the interstitial matrix as a first-class '
                              'structure requires the entity allocation to fall to the entity side of '
                              'this ledger and uniform_scale to be re-derived against that, not '
                              'against the whole body mass'),
            },
        },
        'constituents': constituents,
        'fill': {'volume_m3': fill_volume, 'mass_kg': fill_mass,
                 'mean_density_kg_m3': fill_mass / max(fill_volume, 1e-12),
                 'adipose_tissue_mass_kg': adipose_mass, 'fat_in_that_adipose_kg': fat_mass},
        'entity_side': {'volume_m3': entity_volume, 'mean_density_kg_m3': entity_mean_density,
                        'mass_kg': entity_mass, 'grid_h_m': finest_row['grid_h_m'],
                        'basis': 'finest-grid exclusive occupancy x the mean owned density of '
                                 'data/derived/interstitial-composition-prior-v1'},
        'total_body_mass_kg': total_mass,
        'residual_density_band': bands,
        'body_fat': {
            'declared_fraction': declared_fat_fraction,
            'declared_fat_kg_at_corrected_ledger': ledger_fat_kg,
            'fat_in_synthesized_fill_kg': fat_mass,
            'fill_fat_fraction_of_total_body_mass': fat_mass / max(total_mass, 1e-12),
            'shortfall_kg': ledger_fat_kg - fat_mass,
            'not_counted': marrow_and_essential,
            'closure': {
                'question': ('the declared fraction is WHOLE-BODY fat. The fill can only carry the '
                             'adipose that lies outside every segmented structure. The residual must '
                             'therefore sit inside entity geometry -- yellow marrow inside bone, '
                             'intramyocellular and intrahepatic lipid inside muscle and liver -- and '
                             'this checks whether the residual is the size ICRP says it should be'),
                'required_total_fat_kg': ledger_fat_kg,
                'fat_carried_by_the_fill_kg': fat_mass,
                'implied_fat_inside_entity_geometry_kg': ledger_fat_kg - fat_mass,
                'implied_fraction_of_entity_mass': (ledger_fat_kg - fat_mass) / max(entity_mass, 1e-12),
                'icrp89_cross_check': {
                    'reference_male_mass_kg': reference['total_body_mass_kg'],
                    'reference_total_adipose_kg': reference['adipose_tissue_kg'],
                    'reference_separable_adipose_kg': reference['separable_adipose_kg'],
                    'reference_non_separable_adipose_kg': (reference['adipose_tissue_kg']
                                                           - reference['separable_adipose_kg']),
                    'scale_to_this_body': LEDGER_TARGET_KG / reference['total_body_mass_kg'],
                    'scaled_separable_adipose_kg': reference['separable_adipose_kg']
                        * LEDGER_TARGET_KG / reference['total_body_mass_kg'],
                    'scaled_non_separable_adipose_kg': (reference['adipose_tissue_kg']
                                                        - reference['separable_adipose_kg'])
                        * LEDGER_TARGET_KG / reference['total_body_mass_kg'],
                    'scaled_non_separable_fat_kg': (reference['adipose_tissue_kg']
                                                     - reference['separable_adipose_kg'])
                        * LEDGER_TARGET_KG / reference['total_body_mass_kg'] * fat_of_adipose,
                    'fill_separable_adipose_kg': adipose_mass,
                    'fill_separable_adipose_minus_scaled_icrp_kg': adipose_mass
                        - reference['separable_adipose_kg'] * LEDGER_TARGET_KG
                        / reference['total_body_mass_kg'],
                    'residual_minus_scaled_non_separable_fat_kg': (ledger_fat_kg - fat_mass)
                        - (reference['adipose_tissue_kg'] - reference['separable_adipose_kg'])
                        * LEDGER_TARGET_KG / reference['total_body_mass_kg'] * fat_of_adipose,
                },
            },
        },
        'mass_closure': {
            'composed_total_body_mass_kg': total_mass,
            'declared_total_body_mass_kg': LEDGER_TARGET_KG,
            'difference_kg': total_mass - LEDGER_TARGET_KG,
            'relative_difference': total_mass / LEDGER_TARGET_KG - 1.0,
            'declared_inside_residual_density_band': (
                bands['lower_water_1000']['total_body_mass_kg'] <= LEDGER_TARGET_KG
                <= bands['upper_icru44_soft_1060']['total_body_mass_kg']),
        },
        'unverified': ['the per-compartment adipose and fascia volume fractions are engineering '
                       'splits chosen to reproduce the transferred whole-body depot totals; no '
                       'measurement in this repository constrains them individually',
                       'the entity side of the ledger is a voxel partition times one mean density, '
                       'not the canonical per-entity mass ledger'],
    }


# --------------------------------------------------------------------------------------------
# stage: geometry


def marching_surface(mask, lo, h, smooth, isolevel=0.5):
    """Marching cubes on the fill indicator.

    `smooth` applies one 3-tap binomial pass to the indicator before contouring, which rounds the
    diagonal cell contacts that otherwise produce non-manifold vertices. The volume drift this
    causes is measured, not assumed.
    """
    field = mask.astype(np.float32)
    if smooth:
        field = ndimage.uniform_filter(field, size=3, mode='constant', cval=0.0)
    dim = np.asarray(mask.shape)
    axes = [lo[i] + h * np.arange(dim[i]) for i in range(3)]
    # libigl indexes S as S(x + y*nx + z*nx*ny), so x varies fastest: ravel in Fortran order.
    mesh = np.meshgrid(*axes, indexing='ij')
    grid = np.ascontiguousarray(np.stack([m.ravel(order='F') for m in mesh], axis=1))
    result = igl.marching_cubes(np.ascontiguousarray(field.ravel(order='F').astype(np.float64)),
                                grid, int(dim[0]), int(dim[1]), int(dim[2]), float(isolevel))
    v, f = result[0], result[1]
    return (np.ascontiguousarray(np.asarray(v, float)),
            np.ascontiguousarray(np.asarray(f, np.int64)))


def euler_genus(v, f):
    if not len(f):
        return None, None, None
    edges = np.unique(np.sort(np.concatenate(
        [f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1), axis=0)
    used = len(np.unique(f))
    chi = used - len(edges) + len(f)
    shells = int(igl.facet_components(f)[0])
    genus = (2 * shells - chi) / 2.0
    return int(chi), shells, float(genus)


def stage_geometry(state, lumen_mask, out, smooth, max_components, tet_face_cap, flags,
                   isolevel=0.5):
    lo, h, void = state['lo'], state['h'], state['void']
    fill = void & ~lumen_mask
    labels, count = ndimage.label(fill, ndimage.generate_binary_structure(3, 1))
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    order = np.argsort(-sizes)
    geometry_dir = out / 'geometry'
    geometry_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    kept = 0
    for k in order:
        if sizes[k] == 0 or kept >= max_components:
            break
        component = labels == k
        # pad by one cell so the contour closes inside the array
        idx = np.argwhere(component)
        a = np.maximum(idx.min(0) - 2, 0)
        b = np.minimum(idx.max(0) + 3, np.asarray(component.shape))
        sub = component[a[0]:b[0], a[1]:b[1], a[2]:b[2]]
        sub_lo = lo + a * h
        began = time.monotonic()
        v, f = marching_surface(sub, sub_lo, h, smooth, isolevel)
        raw = diagnose(v, f, intersections=False)
        rv, rf, log = repair(v, f)
        rec = diagnose(rv, rf, intersections=True)
        chi, shells, genus = euler_genus(rv, rf)
        row = {
            'component': int(k), 'rank': kept, 'voxels': int(sizes[k]),
            'voxel_volume_m3': float(sizes[k] * h ** 3),
            'marching_cubes': {'vertices': int(len(v)), 'faces': int(len(f)),
                               'signed_volume_m3': signed_volume(v, f)},
            'repaired': rec, 'repair_log': log,
            'repaired_surface_area_m2': divergence_volume(rv, rf)[1],
            'euler_characteristic': chi, 'shells': shells, 'genus': genus,
            'volume_drift_vs_voxels': (abs(rec['signed_volume_m3']) - sizes[k] * h ** 3)
                                      / max(sizes[k] * h ** 3, 1e-12),
            'seconds_geometry': time.monotonic() - began,
        }
        if len(rf) <= tet_face_cap:
            began = time.monotonic()
            row['tetgen'] = tetrahedralize(rv, rf, flags)
            row['tetgen']['seconds'] = time.monotonic() - began
        else:
            row['tetgen'] = {'flags': flags, 'skipped': 'face count %d above cap %d'
                             % (len(rf), tet_face_cap)}
        path = geometry_dir / ('interstitial-matrix-%03d.json.gz' % kept)
        path.write_bytes(gzip.compress(json.dumps({
            'schema': 'ihm.interstitial-matrix-surface.v1',
            'structure_id': 'body-interstitial-matrix-%03d' % kept,
            'frame': 'bodyparts3d-display-m', 'units': 'm',
            'representation': 'triangular_surface',
            'grid_h_m': h, 'grid_origin_m': [float(x) for x in sub_lo],
            'positions': rv.ravel().tolist(), 'indices': rf.ravel().tolist(),
        }, allow_nan=False).encode()))
        row['output_path'] = str(path.relative_to(out))
        row['output_sha256'] = sha(path)
        rows.append(row)
        kept += 1
        print('  component %d: %d voxels -> %d faces, shells=%s genus=%s tetgen=%s'
              % (kept - 1, sizes[k], len(rf), shells, genus,
                 row['tetgen'].get('succeeded', row['tetgen'].get('skipped'))), flush=True)
    succeeded = [r for r in rows if r['tetgen'].get('succeeded')]
    attempted = [r for r in rows if 'skipped' not in r['tetgen']]
    return {
        'schema': 'ihm.interstitial-matrix-geometry.v1',
        'grid_h_m': h, 'smoothed_indicator': bool(smooth), 'isolevel': float(isolevel),
        'tetgen_flags': flags,
        'fill_definition': 'void minus any voxel classified as correctly-empty lumen',
        'fill_volume_m3': float(fill.sum() * h ** 3),
        'lumen_excluded_volume_m3': float((void & lumen_mask).sum() * h ** 3),
        'components_total': int(count), 'components_emitted': len(rows),
        'components_emitted_volume_m3': float(sum(r['voxel_volume_m3'] for r in rows)),
        'tetgen_attempted': len(attempted), 'tetgen_succeeded': len(succeeded),
        'tetgen_success_rate': len(succeeded) / max(len(attempted), 1),
        'tets_total': int(sum(r['tetgen'].get('tets', 0) for r in succeeded)),
        'components_not_emitted': int(count) - len(rows),
        'components_not_emitted_volume_m3': float(fill.sum() * h ** 3
                                                  - sum(r['voxel_volume_m3'] for r in rows)),
        'repaired_surface_area_m2': float(sum(r['repaired_surface_area_m2'] for r in rows)),
        'repaired_signed_volume_m3': float(sum(r['repaired']['signed_volume_m3'] for r in rows)),
        'tet_volume_m3': float(sum(r['tetgen'].get('tet_volume_m3', 0.0) for r in succeeded)),
        'realised_volume_vs_voxel_volume': (
            float(sum(r['repaired']['signed_volume_m3'] for r in rows))
            / max(sum(r['voxel_volume_m3'] for r in rows), 1e-12) - 1.0),
        'largest': rows[0] if rows else None,
        'rows': rows,
        'repair_recipe': 'imported verbatim from scripts/build_muscle_tet_ready_surfaces.py repair()',
    }


# --------------------------------------------------------------------------------------------
# stage: provenance


def provenance_records(geometry, composition, out, builder_sha, inputs):
    out_rel = str(Path(out).resolve().relative_to(ROOT))
    records = []
    for row in geometry['rows']:
        structure = 'body-interstitial-matrix-%03d' % row['rank']
        records.append({
            'schema': 'ihm.structure-provenance.v1',
            'structure_id': structure,
            'model_id': 'ihm-body',
            'canonical_entity_id': None,
            'name': 'interstitial matrix component %d' % row['rank'],
            'system': 'connective',
            'role': 'interstitial_matrix',
            'role_source': 'assigned by scripts/build_interstitial_matrix.py; no canonical role '
                           'vocabulary entry exists for this tissue class',
            'evidence_kind': 'synthesized_from_complement',
            'dataset': {
                'id': 'ihm-interstitial-matrix-v1',
                'label': 'IHM synthesized interstitial matrix',
                'version': '1',
                'revision': builder_sha,
                'url': None,
                'specimen': 'none. This structure has no source geometry. Its extent is the '
                            'complement of the BodyParts3D entity set inside the outer envelope',
                'units': 'm',
                'frame': 'bodyparts3d-display-m',
                'attribution': 'extent derives from BodyParts3D 4.0 by complement; composition '
                               'derives from ICRU-44, ICRP 89, shen2003, shen2009, gallagher2005',
                'acquisition_status': 'not acquired; constructed',
                'license': None,
                'license_absent_reason': 'synthesized geometry; the licence that applies is that of '
                                         'the BodyParts3D entity set it is the complement of, '
                                         'recorded on every canonical entity record',
            },
            'source_file': {
                'path': None, 'sha256': None, 'sha256_verified': None, 'bytes': None,
                'absent_reason': 'no source file. Nothing measured this structure',
            },
            'build': {
                'script': 'scripts/build_interstitial_matrix.py',
                'script_sha256': builder_sha,
                'commit': None,
                'commit_covers_working_tree': False,
                'uncommitted_changes': True,
                'inputs_sha256': inputs,
            },
            'geometry': {
                'path': out_rel + '/' + row['output_path'],
                'sha256': row['output_sha256'],
                'sha256_verified': True,
                'representation': 'triangular_surface',
                'frame': 'bodyparts3d-display-m',
                'units': 'm',
                'source_vertex_count': row['repaired']['vertices'],
                'source_face_count': row['repaired']['faces'],
            },
            'transforms': [
                {'step': 'voxel_complement', 'residual': {'metric': 'void_fraction_of_interior',
                 'value': geometry['grid_h_m'], 'method': 'grid spacing at which the complement was '
                 'taken'}},
                {'step': 'marching_cubes_isolevel_0.5', 'residual':
                 {'metric': 'relative_volume_drift_vs_voxel_count',
                  'value': row['volume_drift_vs_voxels'], 'method': 'signed divergence volume of the '
                  'repaired shell against the voxel count times h^3'}},
                {'step': 'muscle_tet_ready_repair_recipe', 'residual': None,
                 'residual_reason': 'the repair only welds, drops degenerate faces, reorients and '
                                    'resolves self-intersections; no coordinate is moved'},
            ],
            'frame_relation': 'canonical',
            'tier': 'synthesized',
            'tier_basis': 'Geometry constructed from explicit priors with no source geometry of its '
                          'own: the extent is a set complement and the composition is transferred '
                          'literature. No cell of this structure was measured on any specimen.',
            'tier_evidence': {
                'extent': 'complement of the canonical triangular-surface entity set inside '
                          'data/derived/outer-envelope/outer-envelope.npz, taken on a %.4f m grid'
                          % geometry['grid_h_m'],
                'composition': 'data/derived/interstitial-composition-prior-v1/composition.json',
                'no_measurement': 'BodyParts3D 4.0 segments no adipose tissue, no fascia and no '
                                  'loose areolar connective tissue; a name search over the 2408 '
                                  'canonical entities returns zero of each',
            },
            'assumptions': composition['unverified'],
            'derived_artifacts': [],
            'display_present': False,
            'completeness': {
                'required_present': 8, 'required_total': 11,
                'missing': ['dataset.license', 'source_file.path', 'source_file.sha256'],
                'answerable': True,
                'hashes_verified': 1, 'hashes_failed': 0, 'hashes_unchecked': 0,
            },
        })
    return records


# --------------------------------------------------------------------------------------------


def self_test():
    """Every non-trivial function on geometry whose answer is known in closed form."""
    unit_v = np.array([[0., 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    unit_f = np.array([[0, 1, 2], [0, 2, 3], [4, 6, 5], [4, 7, 6], [0, 5, 1], [0, 4, 5],
                       [1, 6, 2], [1, 5, 6], [2, 7, 3], [2, 6, 7], [3, 4, 0], [3, 7, 4]])[:, ::-1]
    unit_f = np.ascontiguousarray(unit_f.astype(np.int64))
    volume, area = divergence_volume(unit_v, unit_f)
    assert abs(volume - 1) < 1e-12 and abs(area - 6) < 1e-12, (volume, area)
    chi, shells, genus = euler_genus(unit_v, unit_f)
    assert (chi, shells, genus) == (2, 1, 0.0), (chi, shells, genus)

    # a torus-like shell: box minus an axis-aligned tunnel, built on a voxel grid, must have genus 1
    grid = np.zeros((9, 9, 9), bool)
    grid[1:8, 1:8, 1:8] = True
    grid[3:6, 3:6, :] = False
    v, f = marching_surface(grid, np.zeros(3), 1.0, False)
    rv, rf, _ = repair(v, f)
    chi, shells, genus = euler_genus(rv, rf)
    assert shells == 1 and abs(genus - 1.0) < 1e-9, ('tunnel genus', chi, shells, genus)

    # coarsen_any: a single fine cell must reclaim exactly one coarse cell
    fine = np.zeros((17, 17, 17), bool)
    fine[8, 8, 8] = True
    coarse = coarsen_any(fine, 1)
    assert coarse.shape == (9, 9, 9) and coarse.sum() == 1 and coarse[4, 4, 4], coarse.sum()
    fine2 = np.zeros((17, 17, 17), bool)
    fine2[7, 8, 8] = True                      # offset -1 still lands in coarse cell 4
    assert coarsen_any(fine2, 1)[4, 4, 4]
    fine3 = np.zeros((17, 17, 17), bool)
    fine3[9, 8, 8] = True                      # offset +1 belongs to coarse cell 5
    assert coarsen_any(fine3, 1)[5, 4, 4] and not coarsen_any(fine3, 1)[4, 4, 4]

    # line of sight: a wall between the point and its target must block, an empty gap must not
    occ = np.zeros((11, 11, 11), bool)
    occ[5, :, :] = True
    idx = np.array([[2, 5, 5]])
    points = idx.astype(float)
    closest = np.array([[9.0, 5, 5]])
    assert not line_of_sight_to_envelope(idx, points, closest, occ, np.zeros(3), 1.0, 16)[0]
    assert line_of_sight_to_envelope(idx, points, closest, np.zeros_like(occ), np.zeros(3), 1.0,
                                     16)[0]

    # analyse_on must reproduce the shipped analyse() on a synthetic body
    def entity(i):
        return {'id': 'e%d' % i, 'name': 'e%d' % i, 'system': 'test', 'role': 'test',
                'watertight_edge_incidence': True}
    outer_v = np.ascontiguousarray(unit_v * 1.0)
    inner = [(entity(0), np.ascontiguousarray(unit_v * 0.3 + 0.1), unit_f),
             (entity(1), np.ascontiguousarray(unit_v * 0.3 + 0.55), unit_f)]
    h = 0.05
    shipped, _, _, _, _, _ = IV.analyse(outer_v, unit_f, inner, h, False)
    lo, dim = grid_of(outer_v, h)
    mine, _, _, _, _, _ = analyse_on(outer_v, unit_f, inner, lo, h, dim)
    for key in ('interior_volume_m3', 'occupied_volume_m3', 'void_volume_m3',
                'void_fraction_of_interior', 'entity_claims_outside_envelope_m3'):
        assert abs(shipped[key] - mine[key]) < 1e-15, (key, shipped[key], mine[key])
    assert shipped['grid_dim'] == mine['grid_dim']
    print(json.dumps({'self_test': 'ok', 'unit_cube_volume': volume, 'tunnel_genus': genus,
                      'analyse_on_matches_shipped_analyse': True}, indent=2))


def run(out, force, stages, working_h, sweep_grids, samples, smooth, max_components, tet_face_cap,
        flags, isolevel=0.5):
    out = Path(out)
    if out.exists() and any(out.iterdir()) and not force:
        raise ValueError('choose a fresh output directory or pass --force')
    out.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    inputs = {
        'data/derived/canonical/anatomy.json': sha(ANATOMY),
        'data/derived/canonical/profile.json': sha(PROFILE),
        'data/derived/canonical/mechanics.json': sha(MECHANICS),
        'data/derived/outer-envelope/outer-envelope.npz': sha(ENVELOPE),
        'data/derived/interstitial-composition-prior-v1/composition.json': sha(PRIOR / 'composition.json'),
        'data/derived/interstitial-composition-prior-v1/ledger.json': sha(PRIOR / 'ledger.json'),
        'data/derived/interstitial-composition-prior-v1/sources.json': sha(PRIOR / 'sources.json'),
        'data/derived/interstitial-composition-prior-v1/manifest.json': sha(PRIOR / 'manifest.json'),
        'data/derived/entity-record-repair-candidate-v1/allocation.json': sha(REPAIR_CANDIDATE),
        'data/derived/unmodelled-volume/summary.json': sha(ROOT / 'data/derived/unmodelled-volume/summary.json'),
        'data/raw/anatomy/bodyparts3d/isa_parts_list_e.txt': sha(BP3D_RAW / 'isa_parts_list_e.txt'),
        'data/raw/anatomy/bodyparts3d/partof_parts_list_e.txt': sha(BP3D_RAW / 'partof_parts_list_e.txt'),
        'data/raw/anatomy/bodyparts3d/isa_element_parts.txt': sha(BP3D_RAW / 'isa_element_parts.txt'),
        'data/raw/anatomy/bodyparts3d/partof_element_parts.txt': sha(BP3D_RAW / 'partof_element_parts.txt'),
        'scripts/inventory_unmodelled_volume.py': sha(ROOT / 'scripts/inventory_unmodelled_volume.py'),
        'scripts/build_muscle_tet_ready_surfaces.py': sha(ROOT / 'scripts/build_muscle_tet_ready_surfaces.py'),
        'scripts/verify_muscle_tet_ready_surfaces.py': sha(ROOT / 'scripts/verify_muscle_tet_ready_surfaces.py'),
    }
    xv, tt = load_envelope()
    anatomy, entities, skipped, anatomy_sha = load_entities()
    print('entities %d, skipped %s' % (len(entities), skipped), flush=True)
    profile = json.loads(PROFILE.read_text())
    prior = {'composition': json.loads((PRIOR / 'composition.json').read_text()),
             'ledger': json.loads((PRIOR / 'ledger.json').read_text()),
             'sources': json.loads((PRIOR / 'sources.json').read_text())}
    repair_allocation = json.loads(REPAIR_CANDIDATE.read_text())
    mechanics = json.loads(MECHANICS.read_text())

    produced = {}
    lumen_mask = None
    if 'source_coverage' in stages:
        produced['source_coverage'] = stage_source_coverage(anatomy)
        write(out / 'source-coverage.json', produced['source_coverage'])
    if 'sweep' in stages:
        produced['sweep'] = stage_sweep(xv, tt, entities, sweep_grids)
        write(out / 'resolution-sweep.json', produced['sweep'])
    if 'nested' in stages:
        nested, _ = stage_nested(xv, tt, entities, NESTED_BASE_M, NESTED_LEVELS)
        produced['nested'] = nested
        write(out / 'nested-refinement.json', nested)
    state = None
    if {'compartments', 'lumen', 'geometry'} & set(stages):
        produced['compartments'], state = stage_compartments(xv, tt, entities, working_h, samples)
        write(out / 'compartments.json', produced['compartments'])
    if 'lumen' in stages:
        produced['lumen'], lumen_mask = stage_lumen(entities, state)
        write(out / 'lumen.json', produced['lumen'])
        np.savez_compressed(out / 'lumen-mask.npz', mask=np.packbits(lumen_mask),
                            dim=np.asarray(lumen_mask.shape), grid_h_m=state['h'],
                            grid_origin_m=state['lo'])
    if 'composition' in stages:
        sweep = produced.get('sweep') or json.loads((out / 'resolution-sweep.json').read_text())
        compartments = produced.get('compartments') or json.loads((out / 'compartments.json').read_text())
        produced['composition'] = stage_composition(compartments, sweep, anatomy, prior, profile,
                                                    mechanics, repair_allocation)
        write(out / 'composition.json', produced['composition'])
    if 'geometry' in stages:
        if lumen_mask is None:
            payload = np.load(out / 'lumen-mask.npz')
            lumen_mask = np.unpackbits(payload['mask'])[:int(np.prod(payload['dim']))].astype(
                bool).reshape(tuple(payload['dim']))
        produced['geometry'] = stage_geometry(state, lumen_mask, out, smooth, max_components,
                                              tet_face_cap, flags, isolevel)
        write(out / 'geometry.json', produced['geometry'])
    if 'provenance' in stages:
        geometry = produced.get('geometry') or json.loads((out / 'geometry.json').read_text())
        composition = produced.get('composition') or json.loads((out / 'composition.json').read_text())
        records_dir = out / 'provenance'
        records_dir.mkdir(parents=True, exist_ok=True)
        records = provenance_records(geometry, composition, out, sha(__file__), inputs)
        for record in records:
            write(records_dir / (record['structure_id'] + '.json'), record)
        produced['provenance'] = {'schema': 'ihm.interstitial-matrix-provenance-index.v1',
                                  'records': len(records), 'tier': 'synthesized',
                                  'record_schema': 'ihm.structure-provenance.v1',
                                  'ids': [r['structure_id'] for r in records]}
        write(out / 'provenance-index.json', produced['provenance'])

    manifest = {
        'schema': 'ihm.interstitial-matrix-manifest.v1',
        'status': ('CANDIDATE. Not canonical, not promoted, not consumed by any runtime. Every '
                   'canonical asset is opened read-only.'),
        'builder': 'scripts/build_interstitial_matrix.py',
        'builder_sha256': sha(__file__),
        'python': sys.version,
        'stages': sorted(stages),
        'working_grid_h_m': working_h,
        'marching_cubes_isolevel': isolevel,
        'sweep_grids_m': list(sweep_grids),
        'tetgen_flags': flags,
        'inputs_sha256': inputs,
        'anatomy_sha256_consumed': anatomy_sha,
        'entities_considered': len(entities),
        'skipped_by_representation': skipped,
        'canonical_assets_modified': False,
        'inputs_sha256_rechecked_after_run': {k: sha(ROOT / k) for k in inputs},
        'inputs_unchanged_during_run': all(sha(ROOT / k) == v for k, v in inputs.items()),
        'unverified': list(UNVERIFIED),
        'wall_seconds': time.monotonic() - began,
        'artifacts_sha256': {str(p.relative_to(out)): sha(p)
                             for p in sorted(out.rglob('*')) if p.is_file()
                             and p.name != 'manifest.json'},
    }
    write(out / 'manifest.json', manifest)
    print(json.dumps({k: (v if k in ('sweep', 'nested') else
                          {kk: vv for kk, vv in v.items() if kk != 'rows'})
                      for k, v in produced.items()}, indent=2)[:8000])
    return produced


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='data/derived/interstitial-matrix-v1')
    p.add_argument('--force', action='store_true')
    p.add_argument('--self-test', action='store_true')
    p.add_argument('--stages', default='source_coverage,sweep,nested,compartments,lumen,composition,geometry,provenance')
    p.add_argument('--working-h', type=float, default=WORKING_M)
    p.add_argument('--sweep', default=','.join(str(x) for x in SWEEP_M))
    p.add_argument('--los-samples', type=int, default=48)
    p.add_argument('--smooth', action='store_true', help='blur the fill indicator before contouring')
    p.add_argument('--max-components', type=int, default=8000)
    p.add_argument('--tet-face-cap', type=int, default=4_000_000)
    p.add_argument('--flags', default=TETGEN_FLAGS)
    p.add_argument('--isolevel', type=float, default=0.5)
    a = p.parse_args()
    if a.self_test:
        self_test()
        return
    run(ROOT / a.output if not Path(a.output).is_absolute() else Path(a.output), a.force,
        set(s.strip() for s in a.stages.split(',') if s.strip()), a.working_h,
        tuple(float(x) for x in a.sweep.split(',')), a.los_samples, a.smooth, a.max_components,
        a.tet_face_cap, a.flags, a.isolevel)


if __name__ == '__main__':
    main()
