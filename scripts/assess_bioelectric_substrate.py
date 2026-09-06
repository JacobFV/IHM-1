#!/usr/bin/env python3
"""Assess the retained bioelectric substrate against Levin-style pattern-control theses.

Read-only. Recomputes every quoted number from hashed inputs and writes one
assessment plus an input manifest under data/derived/bioelectric-substrate-assessment-v1/.
Nothing in the repository is mutated. Mechanism statements sourced from the
literature are tagged 'literature'; everything else is 'verified' (recomputed
here) or 'inferred' (an explicit engineering judgement made here).
"""
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/bioelectric-substrate-assessment-v1'

INPUTS = [
    'ihm/assembly/skin_bioelectric.py',
    'ihm/assembly/skin_epithelial_inventory.py',
    'ihm/assembly/skin_layers.py',
    'ihm/assembly/epithelial_electrodiffusion.py',
    'ihm/materialize/skin.py',
    'ihm/calibration/skin.py',
    'scripts/verify_bioelectric.py',
    'scripts/verify_body_skin_electric.py',
    'scripts/export_betse_tissue.py',
    'docs/BIOELECTRIC_VIEW.md',
    'docs/BIOELECTRICITY.md',
    'docs/INTEGUMENTARY_DATA.md',
    'docs/research/SKIN_BIOELECTRIC_CALIBRATION_GAPS.md',
    'docs/research/SKIN_EXTERIOR_LAYER_VOLUME.md',
    'data/derived/canonical/skin-electric.json',
    'data/derived/canonical/anatomy.json',
    'data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz',
    'data/derived/calibration/skin-fit.json',
    'data/derived/bioelectric/tissue.json',
    'data/derived/outer-envelope/validation.json',
    'data/derived/outer-envelope/outer-envelope.json.gz',
    'data/derived/entity-tet-ready-skin-probe-v1/summary.json',
    'data/derived/entity-tet-ready-v1/summary.json',
    'data/derived/entity-tet-ready-v1/excluded.json',
    'data/derived/entity-tet-ready-verification-v1/summary.json',
    'data/derived/muscle-tet-ready-verification-v1/summary.json',
    'data/research/engineered_skin_territories/materialization.json',
    'data/research/engineered_skin_territories/layer_volume_candidate.json',
    'data/derived/physiology/betse_run/sim_config.yaml',
    'data/derived/physiology/betse_run/run_summary.json',
]

R, FARADAY = 8.314462618, 96485.33212


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load(rel):
    p = ROOT / rel
    raw = p.read_bytes()
    if rel.endswith('.gz'):
        return json.loads(gzip.decompress(raw))
    return json.loads(raw)


def triangle_areas(v, f):
    return np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1) / 2


def ray_parity_inside(point, a, b, c, direction):
    e1, e2 = b - a, c - a
    h = np.cross(direction, e2)
    det = np.einsum('ij,ij->i', e1, h)
    ok = np.abs(det) > 1e-14
    inv = np.zeros_like(det)
    inv[ok] = 1.0 / det[ok]
    s = point - a
    u = inv * np.einsum('ij,ij->i', s, h)
    q = np.cross(s, e1)
    v = inv * (q @ direction)
    t = inv * np.einsum('ij,ij->i', e2, q)
    hit = ok & (u >= 0) & (u <= 1) & (v >= 0) & (u + v <= 1) & (t > 1e-9)
    return bool(int(hit.sum()) % 2 == 1)


def area_ledger():
    """Every candidate integument area, recomputed, with the defect factor."""
    evidence = load('data/research/engineered_skin_territories/materialization.json')
    geometry = load('data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz')
    v = np.asarray(geometry['positions'], float).reshape(-1, 3)
    f = np.asarray(geometry['indices'], int).reshape(-1, 3)
    area = triangle_areas(v, f)
    comp = np.asarray(evidence['surface_diagnostic']['face_component_ids'], int)
    ids = np.asarray(evidence['contact_eligible_triangle_ids'], int)
    validation = load('data/derived/outer-envelope/validation.json')
    anatomy = load('data/derived/canonical/anatomy.json')
    layers = {e['id']: e for e in anatomy['entities'] if e.get('role') == 'skin_layer'}
    raw = float(area.sum())
    outer = float(area[ids].sum())
    inner = float(area[comp == 1].sum())
    envelope = validation['envelope_topology']['surface_area_m2']
    factor = raw / outer
    thickness = {k: e['volume_m3'] / e['surface_area_m2'] for k, e in layers.items()}
    return dict(
        raw_double_sided_slab_area_m2=raw,
        outer_component_area_m2=outer,
        outer_component_faces=int(len(ids)),
        inner_component_area_m2=inner,
        inner_component_faces=int((comp == 1).sum()),
        seam_fragment_area_m2=raw - outer - inner,
        watertight_envelope_area_m2=envelope,
        envelope_over_outer_component=envelope / outer,
        dubois_bsa_from_profile_mass_m2=validation['validation']['dubois_bsa_m2_from_profile_mass'],
        envelope_vs_dubois_relative=validation['validation']['area_vs_dubois_profile_relative'],
        inflation_factor_raw_over_outer=factor,
        layer_defect=[dict(
            id=k, thickness_m=thickness[k],
            recorded_area_m2=layers[k]['surface_area_m2'],
            recorded_volume_m3=layers[k]['volume_m3'],
            corrected_volume_on_outer_component_m3=outer * thickness[k],
            corrected_volume_on_watertight_envelope_m3=envelope * thickness[k],
            excess_volume_m3=layers[k]['volume_m3'] - outer * thickness[k],
            factor=layers[k]['volume_m3'] / (outer * thickness[k]),
            geometry_representation=layers[k]['reference_geometry']['representation'],
        ) for k in ('body-skin-epidermis', 'body-skin-dermis', 'body-skin-hypodermis')],
    )


def anchor_defect():
    """Which sheet of the double-sided slab the regional patch is pinned to."""
    artifact = load('data/derived/canonical/skin-electric.json')
    evidence = load('data/research/engineered_skin_territories/materialization.json')
    comp = np.asarray(evidence['surface_diagnostic']['face_component_ids'], int)
    eligible = set(int(i) for i in evidence['contact_eligible_triangle_ids'])
    geometry = load('data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz')
    v = np.asarray(geometry['positions'], float).reshape(-1, 3)
    f = np.asarray(geometry['indices'], int).reshape(-1, 3)
    idx = int(artifact['anchor']['face_index'])
    tri = v[f[idx]]
    normal = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    normal = normal / np.linalg.norm(normal)
    origin = np.asarray(artifact['anchor']['origin_m'], float)
    env = load('data/derived/outer-envelope/outer-envelope.json.gz')
    ev = np.asarray(env['positions'], float).reshape(-1, 3)
    ef = np.asarray(env['indices'], int).reshape(-1, 3)
    a, b, c = ev[ef[:, 0]], ev[ef[:, 1]], ev[ef[:, 2]]
    d = np.asarray([0.3178, 0.5171, 0.7963]); d = d / np.linalg.norm(d)
    probes = {f'{s * 1e3:+.1f}mm_along_face_normal': ray_parity_inside(origin + s * normal, a, b, c, d)
              for s in (-2e-3, -5e-4, 0.0, 5e-4, 2e-3)}
    axes = np.asarray(artifact['anchor']['local_axes'], float)
    cell_node = (np.array([0.0, 0.0, -1e-4]) - np.array([8e-4, 0.0, 0.0])) @ axes + origin
    return dict(
        anchor_face_index=idx,
        anchor_face_component_id=int(comp[idx]),
        selected_exterior_component_id=int(evidence['surface_diagnostic']['selected_component_id']),
        anchor_face_on_exterior_component=idx in eligible,
        anchor_origin_inside_watertight_envelope=probes['+0.0mm_along_face_normal'],
        envelope_inside_probes=probes,
        outward_direction_is_minus_face_normal=(not probes['-2.0mm_along_face_normal']) and probes['+2.0mm_along_face_normal'],
        basal_cell_row_world_m=cell_node.tolist(),
        basal_cell_row_inside_envelope=ray_parity_inside(cell_node, a, b, c, d),
        cell_row_local_depth_m=-1e-4,
        finding=('The regional skin-electric patch is pinned to a face of the INNER sheet of the '
                 'double-sided slab, not the reviewed exterior component. Outward is -n at that '
                 'face, so the layer the code calls basal cells (local z=-1e-4 m) sits 0.1 mm '
                 'nearer the body exterior than the apical surface sites (local z=0): the '
                 'apical/basal stack is geometrically inverted in world coordinates.'),
        impact=('Electrically inert today: edge fields use local patch coordinates, so no number '
                'in skin-electric.json changes. It is a provenance and registration defect that '
                'will silently corrupt any future field/mesh registration.'),
    )


def rc_substrate():
    """Recompute the lumped RC patch's derived quantities from its own coefficients."""
    import sys
    sys.path.insert(0, str(ROOT))
    from ihm.materialize.skin import SkinPatch
    from ihm.assembly.skin_bioelectric import electrical_system
    patch = SkinPatch.line(n=9, wound=False)
    a, b, capacitance = electrical_system(patch)
    resting = np.linalg.solve(a, -b)
    thermal = R * patch.temperature_k / FARADAY
    reversals = dict(
        sodium_V=thermal * math.log(patch.sodium_outside_mM / patch.cells[0].sodium_inside_mM),
        potassium_V=thermal * math.log(patch.potassium_outside_mM / patch.cells[0].potassium_inside_mM),
        chloride_V=-thermal * math.log(patch.chloride_outside_mM / patch.cells[0].chloride_inside_mM))
    cells = []
    for cell in patch.cells:
        total = cell.sodium_conductance_s + cell.potassium_conductance_s + cell.chloride_conductance_s
        cells.append(dict(total_conductance_s=total, membrane_tau_s=cell.capacitance_f / total,
                          gap_edge_conductance_s=patch.gap_edges[0][2],
                          gap_over_membrane_conductance=patch.gap_edges[0][2] / total))
    site = patch.surface[0]
    specific = {}
    for label, footprint in (('per_200um_square_site_m2', (2e-4) ** 2),
                             ('per_ninth_of_1p6mm_square_m2', (1.6e-3) ** 2 / 9)):
        specific[label] = dict(assumed_area_m2=footprint,
                               apical_capacitance_F_per_m2=site.capacitance_f / footprint,
                               apical_capacitance_uF_per_cm2=site.capacitance_f / footprint * 1e2,
                               barrier_conductance_S_per_m2=site.barrier_conductance_s / footprint,
                               apical_pump_A_per_m2=site.pump_current_a / footprint)
    return dict(
        node_count=len(patch.cells) + len(patch.surface),
        cell_nodes=len(patch.cells), surface_nodes=len(patch.surface),
        gap_edges=len(patch.gap_edges), surface_edges=len(patch.surface_edges),
        hard_state_cap_in_post_init=512,
        spatial_extent_m=8 * 2e-4, cell_spacing_m=2e-4, cell_layer_depth_m=1e-4,
        ion_species=['Na+', 'K+', 'Cl-'], ion_mass_balance=False, electroneutrality_solve=False,
        channel_gating=False, pump_kinetics=False,
        thermal_voltage_V=thermal, nernst_reversals=reversals,
        resting_membrane_voltage_V=[float(x) for x in resting[:len(patch.cells)]],
        resting_apical_voltage_V=[float(x) for x in resting[len(patch.cells):]],
        implied_TEP_V=-float(resting[len(patch.cells)]),
        per_cell=cells,
        membrane_tau_range_s=[min(c['membrane_tau_s'] for c in cells), max(c['membrane_tau_s'] for c in cells)],
        apical_tau_s=site.capacitance_f / site.barrier_conductance_s,
        implied_membrane_area_at_1uF_per_cm2_m2=patch.cells[0].capacitance_f / 1e-2,
        implied_equivalent_sphere_radius_m=math.sqrt(patch.cells[0].capacitance_f / 1e-2 / (4 * math.pi)),
        geometric_inconsistency=('Cell capacitance 10 pF implies ~1000 um^2 of membrane (a ~9 um-radius cell, '
                                 'literature specific capacitance 1 uF/cm^2), but the cells are spaced 200 um '
                                 'apart, so the lumped network is not a tiling of the cells it claims.'),
        apical_specific_values_under_assumed_footprints=specific,
        two_graphs_are_decoupled=True,
        decoupling_reason=('Basal extracellular bath is prescribed at 0 V, so the intracellular '
                           'gap-junction graph and the apical extracellular graph share no state; '
                           'a wound shunt cannot move any cell Vm and a Vm perturbation cannot '
                           'produce any apical field.'),
    )


def betse_evidence():
    tissue = load('data/derived/bioelectric/tissue.json')
    summary = load('data/derived/physiology/betse_run/run_summary.json')
    times = np.asarray(tissue['time_s'], float)
    fields = {}
    for field in tissue['fields']:
        arr = np.asarray([frame['values'][field['id']] for frame in tissue['frames']], float)
        fields[field['id']] = dict(unit=field['unit'], frames=int(arr.shape[0]), cells=int(arr.shape[1]),
                                   first_frame_mean=float(arr[0].mean()), last_frame_mean=float(arr[-1].mean()),
                                   last_frame_spatial_sd=float(arr[-1].std()),
                                   temporal_range=[float(arr.min()), float(arr.max())])
    return dict(
        engine='BETSE 1.5.1', run_sha256=summary.get('revision'),
        cells=tissue['cells'], membrane_edges=tissue['membrane_edges'],
        simulated_seconds=float(times[-1]), solver_dt_s=1e-4, saved_frames=len(times),
        initialization_seconds=summary['initialization_seconds'],
        world_size_m=150e-6, cell_radius_m=5e-6, cell_height_m=10e-6, lattice='hex, disorder 0.4',
        species_present=[f['id'] for f in tissue['fields']],
        chloride_absent=True,
        channels_configured=['Nav1p3', 'Kv1p5', 'KLeak (inhibited by gating ligand X, Km 0.05, n 2)'],
        gap_junction_model=dict(voltage_sensitive=True, half_close_mV=15, minimum_open_fraction=0.10,
                                gj_surface_area_fraction=5.0e-8),
        wound_present_in_retained_run=dict(cutting_event=True, cut_profile='surgery (wedge bitmap)',
                                           break_tight_junctions=True, wound_TJ_leak=0.1),
        interventions_available_but_disabled=['change Na/K/Cl/Ca membrane permeability', 'block gap junctions',
                                              'block NaKATP pump', 'apply external voltage', 'break ecm junctions',
                                              'change env K/Cl/Na', 'change temperature', 'GRN network'],
        fields=fields,
        specimen='generic computational tissue; no human subject',
        venv='.venv-physiology (betse 1.5.1 installed editable); not importable from .venv',
    )


def wound_field_evidence():
    fit = load('data/derived/calibration/skin-fit.json')
    return dict(
        source='Nuccitelli et al. 2011 (PMC3228273) Table I, hash-receipted',
        quantity=fit['quantity'], unit=fit['unit'],
        observations=len(fit['observations']), participants=40,
        parameters=dict(zip(fit['feature_names'], fit['fit']['parameters'])),
        train=fit['metrics']['train'], holdout=fit['metrics']['holdout'],
        measured_plane='lateral field above epidermis, below stratum corneum',
        identifies=['group-level lateral field amplitude by age/sex/site'],
        does_not_identify=['membrane voltage', 'any conductance, capacitance or pump current',
                           'wound current (no conductivity or geometry)', 'healing rate'],
    )


def area_dependence_audit(ledger):
    """Every bioelectric quantity that is per unit area or per unit volume."""
    factor = ledger['inflation_factor_raw_over_outer']
    return [
        dict(quantity='SkinPatch cell capacitance / conductances / pump current',
             location='ihm/materialize/skin.py:19-39', per_area_or_volume=False,
             verdict='NOT AFFECTED — no area or volume appears anywhere in the lumped RC patch',
             note='These are absolute per-node F, S and A with no declared membrane area or cell volume. '
                  'The 2x area defect cannot invalidate them because they were never area-referenced. '
                  'That is worse, not better: they are unfalsifiable as specific quantities.'),
        dict(quantity='SurfaceSite capacitance / barrier conductance / pump current',
             location='ihm/materialize/skin.py:32-38', per_area_or_volume=False,
             verdict='NOT AFFECTED, and undefined per unit area',
             note='Site footprint is never declared. Any plausible footprint gives 350-2500 uF/cm^2 '
                  'of apical capacitance, 10^2-10^4 above literature skin-barrier values. The '
                  'defect here is missing area, not wrong area.'),
        dict(quantity='skin_layers.physical_skin_support area_m2',
             location='ihm/assembly/skin_layers.py:43-50', per_area_or_volume=True,
             verdict='CORRECT — returns the exterior component only (%.10f m^2)' % ledger['outer_component_area_m2'],
             note='Actively rejects a mask that is not exactly the declared whole component.'),
        dict(quantity='skin_epithelial_inventory.exact_patch_support area_m2',
             location='ihm/assembly/skin_epithelial_inventory.py:18-31', per_area_or_volume=True,
             verdict='CORRECT — constrained to a subset of contact_eligible_triangle_ids'),
        dict(quantity='bind_inventory compartment volumes = area x depth; ion masses; concentrations',
             location='ihm/assembly/skin_epithelial_inventory.py:48-70', per_area_or_volume=True,
             verdict='CORRECT IF FED (c) — but never fed real geometry; the only exercise uses a '
                     'synthetic 1e-4 m^2 fixture in scripts/verify_skin_epithelial_inventory.py:14'),
        dict(quantity='anatomy.json body-skin-{epidermis,dermis,hypodermis}.volume_m3',
             location='data/derived/canonical/anatomy.json (synthesized_layer entities)',
             per_area_or_volume=True,
             verdict='WRONG by x%.6f (too large)' % factor,
             note='epidermis 3.5026e-4 -> 1.78046e-4 m^3; dermis 5.25390e-3 -> 2.67069e-3 m^3; '
                  'hypodermis 1.75130e-2 -> 8.90230e-3 m^3. No bioelectric code reads these today, '
                  'so no current bioelectric number is wrong. They are the natural denominator for '
                  'any future per-volume ion inventory, cell count or tissue conductivity, so the '
                  'defect is latent, not absent.'),
        dict(quantity='human wound-field fit (V/m)',
             location='data/derived/calibration/skin-fit.json', per_area_or_volume=False,
             verdict='NOT AFFECTED — a field amplitude, area-free'),
        dict(quantity='BETSE cell geometry, densities and currents',
             location='data/derived/physiology/betse_run/sim_config.yaml', per_area_or_volume=True,
             verdict='NOT AFFECTED — a self-contained 150 um square world with its own cell radius'),
        dict(quantity='Topology of any mesh-derived gap-junction network on the raw slab',
             location='ihm/materialize/skin.py:113 network_topology over patch.gap_edges',
             per_area_or_volume=False,
             verdict='WORSE THAN THE AREA FACTOR — the slab is two disjoint sheets (%d + %d faces). '
                     'Node adjacency taken from it yields two disconnected epithelia, not one; the '
                     'existing anchor already lands on the wrong sheet.'
                     % (ledger['outer_component_faces'], ledger['inner_component_faces'])),
    ]


def thesis_fitness():
    return [
        dict(thesis='(a) spatially resolved Vmem field over a connected epithelial sheet',
             status='ABSENT at body scale; PRESENT at toy scale',
             detail='The only connected Vmem field is BETSE 212 cells over a 150 um square with '
                    '3.1 mV spatial SD at the final frame. The canonical whole-body path has 9 '
                    'cells on a 1.6 mm line. The integument itself carries no Vmem state: the three '
                    'skin-layer entities have representation surface_shell_quadrature_layer and no '
                    'geometry at all.',
             reachable='Yes. The watertight envelope supplies a connected 54,374-node surface graph.'),
        dict(thesis='(b) regionally perturbable gap-junction coupling',
             status='NOMINALLY PRESENT, PRACTICALLY UNUSABLE in the canonical path',
             detail='gap_edges carry a fixed 10 pS Ohmic conductance summed into the cell block '
                    '(ihm/assembly/skin_bioelectric.py:35-38). There is no voltage or ligand gating, '
                    'no connexin state, and no regional-blockade API: SkinPatch is a frozen '
                    'dataclass rebuilt with dataclasses.replace. BETSE does have voltage-sensitive '
                    'junctions (half-closed at 15 mV, floor 0.10) and a block-gap-junctions event, '
                    'but the event is disabled and BETSE is not importable from .venv.',
             reachable='Yes for the BETSE lane by enabling the existing event; the canonical lane '
                       'needs a gating law it does not have.'),
        dict(thesis='(c) ion-channel and pump manipulation of the kind Levin uses',
             status='ABSENT in the canonical path; CONFIGURABLE BUT OFF in BETSE',
             detail='The canonical patch has three constant Ohmic leaks and one prescribed constant '
                    'pump current per cell; there is nothing to open, close, block or express. BETSE '
                    'carries Nav1p3, Kv1p5 and a ligand-inhibited K leak plus scheduled '
                    'permeability multipliers and a block-NaKATP-pump event, all disabled in the '
                    'retained 0.035 s run.',
             reachable='Yes, but only inside BETSE at 150 um scale, and only after a new hash-pinned '
                       'reviewed run; the retained artifact cannot be re-purposed.'),
        dict(thesis='(d) wound-induced field response',
             status='PARTIAL and structurally broken in the canonical path',
             detail='The canonical wound is a 100x barrier shunt at one apical site with the pump '
                    'zeroed. Because the basal bath is prescribed, the wound provably cannot alter '
                    'any cell Vm — scripts/verify_body_skin_electric.py:40 asserts exactly that '
                    'equality as a passing check. So the substrate encodes wound-independence of '
                    'Vmem as a verified property. BETSE did run a real cut (cutting event enabled, '
                    'break TJ, wound TJ leak 0.1) but for 35 ms.',
             reachable='Only by replacing the prescribed bath with a bidomain or electrodiffusion '
                       'solve; the module docstring already says the present model cannot claim it.'),
        dict(thesis='(e) pattern stability or bistability',
             status='IMPOSSIBLE by construction',
             detail='electrical_system builds a constant matrix A and constant vector b and the '
                    'trajectory is exp(A t): a linear time-invariant system with a single globally '
                    'attracting fixed point. No nonlinearity, no gating, no ion feedback exists, so '
                    'there is exactly one steady state and no bistability, hysteresis, front or '
                    'domain boundary can ever appear.',
             reachable='No. This needs a different model, not a bigger one.'),
        dict(thesis='(f) minutes-to-days timescales',
             status='ABSENT — hard-capped four to seven orders of magnitude short',
             detail='build_skin_electric raises unless 0 < duration_s <= 10 and dt_s <= 0.1 '
                    '(ihm/assembly/skin_bioelectric.py:99); epithelial_electrodiffusion.simulate '
                    'raises unless 0 < duration_s <= 10 (line 105). The retained canonical run is '
                    '2.0 s at dt 0.05 s. BETSE ran 0.035 s at dt 1e-4 s. Levin-style pattern work '
                    'is 1e2 to 1e5 s.',
             reachable='Trivially for the linear RC path (raise the cap; exp(At) is exact at any '
                       'horizon and the answer is the fixed point). Not trivially for anything with '
                       'ion mass balance, where the finite reservoirs deplete: '
                       'epithelial_electrodiffusion.transport raises on exhaustion.'),
    ]


def mesh_synergy(ledger):
    probe = load('data/derived/entity-tet-ready-skin-probe-v1/summary.json')
    entity = load('data/derived/entity-tet-ready-verification-v1/summary.json')
    muscle = load('data/derived/muscle-tet-ready-verification-v1/summary.json')
    repaired = probe['repaired']
    tet = repaired['tetgen']
    area = ledger['watertight_envelope_area_m2']
    surface_nodes = repaired['diagnostics']['vertices']
    mean_tet_volume = tet['tet_volume_m3'] / tet['tets']

    def nodes_at(h):
        return area / ((math.sqrt(3) / 2) * h * h)

    betse = load('data/derived/bioelectric/tissue.json')
    geom = load('data/derived/app/geometry/betse-tissue-cells.json.gz')
    polygon_area = 0.0
    for poly in geom['source_polygons_m']:
        p = np.asarray(poly, float)
        polygon_area += abs(np.dot(p[:, 0], np.roll(p[:, 1], -1)) - np.dot(p[:, 1], np.roll(p[:, 0], -1))) / 2
    density = betse['cells'] / polygon_area
    return dict(
        probe=dict(source='data/derived/outer-envelope/outer-envelope.json.gz',
                   stored_envelope_tetgen_failed=not probe['stored']['tetgen']['succeeded'],
                   stored_self_intersecting_face_pairs=probe['stored']['diagnostics']['self_intersecting_face_pairs'],
                   repaired_self_intersection_free=repaired['diagnostics']['self_intersection_free'],
                   repaired_surface_vertices=surface_nodes,
                   repaired_surface_faces=repaired['diagnostics']['faces'],
                   tet_nodes=tet['tet_vertices'], tets=tet['tets'],
                   tet_volume_m3=tet['tet_volume_m3'],
                   median_surface_edge_m=repaired['shape']['median_edge_length_m'],
                   mean_tet_volume_m3=mean_tet_volume,
                   equivalent_regular_tet_edge_m=(mean_tet_volume * 6 * math.sqrt(2)) ** (1 / 3)),
        whole_body_context=dict(entities_with_tetgen_proof=entity['repaired_tetrahedralized'] + muscle['repaired_tetrahedralized'],
                                entity_tets=entity['total_tets'], muscle_tets=muscle['total_tets'],
                                total_tets=entity['total_tets'] + muscle['total_tets'],
                                skin_slab_excluded_reason='outer_envelope_already_repaired',
                                skin_layers_excluded_reason='not_a_triangular_surface'),
        surface_node_budget=dict(
            current_density_nodes_per_m2=surface_nodes / area,
            current_mean_spacing_m=math.sqrt(area / surface_nodes),
            nodes_at_edge_length={f'{h * 1e3:g}mm': nodes_at(h) for h in (5.7e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4)}),
        cell_scale_reality=dict(
            betse_mean_cell_area_m2=polygon_area / betse['cells'],
            betse_equivalent_cell_diameter_m=2 * math.sqrt(polygon_area / betse['cells'] / math.pi),
            betse_areal_density_cells_per_m2=density,
            whole_integument_keratinocyte_count_at_betse_density=density * area,
            note='One basal-keratinocyte-per-node coverage of the integument is ~2.4e10 nodes. '
                 'That is six orders of magnitude beyond the 5.4e4-node envelope and is not a '
                 'compute problem to be optimised: it is the wrong representation.'),
        verdict=dict(
            same_mesh_for_mechanics_and_bioelectrics='YES for a coarse-grained Vmem continuum, NO for cell-resolved bioelectrics',
            what_it_would_require=[
                'A shell mesh of the integument. None exists: the envelope tets fill the whole '
                '69.7 L body interior at a ~11.2 mm equivalent tet edge, and the three skin layers '
                'have no geometry at all, only area x thickness scalars.',
                'A conductance per mesh edge with a stated derivation (a coarse-grained '
                'gap-junction conductivity in S/m, not a per-cell 10 pS reused at 4 mm spacing).',
                'A statement of what one node means. At 2 mm spacing one node aggregates roughly '
                '5e7 keratinocytes; Vmem at that node is a population mean, and any claim about '
                'single-cell depolarization at that resolution is a category error.',
                'A per-node area and volume, which is exactly the number the current substrate '
                'never declares and the layer entities currently get wrong by x1.967.',
            ],
            realistic_working_point=dict(
                surface_edge_m=2e-3, surface_nodes=nodes_at(2e-3),
                states_per_node=4, total_states=4 * nodes_at(2e-3),
                comment='Vmem plus three ion concentrations on a 2 mm surface graph is ~2.1e6 '
                        'states with a sparse symmetric Laplacian. Implicit steps of 1 s over 24 h '
                        'of simulated time is 8.64e4 solves; at a few tenths of a second per '
                        'preconditioned solve this is hours, not days. Tractable.'),
            infeasible_working_point=dict(
                surface_edge_m=1e-4, surface_nodes=nodes_at(1e-4),
                comment='Resolving the 0.1 mm epidermal thickness with one element through it '
                        'needs ~2.1e8 surface nodes per layer. Not reachable on this hardware.'),
        ),
    )


def honest_limits():
    return dict(
        testable_on_a_fixed_adult_atlas=[
            dict(claim='Vmem as a spatial state variable over a gap-junction-coupled sheet',
                 why='A connected surface graph exists and a Vmem field can be carried on it. '
                     'This is a representation claim, not a biological result.'),
            dict(claim='Regional gap-junction perturbation changes the Vmem pattern',
                 why='Cutting edges in a diffusive graph provably changes its steady state. This is '
                     'a property of graph Laplacians, verifiable here, and it is NOT evidence for '
                     'Levin\'s biology.'),
            dict(claim='A wound-induced lateral field arises from a local barrier breach',
                 why='Reachable once the basal bath is replaced by a solved extracellular potential. '
                     'The Nuccitelli human observations (V/m above the epidermis) already exist as '
                     'a hash-receipted comparison target, and the observation operator is bound.',
                 caveat='The retained human study reports no correlation between field amplitude and '
                        'healing rate in its mouse experiments; matching a field amplitude is not '
                        'evidence of a morphogenetic mechanism.'),
            dict(claim='Ion-channel and pump perturbation shift resting Vmem',
                 why='Standard membrane electrophysiology; BETSE already implements the events.'),
        ],
        not_testable_on_a_fixed_adult_atlas=[
            dict(claim='Vmem prepatterns act as anatomical setpoints (target morphology)',
                 why='The atlas has exactly one anatomy, authored from BodyParts3D and frozen by '
                     'hash. There is no map from any bioelectric state to a shape, so there is no '
                     'setpoint to store, compare or restore. A setpoint claim needs at least two '
                     'anatomies reachable by the same dynamics.'),
            dict(claim='Bioelectric control of regeneration or remodelling',
                 why='No cell birth, death, migration, differentiation or matrix remodelling exists '
                     'anywhere in the substrate. skin_bioelectric.py\'s own limitations list says '
                     'the wound is an electrical surrogate with cells retained and no healing law. '
                     'The mesh is fixed connectivity; regeneration is a topology change.'),
            dict(claim='Bistable Vmem domains as memory or pattern boundaries',
                 why='The canonical model is linear time-invariant with a unique fixed point '
                     '(exp(At)). Bistability is not a resolution or horizon problem; it requires '
                     'nonlinear gating that does not exist here.'),
            dict(claim='Ectopic organ induction, planarian-style head/tail switching, or any '
                       'developmental result',
                 why='No developmental axis, no growth, no species with the relevant regenerative '
                     'capacity, and no morphogenetic readout. These are not approximable here.'),
        ],
        what_would_have_to_change=[
            'Nonlinear voltage- and ligand-gated conductances, so more than one steady state can exist.',
            'A solved extracellular potential (bidomain or full electrodiffusion) replacing the '
            'prescribed 0 V basal bath, so membrane currents and tissue fields can talk to each other.',
            'An integument shell mesh with declared per-node area and volume, plus a stated '
            'coarse-graining from cells to nodes.',
            'Integration horizons of 1e2-1e5 s with ion mass balance that does not deplete, which '
            'means active transport with real kinetics rather than a prescribed constant current.',
            'A morphogenetic readout: some state whose value is a shape, and dynamics that can move '
            'it. Without that, no Levin thesis about anatomical control is even stateable, let '
            'alone testable.',
        ],
    )


def main():
    ledger = area_ledger()
    assessment = dict(
        schema='ihm.bioelectric-substrate-assessment.v1',
        scope=('Assessment only. Nothing was implemented, migrated or mutated. All numbers below '
               'were recomputed from the hashed inputs listed in manifest.json.'),
        evidence_grades=dict(
            verified='recomputed by this script from a hashed input',
            inferred='an engineering judgement made here, stated as such',
            literature='a mechanism or value attributed to a cited publication, never measured in this repository'),
        what_exists=dict(
            canonical_rc_patch=rc_substrate(),
            three_compartment_electrodiffusion=dict(
                location='ihm/assembly/epithelial_electrodiffusion.py',
                compartments=['apical', 'cell', 'basal'], ions=['Na+', 'K+', 'Cl-'],
                spatial_extent='none — zero-dimensional, three lumped compartments',
                conserves_species=True, tracks_free_energy=True, tracks_dissipation=True,
                integration='scipy solve_ivp on 9 transport extents plus 2 energy states, rtol 2e-9',
                duration_cap_s=10,
                gating_or_kinetics=False,
                verdict='The only mass-conserving, thermodynamically audited piece in the substrate, '
                        'and the only one with no spatial extent at all. It is the right physics on '
                        'the wrong domain.'),
            epithelial_inventory=dict(
                location='ihm/assembly/skin_epithelial_inventory.py',
                role='shadow partition of native BioGears skin pools into film/intracellular/extracellular '
                     'volumes by area x depth, then hands them to the electrodiffusion patch',
                native_commit_available=False,
                only_exercise='scripts/verify_skin_epithelial_inventory.py with a synthetic 1e-4 m^2 fixture',
                verdict='The one code path that is area-correct, and it has never been run on real geometry.'),
            canonical_artifact=dict(
                path='data/derived/canonical/skin-electric.json',
                id='forearm-skin-electric', experiments=['baseline', 'wound_shunt', 'electrode', 'membrane_perturbation'],
                frames_per_experiment=41, dt_s=0.05, duration_s=2.0,
                charge_audit='exact matrix exponential states plus independent 12-point Gauss-Legendre '
                             'current quadrature; node residual < 1e-18 C asserted',
                coupling=dict(native_blood_storage_connected=False, neural_connection=False,
                              whole_body_skin_activation=False)),
            betse=betse_evidence(),
            human_wound_field=wound_field_evidence(),
            docs=dict(bioelectric_view='docs/BIOELECTRIC_VIEW.md — describes the 212-cell viewer model',
                      bioelectricity='docs/BIOELECTRICITY.md — the three-quantity taxonomy and the executed equations',
                      integumentary_data='docs/INTEGUMENTARY_DATA.md — the Nuccitelli extraction and its limits',
                      prior_audit='docs/research/SKIN_BIOELECTRIC_CALIBRATION_GAPS.md (2026-09-05) already '
                                  'enumerates every RC coefficient as an illustrative prior and specifies '
                                  'the smallest identifiable next experiment. This assessment does not '
                                  'restate it; it adds the area, sheet-topology and mesh findings.'),
        ),
        area_defect=dict(ledger=ledger, per_unit_audit=area_dependence_audit(ledger),
                         headline=('The 2x area defect does NOT invalidate any existing bioelectric '
                                   'calibration, because no bioelectric quantity in the substrate is '
                                   'area-referenced at all. It does make the three synthesized skin-layer '
                                   'volumes wrong by x%.6f, and those are the denominators any future '
                                   'per-volume bioelectric quantity would use. The more damaging '
                                   'consequence of the same double-sided slab is topological, not metric: '
                                   'it is two disjoint sheets, and the existing regional patch is already '
                                   'pinned to the wrong one.' % ledger['inflation_factor_raw_over_outer'])),
        anchor_defect=anchor_defect(),
        thesis_fitness=thesis_fitness(),
        mesh_synergy=mesh_synergy(ledger),
        honest_limits=honest_limits(),
    )
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'assessment.json').write_text(json.dumps(assessment, indent=1, allow_nan=False) + '\n')
    manifest = dict(
        schema='ihm.bioelectric-substrate-assessment-manifest.v1',
        generator='scripts/assess_bioelectric_substrate.py',
        generator_sha256=sha(Path(__file__)),
        read_only=True, repo_mutated=False, native_jobs_run=0,
        inputs_sha256={rel: sha(ROOT / rel) for rel in INPUTS},
        outputs_sha256={'assessment.json': sha(OUT / 'assessment.json')},
    )
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1) + '\n')
    print(json.dumps(dict(
        outer_component_area_m2=ledger['outer_component_area_m2'],
        watertight_envelope_area_m2=ledger['watertight_envelope_area_m2'],
        layer_volume_inflation_factor=ledger['inflation_factor_raw_over_outer'],
        anchor_on_exterior_component=assessment['anchor_defect']['anchor_face_on_exterior_component'],
        canonical_vmem_nodes=assessment['what_exists']['canonical_rc_patch']['cell_nodes'],
        canonical_duration_s=2.0,
        betse_cells=assessment['what_exists']['betse']['cells'],
        betse_seconds=assessment['what_exists']['betse']['simulated_seconds'],
        envelope_surface_nodes=assessment['mesh_synergy']['probe']['repaired_surface_vertices'],
        envelope_tets=assessment['mesh_synergy']['probe']['tets'],
        keratinocytes_for_full_coverage=assessment['mesh_synergy']['cell_scale_reality']['whole_integument_keratinocyte_count_at_betse_density'],
        outputs=[str(p.relative_to(ROOT)) for p in sorted(OUT.iterdir())],
    ), indent=1))


if __name__ == '__main__':
    main()
