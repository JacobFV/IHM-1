#!/usr/bin/env python3
"""Resampling and intersection: the two geometry facts the wardrobe now rests on.

Neither needs libigl, so both are checked here in the main venv against surfaces
whose answers are known by construction rather than against the wardrobe run.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ihm.assembly.garment_remesh import quality, remesh, split_components
from ihm.assembly.garment_wardrobe import boundary_edge_count, face_components, intersection_report

ROOT = Path(__file__).resolve().parents[1]


def sheet(divisions, size=1.0):
    """A flat square, and the projectors that hold it flat and square."""
    axis = np.linspace(0.0, size, divisions + 1)
    positions = np.array([[axis[i], axis[j], 0.0] for i in range(divisions + 1) for j in range(divisions + 1)])
    at = lambda i, j: i * (divisions + 1) + j
    triangles = np.array([t for i in range(divisions) for j in range(divisions)
                          for t in ((at(i, j), at(i + 1, j), at(i + 1, j + 1)),
                                    (at(i, j), at(i + 1, j + 1), at(i, j + 1)))], np.int64)

    def onto_plane(points):
        return np.column_stack([np.clip(points[:, 0], 0.0, size), np.clip(points[:, 1], 0.0, size),
                                np.zeros(len(points))])

    def onto_edge(points):
        held = onto_plane(points)
        distance = np.column_stack([held[:, 0], size - held[:, 0], held[:, 1], size - held[:, 1]])
        side = np.argmin(distance, axis=1)
        held[side == 0, 0] = 0.0
        held[side == 1, 0] = size
        held[side == 2, 1] = 0.0
        held[side == 3, 1] = size
        return held

    return positions, triangles, onto_plane, onto_edge


def sleeve(around=40, rings=2, radius=0.05, tall=0.40):
    """An open tube meshed the way the acquired garments are: many segments
    around, almost no rings along, so every triangle is a sliver."""
    theta = np.linspace(0.0, 2 * np.pi, around, endpoint=False)
    positions = np.array([[radius * np.cos(a), tall * r / rings, radius * np.sin(a)]
                          for r in range(rings + 1) for a in theta])
    triangles = np.array(
        [[r * around + j, r * around + (j + 1) % around, (r + 1) * around + (j + 1) % around]
         for r in range(rings) for j in range(around)]
        + [[r * around + j, (r + 1) * around + (j + 1) % around, (r + 1) * around + j]
           for r in range(rings) for j in range(around)], np.int64)

    def onto_wall(points):
        span = np.linalg.norm(points[:, [0, 2]], axis=1)
        span = np.where(span > 0.0, span, 1.0)
        return np.column_stack([radius * points[:, 0] / span, np.clip(points[:, 1], 0.0, tall),
                                radius * points[:, 2] / span])

    def onto_rim(points):
        held = onto_wall(points)
        held[:, 1] = np.where(np.abs(held[:, 1]) < np.abs(held[:, 1] - tall), 0.0, tall)
        return held

    return positions, triangles, onto_wall, onto_rim


def area(positions, triangles):
    p = positions[triangles]
    return float(0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum())


def main():
    report = {}

    # ------------------------------------------------------- intersection test
    flat = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0]])
    spike = np.array([[0.2, 0.2, -1.], [0.2, 0.2, 1.], [0.6, 0.3, 1.]])
    one = np.array([[0, 1, 2]], np.int64)
    two = np.array([[0, 1, 2], [3, 4, 5]], np.int64)
    assert intersection_report(np.vstack([flat, spike]), two)['intersecting_face_pairs'] == 1
    assert intersection_report(np.vstack([flat, spike + [10.0, 0, 0]]), two)['intersecting_face_pairs'] == 0
    assert intersection_report(flat, one, spike, one)['intersecting_face_pairs'] == 1
    assert intersection_report(flat, one, spike + [10.0, 0, 0], one)['intersecting_face_pairs'] == 0
    # Faces meeting along a shared edge touch by construction; reporting them as
    # crossing would make every closed mesh look self-intersecting.
    hinge = np.array([[0., 0, 0], [1., 0, 0], [0., 1, 0], [1., 1, 0]])
    assert intersection_report(hinge, np.array([[0, 1, 2], [1, 3, 2]], np.int64))['intersecting_face_pairs'] == 0
    # A face too thin to carry a trustworthy normal is excluded and counted.
    needle = np.array([[0., 0, 0], [1., 0, 0], [0.5, 1e-13, 0]])
    assert intersection_report(np.vstack([needle, spike]), two)['degenerate_faces_excluded'] >= 1
    # Two tubes one inside the other cross nowhere; pushed together they do.
    inner_v, inner_f, _, _ = sleeve(radius=0.05)
    outer_v, outer_f, _, _ = sleeve(radius=0.07)
    assert intersection_report(inner_v, inner_f, outer_v, outer_f)['intersecting_face_pairs'] == 0
    tilted = outer_v + [0.035, 0.0, 0.0]
    assert intersection_report(inner_v, inner_f, tilted, outer_f)['intersecting_face_pairs'] > 0
    report['intersection'] = 'exact on skewered, disjoint, shared-edge, degenerate and nested-tube cases'

    # ------------------------------------------------------------ flat resample
    positions, triangles, onto_plane, onto_edge = sheet(4)
    receipt = {}
    fine, faces = remesh(positions, triangles, 0.10, passes=8,
                         project=onto_plane, project_boundary=onto_edge, report=receipt)
    fine_quality = quality(fine, faces)
    assert boundary_edge_count(faces)[1] == 0, 'resample must stay manifold'
    assert face_components(len(fine), faces)[0] == 1
    assert receipt['boundary_edges_after'] >= receipt['boundary_edges_before']
    assert np.abs(fine[:, 2]).max() < 1e-12, 'a flat sheet must stay flat'
    assert fine[:, :2].min() > -1e-9 and fine[:, :2].max() < 1.0 + 1e-9, 'the outline must not grow'
    assert abs(area(fine, faces) - area(positions, triangles)) < 1e-9, 'a flat resample must preserve area exactly'
    assert 0.8 * 100.0 <= fine_quality['edge_mm_median'] <= 1.34 * 100.0, fine_quality
    report['flat_sheet'] = {'faces': fine_quality['faces'], 'edge_mm_median': fine_quality['edge_mm_median'],
                            'min_corner_angle_deg': fine_quality['min_corner_angle_deg']}

    # ---------------------------------------------------------- sleeve resample
    positions, triangles, onto_wall, onto_rim = sleeve()
    coarse = quality(positions, triangles)
    assert coarse['min_corner_angle_deg'] < 5.0, coarse
    assert coarse['faces_under_10_deg'] == coarse['faces'], 'the fixture must start as all slivers'
    receipt = {}
    fine, faces = remesh(positions, triangles, 0.012, passes=10,
                         project=onto_wall, project_boundary=onto_rim, report=receipt)
    fine_quality = quality(fine, faces)
    assert boundary_edge_count(faces)[1] == 0
    assert face_components(len(fine), faces)[0] == 1
    assert receipt['boundary_edges_after'] == receipt['boundary_edges_before'], 'both rims must survive'
    # Ten degrees is the sliver line the acquired meshes fail: no face may sit
    # below it, and the worst corner in the whole mesh must clear it too.
    assert fine_quality['faces_under_10_deg'] == 0, fine_quality
    assert fine_quality['min_corner_angle_deg'] > 10.0, fine_quality
    assert 0.8 * 12.0 <= fine_quality['edge_mm_median'] <= 1.34 * 12.0, fine_quality
    radius = np.linalg.norm(fine[:, [0, 2]], axis=1)
    assert np.abs(radius - 0.05).max() < 1e-9, 'resampled vertices must lie on the source surface'
    assert fine[:, 1].min() > -1e-9 and fine[:, 1].max() < 0.40 + 1e-9, 'rims must not migrate along the tube'
    assert abs(area(fine, faces) / area(positions, triangles) - 1.0) < 0.01
    # A rim vertex must still have exactly two rim neighbours: an open cylinder
    # has two closed boundary loops and a resample may not open or pinch either.
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    unique, count = np.unique(edges, axis=0, return_counts=True)
    border = unique[count == 1]
    degree = np.bincount(border.ravel())
    assert set(np.unique(degree[degree > 0]).tolist()) == {2}, 'every rim vertex needs exactly two rim edges'
    report['sleeve'] = {'faces': f"{coarse['faces']} -> {fine_quality['faces']}",
                        'edge_mm_median': f"{coarse['edge_mm_median']:.1f} -> {fine_quality['edge_mm_median']:.1f}",
                        'min_corner_angle_deg': f"{coarse['min_corner_angle_deg']:.2f} -> "
                                                f"{fine_quality['min_corner_angle_deg']:.2f}",
                        'faces_under_10_deg': f"{coarse['faces_under_10_deg']} -> "
                                              f"{fine_quality['faces_under_10_deg']}"}

    # A resample that is handed no projector smooths the surface instead of
    # resampling it, which is why the build always supplies one. Stated here so
    # the difference is a measured fact rather than a claim in a docstring.
    drifted, drifted_faces = remesh(positions, triangles, 0.012, passes=10)
    drift = np.abs(np.linalg.norm(drifted[:, [0, 2]], axis=1) - 0.05).max()
    assert drift > 1e-6, 'without projection the surface must be seen to move'
    report['projection_matters'] = {'unprojected_radius_drift_mm': float(drift * 1e3)}

    # Several garments are not one sheet, and a component small enough to be
    # collapsed out of existence has to survive the split as its own piece.
    sheet_v, sheet_f, _, _ = sheet(4)
    scrap_v = np.array([[9., 9, 9], [9.01, 9, 9], [9., 9.01, 9], [9., 9, 9.01]])
    scrap_f = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]], np.int64)
    mixed_v = np.vstack([sheet_v, scrap_v])
    mixed_f = np.vstack([sheet_f, scrap_f + len(sheet_v)])
    parts = split_components(mixed_v, mixed_f)
    assert len(parts) == 2, f'expected two components, got {len(parts)}'
    assert sorted(len(f) for _, f in parts) == sorted([len(sheet_f), len(scrap_f)])
    for part_v, part_f in parts:
        assert part_f.min() == 0 and part_f.max() == len(part_v) - 1, 'each piece must be renumbered onto its own'
    # Splitting partitions; it must not also clean. A component whose faces are
    # all degenerate has to come back as a component, not vanish.
    flatscrap_v = np.array([[3., 0, 0], [3., 0, 0], [3., 0, 0]])
    flatscrap_f = np.array([[0, 1, 2]], np.int64)
    both = split_components(np.vstack([sheet_v, flatscrap_v]),
                            np.vstack([sheet_f, flatscrap_f + len(sheet_v)]))
    assert len(both) == 2, 'a fully degenerate component must survive the split'
    report['components'] = 'split preserves piece count, renumbering, and degenerate pieces'

    print(json.dumps({'passed': True, **report}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
