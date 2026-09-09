"""Isotropic remeshing of a garment surface, with its boundary held.

The acquired garment meshes are authored for a different purpose than fitting.
`cargo-trousers` covers two legs in 392 triangles with a median edge of 76 mm,
`long-skirt` carries a single 307 mm edge, and nearly every garment holds
triangles whose smallest corner is under a tenth of a degree. Two consequences
follow, and both are defects the fitter cannot repair on its own:

  * a triangle wider than the body feature it spans passes through that feature
    while all three of its corners stand clear, so the garment penetrates the
    body in the space between its vertices; and
  * a shrinkwrap has only the vertices it is given, so pulling a coarse garment
    onto the body stretches single edges by up to fifty times rather than
    distributing the deformation.

So the surface is rebuilt at a uniform target edge length before it is ever
fitted, by the standard split/collapse/flip/relax iteration (Botsch and
Kobbelt, "A Remeshing Approach to Multiresolution Modeling", 2004). The
rebuild is resampling, not redesign: every vertex is projected back onto the
original surface each pass, and the boundary — a hem, a neckline, an armhole —
is resampled strictly along its own polyline, so the garment's outline is the
one that was authored.

Topology is numpy and scipy only, so this module imports in both the main venv
and the out-of-process libigl venv. Projection onto the original surface is the
caller's, passed in as `project`, because that is the one step that wants an
AABB tree.
"""
from __future__ import annotations

import numpy as np

from ihm.assembly.garment_wardrobe import boundary_edge_count, compact, edge_table, face_components

# Botsch and Kobbelt's thresholds. Splitting above 4/3 of the target and
# collapsing below 4/5 is the pair that does not fight itself: an edge created
# by a split is never immediately a candidate for collapse, and vice versa.
SPLIT_RATIO = 4.0 / 3.0
COLLAPSE_RATIO = 4.0 / 5.0
DEFAULT_PASSES = 10
RELAX_WEIGHT = 0.6
# A collapse that turns a face over is a fold, not a simplification. Faces are
# small and nearly coplanar with their neighbours here, so requiring the normal
# to stay within a right angle of where it was rejects folds without rejecting
# the ordinary movement a collapse causes.
MIN_NORMAL_AGREEMENT = 0.0
DEGENERATE_AREA_M2 = 1e-14


def _edges_with_faces(triangles):
    """Unique undirected edges, and for each the faces carrying it."""
    corner = np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]])
    owner = np.tile(np.arange(len(triangles)), 3)
    key = np.sort(corner, axis=1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key, owner = key[order], owner[order]
    first = np.ones(len(key), bool)
    first[1:] = (key[1:] != key[:-1]).any(1)
    index = np.cumsum(first) - 1
    edges = key[first]
    counts = np.bincount(index, minlength=len(edges))
    faces = [[] for _ in range(len(edges))]
    for slot, face in zip(index, owner):
        faces[slot].append(int(face))
    return edges, faces, counts


def boundary_vertices(triangles):
    edges = np.sort(np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]), axis=1)
    unique, count = np.unique(edges, axis=0, return_counts=True)
    border = unique[count == 1]
    return set(border.ravel().tolist()), {tuple(e) for e in border.tolist()}


class Surface:
    """A triangle soup that can be split, collapsed and flipped in place.

    Deleted faces are tombstoned rather than removed, because renumbering after
    every operator is what makes a remesher slow and easy to get wrong; `compact`
    is applied once at the end.
    """

    def __init__(self, positions, triangles):
        self.v = [np.asarray(p, float) for p in positions]
        self.f = [list(map(int, t)) for t in triangles]
        self.alive = [True] * len(self.f)
        self.adjacent = [set() for _ in self.v]
        for index, face in enumerate(self.f):
            for corner in face:
                self.adjacent[corner].add(index)

    # ------------------------------------------------------------- accessors
    def faces(self):
        return np.array([f for f, live in zip(self.f, self.alive) if live], np.int64)

    def arrays(self):
        return compact(np.array(self.v, float), self.faces())

    def _ring(self, vertex):
        return {c for face in self.adjacent[vertex] if self.alive[face] for c in self.f[face]} - {vertex}

    def _normal(self, face):
        a, b, c = (self.v[i] for i in self.f[face])
        return np.cross(b - a, c - a)

    def _add_vertex(self, point):
        self.v.append(np.asarray(point, float))
        self.adjacent.append(set())
        return len(self.v) - 1

    def _add_face(self, a, b, c):
        if a == b or b == c or a == c:
            return -1
        self.f.append([a, b, c])
        self.alive.append(True)
        index = len(self.f) - 1
        for corner in (a, b, c):
            self.adjacent[corner].add(index)
        return index

    def _kill_face(self, face):
        self.alive[face] = False
        for corner in self.f[face]:
            self.adjacent[corner].discard(face)

    def _faces_on(self, u, w):
        return [face for face in self.adjacent[u] & self.adjacent[w] if self.alive[face]]

    # -------------------------------------------------------------- operators
    def split_long(self, high):
        """Insert a midpoint on every edge above `high`, then retriangulate each
        touched face by how many of its edges were split. Doing the whole pass at
        once keeps the mesh manifold at every step: a face and its neighbour
        always agree about the midpoint on the edge they share."""
        faces = self.faces()
        if not len(faces):
            return 0
        edges = edge_table(faces)
        length = np.linalg.norm(np.array(self.v)[edges[:, 0]] - np.array(self.v)[edges[:, 1]], axis=1)
        long_edges = edges[length > high]
        if not len(long_edges):
            return 0
        midpoint = {}
        for u, w in long_edges.tolist():
            key = (min(u, w), max(u, w))
            midpoint[key] = self._add_vertex(0.5 * (self.v[u] + self.v[w]))
        touched = {face for key in midpoint for face in self._faces_on(*key)}
        for face in touched:
            a, b, c = self.f[face]
            cut = [midpoint.get((min(x, y), max(x, y)), -1) for x, y in ((a, b), (b, c), (c, a))]
            self._kill_face(face)
            ab, bc, ca = cut
            if ab >= 0 and bc >= 0 and ca >= 0:
                self._add_face(a, ab, ca); self._add_face(ab, b, bc)
                self._add_face(ca, bc, c); self._add_face(ab, bc, ca)
            elif sum(x >= 0 for x in cut) == 2:
                # Cut the two split edges and bridge the quad on its shorter
                # diagonal, so the pair of new faces is as near equilateral as
                # this face allows.
                rot = [(a, b, c, ab, bc, ca), (b, c, a, bc, ca, ab), (c, a, b, ca, ab, bc)]
                for p, q, r, pq, qr, rp in rot:
                    if pq >= 0 and qr >= 0 and rp < 0:
                        self._add_face(q, qr, pq)
                        if np.linalg.norm(self.v[pq] - self.v[r]) < np.linalg.norm(self.v[qr] - self.v[p]):
                            self._add_face(p, pq, r); self._add_face(pq, qr, r)
                        else:
                            self._add_face(p, pq, qr); self._add_face(p, qr, r)
                        break
            else:
                for p, q, r, pq in ((a, b, c, ab), (b, c, a, bc), (c, a, b, ca)):
                    if pq >= 0:
                        self._add_face(p, pq, r); self._add_face(pq, q, r)
                        break
        return len(midpoint)

    def collapse_short(self, low, high, border, border_edges):
        """Collapse edges below `low`, one at a time, refusing any collapse that
        would change the surface rather than simplify it.

        Four refusals matter. The link condition keeps the result manifold. A
        boundary vertex is never pulled into the interior, so a hem stays a hem.
        An edge across the interior joining two boundary vertices is left alone,
        because collapsing it pinches the sheet. And a collapse that would turn
        a neighbouring face over, or stretch some other edge past `high`, buys
        nothing."""
        done = 0
        for u, w in sorted({(min(int(a), int(b)), max(int(a), int(b)))
                            for a, b in edge_table(self.faces()).tolist()},
                           key=lambda e: np.linalg.norm(self.v[e[0]] - self.v[e[1]])):
            if not self.adjacent[u] or not self.adjacent[w]:
                continue
            shared = self._faces_on(u, w)
            if not shared:
                continue
            if np.linalg.norm(self.v[u] - self.v[w]) >= low:
                break
            on_border = (min(u, w), max(u, w)) in border_edges
            if u in border and w in border and not on_border:
                continue
            if u in border and w in border and on_border:
                target = 0.5 * (self.v[u] + self.v[w])
            elif u in border:
                u, w = w, u                      # collapse the interior vertex into the boundary one
                target = self.v[w].copy()
            elif w in border:
                target = self.v[w].copy()
            else:
                target = 0.5 * (self.v[u] + self.v[w])
            ring_u, ring_w = self._ring(u), self._ring(w)
            expected = 1 if on_border else 2
            if len(ring_u & ring_w) != expected:
                continue                          # link condition
            moving = [face for face in self.adjacent[u] if self.alive[face] and face not in shared]
            before = [self._normal(face) for face in moving]
            saved = self.v[u].copy()
            self.v[u] = target
            ok = True
            for face, was in zip(moving, before):
                now = self._normal(face)
                if np.linalg.norm(now) <= DEGENERATE_AREA_M2 or float(was @ now) <= MIN_NORMAL_AGREEMENT:
                    ok = False
                    break
            if ok:
                for other in (ring_u | ring_w) - {u, w}:
                    if np.linalg.norm(target - self.v[other]) > high:
                        ok = False
                        break
            if not ok:
                self.v[u] = saved
                continue
            self.v[w] = target
            for face in shared:
                self._kill_face(face)
            for face in list(self.adjacent[u]):
                if not self.alive[face]:
                    continue
                replaced = [w if c == u else c for c in self.f[face]]
                self._kill_face(face)
                self._add_face(*replaced)
            self.adjacent[u].clear()
            done += 1
        return done

    def flip_to_valence(self, border, border_edges):
        """Flip an interior edge when doing so brings the four vertices closer to
        their ideal valence — six inside, four on the boundary. This is what
        turns the long thin triangles left by splitting and collapsing into
        near-equilateral ones."""
        target = lambda vertex: 4 if vertex in border else 6
        done = 0
        for u, w in edge_table(self.faces()).tolist():
            if (min(u, w), max(u, w)) in border_edges:
                continue
            shared = self._faces_on(u, w)
            if len(shared) != 2:
                continue
            left = next(c for c in self.f[shared[0]] if c not in (u, w))
            right = next(c for c in self.f[shared[1]] if c not in (u, w))
            if right in self._ring(left):
                continue                          # the flipped edge already exists
            valence = {x: len([c for c in self._ring(x)]) for x in (u, w, left, right)}
            now = sum(abs(valence[x] - target(x)) for x in (u, w, left, right))
            after = (abs(valence[u] - 1 - target(u)) + abs(valence[w] - 1 - target(w))
                     + abs(valence[left] + 1 - target(left)) + abs(valence[right] + 1 - target(right)))
            if after >= now:
                continue
            was = [self._normal(face) for face in shared]
            self._kill_face(shared[0]); self._kill_face(shared[1])
            first = self._add_face(left, u, right)
            second = self._add_face(right, w, left)
            made = [f for f in (first, second) if f >= 0]
            good = len(made) == 2 and all(
                np.linalg.norm(self._normal(f)) > DEGENERATE_AREA_M2
                and float(max(w0 @ self._normal(f) for w0 in was)) > 0 for f in made)
            if good:
                done += 1
                continue
            for f in made:
                self._kill_face(f)
            self._add_face(*self.f[shared[0]])
            self._add_face(*self.f[shared[1]])
        return done


def _relax(positions, triangles, border, border_edges, weight=RELAX_WEIGHT):
    """Move each vertex toward the centroid of its neighbours, tangentially.

    An interior vertex drops the component along its own normal, so relaxation
    equalises the sampling without shrinking the surface. A boundary vertex is
    moved only toward the midpoint of its two boundary neighbours, which slides
    it along the hem and never off it."""
    n = len(positions)
    total = np.zeros((n, 3))
    count = np.zeros(n)
    for a, b in edge_table(triangles):
        total[a] += positions[b]; count[a] += 1
        total[b] += positions[a]; count[b] += 1
    normal = np.zeros((n, 3))
    p = positions[triangles]
    face_normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    for corner in range(3):
        np.add.at(normal, triangles[:, corner], face_normal)
    length = np.linalg.norm(normal, axis=1)
    normal = normal / np.where(length[:, None] > 0, length[:, None], 1.0)

    moved = positions.copy()
    interior = np.ones(n, bool)
    interior[list(border)] = False
    live = interior & (count > 0)
    step = np.zeros((n, 3))
    step[live] = total[live] / count[live][:, None] - positions[live]
    step[live] -= normal[live] * np.einsum('ij,ij->i', step[live], normal[live])[:, None]
    moved[live] += weight * step[live]

    along = {v: [] for v in border}
    for a, b in border_edges:
        along[a].append(b); along[b].append(a)
    for vertex, neighbours in along.items():
        if len(neighbours) != 2:
            continue                              # a corner where hems meet: leave it
        midpoint = 0.5 * (positions[neighbours[0]] + positions[neighbours[1]])
        tangent = positions[neighbours[1]] - positions[neighbours[0]]
        norm = np.linalg.norm(tangent)
        if norm <= 0:
            continue
        tangent = tangent / norm
        moved[vertex] = positions[vertex] + weight * tangent * float((midpoint - positions[vertex]) @ tangent)
    return moved


def split_components(positions, triangles):
    """One (positions, triangles) pair per connected surface component.

    Several acquired garments are not one sheet: the polo shirt carries five
    face components, three of them scraps of eight vertices. Resampling them as
    a single mesh lets a scrap be collapsed out of existence or merged into its
    neighbour, so each component is resampled against its own surface and its
    own boundary instead, and can then be checked for having survived.
    """
    _, label = face_components(len(positions), triangles)
    parts = []
    for key in np.unique(label[triangles[:, 0]]):
        keep = triangles[label[triangles[:, 0]] == key]
        # Renumber onto the vertices this component actually uses. Splitting
        # partitions the garment and must not also clean it: the polo shirt
        # carries a component whose faces are all degenerate, and dropping them
        # here would delete the component instead of reporting that it is too
        # small to resample.
        used, inverse = np.unique(keep, return_inverse=True)
        parts.append((positions[used], inverse.reshape(-1, 3).astype(np.int64)))
    return parts


def remesh(positions, triangles, target_edge_m, *, passes=DEFAULT_PASSES, project=None,
           project_boundary=None, report=None):
    """Resample a garment surface at a uniform `target_edge_m`.

    `project` maps points back onto the original surface and `project_boundary`
    maps boundary points back onto the original boundary polyline. Both are the
    caller's, and without them the result is a smoothed surface rather than a
    resampled one, so the build always supplies them.
    """
    positions = np.asarray(positions, float)
    triangles = np.asarray(triangles, np.int64)
    high, low = SPLIT_RATIO * target_edge_m, COLLAPSE_RATIO * target_edge_m
    original_loops = boundary_edge_count(triangles)[0]
    history = []
    surface = Surface(positions, triangles)
    for step in range(passes):
        border, border_edges = boundary_vertices(surface.faces())
        split = surface.split_long(high)
        border, border_edges = boundary_vertices(surface.faces())
        collapsed = surface.collapse_short(low, high, border, border_edges)
        flipped = surface.flip_to_valence(*boundary_vertices(surface.faces()))
        v, f = surface.arrays()
        border, border_edges = boundary_vertices(f)
        v = _relax(v, f, border, border_edges)
        if project is not None:
            inside = np.ones(len(v), bool)
            inside[list(border)] = False
            if inside.any():
                v[inside] = project(v[inside])
            if project_boundary is not None and (~inside).any():
                v[~inside] = project_boundary(v[~inside])
        surface = Surface(v, f)
        history.append({'pass': step + 1, 'split_edges': int(split), 'collapsed_edges': int(collapsed),
                        'flipped_edges': int(flipped), 'vertices': int(len(v)), 'faces': int(len(f))})
    v, f = surface.arrays()
    loops, nonmanifold = boundary_edge_count(f)
    if nonmanifold:
        raise ValueError(f'remesh produced {nonmanifold} non-manifold edges')
    if face_components(len(v), f)[0] != face_components(len(positions), triangles)[0]:
        raise ValueError('remesh changed the number of surface components')
    if report is not None:
        report.update({'target_edge_mm': target_edge_m * 1e3, 'passes': passes, 'history': history,
                       'boundary_edges_before': original_loops, 'boundary_edges_after': loops,
                       'source_vertices': int(len(positions)), 'source_faces': int(len(triangles)),
                       'vertices': int(len(v)), 'faces': int(len(f))})
    return v, f


def quality(positions, triangles):
    """Edge-length and corner-angle statistics: the numbers `crude` meant."""
    edges = edge_table(triangles)
    length = np.linalg.norm(positions[edges[:, 0]] - positions[edges[:, 1]], axis=1)
    p = positions[triangles]
    angles = []
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        u, w = p[:, b] - p[:, a], p[:, c] - p[:, a]
        nu, nw = np.linalg.norm(u, axis=1), np.linalg.norm(w, axis=1)
        cosine = np.einsum('ij,ij->i', u, w) / np.where(nu * nw > 0, nu * nw, 1.0)
        angles.append(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
    smallest = np.min(angles, axis=0)
    return {
        'vertices': int(len(positions)), 'faces': int(len(triangles)), 'edges': int(len(edges)),
        'edge_mm_median': float(np.median(length) * 1e3), 'edge_mm_max': float(length.max() * 1e3),
        'edge_mm_min': float(length.min() * 1e3),
        'edge_mm_percentiles': {str(q): float(np.percentile(length, q) * 1e3) for q in (1, 5, 50, 95, 99)},
        'min_corner_angle_deg': float(smallest.min()),
        'faces_under_10_deg': int((smallest < 10.0).sum()),
        'faces_under_1_deg': int((smallest < 1.0).sum()),
    }
