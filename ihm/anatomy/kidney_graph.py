"""Strict ASCII Amira reader and descriptive kidney-subgraph measurements.

Coordinates stay in the file's local frame. Thickness is deliberately not
renamed radius: the published example does not declare that semantic.
"""
from pathlib import Path
import re
import numpy as np


def read_amira(path: str | Path, *, coordinate_unit: str | None = None) -> dict:
    text = Path(path).read_bytes().decode("latin1")
    if not text.startswith("# Avizo 3D ASCII"):
        raise ValueError("Only ASCII Avizo spatial graphs are supported")
    declared = re.search(r'Coordinates\s+"([^"]+)"', text)
    unit = coordinate_unit or (declared.group(1) if declared else None)
    scale = {"µm": 0.001, "um": 0.001, "mm": 1.0}.get(unit)
    if scale is None:
        raise ValueError(f"Unknown coordinate units: {unit!r}")
    counts = {k: int(v) for k, v in re.findall(r"define (VERTEX|EDGE|POINT) (\d+)", text)}
    blocks = re.split(r"(?m)^@(\d+)\s*$", text)
    data = {int(blocks[i]): blocks[i+1] for i in range(1, len(blocks), 2)}
    fields = {}
    for entity, dtype, width, name, block in re.findall(
        r"(VERTEX|EDGE|POINT)\s*\{\s*(float|int)(?:\[(\d+)\])?\s+(\w+)\s*\}\s*@(\d+)", text
    ):
        values = np.fromstring(data[int(block)], sep=" ", dtype=float)
        if dtype == "int":
            if not np.isfinite(values).all() or np.any(values != np.floor(values)):
                raise ValueError(f"Noninteger values in {name}")
            values = values.astype(np.int64)
        width = int(width or 1)
        if values.size != counts[entity] * width:
            raise ValueError(f"Count mismatch in {name}")
        fields[name] = values.reshape(-1, width) if width > 1 else values
    graph = {
        "nodes_mm": fields["VertexCoordinates"] * scale,
        "edges": fields["EdgeConnectivity"],
        "point_offsets": np.r_[0, np.cumsum(fields["NumEdgePoints"])],
        "points_mm": fields["EdgePointCoordinates"] * scale,
        "thickness_native": fields["thickness"],
    }
    for name in ("strahler", "topo"):
        if name in fields:
            graph[name] = fields[name]
    matrix = re.search(r"(?m)^\s*TransformationMatrix ([^,\n]+)", text)
    if matrix:
        graph["source_transform"] = np.fromstring(matrix.group(1), sep=" ").reshape(4, 4)
    validate_graph(graph)
    return graph


def validate_graph(graph: dict) -> None:
    nodes, edges, points = (graph[x] for x in ("nodes_mm", "edges", "points_mm"))
    offsets = graph["point_offsets"]
    if edges.ndim != 2 or edges.shape[1] != 2 or np.any(edges < 0) or np.any(edges >= len(nodes)):
        raise ValueError("Invalid edge connectivity")
    if len(offsets) != len(edges)+1 or offsets[0] != 0 or offsets[-1] != len(points) or np.any(np.diff(offsets) < 2):
        raise ValueError("Invalid polyline offsets")
    if len(graph["thickness_native"]) != len(points):
        raise ValueError("Thickness/point mismatch")
    if any(not np.isfinite(graph[k]).all() for k in ("nodes_mm", "points_mm", "thickness_native")):
        raise ValueError("Nonfinite measurement")
    for i, (a, b) in enumerate(edges):
        ends = points[[offsets[i], offsets[i+1]-1]]
        if not (np.allclose(ends, nodes[[a,b]], atol=1e-5) or np.allclose(ends, nodes[[b,a]], atol=1e-5)):
            raise ValueError(f"Polyline {i} endpoints disagree with nodes")


def describe(values) -> dict:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"count": 0}
    return {"count": len(a), "mean": float(a.mean()), "std_population": float(a.std()),
            "min": float(a.min()), "q05": float(np.quantile(a,.05)),
            "median": float(np.median(a)), "q95": float(np.quantile(a,.95)), "max": float(a.max())}


def graph_statistics(graph: dict) -> dict:
    """Geometry statistics; radius/Murray remain unavailable without semantics."""
    validate_graph(graph)
    offsets, points, edges = (graph[k] for k in ("point_offsets", "points_mm", "edges"))
    lengths = np.array([np.linalg.norm(np.diff(points[a:b],axis=0),axis=1).sum() for a,b in zip(offsets[:-1],offsets[1:])])
    chords = np.linalg.norm(graph["nodes_mm"][edges[:,0]]-graph["nodes_mm"][edges[:,1]],axis=1)
    tortuosity = np.divide(lengths,chords,out=np.full_like(lengths,np.nan),where=chords>0)
    thickness = np.array([graph["thickness_native"][a:b].mean() for a,b in zip(offsets[:-1],offsets[1:])])
    degrees = np.bincount(edges.ravel(),minlength=len(graph["nodes_mm"]))
    metrics = {"segment_length_mm": lengths,"path_over_chord_tortuosity":tortuosity,"mean_point_thickness_native":thickness}
    result = {"node_count":len(degrees),"edge_count":len(edges),"point_count":len(points),
              "node_degree_histogram":{str(int(k)):int((degrees==k).sum()) for k in np.unique(degrees)},
              "global":{k:describe(v) for k,v in metrics.items()},"conditioned":{},
              "radius_statistics":None,"murray_residual_statistics":None,
              "radius_unavailable_reason":"Example thickness field has no verified radius/diameter semantics or physical units; no correction status established."}
    for order in ("strahler","topo"):
        if order in graph:
            result["conditioned"][order] = {str(int(n)):{k:describe(v[graph[order]==n]) for k,v in metrics.items()} for n in np.unique(graph[order])}
    return result
