"""FEBio 4.13 as an independent check on a prescribed-displacement quasistatic solve, and a
reader for fTetWild's MSH output.

The same compressible neo-Hookean law (FEBio "neo-Hookean": W = mu/2 (I1-3) - mu ln J +
lambda/2 (ln J)^2) on the same linear tetrahedra and the same boundary conditions, solved by
Newton iterations in FEBio instead of L-BFGS-B energy minimisation in DeformableRegion. Spec
4.0 syntax, verified on a known answer: two nodes given per-node x displacements of 0.01 and
0.05 by a NodeData map come back as exactly 0.01 and 0.05.
"""
import os, subprocess
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FEBIO = ROOT / "data/derived/mechanics-reference/build/bin/febio4"


def rcm_order(n_nodes, T):
    """reverse Cuthill-McKee node order for the tets' adjacency graph: FEBio's skyline solver stores
    the whole matrix envelope, and an unstructured tet mesh in arbitrary numbering makes that
    envelope enormous (57,453 equations -> 870 million stored entries on a breast)."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import reverse_cuthill_mckee
    i = np.repeat(T, 4, axis=1).ravel(); j = np.tile(T, (1, 4)).ravel()
    A = sp.coo_matrix((np.ones(len(i), np.int8), (i, j)), shape=(n_nodes, n_nodes)).tocsr()
    return np.asarray(reverse_cuthill_mckee(A, symmetric_mode=True))


def write_feb(path, X, T, mask, disp, *, young_pa, nu, steps=10, auto_stepper=False, reorder=False):
    """auto_stepper=False reproduces the validated block runs exactly; True lets FEBio cut a step
    back when Newton fails, for large prescribed displacements (it changes step control, not the
    converged solution)."""
    X = np.asarray(X, float); T = np.asarray(T, int); mask = np.asarray(mask, bool); disp = np.asarray(disp, float)
    order = rcm_order(len(X), T) if reorder else np.arange(len(X))   # file node k is original node order[k]
    rank = np.empty_like(order); rank[order] = np.arange(len(order))
    X, mask, disp, T = X[order], mask[order], disp[order], rank[T]
    np.save(Path(path).with_suffix(".order.npy"), order)
    out = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">', '<Module type="solid"/>',
           f'<Control><analysis>STATIC</analysis><time_steps>{steps}</time_steps><step_size>{1.0/steps}</step_size>'
           + (f'<time_stepper type="default"><max_retries>10</max_retries><opt_iter>10</opt_iter><dtmin>{1e-4/steps}</dtmin>'
              f'<dtmax>{1.0/steps}</dtmax><aggressiveness>0</aggressiveness><cutback>0.5</cutback></time_stepper>' if auto_stepper else '')
           + '<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness></solver></Control>',
           f'<Material><material id="1" name="m" type="neo-Hookean"><density>1</density><E>{young_pa!r}</E><v>{nu!r}</v></material></Material>',
           '<Mesh><Nodes name="all">']
    out += [f'<node id="{i+1}">{x:.12g},{y:.12g},{z:.12g}</node>' for i, (x, y, z) in enumerate(X)]
    out += ['</Nodes><Elements type="tet4" name="p">']
    out += [f'<elem id="{i+1}">{a+1},{b+1},{c+1},{d+1}</elem>' for i, (a, b, c, d) in enumerate(T)]
    out += ['</Elements>']
    sets = {}
    for c in range(3):
        sel = np.flatnonzero(mask[:, c])
        if len(sel): sets[c] = sel; out.append(f'<NodeSet name="p{c}">' + ",".join(str(s + 1) for s in sel) + '</NodeSet>')
    out += ['</Mesh>', '<MeshDomains><SolidDomain name="p" mat="m"/></MeshDomains>', '<MeshData>']
    for c, sel in sets.items():
        out.append(f'<NodeData name="u{c}" node_set="p{c}" data_type="scalar">'
                   + "".join(f'<node lid="{k+1}">{disp[s, c]:.12g}</node>' for k, s in enumerate(sel)) + '</NodeData>')
    out += ['</MeshData>', '<Boundary>']
    for c in sets:
        out.append(f'<bc name="b{c}" node_set="p{c}" type="prescribed displacement"><dof>{"xyz"[c]}</dof>'
                   f'<value lc="1" type="map">u{c}</value><relative>0</relative></bc>')
    out += ['</Boundary>',
            '<LoadData><load_controller id="1" name="lc1" type="loadcurve"><interpolate>LINEAR</interpolate>'
            '<points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>',
            '<Output><logfile><node_data data="ux;uy;uz" file="nodes.txt" delim=","></node_data></logfile></Output>',
            '</febio_spec>']
    Path(path).write_text("\n".join(out) + "\n")


def run_febio(path, *, threads=8, timeout=7200):
    path = Path(path); log = path.with_suffix(".log")
    with open(log, "w") as h:
        proc = subprocess.run([str(FEBIO), "-i", path.name, "-nosplash"], cwd=path.parent, stdout=h,
                              stderr=subprocess.STDOUT, timeout=timeout, env={**os.environ, "OMP_NUM_THREADS": str(threads)})
    text = log.read_text(errors="replace")
    blocks, t = [], None
    for line in (path.parent / "nodes.txt").read_text().splitlines():
        if line.startswith("*Time"): t = float(line.split("=")[1]); blocks.append((t, []))
        elif line and not line.startswith("*") and blocks: blocks[-1][1].append([float(v) for v in line.split(",")])
    if not blocks: raise RuntimeError(f"FEBio wrote no displacements (exit {proc.returncode}); see {log}")
    t_last, rows = blocks[-1]; rows = np.asarray(rows)
    u_file = np.zeros((int(rows[:, 0].max()), 3)); u_file[rows[:, 0].astype(int) - 1] = rows[:, 1:4]
    order_path = path.with_suffix(".order.npy")
    order = np.load(order_path) if order_path.exists() else np.arange(len(u_file))
    u = np.zeros_like(u_file); u[order] = u_file              # back to the caller's node numbering
    normal = "N O R M A L   T E R M I N A T I O N" in text
    return dict(displacement=u, final_time=t_last, normal_termination=normal, returncode=proc.returncode, log=str(log))


def read_msh_tets(path):
    """fTetWild ASCII .msh, format 2.2 or 4.1 -> (vertices (N,3), tets (M,4), 0-based)."""
    lines = Path(path).read_text().split("\n"); i = 0; V = T = None
    version = None
    while i < len(lines):
        s = lines[i].strip()
        if s == "$MeshFormat": version = float(lines[i + 1].split()[0])
        elif s == "$Nodes":
            if version < 4:
                n = int(lines[i + 1]); arr = np.array([l.split() for l in lines[i + 2:i + 2 + n]], float)
                ids = arr[:, 0].astype(int); V = np.zeros((ids.max(), 3)); V[ids - 1] = arr[:, 1:4]; idmap = None
            else:
                nb, nn = map(int, lines[i + 1].split()[:2]); k = i + 2; ids, xyz = [], []
                for _ in range(nb):
                    m = int(lines[k].split()[3]); k += 1
                    ids += [int(x) for x in lines[k:k + m]]; k += m
                    xyz += [list(map(float, l.split()[:3])) for l in lines[k:k + m]]; k += m
                ids = np.array(ids); V = np.zeros((ids.max(), 3)); V[ids - 1] = np.array(xyz)
        elif s == "$Elements":
            tets = []
            if version < 4:
                n = int(lines[i + 1])
                for l in lines[i + 2:i + 2 + n]:
                    p = l.split()
                    if p[1] == "4": nt = int(p[2]); tets.append([int(x) for x in p[3 + nt:7 + nt]])
            else:
                nb = int(lines[i + 1].split()[0]); k = i + 2
                for _ in range(nb):
                    _, _, et, m = map(int, lines[k].split()[:4]); k += 1
                    if et == 4: tets += [[int(x) for x in l.split()[1:5]] for l in lines[k:k + m]]
                    k += m
            T = np.array(tets, int) - 1
        i += 1
    return V, T
