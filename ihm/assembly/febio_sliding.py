"""FEBio 4.13 counterpart of the sliding base (docs/BODY_PARAMETERS.md, "The sliding base").

The tissue's base slides without friction on a RIGID bed. The bed is a shell of rigid elements,
fixed in space, with its normals facing the tissue. Contact is FEBio's facet-to-facet
"sliding-elastic" (penalty): on the HELD faces -- base faces whose three nodes all start behind
the bed -- tension is allowed, so they stay on the bed both ways; on the other base faces it is
not, so they may not enter it. FEBio's contact is defined on faces while the in-repo solver's is
on nodes, so a held node that sits in no all-held face is held only one way here; writers report
how many.

The held faces start PENETRATING the bed, up to tens of millimetres on a breast. A penalty
contact at full strength would resolve that in one step, so the penalty is ramped by a load curve
from 1e-3 of its value to full over the first `ramp_until` of the run and then held: the overlap is
pushed out gradually, the same continuation the in-repo solver gets from ramping its held gaps.
THE BED IS ADVANCED, NOT BURIED. A breast starts up to 36 mm inside the bed, and resolving that
overlap in place defeated FEBio: 582 s without finishing one time step. So the bed is written
RETRACTED (written at bedV - retract_m) and given a rigid displacement of +retract_m back to its
true position over the run -- so retract_m points FROM the retracted position TO the true one -- with a
search radius smaller than the retraction so tension cannot grab across the opening gap at t = 0.
Contact then engages progressively -- the counterpart of the in-repo solver's ramped gap closure --
and the final configuration is the true one: bed in place, held faces on it.

The residual gap of a penalty contact scales inversely with the penalty: at FEBio's automatic
penalty (factor 1) the cylinder block's held faces still penetrated a median 0.031 mm and worst
0.326 mm, which alone put the two solvers 10.4% apart against a 5% gate. Augmented Lagrangian
closed that (1.47%) but could not finish: with the penalty ramped from a thousandth, enforcing a
0.01 mm gap took ~141 equilibrium iterations per step and exhausted the step retries. So the
penalty is simply STIFF (default factor 50, ramped from 0.02 of it over the first half of the run,
which still eases the initial penetration out). The gap that remains is reported, not judged.
"""
from pathlib import Path
import numpy as np


def write_sliding_feb(path, X, T, held_faces, free_faces, bedV, bedF, *, young_pa, nu, pins=(),
                      steps=20, penalty=1.0, ramp_until=0.5, search_radius=0.05, reorder=None,
                      augment=False, gaptol_m=1e-5, aug_tol=0.05, maxaug=30, ramp_from=0.02, retract_m=(0., 0., 0.)):
    """pins: (node, axis) fixed at zero displacement. reorder: optional node permutation (file
    node k is original node reorder[k]); results are read back through run_febio's order file."""
    X = np.asarray(X, float); T = np.asarray(T, int); N = len(X); nb = len(bedV); ne = len(T)
    order = np.arange(N) if reorder is None else np.asarray(reorder)
    rank = np.empty_like(order); rank[order] = np.arange(N)
    Xf, Tf = X[order], rank[T]
    hf = rank[np.asarray(held_faces, int).reshape(-1, 3)]; ff = rank[np.asarray(free_faces, int).reshape(-1, 3)]
    np.save(Path(path).with_suffix(".order.npy"), order)
    tri = lambda F, off: "".join(f'<tri3 id="{k+1}">{a+off},{b+off},{c+off}</tri3>' for k, (a, b, c) in enumerate(F))
    cp = lambda tension: (f'<laugon>{int(bool(augment))}</laugon><tolerance>{float(aug_tol):.12g}</tolerance><gaptol>{float(gaptol_m):.12g}</gaptol>'
                          f'<minaug>1</minaug><maxaug>{maxaug}</maxaug>'
                          f'<penalty lc="2">{float(penalty):.12g}</penalty><auto_penalty>1</auto_penalty>'
                          f'<two_pass>0</two_pass><search_tol>0.01</search_tol><search_radius>{float(search_radius):.12g}</search_radius>'
                          f'<symmetric_stiffness>1</symmetric_stiffness><tension>{tension}</tension><fric_coeff>0</fric_coeff>')
    L = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">', '<Module type="solid"/>',
         f'<Control><analysis>STATIC</analysis><time_steps>{steps}</time_steps><step_size>{1.0/steps:.12g}</step_size>'
         f'<time_stepper type="default"><max_retries>10</max_retries><opt_iter>10</opt_iter><dtmin>{1e-4/steps:.12g}</dtmin>'
         f'<dtmax>{1.0/steps:.12g}</dtmax><aggressiveness>0</aggressiveness><cutback>0.5</cutback></time_stepper>'
         '<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness></solver></Control>',
         f'<Material><material id="1" name="tissue" type="neo-Hookean"><density>1</density><E>{float(young_pa):.12g}</E><v>{float(nu):.12g}</v></material>',
         '<material id="2" name="bed" type="rigid body"><density>1</density><center_of_mass>'
         + ",".join(f"{v:.9g}" for v in np.asarray(bedV).mean(0)) + '</center_of_mass></material></Material>',
         '<Mesh><Nodes name="tissue">'] + [f'<node id="{i+1}">{x:.12g},{y:.12g},{z:.12g}</node>' for i, (x, y, z) in enumerate(Xf)]
    retract = np.asarray(retract_m, float)
    bedV_written = np.asarray(bedV, float) - retract          # the bed starts retracted and is advanced back
    L += ['</Nodes><Nodes name="bed_n">'] + [f'<node id="{N+i+1}">{x:.12g},{y:.12g},{z:.12g}</node>' for i, (x, y, z) in enumerate(bedV_written)] + ['</Nodes>']
    L += ['<Elements type="tet4" name="tissue_e">'] + [f'<elem id="{i+1}">{a+1},{b+1},{c+1},{d+1}</elem>' for i, (a, b, c, d) in enumerate(Tf)] + ['</Elements>']
    L += ['<Elements type="tri3" name="bed_e">'] + [f'<elem id="{ne+k+1}">{a+N+1},{b+N+1},{c+N+1}</elem>' for k, (a, b, c) in enumerate(bedF)] + ['</Elements>']
    pairs = []
    if len(hf): L.append(f'<Surface name="held">{tri(hf, 1)}</Surface>'); pairs.append(("held", 1))
    if len(ff): L.append(f'<Surface name="free">{tri(ff, 1)}</Surface>'); pairs.append(("free", 0))
    L.append(f'<Surface name="bed_s">{tri(np.asarray(bedF), N + 1)}</Surface>')
    by_axis = {a: sorted(rank[n] + 1 for n, ax in pins if ax == a) for a in range(3)}
    for a, ids in by_axis.items():
        if ids: L.append(f'<NodeSet name="pin{a}">' + ",".join(map(str, ids)) + '</NodeSet>')
    L += [f'<SurfacePair name="{s}_pair"><primary>{s}</primary><secondary>bed_s</secondary></SurfacePair>' for s, _ in pairs] + ['</Mesh>']
    L += ['<MeshDomains><SolidDomain name="tissue_e" mat="tissue"/><ShellDomain name="bed_e" mat="bed"><shell_thickness>0.001</shell_thickness></ShellDomain></MeshDomains>']
    L += ['<Boundary>'] + [f'<bc name="p{a}" node_set="pin{a}" type="zero displacement"><x_dof>{int(a==0)}</x_dof><y_dof>{int(a==1)}</y_dof><z_dof>{int(a==2)}</z_dof></bc>'
                           for a, ids in by_axis.items() if ids] + ['</Boundary>']
    moving = [a for a in range(3) if retract[a] != 0.0]
    L += ['<Rigid>', '<rigid_bc name="fix_bed" type="rigid_fixed"><rb>bed</rb>'
          + "".join(f'<R{"xyz"[a]}_dof>{0 if a in moving else 1}</R{"xyz"[a]}_dof>' for a in range(3))
          + '<Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof></rigid_bc>']
    L += [f'<rigid_bc name="advance_{"xyz"[a]}" type="rigid_displacement"><rb>bed</rb><dof>{"xyz"[a]}</dof>'
          f'<value lc="1">{float(retract[a]):.12g}</value><relative>0</relative></rigid_bc>' for a in moving] + ['</Rigid>']
    L += ['<Contact>'] + [f'<contact type="sliding-elastic" name="c_{s}" surface_pair="{s}_pair">{cp(t)}</contact>' for s, t in pairs] + ['</Contact>']
    L += ['<LoadData><load_controller id="1" name="lc1" type="loadcurve"><interpolate>LINEAR</interpolate><points><pt>0,0</pt><pt>1,1</pt></points></load_controller>',
          f'<load_controller id="2" name="penalty_ramp" type="loadcurve"><interpolate>LINEAR</interpolate><points><pt>0,{float(ramp_from):.12g}</pt><pt>{float(ramp_until):.12g},1</pt><pt>1,1</pt></points></load_controller></LoadData>',
          '<Output><logfile><node_data data="ux;uy;uz" file="nodes.txt" delim="," node_set="tissue"></node_data></logfile></Output>', '</febio_spec>']
    Path(path).write_text("\n".join(L) + "\n")


def split_base_faces(boundary_faces, base, held_nodes):
    """Base faces (all three nodes in the base) split into held (all three held) and free. Returns
    (held_faces, free_faces, held nodes that sit in no all-held face and so are one-way in FEBio)."""
    isb = np.zeros(int(max(boundary_faces.max(), np.max(base))) + 1, bool); isb[base] = True
    ish = np.zeros_like(isb); ish[held_nodes] = True
    F = boundary_faces[isb[boundary_faces].all(1)]
    h = ish[F].all(1)
    covered = np.zeros_like(isb); covered[F[h].ravel()] = True
    return F[h], F[~h], np.flatnonzero(ish & ~covered)
