"""export an illustrative skin network result; requires matplotlib."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ihm.materialize.skin import SkinPatch, skin_model, electric_field

out = Path('artifacts'); out.mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
for wound, label, color in ((False, 'Intact barrier', '#2176ae'), (True, 'Wound shunt', '#d45634')):
    patch = SkinPatch.line();
    if not wound:
        patch = SkinPatch.line(wound=False)
    model = skin_model(patch); model.advance(5)
    surface_idx = [model.components.index(f'skin.surface_{i}.potential') for i in range(len(patch.surface))]
    x = np.array([s.position_m[0] for s in patch.surface])*1000
    axes[0].plot(x, model.mean[surface_idx]*1000, 'o-', color=color, label=label)
    e = electric_field(model, patch)
    mid = [(x[a]+x[b])/2 for a, b, _ in patch.surface_edges]
    axes[1].plot(mid, e['V_per_m'], 'o-', color=color, label=label)
    cells = [model.components.index(f'skin.cell_{i}.membrane_voltage') for i in range(len(patch.cells))]
    cx = [c.position_m[0]*1000 for c in patch.cells]
    axes[2].scatter(cx, model.mean[cells]*1000, color=color, label=label)
    for a, b, _ in patch.gap_edges:
        axes[2].plot([cx[a],cx[b]], model.mean[[cells[a],cells[b]]]*1000, color=color)
axes[0].set_ylabel('Apical potential (mV; basal bath = 0)')
axes[1].set_ylabel('Signed lateral field (V/m)')
axes[2].set_ylabel('Membrane voltage (mV)')
for axis, title in zip(axes, ('Extracellular network', 'Field from surface potential', 'Separate cell-contact network')):
    axis.set_title(title); axis.set_xlabel('Position (mm)'); axis.grid(alpha=.2)
axes[0].legend(frameon=False)
fig.suptitle('IHM-1 skin materialization • illustrative parameters, not human validation', fontsize=13)
fig.tight_layout(); fig.savefig(out/'skin-bioelectricity.png', dpi=180); fig.savefig(out/'skin-bioelectricity.svg')
print(out/'skin-bioelectricity.png')
