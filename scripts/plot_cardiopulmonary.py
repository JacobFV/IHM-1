"""export actual simulated beat and breath trajectories, not drawn waveforms."""
from pathlib import Path
import json
import numpy as np
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ihm.materialize.cardiopulmonary import cardiopulmonary_model

out = Path('artifacts'); out.mkdir(exist_ok=True)
m = cardiopulmonary_model(); m.advance(30)
result = m.forecast(np.linspace(30, 45, 1501))
t = result['time']-30; s = result['signals']
fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
axes[0].plot(t, s['systemic_arterial_pressure_mmHg'], color='#d65342', label='Systemic arterial pressure')
axes[0].set_ylabel('Pressure (mmHg)'); axes[0].legend(loc='upper right', frameon=False)
axes[1].plot(t, s['aortic_flow_mL_s'], color='#973d77', label='Aortic valve flow')
axes[1].set_ylabel('Flow (mL/s)'); axes[1].legend(loc='upper right', frameon=False)
axes[2].plot(t, s['lung_volume_L'], color='#247ca7', label='Lung volume')
axes[2].set_ylabel('Lung volume (L)'); axes[2].set_xlabel('Seconds after 30-second settling period')
axes[2].legend(loc='upper right', frameon=False)
for axis in axes:
    axis.grid(alpha=.18); axis.spines[['top','right']].set_visible(False)
fig.suptitle('IHM-1 • beating heart and tidal breathing, lying supine\nIllustrative mechanistic parameters; not a calibrated human prediction', fontsize=14)
fig.tight_layout(); fig.savefig(out/'cardiopulmonary.png', dpi=180); fig.savefig(out/'cardiopulmonary.svg')
beats, _ = find_peaks(s['aortic_flow_mL_s'], prominence=10)
breaths, _ = find_peaks(s['lung_volume_L'], prominence=.1)
metrics = {'simulated_heart_rate_bpm': float(60/np.mean(np.diff(result['time'][beats]))),
           'simulated_respiratory_rate_per_min': float(60/np.mean(np.diff(result['time'][breaths]))),
           'tidal_volume_L': float(np.ptp(s['lung_volume_L'])),
           'mean_arterial_pressure_mmHg': float(np.mean(s['systemic_arterial_pressure_mmHg'])),
           'mean_cardiac_output_L_min': float(np.mean(s['aortic_flow_mL_s'])*.06),
           'validated_biology': False}
(out/'cardiopulmonary-metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
print(json.dumps(metrics, indent=2))
