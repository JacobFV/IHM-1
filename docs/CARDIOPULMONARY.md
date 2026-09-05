# A resting, supine cardiopulmonary model

This view generates individual heartbeats and breaths. It is a nonlinear explicit
materialization of registered physical state and processes, distinct from the
coarse heart-rate and ventilation components.

The blood loop is:

```
right atrium -> tricuspid valve -> right ventricle -> pulmonary valve
-> pulmonary artery -> pulmonary vascular resistance -> pulmonary vein
-> left atrium -> mitral valve -> left ventricle -> aortic valve
-> systemic artery -> systemic vascular resistance -> systemic vein
-> right atrium
```

Each of the eight blood compartments has a volume. Every edge transfers equal
volume out of one compartment and into the next, preserving total blood volume.
Heart chambers use `P = E(phase) (V - V_unstressed) + P_pleural`. Vessels have
linear compliance; intrathoracic vascular compartments also see pleural pressure.
The four ideal valves open when upstream pressure exceeds downstream pressure.
Pacing phase advances continuously and ventricular/atrial activation changes
elastance, producing contraction, filling and pulsatile ejection.

Respiratory phase drives pleural pressure. Lung volume follows the pressure
difference between the mouth and alveolus through an airway resistance:

```
P_alveolar = P_pleural + (V_lung - V_relaxed) / C_lung
dV_lung/dt = (P_mouth - P_alveolar) / R_airway
```

This produces inspiration and expiration. Pleural pressure simultaneously affects
the heart and pulmonary vessels, so the traces are mechanically coupled in one direction: breathing affects circulation, with no
circulatory feedback into respiratory effort.
Supine posture means no imposed standing hydrostatic gradients or locomotor load.

## Run

```bash
.venv/bin/ihm build cardiopulmonary --subject supine-example \
  --time 30 --output artifacts/cardiopulmonary.npz
.venv/bin/ihm forecast artifacts/cardiopulmonary.npz --times 30 30.1 30.2 30.3
.venv/bin/python scripts/plot_cardiopulmonary.py
```

```python
import numpy as np
from ihm.materialize.cardiopulmonary import cardiopulmonary_model
from ihm.processes.cardiopulmonary import Parameters

model = cardiopulmonary_model('example', Parameters(heart_rate_bpm=70,
                                                  respiratory_rate_per_min=12))
model.advance(30)  # allow the illustrative initial condition to settle
trace = model.forecast(np.linspace(30, 45, 1501))
pressure = trace['signals']['systemic_arterial_pressure_mmHg']
airflow = trace['signals']['airflow_L_s']
```

`--parameters params.json` accepts the `Parameters` schema. Forecasts are bounded
to 600 seconds and 100,000 output points per call; sampling times do not determine
the adaptive solver step. Save/load preserves phase and state. The nonlinear
executor deliberately rejects unsupported data assimilation and hard clamps.

## Scope and evidence

Structural inspiration: [Physiome Windkessel model](https://models.cellml.org/workspace/cardiovascular_circulation_windkessel/file/5431a35ead45c625c50e08be47fdf023b492ca70/MainWindkessel.cellml),
[coupled cardiovascular/respiratory compartment modeling](https://pmc.ncbi.nlm.nih.gov/articles/PMC10581768/)
and [spontaneous respiratory mechanics](https://pmc.ncbi.nlm.nih.gov/articles/PMC4303416/).
This is an original simplified implementation of those classes of equations,
not a reproduction of any paper's fitted model. Default parameter values are
illustrative and cannot establish individual predictive validity.

The pacemakers are prescribed phase oscillators. ECG waveforms, cardiac ion-channel
electrophysiology, arrhythmias, baroreflexes, chemoreflexes, gas exchange, oxygen
consumption and carbon-dioxide chemistry are not implemented in this view. Its
blood pressures and respiratory volumes are mechanical predictions, not evidence
that the whole body is metabolically self-sustaining. It does not yet assimilate
measurements or propagate quantified uncertainty. Other materializations retain
Gaussian evidence conditioning.

Movement would require a posture-dependent vascular geometry, musculoskeletal
mechanics, work-dependent metabolism and neural/reflex feedback. Holding posture
supine is a useful scope reduction while the resting system is developed.
