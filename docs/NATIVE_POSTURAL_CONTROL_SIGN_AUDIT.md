# Native postural control sign correction

The earlier analytical inference that a positive lateral COM error should request positive right/negative left hip-adduction torque is **rejected by native pulse evidence**. It conflated pelvis orientation with whole-body translation under foot contact.

In the retained matched-checkpoint experiment, right hip-adduction +20 Nm and left −20 Nm increased COM velocity toward positive source Z by approximately 0.01536 m/s at 50 ms and 0.03762 m/s at 150 ms, while pelvis list decreased. That input therefore reinforces a positive lateral COM displacement rather than correcting it.

The controller now uses `lateral_torque = -Kz * (COMz - target_z) - Dz * COMvz`, distributed as half to the right hip and minus half to the left. Pelvis-list orientation feedback retains its separate sign and uses error from the initial reference. Both lateral COM gains default to zero; a gain value requires separate closed-loop evaluation.

Evidence: `data/research/locomotion_control/native_lateral_pulse_hhtz1q9x/report.json`, including retained experiment/controller source snapshots. Regression coverage: `scripts/verify_moment_arm_control.py`. A local pulse response is not evidence of sustained balance, walking, or cortical control.
