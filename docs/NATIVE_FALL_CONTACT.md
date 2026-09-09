# Native upright fall contact

The retained walking source supplies twelve foot contact spheres. Until this
change, its torso and limbs could pass through the floor during an uncontrolled
fall. In a bare 70 kg run the pelvis reached -0.99449 m at 1.02 s and subsequently
oscillated beneath the floor. That run continued, unlike the previously reported
1.026 s integrator failure at another mass/actuation configuration. Reaching a
later timestamp alone is therefore insufficient evidence of physical stability.

`scripts/native_mechanical_stream.cpp` now adds a unilateral floor contact sphere
at each non-foot segment's source center of mass. Its radius is the smallest
ellipsoid semi-axis inferred from that segment's mass and diagonal inertia:

`r² = min_k 5 (I_(k+1) + I_(k+2) - I_k) / (2 m)`.

The source `SmoothSphereHalfSpaceForce` law and floor geometry are retained. The
foot assembly (talus, calcaneus, toes) keeps its original contact representation;
the talus's small placeholder inertia is unsuitable for geometric inference.
There is no new pose constraint, tether, gravity compensation, artificial
root force, or integration tolerance relaxation. Force reactions enter the
existing momentum audit and the actual articulated integration.

These are **uncalibrated engineering collision proxies**, not reconstructed skin,
full segment capsules, self-contact, or a balance policy. They catch a fall;
they do not make the body stand or walk. Native frames identify the contact model
as `source_feet_and_inertia_inscribed_body_spheres`, and `execution.json` discloses
the support scope.

Run the real native regression with:

```sh
.venv/bin/python scripts/build_native_mechanical_stream.py
.venv/bin/python scripts/verify_native_fall_contact.py
.venv/bin/python scripts/verify_native_fall_contact.py \
  --augmented-registration data/derived/mechanics/whole_body_arm26_v2/registration.json
```

The test advances five simulated seconds, checks that the pelvis remains above
the floor, requires non-foot contact reaction, and checks total force balance
below 1e-5 N. Each run saves its sampled trajectory and build identity below
`data/derived/native-fall-contact-*`. The standing source pose must not initially
overlap the new proxies with the floor. This test is deliberately independent of
physiology and cortical control; those integrations require separate evidence.

The final body-contact dynamics passed the bare 80-muscle five-second run at
70.7713 kg (`data/derived/native-fall-contact-iqw7x29k/report.json`). At 4.02 s,
pelvis height was 0.01820 m, kinetic energy 0.00898 J, non-foot reaction 651.54 N,
and momentum residual 3.41e-13 N. Over the full run, minimum pelvis height was
0.016970 m and maximum momentum residual was 3.28e-12 N. Peak non-foot impact
reaction was 5749 N. At five seconds kinetic energy was 0.004526 J. The body had
fallen and settled. This is not an upright stance result.

The augmented 92-muscle model also completed five seconds at the same mass:
`data/derived/native-fall-contact-cdd1xjq6/report.json`. It passed the same floor
and momentum checks. Its run took 619.7 wall seconds with concurrent native jobs;
this does not establish real-time execution.

The native stream also exposes a selective read-only
`moment_arms(muscles=[...], coordinates=[...])` query for rotational coordinates.
It returns the source muscle path derivatives at the current state. The native
test `scripts/verify_native_moment_arms.py` measured opposite ankle moment arms
for soleus (-0.049708 m) and tibialis anterior (+0.048273 m), and verified exact
subsequent state/trajectory equality with a branch that omitted the query.
