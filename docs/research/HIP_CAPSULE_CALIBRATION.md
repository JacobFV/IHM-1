# Hip capsule evidence and bounded research calibration

No production passive force is justified yet. The retained evidence supports
an opt-in **fixed-pose scalar research interpolation** at zero hip flexion and
zero ab/adduction, within 5 Nm aggregate restraint. It does not identify a
three-DOF hip potential, a patient-specific spring, or a force onset at the
existing ±0.6981317 rad coordinate bounds. No native model or ROM was changed.

## Retained primary evidence

Original XML/PDF/metadata bytes, licenses, specimen conditions and SHA256
receipts are in `data/research/hip_capsule/sources.json`. The source cards
separate these studies and do not pool different preparations.

- [Van Arkel, Amis & Jeffers 2015](https://doi.org/10.1016/j.jbiomech.2015.09.002)
  supplies the calibration evidence: eight analyzed rotation specimens after
  morphology exclusion and capsule rupture, older fresh-frozen left hips,
  110 N compression directed 20° medially/proximally, room temperature, and
  the second 10-second sinusoidal cycle. It measures pose-dependent slack,
  transition-to-5 Nm travel and tangent stiffness. Published figures provide
  means/CIs; specimen curves were not recovered. The measured restraint also
  includes labral effects. License: CC BY 4.0.
- [Van Arkel et al. 2015 resection study](https://doi.org/10.1302/0301-620X.97B4.34638)
  isolates relative contributions by sequential cutting in nine cadaveric hips
  at original intact ±5 Nm angles. It supports capsule dominance while
  preserving labral/teres contributions. Superposition results do not supply
  complete elastic constitutive laws. License: CC BY, version unspecified in
  retained permissions.
- [Karunaseelan et al. 2021](https://doi.org/10.1302/2046-3758.109.BJR-2020-0536.R1)
  combines those experiments with Twente enthesis centroids and ArtiSynth
  wrapping to estimate ligament moment arms and forces. It reuses earlier
  specimens. Executable MATLAB/ArtiSynth models were not recovered. Its
  CC BY-NC-ND 4.0 article is retained unchanged; that license is not a license
  for our original implementation.
- [Myers et al. 2020](https://doi.org/10.1080/10255842.2020.1764543)
  describes a calibrated probabilistic FE capsule. Only primary bibliographic
  metadata/abstract were acquired; the provider lists subscription access.
  No executable model or reuse license was established.
- [Anantha-Krishnan et al. 2024](https://doi.org/10.3390/bioengineering11010037)
  uses six FE membrane sectors with tension-only springs for an implanted hip,
  with validation against five THA specimens. Data are available upon request;
  no regression weights or calibrated model files were recovered. THA-specific
  geometry/pre-strain cannot be transplanted into a native hip. CC BY 4.0.
- [Pieroh et al. 2016](https://doi.org/10.1371/journal.pone.0163306)
  reports excised ligament-strip material measurements after ethanol-glycerin
  embalming. These preparation-sensitive data cannot directly determine a
  fresh-frozen whole-joint torsional law. The source has inconsistent counts;
  cards preserve that issue. License: CC BY, version unspecified in permissions.

## Concrete scalar law

`neutral_slice.json` records a provisional F0/A0 mean reconstruction. Figure 4
supplies slack width ≈27.5°; Figure 6 supplies transition travel ≈14.5° and
5 Nm tangent ≈0.82 Nm/degree. Figure 5's published neutral regression gives
mid-slack `c=5.7°` internal rotation at zero flexion. These are plot-derived
summaries and a regression intercept, not raw specimen observations.

The retained PDF and unmodified rendered pages have pixel coordinates for every
manual selection. Approximate published 95% CIs are 11.7–42.6° for slack width,
11.6–17.2° for travel and 0.67–0.97 Nm/degree for tangent. Plot-reading uncertainty
is separately recorded as 0.8°, 0.4° and 0.02 Nm/degree. These are CIs of means,
not individual prediction intervals. Their covariance and uncertainty in the
regression intercept are unavailable; corner combinations are sensitivity
scenarios, not a joint confidence envelope.

Let theta be ISB internal rotation, h half the slack width, d transition travel,
and K the tangent converted to Nm/radian. Set

```
x = max(0, abs(theta-c)-h)
p = K*d/5                         # nominal 2.378
V = 5*d/(p+1) * (x/d)^(p+1)       # joules
Q = -sign(theta-c)*5*(x/d)^p       # Nm, Q = -dV/dtheta
```

All angles inside this law are radians. The nominal potential is nonnegative,
with zero torque/stiffness in slack and continuous onset because p>1. At x=d
it matches 5 Nm and the digitized tangent (≈46.98 Nm/rad). Its nonlinear
interior is an engineering interpolation, not independently validated data.
The experimental onset used a 0.03 Nm/degree gradient threshold; zero modeled
onset is an explicit approximation. Equal direction-specific toe laws use
pooled means. No damping or hysteresis coefficient is inferred.

`ihm/assembly/hip_capsule.py` rejects changed loading conditions, nonzero flexion/adduction and rotations
outside `c±(h+d)`, nominally **−22.55° to +33.95°**. Thus both ±40° native ROM
bounds are outside this calibration. Rejection is deliberate: extrapolating
endpoint stiffness to stop the support solver would exceed measured evidence.
The implementation returns potential, restoring torque and positive tangent;
it supplies no native hook or registration.

## Multi-DOF and model integration requirements

A family of rotation slices is insufficient to identify V(flexion,adduction,
rotation). Even if slice parameters were interpolated across poses, a common
potential would generate flexion and adduction derivatives. Ignoring those
terms generally breaks conservative work balance; adding them without
cross-axis data invents uncalibrated forces. A pose-dependent integration
constant is also unidentifiable from rotational measurements alone.

Before any source-to-native mapping, establish the target hip joint convention,
left/right signs, neutral alignment and virtual-work Jacobian against ISB
anatomical axes. No sign or frame equivalence follows from matching coordinate
names. Acquire specimen torque-angle paths/cross-axis data or a reusable
validated ligament model with entheses, wrapping and slack lengths. Calibrate
against the intended anatomy and loading conditions, account for overlap with
labrum/teres and existing passive tissues, and validate independent poses.
Only then assess whether the physiological model explains the support residual.
The support result must not determine the desired spring constant.

## Verification

Under a 1 GiB address-space limit, five tests check endpoint torque and radian
units, energy-gradient/passivity, slack continuity, unmeasured-pose and ±40°
rejection, and all retained source hashes. Run:

```
OPENBLAS_NUM_THREADS=1 prlimit --as=1073741824 nice -n 10 .venv/bin/python -m scripts.verify_hip_capsule
.venv/bin/python data/research/hip_capsule/retain_sources.py --check
```
