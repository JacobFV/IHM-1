# Acquired human tissue constitutive evidence

The institutional repository full text of Khorshidi et al. (2024),
[DOI 10.1016/j.actbio.2024.06.035](https://repository.rcsi.com/articles/journal_contribution/26563843),
is now held with metadata, repository version/file identity, MD5, SHA-256 and
layout-preserving text extraction. The preceding 2023 experimental-design
paper is also held. `scripts/acquire_penile_mechanics.py` verifies retained
bytes and does not replace a previous acquisition.

The 2024 study's Table 1 contains actual inverse-FE point estimates. The
transcription preserves the cohort and layer-specific evidence distinctions:
three human donors aged75–89, frozen postmortem tissue, tests in37±2°C PBS.
Whole-organ samples, individual tunica tension and cavernosum compression
constrain different parts of the inference. Fascia lacks separate specimen
data, and the CS fit represents the whole CS including urethra. These are
neither measured glans coefficients nor a calibrated generic living adult.

| Layer | Executable law | Initial matrix shear modulus | Initial bulk modulus |
|---|---|---:|---:|
| Cavernosum | Two-term compressible Ogden |22 Pa|1111.11 Pa|
| Spongiosum including urethra | Two-term compressible Ogden |600 Pa|20000 Pa|
| Fascia | Split neo-Hookean |50 Pa|2500 Pa|
| Tunica albuginea | HGO with explicit reference fiber |200000 Pa|200000000 Pa|

`ihm/calibration/penile.py` implements the published energy functions and
analytic first-Piola stress, with SI conversions from kPa and kPa⁻¹.
`energy_piola(F)` returns reference-volume energy density and its gradient
with respect to the deformation gradient. HGO retains its distinct
volumetric term and tension-only fiber contribution; assigning tunica requires
an explicit unit reference-space fiber direction. The scalar shear entry for
tunica excludes its anisotropic fiber tangent.

The paper prints `D2=0` while its general Ogden expression contains `1/D2`.
The implementation explicitly interprets zero as omitting that higher-order
volumetric term. Native Abaqus execution parity has not been established.
This qualification stays attached to every materialization.

```python
from ihm.human import ImplicitHuman
law = ImplicitHuman.open().materialize(
    'penile-constitutive', tissue='corpus_cavernosum')
energy_density, first_piola = law.energy_piola(deformation_gradient)
```

The implicit body binds the full text and parameter table by content hash;
saved manifests support both JSON and binary evidence. Unsupported tissue
names and incomplete fiber assignments fail. Existing generic-prior material
and garment experiments retain their original identities and results.

`scripts/verify_penile_constitutive.py` checks zero reference stress, frame
indifference, analytic stress against independent energy differences, pure
dilation, Table3 initial moduli, fiber anisotropy and scope rejection. These
are constitutive implementation checks, not a new reproduction of the paper's
full-organ fit or proof of clothing containment. Geometry, layer interfaces,
vascular pressure, viscosity, damage and population discrepancy remain separate.
