# Native reference mechanics benchmark

The reference backend is a repository-local, unmodified FEBio 4.13 build pinned to commit `32ae206ff4881dfb54f62296cd1558e58ed9fcc6`. This host is Linux AArch64 and already provides CMake, GNU C++, and OpenMP. No system package installation is required for the selected configuration. The built-in skyline solver provides the symmetric direct solve; MKL, HYPRE, SuperLU_MT, MMG, LEVMAR, NLOPT, PDL, FFTW, and ZLIB are explicitly disabled. These disabled features are not tested or claimed available.

The native build and compression benchmark passed. Maximum error across the initial state and all ten increments was `4.4448711378208827e-10 Pa` for normal stress and `1.78792204428e-16 m` for displacement. FEBio reported normal termination and selected skyline with 20 equations. Final axial stress was `-105.555555556 Pa` in both elements and interior axial displacement was `-0.05 m` at all four middle-layer nodes.

The first build exposed an upstream no-MKL Linux link omission: OpenMP-compiled modules had unresolved `omp_*` and `GOMP_*` symbols. The retained `initial-openmp-link-failure.log` documents it. The builder supplies `-fopenmp` through CMake executable and shared-library linker flags; the original FEBio source remains unchanged. Rebuilding with those flags and running the native benchmark succeeded.

The [official build guide](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/BUILD.md) explains optional dependencies. The pinned [dependency discovery](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/FindDependencies.cmake) and [solver initialization](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/NumCore/NumCore.cpp) establish the no-MKL configuration. This is an actual source build, not evidence inferred from a download page.

## Reproduction and retained evidence

From the repository root:

```sh
python3 scripts/build_reference_mechanics.py --jobs 8
python3 scripts/verify_reference_mechanics.py
```

The first command clones the public source if absent, rejects a mismatched or modified source tree, archives the exact commit, and configures/builds inside the repository. `data/raw/mechanics/FEBio-4.13/` retains source and its upstream license; `data/raw/mechanics/FEBio-4.13.tar` retains an exact source archive. The builder records its SHA-256, commands, exit codes, platform, executable hash, shared-library hashes, and linked-library inventory in `data/derived/mechanics-reference/`. No `cmake --install` or global writes occur. Reproduction on another machine requires its own compiler toolchain and produces its own hashes.

The second command checks the native executable and library hashes before executing a generated FEBio model. Its `compression/` directory retains input XML, native stdout/log, nodal displacements/reactions, element stresses, XPLT, and a `benchmark.json` with numerical comparisons and output hashes. Failures raise a nonzero exit status; a build manifest alone never means the benchmark passed.

## Independent expected solution

Two stacked eight-node hexahedra fill a unit cube. Young's modulus is 1000 Pa and Poisson's ratio is zero. The bottom is fixed axially, the top moves downward by 0.1 m over ten quasi-static increments, and symmetry constraints fix normal displacement at the x=0 and y=0 planes. Other lateral displacements and all middle-layer axial displacements remain unknowns. There is no gravity or contact.

For the compressible neo-Hookean energy `W = mu/2 (I1 - 3) - mu log J + Lambda/2 (log J)^2`, `mu=500 Pa` and the first Lamé parameter `Lambda=0`. Homogeneous deformation has `F=diag(1,1,s)` and `J=s`, with stretch `s=1-0.1t`. Differentiating the energy yields Cauchy stress `sigma = mu/J (FF^T-I)`. Thus lateral stresses vanish, `sigma_zz=500(s^2-1)/s`, and `u_z=(s-1)Z`. At the final increment, `sigma_zz=-105.5555556 Pa` and the middle layer displacement is `-0.05 m`. The formula is evaluated directly in the checker, without using native material routines or fitting output. The upstream [material definition](https://github.com/febiosoftware/FEBio/blob/32ae206ff4881dfb54f62296cd1558e58ed9fcc6/FEBioMech/FENeoHookean.cpp) confirms the constitutive convention.

Every reported increment is checked for finite values, both elements' three normal stresses, all twelve nodes' three displacements, and the final simulation time. Absolute tolerances are 0.002 Pa and 0.000002 m, allowing for FEBio's text-output precision. Native normal termination and successful process exit are also required.

## Scope

This is a homogeneous finite-strain constitutive and equilibrium reference. It does not validate collision/contact, heterogeneous anatomy, articulated joints, bed equilibrium, all-organ interaction, or patient material calibration. The application's affine/reduced mechanics has not been compared because an equivalent constitutive law, geometry, loading and boundary-condition endpoint has not been established. Installing this reference does not replace that model or imply its validation. The separate [rigid-wall contact benchmark](REFERENCE_CONTACT_BENCHMARK.md) now passes a confined analytical limit and penalty refinement; nonuniform anatomical contact remains a separate gate.
