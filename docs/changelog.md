# Changes on the OSG branch

Changes since commit `16e6d38` (see `git log` for details).

## Physics fixes
- Neutron CFFs H and E: flavour form-factor normalisation (isospin mirror of the proton);
  the neutron beam-spin asymmetry changed sign.
- DVCS²: Ẽ Ẽ* used Re H̃; spurious helicity-dependent s₁ terms removed (BMK eqs. 44, 47).
- GPD table interpolation: off-by-one in the ξ step (CFFs were evaluated at ξ about 5% too high).
- Jacobian of the event-dependent t range in the keep/reject (high-x_B events were over-produced).
- Fermi motion: isotropic nucleon direction (was flat in θ); nitrogen sampler took the random
  generator by copy (momentum and angle correlated); nitrogen envelope below the peak.
- Meson momentum on the Fermi-motion path (missing square root; no effect on the output).
- DVCS photon angle cuts `--gamma_beam_min`, `--gamma_e_min` (default 1°) removing the BH
  collinear poles.

## Output
- Header target A/Z per event (was always 2/1); real ³He/triton codes for ⁴He (the fake codes
  1234/4321 made gemc exit).

## Robustness (grid production)
- Output file opened once per run (one file handle leaked per event; jobs crashed after ~1000
  events).
- Seeds: no more time-of-day collisions; JETSET seeded; range [1, 2³¹−1].
- Raster vertex uniform in the ellipse and reproducible.
- Strict command-line validation; no endless loop on empty phase space; `--max_trials`;
  capped log messages; meaningful exit codes; `--output`; GPD table found next to the
  executable; output opened only after all checks.
- Reports keep/reject weights above the maximum.

## Repository
- Removed tracked object files and libraries, the duplicate `jetset/` and `src/gpd_table.dat`.
- Makefile builds unmodified on Linux and macOS (`root-config` from `$ROOTSYS` or the PATH).
- Tests (`tests/`), GitHub Actions CI and this documentation.
