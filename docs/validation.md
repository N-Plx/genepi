# Validation and known limitations

## Tests

See `tests/README.md`. Every push and pull request to `OSG` runs in GitHub Actions
(`.github/workflows/ci.yml`, ROOT container):

1. build with the repository Makefile, compiler warnings in the job summary;
2. smoke tests: every target and process, each output validated (format, E²−p² = m²,
   momentum conservation, PIDs, header A/Z vs remnant, photon cuts, vertex), option checks
   and exit codes, seed reproducibility;
3. cross-section regression at 240 fixed points against `tests/ref/xsec_points.ref`;
4. distribution regression of short samples against high-statistics references.

An intended physics change must regenerate the references (`tests/make_references.sh`) in the
same commit. Detector simulation and reconstruction are not part of these tests.

## Known limitations

| Topic | Status |
|---|---|
| Keep/reject maximum weight | fixed w_max; ~25–30% of the DVCS cross-section clipped (reported in the log) |
| BH | massless-electron formulas; collinear poles removed by the photon angle cuts |
| GPD model | ξ-only tables, no Q² evolution; x_B > 0.0198; overall sign convention of the interference not yet cross-checked against data |
| Polarization | beam and target helicities fixed to 0 |
| Radiative corrections | none |
| Meson model | no ε dependence, no σ_TT/σ_LT; free-neutron meson output writes the π⁰ with only one photon; η intermediate π⁰ written with their photons |
| Nuclear effects | on-shell nucleons, no flux factor, no coherent channel, no target dilution |
| GPD models 5/6 | disabled; known bugs |
