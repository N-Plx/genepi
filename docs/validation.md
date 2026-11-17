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