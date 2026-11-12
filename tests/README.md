# GENEPI tests

Run from the top directory after building:

```
make                    # builds ./genepi
make -C tests           # builds tests/xsec_points
tests/run_smoke.sh      # ~30 s
tests/run_regression.sh # ~2 min
```

The same steps run in GitHub Actions on every push and pull request to `OSG`
(`.github/workflows/ci.yml`).

| Script | What it checks |
|---|---|
| `run_smoke.sh` | Short fixed-seed runs for every target (p, n, d, He4, N14, NH3, ND3) and process (DVCS, pi0, eta), each validated by `check_lund.py`; seed reproducibility; command-line validation and exit codes (bad options, inverted ranges, invalid seed/target, empty phase space, `--max_trials`, missing GPD table). |
| `check_lund.py` | One LUND file: format, finite values, E²−p² = m², momentum conservation, allowed PIDs, header target A/Z consistent with the remnant, DVCS photon angle cuts, vertex inside the raster ellipse. |
| `xsec_points` | DVCS/BH/interference cross-section of `src/dvcs.cc` at 240 fixed points (proton/neutron, 4 kinematic points, 5 φ values, beam and target helicities) compared with `ref/xsec_points.ref` (relative tolerance 1e-6). |
| `check_distributions.py` | Means of x_B, Q², −t, cos φ, cos 2φ, near-beam photon fraction and Fermi-motion isotropy of short samples (`regression_samples.txt`) compared with high-statistics references in `ref/distributions.json` (fails beyond 5 standard errors). |

## After an intended physics change

A change of the cross-section or of the sampling makes the regression tests fail on purpose.
Regenerate the references and commit them together with the change, explaining why:

```
tests/make_references.sh        # 20000 events per sample, several minutes
git add tests/ref && git commit
```

Not covered yet (known issues): free-neutron meson output (pi0 written with one decay photon),
and eta decays (intermediate pi0 written together with their photons). Detector simulation
(gemc) and reconstruction (coatjava) are deliberately not part of these tests.
