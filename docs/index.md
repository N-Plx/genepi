# GENEPI

GENEPI is a Monte Carlo event generator for exclusive electroproduction off free and bound
nucleons, used by CLAS12 (Jefferson Lab) and run in grid productions on OSG:

- **DVCS + Bethe-Heitler** (e N → e N γ), including their interference;
- **π⁰ / η electroproduction** (e N → e N π⁰/η), decayed with JETSET.

Targets: free proton and neutron, deuterium, ⁴He, ¹⁴N, and the polarized-target materials
NH₃ and ND₃ (the nucleus is chosen event by event). Events are unweighted (keep/reject on the
cross-section) and written in the LUND format read by gemc.

Original documentation: CLAS note 2009-024
(<https://misportal.jlab.org/ul/Physics/Hall-B/clas/viewFile.cfm/2009-024.pdf?documentId=554>).

## Quick start

```
make                                  # needs ROOT (root-config) and gfortran
./genepi --trig 10000 --targ_A 20 --targ_Z 10 --seed 12345 --output nd3_dvcs.lund
```

`gpd_table.dat` is found next to the executable (or in `$GENEPI`). Run `make -C tests` and
`tests/run_smoke.sh` to check a build.

## Contents

| Page | |
|---|---|
| [Physics model](physics.md) | Processes, cross-sections, GPD/CFF tables and their limits, nuclear targets |
| [Sampling and weights](sampling.md) | What is sampled, Jacobians, keep/reject, maximum weight, photon angle cuts |
| [Options](options.md) | Command line, seeds, exit codes |
| [Output format](output.md) | LUND columns, particle IDs, gemc and coatjava compatibility |
| [Validation and limitations](validation.md) | Tests, CI, known limitations |
| [Changes](changelog.md) | History of the OSG branch |
