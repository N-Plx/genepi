# Options

```
genepi [--option value ...]
```

| Option | Default | Meaning |
|---|---|---|
| `--trig` | 10000 | number of events to write |
| `--seed` | generated | random seed, integer in [1, 2147483647]; printed in the log |
| `--output` | genepi.dat | output LUND file |
| `--Ebeam` | 10.6 | beam energy [GeV] |
| `--process` | 0 | 0 = DVCS + BH, 1 = meson production |
| `--meson` | 0 | 0 = π⁰, 1 = η |
| `--targ_A`, `--targ_Z` | 1, 1 | target, see [Physics](physics.md#nuclear-targets) |
| `--x_min`, `--x_max` | 0.001, 1 | x_B range |
| `--y_min`, `--y_max` | 0.001, 1 | y = P·q / P·k range |
| `--Q2_min`, `--Q2_max` | 1, 20 | Q² range [GeV²] |
| `--W2_min`, `--W2_max` | 4, 50 | W² range [GeV²] |
| `--nu_min`, `--nu_max` | 0.3, 11 | ν range [GeV] |
| `--t_min`, `--t_max` | −1.2, 0 | t range [GeV²] |
| `--ycol_min`, `--ycol_max` | −999, 0.025 | collider y range |
| `--gamma_beam_min` | 1 | DVCS: minimum photon–beam angle [deg] |
| `--gamma_e_min` | 1 | DVCS: minimum photon–scattered-electron angle [deg] |
| `--vx`, `--vy`, `--vz` | 0, 0, −3 | vertex centre [cm] |
| `--raster_x`, `--raster_y` | 0, 0 | raster half-axes [cm] |
| `--max_trials` | 0 (no limit) | stop after this many trials; also replaces the 10⁷-trial empty-phase-space limit |
| `--docker` | — | accepted and ignored (grid scripts) |

Invalid input (unknown option, missing or non-numeric value, min ≥ max, unsupported target,
…) stops the program with a message. Unique abbreviations of option names are accepted.

## Exit status

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | invalid input, GPD table not found, output cannot be opened, or fatal JETSET error |
| 2 | no event accepted after 10⁷ trials (or `--max_trials`): empty phase space |
| 3 | `--max_trials` reached after some events were kept (partial output written) |
| 4 | error while writing the output (e.g. disk full) |

## Files

`gpd_table.dat` is looked up in `$GENEPI` if set, otherwise next to the executable, then in the
current directory. The output file is opened only after all checks, so a failing job does not
overwrite a previous output.
