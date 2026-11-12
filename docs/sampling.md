# Sampling and weights

## Phase space

For each trial GENEPI samples Q² and x_B flat in the requested ranges, the photon/meson
azimuth φ flat in [0, 2π), and t flat in the **event-dependent** physical range
[max(t_kin,min, t_min), min(t_kin,max, t_max)]. For bound nucleons the Fermi momentum is drawn
first (magnitude from the nucleus distribution, isotropic direction) and t is generated in the
γ*N centre-of-mass frame (`src/fmotion.cc`). Cuts on y, W², ν and y_col are then applied.

Because the t range depends on the event, the keep/reject weight is

    w = dσ × Δt / (t_max − t_min)        (Δt = width of the event's t range)

so that kept events follow dσ. The weight written in the LUND header is dσ alone.

## Keep / reject and the maximum weight

An event is kept if w > u·w_max with u uniform in [0,1) and a fixed **w_max = 2** (DVCS) or
**800** (mesons). Trials with w > w_max are always kept, so the region where the cross-section
exceeds w_max is under-represented. GENEPI reports at the end of the run:

```
Max keep/reject weight seen   ...
Trials with weight above max   ...
Fraction of total cross-section lost to max-weight clipping   ...
```

For mesons nothing is clipped. For DVCS at 10.6 GeV with the default cuts and the 1° photon
cuts, about 25–30% of the sampled cross-section is clipped (it was ~70% without the photon
cuts), mostly near the remaining BH peaks and at small −t. A complete fix (weighted events or
an adaptive maximum, e.g. FOAM as in EpIC) is still to be decided.

## Photon angle cuts (DVCS)

`--gamma_beam_min` and `--gamma_e_min` (degrees, default 1 each) reject events whose photon is
closer than this to the beam or to the scattered electron, before the cross-section is
evaluated. They remove the BH collinear poles: without them, ~55% of kept DVCS events at
10.6 GeV had the photon within 2.5° of the beam and keep/reject weights reached ~4000. Setting
both to 0 restores the previous output exactly.

## Random numbers

One seed initializes ROOT's TRandom1 (kinematics), gRandom (vertex) and JETSET (decays).
Without `--seed`, a seed is generated from the system random device, the clock, the process id
and the host name, and printed. Valid seeds: 1 to 2147483647 (TRandom1 uses 31 bits).

## Vertex

x, y uniform inside the ellipse of half-axes (raster_x, raster_y) centred on (vx, vy); z = vz.
gemc can override it (`RANDOMIZE_LUND_VZ`, `RASTER_VERTEX` in the production scripts).
