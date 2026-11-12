# Physics model

## DVCS and Bethe-Heitler (`--process 0`)

The cross-section dσ/(dx_B dQ² dt dφ) (nb/GeV⁴, `x_sec_type = 3`) is the sum of the
Bethe-Heitler (BH), DVCS and interference terms in the twist-2 formalism of
Belitsky, Müller and Kirchner (BMK, hep-ph/0112108), implemented in `src/dvcs.cc`:

- BH: Fourier coefficients c₀, c₁, c₂ (unpolarized) and the longitudinally polarized target
  terms; propagators P₁, P₂ computed with a **massless electron**.
- DVCS²: c₀ and c₁ with the effective CFF approximation F_eff = −2ξ/(1+ξ) F. In that
  approximation the helicity-dependent s₁ terms vanish (BMK eqs. 44, 47) and are set to zero.
- Interference: c₀…c₂ and s₁, s₂, including the polarized-target terms; sign for an e⁻ beam
  (`beamcharge = −1`).
- φ is the Trento angle of the photon; φ_BMK = π − φ_Trento is used internally.
- Beam and target helicities are currently fixed to 0 (unpolarized samples).

### Compton form factors (CFFs)

The CFFs H, E, H̃, Ẽ are built from flavour CFFs read from `gpd_table.dat`:

- proton: (4/9)·(F₁ᵘ/2)·H_u + (1/9)·F₁ᵈ·H_d (and the same with F₂ for E), where F₁ᵘ, F₁ᵈ are
  the flavour form factors from the Dirac/Pauli form factors;
- neutron: isospin mirror, (4/9)·F₁ᵈ·H_d + (1/9)·(F₁ᵘ/2)·H_u;
- H̃ from the axial form factor (dipole), Ẽ from the pion pole.

| Table | Grid | Used by |
|---|---|---|
| Re/Im of the 8 flavour CFFs, 5 models | ξ only: 51 log-spaced points, **0.01 ≤ ξ ≤ 1** | `gpd_type` 0–4 (default 1, hard-coded) |
| Re/Im of H_u, H_d, 2 models | ξ (50 points) × \|t\| (21 points, 0.01–1.2 GeV²) | `gpd_type` 5, 6 (disabled) |

**Kinematic limits.** ξ = x_B/(2 − x_B) must satisfy 0.01 < ξ < 1, i.e. **x_B > 0.0198**.
Events outside are rejected and counted as "Out of xsec limits" (no extrapolation). At 10.6 GeV
with Q² > 1 GeV² the smallest x_B is 0.050, so the limit is not reached; it is reached for
Q²_min ≲ 0.5 GeV² at JLab energies or at higher beam energies. The CFFs have **no Q²
evolution** and the t dependence comes only from the analytic form factors. Models 5/6 have
known bugs (t interpolation index, out-of-range t not guarded) and must not be enabled as is.

### Bethe-Heitler collinear poles

With a massless electron, BH diverges when the real photon is collinear with the incoming
electron (P₁ → 0) or with the scattered electron (P₂ → 0), e.g. at φ = 0 for specific (x_B,
Q², t). Physically the electron mass makes the peak finite but extremely narrow. These photons
go down the beam pipe or merge with the electron; they are removed by the photon angle cuts
(see [Sampling](sampling.md)).

## Meson production (`--process 1`)

`src/ms_xsec.cc`: a parametrized π⁰/η cross-section (Q⁻⁴(Q²+M²)⁻² scaling, exponential t slope
of 1 GeV⁻², simple φ dependence). It has no virtual-photon flux / ε dependence and no
σ_TT, σ_LT terms, so it is a phase-space-like model. Mesons are decayed by JETSET (all
channels; π⁰ → γγ or Dalitz e⁺e⁻γ; η → γγ, 3π⁰, π⁺π⁻π⁰, …).

## Nuclear targets

| Target | `--targ_A --targ_Z` | Struck nucleon | Fermi momentum |
|---|---|---|---|
| p, n | 1 1, 1 0 | at rest | — |
| d | 2 1 | p or n (50/50) | AV18 deuteron momentum distribution |
| ⁴He | 4 4 | p or n | ⁴He distribution (up to 430 MeV/c) |
| ¹⁴N | 14 7 | p or n | C/O-interpolated parametrization (nucl-th/9507024) |
| NH₃ | 17 10 | H (3/17) or N (14/17) by nucleon count | as above |
| ND₃ | 20 10 | D (6/20) or N (14/20) | as above |

The direction of the Fermi momentum is isotropic. The struck nucleon and the spectator system
are on shell (no binding energy, no off-shell or flux correction), so the initial state does
not conserve energy exactly. The cross-section is the free-nucleon one evaluated at the
event kinematics. Coherent scattering and the polarized-target dilution (He bath, packing
fraction) are not modelled.
