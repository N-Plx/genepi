# Output format (LUND)

One header line per event followed by one line per particle.

**Header** (10 columns):

| # | Content |
|---|---|
| 1 | number of particles |
| 2, 3 | target A, Z of the nucleus hit in this event (1/1 H, 1/0 n, 2/1 D, 4/2 ⁴He, 14/7 ¹⁴N) |
| 4, 5 | target and beam polarization (0) |
| 6 | beam particle (11) |
| 7 | beam energy [GeV] |
| 8 | struck nucleon (2212 or 2112) |
| 9 | process (0 DVCS, 1 meson) |
| 10 | dσ of the event [nb/GeV⁴] (without the t Jacobian) |

**Particles** (14 columns): index, lifetime (−1), type (1 = active), PID, parent (0),
daughter (0), px, py, pz, E, m [GeV], vx, vy, vz [cm].

## Particles written

| Process / target | Particles |
|---|---|
| DVCS, p or n | e⁻, recoil nucleon, γ |
| DVCS, d | e⁻, spectator nucleon, recoil nucleon, γ |
| DVCS, ⁴He | e⁻, spectator ³He (1000020030, struck n) or triton (1000010030, struck p), recoil nucleon, γ |
| DVCS, ¹⁴N (and N in NH₃/ND₃) | e⁻, remnant written as **PID 12** (see below), recoil nucleon, γ |
| mesons | as above with the meson decay products from JETSET instead of the γ |

### The ¹⁴N remnant (PID 12)

The ¹³C/¹³N remnant is written with the ν_e code 12, active. gemc cannot track a generic ion
from a LUND file: an active line with the real code (1000060130 / 1000070130) makes gemc exit
("Particle id … not found in G4 table"), and gemc's special code 80000 aborts. As a neutrino the
remnant appears in MC::Particle without interacting, and nitrogen events can be identified by
PID 12 in MC::Particle or by `atarget = 14` in MC::Event. The remnant mass (column 11) is the
sum of the nucleon masses.

## gemc and coatjava

Tested with gemc 5.14 / 6.1 and coatjava 10.0.9, 11.1.0, 14.1.3:

- gemc copies the header into MC::Event (atarget, ztarget, targetid, processid, weight) and
  all LUND lines into MC::Lund. MC::Particle holds the tracked particles; gemc writes the light
  ions as 45 (d), 46 (t), 47 (α), 49 (³He) there, whatever code the generator used, while
  MC::True keeps full PDG codes.
- **⁴He with truth matching**: coatjava's TruthMatch looks every MC::Particle PID up in
  PDGDatabase, which (up to at least 14.2) does not know 46/47/49: the engine fails on those
  events and MC::GenMatch/MC::RecMatch stay empty. Fixed in coatjava by adding the gemc codes
  to PDGDatabase and guarding the lookup; until that is released, ⁴He samples reconstructed
  with truth matching lose the match banks (the rest of the reconstruction is unaffected).
