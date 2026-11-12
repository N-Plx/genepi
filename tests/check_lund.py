#!/usr/bin/env python3
"""Consistency checks of a GENEPI LUND output file (standard library only).

usage: check_lund.py FILE --target {p,n,d,he4,n14,nh3,nd3} [--gamma-beam-min DEG]
                     [--gamma-e-min DEG] [--raster RX RY] [--vertex VX VY VZ] [--ebeam E]

Checks, per event:
  - header has 10 columns and the announced number of particle lines (14 columns each)
  - all numbers are finite; particle type is 1 (active)
  - E^2 - p^2 = m^2 for each particle
  - 3-momentum conservation: beam = sum of final particles (target at rest); for meson
    events the decay photons are used instead of the meson line, if both are present
  - PIDs belong to the allowed set; the struck nucleon (header col. 8) is p or n
  - header target A/Z is consistent with the target and with the remnant present
  - DVCS (process 0): photon angle to the beam and to the scattered electron >= cuts
  - vertex inside the raster ellipse
Exits 1 if any check fails, printing the first failures.
"""
import argparse, math, sys

ALLOWED = {11, -11, 22, 2212, 2112, 12, 111, 221, 211, -211, 1000010030, 1000020030}
HEAVY = {12: "N14 remnant", 1000010030: "triton", 1000020030: "He3"}
TARGET_AZ = {"p": {(1, 1)}, "n": {(1, 0)}, "d": {(2, 1)}, "he4": {(4, 2)}, "n14": {(14, 7)},
             "nh3": {(1, 1), (14, 7)}, "nd3": {(2, 1), (14, 7)}}


def angle(a, b):
    na, nb = math.sqrt(sum(x * x for x in a)), math.sqrt(sum(x * x for x in b))
    c = sum(x * y for x, y in zip(a, b)) / (na * nb)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file")
    ap.add_argument("--target", required=True, choices=sorted(TARGET_AZ))
    ap.add_argument("--gamma-beam-min", type=float, default=1.0)
    ap.add_argument("--gamma-e-min", type=float, default=1.0)
    ap.add_argument("--raster", type=float, nargs=2, default=(0.0, 0.0))
    ap.add_argument("--vertex", type=float, nargs=3, default=(0.0, 0.0, -3.0))
    ap.add_argument("--ebeam", type=float, default=10.6)
    a = ap.parse_args()

    fails, nev = [], 0
    def fail(ev, msg):
        fails.append(f"event {ev}: {msg}")

    lines = open(a.file).read().splitlines()
    i = 0
    while i < len(lines):
        h = lines[i].split()
        nev += 1
        if len(h) != 10:
            fail(nev, f"header has {len(h)} columns"); break
        n = int(h[0]); A, Z = int(h[1]), int(h[2]); eb = float(h[6]); struck = int(h[7]); proc = int(h[8])
        parts = [l.split() for l in lines[i + 1:i + 1 + n]]
        i += 1 + n
        if len(parts) != n or any(len(p) != 14 for p in parts):
            fail(nev, "wrong number of particle lines or columns"); continue
        try:
            vals = [[float(x) for x in p] for p in parts]
        except ValueError:
            fail(nev, "non-numeric value"); continue
        if not all(math.isfinite(x) for p in vals for x in p) or not math.isfinite(float(h[9])):
            fail(nev, "non-finite number"); continue
        pids = [int(p[3]) for p in vals]
        if abs(eb - a.ebeam) > 1e-6:
            fail(nev, f"beam energy {eb} != {a.ebeam}")
        if any(int(p[2]) != 1 for p in vals):
            fail(nev, "particle with type != 1")
        bad = set(pids) - ALLOWED
        if bad:
            fail(nev, f"unexpected PIDs {sorted(bad)}")
        if struck not in (2212, 2112):
            fail(nev, f"struck nucleon id {struck}")
        # E^2 - p^2 = m^2
        for p in vals:
            px, py, pz, E, m = p[6:11]
            if abs(E * E - (px * px + py * py + pz * pz) - m * m) > 1e-4 * E * E + 1e-4:
                fail(nev, f"pid {int(p[3])}: E^2-p^2 != m^2"); break
        # momentum conservation (target at rest)
        has_g = 22 in pids
        use = [p for p in vals if not (has_g and int(p[3]) in (111, 221))]
        s = [sum(p[k] for p in use) for k in (6, 7, 8)]
        if math.sqrt(s[0] ** 2 + s[1] ** 2 + (s[2] - eb) ** 2) > 1e-3:
            fail(nev, "3-momentum not conserved")
        # target A/Z and remnant
        if (A, Z) not in TARGET_AZ[a.target]:
            fail(nev, f"header A/Z = {A}/{Z} not allowed for target {a.target}")
        heavy = [x for x in pids if x in HEAVY]
        if (A, Z) == (14, 7) and heavy != [12]:
            fail(nev, "A=14 event without exactly one N14 remnant (pid 12)")
        if (A, Z) == (4, 2) and len([x for x in heavy if x in (1000010030, 1000020030)]) != 1:
            fail(nev, "He4 event without exactly one He3/triton")
        if (A, Z) in ((1, 1), (1, 0), (2, 1)) and heavy:
            fail(nev, f"A/Z = {A}/{Z} event with a heavy remnant")
        if (A, Z) == (2, 1) and sorted(x for x in pids if x in (2212, 2112)) != [2112, 2212]:
            fail(nev, "deuteron event without one proton and one neutron")
        # DVCS photon angle cuts
        if proc == 0:
            e = [p for p in vals if int(p[3]) == 11]
            g = [p for p in vals if int(p[3]) == 22]
            if len(e) != 1 or len(g) != 1:
                fail(nev, "DVCS event without exactly one electron and one photon")
            else:
                ab = angle(g[0][6:9], (0.0, 0.0, 1.0))
                ae = angle(g[0][6:9], e[0][6:9])
                if ab < a.gamma_beam_min - 1e-6 or ae < a.gamma_e_min - 1e-6:
                    fail(nev, f"photon at {ab:.3f} deg from beam / {ae:.3f} deg from e'")
        # vertex inside the raster ellipse, same for all particles
        vx, vy, vz = a.vertex; rx, ry = a.raster
        for p in vals:
            x, y, z = p[11:14]
            u = (((x - vx) / rx) ** 2 if rx > 0 else (0.0 if abs(x - vx) < 1e-6 else 2.0)) + \
                (((y - vy) / ry) ** 2 if ry > 0 else (0.0 if abs(y - vy) < 1e-6 else 2.0))
            if u > 1 + 1e-6 or abs(z - vz) > 1e-6:
                fail(nev, "vertex outside the raster ellipse"); break
        if len(fails) > 20:
            break

    if nev == 0:
        fails.append("file is empty")
    print(f"{a.file}: {nev} events checked, {len(fails)} failure(s)")
    for f in fails[:20]:
        print("  FAIL", f)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
