#!/usr/bin/env python3
"""Statistical regression of GENEPI distributions (standard library only).

usage: check_distributions.py NAME FILE REF_JSON            compare FILE with REF_JSON[NAME]
       check_distributions.py NAME FILE REF_JSON --write    store FILE's statistics as REF_JSON[NAME]

For each observable the mean over the events is compared with the reference mean; the test
fails if they differ by more than NSIGMA standard errors of the tested sample. References
should be produced with much larger samples than the test (see tests/make_references.sh).
Observables: x_B, Q2, -t, <cos phi>, <cos 2phi> (Trento phi of the photon / photon pair),
fraction of photons within 2.5 deg of the beam (DVCS), and for bound nucleons
P(|cos theta(P_Fermi, q)| > 0.9) (= 0.1 for an isotropic Fermi motion).
"""
import json, math, os, sys

NSIGMA = 5.0
MP, EB = 0.938272, 10.6


def events(path):
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        n = int(lines[i].split()[0])
        yield int(lines[i].split()[8]), [l.split() for l in lines[i + 1:i + 1 + n]]
        i += 1 + n


def vec(p): return [float(p[6]), float(p[7]), float(p[8])]
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def norm(a): return math.sqrt(dot(a, a))


def observables(proc, parts):
    o = {}
    e = [p for p in parts if p[3] == "11" and p[0] == "1"][0]
    g = [p for p in parts if p[3] == "22"]
    if not g:
        return o
    pe = vec(e); Ee = float(e[9])
    qp = [sum(float(p[k]) for p in g) for k in (6, 7, 8)]; Eqp = sum(float(p[9]) for p in g)
    k = [0.0, 0.0, EB]; q = [k[j] - pe[j] for j in range(3)]; nu = EB - Ee
    Q2 = dot(q, q) - nu * nu
    o["xB"] = Q2 / (2 * MP * nu); o["Q2"] = Q2
    d = [q[j] - qp[j] for j in range(3)]
    o["mt"] = -((nu - Eqp) ** 2 - dot(d, d))
    nl, nh = cross(q, k), cross(q, qp)
    phi = math.atan2(dot(cross(nl, nh), q) / norm(q), dot(nl, nh))
    o["cosphi"] = math.cos(phi); o["cos2phi"] = math.cos(2 * phi)
    if proc == 0:   # DVCS only (a single photon in a meson event is a Dalitz decay)
        o["f_gbeam_2.5"] = 1.0 if math.degrees(math.acos(qp[2] / norm(qp))) < 2.5 else 0.0
    spec = [p for p in parts if p[3] in ("12", "1000010030", "1000020030")]
    if not spec:   # deuteron: spectator is the nucleon that is not the struck one (header col. 8)
        nuc = [p for p in parts if p[3] in ("2212", "2112")]
        if len(nuc) == 2:
            spec = [min(nuc, key=lambda p: norm(vec(p)))]
    if spec:
        pf = [-x for x in vec(spec[0])]
        if norm(pf) > 0:
            o["f_fermi_cos0.9"] = 1.0 if abs(dot(pf, q) / (norm(pf) * norm(q))) > 0.9 else 0.0
    return o


def stats(path):
    acc = {}
    for proc, parts in events(path):
        for key, v in observables(proc, parts).items():
            acc.setdefault(key, []).append(v)
    out = {}
    for key, v in acc.items():
        n = len(v); m = sum(v) / n
        sd = math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
        out[key] = {"n": n, "mean": m, "sd": sd}
    return out


def main():
    if len(sys.argv) < 4:
        print(__doc__); sys.exit(2)
    name, path, ref = sys.argv[1:4]
    s = stats(path)
    db = json.load(open(ref)) if os.path.exists(ref) else {}
    if len(sys.argv) > 4 and sys.argv[4] == "--write":
        db[name] = s
        json.dump(db, open(ref, "w"), indent=1, sort_keys=True)
        print(f"wrote reference '{name}' ({max(v['n'] for v in s.values())} events)")
        return
    if name not in db:
        print(f"no reference '{name}' in {ref}"); sys.exit(2)
    nbad = 0
    for key, r in sorted(db[name].items()):
        if key not in s:
            print(f"FAIL  {name} {key}: missing in test sample"); nbad += 1; continue
        t = s[key]
        err = max(t["sd"], r["sd"], 1e-12) / math.sqrt(t["n"])
        pull = (t["mean"] - r["mean"]) / err
        flag = "FAIL" if abs(pull) > NSIGMA else "ok  "
        nbad += abs(pull) > NSIGMA
        print(f"{flag}  {name:10s} {key:15s} test {t['mean']:+.4f} +- {err:.4f}   ref {r['mean']:+.4f}   pull {pull:+.1f}")
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
