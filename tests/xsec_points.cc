// Regression test of the DVCS/BH cross-section code (src/dvcs.cc) at fixed kinematics.
//
// usage: xsec_points GPD_TABLE REF_FILE          compare with the reference (exit 1 on mismatch)
//        xsec_points GPD_TABLE REF_FILE --write  (re)write the reference after an intended change
//
// For each point: target (1 = proton, 0 = neutron), x_B, Q2, t, phi_Trento, beam and target
// helicities; values: BH, DVCS^2, interference (as entering the total) and total, in nb/GeV^4.
#include "genepi.h"
#include "dvcs_vars.h"
#include "inl_funcs.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>

struct Val { double bh, dvcs, inter, tot; };

// Same assembly of the helicity-dependent terms as in genepi.cc (process 0)
static bool eval(DVCS* dv, ReadOptFile* ro, NuclFF* ff, double x, double y, double Q2, double t, double phi, Val& v) {
  dv->init_dvcs();
  if (dv->phot_xsec(ro, ff, x, y, Q2, t, phi) != 0) return false;
  double bh = ro->get_fBheli(), th = ro->get_fTheli(), c = pi(1) - phi;
  v.bh = dv->hc0_BH + dv->hc1_BH * cos(c) + dv->hc2_BH * cos(2 * c) + bh * th * (dv->hc0_BH_LP + dv->hc1_BH_LP * cos(c));
  v.dvcs = dv->hc0_DVCS + dv->hc1_DVCS * cos(c) + bh * dv->hs1_DVCS * sin(c) +
           bh * th * (dv->hc0_DVCS_LP + dv->hc1_DVCS_LP * cos(c)) + th * dv->hs1_DVCS_LP * sin(c);
  double in = dv->hc0_Int + dv->hc1_Int * cos(c) + dv->hc2_Int * cos(2 * c) +
              bh * (dv->hs1_Int * sin(c) + dv->hs2_Int * sin(2 * c)) +
              bh * th * (dv->hc0_Int_LP + dv->hc1_Int_LP * cos(c) + dv->hc2_Int_LP * cos(2 * c)) +
              th * (dv->hs1_Int_LP * sin(c) + dv->hs2_Int_LP * sin(2 * c));
  v.inter = -ro->get_fBchg() * in;
  v.tot = v.bh + v.dvcs + v.inter;
  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) { fprintf(stderr, "usage: xsec_points GPD_TABLE REF_FILE [--write]\n"); return 2; }
  bool write = argc > 3 && strcmp(argv[3], "--write") == 0;
  ReadOptFile* ro = new ReadOptFile(); NuclFF* ff = new NuclFF(); DVCS* dv = new DVCS();
  ro->set_fMgpd(1); ro->set_fVol(3); ro->set_fBchg(-1); ro->set_fEb(10.6);
  dv->read_gpds(argv[1]);

  const double K[][3] = {{0.15, 1.5, -0.15}, {0.20, 2.0, -0.30}, {0.35, 3.0, -0.40}, {0.50, 5.0, -0.60}};
  const double PHI[] = {30., 90., 150., 210., 330.};   // away from the BH collinear poles at phi ~ 0
  std::vector<double> out;
  char line[512];
  std::vector<std::string> rows;
  for (int ipn = 1; ipn >= 0; ipn--) {
    ro->set_fIpn(ipn);
    double M = m_targ(ipn, 1.);
    for (auto& k : K) {
      double x = k[0], Q2 = k[1], t = k[2], y = Q2 / (2 * M * 10.6 * x);
      for (double phd : PHI)
        for (int bh = -1; bh <= 1; bh += 2)
          for (int th = -1; th <= 1; th++) {
            ro->set_fBheli(bh); ro->set_fTheli(th);
            Val v;
            if (!eval(dv, ro, ff, x, y, Q2, t, phd * pi(1) / 180., v)) { fprintf(stderr, "point outside limits\n"); return 2; }
            snprintf(line, sizeof line, "%d %.2f %.2f %.2f %.0f %+d %+d %.10e %.10e %.10e %.10e",
                     ipn, x, Q2, t, phd, bh, th, v.bh, v.dvcs, v.inter, v.tot);
            rows.push_back(line);
          }
    }
  }

  if (write) {
    FILE* f = fopen(argv[2], "w");
    fprintf(f, "# ipn x Q2 t phi bheli theli BH DVCS INT TOT  [nb/GeV^4], E=10.6 GeV, e- beam, gpd_type 1\n");
    for (auto& r : rows) fprintf(f, "%s\n", r.c_str());
    fclose(f);
    printf("wrote %zu reference points to %s\n", rows.size(), argv[2]);
    return 0;
  }

  FILE* f = fopen(argv[2], "r");
  if (!f) { fprintf(stderr, "cannot open reference %s\n", argv[2]); return 2; }
  int nbad = 0; size_t n = 0;
  while (fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    if (n >= rows.size()) { nbad++; break; }
    int ipn, bh, th; double x, Q2, t, phi, r[4], c[4];
    sscanf(line, "%d %lf %lf %lf %lf %d %d %lf %lf %lf %lf", &ipn, &x, &Q2, &t, &phi, &bh, &th, &r[0], &r[1], &r[2], &r[3]);
    sscanf(rows[n].c_str(), "%*d %*f %*f %*f %*f %*d %*d %lf %lf %lf %lf", &c[0], &c[1], &c[2], &c[3]);
    const char* name[4] = {"BH", "DVCS", "INT", "TOT"};
    for (int i = 0; i < 4; i++)
      if (fabs(c[i] - r[i]) > 1e-6 * fabs(r[i]) + 1e-12) {
        if (nbad < 20) printf("MISMATCH %s ipn=%d x=%.2f Q2=%.2f t=%.2f phi=%.0f bheli=%+d theli=%+d: %.6e (ref %.6e)\n",
                              name[i], ipn, x, Q2, t, phi, bh, th, c[i], r[i]);
        nbad++;
      }
    n++;
  }
  fclose(f);
  if (n != rows.size()) { printf("reference has %zu points, expected %zu\n", n, rows.size()); nbad++; }
  printf("cross-section regression: %zu points, %d mismatch(es)\n", rows.size(), nbad);
  return nbad ? 1 : 0;
}
