#define MAXT 50

typedef struct TRACK{
  int    Evt_ID;
  int    Ntracks;
  double Eb;
  int    Bheli;
  int    TarA;
  int    TarZ;
  int    Theli;
  int    Struck_Nucl;
  int    Process;
  int    ExcMS;

  double xbj;
  double y;
  double Q2;
  double W2;
  double nu;
  double t;
  double ycol;

  double bh_xsec;
  double dvcs_xsec;
  double int_xsec;
  double tot_xsec;
  double ms_xsec;
  double kr;
  int    isN;

  int    Type[MAXT];
  int    Charge[MAXT];
  double Px[MAXT];
  double Py[MAXT];
  double Pz[MAXT];
  double P[MAXT];
  double E[MAXT];
  double Theta[MAXT];
  double Phi[MAXT];
};
extern TRACK trk;
