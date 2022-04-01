#include "track_vars.h"

void init_track()
{
  trk.Struck_Nucl   = -1000;

  trk.xbj           = -1000.;
  trk.y             = -1000.;
  trk.Q2            = -1000.;
  trk.W2            = -1000.;
  trk.nu            = -1000.;
  trk.t             = -1000.;
  trk.ycol          = -1000.;

  trk.bh_xsec       = -1000.;
  trk.dvcs_xsec     = -1000.;
  trk.int_xsec      = -1000.;
  trk.tot_xsec      = -1000.;
  trk.ms_xsec       = -1000.;
  trk.kr            = -1000.;
  trk.isN           = -1000;

  for(int i=0; i<MAXT; i++)
  {
    trk.Type[i]    = -1000;
    trk.Charge[i]  = -1000;
    trk.Px[i]      = -1000.;
    trk.Py[i]      = -1000.;
    trk.Pz[i]      = -1000.;
    trk.P[i]       = -1000.;
    trk.E[i]       = -1000.;
    trk.Theta[i]   = -1000.;
    trk.Phi[i]     = -1000.;
  }
}
