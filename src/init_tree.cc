#include "cpp.h"
#include "root.h"
#include "track_vars.h"

void init_tree(TTree *tr)
{
  cout<<"init_tree beam heliciy = "<<trk.Bheli<<endl;

  tr->Branch("Evt_ID",       &trk.Evt_ID,       "Evt_ID/I");
  tr->Branch("Ntracks",      &trk.Ntracks,      "Ntracks/I");
  tr->Branch("Eb",           &trk.Eb,           "Eb/D");
  tr->Branch("Bheli",        &trk.Bheli,        "Bheli/I");
  tr->Branch("TarA",         &trk.TarA,         "TarA/I");
  tr->Branch("TarZ",         &trk.TarZ,         "TarZ/I");
  tr->Branch("Theli",        &trk.Theli,        "Theli/I");
  tr->Branch("Struck_Nucl",  &trk.Struck_Nucl,  "Struck_Nucl/I");
  tr->Branch("Process",      &trk.Process,      "Process/I");
  tr->Branch("ExcMS",        &trk.ExcMS,        "ExcMS/I");
  tr->Branch("xbj",          &trk.xbj,          "xbj/D");
  tr->Branch("y",            &trk.y,            "y/D");
  tr->Branch("Q2",           &trk.Q2,           "Q2/D");
  tr->Branch("W2",           &trk.W2,           "W2/D");
  tr->Branch("nu",           &trk.nu,           "nu/D");
  tr->Branch("t",            &trk.t,            "t/D");
  tr->Branch("ycol",         &trk.ycol,         "ycol/D");

  tr->Branch("bh_xsec",      &trk.bh_xsec,       "bh_xsec/D");
  tr->Branch("dvcs_xsec",    &trk.dvcs_xsec,     "dvcs_xsec/D");
  tr->Branch("int_xsec",     &trk.int_xsec,      "int_xsec/D");
  tr->Branch("tot_xsec",     &trk.tot_xsec,      "tot_xsec/D");
  tr->Branch("ms_xsec",      &trk.ms_xsec,       "ms_xsec/D");
  tr->Branch("keep_reject",  &trk.kr,            "kr/D");
  tr->Branch("isN",          &trk.isN,           "isN/I");

  tr->Branch("Type",        trk.Type,         "Type[Ntracks]/I");
  tr->Branch("Charge",      trk.Charge,       "Charge[Ntracks]/I");
  tr->Branch("Px",          trk.Px,           "Px[Ntracks]/D");
  tr->Branch("Py",          trk.Py,           "Py[Ntracks]/D");
  tr->Branch("Pz",          trk.Pz,           "Pz[Ntracks]/D");
  tr->Branch("P",           trk.P,            "P[Ntracks]/D");
  tr->Branch("E",           trk.E,            "E[Ntracks]/D");
  tr->Branch("Theta",       trk.Theta,        "Theta[Ntracks]/D");
  tr->Branch("Phi",         trk.Phi,          "Phi[Ntracks]/D");
}
