/*
  Author Ahmed El Alaoui, LPSC (Grenoble)
  Ref: BMK, hep-ph/0112108

  Additional features:
  Fermi motion on N : Mathieu Erhart, IJCLab
  NH3/ND3 targets : Noémie Pilleux, IJCLab
*/

#include "genepi.h"
#include "inl_funcs.h"
#include "lujets_cc.h"
#include "jetset.h"
#include "hepevt.h"

#include "track_vars.h"
#include "dvcs_vars.h"
#include "vect.h"

#include "root.h"

//event_vars.h
//beam 0 =============================================================
vector<double> vk;
double k;
double Ee;

//scat. elec. 1 ======================================================
vector<double> vkp;
double kp;
double Eep;
double Thetakkp ,Phikkp;

//q(virtual photon) 2 ================================================
vector<double> vq;
double q;
double nu;
double Thetakq, Phikq;

vector<double> vq2;

//initial nucleon P(Ep;P) 3 ==========================================
vector<double> vP;
double P;
double E;
double ThetakP ,PhikP;

//initial nucleon P(Ep;P) 3 ==========================================
vector<double> vP1;
double ThetaqP, PhiqP;

//target remnant Pr(EPr,Pr) 4 ========================================
vector<double> vPr;
double Pr;
double Er;
double ThetakPr, PhikPr; 

//final nucleon P'(EP',P') wrt virt. phot. frame =====================
vector<double> vPp1;
double Pp;
double Ep;
double ThetaqPp, PhiqPp; 

//final nucleon P'(EP',P') wrt inc. elec. frame 5 ====================
vector<double> vPp;
double ThetakPp, PhikPp;

//qp(real photon) wrt virt. phot. frame ==============================
vector<double> vqp1;
double qp;
double nup;
double Thetaqqp, Phiqqp;

//qp(real photon) wrt inc. elec. frame 6 =============================
vector<double> vqp;
double Thetakqp, Phikqp;

//Pm(meson) wrt inc. elec. frame 6 =============================
vector<double> vPm;
double Pm;
double Em;
double ThetakPm, PhikPm;

//polar and azimuthal angles between scat. elec. and real photon =====
double Thetakpqp, Phikpqp;

//kinematics =========================================================
double Q2;   //photon virtuality
double xbj;  //xbjorken variable
double y;    //elec. energy fraction nu/E
double W2;   //invariant mass of proton + virtual photon system
double t;
double ycol;

//Cross Section ======================================================
double ds_bh;
double ds_dvcs;
double ds_int;
double ds_tot;
double ds_ms;
double kr;
int isN;

HEPEVT hepevt;
TRACK trk;

char fdump[1000];
FILE *ptr;

int main(int argc, char*argv[])
{

  if(argc != 2)
  {
    cout<<"Error"<<endl;
    cout<<"Usage:"<<endl;
    cout<<"./genepi.exe <Option_file>"<<endl;
    return 0;
  }

  //============Handling input file and output location======================//

  ReadOptFile *ro = new ReadOptFile();
  DVCS *dv        = new DVCS();
  NuclFF *ff      = new NuclFF();
  VECT *vect      = new VECT();

  string InputFile = argv[1];
  if(fexist(InputFile.c_str()))
  {
    cout<<"File "<<InputFile.c_str()<<" does not exist"<<endl;
    return 0;
  }

  if(ro->ReadInputFile(InputFile.c_str()) == -1)
  {
    cout<<"Unable to read input file"<<endl;
    return 0;
  }

  string fDir(ro->get_fDir());
  if(ro->get_fDate() == 1)
  {
    fDir += get_date();
    cout<<"output directory : "<<fDir<<endl;
  }

  prod_dir(fDir);

  string cmd("cp "+InputFile+" "+fDir);
  system(cmd.c_str());

  string fout0;
  fout0 = fDir+"/summary.dat";
  ofstream out0(fout0.c_str());
  out0<<"This is the summary file"<<endl;

  int seed;
  if(ro->get_fSeed() == 1)
  {
    seed = get_seed();
    cout<<"SEED: "<< seed <<endl;
    out0<<"SEED: "<< seed <<endl;
  }
  else
  {
    seed = ro->get_fSeed();
  }

  //=================Seed for generation of kin variables=========================//
  TRandom1 rndm; 
  rndm.SetSeed(seed);
  gRandom->SetSeed(seed);

  //Read input A and Z
  int iApZ = ro->get_fAt() + ro->get_fZt();
  //We store this value in case of a molecular target where we will change iApZ when randomly selecting the target
  int iApZ_stored = iApZ; 
  //Some counters for testing/debugging
  int count_N = 0;
  int count_p = 0;
  int counter_kept = 0;
  int counter_rejected = 0;

  //Input run number
  int irunnum = ro->get_fRunnum();

  //Selecting the target
  string target;
  if(iApZ == 1)
  {
    target = "neut";
  }
  else if(iApZ == 2)
  {
    target = "prot";
  }
  else if(iApZ == 3)
  {
    target = "deut";
  }
  else if(iApZ == 8)
  {
    target = "hel4";
  }
  else if(iApZ == 21) //A+Z=14+7
    {
      target = "nitro";
    }
  else if(iApZ == 27) //NH3 A+Z = 17 + 10
    {
      target = "NH3";
    }
  else if(iApZ == 30) //ND3 A+Z = 20 + 10                                                                                                                                                                         
    {
      target = "ND3";
    }

  out0<<"Target:  "<< target<<endl;

  //Needed for root file output
  int ntp_cnt(0);
  char root_file[1000];
  TFile *f;
  TTree *tree;

  if(ro->get_fNtp())
  {
    sprintf(root_file, "%s/ntup_%s_%3.2fgev_%d.root", 
            fDir.c_str(), target.c_str(), ro->get_fEb(), irunnum);
    f    = new TFile(root_file,"RECREATE");
    tree = new TTree("DVCS", "/dvcs_tree");
    init_tree(tree);
    //TTree::SetMaxTreeSize(200000000);
    cout<<"Filling ntuple "<<root_file<<endl;
    out0<<"Filling ntuple "<<root_file<<endl;
  }

  //read GPDs from gpd_table.dat
  string gpd_tbl=("gpd_table.dat");
  if(fexist(gpd_tbl.c_str()))
  {
    cout<<"File "<<gpd_tbl.c_str()<<" does not exist"<<endl;
    return 0;
  }

  dv->read_gpds(gpd_tbl);

  //Kinematic variables initialization
  double xmin,   xmax;
  double numin,  numax;
  double Q2min,  Q2max;
  double tmin,   tmax;
  double nupmin, nupmax;
  double Epmin,  Epmax;

  double Xmin(1000.),   Xmax(-1000.);
  double Ymin(1000.),   Ymax(-1000.);
  double NUmin(1000.),  NUmax(-1000.);
  double QQ2min(1000.), QQ2max(-1000.);
  double WW2min(1000.), WW2max(-1000.);
  double Tmin(1000.),   Tmax(-1000.);
  double NUPmin(1000.), NUPmax(-1000.);
  double EPmin(1000.),  EPmax(-1000.);

  int    Nhit_neut(0);
  int    Nhit_prot(0);
  int    run(0);
  int    Nevts_per_Ntup(0);
  int    outside_klim(0);
  int    outside_klim_t(0);
  int    outside_klim_nu(0);
  int    outside_klim_t2(0);
  int    outside_klim_x(0);
  int    outside_klim_x2(0);
  int    outside_klim_W2(0);

  //Recoil, targer
  double M_RECO(-1000.), M_RECO2(-1000.), M_TARG, M_TARG2;
  //Ipn: proton or neutron, Ims: identify meson
  int    Ipn, Ims, RECO_ID(-1000), RECO_CH(-1000);
  double cosThetakkp, cosThetaqqp, cosThetakpqp;
  //zero 3-vector 
  vector<double> v0(3,0);

  //output variables
  lujets_cc.N = 0;
  lujets_.n   = 0;
  for(int ii=0; ii<5; ii++)
  {
    for(int ij=0; ij<4000; ij++)
    {
      lujets_cc.K[ii][ij] = 0;
      lujets_cc.P[ii][ij] = 0;
      lujets_cc.V[ii][ij] = 0;

      lujets_.k[ii][ij] = 0;
      lujets_.p[ii][ij] = 0;
      lujets_.v[ii][ij] = 0;
    }
  }

  hepevt.NEVHEP = 0;
  hepevt.NHEP   = 0;

  for(int ii=0; ii<NMXHEP; ii++)
  {
    hepevt.ISTHEP[ii] = 0;
    hepevt.IDHEP[ii]  = 0;
    for(int jj=0; jj<2; jj++) hepevt.JMOHEP[ii][jj] = 0;
    for(int jj=0; jj<2; jj++) hepevt.JDAHEP[ii][jj] = 0;
    for(int jj=0; jj<5; jj++) hepevt.PHEP[ii][jj]   = 0;
    for(int jj=0; jj<4; jj++) hepevt.VHEP[ii][jj]   = 0;
  }

  //=====================Start event loop===================//
  for(long long ievt=0; ievt<ro->get_fNevts(); ievt++)
  {
    //Get event info from input
    trk.Evt_ID  = ievt;
    trk.Bheli   = ro->get_fBheli();
    trk.Theli   = ro->get_fTheli();
    trk.Process = ro->get_fProc();
    trk.ExcMS   = ro->get_fIms();
    Nevts_per_Ntup++;

    if(ievt%ro->get_fPrint() == 0)
    {
      cout<<"Processing Event: "<<ievt<<"/"<<ro->get_fNevts()<<endl;
      out0<<"Processing Event: "<<ievt<<"/"<<ro->get_fNevts()<<endl;
    }

    if(ro->get_fAscii() != 0 && ievt%ro->get_fNevtsPerFile() == 0)
    {
      if(run<10)
      {
        sprintf(fdump, "%s/events_%s_run%d.dat", 
	        fDir.c_str(), target.c_str(), irunnum);
      }
      else if(run<100)
      {
        sprintf(fdump, "%s/events_%s_run%d.dat", 
	        fDir.c_str(), target.c_str(), irunnum);
      }
      else
      {
        sprintf(fdump, "%s/events_%s_run%d.dat", 
	        fDir.c_str(), target.c_str(), irunnum);
      }
      run++;
      ptr = fopen(fdump,"w+");
    }

    //Inititalization
    dv->init_dvcs();
    init_event();
    init_track();

    //Beam electron
    vk   = vect->init_vect();
    //Scattered electron
    vkp  = vect->init_vect();
    //Virtual photon
    vq   = vect->init_vect();
    vq2  = vect->init_vect();
    //initial nucleon
    vP   = vect->init_vect();
    vP1  = vect->init_vect();
    //Target remnant
    vPr  = vect->init_vect();
    //Real photon, inc.elec. frame
    vqp  = vect->init_vect();
    //Real photon, virt.phot. frame
    vqp1 = vect->init_vect();
    //Final nucleon, inc.elec. frame
    vPp  = vect->init_vect();
    //Final nucleon, virt.phot. frame
    vPp1 = vect->init_vect();

    //In case of a molecular target, we select the nuclei with which to interact
    if(iApZ_stored == 27) //NH3 : 3/17 chances to interact with a p from H3                                                                                                                                   
      {                                                                                                                                                       
	double rndmtarg  = 17*rndm.Rndm();
	if(rndmtarg < 3)
	  {
	    isN=0;
	    count_p++;
	    iApZ = 2; //we change the value of iApZ to consider the right target for fermi motion                                                                                                                 
	  }
	else
	  {
	    isN=1;
	    count_N++;
	    iApZ = 21;                                                                                                                    
	  }
      }
    
    if(iApZ_stored == 30) //ND3 6/20 chances to interact with deuterium                                                                                                                                                  
      {
	double rndmtarg  = 20*rndm.Rndm();
	if(rndmtarg < 6)
	  {
	    isN=0;
	    count_p++;
	    iApZ=3;
	  }
	else
	  {
	    isN=1;
	    count_N++;
	    iApZ = 21;
	  }
      }


    if(iApZ == 1) //neutron
    {
      Ipn     = 0;
      RECO_ID = -1000;
      RECO_CH = -1000;
      M_RECO  = -1000;
    }
    else if(iApZ == 2) //proton
    {
      Ipn     = 1;
      RECO_ID = -1000;
      RECO_CH = -1000;
      M_RECO  = -1000;
    }
    else if(iApZ == 3) //deuteron
    {
      //select which nucleon to interact with in case of nuclei target
      double rndmnucl  = rndm.Rndm();
      if(rndmnucl >= 0.5)
      {
        Ipn     = 0;
        RECO_ID = prot_id();
        RECO_CH = 1;
	M_RECO  = m_prot(1);
      }
      else if(rndmnucl < 0.5)
      {
        Ipn     = 1;
        RECO_ID = neut_id();
        RECO_CH = 0;
	M_RECO  = m_neut(1);
      }
    }
    else if(iApZ == 8) //He4
    {
      //select which nucleon to interact with in case of nuclei target
      double rndmnucl  = rndm.Rndm();
      if(rndmnucl >= 0.5)
      {
        Ipn     = 0;
        RECO_ID = hel3_id();
        RECO_CH = 2;
	M_RECO  = 2.*m_prot(1) + m_neut(1); //add binding energy
      }
      else if(rndmnucl < 0.5)
      {
        Ipn     = 1;
        RECO_ID = trit_id();
        RECO_CH = 1;
	M_RECO  = 2.*m_neut(1) + m_prot(1); //add binding energy
      }
    }
    else if(iApZ == 21) //Nitro     
      {       
	//select which nucleon to interact with in case of nuclei target
	double rndmnucl  = rndm.Rndm();
	if(rndmnucl >= 0.5)
	  {
	    Ipn     = 0;
	    RECO_ID = remnant_id();
	    RECO_CH = 2;
	    M_RECO  = 7.*m_prot(1) + 6.*m_neut(1); //add binding energy
	  }       
	else if(rndmnucl < 0.5)
	  {
	    Ipn     = 1;
	    RECO_ID = remnant_id();
	    RECO_CH = 1;
	    M_RECO  = 6.*m_neut(1) + 7.*m_prot(1); //add binding energy
	  }     
      }

    //proton or neutron
    trk.Struck_Nucl = Ipn;
    trk.isN=isN;
    M_TARG          = m_targ(Ipn,1.);
    M_TARG2         = m_targ(Ipn,2.);
    if(iApZ > 2 ) M_RECO2 = sqr(M_RECO); 
    ro->set_fIpn(Ipn);
    Ims = ro->get_fIms();

    //Kinem limits computed from input file
    get_kinem_limits(ro);

    //generate Q2 and Xbj
    double rndmQ2, rndmXbj;
    rndmQ2  = rndm.Rndm();
    Q2      = ro->get_fQ2min() + (ro->get_fQ2max() - ro->get_fQ2min())*rndmQ2;
    rndmXbj = rndm.Rndm();
    xbj     = ro->get_fXbjmin() + (ro->get_fXbjmax() - ro->get_fXbjmin())*rndmXbj;
    nu      = Q2/(2.*M_TARG*xbj);
    y       = nu/ro->get_fEb();
    W2      = M_TARG2 - Q2 + Q2/xbj;
    q       = sqrt(Q2+sqr(nu));

    numin   = ro->get_fNumin();
    numax   = ro->get_fEb() - Q2/(4.*ro->get_fEb()); 
    numax   = min(numax, ro->get_fNumax());
   
    if(nu > numax  || nu < numin)
    {
      //cout<<"Event = "<<ievt<<
      //  "; nu out of range "<<numin<<" "<<nu<<" "<<numax<<endl;
      //out0<<"Event = "<<ievt<<
      //  "; nu out of range "<<numin<<" "<<nu<<" "<<numax<<endl;

      outside_klim++;
      outside_klim_nu++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

    Q2min   = ro->get_fQ2min();
    Q2max   = 4.*ro->get_fEb()*(ro->get_fEb() -  nu);
    Q2max   = min(Q2max, ro->get_fQ2max());
   
    if(Q2 > Q2max || Q2 < Q2min)
    {
      cout<<"Event = "<<ievt<<"; Q2 out of range "<<
        Q2min<<" "<<Q2<<" "<<Q2max<<endl;
      out0<<"Event = "<<ievt<<"; Q2 out of range "<<
        Q2min<<" "<<Q2<<" "<<Q2max<<endl;

      outside_klim++;
      //outside_klim_Q2++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

    xmin  = (2.*ro->get_fEb()*Q2)/(M_TARG*(4.*sqr(ro->get_fEb()) - Q2));
    xmin  = max(xmin, ro->get_fXbjmin());
    xmax  = ro->get_fEb()*(q - nu)/(M_TARG*nu);
    xmax  = min(xmax, ro->get_fXbjmax());
   
    if(xbj > xmax || xbj < xmin)
    {
      //cout<<"event "<<ievt<<"; xbj out of range "<<
      //  xmin<<" "<<xbj<<" "<<xmax<<endl;
      //out0<<"event "<<ievt<<"; xbj out of range "<<
      //  xmin<<" "<<xbj<<" "<<xmax<<endl;

      outside_klim++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

    //inc. elec.
    k  = sqrt(sqr(ro->get_fEb()) - m_elec(2));
    Ee = ro->get_fEb();
    vk = vect->fill_vect(k, 0., 0.);
    
    //Target nucleon
    P = 0.0;
    E = m_targ(Ipn,1);
    vP = vect->fill_vect(P, 0., 0.);

    //polar angle between inc. and scat. elec.
    cosThetakkp = 1. - Q2/(2.*Ee*(Ee-nu));
    if(abs(cosThetakkp) >= 1.)
    {
      cout<<"event "<<ievt<<
        "; cosThetakkp larger than one; cosThetakkp = "<<cosThetakkp<< endl;
      cosThetakkp = sign(cosThetakkp);
    }

    Thetakkp =  acos(cosThetakkp);//[0,PI]

    if(Thetakkp > pi(1)) 
      cout<<"Thetakkp larger than pi; Thetakkp = "<<cosThetakkp<< endl;

    //generate the scat. elect. azimuthal angle [0,2PI]
    Phikkp = 2.*pi(1)*rndm.Rndm();

    //scat. elec. energy and momentum wrt inc. elec. frame
    Eep = Ee - nu;
    kp  = sqrt(sqr(Eep) - m_elec(2));
    vkp = vect->fill_vect(kp, Thetakkp, Phikkp);

    //virtual photon 3-momentum vector 
    vq = vect->diff_vect(vk, vkp);

    //polar angle of virt. photon
    Thetakq     = get_theta(vq.at(2), q);
    
    //azimuthal angle of virt. photon
    Phikq       = get_phi(vq.at(0), vq.at(1));

    //upper and lower limits of nup, Ep and t
    nupmin = (m_targ(Ipn,1)*nu - Q2/2.)/(m_targ(Ipn,1) + nu + q);
    nupmax = (m_targ(Ipn,1)*nu - Q2/2.)/(m_targ(Ipn,1) + nu - q);

    Epmin = m_targ(Ipn,1) + nu - nupmax;
    Epmax = m_targ(Ipn,1) + nu - nupmin;

    tmin   = 2.*m_targ(Ipn,1)*(nupmin - nu);
    tmin   = max(tmin, ro->get_ftmin());
    tmax   = 2.*m_targ(Ipn,1)*(nupmax - nu);
    tmax   = min(tmax, ro->get_ftmax());
   
    //generate t
    double rndmt;
    rndmt = rndm.Rndm();
    t     = tmin + (tmax - tmin)*rndmt;
    if(t > tmax || t < tmin)
    {
      //cout<<"event "<<ievt<<"; t out of range "<<tmin<<" "<<t<<" "<<tmax<<endl;
      //out0<<"event "<<ievt<<"; t out of range "<<tmin<<" "<<t<<" "<<tmax<<endl;

      outside_klim++;
      outside_klim_t++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

    //real photon energy and momentum
    nup = nu + t/(2.*m_targ(Ipn,1));
    qp  = nup;
    if(nup > nupmax || nup < nupmin)
    {
      cout<<"event "<<ievt<<"; nup out of range "<<
        nupmin<<" "<<nup<<" "<<nupmax<<endl;
      out0<<"event "<<ievt<<"; nup out of range "<<
        nupmin<<" "<<nup<<" "<<nupmax<<endl;

      outside_klim++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

    //final nucleon energy and momentum
    Ep = m_targ(Ipn,1) - t/(2.*m_targ(Ipn,1));
    Pp = sqrt(sqr(Ep) - m_targ(Ipn,2));
    if(Ep > Epmax || Ep < Epmin)
    {
      cout<<"event "<<ievt<<"; Ep out of range "<<
        Epmin<<" "<<Ep<<" "<<Epmax<<endl;
      out0<<"event "<<ievt<<"; Ep out of range "<<
        Epmin<<" "<<Ep<<" "<<Epmax<<endl;

      outside_klim++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

    //polar angle of real photon wrt virt. phot. frame
    cosThetaqqp = (t + 2.*nu*nup + Q2)/(2.*q*nup);
    if(abs(cosThetaqqp) >= 1.)
    {
      cout<<"event "<<ievt<<
        "; cosThetaqqp larger than one; cosThetaqqp = "<<cosThetaqqp<< endl;
      cosThetaqqp = sign(cosThetaqqp);
    }

    Thetaqqp =  acos(cosThetaqqp);//[0,PI]

    if(Thetaqqp > pi(1)) 
      cout<<"Thetaqqp larger than pi; Thetaqqp = "<<cosThetaqqp<< endl;

    //generate real photon azimuthal angle wrt virt. phot. frame
    double rndmphi = rndm.Rndm();
    Phiqqp = 2.*pi(1)*rndmphi;

    //real photon 3-momentum vector wrt virt. phot. frame
    vqp1 = vect->fill_vect(qp, Thetaqqp, Phiqqp);

    //virt. photon 3-momentum vector wrt virt. phot. frame
    vq2 = vect->fill_vect(q, 0., 0.);

    //final nucleon 3-momentum vector wrt virt. phot. frame
    vPp1 = vect->diff_vect(vq2, vqp1);

    //polar angle of final nucleon wrt virt. phot. frame
    ThetaqPp     = get_theta(vPp1.at(2), Pp);

    //azimuthal angle of final nucleon wrt virt. phot. frame
    PhiqPp       = get_phi(vPp1.at(0), vPp1.at(1));

    //rotate the real photon 3-momentum vector from 
    //the virt. phot. frame to the inc. elec. frame
    vqp = vect->invrot_vect(vqp1, Thetakq, Phikq);

    //polar angle of real photon wrt inc. elec. frame
    Thetakqp     = get_theta(vqp.at(2), qp);

    //azimuthal angle of real photon wrt inc. elec. frame
    Phikqp   = get_phi(vqp.at(0), vqp.at(1));

    //final nucleon 3-momentum vector wrt inc. elec. frame
    vPp  = vect->diff_vect(vect->add_vect(vP,vq), vqp);

    //polar angle of final nucleon wrt inc. elec. frame
    ThetakPp     = get_theta(vPp.at(2), Pp);

    //azimuthal angle of final nucleon wrt inc. elec. frame
    PhikPp       = get_phi(vPp.at(0), vPp.at(1));


    //Angle between scat. elec. and real photon
    cosThetakpqp = cos(Phikqp)*sin(Thetakqp)*cos(Phikkp)*sin(Thetakkp) + 
                   sin(Phikqp)*sin(Thetakqp)*sin(Phikkp)*sin(Thetakkp) +
                   cos(Thetakqp)*cos(Thetakkp);
    Thetakpqp    = acos(cosThetakpqp);

    //======================START FERMI MOTION PART (for incoherent scattering off nucleus target)
    if(iApZ == 3 || iApZ == 8 || iApZ == 21)
    {
      //generate initial nucleon momentum wrt the inc. elec. frame
      int ifermi(1);
      double Output[7];
      double rnd = rndm.Rndm();

      //Call fermi motion corresponding to the nuclei which was hit
      if(iApZ == 3)
      {
        P = fer_mom_deut(ifermi, rnd);
      }
      else if(iApZ == 8)
      {
        P = fer_mom_hel4(rnd);
      }
      else if(iApZ == 21)       
	{
	  P = fer_mom_nitro(rndm);
	}
      E = sqrt(sqr(P) + M_TARG2);
      
      //generate polar angle of initial nucleon wrt the inc. elec. frame
      ThetakP = pi(1)*rndm.Rndm();
      
      //generate azimuthal angle of initial nucleon wrt the inc. elec. frame
      PhikP = 2.*pi(1)*rndm.Rndm();

      //initial nucleon 3-momentum vector wrt inc. elec. frame
      vP = vect->fill_vect(P, ThetakP, PhikP);

      //recoil nucleus 3-momentum vector wrt inc. elec. frame
      Pr  = P;
      Er  = sqrt(sqr(P) + M_RECO2);
      vPr = vect->diff_vect(v0, vP);

      //recoil nucleus polar angle wrt the inc. elec. frame
      ThetakPr     = get_theta(vPr.at(2),Pr);

      //recoil nucleus azimuthal angle wrt the inc. elec. frame
      PhikPr       = get_phi(vPr.at(0),vPr.at(1));

      //rotate the initial nucleon 3-momentum vector from 
      //the inc. elec. frame to the virt. phot. frame
      vP1 = vect->rot_vect(vP, Thetakq, Phikq);

      //polar angle of initial nucleon in virt. phot. frame
      ThetaqP = get_theta(vP1.at(2),P);

      //azimuthal angle of initial nucleon in virt. phot. frame
      PhiqP = get_phi(vP1.at(0),vP1.at(1));

      double Pq = E*nu - P*q*cos(ThetaqP);
      double Pk = E*Ee - P*k*cos(ThetakP);

      xbj =  Q2/(2.*Pq);
      y   = (Pq)/(Pk);
      W2  = -Q2 + M_TARG2 + 2*Pq;

      double W2min_fm = -ro->get_fQ2max() + M_TARG2 + 2.*(E*nu - P*q);
      double W2max_fm = -ro->get_fQ2min() + M_TARG2 + 2.*(E*nu + P*q);
      if(W2 > W2max_fm || W2 < W2min_fm || W2max_fm <= W2min_fm || W2 <= 0)
      {
        //cout<<"event "<<ievt<<" W2 out of range "<<
	// W2min_fm<<" "<<W2<<" "<<W2max_fm<<endl;
        //out0<<"event "<<ievt<<" W2 out of range "<<
        //  W2min_fm<<" "<<W2<<" "<<W2max_fm<<endl;

        outside_klim++;
	outside_klim_W2++;
        dv->init_dvcs();
        init_event();
        init_track();
        continue;
      }

      double xmin_fm = Q2/(2.*(E*nu + P*q));
      double xmax_fm = Q2/(2.*(E*nu - P*q));
      xmin_fm = max(ro->get_fXbjmin(), xmin_fm);
      xmax_fm = min(ro->get_fXbjmax(), xmax_fm);
      if(xbj > xmax_fm || xbj < xmin_fm || xmax_fm <= xmin_fm)
      {
        //cout<<"event "<<ievt<<" xbj out of range "<<
        //  xmin_fm<<" "<<xbj<<" "<<xmax_fm<<endl;
        //out0<<"event "<<ievt<<" xbj out of range "<<
        //  xmin_fm<<" "<<xbj<<" "<<xmax_fm<<endl;

	outside_klim++;
	outside_klim_x2++;
        dv->init_dvcs();
        init_event();
        init_track();
        continue;
      }

      //double rndmt2   = rndm.Rndm();
      double rndmphi2 = rndm.Rndm();
      if(fmotion(rndmt, rndmphi2, Ipn, q, nu, W2, P, ThetakP, PhikP, 
                 Thetakq, Phikq, Output) == -1)
      {
        //cout<<"Error in Fermi motion subroutine"<<endl;
        //out0<<"Error in Fermi motion subroutine"<<endl;

        outside_klim++;
        dv->init_dvcs();
        init_event();
        init_track();
        continue;
      }

      t        = Output[0];
      nup      = Output[1];
      Ep       = Output[2];
      Thetakqp = Output[3];
      Phikqp   = Output[4];
      ThetakPp = Output[5];
      PhikPp   = Output[6];

    if(t > tmax || t < tmin)
    {
      //cout << "second one" << endl;
      //      cout<<"event "<<ievt<<"; t out of range "<<
      //tmin<<" "<<t<<" "<<tmax<<endl;
      //out0<<"event "<<ievt<<"; t out of range "<<
      //tmin<<" "<<t<<" "<<tmax<<endl;

      outside_klim++;
      outside_klim_t2++;
      dv->init_dvcs();
      init_event();
      init_track();
      continue;
    }

      qp = nup;
      Pp = sqrt(sqr(Ep) - M_TARG2);

      vqp = vect->fill_vect(qp, Thetakqp, Phikqp);

      vPp = vect->fill_vect(Pp, ThetakPp, PhikPp);

      //rotate the real photon 3-momentum vector from 
      //the inc. elec. frame to the virt. phot. frame
      vqp1 = vect->rot_vect(vqp, Thetakq, Phikq);

      //polar angle of real photon wrt virt. phot. frame
      Thetaqqp     = get_theta(vqp1.at(2),qp);

      //azimuthal angle of real photon wrt virt. phot. frame
      Phiqqp       = get_phi(vqp1.at(0),vqp1.at(1));

      //final nucleon 3-momentum vector wrt virt. phot. frame
      vPp1 = vect->add_vect(vP1, vect->diff_vect(vq2,vqp1));

      //polar angle of final nucleon wrt virt. phot. frame
      ThetaqPp     = get_theta(vPp1.at(2),Pp);

      //azimuthal angle of final nucleon wrt virt. phot. frame
      PhiqPp       = get_phi(vPp1.at(0), vPp1.at(1));

      //Angle between scat. elec. and real photon
      cosThetakpqp = cos(Phikqp)*sin(Thetakqp)*cos(Phikkp)*sin(Thetakkp) + 
                     sin(Phikqp)*sin(Thetakqp)*sin(Phikkp)*sin(Thetakkp) +
                     cos(Thetakqp)*cos(Thetakkp);
      Thetakpqp    = acos(cosThetakpqp);
    }
//=======================END FERMI MOTION PART

    ycol  = get_ycol(t, xbj, Q2);
    
    lujets_cc.N = 0;
    lujets_.n   = 0;
    for(int ii=0; ii<5; ii++)
    {
      for(int ij=0; ij<30; ij++)
      {
        lujets_cc.K[ii][ij] = 0;
        lujets_cc.P[ii][ij] = 0;
        lujets_cc.V[ii][ij] = 0;

        lujets_.k[ii][ij] = 0;
        lujets_.p[ii][ij] = 0;
        lujets_.v[ii][ij] = 0;
      }
    }

    hepevt.NEVHEP = 0;
    hepevt.NHEP   = 0;
    for(int ii=0; ii<NMXHEP; ii++)
    {
      hepevt.ISTHEP[ii] = 0;
      hepevt.IDHEP[ii]  = 0;
      for(int jj=0; jj<2; jj++) hepevt.JMOHEP[ii][jj] = 0;
      for(int jj=0; jj<2; jj++) hepevt.JDAHEP[ii][jj] = 0;
      for(int jj=0; jj<5; jj++) hepevt.PHEP[ii][jj]   = 0;
      for(int jj=0; jj<4; jj++) hepevt.VHEP[ii][jj]   = 0;
    }

    //Cross sections
    vPm = vect->init_vect();
    double Phi_b = Phiqqp;
    double xsec;
    if(trk.Process == 0)
      { 
      if(dv->phot_xsec(ro, ff, xbj, y, Q2, t, Phi_b) != 0 )
      {
        //cout<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;
        //out0<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;

        outside_klim++;
        dv->init_dvcs();
        init_event();
        init_track();
	continue;
      }
      ds_bh   = dv->hc0_BH + dv->hc1_BH*cos(pi(1) - Phi_b) + dv->hc2_BH*cos(2*(pi(1) - Phi_b)) + 
	        ro->get_fBheli()*ro->get_fTheli()*(dv->hc0_BH_LP + dv->hc1_BH_LP*cos(pi(1) - Phi_b));

      ds_dvcs = dv->hc0_DVCS + dv->hc1_DVCS*cos(pi(1) - Phi_b) + ro->get_fBheli()*dv->hs1_DVCS*sin(pi(1) - Phi_b) + 
	        ro->get_fBheli()*ro->get_fTheli()*(dv->hc0_DVCS_LP + dv->hc1_DVCS_LP*cos(pi(1) - Phi_b)) + 
		ro->get_fTheli()*dv->hs1_DVCS_LP*sin(pi(1) - Phi_b);

      ds_int  = dv->hc0_Int + dv->hc1_Int*cos(pi(1) - Phi_b) + dv->hc2_Int*cos(2*(pi(1) - Phi_b)) + 
	        ro->get_fBheli()*(dv->hs1_Int*sin(pi(1) - Phi_b) + dv->hs2_Int*sin(2*(pi(1) - Phi_b))) + 
		ro->get_fBheli()*ro->get_fTheli()*(
		dv->hc0_Int_LP + dv->hc1_Int_LP*cos(pi(1) - Phi_b) + dv->hc2_Int_LP*cos(2*(pi(1) - Phi_b))) + 
		ro->get_fTheli()*(dv->hs1_Int_LP*sin(pi(1) - Phi_b) + dv->hs2_Int_LP*sin(2*(pi(1) - Phi_b)));

      ds_tot  = ds_bh + ds_dvcs - ro->get_fBchg()*ds_int;
      xsec = ds_tot;
      }
    else if(trk.Process == 1)
    {
      ds_ms = ms_xsec(ro, xbj, Q2, t, Phi_b);
      xsec =  ds_ms;
      
      if(get_ms(ro, xbj, Q2, nu, nup, t, Phikkp, Phiqqp) != 0)
      {
        cout<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;
	out0<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;

        outside_klim++;
        dv->init_dvcs();
        init_event();
        init_track();
	continue;
      }

      for(int i =0; i<3; i++) vPm.at(i) = lujets_.p[i][0];
      Pm       = vect->mom_vect(vPm);
      Em       = sqrt(sqr(Pm) + m_ms(Ims,2));
      ThetakPm = get_theta(vPm.at(2),Pm);
      PhikPm   = get_phi(vPm.at(0),vPm.at(1)); 

      //correct the momentum and energy of final nucleon
      vPp      = vect->diff_vect(vect->add_vect(vq, vP), vPm);
      Pp       = vect->mom_vect(vPp);
      Ep       = sqrt(sqr(Pp) + m_targ(Ipn,2));
      ThetakPp = get_theta(vPp.at(2),Pp);
      PhikPp   = get_phi(vPp.at(0),vPp.at(1));
    }
    else
    {
      out0<<"You Should select a Process"<<endl;
    }

    //fill track
    //be carefull here
    //get_ms() returns the particles coming out from the pi0/eta decay
    //the number of particles, stocked in lujets_.n, is not constant

    trk.Ntracks       = 7;
    trk.TarA          = ro->get_fAt();
    trk.TarZ          = ro->get_fZt();
    trk.Eb            = Ee;

    trk.xbj           = xbj;
    trk.y             = y;
    trk.nu            = nu;
    trk.Q2            = Q2;
    trk.W2            = W2;
    trk.t             = t;
    trk.ycol          = ycol;

    trk.bh_xsec       = ds_bh;
    trk.dvcs_xsec     = ds_dvcs;
    trk.int_xsec      = ds_int;
    trk.tot_xsec      = ds_tot;
    trk.ms_xsec       = ds_ms;

    //inc. elec.
    trk.Type[0]   = elec_id();
    trk.Charge[0] = elec_ch();
    trk.Px[0]      = vk.at(0);
    trk.Py[0]      = vk.at(1);
    trk.Pz[0]      = vk.at(2);
    trk.P[0]       = k;
    trk.E[0]       = Ee;
    trk.Theta[0]   = 0.;
    trk.Phi[0]     = 0.;
    
    //target nucleon
    trk.Type[1]    = targ_id(Ipn);
    trk.Charge[1]  = targ_ch(Ipn);
    trk.Px[1]      = vP.at(0);
    trk.Py[1]      = vP.at(1);
    trk.Pz[1]      = vP.at(2);
    trk.P[1]       = P;
    trk.E[1]       = E;
    trk.Theta[1]   = ThetakP;
    trk.Phi[1]     = PhikP;

    //scat. elec.
    trk.Type[2]    = elec_id();
    trk.Charge[2]  = elec_ch();
    trk.Px[2]      = vkp.at(0);
    trk.Py[2]      = vkp.at(1);
    trk.Pz[2]      = vkp.at(2);
    trk.P[2]       = kp;
    trk.E[2]       = Eep;
    trk.Theta[2]   = Thetakkp;
    trk.Phi[2]     = Phikkp;

    //virt. photon
    trk.Type[3]    = phot_id();
    trk.Charge[3]  = phot_ch();
    trk.Px[3]      = vq.at(0);
    trk.Py[3]      = vq.at(1);
    trk.Pz[3]      = vq.at(2);
    trk.P[3]       = q;
    trk.E[3]       = nu;
    trk.Theta[3]   = Thetakq;
    trk.Phi[3]     = Phikq;

    //recoil nucleon/nucleus
    trk.Type[4]    = RECO_ID;
    trk.Charge[4]  = RECO_CH;
    trk.Px[4]      = vPr.at(0);
    trk.Py[4]      = vPr.at(1);
    trk.Pz[4]      = vPr.at(2);
    trk.P[4]       = Pr;
    trk.E[4]       = Er;
    trk.Theta[4]   = ThetakPr;
    trk.Phi[4]     = PhikPr;

    //scat. nucleon
    trk.Type[5]    = targ_id(Ipn);
    trk.Charge[5]  = targ_ch(Ipn);
    trk.Px[5]      = vPp.at(0);
    trk.Py[5]      = vPp.at(1);
    trk.Pz[5]      = vPp.at(2);
    trk.P[5]       = Pp;
    trk.E[5]       = Ep;
    trk.Theta[5]   = ThetakPp;
    trk.Phi[5]     = PhikPp;

    int NNT = trk.Ntracks - 1;

    if(trk.Process == 0)
    {
      //out. photon
      trk.Type[6]    = phot_id();
      trk.Charge[6]  = phot_ch();
      trk.Px[6]      = vqp.at(0);
      trk.Py[6]      = vqp.at(1);
      trk.Pz[6]      = vqp.at(2);
      trk.P[6]       = qp;
      trk.E[6]       = nup;
      trk.Theta[6]   = Thetakqp;
      trk.Phi[6]     = Phikqp;
    }
    else if(trk.Process == 1)
    {
      trk.Ntracks  += lujets_.n - 1;
      for(int in=0; in<lujets_.n; in++)
      {
        int nnt = NNT + in;
	trk.Type[nnt]   = lujets_.k[1][in];
        if(trk.Type[nnt] == neut_id() || trk.Type[nnt] == eta_id() 
	|| trk.Type[nnt] == pi0_id() || trk.Type[nnt] == phot_id())
	{
	  trk.Charge[nnt] = 0;
	}
	else if(abs(trk.Type[nnt]) == elec_id())
	{
	  trk.Charge[nnt] = -sign(trk.Type[nnt]);
	}
	else
	{
	  trk.Charge[nnt] = sign(trk.Type[nnt]);
	}
        trk.Px[nnt]    = lujets_.p[0][in];
        trk.Py[nnt]    = lujets_.p[1][in];
        trk.Pz[nnt]    = lujets_.p[2][in];
        trk.P[nnt]     = sqrt(sqr(trk.Px[nnt]) + sqr(trk.Py[nnt]) + sqr(trk.Pz[nnt]));
        trk.E[nnt]     = lujets_.p[3][in];
        trk.Theta[nnt] = get_theta(trk.Pz[nnt], trk.P[nnt]);
        trk.Phi[nnt]   = get_phi(trk.Px[nnt], trk.Py[nnt]);
      }
    }
    
    if(ro->get_fAscii() != 0)
    {
      //Fill lujets_cc common block
      lujets_cc.N = trk.Ntracks;
      for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[0][ii] = 1;

      lujets_cc.K[1][0] = elec_id();    //inc. electron
      lujets_cc.K[1][1] = targ_id(Ipn); //proton target
      lujets_cc.K[1][2] = elec_id();    //scat. electron
      lujets_cc.K[1][3] = phot_id();    //virtual photon
      lujets_cc.K[1][4] = RECO_ID;      //recoil proton
      lujets_cc.K[1][5] = targ_id(Ipn); //scat. proton

      for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[2][ii] = 0;
      for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[3][ii] = 0;
      for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[4][ii] = 0;
 
      if(trk.Process == 0)
      {
        lujets_cc.K[0][6] = 1;
        lujets_cc.K[1][6] = phot_id();   //real photon
	lujets_cc.K[2][6] = 0;
	lujets_cc.K[3][6] = 0;
	lujets_cc.K[4][6] = 0;
      }
      else if(trk.Process == 1)
      {
        for(int in=0; in<lujets_.n; in++)
          for(int jn=0; jn<5; jn++)
	    lujets_cc.K[jn][NNT + in] = lujets_.k[jn][in];
      }

      lujets_cc.P[0][0] = vk.at(0);
      lujets_cc.P[0][1] = vP.at(0);
      lujets_cc.P[0][2] = vkp.at(0);
      lujets_cc.P[0][3] = vq.at(0);
      lujets_cc.P[0][4] = vPr.at(0);
      lujets_cc.P[0][5] = vPp.at(0);

      lujets_cc.P[1][0] = vk.at(1);
      lujets_cc.P[1][1] = vP.at(1);
      lujets_cc.P[1][2] = vkp.at(1);
      lujets_cc.P[1][3] = vq.at(1);
      lujets_cc.P[1][4] = vPr.at(1);
      lujets_cc.P[1][5] = vPp.at(1);

      lujets_cc.P[2][0] = vk.at(2);
      lujets_cc.P[2][1] = vP.at(2);
      lujets_cc.P[2][2] = vkp.at(2);
      lujets_cc.P[2][3] = vq.at(2);
      lujets_cc.P[2][4] = vPr.at(2);
      lujets_cc.P[2][5] = vPp.at(2);

      lujets_cc.P[3][0] = Ee;
      lujets_cc.P[3][1] = E;
      lujets_cc.P[3][2] = Eep;
      lujets_cc.P[3][3] = nu;
      lujets_cc.P[3][4] = Er;
      lujets_cc.P[3][5] = Ep;

      lujets_cc.P[4][0] = m_elec(1);
      lujets_cc.P[4][1] = m_targ(Ipn,1);
      lujets_cc.P[4][2] = m_elec(1);
      lujets_cc.P[4][3] = Q2;
      lujets_cc.P[4][4] = M_RECO;
      lujets_cc.P[4][5] = m_targ(Ipn,1);

      if(trk.Process == 0)
      {
	lujets_cc.P[0][6] = vqp.at(0);
        lujets_cc.P[1][6] = vqp.at(1);
        lujets_cc.P[2][6] = vqp.at(2);
        lujets_cc.P[3][6] = nup;
        lujets_cc.P[4][6] = m_phot(1);
      }
      else if(trk.Process == 1)
      {
        for(int in=0; in<lujets_.n; in++)
          for(int jn=0; jn<5; jn++)
	    lujets_cc.P[jn][NNT + in] = lujets_.p[jn][in];
      }

      plu[0]            = elec_ch();    
      plu[1]            = targ_ch(Ipn);
      plu[2]            = elec_ch();
      plu[3]            = phot_ch();
      plu[4]            = RECO_CH;
      plu[5]            = targ_ch(Ipn);

      if(trk.Process == 0)
      {
        plu[6]  = phot_ch();
      }
      else if(trk.Process == 1)
      {
        for(int in=0; in<lujets_.n; in++)
          plu[in] = trk.Charge[NNT + in];
      }

      for(int ii=0; ii<5; ii++)
        for(int ij=0; ij<lujets_cc.N; ij++) lujets_cc.V[ii][ij] = 0;

      //Fill hepevt common block
      hepevt.NHEP = lujets_cc.N;

      for(int ii=0; ii<hepevt.NHEP; ii++)
      {
        hepevt.ISTHEP[ii] = 0;
        if(lujets_cc.K[0][ii] >=  1 && lujets_cc.K[0][ii] <= 10)
	  hepevt.ISTHEP[ii] = 1;
        if(lujets_cc.K[0][ii] >= 11 && lujets_cc.K[0][ii] <= 20)
	  hepevt.ISTHEP[ii] = 2;
        if(lujets_cc.K[0][ii] >= 21 && lujets_cc.K[0][ii] <= 30)
	  hepevt.ISTHEP[ii] = 3;
        if(lujets_cc.K[0][ii] >= 31 && lujets_cc.K[0][ii] <= 100)
	  hepevt.ISTHEP[ii] = lujets_cc.K[0][ii];
        hepevt.IDHEP[ii]     = lujets_cc.K[1][ii];
        hepevt.JMOHEP[ii][0] = lujets_cc.K[2][ii];
        hepevt.JMOHEP[ii][1] = 0;
        if(lujets_cc.K[0][ii] != 3 && lujets_cc.K[0][ii] != 13 && 
	   lujets_cc.K[0][ii] != 14)
	{
          hepevt.JDAHEP[ii][0] = lujets_cc.K[3][ii];
          hepevt.JDAHEP[ii][1] = lujets_cc.K[4][ii];
        }
        else
	{
          hepevt.JDAHEP[ii][0] = 0;
          hepevt.JDAHEP[ii][1] = 0;
        }
        for(int ij=0; ij<5; ij++) hepevt.PHEP[ii][ij] = lujets_cc.P[ij][ii];
        for(int ij=0; ij<4; ij++) hepevt.VHEP[ii][ij] = lujets_cc.V[ij][ii];
      }

      //Write the output file      
      //dump_file(ro->get_fMode(), xsec, ptr);

      //=====================================KEEP REJECT======================================//
      //Output file is written according to probability distrib.      
      //Limits for this keep/reject have been tuned by Rong Wang, Ipno
      if((trk.Process == 0 && xsec>rndm.Rndm()*2)||(trk.Process == 1 && xsec>rndm.Rndm()*800))
        {
          kr = 1;
          counter_kept++;
          dump_file(ro->get_fMode(), ro->get_fProc(),xsec);
        }
      else
        {
          counter_rejected++;
          kr= 0;
        }
      trk.kr = kr;

    }
    if(ro->get_fNtp() && trk.kr==1) 
      {
	tree->Fill();
      }

    if(trk.Struck_Nucl == 0)
    {
      Nhit_neut++;
    }
    else if(trk.Struck_Nucl == 1)
    {
      Nhit_prot++;
    }

    if(ro->get_fNtp() && Nevts_per_Ntup > ro->get_fNevtsPerNtup())
    {
      Nevts_per_Ntup=0;
      ntp_cnt++;
      f->Write();
      f->Close();
      if(ntp_cnt<10)
      {
        sprintf(root_file, "%s/ntup_%s_%3.2fgev_00%d.root", 
	        fDir.c_str(), target.c_str(), Ee, ntp_cnt);
      }
      else if(ntp_cnt<100)
      {
        sprintf(root_file, "%s/ntup_%s_%3.2fgev_%0d.root", 
	        fDir.c_str(), target.c_str(), Ee, ntp_cnt);
      }
      else
      {
        sprintf(root_file, "%s/ntup_%s_%3.2fgev_%d.root", 
	        fDir.c_str(), target.c_str(), Ee, ntp_cnt);
      }
      cout<<"Start Filling ntuple "<<root_file<<endl;
      out0<<"Start Filling ntuple"<<endl;
      f    = new TFile(root_file,"RECREATE");
      tree = new TTree("DVCS", "/dvcs_tree");
      init_tree(tree);
      //TTree::SetMaxTreeSize(2000000000);
    }

    Xmin    = min(Xmin, xbj);
    Xmax    = max(Xmax, xbj);
    Ymin    = min(Ymin, y);
    Ymax    = max(Ymax, y);
    NUmin   = min(NUmin, nu);
    NUmax   = max(NUmax, nu);
    QQ2min  = min(QQ2min, Q2);
    QQ2max  = max(QQ2max,Q2);
    WW2min  = min(WW2min, W2);
    WW2max  = max(WW2max,W2);
    Tmin    = min(Tmin, t);
    Tmax    = max(Tmax, t);
    NUPmin  = min(NUPmin, nup);
    NUPmax  = max(NUPmax, nup);
    EPmin   = min(EPmin, Ep);
    EPmax   = max(EPmax, Ep);

  }
//=======================loop ends here

  if(ro->get_fNtp())
  {
    f->Write();
    f->Close();
  }

  cout<<"initial Number of events              "<<ro->get_fNevts()<<endl;
  cout<<"Number of hit neutrons                "<<Nhit_neut       <<endl;
  cout<<"Number of hit protons                 "<<Nhit_prot       <<endl;
  if (iApZ_stored == 30 || iApZ_stored == 27)
    {
      cout<<"Number of events on N             "<<count_N      <<endl;
      cout<<"Number of events on H or D        "<<count_p      <<endl;
    }
  cout<<"Number of events kept   " << counter_kept   <<endl;
  cout<<"Number of events rejected   "<< counter_rejected   <<endl;
  cout<<"Number of events Out of kin. limits   "<< outside_klim   <<endl;
  cout<<"Number of events Out of t limits   "<< outside_klim_t   <<endl;
  cout<<"Number of events Out of t limits 2  "<< outside_klim_t2   <<endl;
  cout<<"Number of events Out of nu limits   "<< outside_klim_nu   <<endl;
  cout<<"Number of events Out of x limits 2   "<< outside_klim_x2   <<endl;
  cout<<"Number of events Out of W2 limits   "<< outside_klim_W2   <<endl;
  cout<<"kinem. limits :"<<endl;
  cout<<"xbjmin   = "<<setw(10)<<Xmin  <<";    xbjmax   = "<<setw(10)<<Xmax  <<endl;
  cout<<"ymin     = "<<setw(10)<<Ymin  <<";    ymax     = "<<setw(10)<<Ymax  <<endl;
  cout<<"numin    = "<<setw(10)<<NUmin <<";    numax    = "<<setw(10)<<NUmax <<endl;
  cout<<"Q2min    = "<<setw(10)<<QQ2min<<";    Q2max    = "<<setw(10)<<QQ2max<<endl;
  cout<<"W2min    = "<<setw(10)<<WW2min<<";    W2max    = "<<setw(10)<<WW2max<<endl;
  cout<<"tmin     = "<<setw(10)<<Tmin  <<";    tmax     = "<<setw(10)<<Tmax  <<endl;
  cout<<"nupmin   = "<<setw(10)<<NUPmin<<";    nupmax   = "<<setw(10)<<NUPmax<<endl;
  cout<<"Epmin    = "<<setw(10)<<EPmin <<";    Epmax    = "<<setw(10)<<EPmax <<endl;
  cout<<endl;

  //Output summary file
  out0<<"initial Number of events              "<<ro->get_fNevts()<<endl;
  out0<<"Number of hit neutrons                "<<Nhit_neut       <<endl;
  out0<<"Number of hit protons                 "<<Nhit_prot       <<endl;
  if (iApZ_stored == 30 || iApZ_stored == 27)
    {
      cout<<"Number of events on N             "<<count_N      <<endl;
      cout<<"Number of events on H or D        "<<count_p      <<endl;
    }
  out0<<"Number of events kept   " << counter_kept   <<endl;
  out0<<"Number of events rejected   "<< counter_rejected   <<endl;
  out0<<"Number of events Out of kin. limits   "<< outside_klim   <<endl;
  out0<<"Number of events Out of t limits   "<< outside_klim_t   <<endl;
  out0<<"Number of events Out of t limits 2  "<< outside_klim_t2   <<endl;
  out0<<"Number of events Out of nu limits   "<< outside_klim_nu   <<endl;
  out0<<"Number of events Out of x limits 2   "<< outside_klim_x2   <<endl;
  out0<<"Number of events Out of W2 limits   "<< outside_klim_W2   <<endl;
  out0<<"kinem. limits :"<<endl;
  out0<<"xbjmin   = "<<setw(10)<<Xmin  <<";    xbjmax   = "<<setw(10)<<Xmax  <<endl;
  out0<<"ymin     = "<<setw(10)<<Ymin  <<";    ymax     = "<<setw(10)<<Ymax  <<endl;
  out0<<"numin    = "<<setw(10)<<NUmin <<";    numax    = "<<setw(10)<<NUmax <<endl;
  out0<<"Q2min    = "<<setw(10)<<QQ2min<<";    Q2max    = "<<setw(10)<<QQ2max<<endl;
  out0<<"W2min    = "<<setw(10)<<WW2min<<";    W2max    = "<<setw(10)<<WW2max<<endl;
  out0<<"tmin     = "<<setw(10)<<Tmin  <<";    tmax     = "<<setw(10)<<Tmax  <<endl;
  out0<<"nupmin   = "<<setw(10)<<NUPmin<<";    nupmax   = "<<setw(10)<<NUPmax<<endl;
  out0<<"Epmin    = "<<setw(10)<<EPmin <<";    Epmax    = "<<setw(10)<<EPmax <<endl;
  out0<<endl;

  return 0;
}

