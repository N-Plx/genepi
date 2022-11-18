/*
  Author Ahmed El Alaoui, LPSC (Grenoble)
  Ref: BMK, hep-ph/0112108
  
  Additional features:
  Fermi motion on N : Mathieu Erhart, IJCLab
  NH3/ND3 targets : NoÅÈmie Pilleux, IJCLab
  Restructuring / kinematic bug-fixes, phi-meson production : Daria Sokhan, CEA Saclay / Glasgow
*/
#include <cstdlib>

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
// ==== These variables don't seem to be used anywhere, so have been removed. The meson info is in the same variables as the photon info.

//polar and azimuthal angles between scat. elec. and real photon =====
double Thetakpqp, Phikpqp;

//kinematics =========================================================
double Q2;   //photon virtuality
double xbj;  //xbjorken variable
double y;    //elec. energy fraction 
double W2;   //invariant mass of nucleon + virtual photon system
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
  //  seed = 333;
  
  //=================Seed for generation of kin variables=========================//
  TRandom1 rndm;           // changed from the original TRandom1, which is very slow, to TRandom3 
  rndm.SetSeed(seed);
  gRandom->SetSeed(seed);
  
  double rndmt;           // for generating t later on

  //Read input A and Z
  int iApZ = ro->get_fAt() + ro->get_fZt();
  //We store this value in case of a molecular target where we will change iApZ when randomly selecting the target
  int iApZ_stored = iApZ; 
  //Some counters for testing/debugging
  int count_N = 0;
  int count_p = 0;
  int counter_kept = 0;
  int counter_kept_neut = 0;
  int counter_kept_prot = 0;
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
  std::string gpd_tbl = "gpd_table.dat";
  if (std::getenv("GENEPI") != NULL) {
      std::string absolute_path = std::getenv("GENEPI");
      absolute_path += "/";
      absolute_path += gpd_tbl;
      gpd_tbl = absolute_path;
  }
  if (fexist(gpd_tbl.c_str()))
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
  
  double Pq;   // P.q 
  
  int    Nhit_neut(0);
  int    Nhit_prot(0);
  int    run(0);
  int    Nevts_per_Ntup(0);
  int    outside_klim(0);
  int    outside_klim_t(0);
  int    outside_klim_nu(0);
  int    outside_klim_t2(0);
  int    outside_klim_tmin(0);
  int    outside_klim_x(0);
  int    outside_klim_W2(0);
  int    outside_klim_Q2(0);
  int    outside_klim_xb(0);
  int    outside_klim_y(0);
  int    outside_klim_nup(0);
  int    outside_klim_Ep(0);
  int    outside_klim_fm(0);
  int    outside_klim_xsec(0);  

  //Recoil, target
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
  
  long long ievt = 0;

  //while(counter_kept_neut < ro->get_fNevts())
  //{
      for(long long ievt=0; ievt<ro->get_fNevts(); ievt++)
        {
      //Get event info from input
      trk.Evt_ID  = ievt;
      trk.Bheli   = ro->get_fBheli();
      trk.Theli   = ro->get_fTheli();
      trk.Process = ro->get_fProc();
      trk.ExcMS   = ro->get_fIms();
      //Nevts_per_Ntup++;
      
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
      
      else if(iApZ_stored == 30) //ND3 6/20 chances to interact with deuterium
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
      
      
      /////////////// Start of event kinematics /////////////////
      
      //Kinem limits computed from input file
      get_kinem_limits(ro);  // Don't do this!! Function uses collider kinematics, neglects target mass!
      
      
      // Generate Q2, Xbj:
      double rndmQ2, rndmXbj;
      rndmQ2  = rndm.Rndm();
      Q2      = ro->get_fQ2min() + (ro->get_fQ2max() - ro->get_fQ2min())*rndmQ2;
      rndmXbj = rndm.Rndm();
      xbj     = ro->get_fXbjmin() + (ro->get_fXbjmax() - ro->get_fXbjmin())*rndmXbj;      
      
      //inc. elec.
      Ee = ro->get_fEb();
      k  = sqrt(sqr(Ee) - m_elec(2));
      vk = vect->fill_vect(k, 0., 0.);
      
      // On the basis of the beam electron and generated Q2, one can calculate the max possible nu
      // this will happen for the max scattering angle (i.e.: where theta of scattered electron = 180 deg in lab frame):
      
      double numax_kin = (k * sqrt(Q2*(4.*m_elec(2) + Q2)) - Ee*Q2) / (2.*m_elec(2));    // this is the full formula
      
      // The approximation which neglects electron mass is also valid at these kinematics:
      // numax_kin = Ee - Q2/(4.*Ee);

      
      // If the initial nucleon isn't stationary in the lab frame:
      // Calculates E,P and ThetaqP for initial nucleon, and nu,q for virtual photon.
      if (iApZ == 3 || iApZ == 8 || iApZ == 21)
	{
	  //generate initial nucleon momentum wrt the inc. elec. frame
	  int ifermi(1);
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
	  
	  //generate polar angle of initial nucleon wrt the q-direction: needed to calculate nu
	  ThetaqP = pi(1)*rndm.Rndm();
	  
	  // azimuthal angle of initial nucleon wrt the q-direction:
	  double PhiqP = 2.*pi(1)*rndm.Rndm();
	  
	  // construct the initial nucleon vector wrt q-axis:
	  vP1 = vect->fill_vect(P, ThetaqP, PhiqP);
	  
	  // nu is obtained from solving the quadratic equation which arises from the full expression for xbj = Q2/(2P.q)
	  // this is the "sqrt(b^2 - 4ac)" part of the quadratic solution:
	  double nupart  = Q2*P*cos(ThetaqP) * sqrt((1/sqr(xbj)) + 4*(M_TARG2 + sqr(P*sin(ThetaqP)))/Q2);
	  
	  // the two possible solutions for nu:
	  double nuplus  = (E*Q2/xbj + nupart) / (2*(M_TARG2 + sqr(P*sin(ThetaqP))));
	  double numinus = (E*Q2/xbj - nupart) / (2*(M_TARG2 + sqr(P*sin(ThetaqP))));
	  
	  // Only one of them will work, check which one -- this should give you Q2:
	  double test1 = xbj * 2. * (E*nuplus - P*cos(ThetaqP)*sqrt(sqr(nuplus) + Q2));
	  double test2 = xbj * 2. * (E*numinus - P*cos(ThetaqP)*sqrt(sqr(numinus) + Q2));
	 
	  // It's always the + solution, but test just in case anyway: 
	  if (fabs(test1 - Q2) < fabs(test2 - Q2)) nu = nuplus;
	  else if (fabs(test2 - Q2) < fabs(test1 - Q2)) nu = numinus;
	  else 
	    {
	      cout << "Ambiguous calculation of nu. Aborting event." << endl;
	      cout << fabs(test1 - Q2) << ", " << fabs(test2 - Q2) << endl;
	      outside_klim++;
              outside_klim_nu++;
              dv->init_dvcs();
              init_event();
              init_track();
              continue;
	    }
	  
	  if (nu < 0 || nu > numax_kin)    // max value of nu calculated above, where the beam is defined
	    {
	      outside_klim++;
	      outside_klim_nu++;
	      dv->init_dvcs();
	      init_event();
	      init_track();
	      continue;
	    }
	  
	  q  = sqrt(Q2+sqr(nu));
	  
	  // Calculate P.q, this will come in useful later:	  
	  Pq = E*nu - P*q*cos(ThetaqP);
	  
	}  // Fermi-target bit closed
      
      else   // standard stationary target kinematics, get the full initial nucleon kinematics, nu, q and y:
	{
	  //Target nucleon -- assumes stationary!
	  P = 0.0;
          E = M_TARG;
          vP = vect->fill_vect(P, 0., 0.);
	  
	  nu      = Q2/(2.*M_TARG*xbj);   // only true if stationary target!
	  y       = nu/ro->get_fEb();     // only true for stationary target!
	  q       = sqrt(Q2+sqr(nu));
	  Pq      = M_TARG * nu;          // P.q, expression only true for a stationary target!
	}
      
      W2      = M_TARG2 - Q2 + Q2/xbj;   // also true for a moving target. No problem here.      
      
      // Check that you're within "good" kinematic limits for nu and Q2:
      
      numin   = ro->get_fNumin();  // Given in the set-up file 
      // numax_kin is defined above, where beam electron parameters are defined 
      numax   = min(numax_kin, ro->get_fNumax());
      
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
     
      // This is a sanity check -- nu used to calculate Q2max was calculated with Q2 generated, so Q2 should of course be less than this Q2max.
      // If any events don't meet this condition, there's something very wrong in the code. 
      Q2min   = ro->get_fQ2min();
      Q2max = 2. * (k*k - Ee*nu + k * sqrt(k*k + nu*nu - 2.*Ee*nu));    // this is the full expression, not neglecting electron mass
      // Approximating electron mass to zero holds fine in these kinematics:
      // Q2max   = 4.*Ee*(Ee -  nu);   // neglects electron mass
      Q2max   = min(Q2max, ro->get_fQ2max());
      
      if(Q2 > Q2max || Q2 < Q2min)
	{
	  cout<<"Event = "<<ievt<<"; Q2 out of range "<<
	    Q2min<<" "<<Q2<<" "<<Q2max<<endl;
	  out0<<"Event = "<<ievt<<"; Q2 out of range "<<
	    Q2min<<" "<<Q2<<" "<<Q2max<<endl;
	  
	  outside_klim++;
	  outside_klim_Q2++;
	  dv->init_dvcs();
	  init_event();
	  init_track();
	  continue;
	}
      
      
      // Now for the xbj limits, which is again a sanity check as it depends on quantities which were calculated using xbj that was already generated.:
      // This is less straight-forward. xb = Q2 / (P.q). Can get min xb either by using min Q2 (and then calculated nu) or max P.q (and the generated Q2).
      // Check for both bases.
      // Using max P.q, which you get when nu is max and (for non-stationary target) cosThetaqp = -1: 
      xmin = Q2 / (2 * (E*numax + P * sqrt(pow(numax,2) + Q2)));    // works for moving and stationary target, also for real electron mass
      xmin  = max(xmin, ro->get_fXbjmin());
      // Using min Q2 and calculated P.q:
      double xmin_other = Q2min / (2*Pq);
      xmin  = max(xmin, xmin_other);      
      
      // max value of xbj is obtained when Q2 is at its max:
      xmax = Q2max / (2*Pq);
      xmax  = min(xmax, ro->get_fXbjmax());
      // or when P.q is at its min (when nu is min and, for non-stationary target, cosThetaqp = 1):
      double xmax_other = Q2 / (2 * (E*numin - P * sqrt(pow(numin,2) + Q2)));
      if (xmax_other > 0.) xmax  = min(xmax, xmax_other);    // for some cases, P.q could be negative above, which then messes up xmax
      
      if(xbj > xmax || xbj < xmin)
	{
	  cout<<"event "<<ievt<<"; xbj out of range "<<
	    xmin<<" "<<xbj<<" "<<xmax<<endl;
	  out0<<"event "<<ievt<<"; xbj out of range "<<
	    xmin<<" "<<xbj<<" "<<xmax<<endl;
	  
	  outside_klim++;
	  outside_klim_xb++;
	  dv->init_dvcs();
	  init_event();
	  init_track();
	  continue;
	}	  
      
      // Also check the W2, just using the Q2min and Q2max values:
      
      double W2min = M_TARG2 - Q2max  + 2.*Pq;     // in the stationary target case, P.q = Pq = M_TARG * nu
      double W2max = M_TARG2 - Q2min  + 2.*Pq;
      
      if (iApZ > 2 )    // moving target, recalculate W2 max and min values:
	{ 
	  W2min = M_TARG2 - Q2max  + 2.*(E*nu - P*q);
	  W2max = M_TARG2 - Q2min  + 2.*(E*nu + P*q);
	}    
     
      W2max   = min(W2max, ro->get_fW2max());
      W2min   = max(W2min, ro->get_fW2min());
 
      if(W2 > W2max || W2 < W2min || W2max <= W2min || W2 <= 0)
	{
	  //	  cout<<"event "<<ievt<<" W2 out of range "<<                                                                                        
	  // W2min<<" "<<W2<<" "<<W2max<<endl;                                                                                                                      
	  // out0<<"event "<<ievt<<" W2 out of range "<<                                                                                                                   
	  // W2min<<" "<<W2<<" "<<W2max<<endl;                                                                          
	  
	  outside_klim++;
	  outside_klim_W2++;
	  dv->init_dvcs();
	  init_event();
	  init_track();
	  continue;
	}
      
      // Down to here Q2, xbj, nu and W2 checked to be within limits. Can proceed.
      
      //polar angle between inc. and scat. elec.
      cosThetakkp = (k*k - Ee*nu - Q2/2.) / (k*sqrt(k*k + nu*nu - 2.*Ee*nu));
      if(abs(cosThetakkp) > 1.)
	{
	  cout<<"event "<<ievt<<
	    "; cosThetakkp larger than one; cosThetakkp = "<<cosThetakkp<< endl;
	  cosThetakkp = sign(cosThetakkp);
	}
      
      Thetakkp =  acos(cosThetakkp);//[0,PI]
      
      // this check here seems redundant, since cosThetakkp is forced to be +/-1 at most:
      if(Thetakkp > pi(1)) 
	cout<<"Thetakkp larger than pi; cosThetakkp = "<<cosThetakkp<< endl;
      
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
      
      
      // Next take care of the moving initial nucleon, to construct its full four-momentum:
      
      if (iApZ > 2 )     // moving target
	{ 
	  // An angle of initial nucleon wrt q-vector has already been defined: ThetaqP. Pick a random phi angle about the same axis:
	  double PhiqP = 2.*pi(1)*rndm.Rndm();
	  
	  // Construct a vector wrt q-axis:
	  vP1 = vect->fill_vect(P, ThetaqP, PhiqP);
	  
	  // Now rotate it to align it with the z-axis, given by the beam vector k:
	  vP = vect->invrot_vect(vP1, Thetakq, Phikq);
	  
	  //polar angle of initial nucleon in beam electron lab-frame: 
	  ThetakP = get_theta(vP.at(2),P);
	  
	  //azimuthal angle of initial nucleon in beam electron lab-frame:
	  PhikP = get_phi(vP.at(0),vP.at(1));
	  
	  double Pk = E*Ee - P*k*cos(ThetakP);    // P.k dot-product	  
	  
	  y  = (Pq)/(Pk);    // this has already been defined for a stationary target scenario
	  
	  // The next part constructs the spectator nucleon in the nucleus (it's called "recoil" here):
	  //recoil nucleus 3-momentum vector wrt inc. elec. frame  (usually called spectator!)
	  Pr  = P;
          Er  = sqrt(sqr(P) + M_RECO2);
          vPr = vect->diff_vect(v0, vP);      // v0 is a null vector, this is just to obtain vPr which is equal and opposite to vP
	  
	  //recoil nucleus polar angle wrt the inc. elec. frame
	  ThetakPr     = get_theta(vPr.at(2),Pr);
	  
	  //recoil nucleus azimuthal angle wrt the inc. elec. frame
	  PhikPr       = get_phi(vPr.at(0),vPr.at(1));
	}
      
      // Now y has been calculated regardless of whether target is stationary or not.
      
      if(y > ro->get_fYmax() || y < ro->get_fYmin())
	{
	  cout<<"event "<<ievt<<"; y out of range "<<
	    ro->get_fYmin()<<" "<<y<<" "<<ro->get_fYmax()<<endl;
	  out0<<"event "<<ievt<<"; y out of range "<<
	    ro->get_fYmin()<<" "<<y<<" "<<ro->get_fYmax()<<endl;
	  
	  outside_klim++;
	  outside_klim_y++;
	  dv->init_dvcs();
	  init_event();
	  init_track();
	  continue;
	}
      
      
      // By this point have four-momenta of: beam e, scattered e, q, target p, spectator p if it exists. The whole initial state. 
      
      // Get the final state: produced particle (photon or meson) and the scattered nucleon (either same as target or different, if a charged meson was produced)
      
      // Set the square mass of the scattering nucleon in the final state (used for the meson-production case)  
      // This should be set somewhere at the start!! For now hard-code it here and use target mass.
      // It's equal to the square target mass only for produced neutrals (photon, pi0, etc).  
      double M_scat2 = M_TARG2;

      if (iApZ < 3 )   // stationary target scenario
        {
          // Upper and lower limits of nup, Ep and t
	  
          if  (trk.Process != 1)     // DVCS
            {
              // Assuming stationary target and produced photon:
              // derived by calculating the energy of the photon from energy-momentum conservation
              
              nupmin = (m_targ(Ipn,1)*nu - Q2/2.)/(m_targ(Ipn,1) + nu + q);   // assumes stationary target and real photon in final state!
              nupmax = (m_targ(Ipn,1)*nu - Q2/2.)/(m_targ(Ipn,1) + nu - q);   // assumes stationary target and real photon in final state!
            }
          
          else if (trk.Process == 1)     // meson-production
            {
              // If a meson is produced, the expression for energy of produced meson is a horrid quadratic.
              // Assumes stationary target but makes no assumption about produced particle or mass of the scattering nucleon (can be same as target or different):
              
              double U = pow(nu + M_TARG,2) - q*q;   // group some terms to make the expression below neater
              
              // Set the square mass of the produced meson:
              double M_mes2 = m_ms(Ims,2);
              
              nupmin = ((nu + M_TARG)*(U + M_mes2 - M_scat2) - q * sqrt(pow(U + M_mes2 - M_scat2,2) - 4.*M_mes2*U)) / (2.*U);  // stationary target, produced meson
              nupmax = ((nu + M_TARG)*(U + M_mes2 - M_scat2) + q * sqrt(pow(U + M_mes2 - M_scat2,2) - 4.*M_mes2*U)) / (2.*U);  // stationary target, produced meson
            }
          
          Epmin = M_TARG + nu - nupmax;
          Epmax = M_TARG + nu - nupmin;
          
          if (trk.Process != 1)     // DVCS
	    {
              tmin   = 2.* M_TARG * (nupmin - nu);   // assumes stationary target!
              tmax   = 2.* M_TARG * (nupmax - nu);
            }
          else if (trk.Process == 1)     // meson-production, general case (works for charged mesons)
            {
              tmin   = M_TARG2 + M_scat2 - 2.* M_TARG * (M_TARG + nu - nupmin);   // assumes stationary target!
              tmax   = M_TARG2 + M_scat2 - 2.* M_TARG * (M_TARG + nu - nupmax);
            }
          
	  //cout << "User-defined t-range: " << ro->get_ftmin() << ", " << ro->get_ftmax() << endl;
	  //cout << "Calculated t-range: " << tmin << ", " << tmax << endl;

          tmin   = max(tmin, ro->get_ftmin());
          tmax   = min(tmax, ro->get_ftmax());
          
          //generate t:
          rndmt = rndm.Rndm();
          t     = tmin + (tmax - tmin)*rndmt;
          if(t > tmax || t < tmin)
            {
	      // cout<<"event "<<ievt<<"; t out of range "<<tmin<<" "<<t<<" "<<tmax<<endl;
              //out0<<"event "<<ievt<<"; t out of range "<<tmin<<" "<<t<<" "<<tmax<<endl;
	      
              outside_klim++;
              outside_klim_t++;
              dv->init_dvcs();
              init_event();
              init_track();
              continue;
            }
          
          // Now get the info on the final state particles (scattered nucleon and produced photon or meson):
          
          if (trk.Process != 1)     // DVCS
            {
              // real photon energy and momentum
              nup = nu + t/(2.*M_TARG);           // assumes stationary target, neutral produced particle
              qp  = nup;                          // assumes real photon of zero mass!
              
              //final nucleon energy
              Ep = M_TARG - t/(2.*M_TARG);   // assumes stationary target, produced neutral particle
            }
          
          else if (trk.Process == 1)     // meson-production
            {
              // meson energy and momentum
              nup = M_TARG + nu - (M_TARG2 + M_scat2 - t)/(2.*M_TARG);
              qp = sqrt(sqr(nup) - m_ms(Ims,2));
              
              //final nucleon energy
              Ep = (M_TARG2 + M_scat2 - t)/(2.*M_TARG);   // assumes stationary target!
            }
          
          // Check that the produced particle is in a possible range:
          if(nup > nupmax || nup < nupmin)
            {
              cout<<"event "<<ievt<<"; nup out of range "<<
		nupmin<<" "<<nup<<" "<<nupmax<<endl;
              out0<<"event "<<ievt<<"; nup out of range "<<
		nupmin<<" "<<nup<<" "<<nupmax<<endl;
              
              outside_klim++;
	      outside_klim_nup++;              
	      dv->init_dvcs();
              init_event();
              init_track();
              continue;
            }
          
          //final nucleon momentum
          Pp = sqrt(sqr(Ep) - M_scat2);
          if(Ep > Epmax || Ep < Epmin)
            {
              cout<<"event "<<ievt<<"; Ep out of range "<<
		Epmin<<" "<<Ep<<" "<<Epmax<<endl;
              out0<<"event "<<ievt<<"; Ep out of range "<<
		Epmin<<" "<<Ep<<" "<<Epmax<<endl;
              
              outside_klim++;
	      outside_klim_Ep++;
              dv->init_dvcs();
              init_event();
              init_track();
              continue;
            }
          
          //polar angle of produced particle wrt virt. phot. frame
          if (trk.Process != 1)     // DVCS
            {
              cosThetaqqp = (t + 2.*nu*nup + Q2)/(2.*q*nup);       // assumes stationary target and real photon in final state!
            }
          else if (trk.Process == 1)     // meson-production
            {
              cosThetaqqp = (t + Q2 - m_ms(Ims,2) + 2.*nu*nup) / (2.*q*qp);   // assumes stationary target and meson in final state!
            }
          
          if(abs(cosThetaqqp) >= 1.)
            {
              cout<<"event "<<ievt<<
		"; cosThetaqqp larger than one; cosThetaqqp = "<<cosThetaqqp<< endl;
              cosThetaqqp = sign(cosThetaqqp);
            }
          
          Thetaqqp =  acos(cosThetaqqp);//[0,PI]
          
          if(Thetaqqp > pi(1))
            cout<<"Thetaqqp larger than pi; cosThetaqqp = "<<cosThetaqqp<< endl;
          
          //generate produced particle azimuthal angle wrt virt. phot. frame
          double rndmphi = rndm.Rndm();
          Phiqqp = 2.*pi(1)*rndmphi;
          
          //produced particle 3-momentum vector wrt virt. phot. frame
          vqp1 = vect->fill_vect(qp, Thetaqqp, Phiqqp);
          
          //virt. photon 3-momentum vector wrt virt. phot. frame
          vq2 = vect->fill_vect(q, 0., 0.);
          
          //final nucleon 3-momentum vector wrt virt. phot. frame
          vPp1 = vect->diff_vect(vq2, vqp1);                         // assumes stationary target!
          
          //polar angle of final nucleon wrt virt. phot. frame
          ThetaqPp     = get_theta(vPp1.at(2), Pp);
          
          //azimuthal angle of final nucleon wrt virt. phot. frame
          PhiqPp       = get_phi(vPp1.at(0), vPp1.at(1));
          
          //rotate the produced particle 3-momentum vector from
          //the virt. phot. frame to the inc. elec. frame
          vqp = vect->invrot_vect(vqp1, Thetakq, Phikq);
          
          //polar angle of produced particle wrt inc. elec. frame
          Thetakqp     = get_theta(vqp.at(2), qp);
          
          //azimuthal angle of produced particle wrt inc. elec. frame
          Phikqp   = get_phi(vqp.at(0), vqp.at(1));
          
          //final nucleon 3-momentum vector wrt inc. elec. frame
          vPp  = vect->diff_vect(vect->add_vect(vP,vq), vqp);
          
          //polar angle of final nucleon wrt inc. elec. frame
          ThetakPp     = get_theta(vPp.at(2), Pp);
          
          //azimuthal angle of final nucleon wrt inc. elec. frame
          PhikPp       = get_phi(vPp.at(0), vPp.at(1));
          
          //Angle between scat. elec. and produced particle (from dot product of kp and qp 3mom vectors)
          cosThetakpqp = cos(Phikqp)*sin(Thetakqp)*cos(Phikkp)*sin(Thetakkp) +
            sin(Phikqp)*sin(Thetakqp)*sin(Phikkp)*sin(Thetakkp) +
            cos(Thetakqp)*cos(Thetakkp);
          Thetakpqp    = acos(cosThetakpqp);
          
        }
      
      // In the case of a moving target nucleon, a separate function (fmotion) calculates the final state particles and other output:
      // FERMI MOTION ROUTINE TO CALCULATE FINAL STATE PARTICLES FOR A MOVING TARGET

      if (iApZ == 3 || iApZ == 8 || iApZ == 21)
        {

	  double Output[7];
	  
	  rndmt = rndm.Rndm();
          double rndmphi2 = rndm.Rndm();
          int ims_flag = -1;  // default value, set that for DVCS (produced real photon)                                                                                      
          if(trk.Process == 1) ims_flag = Ims;   // will carry the ID value for the meson if it's a DVMP process                                                              
          if(fmotion(rndmt, rndmphi2, Ipn, ims_flag, q, nu, W2, P, ThetakP, PhikP, Thetakq, Phikq, ro->get_ftmin(), ro->get_ftmax(), Output) == -1)
            {
              //cout<<"Error in Fermi motion subroutine"<<endl;
              //out0<<"Error in Fermi motion subroutine"<<endl;
	      
              outside_klim++;
              outside_klim_fm++;
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
	  
	  // The calculation on basis of max and min possible energy of produced particle, which is used in the stationary target case,
	  // doesn't apply here -- but t was generated in the range possible, given the energy of the photon or meson, which is well constrained in the CM frame. 
	  // So just check that it also falls within the t-range selected for the file (but this condition was also applied when it was generated!):
	  
          if(t > ro->get_ftmax() || t < ro->get_ftmin())
            {
	      cout << "Non-stationary target " << endl;
	      cout<<"event "<<ievt<<"; t out of range "<<
	      	ro->get_ftmin()<<" "<<t<<" "<< ro->get_ftmax() <<endl;
	      out0<<"event "<<ievt<<"; t out of range "<<
	      	ro->get_ftmax() <<" "<<t<<" "<< ro->get_ftmax() <<endl;
	      
              outside_klim++;
              outside_klim_t2++;
              dv->init_dvcs();
              init_event();
              init_track();
              continue;
            }
	  
	  if (trk.Process != 1) qp = nup;     // DVCS, assumes produced particle is a photon
	  else if (trk.Process == 1) qp = (nup*nup - m_ms(Ims,2));   // produced particle is a meson
	  
	  // The fmotion routine assumes that the initial and scattered nucleons are the same. This is easy to change -- you need to edit only the lines which calculate the scattered (recoil) nucleon energy at the end of the routine (two lines there, search for M_TARG).
	  Pp = sqrt(sqr(Ep) - M_scat2);            // Scattered nucleon might not be the same as the target one!
	  
          vqp = vect->fill_vect(qp, Thetakqp, Phikqp);
	  
          vPp = vect->fill_vect(Pp, ThetakPp, PhikPp);
	  
          //rotate the produced particle 3-momentum vector from
          //the inc. elec. frame to the virt. phot. frame                                                                                                                     
          vqp1 = vect->rot_vect(vqp, Thetakq, Phikq);	  
	  
          //polar angle of produced particle wrt virt. phot. frame
          Thetaqqp     = get_theta(vqp1.at(2),qp);
	  
          //azimuthal angle of produced particle wrt virt. phot. frame
          Phiqqp       = get_phi(vqp1.at(0),vqp1.at(1));
	  
          //final nucleon 3-momentum vector wrt virt. phot. frame:
	  // rotate the nucleon vector from the inc. electron. frame to the virt. phot. frame:
	  vPp1 = vect->rot_vect(vPp, Thetakq, Phikq);
	  
          //polar angle of final nucleon wrt virt. phot. frame                                                                                                                
          ThetaqPp     = get_theta(vPp1.at(2),Pp);
	  
          //azimuthal angle of final nucleon wrt virt. phot. frame                                                                                                            
          PhiqPp       = get_phi(vPp1.at(0), vPp1.at(1));
	  
          //Angle between scat. elec. and produced particle
          cosThetakpqp = cos(Phikqp)*sin(Thetakqp)*cos(Phikkp)*sin(Thetakkp) +
            sin(Phikqp)*sin(Thetakqp)*sin(Phikkp)*sin(Thetakkp) +
            cos(Thetakqp)*cos(Thetakkp);
          Thetakpqp    = acos(cosThetakpqp);
	  
	} // END OF FERMI MOTION SECTION FOR A MOVING TARGET
      
      ycol  = get_ycol(t, xbj, Q2);

      /**** Caluclate the kinematic tmin here and apply a final check on that: ****/
      // Calculate xi:
      vector<double> vPdiff = vect->init_vect();
      vector<double> vPsum = vect->init_vect();
      for (int i=0; i<3; i++){
	vPdiff.at(i) = vP.at(i) - vPp.at(i);     // p - p', 3mom
	vPsum.at(i) = vP.at(i) + vPp.at(i);      // p + p', 3mom. Gives direction of axis for light-cone co-ordinates
      }
      double Pdiff = sqrt(pow(vPdiff.at(0),2) + pow(vPdiff.at(1),2) + pow(vPdiff.at(2),2));  // magnitude of p - p' 3mom
      double Psum = sqrt(pow(vPsum.at(0),2) + pow(vPsum.at(1),2) + pow(vPsum.at(2),2));  // magnitude of p - p' 3mom
      double cosThetaPdiffPsum = (vPdiff.at(0)*vPsum.at(0) + vPdiff.at(1)*vPsum.at(1) + vPdiff.at(2)*vPsum.at(2)) / (Pdiff*Psum);   // from dot product of vPdiff and vPsum
      double PdiffZ = Pdiff * cosThetaPdiffPsum;  // Z-component of Pdiff along the Psum axis

      double xi = 0.; // xi = (p-p')^+ / (p+p')^+

      // To find whether you take the plus or minus component in calculation of xi, check whether p^+ or p'^+ is greater:
      double cosThetaPPsum = (vP.at(0)*vPsum.at(0) + vP.at(1)*vPsum.at(1) + vP.at(2)*vPsum.at(2)) / (P*Psum);
      double cosThetaPpPsum = (vPp.at(0)*vPsum.at(0) + vPp.at(1)*vPsum.at(1) + vPp.at(2)*vPsum.at(2)) / (Pp*Psum);
      if ((E + P*cosThetaPPsum) > (Ep + Pp*cosThetaPpPsum)){
	xi = (E-Ep + PdiffZ) / (E+Ep + vPsum.at(2));   // + component is E + pz.
      }
      else xi = (E-Ep - PdiffZ) / (E+Ep - vPsum.at(2));

      tmin = (-4 * xi*xi * M_TARG2) / (1 - xi*xi);   // kinematic tmin, actually tmax given that t is negative

      //      cout << "tmin calculated: " << tmin << " and xi " << xi << " and P " << P << " and Pp " << Pp << " and t: " << t << endl;
      
      if(t > tmin)    // t is negative, tmin is actually a max (min absolute value)
	{
	  //  cout<<"event "<<ievt<<"; t out of range "<<
	  //    tmin<<" "<<t<<" "<<endl;
	  //  out0<<"event "<<ievt<<"; t out of range "<<
	  //   tmin <<" "<<t<<" "<<endl;

	  outside_klim++;
	  outside_klim_tmin++;
	  dv->init_dvcs();
	  init_event();
	  init_track();
	  continue;
	}
      
      /**************************/
      
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
      if(trk.Process == 0)  // DVCS 
	{ 
	  if(dv->phot_xsec(ro, ff, xbj, y, Q2, t, Phi_b) != 0 )
	    {
	      //cout<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;
	      //out0<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;
	      
	      outside_klim++;
	      outside_klim_xsec++;
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
      else if(trk.Process == 1)   // meson-production
	{
	  ds_ms = ms_xsec(ro, xbj, Q2, t, Phi_b, Ims, tmin);    // this uses the kinematic tmin! Calculated higher-up.
	  xsec =  ds_ms;
	  
	  if(get_ms(ro, xbj, Q2, nu, nup, t, Phikkp, Phiqqp) != 0)
	    {
	      cout<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;
	      out0<<"event "<<ievt<<"; Error in xsec calculation. "<<endl;
	      
	      outside_klim++;
	      outside_klim_xsec++;
	      dv->init_dvcs();
	      init_event();
	      init_track();
	      continue;
	    }
	}
      else
	{
	  out0<<"You Should select a Process"<<endl;
	}
      
      //fill track
      //be carefull here
      //get_ms() returns the particles coming out from the pi0/eta decay
      //the number of particles, stocked in lujets_.n, is not constant
      
      int id_inc_el = 0;
      int id_targ_nuc = 1;
      int id_scat_el = 2;
      int id_virt_ph = 3;
      int id_recoil = 4;
      int id_scat_nuc = 5;
      int id_out_ph = 6;
      trk.Ntracks   = 7;
      
      //DVCS on H : no recoil
      if(iApZ == 2) 
	{
	  trk.Ntracks = 6;
	  id_scat_nuc = 4;
	  id_out_ph = 5;
	} 

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
      trk.Type[id_inc_el]   = elec_id();
      trk.Charge[id_inc_el] = elec_ch();
      trk.Px[id_inc_el]      = vk.at(0);
      trk.Py[id_inc_el]      = vk.at(1);
      trk.Pz[id_inc_el]      = vk.at(2);
      trk.P[id_inc_el]       = k;
      trk.E[id_inc_el]       = Ee;
      trk.Theta[id_inc_el]   = 0.;
      trk.Phi[id_inc_el]     = 0.;
    
      //target nucleon
      trk.Type[id_targ_nuc]    = targ_id(Ipn);
      trk.Charge[id_targ_nuc]  = targ_ch(Ipn);
      trk.Px[id_targ_nuc]      = vP.at(0);
      trk.Py[id_targ_nuc]      = vP.at(1);
      trk.Pz[id_targ_nuc]      = vP.at(2);
      trk.P[id_targ_nuc]       = P;
      trk.E[id_targ_nuc]       = E;
      trk.Theta[id_targ_nuc]   = ThetakP;
      trk.Phi[id_targ_nuc]     = PhikP;

      //scat. elec.
      trk.Type[id_scat_el]    = elec_id();
      trk.Charge[id_scat_el]  = elec_ch();
      trk.Px[id_scat_el]      = vkp.at(0);
      trk.Py[id_scat_el]      = vkp.at(1);
      trk.Pz[id_scat_el]      = vkp.at(2);
      trk.P[id_scat_el]       = kp;
      trk.E[id_scat_el]       = Eep;
      trk.Theta[id_scat_el]   = Thetakkp;
      trk.Phi[id_scat_el]     = Phikkp;

      //virt. photon
      trk.Type[id_virt_ph]    = phot_id();
      trk.Charge[id_virt_ph]  = phot_ch();
      trk.Px[id_virt_ph]      = vq.at(0);
      trk.Py[id_virt_ph]      = vq.at(1);
      trk.Pz[id_virt_ph]      = vq.at(2);
      trk.P[id_virt_ph]       = q;
      trk.E[id_virt_ph]       = nu;
      trk.Theta[id_virt_ph]   = Thetakq;
      trk.Phi[id_virt_ph]     = Phikq;
    
      if(iApZ!=2)//if(iApZ_stored!=2)
	{
	  //recoil nucleon/nucleus
	  trk.Type[id_recoil]    = RECO_ID;
	  trk.Charge[id_recoil]  = RECO_CH;
	  trk.Px[id_recoil]      = vPr.at(0);
	  trk.Py[id_recoil]      = vPr.at(1);
	  trk.Pz[id_recoil]      = vPr.at(2);
	  trk.P[id_recoil]       = Pr;
	  trk.E[id_recoil]       = Er;
	  trk.Theta[id_recoil]   = ThetakPr;
	  trk.Phi[id_recoil]     = PhikPr;
	}

      //scat. nucleon
      trk.Type[id_scat_nuc]    = targ_id(Ipn);
      trk.Charge[id_scat_nuc]  = targ_ch(Ipn);
      trk.Px[id_scat_nuc]      = vPp.at(0);
      trk.Py[id_scat_nuc]      = vPp.at(1);
      trk.Pz[id_scat_nuc]      = vPp.at(2);
      trk.P[id_scat_nuc]       = Pp;
      trk.E[id_scat_nuc]       = Ep;
      trk.Theta[id_scat_nuc]   = ThetakPp;
      trk.Phi[id_scat_nuc]     = PhikPp;

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

	  lujets_cc.K[1][id_inc_el] = elec_id();    //inc. electron
	  lujets_cc.K[1][id_targ_nuc] = targ_id(Ipn); //proton target
	  lujets_cc.K[1][id_scat_el] = elec_id();    //scat. electron
	  lujets_cc.K[1][id_virt_ph] = phot_id();    //virtual photon
	  //if(trk.Process != 0 || iApZ_stored!=2) lujets_cc.K[1][4] = RECO_ID;      //recoil proton
	  if(iApZ!=2) lujets_cc.K[1][4] = RECO_ID;
	  lujets_cc.K[1][id_scat_nuc] = targ_id(Ipn); //scat. proton

	  for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[2][ii] = 0;
	  for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[3][ii] = 0;
	  for(int ii=0; ii<lujets_cc.N; ii++) lujets_cc.K[4][ii] = 0;
 
	  if(trk.Process == 0)
	    {
	      lujets_cc.K[0][id_out_ph] = 1;
	      lujets_cc.K[1][id_out_ph] = phot_id();   //real photon
	      lujets_cc.K[2][id_out_ph] = 0;
	      lujets_cc.K[3][id_out_ph] = 0;
	      lujets_cc.K[4][id_out_ph] = 0;
	    }
	  else if(trk.Process == 1)
	    {
	      for(int in=0; in<lujets_.n; in++)
		{
		  for(int jn=0; jn<5; jn++)
		    {
		      lujets_cc.K[jn][NNT + in] = lujets_.k[jn][in];
		    }
		}
	    }

	  lujets_cc.P[0][id_inc_el] = vk.at(0);
	  lujets_cc.P[0][id_targ_nuc] = vP.at(0);
	  lujets_cc.P[0][id_scat_el] = vkp.at(0);
	  lujets_cc.P[0][id_virt_ph] = vq.at(0);
	  //if(trk.Process != 0 || iApZ_stored!=2)lujets_cc.P[0][4] = vPr.at(0);
	  if (iApZ != 2) lujets_cc.P[0][4] = vPr.at(0);
	  lujets_cc.P[0][id_scat_nuc] = vPp.at(0);

	  lujets_cc.P[1][0] = vk.at(1);
	  lujets_cc.P[1][1] = vP.at(1);
	  lujets_cc.P[1][id_scat_el] = vkp.at(1);
	  lujets_cc.P[1][id_virt_ph] = vq.at(1);
	  //if(trk.Process != 0 || iApZ_stored!=2) lujets_cc.P[1][4] = vPr.at(1);
	  if (iApZ != 2) lujets_cc.P[1][4] = vPr.at(1);
	  lujets_cc.P[1][id_scat_nuc] = vPp.at(1);

	  lujets_cc.P[2][0] = vk.at(2);
	  lujets_cc.P[2][1] = vP.at(2);
	  lujets_cc.P[2][id_scat_el] = vkp.at(2);
	  lujets_cc.P[2][id_virt_ph] = vq.at(2);
	  //if(trk.Process != 0 || iApZ_stored!=2)lujets_cc.P[2][4] = vPr.at(2);
	  if (iApZ != 2) lujets_cc.P[2][4] = vPr.at(2);
	  lujets_cc.P[2][id_scat_nuc] = vPp.at(2);

	  lujets_cc.P[3][0] = Ee;
	  lujets_cc.P[3][1] = E;
	  lujets_cc.P[3][id_scat_el] = Eep;
	  lujets_cc.P[3][id_virt_ph] = nu;
	  //if(trk.Process != 0 || iApZ_stored!=2)lujets_cc.P[3][4] = Er;
	  if (iApZ != 2) lujets_cc.P[3][4] = Er;
	  lujets_cc.P[3][id_scat_nuc] = Ep;

	  lujets_cc.P[4][0] = m_elec(1);
	  lujets_cc.P[4][1] = m_targ(Ipn,1);
	  lujets_cc.P[4][id_scat_el] = m_elec(1);
	  lujets_cc.P[4][id_virt_ph] = Q2;
	  //if(trk.Process != 0 || iApZ_stored!=2)lujets_cc.P[4][4] = M_RECO;
	  if (iApZ != 2) lujets_cc.P[4][4] = M_RECO;
	  lujets_cc.P[4][id_scat_nuc] = m_targ(Ipn,1);

	  if(trk.Process == 0)
	    {
	      lujets_cc.P[0][id_out_ph] = vqp.at(0);
	      lujets_cc.P[1][id_out_ph] = vqp.at(1);
	      lujets_cc.P[2][id_out_ph] = vqp.at(2);
	      lujets_cc.P[3][id_out_ph] = nup;
	      lujets_cc.P[4][id_out_ph] = m_phot(1);
	    }
	  else if(trk.Process == 1)
	    {
	      for(int in=0; in<lujets_.n; in++)
		{
		  for(int jn=0; jn<5; jn++)
		    {
		      lujets_cc.P[jn][NNT + in] = lujets_.p[jn][in];
		    }
		}
	    }
	  
	  plu[id_inc_el]              = elec_ch();    
	  plu[id_targ_nuc]            = targ_ch(Ipn);
	  plu[id_scat_el]             = elec_ch();
	  plu[id_virt_ph]             = phot_ch();
	  if(iApZ!=2) plu[id_recoil]  = RECO_CH;
	  plu[id_scat_nuc]            = targ_ch(Ipn);

	  if(trk.Process == 0)
	    {
	      plu[id_out_ph]  = phot_ch();
	    }
	  else if(trk.Process == 1)
	    {
	      for(int in=0; in<lujets_.n; in++)
		{
		  plu[in] = trk.Charge[NNT + in];
		}
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
	      dump_file(ro->get_fMode(), ro->get_fProc(),xsec,iApZ);
	    }
	  else
	    {
	      counter_rejected++;
	      kr= 0;
	    }
	  trk.kr = kr;

	}
      if(ro->get_fNtp())// && trk.kr==1) 
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

      ievt++;
      trk.Evt_ID  = ievt;
      Nevts_per_Ntup++;

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
  cout<<"Number of hit neutron events kept   " << counter_kept_neut   <<endl; 
  cout<<"Number of hit proton events kept   " << counter_kept_prot   <<endl;
  cout<<"Number of events Out of kin. limits   "<< outside_klim   <<endl;
  cout<<"Number of events Out of t limits   "<< outside_klim_t   <<endl;
  cout<<"Number of events Out of t limits 2  "<< outside_klim_t2   <<endl;
  cout<<"Number of events Out of t limits (tmin)  "<< outside_klim_tmin   <<endl;
  cout<<"Number of events Out of nu limits   "<< outside_klim_nu   <<endl;
  cout<<"Number of events Out of W2 limits   "<< outside_klim_W2   <<endl;
  cout<<"Number of events Out of Q2 limits   "<< outside_klim_Q2   <<endl;
  cout<<"Number of events Out of xb limits   "<< outside_klim_xb   <<endl;
  cout<<"Number of events Out of y limits   "<< outside_klim_y   <<endl;
  cout<<"Number of events Out of nup limits   "<< outside_klim_nup   <<endl;
  cout<<"Number of events Out of Ep limits   "<< outside_klim_Ep   <<endl;
  cout<<"Number of events Out of limits within the fmotion routine  "<< outside_klim_fm   <<endl;
  cout<<"Number of events Out of xsec limits   "<< outside_klim_xsec   <<endl;
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
  out0<<"Number of hit neutron events kept   " << counter_kept_neut   <<endl;
  out0<<"Number of hit proton events kept   " << counter_kept_prot   <<endl;
  out0<<"Number of events Out of kin. limits   "<< outside_klim   <<endl;
  out0<<"Number of events Out of t limits   "<< outside_klim_t   <<endl;
  out0<<"Number of events Out of t limits 2  "<< outside_klim_t2   <<endl;
  out0<<"Number of events Out of t limits (tmin)  "<< outside_klim_tmin   <<endl;
  out0<<"Number of events Out of nu limits   "<< outside_klim_nu   <<endl;
  out0<<"Number of events Out of W2 limits   "<< outside_klim_W2   <<endl;
  out0<<"Number of events Out of Q2 limits   "<< outside_klim_Q2   <<endl;
  out0<<"Number of events Out of xb limits   "<< outside_klim_xb   <<endl;
  out0<<"Number of events Out of y limits   "<< outside_klim_y   <<endl;
  out0<<"Number of events Out of nup limits   "<< outside_klim_nup   <<endl;
  out0<<"Number of events Out of Ep limits   "<< outside_klim_Ep   <<endl;
  out0<<"Number of events Out of limits within the fmotion routine  "<< outside_klim_fm   <<endl;
  out0<<"Number of events Out of xsec limits   "<< outside_klim_xsec   <<endl;
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

