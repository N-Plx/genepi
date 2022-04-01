#include"dvcs_vars.h"
#include "genepi.h"
#include "inl_funcs.h"

double (NuclFF::*FF1[2])(double, int) = {NULL};
double (NuclFF::*FF2[2])(double, int) = {NULL};


DVCS::DVCS()
{};


DVCS::~DVCS()
{};


void DVCS::init_dvcs()
{
  Re_E    = -1000.;
  Im_E    = -1000.;
  Re_Et   = -1000.;
  Im_Et   = -1000.;

  Re_H    = -1000.;
  Im_H    = -1000.;
  Re_Ht   = -1000.;
  Im_Ht   = -1000.;

  hc0_BH      = 0.;
  hc1_BH      = 0.; 
  hc2_BH      = 0.;

  hc0_BH_LP   = 0.;
  hc1_BH_LP   = 0.; 

  hc0_DVCS    = 0.;
  hc1_DVCS    = 0.;
  hs1_DVCS    = 0.;

  hc0_DVCS_LP = 0.;
  hc1_DVCS_LP = 0.;
  hs1_DVCS_LP = 0.;

  hc0_Int     = 0.;
  hc1_Int     = 0.;
  hc2_Int     = 0.;
  hs1_Int     = 0.;
  hs2_Int     = 0.;

  hc0_Int_LP  = 0.;
  hc1_Int_LP  = 0.;
  hc2_Int_LP  = 0.;
  hs1_Int_LP  = 0.;
  hs2_Int_LP  = 0.;
}


int DVCS::read_gpds(string gpd_tbl)
{
  ifstream in(gpd_tbl.c_str());
  
  //real part
  for(int j=0; j<5; j++)
    {
      for(int i=0; i<MAXSK; i++)
	{
	  in>>re_Eu[i][j]>>re_Ed[i][j]>>re_Etu[i][j]>>re_Etd[i][j]>>
	    re_Hu[i][j]>>re_Hd[i][j]>>re_Htu[i][j]>>re_Htd[i][j];
	}
    }
  
  //imaginary part
  for(int j=0; j<5; j++)
    {
      for(int i=0; i<MAXSK; i++)
	{
	  in>>im_Eu[i][j]>>im_Ed[i][j]>>im_Etu[i][j]>>im_Etd[i][j]>>
	    im_Hu[i][j]>>im_Hd[i][j]>>im_Htu[i][j]>>im_Htd[i][j];
	}
    }
  
  //real part
  for(int k=0; k<2; k++)
    {
      for(int j=0; j<21; j++)
	{
	  for(int i=0; i<MAXSK-1; i++)
	    {
	      in>>re_HuNF[i][j][k]>>re_HdNF[i][j][k];  
	    }
	}
    }
  
  //imaginary part
  for(int k=0; k<2; k++)
    {
      for(int j=0; j<21; j++)
	{
	  for(int i=0; i<MAXSK-1; i++)
	    {
	      in>>im_HuNF[i][j][k]>>im_HdNF[i][j][k];  
	    }
	}
    }
  
  return 0;
}


int DVCS::eval_cffs(ReadOptFile *ro, NuclFF *ff, double skew, double t)
{
  double sk, dsk, skD;
  double tt, dtt, ttD;
  
  double F1u, F1d, F2u, F2d;
  double gA, gA0, gAu, gAd, Fpi;
  
  double Re_Hu, Re_Htu, Re_Eu, Re_Etu;
  double Re_Hd, Re_Htd, Re_Ed, Re_Etd;

  double Im_Hu, Im_Htu, Im_Eu, Im_Etu;
  double Im_Hd, Im_Htd, Im_Ed, Im_Etd;

  int    isk;
  int    itt;
  int    iGPD1(-1), iGPD2(-1);
  int    iGPD0(ro->get_fMgpd()), Ipn(ro->get_fIpn());
  double M_TARG2(m_targ(Ipn,2.));

  //skewmin = 0.01; skewmin = 1.0
  dsk  = (log10(1.00) - log10(0.01))/MAXSK;
  sk   = (log10(skew) - log10(0.01))/dsk;
  isk  = (int)sk;
  skD  = sk - (double)isk;
  //if(isk <= 10) cout<<"sk = "<<sk<<"; skew = "<<skew <<"  "<<isk<<endl;
  
  //|t_min| = 0.01; |t_max| = 1.20;
  //only needed if iGPD0=5,6 H in proton
  dtt  = (log10(1.20) - log10(0.01))/20.;
  tt   = (log10(abs(t))- log10(0.01))/dtt;
  itt  = (int)tt;
  ttD  = tt - (double)itt;
  //cout<<"tt = "<<tt<<"; t = "<<t <<" "<<itt<<endl;
  
  iGPD1 = iGPD0;
  if( iGPD0 == 5 || iGPD0 == 6 )
    {
      iGPD1 = iGPD0 - 4;
      iGPD2 = iGPD0 - 5;
    }
  else if(iGPD0 > 6)
    {
      cout<<"Unknown GPD Model"<<endl;
      return -1;
    }
  
  Re_Hu  = re_Hu[isk][iGPD1]*(1. - skD) + re_Hu[isk+1][iGPD1]*skD;
  Im_Hu  = im_Hu[isk][iGPD1]*(1. - skD) + im_Hu[isk+1][iGPD1]*skD;
  
  Re_Hd  = re_Hd[isk][iGPD1]*(1. - skD) + re_Hd[isk+1][iGPD1]*skD;
  Im_Hd  = im_Hd[isk][iGPD1]*(1. - skD) + im_Hd[isk+1][iGPD1]*skD;

  Re_Htu = re_Htu[isk][iGPD1]*(1. - skD) + re_Htu[isk+1][iGPD1]*skD;
  Im_Htu = im_Htu[isk][iGPD1]*(1. - skD) + im_Htu[isk+1][iGPD1]*skD;

  Re_Htd = re_Htd[isk][iGPD1]*(1. - skD) + re_Htd[isk+1][iGPD1]*skD;
  Im_Htd = im_Htd[isk][iGPD1]*(1. - skD) + im_Htd[isk+1][iGPD1]*skD;

  Re_Eu  = re_Eu[isk][iGPD1]*(1. - skD)  + re_Eu[isk+1][iGPD1]*skD;
  Im_Eu  = im_Eu[isk][iGPD1]*(1. - skD)  + im_Eu[isk+1][iGPD1]*skD;

  Re_Ed  = re_Ed[isk][iGPD1]*(1. - skD)  + re_Ed[isk+1][iGPD1]*skD;
  Im_Ed  = im_Ed[isk][iGPD1]*(1. - skD)  + im_Ed[isk+1][iGPD1]*skD;

  Re_Etu = re_Etu[isk][iGPD1]*(1. - skD) + re_Etu[isk+1][iGPD1]*skD;
  Im_Etu = im_Etu[isk][iGPD1]*(1. - skD) + im_Etu[isk+1][iGPD1]*skD;

  Re_Etd = re_Etd[isk][iGPD1]*(1. - skD) + re_Etd[isk+1][iGPD1]*skD;
  Im_Etd = im_Etd[isk][iGPD1]*(1. - skD) + im_Etd[isk+1][iGPD1]*skD;

  if(iGPD0 > 4)
    {
      Re_Hu = (re_HuNF[isk][itt][iGPD2]*(1. - skD)   + re_HuNF[isk+1][itt][iGPD2]*skD)*(1. - ttD) + 
	(re_HuNF[isk][itt+1][iGPD2]*(1. - skD) + re_HuNF[isk+1][itt][iGPD2]*skD)*ttD;
      
      Im_Hu = (im_HuNF[isk][itt][iGPD2]*(1. - skD)   + im_HuNF[isk+1][itt][iGPD2]*skD)*(1. - ttD) + 
	(im_HuNF[isk][itt+1][iGPD2]*(1. - skD) + im_HuNF[isk+1][itt][iGPD2]*skD)*ttD;
      
      Re_Hd = (re_HdNF[isk][itt][iGPD2]*(1. - skD)   + re_HdNF[isk+1][itt][iGPD2]*skD)*(1. - ttD) + 
	(re_HdNF[isk][itt+1][iGPD2]*(1. - skD) + re_HdNF[isk+1][itt][iGPD2]*skD)*ttD;
      
      Im_Hd = (im_HdNF[isk][itt][iGPD2]*(1. - skD)   + im_HdNF[isk+1][itt][iGPD2]*skD)*(1. - ttD) + 
	(im_HdNF[isk][itt+1][iGPD2]*(1. - skD) + im_HdNF[isk+1][itt][iGPD2]*skD)*ttD;
    }
  
  F1u = 1.*ff->get_FF1(0) + 2.*ff->get_FF1(1);
  F1d = 2.*ff->get_FF1(0) + 1.*ff->get_FF1(1);
  F2u = 1.*ff->get_FF2(0) + 2.*ff->get_FF2(1);
  F2d = 2.*ff->get_FF2(0) + 1.*ff->get_FF2(1);
  
  gA  = 1.267/sqr(1. - t/0.84);
  gA0 = 0.6*gA;
  gAu = 0.5*(gA + gA0)/(0.8*1.267);
  gAd = 0.5*(gA - gA0)/(0.2*1.267);
  Fpi = 4*1.267*M_TARG2/(m_pion(2) - t);
  
  if(Ipn == 1)
    {
      Re_H  = (4./9.)*(F1u/2.)*Re_Hu + (1./9.)*F1d*Re_Hd;
      Im_H  = (4./9.)*(F1u/2.)*Im_Hu + (1./9.)*F1d*Im_Hd;
      
      Re_Ht = (4./9.)*gAu*Re_Htu     + (1./9.)*gAu*Re_Htd;
      Im_Ht = (4./9.)*gAd*Im_Htu     + (1./9.)*gAd*Im_Htd;
      
      Re_E  = (4./9.)*(F2u/2.)*Re_Eu + (1./9.)*F2d*Re_Ed;
      Im_E  = (4./9.)*(F2u/2.)*Im_Eu + (1./9.)*F2d*Im_Ed;
      
      Re_Et = (4./9.)*Fpi*Re_Etu     + (1./9.)*Fpi*Re_Etd;
      Im_Et = (4./9.)*Fpi*Im_Etu     + (1./9.)*Fpi*Im_Etd;
      
      if(iGPD0 > 4)
	{
	  Re_H = (4./9.)*Re_Hu + (4./9.)*Re_Hd;
	  Im_H = (4./9.)*Im_Hu + (4./9.)*Im_Hd;
	}
    }
  else if(Ipn == 0)
    {
      Re_H  = (4./9.)*(F1d/2.)*Re_Hd + (1./9.)*F1u*Re_Hu;
      Im_H  = (4./9.)*(F1d/2.)*Im_Hd + (1./9.)*F1u*Im_Hu;
      
      Re_Ht = (4./9.)*gAd*Re_Htd     + (1./9.)*gAu*Re_Htu;
      Im_Ht = (4./9.)*gAd*Im_Htd     + (1./9.)*gAu*Im_Htu;
      
      Re_E  = (4./9.)*(F2d/2.)*Re_Ed + (1./9.)*F2u*Re_Eu;
      Im_E  = (4./9.)*(F2d/2.)*Im_Ed + (1./9.)*F2u*Im_Eu;
      
      Re_Et = (4./9.)*Fpi*Re_Etd     + (1./9.)*Fpi*Re_Etu;
      Im_Et = (4./9.)*Fpi*Im_Etd     + (1./9.)*Fpi*Im_Etu;
    }
  
  return 0;
}


/*==============================================================================
  Calculation of the BH xsec off the nucleon

  inputs:
   xbj    Bjorken x defined on the nucleon at rest
   y
   Q2     virtual photon mass with sign Q2 > 0 (GeV2 units)
   t      4-momentum transfer to the nucleon with sign t < 0 (GeV2 units)
   Phi_b   out-of-plane angle in the Belitsky frame (rd unit)

  output:
   BHxsec BH xsection in the given kin. variables and phase space volume index
==============================================================================*/
int DVCS::bh_xsec(ReadOptFile *ro, double *FF, double x, double y, double Q2, double t, double Phi_b)
{
  int    Ipn = ro->get_fIpn();
  double M_TARG(m_targ(Ipn,1.));
  double M_TARG2(m_targ(Ipn,2.));
  double Phi_bs = pi(1) - Phi_b;
  double K, K2, J, P1, P2;
  double x2     = x*x;
  double y2     = y*y;
  double eps    = 2.*x*M_TARG/sqrt(Q2);
  double eps2   = eps*eps;
  double tau    = t/(4.*M_TARG2);
  double BHfact, XSfact;
 
  K2 = get_K2(Ipn, x, y, Q2, t);
  J  = get_J(Ipn, x, y, Q2, t);
  K  = sqrt(K2);
  P1 = -(J + 2.*K*cos(Phi_bs))/(y*(1.+eps2));
  P2 = 1.+t/Q2 - P1;

  double F1         = FF[0];//FF1[Ipn];
  double F2         = FF[1];//FF2[Ipn];
  double F12mtauF22 = F1*F1-tau*F2*F2;

  //unpol target
  double c00 = 8.*K2*((Q2/t)*(2.+3.*eps2)*F12mtauF22 + 2.*x2*sqr(F1+F2));
  double c01 = sqr(2.-y)*(2.+eps2)*((eps2*Q2/t)*sqr(1.+t/Q2) + 4.*(1.-x)*(1.+x*t/Q2))*F12mtauF22;
  double c02 = 4.*sqr(2.-y)*x2*(x + (1.-x+eps2/2.)*sqr(1.-t/Q2) - x*(1.-2.*x)*sqr(t/Q2))*sqr(F1+F2);
  double c03 = 8.*(1.+eps2)*(1.-y-y2*eps2/4.)*(2.*eps2*(1.-tau)*F12mtauF22 - x2*sqr(1.-t/Q2)*sqr(F1+F2));
  hc0_BH = c00 + c01 + c02 + c03;

  hc1_BH = 8.*K*(2.-y)*((4.*x2*M_TARG2/t-2.*x-eps2)*F12mtauF22 + 2.*x2*(1.-(1.-2.*x)*t/Q2)*sqr(F1+F2));

  hc2_BH = 8.*x2*K2*(F12mtauF22/tau + 2.*sqr(F1+F2));

  //long. pol. target
  double c00_LP = 8.*x*y*(2.-y)*sqrt(1.+eps2)*(F1+F2)/(1.-tau);
  double c01_LP = (x*(1.-t/Q2)/4.-tau/2.)*(2.-x-2.*t*sqr(1.-x)/Q2+eps2*(1.-t/Q2)-x*(1.-2.*x)*sqr(t/Q2))*(F1+F2);
  double c02_LP = (1.-(1.-x)*t/Q2)*(x2*M_TARG2*sqr(1.+t/Q2)/t+(1.-x)*(1.+x*t/Q2))*(F1+tau*F2);
  hc0_BH_LP = c00_LP*(c01_LP + c02_LP);
  
  double c10_LP = -8.*x*y*K*sqrt(1.+eps2)*(F1+F2)/(1.-tau);
  double c11_LP = (t/(2.*M_TARG2)-x*(1.-t/Q2))*(1.-x*(1.-t/Q2))*(F1+F2);
  double c12_LP = (1.+x-(3.-2.*x)*(1.+x*t/Q2)-x2*(1.+sqr(t/Q2))/tau)*(F1+tau*F2);
  hc1_BH_LP = c10_LP*(c11_LP + c12_LP);

  XSfact = alpha(3)*x*y/(16.*pi(2)*Q2*sqrt(1.+eps2));
  BHfact = 1./(x2*y2*sqr(1.+eps2)*t*P1*P2);

  double PHAS;
  if(ro->get_fVol() == 1)
    {
      // Differential cross section in the lab frame in nb.GeV-2.rd-2
      // d5sigma / (dx_B dy dt dPhi_e dPhi_g)
      PHAS  = 1.;
    }
  else if(ro->get_fVol() == 2)
    {
      // Differential cross section in the lab frame in nb.GeV-1.sr-2
      // d5sigma / (dk_e dOmega_e dOmega_g)
      double nu   = y*ro->get_fEb();
      double qmom = sqrt(Q2 + nu*nu);// Virtual photon momentum
      // Jacobian dx dy to scattered electron (GeV-1)
      PHAS = (ro->get_fEb()-nu)/(M_TARG*nu);
      // Jacobian momentum transfert to photon angle (GeV+2)
      PHAS *= sqr(Q2 + x*t)*qmom/(M_TARG*x*(1.-x)*Q2);
    }
  else
    {
      // Differential cross section in the lab frame in nb.GeV-4.rd-2
      // d5sigma / (dQ2 dx_B dt dPhi_e dPhi_g)
      // Jacobian dy to virtual photon mass (GeV-2)
      PHAS = y/Q2;
    }
  
  //xsec in nb
  hc0_BH *= 1.e7*hbarc(2)*BHfact*XSfact*PHAS;
  hc1_BH *= 1.e7*hbarc(2)*BHfact*XSfact*PHAS;
  hc2_BH *= 1.e7*hbarc(2)*BHfact*XSfact*PHAS;

  //long. pol. target
  hc0_BH_LP *= 1.e7*hbarc(2)*BHfact*XSfact*PHAS;
  hc1_BH_LP *= 1.e7*hbarc(2)*BHfact*XSfact*PHAS;

  return 0;
}


/*
==============================================================================
  Evaluate the Fourier coefficients of the DVCS term
==============================================================================
*/
int DVCS::dvcs_xsec(ReadOptFile *ro, double x, double y, double Q2, double t, double skew)
{
  int    Ipn = ro->get_fIpn();
  double M_TARG(m_targ(Ipn,1.));
  double M_TARG2(m_targ(Ipn,2.));
  double K, K2;
  double x2     = x*x;
  double y2     = y*y;
  double eps    = 2.*x*M_TARG/sqrt(Q2);
  double eps2   = eps*eps;
  double tau    = t/(4.*M_TARG2);
  double sskew  = -2.*skew/(1.+skew);
  double C_DVCS, C_DVCS_eff;
  double C_DVCS_LP, C_DVCS_eff_LP;
  double DVCSfact, XSfact;

  K2 = get_K2(Ipn, x, y, Q2, t);
  K  = sqrt(K2);
  
//  DVCS
  double HHconj            = sqr(Re_H)       + sqr(Im_H);    //HH*
  double EEconj            = sqr(Re_E)       + sqr(Im_E);    //EE*
  double HtHtconj          = sqr(Re_Ht)      + sqr(Im_Ht);   //H~H~*
  double EtEtconj          = sqr(Re_Ht)      + sqr(Im_Et);   //E~E~*
  double HEconjpEHconj     = 2.*(Re_H*Re_E   + Im_H*Im_E);   //HE* + EH*
  double HtEtconjpEtHtconj = 2.*(Re_Ht*Re_Et + Im_Ht*Im_Et); //H~E~* + E~H~*
  double HHtconjpHtHconj   = 2.*(Re_H*Re_Ht  + Im_H*Im_Ht);  //HH~* + H~H*  
  double EEtconjpEtEconj   = 2.*(Re_E*Re_Et  + Im_E*Im_Et);  //EE~* + E~E*  
  double HEtconjpEtHconj   = 2.*(Re_H*Re_Et  + Im_H*Im_Et);  //HE~* + E~H*  
  double HtEconjpEHtconj   = 2.*(Re_Ht*Re_E  + Im_Ht*Im_E);  //H~E* + EH~*  

  //eq. 66 (real part)
  C_DVCS        = (4.*(1.-x)*(HHconj+HtHtconj) - x2*tau*EtEtconj - (x2+sqr(2.-x)*tau)*EEconj 
                   - x2*(HEconjpEHconj+HtEtconjpEtHtconj))/sqr(2.-x);

  //C_DVCS_eff  = -x*C_DVCS; // approximation of eq. 52
  C_DVCS_eff    = sskew*C_DVCS; // approximation of eq. 52

  hc0_DVCS      =  2.*(1.+sqr(1.-y))*C_DVCS;
  hc1_DVCS      =  8.*K*(2.-y)/(2.-x)*C_DVCS_eff;
  hs1_DVCS      = -8.*K*y/(2.-x)*C_DVCS_eff;

  //long. pol. target
  C_DVCS_LP     = (4.*(1.-x)*HHtconjpHtHconj - x*(x2/2. + (2.-x)*tau)*EEtconjpEtEconj 
                   - x2*(HEtconjpEtHconj+HtEconjpEHtconj))/sqr(2.-x);

  C_DVCS_eff_LP = sskew*C_DVCS_LP;

  hc0_DVCS_LP   = 2.*y*(2. - y)*C_DVCS_LP;
  hc1_DVCS_LP   = 8.*K*y/(2. - x)*C_DVCS_eff_LP;
  hs1_DVCS_LP   = -8.*K*(2.-y)/(2. - x)*C_DVCS_eff_LP;

  XSfact     = alpha(3)*x*y/(16.*pi(2)*Q2*sqrt(1. + eps2));
  DVCSfact   = 1./(y2*Q2);

  double PHAS;
  if(ro->get_fVol() == 1)
    {
      //Differential cross section in the lab frame in nb.GeV-2.rd-2
      //d5sigma / (dx_B dy dt dPhi_e dPhi_g)
      PHAS  = 1.;
    }
  else if(ro->get_fVol() == 2)
    {
      // Differential cross section in the lab frame in nb.GeV-1.sr-2
      // d5sigma / (dk_e dOmega_e dOmega_g)
      double nu   = y*ro->get_fEb();
      double qmom = sqrt(Q2 + nu*nu);// Virtual photon momentum
      // Jacobian dx dy to scattered electron (GeV-1)
      PHAS = (ro->get_fEb() - nu)/(M_TARG*nu);
      // Jacobian momentum transfert to photon angle (GeV+2)
      PHAS *= sqr(Q2 + x*t)*qmom/(M_TARG*x*(1. - x)*Q2);
    }
  else
    {
      // Differential cross section in the lab frame in nb.GeV-4.rd-2
      // d5sigma / (dQ2 dx_B dt dPhi_e dPhi_g)
      // Jacobian dy to virtual photon mass (GeV-2)
      PHAS = y/Q2;
    }
  
  hc0_DVCS *= 1.e7*hbarc(2)*DVCSfact*XSfact*PHAS;
  hc1_DVCS *= 1.e7*hbarc(2)*DVCSfact*XSfact*PHAS;
  hs1_DVCS *= 1.e7*hbarc(2)*DVCSfact*XSfact*PHAS;

  //long. pol. target
  hc0_DVCS_LP *= 1.e7*hbarc(2)*DVCSfact*XSfact*PHAS;
  hc1_DVCS_LP *= 1.e7*hbarc(2)*DVCSfact*XSfact*PHAS;
  hs1_DVCS_LP *= 1.e7*hbarc(2)*DVCSfact*XSfact*PHAS;

  return 0;
}


/*
==============================================================================
  Evaluate the Fourier coefficients of the Interference term
==============================================================================
*/
int DVCS::int_xsec(ReadOptFile *ro, double *FF, double x, double y, double Q2, double t, double skew, double Phi_b)
{
  int    Ipn = ro->get_fIpn();
  double M_TARG(m_targ(Ipn,1.));
  double M_TARG2(m_targ(Ipn,2.));
  double Phi_bs = pi(1) - Phi_b;
  double K, K2, J, P1, P2;
  double y3     = y*y*y;
  double eps    = 2.*x*M_TARG/sqrt(Q2);
  double eps2   = eps*eps;
  double tau    = t/(4.*M_TARG2);
  double sskew  = -2.*skew/(1.+skew);
  double Re_C_I, Re_DelC_I, Re_C_I_eff;
  double Im_C_I, Im_DelC_I, Im_C_I_eff;
  double Re_C_I_LP, Re_DelC_I_LP, Re_C_I_eff_LP;
  double Im_C_I_LP, Im_DelC_I_LP, Im_C_I_eff_LP;
  double INTERfact, XSfact;

  double F1 = FF[0];
  double F2 = FF[1];

  K2 = get_K2(Ipn, x, y, Q2, t);
  J  = get_J(Ipn, x, y, Q2, t);
  K  = sqrt(K2);
  P1 = -(J + 2.*K*cos(Phi_bs))/(y*(1. + eps2));
  P2 = 1. + t/Q2 - P1;

  // INTERF
  Re_C_I     = F1*Re_H + x*(F1+F2)*Re_Ht/(2.-x) - tau*F2*Re_E;
  Im_C_I     = F1*Im_H + x*(F1+F2)*Im_Ht/(2.-x) - tau*F2*Im_E;

  Re_DelC_I  = -x*(F1+F2)*(x*(Re_H+Re_E)/(2.-x) + Re_Ht)/(2.-x);
  Im_DelC_I  = -x*(F1+F2)*(x*(Im_H+Im_E)/(2.-x) + Im_Ht)/(2.-x);

  Re_C_I_eff = sskew*Re_C_I;
  Im_C_I_eff = sskew*Im_C_I;

  hc0_Int    = -8.*(2.-y)*(sqr(2.-y)*K2*Re_C_I/(1.-y) + (t/Q2)*(1.-y)*(2.-x)*(Re_C_I+Re_DelC_I));
  hc1_Int    = -8.*K*(1.+sqr(1.-y))*Re_C_I;
  hc2_Int    = -16.*K2*(2.-y)*Re_C_I_eff/(2.-x);
  hs1_Int    =  8.*K*y*(2.-y)*Im_C_I;
  hs2_Int    =  16.*K2*y*Im_C_I_eff/(2.-x);

  //long. pol. target
  Re_C_I_LP     = x*(F1+F2)*(Re_H+(x/2.)*Re_E)/(2.-x) + F1*Re_Ht - x*((x/2.)*F1+tau*F2)*Re_Et/(2.-x);
  Im_C_I_LP     = x*(F1+F2)*(Im_H+(x/2.)*Im_E)/(2.-x) + F1*Im_Ht - x*((x/2.)*F1+tau*F2)*Im_Et/(2.-x);

  Re_DelC_I_LP  = -x*(F1+F2)*(Re_H + (x/2.)*Re_E + x*(Re_Ht+(x/2.)*Re_Et)/(2.-x))/(2.-x);
  Im_DelC_I_LP  = -x*(F1+F2)*(Im_H + (x/2.)*Im_E + x*(Im_Ht+(x/2.)*Im_Et)/(2.-x))/(2.-x);

  Re_C_I_eff_LP = sskew*Re_C_I_LP;
  Im_C_I_eff_LP = sskew*Im_C_I_LP;

  hc0_Int_LP    = -8.*y*((sqr(2.-y)/(1.-y)+2.)*K2*Re_C_I_LP + (t/Q2)*(1.-y)*(2.-x)*(Re_C_I_LP+Re_DelC_I_LP));
  hc1_Int_LP    = -8.*K*y*(2.-y)*Re_C_I_LP;
  hc2_Int_LP    = -16.*K2*y*Re_C_I_eff_LP/(2.-x);
  hs1_Int_LP    =  8.*K*(1.+sqr(1.-y))*Im_C_I_LP;
  hs2_Int_LP    =  16.*K2*(2.-y)*Im_C_I_eff_LP/(2.-x);

  XSfact     = alpha(3)*x*y/(16.*pi(2)*Q2*sqrt(1. + eps2));
  INTERfact  = 1./(x*y3*t*P1*P2);

  double PHAS;
  if(ro->get_fVol() == 1)
    {
      // Differential cross section in the lab frame in nb.GeV-2.rd-2
      // d5sigma / (dx_B dy dt dPhi_e dPhi_g)
      PHAS  = 1.;
    }
  else if(ro->get_fVol() == 2)
    {
      // Differential cross section in the lab frame in nb.GeV-1.sr-2
      // d5sigma / (dk_e dOmega_e dOmega_g)
      double nu   = y*ro->get_fEb();
      double qmom = sqrt(Q2 + nu*nu);// Virtual photon momentum
      // Jacobian dx dy to scattered electron (GeV-1)
      PHAS = (ro->get_fEb() - nu)/(M_TARG*nu);
      // Jacobian momentum transfert to photon angle (GeV+2)
      PHAS *= sqr(Q2 + x*t)*qmom/(M_TARG*x*(1. - x)*Q2);
    }
  else
    {
      // Differential cross section in the lab frame in nb.GeV-4.rd-2
      // d5sigma / (dQ2 dx_B dt dPhi_e dPhi_g)
      // Jacobian dy to virtual photon mass (GeV-2)
      PHAS = y/Q2;
    }
  
  hc0_Int  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hc1_Int  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hc2_Int  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hs1_Int  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hs2_Int  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;

  //long. pol. target
  hc0_Int_LP  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hc1_Int_LP  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hc2_Int_LP  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hs1_Int_LP  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;
  hs2_Int_LP  *= 1.e7*hbarc(2)*INTERfact*XSfact*PHAS;

  return 0;
}


int DVCS::phot_xsec(ReadOptFile *ro, NuclFF *ff, double x, double y, double Q2, double t, double Phi_g)
{
  //Nucleon form factors
  FF1[0] = FF1[1] = &NuclFF::dirac_ff;
  FF2[0] = FF2[1] = &NuclFF::pauli_ff;
  NuclFF* instance = new NuclFF;
  for(int i=0; i<2; i++) ff->set_FF1(i, (instance->*FF1[i])(-t,i));
  for(int i=0; i<2; i++) ff->set_FF2(i, (instance->*FF2[i])(-t,i));
  delete instance;

  double FF[2];
  FF[0] = ff->get_FF1(ro->get_fIpn());
  FF[1] = ff->get_FF2(ro->get_fIpn());

  double skew_old = x/(2.-x);
  //double skew_new = x*(1.+t/(2.*Q2))/(2.-x*(1.- t/Q2));

  double skew=skew_old;
  if(skew <= 0.01 || skew >= 1.0)
    {
      //cout<<"skewdness out of limit "<<skew<<" "<<skew_old<<endl;
      return -1;
    }
  
  eval_cffs(ro, ff, skew, t);
  
  bh_xsec(ro, FF, x, y, Q2, t, Phi_g);
  dvcs_xsec(ro, x, y, Q2, t, skew);
  int_xsec(ro, FF, x, y, Q2, t, skew, Phi_g);

  return 0;
}
