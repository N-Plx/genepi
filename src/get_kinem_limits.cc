#include "cpp.h"
#include "read_optFile.h"
#include "genepi.h"
#include "inl_funcs.h"

void get_kinem_limits(ReadOptFile *ro)
{
  int    Ipn = ro->get_fIpn();
  double M_TARG(m_targ(Ipn,1.));
  double M_TARG2(m_targ(Ipn,2.));
  double  s = 2.*ro->get_fEb()*M_TARG;//s = parl21 = 2kP com Energy
  double fXbjmin(ro->get_fXbjmin()), fXbjmax(ro->get_fXbjmax());
  double fYmin(ro->get_fYmin()), fYmax(ro->get_fYmax());
  double fQ2min(ro->get_fQ2min()), fQ2max(ro->get_fQ2max());
  double fNumin(ro->get_fNumin()), fNumax(ro->get_fNumax());
  double fW2min(ro->get_fW2min()), fW2max(ro->get_fW2max());
  /*
  for(int j=0;j<2;j++)
  {
    fXbjmin = max(max(fXbjmin,fQ2min/(s*fYmax)), 
              max(fQ2min/(2.*M_TARG*fNumax), 1.-(fW2max-M_TARG2)/max(s*fYmin,1.e-22)));
    fXbjmin = max(fXbjmin, 
                  1.-(fW2max-M_TARG2)/max(2.*M_TARG*fNumin,1.E-22));

    fXbjmax = min(min(fXbjmax,fQ2max/max(s*fYmin,1.E-22)),
                  min(fQ2max/max(2.*M_TARG*fNumin,1.E-22),1.-(fW2min-M_TARG2)/(s*fYmax)));
    fXbjmax = min(fXbjmax,
                  1.-(fW2min-M_TARG2)/(2.*M_TARG*fNumax));

    fYmin = max(max(fYmin,fQ2min/(s*fXbjmax)),
                max((fW2min-M_TARG2)/(s*(1.-fXbjmin)),(fW2min-M_TARG2+fQ2min)/s));
    fYmin = max(fYmin,
                2.*M_TARG*fNumin/s);

    fYmax = min(min(fYmax,fQ2max/max(s*fXbjmin,1.E-22)),
                min((fW2max-M_TARG2)/max(s*(1.-fXbjmax),1.E-22),(fW2max-M_TARG2+fQ2max)/s));
    fYmax = min(fYmax,
                2.*M_TARG*fNumax/s);

    fQ2min = max(max(fQ2min,s*fXbjmin*fYmin),
                 max(s*fYmin-fW2max+M_TARG2,2.*M_TARG*fNumin*fXbjmin));
    fQ2min = max(fQ2min,
                 fXbjmin*(fW2min-M_TARG2)/(1.-fXbjmin));

    fQ2max = min(min(fQ2max,s*fXbjmax*fYmax),
                 min(s*fYmax-fW2min+M_TARG2,2.*M_TARG*fNumax*fXbjmax));
    fQ2max = min(fQ2max,
                 fXbjmax*(fW2max-M_TARG2)/max(1.-fXbjmax,1.E-22));

    fW2min = max(max(fW2min,s*(1.-fXbjmax)*fYmin+M_TARG2),
                 max(fQ2min*(1.-fXbjmax)/fXbjmax+M_TARG2,s*fYmin-fQ2max+M_TARG2));
    fW2min = max(fW2min,
                 2.*M_TARG*fNumin*(1.-fXbjmax)+M_TARG2);

    fW2max = min(min(fW2max,s*(1.-fXbjmin)*fYmax+M_TARG2),
                 min(fQ2min*(1.-fXbjmin)/max(fXbjmin,1.E-22)+M_TARG2,s*fYmax-fQ2min+M_TARG2));
    fW2max = min(fW2max,
                 2.*M_TARG*fNumax*(1.-fXbjmin)+M_TARG2);

    fNumin = max(fNumin,
                 (fW2min+fQ2min-M_TARG2)/(2.*M_TARG));
    
    fNumax = min(fNumax,
                 min(ro->get_fEb(),(fW2max+fQ2max-M_TARG2)/(2.*M_TARG)));
  }
  */
  ro->set_fXbjmin(fXbjmin); ro->set_fXbjmax(fXbjmax);
  ro->set_fYmin(fYmin);     ro->set_fYmax(fYmax);
  ro->set_fQ2min(fQ2min);   ro->set_fQ2max(fQ2max);
  ro->set_fNumin(fNumin);   ro->set_fNumax(fNumax);
  ro->set_fW2min(fW2min);   ro->set_fW2max(fW2max);
}


double get_tmin(int Ipn, double x, double Q2)
{
  double eps    = 2.*x*m_targ(Ipn,1.)/sqrt(Q2);
  double eps2   = sqr(eps);
  double tmin   = -Q2*(2.*(1.-x)*(1.-sqrt(1.+eps2)) + eps2)/(4.*x*(1.-x) + eps2);

  return tmin;
}


double get_K2(int Ipn, double x, double y, double Q2, double t)
{
  double eps     = 2.*x*m_targ(Ipn,1.)/sqrt(Q2);
  double eps2    = sqr(eps);
  double tmin    = get_tmin(Ipn, x, Q2);
  double ttminQ2 = (t - tmin)/Q2;
  double K2      = -ttminQ2*(1. - x)*(1. - y - y*y*eps2/4.);
  K2            *= (sqrt(1. + eps2)+(4.*x*(1. - x)+eps2)*ttminQ2/(4.*(1. - x)));
  K2             = max(0., K2);

  return K2;
}


double get_J(int Ipn, double x, double y, double Q2, double t)
{
  double eps  = 2.*x*m_targ(Ipn,1.)/sqrt(Q2);
  double eps2 = sqr(eps);
  double J    = (1.-y - y*eps2/2.)*(1.+t/Q2) - (1.-x)*(2.-y)*t/Q2;

  return J;
}


double get_ycol(double t, double x, double Q2)
{
  double ycol = (Q2+t)/(Q2+x*t);

  return ycol;
}
