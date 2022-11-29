#include "genepi.h"
#include "dvcs_vars.h"
#include "inl_funcs.h"

double ms_xsec(ReadOptFile *ro, double x, double Q2, double t, double Phi_g, int Ims, double tmin)
{

  // Ims: specifies which meson

  int Ipn(ro->get_fIpn());
  double M_TARG2(m_targ(Ipn,2.));
  double dd;
  double x2 = x*x;
  double x3 = x*x*x;
  double xsec, dmsunp, dmspol;
  
  double delu = ups(x) - ums(x);
  double deld = dps(x) - dms(x);
  //  double tmin = get_tmin(Ipn, x, Q2);   // assumes zero meson mass! Just lovely for a function specifically aimed at mesons.
  
  double M_mes2 = m_ms(Ims,2);
 
  if(ro->get_fIms() == 0) //pi0
  {
    if(Ipn == 0) //neutron
    {
      dd = sqr(2.*deld + delu);
    }
    else if(Ipn == 1) //proton
    {
      dd = sqr(2.*delu + deld);
    }
  }
  if(ro->get_fIms() == 1) //eta
  {
    if(Ipn == 0) //neutron
    {
      dd = sqr(2.*deld - delu)/3.;
    }
    else if(Ipn == 1) //proton
    {
      dd = sqr(2.*delu - deld)/3.;
    }
  }

  dmsunp = (25.2*dd*x3*(1. - x))/(sqr(Q2*(Q2 + M_TARG2)))*
            (1. + 2.*ro->get_fBheli()*x*pow((1. - x),5)*sin(Phi_g))*exp((t - tmin));

  dmspol = 6*x2*pow((1. - x),5);

  xsec  = xsecpi0()*dmsunp*(1. + ro->get_fTheli()*dmspol);

  return xsec;
}
