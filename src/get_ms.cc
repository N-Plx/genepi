#include "genepi.h"
#include "vect.h"
#include "inl_funcs.h"
#include "jetset.h"
#include "read_optFile.h"

// this recalculates the produced particle energy and momentum assuming a stationary target but a non-zero produced particle mass, writes it to lujets.

int get_ms(ReadOptFile *ro, double x, double Q2, double nu, double Ems, double t, double Phikkp, double Phiqms)
{
  VECT   *vect = new VECT();
  double k;
  double Ep, kp, Thetakkp, cosThetakkp;
  double q, Thetakq, Phikq;
  double Pms, Pms2;
  double Thetaqms, cosThetaqms;
  double Thetakms, Phikms;  //wrt inc. lepton frame
  static int  nlist(0);     //counter to print lulist()
  
  int    MS_ID  = ms_id(ro->get_fIms());
  double M_MS2  = m_ms(ro->get_fIms(),2.);
  
  q  = sqrt(nu*nu + Q2);
  Ep = ro->get_fEb() - nu;   // scattered electron energy. Unfortunate naming.
  
  // beam electron 3-momentum:
  k = sqrt(sqr(ro->get_fEb()) - m_elec(2)); 
  vector<double> vk = vect->fill_vect(k,0.,0.);
  
  // angle between scattered electron and beam electron. Doesn't neglect electron mass:
  cosThetakkp = (k*k - ro->get_fEb()*nu - Q2/2.) / (k*sqrt(k*k + nu*nu - 2.*ro->get_fEb()*nu));
  
  if(abs(cosThetakkp) >= 1.)
    {
      //cout<<"get_ms() : cosThetakkp larger than one; cosThetakkp = "<<cosThetakkp<< endl;
      cosThetakkp = sign(cosThetakkp);
    }
  
  Thetakkp =  acos(cosThetakkp);//[0,PI]
  if(Thetakkp > pi(1)) 
    cout<<"Thetakkp larger than pi; cosThetakkp = "<<cosThetakkp<< endl;
  
  
  kp = sqrt(Ep*Ep - m_elec(2));
  vector<double> vkp = vect->fill_vect(kp, Thetakkp, Phikkp); // 3mom of the scattered electron
  
  vector<double> vq =  vect->diff_vect(vk, vkp);
  
  Thetakq = get_theta(vq.at(2), q);
  Phikq   = get_phi(vq.at(0), vq.at(1));
  
  Pms2 = sqr(Ems) - M_MS2;   // squared 3mom of produced meson
  if(Pms2 <= 0) 
    {
      //cout<<"Negative Pms2 "<<Pms2<<endl;
      return -1;
    }
  Pms = sqrt(Pms2);
  
  cosThetaqms = (t + Q2 + 2.*Ems*nu - M_MS2)/(2*Pms*q);
  if(abs(cosThetaqms) >= 1.)
    {
      //cout<<"get_ms() : cosThetaqms larger than one; cosThetaqms = "<<cosThetaqms<< endl;
      return -1;
    }
  
  Thetaqms =  acos(cosThetaqms);//[0,PI]
  if(Thetaqms > pi(1)) 
    cout<<"Thetaqms larger than pi; Thetaqms = "<<cosThetaqms<< endl;
  
  //meson momentum in virt. phot. frame
  vector<double> vPms1 = vect->fill_vect(Pms, Thetaqms, Phiqms);
  
  //meson momentum in inc. elec. frame
  vector<double> vPms = vect->invrot_vect(vPms1, Thetakq, Phikq);
  
  Thetakms = get_theta(vPms.at(2), Pms);
  Phikms   = get_phi(vPms.at(0), vPms.at(1));

  int   i1(1);
  int   i2(2);

  float fl_Ems      = (float) Ems;
  float fl_Thetakms = (float) Thetakms;
  float fl_Phikms   = (float) Phikms;

  lu1ent_(i1, MS_ID, fl_Ems, fl_Thetakms, fl_Phikms);

  luexec_();

  if(nlist<5)
  {
    lulist_(i2);
    nlist++;
  }

  return 0;
}
