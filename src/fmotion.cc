#include "genepi.h"
#include "dvcs_vars.h"
#include "inl_funcs.h"

/*
input:
  rndm_t random variable to generate the squared momentum transfer t
  rndm_phi azimutal angle of the produced particle wrt virtual photon in CM
  q (nu) Energy (momentum) of the virtual photon
  Ep (Pp) Energy (momentum) of the incoming nucleon
  thetap (phip) polar (azimuthal) angle of the incoming nucleon wrt inc. electron (LAB)
  thetag (phig) polar (azimuthal) angle of the virtual photon wrt inc. electron (LAB)
  ipn specifies the type of target nucleon (p or n)
  ims specifies the produced meson (pi0 or eta). If it's a DVCS process,
      the value -1 is passed to the fmotion function in the place of ims
output:
  t  the squared momentum transfer
  nup Energy of the outgoing produced particle (photon or meson)
  Epp Energy of the outgoing nucleon
  thetaqp (Phiqp) polar (azimuthal) angle of the outgoing produced particle wrt inc. electron (LAB)
  thetapp (Phipp) polar (azimuthal) angle of the outgoing nucleon wrt inc. electron (LAB)
*/
int fmotion(double rndm_t,double rndm_phi,int ipn,int ims,double q,double nu,double W2, double Pp,
            double thetap,double phip,double thetag,double phig,double * out)
{
  double M_TARG2(m_targ(ipn,2.));
  double Ep = sqrt(M_TARG2 + sqr(Pp));
  double t;
  double W, Pcm, Pcm2;//\vec{Pcm} = \vec{q} + \vec{P}
  double nup, qp, thetaqp, Phiqp;  //real photon
  double Epp, Ppp, thetapp, Phipp; //scat nucleon
  double Q2    = q*q - nu*nu;
  double cthp  = cos(thetap);
  double sthp  = sin(thetap);
  double cphip = cos(phip);
  double sphip = sin(phip);
  double cthg  = cos(thetag);
  double sthg  = sin(thetag);
  double cphig = cos(phig);
  double sphig = sin(phig);
  double cthpq = sthg*cphig*sthp*cphip + sthg*sphig*sthp*sphip + cthg*cthp;
  //save nucleon's variables
  double Ppf    = Pp;
  double Epf    = Ep;

  double qp0, qp02, nup0; //photon momentum and Energy in CM frame
  double btcm, gacm;
  double cthcm, sthcm, cphicm, sphicm;
  double cthgcm, sthgcm, cphigcm, sphigcm;
  double ql, qt, nu0, q0l, q0t, q0;
  double cthg0, sthg0, cphig0, sphig0;
  double cthpfcm, sthpfcm, Ppfl, Ppft, Epf0, Ppf0l, Ppf0t, Ppf0;
  double tmin, tmax, xabs;
  double cthqqp0, sthqqp0, phiqqp0, cphiqqp0, sphiqqp0;
  double cthqp0, sthqp0, cphiqp0, sphiqp0;
  double qpl, qpt;
  double cthqpcm, sthqpcm, cphiqpcm, sphiqpcm;
  double cthqp, sthqp, cphiqp, sphiqp, cthqqp;
  double Ep0, Ppl, Ppt;
  double cthpcm, sthpcm, cphipcm, sphipcm;

  W = sqrt(W2);
  Pcm2 = Pp*Pp + q*q + 2.*Pp*q*cthpq;

  if(Pcm2 < 0.)
    {
      cout<<"negative Pcm2 "<<Pcm2<<endl;
      return -1;
    }
  Pcm  = sqrt(Pcm2);
  
  // ims 0 = pi0, 1 = eta. -1 is the value passed in place of ims if it's a real photon
  if (ims == -1) qp0  = (W2 - M_TARG2)/(2.*W);
  else qp0  = sqrt(0.25 * sqr(W2 - M_TARG2 - m_ms(ims,2)) - M_TARG2 * m_ms(ims,2))/W;
  if(qp0 < 0.)
    {
      cout<<"negative qp0 "<<qp0<<"; W2 = "<<W2<<endl;
      return -1;
    }
  qp02 = sqr(qp0);
  if (ims == -1) nup0 = qp0;
  else nup0 = sqrt(qp02 + m_ms(ims,2));
  
  // beta and gamma of the CM
  btcm = Pcm/(nu + Ep);
  gacm = 1./sqrt(1. - btcm*btcm);//(nu + Ep)/W;
  
  // Calculate the CM polar angle in the incident electron frame (Lab)
  cthcm = (q*cthg + Pp*cthp)/Pcm;
  if(abs(cthcm) > 1.0) cthcm = cthcm/abs(cthcm);  // presumably to fix issues from rounding errors, if it goes a tiny bit over 1.
  sthcm = sqrt(1. - cthcm*cthcm);
  
  // Calculate the CM azimutal angle in the incident electron frame (Lab)
  if(sthcm != 0.0)
    {
      cphicm = (q*sthg*cphig + Pp*sthp*cphip)/(Pcm*sthcm);
      sphicm = (q*sthg*sphig + Pp*sthp*sphip)/(Pcm*sthcm);
    }
  else
    {  // if theta_cm is 0 or 180deg, Pcm is collinear with incident electron. In that frame phi doesn't matter, so can set it to anything.
      cout<<"sthcm == 0.0"<<endl;
      cphicm = cphig;
      sphicm = sphig;
    }
  
  if(abs(cphicm) > 1)
    {
      cphicm = cphicm/abs(cphicm);
      sphicm = 0.0;
    }
  else if(abs(sphicm) > 1)
    {
      sphicm = sphicm/abs(sphicm);
      cphicm = 0.0;
    }
  // Now you have the full three-vector Pcm, the momentum of the CM in the incident electron Lab frame.
    
  // Transform (nu,q) to the CM (nu0, q0)
  // Calculate photon angles with CM in the Lab frame
  cthgcm = (q + Ppf*cthpq)/Pcm;   // q + Ppf*cthpq gives total 3mom along direction of q, then apply trig. to get the angle this makes to the Pcm vector.
  if(abs(cthgcm) > 1.) cthgcm = cthgcm/abs(cthgcm);
  sthgcm = sqrt(1. - cthgcm*cthgcm);
  
  // For the phi angles, obtain them from calculating the rotation matrix:
  // First, find matrix for rotation by -phi_cm around lab z-axis, to align lab x-axis with Pcm projection into xy plane.
  // Next, find matrix for rotation by -theta_cm around lab y-axis, to align lab z-axis with Pcm.
  // Matrix resulting from multiplication of these two will rotate q into frame where z-axis is aligned with Pcm.
  if(sthgcm != 0.0)
    {
      cphigcm = (cthcm*cphicm*sthg*cphig + cthcm*sphicm*sthg*sphig - sthcm*cthg)/sthgcm; // from the x-component of the rotated q-vector
      sphigcm = (-sphicm*sthg*cphig + cphicm*sthg*sphig)/sthgcm;  // from the y-component of the rotated q-vector
    }
  else
    {
      cphigcm = 1.0;     // either aligned or anti-alligned with Pcm, so phi_gcm is arbitrary
      sphigcm = 0.0;
    }
  
  ql  = q*cthgcm;   // longitudinal component of rotated q-vector (for the application of the boost)
  qt  = q*sthgcm;   // transverse component of the rotated q-vector
  nu0 = gacm*(nu - btcm*ql);
  if(nu0<0)
    {
      //cout<<"neg. nu0 : "<<nu0<<" "<<gacm<<"  "<<nu<<"  "<<btcm<<"  "<<
      //ql<<"  "<<q<<"  "<<cthgcm<<"  "<<sthgcm<<" "<<Q2<<" "<<W2<<endl;
      return -1;
    }
  q0l = gacm*(ql - btcm*nu);        // q0 || to P_cm direction
  q0t = qt;                         // q0 perp. to P_cm direction
  q0  = sqrt(q0l*q0l + q0t*q0t);    // q0
    
  // Now you have the virtual photon four-momentum in the CM frame

  // Calculate virtual photon angles with the CM frame
  cthg0 = q0l/q0;
  sthg0 = sqrt(1. - cthg0*cthg0);
  if(sthg0 != 0.0)  // transverse components unaffected by boost
    {
      cphig0 = cphigcm;
      sphig0 = sphigcm;
    }
  else
    {
      cphig0 = 1.0;   // as usual, if momentum is aligned with z in the CM frame, set phi arbitrarily to 0.
      sphig0 = 0.0;
    }
  
  // Transform (Epf,Ppf) to CM (Epf0,Ppf0)
  // Polar angle of the incident nucleon with CM direction
  cthpfcm = (Ppf + q*cthpq)/Pcm;         // q*cthpq = projection of q along Ppf. Then apply trig. to right-angled triangle of (Ppf + q*cthpq) and Pcm
  sthpfcm = sqrt(1. - cthpfcm*cthpfcm);
  Ppfl    = Ppf*cthpfcm;
  Ppft    = Ppf*sthpfcm;
  Epf0  = gacm*(Epf  - btcm*Ppfl);
  Ppf0l = gacm*(Ppfl - btcm*Ep);
  Ppf0t = Ppft;
  Ppf0  = sqrt(Ppf0l*Ppf0l + Ppf0t*Ppf0t);
 
  // Now have the incident nucleuon energy, Epf0, and momentum, Ppf0, in the CM frame. The angles are opposite to the q0 vector.

  // Calculate the values of tmin and tmax using resp. thqqp0 =0 and PI
  // ims 0 = pi0, 1 = eta. -1 is the value passed in place of ims if it's a real photon
  if (ims == -1)
    {
      tmax = -Q2 - 2.*nu0*nup0 + 2.*q0*qp0;
      tmin = -Q2 - 2.*nu0*nup0 - 2.*q0*qp0;
    }
  else // produced particle is a meson
    {
      tmax = -Q2 + m_ms(ims,2) - 2.*nu0*nup0 + 2.*q0*qp0;
      tmin = -Q2 + m_ms(ims,2) - 2.*nu0*nup0 - 2.*q0*qp0;
    }
        
  // Generate the Momentum transfer t (equivalent to generating a random phi angle for the produced particle)
  t = tmin + (tmax - tmin)*rndm_t;

  // Calculate the produced particle polar angle with virtual photon in the CM
  if (ims == -1) cthqqp0 = (t + Q2 + 2.*nup0*nu0)/(2.*q0*qp0);
  else cthqqp0 = (t + Q2 - m_ms(ims,2) + 2.*nup0*nu0)/(2.*q0*qp0);
  xabs = abs(cthqqp0);
  if(xabs > 1.0)
    {
      cout<<"bad generated t "<<endl;
      return -1;
    }
  sthqqp0 = sqrt(1.0 - cthqqp0*cthqqp0);

  // Generate the azimutal angle of the produced particle with the photon in CM : 0->2PI
  phiqqp0 = 2.*pi(1)*rndm_phi; // CM is not || to photon
  cphiqqp0 = cos(phiqqp0);
  sphiqqp0 = sin(phiqqp0);

  // Transform from frame aligned with the virtual photon to CM frame : simple rotation
  // You have a vector w.r.t. q0. If you rotate it by whatever you need to rotate q0 by to get it aligned with CM axes, you'll have qp0 w.r.t. CM axes.
  // To rotate q0 to align it with CM axes, first rotate by theta_g0 around CM y, then by phi_g0 around CM z: matrix is multiplication of the two.
  // Polar angle
  cthqp0 = -sthg0*sthqqp0*cphiqqp0 + cthg0*cthqqp0;
  sthqp0 = sqrt(1. - cthqp0*cthqp0);
  // Azimuthal angle
  cphiqp0 = (cthg0*cphig0*sthqqp0*cphiqqp0 - sphig0*sthqqp0*sphiqqp0 + sthg0*cphig0*cthqqp0)/sthqp0;
  sphiqp0 = (cthg0*sphig0*sthqqp0*cphiqqp0 + cphig0*sthqqp0*sphiqqp0 + sthg0*sphig0*cthqqp0)/sthqp0;

  // Transform from CM to Lab system : Lorentz transform
  qpl = gacm *(btcm*nup0 + qp0*cthqp0);
  qpt = qp0 *sthqp0;
  qp  = sqrt(qpl*qpl + qpt*qpt);
  if (ims == -1) nup = qp;
  else nup = gacm *(nup0 + btcm*qp0*cthqp0);

  // Polar angle of the produced particle with the CM momentum direction in the lab system
  cthqpcm = qpl/qp;
  sthqpcm = qpt/qp;

  // Azimutal angle of the produced particle with the CM direction: 0->2PI
  cphiqpcm = cphiqp0;  // it's the same because it is transverse to the boost
  sphiqpcm = sphiqp0;

  // Angles of the produced particle with the incident electron (z-axis) in the lab (same principle of rotation as for the qp0 vector w.r.t. CM frame above):
  // Polar angle of produced particle (qp) in the lab:
  cthqp = -sthqpcm*cphiqpcm*sthcm + cthqpcm*cthcm;
  xabs = abs(cthqp);
  if(xabs > 1.0) cthqp = cthqp/xabs;
  sthqp = sqrt(1.0 -cthqp*cthqp);
  thetaqp = acos(cthqp);

  // Calculate the azimutal angle of the produced particle with the lab axes:
  cphiqp = (sthqpcm*cphiqpcm*cthcm*cphicm - sthqpcm*sphiqpcm*sphicm + cthqpcm*sthcm*cphicm)/sthqp;
  sphiqp = (sthqpcm*cphiqpcm*cthcm*sphicm + sthqpcm*sphiqpcm*cphicm + cthqpcm*sthcm*sphicm)/sthqp;

  xabs = abs(cphiqp);
  if(xabs > 1.0)
    {
      cphiqp = cphiqp/xabs;
      sphiqp = 0.0;
    }
  Phiqp = acos(cphiqp);  // between 0 and PI --> next check the sin
  if(sphiqp < 0.0) Phiqp = 2.*pi(1) - Phiqp; // to be between 0 and 2PI
    
  // Now have the produced particle energy, momentum and polar and azimuthal angles in the Lab.

  // produced particle and virtual photon relative angle to check t
  cthqqp = (sthqp*cphiqp)*(sthg*cphig) + (sthqp*sphiqp)*(sthg*sphig) + cthqp*cthg;  // via dot product of qp and q
  double ttmp = 0.;
  if (ims == -1) ttmp = -Q2 -2.*nu*nup + 2.*q*qp*cthqqp;
  else ttmp = -Q2 + m_ms(ims,2) - 2.*nu*nup + 2.*q*qp*cthqqp;
  if(abs(ttmp-t)>tiny())
    {
      cout<<"t check "<<endl;
      cout<<nu <<"  "<<nup <<"  "<<qp <<"  "<<cthqqp <<"  "<<ttmp <<endl;
      cout<<nu0<<"  "<<nup0<<"  "<<qp0<<"  "<<cthqqp0<<"  "<<t <<endl;
      cout<<sthqp<<" "<<cphiqp<<" "<<sthg<<" "<<cphig<<" "<<sthqp<<" "<<
	sphiqp<<" "<<sthg<<" "<<sphig<<" "<<cthqp<<" "<<cthg<<endl<<endl;
    }

  // Calculate the momentum and angle of the recoil nucleon
  Ep0 = sqrt(M_TARG2 + qp02);
  Ppl = gacm *(btcm*Ep0 - qp0*cthqp0);   // momentum of recoil is equal and opposite to momentum of produced particle in CM
  Ppt = qp0*sthqp0;    // no - sign add PI to phi
  Ppp = sqrt(Ppl*Ppl + Ppt*Ppt);   // in Lab
  Epp = sqrt(M_TARG2 + Ppp*Ppp);   // in Lab

  // Polar angle of the recoil nucleon with the CM direction in the lab
  cthpcm = Ppl/Ppp;
  sthpcm = Ppt/Ppp;

  // Azimuthal angle of the recoil nucleon with the CM in the lab
  cphipcm = -cphiqpcm;
  sphipcm = -sphiqpcm;

  // Polar angle of the recoil nucleon with the incident electron : beam axis (same rotation principle as for qp above)
  cthp    = -sthpcm*cphipcm*sthcm + cthpcm*cthcm;
  sthp    = sqrt(1.0 - cthp*cthp);
  thetapp = acos(cthp);

  // Azimutal angle of the nucleon with the incident electron : beam axis (same rotation principle as for qp above)
  cphip = (sthpcm*cphipcm*cthcm*cphicm - sthpcm*sphipcm*sphicm + cthpcm*sthcm*cphicm)/sthp;
  sphip = (sthpcm*cphipcm*cthcm*sphicm + sthpcm*sphipcm*cphicm + cthpcm*sthcm*sphicm)/sthp;
  xabs  = abs(cphip);
  if(xabs > 1.0)
    {
      cphip = cphip/xabs;
      sphip = 0.0;
    }

  Phipp = acos(cphip);  // between 0 and PI --> check sin next
  if(sphip < 0.0) Phipp = 2.*pi(1) - Phipp; // to be between 0 and 2PI
 
  out[0] = t;
  out[1] = nup;
  out[2] = Epp;
  out[3] = thetaqp;
  out[4] = Phiqp;
  out[5] = thetapp;
  out[6] = Phipp;

  return 0;
}
