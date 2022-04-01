#include "root.h"

//beam 0 =============================================================
extern double k;
extern double Ee;

//scat. elec. 1 ======================================================
extern double kp;
extern double Eep;
extern double Thetakkp ,Phikkp;

//q(virtual photon) 2 ================================================
extern double q;
extern double nu;
extern double Thetakq ,Phikq;

//initial nucleon P(Ep;P) 3 ==========================================
extern double P;
extern double E;
extern double ThetakP, PhikP;

//initial nucleon P(Ep;P) 3 ==========================================
extern double ThetaqP, PhiqP;

//target remnant Pr(EPr,Pr) 4 ========================================
extern double Pr;
extern double Er;
extern double ThetakPr, PhikPr; 

//final nucleon P'(EP',P') wrt virt. phot. frame =====================
extern double Ep;
extern double Pp;
extern double ThetaqPp, PhiqPp;

//final nucleon P'(EP',P') wrt inc. elec. frame 5 ====================
extern double ThetakPp, PhikPp;

//qp(real photon) wrt virt. phot. frame ==============================
extern double qp;
extern double nup;
extern double Thetaqqp, Phiqqp;

//qp(real photon) wrt inc. elec. frame 6 =============================
extern double Thetakqp, Phikqp;

//Pm(meson) wrt inc. elec. frame 6 =============================
extern double Em;
extern double Pm;
extern double ThetakPm, PhikPm;

//polar and azimuthal angles between scat. elec. and real photon =====
extern double Thetakpqp, Phikpqp;

// kinematics ========================================================
extern double Q2;   // photon virtuality
extern double xbj;  // xbjorken variable
extern double y;    // e- energy loss  
extern double W2;   // invariant mass
extern double t;
extern double ycol;

//Cross Section ======================================================
extern double ds_bh;
extern double ds_dvcs;
extern double ds_int;
extern double ds_tot;
extern double ds_ms;
