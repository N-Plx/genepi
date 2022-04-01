#include "genepi.h"
#include "dvcs_vars.h"
#include "inl_funcs.h"


NuclFF::NuclFF()
{};


NuclFF::~NuclFF()
{};


double NuclFF::dirac_ff(double t, int ic)
{
/*
  t input  var. 4-momentum transfer (GeV2)
  ic nucleon index (0: neut, 1: prot)
  Calculation of the Dirac form factors according to the 
  platchkov parametrization.
*/
  double Gm0[2]={-1.9130427, 2.792847337};//{\mu_n, \mu_p}
  double mv(0.843);
  double mv2=sqr(mv);
  double tau = t/(4.*m_targ(ic,2));
  double shape = 1./sqr(1. + t/mv2);
  double F;

  //Nucleon electric form factor
  double Ge(shape);
  if(ic == 0)
  {
    Ge = -1.25*Gm0[ic]*shape*tau/(1. + 18.3*tau);
  }

  //Nucleon magnetic form factor
  double Gm = Gm0[ic]*shape;

  //Dirac nucleon form factors
  F = Gm - (Gm - Ge)/(1. + tau);

  return F;
}


double NuclFF::pauli_ff(double t, int ic)
{
/*
  t input  var. 4-momentum transfer (GeV2)
  ic nucleon index (0: neut, 1: prot)
  Calculation of the Pauli form factors according to the 
  platchkov parametrization.
*/
  double Gm0[2]={-1.9130427, 2.792847337};//{\mu_n, \mu_p}
  double mv(0.843);
  double mv2=sqr(mv);
  double tau = t/(4.*m_targ(ic,2));
  double shape = 1./sqr(1. + t/mv2);
  double F;

  //Nucleon electric form factor
  double Ge(shape);
  if(ic == 0)
  {
    Ge = -1.25*Gm0[ic]*shape*tau/(1. + 18.3*tau);
  }

  //Nucleon magnetic form factor
  double Gm = Gm0[ic]*shape;

  //Pauli nucleon form factors
  F = (Gm - Ge)/(1. + tau);

  return F;
}
