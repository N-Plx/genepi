#include "genepi.h" 

//REFERENCE arXiv:nucl-th/9507024

double nitro_fermi_dist(double x) 
{

  //parametrization of the nucleon momentum distributions
  //see reference appendix
  double xx = pow((x*5.07),2);
  //coeff for nitrogen are a linear interpolation of carbon and oxygen coeff
  const double coeff[7] = {2.675, 2.995, 5.1, 0.376, 1.5, 0.025, 0.22};
  double n0 = coeff[0]*exp(-coeff[1]*xx)*(1+coeff[2]*xx);
  double n1 = coeff[3]*exp(-coeff[4]*xx)+coeff[5]*exp(-coeff[6]*xx);
  return (n0+n1)*xx; 
} 

// The generator is taken by reference so that the caller's random stream advances.
double fer_mom_nitro(TRandom1& randN) 
{   
  // Envelope of the accept-reject: must be >= max of nitro_fermi_dist, which is
  // 1.1854 at p = 0.152 GeV (1.10522 clipped the peak, 129-175 MeV/c, by 1.24%)
  const double fmax = 1.19;
  const double xmin = 0.; 
  const double xmax = 4/5.07;
  while(true)
    {
      const double x = randN.Uniform(xmin, xmax);
      const double fx_test = randN.Uniform(0,fmax);
      if(fx_test <= nitro_fermi_dist(x)) return x;
    }
}

