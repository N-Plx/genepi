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

double fer_mom_nitro(TRandom1 randN) 
{   
  const double fmax = 1.10522;
  const double xmin = 0.; 
  const double xmax = 4/5.07;
  const double x = randN.Uniform(xmin, xmax);
  const double fx_test = randN.Uniform(0,fmax);
  const double fx = nitro_fermi_dist(x);
  if(fx_test > fx) return fer_mom_nitro(randN);  
  return x; 
}

