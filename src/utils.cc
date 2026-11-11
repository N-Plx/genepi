#include "genepi.h"
#include "inl_funcs.h"
#include <random>
#include <chrono>
#include <functional>
#include <unistd.h>

int sign(int x)
{
  int sign(0);
  if(x>=0) {
    sign = 1;
  }
  else{
    sign = -1;
  }
  return sign;
}


double sign(double x)
{
  double sign(0);
  if(x>=0) {
    sign = 1;
  }
  else{
    sign = -1;
  }
  return sign;
}


double sqr(double x)
{
  return x*x;
}


double get_theta(double z, double norm)
{
  double theta;
  if(norm != 0)
  {
    theta = acos(z/norm);
  }
  else
  {
    theta = -1000;
  }

  return theta;
}


double get_phi(double x, double y)
{
  double phi;
  if(x != 0 && y != 0)
  {
    phi   = atan(abs(y)/abs(x));
    if(x > 0. && y > 0.)
    {
      phi = phi;
    }
    else if(x < 0. && y > 0.)
    {
      phi   = pi(1) - phi;
    }
    else if(x < 0. && y < 0.)    
    {
      phi   = pi(1) + phi;
    }
    else if(x > 0. && y < 0.) 
    {
      phi   = 2.*pi(1) - phi;
    }
  }
  else if(x == 0 && y != 0)
  {
    if(y > 0.)
    {
      phi = pi(1)/2.;
    }
    else if(y < 0.)
    {
      phi = 3.*pi(1)/2.;
    }
  }
  else if(x != 0 && y == 0)
  {
    if(x > 0.)
    {
      phi = 0.;
    }
    else if(x < 0.)
    {
      phi = pi(1);
    }
  }
  else if(x == 0 && y == 0)
  {
    //phi = -1000.;
    phi = 0.;
  }

  return phi;
}


int fexist(const char *filename)
{
  struct stat buffer;
  if(stat(filename, &buffer)){
    return 1;
  } 
  else{
    return 0;
  }
}


template <typename T>
string to_str(T const& value)
{
  stringstream sstr;
  sstr << value;
  return sstr.str();
}


string get_date()
{
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[100];

  time (&rawtime);
  timeinfo = localtime (&rawtime);

  strftime(buffer, 20, "_%y%m%d_%H%M%S", timeinfo);
  string str = to_str(buffer);

  return str;
}


// Default seed when --seed is not given. The previous time-of-day seed (HHMMSS)
// gave identical events for jobs started in the same second, e.g. on a grid.
// Mix std::random_device, the high-resolution clock, the process id and the host
// name instead. Result in [1, 900000000] (also a valid JETSET seed).
unsigned int get_seed()
{
  const unsigned long long k = 0x9E3779B97F4A7C15ULL;
  std::random_device rd;
  unsigned long long h = rd();
  h = h*k ^ (unsigned long long)std::chrono::high_resolution_clock::now().time_since_epoch().count();
  h = h*k ^ (unsigned long long)getpid();
  char host[256] = {0};
  gethostname(host, sizeof(host) - 1);
  h = h*k ^ (unsigned long long)std::hash<std::string>()(std::string(host));
  h ^= h >> 31;
  return (unsigned int)(h % 900000000ULL) + 1;
}


void prod_dir(string dir)
{
  string cmd1("mkdir -p "+dir);
  system(cmd1.c_str());
}


//reference:
//QCD Constrains on the Shape of Polarized Quark and Gluon Distributions
//Stanley J. Brodsky, Matthias Burkardt, Ivan Schmidt
//Nucl.Phys. B441 (1995) 197-214 hep-ph/9401328
double ups(double x)
{
  double Uplus = pow(1./x,1.12)*(3.784*pow((1. - x),3) - 3.672*pow((1. - x),4));
  return Uplus;
}


double ums(double x)
{
  double Umnus = pow(1./x,1.12)*(2.004*pow((1. - x),5) - 1.892*pow((1. - x),6));
  return Umnus;
}


double dps(double x)
{
  double Dplus = pow(1./x,1.12)*(0.757*pow((1. - x),3) - 0.645*pow((1. - x),4));
  return Dplus;
}


double dms(double x)
{
  double Dmnus = pow(1./x,1.12)*(3.23*pow((1. - x),5) - 3.118*pow((1. - x),6));
  return Dmnus;
}


double m_targ(int i, double y)
{
  double mt;
  if(i == 0) //neutron
  {
    mt = m_neut(y);
  }
  else if(i == 1) //proton
  {
    mt =  m_prot(y);
  }
  return mt;
}


int targ_id(int i)
{
  int id;
  if(i == 0) //neutron
  {
    id = neut_id();
  }
  else if(i == 1) //proton
  {
    id = prot_id();
  }
  return id;
}


int targ_ch(int i)
{
  int ich;
  if(i == 0) //neutron
  {
    ich = neut_ch();
  }
  else if(i == 1) //proton
  {
    ich = prot_ch();
  }
  return ich;
}


double m_ms(int i, double y)
{
  double ms;
  if(i == 0) //pi0
  {
    ms = m_pi0(y);
  }
  else if(i == 1) //eta
  {
    ms = m_eta(y);
  }
  return ms;
}


int ms_id(int i)
{
  int id;
  if(i == 0) //pi0
  {
    id = pi0_id();
  }
  else if(i == 1) //eta
  {
    id = eta_id();
  }
  return id;
}
