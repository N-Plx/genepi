#include <cpp.h>
#include "vect.h"

VECT::VECT()
{;}


VECT::~VECT()
{;}


vector<double> VECT::init_vect()
{
  vector<double> v;
  v.resize(3);
  for(unsigned int i=0; i<3; i++) v.at(i) = -1000.;

  return v;
}


double VECT::mom_vect(vector<double> v)
{
  double mom = 0;
  for(unsigned int i=0; i<3; i++) mom += v.at(i)*v.at(i);

  return sqrt(mom);
}


vector<double> VECT::fill_vect(double mom, double Theta, double Phi)
{
  vector<double> v;
  v.resize(3);
  v.at(0) = mom*sin(Theta)*cos(Phi);
  v.at(1) = mom*sin(Theta)*sin(Phi);
  v.at(2) = mom*cos(Theta);

  return v;
}


vector<double> VECT::add_vect(vector<double> v1, vector<double> v2)
{
  vector<double> v;
  v.resize(3);
  for(unsigned int i=0; i<v1.size(); i++) v.at(i) = v1.at(i) + v2.at(i);

  return v;
}


vector<double> VECT::diff_vect(vector<double> v1, vector<double> v2)
{
  vector<double> v;
  v.resize(3);
  for(unsigned int i=0; i<v1.size(); i++) v.at(i) = v1.at(i) - v2.at(i);

  return v;
}


vector<double> VECT::fill_matr(double Theta, double Phi)
{
  vector<double> matr;
  matr.resize(9);

  matr.at(0) = -cos(Theta)*cos(Phi); //(1,1)
  matr.at(1) = -cos(Theta)*sin(Phi); //(1,2)
  matr.at(2) =  sin(Theta);          //(1,3)
  matr.at(3) =  sin(Phi);            //(2,1)
  matr.at(4) = -cos(Phi);            //(2,2)
  matr.at(5) =  0.0;                 //(2,3)
  matr.at(6) =  sin(Theta)*cos(Phi); //(3,1)
  matr.at(7) =  sin(Theta)*sin(Phi); //(3,2)
  matr.at(8) =  cos(Theta);          //(3,3)

  return matr;
}


vector<double> VECT::rot_vect(vector<double> v, double Theta, double Phi)
{
  vector<double> matr = fill_matr(Theta, Phi);
  vector<vector<double> > mymatr(3);
  unsigned int i, j;
  for(i=0; i<mymatr.size(); i++) mymatr[i].resize(3);
  for(i=0; i<mymatr.size(); i++)
    for(j=0; j<mymatr[i].size(); j++) mymatr[i][j] = matr.at(3*i+j);

  vector<double> vv;
  vv.resize(3);
  for(i=0; i<v.size(); i++)
  {
    double sum = 0.;
    for(j=0; j<mymatr[i].size(); j++) sum += mymatr[i][j]*v.at(j);
    vv.at(i) = sum;
  }

  return vv;
}


vector<double> VECT::fill_invmatr(double Theta, double Phi)
{
  vector<double> matr;
  matr.resize(9);

  matr.at(0) =  -cos(Theta)*cos(Phi); //(1,1)
  matr.at(1) =   sin(Phi);           //(1,2)
  matr.at(2) =   sin(Theta)*cos(Phi); //(1,3)
  matr.at(3) =  -cos(Theta)*sin(Phi); //(2,1)
  matr.at(4) =  -cos(Phi);            //(2,2)
  matr.at(5) =   sin(Theta)*sin(Phi); //(2,3)
  matr.at(6) =   sin(Theta);         //(3,1)
  matr.at(7) =   0.0;                 //(3,2)
  matr.at(8) =   cos(Theta);          //(3,3)

  return matr;
}


vector<double> VECT::invrot_vect(vector<double> v, double Theta, double Phi)
{
  vector<double> matr = fill_invmatr(Theta, Phi);
  vector<vector<double> > mymatr(3);
  unsigned int i, j;
  for(i=0; i<mymatr.size(); i++) mymatr[i].resize(3);
  for(i=0; i<mymatr.size(); i++)
    for(j=0; j<mymatr[i].size(); j++) mymatr[i][j] = matr.at(3*i+j);

  vector<double> vv;
  vv.resize(3);
  for(i=0; i<v.size(); i++)
  {
    double sum = 0.;
    for(j=0; j<mymatr[i].size(); j++) sum += mymatr[i][j]*v.at(j);
    vv.at(i) = sum;
  }

  return vv;
}
