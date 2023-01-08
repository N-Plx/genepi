#include"vertex.h"
#include <cmath>
TF1 fexp("fexp","exp(-x/10.0)",-10,10);

// vertex in unit of cm.
double vertex_x(){
  return 0.0;
//  return fexp.GetRandom();
}

double vertex_y(){
  return 0.0;
//  return fexp.GetRandom();
}

double vertex_z(){
  return -3.0;
//  return fexp.GetRandom();
}

void vertex_xyz(double* x, double* y, double* z, double vx, double vy, double vz, double raster_x, double raster_y){
  
  TRandom1 random, random2;                                                                                                                                                                                       
  double rx, ry, theta;
  rx = raster_x*random.Rndm();
  ry = raster_y*random.Rndm();
  theta = 2*M_PI*random.Rndm();
  *x = vx + rx*cos(theta);                                                                                                                                                                                        
  *y = vy + ry*sin(theta);                                                                                                                                                                                        
  *z = vz;
  //Raster position
  /*r=0.9;
  theta=M_PI/3.;
  *x = r*cos(theta);
  *y = r*sin(theta);
  *x = vertex_x();
  *y = vertex_y();
  *z = vertex_z();*/

}
