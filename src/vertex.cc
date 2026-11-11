#include"vertex.h"
#include <cmath>
#include "TRandom.h"
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

// Vertex uniformly distributed inside the ellipse of half-axes (raster_x, raster_y)
// centred on (vx, vy), at z = vz. Units: cm.
// - r = sqrt(u) gives a uniform density in area (r = u piles events up at the centre).
// - Uses gRandom, seeded from --seed in main(): a new TRandom1 per call ignored the seed.
void vertex_xyz(double* x, double* y, double* z, double vx, double vy, double vz, double raster_x, double raster_y){

  double r     = sqrt(gRandom->Rndm());
  double theta = 2*M_PI*gRandom->Rndm();
  *x = vx + raster_x*r*cos(theta);
  *y = vy + raster_y*r*sin(theta);
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
