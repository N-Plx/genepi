typedef struct lujets{
  int   n;
  int   k[5][4000];
  float p[5][4000];
  float v[5][4000];    
};
extern "C" lujets lujets_;

typedef struct ludat1{
  int mstu[200];
  int mstj[200];
  float paru[200];
  float parj[200];
};
extern "C" ludat1 ludat1_;

extern "C" void luexec_();
extern "C" void lulist_(int&);
extern "C" void lu1ent_(int&,int&,float&,float&,float&);
