#include "nucl_ff.h"
#include "read_optFile.h"
#define MAXSK   51

class DVCS
{
public:
  DVCS();
  ~DVCS();
  void init_dvcs();
  int  read_gpds(string);
  int  eval_cffs(ReadOptFile*, NuclFF*, double, double);
  int  bh_xsec(ReadOptFile*, double*, double, double, double, double, double);
  int  dvcs_xsec(ReadOptFile*, double, double, double, double, double);
  int  int_xsec(ReadOptFile*, double*, double, double, double, double, double, double);
  int  phot_xsec(ReadOptFile*, NuclFF*, double, double, double, double, double);

private:
  double re_Eu[MAXSK][5];
  double re_Ed[MAXSK][5];
  double re_Etu[MAXSK][5];
  double re_Etd[MAXSK][5];
  double re_Hu[MAXSK][5];
  double re_Hd[MAXSK][5];
  double re_Htu[MAXSK][5];
  double re_Htd[MAXSK][5];
  double im_Eu[MAXSK][5];
  double im_Ed[MAXSK][5];
  double im_Etu[MAXSK][5];
  double im_Etd[MAXSK][5]; 
  double im_Hu[MAXSK][5];
  double im_Hd[MAXSK][5];
  double im_Htu[MAXSK][5];
  double im_Htd[MAXSK][5]; 
  double re_HuNF[MAXSK-1][21][2];
  double re_HdNF[MAXSK-1][21][2];
  double im_HuNF[MAXSK-1][21][2];
  double im_HdNF[MAXSK-1][21][2];

public:
  double Re_E;
  double Im_E;
  double Re_Et;
  double Im_Et;

  double Re_H;
  double Im_H;
  double Re_Ht;
  double Im_Ht;

  double hc0_BH;
  double hc1_BH;
  double hc2_BH;

  double hc0_BH_LP;
  double hc1_BH_LP;

  double hc0_DVCS;
  double hc1_DVCS;
  double hs1_DVCS;

  double hc0_DVCS_LP;
  double hc1_DVCS_LP;
  double hs1_DVCS_LP;

  double hc0_Int;
  double hc1_Int;
  double hc2_Int;
  double hs1_Int;
  double hs2_Int;

  double hc0_Int_LP;
  double hc1_Int_LP;
  double hc2_Int_LP;
  double hs1_Int_LP;
  double hs2_Int_LP;

  double get_re_Eu(int i,int j)  {return re_Eu[i][j];}
  double get_re_Ed(int i,int j)  {return re_Ed[i][j];}
  double get_re_Etu(int i,int j) {return re_Etu[i][j];}
  double get_re_Etd(int i,int j) {return re_Etd[i][j];}
  double get_re_Hu(int i,int j)  {return re_Hu[i][j];}
  double get_re_Hd(int i,int j)  {return re_Hd[i][j];}
  double get_re_Htu(int i,int j) {return re_Htu[i][j];}
  double get_re_Htd(int i,int j) {return re_Htd[i][j];}
  double get_im_Eu(int i,int j)  {return im_Eu[i][j];}
  double get_im_Ed(int i,int j)  {return im_Ed[i][j];}
  double get_im_Etu(int i,int j) {return im_Etu[i][j];}
  double get_im_Etd(int i,int j) {return im_Etd[i][j];}
  double get_im_Hu(int i,int j)  {return im_Hu[i][j];}
  double get_im_Hd(int i,int j)  {return im_Hd[i][j];}
  double get_im_Htu(int i,int j) {return im_Htu[i][j];}
  double get_im_Htd(int i,int j) {return im_Htd[i][j];}
  double get_re_HuNF(int i,int j,int k) {return re_HuNF[i][j][k];}
  double get_re_HdNF(int i,int j,int k) {return re_HdNF[i][j][k];}
  double get_im_HuNF(int i,int j,int k) {return im_HuNF[i][j][k];}
  double get_im_HdNF(int i,int j,int k) {return im_HdNF[i][j][k];}
};
