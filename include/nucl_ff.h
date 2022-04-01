class NuclFF
{
public:
  NuclFF();
  ~NuclFF();
  double dirac_ff(double, int);
  double pauli_ff(double, int);

private:
  int    Ipn; //(0: neutron; 1: proton)
  double M_TARG;
  double M_TARG2;

public:
  double FF1[2];
  double FF2[2];
  int    get_Ipn()      {return Ipn;}
  double get_FF1(int i) {return FF1[i];}
  double get_FF2(int i) {return FF2[i];}
  void   set_Ipn(int i) {Ipn = i;}
  void   set_FF1(int i, double ff1) {FF1[i] = ff1;}
  void   set_FF2(int i, double ff2) {FF2[i] = ff2;}
};
//FF1[0] //Dirac neutron
//FF1[1] //Dirac proton
//FF2[0] //Pauli neutron
//FF2[1] //Pauli proton
