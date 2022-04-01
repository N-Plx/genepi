#include <vector>
using namespace std;

class VECT
{
public:
  VECT();
  ~VECT();
  vector<double> init_vect();
  double mom_vect(vector<double>);
  vector<double> fill_vect(double, double, double);
  vector<double> add_vect(vector<double>, vector<double>);
  vector<double> diff_vect(vector<double>, vector<double>);
  vector<double> fill_matr(double, double);
  vector<double> rot_vect(vector<double>, double, double);
  vector<double> fill_invmatr(double, double);
  vector<double> invrot_vect(vector<double>, double, double);
};
