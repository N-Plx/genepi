#include "cpp.h"
#include "root.h"

double sqr(double);
double sign(double);
int    sign(int);

int    fexist(const char*);
void   init_tree(TTree *t);
void   init_track();
void   init_event();

double get_theta(double, double);
double get_phi(double, double);
double get_ycol(double, double, double);
void   get_matr(double, double, double matr[3][3]);
void   dump_file(int, int, double);
//void   dump_file(int, double, FILE*);
int    get_seed();
string get_date();
void   prod_dir(string);

double fer_mom_deut(int, double);
double fer_mom_hel4(double);
double fer_mom_nitro(TRandom1);
int    fmotion(double, double, int, int, double, double, double, double, double, double, double, double, double*);

double ups(double);
double ums(double);
double dps(double);
double dms(double);
double get_tmin(int, double, double);
double get_K2(int, double, double, double, double);
double get_J(int, double, double, double, double);

double m_targ(int, double);
int    targ_id(int);
int    targ_ch(int);

double m_ms(int, double);
int    ms_id(int);
