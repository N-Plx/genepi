#include "cpp.h"

class ReadOptFile
{
public:
  ReadOptFile();
  ~ReadOptFile();
  int ReadInputFile(const char* InpFile);

private:
  int    fNevts;
  double fEb;
  int    fBchg;
  int    fBheli;
  int    fAt;
  int    fZt;
  int    fTheli;
  double fXbjmin;
  double fXbjmax;
  double fYmin;
  double fYmax;
  double fQ2min;
  double fQ2max;
  double fW2min;
  double fW2max;
  double fNumin;
  double fNumax;
  double ftmin;
  double ftmax;
  double fYcolmin;
  double fYcolmax;
  int    fVol;
  string fDir;
  int    fDate;
  int    fNtp;
  int    fNevtsPerNtup;
  int    fAscii;
  int    fMode;
  int    fNevtsPerFile;
  int    fPrint;
  int    fSeed;
  int    fProc;
  int    fMgpd;
  int    fIms;
  int    fIpn;//does not exist in input file
  int    fRunnum;

public:
  int    get_fNevts()   {return fNevts;}
  double get_fEb()      {return fEb;}
  int    get_fBchg()    {return fBchg;}
  int    get_fBheli()   {return fBheli;}
  int    get_fAt()      {return fAt;}
  int    get_fZt()      {return fZt;}
  int    get_fTheli()   {return fTheli;}
  double get_fXbjmin()  {return fXbjmin;}
  double get_fXbjmax()  {return fXbjmax;}
  double get_fYmin()    {return fYmin;}
  double get_fYmax()    {return fYmax;}
  double get_fQ2min()   {return fQ2min;}
  double get_fQ2max()   {return fQ2max;}
  double get_fW2min()   {return fW2min;}
  double get_fW2max()   {return fW2max;}
  double get_fNumin()   {return fNumin;}
  double get_fNumax()   {return fNumax;}
  double get_ftmin()    {return ftmin;}
  double get_ftmax()    {return ftmax;}
  double get_fYcolmin() {return fYcolmin;}
  double get_fYcolmax() {return fYcolmax;}
  double get_fVol()     {return fVol;}
  string get_fDir()     {return fDir;}
  int    get_fDate()    {return fDate;}
  int    get_fNtp()     {return fNtp;}
  int    get_fNevtsPerNtup() {return fNevtsPerNtup;}
  int    get_fAscii()   {return fAscii;}
  int    get_fMode()    {return fMode;}
  int    get_fNevtsPerFile() {return fNevtsPerFile;}
  int    get_fPrint()   {return fPrint;}
  int    get_fSeed()    {return fSeed;}
  int    get_fProc()    {return fProc;}
  int    get_fMgpd()    {return fMgpd;}
  int    get_fIms()     {return fIms;}
  int    get_fIpn()     {return fIpn;}
  int    get_fRunnum()     {return fRunnum;}

  void set_fXbjmin(double Xbjmin) {fXbjmin = Xbjmin;}
  void set_fXbjmax(double Xbjmax) {fXbjmax = Xbjmax;}
  void set_fYmin(double Ymin)     {fYmin = Ymin;}
  void set_fYmax(double Ymax)     {fYmax = Ymax;}
  void set_fQ2min(double Q2min)   {fQ2min = Q2min;}
  void set_fQ2max(double Q2max)   {fQ2max = Q2max;}
  void set_fW2min(double W2min)   {fW2min = W2min;}
  void set_fW2max(double W2max)   {fW2max = W2max;}
  void set_fNumin(double Numin)   {fNumin = Numin;}
  void set_fNumax(double Numax)   {fNumax = Numax;}
  void set_fIpn(int Ipn)          {fIpn=Ipn;}
};

void   get_kinem_limits(ReadOptFile*);
double ms_xsec(ReadOptFile*, double, double, double, double);
int    get_ms(ReadOptFile*, double, double, double, double, double, double, double);
