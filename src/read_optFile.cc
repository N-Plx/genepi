#include "cpp.h"
#include "options.h"
#include "read_optFile.h"


ReadOptFile::ReadOptFile()
{};


ReadOptFile::~ReadOptFile()
{};


int ReadOptFile::ReadInputFile(const char* inpfile)
{
  cout<<endl;
  cout<<"------ reading input file : "<<inpfile<<endl;
  parse_file(inpfile);
 
  get_opt("GENE","EVENTS", fNevts);
  cout<< "number of events to be generated: " << fNevts << endl;

  get_opt("BEAM","ENERGY", fEb);
  cout<< "beam energy: " << fEb <<" GeV"<<endl;

  get_opt("BEAM","CHRG", fBchg);
  cout<< "beam charge: " << fBchg <<endl;

  get_opt("BEAM","HELI", fBheli);
  if(abs(fBheli) > 1)
  {
    cout<<"unknown beam helicity. It should be -1, 0 or +1"<<endl;
    return -1;
  }
  cout<<"beam helicity: " << fBheli <<endl;

  get_opt("TARGET","A", fAt);

  get_opt("TARGET","Z", fZt);

  int iApZ = (int)(fAt + fZt);
  string target;
  if(iApZ == 1)
  {
    target = "neut";
  }
  else if(iApZ == 2)
  {
    target = "prot";
  }
  else if(iApZ == 3)
  {
    target = "deut";
  }
  else if(iApZ == 8)
  {
    target = "hel4";
  }
  else if(iApZ == 21)   
    {     
      target = "nitro";
    }
  else if(iApZ == 27)
    {
      target = "NH3";
    }
  else if(iApZ == 30)
    {
      target = "ND3";
    }
  else
  {
    cout<<"invalid target"<<endl;
    return -1;
  }
  cout<<"target type :  "<< target<<endl;

  get_opt("TARGET","HELI", fTheli);
  if(abs(fTheli) > 1)
  {
    cout<<"unknown target helicity. It should be -1, 0 or +1"<<endl;
    return -1;
  }
  cout<< "target helicity: " << fTheli <<endl;

  vector<double> ftmp;
  get_opt("KINEMATIC","XBJ", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad xbj range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  fXbjmin = ftmp[0]; fXbjmax = ftmp[1];
  ftmp.clear();

  get_opt("KINEMATIC","Y", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad y range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  fYmin = ftmp[0]; fYmax = ftmp[1];
  ftmp.clear();

  get_opt("KINEMATIC","Q2", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad Q2 range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  fQ2min = ftmp[0]; fQ2max = ftmp[1];
  ftmp.clear();

  get_opt("KINEMATIC","W2", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad W2 range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  fW2min = ftmp[0]; fW2max = ftmp[1];
  ftmp.clear();

  get_opt("KINEMATIC","NU", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad nu range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  fNumin = ftmp[0]; fNumax = ftmp[1];
  ftmp.clear();

  get_opt("KINEMATIC","T", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad t range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  ftmin = ftmp[0]; ftmax = ftmp[1];
  ftmp.clear();

  get_opt("KINEMATIC","YCOL", ftmp);
  if(ftmp.size() != 2)
  {
    cerr<<"bad ycol range"<<endl;
    ftmp.clear();
    ftmp.push_back(0);
    ftmp.push_back(100);
  }
  fYcolmin = ftmp[0]; fYcolmax = ftmp[1];
  ftmp.clear();

  get_opt("XSEC","DIFF", fVol);
  if(fVol == 1)
  {
    cout<<"xsec Elementary volume : d5sigma /(dx_B dy dt dPhi_e dPhi_g)"<<endl;
  }
  else if(fVol == 2)
  {
    cout<<"xsec Elementary volume : d5sigma /(dk_e dOmega_e dOmega_g)"<<endl;
  }
  else
  {
    cout<<"xsec Elementary volume : d5sigma /(dQ2 dx_B dt dPhi_e dPhi_g)"<<endl;
  }

  get_opt("PROD",   "DIR",          fDir);
  cout<<"output directory "<<fDir<<endl;

  get_opt("APPEND", "DATE",         fDate);
  if(fDate == 1) cout<<"append time string to the output directory : yes"<<endl;
  if(fDate == 0) cout<<"append time string to the output directory : no"<<endl;

  get_opt("NTUP",   "FILL",         fNtp);

  get_opt("NTUP",   "EVTSPERNTUP",  fNevtsPerNtup);

  if(fNtp == 0) cout<<"producing root ntuples : no"<<endl;
  if(fNtp == 1) 
  {
    cout<<"producing root ntuples : yes"<<endl;
    cout<<"number of events per ntuple : "<<fNevtsPerNtup<<endl;
  }

  get_opt("PROD",   "ASCII",    fAscii);

  get_opt("PROD",   "FORMAT",   fMode);

  get_opt("OUTPUT","EVTSPERFILE", fNevtsPerFile);

  if(fAscii == 1)
  {
    string sMode;
    if(fMode == 0)
    {
      sMode = "HEPEVT";
    }
    else
    {
      sMode = "LUJETS";
    }
    cout<<"producing text file in "<<sMode<<" format : yes"<<endl;
    cout<<"number of events per ASCII file : "<<fNevtsPerFile<<endl;
  }
  if(fAscii == 0) cout<<"producing text file : no"<<endl;

  get_opt("OUTPUT","PRINT",       fPrint);
  cout<<"print event info each "<<fPrint<<" events"<<endl;

  get_opt("SEED",  "APPLY",       fSeed);
  if(fSeed == 1)
  {
    cout<<"seed is set to the current  machine clock"<<endl;
  }
  else
  {
    cout<<"seed ="<<fSeed<<endl;
  }

  get_opt("PROC",  "TYPE",        fProc);

  get_opt("GPD",  "TYPE",         fMgpd);

  get_opt("MESON", "TYPE",        fIms);

  if(fProc == 0)
  {
    cout<<"DVCS process"<<endl;
    if(fMgpd >= 0 && fMgpd <= 4)
    {
      cout<<"GPD model "<<fMgpd<<endl;
    }
    else
    {
      cout<<"unknown GPD model"<<endl;
      return -1;
    }
  }
  else if(fProc == 1)
  {
    if(fIms == 0)
    {
      cout<<"pi0 electroproduction process"<<endl;
    }
    else if(fIms == 1)
    {
      cout<<"eta electroproduction process"<<endl;
    }
    else if(fIms == 2)
      {
	cout<<"phi electroproduction process"<<endl;
      }
    else
    {
      cout<<"unknown meson"<<endl;
      return -1;
    }
  }
  else
  {
    cout<<"unknown process"<<endl;
    return -1;
  }

  get_opt("RUN","NUMBER", fRunnum);
  cout<< "Run number: " << fRunnum << endl;

  get_opt("VERTEX","X", fVx);
  cout<< "Vertex x : " << fVx << endl;
  get_opt("VERTEX","Y", fVy);
  cout<< "Vertex y : " << fVy << endl;
  get_opt("VERTEX","Z", fVz);
  cout<< "Vertex z : " << fVz << endl;
  get_opt("RASTER","X", fRasterx);
  cout<< "Raster x : " << fRasterx << endl;
  get_opt("RASTER","Y", fRastery);
  cout<< "Raster y : " << fRastery << endl;

  cout<<"------ reading input file : DONE WITH SUCCESS ---------"<<endl;

  cout<<endl;

  return 0;
}
