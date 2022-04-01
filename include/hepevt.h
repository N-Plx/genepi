//geant4/source/event/include/G4HEPEvtInterface.hh
#define NMXHEP 2000
typedef struct HEPEVT
{
  int    NEVHEP;
  int    NHEP;
  int    ISTHEP[NMXHEP];
  int    IDHEP[NMXHEP];
  int    JMOHEP[NMXHEP][2];
  int    JDAHEP[NMXHEP][2];
  double PHEP[NMXHEP][5];
  double VHEP[NMXHEP][4];
};
extern HEPEVT hepevt;
