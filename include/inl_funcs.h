#include<cmath>

inline double m_neut(double);
inline double m_neut(double y) {return pow(0.93957,y);}

inline double m_prot(double);
inline double m_prot(double y) {return pow(0.93827,y);}

inline double m_eta(double);
inline double m_eta(double y) {return pow(0.54775,y);}

inline double m_kaon(double);
inline double m_kaon(double y) {return pow(0.49368,y);}

inline double m_pion(double);
inline double m_pion(double y) {return pow(0.13957,y);}

inline double m_pi0(double);
inline double m_pi0(double y) {return pow(0.13498,y);}

inline double m_elec(double);
inline double m_elec(double y) {return pow(0.00051,y);}

inline double m_phot(double);
inline double m_phot(double y) {return pow(0.00000,y);}

int inline neut_id();
int inline neut_id() {return 2112;}

int inline prot_id();
int inline prot_id() {return 2212;}

int inline eta_id();
int inline eta_id() {return 221;}

int inline kaon_id();
int inline kaon_id() {return 321;}

int inline ka0_id();
int inline ka0_id() {return 310;}

int inline pion_id();
int inline pion_id() {return 211;}

int inline pi0_id();
int inline pi0_id() {return 111;}

int inline elec_id();
int inline elec_id() {return 11;}

int inline phot_id();
int inline phot_id() {return 22;}

int inline hel3_id();
int inline hel3_id() {return 1234;}//FAKE

int inline trit_id();
int inline trit_id() {return 4321;}//FAKE

int inline remnant_id(); 
int inline remnant_id() {return 12;}//id of a neutrino

int inline neut_ch();
int inline neut_ch() {return 0;}

int inline prot_ch();
int inline prot_ch() {return 1;}

int inline elec_ch();
int inline elec_ch() {return -1;}

int inline phot_ch();
int inline phot_ch() {return 0;}

inline double alpha(double);
inline double alpha(double y) {return pow(1./137.0359998, y);}

inline double hbarc(double);
inline double hbarc(double y) {return pow(0.1973269602,y);}//GeV fm

inline double tiny();
inline double tiny() {return 1e-4;}

inline double verytiny();
inline double verytiny() {return 1e-15;}

inline double pi(double);
inline double pi(double y) {return pow(3.14159265,y);}

inline double torad();
inline double torad() {return 0.01745;}

inline double todeg();
inline double todeg() {return 57.299578;}

inline double xsecpi0();
inline double xsecpi0() {return 394.73;}
