#include "Globals.h"

/*************************************************/
/****            Global variables             ****/
/*************************************************/

/*------ GRETINA DATA STRUCTURES ------*/

globalHeader gHeader;
globalHeader gHeaderOUT;

/* Waveforms...easier as globals */
GRETINAWF *gWf;

Float_t WFbaseline;
Float_t WFrunningBaseline;
Int_t WFid;
Float_t WFenergy;

GRETINA *gret;

Track *track;

/*------ CHICO DATA STRUCTURES ------*/
#ifdef WITH_CHICO
CHICOFull *chico;
#endif

/*------ GODDESS DATA STRUCTURES -------*/
#ifdef WITH_GOD
goddessFull *goddess;
#endif

/*------ DFMA DATA STRUCTURES -------*/
#ifdef WITH_DFMA
DFMAFull *dfma;
#endif

/*------ PWALL DATA STRUCTURES ------*/
#ifdef WITH_PWALL
phosWallFull *phosWall;
#endif

/*------ S800 DATA STRUCTURES ------*/

#ifdef WITH_S800
S800Full *s800;
S800Scaler *s800Scaler;
#endif

/*------ LENDA DATA STRUCTURES ------*/

#ifdef WITH_LENDA
ddasEvent *ddasEv;
lendaEvent *lendaEv;
lendaPacker *lendaPack;
lendaSettings *lendaSet;
#endif

/*------ ROOT TREES ------*/
TTree *teb; 
TTree *wave; 
TTree *scaler;
