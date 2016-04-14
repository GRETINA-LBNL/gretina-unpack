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

/*------ PWALL DATA STRUCTURES ------*/
#ifdef WITH_PWALL
phosWallFull *phosWall;
#endif

/*------ S800 DATA STRUCTURES ------*/

#ifdef WITH_S800
S800Full *s800;
S800Scaler *s800Scaler;
#endif

/*------ ROOT TREES ------*/
TTree *teb; 
TTree *wave; 
TTree *scaler;
