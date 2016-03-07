#include "GlobalsAGH.h"

/*************************************************/
/****            Global variables             ****/
/*************************************************/

/*------ GRETINA DATA STRUCTURES ------*/

globalHeader gHeader;
globalHeader gHeaderOUT;

cloverPacket clover;
cloverEvent *cloverEventOUT;

GRawPacket gMode3;
GRawEvent *gMode3Event;

/* Waveforms...easier as globals */
GRETINAWF *gWf;

Float_t WFbaseline;
Float_t WFrunningBaseline;
Int_t WFid;
Float_t WFenergy;

mode2StructOld gMode2reallyOld;
mode2Struct gMode2old;
mode2StructNew gMode2;

trackedGamma gMode1;

g4Sim_abcd1234 g4sim;
g4SimOUT *gSimOUT;

/*------ CHICO DATA STRUCTURES ------*/
#ifdef WITH_CHICO
CHICOFull *chico;
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
