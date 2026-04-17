#ifndef __GLOBALS_H
#define __GLOBALS_H

#include "colors.h"

#include "TTree.h"

#include "GRETINA.h"
#include "GRETINAWavefunction.h"

/*************************************************/
/****            Global variables             ****/
/*************************************************/

/*------ GRETINA DATA STRUCTURES ------*/

extern globalHeader gHeader;
extern globalHeader gHeaderOUT;

extern cloverPacket clover;
extern cloverEvent *cloverEventOUT;

extern GRawPacket gMode3;
extern GRawEvent *gMode3Event;

extern mode2StructOld gMode2reallyOld;
extern mode2Struct gMode2old;
extern mode2StructNew gMode2;

extern trackedGamma gMode1;

/* Waveforms...easier as globals */
extern GRETINAWF *gWf;

extern Float_t WFbaseline;
extern Float_t WFrunningBaseline;
extern Int_t WFid;
extern Float_t WFenergy;

extern GRETINA *gret;

extern g4Sim_abcd1234 g4sim;
extern g4SimOUT *gSimOUT;

/*------ ROOT TREES ------*/
extern TTree *teb;
extern TTree *wave;
extern TTree *scaler;

#endif
