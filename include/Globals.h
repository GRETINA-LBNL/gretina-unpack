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

/* Waveforms...easier as globals */
extern GRETINAWF *gWf;

extern Float_t WFbaseline;
extern Float_t WFrunningBaseline;
extern Int_t WFid;
extern Float_t WFenergy;

extern GRETINA *gret;

/*------ ROOT TREES ------*/
extern TTree *teb;
extern TTree *wave;
extern TTree *scaler;

#endif
