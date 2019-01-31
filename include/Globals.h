#ifndef __GLOBALS_H
#define __GLOBALS_H

#include "colors.h"

#include "TTree.h"

#include "GRETINA.h"
#include "GRETINAWavefunction.h"

#include "Track.h"

#ifdef WITH_CHICO
#include "CHICO.h"
#include "CHICODefinitions.h"
#endif

#ifdef WITH_GOD
#include "GODDESS.h"
#endif

#ifdef WITH_DFMA
#include "DFMA.h"
#endif

#ifdef WITH_PWALL
#include "PhosWall.h"
#endif

#ifdef WITH_S800
#include "S800Parameters.h"
#include "S800Definitions.h"
#endif

#ifdef WITH_LENDA
#include "LENDA-DDAS.h"
#include "LENDA-Controls.h"
#endif

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

extern Track *track;

/*------ CHICO DATA STRUCTURES ------*/
#ifdef WITH_CHICO
extern CHICOFull *chico;
#endif

/*------ GODDESS DATA STRUCTURES -------*/
#ifdef WITH_GOD
extern goddessFull *goddess;
#endif

/*------ DFMA DATA STRUCTURES -------*/
#ifdef WITH_DFMA
extern DFMAFull *dfma;
#endif

/*------ PWALL DATA STRUCTURES ------*/
#ifdef WITH_PWALL
extern phosWallFull *phosWall;
#endif

/*------ S800 DATA STRUCTURES ------*/
#ifdef WITH_S800
extern S800Full *s800;
extern S800Scaler *s800Scaler;
#endif

/*------ LENDA DATA STRUCTURES ------*/
#ifdef WITH_LENDA
extern ddasEvent *ddasEv;
extern lendaEvent *lendaEv;
extern lendaPacker *lendaPack;
extern lendaSettings *lendaSet;
#endif

/*------ ROOT TREES ------*/
extern TTree *teb;
extern TTree *wave;
extern TTree *scaler;

#endif
