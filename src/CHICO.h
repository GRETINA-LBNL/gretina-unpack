/*! \file CHICO.h
    \brief Parameter and function definitions for CHICO analysis.

    This file provides the data structures and function prototypes for CHICO
    analysis, based on an analysis code CheckChico.cc from S. Zhu (ANL).

    Author: H.L. Crawford
    Date: April 29, 2014
*/

#ifndef __CHICO_H
#define __CHICO_H

#include "CHICODefinitions.h"
#include "colors.h"

#include "TObject.h"
#include "TMath.h"
#include "TString.h"
#include <cstring>
#include <vector>
#include <iostream>
#include <RConfigure.h>

#include "TFile.h"
#include "TCutG.h"

class CHICOFull;

class CHICORaw : public TObject {
 private:
  CHICOFull* m_top;
  
 public:
  unsigned long long int LEDts;
  UShort_t cathode_tdc_num;
  Int_t cathode_tdc_val[128];
  UShort_t cathode_tdc_ch[128];
  UShort_t anode_tdc_num;
  Int_t anode_tdc_val[128];
  UShort_t anode_tdc_ch[128];
  UShort_t anode_qdc_num;
  Int_t anode_qdc_val[32];
  Short_t anode_qdc_ch[32];
  UInt_t status;
  Int_t RF;

 public:
  void Initialize();
  void Reset();
  
  ClassDef(CHICORaw, 1);
};

class CHICOParticle : public TObject {
 private:
  CHICOFull* m_top;

 public:
  unsigned long long int LEDts; /*2018-05-10 CMC added*/

  Int_t id;
  ULong64_t t;
  Double_t tof; /* left - right */
  Double_t mass;

  Int_t thetaL;  Int_t phiL;
  Int_t thetaR;  Int_t phiR;
  /*
  Float_t fThetaL;  Float_t fPhiL;
  Float_t fThetaR;  Float_t fPhiR;
  */
  Double_t fThetaL; Double_t  fPhiL;
  Double_t fThetaR; Double_t  fPhiR;

  Int_t eL;  Int_t eR;

  Double_t rf;

  /*  Float_t pgCosL, pgCosR;*/
  Double_t pgCosT, pgCosP;

  Double_t dL, dR;
  Double_t massT, massP, betaT, betaP, Qvalue;

 public:
  void Initialize();
  void Reset();
  
  ClassDef(CHICOParticle, 1);
};

class CHICOFull : public TObject {

 public:
  CHICORaw raw;
  CHICOParticle particle;

  UInt_t gotParticle;
  UInt_t idParticle;
  UInt_t multiAnodeTDCNum;
  UInt_t multiCathodeTDCNum;

  /* Calibration parameters */
  Float_t thetaOffset[PPAC_NUM];  Float_t thetaGain[PPAC_NUM];
  Float_t thetaQuad[PPAC_NUM];
  Float_t phiOffset[PPAC_NUM];  Float_t phiGain[PPAC_NUM];
  Float_t betaT[160]; Float_t betaP[160];
  // Float_t at[ANGNUM], bt[ANGNUM], ct[ANGNUM];
  // Float_t am[ANGNUM], bm[ANGNUM], cm[ANGNUM];
  Float_t betaP0,betaP1,betaP2,betaP3,betaP4,betaP5;
  Float_t betaT0,betaT1,betaT2,betaT3,betaT4,betaT5;
  Float_t betaB[180];

  /* CHICO data */
  UInt_t rawCHICO[8192];

  /* GRETINA target offset with CHICO */
  Float_t offsetTarget;

  /* CHICO cuts */
  TCutG *massLCut,*massRCut;
  TCutG *etCut,*etCutTK;
  TCutG *GT_TKCut;
  TCutG *ThetaLRCut;

 public:
  CHICOFull();
  
 public:
  void Initialize();
  void Reset();
  void InitializeCHICOVariables(TString thetaCalFile, TString phiCalFile,
				TString betaCalFile);
  void ReadThetaCalibration(TString thetaCalFile);
  void ReadPhiCalibration(TString phiCalFile);
  void ReadBetaCalibration(TString betaCalFile);
  Int_t getAndUnpackCHICO(FILE *inf, Int_t length);
  void GetParticle();
  Float_t calcCos(Float_t pTheta, Float_t pPhi, Float_t gTheta, Float_t gPhi);
  
  void setupCut();
  Int_t IDparticles();
  Float_t GetMass(int LR, int PT);
  Float_t GetBeta(int LR, int PT, float MASS);
  Float_t GetQval(int LR);

 private:
  
  ClassDef(CHICOFull, 1);
};

#endif
