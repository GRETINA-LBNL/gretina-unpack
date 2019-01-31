#ifndef __HISTOS_H
#define __HISTOS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "Globals.h"
#include "Unpack.h"
#include "GRETINA.h"
#include "SortingStructures.h"

class Histos : public TObject {
 public:
  /* ROOT histograms */
  TH1F *cc1cal[MAXCRYSTALS];  TH1F *cc2cal[MAXCRYSTALS];
  TH1F *cc3cal[MAXCRYSTALS];  TH1F *cc4cal[MAXCRYSTALS];
  TH1F *ccTot[4];  TH1F *ccDop[4];
  
  TH2F *ccSummary;  TH2F *ccDopSummary;
  
  TH1F *mult;
  TH1F *ccHit;
  
  TH2F *ggMat;  TH2F *ggMatDop;
  
  TH2F *baseVsE[MAXCRYSTALS];
  TH1I *bankTSdiff[MAXCRYSTALS], *bankTSlow[MAXCRYSTALS], *bankTShigh[MAXCRYSTALS];
  
  TH2F *xyProj, *xzProj, *yzProj;
  TH3F *xyz;
  TH2F *thetaPhi, *thetaE, *thetaEDop;
  
  TH2F *phiVsE[MAXCRYSTALS];
  
  TH2F *chi2VsXtal;
  TH1F *deltaT;
  
  TH1F *test;
  
  /* Energy calibration histograms */
  TH1F *eraw[MAXCHANNELS], *ecal[MAXCHANNELS];
  TH1F *cc1All, *cc2All, *cc3All, *cc4All;
  TH2F *cc1Summary, *cc2Summary, *cc3Summary, *cc4Summary;
  TH2F *xtalSummary[MAXCRYSTALS];
  TH2F *smoke[112];
  
  /* XTalk calibration histograms */
  TH1F *erawXT[40], *ecalXT[40];
  TH1F *mult1Singles[40];
  TH1F *mult2Pair[40][40]; TH2F *mult2Pair_2D[40][40];
  TH2F *mult1CC[40];
  TH1F *mult1, *mult2, *mult3, *mult4, *mult5;
  TH1F *mult1cal, *mult2cal, *mult3cal, *mult4cal, *mult5cal;
  TH1F *segMult, *segSum, *segSumCal;

 public:
  Histos() { ; }
  ~Histos() { ; }
  void InitializeHist(controlVariables *ctrl);
  void WriteHistos(controlVariables *ctrl);
  
 private:
  ClassDef(Histos, 1);
};

#endif
