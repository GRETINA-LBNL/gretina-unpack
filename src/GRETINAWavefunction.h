/**********************************************************************/
/* File: GRETINAWavefunction.h                                        */
/* Description: Parameter and function definitions for GRETINA        */
/*              waveform analysis                                     */
/* Author: H. Crawford                                                */
/* Date: January 2013                                                 */
/**********************************************************************/

#ifndef __GRETINAWF_H
#define __GRETINAWF_H

#include <stdio.h>
#include <stdlib.h>

#include "TObject.h"

#include "Defines.h"

using std::vector;

class GRETINAWF : public TObject {
 public:
  /* The waveform we usually look at */
  vector<Short_t> waveform2Out[MAXCRYSTALS][40];
  vector<Short_t> waveform2Print;
  vector<Float_t> pzWaveform;
  vector<Short_t> scfd;
  Int_t tracelength;

  /* Some wavefunction analysis variables... */
  Float_t minCrossTime, maxCrossTime;
  Int_t eFilterK, eFilterM;
  Int_t puFilterK, puFilterM, puFilterThresh;

  /* Waveform timing information. */
  Float_t LEDcrossing[MAXCRYSTALS][40];

  /* CC baseline information */
  Float_t cc1BL[MAXCRYSTALS], cc2BL[MAXCRYSTALS], cc3BL[MAXCRYSTALS], cc4BL[MAXCRYSTALS];

  /* Resting baseline values */
  Float_t restingBase[MAXCHANNELS];
  Int_t baseSamples[MAXCHANNELS];
  
  Float_t tau[MAXCHANNELS];
  Float_t baseline[MAXCHANNELS];

  /* Arrays for the RadfordBaseline algorithms. */
  Double_t sa[MAXCHANNELS], sb[MAXCHANNELS], sc[MAXCHANNELS], sd[MAXCHANNELS], se[MAXCHANNELS];
  Int_t his[MAXCHANNELS][3][4096];

  /* Segment threshold value for summing */
  Float_t segmentThresh;

  /* Timing stuff */
  Int_t minT; Int_t maxT;

  /* Ken's data structures */
  struct peak{
    Int_t pileup;
    Double_t amp;
    Double_t cen;
    Double_t base;
  };
  
  struct linfit{
    Double_t m;
    Double_t b;
    Double_t s;
  };

  /* These are vectors to hold waveforms for output tree. */
  vector<Short_t> waveformConCat;
  vector<Double_t> waveE;
  vector<Double_t> waveCE;

 public:
  GRETINAWF();
  ~GRETINAWF() { ; }
  
 public:
  void Initialize();
  void Reset();
  Int_t ANLAmplitude(Int_t T, Int_t cNum, Int_t chNum);
  Int_t ANLBaseCal(Int_t T, Int_t cNum, Int_t chNum);
  Double_t ANLEnergy1(Int_t T, Int_t id, Int_t cNum, Int_t chNum);
  Double_t ANLEnergy2(Int_t T, Int_t id, Int_t cNum, Int_t chNum);
  Float_t ANLEnergy(Int_t T, Int_t id, Int_t cNum, Int_t chNum);
  Float_t ANLRiseIntegral(Int_t T, Int_t id, Int_t cNum, Int_t chNum);
  Int_t LEDLevel(Int_t index, Int_t thresh, Int_t cNum, Int_t chNum);
  Int_t LED(Int_t index, Int_t thresh, Int_t cNum, Int_t chNum);
  Double_t ComputeChiSquare(Int_t startIndex, Int_t nPoints, 
			    Double_t tau, Double_t Nstart, 
			    Double_t baseline, Int_t cNum, Int_t chNum);
  Int_t LSFitExponential(Int_t startIndex, Int_t nPoints, 
			 Double_t Nstart, Double_t tau, 
			 Double_t &NstartFit, Double_t &baselineFit,
			 Double_t &NstartFitError, 
			 Double_t &baselineFitError, Int_t cNum, Int_t chNum);  
  Double_t LSFitLinear(Int_t startIndex, Int_t nPoints,
		       Double_t &slopeFit, Double_t &offsetFit);
  Float_t Baseline2Flatten(Int_t id, Float_t tau, 
			   Int_t startIndexBase, Int_t lengthBase,
			   Int_t startIndexTop, Int_t lengthTop,
			   Float_t smallStep, Int_t cNum, Int_t chNum);
  Int_t FlattenBaselineToZero(Float_t baseStart, Float_t ampStart, 
			      Float_t itau, Int_t startIndex, 
			      Int_t nPoints, Float_t &baseOut,
			      Float_t &ampOut, Int_t cNum, Int_t chNum);
  void RadfordBaselineFill(Int_t id, Int_t cNum, Int_t chNum);
  void RadfordBaselineFit(const char *file);
  Float_t RadfordEnergy(Int_t id, Int_t cNum, Int_t chNum);
  void PoleZeroTrace(Int_t id, Int_t cNum, Int_t chNum);
  Float_t FPGAFilter(Int_t id, Int_t cNum, Int_t chNum);
  Float_t SimpleEnergy(Int_t cNum, Int_t chNum);
  linfit LinReg(Double_t *trace, Int_t startReg, Int_t stopReg);
  peak GregorichTrapFilter(Int_t id, Int_t cNum, Int_t chNum);
  Float_t CalcT0(Int_t startSample, Int_t cNum, Int_t chNum);
  Int_t Look4Pileup(Int_t cNum, Int_t chNum);
  Int_t TimeSeg(Int_t len, Int_t cNum, Int_t chNum);

 private:
  ClassDef(GRETINAWF, 1);
};

#endif
