#ifndef __GRETINA_H
#define __GRETINA_H

#include "Riostream.h"
#include "Rtypes.h"

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TClass.h"

#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include "SortingStructures.h"
#include "Defines.h"
#include "colors.h"

using std::vector;
using std::cout;  using std::endl;
using std::swap;

class GRETINA;

/******** Structures to store raw data ********/

struct globalHeader {
  Int_t type;
  Int_t length;
  long long timestamp;
};

/*--------------------------------------------------------*/
/*---       GRETINA MODE 3 - RAW DATA STRUCTURES       ---*/
/*--------------------------------------------------------*/

struct mode3DataPacket {
  unsigned short aahdr[2];
  unsigned short hdr[14];
  /* with traces... */
  unsigned short waveform[MAX_TRACE_LENGTH];
};

class g3Waveform : public TObject {

 public:
  vector<Short_t> raw;
  vector<Float_t> pz;
  
 public:
  void Clear();
  Int_t LEDLevel(Int_t index, Int_t thresh, Int_t filterK);
  Int_t LED(Int_t index, Int_t thresh, Int_t filterK);
  Float_t CFD(Int_t startSample);
  Float_t BL(Int_t start, Int_t end);
  Int_t Look4Pileup();
  void calcPZ(Float_t base, Float_t tau);
  Float_t simpleE();
  Float_t riseTime(Float_t fLow, Float_t fHigh);
  Double_t computeChiSquare(Int_t startIndex, Int_t nPoints,
			    Double_t tau, Double_t nStart, Double_t baseline);
  Double_t LSFitExpo(Int_t startIndex, Int_t nPoints,
		     Double_t nStart, Double_t tau,
		     Double_t &nStartFit, Double_t &baselineFit,
		     Double_t &nStartFitError,
		     Double_t &baselineFitError);
  Double_t LSFitLinear(Int_t startIndex, Int_t nPoints,
		       Double_t &slopeFit, Double_t &offsetFit);
  Float_t baseline2Flatten(Float_t tau, Int_t startIndexBase,
			   Int_t lengthBase, Int_t startIndexTop,
			   Int_t lengthTop, Float_t smallStep);
  
  ClassDef(g3Waveform, 1);
};

class g3ChannelEvent : public TObject {

 public:
  UShort_t hdr0, hdr1, hdr7;
  Int_t ID;
  Int_t segNum;
  Float_t eRaw, eCal;
  /* These are for event history, available in firmware 2.00_006f and beyond. */
  Float_t eCalPO, prevE1, prevE2;
  UShort_t deltaT1, deltaT2;
  UShort_t PZrollover;

  Float_t eCalc;
  Float_t calcTime;
  long long int timestamp;
  long long int CFDtimestamp;
  Float_t baseline;
  Float_t riseTime;
  g3Waveform wf;

 public:
  Int_t module() { return (hdr1 >> 4); }
  Int_t boardID() { return (hdr0 >> 11); }
  Int_t chanID() { return (hdr1 & 0xf); }
  Int_t tracelength() { return (((hdr0 & 0x7ff)*2) - 14); }
  Int_t chnNum() { return (ID%40); }
  Bool_t sign() { return (hdr7 & 0x0100); }
  Int_t pileUp() { return ((hdr7 & 0x8000) >> 15); }    

 public:

  ClassDef(g3ChannelEvent, 3);
};

class g3CrystalEvent : public TObject {
 
 public:
  Int_t module;
  Int_t bankNum;
  Int_t runNumber;
  Int_t traceLength;
  vector<g3ChannelEvent> chn;

 public:
  void OrderChannels();
  long long int LEDLow();
  long long int LEDHigh();
  long long int LEDRange();

 private:
  ClassDef(g3CrystalEvent, 2);
};

class g3OUT : public TObject {

 public:
  vector <g3CrystalEvent> xtals;

 public:
  g3OUT() { ; }
  ~g3OUT() { ; }
  void Reset();
  UInt_t bankMult();

 private:
  ClassDef(g3OUT, 1);
};

struct mode3HistoryPacket {
  unsigned short aahdr[2];
  unsigned short hdr[8];
  unsigned short data[MAX_TRACE_LENGTH];
};

struct historyEvent{
  Float_t energy;
  long long int TS;
  Int_t BLpreSum;
  // 2016-07-23 CMC added module to gH to differentiate between digitizers
  Int_t module;
};

class g3HistoryEvent : public TObject {
  
 public:
  std::vector<historyEvent> past;
  
 public:
  g3HistoryEvent() { ; }
  ~g3HistoryEvent() { ; }
  UInt_t lengthOfPast() { return past.size(); }

 private: 
  ClassDef(g3HistoryEvent, 1);
};

class Bank88 : public TObject {
 
 public:
  vector<g3ChannelEvent> chn;
  long long int timestamp;
  long long int wfCFD;

 public:
  Bank88() { ; }
  ~Bank88() { ; }
 
  ClassDef(Bank88, 2);
};

class GRETINA : public TObject {

 public:
  g3CrystalEvent g3X; g3ChannelEvent g3ch;
  vector<g3ChannelEvent> g3Temp;
  historyEvent gH; 

  unsigned char gBuf[32*32*1024];

  g3HistoryEvent g3H;
  g3OUT g3out;
  
  Bank88 b88;

 public:
  Float_t wfMinCrossTime, wfMaxCrossTime;
  
  UShort_t scalerBuf[8192];

 public:
  GRETINA() { ; }
  ~GRETINA() { ; }
  void Initialize();
  void Reset();
  Int_t getMode3(FILE *inf, Int_t evtLength, counterVariables *cnt,
		 controlVariables *ctrl);
  Int_t getMode3History(FILE *inf, Int_t evtLength, long long int hTS, counterVariables *cnt);
  Int_t getBank88(FILE *inf, Int_t evtLength, counterVariables *cnt);

  void analyzeMode3(controlVariables *ctrl);
  void calibrateMode3(g3ChannelEvent *g3);

 public:
  ClassDef(GRETINA, 1);

};

#endif
