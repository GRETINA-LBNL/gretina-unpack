#ifndef __SORTINGSTRUCTURES_H
#define __SORTINGSTRUCTURES_H

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TString.h"
#include <vector>
#include <stdint.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include "Defines.h"

using std::cout; using std::endl;
using std::cin; using std::cerr;

/* Analysis control variables...   */

class controlVariables : public TObject {
  
 public:
  Int_t doTRACK;
  Int_t withWAVE;
  Int_t withSEG;

  Int_t withHISTOS;
  Int_t withTREE;  

  Int_t dopplerSimple;

  Int_t startRun;
  TString fileType;
  TString fileName;
  TString outputSuffix;
  TString outfileName;
  TString directory;
  Int_t outputON;
  Int_t outputName;
  TString outputFileName;
  Int_t compressedFile;
  Int_t compressedFileB;
  Int_t noHFC;
  Int_t suppressTS;
  Int_t pgh;
  Int_t noEB;
  Int_t calibration;
  Int_t specifyCalibration;
  TString calibrationFile;
  Int_t xtalkAnalysis;
  Int_t xtalkID; 
  Int_t xtLowE, xtHighE;
  
  /* Superpulse information */
  Bool_t superPulse;
  Float_t spLowE, spHighE;
  TString spXtalkFile;
  char* spCalibDirectory;

  /* INL correction information */
  Bool_t INLcorrection;
  Bool_t INLCConly;
  TString digMapFileName;

  Bool_t gateTree;

  Bool_t mode2Old;

  /* GRETINA waveform analysis flags. */
  Bool_t WITH_TRACETREE;
  Bool_t CHECK_PILEUP;
  Bool_t FALLON_TIME;
  Bool_t LEDCROSSING;
  Bool_t FIT_BASELINE;
  Bool_t PZ_BASELINE;
  Bool_t PZ_ZERO_BASELINE;
  Bool_t RADFORD_BASELINE;
  Bool_t RADFORD_ENERGY;
  Bool_t GREGORICH_ENERGY;
  Bool_t POLEZERO_TRACE;
  Bool_t FPGA_ENERGY;
  Bool_t BASIC_ENERGY;
  Bool_t INL_CORRECT;

  /* Analysis of Mode2 and Mode3 together */
  Bool_t analyze2AND3;
  TString fileName2;

  //#ifdef S800
  Int_t s800File;
  TString s800ControlFile;
  TString s800VariableFile;
  
  Bool_t S800_DIAG;
  Bool_t E1_RAW, E1_CAL, E2_RAW, E2_CAL, E3_RAW, E3_CAL;
  Bool_t IC_RAW, IC_CAL, IC_SUMS;
  Bool_t CRDC1_RAW_PADS, CRDC1_RAW_CALC, CRDC1_CALC;
  Bool_t CRDC2_RAW_PADS, CRDC2_RAW_CALC, CRDC2_CALC;
  Bool_t FP_TRACK_RAW, FP_TRACK_COR;
  Bool_t HODO_RAW, HODO_CAL;
  Bool_t TARGET_PPAC_RAW, TARGET_PPAC_CALC;
  Bool_t TARGET_PIN1_RAW, TARGET_PIN1_CAL;
  Bool_t TARGET_PIN2_RAW, TARGET_PIN2_CAL;
  Bool_t TARGET_TOTAL;
  Bool_t IMAGE_CALC;
  Bool_t IMAGE_TPPAC1_RAW, IMAGE_TPPAC1_CALC;
  Bool_t IMAGE_TPPAC2_RAW, IMAGE_TPPAC2_CALC;
  Bool_t IMAGE_TRACK;
  Bool_t IMAGE_PPAC1_RAW, IMAGE_PPAC1_CALC;
  Bool_t IMAGE_PPAC2_RAW, IMAGE_PPAC2_CALC;
  Bool_t OBJECT_PIN_RAW, OBJECT_PIN_CAL;
  Bool_t TRIGGER, S800_TIMESTAMP;
  Bool_t TOF;
  //#endif

 public:
  controlVariables();
  ~controlVariables() { ; }
  void Initialize();
  Int_t InterpretCommandLine(int argc, char *argv[]);  
  Int_t ReportRunFlags();
  void SetS800Controls(TString fileName);
  
  private:  
    ClassDef(controlVariables, 1);
};

/* Analysis counters... */

class counterVariables : public TObject {

 public:
  Int_t event;

  /* Run number */
  Int_t runNum;

  /* TS for run length info */
  long long int TSFirst, TSLast;

  /* GEB header counters */
  Int_t headerType[100];

  Int_t TSerrors;

  /* Overall sorting variables */
  long long int bytes_read;
  long long int bytes_read_since_last_time;
  UInt_t MBread;

  long long int treeWrites;

  /* Mode3 analysis variables */
  UInt_t eoBuffer;
  Int_t eofInBuffer;
  UInt_t eofPosInBuffer;
  UInt_t mode3i;
  UInt_t b88i;

  long long int old3Bytes;

  /* Waveform analysis counters */
  Int_t badEvent;
  Int_t badSegment;
  Int_t badCC1;
  Int_t badCC2;
  Int_t goodTraceE_PU[MAXCHANNELS];
  Int_t goodTraceE[MAXCHANNELS];
  Int_t badFPGA_zero_PU[MAXCHANNELS];
  Int_t badFPGA_zero[MAXCHANNELS];
  Int_t badFPGA_neg_PU[MAXCHANNELS];
  Int_t badFPGA_neg[MAXCHANNELS];
  Int_t crystalBuildEvent;
  Int_t crystalBuildEventXtal[MAXCRYSTALS];
  Int_t totalCrystalEvent;
  Int_t totalCrystalEventXtal[MAXCRYSTALS];
  
  Int_t tossed4Time;

  long long int lastBdTS[MAXCRYSTALS*4];

  /* Mode2 analysis variables */
  Int_t nGammasThisHeader;
  Int_t nGammasRead;

  /* Mode2+3 counters */
  Int_t nMode3Skipped[MAXCRYSTALS][41];
  Int_t nMode3SkippedAtEnd[MAXCRYSTALS][41];
  Int_t nMode2NoTraces[MAXCRYSTALS];
  Int_t nMode3NoMode2[MAXCRYSTALS];
  
 public:
  counterVariables() { ; }
  ~counterVariables() { ; }
  void Initialize();
  void ResetRunCounters();
  void PrintRunStatistics(Int_t pgh, Int_t withWAVE, Int_t superPulse, Int_t sort2and3);
  void Increment(Int_t bytes);
  void setEventBit(Int_t bit);
  Int_t getEventBit(Int_t bit);

 private:  
  ClassDef(counterVariables, 1);
}; 

#endif
