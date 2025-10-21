/**********************************************************************/
/* File: SortingStructures.C                                          */
/* Description: Functions for GRETINA/Aux control variables and       */
/*              counting variables in analysis                        */
/* Author: H. Crawford                                                */
/* Date: January 2013                                                 */
/**********************************************************************/

#include <stdlib.h>

#include "SortingStructures.h"

ClassImp(controlVariables);
ClassImp(counterVariables);

/****************************************************/
/* Control variable functions...                    */
/****************************************************/

controlVariables::controlVariables() { 
  doTRACK = 0;
  withWAVE = 0;
  withSEG = 1;
  withHISTOS = 1;
  withTREE = 1;
  dopplerSimple = 0;
  
  specifyCalibration = 0;
  calibrationFile = "";
  
  startRun = 0;
  fileType = "";
  directory = "";
  outputSuffix = "";
  outputON = 0;
  outputName = 0;
  outfileName = "";
  outputFileName = "";
  compressedFile = 0;
  compressedFileB = 0;
  noHFC = 0;
  suppressTS = 0;
  pgh = 0;
  noEB = 0;
  calibration = 0;
  superPulse = 0;
  xtalkAnalysis = 0;

  INLcorrection = 0;
  INLCConly = 0;
  digMapFileName = "";

  gateTree = 0;

  mode2Old = 0;

  analyze2AND3 = 0;
  fileName = "";
  
  WITH_TRACETREE = 0;
  CHECK_PILEUP = 0;
  FALLON_TIME = 0;
  LEDCROSSING = 0;
  FIT_BASELINE = 0;
  PZ_BASELINE = 0;
  PZ_ZERO_BASELINE = 0; 
  RADFORD_BASELINE = 0;
  RADFORD_ENERGY = 0;
  GREGORICH_ENERGY = 0;
  POLEZERO_TRACE = 0;
  FPGA_ENERGY = 0;
  BASIC_ENERGY = 0;
  INL_CORRECT = 0;
  
  //#ifdef WITH_S800
  s800File = 0;
  s800ControlFile = "";
  
  S800_DIAG=0;
  E1_RAW=0; E1_CAL=0; E2_RAW=0; E2_CAL=0; E3_RAW=0; E3_CAL=0;
  IC_RAW=0; IC_CAL=0; IC_SUMS=0;
  CRDC1_RAW_PADS=0; CRDC1_RAW_CALC=0; CRDC1_CALC=0;
  CRDC2_RAW_PADS=0; CRDC2_RAW_CALC=0; CRDC2_CALC=0;
  FP_TRACK_RAW=0; FP_TRACK_COR=0;
  HODO_RAW=0; HODO_CAL=0;
  TARGET_PPAC_RAW=0; TARGET_PPAC_CALC=0;
  TARGET_PIN1_RAW=0; TARGET_PIN1_CAL=0;
  TARGET_PIN2_RAW=0; TARGET_PIN2_CAL=0;
  TARGET_TOTAL=0;
  IMAGE_CALC=0;
  IMAGE_TPPAC1_RAW=0; IMAGE_TPPAC1_CALC=0;
  IMAGE_TPPAC2_RAW=0; IMAGE_TPPAC2_CALC=0;
  IMAGE_TRACK=0;
  IMAGE_PPAC1_RAW=0; IMAGE_PPAC1_CALC=0;
  IMAGE_PPAC2_RAW=0; IMAGE_PPAC2_CALC=0;
  OBJECT_PIN_RAW=0; OBJECT_PIN_CAL=0;
  TRIGGER=0; S800_TIMESTAMP=0;
  TOF=0;
  //#endif
}

/****************************************************/

void controlVariables::Initialize() {  
  doTRACK = 0;
  withWAVE = 0;
  withSEG = 1;
  withHISTOS = 0;
  withTREE = 1;
  dopplerSimple = 0;
  
  specifyCalibration = 0;
  calibrationFile = "";
  
  startRun = 0;
  fileType = "";
  directory = "";
  outputON = 0;
  outputName = 0;
  outputFileName = "";
  outputSuffix = "";
  outfileName = "";
  compressedFile = 0;
  compressedFileB = 0;
  noHFC = 0;
  suppressTS = 0;
  pgh = 0;
  noEB = 0;
  calibration = 0;

  INLcorrection = 0;
  digMapFileName = "";

  gateTree = 0;

  mode2Old = 0;
  
  analyze2AND3 = 0;
  fileName2 = "";

  WITH_TRACETREE = 0;
  CHECK_PILEUP = 0;
  FALLON_TIME = 0;
  LEDCROSSING = 0;
  FIT_BASELINE = 0;
  PZ_BASELINE = 0;
  PZ_ZERO_BASELINE = 0; 
  RADFORD_BASELINE = 0;
  RADFORD_ENERGY = 0;
  GREGORICH_ENERGY = 0;
  POLEZERO_TRACE = 0;
  FPGA_ENERGY = 0;
  BASIC_ENERGY = 0;
  INL_CORRECT = 0;
  
  //#ifdef WITH_S800
  s800File = 0;
  s800ControlFile = "";
  S800_DIAG=0;
  E1_RAW=0; E1_CAL=0; E2_RAW=0; E2_CAL=0; E3_RAW=0; E3_CAL=0;
  IC_RAW=0; IC_CAL=0; IC_SUMS=0;
  CRDC1_RAW_PADS=0; CRDC1_RAW_CALC=0; CRDC1_CALC=0;
  CRDC2_RAW_PADS=0; CRDC2_RAW_CALC=0; CRDC2_CALC=0;
  FP_TRACK_RAW=0; FP_TRACK_COR=0;
  HODO_RAW=0; HODO_CAL=0;
  TARGET_PPAC_RAW=0; TARGET_PPAC_CALC=0;
  TARGET_PIN1_RAW=0; TARGET_PIN1_CAL=0;
  TARGET_PIN2_RAW=0; TARGET_PIN2_CAL=0;
  TARGET_TOTAL=0;
  IMAGE_CALC=0;
  IMAGE_TPPAC1_RAW=0; IMAGE_TPPAC1_CALC=0;
  IMAGE_TPPAC2_RAW=0; IMAGE_TPPAC2_CALC=0;
  IMAGE_TRACK=0;
  IMAGE_PPAC1_RAW=0; IMAGE_PPAC1_CALC=0;
  IMAGE_PPAC2_RAW=0; IMAGE_PPAC2_CALC=0;
  OBJECT_PIN_RAW=0; OBJECT_PIN_CAL=0;
  TRIGGER=0; S800_TIMESTAMP=0;
  TOF=0;
  //#endif
  
}

/****************************************************/

Int_t controlVariables::InterpretCommandLine(int argc, char *argv[]) {
  Int_t i=1;
  while (i < argc && startRun == 0) {
    if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "-gr") == 0 || 
	strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-f2") == 0 ||
	strcmp(argv[i], "-f1") == 0) {
      fileType = argv[i];
      fileType.ReplaceAll("-","");
      if (fileType == "f" || fileType == "f1") {
	fileName = argv[i+1];
	i++;
      }
      if (fileType == "f2") {
	fileName2 = argv[i+1];
	i++;
      }
      i++;
    }
    else if (strcmp(argv[i], "-wf") == 0) {
      withWAVE = 1;
      WITH_TRACETREE = 1;
      cout << "Waveform analysis enabled. " << endl;
      i++;
    }
    else if (strcmp(argv[i], "-rootName") == 0) {
      i++;
      outfileName = argv[i]; i++;
    }
    else if (strcmp(argv[i], "-zip") == 0) {
      compressedFile = 1;
      i++;
    }
    else if (strcmp(argv[i], "-bzip") == 0) {
      compressedFileB = 1;
      i++;
    }
    else if (strcmp(argv[i], "-d") == 0) {
      directory = argv[i+1];
      i += 2;
    }
    else if (strcmp(argv[i], "-run") == 0) {
      startRun = i+1;
    }
    else if (strcmp(argv[i], "-noHFC") == 0) {
      noHFC = 1;
      i++;
    }
    else if (strcmp(argv[i], "-noEB") == 0) {
      noEB = 1;
      cout << "Event building turned off." << endl;
      i++;
    } else {
      cout << "Error -- unrecognized input flag: " << argv[i] << endl;
      return -1;
    }
  }
  return 1;
}

/****************************************************/

Int_t controlVariables::ReportRunFlags() {
  cout << "  Analysis conditions: " << endl;
  cout << "     Expecting ";
  if (compressedFile) {
    cout << " compressed ";
  } 
  if (fileType == "g" && !pgh) {
    cout << "GRETINA + auxiliary file WITH global headers -- Global.dat";
  } else if (fileType == "gr" && !pgh) {
    cout << "GRETINA mode3 file WITH global headers -- GlobalRaw.dat";
  }
  if (compressedFile) {
    cout << ".gz " << endl;
  }  else {
    cout << endl;
  }
  if (withTREE) {
    cout << "     Will write out a ROOT tree, but no histograms. " << endl;
  }
  if (noHFC) {
    cout << "     Will NOT use Dirk's GEB_HFC resorter code. " << endl;
  }
  if (!noHFC) {
    cout << "     Will use Dirk's GEB_HFC resorter code. " << endl;
  }
  if (withWAVE) {
    cout << "     Will do some waveform analysis on GRETINA waveforms, as specified with flags in Unpack.h. " << endl;
  }
  return(1);
}

/****************************************************/
/* Counter variable functions...                    */
/****************************************************/

void counterVariables::Initialize() {
  event = 0x0000;

  TSFirst = 0; TSLast = 0;

  for (Int_t i=0; i<100; i++) { headerType[i] = 0; }
  
  TSerrors = 0;
  
  bytes_read = 0; bytes_read_since_last_time = 0;
  MBread = 0;
  
  eoBuffer = 0; eofInBuffer = 0; eofPosInBuffer = 0;
  mode3i = 0; old3Bytes = 0;
  b88i = 0;

  treeWrites = 0;
  
  badEvent = 0; badSegment = 0; badCC1 = 0; badCC2 = 0;
  for (Int_t i=0; i<(MAXCHANNELS); i++) {
    goodTraceE_PU[i] = 0; goodTraceE[i] = 0;
    badFPGA_zero_PU[i] = 0; badFPGA_zero[i] = 0;
    badFPGA_neg_PU[i] = 0; badFPGA_neg[i] = 0;
  }
  crystalBuildEvent = 0; totalCrystalEvent = 0;
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    crystalBuildEventXtal[i] = 0; totalCrystalEventXtal[i] = 0;
  }
  tossed4Time = 0;
  
  for (Int_t i=0; i<4*MAXCRYSTALS; i++) { /* 4 boards/crystal */
    lastBdTS[i] = 0;
  }
  
  nGammasRead = 0; nGammasThisHeader = 0;

  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<41; j++) {
      nMode3Skipped[i][j] = 0;
      nMode3SkippedAtEnd[i][j] = 0;
    }
    nMode2NoTraces[i] = 0;
    nMode3NoMode2[i] = 0;
  }
  
}

/****************************************************/

void counterVariables::ResetRunCounters() {
  event = 0x0000;
  for (Int_t i=0; i<100; i++) { headerType[i] = 0; }
  
  TSFirst = 0; TSLast = 0;

  TSerrors = 0;
  
  bytes_read = 0; bytes_read_since_last_time = 0;
  MBread = 0;

  treeWrites = 0;
  
  eoBuffer = 0; eofInBuffer = 0; eofPosInBuffer = 0;
  mode3i = 0; old3Bytes = 0;
  b88i = 0;

  badEvent = 0; badSegment = 0; badCC1 = 0; badCC2 = 0;
  for (Int_t i=0; i<(MAXCHANNELS); i++) {
    goodTraceE_PU[i] = 0; goodTraceE[i] = 0;
    badFPGA_zero_PU[i] = 0; badFPGA_zero[i] = 0;
    badFPGA_neg_PU[i] = 0; badFPGA_neg[i] = 0;
  }
  crystalBuildEvent = 0; totalCrystalEvent = 0;
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    crystalBuildEventXtal[i] = 0; totalCrystalEventXtal[i] = 0;
  }
  tossed4Time = 0;
  
  for (Int_t i=0; i<(MAXCRYSTALS*4); i++) {
    lastBdTS[i] = 0;
  }
  
  nGammasRead = 0; nGammasThisHeader = 0;

  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<41; j++) {
      nMode3Skipped[i][j] = 0;
      nMode3SkippedAtEnd[i][j] = 0;
    }
    nMode2NoTraces[i] = 0;
    nMode3NoMode2[i] = 0;
  }
  
}

/****************************************************/

void counterVariables::PrintRunStatistics(Int_t pgh, Int_t withWAVE, Int_t superPulse, 
					  Int_t sort2and3) {
  printf("-----------------------------------------------------------\n");
  if (!pgh) {
    printf(" Run Statistics:\n");
    printf(" GRETINA-related data...\n");
    if (headerType[TRACK] > 0) 
      printf("  Mode1 GRETINA headers:      %d\n", headerType[TRACK]);
    if (headerType[DECOMP] > 0) 
      printf("  Mode2 GRETINA headers:      %d\n", headerType[DECOMP]);
    if (headerType[RAW] > 0)
      printf("  Mode3 GRETINA headers:      %d\n", headerType[RAW]);
    if (headerType[G4SIM] > 0) 
      printf("  Simulation GRETINA headers: %d\n", headerType[G4SIM]);
    
    if (headerType[GRETSCALER] > 0)
      printf("  Scaler GRETINA headers:     %d\n", headerType[GRETSCALER]);
    if (headerType[BANK88] > 0)
      printf("  Bank88 GRETINA headers:     %d\n", headerType[BANK88]);
   
  }
  printf("\n TS errors:   %d\n", TSerrors);
  printf("\n\n Run time (from TS): %0.3f seconds\n", (TSLast - TSFirst)*1e-8);

  cout << "--------------------------------------------------" << endl;

}

void counterVariables::Increment(Int_t bytes) {
  bytes_read += bytes;
  bytes_read_since_last_time += bytes;
}

void counterVariables::setEventBit(Int_t bit) {
  event |= (0x1 << bit);
}

Int_t counterVariables::getEventBit(Int_t bit) {
  return ((event >> bit) & 0x1);
}
