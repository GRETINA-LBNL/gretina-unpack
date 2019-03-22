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
	strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "-m2") == 0 || 
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
    else if (strcmp(argv[i], "-track") == 0) {
      doTRACK = 1;
      cout << "Tracking enabled. " << endl;
      i++;
    }
    else if (strcmp(argv[i], "-wf") == 0) {
      withWAVE = 1;
      WITH_TRACETREE = 1;
      cout << "Waveform analysis enabled. " << endl;
      i++;
    }
    else if (strcmp(argv[i], "-noSeg") == 0) {
      withSEG = 0;
      cout << "Segment analysis disabled. " << endl;
      i++;
    }
    else if (strcmp(argv[i], "-outputON") == 0) {
      outputON = 1;
      cout << "Will write S800Physics + GRETINA output according to gate conditions. " << endl;
      i++;
    }
    else if (strcmp(argv[i], "-outputName") == 0) {
      outputName = 1; outputON = 1;
      cout << "Will write output called ";
      i++;
      outputFileName = argv[i]; i++;
      cout << outputFileName << endl;
    }
    else if (strcmp(argv[i], "-outputSuffix") == 0) {
      i++;
      outputSuffix = argv[i];
      cout << "Will include " << outputSuffix.Data()
	   << " on ROOT files." << endl;
    }
    else if (strcmp(argv[i], "-rootName") == 0) {
      i++;
      outfileName = argv[i]; i++;
    }
    else if (strcmp(argv[i], "-analyze2and3") == 0) {
      i++;
      analyze2AND3 = 1;
      withWAVE = 1;
    }
    else if (strcmp(argv[i], "-gateTree") == 0) {
      gateTree = 1;
      i++;
    }
    else if (strcmp(argv[i], "-zip") == 0) {
      compressedFile = 1;
      i++;
    }
    else if (strcmp(argv[i], "-bzip") == 0) {
      compressedFileB = 1;
      i++;
    }
    else if (strcmp(argv[i], "-withHistos") == 0) {
      withHISTOS = 1;
      i++;
    }
    else if (strcmp(argv[i], "-noTree") == 0) {
      withTREE = 0;
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
    else if (strcmp(argv[i], "-dopplerSimple") == 0) {
      dopplerSimple = 1;
      i++;
    }
    else if (strcmp(argv[i], "-suppressTS") == 0) {
      suppressTS = 1;
      i++;
    }
    else if (strcmp(argv[i], "-s800File") == 0) {
      s800File = 1;
      s800ControlFile = argv[i+1];
      i+=2;
    }
    else if (strcmp(argv[i], "-preGH2") == 0) {
      pgh = 2;
      noHFC = 1;
      i++;
    }
    else if (strcmp(argv[i], "-preGH3") == 0) {
      pgh = 3;
      noHFC = 1;
      i++;
    }
    else if (strcmp(argv[i], "-noEB") == 0) {
      noEB = 1;
      cout << "Event building turned off." << endl;
      i++;
    }
    else if (strcmp(argv[i], "-calibrationRun") == 0) {
      calibration = 1; noEB = 1; withHISTOS = 1; withTREE = 0;
      cout << "Calibration run. Disabled event building too. " << endl;
      i++;
    }
    else if (strcmp(argv[i], "-readCal") == 0) {
      specifyCalibration = 1; calibrationFile = argv[i+1];
      i+=2;
    }
    else if (strcmp(argv[i], "-superPulse") == 0) {
      superPulse = 1; i++;
      withWAVE = 1;
      spLowE = atoi(argv[i]); i++; 
      spHighE = atoi(argv[i]); i++;
      spCalibDirectory = argv[i]; i++;
      spXtalkFile = argv[i]; i++;
      withHISTOS = 0; withTREE = 0;
    }
    else if (strcmp(argv[i], "-xtalkAnalysis") == 0) {
      xtalkAnalysis = 1; i++;
      xtalkID = atoi(argv[i]); i++;
      xtLowE = atoi(argv[i]); i++;
      xtHighE = atoi(argv[i]); i++;
      withTREE = 0;
    } else if (strcmp(argv[i], "-INLcorrect") == 0) {
      INLcorrection = 1; i++;
      INLCConly = atoi(argv[i]); i++;
      digMapFileName = argv[i]; i++; 
      withWAVE = 1;
      WITH_TRACETREE = 1;
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
  } else if (fileType == "g" && pgh == 2) {
    cout << "old GRETINA mode2 only file WITHOUT global headers -- Global.dat";
  } else if (fileType == "g" && pgh == 3) {
    cout << "old GRETINA mode3 only file WITHOUT global headers -- Global.dat";
  } else if (fileType == "gr" && !pgh) {
    cout << "GRETINA mode3 file WITH global headers -- GlobalRaw.dat";
  } else if (fileType == "m") {
    cout << "merged GRETINA mode3 + auxiliary WITH global headers -- Merged.dat";
  } else if (fileType == "m2") {
    cout << "merged GRETINA mode2 + auxiliary WITH global headers -- Merged.Mode2.dat";
  }
  if (compressedFile) {
    cout << ".gz " << endl;
  }  else {
    cout << endl;
  }
  if (withTREE && withHISTOS ) {
    cout << "     Will write out a ROOT tree and defined histograms." << endl;
  }
  if (withTREE && !withHISTOS) {
    cout << "     Will write out a ROOT tree, but no histograms. " << endl;
  }
  if (!withTREE && withHISTOS) {
    cout << "     Will write out defined histograms, but no ROOT tree. " << endl;
  }
  if (!withTREE && !withHISTOS && !superPulse) {
    cout << "     You aren't planning on writing histograms or a tree...are you sure you want to proceed? ('yes' or 'no')" << endl;
    char answer[50];
    cin >> answer;
    if (strcmp(answer, "yes") != 0) {
      cout << "Aborting now. " << endl;
      return(-1);
    }
  }
  if (noHFC && suppressTS) {
    cout << "     Will NOT use Dirk's GEB_HFC resorter code, and will ignore TS order errors. " << endl;
  }
  if (noHFC && !suppressTS) {
    cout << "     Will NOT use Dirk's GEB_HFC resorter code, and will print any observed TS order errors. " << endl;
  }
  if (!noHFC && suppressTS) {
    cout << "     Will use Dirk's GEB_HFC resorter code, and will print any observed TS order errors. " << endl;
  }
  if (!noHFC && !suppressTS) {
    cout << "     Will use Dirk's GEB_HFC resorter code, and will ignore TS order errors. " << endl;
  }
  if (doTRACK) { 
    cout << "     Will try really hard to do tracking, using Torben's code. " << endl;
  }
  if (withWAVE) {
    cout << "     Will do some waveform analysis on GRETINA waveforms, as specified with flags in Unpack.h. " << endl;
  }
  if (outputON) {
    cout << "     Will generate combined output data file with processed aux data where appropriate (i.e. mode9 S800Physics) " << endl;
  }
  return(1);
}

void controlVariables::SetS800Controls(TString fileName) {
  FILE *controlINPUT;
  if ( (controlINPUT = fopen(fileName.Data(), "r")) == NULL ) {
    cerr << "S800 control file " << fileName.Data() << " could not be opened." << endl;
  }
  cout << "Setting S800 controls based on " << fileName.Data() << endl;
  
  char line[300];
  char junk[300];
  char filename[300];
  int value;

  while ( !feof(controlINPUT) ) {
    char* chr = fgets(line, 300, controlINPUT);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "VariableFile", 12) == 0) {
      sscanf(line, "%s %s", junk, filename);
      s800VariableFile = filename;
    }
    if (strncmp(line, "S800_DIAG", 9) == 0) {
      sscanf(line, "%s %d", junk, &value);
      S800_DIAG = value;
    }
    if (strncmp(line, "E1_RAW", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      E1_RAW = value;
    }
    if (strncmp(line, "E1_CAL", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      E1_CAL = value;
    }  
    if (strncmp(line, "E2_RAW", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      E2_RAW = value;
    }
    if (strncmp(line, "E2_CAL", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      E2_CAL = value;
    }    
    if (strncmp(line, "E3_RAW", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      E3_RAW = value;
    }
    if (strncmp(line, "E3_CAL", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      E3_CAL = value;
    }    
    if (strncmp(line, "IC_RAW", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IC_RAW = value;
    }
    if (strncmp(line, "IC_CAL", 6) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IC_CAL = value;
    } 
    if (strncmp(line, "CRDC1_RAW_PADS", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      CRDC1_RAW_PADS = value;
    } 
    if (strncmp(line, "CRDC1_RAW_CALC", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      CRDC1_RAW_CALC = value;
    } 
    if (strncmp(line, "CRDC1_CALC", 10) == 0) {
      sscanf(line, "%s %d", junk, &value);
      CRDC1_CALC = value;
    } 
    if (strncmp(line, "CRDC2_RAW_PADS", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      CRDC2_RAW_PADS = value;
    } 
    if (strncmp(line, "CRDC2_RAW_CALC", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      CRDC2_RAW_CALC = value;
    } 
    if (strncmp(line, "CRDC2_CALC", 10) == 0) {
      sscanf(line, "%s %d", junk, &value);
      CRDC2_CALC = value;
    } 
    if (strncmp(line, "S800_TIMESTAMP", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      S800_TIMESTAMP = value;
    }
    if (strncmp(line, "FP_TRACK_RAW", 12) == 0) {
      sscanf(line, "%s %d", junk, &value);
      FP_TRACK_RAW = value;
    } 
    if (strncmp(line, "FP_TRACK_COR", 12) == 0) {
      sscanf(line, "%s %d", junk, &value);
      FP_TRACK_COR = value;
    } 
    if (strncmp(line, "HODO_RAW", 8) == 0) {
      sscanf(line, "%s %d", junk, &value);
      HODO_RAW = value;
    } 
    if (strncmp(line, "HODO_CAL", 8) == 0) {
      sscanf(line, "%s %d", junk, &value);
      HODO_CAL = value;
    } 
    if (strncmp(line, "TARGET_PPAC_RAW", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_PPAC_RAW = value;
    } 
    if (strncmp(line, "TARGET_PPAC_CALC", 16) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_PPAC_CALC = value;
    } 
    if (strncmp(line, "TARGET_PIN1_RAW", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_PIN1_RAW = value;
    } 
    if (strncmp(line, "TARGET_PIN1_CAL", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_PIN1_CAL = value;
    } 
    if (strncmp(line, "TARGET_PIN2_RAW", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_PIN2_RAW = value;
    } 
    if (strncmp(line, "TARGET_PIN2_CAL", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_PIN2_CAL = value;
    } 
    if (strncmp(line, "TARGET_TOTAL", 12) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TARGET_TOTAL = value;
    } 
    if (strncmp(line, "IMAGE_CALC", 10) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_CALC = value;
    } 
    if (strncmp(line, "IMAGE_TPPAC1_RAW", 16) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_TPPAC1_RAW = value;
    } 
    if (strncmp(line, "IMAGE_TPPAC1_CALC", 17) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_TPPAC1_CALC = value;
    } 
    if (strncmp(line, "IMAGE_TPPAC2_RAW", 16) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_TPPAC2_RAW = value;
    } 
    if (strncmp(line, "IMAGE_TPPAC2_CALC", 17) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_TPPAC2_CALC = value;
    } 
    if (strncmp(line, "IMAGE_TRACK", 11) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_TRACK = value;
    } 
    if (strncmp(line, "IMAGE_PPAC1_RAW", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_PPAC1_RAW = value;
    } 
    if (strncmp(line, "IMAGE_PPAC1_CALC", 16) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_PPAC1_CALC = value;
    } 
    if (strncmp(line, "IMAGE_PPAC2_RAW", 15) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_PPAC2_RAW = value;
    } 
    if (strncmp(line, "IMAGE_PPAC2_CALC", 16) == 0) {
      sscanf(line, "%s %d", junk, &value);
      IMAGE_PPAC2_CALC = value;
    } 
    if (strncmp(line, "OBJECT_PIN_RAW", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      OBJECT_PIN_RAW = value;
    } 
    if (strncmp(line, "OBJECT_PIN_CAL", 14) == 0) {
      sscanf(line, "%s %d", junk, &value);
      OBJECT_PIN_CAL = value;
    } 
    if (strncmp(line, "TRIGGER", 7) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TRIGGER = value;
    } 
    if (strncmp(line, "TOF", 3) == 0) {
      sscanf(line, "%s %d", junk, &value);
      TOF = value;
    } 
  }
  fclose(controlINPUT);
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
    
    printf("\n Auxiliary detectors...\n");
    if (headerType[BGS] > 0) 
      printf("  BGS headers:           %d\n", headerType[BGS]);
    if (headerType[CHICO] > 0) 
      printf("  CHICO headers:         %d\n", headerType[CHICO]);
    if (headerType[DFMA] > 0)
      printf("  DFMA headers:          %d\n", headerType[DFMA]);
    if (headerType[PWALL] > 0) 
      printf("  PhosWall headers:      %d\n", headerType[PWALL]);
    if (headerType[PWALLAUX] > 0) 
      printf("  PhosWall Aux. headers: %d\n", headerType[PWALL]);
    if (headerType[S800] > 0)
      printf("  S800 headers:          %d\n", headerType[S800]);
    if (headerType[S800PHYSICS] > 0)
      printf("  S800Physics headers:   %d\n", headerType[S800PHYSICS]);
    if (headerType[S800AUX] > 0 || headerType[S800AUX_TS] > 0)
      printf("  S800 Aux. headers:     %d\n", (headerType[S800AUX]+headerType[S800AUX_TS]));
    if (headerType[LENDA] > 0) 
      printf("  LENDA headers:   %d\n", headerType[LENDA]);
    if (headerType[GODDESS] > 0) 
      printf("  GODDESS headers: %d\n", headerType[GODDESS]);
  }
  printf("\n TS errors:   %d\n", TSerrors);
  printf("\n\n Run time (from TS): %0.3f seconds\n", (TSLast - TSFirst)*1e-8);

  if (sort2and3) {
    cout << endl << "Mode2 + 3 Analysis concerns..." << endl;
    cout << endl;
    cout << "  Mode3 events skipped: " << endl;
    for (Int_t j=0; j<41; j++) {
      for (Int_t i=1; i<MAXCRYSTALS; i++) {
	nMode3Skipped[0][j]+=nMode3Skipped[i][j];
      }
      if (nMode3Skipped[0][j] > 0) {
	if (j<10) {
	  cout << "      " << j << " channels --> " << nMode3Skipped[0][j] << endl;
	} else {
	  cout << "     " << j << " channels --> " << nMode3Skipped[0][j] << endl;
	}
      }
    }
    cout << endl;
    cout << "  Mode3 events skipped at end of file: " << endl;
    for (Int_t j=0; j<41; j++) {
      for (Int_t i=1; i<MAXCRYSTALS; i++) {
	nMode3SkippedAtEnd[0][j]+=nMode3Skipped[i][j];
      }
      if (nMode3SkippedAtEnd[0][j] > 0) {
	if (j<10) {
	  cout << "      " << j << " channels --> " << nMode3SkippedAtEnd[0][j] << endl;
	} else {
	  cout << "     " << j << " channels --> " << nMode3SkippedAtEnd[0][j] << endl;
	}
      }
    }
    cout << endl;
    
    cout <<  "Mode2 without waveforms: ";
    for (Int_t i=1; i<MAXCRYSTALS; i++) {
      nMode2NoTraces[0]+=nMode2NoTraces[i];
    }
    cout << nMode2NoTraces[0] << endl << endl;
    
    cout <<  "Mode3 without Mode2: ";
    for (Int_t i=1; i<MAXCRYSTALS; i++) {
      nMode3NoMode2[0]+=nMode3NoMode2[i];
    }
    cout << nMode3NoMode2[0] << endl << endl; 
  }    

  /* This is a little antiquated... from the days of missing FPGA energies, 
     etc.  I'll keep it around for now, to remember how much the early data sucked. */
  
  if (withWAVE && !superPulse) {
    cout << "Segments tossed from sums due to bad LED timing: " << tossed4Time << "\n" << endl;
    
    cout << "Good crystal build events: " << crystalBuildEvent << "/" << totalCrystalEvent 
	 << " = " << (float)(100*crystalBuildEvent)/totalCrystalEvent << "%" << endl;
    cout << "  # with a bad channel: " << badEvent << " (" 
	 << (float)(100*badEvent)/crystalBuildEvent << "%)" << endl;
    cout << "  # with bad segment:   " << badSegment << " (" 
	 << (float)(100*badSegment)/crystalBuildEvent << "%)" <<  endl;
    cout << "  # with bad cc1:       " << badCC1 << " (" 
	 << (float)(100*badCC1)/crystalBuildEvent << "%)" << endl;
    cout << "  # with bad cc2:       " << badCC2 << " (" 
	 << (float)(100*badCC2)/crystalBuildEvent << "%)" << endl; 
    cout << endl;
    
    Int_t tempSum1 = 0; Int_t tempSum2 = 0;
    Int_t tempSum1_pu = 0; Int_t tempSum2_pu = 0;
    Int_t tempSum2_neg = 0; Int_t tempSum2_neg_pu = 0;
    
    for (Int_t i=0; i<(MAXCHANNELS); i++) {
      tempSum1 += goodTraceE[i]; tempSum1_pu += goodTraceE_PU[i];
    }
    cout << "  # good channel trace E:  " << tempSum1+tempSum1_pu 
	 << " (" << (float)tempSum1_pu/(tempSum1_pu+tempSum1) << "% PU flag)" << endl;
    
    tempSum2 = 0;  tempSum2_pu = 0;
    tempSum2_neg = 0;   tempSum2_neg_pu = 0;
    for (Int_t i=0; i<(MAXCHANNELS); i++) {
      tempSum2 += (badFPGA_zero[i]); tempSum2_pu += (badFPGA_zero_PU[i]);
      tempSum2_neg += (badFPGA_neg[i]); tempSum2_neg_pu += (badFPGA_neg_PU[i]);
    }
    cout << "  # bad FPGA energies:     " << tempSum2+tempSum2_pu+tempSum2_neg+tempSum2_neg_pu 
	 << " (" 
	 << (float)100*(tempSum2+tempSum2_pu+tempSum2_neg+tempSum2_neg_pu)/(tempSum1+tempSum1_pu) 
	 << "%)" << endl;
    cout << "            zero FPGA:     " << tempSum2+tempSum2_pu << " (" 
	 << (float)100*(tempSum2+tempSum2_pu)/(tempSum1+tempSum1_pu) << "%)" << endl;
    cout << "        negative FPGA:     " << tempSum2_neg+tempSum2_neg_pu << " (" 
	 << (float)100*(tempSum2_neg+tempSum2_neg_pu)/(tempSum1+tempSum1_pu) << "%)" << endl;
    
    tempSum1 = 0;  tempSum1_pu = 0;
    tempSum2 = 0;  tempSum2_pu = 0;
    tempSum2_neg = 0;   tempSum2_neg_pu = 0;
    for (Int_t i=0; i<(MAXCHANNELS); i++) {
      if (i%40 != 39 && i%40 != 29) {  
	tempSum1 += goodTraceE[i]; tempSum1_pu += goodTraceE_PU[i];
      }
      if (i%40 != 39 && i%40 != 29) { 
	tempSum2 += badFPGA_zero[i]; tempSum2_pu += badFPGA_zero_PU[i];
	tempSum2_neg += (badFPGA_neg[i]); tempSum2_neg_pu += (badFPGA_neg_PU[i]);
      }
    }
    cout << "  ---------------------------" << endl;
    cout << "  % segments missing: " 
	 << (float)100*(tempSum2+tempSum2_pu+tempSum2_neg+tempSum2_neg_pu)/(tempSum1+tempSum1_pu) 
	 << "%" << endl;
    cout << "     without PU flag: " << (float)100*(tempSum2)/(tempSum1) << "% zeros/ " 
	 << (float)100*(tempSum2_neg)/(tempSum1) << "% negative" << endl;
    cout << "        with PU flag: " << (float)100*(tempSum2_pu)/(tempSum1_pu) << "% zeros/ "
	 << (float)100*(tempSum2_neg_pu)/(tempSum1_pu) << "% negative" << endl;
    
    tempSum1 = 0;  tempSum1_pu = 0;
    tempSum2 = 0;  tempSum2_pu = 0;
    tempSum2_neg = 0;   tempSum2_neg_pu = 0;
    for (Int_t i=0; i<(MAXCHANNELS); i++) {
      if (i%40 == 29) {  
	tempSum1 += goodTraceE[i]; tempSum1_pu += goodTraceE_PU[i];
      }
      if (i%40 == 29) { 
	tempSum2 += badFPGA_zero[i]; tempSum2_pu += badFPGA_zero_PU[i];
	tempSum2_neg += (badFPGA_neg[i]); tempSum2_neg_pu += (badFPGA_neg_PU[i]);	
      }
    }
    cout << "  ---------------------------" << endl;
    cout << "       % cc1 missing: " 
	 << (float)100*(tempSum2+tempSum2_pu+tempSum2_neg+tempSum2_neg_pu)/(tempSum1+tempSum1_pu) 
	 << "%" << endl;
    cout << "     without PU flag: " << (float)100*(tempSum2)/(tempSum1) << "% zeros/ " 
	 << (float)100*(tempSum2_neg)/(tempSum1) << "% negative" << endl;
    cout << "        with PU flag: " << (float)100*(tempSum2_pu)/(tempSum1_pu) << "% zeros/ "
	 << (float)100*(tempSum2_neg_pu)/(tempSum1_pu) << "% negative" << endl;
    
    tempSum1 = 0;  tempSum1_pu = 0;
    tempSum2 = 0;  tempSum2_pu = 0;
    tempSum2_neg = 0;   tempSum2_neg_pu = 0;
    for (Int_t i=0; i<(MAXCHANNELS); i++) {
      if (i%40 == 39) {  
	tempSum1 += goodTraceE[i]; tempSum1_pu += goodTraceE_PU[i];
      }
      if (i%40 == 39) { 
	tempSum2 += badFPGA_zero[i]; tempSum2_pu += badFPGA_zero_PU[i];
	tempSum2_neg += (badFPGA_neg[i]); tempSum2_neg_pu += (badFPGA_neg_PU[i]);	
      }
    }
    cout << "  ---------------------------" << endl;
    cout << "       % cc2 missing: " 
	 << (float)100*(tempSum2+tempSum2_pu+tempSum2_neg+tempSum2_neg_pu)/(tempSum1+tempSum1_pu) 
	 << "%" << endl;
    cout << "     without PU flag: " << (float)100*(tempSum2)/(tempSum1) << "% zeros/ " 
	 << (float)100*(tempSum2_neg)/(tempSum1) << "% negative" << endl;
    cout << "        with PU flag: " << (float)100*(tempSum2_pu)/(tempSum1_pu) << "% zeros/ "
	 << (float)100*(tempSum2_neg_pu)/(tempSum1_pu) << "% negative" << endl;
    cout << "  ---------------------------" << endl;
  }

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
