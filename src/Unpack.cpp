/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Stefanos Paschalis & Heather Crawford                                    */
/*                                                                          */
/* Unpacking (and optionally event building) GRETINA decomposed data        */
/* and storing them in a ROOT tree, or creating and saving histograms.      */
/*                                                                          */
/* Usage:                                                                   */
/*  Unpack <g, m, or t> <Directory, i.e. CR-5> <Run#s, separated by spaces> */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

/* Standard library includes */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iomanip>
#include "Riostream.h"
#include <vector>
#include <signal.h>
#include <deque>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <fstream>

/* ROOT includes */
#include "TString.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h" 
#include "TSortedList.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCutG.h"

/* Program header files */
#include "Globals.h"
#include "Defines.h"
#include "Unpack.h"

/* GRETINA-specific header files */
#include "GRETINA.h"
#include "SortingStructures.h"
#include "INLCorrection.h"

#include "UnpackUtilities.h"

/* Tracking... */
#include "Track.h"

#include "RdGeCalFile.c"

/* Phoswich Wall header files */
#ifdef WITH_PWALL
#include "PhosWall.h"
#endif

/* CHICO header files */
#ifdef WITH_CHICO
#include "CHICO.h"
#endif

/* S800 header files */
#ifdef WITH_S800
#include "S800Parameters.h"
#include "S800Definitions.h"
#include "S800Functions.h"
#endif

/* BGS header files */
#ifdef WITH_BGS
#include "BGSParamList.h"
#include "BGSConditions.h"
#include "BGSFPPositions.h"
#include "BGSCalibrations.h"
#include "BGSThresholds.h"
#include "BGSAnalyze.h"
#endif

#include "Tree.h"

#define DEBUG2AND3 0

/****************************************************/

void PrintHelpInformation();
void PrintConditions();

void GetData(FILE* inf, controlVariables* ctrl, counterVariables* cnt,
	     GRETINAVariables* gVar, INLCorrection *inlCor, 
	     UShort_t junk[]);

void GetBank29(FILE* inf, UShort_t buf[]);
void ReadMario(FILE* inf);
void SkipData(FILE* inf, UShort_t junk[]);

/****************************************************/

Int_t gotsignal;
void breakhandler(int dummy) {
  cout << "Got break signal.  Aborting sort cleanly..." << endl;
  gotsignal = 1;
}

/****************************************************/

int main(int argc, char *argv[]) {
  
  /* Some CTRL-C interrupt handling stuff... */
  gotsignal = 0;
  signal(SIGINT, breakhandler);
  
  /* When not enough arguments, print the help information */
  if ( (argc<3) ) { PrintHelpInformation(); exit(1); }
  
  /* Initialize analysis control flags to default values,
     then read in the command line arguments. */
  controlVariables *ctrl = new controlVariables();
  ctrl->Initialize();
  Int_t good2Go = ctrl->InterpretCommandLine(argc, argv);
  if (good2Go != 1) { exit(-1); }
  PrintConditions();
  good2Go = ctrl->ReportRunFlags();
  if (good2Go != 1) { exit(-2); }
  printf("\n");

  /* Declare and then initialize GRETINA variables, i.e. geometry
     stuff, and calibration parameters, etc. */
  GRETINAVariables*  gVar = new GRETINAVariables();
  gVar->Initialize();
  gVar->InitializeGRETINAVariables("gretina.set");
  RdGeDINOCalFile(37, "Co60Pairs-dino.dat", gVar);

  gret = new GRETINA();
  gret->Initialize();

  counterVariables *cnt = new counterVariables();

  /* Initialize the GRETINA data structures. */

  /* Superpulse analysis */
  gret->sp.Initialize(ctrl, gVar);
  
  /* INL correction parameters */
  INLCorrection *inlCor = new INLCorrection();
  inlCor->Initialize(ctrl, gVar);

  /* GRETINA Raw data - mode 3... */
  // gMode3Event = new GRawEvent();
   
  cloverEventOUT = new cloverEvent();

  /* Initialize tracking stuff. */
  if (ctrl->doTRACK) {
    gret->track.Initialize();
  }

  /* And data arrays... */
  /* Throw-away/skip data */
  UShort_t junk[8192];

  /* Get the parameters for mapping from crystal coordinate frame
     to the lab frame, for Doppler correction. */
  gret->rot.ReadMatrix("crmat.dat");
  
  /* Get the calibration parameters. */
  if (ctrl->specifyCalibration) {
    RdGeCalFile(ctrl->calibrationFile.Data(), gVar);
  } else {
    cout << "Using default GRETINA energy calibration file. " << endl;
    RdGeCalFile("EhiGainCor-Run0296Scanning.dat", gVar);
  }

  /* Read in baseline values for all segments. These are only
     used with the option "RADFORD_ENERGY". */
  //if (ctrl->withWAVE) {
  // if (ctrl->RADFORD_ENERGY) {
  //  RdGeBaseLineFile("Baselines.out");
  //  for (Int_t i=0; i<MAXCHANNELS; i++) {
  //	if (i%40 == 9 || i%40 == 19) {
  //	  gWf->tau[i] = gWf->tau[i+20];
  //	}
  //  }
  //}
  //}
  cout << endl;

  /* And now the auxiliary detector system data structures.  In alphabetical order :) */

  /* CHICO */
#ifdef WITH_CHICO
  chico = new CHICOFull();
  chico->Initialize();
  chico->InitializeCHICOVariables("ppacTheta.cal", "ppacPhi.cal", 
                                  "beta_76Ge_May14.dat");
  chico->offsetTarget = 15.8; 
#endif
  
  /* Phoswich Wall */
#ifdef WITH_PWALL
  phosWall = new phosWallFull();
  phosWall->InitializeParameters("phosParam.set");
  phosWall->Initialize();
  phosWall->SetPwallPositions();
#endif

  /* S800 */
#ifdef WITH_S800
  s800 = new S800Full();
  s800->Initialize();

  if (ctrl->s800File) {
    ctrl->SetS800Controls(ctrl->s800ControlFile);
    s800->InitializeS800Variables(ctrl->s800VariableFile);
  } else { s800->InitializeS800Variables("s800.set"); }
  
  s800->fp.ic.BuildLookUp();
  s800->fp.track.map.LoadInverseMap(s800->fp.track.map.mapFilename);
  s800->fp.crdc1.calc.CheckBadPads();
  s800->fp.crdc2.calc.CheckBadPads();
  s800->fp.crdc1.pad.BuildLookUp();
  s800->fp.crdc2.pad.BuildLookUp();
#ifdef S800_LINK_IMAGE_TRACK
  s800->im.track.map.LoadInverseMap(s800->im.track.map.imMapFilename);
  s800->im.track.tppac1.CheckBadStrips();
  s800->im.track.tppac2.CheckBadStrips();
#endif

  s800Scaler = new S800Scaler();
  s800Scaler->Initialize();

  s800Scaler->InitializeS800ScalerParameters("S800Scaler.definition");
  cout << endl;
#endif /* WITH_S800 */

  FILE *inf;
  
  FILE *generalOut = NULL;
  if (ctrl->outputON) {
    if (ctrl->outputName) {
      generalOut = fopen(ctrl->outputFileName.Data(), "wb");
      if (!generalOut) {
	printf("Cannot open general output file: %s \n", ctrl->outputFileName.Data());
	exit(2); 
      }
    } else {
      generalOut = fopen("GeneralFile.out", "wb");
      if (!generalOut) {
	printf("Cannot open general output file: GeneralFile.out \n");
	exit(2); 
      }
    }
  }

#ifdef WITH_S800
  /* Here, read in ROOT gates (PID gates) that will be used to make 
     filtered files, or just cut the crap we don't need out of the 
     ROOT trees. */
  TCutG *incomingBeam = new TCutG();
  TCutG *outgoingBeam = new TCutG();

  TCutG *incomingAll = new TCutG("incomingAll", 4);
  incomingAll->SetVarX("tof.xfp");
  incomingAll->SetVarY("tof.obj");
  incomingAll->SetPoint(0, -100000, 100000);
  incomingAll->SetPoint(1, 100000, 100000);
  incomingAll->SetPoint(2, 100000, -100000);
  incomingAll->SetPoint(3, -100000, -100000);
  
  TCutG *outgoingAll = new TCutG("outgoingAll", 4);
  outgoingAll->SetVarX("tof.obje1");
  outgoingAll->SetVarY("fp.ic.de");
  outgoingAll->SetPoint(0, -100000, 100000);
  outgoingAll->SetPoint(1, -100000, -100000);
  outgoingAll->SetPoint(2, 100000, -100000);
  outgoingAll->SetPoint(3, 100000, 100000);

  incomingBeam = incomingAll;
  outgoingBeam = outgoingAll;

  /* For now, include the user gates here...eventually we need to pull this
     out of the compiled code.  I need to think about this one... */
#include "ROOTGates.var"

#endif /* WITH_S800 */

  TStopwatch timer;

  /* Loop over each run given at the command line. */
  if (ctrl->fileType == "f1" || ctrl->fileType == "f2" || 
      ctrl->fileType == "f") {
    ctrl->startRun = 0; argc = 1; 
  } /* For specific file name, once through loop only. */
  
  for (Int_t mm = ctrl->startRun; mm<argc; mm++) {
    if (!gotsignal) { /* We haven't aborted for some reason. */
      
      timer.Reset(); timer.Start();
      cnt->ResetRunCounters();
      
      TString runNumber = argv[mm];
      cnt->runNum = atoi(argv[mm]);
      
      Int_t fileOK = OpenInputFile(&inf, ctrl, runNumber);
      if (fileOK != 0) { exit(2); }
      
#ifdef WITH_S800
      if (ctrl->fileType != "f") {
	TString runVariableFileName = ctrl->directory + "Run" + runNumber + "/Run" + runNumber + ".var";
	cout << "Looking for run variable file " << runVariableFileName.Data() << endl;
	s800->UpdateS800RunVariables(runVariableFileName);
      }
#endif
      
      /* Open output file, set up tree and/or histograms. */
      TFile *fout_root = NULL;
      if (ctrl->withTREE || ctrl->withHISTOS) {
	fout_root = new TFile(ctrl->outfileName.Data(), "RECREATE");
	fout_root->SetCompressionAlgorithm(1);
	fout_root->SetCompressionLevel(2);
	cout << "Output file: " << ctrl->outfileName << " (compression " 
	     << fout_root->GetCompressionLevel() << ")" << endl;
      } else {
	cout << "No ROOT output requested -- no histos or trees. " << endl;
      }
      
      if (ctrl->withTREE) { InitializeTree(ctrl); }
      
      // teb->SetMaxTreeSize(1000000000LL); /* Max tree size is 1GB */

      cout << "********************************************************" << endl;
      cout << endl;
      
      /* Reset variables needed for unpacking, histogramming, etc. */
      cnt->ResetRunCounters();
      Bool_t GO_FOR_BUILD = 1;
      
      Int_t TSerrors = 0;
      long long int currTS = 0;  long long int lastTS = 0;
      long long int deltaEvent = 0;        
      Int_t siz = 0; 
      Int_t remaining = 0;
      Int_t s800Length = 0;
      Int_t timeToOptimize = 0;
      Int_t atSTARTFile2 = 1; Int_t BonusMode3 = 0;
      
      /********************************************************/
      /*  THE MAIN EVENT -- SORTING LOOP                      */
      /********************************************************/
      
      /* Loop over file, reading data, and building events... */
      if (ctrl->pgh == 0) { /* We expect global headers, so read one */
	siz = fread(&gHeader, sizeof(struct globalHeader), 1, inf);
      }
      
      while (siz && !gotsignal) {

	if (cnt->TSFirst == 0 && gHeader.timestamp > 0) { cnt->TSFirst = gHeader.timestamp; }
	cnt->TSLast = gHeader.timestamp;

	gHeader.timestamp++; /* For simulation - first TS is 0, so just
				adding 1 solves a bunch of problems... */

	if (ctrl->noEB) { /* Just get the data, don't event build. */

	  GetData(inf, ctrl, cnt, gVar, inlCor, junk);

	  /* Fill singles spectra as appropriate */
	  if (ctrl->withHISTOS && ctrl->calibration && !ctrl->xtalkAnalysis) { 
	    gret->fillHistos(1); 
	  } 
	  ResetEvent(ctrl, cnt);
	  cnt->event = 0.0;

	} else { /* !ctrl->noEB */
	  
	  /* We are building events (i.e. grouping by timestamp) -- 
	     Check against timestamps in file being out of order... */
	  
	  if (gHeader.timestamp < lastTS) {
	    if (!ctrl->suppressTS) {
	      cout << ALERTTEXT;
	      printf("Unpack(): TS out of order: lastTS %lld, current %lld\n",
		     lastTS, gHeader.timestamp);
	      cout << RESET_COLOR;  fflush(stdout);
	    }
	    TSerrors++;
	  }
	  lastTS = gHeader.timestamp;
	  
	  if ( (gHeader.timestamp != 0) && (currTS == 0) ) {
	    currTS = gHeader.timestamp;
	  }
	  
	  /* Just throw out stupid events with no valid TS. */
	  if (gHeader.timestamp != 0 && remaining != 5) { 
	    GO_FOR_BUILD = 1; 
	  } else { GO_FOR_BUILD = 0; }
	  
	  if (GO_FOR_BUILD) {
	    
	    /* Build events based on timestamp differences between the
	       start of an event and current timestamp. */
	    deltaEvent = (Float_t)(gHeader.timestamp - currTS);
	    
	    if (abs(deltaEvent) < EB_DIFF_TIME) {
	      
	      GetData(inf, ctrl, cnt, gVar, inlCor, junk);
	      
	      if (ctrl->superPulse) {
		if (gHeader.type == RAW) {
		  gret->sp.trLength = gret->g3Temp[0].wf.raw.size();
		}
	      }

	      /* Add to event, as necessary -- trying to get rid of this... */
	      switch(gHeader.type) {
	      case 0xC10:
		cloverEventOUT->obj.push_back(clover);
		break;
	      default:  break;
	      }
	      
	    } else { /* Time difference is big...old event should be closed. */
	      
	      /* We need to be careful of tracelengths in superPulse analysis... */
	      if (ctrl->superPulse) {
		if (gHeader.type == RAW) {
		  gret->sp.trLength = gret->g3Temp[0].wf.raw.size();
		}
	      }

	      if (ctrl->gateTree) {
#ifdef WITH_S800
		Int_t pidOK = CheckS800PIDGates(incomingBeam, outgoingBeam);
		if (pidOK) { 
		  int evtOK = ProcessEvent(currTS, ctrl, cnt, gVar);
		  if (evtOK < 0) { raise(SIGINT); }
		}
#endif
	      } else {
		Int_t evtOK = ProcessEvent(currTS, ctrl, cnt, gVar);
		if (evtOK < 0) { raise(SIGINT); }
	      }
	      
	      /* Check on gate conditions...do we write filtered output? */
	      if (ctrl->outputON) {
#ifdef WITH_S800
		Int_t writeOK = CheckS800PIDGates(incomingBeam, outgoingBeam);
		if (writeOK) { 
		  // WriteS800PhysicsFile(generalOut); 
		}
#endif
	      }
	      ResetEvent(ctrl, cnt);
	      cnt->event = 0.0;
	      
	      /* Update timestamp information for beginning of new event, 
		 and start filling new event. */
	      currTS = gHeader.timestamp;
	      deltaEvent = (Float_t)(gHeader.timestamp - currTS);
	      
	      GetData(inf, ctrl, cnt, gVar, inlCor, junk);
	      
	      switch(gHeader.type) {
	      case 0xC10:
		cloverEventOUT->obj.push_back(clover);
		break;
	      default:  break;
	      }
	    }
	  } else { /* End of "if (GO_FOR_BUILD)" */	
	    SkipData(inf, junk);
	    cnt->Increment(gHeader.length);
	  }
	} /* End of !ctrl->noEB */	  

	if (cnt->bytes_read_since_last_time > 2*1024*1024) { 
	  if (cnt->bytes_read > 1024*1024*1024) {
	    cerr << "Processing " << (Float_t)cnt->bytes_read/(1024.*1024.*1024.) << " GB         " <<"\r";
	  } else { 
	    cerr << "Processing " << cnt->bytes_read/(1024*1024) << " MB    " <<"\r";
	  }
	  cnt->MBread+=2;
	  cnt->bytes_read_since_last_time = 0;
	}
	
	siz = fread(&gHeader, sizeof(struct globalHeader), 1, inf);
	
      } /* End of "while we still have data and no interrupt signal" */
      
      if (ctrl->superPulse) {
	gret->checkSPIntegrity();
	gret->sp.MakeSuperPulses();
	gret->sp.FinishSuperPulses();
	gret->sp.WriteSuperPulses();
      }
      
      if (ctrl->outputON) {
#ifdef WITH_S800
	Int_t writeOK = CheckS800PIDGates(incomingBeam, outgoingBeam);
	if (writeOK) { 
	  //WriteS800PhysicsFile(generalOut); 
	  //WriteS800GatedFile(generalOut, rawS800, s800Length); 
	}
#endif
      }
      
      /* Write the last event... */
      if (ctrl->gateTree) {
#ifdef WITH_S800
	Int_t pidOK = CheckS800PIDGates(incomingBeam, outgoingBeam);
	if (pidOK && ctrl->withTREE) { teb->Fill(); }	
#endif
      } else {
	if (ctrl->withTREE) { teb->Fill(); }
      }
      
      /* Do David's baseline fitting routine, if that's what we're looking at. */    
      //if (ctrl->RADFORD_BASELINE) {
      //	gWf->RadfordBaselineFit("Baselines.out");
      //}
      
      timer.Stop();
      
      cout << "\n CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
      cout << " Average processing speed: " << (cnt->bytes_read/(1024*1024))/timer.RealTime() 
	   << "MB/s -- File size was " << cnt->bytes_read/(1024*1024) << " MB \n" << endl;
      
      cnt->PrintRunStatistics(ctrl->pgh, ctrl->withWAVE, ctrl->superPulse, ctrl->analyze2AND3);
      cnt->ResetRunCounters();
      
      if (ctrl->withTREE) { 
	cout << endl << "Writing ROOT tree..." << endl;
	teb->Write(); 
#ifdef WITH_S800
	//scaler->Write();
#endif
	if (ctrl->withWAVE) {
	  if (ctrl->WITH_TRACETREE) {
	    wave->Write();
	  }
	}
      }
      
      if (ctrl->withHISTOS) { 
	cout << "Writing histograms..." << endl;
	gret->gHist.writeHistos(1); 
      }
      if (ctrl->withHISTOS || ctrl->withTREE) {
	printf("ROOT file \"%s\" closing...\n", ctrl->outfileName.Data());
	//fout_root->Write();
	fout_root->Close();
      }
      
      cout << "*******************************************************" << endl;
      cout << endl;
      
    } /* End of if !gotsignal */
  } /* End of iterating over different run # */

  if (ctrl->outputON) {
    cout << endl;
    cout << "Closing output file...";
    fclose(generalOut);
    cout << "Done. " << endl;    
  }
  
  /* Declare victory!!! */
  
  cout << endl;
  timer.Delete();
  cout << "We finished!! Yay us! :)" << endl;
  cout << endl;
  
  return 1;
}

/****************************************************/

void GetData(FILE* inf, controlVariables* ctrl, counterVariables* cnt,
	     GRETINAVariables* gVar, INLCorrection *inlCor, UShort_t junk[]) {
  
  cnt->Increment(sizeof(struct globalHeader));
  
  switch(gHeader.type) {
    
  case DECOMP:
    { gret->getMode2(inf, gHeader.length, gVar, cnt); }
    break;
  case TRACK:
    { gret->getMode1(inf, gVar, cnt); }
    break;
  case RAW:
    { gret->getMode3(inf, gHeader.length, cnt, ctrl, gVar); }
    break;
#ifdef WITH_BGS
  case BGS:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
#else 
  case BGS:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
#endif
#ifdef WITH_S800
  case S800:
    { s800->getAndProcessS800(inf, gHeader.length);  cnt->Increment(gHeader.length); }
    break;
  case S800AUX:
    {
      Bool_t scalerPacket = s800Scaler->getAndProcessS800Aux(inf, gHeader.length, gHeader.timestamp);
      if (scalerPacket && ctrl->withTREE) { /* scaler->Fill(); */ }
      cnt->Increment(gHeader.length);
    }
    break;
  case S800AUX_TS:
    {
      Bool_t scalerPacket = s800Scaler->getAndProcessS800Aux(inf, gHeader.length, gHeader.timestamp);
      if (scalerPacket && ctrl->withTREE) { /* scaler->Fill(); */ }
      cnt->Increment(gHeader.length);
    }
    break;
  case S800PHYSICS:
    {
      SkipData(inf, junk);
      //printf("%0.3f\n", s800->fp.crdc1.x);
      cnt->Increment(gHeader.length);
    }
    break; 
#else 
  case S800:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case S800AUX:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case S800AUX_TS:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case S800PHYSICS:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
#endif
  case BANK29:
    { GetBank29(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case GRETSCALER:
    { gret->getScaler(inf, gHeader.length);  cnt->Increment(gHeader.length); }
    break;
  case G4SIM:
    { gret->getSimulated(inf);  cnt->Increment(gHeader.length); }
    break;
#ifdef WITH_CHICO
  case CHICO: { chico->getAndUnpackCHICO(inf, gHeader.length);  cnt->Increment(gHeader.length); }
    break;
#else
  case CHICO:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
#endif
#ifdef WITH_PWALL
  case PWALL:
    {
      phosWall->getAndUnpackPhosWall(inf, gHeader.length);
      phosWall->ProcessPhosWall();
      phosWall->timestamp = gHeader.timestamp;
      cnt->Increment(gHeader.length);
    }
    break;
  case PWALLAUX:
    {
      phosWall->getAndUnpackPhosWallAux(inf, gHeader.length);
      phosWall->ProcessPhosWallAux();
      cnt->Increment(gHeader.length);
    }
    break;
#else 
  case PWALL:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case PWALLAUX:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
#endif
  case 0xC10:
    { 
      Int_t siz = fread(&clover, 1, gHeader.length, inf);
      if (siz != gHeader.length) {
	cerr << "GetClover failed in bytes read." << endl;
      }
      cnt->Increment(gHeader.length); 
    }
    break;
  case 0:
    {
      cout << "GlobalHeader type = 0.  Ignoring." << endl;
      cnt->Increment(gHeader.length);
    }
    break;
  case 100:
    {
      cout << "GlobalHeader Mario found. " << endl;
      ReadMario(inf);
      cnt->Increment(gHeader.length);
    }
    break;
  default:
    {
      cout << "GlobalHeader type not recognized: " << gHeader.type << endl;
      SkipData(inf, junk);
      cnt->Increment(gHeader.length);
    }
    break;
  }   

  if (gHeader.type < 50) {
    cnt->headerType[gHeader.type]++;
    cnt->setEventBit(gHeader.type);
  }
}

/****************************************************/

void GetBank29(FILE* inf, UShort_t buf[]) {
  gret->b29.timestamp = gHeader.timestamp;

  /* Skip the actual data for now. Maybe we'll do 
     something like CFD down the road. */
  Int_t siz = fread(buf, 1, gHeader.length, inf);
  if (siz != gHeader.length) {
    cerr << "GetBank29 failed in bytes read." << endl;
  }
}

/****************************************************/

struct out4Mario {
  float ccEnergy;
  float segEnergy[36];
  float pad;
  short wf[37][300];
};

void ReadMario(FILE* inf) {
  out4Mario mario;
  fread(&mario, 1, sizeof(struct out4Mario), inf);
  cout << mario.ccEnergy << endl;
  for (Int_t i=0; i<36; i++) {
    cout << mario.segEnergy[i] << endl;
  }
}

/****************************************************/

void SkipData(FILE* inf, UShort_t junk[]) {
  Int_t siz = fread(junk, 1, gHeader.length, inf);
  if (siz != gHeader.length) {
    cout << ALERTTEXT;
    printf("SkipData(): Failed.\n");
    cout << RESET_COLOR;  fflush(stdout);
  } 
}

/****************************************************/

void PrintHelpInformation() {
  printf("\n");
  printf("Usage: ./Unpack <File Type> <Usage flags> -d <Subdirectory, i.e. CR-5> -run <input run ### (separate multiple files by a space)>\n");
  printf("    Valid file types: -g  (Global.dat)\n");
  printf("                      -gr (GlobalRaw.dat)\n");
  printf("                      -m  (Merged.dat, externally merged file, usually GRETINA mode 3)\n");
  printf("                      -m2 (Merged.Mode2.dat, GRETINA mode2 externally merged file)\n");
  printf("                      -f  <FILENAME> (specific file -- do not require subdirectory or run #)\n");
  printf("    Valid usage flags: -preGH2 (file is Mode2, from the pre global header (pGH) era) \n");
  printf("                       -preGH3 (file is Mode3, from the pre global header (pGH) era)\n");
  printf("                       -wf (do waveform analysis; algorithms to run flagged in Unpack.h)\n");
  printf("                       -noEB (DISABLE event building; default is to BUILD events)\n");
  printf("                       -calibrationRun (DISABLE event building, and fill calibration histograms)\n");
  printf("                       -superPulse <lowE> <highE> <detMapFile directory> <XtalkFile>\n");
  printf("                               (turns off tree and histograms, builds superpulse .spn files)\n");
  printf("                       -noSeg (DISABLE segment analysis, i.e. segment summing; default is ON \n");   
#ifdef WITH_S800
  printf("                       -s800File (define the s800 control file; without default is s800.set\n                                  for S800 parameters, and NO s800 included in ROOT tree)\n");
#endif
  printf("                       -track (do tracking, such as it is -- options specified in track.chat)\n");
#ifdef WITH_S800
  printf("                       -outputON (write output file, with S800Physics, gated on PID)\n");
#endif
  printf("                       -zip (compressed data file, will append .gz to filename)\n");
  printf("                       -bzip (compressed data file, will append .bz2 to filename)\n");
  printf("                       -noHistos (turn OFF histograms, as defined in Histos.h, default is ON)\n");
  printf("                       -noTree (turn OFF root tree, default is ON)\n");
  printf("                       -noHFC (turn OFF HFC presort piping; default is ON)\n");
  printf("                       -dopplerSimple (use simple GRETINA decomp position for Doppler correction)\n");
  printf("                       -suppressTS (suppress TS error warnings in analysis)\n");
  /* 2015-04-21 CMC added command-line flag descriptions, as I understand them, feel free to correct or update */
  printf("                       -outputName <FILENAME> (set the PID-gated output event file name)\n");
  printf("                       -rootName <FILENAME> (set the output ROOT file name)\n");
  printf("                       -analyze2and3 (analyze Mode2 and Mode3, matching by timestamps)\n");
  printf("                       -gateTree (gates tree and histogramm by a PID gate)\n");
  printf("                       -readCal <FILENAME> (read in a calibration file)\n");
  printf("                       -xtalkAnalysis <ID> <lowE> <highE> \n");
  printf("                               (turns off tree, builds histograms for segment cross-talk analysis)\n");
  printf("                       -INLcorrect <INLCConly> <digMapFileName> \n");
  printf("                               (turns on waveform analysis, attempts to correct energy filter INL)\n");
  printf("\n");
}

/****************************************************/

void PrintConditions() {
  printf("\n***************************************************************************\n\n");
  printf("    Initializing -- GRETINA ");
#ifdef WITH_BGS
  printf("+ BGS ");
#endif
#ifdef WITH_CHICO
  printf("+ CHICO2 ");
#endif
#ifdef WITH_PWALL
  printf("+ PHOSWALL ");
#endif
#ifdef WITH_S800
  printf("+ S800 ");
#endif
#ifndef WITH_BGS
#ifndef WITH_CHICO
#ifndef WITH_PWALL
#ifndef WITH_S800
  printf("only ");
#endif
#endif
#endif
#endif
  printf("sort...\n");

}
