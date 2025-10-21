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

#include "UnpackUtilities.h"

#include "Tree.h"

/****************************************************/

void PrintHelpInformation();
void GetData(FILE* inf, controlVariables* ctrl, counterVariables* cnt,
	     UShort_t junk[]);

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
  good2Go = ctrl->ReportRunFlags();
  if (good2Go != 1) { exit(-2); }
  printf("\n");

  gret = new GRETINA();
  gret->Initialize();

  counterVariables *cnt = new counterVariables();

  /* Initialize the GRETINA data structures. */

  /* And data arrays... */
  /* Throw-away/skip data */
  UShort_t junk[8192];


  FILE *inf;
  
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
      
      if (ctrl->withTREE) { 
	InitializeTree(); 
      }

      UInt_t counter = 0;
      
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


	//printf("%lld %d %d\n", gHeader.timestamp, gHeader.type, gHeader.length);
	
	if (ctrl->noEB) { /* Just get the data, don't event build. */

	  GetData(inf, ctrl, cnt, junk);
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
	      
	      GetData(inf, ctrl, cnt, junk);
	      	      
	    } else { /* Time difference is big...old event should be closed. */

	      counter++;
	     
	      Int_t evtOK = ProcessEvent(currTS, ctrl, cnt);
	      if (evtOK < 0) { raise(SIGINT); }
	      
	      ResetEvent(ctrl, cnt);
	      cnt->event = 0.0;
	      
	      /* Update timestamp information for beginning of new event, 
		 and start filling new event. */
	      currTS = gHeader.timestamp;
	      deltaEvent = (Float_t)(gHeader.timestamp - currTS);
	      
	      GetData(inf, ctrl, cnt, junk);
	      
	    }
	  } else { /* End of "if (GO_FOR_BUILD)" */	
	    SkipData(inf, junk);
	    cnt->Increment(gHeader.length);
	  }
	} /* End of !ctrl->noEB */	  

	//printf("%d \n", cnt->bytes_read_since_last_time);
	
	if (cnt->bytes_read_since_last_time > 2*1024*1024) { 
	  if (cnt->bytes_read > 1024*1024*1024) {
	    cout << "Processing " << (Float_t)cnt->bytes_read/(1024.*1024.*1024.) << " GB         " <<"\r";
	  } else { 
	    cout << "Processing " << cnt->bytes_read/(1024*1024) << " MB    " <<"\r";
	  }
	  cnt->MBread+=2;
	  cnt->bytes_read_since_last_time = 0;
	}
	
	siz = fread(&gHeader, sizeof(struct globalHeader), 1, inf);
	
      } /* End of "while we still have data and no interrupt signal" */
      
      
      /* Write the last event... */
      if (ctrl->withTREE) { teb->Fill();  cnt->treeWrites++; }
      
      timer.Stop();
      
      cout << "\n CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
      cout << " Average processing speed: " << (cnt->bytes_read/(1024*1024))/timer.RealTime() 
	   << "MB/s -- File size was " << cnt->bytes_read/(1024*1024) << " MB \n" << endl;
      
      cnt->PrintRunStatistics(ctrl->pgh, ctrl->withWAVE, ctrl->superPulse, ctrl->analyze2AND3);
      cnt->ResetRunCounters();
      
      if (ctrl->withTREE) { 
	cout << endl << "Writing ROOT tree..." << endl;
	teb->Write(); 
	if (ctrl->withWAVE) {
	  if (ctrl->WITH_TRACETREE) {
	    wave->Write();
	  }
	}
      }
      
      if (ctrl->withTREE) {
	printf("ROOT file \"%s\" closing...\n", ctrl->outfileName.Data());
	//fout_root->Write();
	fout_root->Close();
      }
      
      cout << "*******************************************************" << endl;
      cout << endl;
      
    } /* End of if !gotsignal */
  } /* End of iterating over different run # */

 
  /* Declare victory!!! */
  cout << endl;
  timer.Delete();

  return 1;
}

/****************************************************/

void GetData(FILE* inf, controlVariables* ctrl, counterVariables* cnt,
	     UShort_t junk[]) {
  
  cnt->Increment(sizeof(struct globalHeader));

  switch(gHeader.type) {
    
  case DECOMP:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case TRACK:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case RAW:
    { 
      if (cnt->headerType[RAW] == 0 && ctrl->withTREE) {
	InitializeTreeMode3();
	for (Int_t i=0; i<cnt->treeWrites; i++) {
	  teb->FindBranch("g3")->Fill();
	}
      }
      gret->getMode3(inf, gHeader.length, cnt, ctrl); 
    }
    break;
  case RAWHISTORY:
    { 
      if (cnt->headerType[RAWHISTORY] == 0 && ctrl->withTREE) {
	InitializeTreeHistory();
	for (Int_t i=0; i<cnt->treeWrites; i++) {
	  teb->FindBranch("g3H")->Fill();
	}
      }
      gret->getMode3History(inf, gHeader.length, gHeader.timestamp, cnt); 
    }
    break;
  case BGS:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
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
  case BANK88:
    { 
      if (cnt->headerType[BANK88] == 0 && ctrl->withTREE) {
	InitializeTreeBank88();
	for (Int_t i=0; i<cnt->treeWrites; i++) {
	  teb->FindBranch("b88")->Fill();
	}
      }
      gret->getBank88(inf, gHeader.length, cnt);  cnt->Increment(gHeader.length); 
    }
    break;
  case GRETSCALER:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case G4SIM:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); } 
    break;
  case CHICO:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case PWALL:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case PWALLAUX:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case GODDESS:
    { SkipData(inf, junk);  cnt->Increment(gHeader.length); }
    break;
  case LENDA:
    { SkipData(inf, junk); cnt->Increment(gHeader.length); }
    break;
  case DFMA:
    { SkipData(inf, junk); cnt->Increment(gHeader.length); }
    break;
  case 0:
    {
      cout << "GlobalHeader type = 0.  Ignoring." << endl;
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
  if (gHeader.type == RAWHISTORY) {
    cnt->setEventBit(RAW);
    cnt->headerType[RAW]++;
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
  printf("                      -f  <FILENAME> (specific file -- do not require subdirectory or run #)\n");
  printf("    Valid usage flags: -wf (do waveform analysis; algorithms to run flagged in Unpack.h)\n");
  printf("                       -noEB (DISABLE event building; default is to BUILD events)\n");
  printf("                       -zip (compressed data file, will append .gz to filename)\n");
  printf("                       -bzip (compressed data file, will append .bz2 to filename)\n");
  printf("                       -noHFC (turn OFF HFC presort piping; default is ON)\n");
  printf("                       -suppressTS (suppress TS error warnings in analysis)\n");
  printf("                       -rootName <FILENAME> (set the output ROOT file name)\n");
  printf("\n");
}

