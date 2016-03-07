/****************************************************************************/
/****************************************************************************/
/*                                                                          */
/* Stefanos Paschalis & Heather Crawford                                    */
/*                                                                          */
/* Adding global headers to old data for ease of analysis, signal decomp,   */
/* and other stuff.  This is to make life easier with old data.             */
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
#include "GlobalsAGH.h"
#include "Unpack.h"

/* GRETINA-specific header files */
#include "GRETINA.h"
#include "SortingStructures.h"

/****************************************************/

void PrintHelpInformation();

int GetGRETINA(FILE* inf, FILE *outfile, unsigned char gBuf[], int sizegBuf,
	       int type, int remaining, controlVariables* ctrl, 
	       counterVariables* cnt);

/****************************************************/

Int_t gotsignal = 0;

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
  if ( (argc<3) ) {
    PrintHelpInformation();
    exit(1);
  }
  
  /* Initialize analysis control flags to default values,
     then read in the command line arguments. */
  controlVariables *ctrl = new controlVariables();
  ctrl->Initialize();
  Int_t good2Go = ctrl->InterpretCommandLine(argc, argv);
  if (good2Go != 1) { exit(-1); }

  cout << endl;
  cout << "********************************************************" << endl;
  cout << endl;
  
  cout << endl;
  
  good2Go = ctrl->ReportRunFlags();
  if (good2Go != 1) { exit(-2); }
  
  cout << endl;

  TStopwatch timer;

  /* Declare and then initialize GRETINA variables, i.e. geometry
     stuff, and calibration parameters, etc. */
  counterVariables *cnt = new counterVariables();

  /* Initialize the GRETINA data structures. */

  /* And data arrays... */
  /* GRETINA */
  unsigned char gBuf[32*32*1024];

  FILE *inf = NULL;

  FILE *generalOut = NULL;
  
  /* Loop over each run given at the command line. */
  if (ctrl->fileType == "f") { ctrl->startRun = 0; argc = 1; } /* For specific file name, 
								  once through loop only. */

  for (Int_t mm = ctrl->startRun; mm<argc; mm++) {
    
    timer.Reset();
    timer.Start();
    
    cnt->ResetRunCounters();

    TString runNumber = argv[mm];
    cnt->runNum = atoi(argv[mm]);
    
    /* Figure out the filename, etc. */
    if (ctrl->fileType == "g") {
      ctrl->fileName = ctrl->directory + "/Run" + runNumber + "/Global.dat";
      ctrl->outfileName = "./ROOTFiles/Run" + runNumber + "wGH.root";
    } else if (ctrl->fileType == "f") {
      ctrl->outfileName = "./ROOTFiles/test-wGH.root";
    } else {
      cerr << "WHAT???" << endl;
      return 0;
    }

    cout << "********************************************************" << endl;
    cout << endl;

    /* Open the input data file... */
    if (ctrl->compressedFile) {
      ctrl->fileName = "zcat " + ctrl->fileName;
      if (ctrl->fileName.EndsWith(".gz")) { }
      else { ctrl->fileName = ctrl->fileName + ".gz"; }
      inf = popen(ctrl->fileName.Data(), "r");
    } else if (ctrl->compressedFileB) {
      ctrl->fileName = "bzcat " + ctrl->fileName;
      if (ctrl->fileName.EndsWith(".bz2")) { }
      else { ctrl->fileName = ctrl->fileName + ".bz2"; }
      inf = popen(ctrl->fileName.Data(), "r");
    } else {
      inf = fopen(ctrl->fileName.Data(), "r");
    }

    ctrl->outputON = 1;
    
    if (ctrl->outputON) {
      if (ctrl->outputName) {
	generalOut = fopen(ctrl->outputFileName.Data(), "wb");
      } else if (!ctrl->outputName) {
	ctrl->outputFileName = ctrl->fileName;
	if (ctrl->fileName.EndsWith(".bz2")) { 
	  ctrl->outputFileName.ReplaceAll(".dat.bz2","GH.dat");
	} else if (ctrl->fileName.EndsWith(".gz")) {
	  ctrl->outputFileName.ReplaceAll(".dat.gz","GH.dat");
	} else {
	  ctrl->outputFileName.ReplaceAll(".dat","GH.dat");
	}
	ctrl->outputFileName.ReplaceAll("zcat ", "");
	generalOut = fopen(ctrl->outputFileName.Data(), "wb");
      }
      if (!generalOut) {
	cout << "Could not open output file " << ctrl->outputFileName.Data() << endl;
	exit(2); 
      } else {
	cout << "Opened output file " << ctrl->outputFileName.Data() << endl;
      } 
    } else {
      generalOut = NULL;
    }
    
    /* Reset variables needed for unpacking, histogramming, etc. */
    cnt->ResetRunCounters();
    
    Int_t TSerrors = 0;
    long long int currTS = 0;  long long int lastTS = 0;
    
    Int_t remaining = 0;
    
    /********************************************************/
    /*  THE MAIN EVENT -- SORTING LOOP                      */
    /********************************************************/

    /* Loop over file, reading data, and building events... */
    while (!gotsignal) {

      if (ctrl->pgh == 2) {
	
	gHeader.type = DECOMP;
	remaining = GetGRETINA(inf, generalOut, gBuf, sizeof(gBuf), gHeader.type, 0, ctrl, cnt);
	gHeader.timestamp = gMode2.timestamp;

      } else if (ctrl->pgh == 3) {

	gHeader.type = RAW;
	remaining = GetGRETINA(inf, generalOut, gBuf, sizeof(gBuf), gHeader.type, 0, ctrl, cnt);
	gHeader.timestamp = gMode3.led_timestamp;

      }

      /* Check against timestamps in file being out of order... */
      if (gHeader.timestamp < lastTS) {
	if (!ctrl->suppressTS)
	  cout << "TS out of order : lastTS " << lastTS << " current " 
	       << gHeader.timestamp << endl;
	TSerrors++;
      }
      lastTS = gHeader.timestamp;
	
      if ( (gHeader.timestamp != 0) && (currTS == 0) ) {
	currTS = gHeader.timestamp;
      }
      
      if (cnt->bytes_read_since_last_time > 2*1024*1024) {    
	cerr << "Processing " << cnt->bytes_read/(1024*1024) << " MB" <<"\r";
	cnt->MBread+=2;
	cnt->bytes_read_since_last_time = 0;
      }
      
      /**********************************************************/
      /* More GRETINA information before another global header. */
      /**********************************************************/
      while (remaining && !gotsignal) { 

	remaining = GetGRETINA(inf, generalOut, gBuf, sizeof(gBuf), gHeader.type, 
			       remaining, ctrl, cnt);
	if (ctrl->pgh == 2) { gHeader.timestamp = gMode2.timestamp; }
	else if (ctrl->pgh == 3) { gHeader.timestamp = gMode3.led_timestamp; }

	if (cnt->bytes_read_since_last_time > 2*1024*1024) {    
	  cerr << "Processing " << cnt->bytes_read/(1024*1024) << " MB" << "\r";
	  cnt->MBread+=2;
	  cnt->bytes_read_since_last_time = 0;
	}
      } /* end of ' while (remaining) ' */

    } /* End of "while we still have data and no interrupt signal" */
    
    timer.Stop();
    
    cout << "\n CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;
    cout << " Average processing speed: " << (cnt->bytes_read/(1024*1024))/timer.RealTime() 
	 << "MB/s -- File size was " << cnt->bytes_read/(1024*1024) << " MB \n" << endl;
    
    cnt->PrintRunStatistics(ctrl->pgh, ctrl->withWAVE, ctrl->superPulse, ctrl->analyze2AND3);
    cnt->ResetRunCounters();

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
  cout << "We finished without crashing!! Yay us! :)" << endl;
  cout << endl;

  return 1;
}

/****************************************************/

int GetGRETINA(FILE* inf, FILE* outfile, unsigned char gBuf[], int sizegBuf, 
	       int type, int remaining, controlVariables* ctrl, 
	       counterVariables* cnt) {
  
  Int_t j, siz;

  Int_t debugAdd = 0;
  unsigned char gOutBuf[32*32*1024];

  if (type == RAW) {
    
    mode3DataPacket *dp;
    if ( !(dp = 
	   (mode3DataPacket*)malloc(sizeof(dp->aahdr) +
				    sizeof(dp->hdr) +
				    sizeof(dp->waveform))) ) {
      exit (-1);
    }    
    
    if (ctrl->pgh) { /* ctrl->pgh == pre global headers */
      
      /* No global headers, old school unpacking! */
      
      if (remaining == 0) { /* We need to read in more data. */
	if (debugAdd) { cout << "Remaining = 0, getting data... " << endl; }

	siz = fread(gBuf, 1, sizegBuf, inf);
	cnt->eoBuffer = siz;
	cnt->eofPosInBuffer = 0;
	
	/* Check if we're running into the end of the file. */
	if (feof(inf)) {
	  cout << "Reached end of mode3 (pre global-header) file." << endl;
	  cnt->eofInBuffer = 1;
	  cnt->eofPosInBuffer = siz;
	  cout << "EOF position in this buffer = " << cnt->eofPosInBuffer << endl;
	}
	
	for (j=0; j<siz; j+=2) {
	  swap(*(gBuf + j), *(gBuf + j + 1));
	}

	cnt->mode3i = 0;
	
      } else if (remaining > 0) { } /* Still more data in the buffer, just 
				       get the next packet... */
      
      /* Okay, we've got data in the buffer, let's go... */
      
      unsigned char *tmp = (gBuf + cnt->mode3i*2);
      
      /* Memory allocation for data packet... */
      memset(dp->waveform, 1, gMode3.tracelength*sizeof(UShort_t));
      
      /* Check there is enough data in the buffer to get the 'AAAA header'...
	 the buffer was arbitrary length, so there might not be. */
      if ( (cnt->eoBuffer - cnt->mode3i*2) < ((sizeof(dp->aahdr) + sizeof(dp->hdr))/2) && 
	   ( (cnt->eofInBuffer != 1) ) ) { 

	if (debugAdd) {
	  cout << "Failed at sufficient data for 'AAAA' -- getting more. " << endl;
	}

	/* Not enough data left in this buffer...get more. */
	int remainingBytes = (cnt->eoBuffer - cnt->mode3i*2);

	/* Move remaining data from the end of the buffer, to the top. */
	tmp = (gBuf + cnt->mode3i*2);
	memmove(gBuf, tmp, remainingBytes);

	for (j=0; j<remainingBytes+10; j=j+2) {
	  swap(*(gBuf + j), *(gBuf + j + 1));
	}
	
	siz = fread(gBuf + remainingBytes, 1, sizegBuf-(remainingBytes), inf);
	cnt->eoBuffer = siz + remainingBytes;
	
	if (feof(inf)) {
	  cout << "Reached end of mode3 (pre global-header) file." << endl;
	  cnt->eofInBuffer = 1;
	  cnt->eofPosInBuffer = siz;
	  cout << " EOF position in this buffer = " << cnt->eofPosInBuffer << endl;
	}

	for (j=0; j<(siz + remainingBytes); j=j+2) {
	  swap(*(gBuf + j), *(gBuf + j + 1));
	}
	
	cnt->mode3i = 0;
	
      } else if ( (cnt->eofInBuffer !=1) || 
		  (cnt->eofInBuffer == 1 && 
		   cnt->eofPosInBuffer >= cnt->mode3i*2 + ((sizeof(dp->aahdr) + sizeof(dp->hdr))) ) ) {
	if (debugAdd) { cout << "Enough data to get 'AAAA'" << endl; }
      }
      
      /* We have data, go ahead and get header information. */
      tmp = (gBuf + cnt->mode3i*2);
      memmove(&dp->aahdr[0], tmp, (sizeof(dp->hdr) + sizeof(dp->aahdr)));

      if ((dp->aahdr[0] != 43690) || (dp->aahdr[1] != 43690)) {
	cerr << "Didn't find 'AAAA' header as expected!" << endl;
	cout << hex << "aahdr[0] " << dp->aahdr[0] << " [1] " 
	     << dp->aahdr[1] << dec << endl;
      } else {
	if (debugAdd) { cout << "  Got 'AAAA' header okay" << endl; }
      }
      
      cnt->mode3i += (sizeof(dp->aahdr) + sizeof(dp->hdr))/2;
      cnt->old3Bytes += (sizeof(dp->aahdr) + sizeof(dp->hdr))/2;
      cnt->bytes_read += (sizeof(dp->aahdr) + sizeof(dp->hdr))/2;
      cnt->bytes_read_since_last_time += (sizeof(dp->aahdr) + sizeof(dp->hdr))/2;
      
      gMode3.tracelength = (((dp->hdr[0] & 0x7ff))*2) - 14;
      tmp = (gBuf + cnt->mode3i*2);

      if (debugAdd) {
	cout << "tracelength " << gMode3.tracelength << endl;
      }
      
      /* Check there is enough data in the buffer to get the waveform...
	 the buffer was arbitrary length, so there might not be. */
      if ( (cnt->eoBuffer - cnt->mode3i*2) < gMode3.tracelength*sizeof(UShort_t) &&
	   ( (cnt->eofInBuffer != 1) ) ) {

	if (debugAdd) {
	  cout << "Failed at sufficient data for waveform -- getting more. " << endl;
	  cout << " cnt->eoBuffer " << cnt->eoBuffer << ", cnt->mode3i*2 " << cnt->mode3i*2 << endl;
	}
	
	/* Not enough data left in this buffer...get more. */
	int remainingBytes = (cnt->eoBuffer - cnt->mode3i*2);

	/* Move remaining data from the end of the buffer, to the top. */
	tmp = (gBuf + cnt->mode3i*2);
	memmove(gBuf, tmp, remainingBytes);

	for (j=0; j<remainingBytes+10; j=j+2) {
	  swap(*(gBuf + j), *(gBuf + j + 1));
	}
	
	siz = fread(gBuf + remainingBytes, 1, sizegBuf-(remainingBytes), inf);
	cnt->eoBuffer = siz + remainingBytes;
	
	if (feof(inf)) {
	  cout << "Reached end of mode3 (pre global-header) file." << endl;
	  cnt->eofInBuffer = 1;
	  cnt->eofPosInBuffer = siz;
	  cout << " EOF position in this buffer = " << cnt->eofPosInBuffer << endl;
	}
	
	for (j=0; j<(siz + remainingBytes); j=j+2) {
	  swap(*(gBuf + j), *(gBuf + j + 1));
	}
	
	cnt->mode3i = 0;

      } else if ( (cnt->eofInBuffer != 1) || 
		  (cnt->eofInBuffer == 1 && 
		   cnt->eofPosInBuffer >= (cnt->mode3i*2 + gMode3.tracelength*sizeof(UShort_t) ) ) ) {
      }
      
      /* Now we have more data, go ahead and get waveform. */
      tmp = (gBuf + cnt->mode3i*2);
      memmove(&dp->waveform[0], tmp, gMode3.tracelength*sizeof(UShort_t));
      
      
      cnt->mode3i += (gMode3.tracelength*sizeof(UShort_t))/2;
      cnt->old3Bytes += (gMode3.tracelength*sizeof(UShort_t))/2;
      cnt->bytes_read += (gMode3.tracelength*sizeof(UShort_t))/2;
      cnt->bytes_read_since_last_time += (gMode3.tracelength*sizeof(UShort_t))/2;
      tmp = (gBuf + cnt->mode3i*2); 
    }

    /* One way or another, we got a mode3 channel event, now interpret it. */

    /* Interpret the header information... */
    gMode3.led_timestamp = 0;
    gMode3.led_timestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[3])) +
					((ULong64_t)(dp->hdr[2]) << 16) +
					((ULong64_t)(dp->hdr[5]) << 32) );
    
    /* Now move data to another buffer for output with a global header. */
    gHeader.length = (sizeof(dp->hdr) + sizeof(dp->aahdr) + 
		      gMode3.tracelength*sizeof(UShort_t));  
    unsigned char *tmpOut = (gOutBuf);
    memcpy(tmpOut, &dp->aahdr[0], (sizeof(dp->hdr) + sizeof(dp->aahdr)));
    tmpOut = (gOutBuf + (sizeof(dp->hdr) + sizeof(dp->aahdr)));
    memcpy(tmpOut, &dp->waveform[0], gMode3.tracelength*sizeof(UShort_t));
    
    for (j=0; j<gHeader.length; j=j+2) {
      swap(*(gOutBuf + j), *(gOutBuf + j + 1));
    }
    
    /* And write out with new global header */
    siz = fwrite(&gHeader, 1, sizeof(struct globalHeader), outfile);
    siz = fwrite(&gOutBuf, 1, (sizeof(dp->hdr) + sizeof(dp->aahdr) + 
			       gMode3.tracelength*sizeof(UShort_t)), outfile);
    
    free(dp);
    
    /* No global header...will we need more data next time?  
       Or are we all done? */
    
    if ( (UInt_t)(cnt->mode3i*2) == cnt->eoBuffer && cnt->eofInBuffer ) { gotsignal = 1; }
    if ( (UInt_t)(cnt->mode3i*2) == cnt->eoBuffer ) { return 0; }
    else if ( (UInt_t)(cnt->mode3i*2) < cnt->eoBuffer ) { return 1; }
    
  } else if (type == DECOMP) { /* Mode 2 data */

    if (ctrl->pgh) { /* ctrl->pgh = pre global-header */
      
      /* We're going to assume if the data is old, it's probably
	 a really old mode 2 data format too :) */
      siz = fread(&gMode2reallyOld, 1, sizeof(struct mode2StructOld),  inf);
      if (siz != sizeof(struct mode2StructOld)) {
      	cerr << "GetGRETINA(MODE2) failed in bytes read (ctrl->pgh)." << endl;
      } else {
	cnt->bytes_read += sizeof(struct mode2StructOld);
	cnt->bytes_read_since_last_time += sizeof(struct mode2StructOld);
      }
      
      /* To avoid problems with multiple mode2Event containers, I'm going
	 to repack all old mode2 data formats into the newest.  Variables that
	 don't exist in old formats will be set to -282. Just because. */
      
      gMode2.type = 0xabcd5678;
      gMode2.crystal_id = gMode2reallyOld.crystal_id;
      gMode2.num = gMode2reallyOld.num;
      gMode2.tot_e = gMode2reallyOld.tot_e;
      for (Int_t i=0; i<4; i++) {
	gMode2.core_e[i] = -282;
      }
      gMode2.timestamp = gMode2reallyOld.timestamp;
      gMode2.trig_time = -282;
      gMode2.t0 = gMode2reallyOld.t0;
      gMode2.cfd = -282;
      gMode2.chisq = gMode2reallyOld.chisq;
      gMode2.norm_chisq = gMode2reallyOld.norm_chisq;
      gMode2.baseline = -282;
      gMode2.prestep = -282; gMode2.poststep = -282;
      gMode2.pad = -282;
      for (Int_t i=0; i<MAX_INTPTS; i++) {
	gMode2.intpts[i].x = gMode2reallyOld.intpts[i].x;
	gMode2.intpts[i].y = gMode2reallyOld.intpts[i].y;
	gMode2.intpts[i].z = gMode2reallyOld.intpts[i].z;
	gMode2.intpts[i].e = gMode2reallyOld.intpts[i].e;
	gMode2.intpts[i].seg = -282;
	gMode2.intpts[i].seg_energy = -282;
      }

      if ((UInt_t)siz < sizeof(struct mode2StructOld)) {
	gotsignal = 1;
      }

      /* And repack for writing out with new global header */
      if (!gotsignal) {
	gHeader.type = 1;
	gHeader.timestamp = gMode2.timestamp;
	gHeader.length = sizeof(struct mode2StructNew);
	siz = fwrite(&gHeader, 1, sizeof(struct globalHeader), outfile);
	siz = fwrite(&gMode2, 1, sizeof(struct mode2StructNew), outfile);
      }
      
      return 0;

    }
    
  } else if (type == TRACK) {
  
    if (ctrl->pgh) { /* ctrl->pgh = pre global-headers */
      cerr << "Sorry, there is no option for analysis of pre global-header mode1 data. " << endl;
      cerr << " (I don't think mode1 was very well-defined at that time?) " << endl;
    }
  }
  
  return 0;
}

/****************************************************/

void PrintHelpInformation() {
  printf("\n");
  printf("Usage: ./AddGlobalHeaders <File Type> <Usage flags> -d <Subdirectory, i.e. CR-5> -run <input run ### (separate multiple files by a space)>\n");
  printf("    Valid file types: -g  (Global.dat)\n");
  printf("                      -f  <FILENAME> (specific file -- do not require subdirectory or run #)\n");
  printf("    Valid usage flags: -preGH2 (file is Mode2, from the pre global header (pGH) era) \n");
  printf("                       -preGH3 (file is Mode3, from the pre global header (pGH) era)\n");
  printf("                       -outputName <OUTPUT FILENAME> (write output file with user defined name)\n");
  printf("                       -zip (compressed data file, will append .gz to filename)\n");
  printf("                       -bzip (compressed data file, will append .bz2 to filename)\n");
  printf("                       -suppressTS (suppress TS error warnings in analysis)\n\n");

}
