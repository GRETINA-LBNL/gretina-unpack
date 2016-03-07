#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <signal.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <time.h>

#include "Structures.h"

FILE *mergeOutFile;

Int_t gotsignal;
void breakhandler(int dummy) {
  gotsignal = 1; 
}

/**************************************/
int main(int argc, char *argv[]) {

  gotsignal = 0;
  signal(SIGINT, breakhandler);
  signal(SIGPIPE, breakhandler);

  TString inputFile1, inputFile2, outputFile;
  Bool_t haveFile1=0, haveFile2=0, haveOutput=0;
  Int_t file1zipped = -1; Int_t file2zipped = -1;
  Int_t file1Events=0, file2Events=0;
  Int_t dataWritten1=0, dataWritten2=0;

  TFile *mOut = new TFile("merge.root", "RECREATE");
  TH1F *tac = new TH1F("tac", "tac", 1000, -500, 500);

  if (argc < 7) {
    printf("\nMinimum usage: %s -f1 <Filename1> -f2 <Filename2> -fOut <OutputFilename>\n", argv[0]);
    exit(1);
  } else if ((strcmp(argv[1], "-h") == 0)) {
    printf("\nMinimum usage: %s -f1 <Filename1> -f2 <Filename2> -fOut <OutputFilename>\n", argv[0]);
    exit(1);
  } else {
    Int_t i=1;
    while( i<argc ) {
      if (strcmp(argv[i], "-f1") == 0) {
	haveFile1 = 1; i++;
	inputFile1 = argv[i]; i++;
	file1zipped = 0;
      } else if (strcmp(argv[i], "-f1zip") == 0) {
	haveFile1 = 1; i++;
	inputFile1 = argv[i]; i++;
	file1zipped = 1;
      } else if (strcmp(argv[i], "-f2") == 0) {
	haveFile2 = 1; i++;
	inputFile2 = argv[i]; i++;
	file2zipped = 0;
      } else if (strcmp(argv[i], "-f2zip") == 0) {
	haveFile2 = 1; i++;
	inputFile2 = argv[i]; i++;
	file2zipped = 1;
      } else if (strcmp(argv[i], "-fOut") == 0) {
	haveOutput = 1; i++;
	outputFile = argv[i]; i++;
      } else if (strcmp(argv[i], "-f1ExcludeCrysID") == 0) {
	i++;
	nFile1Exclude = atoi(argv[i]); i++;
	for (Int_t n=0; n<nFile1Exclude; n++) {
	  file1Exclude[n] = atoi(argv[i]); i++; 
	}
      } else if (strcmp(argv[i], "-h") == 0) {
	 printf("\nMinimum usage: %s -f1 <Filename1> -f2 <Filename2> -fOut <OutputFilename>\n", argv[0]);
	 printf(" Option flags: -f1ExcludeCrysID <N> <Value1> ... <ValueN>\n");
	 printf("                     This excludes certain mode2 crystal_id values in File1 from\n");
	 printf("                     being included in the merged output.\n");
	 printf("\n");
	 exit(1);
      }
    }
  }

  TStopwatch timer;
  timer.Start();
  
  if ( !haveFile1 || !haveFile2 || !haveOutput ) {
    cout << "Problem: we don't have all the file information required." << endl;
    printf("Minimum usage: %s -f1 <Filename1> -f2 <Filename2> -fOut <OutputFilename>", argv[0]);
    exit(1);
  }

  /***********************************************************/
  /* MERGED OUTPUT FILE NAME                                 */
  /***********************************************************/

  int pipeFlag = 0;

  /* Set-up merged output file */
  FILE *mergeOutFile;
  if (!strcmp(outputFile.Data(), "pipe")) {
    mergeOutFile = stdout;
    pipeFlag = 1;
  } else {
    mergeOutFile = fopen64(outputFile.Data(), "w");
  }
  if (!mergeOutFile) {
    printf("Cannot open \"%s\" for reading\n", outputFile.Data());
    exit(3);
  } else {
    if (!pipeFlag) 
      printf("Output file \"%s\" opened\n", outputFile.Data());
  }

  /***********************************************************/
  /* INPUT FILE #1 FOR MERGING                               */
  /***********************************************************/

  if (file1zipped) {  
    inputFile1 = "./GEB_HFC -z -p " + inputFile1;
  } else {
    inputFile1 = "./GEB_HFC -p " + inputFile1;
  }
  inFile1 = popen(inputFile1.Data(), "r");
  if (!inFile1) {
    printf("Cannot open %s for reading\n", inputFile1.Data());
    exit(2);
  } else {
    if (!pipeFlag)
      printf("File #1 -- \"%s\" opened for reading\n", inputFile1.Data());
  }

  /***********************************************************/
  /* INPUT FILE #2 FOR MERGING                               */
  /***********************************************************/

  if (file2zipped) {
    inputFile2 = "./GEB_HFC -z -p " + inputFile2;
  } else {
    inputFile2 = "./GEB_HFC -p " + inputFile2;
  }
  inFile2 = popen(inputFile2.Data(), "r");
  if (!inFile2) {
    printf("Cannot open %s for reading\n", inputFile2.Data());
    exit(2);
  } else {
    if (!pipeFlag) 
      printf("File #2 -- \"%s\" opened for reading\n", inputFile2.Data());
  }
 
  /***********************************************************/
  /* Loop and merge files...                                 */
  /***********************************************************/

long long int lastTS = 0;  int lastTSfile = 0;

  while ( !doneFile2  && ! gotsignal) {
 
    GetNextEvent2(dataWritten2);
   
    if ( !doneFile1 ) {
      while ( GetNextEvent1(dataWritten1) && !doneFile1  && !gotsignal) {

	if ((fileBuf1.header.timestamp) < 
	    (fileBuf2.header.timestamp)) {

//cout << fileBuf1.header.timestamp << " " << fileBuf2.header.timestamp << " -- 1 lower " << endl;

if (fileBuf1.header.timestamp - lastTS < 1000) {
if (lastTSfile == 2) {tac->Fill(fileBuf1.header.timestamp - lastTS);}
} else {
   lastTS = fileBuf1.header.timestamp;  lastTSfile = 1;
}

	  /* Write the file 1 event to disk, it's go the lower TS. */
	  fwrite(&fileBuf1, 
		 fileBuf1.header.length + sizeof(struct globalHeader),
		 1, mergeOutFile);
	  file1Events++;
	  dataWritten1 = 1;
	  dataWritten2 = 0;
	} else {

//cout << fileBuf1.header.timestamp << " " << fileBuf2.header.timestamp << " -- 2 lower " << endl;
//cin.get();
	  /* Write the file 2 event to disk, it's go the lower TS. */
	  fwrite(&fileBuf2, 
		 fileBuf2.header.length + sizeof(struct globalHeader),
		 1, mergeOutFile);	  
	  file2Events++;
	  dataWritten1 = 0;
	  dataWritten2 = 1;
	  break;
	}
     
      }
	
    } else { /* if (doneFile1) */
      /* File 1 is done, just empty out file 2 now */
      fwrite(&fileBuf2, fileBuf2.header.length + sizeof(struct globalHeader),
	     1, mergeOutFile);
      file2Events++;    
      dataWritten2 = 1; 
    }

  }

  fclose(inFile2);
  
  if (!doneFile1) {
    /* Out of file2 events, dump the remaining file1 events, if any */
    while ( GetNextEvent1(dataWritten1) && !doneFile1) {
      fwrite(&fileBuf1, fileBuf1.header.length + sizeof(struct globalHeader), 
	     1, mergeOutFile);
      file1Events++;
      dataWritten1 = 1;
    }
  }
  
  if (!pipeFlag) 
    cout << "File 1 events " << file1Events << " " 
	 << "File 2 events " << file2Events << endl;
    
  fclose(mergeOutFile);

	tac->Write();
mOut->Write();
mOut->Close();
  return 0;
}

