#include "Riostream.h"
#include "Rtypes.h"
#include <iomanip>
#include <vector>
#include <deque>
#include <assert.h>
#include <cmath>
#include <iostream>

#include "TObject.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSortedList.h"
#include "TStopwatch.h"

#include "RdGeCalFile.c"

/******** Definitions ********/

#define DECOMP 1
#define RAW    2
#define TRACK  3
#define BGS    4
#define S800   5
#define NSCLNONEVENT 6
#define GRETSCALER   7
#define BANK29 8
#define S800PHYSICS  9
#define TSNSCLNONEVENT 10
#define G4SIMULATION 11
#define CHICORAW     12
#define MYSTERY      13
#define DGS    14
#define DGSTRIG 15
#define DFMA 16
#define PHOSWALL 17
#define PHOSWALLAUX 18

#define maxGRETINAwords 4096*32*32
#define MAX_INTPTS 16

/******** Data structures ********/

struct globalHeader {
  int type;
  int length;
  long long timestamp;
};

struct fileData {
  globalHeader header;
  unsigned short data[maxGRETINAwords];
};

struct ip {
  float x, y, z, e;
  int seg;
  float seg_energy;
};

struct mode2StructABCD1234 {
  int type;
  int crystal_id;
  int num;
  float tot_e;
  long long int timestamp;
  long long trig_time;
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  int pad;
  ip intpts[MAX_INTPTS];
};

struct mode2StructABCD5678 {
  int type;
  int crystal_id;
  int num;
  float tot_e;
  long long int timestamp;
  long long trig_time;
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  int pad;
  ip intpts[MAX_INTPTS];
};

/******** GRETINA Variables ********/

FILE *inFile1, *inFile2;

fileData fileBuf1;
fileData fileBuf2;
fileData fileTemp;

mode2StructABCD5678 g2;

Bool_t doneFile1=0, doneFile2=0;
Int_t dataRead1 = 0, dataRead2 = 0;
long long int ktn = 0, ktn2 = 0;
int MBread1 = 0, MBread2 = 0;

long long int tsLAST;

Int_t nFile1Exclude = 0;
Int_t file1Exclude[2048] = {-1};

/******** GRETINA Functions ********/

int GetNextEvent1(int written) {
  if (written == 1) {

    Int_t excludeThis = 1;
    
    while (excludeThis == 1) {
      
      if ( (dataRead1 = fread(&fileBuf1.header, sizeof(struct globalHeader), 
			      1, inFile1)) ) {

	excludeThis = 0;
	ktn += dataRead1*sizeof(struct globalHeader);
	// cout << " Type " << fileBuf1.header.type
	//      << " length " << fileBuf1.header.length
	//      << " TS " << fileBuf1.header.timestamp << endl;
	
	dataRead1 = fread(&fileBuf1.data, fileBuf1.header.length, 
			  1, inFile1);
	
	ktn += dataRead1*fileBuf1.header.length;
	
	if (ktn > 1024*1024*1024) {
	  MBread1++;	 
	  cout << "Read " << MBread1 << " GB from file #1" << endl;
	  ktn -= 1024*1024*1024;
	}
	
	if (dataRead1 != 1 && fileBuf1.header.length > 0) {
	  cerr << "Mismatch between event length and read - file #1." << endl;
	  doneFile1 = 1;
	  return 1;
	}
	
	/* Now, check if this is data we want to use for the merging. 
	   Here, we exclude any Mode2 crystal_id values dictated in the
	   command line input information. */
	if (fileBuf1.header.type == DECOMP) { /* Mode 2 data */
	  excludeThis = 0;
	  memcpy(&g2, &fileBuf1.data, sizeof(struct mode2StructABCD5678));
	  
	  for (Int_t j=0; j<nFile1Exclude; j++) {
	    if (g2.crystal_id == file1Exclude[j]) {
	      excludeThis = 1;
	    }
	  }
	}
      } else {
	doneFile1 = 1;
	return 0;
      }
    }
    
    doneFile1 = 0;
    return 1;   
    
  } else {
    return 1;
  }
}

int GetNextEvent2(int written) {
  if (written == 1) {

    Int_t excludeThis = 1;
    
    while (excludeThis == 1) {
      
      if ( (dataRead2 = fread(&fileBuf2.header, sizeof(struct globalHeader), 
			      1, inFile2)) ) {

	excludeThis = 0;
	ktn2 += dataRead2*sizeof(struct globalHeader);
	// cout << " Type " << fileBuf1.header.type
	//      << " length " << fileBuf1.header.length
	//      << " TS " << fileBuf1.header.timestamp << endl;
	
	dataRead2 = fread(&fileBuf2.data, fileBuf2.header.length, 
			  1, inFile2);
	
	ktn2 += dataRead2*fileBuf2.header.length;
	
	if (ktn2 > 1024*1024*1024) {
	  MBread2++;	 
	  cout << "Read " << MBread2 << " GB from file #2" << endl;
	  ktn2 -= 1024*1024*1024;
	}
	
	if (dataRead2 != 1 && fileBuf2.header.length > 0) {
	  cerr << "Mismatch between event length and read - file #2." << endl;
	  doneFile2 = 1;
	  return 1;
	}
	
      } else {
	doneFile2 = 1;
	return 0;
      }
    }
    
    doneFile2 = 0;
    return 1;   
    
  } else {
    return 1;
  }
}


int GetNextEvent2Old(int written) {
  if (written == 1) {
    if ( (dataRead2 = fread(&fileBuf2.header, sizeof(struct globalHeader), 
			    1, inFile2)) ) {
      ktn2 += dataRead2*sizeof(struct globalHeader);
      // cout << " Type " << fileBuf2.header.type
      //      << " length " << fileBuf2.header.length
      //      << " TS " << fileBuf2.header.timestamp << endl;
      

      dataRead2 = fread(&fileBuf2.data, fileBuf2.header.length, 
			1, inFile2);
      
      ktn2 += dataRead2*fileBuf2.header.length;
      
	if (ktn2 > 1024*1024*1024) {
	  MBread2++;	 
	  cout << "Read " << MBread2 << " GB from file #2" << endl;
	  ktn2 -= 1024*1024*1024;
	}

	if (dataRead2 != 1) {
	  cerr << "Mismatch between GRETINA event length and read. " << endl;
	  doneFile2 = 1;
	  return 1;
	}

	doneFile2 = 0;
	return 1;   
    } else {
      doneFile2 = 1;
      return 0;
    }
  } else {
    return 1;
  }
}

