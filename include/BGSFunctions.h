#ifndef __BGSFUNCTIONS_H
#define __BGSFUNCTIONS_H

/* BGS Analysis header file */

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

#include "Unpack.h"
#include "Structures.h"

#include "BGSParamList.h"
#include "BGSConditions.h"
#include "BGSFPPositions.h"
#include "BGSThresholds.h"

void InitializeBGSSettings(TString filename) {
  FILE *input; 
  if ( (input = fopen(filename.Data(), "r")) == NULL ) {
    cerr << "BGS settings file " << filename.Data() << " could not be opened. " << endl;
  }
  cout << "BGS settings based on " << filename.Data() << endl;
  
  char line[300];
  char junk[300];
  char* retVal;

  while (!feof(input)) {
    retVal = fgets(line, 300, in);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "Calibration", 11) == 0) { 
      
    }


}

void ProcessBGS();
void DoStripSearch();

void ProcessBGS() { /* Main BGS analysis function for the BGS */

  /* Here we process the BGS data. */
  for (int i=0; i<maxBGSwords; i++) {
    data[i] = 0;
  }

  BGSevent->ClearBGS();

  /* Skip any empty events */
  if (gHeader.length/sizeof(struct datavalue) == 0) { return; }
  
  /* Sort out the BGS data into a nice big easy to manipulate array. */
  for (UInt_t i=0; i<(gHeader.length/sizeof(struct datavalue)); i++) {
    data[(int)rawBGS.data[i].param] = rawBGS.data[i].value; 
#ifdef WITH_HISTOS
    bgsparhit->Fill((int)rawBGS.data[i].param);
#endif
  }
  
  /* Check error words, skip events with flags... */
  if (data[ierr1] || data[ierr2]) { return; }

  /* Calculate times from the different BGS scalers... */

  /* 10 MHz pulser... */
  BGSevent->Time10MHz = ((double)((unsigned long)data[iusec] + 
 				  ((unsigned long)data[iUSEC])*65536)/10000000.
			 + (double)(BGSevent->Overflows10MHz*429.4967296));

  while (BGSevent->Time10MHz < BGSevent->Last10MHzTime) { // Check for overflow
    BGSevent->Time10MHz += 429.4967296; // Adjust for overflow, and increment
    BGSevent->Overflows10MHz++;
  }
  BGSevent->Last10MHzTime = BGSevent->Time10MHz;

  /* GRETINA 50 MHz clock (just use from header, then don't need
     to count the overflows for this clock) */
  BGSevent->timestamp = gHeader.timestamp;
  BGSevent->time = BGSevent->timestamp*0.00000001; // Calculate time in seconds

  /* Do the hit search and calibrate... */
  DoStripSearch();

  /* Get the gamma information for the signals plugged into the BGS
     electronics -- param 352 - 368. */
  for (int m=0; m<16; m++) {
    BGSevent->egam[m] = data[352 + m];
  }

  /* Only sort event with valid pixel IDs... */
  if ( (BGSevent->ORchip == -1) || (BGSevent->ORfront == -1) ||
       (BGSevent->ORback == -1) ) {
    return;
  }
  
  /* This is a silly correction to the energies based on 
     the observed positions in the 176Pt alpha peak. */
  if ( (BGSevent->ORchip == 0) && (BGSevent->ALenergy > 0) && 
       (BGSevent->TEchip == -1) ) {
    BGSevent->ALenergy += 47.;
  } else if ( (BGSevent->ORchip == 1) && (BGSevent->ALenergy > 0) && 
	      (BGSevent->TEchip == -1) ) {
    BGSevent->ALenergy += 95.;
  } else if ( (BGSevent->ORchip == 2) && (BGSevent->ALenergy > 0) && 
	      (BGSevent->TEchip == -1) ) {
    BGSevent->ALenergy += 73.;
  } else if ( (BGSevent->ORchip == 3) && (BGSevent->ALenergy > 0) && 
	      (BGSevent->TEchip == -1) ) {
    BGSevent->ALenergy += 66.;
  }
  if ( (BGSevent->TEchip == 0) && (BGSevent->ALenergy > 0) ) {
    BGSevent->ALenergy += 47.;
  } else if ( (BGSevent->TEchip == 1) && (BGSevent->ALenergy > 0) ) {
    BGSevent->ALenergy += 95.;
  } else if ( (BGSevent->TEchip == 2) && (BGSevent->ALenergy > 0) ) {
    BGSevent->ALenergy += 73.;
  } else if ( (BGSevent->TEchip == 3) && (BGSevent->ALenergy > 0) ) {
    BGSevent->ALenergy += 66.;
  }

  /* There are no gamma energies at the FP right now... otherwise put
     analysis for these here. */

  /* Fill in MWPC data. */
  if (data[imae]>0 && data[imae]<4095) { BGSevent->eMWPC = data[imae]; }
  if (data[imat]>0 && data[imat]<4095) { BGSevent->tMWPC = data[imat]; }

  /* Define TOF condition */
#define tofCondition  (data[imat] >= mwa_fp_tdc_min && \
		       data[imat] <= mwa_fp_tdc_max)

  /* Define MWPC condition */
#define mwpcCondition  (data[imae] >= mwa_min && \
			data[imae] <= mwa_max)

  /* Define the recoil condition... good energy, good pixel
     ID and MWPC hit */
#define recoilCondition  (BGSevent->SIenergy >= r_mine &&  \
			  BGSevent->SIenergy <= r_maxe &&  \
			  BGSevent->ORchip > -1 &&   \
			  (BGSevent->ORfront > -1 && BGSevent->ORback > -1) && \
			  (tofCondition || mwpcCondition))

  /* Define the parent c.e. condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define peCondition  (BGSevent->CEenergy >= pe_mine &&  \
		      BGSevent->CEenergy <= pe_maxe &&  \
		      BGSevent->ORchip > -1 &&  \
		      BGSevent->ORfront > -1 &&  \
		      BGSevent->ORback > -1 &&  \
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the daughter c.e. condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define deCondition  (BGSevent->CEenergy >= de_mine &&	\
		      BGSevent->CEenergy <= de_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&   \
		      mwpcCondition == 0 &&    \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the granddaughter c.e. condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define geCondition  (BGSevent->CEenergy >= ge_mine &&	\
		      BGSevent->CEenergy <= ge_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the parent alpha condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define paCondition  (BGSevent->ALenergy >= pa_mine &&	\
		      BGSevent->ALenergy <= pa_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the daughter alpha condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define daCondition  (BGSevent->ALenergy >= da_mine &&	\
		      BGSevent->ALenergy <= da_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the granddaughter alpha condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define gaCondition  (BGSevent->ALenergy >= ga_mine &&	\
		      BGSevent->ALenergy <= ga_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the granddaughter alpha condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define gaCondition  (BGSevent->ALenergy >= ga_mine &&	\
		      BGSevent->ALenergy <= ga_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the parent escape condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define pxCondition  (BGSevent->SIenergy >= px_mine &&	\
		      BGSevent->SIenergy <= px_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the daughter escape condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define dxCondition  (BGSevent->SIenergy >= dx_mine &&	\
		      BGSevent->SIenergy <= dx_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the granddaughter escape condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define gxCondition  (BGSevent->SIenergy >= gx_mine &&	\
		      BGSevent->SIenergy <= gx_maxe &&	\
		      BGSevent->ORchip > -1 &&	\
		      BGSevent->ORfront > -1 &&	 \
		      BGSevent->ORback > -1 &&	\
		      tofCondition == 0 &&  \
		      mwpcCondition == 0 &&  \
		      evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  /* Define the fission condition... good energy,
     good pixel ID, no TOF, no MWPC, at least one recoil 
     in pixel so far */
#define fissionCondition  (BGSevent->FIenergy >= f_mine &&	\
			   BGSevent->FIenergy <= f_maxe &&	\
			   BGSevent->ORchip > -1 &&		\
			   BGSevent->ORfront > -1 &&		\
			   BGSevent->ORback > -1 &&		\
			   tofCondition == 0 &&			\
			   mwpcCondition == 0 &&			\
			   evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback] >= 1)

  if (recoilCondition) { BGSevent->EVR = 1; }
  else if (peCondition || deCondition || geCondition) { BGSevent->CONV = 1; }
  else if (paCondition || daCondition || gaCondition) { BGSevent->ALPHA = 1; }
  else if (pxCondition || dxCondition || gxCondition) { BGSevent->ESCAPE = 1; }
  else if (fissionCondition) { BGSevent->FISSION = 1; }

  /* Correlation analysis now... */
  /* Recoil events first. */
  if (recoilCondition) {
    evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].implanted = 1;
    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->ORenergy;
    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].TOF = BGSevent->TOF;
    // Gammas from GRETINA get added in the main analysis code after events are built.

    BGSevent->EVENTTYPE = R;
  }

  /* Parent C.E. events. */
  if (peCondition) {
    parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
    parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->CEenergy;
    parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
    parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
    /* If the c.e. event was reconstructed from two signals in the FP, the origin vs. terminus is
       ambiguous.  Put the event in the buffer again with origin and terminus reversed. */
    if (BGSevent->TEchip > -1) {
      parentCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]++;
      parentCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][parentCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].energy = BGSevent->CEenergy;
      parentCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][parentCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].timestamp = BGSevent->timestamp;
      parentCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][parentCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].time = BGSevent->time;
    }
  }

 if (deCondition) {
    daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
    daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->CEenergy;
    daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
    daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
    /* If the c.e. event was reconstructed from two signals in the FP, the origin vs. terminus is
       ambiguous.  Put the event in the buffer again with origin and terminus reversed. */
    if (BGSevent->TEchip > -1) {
      daughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]++;
      daughterCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][daughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].energy = BGSevent->CEenergy;
      daughterCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][daughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].timestamp = BGSevent->timestamp;
      daughterCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][daughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].time = BGSevent->time;
    }
  }

 if (geCondition) {
    granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
    granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->CEenergy;
    granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
    granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
    /* If the c.e. event was reconstructed from two signals in the FP, the origin vs. terminus is
       ambiguous.  Put the event in the buffer again with origin and terminus reversed. */
    if (BGSevent->TEchip > -1) {
      granddaughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]++;
      granddaughterCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][granddaughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].energy = BGSevent->CEenergy;
      granddaughterCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][granddaughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].timestamp = BGSevent->timestamp;
      granddaughterCECorr[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback][granddaughterCENum[BGSevent->TEchip][BGSevent->TEfront][BGSevent->TEback]%10].time = BGSevent->time;
    }
  }
 
 if (paCondition) {
   parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->ALenergy;
   parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }

 if (daCondition) {
   daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->ALenergy;
   daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }

 if (gaCondition) {
   granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->ALenergy;
   granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }
  
 if (pxCondition) {
   parentXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   parentXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->SIenergy;
   parentXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   parentXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][parentXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }

 if (dxCondition) {
   daughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   daughterXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->SIenergy;
   daughterXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   daughterXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][daughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }
 
 if (gxCondition) {
   granddaughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   granddaughterXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->SIenergy;
   granddaughterXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   granddaughterXCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][granddaughterXNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }

 if (fissionCondition) {
   fissionNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]++;
   fissionCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][fissionNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].energy = BGSevent->FIenergy;
   fissionCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][fissionNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].timestamp = BGSevent->timestamp;
   fissionCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][fissionNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]%10].time = BGSevent->time;
 }

 double deltaT1 = 0, deltaT2 = 0, deltaT3 = 0;

 /* Now the actual correlation search! */
 /* Recoil - parent C.E. search */
 if (peCondition) {
   for (int r=0; r<10; r++) {
     deltaT1 = (double)(parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
     if (deltaT1 > pe_maxt) { break; }
     if (deltaT1 >= pe_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

       BGSevent->EVENTTYPE = RE;
       
       BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
       BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
       BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
       BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
       if (BGSevent->ievrGamma.vgam.size() == 0) {
	 for (int i=0; i<100; i++) {
	   if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
	     gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
	     gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);
	     gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
	     BGSevent->ievrGamma.vgam.push_back(gammaTemp);
	   }
	 }
       }
       
       BGSevent->dparentCE.time = parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
       BGSevent->dparentCE.timestamp = parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
       BGSevent->dparentCE.energy = parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;

     } // End of time checks
   } // End of search over implants
 } // End of peCondition

 /* Recoil - parent C.E. - daughter C.E. search */
 if (deCondition) {
   for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
     if (deltaT1 > de_maxt) { break; }
     if (deltaT1 >= de_mint && parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */

       for (int r=0; r<10; r++) { // Step through recoil buffer to find correlation. 
	 deltaT2 = (double)(parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
			    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	 if (deltaT2 > pe_maxt) { break; }
	 if (deltaT2 >= pe_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

	   BGSevent->EVENTTYPE = REE;
	  	 	   
	   BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	   BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	   BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	   BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	   if (BGSevent->ievrGamma.vgam.size() == 0) {
	     for (int i=0; i<100; i++) {
	       if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		 gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		 gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);	      
		 gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		 BGSevent->ievrGamma.vgam.push_back(gammaTemp);
	       }
	     }
	   }
	   
	   BGSevent->dparentCE.time = 
	     parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	   BGSevent->dparentCE.timestamp = 
	     parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	   BGSevent->dparentCE.energy = 
	     parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;
	   
	   BGSevent->ddaughCE.time = 
	     daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	   BGSevent->ddaughCE.timestamp = 
	     daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	   BGSevent->ddaughCE.energy = 
	     daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;
	   
	 } // End of time checks (parent-recoil)
       } // End of search over implants
     } // End of time checks (daughter-parent)
   } // End of search over parent decays
 } // End of deCondition
  
 /* Recoil - parent C.E. - daughter C.E. - granddaughter C.E. search */
 if (geCondition) {
   for (int d=0; d<10; d++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time);
     if (deltaT1 > ge_maxt) { break; }
     if (deltaT1 >= ge_mint && daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d > 0) { /* Valid correlation! */

       for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
	 deltaT2 = (double)(daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time - 
			    parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
	 if (deltaT2 > de_maxt) { break; }
	 if (deltaT2 >= de_mint && parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */
	   
	   for (int r=0; r<10; r++) { // Step through recoil buffer to find correlation. 
	     deltaT3 = (double)(parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
				evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	     if (deltaT3 > pe_maxt) { break; }
	     if (deltaT3 >= pe_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

	       BGSevent->EVENTTYPE = REEE;
	       
	       BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	       BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	       BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	       BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	       if (BGSevent->ievrGamma.vgam.size() == 0) {	    
		 for (int i=0; i<100; i++) {
		   if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		     gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		     gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);		  
		     gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		     BGSevent->ievrGamma.vgam.push_back(gammaTemp);
		   }
		 }
	       }
	       
	       BGSevent->dparentCE.time = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	       BGSevent->dparentCE.timestamp = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	       BGSevent->dparentCE.energy = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;

	       BGSevent->ddaughCE.time = 
		 daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time;
	       BGSevent->ddaughCE.timestamp = 
		 daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].timestamp;
	       BGSevent->ddaughCE.energy = 
		 daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].energy;

	       BGSevent->dgdaughCE.time = 
		 granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	       BGSevent->dgdaughCE.timestamp = 
		 granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	       BGSevent->dgdaughCE.energy = 
		 granddaughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;
	       
	     } // End of time checks (parent-recoil)
	   } // End of search over implants
	 } // End of time checks (daughter-parent)
       } // End of search over parent decays
     } // End of time checks (granddaughter-daughter) 
   } // End of search over daughter decays
 } // End of geCondition

 /* Recoil - parent alpha search */
 if (paCondition) {
   for (int r=0; r<10; r++) {
     deltaT1 = (double)(parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
     if (deltaT1 > pa_maxt) { break; }
     if (deltaT1 >= pa_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

       BGSevent->EVENTTYPE = RA;
       
       BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
       BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
       BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
       BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
       if (BGSevent->ievrGamma.vgam.size() == 0) {
	 for (int i=0; i<100; i++) {
	   if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
	     gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
	     gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);	  
	     gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
	     BGSevent->ievrGamma.vgam.push_back(gammaTemp);
	   }
	 }
       }

       BGSevent->dparentAL.time = 
	 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
       BGSevent->dparentAL.timestamp = 
	 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
       BGSevent->dparentAL.energy = 
	 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;

     } // End of time checks
   } // End of search over implants

   if (BGSevent->ievr.energy>7000 && BGSevent->ievr.energy<13000) {
     if (BGSevent->dparentAL.energy>8000 && BGSevent->dparentAL.energy<8150) {
       if ((BGSevent->dparentAL.time - BGSevent->ievr.time)<300) {
	 // cout << BGSevent->ievr.timestamp/2 << "  " << BGSevent->dparentAL.timestamp/2 << endl;
       }
     }
   }

 } // End of paCondition

 /* Recoil - parent C.E. - parent alpha search */
 if (paCondition) {
   for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
     if (deltaT1 > pa_maxt) { break; }
     if (deltaT1 >= pa_mint && parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */
       for (int r=0; r<10; r++) {
	 deltaT2 = (double)(parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
			    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	 if (deltaT2 > pe_maxt) { break; }
	 if (deltaT2 >= pe_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */
	  
	   BGSevent->EVENTTYPE = REA;

	   BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	   BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	   BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	   BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	   if (BGSevent->ievrGamma.vgam.size() == 0) {	
	     for (int i=0; i<100; i++) {
	       if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		 gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		 gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);	      
		 gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		 BGSevent->ievrGamma.vgam.push_back(gammaTemp);
	       }
	     }
	   }

	   BGSevent->dparentCE.time = 
	     parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	   BGSevent->dparentCE.timestamp = 
	     parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	   BGSevent->dparentCE.energy = 
	     parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;

	   BGSevent->dparentAL.time = 
	     parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	   BGSevent->dparentAL.timestamp = 
	     parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	   BGSevent->dparentAL.energy = 
	     parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;
	   
	 } // End of time checks (recoil-parent CE)
       } // End of search over implants
     } // End of time checks (parent CE-parent AL)
   } // End of search over parent CE
 } // End of paCondition

 /* Recoil - parent C.E. - daughter C.E. - parent alpha search */
 if (paCondition) {
   for (int d=0; d<10; d++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time);
     if (deltaT1 > pa_maxt) { break; }
     if (deltaT1 >= pa_mint && daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d > 0) { /* Valid correlation! */
       
       for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
	 deltaT2 = (double)(daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time - 
			    parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
	 if (deltaT2 > de_maxt) { break; }
	 if (deltaT2 >= de_mint && parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */
	   
	   for (int r=0; r<10; r++) { // Step through recoil buffer to find correlation. 
	     deltaT3 = (double)(parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
				evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	     if (deltaT3 > pe_maxt) { break; }
	     if (deltaT3 >= pe_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

	       BGSevent->EVENTTYPE = REEA;
	       
	       BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	       BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	       BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	       BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	       if (BGSevent->ievrGamma.vgam.size() == 0) {	    
		 for (int i=0; i<100; i++) {
		   if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		     gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		     gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);	
		     gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		     BGSevent->ievrGamma.vgam.push_back(gammaTemp);
		   }
		 }
	       }

	       BGSevent->dparentCE.time = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	       BGSevent->dparentCE.timestamp = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	       BGSevent->dparentCE.energy = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;

	       BGSevent->ddaughCE.time = 
		 daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time;
	       BGSevent->ddaughCE.timestamp = 
		 daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].timestamp;
	       BGSevent->ddaughCE.energy = 
		 daughterCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].energy;

	       BGSevent->dparentAL.time = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	       BGSevent->dparentAL.timestamp = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	       BGSevent->dparentAL.energy = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;
	       
	     } // End of time checks (parent-recoil)
	   } // End of search over implants
	 } // End of time checks (daughter-parent)
       } // End of search over parent decays
     } // End of time checks (alpha-daughter CE) 
   } // End of search over daughter decays
 } // End of paCondition

 /* Recoil - parent alpha - daughter alpha search */
 if (daCondition) {
   for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
     if (deltaT1 > da_maxt) { break; }
     if (deltaT1 >= da_mint && parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */
       
       for (int r=0; r<10; r++) { // Step through recoil buffer to find correlation. 
	 deltaT2 = (double)(parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
			    evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	 if (deltaT2 > pa_maxt) { break; }
	 if (deltaT2 >= pa_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

	   BGSevent->EVENTTYPE = RAA;
	   
	   BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	   BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	   BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	   BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	   if (BGSevent->ievrGamma.vgam.size() == 0) {
	     for (int i=0; i<100; i++) {
	       if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		 gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		 gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);	     
		 gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		 BGSevent->ievrGamma.vgam.push_back(gammaTemp);
	       }
	     }
	   }

	   BGSevent->dparentAL.time = 
	     parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	   BGSevent->dparentAL.timestamp = 
	     parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	   BGSevent->dparentAL.energy = 
	     parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;

	   BGSevent->ddaughAL.time = 
	     daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	   BGSevent->ddaughAL.timestamp = 
	     daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	   BGSevent->ddaughAL.energy = 
	     daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;
       
	 } // End of time checks (parent-recoil)
       } // End of search over implants
     } // End of time checks (daughter-parent)
   } // End of search over parent decays
 } // End of daCondition

 /* Recoil - parent C.E. - parent alpha - daughter alpha search */
 if (daCondition) {
   for (int d=0; d<10; d++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time);
     if (deltaT1 > da_maxt) { break; }
     if (deltaT1 >= da_mint && parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d > 0) { /* Valid correlation! */

       for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
	 deltaT2 = (double)(parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time - 
			    parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
	 if (deltaT2 > pa_maxt) { break; }
	 if (deltaT2 >= pa_mint && parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */
	   
	   for (int r=0; r<10; r++) { // Step through recoil buffer to find correlation. 
	     deltaT3 = (double)(parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
				evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	     if (deltaT3 > pe_maxt) { break; }
	     if (deltaT3 >= pe_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

	       BGSevent->EVENTTYPE = REAA;
	       
	       BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	       BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	       BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	       BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	       if (BGSevent->ievrGamma.vgam.size() == 0) {	   
		 for (int i=0; i<100; i++) {
		   if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		     gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		     gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);		 
		     gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		     BGSevent->ievrGamma.vgam.push_back(gammaTemp);
		   }
		 }
	       }

	       BGSevent->dparentCE.time = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	       BGSevent->dparentCE.timestamp = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	       BGSevent->dparentCE.energy = 
		 parentCECorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentCENum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;

	       BGSevent->dparentAL.time = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time;
	       BGSevent->dparentAL.timestamp = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].timestamp;
	       BGSevent->dparentAL.energy = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].energy;
	       
	       BGSevent->ddaughAL.time = 
		 daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	       BGSevent->ddaughAL.timestamp = 
		 daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	       BGSevent->ddaughAL.energy = 
		 daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;

	     } // End of time checks (parent-recoil)
	   } // End of search over implants
	 } // End of time checks (parent CE-parent alpha)
       } // End of search over parent CE decays
     } // End of time checks (parent alpha-daughter alpha) 
   } // End of search over parent alpha decays
 } // End of daCondition

 /* Recoil - parent alpha - daughter alpha - granddaughter alpha search */
 if (gaCondition) {
   for (int d=0; d<10; d++) { // Step through parent buffer to find correlation. 
     deltaT1 = (double)(granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time - 
			daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time);
     if (deltaT1 > ga_maxt) { break; }
     if (deltaT1 >= ga_mint && daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d > 0) { /* Valid correlation! */
       
       for (int p=0; p<10; p++) { // Step through parent buffer to find correlation. 
	 deltaT2 = (double)(daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time - 
			    parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time);
	 if (deltaT2 > da_maxt) { break; }
	 if (deltaT2 >= da_mint && parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p > 0) { /* Valid correlation! */
	   
	   for (int r=0; r<10; r++) { // Step through recoil buffer to find correlation. 
	     deltaT3 = (double)(parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time - 
				evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time);
	     if (deltaT3 > pa_maxt) { break; }
	     if (deltaT3 >= pa_mint && evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r > 0) { /* Valid correlation! */

	       BGSevent->EVENTTYPE = RAAA;
	       
	       BGSevent->ievr.time = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].time;
	       BGSevent->ievr.timestamp = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].timestamp;
	       BGSevent->ievr.energy = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].energy;
	       BGSevent->ievr.TOF = evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].TOF;
	       if (BGSevent->ievrGamma.vgam.size() == 0) {	      
		 for (int i=0; i<100; i++) {
		   if (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i] > 0) {
		     gammaTemp.energy = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].egamma[i]);
		     gammaTemp.segmentSum = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].segegamma[i]);		 
		     gammaTemp.time = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].bgsTgamma[i]);
	     gammaTemp.fom = (evrCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(evrNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-r)%10].fomgamma[i]);
		     BGSevent->ievrGamma.vgam.push_back(gammaTemp);
		   }
		 }
	       }
	       
	       BGSevent->dparentAL.time = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].time;
	       BGSevent->dparentAL.timestamp = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].timestamp;
	       BGSevent->dparentAL.energy = 
		 parentALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(parentALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-p)%10].energy;
	       
	       BGSevent->ddaughAL.time = 
		 daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].time;
	       BGSevent->ddaughAL.timestamp = 
		 daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].timestamp;
	       BGSevent->ddaughAL.energy = 
		 daughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(daughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback]-d)%10].energy;

	       BGSevent->dgdaughAL.time = 
		 granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].time;
	       BGSevent->dgdaughAL.timestamp = 
		 granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].timestamp;
	       BGSevent->dgdaughAL.energy = 
		 granddaughterALCorr[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback][(granddaughterALNum[BGSevent->ORchip][BGSevent->ORfront][BGSevent->ORback])%10].energy;
	       
	     } // End of time checks (parent-recoil)
	   } // End of search over implants
	 } // End of time checks (daughter-parent)
       } // End of search over parent decays
     } // End of time checks (granddaughter-daughter) 
   } // End of search over daughter decays
 } // End of geCondition

 /* Come back to finish off the escapes and fission correlations later. */
 
}

void DoStripSearch() {

  /* Stolen from Ken...goes through the fronts of FP chips, remembers
     the two highest-energy hits, as long as they're in different hits. 
     Then two highest-energy hits on the back of FP chips, highest energy 
     hit from each of upstream (US) and punch-through (PT) strips.  Then
     does the reconstruction of fissions, alphas and conversion electrons. */
  
  int chip; // loop counter
  int st=0;   // strip counter
  int idxl, idxh; // 0th strip index for lo- and hi-energy signals
  double energy = 0; // energy in Si
  double elin = 0;   // energy from linear high-energy calibration
  double eexp = 0;   // energy from exponential calibration
  double chan = 0;  

  /* Step through the strips in the 3 FP fronts, and remember
     the largest two above-threshold hits. Must have a good prompt
     TDC signal. */
  for (chip=0; chip<3; chip++) { 

    if ( (chip==0 && ((data[iftt] > fttdc_min) && (data[iftt] < fttdc_max))) || 
	 (chip==1 && ((data[ifet] > fetdc_min) && (data[ifet] < fetdc_max))) ||
	 (chip==2 && ((data[ifwt] > fwtdc_min) && (data[ifwt] < fwtdc_max))) ) { 
      /* Here, assign start param numbers for given chip. */
      if (chip==0) { 
	idxl = iftl;
	idxh = ifth;
      } else if (chip==1) {
	idxl = ifel;
	idxh = ifeh;
      } else if (chip==2) {
	idxl = ifwl;
	idxh = ifwh;
      }

      for (st=0; st<nfp; st++) { // Step through strips

	if (data[idxl+st] > FPLMinChan[chip][st] || data[idxh+st]>0) {
	  
	  if (data[idxl+st] <= FPLMaxChan[chip][st]) { // Low energy branch
	    chan = (double)data[idxl+st] + (double)rand()/(double)RAND_MAX-0.5;
	    energy = (chan-FPLa[chip][st])*FPLm[chip][st];
	    BGSevent->fCal[chip][st] = energy;
	  } else { // High energy branch
	    chan = (double)data[idxh+st] + (double)rand()/(double)RAND_MAX-0.5;
	    elin = (chan-FPHa[chip][st])*FPHm[chip][st];
	    eexp = FPHc[chip][st]*exp(chan*FPHd[chip][st]);
	    if (elin >= eexp) { energy = elin; }
	    else { energy = eexp; }
	    BGSevent->fCal[chip][st] = energy;
	  }
	  
	  /* Now, remember the two largest hits from the front
	     of the focal plane chips. */
	  if (energy > BGSevent->Fenergy1) {
	    BGSevent->Fenergy2 = BGSevent->Fenergy1;
	    BGSevent->Fchip2 = BGSevent->Fchip1;
	    BGSevent->Fstrip2 = BGSevent->Fstrip1;
	    BGSevent->Fenergy1 = energy;
	    BGSevent->Fchip1 = chip;
	    BGSevent->Fstrip1 = st;
	  } else if (energy > BGSevent->Fenergy2) {
	    BGSevent->Fenergy2 = energy;
	    BGSevent->Fchip2 = chip;
	    BGSevent->Fstrip2 = st;
	  }
	} // Done with checking above threshold
      } // Done with stepping through strips
    } // Done with prompt TDC requirement
  } // Done with stepping through chips

  /* If there is a MWPC hit, calculate the TOF from MWPC to the FP */
  if (BGSevent->Fchip1==0 && data[imat]>=mwa_fp_tdc_min && data[imat]<=mwa_fp_tdc_max) {
    BGSevent->TOF = (data[iftt]-FPTa[0])*FPTm[0] - (data[imat]-FPTa[5])*FPTm[5];
  } else if (BGSevent->Fchip1==1 && data[imat]>=mwa_fp_tdc_min && data[imat]<=mwa_fp_tdc_max) {
    BGSevent->TOF = (data[iftt]-FPTa[1])*FPTm[1] - (data[imat]-FPTa[5])*FPTm[5];
  } else if (BGSevent->Fchip1==2 && data[imat]>=mwa_fp_tdc_min && data[imat]<=mwa_fp_tdc_max) {
    BGSevent->TOF = (data[iftt]-FPTa[2])*FPTm[2] - (data[imat]-FPTa[5])*FPTm[5];
  }

  /* For now, require the two highest-energy FP hits to be in
     different chips. */
  if (BGSevent->Fchip1 == BGSevent->Fchip2) {
    BGSevent->Fenergy2 = 0.;
    BGSevent->Fchip2 = -1;
    BGSevent->Fstrip2 = -1;
  }

  energy = 0;

  /* Step through the back strips and remember the largest two 
     above-threshold hits. */
  if ( (data[ifbt] > fbtdc_min) && (data[ifbt] < fbtdc_max) ) {
    /* Here, assign start param numbers. */
    idxl = ifbl;
    idxh = ifbh;
    chip = 3;
    
    for (st=0; st<nfp; st++) { // Step through strips
     
      energy = 0;

      /* This is a change because strips 28 and 29 are shorted together. */
      if (st==29) {
	data[idxh+29] = 0;
	data[idxl+29] = 0;
      }
      if (st==28) {
	data[idxh+28] += data[idxh+29];
	data[idxl+28] += data[idxl+29];
      }
 
      if (data[idxl+st] > FPLMinChan[chip][st] || data[idxh+st]>0) {
	
	if (data[idxl+st] <= FPLMaxChan[chip][st]) { // Low energy branch
	  chan = (double)data[idxl+st] + (double)rand()/(double)RAND_MAX-0.5;
	  energy = (chan-FPLa[chip][st])*FPLm[chip][st];
	  BGSevent->bCal[st] = energy;
	} else {
	  chan = (double)data[idxh+st] + (double)rand()/(double)RAND_MAX-0.5;
	  elin = (chan-FPHa[chip][st])*FPHm[chip][st];
	  eexp = FPHc[chip][st]*exp(chan*FPHd[chip][st]);
	  if (elin >= eexp) { energy = elin; }
	  else { energy = eexp; }
	  BGSevent->bCal[st] = energy;
	}
	
	/* Now, remember the two largest hits from the front
	   of the focal plane chips. */
	if (energy > BGSevent->Benergy1) {
	  BGSevent->Benergy2 = BGSevent->Benergy1;
	  BGSevent->Bstrip2 = BGSevent->Bstrip1;
	  BGSevent->Benergy1 = energy;
	  BGSevent->Bstrip1 = st;
	} else if (energy > BGSevent->Benergy2) {
	  BGSevent->Benergy2 = energy;
	  BGSevent->Bstrip2 = st;
	}
      } // Done with checking above threshold
    } // Done with stepping through strips
  } // Done with back strips
	    
  energy = 0;

  /* Step through the upstream strips and find the largest above-threshold
     hit (only a high-energy branch due to lack of ADCs). */
  if ( (data[iust] >= ustdc_min) && (data[iust] < ustdc_max) ) {
    for (st=0; st<nus; st++) {
      if (data[iush+st] > USMinChan[st]) {
	chan = (double)data[iush+st]+(double)rand()/(double)RAND_MAX-0.5;
	elin = (chan-USHa[st])*USHm[st];
	eexp = USHc[st]*exp(chan*USHd[st]);
	if (elin >= eexp) { energy = elin; }
	else { energy = eexp; }

	if (energy > BGSevent->Uenergy1) {
	  BGSevent->Uenergy1 = energy;
	  BGSevent->Ustrip1 = st;
	}
      } // End of checking above threshhold
    } // Done stepping through upstream strips
  } // End of upstreams
  
  energy = 0;

  /* Step through the punchthru detectors and find the largest 
     above-threshold hit */
  if ( (data[iptt+st] >= pttdc_min) && (data[iptt] < pttdc_max) ) {
    for (st=0; st<npt; st++) {
      if ((data[iptl+st]>PTMinChan[st])) {
	chan = (double)data[iptl+st] + (double)rand()/(double)RAND_MAX-0.5;
	energy = (chan-PTLa[st])*PTLm[st];
	
	if (energy > BGSevent->Penergy1) {
	  BGSevent->Penergy1 = energy;
	  BGSevent->Pstrip1 = st;
	}
      } // End of thresholds check
    } // End of stepping through strips
  } // End of punch-throughs

  energy = 0;

  /* Step through the wing detector strips and find the largest 
     above-threshold hit (only a high-energy branch dues to lack of ADCs). */
  for (st=0; st<nws; st++) {
    if (data[iwdh+st] > WMinChan[st]) {
      chan = (double)data[iwdh+st] + (double)rand()/(double)RAND_MAX-0.5;
      elin = (chan-WHa[st])*WHm[st];
      eexp = WHc[st]*exp(chan*WHd[st]);
      if (elin >= eexp) { energy = elin; }
      else { energy = eexp; }
      
      if (energy > BGSevent->Wenergy1) {
	BGSevent->Wenergy1 = energy;
	BGSevent->Wstrip1 = st;
      }
    } // End of threshold check
  }  // End of wing detectors
  
  /* Try to reconstruct a fission with upstream detector. */
  if ( (recon_fiss & 4) && (BGSevent->Benergy1 > 0.) && 
       (BGSevent->Fenergy1 > 0.) && (BGSevent->Uenergy1 > 0.) && 
       (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1)< fiss_match) &&
       (BGSevent->Fenergy1 + BGSevent->Uenergy1 >= fiss_mine) &&
       (BGSevent->Fenergy1 + BGSevent->Uenergy1 <= fiss_maxe) &&
       (BGSevent->Uenergy1 < BGSevent->Fenergy1) &&
       (BGSevent->Uenergy1 > BGSevent->Fenergy2) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->TEenergy = BGSevent->Fenergy2;
    BGSevent->FIenergy = BGSevent->ORenergy + BGSevent->TEenergy;
    BGSevent->SIenergy = BGSevent->ORenergy + BGSevent->TEenergy;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);

    if ( (BGSevent->Ustrip1>=0) && (BGSevent->Ustrip1<4) ) { // 1 o'clock chip
      BGSevent->TEchip = 5;
      BGSevent->TEfront = BGSevent->Ustrip1;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.x - t_later.x) + (t_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.y - t_later.y) + (t_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.z - t_later.z) + (t_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=4) && (BGSevent->Ustrip1<8) ) { // 3 o'clock chip
      BGSevent->TEchip = 6;
      BGSevent->TEfront = BGSevent->Ustrip1-4;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.x - w_early.x) + (w_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.y - w_early.y) + (w_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.z - w_early.z) + (w_early.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=8) && (BGSevent->Ustrip1<12) ) { // 5 o'clock chip
      BGSevent->TEchip = 7;
      BGSevent->TEfront = BGSevent->Ustrip1-8;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.x - w_later.x) + (w_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.y - w_later.y) + (w_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.z - w_later.z) + (w_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=12) && (BGSevent->Ustrip1<16) ) { // 7 o'clock strip
      BGSevent->TEchip = 8;
      BGSevent->TEfront = BGSevent->Ustrip1-12;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.x - e_early.x) + (e_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.y - e_early.y) + (e_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.z - e_early.z) + (e_early.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=16) && (BGSevent->Ustrip1<20) ) { // 9 o'clock strip
      BGSevent->TEchip = 9;
      BGSevent->TEfront = BGSevent->Ustrip1-16;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.x - e_later.x) + (e_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.y - e_later.y) + (e_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.z - e_later.z) + (e_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=20) && (BGSevent->Ustrip1<24) ) { // 11 o'clock strip
      BGSevent->TEchip = 10;
      BGSevent->TEfront = BGSevent->Ustrip1-20;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.x - t_early.x) + (t_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.y - t_early.y) + (t_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.z - t_early.z) + (t_early.z) - 1.;
    }
  }

  /* Now, try to reconstruct fission with another FP chip, two different 
     back strips. */
  else if ( (recon_fiss & 1) && 
	    (BGSevent->Fenergy1 > 0.) && (BGSevent->Fenergy2 > 0.) &&
	    (BGSevent->Benergy1 > 0.) &&
	    ( (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1)<fiss_match &&
	       fabs(BGSevent->Fenergy2 - BGSevent->Benergy2)<fiss_match &&
	       BGSevent->Benergy2>10000.) ||
	      (fabs(BGSevent->Fenergy1 + BGSevent->Fenergy2 - BGSevent->Benergy1)<fiss_match && 
	       (recon_fiss&2)>1 && BGSevent->Benergy2<=10000.) ) &&
	    (BGSevent->Fenergy2 > BGSevent->Uenergy1) &&
	    (BGSevent->Fenergy1 + BGSevent->Fenergy2 >= fiss_mine) &&
	    (BGSevent->Fenergy1 + BGSevent->Fenergy2 <= fiss_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->TEenergy = BGSevent->Fenergy2;
    BGSevent->FIenergy = BGSevent->ORenergy + BGSevent->TEenergy;
    BGSevent->SIenergy = BGSevent->ORenergy + BGSevent->TEenergy;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
    BGSevent->TEchip = BGSevent->Fchip2;
    BGSevent->TEfront = BGSevent->Fstrip2;
    if (BGSevent->Benergy2 > 0.) { BGSevent->TEback = BGSevent->Bstrip2; }
    else { BGSevent->TEback = BGSevent->Bstrip1; }
    BGSevent->TEpos = xyz(BGSevent->TEchip, BGSevent->TEfront, 
			  BGSevent->TEback, 1);
  }			       

  /* Now, reconstruct an alpha with upstream detector. */
  else if ( (recon_alph & 4) && 
	    (BGSevent->Benergy1 > 0.) && (BGSevent->Fenergy1 > 0.) &&
	    (BGSevent->Uenergy1 > 0.) && 
	    (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1)<alph_match) &&
	    (BGSevent->Uenergy1 >= alph_minterme) &&
	    (BGSevent->Uenergy1 > BGSevent->Fenergy2) &&
	    (BGSevent->Fenergy1 + BGSevent->Uenergy1 >= alph_mine) &&
	    (BGSevent->Fenergy1 + BGSevent->Uenergy1 <= alph_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->TEenergy = BGSevent->Uenergy1;
    BGSevent->ALenergy = (BGSevent->ORenergy + BGSevent->TEenergy)*ADRC;
    BGSevent->SIenergy = BGSevent->ORenergy + BGSevent->TEenergy;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
  
    if ( (BGSevent->Ustrip1>=0) && (BGSevent->Ustrip1<4) ) { // 1 o'clock chip
      BGSevent->TEchip = 5;
      BGSevent->TEfront = BGSevent->Ustrip1;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.x - t_later.x) + (t_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.y - t_later.y) + (t_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.z - t_later.z) + (t_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=4) && (BGSevent->Ustrip1<8) ) { // 3 o'clock chip
      BGSevent->TEchip = 6;
      BGSevent->TEfront = BGSevent->Ustrip1-4;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.x - w_early.x) + (w_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.y - w_early.y) + (w_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.z - w_early.z) + (w_early.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=8) && (BGSevent->Ustrip1<12) ) { // 5 o'clock chip
      BGSevent->TEchip = 7;
      BGSevent->TEfront = BGSevent->Ustrip1-8;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.x - w_later.x) + (w_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.y - w_later.y) + (w_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.z - w_later.z) + (w_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=12) && (BGSevent->Ustrip1<16) ) { // 7 o'clock strip
      BGSevent->TEchip = 8;
      BGSevent->TEfront = BGSevent->Ustrip1-12;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.x - e_early.x) + (e_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.y - e_early.y) + (e_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.z - e_early.z) + (e_early.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=16) && (BGSevent->Ustrip1<20) ) { // 9 o'clock strip
      BGSevent->TEchip = 9;
      BGSevent->TEfront = BGSevent->Ustrip1-16;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.x - e_later.x) + (e_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.y - e_later.y) + (e_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.z - e_later.z) + (e_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=20) && (BGSevent->Ustrip1<24) ) { // 11 o'clock strip
      BGSevent->TEchip = 10;
      BGSevent->TEfront = BGSevent->Ustrip1-20;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.x - t_early.x) + (t_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.y - t_early.y) + (t_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.z - t_early.z) + (t_early.z) - 1.;
    }
    
    /* Correct the alpha energy for energy loss in detector dead layers. */
    ALpathlength = sqrt( ((BGSevent->TEpos.x - BGSevent->ORpos.x)*
			  (BGSevent->TEpos.x - BGSevent->ORpos.x)) + 
			 ((BGSevent->TEpos.y - BGSevent->ORpos.y)*
			  (BGSevent->TEpos.y - BGSevent->ORpos.y)) + 
			 ((BGSevent->TEpos.z - BGSevent->ORpos.z)*
			  (BGSevent->TEpos.z - BGSevent->ORpos.z)) );
    amx = (BGSevent->TEpos.x - BGSevent->ORpos.x)/ALpathlength;
    amy = (BGSevent->TEpos.y - BGSevent->ORpos.y)/ALpathlength;
    amz = (BGSevent->TEpos.z - BGSevent->ORpos.z)/ALpathlength;
    
    if (BGSevent->ORchip == 0) {
      ORdeadlength = deadlayer/(fabs(amx*tmx + amy*tmy + amz*tmz)/
				sqrt(tmx*tmx + tmy*tmy + tmz*tmz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->ORchip == 1) {
      ORdeadlength = deadlayer/(fabs(amx*emx + amy*emy + amz*emz)/
				sqrt(emx*emx + emy*emy + emz*emz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->ORchip == 2) {
      ORdeadlength = deadlayer/(fabs(amx*wmx + amy*wmy + amz*wmz)/
				sqrt(wmx*wmx + wmy*wmy + wmz*wmz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    }
    if (BGSevent->TEchip == 5) {
      TEdeadlength = deadlayer/(fabs(amx*umx01 + amy*umy01 + amz*umz01)/
				sqrt(umx01*umx01 + umy01*umy01 + umz01*umz01)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 6) {
      TEdeadlength = deadlayer/(fabs(amx*umx03 + amy*umy03 + amz*umz03)/
				sqrt(umx03*umx03 + umy03*umy03 + umz03*umz03)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 7) {
      TEdeadlength = deadlayer/(fabs(amx*umx05 + amy*umy05 + amz*umz05)/
				sqrt(umx05*umx05 + umy05*umy05 + umz05*umz05)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 8) {
      TEdeadlength = deadlayer/(fabs(amx*umx07 + amy*umy07 + amz*umz07)/
				sqrt(umx07*umx07 + umy07*umy07 + umz07*umz07)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 9) {
      TEdeadlength = deadlayer/(fabs(amx*umx09 + amy*umy09 + amz*umz09)/
				sqrt(umx09*umx09 + umy09*umy09 + umz09*umz09)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 10) {
      TEdeadlength = deadlayer/(fabs(amx*umx11 + amy*umy11 + amz*umz11)/
				sqrt(umx11*umx11 + umy11*umy11 + umz11*umz11)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    }

    /* Micron detector have 0.3 um Si and 0.3 um Al dead layers.  
       Approximate with 0.6 um Si equivalent.  dE/dx in Si was
       approximated fitting a line between SRIM results from 6-10MeV. */
    BGSevent->ALenergy += ((ORdeadlength + TEdeadlength)*
			   (178.92 - 0.00889*BGSevent->ALenergy));
  }
  
  /* And if that didn't happen, try to reconstruct an alpha with another
     FP chip. */
  else if ( (recon_alph & 1) && 
	    (BGSevent->Fenergy1 > 0.) && (BGSevent->Fenergy2 > 0.) &&
	    (BGSevent->Benergy1 > 0.) &&
	    (BGSevent->Fenergy1 > alph_mintee) && 
	    (BGSevent->Fenergy2 > alph_minore) &&
	    ( ((fabs(BGSevent->Fenergy1 - BGSevent->Benergy1) < alph_match) &&
	       (fabs(BGSevent->Fenergy2 - BGSevent->Benergy2) < alph_match) &&
	       BGSevent->Benergy2 > 0.) ||
	      ((fabs(BGSevent->Fenergy1 + BGSevent->Fenergy2 - BGSevent->Benergy1) < alph_match) 
	       && (recon_alph & 2) && (BGSevent->Benergy2 == 0.)) ) &&
	    (BGSevent->Fenergy2 >= BGSevent->Uenergy1) &&
	    (BGSevent->Fenergy1 + BGSevent->Fenergy2 >= alph_mine) &&
	    (BGSevent->Fenergy1 + BGSevent->Fenergy2 <= alph_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy2;
    BGSevent->TEenergy = BGSevent->Fenergy1;
    BGSevent->ALenergy = (BGSevent->ORenergy + BGSevent->TEenergy)*ADRC;
    BGSevent->SIenergy = BGSevent->ORenergy + BGSevent->TEenergy;
    BGSevent->ORchip = BGSevent->Fchip2;
    BGSevent->ORfront = BGSevent->Fstrip2;
    if (BGSevent->Benergy2 > 0.) { BGSevent->ORback = BGSevent->Bstrip2; }
    else { BGSevent->ORback = BGSevent->Bstrip1; }
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
    BGSevent->TEchip = BGSevent->Fchip1;
    BGSevent->TEfront = BGSevent->Fstrip1;
    BGSevent->TEback = BGSevent->Bstrip1;
    BGSevent->TEpos = xyz(BGSevent->TEchip, BGSevent->TEfront, 
			  BGSevent->TEback, 1);

    /* Now correct the alpha energy for energy loss in the dead layers. */
    ALpathlength = sqrt( ((BGSevent->TEpos.x - BGSevent->ORpos.x)*
			  (BGSevent->TEpos.x - BGSevent->ORpos.x)) + 
			 ((BGSevent->TEpos.y - BGSevent->ORpos.y)*
			  (BGSevent->TEpos.y - BGSevent->ORpos.y)) + 
			 ((BGSevent->TEpos.z - BGSevent->ORpos.z)*
			  (BGSevent->TEpos.z - BGSevent->ORpos.z)) );
    amx = (BGSevent->TEpos.x - BGSevent->ORpos.x)/ALpathlength;
    amy = (BGSevent->TEpos.y - BGSevent->ORpos.y)/ALpathlength;
    amz = (BGSevent->TEpos.z - BGSevent->ORpos.z)/ALpathlength;
    
    if (BGSevent->ORchip == 0) {
      ORdeadlength = deadlayer/(fabs(amx*tmx + amy*tmy + amz*tmz)/
				sqrt(tmx*tmx + tmy*tmy + tmz*tmz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->ORchip == 1) {
      ORdeadlength = deadlayer/(fabs(amx*emx + amy*emy + amz*emz)/
				sqrt(emx*emx + emy*emy + emz*emz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->ORchip == 2) {
      ORdeadlength = deadlayer/(fabs(amx*wmx + amy*wmy + amz*wmz)/
				sqrt(wmx*wmx + wmy*wmy + wmz*wmz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    }
    if (BGSevent->TEchip == 0) {
      TEdeadlength = deadlayer/(fabs(amx*tmx + amy*tmy + amz*tmz)/
				sqrt(tmx*tmx + tmy*tmy + tmz*tmz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 1) {
      TEdeadlength = deadlayer/(fabs(amx*emx + amy*emy + amz*emz)/
				sqrt(emx*emx + emy*emy + emz*emz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    } else if (BGSevent->TEchip == 2) {
      TEdeadlength = deadlayer/(fabs(amx*wmx + amy*wmy + amz*wmz)/
				sqrt(wmx*wmx + wmy*wmy + wmz*wmz)/
				sqrt(amx*amx + amy*amy + amz*amz) );
    }

    /* Micron detector have 0.3 um Si and 0.3 um Al dead layers.  
       Approximate with 0.6 um Si equivalent.  dE/dx in Si was
       approximated fitting a line between SRIM results from 6-10MeV. */
    BGSevent->ALenergy += ((ORdeadlength + TEdeadlength)*
			   (178.92 - 0.00889*BGSevent->ALenergy));
  }
  
  /* Still more: try to reconstruct a beta (c.e.) energy between FP and UP */
  else if ( (recon_beta & 4) && 
	    (BGSevent->Benergy1 > 0.) && (BGSevent->Fenergy1 > 0.) &&
	    (BGSevent->Uenergy1 > 0.) && 
	    (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1) < beta_match) &&
	    (BGSevent->Uenergy1 >= BGSevent->Fenergy2) && 
	    (BGSevent->Fenergy1 + BGSevent->Uenergy1 >= beta_mine) &&
	    (BGSevent->Fenergy1 + BGSevent->Uenergy1 <= beta_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->TEenergy = BGSevent->Uenergy1;
    BGSevent->CEenergy = (BGSevent->ORenergy + BGSevent->TEenergy);
    BGSevent->SIenergy = (BGSevent->ORenergy + BGSevent->TEenergy);
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
    
    if ( (BGSevent->Ustrip1>=0) && (BGSevent->Ustrip1<4) ) { // 1 o'clock chip
      BGSevent->TEchip = 5;
      BGSevent->TEfront = BGSevent->Ustrip1;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.x - t_later.x) + (t_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.y - t_later.y) + (t_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 0. + 0.5)/4.*(t_outer.z - t_later.z) + (t_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=4) && (BGSevent->Ustrip1<8) ) { // 3 o'clock chip
      BGSevent->TEchip = 6;
      BGSevent->TEfront = BGSevent->Ustrip1-4;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.x - w_early.x) + (w_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.y - w_early.y) + (w_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 4. + 0.5)/4.*(w_outer.z - w_early.z) + (w_early.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=8) && (BGSevent->Ustrip1<12) ) { // 5 o'clock chip
      BGSevent->TEchip = 7;
      BGSevent->TEfront = BGSevent->Ustrip1-8;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.x - w_later.x) + (w_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.y - w_later.y) + (w_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 8. + 0.5)/4.*(w_outer.z - w_later.z) + (w_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=12) && (BGSevent->Ustrip1<16) ) { // 7 o'clock strip
      BGSevent->TEchip = 8;
      BGSevent->TEfront = BGSevent->Ustrip1-12;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.x - e_early.x) + (e_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.y - e_early.y) + (e_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 12. + 0.5)/4.*(e_outer.z - e_early.z) + (e_early.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=16) && (BGSevent->Ustrip1<20) ) { // 9 o'clock strip
      BGSevent->TEchip = 9;
      BGSevent->TEfront = BGSevent->Ustrip1-16;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.x - e_later.x) + (e_later.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.y - e_later.y) + (e_later.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 16. + 0.5)/4.*(e_outer.z - e_later.z) + (e_later.z) - 1.;
      
    } else if ( (BGSevent->Ustrip1>=20) && (BGSevent->Ustrip1<24) ) { // 11 o'clock strip
      BGSevent->TEchip = 10;
      BGSevent->TEfront = BGSevent->Ustrip1-20;
      BGSevent->TEpos.x = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.x - t_early.x) + (t_early.x);
      BGSevent->TEpos.y = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.y - t_early.y) + (t_early.y);
      BGSevent->TEpos.z = ((double)BGSevent->Ustrip1 - 20. + 0.5)/4.*(t_outer.z - t_early.z) + (t_early.z) - 1.;
    }
  }

  /* And if that didn't happen, try to reconstruct a beta with another
     FP chip. */
  else if ( (recon_beta & 1) && 
	    (BGSevent->Fenergy1 > 0.) && (BGSevent->Fenergy2 > 0.) &&
	    (BGSevent->Benergy1 > 0.) &&
	    ( ((fabs(BGSevent->Fenergy1 - BGSevent->Benergy1) < beta_match) &&
	       (fabs(BGSevent->Fenergy2 - BGSevent->Benergy2) < beta_match) &&
	       BGSevent->Benergy2 > 0.) ||
	      ((fabs(BGSevent->Fenergy1 + BGSevent->Fenergy2 - BGSevent->Benergy1) < beta_match) 
	       && (recon_beta & 2) && (BGSevent->Benergy2 == 0.)) ) &&
	    (BGSevent->Fenergy2 >= BGSevent->Uenergy1) &&
	    (BGSevent->Fenergy1 + BGSevent->Fenergy2 >= beta_mine) &&
	    (BGSevent->Fenergy1 + BGSevent->Fenergy2 <= beta_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->TEenergy = BGSevent->Fenergy2;
    BGSevent->CEenergy = (BGSevent->ORenergy + BGSevent->TEenergy);
    BGSevent->SIenergy = (BGSevent->ORenergy + BGSevent->TEenergy);
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
    BGSevent->TEchip = BGSevent->Fchip2;
    BGSevent->TEfront = BGSevent->Fstrip2;
    if (BGSevent->Benergy2 > 0.) { BGSevent->TEback = BGSevent->Bstrip2; }
    else { BGSevent->TEback = BGSevent->Bstrip1; }
    BGSevent->TEpos = xyz(BGSevent->TEchip, BGSevent->TEfront, 
			  BGSevent->TEback, 1);
  }

  /* What if the reconstruct options are off, or single-ended events? */

  /* Fission-like energies */
  else if ( (BGSevent->Fenergy1 > 0.) && (BGSevent->Benergy1 > 0.) &&
	    (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1) < fiss_match) &&
	    (BGSevent->Fenergy1 >= fiss_mine) && 
	    (BGSevent->Fenergy1 <= fiss_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->FIenergy = BGSevent->Fenergy1;
    BGSevent->SIenergy = BGSevent->Fenergy1;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
  }

  /* Alpha-like energies */
  else if ( (BGSevent->Fenergy1 > 0.) && (BGSevent->Benergy1 > 0.) &&
	    ( (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1) < alph_match) ||
	      (BGSevent->Benergy2 > 0.) && (recon_alph & 8) &&
	      (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1 - 
		    BGSevent->Benergy2)< alph_match) &&
	      abs(BGSevent->Bstrip1 - BGSevent->Bstrip2)==1 ) &&
	    (BGSevent->Fenergy1 >= alph_mine) && 
	    (BGSevent->Fenergy1 <= alph_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->ALenergy = BGSevent->Fenergy1*ADRC;
    BGSevent->SIenergy = BGSevent->Fenergy1;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
  }

  /* Beta-like energies */
  else if ( (BGSevent->Fenergy1 > 0.) && (BGSevent->Benergy1 > 0.) &&
	    (fabs(BGSevent->Fenergy1 - BGSevent->Benergy1) < beta_match) &&
	    (BGSevent->Fenergy1 >= beta_mine) && 
	    (BGSevent->Fenergy1 <= beta_maxe) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->CEenergy = BGSevent->Fenergy1;
    BGSevent->SIenergy = BGSevent->Fenergy1;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
  }

  /* Finally, other energies. */
  else if ( (BGSevent->Fenergy1 > 0.) && (BGSevent->Benergy1 > 0.) ) {
    BGSevent->ORenergy = BGSevent->Fenergy1;
    BGSevent->SIenergy = BGSevent->Fenergy1;
    BGSevent->ORchip = BGSevent->Fchip1;
    BGSevent->ORfront = BGSevent->Fstrip1;
    BGSevent->ORback = BGSevent->Bstrip1;
    BGSevent->ORpos = xyz(BGSevent->ORchip, BGSevent->ORfront, 
			  BGSevent->ORback, 1);
  }
  return;

}

void ClearBGSCorrelations() {
  BGSevent->timestamp = 0;
  BGSevent->Time10MHz = 0;
  BGSevent->time = 0;
  BGSevent->Last10MHzTime = 0;
  BGSevent->Overflows10MHz = 0;
  BGSevent->TOF = 0;

  BGSevent->EVENTTYPE = 0;

  for (int i=0; i<nFPstrips; i++) {
    BGSevent->fCal[0][i] = 0;
    BGSevent->fCal[1][i] = 0;
    BGSevent->fCal[2][i] = 0;
    BGSevent->bCal[i] = 0;
  }

  BGSevent->eMWPC = 0;
  BGSevent->tMWPC = 0;
  BGSevent->EVR = 0;
  BGSevent->CONV = 0;
  BGSevent->ALPHA = 0;
  BGSevent->ESCAPE = 0;
  BGSevent->FISSION = 0;
  
  BGSevent->ievr.timestamp = 0;
  BGSevent->ievr.time = 0;
  BGSevent->ievr.TOF = 0;
  BGSevent->ievr.energy = 0;
  BGSevent->ievrGamma.vgam.clear();
  
  BGSevent->dparentCE.timestamp = 0;
  BGSevent->dparentCE.time = 0;
  BGSevent->dparentCE.energy = 0;
  BGSevent->ddaughCE.timestamp = 0;
  BGSevent->ddaughCE.time = 0;
  BGSevent->ddaughCE.energy = 0;
  BGSevent->dgdaughCE.timestamp = 0;
  BGSevent->dgdaughCE.time = 0;
  BGSevent->dgdaughCE.energy = 0;

  BGSevent->dparentAL.timestamp = 0;
  BGSevent->dparentAL.time = 0;
  BGSevent->dparentAL.energy = 0;
  BGSevent->ddaughAL.timestamp = 0;
  BGSevent->ddaughAL.time = 0;
  BGSevent->ddaughAL.energy = 0;
  BGSevent->dgdaughAL.timestamp = 0;
  BGSevent->dgdaughAL.time = 0;
  BGSevent->dgdaughAL.energy = 0;

  BGSevent->dparentX.timestamp = 0;
  BGSevent->dparentX.time = 0;
  BGSevent->dparentX.energy = 0;
  BGSevent->ddaughX.timestamp = 0;
  BGSevent->ddaughX.time = 0;
  BGSevent->ddaughX.energy = 0;
  BGSevent->dgdaughX.timestamp = 0;
  BGSevent->dgdaughX.time = 0;
  BGSevent->dgdaughX.energy = 0;

  BGSevent->dfission.timestamp = 0;
  BGSevent->dfission.time = 0;
  BGSevent->dfission.energy = 0;
  
  BGSevent->Fenergy1 = 0;
  BGSevent->Fchip1 = -1;  
  BGSevent->Fstrip1 = -1; 
  BGSevent->Fenergy2 = 0;
  BGSevent->Fchip2 = -1;  
  BGSevent->Fstrip2 = -1; 
  BGSevent->Benergy1 = 0;
  BGSevent->Bstrip1 = -1; 
  BGSevent->Benergy2 = 0;
  BGSevent->Bstrip2 = -1; 
  BGSevent->Uenergy1 = 0;
  BGSevent->Ustrip1 = -1; 
  BGSevent->Penergy1 = 0;
  BGSevent->Pstrip1 = -1; 
  BGSevent->Wenergy1 = 0;
  BGSevent->Wstrip1 = -1; 
  BGSevent->ORenergy = 0;
  BGSevent->TEenergy = 0;
  BGSevent->SIenergy = 0;
  BGSevent->CEenergy = 0;
  BGSevent->ALenergy = 0;
  BGSevent->FIenergy = 0;
  BGSevent->ORpos.x = 0;
  BGSevent->ORpos.y = 0;
  BGSevent->ORpos.z = 0;
  BGSevent->ORchip = -1;
  BGSevent->ORfront = -1;
  BGSevent->ORback = -1;
  BGSevent->TEpos.x = 0;
  BGSevent->TEpos.y = 0;
  BGSevent->TEpos.z = 0;
  BGSevent->TEchip = -1; 
  BGSevent->TEfront = -1;
  BGSevent->TEback = -1;
  BGSevent->EVRenergy = 0; 

  for (int ch=0; ch<3; ch++) {
    for (int fr=0; fr<nfp; fr++) {
      for (int ba=0; ba<nfp; ba++) {
	
	evrNum[ch][fr][ba] = 0;
	parentCENum[ch][fr][ba] = 0;
	daughterCENum[ch][fr][ba] = 0;
	granddaughterCENum[ch][fr][ba] = 0;
	parentALNum[ch][fr][ba] = 0;
	daughterALNum[ch][fr][ba] = 0;
	granddaughterALNum[ch][fr][ba] = 0;
	parentXNum[ch][fr][ba] = 0;
	daughterXNum[ch][fr][ba] = 0;
	granddaughterXNum[ch][fr][ba] = 0;
	fissionNum[ch][fr][ba] = 0;
	
	for (int de=0; de<10; de++) {
	  evrCorr[ch][fr][ba][de].implanted = 0;
	  evrCorr[ch][fr][ba][de].timestamp = 0;
	  evrCorr[ch][fr][ba][de].time = 0;
	  evrCorr[ch][fr][ba][de].energy = 0;
	  evrCorr[ch][fr][ba][de].TOF = 0;
	  
	  for (int n=0; n<100; n++) {
	    evrCorr[ch][fr][ba][de].egamma[n] = 0;
	    evrCorr[ch][fr][ba][de].segegamma[n] = 0;
	    evrCorr[ch][fr][ba][de].bgsTgamma[n] = 0;	    
	  }
	  
	  parentCECorr[ch][fr][ba][de].timestamp = 0;
	  parentCECorr[ch][fr][ba][de].time = 0;
	  parentCECorr[ch][fr][ba][de].energy = 0;
	  daughterCECorr[ch][fr][ba][de].timestamp = 0;
	  daughterCECorr[ch][fr][ba][de].time = 0;
	  daughterCECorr[ch][fr][ba][de].energy = 0;
	  granddaughterCECorr[ch][fr][ba][de].timestamp = 0;
	  granddaughterCECorr[ch][fr][ba][de].time = 0;
	  granddaughterCECorr[ch][fr][ba][de].energy = 0;
	  
	  parentALCorr[ch][fr][ba][de].timestamp = 0;
	  parentALCorr[ch][fr][ba][de].time = 0;
	  parentALCorr[ch][fr][ba][de].energy = 0;
	  daughterALCorr[ch][fr][ba][de].timestamp = 0;
	  daughterALCorr[ch][fr][ba][de].time = 0;
	  daughterALCorr[ch][fr][ba][de].energy = 0;
	  granddaughterALCorr[ch][fr][ba][de].timestamp = 0;
	  granddaughterALCorr[ch][fr][ba][de].time = 0;
	  granddaughterALCorr[ch][fr][ba][de].energy = 0;
	  
	  parentXCorr[ch][fr][ba][de].timestamp = 0;
	  parentXCorr[ch][fr][ba][de].time = 0;
	  parentXCorr[ch][fr][ba][de].energy = 0;
	  daughterXCorr[ch][fr][ba][de].timestamp = 0;
	  daughterXCorr[ch][fr][ba][de].time = 0;
	  daughterXCorr[ch][fr][ba][de].energy = 0;
	  granddaughterXCorr[ch][fr][ba][de].timestamp = 0;
	  granddaughterXCorr[ch][fr][ba][de].time = 0;
	  granddaughterXCorr[ch][fr][ba][de].energy = 0;
	  
	  fissionCorr[ch][fr][ba][de].timestamp = 0;
	  fissionCorr[ch][fr][ba][de].time = 0;
	  fissionCorr[ch][fr][ba][de].energy = 0;
	}
      }
    }
  }

}

#endif
