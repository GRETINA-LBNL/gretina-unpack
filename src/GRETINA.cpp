#include "GRETINA.h"

ClassImp(rotationMatrix);
ClassImp(globalHeader);

ClassImp(cloverPacket);
ClassImp(cloverEvent);

ClassImp(mode3DataPacket);

ClassImp(ipOLD);
ClassImp(ip);
ClassImp(mode2Old);
ClassImp(mode2ABCD1234);
ClassImp(mode2ABCD5678);

ClassImp(g1GammaEvent);
ClassImp(g1OUT);

ClassImp(GRETINAVariables);

/**************************************************************/
/* rotationMatrix class functions *****************************/
/**************************************************************/

/*! Loads the rotation matrix that maps from crystal coordinates 
    to world coordinates based on the hole position of a detector.
    This function reads the text format file, NOT the binary file.

    \param file A string value, the filename of the text format 
           rotation matrix (typically crmat.dat)
    \return Returns 0 if successful, and 1 if the file couldn't be 
            opened 
*/

Int_t rotationMatrix::ReadMatrix(TString file) {
  FILE *fp;

  Float_t f1, f2, f3, f4;
  Int_t pos, xtal;
  Int_t nn = 0;
  char *st, str[256];
  
  fp = fopen(file.Data(), "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", file.Data());
    exit(1);
  } else {
    printf("\"%s\" open...", file.Data());
  }
  
  nn = 0;
  st = fgets(str, 256, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%i %i", &pos, &xtal);
      for (Int_t i=0; i<4; i++) {
	st = fgets(str, 256, fp);
	sscanf(str, "%f %f %f %f", &f1, &f2, &f3, &f4);
	crmat[pos-1][xtal][i][0] = f1;
	crmat[pos-1][xtal][i][1] = f2;
	crmat[pos-1][xtal][i][2] = f3;
	crmat[pos-1][xtal][i][3] = f4;
      }
      nn++;
    }
    /* Attempt to read the next line */
    st = fgets(str, 256, fp);
  }
  
  printf("Read %i rotation matrix coefficients.\n", nn);

  /* Done! */
  fclose(fp);
  return (0);
}

/**************************************************************/

/*! Calculates the world position for a point in crystal coordinates
    based on the rotation matrix provided and the hole number (position)
    of the crystal

    \param crystalID The integer crystal identification - the hole 
           number comes from crystalID/4, and the crystal number from
	   crystalID%4
    \param xyz TVector3 value of the position to be transformed into 
           world coordinates space
    \return Returns the TVector3 corresponding to the point in world
            coordinate space
*/

TVector3 rotationMatrix::crys2Lab(Int_t crystalID, TVector3 xyz) {

  Int_t detectorPosition = ((crystalID & 0xfffc)>>2);
  Int_t crystalNumber = (crystalID & 0x0003);

  detectorPosition -= 1;

  TVector3 xyzLab;
  xyzLab.SetX((Double_t)((crmat[detectorPosition][crystalNumber][0][0] * xyz.X()) +
			 (crmat[detectorPosition][crystalNumber][0][1] * xyz.Y()) +
			 (crmat[detectorPosition][crystalNumber][0][2] * xyz.Z()) +
			 (crmat[detectorPosition][crystalNumber][0][3]) ));
  xyzLab.SetY((Double_t)((crmat[detectorPosition][crystalNumber][1][0] * xyz.X()) +
			 (crmat[detectorPosition][crystalNumber][1][1] * xyz.Y()) +
			 (crmat[detectorPosition][crystalNumber][1][2] * xyz.Z()) +
			 (crmat[detectorPosition][crystalNumber][1][3]) ));
  xyzLab.SetZ((Double_t)((crmat[detectorPosition][crystalNumber][2][0] * xyz.X()) +
			 (crmat[detectorPosition][crystalNumber][2][1] * xyz.Y()) +
			 (crmat[detectorPosition][crystalNumber][2][2] * xyz.Z()) +
			 (crmat[detectorPosition][crystalNumber][2][3]) ));
  
  return xyzLab;
}

/**************************************************************/
/* g3CrystalEvent Class Functions *****************************/
/**************************************************************/

/*! Puts channels from a mode3 crystal event in numerical order
    within the vector chn<g3ChannelEvent>

    \return No return value -- directly alters class structures
*/

void g3CrystalEvent::OrderChannels() {
  g3ChannelEvent temp;
  Bool_t finished = 0;
  while (!finished) {
    finished = 1;
    for (UInt_t ui=0; ui<chn.size()-1; ui++) {
      if (chn[ui].chnNum() > chn[ui].chnNum()) {
	temp = chn[ui];
	chn[ui] = chn[ui+1];
	chn[ui+1] = temp;
	finished = 0;
      }
    }
  }
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC1 --> electronics channel 9 for a given crystal) 

    \return Returns the float value of the calibrated CC1 energy
*/

Float_t g3CrystalEvent::cc1() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 9) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC2 --> electronics channel 19 for a given crystal) 

    \return Returns the float value of the calibrated CC2 energy
*/

Float_t g3CrystalEvent::cc2() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 19) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC3 --> electronics channel 29 for a given crystal) 

    \return Returns the float value of the calibrated CC3 energy
*/

Float_t g3CrystalEvent::cc3() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 29) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the calibrated first central contact energy from mode3
    data (CC4 --> electronics channel 39 for a given crystal) 

    \return Returns the float value of the calibrated CC4 energy
*/

Float_t g3CrystalEvent::cc4() { 
 for (UInt_t ui=0; ui<chn.size(); ui++) {
   if (chn[ui].chnNum() == 39) { return chn[ui].eCal; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC1 --> electronics channel 9 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC1 energy
*/

Float_t g3CrystalEvent::cc1Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 9) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC2 --> electronics channel 19 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC2 energy
*/

Float_t g3CrystalEvent::cc2Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 19) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC3 --> electronics channel 29 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC3 energy
*/

Float_t g3CrystalEvent::cc3Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 29) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the uncalibrated first central contact energy from mode3
    data (CC4 --> electronics channel 39 for a given crystal) 

    \return Returns the float value of the raw (uncalibrated) CC4 energy
*/

Float_t g3CrystalEvent::cc4Raw() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 39) { return chn[ui].eRaw; }
  }
  return -1.0;
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC1 --> electronics channel 9 for a given crystal) 

    \return Returns the float value of the waveform calculated CC1 energy
*/

Float_t g3CrystalEvent::cc1Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 9) { return chn[ui].eCalc; }
  }
  return -1.0; 
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC2 --> electronics channel 19 for a given crystal) 

    \return Returns the float value of the waveform calculated CC2 energy
*/

Float_t g3CrystalEvent::cc2Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 19) { return chn[ui].eCalc; }
  }
  return -1.0;  
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC3 --> electronics channel 29 for a given crystal) 

    \return Returns the float value of the waveform calculated CC3 energy
*/

Float_t g3CrystalEvent::cc3Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 29) { return chn[ui].eCalc; }
  }
  return -1.0;  
}

/**************************************************************/

/*! Extracts the first central contact energy calculated from the waveform
    in  mode3 data (CC4 --> electronics channel 39 for a given crystal) 

    \return Returns the float value of the waveform calculated CC4 energy
*/

Float_t g3CrystalEvent::cc4Calc() { 
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].chnNum() == 39) { return chn[ui].eCalc; }
  }
  return -1.0;  
}

/**************************************************************/

/*! Calculates a raw (NOT Dino corrected) segment sum for the crystal, 
    with no threshold on the segments included

    \return Returns the float value of the raw segment sum for
            the crystal 
*/

Float_t g3CrystalEvent::segSumRaw() {
  return segSumRawThresh(0.0);
}

/**************************************************************/

/*! Calculates a raw (NOT Dino corrected) segment sum for the crystal,
    including segments that are above a threshold value

    \param thresh The float energy threshold; any segments below
           this value are not included in the sum
    \return Returns the float value of the raw segment sum for the crystal
*/

Float_t g3CrystalEvent::segSumRawThresh(Float_t thresh) {
  Float_t segSum = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > thresh && chn[ui].segNum < 36) {
      segSum += chn[ui].eCal;
    }
  }
  return segSum;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumRawTimed(Float_t minT, Float_t maxT) {
  Float_t segSum = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 36) {
      if (chn[ui].calcTime > minT &&
	  chn[ui].calcTime > maxT) {
	segSum += chn[ui].eCal;
      }
    }
  }
  return segSum;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumRawThreshTimed(Float_t thresh, Float_t minT, Float_t maxT) {
  Float_t segSum = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > thresh && chn[ui].segNum < 36) {
      if (chn[ui].calcTime > minT &&
	  chn[ui].calcTime > maxT) {
	segSum += chn[ui].eCal;
      }
    }
  }
  return segSum;
}

/**************************************************************/

/*! Determines the number of segments hit in the crystal above a given
    energy threshold

    \param thresh The float energy threshold -- segments with energy
           depositions above threshold are counted as hit
    \return Returns integer value of the number of segments hit above 
            the energy threshold
*/

Int_t g3CrystalEvent::segsHit(Float_t thresh) {
  Int_t hit = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > thresh && chn[ui].segNum < 36) {
      hit++;
    }
  }
  return hit;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumCalc() {
  Float_t segSum = 0.;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 36) {
      segSum += chn[ui].eCalc;
    }
  }
  return segSum;
}

/**************************************************************/

Float_t g3CrystalEvent::segSumCalcThresh(Float_t thresh) {
  Float_t segSum = 0.;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 36 && chn[ui].eCalc < thresh) {
      segSum += chn[ui].eCalc;
    }
  }
  return segSum;
}

/**************************************************************/

Int_t g3CrystalEvent::segsHitCalc(Float_t thresh) {
  Int_t hit = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCalc > thresh && chn[ui].segNum < 36) {
      hit++;
    }
  }
  return hit;
}

/**************************************************************/

Float_t g3CrystalEvent::maxSegE() {
  Float_t maxE = 0.;  Int_t max = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCal > maxE && chn[ui].segNum < 36) {
      max = ui;
      maxE = chn[ui].eCal;
    }
  }
  return maxE;
}

/**************************************************************/

Float_t g3CrystalEvent::maxSegECalc() {
  Float_t maxE = 0.;  Int_t max = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCalc > maxE && chn[ui].segNum < 36) {
      max = ui;
      maxE = chn[ui].eCalc;
    }
  }
  return maxE;
}

/**************************************************************/

Int_t g3CrystalEvent::maxSegNum() {
  Float_t maxE = 0.;  Int_t max = -1;
 
  for (Int_t u=0; u<chn.size(); u++) {
    if (chn[u].eCal > maxE && chn[u].segNum < 36) {
      max = u;
      maxE = chn[u].eCal;
    }
  }

  return chn[max].segNum;
}

/**************************************************************/

Int_t g3CrystalEvent::secondSegNum() {
  Float_t maxE = 0.;  Int_t max = -1;
  
  for (Int_t u=0; u<chn.size(); u++) {
    if (chn[u].eCal > maxE && chn[u].segNum < 36 && chn[u].eCal != maxSegE()) {
      max = u;
      maxE = chn[u].eCal;
    }
  }

  return chn[max].segNum;
}

/**************************************************************/

Int_t g3CrystalEvent::maxSegNumCalc() {
  Float_t maxE = 0.;  UInt_t max = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].eCalc > maxE && chn[ui].segNum < 36) {
      max = ui;
      maxE = chn[ui].eCalc;
    }
  }
  return chn[max].segNum;
}

/**************************************************************/

Int_t g3CrystalEvent::deepRingHit(Float_t thresh) {
  Int_t deep = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (chn[ui].segNum < 6 && chn[ui].eCal > thresh && deep < 1) {
      deep = 1; 
    } else if (chn[ui].segNum < 12 && chn[ui].segNum > 5 && 
	       chn[ui].eCal > thresh && deep < 2) {
      deep = 2;
    } else if (chn[ui].segNum < 18 && chn[ui].segNum > 11 && 
	       chn[ui].eCal > thresh && deep < 3) {
      deep = 3;
    } else if (chn[ui].segNum < 24 && chn[ui].segNum > 17 && 
	       chn[ui].eCal > thresh && deep < 4) {
      deep = 4;
    } else if (chn[ui].segNum < 30 && chn[ui].segNum > 23 && 
	       chn[ui].eCal > thresh && deep < 5) {
      deep = 5;
    } else if (chn[ui].segNum < 36 && chn[ui].segNum > 29 && 
	       chn[ui].eCal > thresh && deep < 1) {
      deep = 6;
    }
  }
  return deep;
}

/**************************************************************/

long long int g3CrystalEvent::LEDLow() {
  long long int ledLow = -1;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (ledLow == -1) { ledLow = chn[ui].timestamp; }
    if (ledLow > chn[ui].timestamp) {
      ledLow = chn[ui].timestamp;
    }
  }
  return ledLow;
}

/**************************************************************/

long long int g3CrystalEvent::LEDHigh() {
  long long int ledHigh = 0;
  for (UInt_t ui=0; ui<chn.size(); ui++) {
    if (ledHigh < chn[ui].timestamp) {
      ledHigh = chn[ui].timestamp;
    }
  }
  return ledHigh;
}

/**************************************************************/

UInt_t g3CrystalEvent::crystalBuild() { return chn.size(); }

/**************************************************************/

/*! Function calculates the difference between the high and low 
    LED timestamps for channel events within the mode3 crystal event,
    essentially the 'width' of the event for the crystal

    \return Returns the difference between the first and last LED
            timestamps as a long long int
*/

long long int g3CrystalEvent::LEDRange() {
  return (LEDHigh() - LEDLow());
}

/**************************************************************/

void g3CrystalEvent::mode3ToFile(FILE *outFile) {

  mode3DataPacket *dp;
  dp = (mode3DataPacket*)malloc(sizeof(dp->aahdr) + 
				sizeof(dp->hdr) + 
				sizeof(dp->waveform));
  memset(dp->waveform, 1, MAX_TRACE_LENGTH * sizeof(UShort_t));

  globalHeader *gHeader;
  gHeader = (globalHeader*)malloc(sizeof(gHeader));

  unsigned char outBuf[32*1024];
  unsigned char *tmp = (outBuf);

  dp->aahdr[0] = 0xAAAA;
  dp->aahdr[1] = 0xAAAA;
  
  for (Int_t i=0; i<chn.size(); i++) {

    gHeader->type = RAW;
    gHeader->timestamp = chn[i].timestamp;
    Int_t evtLength = (sizeof(dp->aahdr) + sizeof(dp->hdr) + chn[i].tracelength()*sizeof(UShort_t));
    gHeader->length = evtLength;

    dp->hdr[0] = chn[i].hdr0;
    dp->hdr[1] = chn[i].hdr1;
    dp->hdr[2] = ((chn[i].timestamp & 0xffff0000) >> 16);
    dp->hdr[3] = ((chn[i].timestamp & 0xffff));
    dp->hdr[5] = ((chn[i].timestamp & 0xffff00000000) >> 32);
    dp->hdr[6] = ((chn[i].CFDtimestamp & 0xffff));
    dp->hdr[8] = ((chn[i].CFDtimestamp & 0xffff00000000) >> 32);
    dp->hdr[9] = ((chn[i].CFDtimestamp & 0xffff0000) >> 16);

    Bool_t sign = 0;
    Int_t tmpIntEnergy = (Int_t)(32*chn[i].eRaw);
    UInt_t tmpEnergy = 0;
    if (chn[i].chanID()%10 == 9) { /* CC */
      if (chn[i].eRaw < 0) {
	tmpEnergy = (tmpIntEnergy + (UInt_t)(0x01000000));
	sign = 1;
      } else {
	tmpEnergy = (UInt_t)(tmpIntEnergy);
      }
    } else {
      if (chn[i].eRaw < 0) {
	tmpEnergy = (UInt_t)(-tmpIntEnergy);
      } else {
	tmpEnergy = (UInt_t)((0x01000000)-tmpIntEnergy);
	sign = 1;
      }
    }
    dp->hdr[4] = (tmpEnergy & 0xffff);
    dp->hdr[7] = ((tmpEnergy & 0x00ff0000) >> 16);
    if (sign) { dp->hdr[7] = (dp->hdr[7] & 0x0100); }
    dp->hdr[7] += (chn[i].pileUp() << 15);

    for (Int_t j=0; j<chn[i].tracelength()+1; j=j+2) {
      if (chn[i].wf.raw[j] < 0) {
    	dp->waveform[j+1] = (UShort_t)(chn[i].wf.raw[j] + std::numeric_limits<unsigned int>::max());
      } else {
    	dp->waveform[j+1] = (UShort_t)(chn[i].wf.raw[j]);
      }
      if (chn[i].wf.raw[j+1] < 0) {
    	dp->waveform[j] = (UShort_t)(chn[i].wf.raw[j+1] + std::numeric_limits<unsigned int>::max());
      } else {
    	dp->waveform[j] = (UShort_t)(chn[i].wf.raw[j+1]);
      }
    }

    memmove(tmp, &dp->aahdr[0], evtLength);
    for (Int_t j=0; j<evtLength; j=j+2) {
      swap(*(outBuf + j), *(outBuf + j + 1));
    }
    
    fwrite(gHeader, sizeof(struct globalHeader), 1, outFile);
    fwrite(tmp, evtLength, 1, outFile);
    
  }

}

/**************************************************************/
/* g3OUT Class Functions **************************************/
/**************************************************************/

/*! Resets the mode3 output class -- clears channel vector

    \return No return value -- clears data structures in class
            directly
*/

void g3OUT::Reset() {
  UInt_t crystals = crystalMult();
  for (UInt_t ui=0; ui<crystals; ui++) {
    xtals[ui].chn.clear();
  }
  xtals.clear();
}

/**************************************************************/

/*! Extracts the number of crystals hit within a mode3 event, namely
    the size of the xtals vector

    \return Returns the unsigned integer value of the size of the
            xtals array
*/

UInt_t g3OUT::crystalMult() { return xtals.size(); }

/**************************************************************/

Float_t g3OUT::calorimeterE() {
  Float_t sum = 0.;
  for (UInt_t ui=0; ui<crystalMult(); ui++) {
    sum += xtals[ui].cc1();
  }
  return sum;
}

/**************************************************************/
/* g2CrystalEvent Class Functions *****************************/
/**************************************************************/

void g2CrystalEvent::Reset() { 
  intpts.clear();
  crystalNum = -1;  quadNum = -1;
  timestamp = 0;
  t0 = -1.0; chiSq = -1.0;
  bl = -9999; /* 2015-09-20 CMC: added for wf baseline and PZ checks */
  error = -1;
  cc = 0.0; segSum = 0.0;
  doppler = 0.0; dopplerSeg = 0.0; dopplerCrystal = 0.0;
}

/**************************************************************/

Int_t g2CrystalEvent::numIntPts() {  return intpts.size();  }

/**************************************************************/

Int_t g2CrystalEvent::numSegHit() {
  Int_t hit[40] = {0};
  Int_t segsHit = 0;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (hit[intpts[i].segNum] == 0) {
      hit[intpts[i].segNum] = 1;
      segsHit++;
    }
  }
  return segsHit;
}

/**************************************************************/

Int_t g2CrystalEvent::segsHitAboveThresh(Float_t thresh) {
  Int_t hit[40] = {0};
  Int_t segsHit = 0;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (hit[intpts[i].segNum] == 0 && intpts[i].segE > thresh) {
      hit[intpts[i].segNum] = 1;
      segsHit++;
    }
  }
  return segsHit;
}

/**************************************************************/

Float_t g2CrystalEvent::averageIntPtsPerSeg() {
  if (numSegHit() > 0) {
    return ((Float_t)intpts.size() / (Float_t)(numSegHit()));
  } else { return 0.; }
} 

/**************************************************************/

Float_t g2CrystalEvent::segSumRaw() {
  Float_t sum = 0;
  Int_t used[40] = {0};
  for (Int_t i=0; i<intpts.size(); i++) {
    if (used[intpts[i].segNum] == 0) {
      sum += intpts[i].segE;
      used[intpts[i].segNum] = 1;
    }
  }
  return sum;
}

/**************************************************************/

Float_t g2CrystalEvent::segSumRawThresh(Float_t thresh) {
  Float_t sum = 0;
  Int_t used[40] = {0};
  for (Int_t i=0; i<intpts.size(); i++) {
    if (used[intpts[i].segNum] == 0 && intpts[i].segE > thresh) {
      sum += intpts[i].segE;
      used[intpts[i].segNum] = 1;
    }
  }
  return sum;
}

/**************************************************************/

TVector3 g2CrystalEvent::maxIntPtXYZ() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].xyz; }
  else { return TVector3(0., 0., 0.); }
}

/**************************************************************/

TVector3 g2CrystalEvent::maxIntPtXYZLab() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].xyzLab; }
  else { return TVector3(0., 0., 0.); }
}

/**************************************************************/

TVector3 g2CrystalEvent::maxIntPtXYZLabSeg() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].xyzLabSeg; }
  else { return TVector3(0., 0., 0.); }
}

/**************************************************************/

TVector3 g2CrystalEvent::maxIntPtXYZLabCrystal() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].xyzLabCrys; }
  else { return TVector3(0., 0., 0.); }
}

/**************************************************************/

Float_t g2CrystalEvent::maxIntPtE() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].e; }
  else { return 0.; }
}

/**************************************************************/

Int_t g2CrystalEvent::maxIntPtSegNum() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].segNum; }
  else { return -1; }
}

/**************************************************************/

Float_t g2CrystalEvent::maxIntPtSegE() {
  Float_t maxE = 0; Int_t max = -1;
  for (Int_t i=0; i<intpts.size(); i++) {
    if (intpts[i].e > maxE) {
      maxE = intpts[i].e;
      max = i;
    }
  }
  if (max > -1) { return intpts[max].segE; }
  else { return 0.0; }
}

/**************************************************************/

Float_t g2CrystalEvent::gTheta() {
  if (numIntPts() > 0) {
    return maxIntPtXYZLab().Theta();
  } else { return 0.0; }
}

/**************************************************************/

Float_t g2CrystalEvent::gPhi() {
  if (numIntPts() > 0) {
    if (maxIntPtXYZLab().Phi() < 0) {
      return (maxIntPtXYZLab().Phi() + TMath::TwoPi());
    } else {
      return maxIntPtXYZLab().Phi();
    }
  } else { return 0.0; }
}

/**************************************************************/

Float_t g2CrystalEvent::CosDop() {
  if (numIntPts() > 0) {
    return maxIntPtXYZLab().CosTheta();
  } else {
    return 0.0;
  }
}

/**************************************************************/

Float_t g2CrystalEvent::cDoppler(Float_t beta) {
  if (numIntPts() > 0) {
    Float_t cosDop = maxIntPtXYZLab().CosTheta();
    Float_t gamma = 1/TMath::Sqrt(1. - beta*beta);
    return (gamma*(1. - beta*cosDop));
  } else {
    return 0.0;
  }
}

/**************************************************************/

Float_t g2CrystalEvent::cDopplerSeg(Float_t beta) {
  if (numIntPts() > 0) {
    Float_t cosDop = maxIntPtXYZLabSeg().CosTheta();
    Float_t gamma = 1/TMath::Sqrt(1. - beta*beta);
    return (gamma*(1. - beta*cosDop));
  } else {
    return 0.0;
  }
}

/**************************************************************/

Float_t g2CrystalEvent::cDopplerCrystal(Float_t beta) {
  if (numIntPts() > 0) {
    Float_t cosDop = maxIntPtXYZLabCrystal().CosTheta();
    Float_t gamma = 1/TMath::Sqrt(1. - beta*beta);
    return (gamma*(1. - beta*cosDop));
  } else {
    return 0.0;
  }
}


/**************************************************************/
/* g2OUT Class Functions **************************************/
/**************************************************************/

void g2OUT::Reset() {
  runNumber = -1;
  for (UInt_t ui=0; ui<crystalMult(); ui++) {
    xtals[ui].intpts.clear();
    xtals[ui].waveAll.clear();
  }
  xtals.clear();
}

/**************************************************************/

Float_t g2OUT::calorimeterDoppler(Float_t beta) {
  Float_t maxE = 0;
  Int_t max = -1;
  for (UInt_t ui=0; ui<crystalMult(); ui++) {
    if (xtals[ui].cc > maxE) { maxE = xtals[ui].cc;  max = ui; }
  }
  if (max >= 0) {
    if (beta > 0) {
      return (xtals[max].cDoppler(beta));
    } else { return (xtals[max].doppler); }
  } else { return 0.; }
}

/**************************************************************/

/*! Determines the calorimeter mode energy for a mode2 event, 
    subject to an energy threshold for inclusion of a given crystal
    energy.  Calculated based on the CC energy values.
    
    \param thresh The float value of the energy threshold that is 
           checked for each crystal to be included in the calorimeter
           sum
    \return Returns the float value of the calorimeter mode energy sum
            for the mode2 event
*/

Float_t g2OUT::calorimeterE(Float_t thresh) {
  Float_t sum = 0.;
  for (UInt_t ui=0; ui<crystalMult(); ui++) {
    if (xtals[ui].cc3 > thresh) { sum += xtals[ui].cc3; }
  }
  return sum;
}

/**************************************************************/

Float_t g2OUT::calorimeterE() {
  return calorimeterE(0.0);
}

/**************************************************************/

UInt_t g2OUT::crystalMult() {
  return xtals.size();
}

/**************************************************************/
/* g1GammaEvent Class Functions *******************************/
/**************************************************************/

Float_t g1GammaEvent::gTheta() {
  return xyzLab1.Theta();
}

/**************************************************************/

Float_t g1GammaEvent::gPhi() {
  return xyzLab1.Phi();
}

/**************************************************************/

TVector3 g1GammaEvent::normalVector1() {
  TVector3 xyz;
  xyz.SetXYZ(1.0*xyzLab1.Y(), -1.0*xyzLab1.X(), 0);
  return xyz.Unit();
}

/**************************************************************/

TVector3 g1GammaEvent::normalVector2() {
  TVector3 xyz;
  xyz.SetXYZ((xyzLab1.Y()*(xyzLab2.Z()-xyzLab1.Z())-
	      xyzLab1.Z()*(xyzLab2.Y()-xyzLab1.Y())),
	     (xyzLab1.Z()*(xyzLab2.X()-xyzLab1.X())-
	      xyzLab1.X()*(xyzLab2.Z()-xyzLab1.Z())),
	     (xyzLab1.X()*(xyzLab2.Y()-xyzLab1.Y())-
	      xyzLab1.Y()*(xyzLab2.X()-xyzLab1.X())));	
  return xyz.Unit();
}

/**************************************************************/

Float_t g1GammaEvent::phiCompton() {
  Float_t angle;
  Double_t dotProd = normalVector1().Dot(normalVector2());
  if (normalVector2().Z() > 0) {
    angle = (TMath::ACos(dotProd / 
			 (normalVector1().Mag()*normalVector2().Mag())) + 
	     TMath::Pi()) * TMath::RadToDeg();
  } else if (normalVector2().Z() <= 0) {
    angle = ((TMath::ACos(dotProd / 
			  (normalVector1().Mag()*normalVector2().Mag()))) * 
	     TMath::RadToDeg());
  }
  return angle;
}

/**************************************************************/

Float_t g1GammaEvent::intPtDistance() {
  return ((normalVector2() - normalVector1()).Mag());
}

/**************************************************************/

Float_t g1GammaEvent::scatterAngle() {
  TVector3 delta = normalVector1() - normalVector2();
  Float_t dotProdScat = delta.Dot(normalVector1());
  return (TMath::ACos(dotProdScat / xyzLab1.Mag() * delta.Mag()) * 
	  TMath::RadToDeg());
}

/**************************************************************/

Float_t g1GammaEvent::cDoppler(Float_t beta) {
  Float_t cosDop = xyzLab1.CosTheta();
  Float_t gamma = 1/TMath::Sqrt(1. - beta*beta);
  return (gamma*(1. - beta*cosDop));
}

/**************************************************************/
/* g1OUT Class Functions **************************************/
/**************************************************************/

void g1OUT::Reset() {
  gammas.clear();
}

/**************************************************************/

Float_t g1OUT::calorimeterE() {
  Float_t sum = 0.;
  for (UInt_t ui=0; ui<gammaMult(); ui++) {
    sum += gammas[ui].cc;
  }
  return sum;
}

/**************************************************************/

UInt_t g1OUT::gammaMult() {
  return gammas.size();
}

/**************************************************************/
/* GRETINAVariable Class Functions ****************************/
/**************************************************************/

GRETINAVariables::GRETINAVariables() { 
  beta = 0.0;
  targetXYZ.SetXYZ(0.0, 0.0, 0.0);
  ataCor = 0.0; btaCor = 0.0;
  for (Int_t i=0; i<MAXCRYSTALS; i++) { radiusCor[i] = 1.0; }
  for (Int_t i=0; i<MAXQUADS; i++) { hole[i] = -1; electronicsOrder[i] = -1; }
  
  for (Int_t i=0; i<MAXCHANNELS; i++) {
    ehiGeOffset[i] = 0.0;
    ehiGeGain[i] = 1.0;
  }
  
  for (Int_t i=0; i<MAXCHANNELS; i++) {
    for (Int_t j=0; j<5000; j++) {
      dnlLU[i][j] = 0.0;
    }
  }
  
  for (Int_t i=0; i<(MAXCRYSTALS); i++) {
    for (Int_t j=0; j<40; j++) {
      for (Int_t k=0; k<40; k++) {
	dinoFactor[i][j][k] = 0.0;
      }
    }
  }
  
  for (Int_t i=0; i<36; i++) {
    for (Int_t j=0; j<4; j++) {
      Q1Special[j][i] = i;
      Q2Special[j][i] = i;
      QSegmentOrder[j][i] = i;
    }
  }
  
  for (Int_t i=0; i<3; i+=2) {
    Q1Special[i][ 0] =  1; Q1Special[i][ 1] =  8; Q1Special[i][ 2] =  3; Q1Special[i][ 3] =  4;
    Q1Special[i][ 4] = 35; Q1Special[i][ 5] =  0; Q1Special[i][ 6] =  7; Q1Special[i][ 7] = 14;
    Q1Special[i][ 8] =  9; Q1Special[i][ 9] = 10; Q1Special[i][10] = 29; Q1Special[i][11] = 30;
    Q1Special[i][12] = 13; Q1Special[i][13] = 20; Q1Special[i][14] = 15; Q1Special[i][15] = 16;
    Q1Special[i][16] = 23; Q1Special[i][17] = 24; 
  }
  
  Q1Special[0][18] = 19; Q1Special[0][19] = 26; Q1Special[0][20] = 21; Q1Special[0][21] = 22; 
  Q1Special[0][22] = 17; Q1Special[0][23] = 18; Q1Special[0][24] = 25; Q1Special[0][25] = 32; 
  Q1Special[0][26] = 27; Q1Special[0][27] = 28; Q1Special[0][28] = 11; Q1Special[0][29] = 12;
  Q1Special[0][30] = 31; Q1Special[0][31] =  2; Q1Special[0][32] = 33; Q1Special[0][33] = 34; 
  Q1Special[0][34] =  5; Q1Special[0][35] =  6;
  
  Q1Special[1][18] = 27; Q1Special[1][19] = 28; Q1Special[1][20] = 29; Q1Special[1][21] = 30;
  Q1Special[1][22] = 31; Q1Special[1][23] = 32; Q1Special[1][24] = 33; Q1Special[1][25] = 34;
  Q1Special[1][26] = 35; Q1Special[1][27] = 18; Q1Special[1][28] = 19; Q1Special[1][29] = 20;
  Q1Special[1][30] = 21; Q1Special[1][31] = 22; Q1Special[1][32] = 23; Q1Special[1][33] = 24;
  Q1Special[1][34] = 25; Q1Special[1][35] = 26;
  
  Q1Special[2][18] = 28; Q1Special[2][19] = 11; Q1Special[2][20] = 12; Q1Special[2][21] = 31; 
  Q1Special[2][22] =  2; Q1Special[2][23] = 33; Q1Special[2][24] = 34; Q1Special[2][25] =  5; 
  Q1Special[2][26] =  6; Q1Special[2][27] = 19; Q1Special[2][28] = 26; Q1Special[2][29] = 21;
  Q1Special[2][30] = 22; Q1Special[2][31] = 17; Q1Special[2][32] = 18; Q1Special[2][33] = 25; 
  Q1Special[2][34] = 32; Q1Special[2][35] = 27;
  
  Q2Special[0][16] = 22; Q2Special[0][22] = 16;
  Q2Special[2][16] = 22; Q2Special[2][22] = 16;
}

/**************************************************************/

void GRETINAVariables::Initialize() {
  beta = 0.0;
  targetXYZ.SetXYZ(0.0, 0.0, 0.0);
  ataCor = 0.0; btaCor = 0.0;
  for (Int_t i=0; i<MAXCRYSTALS; i++) { radiusCor[i] = 1.0; }
  for (Int_t i=0; i<MAXQUADS; i++) { hole[i] = -1; electronicsOrder[i] = -1; }
  
  minCrossTime = 70; maxCrossTime = 90;
}

/**************************************************************/

Double_t GRETINAVariables::GetNextValue(FILE *file) {
 
  /* This function assumes the line is a string descriptor
     and a value -- returns the value. */
  char line[100];
  char junk[300];
  double value;

  char* retVal;
  
  retVal = fgets(line, 100, file);
  sscanf(line, "%s %lf", junk, &value);
  return value;
}

/**************************************************************/

void GRETINAVariables::InitializeGRETINAVariables(TString inputFilename) {
  FILE *in;
  if ( (in = fopen(inputFilename.Data(), "r")) == NULL ) {
    perror("Input variable file could not be opened.");
  }

  cout << "Reading variable initialization file: "
       << inputFilename.Data() << endl;

  char line[300];
  double value[20];
  int intvalue[5];
  int num = 0;

  char* retVal;

  while ( !feof(in) ) {
    retVal = fgets(line, 300, in);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, "Doppler Correction", 18) == 0) {
      beta = GetNextValue(in);
      ataCor = GetNextValue(in);
      btaCor = GetNextValue(in);
      targetXYZ.SetX(GetNextValue(in));
      targetXYZ.SetY(GetNextValue(in));
      targetXYZ.SetZ(GetNextValue(in));
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 100, in);
      sscanf(line, "%d", &num);
      for (int i=0; i<num; i++) {
    	retVal = fgets(line, 300, in);
    	sscanf(line, "%d %lf", &intvalue[0], &value[0]);
    	radiusCor[(intvalue[0]-1)] = value[0];
      } 
    }
    if (strncmp(line, "Geometry", 8) == 0) {
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 100, in); /* # of quads identified */
      sscanf(line, "%d", &num);
      for (int i=0; i<num; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%d %d", &intvalue[0], &intvalue[1]);
	hole[intvalue[0]-1] = intvalue[1];
      }
    }
    if (strncmp(line, "Electronics Order", 17) == 0) {
      retVal = fgets(line, 300, in); /* Comment */
      retVal = fgets(line, 100, in); /* # of quads identified */
      sscanf(line, "%d", &num);
      for (int i=0; i<num; i++) {
	retVal = fgets(line, 300, in);
	sscanf(line, "%d %d", &intvalue[0], &intvalue[1]);
	electronicsOrder[intvalue[0]-1] = intvalue[1];
      }
    }
  }

  fclose(in);
}

/**************************************************************/

void GRETINAVariables::Reset() {
  beta = 0.0;
  targetXYZ.SetXYZ(0.0, 0.0, 0.0);
  ataCor = 0.0; btaCor = 0.0;
  for (Int_t i=0; i<MAXCRYSTALS; i++) { radiusCor[i] = 1.0; }
  for (Int_t i=0; i<MAXQUADS; i++) { hole[i] = -1; electronicsOrder[i] = -1; }
}

/**************************************************************/
/* GRETINA Class Functions ************************************/
/**************************************************************/

void GRETINA::Initialize() {

  wfMinCrossTime = 70.; wfMaxCrossTime = 90.;
  
  cloverModuleNumber = 99999; /* No clover in setup; 
				 all GlobalID should be available. */

  // cloverModuleNumber = 37;
  /* Old clover calibration */
  cloverSlope[30] = 0.0980; cloverOffset[30] = 0.2323;
  cloverSlope[31] = 0.1331; cloverOffset[31] = 1.1686;
  cloverSlope[32] = 0.1299; cloverOffset[32] = 0.8119;
  cloverSlope[33] = 0.0973; cloverOffset[33] = -1.1761;
  cloverSlope[35] = 0.1342; cloverOffset[35] = 0.1717;
  cloverSlope[36] = 0.1299; cloverOffset[36] = 0.1809;
  cloverSlope[37] = 0.1319; cloverOffset[37] = 1.0462;
  cloverSlope[38] = 0.1313; cloverOffset[38] = 1.3180;

  /* Run0296 Calibration */
  //cloverSlope[30] = 0.0980; cloverOffset[30] = 0.2323;
  //cloverSlope[31] = 0.1331; cloverOffset[31] = 1.1686;
  //cloverSlope[32] = 0.1299; cloverOffset[32] = 0.8119;
  //cloverSlope[33] = 0.0973; cloverOffset[33] = -1.1761;
  //cloverSlope[35] = 0.1342; cloverOffset[35] = 0.1717;
  //cloverSlope[36] = 0.1299; cloverOffset[36] = 0.1809;
  //cloverSlope[37] = 0.1319; cloverOffset[37] = 1.0462;
  //cloverSlope[38] = 0.1313; cloverOffset[38] = 1.3180;

    /* Segment positions for Doppler corrections */
    xA[0] = 20.6;    yA[0] = 0.0;     zA[0] = 4.1;
    xA[1] = 9.7;     yA[1] = -17.2;   zA[1] = 4.1;
    xA[2] = -10.3;   yA[2] = -15.7;   zA[2] = 4.1;
    xA[3] = -18.3;   yA[3] = 0.9;     zA[3] = 4.1;
    xA[4] = -8.4;    yA[4] = 16.2;    zA[4] = 4.1;
    xA[5] = 11.1;    yA[5] = 16.6;    zA[5] = 4.1;
    xA[6] = 22.2;    yA[6] = 0.0;     zA[6] = 15.1;
    xA[7] = 10.4;    yA[7] = -18.7;   zA[7] = 15.1;
    xA[8] = -11.2;   yA[8] = -17.0;   zA[8] = 15.1;
    xA[9] = -20.0;   yA[9] = 0.9;     zA[9] = 15.1;
    xA[10] = -9.2;   yA[10] = 17.6;   zA[10] = 15.1;
    xA[11] = 12.0;   yA[11] = 18.0;   zA[11] = 15.1;
    xA[12] = 24.0;   yA[12] = 0.0;    zA[12] = 30.2;
    xA[13] = 11.2;   yA[13] = -20.4;  zA[13] = 30.2;
    xA[14] = -12.4;  yA[14] = -18.7;  zA[14] = 30.2;
    xA[15] = -21.9;  yA[15] = 1.1;    zA[15] = 30.2;
    xA[16] = -10.1;  yA[16] = 19.3;   zA[16] = 30.2;
    xA[17] = 12.9;   yA[17] = 19.7;   zA[17] = 30.2;
    xA[18] = 25.3;   yA[18] = 0.0;    zA[18] = 46.7;
    xA[19] = 11.7;   yA[19] = -21.8;  zA[19] = 46.7;
    xA[20] = -13.4;  yA[20] = -20.2;  zA[20] = 46.7;
    xA[21] = -23.7;  yA[21] = 1.2;    zA[21] = 46.7;
    xA[22] = -11.0;  yA[22] = 20.9;   zA[22] = 46.7;
    xA[23] = 13.6;   yA[23] = 21.1;   zA[23] = 46.7;
    xA[24] = 25.9;   yA[24] = -0.3;   zA[24] = 65.6;
    xA[25] = 12.1;   yA[25] = -22.7;  zA[25] = 65.6;
    xA[26] = -14.2;  yA[26] = -21.3;  zA[26] = 65.7;
    xA[27] = -25.5;  yA[27] = 1.2;    zA[27] = 65.7;
    xA[28] = -11.9;  yA[28] = 22.4;   zA[28] = 65.7;
    xA[29] = 13.7;   yA[29] = 22.0;   zA[29] = 65.6;
    xA[30] = 26.1;   yA[30] = -0.4;   zA[30] = 82.9;
    xA[31] = 12.1;   yA[31] = -23.1;  zA[31] = 82.9;
    xA[32] = -14.6;  yA[32] = -21.7;  zA[32] = 82.9;
    xA[33] = -26.2;  yA[33] = 1.4;    zA[33] = 82.9;
    xA[34] = -12.3;  yA[34] = 23.1;   zA[34] = 82.9;
    xA[35] = 13.6;   yA[35] = 22.4;   zA[35] = 82.9;
    
    xB[0] = 18.0;    yB[0] = 10.0;    zB[0] = 4.1;
    xB[1] = 18.4;    yB[1] = -9.8;    zB[1] = 4.1;
    xB[2] = 1.8;     yB[2] = -19.9;   zB[2] = 4.1;
    xB[3] = -16.1;   yB[3] = -10.1;   zB[3] = 4.1;
    xB[4] = -17.5;   yB[4] = 9.1;     zB[4] = 4.1;
    xB[5] = -0.3;    yB[5] = 19.8;    zB[5] = 4.1;
    xB[6] = 19.4;    yB[6] = 10.8;    zB[6] = 15.1;
    xB[7] = 19.9;    yB[7] = -10.6;   zB[7] = 15.1;
    xB[8] = 1.9;     yB[8] = -21.5;   zB[8] = 15.1;
    xB[9] = -17.6;   yB[9] = -11.0;   zB[9] = 15.1;
    xB[10] = -19.0;  yB[10] = 9.9;    zB[10] = 15.1;
    xB[11] = -0.5;   yB[11] = 21.5;   zB[11] = 15.1;
    xB[12] = 20.9;   yB[12] = 11.8;   zB[12] = 30.2;
    xB[13] = 21.5;   yB[13] = -11.6;  zB[13] = 30.2;
    xB[14] = 1.9;    yB[14] = -23.4;  zB[14] = 30.2;
    xB[15] = -19.3;  yB[15] = -12.0;  zB[15] = 30.2;
    xB[16] = -20.9;  yB[16] = 10.8;   zB[16] = 30.2;
    xB[17] = -0.7;   yB[17] = 23.4;   zB[17] = 30.2;
    xB[18] = 22.1;   yB[18] = 12.4;   zB[18] = 46.6;
    xB[19] = 22.6;   yB[19] = -12.2;  zB[19] = 46.6;
    xB[20] = 1.9;    yB[20] = -25.0;  zB[20] = 46.7;
    xB[21] = -20.9;  yB[21] = -12.9;  zB[21] = 46.7;
    xB[22] = -22.3;  yB[22] = 11.6;   zB[22] = 46.7;
    xB[23] = -0.8;   yB[23] = 24.9;   zB[23] = 46.7;
    xB[24] = 22.4;   yB[24] = 13.0;   zB[24] = 65.6;
    xB[25] = 22.9;   yB[25] = -12.5;  zB[25] = 65.5;
    xB[26] = 1.5;    yB[26] = -25.8;  zB[26] = 65.6;
    xB[27] = -21.8;  yB[27] = -13.7;  zB[27] = 65.6;
    xB[28] = -23.1;  yB[28] = 11.9;   zB[28] = 65.6;
    xB[29] = -0.9;   yB[29] = 25.8;   zB[29] = 65.6;
    xB[30] = 22.5;   yB[30] = 13.2;   zB[30] = 82.9;
    xB[31] = 23.1;   yB[31] = -12.6;  zB[31] = 82.9;
    xB[32] = 1.3;    yB[32] = -26.1;  zB[32] = 82.9;
    xB[33] = -22.1;  yB[33] = -14.0;  zB[33] = 82.9;
    xB[34] = -23.3;  yB[34] = 12.0;   zB[34] = 82.9;
    xB[35] = -1.0;   yB[35] = 26.0;   zB[35] = 82.9;

    QNormal[0][0] = 0;   QNormal[1][0] = 0;   QNormal[2][0] = 0;   QNormal[3][0] = 0;  
    QNormal[0][1] = 1;   QNormal[1][1] = 1;   QNormal[2][1] = 1;   QNormal[3][1] = 1;  
    QNormal[0][2] = 2;   QNormal[1][2] = 2;   QNormal[2][2] = 2;   QNormal[3][2] = 2;  
    QNormal[0][3] = 3;   QNormal[1][3] = 3;   QNormal[2][3] = 3;   QNormal[3][3] = 3;  
    QNormal[0][4] = 4;   QNormal[1][4] = 4;   QNormal[2][4] = 4;   QNormal[3][4] = 4;  
    QNormal[0][5] = 5;   QNormal[1][5] = 5;   QNormal[2][5] = 5;   QNormal[3][5] = 5;  
    QNormal[0][6] = 6;   QNormal[1][6] = 6;   QNormal[2][6] = 6;   QNormal[3][6] = 6;  
    QNormal[0][7] = 7;   QNormal[1][7] = 7;   QNormal[2][7] = 7;   QNormal[3][7] = 7;  
    QNormal[0][8] = 8;   QNormal[1][8] = 8;   QNormal[2][8] = 8;   QNormal[3][8] = 8;  
    QNormal[0][9] = 36;  QNormal[1][9] = 36;  QNormal[2][9] = 36;  QNormal[3][9] = 36; 
    QNormal[0][10] = 9;  QNormal[1][10] = 9;  QNormal[2][10] = 9;  QNormal[3][10] = 9; 
    QNormal[0][11] = 10; QNormal[1][11] = 10; QNormal[2][11] = 10; QNormal[3][11] = 10;
    QNormal[0][12] = 11; QNormal[1][12] = 11; QNormal[2][12] = 11; QNormal[3][12] = 11;
    QNormal[0][13] = 12; QNormal[1][13] = 12; QNormal[2][13] = 12; QNormal[3][13] = 12;
    QNormal[0][14] = 13; QNormal[1][14] = 13; QNormal[2][14] = 13; QNormal[3][14] = 13;
    QNormal[0][15] = 14; QNormal[1][15] = 14; QNormal[2][15] = 14; QNormal[3][15] = 14;
    QNormal[0][16] = 15; QNormal[1][16] = 15; QNormal[2][16] = 15; QNormal[3][16] = 15;
    QNormal[0][17] = 16; QNormal[1][17] = 16; QNormal[2][17] = 16; QNormal[3][17] = 16;
    QNormal[0][18] = 17; QNormal[1][18] = 17; QNormal[2][18] = 17; QNormal[3][18] = 17;
    QNormal[0][19] = 37; QNormal[1][19] = 37; QNormal[2][19] = 37; QNormal[3][19] = 37;
    QNormal[0][20] = 18; QNormal[1][20] = 18; QNormal[2][20] = 18; QNormal[3][20] = 18;
    QNormal[0][21] = 19; QNormal[1][21] = 19; QNormal[2][21] = 19; QNormal[3][21] = 19;
    QNormal[0][22] = 20; QNormal[1][22] = 20; QNormal[2][22] = 20; QNormal[3][22] = 20;
    QNormal[0][23] = 21; QNormal[1][23] = 21; QNormal[2][23] = 21; QNormal[3][23] = 21;
    QNormal[0][24] = 22; QNormal[1][24] = 22; QNormal[2][24] = 22; QNormal[3][24] = 22;
    QNormal[0][25] = 23; QNormal[1][25] = 23; QNormal[2][25] = 23; QNormal[3][25] = 23;
    QNormal[0][26] = 24; QNormal[1][26] = 24; QNormal[2][26] = 24; QNormal[3][26] = 24;
    QNormal[0][27] = 25; QNormal[1][27] = 25; QNormal[2][27] = 25; QNormal[3][27] = 25;
    QNormal[0][28] = 26; QNormal[1][28] = 26; QNormal[2][28] = 26; QNormal[3][28] = 26;
    QNormal[0][29] = 38; QNormal[1][29] = 38; QNormal[2][29] = 38; QNormal[3][29] = 38;
    QNormal[0][30] = 27; QNormal[1][30] = 27; QNormal[2][30] = 27; QNormal[3][30] = 27;
    QNormal[0][31] = 28; QNormal[1][31] = 28; QNormal[2][31] = 28; QNormal[3][31] = 28;
    QNormal[0][32] = 29; QNormal[1][32] = 29; QNormal[2][32] = 29; QNormal[3][32] = 29;
    QNormal[0][33] = 30; QNormal[1][33] = 30; QNormal[2][33] = 30; QNormal[3][33] = 30;
    QNormal[0][34] = 31; QNormal[1][34] = 31; QNormal[2][34] = 31; QNormal[3][34] = 31;
    QNormal[0][35] = 32; QNormal[1][35] = 32; QNormal[2][35] = 32; QNormal[3][35] = 32;
    QNormal[0][36] = 33; QNormal[1][36] = 33; QNormal[2][36] = 33; QNormal[3][36] = 33;
    QNormal[0][37] = 34; QNormal[1][37] = 34; QNormal[2][37] = 34; QNormal[3][37] = 34;
    QNormal[0][38] = 35; QNormal[1][38] = 35; QNormal[2][38] = 35; QNormal[3][38] = 35;
    QNormal[0][39] = 39; QNormal[1][39] = 39; QNormal[2][39] = 39; QNormal[3][39] = 39;

    Q1Special[0][0] = 1;   Q1Special[1][0] = 0;   Q1Special[2][0] = 1;   Q1Special[3][0] = 0;  
    Q1Special[0][1] = 8;   Q1Special[1][1] = 1;   Q1Special[2][1] = 8;   Q1Special[3][1] = 1;  
    Q1Special[0][2] = 3;   Q1Special[1][2] = 2;   Q1Special[2][2] = 3;   Q1Special[3][2] = 2;  
    Q1Special[0][3] = 4;   Q1Special[1][3] = 3;   Q1Special[2][3] = 4;   Q1Special[3][3] = 3;  
    Q1Special[0][4] = 35;  Q1Special[1][4] = 4;   Q1Special[2][4] = 35;  Q1Special[3][4] = 4;  
    Q1Special[0][5] = 0;   Q1Special[1][5] = 5;   Q1Special[2][5] = 0;   Q1Special[3][5] = 5;  
    Q1Special[0][6] = 7;   Q1Special[1][6] = 6;   Q1Special[2][6] = 7;   Q1Special[3][6] = 6;  
    Q1Special[0][7] = 14;  Q1Special[1][7] = 7;   Q1Special[2][7] = 14;  Q1Special[3][7] = 7;  
    Q1Special[0][8] = 9;   Q1Special[1][8] = 8;   Q1Special[2][8] = 9;   Q1Special[3][8] = 8;  
    Q1Special[0][9] = 36;  Q1Special[1][9] = 36;  Q1Special[2][9] = 36;  Q1Special[3][9] = 36; 
    Q1Special[0][10] = 10; Q1Special[1][10] = 9;  Q1Special[2][10] = 10; Q1Special[3][10] = 9; 
    Q1Special[0][11] = 29; Q1Special[1][11] = 10; Q1Special[2][11] = 29; Q1Special[3][11] = 10;
    Q1Special[0][12] = 30; Q1Special[1][12] = 11; Q1Special[2][12] = 30; Q1Special[3][12] = 11;
    Q1Special[0][13] = 13; Q1Special[1][13] = 12; Q1Special[2][13] = 13; Q1Special[3][13] = 12;
    Q1Special[0][14] = 20; Q1Special[1][14] = 13; Q1Special[2][14] = 20; Q1Special[3][14] = 13;
    Q1Special[0][15] = 15; Q1Special[1][15] = 14; Q1Special[2][15] = 15; Q1Special[3][15] = 14;
    Q1Special[0][16] = 16; Q1Special[1][16] = 15; Q1Special[2][16] = 16; Q1Special[3][16] = 15;
    Q1Special[0][17] = 23; Q1Special[1][17] = 16; Q1Special[2][17] = 23; Q1Special[3][17] = 16;
    Q1Special[0][18] = 24; Q1Special[1][18] = 17; Q1Special[2][18] = 24; Q1Special[3][18] = 17;
    Q1Special[0][19] = 37; Q1Special[1][19] = 37; Q1Special[2][19] = 37; Q1Special[3][19] = 37;
    Q1Special[0][20] = 19; Q1Special[1][20] = 18; Q1Special[2][20] = 19; Q1Special[3][20] = 18;
    Q1Special[0][21] = 26; Q1Special[1][21] = 19; Q1Special[2][21] = 26; Q1Special[3][21] = 19;
    Q1Special[0][22] = 21; Q1Special[1][22] = 20; Q1Special[2][22] = 21; Q1Special[3][22] = 20;
    Q1Special[0][23] = 22; Q1Special[1][23] = 21; Q1Special[2][23] = 22; Q1Special[3][23] = 21;
    Q1Special[0][24] = 17; Q1Special[1][24] = 22; Q1Special[2][24] = 17; Q1Special[3][24] = 22;
    Q1Special[0][25] = 18; Q1Special[1][25] = 23; Q1Special[2][25] = 18; Q1Special[3][25] = 23;
    Q1Special[0][26] = 25; Q1Special[1][26] = 24; Q1Special[2][26] = 25; Q1Special[3][26] = 24;
    Q1Special[0][27] = 32; Q1Special[1][27] = 25; Q1Special[2][27] = 32; Q1Special[3][27] = 25;
    Q1Special[0][28] = 27; Q1Special[1][28] = 26; Q1Special[2][28] = 27; Q1Special[3][28] = 26;
    Q1Special[0][29] = 38; Q1Special[1][29] = 38; Q1Special[2][29] = 38; Q1Special[3][29] = 38;
    Q1Special[0][30] = 28; Q1Special[1][30] = 27; Q1Special[2][30] = 28; Q1Special[3][30] = 27;
    Q1Special[0][31] = 11; Q1Special[1][31] = 28; Q1Special[2][31] = 11; Q1Special[3][31] = 28;
    Q1Special[0][32] = 12; Q1Special[1][32] = 29; Q1Special[2][32] = 12; Q1Special[3][32] = 29;
    Q1Special[0][33] = 31; Q1Special[1][33] = 30; Q1Special[2][33] = 31; Q1Special[3][33] = 30;
    Q1Special[0][34] = 2;  Q1Special[1][34] = 31; Q1Special[2][34] = 2;  Q1Special[3][34] = 31;
    Q1Special[0][35] = 33; Q1Special[1][35] = 32; Q1Special[2][35] = 33; Q1Special[3][35] = 32;
    Q1Special[0][36] = 34; Q1Special[1][36] = 33; Q1Special[2][36] = 34; Q1Special[3][36] = 33;
    Q1Special[0][37] = 5;  Q1Special[1][37] = 34; Q1Special[2][37] = 5;  Q1Special[3][37] = 34;
    Q1Special[0][38] = 6;  Q1Special[1][38] = 35; Q1Special[2][38] = 6;  Q1Special[3][38] = 35;
    Q1Special[0][39] = 39; Q1Special[1][39] = 39; Q1Special[2][39] = 39; Q1Special[3][39] = 39;

    Q2Special[0][0] = 0;   Q2Special[1][0] = 0;   Q2Special[2][0] = 0;   Q2Special[3][0] = 0;  
    Q2Special[0][1] = 1;   Q2Special[1][1] = 1;   Q2Special[2][1] = 1;   Q2Special[3][1] = 1;  
    Q2Special[0][2] = 2;   Q2Special[1][2] = 2;   Q2Special[2][2] = 2;   Q2Special[3][2] = 2;  
    Q2Special[0][3] = 3;   Q2Special[1][3] = 3;   Q2Special[2][3] = 3;   Q2Special[3][3] = 3;  
    Q2Special[0][4] = 4;   Q2Special[1][4] = 4;   Q2Special[2][4] = 4;   Q2Special[3][4] = 4;  
    Q2Special[0][5] = 5;   Q2Special[1][5] = 5;   Q2Special[2][5] = 5;   Q2Special[3][5] = 5;  
    Q2Special[0][6] = 6;   Q2Special[1][6] = 6;   Q2Special[2][6] = 6;   Q2Special[3][6] = 6;  
    Q2Special[0][7] = 7;   Q2Special[1][7] = 7;   Q2Special[2][7] = 7;   Q2Special[3][7] = 7;  
    Q2Special[0][8] = 8;   Q2Special[1][8] = 8;   Q2Special[2][8] = 8;   Q2Special[3][8] = 8;  
    Q2Special[0][9] = 36;  Q2Special[1][9] = 36;  Q2Special[2][9] = 36;  Q2Special[3][9] = 36; 
    Q2Special[0][10] = 9;  Q2Special[1][10] = 9;  Q2Special[2][10] = 9;  Q2Special[3][10] = 9; 
    Q2Special[0][11] = 10; Q2Special[1][11] = 10; Q2Special[2][11] = 10; Q2Special[3][11] = 10;
    Q2Special[0][12] = 11; Q2Special[1][12] = 11; Q2Special[2][12] = 11; Q2Special[3][12] = 11;
    Q2Special[0][13] = 12; Q2Special[1][13] = 12; Q2Special[2][13] = 12; Q2Special[3][13] = 12;
    Q2Special[0][14] = 13; Q2Special[1][14] = 13; Q2Special[2][14] = 13; Q2Special[3][14] = 13;
    Q2Special[0][15] = 14; Q2Special[1][15] = 14; Q2Special[2][15] = 14; Q2Special[3][15] = 14;
    Q2Special[0][16] = 15; Q2Special[1][16] = 15; Q2Special[2][16] = 15; Q2Special[3][16] = 15;
    Q2Special[0][17] = 22; Q2Special[1][17] = 16; Q2Special[2][17] = 22; Q2Special[3][17] = 16;
    Q2Special[0][18] = 17; Q2Special[1][18] = 17; Q2Special[2][18] = 17; Q2Special[3][18] = 17;
    Q2Special[0][19] = 37; Q2Special[1][19] = 37; Q2Special[2][19] = 37; Q2Special[3][19] = 37;
    Q2Special[0][20] = 18; Q2Special[1][20] = 18; Q2Special[2][20] = 18; Q2Special[3][20] = 18;
    Q2Special[0][21] = 19; Q2Special[1][21] = 19; Q2Special[2][21] = 19; Q2Special[3][21] = 19;
    Q2Special[0][22] = 20; Q2Special[1][22] = 20; Q2Special[2][22] = 20; Q2Special[3][22] = 20;
    Q2Special[0][23] = 21; Q2Special[1][23] = 21; Q2Special[2][23] = 21; Q2Special[3][23] = 21;
    Q2Special[0][24] = 16; Q2Special[1][24] = 22; Q2Special[2][24] = 16; Q2Special[3][24] = 22;
    Q2Special[0][25] = 23; Q2Special[1][25] = 23; Q2Special[2][25] = 23; Q2Special[3][25] = 23;
    Q2Special[0][26] = 24; Q2Special[1][26] = 24; Q2Special[2][26] = 24; Q2Special[3][26] = 24;
    Q2Special[0][27] = 25; Q2Special[1][27] = 25; Q2Special[2][27] = 25; Q2Special[3][27] = 25;
    Q2Special[0][28] = 26; Q2Special[1][28] = 26; Q2Special[2][28] = 26; Q2Special[3][28] = 26;
    Q2Special[0][29] = 38; Q2Special[1][29] = 38; Q2Special[2][29] = 38; Q2Special[3][29] = 38;
    Q2Special[0][30] = 27; Q2Special[1][30] = 27; Q2Special[2][30] = 27; Q2Special[3][30] = 27;
    Q2Special[0][31] = 28; Q2Special[1][31] = 28; Q2Special[2][31] = 28; Q2Special[3][31] = 28;
    Q2Special[0][32] = 29; Q2Special[1][32] = 29; Q2Special[2][32] = 29; Q2Special[3][32] = 29;
    Q2Special[0][33] = 30; Q2Special[1][33] = 30; Q2Special[2][33] = 30; Q2Special[3][33] = 30;
    Q2Special[0][34] = 31; Q2Special[1][34] = 31; Q2Special[2][34] = 31; Q2Special[3][34] = 31;
    Q2Special[0][35] = 32; Q2Special[1][35] = 32; Q2Special[2][35] = 32; Q2Special[3][35] = 32;
    Q2Special[0][36] = 33; Q2Special[1][36] = 33; Q2Special[2][36] = 33; Q2Special[3][36] = 33;
    Q2Special[0][37] = 34; Q2Special[1][37] = 34; Q2Special[2][37] = 34; Q2Special[3][37] = 34;
    Q2Special[0][38] = 35; Q2Special[1][38] = 35; Q2Special[2][38] = 35; Q2Special[3][38] = 35;
    Q2Special[0][39] = 39; Q2Special[1][39] = 39; Q2Special[2][39] = 39; Q2Special[3][39] = 39;
}

/**************************************************************/

void GRETINA::Reset() {
  g3Temp.clear();
   
  g1X.Clear(); g2X.Clear(); pt.Clear();
  g3X.Clear();

  g3out.Reset();  g3out.Clear();
  g2out.Reset();  g2out.Clear();
  g1out.Reset();  g1out.Clear();
  gSimOut.Reset();  gSimOut.Clear();
  b29.timestamp = 0;  b29.wfCFD = 0;
  b29.Clear();
};

/**************************************************************/

Int_t GRETINA::getMode1(FILE* inf, GRETINAVariables *gVar,
			counterVariables *cnt) {
  
  Int_t siz = 0;
  trackedGamma g1;
  Int_t remaining = 0;
  
  if (remaining == 0) {
    Int_t stuff[2];
    siz = fread(&stuff, sizeof(Int_t), 2, inf);
    cnt->nGammasThisHeader = stuff[0];
    cnt->nGammasRead = 0;
    cnt->Increment(2*sizeof(Int_t));
    remaining = 1;
  }
  
  while (remaining) {
    
    if (cnt->nGammasThisHeader > 0) {
      siz = fread(&g1, sizeof(struct trackedGamma), 1, inf);
      if (siz != 1) {
	cout << ALERTTEXT;
	printf("GetMode1(): Failed in bytes read.\n");
	raise(SIGINT);
	cout << RESET_COLOR;
      } else {
	cnt->nGammasRead++;
      }
      cnt->Increment(sizeof(struct trackedGamma));
    }
    
    g1X.Clear();
    
    /* Basic variable assignments to output structure */
    g1X.tracked = g1.tracked;
    g1X.FOM = g1.fom;
    g1X.nIntPts = g1.ndet;
    g1X.cc = g1.esum;
    g1X.xyzLab1.SetXYZ(g1.x0, g1.y0, g1.z0);
    g1X.xyzLab2.SetXYZ(g1.x1, g1.y1, g1.z1);
    g1X.e1 = g1.e0*1000.;
    g1X.e2 = g1.e1*1000.;
    g1X.timestamp = g1.timestamp;
    
    /* Doppler correction -- simple only right here */
    g1X.doppler = getDopplerSimple(g1X.xyzLab1, gVar->beta);
    
    g1out.gammas.push_back(g1X);
    
    if (cnt->nGammasRead == cnt->nGammasThisHeader) { remaining = 0; }
    else { remaining = 1; }
    
  }

  return 0;

}

/**************************************************************/

Int_t GRETINA::getMode2(FILE* inf, Int_t evtLength, GRETINAVariables *gVar, 
			counterVariables *cnt) {
 
  Int_t siz = 0; Int_t old = 0;
  mode2ABCD5678 g2; mode2ABCD1234 g2old; 

  Int_t remaining = 1;

  while (remaining) {
    siz = fread(&g2, sizeof(Int_t), 1, inf);
    if (g2.type == (Int_t)0xabcd1234) {
      old = 1;
      g2old.type = g2.type;
    }
    
    if (!old) {
      siz = fread(&g2.crystal_id, (sizeof(struct mode2ABCD5678) - 
				   sizeof(Int_t)), 1, inf);
      if (siz != 1) {
	cout << ALERTTEXT;
	printf("GetMode2(): Failed in bytes read.\n");
	raise(SIGINT);
	cout << RESET_COLOR;
      } else {
	cnt->Increment(sizeof(struct mode2ABCD5678));
	evtLength -= sizeof(struct mode2ABCD5678);
      }
    } else if (old) {
      siz = fread(&g2old.crystal_id, (sizeof(struct mode2ABCD1234) - 
				      sizeof(Int_t)), 1, inf);
      if (siz != 1) {
	cout << ALERTTEXT;
	printf("GetMode2(): Failed in bytes read.\n");
	raise(SIGINT);
	cout << RESET_COLOR;
      } else {
	cnt->Increment(sizeof(struct mode2ABCD1234));
	evtLength -= sizeof(struct mode2ABCD1234);
      }
    }
    
    if (!old) {
      g2X.Reset();
      g2X.crystalID = g2.crystal_id;
      /* Simple data transfers here... */
      g2X.timestamp = g2.timestamp;
      g2X.t0 = g2.t0;
      g2X.chiSq = g2.chisq;
      g2X.normChiSq = g2.norm_chisq;
      g2X.bl = g2.prestep; /* 2015-09-20 CMC: added for wf baseline and PZ checks */
      g2X.error = g2.pad;
      g2X.cc = g2.tot_e;
      
      /* Figure out the basics...what detector is this in terms of quads? 
	 We need to know this for the calibration...*/
      Int_t crystal = -1;
      for (Int_t index=0; index<MAXQUADS; index++) {
	if ((Int_t)(g2X.crystalID/4) == gVar->hole[index]) {
	  crystal = ((gVar->electronicsOrder[index]*4) + 
		     (Int_t)(g2X.crystalID%4));
	}
      }
      
      if (crystal != -1) { calibrateMode2CC(crystal, &g2, &g2X, gVar); }
      
      pt.Clear();
      for (Int_t m=0; m<MAX_INTPTS; m++) {
	if (g2.intpts[m].e != 0) {
	  pt.Clear();
	  pt.segNum = g2.intpts[m].seg;
	  pt.xyz.SetXYZ(g2.intpts[m].x,
			g2.intpts[m].y,
			g2.intpts[m].z);
	  pt.xyzLab = rot.crys2Lab(g2.crystal_id, pt.xyz);
	  pt.xyzLabSeg = rot.crys2Lab(g2.crystal_id, TVector3(xA[pt.segNum], 
							     yA[pt.segNum], 
							     zA[pt.segNum]));
	  pt.xyzLabCrys = rot.crys2Lab(g2.crystal_id, TVector3(0., 0., 0.));
	  pt.e = g2.intpts[m].e;
	  pt.segE = g2.intpts[m].seg_energy;
	  g2X.intpts.push_back(pt);
	}
      }
      
      analyzeMode2(&g2X, gVar);
      g2out.xtals.push_back(g2X);

    } else if (old) {
      g2X.Reset();
      g2X.crystalID = g2old.crystal_id;
      /* Simple data transfers here... */
      g2X.timestamp = g2old.timestamp;
      g2X.t0 = g2old.t0;
      g2X.chiSq = g2old.chisq;
      g2X.normChiSq = g2old.norm_chisq;
      g2X.error = g2old.pad; 

      g2X.cc = g2old.tot_e;
      
      pt.Clear();
      for (Int_t m=0; m<MAX_INTPTS; m++) {
	if (g2old.intpts[m].e != 0) {
	  pt.Clear();
	  pt.segNum = g2old.intpts[m].seg;
	  pt.xyz.SetXYZ(g2old.intpts[m].x,
			g2old.intpts[m].y,
			g2old.intpts[m].z);
	  pt.xyzLab = rot.crys2Lab(g2old.crystal_id, pt.xyz);
	  pt.e = g2old.intpts[m].e;
	  pt.segE = g2old.intpts[m].seg_energy;
	  g2X.intpts.push_back(pt);
	}
      }
      
      analyzeMode2(&g2X, gVar);
      g2out.xtals.push_back(g2X);   
    }
    
    if (evtLength > 0) { remaining = 1; }
    else { remaining = 0; }
  }
    
  return 0;
}

/**************************************************************/

/*! Function to read in a g4Simulated data structure -- this function
    is called from the GetData() function in Unpack.cpp when a G4SIM
    global header type is found
    
    \param inf A FILE* pointer to the data file to read in
    \return Returns an integer value of 0 if successful
*/

Int_t GRETINA::getSimulated(FILE* inf) {
 
  Int_t siz = 0;
  g4Sim_abcd1234 gSim;

  siz = fread(&gSim, sizeof(Int_t), 1, inf);
  if (gSim.type == (Int_t)0xabcd1234) {
    siz = fread(&gSim.num, (sizeof(struct g4Sim_abcd1234) - 
			    sizeof(Int_t)), 1, inf);
    if (siz != 1) {
      cout << ALERTTEXT;
      printf("GetSimulated(): Failed in bytes read.\n");
      raise(SIGINT);
      cout << RESET_COLOR;
    }
    gSimOut.sim.push_back(gSim);
  } else {
    cout << ALERTTEXT;
    printf("getSimulated(): Unexpected version for simulated data!\n");
    raise(SIGINT);
    cout << RESET_COLOR; fflush(stdout);
  }
  return 0;
}

/**************************************************************/

Int_t GRETINA::analyzeMode2(g2CrystalEvent *g2, GRETINAVariables* gVar) {

  /* Figure out the basics...what detector is this in terms of quads? */
  Int_t detectorFound = 0;  Int_t crystal = -1;
  for (Int_t index=0; index<MAXQUADS; index++) {
    if ((Int_t)(g2->crystalID/4) == gVar->hole[index]) {
      detectorFound = 1;
      crystal = ((gVar->electronicsOrder[index]*4) + 
   		 (Int_t)(g2->crystalID%4));
      g2->crystalNum = crystal+1; /* crystal starts from 0; crystalNum goes from 1 */
      g2->quadNum = index+1;
    }
  }

  /* If we didn't find that hole number in our settings, keep processing -- 
     but crystalNum gets set big (over edge of maximum of 120), and quadNum similarly 
     gets set to represent holeNumber + 30 */
  if (detectorFound == 0) {  
    g2->crystalNum = g2->crystalID + 120; 
    crystal = g2->crystalNum - 1;
    g2->quadNum = (g2->crystalNum/4) + 30; /* quadNum = holeNum + 30 */
  }

  /* Segment sum with cross-talk (Dino) correction applied -- can this be moved? */
  Float_t divisor = 0;
  Double_t sum = 0; Double_t netdino = 0;
  Float_t segmentSum = 0;
  
  vector<Int_t> hitID;  vector<Float_t> hitE;
  Int_t hit[40] = {0};
  
  UInt_t vecSize = g2->intpts.size();
  if (vecSize > 0) {
    for (UInt_t uj = 0; uj<vecSize; uj++) {
      if (g2->intpts[uj].segE > 10 && hit[g2->intpts[uj].segNum] == 0) {
	hit[g2->intpts[uj].segNum] = 1;
	hitID.push_back(g2->intpts[uj].segNum);
	hitE.push_back(g2->intpts[uj].segE);
      }
    }
  }
  UInt_t vecSizeID = hitID.size();

  /* Remember crystal goes from 0 */
  if (vecSizeID>0 && crystal >= 0 && crystal < MAXCRYSTALS) {
    for (UInt_t uj=0; uj<vecSizeID; uj++) {
      for (UInt_t uk=uj+1; uk<vecSizeID; uk++) {
	divisor += (gVar->dinoFactor[crystal][hitID[uj]][hitID[uk]]*
		    gVar->dinoFactor[crystal][hitID[uk]][hitID[uj]]);
      }
      netdino = 1;
      for (UInt_t um=0; um<vecSizeID; um++) {
	if (um != uj) {
	  netdino -= gVar->dinoFactor[crystal][hitID[uj]][hitID[um]];
	}
      }  
      sum += (hitE[uj])*netdino;
    }
  }
  g2->segSum = sum/(1-divisor);
  
  /* Doppler correction time -- right now, simple correction only.
     We'll check for particle information later, and recalculate Doppler
     if necessary. */
  
  if (g2->intpts.size() > 0) {
    g2->doppler = getDopplerSimple(g2->maxIntPtXYZLab(), gVar->beta);
    if ((g2->crystalNum)%2 == 0) { // Position 2 or 4: A type (--> crystalNum numbers starting at 1)
      g2->dopplerSeg = getDopplerSimple(rot.crys2Lab(g2->crystalID,
						     TVector3(xA[g2->maxIntPtSegNum()], 
							      yA[g2->maxIntPtSegNum()], 
							      zA[g2->maxIntPtSegNum()])),
					gVar->beta);
    } else if ((g2->crystalNum)%2 == 1) { // Position 1 or 3: B type
      g2->dopplerSeg = getDopplerSimple(rot.crys2Lab(g2->crystalID,
						     TVector3(xB[g2->maxIntPtSegNum()], 
							      yB[g2->maxIntPtSegNum()], 
							      zB[g2->maxIntPtSegNum()])),
					gVar->beta);
    }
    g2->dopplerCrystal = getDopplerSimple (rot.crys2Lab(g2->crystalID,
							TVector3(0., 0., 0.)), gVar->beta);
    
  }
  
  return (0);
}

/**************************************************************/

/*! Calculates the Doppler correction factor based on the velocity
    (beta) and the gamma-ray angle only. 
    
    \param xyz The TVector3 position vector of the first interaction
           point -- used to calculate the angle of gamma-ray emission
    \param beta The float value of the velocity of the moving emitting
           particle, in units of c
    \return Retuns a float value corresponding to the Doppler correction
            factor
*/

Float_t GRETINA::getDopplerSimple(TVector3 xyz, Float_t beta) {
  Float_t cosDop = xyz.CosTheta();
  Float_t gamma = 1/TMath::Sqrt(1. - beta*beta);
  return (gamma*(1. - beta*cosDop));
}

/**************************************************************/

void GRETINA::calibrateMode2CC(Int_t crystal, mode2ABCD5678 *g2, g2CrystalEvent *g2crystal, GRETINAVariables *gVar) {
  TRandom *random = new TRandom1();
  for (Int_t i=0; i<4; i++) {
    Float_t tmpE = (Float_t)g2->core_e[i] + random->Rndm() - 0.5;
    Int_t idNum;
    if (i == 0) { 
      idNum = (crystal)*40 + 9; 
      g2crystal->cc1 = (Float_t)(tmpE*gVar->ehiGeGain[idNum] + gVar->ehiGeOffset[idNum]);
    }  else if (i == 1) { 
      idNum = (crystal)*40 + 19; 
      g2crystal->cc2 = (Float_t)(tmpE*gVar->ehiGeGain[idNum] + gVar->ehiGeOffset[idNum]);
    } else if (i == 2) { 
      idNum = (crystal)*40 + 29; 
      g2crystal->cc3 = (Float_t)(tmpE*gVar->ehiGeGain[idNum] + gVar->ehiGeOffset[idNum]);
    } else if (i == 3) { 
      idNum = (crystal)*40 + 39;  
      g2crystal->cc4 = (Float_t)(tmpE*gVar->ehiGeGain[idNum] + gVar->ehiGeOffset[idNum]);
    }
  }
  random->Delete();
}

/**************************************************************/

void GRETINA::printSegCenters() {
  for (Int_t i=0; i<36; i++) {
    printf("%d:  A: (%0.3f, %0.3f, %0.3f)  B: (%0.3f, %0.3f, %0.3f)\n",
	   i, xA[i], yA[i], zA[i],xB[i], yB[i],zB[i]);
  }
}

/**************************************************************/

Int_t GRETINA::fillShell2Track() {

  /* Transfer mode2 data to shell structure for tracking. */
  Int_t i, j, k, l;
  Float_t sum, segSum[36], ff;

  /* Clear cluster structure... */
  for (j=0; j<MAXCLUSTERHITS; j++) { track.clust[j].Clear(); }
  /* Clear shell structure... */
  for (j=0; j<MAXSHELLHITS; j++) {
    track.shell.t0[j] = -1.0;  
    track.shell.chiSq[j] = -1.0; track.shell.normChiSq[j] = -1.0;
    track.shell.timestamp[j] = 0; track.shell.crystalID[j] = -1;
    track.shell.eDet[j] = -1.0;  track.shell.eSum[j] = -1.0; 
    track.shell.origPos[j] = -1; 
    track.shell.detNum[j] = -1;  track.shell.module[j] = -1; 
    track.shell.crystalType[j] = -1;
  }

  for (i=0; i<g2out.crystalMult(); i++) {
    if (g2out.xtals[i].error != 0) {
      return (8);
    }
  }
  
  /* # of mode2 events */
  track.shell.knownNumClusters = g2out.crystalMult();
  track.shell.nHit = 0;
  k=0;

  /* Loop over all mode2 events */
  for (i=0; i<g2out.crystalMult(); i++) {

    /* Modify energies based on tracking options... */
    if (track.ctrl.useSegEnergy) {
      for (j=0; j<36; j++) { segSum[j] = 0; }
      for (j=0; j<g2out.xtals[i].numIntPts(); j++) {
	if (g2out.xtals[i].intpts[j].e > 0) {
	  segSum[g2out.xtals[i].intpts[j].segNum] += g2out.xtals[i].intpts[j].e;
	}
      }
    } else if (track.ctrl.useCCEnergy) {
      sum = 0;
      for (j=0; j<g2out.xtals[i].numIntPts(); j++) {
	if (g2out.xtals[i].intpts[j].e > 0) {
	  sum += g2out.xtals[i].intpts[j].e;
	}
      }
    }

    /* Loop over interaction points */
    for (j=0; j<g2out.xtals[i].numIntPts(); j++) { 
      track.shell.knownClusterNum[k] = i;
      track.shell.t0[k] = g2out.xtals[i].t0;
      track.shell.chiSq[k] = g2out.xtals[i].chiSq;
      track.shell.normChiSq[k] = g2out.xtals[i].normChiSq;
      track.shell.timestamp[k] = g2out.xtals[i].timestamp;
      track.shell.t0[k] = g2out.xtals[i].t0;
      track.shell.crystalID[k] = g2out.xtals[i].crystalID;
 
      track.shell.xyz[k] = rot.crys2Lab(g2out.xtals[i].crystalID, 
					TVector3(g2out.xtals[i].intpts[j].xyz.X(),
						 g2out.xtals[i].intpts[j].xyz.Y(),
						 g2out.xtals[i].intpts[j].xyz.Z()));
 
      /* Tracking requires cm, not mm */
      track.shell.xyz[k].SetXYZ(track.shell.xyz[k].X()/10., 
				track.shell.xyz[k].Y()/10., 
				track.shell.xyz[k].Z()/10.);
 
      /* Tracking requires energy in MeV, not keV */
      if (track.ctrl.useSegEnergy) {
	ff = g2out.xtals[i].intpts[j].e / segSum[g2out.xtals[i].intpts[j].segNum];
	track.shell.eDet[k] = ((g2out.xtals[i].intpts[j].segE * ff) / 1000.);
      } else if (track.ctrl.useCCEnergy) {
	track.shell.eDet[k] = ((g2out.xtals[i].cc * 
				(g2out.xtals[i].intpts[j].e / sum)) / 1000.);
      }
      track.shell.eSum[k] = g2out.xtals[i].cc / 1000.;
      
      track.shell.origPos[k] = j;
      track.shell.crystalType[k] = (track.shell.crystalID[k] & 0x0003);
      track.shell.module[k] = (track.shell.crystalID[k] & 0xfffc) >> 2;
      track.shell.detNum[k] = ((4 * track.shell.module[k]) + track.shell.crystalType[k]);
      k++;
    }
  }
  track.shell.nHit = k;

  /* Done */
  return (0);
}

/**************************************************************/

void GRETINA::fillMode1(Int_t trackStatus, GRETINAVariables *gVar) {
  Int_t iClust, nmTS;
  Int_t i, j, ng = 0;
  Int_t writeOut[1000];
  long long int mTS[1000];
  
  for (iClust = 0; iClust < track.nClusters; iClust++) {

    /* Clusters need to qualify to be transferred as gammas. */
    writeOut[iClust] = 0;
    mTS[iClust] = (long long int) 0;
    nmTS = 0;
    
    /* Loop and find valid gammas */
    if (track.clust[iClust].valid) {
      track.stat.gammasForOutput++;
      if (track.clust[iClust].tracked &&
	  track.clust[iClust].nDet > 0 &&
	  track.clust[iClust].nDet < MAXNUMDET) {
	
	track.stat.gammasForOutputOK1++;
	
	/* Find the mean time for this gamma ray */
	for (j=0; j<track.clust[iClust].nDet; j++) {
	  mTS[iClust] += track.clust[iClust].intpts[j].timestamp;
	  nmTS++;
	}
	if (nmTS > 1) { mTS[iClust] /= (long long int) nmTS; }
	
	/* Mark for writeout if it qualifies */
	if (nmTS > 0 && mTS[iClust] > (long long int) 0) {
	  writeOut[iClust] = 1;
	  ng++;
	  track.stat.gammasForOutputOK2++;
	}
      }
    }
  }

  /* If there was a tracking problem, suppress the output */
  if (trackStatus != 0) { ng = 0; }
  
  /* Put the tracked gammas into output */
  g1GammaEvent g1;
  if (ng > 0) {
    for (iClust = 0; iClust < track.nClusters; iClust++) {
      if (writeOut[iClust]) {
	g1.cc = track.clust[iClust].eSum * 1000;
	g1.nIntPts = track.clust[iClust].nDet;
	g1.FOM = track.clust[iClust].FOM;
	g1.tracked = track.clust[iClust].tracked;
	/* Find and extract the first and second interactions */
	for (j=0; j<track.clust[iClust].nDet; j++) {
	  if (track.clust[iClust].intpts[j].order == 0) {
	    g1.xyzLab1.SetXYZ(track.clust[iClust].intpts[j].xyz.X()*10,
			      track.clust[iClust].intpts[j].xyz.Y()*10,
			      track.clust[iClust].intpts[j].xyz.Z()*10);
	    g1.e1 = track.clust[iClust].intpts[j].eDet * 1000;
	    g1.timestamp = track.clust[iClust].intpts[j].timestamp;
	    g1.t0 = track.clust[iClust].intpts[j].t0;
	  }
	  if (track.clust[iClust].intpts[j].order == 1) {
	    g1.xyzLab2.SetXYZ(track.clust[iClust].intpts[j].xyz.X()*10,
			      track.clust[iClust].intpts[j].xyz.Y()*10,
			      track.clust[iClust].intpts[j].xyz.Z()*10);
	    g1.e2 = track.clust[iClust].intpts[j].eDet * 1000;
	  }
	}

	/* Doppler correction... simple only here */
	g1.doppler = getDopplerSimple(g1.xyzLab1, gVar->beta);

	g1out.gammas.push_back(g1);
      }
    }
  }

}

/**************************************************************/

void GRETINA::getScaler(FILE *inf, Int_t evtLength) {
  Int_t siz = fread(scalerBuf, 1, evtLength, inf);
  if (siz != evtLength) {
    cout << ALERTTEXT;
    printf("GRETINA: getScaler(): Failed in bytes read.\n");
    cout << RESET_COLOR;  fflush(stdout);
  }
  
  printf("GRETINA: getScaler(): Buffer found:\n");
  for (Int_t i=0; i<(evtLength/4); i++) {
    printf("  --> [%d]: %d\b", i, scalerBuf[i]);
  }
}

/**************************************************************/

void GRETINA::checkSPIntegrity() {
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    sp.netSeg[i] = -1;
    sp.mult[i] = 0;
    sp.crystalBuild[i] = 0;
    sp.ccE[i] = -1;
    sp.segE[i] = -1;
    for (Int_t j=0; j<40; j++) {
      sp.segEventIndex[i][j] = -1;
      sp.segsNet[i][j] = -1;
    }
  }

  UInt_t vecSize = g3Temp.size();
  for (UInt_t i=0; i<vecSize; i++) {
    Int_t xtalNum = (Int_t)(g3Temp[i].ID/40);
    Int_t segNum = sp.map[xtalNum][(Int_t)(g3Temp[i].ID%40)];
  
    sp.segEventIndex[xtalNum][segNum] = i;

    /* The trace was pulled and put into the waves array when
       the GRETINA data was first unpacked. */
    
    /* Adjust the trace baseline offset */
    Int_t s = 0;
    for (Int_t b=0; b<25; b++) {
      s += sp.waves[xtalNum][segNum][b];
    }
    if (s >= 0) { s = (s+7)/25; }
    else { s = (s-7)/25; }
    for (Int_t b=0; b<sp.trLength; b++) {
      sp.waves[xtalNum][segNum][b] -= s;
    }
    
    /* Check for net energy */
    Int_t avg = 0;
    for (Int_t b = sp.trLength-10; b<sp.trLength; b++) {
      avg += TMath::Abs(sp.waves[xtalNum][segNum][b]);
    }
    avg /= 10;
    
    if (avg > 20 && g3Temp[i].eCal > 30 && segNum < 36) { /* Net */
      if (g3Temp[i].eCal > sp.segE[xtalNum]) {
	sp.segE[xtalNum] = g3Temp[i].eCal;
	sp.netSeg[xtalNum] = segNum;
      }
      sp.segs[xtalNum][sp.mult[xtalNum]] = segNum;
      sp.mult[xtalNum]++;
    } else if (avg <= 20 && segNum < 36 && g3Temp[i].eCal > 30) {
      if (0) {
	cout << ALERTTEXT;
	printf("GRETINA: checkSPIntegrity(): FPGA reports net, trace disagrees: %d\n",
	       g3Temp[i].ID);
	printf("          eCal: %0.3f, TraceAvg: %d\n", g3Temp[i].eCal, avg);
	cout << RESET_COLOR;  fflush(stdout);
      }
    }
    
    /* Pull out the CC energy -- use CC on last digitizer, totally arbitrary */
    if (g3Temp[i].ID%40 == 39) {
      sp.ccE[xtalNum] = g3Temp[i].eCal;
    }

    /* Build up the hit pattern */
    sp.crystalBuild[xtalNum] += (segNum + 1);
  }
}

/**************************************************************/

Int_t GRETINA::getMode3(FILE *inf, Int_t evtLength, counterVariables *cnt,
			controlVariables *ctrl, GRETINAVariables *gVar) {

  Int_t siz = 0, remaining = 0;
  mode3DataPacket *dp;
  
  siz = fread(gBuf, evtLength, 1, inf);
  if (siz != 1) {
    cout << ALERTTEXT;
    printf("getMode3(): Error in read attempt (A).  Aborting now...\n");
    cout << RESET_COLOR;  fflush(stdout);
    return 5;
  }
  cnt->Increment(evtLength);
  
  /* Byte swapping, due to little/big endian problem */
  for (Int_t j=0; j<evtLength; j=j+2) {
    swap(*(gBuf + j), *(gBuf + j + 1));
  }
  cnt->mode3i = 0;
  
  remaining = 1;
  
  while (remaining) {
    unsigned char *tmp = (gBuf);
    tmp = (gBuf + cnt->mode3i*2);
    
    /* Allocate memory... */
    if ( !(dp = (mode3DataPacket*)malloc(sizeof(dp->aahdr) +
					 sizeof(dp->hdr) + 
					 sizeof(dp->waveform))) ) {
      cout << ALERTTEXT;
      printf("getMode3(): Failed in memory allocation.\n");
      cout << RESET_COLOR;  fflush(stdout);
      exit(-1);
    }
    memset(dp->waveform, 1, MAX_TRACE_LENGTH * sizeof(UShort_t));
    
    /* Copy 'AAAA' header */
    memmove(&dp->aahdr[0], tmp, sizeof(dp->hdr) + sizeof(dp->aahdr));
    if ((dp->aahdr[0] != 0xAAAA) || (dp->aahdr[1] != 0xAAAA)) {
      cout << ALERTTEXT;
      printf("getMode3(): Didn't get 'AAAA' header as expected!\n");
      printf("getMode3(): Found this instead: %x %x\n", dp->aahdr[0], dp->aahdr[1]);
      cout << RESET_COLOR;  fflush(stdout);
      exit(-1);
    }
    
    /* We've got the data packet, pull out information */
    g3ch.Clear();
    
    /* Interpret the header information */
    g3ch.hdr0 = dp->hdr[0];  g3ch.hdr1 = dp->hdr[1];
    g3ch.hdr7 = dp->hdr[7];

    Int_t module = g3ch.module();
    Int_t channel = g3ch.chanID();
    Int_t sign = g3ch.sign();
    Int_t TL = g3ch.tracelength();

    cnt->mode3i += (sizeof(dp->aahdr) + sizeof(dp->hdr)) / 2;
    tmp = (gBuf + cnt->mode3i*2);
    
    /* Copy the waveform! */
    memmove(&dp->waveform[0], tmp, TL*sizeof(UShort_t));
    
    cnt->mode3i += (TL * sizeof(UShort_t)) / 2;
    tmp = (gBuf + cnt->mode3i*2);
  
    g3ch.ID = -1;
    for (Int_t i=0; i<MAXQUADS; i++) {
      if (module/16 == gVar->hole[i]) {
	g3ch.ID = (module - (gVar->hole[i] - gVar->electronicsOrder[i])*4*4)*10 + channel;
      } 
    }
    if (g3ch.ID < 0) { g3ch.ID = module*10 + channel; } /* Need to think about this... */
    
    Int_t hiEnergy = 0;
    hiEnergy = (dp->hdr[7] & 0x00ff);
    UInt_t tmpEnergy = 0;  Int_t tmpIntEnergy = 0;
    tmpEnergy = ((UInt_t)(hiEnergy) << 16);
    tmpEnergy += dp->hdr[4];
    tmpIntEnergy = (Int_t)tmpEnergy;
    if (sign) {
      tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    } else {
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    }
    if (tmpIntEnergy == 65536) { /* Guard against weird FPGA energy anomoly */
      g3ch.eRaw = 0.; 
    } else { g3ch.eRaw = (Float_t)(tmpIntEnergy/32.); }
    
    hiEnergy = 0;  sign = 0;  tmpEnergy = 0;  tmpIntEnergy = 0;
    hiEnergy = (dp->hdr[11] & 0x00ff);
    sign = (dp->hdr[11] & 0x0100);
    tmpEnergy = ((UInt_t)(hiEnergy) << 16);
    tmpEnergy += dp->hdr[8];
    tmpIntEnergy = (Int_t)(tmpEnergy);
    if (sign) {
      tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    } else {
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    }
    g3ch.eCalPO = (Float_t)(tmpIntEnergy/32.);

    hiEnergy = 0;  sign = 0;  tmpEnergy = 0;  tmpIntEnergy = 0;
    hiEnergy = (dp->hdr[13] & 0x0001);
    sign = (dp->hdr[13] & 0x0002);
    tmpEnergy = ((UInt_t)(hiEnergy) << 23);
    tmpEnergy += ((UInt_t)(dp->hdr[10]) << 7);
    tmpEnergy += ((UInt_t)(dp->hdr[11] & 0xfe00) >> 9);
    tmpIntEnergy = (Int_t)(tmpEnergy);
    if (sign) {
      tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    } else {
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    }
    g3ch.prevE1 = (Float_t)(tmpIntEnergy/32.);
    
    hiEnergy = 0;  sign = 0;  tmpEnergy = 0;  tmpIntEnergy = 0;
    hiEnergy = (dp->hdr[12] & 0x03ff);
    sign = (dp->hdr[12] & 0x0400);
    tmpEnergy = ((UInt_t)(hiEnergy) << 14);
    tmpEnergy += ((UInt_t)(dp->hdr[13] & 0xfffc) >> 2);
    tmpIntEnergy = (Int_t)(tmpEnergy);
    if (sign) {
      tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    } else {
      if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
      }
    }
    g3ch.prevE2 = (Float_t)(tmpIntEnergy/32.);
    g3ch.PZrollover = ((UInt_t)(dp->hdr[12] & 0xf800) >> 11);

    // printf("Real:   HDR 4  7: %#08x %#08x\n", dp->hdr[4], dp->hdr[7]);
    // printf("        Chn: %d Sign: %d, tmpIntEnergy: %d\n", g3ch.chanID(), g3ch.sign(), tmpIntEnergy);
    // printf("        eRaw: %0.3f  pileUp: %d\n", g3ch.eRaw, g3ch.pileUp());

    // Bool_t sign2 = 0;
    // Int_t tmpIntEnergy2 = (Int_t)(32*g3ch.eRaw);
    // UInt_t tmpEnergy2 = 0;
    // if (g3ch.chanID()%10 == 9) { /* CC */
    //   if (g3ch.eRaw <= 0) {
    // 	tmpEnergy2 = (tmpIntEnergy2 + (UInt_t)(0x01000000));
    // 	sign2 = 1;
    //   } else {
    // 	tmpEnergy2 = (UInt_t)(tmpIntEnergy2);
    //   }
    // } else {
    //   if (g3ch.eRaw <= 0) {
    // 	tmpEnergy2 = (UInt_t)(-tmpIntEnergy2);
    //   } else {
    // 	tmpEnergy2 = (UInt_t)((0x01000000)-tmpIntEnergy2);
    // 	sign2 = 1;
    //   }
    // }
    // UShort_t hdr4 = (tmpEnergy & 0xffff);
    // UShort_t hdr7 = ((tmpEnergy & 0x00ff0000) >> 16);
    // if (sign2) { hdr7 = (hdr7 | 0x0100); }
    // hdr7 += (g3ch.pileUp() << 15);

    // printf("Calculated HDR 4  7: %#08x %#08x\n", hdr4, hdr7);
    // printf("           Sign2: %d\n", sign2);
    // cin.get();

    if (!ctrl->superPulse) { calibrateMode3(&g3ch, gVar); }
    else { calibrateMode3SP(&g3ch); }
    
    /* Transform the waveform, if needed */
    if (ctrl->withWAVE || channel%10 == 9 || g3ch.eRaw > 100) {
      g3ch.wf.raw.clear();

      for (Int_t j=0; j<TL+1; j=j+2) {
	if (dp->waveform[j+1] & 0x8000) {
	  g3ch.wf.raw.push_back(dp->waveform[j+1] - std::numeric_limits<unsigned int>::max());
	} else {
	  g3ch.wf.raw.push_back(dp->waveform[j+1]);
	}
	if (dp->waveform[j] & 0x8000) {
	  g3ch.wf.raw.push_back(dp->waveform[j] - std::numeric_limits<unsigned int>::max());
	} else {
	  g3ch.wf.raw.push_back(dp->waveform[j]);
	} 
      }

      if (channel%10 != 9) { 
	for (Int_t j=0; j<g3ch.wf.raw.size(); j++) {
	  g3ch.wf.raw[j] = -1*g3ch.wf.raw[j];
	}
      }
      
      /* For the CC always get a baseline value from the minimum trace, which is 6 samples. */
      if (g3ch.wf.raw.size() >= 6) {  g3ch.baseline = g3ch.wf.BL(0, 6);  }

      g3ch.calcTime = g3ch.wf.CFD(0);

    }
    
    if (ctrl->withWAVE) {
      if (ctrl->superPulse) {
	for (Int_t j=0; j<sp.trLength; j++) {
	  sp.waves[(Int_t)(g3ch.ID/40)][sp.map[(Int_t)(g3ch.ID/40)][(Int_t)(g3ch.ID%40)]][j] = g3ch.wf.raw[j];
	}
      }
    }
    
    g3ch.timestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[3])) + 
				  ((ULong64_t)(dp->hdr[2]) << 16) +
				  ((ULong64_t)(dp->hdr[5]) << 32) );
    g3ch.CFDtimestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[6])) + 
				     ((ULong64_t)(dp->hdr[9]) << 16) +
				     ((ULong64_t)(dp->hdr[8]) << 32) );
    g3ch.deltaT1 = (UShort_t)(dp->hdr[6]);
    g3ch.deltaT2 = (UShort_t)(dp->hdr[9]);
    cnt->lastBdTS[(Int_t)(g3ch.ID/10)] = g3ch.timestamp;

    if (!ctrl->withWAVE) {
      g3ch.wf.Clear();
    }
    
    g3Temp.push_back(g3ch);
    
    free(dp);

    if ( (Int_t)(cnt->mode3i*2) == evtLength ) { remaining = 0; }
    else if ( (Int_t)(cnt->mode3i*2) < evtLength ) { remaining = 1; }
  }
  
  return (0);
  
}

/**************************************************************/

void GRETINA::calibrateMode3(g3ChannelEvent *g3, GRETINAVariables *gVar) {

  Float_t tmpE = (g3->eRaw)*0.25;
  g3->eCal = (tmpE * gVar->ehiGeGain[g3->ID] + gVar->ehiGeOffset[g3->ID]);
  tmpE = (g3->eCalPO)*0.25;
  g3->eCalPO = (tmpE * gVar->ehiGeGain[g3->ID] + gVar->ehiGeOffset[g3->ID]);
  tmpE = (g3->prevE1)*0.25;
  g3->prevE1 = (tmpE * gVar->ehiGeGain[g3->ID] + gVar->ehiGeOffset[g3->ID]);
  tmpE = (g3->prevE2)*0.25;
  g3->prevE2 = (tmpE * gVar->ehiGeGain[g3->ID] + gVar->ehiGeOffset[g3->ID]);

}

/**************************************************************/

void GRETINA::calibrateMode3SP(g3ChannelEvent *g3) {
  
  Float_t tmpE = (g3->eRaw)*0.25;
  g3->eCal = (tmpE * sp.ehiGeGain[g3->ID] + sp.ehiGeOffset[g3->ID]);

}

/**************************************************************/

void GRETINA::analyzeMode3(GRETINAVariables *gVar, controlVariables *ctrl) {
  
  if (g3Temp.size() > 0) {
    
    for (UInt_t ui=0; ui<g3Temp.size(); ui++) {
      
      Int_t newCrystal = 0, found = 0;
      UInt_t xIndex = 1000;

      g3X.Clear();
      
      Int_t xid = g3Temp[ui].ID/40; /* xid numbers from 0 */
      Int_t cid = g3Temp[ui].ID%40;

      if (xid >= MAXCRYSTALS || xid < 0) {
	cout << ALERTTEXT;
	printf("\n -------> Oops, you've got an unexpected crystal ID (hole %i).", 
	       g3Temp[ui].module()/16);
	printf("\n          Did you edit the \"gretina.set\" file?\n\n");
	cout << RESET_COLOR << endl;
	fflush(stdout);
	exit(0);      
      }

      /* Simple assignments first */
      if (g3out.crystalMult() == 0) {
	newCrystal = 1;
	g3X.crystalNum = xid; /* xid still numbering from 0, 
				 so crystalNum at this point starts at 0 */
      } else {
	for (UInt_t uj=0; uj<g3out.crystalMult(); uj++) {
	  if (g3out.xtals[uj].crystalNum == xid) {
	    xIndex = uj;
	    found = 1;
	  }
	}
	if (!found) { newCrystal = 1;  g3X.crystalNum = xid; }
      }
      if (newCrystal) {      
	g3X.module = g3Temp[ui].module();
	g3out.xtals.push_back(g3X);
	xIndex = g3out.crystalMult() - 1;
      }
      g3out.xtals[xIndex].chn.push_back(g3Temp[ui]);
    }

    g3Temp.clear();

  }

  /*------ Crystal-by-crystal standard analysis (the main analysis). ------*/
  
  vector<Int_t> hitID;  vector<Double_t> hitE;
  vector<Int_t> hitIDT;  vector<Double_t> hitET;
  vector<Int_t> hitIDT2;  vector<Double_t> hitET2;

  Double_t ccCalcTime;

  for (UInt_t ui = 0; ui<g3out.crystalMult(); ui++) {

    if (g3out.xtals[ui].module/4 != cloverModuleNumber) {
      g3out.xtals[ui].OrderChannels();
      g3out.xtals[ui].crystalNum += 1; /* And now crystalNum goes from 1 */
      
      for (Int_t qN=0; qN<MAXQUADS; qN++) {
	if ((g3out.xtals[ui].crystalNum-1) < (gVar->electronicsOrder[qN]+1)*4 &&
	    (g3out.xtals[ui].crystalNum-1) >= (gVar->electronicsOrder[qN]*4)) {
	  g3out.xtals[ui].quadNum = qN + 1; /* quadNum goes from 1 as well... */
	}
      }

      /* Something to realize -- anything not corresponding to something in the
	 gretina.set file gets assigned quadNum = 0... fix this maybe? */
      
      /* Assign segment numbers properly */
      if (g3out.xtals[ui].quadNum == 1) {
	for (Int_t i=0; i<g3out.xtals[ui].chn.size(); i++) {
	  g3out.xtals[ui].chn[i].segNum = Q1Special[(Int_t)(g3out.xtals[ui].crystalNum-1)%4][g3out.xtals[ui].chn[i].chnNum()]; 
	}
      } else if (g3out.xtals[ui].quadNum == 2) {
	for (Int_t i=0; i<g3out.xtals[ui].chn.size(); i++) {
	  g3out.xtals[ui].chn[i].segNum = Q2Special[(Int_t)(g3out.xtals[ui].crystalNum-1)%4][g3out.xtals[ui].chn[i].chnNum()]; 
	}
      } else {
	for (Int_t i=0; i<g3out.xtals[ui].chn.size(); i++) {
	  g3out.xtals[ui].chn[i].segNum = QNormal[(Int_t)(g3out.xtals[ui].crystalNum-1)%4][g3out.xtals[ui].chn[i].chnNum()]; 
	  if (g3out.xtals[ui].chn[i].segNum == 36) { ccCalcTime = g3out.xtals[ui].chn[i].calcTime; }
	}
      }
      
      Float_t divisor = 0;  Float_t sum = 0;  Float_t netdino = 0;
      
      /* Segment sum analysis */
      
      hitID.clear();  hitE.clear();
      for (Int_t i=0; i<g3out.xtals[ui].chn.size(); i++) {
	if (g3out.xtals[ui].chn[i].eCal > 25) {
	  if (g3out.xtals[ui].chn[i].chnNum()%10 != 9) {
	    hitID.push_back(g3out.xtals[ui].chn[i].segNum);
	    hitE.push_back(g3out.xtals[ui].chn[i].eCal);
	    if (ctrl->withWAVE) {
	      if (((g3out.xtals[ui].chn[i].calcTime >= wfMinCrossTime) &&
		   (g3out.xtals[ui].chn[i].calcTime <= wfMaxCrossTime)) ||
		  (g3out.xtals[ui].chn[i].calcTime == 1)) {
		hitIDT.push_back(g3out.xtals[ui].chn[i].segNum);
		hitET.push_back(g3out.xtals[ui].chn[i].eCal);
	      } 
	      if ((((ccCalcTime - g3out.xtals[ui].chn[i].calcTime) > -3) &&
		   ((ccCalcTime - g3out.xtals[ui].chn[i].calcTime) < 5))) {
		hitIDT2.push_back(g3out.xtals[ui].chn[i].segNum);
		hitET2.push_back(g3out.xtals[ui].chn[i].eCal);
	      } 
	    }
	  } /* end of segment requirement */
	} /* end of threshold check */
      } /* end of loop over channels */     
      
      /*--- Basic sum, no cuts. ---*/
      Float_t maxE = 0;
      divisor = 0.0; sum = 0.0; netdino = 0.0;
      UInt_t vecSizeID = hitID.size();      
      if (vecSizeID>0) {
	for (UInt_t uj=0; uj<vecSizeID; uj++) {
	  for (UInt_t uk=uj+1; uk<vecSizeID; uk++) {
	    divisor += (gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitID[uj]][hitID[uk]]*
			gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitID[uk]][hitID[uj]]);
	  }
	  netdino = 1.;
	  for (UInt_t um=0; um<vecSizeID; um++) {
	    if (um != uj) { netdino -= gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitID[uj]][hitID[um]]; }
	  }  
	  sum += (hitE[uj])*netdino;
	}
      }
      g3out.xtals[ui].segSum = sum / (1-divisor);
      
      hitID.clear();
      hitE.clear();
      
      /*--- Sum with segment timing cut. ---*/
      if (ctrl->withWAVE) {
	divisor = 0.0; sum = 0.0; netdino = 0.0;  
	UInt_t vecSizeIDT = hitIDT.size();
	if (vecSizeIDT>0) {
	  for (UInt_t uj=0; uj<vecSizeIDT; uj++) {
	    for (UInt_t uk=uj+1; uk<vecSizeIDT; uk++) {
	      divisor += (gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitIDT[uj]][hitIDT[uk]]*
			  gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitIDT[uk]][hitIDT[uj]]);
	    }
	    netdino = 1;
	    for (UInt_t um=0; um<vecSizeIDT; um++) {
	      if (um != uj) { netdino -= gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitIDT[uj]][hitIDT[um]]; }
	    }  
	    sum += (hitET[uj])*netdino;
	  }
	}
      }
      g3out.xtals[ui].segSumTimed = sum / (1-divisor);
      
      hitIDT.clear();
      hitET.clear();
          
      /*--- Sum with segment timing cut 2. ---*/
      if (ctrl->withWAVE) {
	divisor = 0.0; sum = 0.0; netdino = 0.0;  
	UInt_t vecSizeIDT2 = hitIDT2.size();
	if (vecSizeIDT2>0) {
	  for (UInt_t uj=0; uj<vecSizeIDT2; uj++) {
	    for (UInt_t uk=uj+1; uk<vecSizeIDT2; uk++) {
	      divisor += (gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitIDT2[uj]][hitIDT2[uk]]*
			  gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitIDT2[uk]][hitIDT2[uj]]);
	    }
	    netdino = 1;
	    for (UInt_t um=0; um<vecSizeIDT2; um++) {
	      if (um != uj) { netdino -= gVar->dinoFactor[(g3out.xtals[ui].crystalNum - 1)][hitIDT2[uj]][hitIDT2[um]]; }
	    }  
	    sum += (hitET2[uj])*netdino;
	  }
	}
      }
      g3out.xtals[ui].segSumTimed2 = sum / (1-divisor);
      
      hitIDT2.clear();
      hitET2.clear();


      /* Now some real physics -- Doppler correction.  We use the center of 
	 the max E hit segment, could easily use center of crystal. */
      TVector3 xyzL;
      if (g3out.xtals[ui].module > 0 && g3out.xtals[ui].module/16 <= 30) {
	if (g3out.xtals[ui].maxSegNum() < 0 || g3out.xtals[ui].maxSegNum() > 35) { /* No segment! use crystal origin */
	  xyzL = rot.crys2Lab(g3out.xtals[ui].module/4, TVector3(0., 0., 0.));
	} else if ((g3out.xtals[ui].crystalNum)%2 == 0) { /* Position 2 or 4: A type 
							     (--> remember crystalNum starts at 1 ) */
	  xyzL = rot.crys2Lab(g3out.xtals[ui].module/4, TVector3(xA[g3out.xtals[ui].maxSegNum()],
								 yA[g3out.xtals[ui].maxSegNum()],
								 zA[g3out.xtals[ui].maxSegNum()]));
	  
	} else if ((g3out.xtals[ui].crystalNum)%2 == 1) { /* Position 1 or 3: B type */
	  xyzL = rot.crys2Lab(g3out.xtals[ui].module/4, TVector3(xB[g3out.xtals[ui].maxSegNum()],
								 yB[g3out.xtals[ui].maxSegNum()],
								 zB[g3out.xtals[ui].maxSegNum()]));
	}
      }
      g3out.xtals[ui].dopplerSeg = getDopplerSimple(xyzL, gVar->beta);
      
      if (g3out.xtals[ui].module > 0 && g3out.xtals[ui].module/4 <= 30) {
	xyzL = rot.crys2Lab(g3out.xtals[ui].module/4, TVector3(0., 0., 0.));
	g3out.xtals[ui].dopplerCrystal = getDopplerSimple(xyzL, gVar->beta);
      }
      
    } else { /* Clover module ID -- relevant to coincidence scanning only at this point */

      cloverOut.Reset();

      /* Map mode3 channel data into clover output structure */
      for (Int_t m=0; m<g3out.xtals[ui].chn.size(); m++) {
	cloverX.Clear();
	if (g3out.xtals[ui].chn[m].eRaw > 0) {
	  cloverX.eRaw = g3out.xtals[ui].chn[m].eRaw;
	  cloverX.eCal = ( (cloverSlope[g3out.xtals[ui].chn[m].chnNum()] * cloverX.eRaw) +
			   cloverOffset[g3out.xtals[ui].chn[m].chnNum()] );
	  cloverX.ID = g3out.xtals[ui].chn[m].chnNum() - 29;
	  if (ctrl->withWAVE) { cloverX.CFD = g3out.xtals[ui].chn[m].wf.CFD(10); }
	  if (cloverX.ID > 4) { cloverX.ID -= 1; }
	  if (cloverX.ID < 0) { cloverX.ID = g3out.xtals[ui].chn[m].chnNum() - 9; }
	  cloverX.TS = g3out.xtals[ui].chn[m].timestamp;
	  cloverOut.xtals.push_back(cloverX);
	}
      }

      g3out.xtals.erase(g3out.xtals.begin() + ui);

    }

  } /* Loop over hit crystals */
    
}

gHistos::gHistos() { ; }

gHistos::~gHistos() { ; }

void gHistos::writeHistos(Int_t ctrl) {
  if (ctrl == 1) { /* calibration for energies */
    for (Int_t m=0; m<MAXCHANNELS; m++) {
      if (eRaw[m]) { eRaw[m]->Write();  eRaw[m]->Delete(); }
      if (eCal[m]) { eCal[m]->Write();  eCal[m]->Delete(); }
    }
  } else if (ctrl == 2) { /* Cross-talk calibrations */
    for (Int_t m=0; m<MAXCHANNELS; m++) {
      if (eRaw[m]) { eRaw[m]->Write();  eRaw[m]->Delete(); }
      if (eCal[m]) { eCal[m]->Write();  eCal[m]->Delete(); }
    }
    for (Int_t m=0; m<36; m++) {
      for (Int_t n=m; n<36; n++) {
	if (mult2Pair[n][m]) { mult2Pair[n][m]->Write(); mult2Pair[n][m]->Delete(); }
      }
    }
  }
  
}

void GRETINA::fillHistos(Int_t ctrl) {

  gHist.xtalkThresh = 20.;

  if (ctrl == 1) { /* ctrl->calibration, so energy histos */
    for (UInt_t um=0; um<g3Temp.size(); um++) {
      if (g3Temp[um].ID < MAXCHANNELS) { 
	if (!gHist.eRaw[g3Temp[um].ID]) {
	  printf("Need histogram with index: %d... building those now.\n", g3Temp[um].ID);
	  char title[300];
	  sprintf(title, "eRaw%d", g3Temp[um].ID);
	  gHist.eRaw[g3Temp[um].ID] = new TH1F(title, title, 30000, 0, 60000);
	  sprintf(title, "eCal%d", g3Temp[um].ID);
	  gHist.eCal[g3Temp[um].ID] = new TH1F(title, title, 5000, 0, 5000);
	} else { ; }
	gHist.eRaw[g3Temp[um].ID]->Fill(g3Temp[um].eRaw*0.25); 
	gHist.eCal[g3Temp[um].ID]->Fill(g3Temp[um].eCal); 
      }
    }
  } else if (ctrl == 2) { /* Cross-talk calibrations */
    vector<Int_t> segNum;
    vector<Float_t> segE;
    float_t ccE;
    for (UInt_t um=0; um<g3out.xtals[0].chn.size(); um++) {
      if (g3out.xtals[0].chn[um].ID < MAXCHANNELS) { 
   	if (!gHist.eRaw[g3out.xtals[0].chn[um].ID]) {
	  printf("Need histogram with index: %d... building those now.\n", g3out.xtals[0].chn[um].ID);
	  char title[300];
	  sprintf(title, "eRaw%d", g3out.xtals[0].chn[um].ID);
	  gHist.eRaw[g3out.xtals[0].chn[um].ID] = new TH1F(title, title, 30000, 0, 60000);
	  sprintf(title, "eCal%d", g3out.xtals[0].chn[um].ID);
	  gHist.eCal[g3out.xtals[0].chn[um].ID] = new TH1F(title, title, 5000, 0, 5000);
	} else { ; }
	gHist.eRaw[g3out.xtals[0].chn[um].ID]->Fill(g3out.xtals[0].chn[um].eRaw*0.25); 
	gHist.eCal[g3out.xtals[0].chn[um].ID]->Fill(g3out.xtals[0].chn[um].eCal); 
      
	if (g3out.xtals[0].chn[um].eCal > gHist.xtalkThresh && g3out.xtals[0].chn[um].segNum<36) {
	  segNum.push_back(g3out.xtals[0].chn[um].segNum);
	  segE.push_back(g3out.xtals[0].chn[um].eCal);
	}
      }
      if (g3out.xtals[0].chn[um].ID == 1459) { ccE = g3out.xtals[0].chn[um].eCal; }
    }
    if (segNum.size() == 2 && ((ccE > 1170 && ccE < 1176) || (ccE > 1330 && ccE < 1335))) { 
      if (segNum[0] < segNum[1]) { 
	Int_t tmpNum = segNum[0]; segNum[0] = segNum[1]; segNum[1] = tmpNum;
	Float_t tmpE = segE[0]; segE[0] = segE[1]; segE[1] = tmpE;
      }
      if (!gHist.mult2Pair[segNum[0]][segNum[1]]) {
	char title[300];
	sprintf(title, "mult2Pair_%d_%d", segNum[0], segNum[1]);
	printf("Need pair histogram: %s\n", title);
	gHist.mult2Pair[segNum[0]][segNum[1]] = new TH1F(title, title, 5000, 0, 5000);
      }
      gHist.mult2Pair[segNum[0]][segNum[1]]->Fill(segE[0]+segE[1]);
    }
  }
}

void cloverCrystalEvent::Clear() {
  eRaw = 0.0; eCal = 0.0; ID = 0; TS = 0;
}

void cloverOUT::Reset() {
  xtals.clear();
};

UInt_t cloverOUT::cloverMult(Float_t thresh) {
  Int_t mult = 0;
  for (UInt_t i=0; i<xtals.size(); i++) {
    if (xtals[i].eCal > thresh && xtals[i].ID < 10) { mult++; }
  }
  return mult;
}

Int_t cloverOUT::maxCloverID() {
  Float_t maxE = 0; Int_t maxID = 0;
  for (UInt_t i=0; i<xtals.size(); i++) {
    if (xtals[i].eCal > maxE && xtals[i].ID < 10) {
      maxE = xtals[i].eCal;
      maxID = xtals[i].ID;
    }
  }
  return maxID;
}

Float_t cloverOUT::maxCloverE() {
  Float_t maxE = 0; Int_t maxID = 0;
  for (UInt_t i=0; i<xtals.size(); i++) {
    if (xtals[i].eCal > maxE && xtals[i].ID < 10) {
      maxE = xtals[i].eCal;
      maxID = xtals[i].ID;
    }
  }
  return maxE;
}

Bool_t cloverOUT::middleA() {
  
  Float_t top = 0., bottom = 0., middle = 0.;

  if (maxCloverID() <= 4) { /* Crystal A! */ 
    for (UInt_t i=0; i<xtals.size(); i++) {
      if (xtals[i].ID == 14) { top = xtals[i].eRaw; }
      if (xtals[i].ID == 15) { middle = xtals[i].eRaw; }
      if (xtals[i].ID == 16) { bottom = xtals[i].eRaw; }
    }
  }
  
  if (middle > top && middle > bottom) { return 1; }
  else { return 0; }

}

Bool_t cloverOUT::middleB() {
  
  Float_t top = 0., bottom = 0., middle = 0.;

  if (maxCloverID() > 4) { /* Crystal B! */ 
    for (UInt_t i=0; i<xtals.size(); i++) {
      if (xtals[i].ID == 11) { top = xtals[i].eRaw; }
      if (xtals[i].ID == 12) { middle = xtals[i].eRaw; }
      if (xtals[i].ID == 13) { bottom = xtals[i].eRaw; }
    }
  }
  
  if (top < 100. && bottom < 100.) { return 1; }
  else { return 0; }

}
