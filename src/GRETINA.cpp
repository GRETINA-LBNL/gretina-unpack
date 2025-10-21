#include "GRETINA.h"

ClassImp(globalHeader);

ClassImp(mode3DataPacket);


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
/* g3OUT Class Functions **************************************/
/**************************************************************/

/*! Resets the mode3 output class -- clears channel vector

    \return No return value -- clears data structures in class
            directly
*/

void g3OUT::Reset() {
  UInt_t banks = bankMult();
  for (UInt_t ui=0; ui<banks; ui++) {
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

UInt_t g3OUT::bankMult() { return xtals.size(); }


/**************************************************************/
/* GRETINA Class Functions ************************************/
/**************************************************************/

void GRETINA::Initialize() {

  wfMinCrossTime = 70.; wfMaxCrossTime = 90.;
     
}

/**************************************************************/

void GRETINA::Reset() {
  g3Temp.clear();
  g3X.Clear();
  g3out.Reset();  g3out.Clear();
  g3H.past.clear(); g3H.Clear();

  b88.chn.clear();
  b88.timestamp = 0;  b88.wfCFD = 0;
  b88.Clear();
};

/**************************************************************/

Int_t GRETINA::getMode3(FILE *inf, Int_t evtLength, counterVariables *cnt,
			controlVariables *ctrl) {

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

  /* Byte swapping, due to little/big endian mismatch */
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
    g3ch.hdr0 = dp->hdr[0];
    g3ch.hdr1 = dp->hdr[1];
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
  
    g3ch.ID = module*10 + channel;  

    /* Extract energy information (peak find) */
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

    /* Pick-off energy extraction */
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
    g3ch.eRawPO = (Float_t)(tmpIntEnergy/32.);

    /* Last previous energy extraction */
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

    /* Second last previous energy extraction */
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

    /* Transform the waveform, if needed */
    if (ctrl->withWAVE) {
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
       
    g3ch.timestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[3])) + 
				  ((ULong64_t)(dp->hdr[2]) << 16) +
				  ((ULong64_t)(dp->hdr[5]) << 32) );
    g3ch.CFDtimestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[6])) + 
				     ((ULong64_t)(dp->hdr[9]) << 16) +
				     ((ULong64_t)(dp->hdr[8]) << 32) );
    g3ch.deltaT1 = (UShort_t)(dp->hdr[6]);
    g3ch.deltaT2 = (UShort_t)(dp->hdr[9]);
    cnt->lastBdTS[(Int_t)(g3ch.ID/10)] = g3ch.timestamp;

    g3Temp.push_back(g3ch);
    
    free(dp);

    if ( (Int_t)(cnt->mode3i*2) == evtLength ) { remaining = 0; }
    else if ( (Int_t)(cnt->mode3i*2) < evtLength ) { remaining = 1; }
  }
  
  return (0);
  
}

/**************************************************************/


void GRETINA::analyzeMode3(controlVariables *ctrl) {
  
  if (g3Temp.size() > 0) {
    
    for (UInt_t ui=0; ui<g3Temp.size(); ui++) {
      
      Int_t newBank = 0, found = 0;
      UInt_t xIndex = 1000;

      g3X.Clear();
      
      Int_t xid = g3Temp[ui].ID/40; /* xid numbers from 0 */
      Int_t cid = g3Temp[ui].ID%40;

      /* Simple assignments first */
      if (g3out.bankMult() == 0) {
	newBank = 1;
	g3X.bankNum = xid; /* xid still numbering from 0, 
			      so crystalNum at this point starts at 0 */
      } else {
	for (UInt_t uj=0; uj<g3out.bankMult(); uj++) {
	  if (g3out.xtals[uj].bankNum == xid) {
	    xIndex = uj;
	    found = 1;
	  }
	}
	if (!found) { newBank = 1;  g3X.bankNum = xid; }
      }
      if (newBank) {      
	g3X.module = g3Temp[ui].module();
	g3out.xtals.push_back(g3X);
	xIndex = g3out.bankMult() - 1;
      }
      g3out.xtals[xIndex].chn.push_back(g3Temp[ui]);
    }

    g3Temp.clear();

  }

  for (UInt_t ui = 0; ui<g3out.bankMult(); ui++) {

    g3out.xtals[ui].OrderChannels();
    g3out.xtals[ui].bankNum += 1; /* And now bankNum goes from 1 */
      
  } /* Loop over hit crystals */
  
}


Int_t GRETINA::getMode3History(FILE *inf, Int_t evtLength, long long int hTS, counterVariables *cnt) {

  Int_t siz = 0, remaining = 0;
  mode3HistoryPacket *dp;

  siz = fread(gBuf, evtLength, 1, inf);
  if (siz != 1) {
    cout << ALERTTEXT;
    printf("getMode3History(): Error in read attempt (A).  Aborting now...\n");
    cout << RESET_COLOR;  fflush(stdout);
    return 5;
  }
  cnt->Increment(evtLength);
  
  /* Byte swapping, due to little/big endian problem */
  for (Int_t j=0; j<evtLength; j=j+2) {
    swap(*(gBuf + j), *(gBuf + j + 1));
  }
  
  unsigned char *tmp = (gBuf);
  
  /* Allocate memory... */
  if ( !(dp = (mode3HistoryPacket*)malloc(sizeof(dp->aahdr) +
					  sizeof(dp->hdr) + 
					  sizeof(dp->data))) ) {
    cout << ALERTTEXT;
    printf("getMode3History(): Failed in memory allocation.\n");
    cout << RESET_COLOR;  fflush(stdout);
    exit(-1);
  }
  memset(dp->data, 1, MAX_TRACE_LENGTH * sizeof(UShort_t));
  
  /* Copy 'AAAA' header */
  memmove(&dp->aahdr[0], tmp, sizeof(dp->aahdr));
  if ((dp->aahdr[0] != 0xAAAA) || (dp->aahdr[1] != 0xAAAA)) {
    cout << ALERTTEXT;
    printf("getMode3History(): Didn't get 'AAAA' header as expected!\n");
    printf("getMode3History(): Found this instead: %x %x\n", dp->aahdr[0], dp->aahdr[1]);
    cout << RESET_COLOR;  fflush(stdout);
    exit(-1);
  }
  
  tmp = (gBuf + (sizeof(dp->aahdr)));
  
  /* Now copy the rest of the event */  
  memmove(&dp->hdr[0], tmp, (evtLength - sizeof(dp->aahdr)));

  /* We've got the data, pull out the information now... */
  if ( (dp->hdr[1] & 0xf) == 0xB) {
    Int_t eventsize = (dp->hdr[0] & 0x3ff);

    gH.energy = 0.;
    gH.TS = 0;
    gH.module = 0;

    // 2016-07-23 CMC added module to gH to differentiate between digitizers
    // unlike energy and TS, module should not be reset within the while loop
    // a mode3 channel event/packet comes from one digitzer, channel 9 by firmware
    // Int_t module() { return (hdr1 >> 4); }
    gH.module = (dp->hdr[1]) >> 4;
    
    gH.TS = (ULong64_t)( ((ULong64_t)(dp->hdr[3])) +  
			 ((ULong64_t)(dp->hdr[2]) << 16) +
			 ((ULong64_t)(dp->hdr[5] & 0x7fff) << 32) +
			 ((ULong64_t)(hTS & 0x800000000000)) );

    long long int headerTS = gH.TS;

    Int_t overflow = (dp->hdr[5] & 0x8000);
    gH.energy = dp->hdr[4] + ((dp->hdr[7] & 0xff)<< 16);
    gH.energy /= 32.;
       
    eventsize -= (sizeof(dp->hdr)/4);

    g3H.past.push_back(gH);
    gH.energy = 0; gH.TS = 0; gH.BLpreSum = 0;
        
    int i=0;
    while (eventsize) {
     
      /**************************************************/
      /* Pairs of 32-bit words:                         */
      /* 0:   31 -- 25   24 -- 0                        */
      /*       TS(6-0)   E(24-0)                        */
      /* 1:    31        30 -- 0                        */
      /*    Overflow     TS(37-7)                       */
      /* 2:   31 -- 18   17 -- 0                        */
      /*         ??     BLpre-sum                       */
      /**************************************************/

      long long int TSintermediate = ((dp->data[i+3]) + ((dp->data[i+2] & 0x7fff) << 16));

      gH.TS = (ULong64_t)( ((ULong64_t)(dp->data[i] & 0xfe00) >> 9) + 
			   ((ULong64_t)(TSintermediate) << 7) + 
			   ((ULong64_t)(headerTS & 0xffc000000000)) );   
      gH.energy = dp->data[i+1] + ((dp->data[i] & 0x1ff) << 16);
      gH.energy /= 32;
      overflow = (dp->data[i+2] & 0x8000);
      gH.BLpreSum = dp->data[i+5];

      i+=6;
      eventsize -= (sizeof(unsigned short)*3/2);
      g3H.past.insert(g3H.past.end(), gH);
      gH.energy = 0; gH.TS = 0; gH.BLpreSum = 0;
    }
  }

  free(dp);

  return (0);
  
}

Int_t GRETINA::getBank88(FILE *inf, Int_t evtLength, counterVariables *cnt) {

  Int_t siz = 0, remaining = 0;
  mode3DataPacket *dp;
  
  siz = fread(gBuf, evtLength, 1, inf);
  if (siz != 1) {
    cout << ALERTTEXT;
    printf("getBank29(): Error in read attempt (A).  Aborting now...\n");
    cout << RESET_COLOR;  fflush(stdout);
    return 5;
  }
  cnt->Increment(evtLength);
  
  /* Byte swapping, due to little/big endian problem */
  for (Int_t j=0; j<evtLength; j=j+2) {
    swap(*(gBuf + j), *(gBuf + j + 1));
  }
  cnt->b88i = 0;
  
  remaining = 1;
  
  while (remaining) {
    unsigned char *tmp = (gBuf);
    tmp = (gBuf + cnt->b88i*2);
    
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
      printf("getBank29(): Didn't get 'AAAA' header as expected!\n");
      printf("getBank29(): Found this instead: %x %x\n", dp->aahdr[0], dp->aahdr[1]);
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

    cnt->b88i += (sizeof(dp->aahdr) + sizeof(dp->hdr)) / 2;
    tmp = (gBuf + cnt->b88i*2);
    
    /* Copy the waveform! */
    memmove(&dp->waveform[0], tmp, TL*sizeof(UShort_t));
    
    cnt->b88i += (TL * sizeof(UShort_t)) / 2;
    tmp = (gBuf + cnt->b88i*2);
  
    // g3ch.ID = channel;

    // Int_t hiEnergy = 0;
    // hiEnergy = (dp->hdr[7] & 0x00ff);
    // UInt_t tmpEnergy = 0;  Int_t tmpIntEnergy = 0;
    // tmpEnergy = ((UInt_t)(hiEnergy) << 16);
    // tmpEnergy += dp->hdr[4];
    // tmpIntEnergy = (Int_t)tmpEnergy;
    // if (sign) {
    //   tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // } else {
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // }
    // if (tmpIntEnergy == 65536) { /* Guard against weird FPGA energy anomoly */
    //   g3ch.eRaw = 0.; 
    // } else { g3ch.eRaw = (Float_t)(tmpIntEnergy/32.); }
    
    // hiEnergy = 0;  sign = 0;  tmpEnergy = 0;  tmpIntEnergy = 0;
    // hiEnergy = (dp->hdr[11] & 0x00ff);
    // sign = (dp->hdr[11] & 0x0100);
    // tmpEnergy = ((UInt_t)(hiEnergy) << 16);
    // tmpEnergy += dp->hdr[8];
    // tmpIntEnergy = (Int_t)(tmpEnergy);
    // if (sign) {
    //   tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // } else {
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // }
    // g3ch.eCalPO = (Float_t)(tmpIntEnergy/32.);

    // hiEnergy = 0;  sign = 0;  tmpEnergy = 0;  tmpIntEnergy = 0;
    // hiEnergy = (dp->hdr[13] & 0x0001);
    // sign = (dp->hdr[13] & 0x0002);
    // tmpEnergy = ((UInt_t)(hiEnergy) << 23);
    // tmpEnergy += ((UInt_t)(dp->hdr[10]) << 7);
    // tmpEnergy += ((UInt_t)(dp->hdr[11] & 0xfe00) >> 9);
    // tmpIntEnergy = (Int_t)(tmpEnergy);
    // if (sign) {
    //   tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // } else {
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // }
    // g3ch.prevE1 = (Float_t)(tmpIntEnergy/32.);
    
    // hiEnergy = 0;  sign = 0;  tmpEnergy = 0;  tmpIntEnergy = 0;
    // hiEnergy = (dp->hdr[12] & 0x03ff);
    // sign = (dp->hdr[12] & 0x0400);
    // tmpEnergy = ((UInt_t)(hiEnergy) << 14);
    // tmpEnergy += ((UInt_t)(dp->hdr[13] & 0xfffc) >> 2);
    // tmpIntEnergy = (Int_t)(tmpEnergy);
    // if (sign) {
    //   tmpIntEnergy = (Int_t)(tmpIntEnergy - (Int_t)0x01000000);
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // } else {
    //   if ( (Int_t)(channel%10) != 9 ) { /* Not a CC */
    // 	tmpIntEnergy = -(Int_t)(tmpIntEnergy);
    //   }
    // }
    // g3ch.prevE2 = (Float_t)(tmpIntEnergy/32.);
    // g3ch.PZrollover = ((UInt_t)(dp->hdr[12] & 0xf800) >> 11);

    // /* Transform the waveform */
    // g3ch.wf.raw.clear();
    // for (Int_t j=0; j<TL+1; j=j+2) {
    //   if (dp->waveform[j+1] & 0x8000) {
    // 	g3ch.wf.raw.push_back(dp->waveform[j+1] - std::numeric_limits<unsigned int>::max());
    //   } else {
    // 	g3ch.wf.raw.push_back(dp->waveform[j+1]);
    //   }
    //   if (dp->waveform[j] & 0x8000) {
    // 	g3ch.wf.raw.push_back(dp->waveform[j] - std::numeric_limits<unsigned int>::max());
    //   } else {
    // 	g3ch.wf.raw.push_back(dp->waveform[j]);
    //   } 
    // }

    // /* For the CC always get a baseline value from the minimum trace, which is 6 samples. */
    // if (g3ch.wf.raw.size() >= 6) {  g3ch.baseline = g3ch.wf.BL(0, 6);  }
    
    // g3ch.calcTime = g3ch.wf.CFD(0);

    g3ch.timestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[3])) + 
				  ((ULong64_t)(dp->hdr[2]) << 16) +
				  ((ULong64_t)(dp->hdr[5]) << 32) );
    g3ch.CFDtimestamp = (ULong64_t)( ((ULong64_t)(dp->hdr[6])) + 
				     ((ULong64_t)(dp->hdr[9]) << 16) +
				     ((ULong64_t)(dp->hdr[8]) << 32) );
    g3ch.deltaT1 = (UShort_t)(dp->hdr[6]);
    g3ch.deltaT2 = (UShort_t)(dp->hdr[9]);
   
    b88.timestamp = g3ch.timestamp; 
    b88.chn.push_back(g3ch);

    free(dp);

    if ( (Int_t)(cnt->b88i*2) == evtLength ) { remaining = 0; }
    else if ( (Int_t)(cnt->b88i*2) < evtLength ) { remaining = 1; }
  }
  
  return (0);
  
}
