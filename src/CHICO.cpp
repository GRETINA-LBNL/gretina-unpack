/*! \file CHICO.cpp
    \brief Parameter and functions for CHICO analysis.

    This file provides the analysis functions for CHICO analysis, based on 
    an analysis code CheckChico.cc from S. Zhu (ANL).

    Author: H.L. Crawford
    Date: April 29, 2014
*/

#include "CHICO.h"

void CHICORaw::Initialize() {
  for (Int_t i=0; i<128; i++) {
    cathode_tdc_val[i] = -1;
    cathode_tdc_ch[i] = 0;
    anode_tdc_val[i] = -1;
    anode_tdc_ch[i] = 0;
    if (i<32) {
      anode_qdc_val[i] = -1;
      anode_qdc_ch[i] = -1;
    }
  }
  LEDts = 0;
  cathode_tdc_num = 0;
  anode_tdc_num = 0;
  anode_qdc_num = 0;
  status = 0;
  RF = -1;

    
}

void CHICORaw::Reset() {
  Initialize();
}

void CHICOParticle::Initialize() {

  /* printf("CHICOParticle::Initialize() particle.fThetaL %f \n"); */

  id = -1;
  t = 0.0;
  tof = 0.0;
  mass = -1.0;
  thetaL = -1; phiL = -1;
  thetaR = -1; phiR = -1;
  fThetaL = -1.0; fPhiL = -1.0;
  fThetaR = -1.0; fPhiR = -1.0;
  eL = -1; eR = -1;
  rf = 0.0;
  /*  pgCosL = 0.0; pgCosR = 0.0;*/
  dL = -1.0; dR = -1.0;
  massT = -1.0; massP = -1.0; betaT = -1.0; betaP = -1.0; Qvalue = -999;
  pgCosT = 0.0; pgCosP = 0.0;
}

void CHICOParticle::Reset() {
  Initialize();
}

CHICOFull::CHICOFull() {
  gotParticle = 0;
  idParticle = 999;
}

void CHICOFull::Initialize() {
  gotParticle = 0;
  
  raw.Initialize();
  particle.Initialize();
}

void CHICOFull::Reset() {
  gotParticle = 0;

  raw.Reset();
  particle.Reset();
}

void CHICOFull::ReadThetaCalibration(TString thetaCalFile) {
  FILE *in;
  if ((in = fopen(thetaCalFile.Data(), "r")) == NULL) {
    printf("Input calibration file %s not found.", thetaCalFile.Data());
    return;
  }

  printf("Reading theta calibration file: %s ...", thetaCalFile.Data());

  char line[500]; 
  char* retval;
  Int_t nn = 0, i1;
  Float_t f1, f2, f3;

  retval = fgets(line, 500, in);
  while (retval != NULL) {
    if (line[0] == 35) { 
      /* '#' comment line, do nothing. */
    } else if (line[0] == 59) {
      /* ';' comment line, do nothing. */
    } else if (line[0] == 10) {
      /* Empty line, do nothing. */
    } else {
      sscanf(line, "%i %f %f %f", &i1, &f1, &f2, &f3);
      thetaOffset[nn] = f1; thetaGain[nn] = f2; thetaQuad[nn] = f3;
      nn++;
    }
    retval = fgets(line, 500, in);
  }
  printf("Done.\n");
  if (nn != PPAC_NUM) {
    perror("Unexpected number of theta calibration parameters read.");
    return;
  } else {
    printf("Theta calibration for %i PPACs loaded.\n", nn);
  }
    
  fclose(in);
}

void CHICOFull::ReadPhiCalibration(TString phiCalFile) {
  FILE *in;
  if ((in = fopen(phiCalFile.Data(), "r")) == NULL) {
    printf("Input calibration file %s not found.", phiCalFile.Data());
    return;
  }

  printf("Reading theta calibration file: %s ...", phiCalFile.Data());

  char line[500]; 
  char* retval;
  Int_t nn = 0, i1;
  Float_t f1, f2, f3;

  retval = fgets(line, 500, in);
  while (retval != NULL) {
    if (line[0] == 35) { 
      /* '#' comment line, do nothing. */
    } else if (line[0] == 59) {
      /* ';' comment line, do nothing. */
    } else if (line[0] == 10) {
      /* Empty line, do nothing. */
    } else {
      sscanf(line, "%i %f %f", &i1, &f1, &f2);
      phiOffset[nn] = f1; phiGain[nn] = f2;
      nn++;
    }
    retval = fgets(line, 500, in);
  }
  printf("Done.\n");
  if (nn != PPAC_NUM) {
    perror("Unexpected number of phi calibration parameters read.");
    return;
  } else {
    printf("Phi calibration for %i PPACs loaded.\n", nn);
  }
    
  fclose(in);

}

void CHICOFull::ReadBetaCalibration(TString betaCalFile) {
  FILE *in;
  if ((in = fopen(betaCalFile.Data(), "r")) == NULL) {
    printf("Input calibration file %s not found.", betaCalFile.Data());
    return;
  }

  printf("Reading beta calibration file: %s ...", betaCalFile.Data());

  char line[500]; 
  char* retval;
  Int_t nn = 0, i1, i;
  Float_t f1, f2, f3, f4, f5, f6;

  for(i=0;i<96;i++)betaB[i]=0.0;
  nn = 96;

  retval = fgets(line, 500, in);

  sscanf(line, "%f %f %f %f %f %f", &f1, &f2, &f3, &f4, &f5, &f6);
  retval = fgets(line, 500, in);

  betaP0=f1; betaP1=f2; betaP2=f3; betaP3=f4; betaP4=f5; betaP5=f6; 

  sscanf(line, "%f %f %f %f %f %f", &f1, &f2, &f3, &f4, &f5, &f6);
  retval = fgets(line, 500, in);

  betaT0=f1; betaT1=f2; betaT2=f3; betaT3=f4; betaT4=f5; betaT5=f6; 


  while (retval != NULL) {
    if (line[0] == 35) { 
      /* '#' comment line, do nothing. */
    } else if (line[0] == 59) {
      /* ';' comment line, do nothing. */
    } else if (line[0] == 10) {
      /* Empty line, do nothing. */
    } else {
      sscanf(line, "%f %f", &f1, &f2);
      betaB[nn] = f2;     
      nn++;
    }
    retval = fgets(line, 500, in);
  }
  printf("Done.\n");
  /* 2018-05-08 CMC replace
  if (nn != 160) {
  */
  if (nn != 169) {
    printf("nn= %d \n",nn );
    perror("Unexpected number of beta calibration parameters read.");
    return;
  } else {
    printf("Beta calibration loaded.\n");
  }

  for(i=169;i<180;i++)betaB[i]=0.0;
    
  fclose(in);
}

void CHICOFull::InitializeCHICOVariables(TString thetaCalFile,
					 TString phiCalFile,
					 TString betaCalFile) {
  ReadThetaCalibration(thetaCalFile);
  ReadPhiCalibration(phiCalFile);
  ReadBetaCalibration(betaCalFile);
  setupCut();
};

Int_t CHICOFull::getAndUnpackCHICO(FILE *inf, Int_t length) {

  Reset(); /* Resets CHICO data structures (raw + processed) */
  memset(rawCHICO, 0, sizeof(rawCHICO)); 

  Int_t siz = fread(rawCHICO, 1, length, inf);
  if (siz != length) {
    std::cout << ALERTTEXT;
    printf("GetCHICO(): Failed in bytes read.\n");
    std::cout << RESET_COLOR;  fflush(stdout);
  }

  UInt_t *p = rawCHICO;
 
  UShort_t evSize;
  UInt_t nextInt;
  UShort_t chan = 0; Int_t val = 0, refval = 0;

  Int_t i, j, k = 0;
  Int_t chCounter = 0;
  Int_t seenTrailer = 0;
  
  evSize = *p/4;
  raw.status = (UShort_t)((*p & 0xffff0000) >> 16); p++;
  raw.LEDts = (unsigned long long int)(*p & 0xffffffff); p++;
  raw.LEDts += ((unsigned long long int)(*p & 0xffff) << 32);

  nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
  p++;

  /***** Anode QDC *****/
  if (nextInt != 0xffffffff) {
    if (((nextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != QDCHEADER) return 0;
    if (((nextInt & QDCGEOMASK) >> QDCGEOSHIFT) != ANODE_E_VSN) return 0;
    chCounter = (nextInt & COUNTMASK) >> COUNTSHIFT;
    k = 0;
    for (i = 0; i<chCounter; i++) {
      nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
      p++;
      if (((nextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != DATA) return 0;
      chan = (UShort_t)((nextInt & QDCCHANMASK) >> QDCCHANSHIFT);
      val = (nextInt & QDCDATAMASK);
      if (chan < PPAC_NUM && val > 0) {
	raw.anode_qdc_ch[k] = chan;
	raw.anode_qdc_val[k] = val;
	k++;
      }
    }
    nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
    p++;
    if (((nextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != QDCTRAILER) return 0;
    nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
    p++;
    if (nextInt != 0xffffffff) {
      while (1) {
	nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
	p++;
	if (nextInt == 0xffffffff) break;
      }
    }
  }
  raw.anode_qdc_num = k;

  chCounter = 0;
  chan = 0; val = 0;
  nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
  p++;

  /***** Anode TDC *****/
  if (nextInt != 0xffffffff) {
    if ((nextInt & TDCTYPEMASK) != TDCHEADER) return 0;
    if ((nextInt & TDCGEOMASK) != ANODE_T_VSN) return 0;
    while (!seenTrailer) {
      nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
      p++;
      switch(nextInt & TDCTYPEMASK) {
      case DATA:
	chan = (UShort_t)((nextInt & TDCCHANMASK) >> TDCCHANSHIFT);
	val = (nextInt & TDCDATAMASK);
	if (chan != ANODE_REFCH && chan != RFCH) {
	  if (chan < PPAC_NUM && chCounter < 128) {
	    raw.anode_tdc_ch[chCounter] = chan;
	    raw.anode_tdc_val[chCounter] = val;
	    chCounter++;
	  } else if (chCounter >= 128) {
	    return 128;
	  }
	} else if (chan == RFCH) {
	  raw.RF = val;
	} else if (chan == ANODE_REFCH) {
	  refval = val;
	}
	break;
      case TDCTRAILER:
	seenTrailer = 1;
	break;
      default:
	break;
      }
    }
    // for (i=0; i<chCounter; i++) {
    //   raw.anode_tdc_val[i] -= refval;
    // }
    raw.RF -= refval;

    nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
    p++;
    if (nextInt != 0xffffffff) {
      multiAnodeTDCNum++;
      while (1) {
	nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
	p++;
	if (nextInt == 0xffffffff) break;
      }
    }
  }
  raw.anode_tdc_num = chCounter;

  chCounter = 0; seenTrailer = 0;
  chan = 0; val = 0; refval = 0;
  nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
  p++;

  /***** Cathode TDC *****/
  if (nextInt != 0xffffffff) {
    if ((nextInt & TDCTYPEMASK) != TDCHEADER) return 0;
    if ((nextInt & TDCGEOMASK) != CATHODE_T_VSN) return 0;
    while (!seenTrailer) {
      nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
      p++;
      switch(nextInt & TDCTYPEMASK) {
      case DATA:
	chan = (UShort_t)((nextInt & TDCCHANMASK) >> TDCCHANSHIFT);
	val = (nextInt & TDCDATAMASK);
	if (chan != CATHODE_REFCH) {
	  if (chan < PPAC_NUM*4 && chCounter < 128) {
	    raw.cathode_tdc_ch[chCounter] = chan;
	    raw.cathode_tdc_val[chCounter] = val;
	    chCounter++;
	  } else if (chCounter >= 128) {
	    return 128;
	  }
	} else if (chan == CATHODE_REFCH) {
	  refval = val;
	}
	break;
      case TDCTRAILER:
	seenTrailer = 1;
	break;
      default:
	break;
      }
    }
    nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
    p++;
    if (nextInt != 0xffffffff) {
      multiCathodeTDCNum++;
      while (1) {
	nextInt = ((*p & 0xffff0000) >> 16) + ((*(p+1) & 0xffff) << 16);
	p++;
	if (nextInt == 0xffffffff) break;
      }
    }
  }
  raw.cathode_tdc_num = chCounter;


  /*
  for (Int_t i=0; i<raw.cathode_tdc_num; i++) { printf("i %d raw.cathode_tdc_val[i] %d raw.cathode_tdc_ch[i] %d \n ", i , raw.cathode_tdc_val[i], raw.cathode_tdc_ch[i]); }
  for (Int_t i=0; i<raw.anode_tdc_num; i++) { printf("i %d raw.anode_tdc_val[i] %d raw.anode_tdc_ch[i] %d \n ", i , raw.anode_tdc_val[i], raw.anode_tdc_ch[i]); }
  //for (Int_t i=0; i<raw.anode_qdc_num; i++) { printf("i %d raw.anode_qdc_val[i] %d raw.anode_qdc_ch[i] %d \n ", i , raw.anode_qdc_val[i], raw.anode_qdc_ch[i]); }
  printf("\n");
  */

  GetParticle();
  idParticle=IDparticles();
  return 1;
}

void CHICOFull::GetParticle() {
  
  Float_t d = (Float_t) (RAND_MAX) + 1.0;
  Int_t id = 160;
  Int_t dT = 0;
  Int_t anodeID = 0; Int_t anodeL = 100, anodeR = 100;
 
  Int_t validL = 0, validR = 0;
  Int_t validTL = 0, validPL = 0;
  Int_t validTR = 0, validPR = 0;
  Int_t validT = 0, validP = 0;

  Int_t thetaL = 0, thetaR = 0;
  Int_t thetaLs, thetaRs;
  Int_t phiL = 0, phiR = 0;
  Int_t phiLs, phiRs;
  Float_t fThetaL = 0., fThetaR = 0.;
  Float_t fPhiL = 0., fPhiR = 0., fPhiL_R = 0., fPhiR_L = 0.;
  Float_t fPhiL_av = 0., fPhiR_av = 0.;;
  
  Float_t tof = 0.;

  UShort_t tmp1 = 0, tmp2 = 0;
  Int_t vtmp1 = 0, vtmp2 = 0, doubleOK;

  doubleOK = 0;
  if (raw.anode_tdc_num > 2) {
    for (Int_t i = 0; i < (raw.anode_tdc_num - 1); i++) {
      for (Int_t j = (i+1); j < raw.anode_tdc_num; j++) {
	if ((raw.anode_tdc_ch[i]%5) == (raw.anode_tdc_ch[j]%5)) {
	  tmp1 = raw.anode_tdc_ch[i];
	  tmp2 = raw.anode_tdc_ch[j];
	  vtmp1 = raw.anode_tdc_val[i];
	  vtmp2 = raw.anode_tdc_val[j];
	  doubleOK = 1;
	  /* printf("doubleOK i %d j %d \n",i,j); */
	}
      }
    }
  }

  if (doubleOK == 1) {
    raw.anode_tdc_num = 2;
    raw.anode_tdc_ch[0] = tmp1;
    raw.anode_tdc_ch[1] = tmp2;
    raw.anode_tdc_val[0] = vtmp1;
    raw.anode_tdc_val[1] = vtmp2;
  }


  if(raw.anode_tdc_num == 2) {
    if ((raw.anode_tdc_ch[0]%5) == (raw.anode_tdc_ch[1]%5)) {
      if (abs(raw.anode_tdc_ch[0]/PPACMATCH - raw.anode_tdc_ch[1]/PPACMATCH) == 1) {
	if ((raw.anode_tdc_ch[0]/PPACMATCH) == 0 &&
	    (raw.anode_tdc_ch[1]/PPACMATCH) == 1) {
	  tof = (Float_t)(raw.anode_tdc_val[0] - raw.anode_tdc_val[1]);
	  anodeL = raw.anode_tdc_ch[0];
	  anodeR = raw.anode_tdc_ch[1];
	  anodeID = raw.anode_tdc_ch[0]%PPACMATCH;
	} else if ((raw.anode_tdc_ch[0]/PPACMATCH) == 1 &&
		 (raw.anode_tdc_ch[1]/PPACMATCH) == 0) {
	  tof = (Float_t)(raw.anode_tdc_val[1] - raw.anode_tdc_val[0]);
	  anodeL = raw.anode_tdc_ch[1];
	  anodeR = raw.anode_tdc_ch[0];
	  anodeID = raw.anode_tdc_ch[1]%PPACMATCH;
	}
	/*
	printf("Ready for TDC diff, anodeID %d anodeL %d anodeR %d \n",anodeID,anodeL,anodeR);
	*/
	for (Int_t i=0; i<(raw.cathode_tdc_num-1); i++) {
	  for (Int_t j=(i+1); j<raw.cathode_tdc_num; j++) {
	    if ((raw.cathode_tdc_ch[i]/4) == (raw.cathode_tdc_ch[j]/4)) {
	      if (raw.cathode_tdc_ch[i]/4 == anodeL) {
		if ((raw.cathode_tdc_ch[i]%4) == 0 &&
		    (raw.cathode_tdc_ch[j]%4) == 1) {
		  thetaL = raw.cathode_tdc_val[i] - raw.cathode_tdc_val[j];
		  thetaLs = raw.cathode_tdc_val[i] + raw.cathode_tdc_val[j];
		  id = raw.cathode_tdc_ch[i]/4;
		  fThetaL = thetaGain[id]*(Float_t)thetaL;
		  fThetaL += thetaQuad[id]*(TMath::Power((Float_t)thetaL, 2));
		  fThetaL += thetaOffset[id];
		  fThetaL += (Float_t)rand()/d - 0.5;
		  validTL = 1;
		} else if ((raw.cathode_tdc_ch[i]%4) == 1 &&
			   (raw.cathode_tdc_ch[j]%4) == 0) {
		  thetaL = raw.cathode_tdc_val[j] - raw.cathode_tdc_val[i];
		  thetaLs = raw.cathode_tdc_val[j] + raw.cathode_tdc_val[i];
		  id = raw.cathode_tdc_ch[j]/4;
		  fThetaL = thetaGain[id]*(Float_t)thetaL;
		  fThetaL += thetaQuad[id]*(TMath::Power((Float_t)thetaL, 2));
		  fThetaL += thetaOffset[id];
		  fThetaL += (Float_t)rand()/d - 0.5;
		  validTL = 1;
		} else if ((raw.cathode_tdc_ch[i]%4) == 2 &&
			   (raw.cathode_tdc_ch[j]%4) == 3) {
		  phiL = raw.cathode_tdc_val[i] - raw.cathode_tdc_val[j];
		  phiLs = raw.cathode_tdc_val[i] + raw.cathode_tdc_val[j];
		  id = raw.cathode_tdc_ch[i]/4;
		  fPhiL = phiGain[id]*(Float_t)phiL;
		  fPhiL += phiOffset[id];
		  fPhiL += 36.0*(Float_t)anodeL;
		  fPhiL += (Float_t)rand()/d - 0.5;
		  validPL = 1;
		} else if ((raw.cathode_tdc_ch[i]%4) == 3 &&
			   (raw.cathode_tdc_ch[j]%4) == 2) {
		  phiL = raw.cathode_tdc_val[j] - raw.cathode_tdc_val[i];
		  phiLs = raw.cathode_tdc_val[j] + raw.cathode_tdc_val[i];
		  id = raw.cathode_tdc_ch[j]/4;
		  fPhiL = phiGain[id]*(Float_t)phiL;
		  fPhiL += phiOffset[id];
		  fPhiL += 36.0*(Float_t)anodeL;
		  fPhiL += (Float_t)rand()/d - 0.5;
		  validPL = 1;
		}
		if (validTL == 1 && validPL == 1) { validL = 1; }
	      } else if (raw.cathode_tdc_ch[i]/4 == anodeR) {
		if ((raw.cathode_tdc_ch[i]%4) == 0 &&
		    (raw.cathode_tdc_ch[j]%4) == 1) {
		  thetaR = raw.cathode_tdc_val[i] - raw.cathode_tdc_val[j];
		  thetaRs = raw.cathode_tdc_val[i] + raw.cathode_tdc_val[j];
		  id = raw.cathode_tdc_ch[i]/4;
		  fThetaR = thetaGain[id]*(Float_t)thetaR;
		  fThetaR += thetaQuad[id]*(TMath::Power((Float_t)thetaR, 2));
		  fThetaR += thetaOffset[id];
		  fThetaR += (Float_t)rand()/d - 0.5;
		  validTR = 1;
		} else if ((raw.cathode_tdc_ch[i]%4) == 1 &&
			   (raw.cathode_tdc_ch[j]%4) == 0) {
		  thetaR = raw.cathode_tdc_val[j] - raw.cathode_tdc_val[i];
		  thetaRs = raw.cathode_tdc_val[j] + raw.cathode_tdc_val[i];
		  id = raw.cathode_tdc_ch[j]/4;
		  fThetaR = thetaGain[id]*(Float_t)thetaR;
		  /* fThetaR = thetaQuad[id]*(TMath::Power((Float_t)thetaR, 2)); replaced 2018-05-18*/
		  fThetaR += thetaQuad[id]*(TMath::Power((Float_t)thetaR, 2));
		  fThetaR += thetaOffset[id];
		  fThetaR += (Float_t)rand()/d - 0.5;
		  validTR = 1;
		} else if ((raw.cathode_tdc_ch[i]%4) == 2 &&
			   (raw.cathode_tdc_ch[j]%4) == 3) {
		  phiR = raw.cathode_tdc_val[i] - raw.cathode_tdc_val[j];
		  phiRs = raw.cathode_tdc_val[i] + raw.cathode_tdc_val[j];
		  id = raw.cathode_tdc_ch[i]/4;
		  fPhiR = phiGain[id]*(Float_t)phiR;
		  fPhiR += phiOffset[id];
		  fPhiR += 36.0*(Float_t)anodeR;
		  fPhiR += (Float_t)rand()/d - 0.5;
		  validPR = 1;
		} else if ((raw.cathode_tdc_ch[i]%4) == 3 &&
			   (raw.cathode_tdc_ch[j]%4) == 2) {
		  phiR = raw.cathode_tdc_val[j] - raw.cathode_tdc_val[i];
		  phiRs = raw.cathode_tdc_val[j] + raw.cathode_tdc_val[i];
		  id = raw.cathode_tdc_ch[j]/4;
		  fPhiR = phiGain[id]*(Float_t)phiR;
		  fPhiR += phiOffset[id];
		  fPhiR += 36.0*(Float_t)anodeR;
		  fPhiR += (Float_t)rand()/d - 0.5;
		  validPR = 1;
		}
		if (validTR == 1 && validPR == 1) { validR = 1; }
	      }
	    }
	  }
	}
      }
    }
  }

  /* 2018-05-08 CMC KLUDGE to see TS without proper CHICO files*/
  /*    particle.t = (Double_t)raw.LEDts; */
  /*
  printf("\n\nin CHICO, raw.LEDts %llu thetaL %d thetaR %d fThetaL %f fThetaR %f \n", raw.LEDts,thetaL,thetaR,fThetaL,fThetaR);
  printf("fPhiL %f fPhiR %f \n",fPhiL,fPhiR);
  printf("validTL %d validTR %d validPL %d validPR %d \n",validTL,validTR,validPL,validPR);
  printf("validT %d validP %d validL %d validR %d \n",validTL,validTR,validPL,validPR);
  */
  particle.LEDts = 0;

  if (validTL == 1 && validTR == 1) { validT = 1; }
  if (validPL == 1 && validPR == 0) { fPhiR = fPhiL + 180.; }
  if (validPL == 0 && validPR == 1) { fPhiL = fPhiR - 180.; }
  if (validPL == 1 && validPR == 1) { validP = 1; }
  if (validT && validP || 1) {
    if (validL && validR || 1) {
      particle.t = (Double_t)raw.LEDts;
      particle.LEDts = raw.LEDts;
      // particle.rf = ((Double_t)raw.RF * CH2NS)*0.1 + (Double_t)rand()/d - 0.5;
      particle.rf = raw.RF;
      particle.id = anodeID;
      particle.tof = tof;
      particle.thetaL = thetaL;
      particle.thetaR = thetaR;
      particle.fThetaL = (fThetaL)*M_PI/180.;
      particle.fThetaR = (fThetaR)*M_PI/180.;
      particle.phiL = phiL;
      particle.phiR = phiR;

      fPhiL_R = fPhiR + 180.;
      if (fPhiL_R > 360.) {
	fPhiL_R = fPhiL_R - 360.;
      }
      fPhiL_av = (fPhiL + fPhiL_R)/2.;
      fPhiR_L = fPhiL + 180.;
      if (fPhiR_L > 360.) {
	fPhiR_L = fPhiR_L - 360.;
      }
      fPhiR_av = (fPhiR + fPhiR_L)/2.;
      
      particle.fPhiL = (fPhiL)*M_PI/180.;
      particle.fPhiR = (fPhiR)*M_PI/180.;
      
      particle.mass = tof;
      gotParticle = 1;

      /* distances for TOF->Mass
	  anodeL = raw.anode_tdc_ch[1];
	  anodeR = raw.anode_tdc_ch[0];
      */
      particle.dL = 12.8/(0.75471*sinf(particle.fThetaL)*cosf(particle.fPhiL-(18.0+(float)(anodeL)*36.0)*M_PI/180.)
        +0.65606*cosf(particle.fThetaL));
      particle.dR = 12.8/(0.75471*sinf(particle.fThetaR)*cosf(particle.fPhiR-(18.0+(float)(anodeR)*36.0)*M_PI/180.)
        +0.65606*cosf(particle.fThetaR));
      /**/


    } else { gotParticle = 0; }
  } else { gotParticle = 0; }
  /*printf("gotParticle %d fThetaL %f particle.fThetaL %f \n",gotParticle, fThetaL, particle.fThetaL);*/

  /*
  if ( (fabs(particle.fThetaR-0.87)<0.01) && (particle.thetaR>200) ) {
    
    printf("\n\nBad Event\n");
    printf("\n if ( (fabs(particle.fThetaR-0.87) %f  <0.01) && (particle.thetaR>200) %d ) \n",fabs(particle.fThetaR-0.87),particle.thetaR);
  printf("\n\nin CHICO, raw.LEDts %llu thetaL %d thetaR %d fThetaL %f fThetaR %f \n", raw.LEDts,particle.thetaL,particle.thetaR,particle.fThetaL,particle.fThetaR);
  printf("fPhiL %f fPhiR %f \n",fPhiL,fPhiR);
  printf("validTL %d validTR %d validPL %d validPR %d \n",validTL,validTR,validPL,validPR);
  printf("validT %d validP %d validL %d validR %d \n",validTL,validTR,validPL,validPR);

  for (Int_t i=0; i<raw.cathode_tdc_num; i++) { printf("i %d raw.cathode_tdc_val[i] %d raw.cathode_tdc_ch[i] %d \n ", i , raw.cathode_tdc_val[i], raw.cathode_tdc_ch[i]); }
  for (Int_t i=0; i<raw.anode_tdc_num; i++) { printf("i %d raw.anode_tdc_val[i] %d raw.anode_tdc_ch[i] %d \n ", i , raw.anode_tdc_val[i], raw.anode_tdc_ch[i]); }

  } else if ( (fabs(particle.fThetaR-0.87)<0.01) && (particle.thetaR<0) ) {
    
    printf("\n\nGood Event\n");
  printf("\n\nin CHICO, raw.LEDts %llu thetaL %d thetaR %d fThetaL %f fThetaR %f \n", raw.LEDts,particle.thetaL,particle.thetaR,particle.fThetaL,particle.fThetaR);
  printf("fPhiL %f fPhiR %f \n",fPhiL,fPhiR);
  printf("validTL %d validTR %d validPL %d validPR %d \n",validTL,validTR,validPL,validPR);
  printf("validT %d validP %d validL %d validR %d \n",validTL,validTR,validPL,validPR);

  for (Int_t i=0; i<raw.cathode_tdc_num; i++) { printf("i %d raw.cathode_tdc_val[i] %d raw.cathode_tdc_ch[i] %d \n ", i , raw.cathode_tdc_val[i], raw.cathode_tdc_ch[i]); }
  for (Int_t i=0; i<raw.anode_tdc_num; i++) { printf("i %d raw.anode_tdc_val[i] %d raw.anode_tdc_ch[i] %d \n ", i , raw.anode_tdc_val[i], raw.anode_tdc_ch[i]); }

  }
  */
}

Float_t CHICOFull::calcCos(Float_t pTheta, Float_t pPhi, 
			   Float_t gTheta, Float_t gPhi) {
  Float_t calcCos;
  calcCos = sinf(pTheta)*sinf(gTheta)*cosf(pPhi-gPhi) + cosf(pTheta)*cosf(gTheta);
  return calcCos;
}



void CHICOFull::setupCut(){
  TFile *fcutg;
  char WinName[64];

  fcutg = new TFile("chicoCalibrations/ChicoCut", "read");

  sprintf(WinName, "massLCut");
  massLCut = (TCutG*)fcutg->Get(WinName);
  if(massLCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  sprintf(WinName, "massRCut");
  massRCut = (TCutG*)fcutg->Get(WinName);
  if(massRCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  sprintf(WinName, "etCut");
  etCut = (TCutG*)fcutg->Get(WinName);
  if(etCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  sprintf(WinName, "etCutTK");
  etCutTK = (TCutG*)fcutg->Get(WinName);
  if(etCutTK == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  sprintf(WinName, "GT_TKCut");
  GT_TKCut = (TCutG*)fcutg->Get(WinName);
  if(GT_TKCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  sprintf(WinName, "ThetaLRCut");
  ThetaLRCut = (TCutG*)fcutg->Get(WinName);
  if(ThetaLRCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  //fcutg->ls();
  fcutg->Close();
}


Int_t CHICOFull::IDparticles() {

  /*
    Return the left/right mapping of the particles
    Define:    0 : fails gating, no assignment
               1 : Target = Right, Projectile = Left
               2 : Target = Left , Projectile = Right
  */

  Int_t mapping = 0;
  Double_t MassT,MassP,Qval,BetaP,BetaT;

  if(massLCut->IsInside(particle.fThetaL/M_PI*180.,particle.tof)){	
    if(ThetaLRCut->IsInside(particle.fThetaL/M_PI*180.,particle.fThetaR/M_PI*180.)){
      mapping = 1;
    MassT = GetMass(LEFT,TARGET);
    MassP = GetMass(LEFT,BEAM);
    Qval = GetQval(LEFT);
    BetaP = GetBeta(LEFT,BEAM,MassP);
    BetaT = GetBeta(LEFT,TARGET,MassT);
    /*
    MassTheta->Fill(MassT,particle.fThetaR/M_PI*180); //Target=Right
    MassTheta->Fill(MassP,particle.fThetaL/M_PI*180);
    QvalTheta->Fill(Qval,particle.fThetaL/M_PI*180);
    MassTQval->Fill(MassT,Qval);
    MassPQval->Fill(MassP,Qval);
    */
    }
  }
  else if(massRCut->IsInside(particle.fThetaR/M_PI*180.,particle.tof)){	
    if(ThetaLRCut->IsInside(particle.fThetaR/M_PI*180.,particle.fThetaL/M_PI*180.)){
      mapping = 2;
    MassT = GetMass(RIGHT,TARGET);
    MassP = GetMass(RIGHT,BEAM);
    Qval = GetQval(RIGHT);
    BetaP = GetBeta(RIGHT,BEAM,MassP);
    BetaT = GetBeta(RIGHT,TARGET,MassT);
    /*
    MassTheta->Fill(MassT,particle.fThetaL/M_PI*180);
    MassTheta->Fill(MassP,particle.fThetaR/M_PI*180);
    QvalTheta->Fill(Qval,particle.fThetaR/M_PI*180);
    MassTQval->Fill(MassT,Qval);
    MassPQval->Fill(MassP,Qval);
    */
    }
  }

  if (mapping == 0) {
    particle.massP=0;
    particle.massT=0;
    particle.Qvalue=0;
    particle.betaP=0;
    particle.betaT=0;
  } else  {
    /*   printf("mapping %d MassT %f MassP %f Qval %f BetaP %f BetaT %f \n",mapping,MassT,MassP,Qval,BetaP,BetaT); */
    particle.massP=MassP;
    particle.massT=MassT;
    particle.Qvalue=Qval;
    particle.betaP=BetaP;
    particle.betaT=BetaT;
  }

  return mapping;
}


//////////////////////////////////////////////////////////////////////////////
Float_t CHICOFull::GetMass(int LR, int PT){
  float theta1=0.0,theta2=0.0;
  float L1=0.0,L2=0.0;
  float Mass;
  float dT=0.0;
  float f1,f2;
  float p1,p2;
  if(LR==LEFT){
    theta1 = particle.fThetaL;
    theta2 = particle.fThetaR;
    L1 = particle.dL;
    L2 = particle.dR;
    //dT = (particle.dT-2.5)*ch2ns;
    dT = (particle.tof)*CH2NS;
  }
  else if(LR==RIGHT){
    theta1 = particle.fThetaR;
    theta2 = particle.fThetaL;
    L1 = particle.dR;
    L2 = particle.dL;
    //dT = -1.*(particle.dT-2.5)*ch2ns;
    dT = -1.*(particle.tof)*CH2NS;
  }

  f1 = 1.0127 - 0.015567/cosf(theta1) - 0.00027722/powf(cosf(theta1),2.);  //1088 MeV;
  f2 = 1.0068 - 0.0040415/cosf(theta2) - 0.0088646/powf(cosf(theta2),2);   //1088 MeV;
  //f1 = 1.0132 - 0.016286/cosf(theta1) - 0.00027716/powf(cosf(theta1),2.);    //1028.5 MeV;
  //f2 = 0.83469 - 0.0037606/cosf(theta2) - 0.0074352/powf(cosf(theta2),2);    //1028.5 MeV;
  p1 = P0*sinf(theta2)/sinf(theta1+theta2)*f1;
  p2 = P0*sinf(theta1)/sinf(theta1+theta2)*f2;
  //p1 = P0*sinf(theta2)/sinf(theta1+theta2);
  //p2 = P0*sinf(theta1)/sinf(theta1+theta2);
  if(PT==BEAM){
    Mass = (0.032206*dT + Mtot*L2/p2)/(L1/p1+L2/p2);
    return Mass;
  }
  else if(PT==TARGET){
    Mass = (-0.032206*dT + Mtot*L1/p1)/(L1/p1+L2/p2);
    return Mass;
  }
  return 0.0;
}
//////////////////////////////////////////////////////////////////////////////
Float_t CHICOFull::GetBeta(int LR, int PT, float MASS){
  float theta1=0.0,theta2=0.0;
  //float e;
  float f1,f2;
  float beta;
  float p1,p2;

  if(LR==1){
    theta1 = particle.fThetaL;
    theta2 = particle.fThetaR;
  }
  else if(LR==2){
    theta1 = particle.fThetaR;
    theta2 = particle.fThetaL;
  }

  f1 = 1.0127 - 0.015567/cosf(theta1) - 0.00027722/powf(cosf(theta1),2.); //1088 MeV
  f2 = 1.0068 - 0.0040415/cosf(theta2) - 0.0088646/powf(cosf(theta2),2); //1088 MeV
  //f1 = 1.0132 - 0.016286/cosf(theta1) - 0.00027716/powf(cosf(theta1),2.);    //1028.5 MeV;
  //f2 = 0.83469 - 0.0037606/cosf(theta2) - 0.0074352/powf(cosf(theta2),2);    //1028.5 MeV;
  //f1 = 0.045129 + 0.083589*cosf(theta1); //beta vs theta1
  //f2 = -0.003948 + 0.10964*cosf(theta2); //beta vs theta2

  p1 = P0*sinf(theta2)/sinf(theta1+theta2)*f1;
  p2 = P0*sinf(theta1)/sinf(theta1+theta2)*f2;
  if(PT==1){
    beta = p1/MASS*1.0658e-3;
    //e = powf(sinf(theta2)/sinf(theta1+theta2),2)*beamE;
    //beta=sqrt(e/beamA)*0.046;
    //beta*=f1*0.864;
    //beta*=f1;
    return beta;
  }
  else if(PT==2){
    beta = p2/MASS*1.0658e-3;
    //e = powf(sinf(theta1)/sinf(theta1+theta2),2)*beamE*beamA/targetA;
    //beta=sqrt(e/targetA)*0.046;
    //beta*=f2*1.136;              
    //beta*=f2;             
    return beta;
  }
  return 0.0;
}
//////////////////////////////////////////////////////////////////////////////
Float_t CHICOFull::GetQval(int LR){
  float theta1=0.0,theta2=0.0;
  float e1,e2;
  float Qval;

  if(LR==1){
    theta1 = particle.fThetaL;
    theta2 = particle.fThetaR;
  }
  else if(LR==2){
    theta1 = particle.fThetaR;
    theta2 = particle.fThetaL;
  }
  e1 = powf(sinf(theta2)/sinf(theta1+theta2),2)*beamE;
  e2 = powf(sinf(theta1)/sinf(theta1+theta2),2)*beamE*beamA/targetA;
  Qval = e1 + e2 - beamE;
  return Qval;
}
//////////////////////////////////////////////////////////////////////
