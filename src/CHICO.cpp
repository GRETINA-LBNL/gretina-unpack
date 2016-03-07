/*! \file CHICO.cpp
    \brief Parameter and functions for CHICO analysis.

    This file provides the analysis functions for CHICO analysis, based on 
    an analysis code CheckChico.cc from S. Zhu (ANL).

    Author: H.L. Crawford
    Date: April 29, 2014
*/

#include "CHICO.h"

void CHICORaw::Initialize() {

}

void CHICORaw::Reset() {
  Initialize();
}

void CHICOParticle::Initialize() {

}

void CHICOParticle::Reset() {
  Initialize();
}

CHICOFull::CHICOFull() {
  gotParticle = 0;
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
  Int_t nn = 0, i1;
  Float_t f1, f2, f3;

  nn = 20;

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
      betaP[nn] = f1; betaT[nn] = f2;
      nn++;
    }
    retval = fgets(line, 500, in);
  }
  printf("Done.\n");
  if (nn != 160) {
    perror("Unexpected number of beta calibration parameters read.");
    return;
  } else {
    printf("Beta calibration loaded.\n");
  }
    
  fclose(in);
}

void CHICOFull::InitializeCHICOVariables(TString thetaCalFile,
					 TString phiCalFile,
					 TString betaCalFile) {
  ReadThetaCalibration(thetaCalFile);
  ReadPhiCalibration(phiCalFile);
  ReadBetaCalibration(betaCalFile);
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

  GetParticle();

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
		  fThetaR = thetaQuad[id]*(TMath::Power((Float_t)thetaR, 2));
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

  if (validTL == 1 && validTR == 1) { validT = 1; }
  if (validPL == 1 && validPR == 0) { fPhiR = fPhiL + 180.; }
  if (validPL == 0 && validPR == 1) { fPhiL = fPhiR - 180.; }
  if (validPL == 1 && validPR == 1) { validP = 1; }
  if (validT && validP) {
    if (validL && validR) {
      particle.t = (Double_t)raw.LEDts;
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
    } else { gotParticle = 0; }
  } else { gotParticle = 0; }

}

Float_t CHICOFull::calcCos(Float_t pTheta, Float_t pPhi, 
			   Float_t gTheta, Float_t gPhi) {
  Float_t calcCos;
  calcCos = sinf(pTheta)*sinf(gTheta)*cosf(pPhi-gPhi) + cosf(pTheta)*cosf(gTheta);
  return calcCos;
}
