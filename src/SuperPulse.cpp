/**********************************************************************/
/* File: SuperPulse.C                                                 */
/* Description: Functions for GRETINA superpulse analysis             */
/* Author: H. Crawford                                                */
/* Date: January 2013                                                 */
/**********************************************************************/

#include "GRETINA.h"

ClassImp(SuperPulse);

void SuperPulse::Initialize(controlVariables* ctrl, GRETINAVariables* gVar) {
  CFD_INT_LEN = 4;
  CFD_DELAY = 4;
  CFD_FRACTION = 4;
  TR_SCALE = 10000.;

  for (Int_t i=0; i<MAXCHANNELS; i++) {
    ehiGeOffset[i] = 0;
    ehiGeGain[i] = 1.0;
  }

  lowE = 0; highE = 0;
  
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    mult[i] = 0;
    crystalBuild[i] = 0;
    ccE[i] = 0.;
    segE[i] = 0.;
    
    netSeg[i] = 0;
    for (Int_t j=0; j<40; j++) {
      segEventIndex[i][j] = 0;
      segsNet[i][j] = 0;
      segs[i][j] = 0;
      if (j<36) { data4net[i][j] = 0; }
      for (Int_t k=0; k<4096; k++) {
	averageTrace[i][j][k] = 0.0;
	averageTraceINT[i][j][k] = 0;
      }
    }
  }

  /* Check if this is a superpulse analysis, if it is, 
     get calibrations, etc. */
  if (ctrl->superPulse) {
    lowE = ctrl->spLowE; highE = ctrl->spHighE;
    ReadDetMaps(ctrl->spCalibDirectory, gVar);
    ReadParams(ctrl->spXtalkFile, "delay1", delay1, 40, gVar);
    for (Int_t i=0; i<MAXCRYSTALS; i++) {
      delay1[i][39] = 0.0f;
    }
  }
  
}

/****************************************************/

Int_t SuperPulse::ReadDetMaps(char *fn, GRETINAVariables* gVar) {

  /***********************************************/
  /* Detmap format:                              */
  /*  BdID(0x3 - 0x6) Chn# 0 SegID Offset Gain   */
  /***********************************************/

  FILE *fp;
  Int_t ccNum, posN, Qnum, fileFlag;
  Int_t bd, ch, i1, seg, nn;
  Float_t off, gn;
  char *st, str[128];
  
  char directory[512];
  char filename[512];

  // for (Qnum = 1; Qnum <= MAXQUADS; Qnum++) {
  for (Qnum = 13; Qnum <= 13; Qnum++) {
    //   for (posN = 1; posN <= 4; posN++) {
    for (posN = 2; posN <= 2; posN++) {
      for (ccNum = 0; ccNum < 4; ccNum++) {
	// for (fileFlag = 0; fileFlag < 2; fileFlag++) {
	for (fileFlag = 0; fileFlag < 1; fileFlag++) {
	  /* Set filename...loop over all CC */
	  strcpy(directory, fn);
	  strcat(directory, "/");
	  if (fileFlag == 0) {
	    strcpy(filename, "detmap_Q");
	  } else if (fileFlag == 1) {
	    strcpy(filename, "tr_gain_Q");
	  }
	  char number[10];
	  sprintf(number, "%dpos", Qnum);
	  strcat(filename, number);
	  sprintf(number, "%d", posN);
	  strcat(filename, number);
	  if (ccNum == 0) { strcat(filename, "_CC9.txt"); }
	  else if (ccNum == 1) { strcat(filename, "_CC19.txt"); }
	  else if (ccNum == 2) { strcat(filename, "_CC29.txt"); }
	  else if (ccNum == 3) { strcat(filename, "_CC39.txt"); }
	  
	  strcat(directory, filename);
	  
	  fp = fopen(directory, "r");
	  if (fp == NULL) {
	    //printf("Could not open \"%s\". --> Skipping! \n", filename);
	  } else {
	    printf("Opened \"%s\".\n", filename);
	    /* Determine the array ID value, based on electronics ordering */
	    Int_t gIDstart = gVar->electronicsOrder[Qnum-1]*160 + (posN-1)*40; 
	    
	    /* Read values */
	    if (fileFlag == 0) {
	      nn = 0;
	      st = fgets(str, 64, fp);
	      while (st != NULL) {
		if (str[0] == 35) {
		  /* '#' comment line, do nothing */
		} else if (str[0] == 59) {
		  /* ';' comment line, do nothing */
		} else if (str[0] == 10) {
		  /* Empty line, do nothing */
		} else {
		  sscanf(str, "%X %d %d %d %f %f", &bd, &ch, &i1, &seg, &off, &gn);
		  Int_t segCh = (bd-3)*10 + ch;
		  map[(Int_t)gIDstart/40][segCh] = seg;
		  if (seg == 36 && ccNum==2) { map[(Int_t)gIDstart/40][segCh] = seg; }
		  else if (seg == 36 && ccNum==0) { map[(Int_t)gIDstart/40][segCh] = 38; }
		  else if (seg == 36) { map[(Int_t)gIDstart/40][segCh] = seg + ccNum; }
		  ehiGeOffset[(Int_t)segCh + gIDstart] = off;
		  ehiGeGain[(Int_t)segCh + gIDstart] = gn;
		  
		  fflush(stdout);
		  nn++;
		}
		
		/* Attempt to read the next line */
		st = fgets(str, 64, fp);
	      }
	    
	      /* Done */
	      fclose(fp);
	    } else if (fileFlag == 1) {
	      
	      nn = 0;
	      st = fgets(str, 64, fp);
	      while (st != NULL) {
		if (str[0] == 35) {
		  /* '#' comment line, do nothing */
		} else if (str[0] == 59) {
		  /* ';' comment line, do nothing */
		} else if (str[0] == 10) {
		  /* Empty line, do nothing */
		} else {
		  sscanf(str, "%f", &gn);
		  if (nn < 36) {
		    trGain[gIDstart + nn] = gn;
		  } else if (nn==36) {
		    trGain[gIDstart + nn + ccNum] = gn;
		  }
		  fflush(stdout);
		  nn++;
		}
		
		/* Attempt to read the next line */
		st = fgets(str, 64, fp);
	      }
	      
	      /* Done */
	      fclose(fp);
	    }
	  }
	}
      }
    }
  }

  printf("Finished detmap / trgain files.\n");

  return (0);
      
}

/****************************************************/

Int_t SuperPulse::ReadParams(TString filename, const char *label, 
			     Float_t x[][40], Int_t len, GRETINAVariables* gVar) {
  FILE *listFile;
  listFile = fopen(filename.Data(), "r");
  if (listFile == NULL) {
    printf("Could not open \"%s\".\n", filename.Data());
    exit(1); 
  }
  
  char *st, str[512];
  Int_t Qnum, posNum;
  Int_t lines = 0;
  
  char temp[512];
  char *fNames[MAXCRYSTALS];
  
  lines = 0;
  st = fgets(str, 512, listFile);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%d %d %s", &Qnum, &posNum, temp);
      fNames[gVar->electronicsOrder[Qnum-1] + posNum - 1] = temp; 
      lines++;
    }

    /* Attempt to read next line */
    st = fgets(str, 512, listFile);
  }
  fclose(listFile);

  FILE *fin;
  Int_t max_length = 120;
  char s[max_length], *tok0, *tok, *u, *s0;
  Int_t nn = 0;
  Float_t a;

  for (Qnum = 0; Qnum <= MAXQUADS; Qnum++) {
    for (posNum = 1; posNum <= 4; posNum++) {
      Int_t xNum = -1;
      xNum = gVar->electronicsOrder[Qnum-1] + posNum - 1; 
      fin = fopen(fNames[xNum], "r");
      if (fin == NULL) {
	//printf("Could not open \"%s\". --> Skipping!\n", fNames[xNum]);
      } else {
	//printf("Opened \"%s\".\n", fNames[xNum]);
	while (fgets(s, max_length, fin) != 0) {
	  tok0 = strtok(s, " ");
	  if ((u = strchr(tok0, atoi(":"))) != 0) {
	    *u = '\0';
	    if (strcmp(label, tok0) == 0) {
	      while(1) {
		s0 = fgets(s, max_length, fin);
		if (s0 == 0 || s[0] == '\n') {
		  fclose(fin);
		  return nn;
		}
		s0 = s;
		while ((tok = strtok(s0, " \t\n")) != 0) {
		  s0 = 0;
		  a = (Float_t) atof(tok);
		  x[xNum][nn++] = a;
		  if (nn >= len) {
		    fclose(fin);
		    return nn;
		  }
		}
	      }
	    }
	  }
	  fclose(fin);
	}
      }
    }
  }
  printf("Finished param files.\n");
  return 0;
}

/****************************************************/
 
void SuperPulse::MakeSuperPulses() {
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    if (mult[i] == 1 && crystalBuild[i] == 820) { /* Mult 1, 40 segments */
      if (segE[i] >= lowE && ccE[i] >= lowE && segE[i] <= highE && ccE[i] <= highE) {
	Int_t t0 = AlignCFD(i);
	if (t0 >= 0) { /* CFD alignment was successful. */
	  data4net[i][netSeg[i]]++;
	  for (Int_t m=0; m<37; m++) { /* 36 segments + 1 CC only */
	    for (Int_t j=0; j<AVG_TR_LENGTH; j++) {
	      /* We fill one giant array with all 37 traces,
		 with small gaps between waveforms. */
	      averageTrace[i][netSeg[i]][m*(AVG_TR_STRIDE) + j] += waves[i][m][j];
	    } /* Loop over waveform samples */
	  } /* Loop over segments */
	} /* if CFD alignment is OK */
      } /* in energy window */
    } /* if mult = 1 and 40 segments */
  } /* loop over MAXCRYSTALS crystals */
}

/****************************************************/

Int_t SuperPulse::AlignCFD(Int_t crystalNum) {
  Int_t i, j, k;
  Float_t t=0, tmin = 1000.0;

  /* Find the time to align to. */
  for (j=0; j<40; j++) {
    if (netSeg[crystalNum] == j || (j ==38)) { /* Net segment or CC */
      t = cfdTime(crystalNum, j);
      if (t < 0.0) { return -1; }
      t -= delay1[crystalNum][j]; /* DCR - there may be a question about - sign */
      if (tmin > t) { tmin = t; }
    }
  }
  i = (Int_t)(t + 0.5) - 16; /* Starting step of the trace to extract for SP */
  if (i == 0) { return 0; }
  
  /* Shift the traces.  Order depends on if this is a left or right shift. */
  if (i < 0) {
    for (j=0; j<40; j++) {
      for (k=trLength - 2; k > (0-i); k--) {
	waves[crystalNum][j][k] = waves[crystalNum][j][k+i];
      }
      for (k=0; k<(0-i); k++) {
	waves[crystalNum][j][k] = 0;
      }
    }
  } else {
    for (j=0; j<40; j++) {
      for (k=0; k<=trLength - i - 2; k++) {
	waves[crystalNum][j][k] = waves[crystalNum][j][k+i];
      }
      for (k=trLength - i - 1; k<trLength; k++) {
	waves[crystalNum][j][k] = waves[crystalNum][j][trLength-2];
      }
    }
  }

  return 0;
}

/****************************************************/

Float_t SuperPulse::cfdTime(Int_t crystalNum, Int_t segNum) {
  /* Ge-style CFD to get time from the signal trace.  
     From DCR. */
  Int_t deriv[2048];
  Int_t cfd[2048];
  Int_t i, imax = 0, max_deriv = 0;
  
  deriv[0] = 0;
  for (i=0; i< CFD_INT_LEN; i++) {
    deriv[0] += (waves[crystalNum][segNum][i+CFD_INT_LEN] -
		 waves[crystalNum][segNum][i]);
  }
  for (i=1; i<trLength - 5 - 2*CFD_INT_LEN; i++) {
    deriv[i] = (deriv[i-1] + 
		waves[crystalNum][segNum][i+2*CFD_INT_LEN] -
		2*waves[crystalNum][segNum][i+CFD_INT_LEN] + 
		waves[crystalNum][segNum][i]);
    if (max_deriv < deriv[i]) {
      max_deriv = deriv[i];
      imax = i; 
    }
  }
  for (i=0; i<trLength - 5 - 2*CFD_INT_LEN - CFD_DELAY; i++) {
    cfd[i] = deriv[i] - deriv[i+CFD_DELAY]/CFD_FRACTION;
  }

  for (i=imax + CFD_DELAY; i>0; i--) {
    if (cfd[i] <= 0 && cfd[i+1] > 0) {
      /* Interpolate zero crossing and return time in steps. */
      return ((Float_t)i) - ((Float_t)cfd[i]) / ((Float_t)(cfd[i+1] - cfd[i]));
      break;
    }
  }
  return -1.0;
}

/****************************************************/

void SuperPulse::FinishSuperPulses() {
  /* Scale the superpulses, and get the trace gains. */

  for (Int_t i=0; i<MAXCRYSTALS; i++) {

    for (Int_t j=0; j<36; j++) {
      gain[i][j] = (((Float_t)averageTrace[i][j][j*AVG_TR_STRIDE + AVG_TR_LENGTH - 1]) /
		    ((Float_t)averageTrace[i][j][36*AVG_TR_STRIDE + AVG_TR_LENGTH - 1]));
    }

    gain[i][36] = 1.0f;

    for (Int_t j=0; j<36; j++) {
      Float_t scaleFactor = TR_SCALE / ((Float_t)averageTrace[i][j][36*AVG_TR_STRIDE + AVG_TR_LENGTH -1]);
      for (Int_t k=0; k<37; k++) {
	if (k<36) {
	  for (Int_t m=0; m<AVG_TR_LENGTH; m++) {
	    averageTrace[i][j][k*AVG_TR_STRIDE + m] = ((Float_t)averageTrace[i][j][k*AVG_TR_STRIDE + m])*scaleFactor/gain[i][k];
	    averageTraceINT[i][j][k*AVG_TR_STRIDE + m] = (Int_t)averageTrace[i][j][k*AVG_TR_STRIDE + m];
	  }
	} else {
	  for (Int_t m=0; m<AVG_TR_LENGTH; m++) {
	    averageTrace[i][j][k*AVG_TR_STRIDE + m] = ((Float_t)averageTrace[i][j][k*AVG_TR_STRIDE + m])*scaleFactor/gain[i][k];
	    averageTraceINT[i][j][k*AVG_TR_STRIDE + m] = (Int_t)averageTrace[i][j][k*AVG_TR_STRIDE + m];
	  }

	}
      }
    }
  }
}

/****************************************************/

void SuperPulse::WriteSuperPulses() {
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    if (data4net[i][0] > 0 || data4net[i][1] > 0) {
      char filenameOut[1024];
      sprintf(filenameOut, "SPCrystal%d_%d.spn", i, (Int_t)((lowE + highE)/2));
      cout << "Writing superpulses to " << filenameOut << endl;
      
      FILE *spOut = fopen(filenameOut, "wb");
      for (Int_t j=0; j<36; j++) {
	fwrite(averageTraceINT[i][j], sizeof(Int_t), 4096, spOut);
      }
      fclose(spOut);
      printf("Superpulse statistics:\n");
      for (Int_t k=0; k<36; k++) {
	printf("%d -- %d\n", k, data4net[i][k]);
      }
    }
  }
}

    
							    

