#include "Track.h"

#include "Globals.h"

ClassImp(TrackClusterIP);
ClassImp(TrackCluster);
ClassImp(TrackShell);
ClassImp(TrackCtrl);
ClassImp(TrackStat);
ClassImp(TrackPerm);
ClassImp(Track);

#define DEBUG 0

/**********************************************************/

Double_t TrackClusterIP::xx() { return xyz.X(); }

Double_t TrackClusterIP::yy() { return xyz.Y(); }

Double_t TrackClusterIP::zz() { return xyz.Z(); }

/**********************************************************/

Int_t TrackCluster::findFOM(TVector3 targetPos) {

  Int_t j, k;
  Double_t thetaTmp;  Float_t thetaComp;
  Float_t pFac;
  TVector3 vA, vB, v1, v2;
  Float_t fomTmp, eg;

  /* Find the first vector, from target to the first interaction point */
  //for (k=0; k<nDet; k++) {
  //  if (intpts[k].order == 0) {
  //    vA = (intpts[k].xyz - targetPos).Unit();
  //  }
  //}
  vA = (intpts[intpts[0].order].xyz - targetPos).Unit();
  eg = eSum;
  fomTmp = 0;

  Int_t b1, b2;

  k = 0;
  while (k < nDet - 1) {
    for (j=0; j<nDet; j++) {
      // if (intpts[j].order == k) { v1 = intpts[j].xyz;  b1 = j; }
      // if (intpts[j].order == k+1) { v2 = intpts[j].xyz;  b2 = j; }
    }
    
    v1 = intpts[intpts[k].order].xyz; b1 = intpts[k].order;
    v2 = intpts[intpts[k+1].order].xyz; b2 = intpts[k+1].order;

    vB = (v2 - v1).Unit();
   
    if (DEBUG) { 
      printf("b1, b2: %d %d\n", b1, b2);
      printf("%d: A(%0.3f, %0.3f, %0.3f)\n", k, v1.X(), v1.Y(), v1.Z());
      printf("%d: B(%0.3f, %0.3f, %0.3f)\n", k, v2.X(), v2.Y(), v2.Z());
    }    

    /* Find the scattering angle */
    thetaTmp = vA.Angle(vB);
    
    /* Sanity checks -- invalidate so nobody else will see it */
    if (thetaTmp < 0) { 
      tracked = 0;  valid = 0;  return (5);
    } else if (thetaTmp > TMath::Pi()) {
      tracked = 0;  valid = 0;  return (6);
    } 

    /* Find compton scattering angle */
    pFac =  findComptonAngle(eg, intpts[b1].eDet, &thetaComp);

    /* Add to Figure of Merit */
    fomTmp += (thetaComp - thetaTmp)*(thetaComp - thetaTmp);

    /* Add penalty if the energies make no sense
       in the Compton formula */
    //fomTmp += (pFac*badTheoAngPenalty);

    /* Shift first vector */
    vA = vB;

    /* Subtract energy lost in the scatter */
    eg = eg - intpts[b1].eDet;
  
    k++;
  } /* Loop over interaction points in the sequence */
  
  /* Done! */
  FOM = fomTmp;
 
  return (0);
}

Float_t TrackCluster::findComptonAngle(Float_t eg, Float_t ee, Float_t *thc) {
  /* eg: incoming gamma energy, ee: interaction energy */
  Float_t egp, cosa;
  
  /* Calculate Compton angle */
  egp = eg - ee;
  cosa = 1. + 0.511 / eg - 0.511 / egp;  
  if (cosa > 1) {
    *thc = acosf(1.0);
    return (cosa - 1.0);
  } else if (cosa < -1) {
    *thc = acosf(-1.0);
    return (-1.0 - cosa);
  } else {
    *thc = acosf(cosa);
    return (0);
  }
  
  /* Done */
  return (0);
}

Int_t TrackCluster::sortSegEnergies() {
  Int_t i, j;  Int_t iTmp; Float_t fTmp; TVector3 vTmp;
  long long int tsTmp;

  /* GEANT simulation show that (a) one average
     we lose half the energy for every interaction and that
     photo-absorption mostly happen at about 200keV.
     As a matter of fact, if the last interaction energy
     is not near 200keV the sequence is probably wrong. */
  
  /* Bubble sort segment energies */
  for (i=0; i<nDet; i++) {
    for (j=i+1; j<nDet; j++) {
      if (intpts[i].eDet < intpts[j].eDet) {
	vTmp = intpts[i].xyz;
	intpts[i].xyz = intpts[j].xyz;	
	intpts[j].xyz = vTmp;

	fTmp = intpts[i].eDet;
	intpts[i].eDet = intpts[j].eDet;  
	intpts[j].eDet = fTmp;
	
	iTmp = intpts[i].order;
	intpts[i].order = intpts[j].order;  
	intpts[j].order = iTmp;

	tsTmp = intpts[i].timestamp;
	intpts[i].timestamp = intpts[j].timestamp;  
	intpts[j].timestamp = tsTmp;

	fTmp = intpts[i].t0;
	intpts[i].t0 = intpts[j].t0;  
	intpts[j].t0 = fTmp;

	iTmp = intpts[i].shellPos;
	intpts[i].shellPos = intpts[j].shellPos;  
	intpts[j].shellPos = iTmp;

	iTmp = intpts[i].detNum;
	intpts[i].detNum = intpts[j].detNum;  
	intpts[j].detNum = iTmp;

	/* Simulated values */
	fTmp = intpts[i].eSumOrig;
	intpts[i].eSumOrig = intpts[j].eSumOrig;  
	intpts[j].eSumOrig = fTmp;
	
	iTmp = intpts[i].origPos;
	intpts[i].origPos = intpts[j].origPos;	
	intpts[j].origPos = iTmp;
      }
    }
  }

  /* Done! */
  return(0);
}

void TrackCluster::Reset() {
  nDet = 0;  eSum = 0;  tracked = 0;  valid = 1;  FOM = MAXFOM;
  for (Int_t k=0; k<MAXNUMDET; k++) { intpts[k].order = -1; }
}

void TrackCluster::Clear() {
  valid = 0; nDet = 0; tracked = 0; FOM = 0.0; eSum = 0.0;
  trackNum = 0; bestPermutation = 0; processed = 0;
  for (Int_t k=0; k<MAXNUMDET; k++) {
    intpts[k].xyz.SetXYZ(0., 0., 0.); intpts[k].eDet = 0.;
    intpts[k].order = 0; intpts[k].timestamp = 0;
    intpts[k].shellPos = 0; intpts[k].detNum = 0;
    intpts[k].eSumOrig = 0; intpts[k].origPos = 0;
  }
}

void TrackCluster::PrintCluster() {
  printf("TrackCluster::PrintCluster() -->\n");
  if (valid) {
    printf("    valid, tracked?: %i, %i\n", valid, tracked); 
    printf("    eSum, FOM: %0.3f, %0.3f\n", eSum, FOM);
    printf("    nDet: %i\n", nDet);
    printf("    perm: %i\n", bestPermutation);
    for (Int_t j=0; j<nDet; j++) {
      if (intpts[j].order == 0) {
	printf("----> IntPt %i:\n", j);
      } else {
	printf("      IntPt %i:\n", j);
      }
      printf("       (x,y,z): (%0.3f, %0.3f, %0.3f)\n", intpts[j].xyz.X(), intpts[j].xyz.Y(), intpts[j].xyz.Z());
      printf("       eDet, order: %0.3f, %i\n", intpts[j].eDet, intpts[j].order);
    }
  }
}

/**********************************************************/

Double_t TrackShell::xx(Int_t i) { return xyz[i].X(); }

Double_t TrackShell::yy(Int_t i) { return xyz[i].Y(); }

Double_t TrackShell::zz(Int_t i) { return xyz[i].Z(); }

Double_t TrackShell::nxx(Int_t i) { return xyz[i].Unit().X(); }

Double_t TrackShell::nyy(Int_t i) { return xyz[i].Unit().Y(); }

Double_t TrackShell::nzz(Int_t i) { return xyz[i].Unit().Z(); }

Double_t TrackShell::relAngle(Int_t i, Int_t j) { return xyz[i].Angle(xyz[j]); }

Double_t TrackShell::polAngle(Int_t i) { return xyz[i].Unit().Theta(); }

Double_t TrackShell::aziAngle(Int_t i) { return xyz[i].Unit().Phi(); }

void TrackShell::Print() {
  printf("TrackShell::Print() -->\n");
  for (Int_t i=0; i<nHit; i++) {
    printf("  IntPt %i: \n", i);
    cout << "    t0, chiSq, normChiSq: " << t0[i] << ", " 
	 << chiSq[i] << ", " << normChiSq[i] << endl;
    cout << "    crystalID, timestamp: " << crystalID[i] << ", " 
	 << timestamp[i] << endl;
    cout << "    (x,y,z), eDet, eSum : (" << xyz[i].X() << ", " 
	 << xyz[i].Y() << ", " << xyz[i].Z() << "), " << eDet[i] 
	 << ", " << eSum[i] << endl;
  }
}

void TrackShell::PrintClusters() {
  printf("TrackShell::PrintClusters() -->\n");
  printf("  NumClusters: %i\n", numClusters);
  for (Int_t i=0; i<nHit; i++) {
    printf("    IntPt %i --> cluster %i\n", i, clusterNum[i]);
  }
}

/**********************************************************/

Int_t TrackCtrl::ReadChatFile(TString fileName) {
  FILE *fp;
  char *p, *pc, str[512] = {'0'}, str1[512] = {'0'}, str2[512] = {'0'};
  Int_t echo = 0, nret = 0, nn = 0, nni = 0, i, i1, i2, st;
  Double_t d1;
  Int_t ng, nt, ok, ndet, j;
  Float_t r1, r2, r3, aa, bb, distToCrystal;
  
  if ((fp = fopen(fileName.Data(), "r")) == NULL) {
    cout << ALERTTEXT;
    printf("TRACK: Error: Could not open chat file <%s>\n", fileName.Data());
    cout << RESET_COLOR; fflush(stdout);
    exit(0);
  }
  printf("TRACK: Opened chat file <%s>\n\n", fileName.Data());
  fflush(stdout);
  
  /* Read content and act */
  pc = fgets(str, 512, fp);
  while (pc != NULL) {
    if (echo) { printf("TRACK: chat.%i: %s", nn, str);  fflush(stdout); }
    
    /* Attempt to interpret the line */
    if ((p = strstr(str, "echo")) != NULL) {
      echo = 1;
      if (echo) printf("TRACK: will echo command lines\n");
      fflush(stdout);
    } else if (str[0] == 35) {
      /* '#' comment line, do nothing */
      nni--; 
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
      nni--;
    } else if (str[0] == 10) {
      /* empty line, do nothing */
      nni--;
    } else if ((p = strstr(str, "trackingstrategy")) != NULL) {
      nret = sscanf(str, "%s %i %i %s", str1, &ndet, &i2, str2);
      if (nret == 3) {
	/* Set the option */
	printf("TRACK: ndet == %i, simple\n", ndet);
	trackOps[ndet] = i2;
	i1 = trackOpt(trackOps[ndet]);
	if (i1 != 0) {
	  printf("TRACK: Bad tracking option, quit...\n");
	  exit(1);
	}
      } else if (nret == 4) {
	printf("TRACK: ndet == %i, complex jumping option\n", ndet);
	if (i2 != 5) { exit(-1); }
	/* Set the option */
	trackOps[ndet] = i2;
	i1 = trackOpt(trackOps[ndet]);
	if (i1 != 0) {
	  printf("TRACK: Bad tracking option, quit...\n");
	  exit(1);
	}
	/* Count g's and t's */
	ng = 0;  nt = 0;  ok = 1;  i = 0;
	while (ok && i < 10) {
	  switch(str2[i]) {
	  case 103:
	    ng++; /* It is a g */
	    break;
	  case 116:
	    nt++; /* It is a t */
	    break;
	  default:
	    ok = 0;
	    break;
	  }
	  i++;
	}
	printf("TRACK:  --> ng = %i, nt = %i\n", ng, nt);
	
	/* Check... */
	if (ndet != (ng + nt)) {
	  printf("TRACK: Bad jump strategy length, quit...\n");
	  exit(1);
	}
	
	/* Set group length */
	jumpGroupLength[ndet] = ng;
	
      } else {
	printf("TRACK: Bad trackingstrategy input\n\n");
	exit(1);
      }
    } else if ((p = strstr(str, "fixdoublehits")) != NULL) {
      nret = sscanf(str, "%s %i", str1, &fixDoubleHitsCode);
      CheckNumArgs(nret, 2, str);
      fixDoubleHits = 1;
    } else if ((p = strstr(str, "fomjump")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &FOMJump);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr(str, "nocrystaltoworldrot")) != NULL) {
      nret = sscanf(str, "%s", str1);
      CheckNumArgs(nret, 1, str);
      printf("TRACK: \"nocrystaltoworldrot\" option disabled.  Ignoring.\n");
    } else if ((p = strstr(str, "badTheoAngPenalty")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &badTheoAngPenalty);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr(str, "ndetElim")) != NULL) {
      nret = sscanf(str, "%s %i %f %f %f", str1, &i1, &r1, &r2, &r3);
      CheckNumArgs(nret, 5, str);
      nDetELim_lo[i1] = r1;
      nDetELim_hi[i1] = r2;
      nDetELim_FOM[i1] = r3;
      printf("TRACK: Will set FOM to %5.2f for ndet == %i if outside %9.3f MeV to %9.3f MeV\n",
	     nDetELim_FOM[i1], i1, nDetELim_lo[i1], nDetELim_hi[i1]);
    } else if ((p = strstr(str, "useCCEnergy")) != NULL) {
      CheckNumArgs(nret, 1, str);
      useCCEnergy = 1;
      printf("TRACK: useCCEnergy = %i\n", useCCEnergy);
    } else if ((p = strstr (str, "useSegEnergy")) != NULL) {
      CheckNumArgs(nret, 1, str);
      useSegEnergy = 1;
      printf ("TRACK: useSegEnergy=%i\n", useSegEnergy);
    } else if ((p = strstr (str, "enabled")) != NULL) {
      nret = sscanf(str, "%s %s", str1, str2);
      CheckNumArgs(nret, 2, str);
      strDecomp(str2, MAXDETNUM + 1, enabled);
    } else if ((p = strstr (str, "singlesfom")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &singlesFOM);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr (str, "fomgoodenough")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &FOMGoodEnough);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr (str, "target_x")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &targetX);
      CheckNumArgs(nret, 2, str);
      targetX /= 10;
    } else if ((p = strstr (str, "target_y")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &targetY);
      CheckNumArgs(nret, 2, str);
      targetY /= 10;
    } else if ((p = strstr (str, "target_z")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &targetZ);
      CheckNumArgs(nret, 2, str);
      targetZ /= 10;
    } else if ((p = strstr (str, "singlehitmaxdepth")) != NULL) {
      nret = sscanf(str, "%s %i %f %f", str1, &i1, 
		    &singleHitRangeFOM, &distToCrystal);
      CheckNumArgs(nret, 4, str);
      
      markSingleHitRange = 1;
      printf ("TRACK: marksinglehitrange = %i\n", markSingleHitRange);
      printf ("TRACK: singlehits outside range will be assigned FOM of %f\n", singleHitRangeFOM);
      printf ("TRACK: distToCrystal= %f (will be added)\n", distToCrystal);
      
      for (i = 0; i < i1; i++) {
	int match = fscanf (fp, "%f %f\n", &r1, &r2);
	i2 = (int) (r1 * 1000 + 0.5);     /* now keV */
	printf ("TRACK:  --> Max range at %6i keV ", i2);  fflush(stdout);
	singleHitRange[i2] = r2 + distToCrystal;
	printf ("is %4.2f cm\n", singleHitRange[i2]); fflush (stdout);
      }
      
      /* Interpolate for the rest */
      for (i = 0; i < LONGLEN; i++)
	if (singleHitRange[i] == 0) {
	  /* Find next up that is filled */	
	  j = i + 1;
	  while (singleHitRange[j] == 0 && j < LONGLEN) { j++; }
	  /* Interpolate */
	  bb = (singleHitRange[i - 1] - singleHitRange[j]) / (i - 1 - j);
	  aa = singleHitRange[j] - bb * j;
	  singleHitRange[i] = aa + bb * i;	
	}
      
      /* A hack to introduce i-yang's singles add-back curve */
      if (0) {
	for (j=0; j<LONGLEN; j++) {
	  singleHitRange[j] = (((-0.0006*j*j + 0.2359*j - 9.4143)/10) + 
			       distToCrystal);
	}
      }
      
      /* Quietly make the value the square range 
	 so we do not have to take square roots later */
      
      for (i = 0; i < LONGLEN; i++) { singleHitRange[i] *= singleHitRange[i]; }
      
    } else if ((p = strstr (str, "splitclusters1")) != NULL) {
      nret = sscanf(str, "%s %f %f %i %i %i", str2, 
		    &splitClusters_kickoutFOM[0], &splitClusters_thresh[0],
		    &splitClusters_minN[0], &splitClusters_maxN[0], 
		    &splitClusters_maxTry[0]);
      CheckNumArgs(nret, 6, str);
      splitClusters[0] = 1;
    } else if ((p = strstr (str, "splitclusters2")) != NULL) {
      nret = sscanf(str, "%s %f %f %i %i %i %f", str2, 
		    &splitClusters_kickoutFOM[1], &splitClusters_thresh[1], 
		    &splitClusters_minN[1], &splitClusters_maxN[1], 
		    &splitClusters_maxTry[1], &goodImproveFraction);
      CheckNumArgs(nret, 7, str);
      splitClusters[1] = 1;
    } else if ((p = strstr (str, "recluster1")) != NULL) {
      nret = sscanf(str, "%s %f %f %i %i %f", str2, &reCluster_kickoutFOM[0], 
		    &reCluster_thresh[0], &reCluster_minN[0], &reCluster_maxTry[0], 
		    &reCluster_reduxFactor[0]);
      CheckNumArgs(nret, 6, str);
      reCluster[0] = 1;
    } else if ((p = strstr (str, "recluster2")) != NULL) {
      nret = sscanf (str, "%s %f %f %i %i %f", str2, &reCluster_kickoutFOM[1], 
		     &reCluster_thresh[1], &reCluster_minN[1], &reCluster_maxTry[1], 
		     &reCluster_reduxFactor[1]);
      CheckNumArgs(nret, 6, str);
      reCluster[1] = 1;
    } else if ((p = strstr (str, "combineclusters")) != NULL) {
      nret = sscanf(str, "%s %f %f %i %f", str2, &combineClusters_kickoutFOM, 
		    &combineClusters_thresh, &combineClusters_maxN, 
		    &combineMaxDist);
      CheckNumArgs(nret, 5, str);
      combineClusters = 1;
      
      /* Quietly make this the square so we don't 
	 have to do squareroots later when we use it. */
      combineMaxDist = combineMaxDist * combineMaxDist; 
    } else if ((p = strstr (str, "matchmaker")) != NULL) {
      nret = sscanf(str, "%s %f %f", str2, &matchMaker_kickoutFOM, 
		    &matchMakerMaxDist);
      CheckNumArgs(nret, 3, str);
      matchMaker = 1;
      
      /* Quietly make this the square so we don't 
	 have to do squareroots later when we use it. */
      matchMakerMaxDist = matchMakerMaxDist * matchMakerMaxDist;
    } else if ((p = strstr (str, "iterations")) != NULL) {
      nret = sscanf(str, "%s %i", str1, &iterations);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr (str, "removeinvalidentries")) != NULL) {
      CheckNumArgs(nret, 1, str);
      printf("TRACK: ReadChatFile(): \"removeinvalidentries\" option removed, no longer needed\n");
    } else if ((p = strstr (str, "maxevents")) != NULL) {
      nret = sscanf(str, "%s %i", str1, &maxEvents);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr (str, "untracked_fom_kickout")) != NULL) {
      nret = sscanf(str, "%s %f", str1, &untracked_kickoutFOM);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr (str, "trackDataStream")) != NULL) {
      nret = sscanf (str, "%s %s", str1, str2);
      CheckNumArgs(nret, 2, str);
      trackDataFileName = str2;
    } else if ((p = strstr (str, "nprint")) != NULL) {
      nret = sscanf(str, "%s %i", str1, &nPrint);
      CheckNumArgs(nret, 2, str);
    } else if ((p = strstr (str, "clusterangle")) != NULL) {
      nret = sscanf(str, "%s %i %f", str1, &i1, &r1);
      CheckNumArgs(nret, 3, str);
      if (i1 >= MAXSHELLHITS) {
	printf("TRACK: Error: max number of interaction points are %i\n", MAXSHELLHITS);
	exit(1);
      }
      alpha[i1] = r1;
      alpha[i1] *= TMath::DegToRad();
    } else if ((p = strstr (str, "dtwin")) != NULL) {
      /* For stand-alone tracking only */
    } else {
      /***** Chat script error point *****/
      cout << ALERTTEXT;
      printf("TRACK: line %2.2i in chat script, option: %s \n  -- not understood\n", nn, str);
      printf("TRACK: %i\n", str[0]);
      printf("TRACK: Aborting...\n");
      cout << RESET_COLOR;  fflush(stdout);
      exit(0);
    }
    
    /* Read next line in chat file */
    nn++;
    nni++;
    pc = fgets(str, 512, fp);
  }

  /* Done */
  fclose(fp);
  printf("TRACK: Chat file <%s> closed\n", fileName.Data());
  printf("TRACK:   --> processed %i track instructions and %i lines\n", nni, nn);
  printf("\n");
  fflush(stdout);

  /* Check the alpha values */
  if (alpha[1] <= 0) {
    cout << ALERTTEXT;
    printf("TRACK: You must at least specify the clustering angle for %i interaction points, cluster angle %i xxx.xxx\n", 1, 1);
    cout << RESET_COLOR;  fflush(stdout);
    exit(1); 
  }
  if (alpha[MAXSHELLHITS-1] <= 0) {
    /* Look for the last non-zero one */
    for (i=(MAXSHELLHITS-1); i>=1; i--) {
      if (alpha[i] > 0) {
	alpha[MAXSHELLHITS-1] = alpha[i];
	break;
      }
    }
  }
  /* Process the alpha clustering angle (interpolate) */
  for (i=1; i<MAXSHELLHITS; i++) {
    if (alpha[i] == 0) {
      /* Find next that is filled */
      j = i+1;
      while (alpha[j] == 0 && j < MAXSHELLHITS) { j++; }
      /* Interpolate */
      bb = (alpha[i-1] - alpha[j])/(i - 1 - j);
      aa = alpha[j] - bb*j;
      alpha[i] = aa + bb*i;
    }
  }

  printf("TRACK: Chat file initialization complete. \n");

  return (0);
}

void TrackCtrl::CheckNumArgs(Int_t required, Int_t actual, char *str) {
  if (required < actual) {
    printf("TRACK: Argument problem with chat option\n");
    printf("TRACK: --> %s\n", str);
    printf("TRACK: Required # arguments: %i\n", required-1);
    printf("TRACK: Provided # arguments: %i\n", actual-1);
    printf("TRACK: Please fix and try again.  Quitting now...\n");
    exit(1);
  }
}

/* Decompose a detector list like: 1-5,7,9-12 into an array */

Int_t TrackCtrl::strDecomp(char *str, Int_t dim, Int_t *yy) {
  Int_t i, j, pos, nn, ok, lo, hi;
  char lstr[100];
  // printf("TRACK: strDecomp: Decomposing string <%s> into array [0...%i]\n",
  //	 str, dim-1);
  
  /* Zero the array */
  for (i=0; i<dim; i++) { yy[i] = 0; }
  
  pos = 0; ok = 1;
  while (ok == 1) {
    /* Search for the next subrange */
    nn = 0;
    while ((str[pos+nn] != ',') && (ok == 1)) {
      nn++;
      if (str[pos+nn] == 0) { ok = 0; }
    }
    
    /* Create the local string */
    j = 0;
    for (i=pos; i<(pos+nn); i++) {
      lstr[j] = str[i];
      j++;
    }
    lstr[j] = 0;
    
    /* Remove the non-numeric characters in this string. */
    for (i=0; i<nn; i++) { 
      if ((lstr[i]<48) || (lstr[i]>57)) { lstr[i] = ' '; }
    }
    
    /* Extract limits */
    i = sscanf(lstr, "%i%i", &lo, &hi);
    if (i==1) { hi = lo; }

    /* Check limits */
    if (lo < 0) { 
      printf("TRACK: strDecomp: %i is less than lower limits of zero\n", lo);
      exit(0);
    }
    if (hi < 0) {
      printf("TRACK: strDecomp: %i is less than lower limits of zero\n", hi);
      exit(0);
    }
    if (lo >= dim) {
      printf("TRACK: strDecomp: %i is higher than limit of %i\n", lo, dim);
      exit(0);
    } 
    if (hi >= dim) {
      printf("TRACK: strDecomp: %i is higher than limit of %i\n", hi, dim);
      exit(0);
    }
    
    /* Fill array */
    for (j=lo; j<=hi; j++) {
      yy[j] = 1;
    }
    
    /* Move position pointer forward */
    pos += nn;
    pos++;
  }

  /* Done */
  return(0);
}

void TrackCtrl::PrintPartial() {
  printf("\nTRACK: ---------------------------------------\n");
  printf("TRACK: Tracking parameters (incomplete list):\n");
  printf("TRACK: \n");
  printf("TRACK: Alpha =  %f\n", alpha[2]);
  printf("TRACK: FOM jump = %f\n", FOMJump);
  printf("TRACK: FOM good enough = %f\n", FOMGoodEnough);
  printf("TRACK: nPrint = %i\n", nPrint);
  printf("TRACK: ---------------------------------------\n\n");
}

Int_t TrackCtrl::trackOpt(Int_t i) {
  switch(i) {
  case 0:
    // printf("TRACK: Full tracking\n");  
    return 0;  break;    
  case 1:
    // printf("TRACK: Largest E first\n");  
    return 0;  break;
  case 2:
    // printf("TRACK: Random point\n");  
    return 0;  break;
  case 3:
    // printf("TRACK: 'Kickout'\n");  
    return 0;  break;
  case 4:
    // printf("TRACK: 'GoodEnough'\n");  
    return 0;  break;
  case 5:
    // printf("TRACK: 'Hop/Skip/Jump'\n");  
    return 0;  break;
  default:
    printf("TRACK: Unknown\n");  return 1;  break;
  }
};

/**********************************************************/

void TrackStat::Zero() {
  nEvents = 0;
  trackGetCalls = 0;
  trackPassCnt = 0;
  trackingCalls = 0;
  firstClusterTrackCalls = 0;
  trackFOMokay = 0;
  singleHitOutOfRange = 0;
  badNDetELim = 0;
  sumFOM = 0.0;
  sumFOM_NT = 0.0;
  nSumFOM_NT = 0;
  nClusters = 0;
  notTracked = 0;
  nPerm = 0;
  doubleSegHits = 0;
  segHits = 0;
  numInputDataLines = 0;
  goodTrack = 0;
  badTrack = 0;
  noTrueClusters = 0;
  nDetTooLarge = 0;
  fullEnergyEv = 0;
  comptonLossEnergyEv = 0;
  gammasForOutput = 0;
  gammasForOutputOK1 = 0;
  gammasForOutputOK2 = 0;
  gammasWritten = 0;
  trackedEventsWritten = 0;

  for (Int_t i=0; i<NRETURNCODES; i++) {
    for (Int_t j=0; j<6; j++) {
      if (j < 2) {
	reClusterReturns[j][i] = 0;
	reClusterTrackCalls[j][i] = 0;
	splitClusterReturns[j][i] = 0;
	splitClusterTrackCalls[j][i] = 0;
      }
      trackError[j][i] = 0;
    }
    combineClusterReturns[i] = 0;
    combineClusterTrackCalls[i] = 0;
    matchMakerReturns[i] = 0;
    matchMakerTrackCalls[i] = 0;
  }

  for (Int_t i=0; i<164; i++) {
    nHit[i] = 0;
  }
  for (Int_t i=0; i<MAXNUMSEG; i++) {
    for (Int_t j=0; j<MAXPERM; j++) {
      permHit[i][j] = 0;
    }
  }
}

/**********************************************************/

void TrackPerm::Initialize() {
 
  /* Precalculate possible # of permutations */

  Int_t i1 = 0;
  for (Int_t j=1; j<MAXNUMSEG; j++) {
    i1 = 1;
    for (Int_t i=1; i<=j; i++) {
      i1 *= i;
      nPerm[j] = i1;
    }
  }

  if (0) {
    for (Int_t j=1; j<MAXNUMSEG; j++) {
      printf("TRACK:DEBUG: TrackPerm::Initialize(): nDet = %2i, #permutations = %10i\n", j, nPerm[j]);
    }
  }

  /* Make permutation lookup tables */
  
  /* Calculate permutation lookup table, currently
     lexographic ordered -- haven't found elegant way
     to do this, but only done once, so not critical */
  
  /* 1: nDet, selects the order,
     2: nPerm, permutation number,
     3: actual permutation for this permutation number */
  
  Int_t ndet = 1;
  lookup[ndet][0][0] = 0;
  ndet = 2;
  Int_t nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) {
	lookup[ndet][nn][0] = i1;
	lookup[ndet][nn][1] = i2;
	nn++;
      }
  ndet = 3; nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) 
	for (Int_t i3=0; i3<ndet; i3++) 
	  if (i3 != i1 && i3 != i2) {
	    lookup[ndet][nn][0] = i1;
	    lookup[ndet][nn][1] = i2;
	    lookup[ndet][nn][2] = i3;
	    nn++;
	  }
  ndet = 4; nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) 
	for (Int_t i3=0; i3<ndet; i3++) 
	  if (i3 != i1 && i3 != i2) 
	    for (Int_t i4=0; i4<ndet; i4++) 
	      if (i4 != i1 && i4 != i2 && i4 != i3) {
		lookup[ndet][nn][0] = i1;
		lookup[ndet][nn][1] = i2;
		lookup[ndet][nn][2] = i3;
		lookup[ndet][nn][3] = i4;
		nn++;
	      }
  ndet = 5; nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) 
	for (Int_t i3=0; i3<ndet; i3++) 
	  if (i3 != i1 && i3 != i2) 
	    for (Int_t i4=0; i4<ndet; i4++) 
	      if (i4 != i1 && i4 != i2 && i4 != i3) 
		for (Int_t i5=0; i5<ndet; i5++) 
		  if (i5 != i1 && i5 != i2 && i5 != i3 && i5 != i4) {
		    lookup[ndet][nn][0] = i1;
		    lookup[ndet][nn][1] = i2;
		    lookup[ndet][nn][2] = i3;
		    lookup[ndet][nn][3] = i4;
		    lookup[ndet][nn][4] = i5;
		    nn++;
		  }
  ndet = 6; nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) 
	for (Int_t i3=0; i3<ndet; i3++) 
	  if (i3 != i1 && i3 != i2) 
	    for (Int_t i4=0; i4<ndet; i4++) 
	      if (i4 != i1 && i4 != i2 && i4 != i3) 
		for (Int_t i5=0; i5<ndet; i5++) 
		  if (i5 != i1 && i5 != i2 && i5 != i3 && i5 != i4) 
		    for (Int_t i6=0; i6<ndet; i6++) 
		      if (i6 != i1 && i6 != i2 && i6 != i3 &&
			  i6 != i4 && i6 != i5) {
			lookup[ndet][nn][0] = i1;
			lookup[ndet][nn][1] = i2;
			lookup[ndet][nn][2] = i3;
			lookup[ndet][nn][3] = i4;
			lookup[ndet][nn][4] = i5;
			lookup[ndet][nn][5] = i6;
			nn++;
		      }
  ndet = 7; nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) 
	for (Int_t i3=0; i3<ndet; i3++) 
	  if (i3 != i1 && i3 != i2) 
	    for (Int_t i4=0; i4<ndet; i4++) 
	      if (i4 != i1 && i4 != i2 && i4 != i3) 
		for (Int_t i5=0; i5<ndet; i5++) 
		  if (i5 != i1 && i5 != i2 && i5 != i3 && i5 != i4) 
		    for (Int_t i6=0; i6<ndet; i6++) 
		      if (i6 != i1 && i6 != i2 && i6 != i3 &&
			  i6 != i4 && i6 != i5) 
			for (Int_t i7=0; i7<ndet; i7++) 
			  if (i7 != i1 && i7 != i2 && i7 != i3 &&
			      i7 != i4 && i7 != i5 && i7 != i6) {
			    lookup[ndet][nn][0] = i1;
			    lookup[ndet][nn][1] = i2;
			    lookup[ndet][nn][2] = i3;
			    lookup[ndet][nn][3] = i4;
			    lookup[ndet][nn][4] = i5;
			    lookup[ndet][nn][5] = i6;
			    lookup[ndet][nn][6] = i7;
			    nn++;
			  }
  ndet = 8; nn = 0;
  for (Int_t i1=0; i1<ndet; i1++) 
    for (Int_t i2=0; i2<ndet; i2++) 
      if (i2 != i1) 
	for (Int_t i3=0; i3<ndet; i3++) 
	  if (i3 != i1 && i3 != i2) 
	    for (Int_t i4=0; i4<ndet; i4++) 
	      if (i4 != i1 && i4 != i2 && i4 != i3) 
		for (Int_t i5=0; i5<ndet; i5++) 
		  if (i5 != i1 && i5 != i2 && i5 != i3 && i5 != i4) 
		    for (Int_t i6=0; i6<ndet; i6++) 
		      if (i6 != i1 && i6 != i2 && i6 != i3 &&
			  i6 != i4 && i6 != i5) 
			for (Int_t i7=0; i7<ndet; i7++) 
			  if (i7 != i1 && i7 != i2 && i7 != i3 &&
			      i7 != i4 && i7 != i5 && i7 != i6) 
			    for (Int_t i8=0; i8<ndet; i8++) 
			      if (i8 != i1 && i8 != i2 && i8 != i3 &&
				  i8 != i4 && i8 != i5 && i8 != i6 &&
				  i8 != i7) {
				lookup[ndet][nn][0] = i1;
				lookup[ndet][nn][1] = i2;
				lookup[ndet][nn][2] = i3;
				lookup[ndet][nn][3] = i4;
				lookup[ndet][nn][4] = i5;
				lookup[ndet][nn][5] = i6;
				lookup[ndet][nn][6] = i7;
				lookup[ndet][nn][7] = i8;
				nn++;
			      }
}

/**********************************************************/


void Track::Initialize() {
  
  {  /* Previously in setupTrack() function */ 

    /* Zero out the statistics for tracking */
    stat.Zero();
    
    {  /* Previously in ctkinit() function */
      ctrl.badTheoAngPenalty = 0.0;
      ctrl.reTrack = 0;
      
      for (Int_t i=0; i<MAXNUMSEG; i++) { 
	ctrl.nDetELim_lo[i] = 0;
	ctrl.nDetELim_hi[i] = 1000.0;
      }
      
      /* Initialize permutation tables */
      perm.Initialize();
    }

    /* Faculty function lookup */
    printf("\nTRACK: Faculty function being calculated.....");
    faculty[0] = 1;
    faculty[1] = 1;
    for (Int_t i=2; i<MAXFACULTY; i++) { 
      faculty[i] = (long long) i * faculty[i-1]; 
    }
    printf("Done.\n\n");
 
  }  
  
  ctrl.ReadChatFile("track.chat");
  
  {  /* Previously in setupTrack() function */   
    
    /* Sanity checks */
    if (ctrl.useCCEnergy == 1 && ctrl.useSegEnergy == 1) {
      cout << ALERTTEXT;
      printf("TRACK: Error: only one of ctrl.useCCEnergy and ctrl.useSegEnergy can be true\n");
      cout << RESET_COLOR;  fflush(stdout);
      exit(1);
    }
    
    /* Initialize parameter needed in cluster -- crap way of doing this
       at the moment. */
    for (Int_t i=0; i<MAXCLUSTERHITS; i++) {
      clust[i].badTheoAngPenalty = ctrl.badTheoAngPenalty;
    }

  /* Print tracking parameters that will be used */
   ctrl.PrintPartial();
  
  }
}

  
Int_t Track::findTargetPos() {
  targetPos.SetXYZ(ctrl.targetX, ctrl.targetY, ctrl.targetZ);
  return (0);
}

Int_t Track::trackEvent() {

  Int_t i, j, i1, st, itNum;
  
  stat.nEvents++;

  if (DEBUG) {
    printf("TRACK: DEBUG: Entered trackEvent(): %d\n", stat.nEvents);
  }

  /* Replace bad double hits in the same sector
     with single hits instead, compensating for 
     problems in the decomposition. */
  if (ctrl.fixDoubleHits) {
    // st = fixDoubleHitsInSeg();
    if (st != 0) {
      return (st);
    }
  }

  /* Clustering! */
  findClusters();

  if (DEBUG) {
    printf("TRACK: DEBUG: After findClusters --> shell:\n");
    shell.Print();
    shell.PrintClusters();
    cin.get();
  }

  transferClusters();

  /* Right now the cluster array only has valid entries, 
     since nothing has been rearranged in reclustering. */
  nClusters = shell.numClusters;

  if (DEBUG) {
    printf("TRACK: DEBUG: After transferClusters\n");
    PrintValidClusters();
    cin.get();
  }

  /* Sort segment energies, so groups can be skipped properly */
  for (i=0; i<nClusters; i++) { clust[i].sortSegEnergies(); }

  /* Track each cluster */
  for (i=0; i<nClusters; i++) {
    if (clust[i].valid) {
      stat.trackGetCalls++;
      
      /* Hit statistics */
      stat.nHit[clust[i].nDet]++;
      stat.trackPassCnt++;

      /* Find best interaction point sequence */
      i1 = -stat.trackingCalls;
      clust[i].tracked = 0;
      if (clust[i].nDet < MAXNUMSEG) {
	if (ctrl.trackOps[clust[i].nDet] < 6 && ctrl.trackOps[clust[i].nDet] >= 0 && ctrl.trackOps[clust[i].nDet] != 2) {
	  st = doTrack(ctrl.trackOps[clust[i].nDet], i);
	  stat.trackError[ctrl.trackOps[clust[i].nDet]][st]++;
	} else {
	  printf("TRACK: ctk: tracking option not known!?, option = %i\n, TRACK: Quitting...\n", ctrl.trackOps[clust[i].nDet]);
	  exit(1);
	}
      }

      i1 += stat.trackingCalls;
      stat.firstClusterTrackCalls += i1;

      /* Count good tracks with clust[i].FOM below cut */
      stat.trackFOMokay++;

      /* What permutations work? */
      stat.permHit[clust[i].nDet][clust[i].bestPermutation]++;
      
      /* Sanity checks */
      if (clust[i].valid) {
	if (clust[i].tracked) {
	  if (clust[i].FOM > MAXFOM) {
	    printf("TRACK: ERROR: FOM makes no sense = %f\n", clust[i].FOM);
	    printf("TRACK: iCluster = %i\n", i);
	    printf("TRACK: ctrl.trackOps[nDet] = %i\n", ctrl.trackOps[clust[i].nDet]);
	  }
	}
      }
    }
  } /* Loop over clusters */
  
  if (DEBUG) {
    printf("TRACK: DEBUG: After firstTracking\n");
    PrintValidClusters();
    cin.get();
  }

  Int_t success = 0;
  Int_t numCurrentClusters = nClusters;

  /* Now we have formed and tracked clusters -- need to
     split and add to see if improvement is possible.
     We will try to iteratively improve. */
  
  for (itNum = 0; itNum < ctrl.iterations; itNum++) {
    /* Look for clusters that are so big that they weren't tracked, 
       and split them. */
    if (ctrl.reCluster[0]) {
      success = 0;
      numCurrentClusters = nClusters;
      for (i=0; i<numCurrentClusters; i++) {
	if (clust[i].valid && !clust[i].tracked && 
	    clust[i].nDet >= ctrl.reCluster_minN[0]) {
	  i1 = -stat.trackingCalls;
	  st = reCluster(i, 0); /* Actual recluster/track step */
	  i1 += stat.trackingCalls;
	  stat.reClusterReturns[0][st]++;
	  stat.reClusterTrackCalls[0][st] += i1;
	  if (st != 8) { success += 1; }
	}
      }
    } /* reCluster[0] active */
  
    /* Try to split monster clusters using splitCluster */
    if (ctrl.splitClusters[0]) {
      success = 0;
      numCurrentClusters = nClusters;
      for (i=0; i<numCurrentClusters; i++) {
	if (clust[i].valid && !clust[i].tracked) {
	  i1 = -stat.trackingCalls;
	  st = splitCluster(i, 0);
	  i1 += stat.trackingCalls;
	  stat.splitClusterReturns[0][st]++;
	  stat.splitClusterTrackCalls[0][st] += i1;
	  if (st != 18) { success = 1; }
	}
      }
    } /* splitCluster[0] active */
      
    /* Now take a stab at clusters that most likely are too small,
       and combine them. 
       We try singles hits no matter what since they have no FOM. */
    if (ctrl.combineClusters) {
      success = 0;
      numCurrentClusters = nClusters;
      for (i=0; i<numCurrentClusters; i++) {
	if (clust[i].valid && clust[i].FOM > ctrl.combineClusters_thresh &&
	    clust[i].nDet <= ctrl.combineClusters_maxN) {
	  i1 = -stat.trackingCalls;
	  st = combineCluster(i);
	  i1 += stat.trackingCalls;
	  stat.combineClusterReturns[st]++;
	  stat.combineClusterTrackCalls[st] += i1;
	  if (st == 0) { success = 1; }
	}
      } 
    } /* combineCluster active */

    /* And another attempt to split big clusters */
    if (ctrl.reCluster[1]) {
      success = 0;
      numCurrentClusters = nClusters;
      for (i=0; i<numCurrentClusters; i++) {
	if (clust[i].valid && !clust[i].tracked && 
	    clust[i].FOM > ctrl.reCluster_thresh[1] &&
	    clust[i].nDet >= ctrl.reCluster_minN[1]) {
	  i1 = -stat.trackingCalls;
	  st = reCluster(i, 1); /* Actual recluster/track step */
	  i1 += stat.trackingCalls;
	  stat.reClusterReturns[1][st]++;
	  stat.reClusterTrackCalls[1][st] += i1;
	  if (st != 8) { success += 1; }
	}
      }
    } /* reCluster[1] active */

    if (ctrl.splitClusters[1]) {
      success = 0;
      numCurrentClusters = nClusters;
      for (i=0; i<numCurrentClusters; i++) {
	if (clust[i].valid && clust[i].FOM > ctrl.splitClusters_thresh[1] &&
	    clust[i].nDet >= ctrl.splitClusters_minN[1] &&
	    clust[i].nDet <= ctrl.splitClusters_maxN[1]) {
	  i1 = -stat.trackingCalls;
	  st = splitCluster(i, 1);
	  i1 += stat.trackingCalls;
	  stat.splitClusterReturns[1][st]++;
	  stat.splitClusterTrackCalls[1][st] += i1;
	  if (st != 18) { success = 1; }
	}
      }
    } /* splitCluster[1] active */

    if (ctrl.matchMaker) {
      /* Look for single hits within a threshold distance and
	 try to combine them. */
      success = 0;
      numCurrentClusters = nClusters;
      for (i=0; i<numCurrentClusters; i++) {
	if (clust[i].nDet == 1 && clust[i].valid) {
	  for (j=i+1; j<numCurrentClusters; j++) {
	    if (clust[j].nDet == 1 && clust[j].valid) {
	      i1 = -stat.trackingCalls;
	      st = matchMaker(i, j);
	      i1 += stat.trackingCalls;
	      stat.matchMakerReturns[st]++;
	      stat.matchMakerTrackCalls[st] += i1;
	      if (st == 0) { success = 1; }
	    }
	  }
	}
      }
    }
  } /* End of iteration loop */

  if (DEBUG) {
    cout << DRED;
    printf("TRACK: DEBUG: After iterationLoop\n");
    PrintValidClusters();
    cout << RESET_COLOR;  fflush(stdout);
    cin.get();
  }

  /* Mark singles hits that are out of range */
  if (ctrl.markSingleHitRange) {
    numCurrentClusters = nClusters;
    for (i=0; i<numCurrentClusters; i++) {
      if (clust[i].valid && clust[i].nDet == 1) {
	i1 = (Int_t)(clust[i].eSum*1000);
	if (i1 > 0 && i1 < LONGLEN) {
	  if (clust[i].intpts[0].xyz.Mag2() > ctrl.singleHitRange[i1]) {
	    clust[i].FOM = ctrl.singleHitRangeFOM;
	    stat.singleHitOutOfRange++;
	  }
	}
      }
    }
  }
  
  /* Mark gamma-rays with a silly number of interaction points */
  numCurrentClusters = nClusters;
  for (i=0; i<numCurrentClusters; i++) {
    if (clust[i].valid) {
      if ((clust[i].eSum < ctrl.nDetELim_lo[clust[i].nDet]) ||
	  (clust[i].eSum > ctrl.nDetELim_hi[clust[i].nDet])) {
	clust[i].FOM = ctrl.nDetELim_FOM[clust[i].nDet];
	stat.badNDetELim++;
      }

      /* More tracking statistics */
      if (clust[i].tracked) {
	stat.sumFOM += clust[i].FOM;
	if (clust[i].nDet > 1) {
	  stat.sumFOM_NT += clust[i].FOM;
	  stat.nSumFOM_NT++;
	}
      }
      stat.nClusters++;
      if (!clust[i].tracked) {
	stat.notTracked++;
      }
    }
  }

  /* Done */
  return (0);
}

Int_t Track::doTrack(Int_t mode, Int_t iClust) {
  cout << DCYAN;

  Float_t FOM, minFOM;
  Int_t curPerm, k;
  Int_t index[MAXNUMSEG];
  Float_t eg, th, thc;
  TVector3 v1, v2;
  Int_t j;
  Float_t pFac;

  Float_t oldFOM = 0; /* For use with mode 3 (kickout) */
  Float_t FOMcutGoodEnough; /* For use with mode 4 (goodEnough) */
  /* For use with mode 5 (jump) */
  Float_t FOMcutJump;
  Int_t distNextGrp, baseJump;
  
  stat.trackingCalls++;

  /* Initialize */
  clust[iClust].FOM = 0.;
  
  /* Trap single hits, which are trivial, but
     we need to mark them as tracked anyway. */
  if (clust[iClust].nDet == 1) {
    clust[iClust].tracked = 1;
    clust[iClust].FOM = ctrl.singlesFOM;
    clust[iClust].bestPermutation = 0;
    clust[iClust].intpts[0].order = 0;
    return (1);
  }

  /* Trap if nDet is higher than we can handle */
  if (clust[iClust].nDet >= MAXNUMSEG) {
    clust[iClust].tracked = 0;
    clust[iClust].FOM = MAXFOM;
    return (2);
  } else {
    clust[iClust].tracked = 1;
  }

  if (mode == 4) { /* goodEnough */
    FOMcutGoodEnough = clust[iClust].nDet * ctrl.FOMGoodEnough;
  }

  /* If mode 1 (maxE), just return the order we have,
     since  we sorted by energy before */
  if (mode == 1) {
    j = 0;
    for (k=0; k<clust[iClust].nDet; k++) {
      clust[iClust].intpts[k].order = j;
      j++;
    }
    clust[iClust].tracked = 1;

    /* Assign a proper FOM for the cluster in energy order */
    clust[iClust].findFOM(targetPos);

    /* We are done with this! */
    return (8);
  } /* mode 1 */

  /* Consider mode 5 now, the jump strategy, define variables */
  if (mode == 5) {
    FOMcutJump = clust[iClust].nDet*ctrl.FOMJump;
    /* Distance to the next group, decided by size of tail */
    baseJump = faculty[clust[iClust].nDet - ctrl.jumpGroupLength[clust[iClust].nDet]];
    distNextGrp = baseJump;
 
    if (DEBUG) {
      printf("TRACK: DEBUG: doTrack(): Mode5: FOMcutJump, baseJump, distNextGrp = %0.3f, %d, %d\n",
	     FOMcutJump, baseJump, distNextGrp);
    }
 }

  /* Find sum of interaction energies */
  clust[iClust].eSum = 0;
  for (k=0; k<clust[iClust].nDet; k++) { 
    clust[iClust].eSum += clust[iClust].intpts[k].eDet; 
  }

  /* Loop over ALL permutations, i.e. full search, not tree search */

  minFOM = FLT_MAX;
  for (curPerm = 0; curPerm < perm.nPerm[clust[iClust].nDet]; curPerm++) {

    if (DEBUG) {
      printf("curPerm, FOM, minFOM: %d, %0.3f, %0.3f\n", curPerm, FOM, minFOM);
      printf(" Order: \n");
      for (k=0; k<clust[iClust].nDet; k++) {
	printf("   %d: %d\n", k, perm.lookup[clust[iClust].nDet][curPerm][k]);
      }
      printf("\n");
    }

    /* Kickout mode */
    if (mode == 3) { if (curPerm > 0) { oldFOM = minFOM; } }
   
    /* Count permutations we make */
    stat.nPerm++;

    eg = clust[iClust].eSum;

    /* Index for FOM calculation */
    for (k=0; k<clust[iClust].nDet; k++) {
      clust[iClust].intpts[k].order = perm.lookup[clust[iClust].nDet][curPerm][k];
    }
   
    clust[iClust].findFOM(targetPos);
    FOM = clust[iClust].FOM;

    /* Keep track of the best permutation so far */
    if (FOM < minFOM) {
      minFOM = FOM;
      clust[iClust].bestPermutation = curPerm;
      clust[iClust].FOM = FOM;
    }

    if (mode == 3) { /* kickout -- if we're doing worse than 
			last permutation, kickout here */
      if (curPerm > 0 && FOM > oldFOM) {
	/* Normalize FOM, and return minFOM */
	minFOM /= (clust[iClust].nDet - 1);
	clust[iClust].FOM = minFOM;
	for (k=0; k<clust[iClust].nDet; k++) {
	  clust[iClust].intpts[perm.lookup[clust[iClust].nDet][clust[iClust].bestPermutation][k]].order = k;
	  //clust[iClust].intpts[k].order = perm.lookup[clust[iClust].nDet][clust[iClust].bestPermutation][k];
	}
	return (3);
      }
    }

    if (mode == 4) { /* goodEnough -- leave if FOM is below cut */
      if (FOM < FOMcutGoodEnough) {
	/* Normalize FOM, and return minFOM */
	minFOM /= (clust[iClust].nDet - 1);
	clust[iClust].FOM = minFOM;
	for (k=0; k<clust[iClust].nDet; k++) {
	  clust[iClust].intpts[perm.lookup[clust[iClust].nDet][clust[iClust].bestPermutation][k]].order = k;
	  // clust[iClust].intpts[k].order = perm.lookup[clust[iClust].nDet][clust[iClust].bestPermutation][k];
	}
	return (4);
      }
    }

    if (mode == 5) { /* jump -- jump ahead if FOM is over FOM cut */

      if (FOM > ctrl.FOMJump) {
	/* Advance curPerm to the next group */
	curPerm += (distNextGrp - 1);
	/* Reset distNextGroup */
	distNextGrp = baseJump;
      } else {
	/* Find the distance to the next group */
	distNextGrp--;
	if (distNextGrp <= 0) { distNextGrp = baseJump; }
      }
    }
  
  } /* Loop over permutations */
  
  /* Normalize FOM, and return minFOM */
  minFOM = sqrtf(minFOM) / (clust[iClust].nDet - 1);
  // minFOM /= (clust[iClust].nDet - 1);
  clust[iClust].FOM = minFOM;
  for (k=0; k<clust[iClust].nDet; k++) {
    clust[iClust].intpts[perm.lookup[clust[iClust].nDet][clust[iClust].bestPermutation][k]].order = k;
    //clust[iClust].intpts[k].order = perm.lookup[clust[iClust].nDet][clust[iClust].bestPermutation][k];
  }

  if (DEBUG) {
    printf("TRACK: DEBUG: Finished doTrack: minFOM = %0.3f\n", minFOM);
    cin.get();
  }

  cout << RESET_COLOR;

  /* Done */
  return (0);
}

Int_t Track::reCluster(Int_t ii, Int_t m) {

  /* This function is similar to splitCluster, but based on doing
     a recluster with a smaller alpha value than original. 
     For some clusters this should be okay, and it will be 
     faster than splitCluster.  But if two/three clusters assumed
     to make up this big cluster are too close, it won't work. */

  Int_t nn, i, j, k, l, opi, opj, i1, i2;
  Int_t clusterPos[MAXSHELLHITS];
  Float_t curAlpha, sumNewFOM, smallestNonTrivial = 0, biggestNonTrivial = 0;
  Int_t trialPos[MAXSHELLHITS], oldTrialPos[MAXSHELLHITS], nDiffs = 0;
  Int_t nTrialClusters = 0;
  Int_t trackRtn[MAXSHELLHITS];
  
  /* Initialize */
  for (i=0; i<clust[ii].nDet; i++) { oldTrialPos[i] = -1; }
  curAlpha = ctrl.alpha[shell.nHit];
  
  /* Loop through smaller cluster angles and
     find a split that makes sense. */
  for (nn=0; nn<ctrl.reCluster_maxTry[m]; nn++) {
    nTrialClusters = 0;
    for (i=0; i<clust[ii].nDet; i++) {
      clusterPos[i] = -1;  /* Mark as not found */
      oldTrialPos[i] = 0;  /* Trivial old position */
    }
    
    /* Cluster angle for this iteration */
    curAlpha *= ctrl.reCluster_reduxFactor[m];

    /* Work on multi-member clusters */
    for (i=0; i<clust[ii].nDet; i++) {
      for (j=0; j<clust[ii].nDet; j++) {
	if (i != j) {
	  if (clust[ii].intpts[i].xyz.Angle(clust[ii].intpts[j].xyz) < curAlpha) {

	    if (clusterPos[i] < 0 && clusterPos[j] < 0) { /* New cluster! */
	      clusterPos[i] = nTrialClusters;
	      clusterPos[j] = nTrialClusters;
	      nTrialClusters++;
	    } else if (clusterPos[i] >= 0 && 
		       clusterPos[j] < 0) { /* Join cluster */
	      clusterPos[j] = clusterPos[i];
	    } else if (clusterPos[i] < 0 && 
		       clusterPos[j] >= 0) { /* Join clusters */
	      clusterPos[i] = clusterPos[j];
	    } else if (clusterPos[i] >= 0 && clusterPos[j] >= 0 &&
		       clusterPos[i] != clusterPos[j]) { /* Clusters collide */
	      /* i1 = common index, i2 will move to i1 */
	      if (clusterPos[i] > clusterPos[j]) {
		i1 = clusterPos[j];  i2 = clusterPos[i];
	      } else { i1 = clusterPos[i];  i2 = clusterPos[j]; }
	      
	      /* Merge to lower numbered cluster */
	      for (k=0; k<clust[ii].nDet; k++) {
		if (clusterPos[k] == i2) { clusterPos[k] = i1; }
	      }
	      for (k=0; k<clust[ii].nDet; k++) {
		if (clusterPos[k] > i2) { clusterPos[k]--; }
	      }
	      nTrialClusters--;
	    }
	  }
	}
      }
    }
  
    /* Assign single hit clusters to the rest that don't have a partner */
    for (k=0; k<clust[ii].nDet; k++) {
      if (clusterPos[k] < 0) {
	clusterPos[k] = nTrialClusters;
	nTrialClusters++;
      }
    }
    
    /* Fail reclustering if we only found single clusters since they 
       will always have FOM = 0.0, which is a problem! */
    if (nTrialClusters == clust[ii].nDet) { return (8); }

    /* See if clustering is different from the previous one. */
    nDiffs = 0;
    for (i=0; i<clust[ii].nDet; i++) {
      if (clusterPos[i] != oldTrialPos[i]) {
	nDiffs++;
      }
    }

    /* Do we have more clusters that are different? */
    if (nDiffs > 0 && nTrialClusters > 0) {
      for (i=0; i<clust[ii].nDet; i++) {
	oldTrialPos[i] = clusterPos[i];
      }
    
      /* Make the trial clusters */
      for (k=0; k<nTrialClusters; k++) {
	trialPos[k] = k + nClusters;
	clust[trialPos[k]].nDet = 0;
	clust[trialPos[k]].eSum = 0;
	clust[trialPos[k]].tracked = 0;
	clust[trialPos[k]].valid = 1;
	clust[trialPos[k]].FOM = 98;
	for (l=0; l<MAXNUMDET; l++) {
	  clust[trialPos[k]].intpts[l].order = -1;
	}
      }
      /* Fill the trial clusters */
      for (l=0; l<clust[ii].nDet; l++) {
	k = clusterPos[l];
	clust[trialPos[k]].eSum += clust[ii].intpts[l].eDet;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].xyz = 
	  clust[ii].intpts[l].xyz;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].eDet = 
	  clust[ii].intpts[l].eDet;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].timestamp = 
	  clust[ii].intpts[l].timestamp;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].t0 = 
	  clust[ii].intpts[l].t0;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].detNum = 
	  clust[ii].intpts[l].detNum;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].eSumOrig = 
	  clust[ii].intpts[l].eSumOrig;
	clust[trialPos[k]].intpts[clust[trialPos[k]].nDet].origPos = 
	  clust[ii].intpts[l].origPos;
	clust[trialPos[k]].nDet++;
      }

      /* Track the trial clusters */
      for (l=0; l<nTrialClusters; l++) {
	k = clusterPos[l];
	clust[trialPos[l]].sortSegEnergies(); 	
	trackRtn[l] = -1;

	Int_t st = -1;
	if (clust[k].nDet < MAXNUMSEG) {
	  if (ctrl.trackOps[clust[k].nDet] < 6 && 
	      ctrl.trackOps[clust[k].nDet] >= 0 && 
	      ctrl.trackOps[clust[k].nDet] != 2) {
	    st = doTrack(ctrl.trackOps[clust[k].nDet], trialPos[l]);
	    stat.trackError[ctrl.trackOps[clust[k].nDet]][st]++;
	  } else {
	    printf("TRACK: reCluster: tracking option not known!?, option = %i\n, TRACK: Quitting...\n", ctrl.trackOps[clust[k].nDet]);
	    exit(1);
	  }
	}
      }

      /* Are we doing better?  Several checks below. */
      sumNewFOM = 0.;
      smallestNonTrivial = 999.;  biggestNonTrivial = 0.;
      for (l=0; l<nTrialClusters; l++) {
	i1 = trialPos[l];
	sumNewFOM += clust[i1].FOM;
	if (clust[i1].nDet > 1 && clust[i1].FOM < smallestNonTrivial) {
	  smallestNonTrivial = clust[i1].FOM;
	} 
	if (clust[i1].nDet > 1 && clust[i1].FOM > biggestNonTrivial) {
	  biggestNonTrivial = clust[i1].FOM;
	}
      }
      if (sumNewFOM < clust[ii].FOM) {
	/* Always require that the sum of new FOM is better than it was, 
	   but that is not all. -- May still need more options to declare 
	   success in this section */
	if (sumNewFOM < ctrl.reCluster_kickoutFOM[m]) {
	  /* Kickout success.  Just update nClusters since
	     the trial clusters are valid and in the right spot. */
	  nClusters += nTrialClusters;
	  /* Invalidate the old cluster */
	  clust[ii].valid = 0;
	  return (0);
	} else if (smallestNonTrivial < ctrl.reCluster_kickoutFOM[m]) {
	  /* At least one good non-trivial cluster is now good,
	     declare victory here */
	  nClusters += nTrialClusters;
	  /* Invalidate the old cluster */
	  clust[ii].valid = 0;
	  return (1);
	} else {
	  /* Invalidate the trial clusters and continue. */
	  for (l=0; l<nTrialClusters; l++) {
	    clust[trialPos[l]].valid = 0;
	  }
	}
      } else {
	/* Invalidate the trial clusters and continue. */
	for (l=0; l<nTrialClusters; l++) {
	  clust[trialPos[l]].valid = 0;
	}
      }
    }
  }

  /* Done -- if we get here, we failed to recluster, but we might
     still be doing better... unclear how to proceed here. */
  return (8);
}

Int_t Track::splitCluster(Int_t ii, Int_t m) {
  /* Split a big cluster into two clusters... */
  Int_t nBaseSplits, i, j, k, l, nWays, i0, i1, i2, i3, i4, i5, i6, nSplit = 0;
  Float_t e1, e2, sumFOM, bestSplitFOM;
  Float_t dst[MAXNUMDET][MAXNUMDET], r1, r2;
  Int_t n1, n2, pNum, try1, try2, bestSplit = -1, best1, best2;
  Int_t nTracksTried = 0;
  Int_t order[MAXPNUM], st;
  PM pm[MAXPNUM], tmppm;

  Float_t goodEnough;
  if (m == 0) { goodEnough = ctrl.untracked_kickoutFOM; }
  else if (m == 1) { goodEnough = clust[ii].FOM * ctrl.goodImproveFraction; }

  /* Check to make sure that we are not out of bounds. */
  if (clust[ii].nDet >= MAXFACULTY) {
    printf("TRACK: splitCluster: nDet = %i is too big to be handled\n", clust[i].nDet);
    return (5);
  }

  /* Initialize */
  bestSplitFOM = 100.;
  bestSplit = -1;
  pNum = 0;
  for (j=0; j<MAXPNUM; j++) { bzero((void*) &pm[j], sizeof(struct PM)); }

  /* For now we only consider splitting into two clusters. */
  nBaseSplits = (Int_t)(clust[ii].nDet / 2);
  pNum = 0;
  for (j=1; j<=nBaseSplits; j++) {
    if (j < MAXNUMSEG && (clust[ii].nDet - j) < MAXNUMSEG) {
      nWays = (Int_t) (faculty[clust[ii].nDet] / faculty[j] / faculty[clust[ii].nDet - j]);
      switch(j) {
      case 1:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
	  for (k=0; k<clust[ii].nDet; k++) {
	    if (k == i0) {
	      n1++;
	      pm[pNum].s1 += clust[ii].intpts[k].eDet;
	      pm[pNum].p1[pm[pNum].np1] = k;
	      pm[pNum].np1++;
	    } else {
	      n1++;
	      pm[pNum].s2 += clust[ii].intpts[k].eDet;
	      pm[pNum].p2[pm[pNum].np2] = k;
	      pm[pNum].np2++;
	    } 
	    pNum++;
	  }
	}
	break;
      case 2:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  for (i1=i0+1; i1<clust[ii].nDet; i1++) {
	    pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
	    for (k=0; k<clust[ii].nDet; k++) {
	      if (k == i0 || k == i1) {
		e1 += clust[ii].intpts[k].eDet;
		n1++;
		pm[pNum].s1 += clust[ii].intpts[k].eDet;
		pm[pNum].p1[pm[pNum].np1] = k;
		pm[pNum].np1++;
	      } else {
		e2 += clust[ii].intpts[k].eDet;
		n1++;
		pm[pNum].s2 += clust[ii].intpts[k].eDet;
		pm[pNum].p2[pm[pNum].np2] = k;
		pm[pNum].np2++;
	      } 
	      pNum++;
	    }
	  }
	}
	break;
      case 3:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  for (i1=i0+1; i1<clust[ii].nDet; i1++) {
	    for (i2=i1+1; i2<clust[ii].nDet; i2++) {
	      pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
	      for (k=0; k<clust[ii].nDet; k++) {
		if (k == i0 || k == i1 || k == i2) {
		  e1 += clust[ii].intpts[k].eDet;
		  n1++;
		  pm[pNum].s1 += clust[ii].intpts[k].eDet;
		  pm[pNum].p1[pm[pNum].np1] = k;
		  pm[pNum].np1++;
		} else {
		  e2 += clust[ii].intpts[k].eDet;
		  n1++;
		  pm[pNum].s2 += clust[ii].intpts[k].eDet;
		  pm[pNum].p2[pm[pNum].np2] = k;
		  pm[pNum].np2++;
		} 
		pNum++;
	      }
	    }
	  }
	}
	break;
      case 4:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  for (i1=i0+1; i1<clust[ii].nDet; i1++) {
	    for (i2=i1+1; i2<clust[ii].nDet; i2++) {
	      for (i3=i2+1; i3<clust[ii].nDet; i3++) {
		pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
		for (k=0; k<clust[ii].nDet; k++) {
		  if (k == i0 || k == i1 || k == i2 ||
		      k == i3) {
		    e1 += clust[ii].intpts[k].eDet;
		    n1++;
		    pm[pNum].s1 += clust[ii].intpts[k].eDet;
		    pm[pNum].p1[pm[pNum].np1] = k;
		    pm[pNum].np1++;
		  } else {
		    e2 += clust[ii].intpts[k].eDet;
		    n1++;
		    pm[pNum].s2 += clust[ii].intpts[k].eDet;
		    pm[pNum].p2[pm[pNum].np2] = k;
		    pm[pNum].np2++;
		  } 
		  pNum++;
		}
	      }
	    }
	  }
	}
	break;
      case 5:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  for (i1=i0+1; i1<clust[ii].nDet; i1++) {
	    for (i2=i1+1; i2<clust[ii].nDet; i2++) {
	      for (i3=i2+1; i3<clust[ii].nDet; i3++) {
		for (i4=i3+1; i4<clust[ii].nDet; i4++) {
		  pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
		  for (k=0; k<clust[ii].nDet; k++) {
		    if (k == i0 || k == i1 || k == i2 ||
			k == i3 || k == i4) {
		      e1 += clust[ii].intpts[k].eDet;
		      n1++;
		      pm[pNum].s1 += clust[ii].intpts[k].eDet;
		      pm[pNum].p1[pm[pNum].np1] = k;
		      pm[pNum].np1++;
		    } else {
		      e2 += clust[ii].intpts[k].eDet;
		      n1++;
		      pm[pNum].s2 += clust[ii].intpts[k].eDet;
		      pm[pNum].p2[pm[pNum].np2] = k;
		      pm[pNum].np2++;
		    } 
		    pNum++;
		  }
		}
	      }
	    }
	  }
	}
	break;
      case 6:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  for (i1=i0+1; i1<clust[ii].nDet; i1++) {
	    for (i2=i1+1; i2<clust[ii].nDet; i2++) {
	      for (i3=i2+1; i3<clust[ii].nDet; i3++) {
		for (i4=i3+1; i4<clust[ii].nDet; i4++) {
		  for (i5=i4+1; i5<clust[ii].nDet; i5++) {	
		    pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
		    for (k=0; k<clust[ii].nDet; k++) {
		      if (k == i0 || k == i1 || k == i2 ||
			  k == i3 || k == i4 || k == i5) {
			e1 += clust[ii].intpts[k].eDet;
			n1++;
			pm[pNum].s1 += clust[ii].intpts[k].eDet;
			pm[pNum].p1[pm[pNum].np1] = k;
			pm[pNum].np1++;
		      } else {
			e2 += clust[ii].intpts[k].eDet;
			n1++;
			pm[pNum].s2 += clust[ii].intpts[k].eDet;
			pm[pNum].p2[pm[pNum].np2] = k;
			pm[pNum].np2++;
		      } 
		      pNum++;
		    }
		  }
		}
	      }
	    }
	  }
	}
	break;
      case 7:
	for (i0=0; i0<clust[ii].nDet; i0++) {
	  for (i1=i0+1; i1<clust[ii].nDet; i1++) {
	    for (i2=i1+1; i2<clust[ii].nDet; i2++) {
	      for (i3=i2+1; i3<clust[ii].nDet; i3++) {
		for (i4=i3+1; i4<clust[ii].nDet; i4++) {
		  for (i5=i4+1; i5<clust[ii].nDet; i5++) {
		    for (i6=i5+1; i6<clust[ii].nDet; i6++) {	
		      pm[pNum].s1 = 0;  n1 = 0;  n2 = 0; pm[pNum].s2 = 0;
		      for (k=0; k<clust[ii].nDet; k++) {
			if (k == i0 || k == i1 || k == i2 ||
			    k == i3 || k == i4 || k == i5 || k == i6) {
			  e1 += clust[ii].intpts[k].eDet;
			  n1++;
			  pm[pNum].s1 += clust[ii].intpts[k].eDet;
			  pm[pNum].p1[pm[pNum].np1] = k;
			  pm[pNum].np1++;
			} else {
			  e2 += clust[ii].intpts[k].eDet;
			  n1++;
			  pm[pNum].s2 += clust[ii].intpts[k].eDet;
			  pm[pNum].p2[pm[pNum].np2] = k;
			  pm[pNum].np2++;
			} 
			pNum++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	break;
      default:
	printf("TRACK: splitCluster: too many base groups, event # %i, continuing\n",
	       stat.nEvents);
	break;
      }
    }
  }

  /* Go through the clusters and evaluate the split that gives the 
     lowest overall FOM compared to the total cluster */
  best1 = nClusters;  best2 = best1 + 1; /* ultimate homes for best clusters */
  clust[best1].valid = 0; clust[best2].valid = 0;
  
  try1 = best2 + 1;  try2 = try1 + 1;  /* homes for trial clusters */
  
  /* Try to find an optimum evaluation order */
  for (j=0; j<pNum; j++) {
    /* The more clustered the fragments are, the more likely they are */
    pm[j].likely1 = 1;
    r1 = 0.0;
    if (pm[j].np1 > 1) {
      for (k=0; k<pm[j].np1; k++) {
	for (l=k+1; l<pm[j].np1; l++) {
	  r1 += ((clust[ii].intpts[pm[j].p1[k]].xyz - 
		  clust[ii].intpts[pm[j].p1[l]].xyz).Mag());
	}
      }
      r1 /= (pm[j].np1 * (pm[j].np1 - 1.0) * 2);
    } else {
      r1 = clust[ii].nDet + 1.0;
    }

    r2 = 0.0;
    if (pm[j].np2 > 1) {
      for (k=0; k<pm[j].np2; k++) {
	for (l=k+1; l<pm[j].np2; l++) {
	  r1 += ((clust[ii].intpts[pm[j].p2[k]].xyz - 
		  clust[ii].intpts[pm[j].p2[l]].xyz).Mag());
	}
      }
      r2 /= (pm[j].np2 * (pm[j].np2 - 1.0) * 2);
    } else {
      r2 = clust[ii].nDet + 1.0;
    }

    pm[j].likely1 = 1.0 / r1 + 1.0 / r2;

    /* Also evaluate if number of hits in the cluster 
       is proper for the energy in the cluster, and 
       evaluate likelihood of split in that way 
       -- not implemented yet */
    pm[j].likely2 = 1.0; 
  }

  /* Sort for this order, so that in the following we try the 
     more likely splits first */
  for (j=0; j<pNum; j++) { order[j] = j; }
  for (j=0; j<pNum; j++) {
    for (k=0; k<pNum; k++) {
      if (pm[k].likely1*pm[k].likely2 > pm[j].likely1*pm[j].likely2) {
	memcpy (&tmppm, &pm[j], sizeof(struct PM));
	memcpy (&pm[j], &pm[k], sizeof(struct PM));
	memcpy (&pm[k], &tmppm, sizeof(struct PM));
	i1 = order[j];  order[j] = order[k];  order[k] = i1;
      }
    }
  }

  /* Go through the splits */

  for (j=0; j<pNum; j++) {
    /* Give up going on? */
    if (j > ctrl.splitClusters_maxTry[m]) {
      /* We have reached maximum number of tries, evaluate
	 where we are and return */
      clust[try1].valid = 0;  clust[try2].valid = 0;
      clust[best1].valid = 0;  clust[best2].valid = 0;

      if (nSplit > 0) {
	if (!clust[ii].tracked || (bestSplitFOM < clust[ii].FOM) ||
	    (bestSplitFOM < goodEnough)) {
	  /* Increase cluster index by two.  The new clusters are there 
	     and marked as valid. */
	  nClusters += 2;
	  clust[best1].valid = 1;
	  clust[best2].valid = 1;
	  clust[ii].valid = 0; /* Invalidate old cluster */
	  return (1);
	}
      }
      /* Here if we failed any attempt */
      return (6);
    }

    if (pm[j].np1 < MAXNUMSEG && pm[j].np2 < MAXNUMSEG) {
      /* Initialize trial clusters */
      clust[try1].Reset();  clust[try2].Reset();

      for (k=0; k<pm[j].np1; k++) {
	/* Fill trial cluster 1 */
	clust[try1].eSum += clust[ii].intpts[pm[j].p1[k]].eDet;
	clust[try1].intpts[clust[try1].nDet].xyz = 
	  clust[ii].intpts[pm[j].p1[k]].xyz;
	clust[try1].intpts[clust[try1].nDet].eDet = 
	  clust[ii].intpts[pm[j].p1[k]].eDet;
	clust[try1].intpts[clust[try1].nDet].timestamp = 
	  clust[ii].intpts[pm[j].p1[k]].timestamp;
	clust[try1].intpts[clust[try1].nDet].t0 = 
	  clust[ii].intpts[pm[j].p1[k]].t0;
	clust[try1].intpts[clust[try1].nDet].detNum = 
	  clust[ii].intpts[pm[j].p1[k]].detNum;
	clust[try1].intpts[clust[try1].nDet].eSumOrig = 
	  clust[ii].intpts[pm[j].p1[k]].eSumOrig;
	clust[try1].intpts[clust[try1].nDet].origPos = 
	  clust[ii].intpts[pm[j].p1[k]].origPos;
	clust[try1].intpts[clust[try1].nDet].shellPos = 
	  clust[ii].intpts[pm[j].p1[k]].shellPos;
	clust[try1].nDet++;
      }
      for (k=0; k<pm[j].np2; k++) {
	/* Fill trial cluster 2 */
	clust[try2].eSum += clust[ii].intpts[pm[j].p2[k]].eDet;
	clust[try2].intpts[clust[try2].nDet].xyz = 
	  clust[ii].intpts[pm[j].p2[k]].xyz;
	clust[try2].intpts[clust[try2].nDet].eDet = 
	  clust[ii].intpts[pm[j].p2[k]].eDet;
	clust[try2].intpts[clust[try2].nDet].timestamp = 
	  clust[ii].intpts[pm[j].p2[k]].timestamp;
	clust[try2].intpts[clust[try2].nDet].t0 = 
	  clust[ii].intpts[pm[j].p2[k]].t0;
	clust[try2].intpts[clust[try2].nDet].detNum = 
	  clust[ii].intpts[pm[j].p2[k]].detNum;
	clust[try2].intpts[clust[try2].nDet].eSumOrig = 
	  clust[ii].intpts[pm[j].p2[k]].eSumOrig;
	clust[try2].intpts[clust[try2].nDet].origPos = 
	  clust[ii].intpts[pm[j].p2[k]].origPos;
	clust[try2].intpts[clust[try2].nDet].shellPos = 
	  clust[ii].intpts[pm[j].p2[k]].shellPos;
	clust[try2].nDet++;
      }

      /* Sort for energy */
      clust[try1].sortSegEnergies(); 
      clust[try2].sortSegEnergies();

      /* Track try1 cluster */
      pm[j].FOM1 = -1;
      clust[try1].tracked = 0;
      if (clust[try1].nDet < MAXNUMSEG) {
	if (ctrl.trackOps[clust[try1].nDet] < 6 && 
	    ctrl.trackOps[clust[try1].nDet] >= 0 && 
	    ctrl.trackOps[clust[try1].nDet] != 2) {
	  st = doTrack(ctrl.trackOps[clust[try1].nDet], try1);
	  stat.trackError[ctrl.trackOps[clust[try1].nDet]][st]++;
	} else {
	  printf("TRACK: splitCluster: tracking option not known!?, option = %i\n, TRACK: Quitting...\n", ctrl.trackOps[clust[try1].nDet]);
	  exit(1);
	}
      }
      pm[j].FOM1 = clust[try1].FOM;
    
      /* Track try2 cluster, if it is reasonable to try */
      pm[j].FOM2 = -1;
      if (clust[try1].FOM < bestSplitFOM && clust[try1].FOM < clust[ii].FOM) {
	nTracksTried++;
	clust[try2].tracked = 0;
	if (clust[try2].nDet < MAXNUMSEG) {
	  if (ctrl.trackOps[clust[try2].nDet] < 6 && 
	      ctrl.trackOps[clust[try2].nDet] >= 0 && 
	      ctrl.trackOps[clust[try2].nDet] != 2) {
	    st = doTrack(ctrl.trackOps[clust[try2].nDet], try2);
	    stat.trackError[ctrl.trackOps[clust[try2].nDet]][st]++;
	  } else {
	    printf("TRACK: splitCluster: tracking option not known!?, option = %i\n, TRACK: Quitting...\n", ctrl.trackOps[clust[try2].nDet]);
	    exit(1);
	  }
	}
      } else { clust[try2].tracked = 0; }

      pm[j].FOM2 = clust[try2].FOM;

      /* At this point both trial clusters are tracked, we can evaluate 
	 if the split is better than the best we've seen so far. */
      if (clust[try1].tracked && clust[try2].tracked) {
	/* Combine FOM of split event */
	sumFOM = clust[try1].FOM + clust[try2].FOM;
	if (sumFOM < bestSplitFOM) {
	  nSplit++;
	  bestSplitFOM = sumFOM;
	  bestSplit = j;
	  
	  /* Transfer to best1 and best2 positions */
	  clust[best1] = clust[try1];  clust[best2] = clust[try2];
	  clust[best1].valid = 1;  clust[best2].valid = 1;

	  /* If this result is better than the splitCluster
	     kickout FOM, return with what we have now */
	  if (bestSplitFOM < ctrl.splitClusters_kickoutFOM[m]) {
	    nClusters += 2;
	    clust[ii].valid = 0;
	    return (2);
	  }

	  /* If this result is goodEnough, return at this point. */
	  if (bestSplitFOM < goodEnough) {
	    nClusters += 2;
	    clust[ii].valid = 0;
	    return (3);
	  } 
	} /* end of if (bestSplit so far) */	  
      }
    }
  }

  /* We didn't manage to split -- clean up and return */
  if (bestSplit < 0) {
    clust[try1].valid = 0;  clust[try2].valid = 0;
    clust[best1].valid = 0;  clust[best2].valid = 0;
    
    return (8);
  }

  if (nSplit > 0) {
    if ( !clust[ii].tracked || (bestSplitFOM < clust[ii].FOM ||
			       bestSplitFOM < goodEnough) ) {
      nClusters += 2;
      clust[ii].valid = 0;
      return (0); /* Success! */
    } else {
      clust[try1].valid = 0;  clust[try2].valid = 0;
      clust[best1].valid = 0;  clust[best2].valid = 0;
      return (8);
    }
  }

  return (0); /* We should never get here! */
	       
}
      
Int_t Track::combineCluster(Int_t focusC) {

  Int_t i, j, i1;
  Float_t p1, p2, p3, q1, q2, q3;
  Float_t dist, bestFOMsoFar = FLT_MAX, val, eSum;
  Int_t nCandidates = 0, nBase, best, try1, bestCombSoFar = -1;

  /* Setup the base of the trial cluster */
  best = nClusters; /* where the best combined cluster will end up */
  clust[best].valid = 0;
  try1 = best + 1; /* where the trial will be */

  clust[try1] = clust[focusC];
  nBase = clust[focusC].nDet;
  eSum = clust[focusC].eSum;

  /* Find the mean position of the cluster we want to combine from */
  p1 = 0; p2 = 0; p3 = 0;
  for (i=0; i<clust[focusC].nDet; i++) {
    p1 += clust[focusC].intpts[i].xyz.X();
    p2 += clust[focusC].intpts[i].xyz.Y();
    p3 += clust[focusC].intpts[i].xyz.Z();
  }
  p1 /= clust[focusC].nDet;  
  p2 /= clust[focusC].nDet;  
  p3 /= clust[focusC].nDet;  
  
  /* Go through and find clusters to include.
     Singles hits are always candidates, since FOM is zero.
     For others, we are looking for clusters with FOM over
     a specified threshold. */
  for (i=0; i<nClusters; i++) {
    if (i!=focusC && clust[i].valid && 
	(clust[i].FOM > ctrl.combineClusters_thresh || clust[i].nDet == 1)) {
      /* Find mean position of candidate cluster */
      q1 = 0;  q2 = 0;  q3 = 0;
      for (j=0; j<clust[i].nDet; j++) {
	q1 += clust[i].intpts[j].xyz.X();
	q2 += clust[i].intpts[j].xyz.Y();
	q3 += clust[i].intpts[j].xyz.Z();
      }
      q1 /= clust[i].nDet;  
      q2 /= clust[i].nDet;  
      q3 /= clust[i].nDet; 
      
      /* Find the distance */
      dist = (q1-p1)*(q1-p1) + (q2-p2)*(q2-p2) + (q3-p3)*(q3-p3);
      
      /* Is it a candidate? */
      if (dist < ctrl.combineMaxDist) {
	/* Add candidate combination cluster to base trial cluster */
	clust[try1].eSum = eSum;
	for (j=0; j<clust[i].nDet; j++) {
	  clust[try1].nDet++;
	  clust[try1].eSum += clust[i].intpts[j].eDet;
	  clust[try1].intpts[j+nBase].eDet = clust[i].intpts[j].eDet;
	  clust[try1].intpts[j+nBase].xyz = clust[i].intpts[j].xyz;
	  clust[try1].intpts[j+nBase].timestamp = clust[i].intpts[j].timestamp;
	  clust[try1].intpts[j+nBase].t0 = clust[i].intpts[j].t0;
	  clust[try1].intpts[j+nBase].detNum = clust[i].intpts[j].detNum;
	  clust[try1].intpts[j+nBase].origPos = clust[i].intpts[j].origPos;
	  clust[try1].intpts[j+nBase].eSumOrig = clust[i].intpts[j].eSumOrig;
	}
	
	/* Initialize trial cluster for tracking */
	clust[try1].Reset();

	/* Sort for energy */
	clust[try1].sortSegEnergies(); 
	
	Int_t st = -1;
	/* Track try1 cluster */
	if (clust[try1].nDet < MAXNUMSEG) {
	  if (ctrl.trackOps[clust[try1].nDet] < 6 && 
	      ctrl.trackOps[clust[try1].nDet] >= 0 && 
	      ctrl.trackOps[clust[try1].nDet] != 2) {
	    st = doTrack(ctrl.trackOps[clust[try1].nDet], try1);
	    stat.trackError[ctrl.trackOps[clust[try1].nDet]][st]++;
	  } else {
	    printf("TRACK: combineCluster: tracking option not known!?, option = %i\n, TRACK: Quitting...\n", ctrl.trackOps[clust[try1].nDet]);
	    exit(1);
	  }
	} else {
	  clust[try1].tracked = 0;
	}

	/* Store as best combination so far? */
	if (clust[try1].FOM < bestFOMsoFar) {
	  nCandidates++;
	  bestFOMsoFar = clust[try1].FOM;
	  bestCombSoFar = i;
	  clust[best] = clust[try1];
	}
      }
    }
  }

  /* We've tried all combinations as determined by parameter
     set.  Now decide whether to accept combination or not. */
  if (nCandidates > 0 && bestCombSoFar >= 0) {
    /* Find the uncombined 'combined' FOM */
    val = clust[focusC].FOM * clust[focusC].nDet;
    val += clust[bestCombSoFar].FOM * clust[bestCombSoFar].nDet;
    i1 = clust[focusC].nDet + clust[bestCombSoFar].nDet;
    val /= i1;
    
    /* Are we better off with the combination? */
    if (clust[best].FOM < val) {
      /* Validate the new combined cluster */
      clust[best].valid = 1; 
      nClusters++;
      
      /* Invalidate the two clusters we combined */
      clust[focusC].valid = 0;  clust[bestCombSoFar].valid = 0;

      /* Return success */
      return (0);
      
    } else {
      return (1);   
    }
  } else {
    return (2);
  }

  /* Done, we shouldn't really get here */
  return (0);
}	      
	      

Int_t Track::matchMaker(Int_t ii, Int_t jj) {
  Float_t d1;
  Double_t dist;
  Int_t try1, j;

  /* Find distance between points, and return if they're too far apart */
  dist = (clust[ii].intpts[0].xyz - clust[ii].intpts[0].xyz).Mag2();
  if (dist > ctrl.matchMakerMaxDist) { return (1); }

  /* Two single hits are close enough to consider as a trial cluster */
  try1 = nClusters;
  clust[try1].nDet = 2;
  clust[try1].eSum = clust[ii].eSum + clust[jj].eSum;
  clust[try1].intpts[0].detNum = clust[ii].intpts[0].detNum;
  clust[try1].intpts[0].xyz = clust[ii].intpts[0].xyz;
  clust[try1].intpts[0].timestamp = clust[ii].intpts[0].timestamp;
  clust[try1].intpts[0].t0 = clust[ii].intpts[0].t0;
  clust[try1].intpts[0].shellPos = clust[ii].intpts[0].shellPos;
  clust[try1].intpts[0].eSumOrig = clust[ii].intpts[0].eSumOrig;
  clust[try1].intpts[0].origPos = clust[ii].intpts[0].origPos;
 
  clust[try1].intpts[1].detNum = clust[jj].intpts[0].detNum;
  clust[try1].intpts[1].xyz = clust[jj].intpts[0].xyz;
  clust[try1].intpts[1].timestamp = clust[jj].intpts[0].timestamp;
  clust[try1].intpts[1].t0 = clust[jj].intpts[0].t0;
  clust[try1].intpts[1].shellPos = clust[jj].intpts[0].shellPos;
  clust[try1].intpts[1].eSumOrig = clust[jj].intpts[0].eSumOrig;
  clust[try1].intpts[1].origPos = clust[jj].intpts[0].origPos;
  
  /* Initialize trial cluster for tracking */
  clust[try1].Reset();

  /* Sort for energy */
  clust[try1].sortSegEnergies(); 
  
  /* Track try1 cluster */
  Int_t st = -1;
  if (clust[try1].nDet < MAXNUMSEG) {
    if (ctrl.trackOps[clust[try1].nDet] < 6 && 
	ctrl.trackOps[clust[try1].nDet] >= 0 && 
	ctrl.trackOps[clust[try1].nDet] != 2) {
      st = doTrack(ctrl.trackOps[clust[try1].nDet], try1);
      stat.trackError[ctrl.trackOps[clust[try1].nDet]][st]++;
    } else {
      printf("TRACK: matchMaker: tracking option not known!?, option = %i\n, TRACK: Quitting...\n", ctrl.trackOps[clust[try1].nDet]);
      exit(1);
    }
  } else {
    clust[try1].tracked = 0;
  }
  
  /* If OK, promote to valid cluster */
  if (clust[try1].FOM < ctrl.matchMaker_kickoutFOM) {
    /* Validate the new cluster */
    clust[try1].valid = 1;
    nClusters++;
    
    /* Invalidate the old single hits */
    clust[ii].valid = 0;  clust[jj].valid = 0;
    
    return (0);
  } else {
    clust[try1].valid = 0;
    return (2);
  }
  
  /* Done, but shouldn't really get here */
  return (0);
  
}

Int_t Track::findClusters() {
  /* This is actually just assigning cluster numbers */
  Int_t i, j, k, i1, i2;
  
  for (k=0; k<shell.nHit; k++) { shell.clusterNum[k] = -1; }

  shell.numClusters = 0;
  
  /* Hunt for multihit clusters */
  for (i=0; i<shell.nHit; i++) {
    for (j=0; j<shell.nHit; j++) {
      if (i != j) {
	/* Single out hits that are close */
	if (shell.relAngle(i,j) < ctrl.alpha[shell.nHit]) {
	  /* We have two interactions close in angle */
	  if (shell.clusterNum[i] < 0 && shell.clusterNum[j] < 0) {
	    /* Neither has been clustered previously, make a new cluster. */
	    shell.clusterNum[i] = shell.numClusters;
	    shell.clusterNum[j] = shell.numClusters;
	    shell.numClusters++;
	  } else if (shell.clusterNum[i] >= 0 && shell.clusterNum[j] < 0) {
	    /* [j] should join [i]'s cluster */
	    shell.clusterNum[j] = shell.clusterNum[i];
	  } else if (shell.clusterNum[i] < 0  && shell.clusterNum[j] >= 0) {
	    /* [i] should join [j]'s cluster */
	    shell.clusterNum[i] = shell.clusterNum[j];
	  } else if (shell.clusterNum[i] >= 0 && shell.clusterNum[j] >= 0 &&
		     shell.clusterNum[i] != shell.clusterNum[j]) {
	    /* Oops, colliding clusters - put them together! 
	       -. i1 is common cluster index, i2 will move to i1 */
	    if (shell.clusterNum[i] > shell.clusterNum[j]) {
	      i1 = shell.clusterNum[j];  i2 = shell.clusterNum[i];
	    } else {
	      i1 = shell.clusterNum[i];  i2 = shell.clusterNum[j];
	    }
	    
	    for (k=0; k<shell.nHit; k++) {
	      if (shell.clusterNum[k] == i2) {
		shell.clusterNum[k] = i1;
	      }
	    }
	    
	    /* Bump everyone above i2 down by one index */
	    for (k=0; k<shell.nHit; k++) {
	      if (shell.clusterNum[k] > i2) {
		shell.clusterNum[k]--;
	      }
	    }
	    
	    shell.numClusters--;
	  }
	}
      }
    }
  }

  /* Assign single hit clusters to the rest that never
     found any partner. */
  for (i=0; i<shell.nHit; i++) {
    if (shell.clusterNum[i] < 0) {
      shell.clusterNum[i] = shell.numClusters;
      shell.numClusters++;
    }
  }

  /* Done! */
  return (0);
 
}

Int_t Track::transferClusters() {
  Int_t i, j, clustNum, insertPos;

  /* Initialize all cluster counters 
     and sum energy, valid flags, etc. */
  for (i=0; i<MAXCLUSTERHITS; i++) {
    clust[i].Reset();
  }

  /* Loop over all shell hits and use the cluster 
     number to fill up the cluster array. */
  for (i=0; i<shell.nHit; i++) {
    
    /* Transfer to proper cluster array (index: clustNum)
       and the right position in this cluster (index: insertPos)
       -. the cluster may already have some members */
    clustNum = shell.clusterNum[i];
    insertPos = clust[clustNum].nDet;

    if (insertPos == 0) { clust[clustNum].valid = 1; }
    /* Check we don't have too many points... */
    if (insertPos >= MAXNUMDET) {
      printf("TRACK: Too many hits...\n");
      return (1);
    }
    
    clust[clustNum].trackNum = clustNum;
    clust[clustNum].intpts[insertPos].shellPos = i; /* Needed later */
    clust[clustNum].intpts[insertPos].xyz = shell.xyz[i];
    clust[clustNum].intpts[insertPos].eDet = shell.eDet[i];
    clust[clustNum].intpts[insertPos].timestamp = shell.timestamp[i];
    clust[clustNum].intpts[insertPos].t0 = shell.t0[i];
    clust[clustNum].intpts[insertPos].detNum = shell.detNum[i];

    /* For now, add up track energies to fill eSum 
       -- will be reassigned later. */
    clust[clustNum].eSum += shell.eDet[i];
    
    clust[clustNum].nDet++;
    if (clust[clustNum].nDet >= MAXNUMDET) {
      return (TOOMANYHITS);
    }
  }

  /* Done! */
  return (0);
}

void Track::PrintValidClusters() {
  printf("Total clusters: %i\n", nClusters);
  for (Int_t i=0; i<nClusters; i++) {
    if (clust[i].valid) {
      printf(" Valid cluster %i: \n", i);
      clust[i].PrintCluster();
    }
  }
}
