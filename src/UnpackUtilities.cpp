#include "UnpackUtilities.h"

Int_t OpenInputFile(FILE** inf, controlVariables* ctrl, TString runNumber) {
  
  if (ctrl->fileType != "f" && ctrl->fileType != "f1" && ctrl->fileType != "f2")  {

    if (!ctrl->analyze2AND3) {
      if (!ctrl->directory.EndsWith("/")) { ctrl->directory = ctrl->directory + "/"; }
      
      if (ctrl->fileType == "g") {
	ctrl->fileName = ctrl->directory + "Run" + runNumber + "/Global.dat";
      } else if (ctrl->fileType == "gr") {
	ctrl->fileName = ctrl->directory + "Run" + runNumber + "/GlobalRaw.dat";
      } else if (ctrl->fileType == "m") {
	ctrl->fileName = ctrl->directory + "Run" + runNumber + "/Merged.dat";
      } else if (ctrl->fileType == "m2") {
	ctrl->fileName = ctrl->directory + "Run" + runNumber + "/Merged.Mode2.dat";
      } else {
	cerr << "WHAT???" << endl;
	return 0;
      }
      
      if (ctrl->compressedFile) {

	if (ctrl->noHFC) {
	  if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) {
	    ctrl->fileName = ctrl->fileName + ".gz";
	  }
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    fclose(*inf);
	    ctrl->fileName = "zcat " + ctrl->fileName;	  
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	  
	} else if (!ctrl->noHFC) {

	  if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) { 
	    ctrl->fileName = ctrl->fileName + ".gz";
	  }
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    fclose(*inf);
	    ctrl->fileName = "./GEB_HFC -z -p " + ctrl->fileName;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	}
	
      } else if (ctrl->compressedFileB) {

	if (ctrl->noHFC) {
	  if (!ctrl->fileName.EndsWith(".bz2")) { 
	    ctrl->fileName = ctrl->fileName + ".bz2";
	  }	 
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    fclose(*inf);
	    ctrl->fileName = "bzcat " + ctrl->fileName;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	  
	} else if (!ctrl->noHFC) {

	  if (!ctrl->fileName.EndsWith(".bz2")) { 
	    ctrl->fileName = ctrl->fileName + ".bz2";
	  }
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    ctrl->fileName = "./GEB_HFC -bz -p " + ctrl->fileName;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	}
	
      } else if (ctrl->noHFC) {

	*inf = fopen(ctrl->fileName.Data(), "r");

      } else {

	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  ctrl->fileName = "./GEB_HFC -p " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      }

    } else if (ctrl->analyze2AND3) {

      if (!ctrl->directory.EndsWith("/")) { ctrl->directory = ctrl->directory + "/"; }
      
      ctrl->fileName = ctrl->directory + "Run" + runNumber + "/GlobalRaw.dat"; /* Mode3 data */
      ctrl->fileName2 = ctrl->directory + "Run" + runNumber + "/Global.dat"; /* Mode2 data */
      
      ctrl->fileName = "./MergeArbFiles -f1 " + ctrl->fileName + " -f2 " + ctrl->fileName2 + " -fOut pipe";
      
      cout << ctrl->fileName.Data() << endl;
      
      *inf = popen(ctrl->fileName.Data(), "r");
    }
    
  } else if (ctrl->fileType == "f" || ctrl->fileType == "f1" || ctrl->fileType == "f2") {

    if (ctrl->compressedFile) {

      if (!ctrl->analyze2AND3) {
	if (ctrl->noHFC) {
	  if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) {
	    ctrl->fileName = ctrl->fileName + ".gz";
	  }
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    fclose(*inf);
	    ctrl->fileName = "zcat " + ctrl->fileName;
	    *inf = NULL;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	  
	} else if (!ctrl->noHFC) {
	  
	  if (!ctrl->fileName.EndsWith(".gz") && !ctrl->fileName.EndsWith(".gzip")) { 
	    ctrl->fileName = ctrl->fileName + ".gz";
	  }
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    fclose(*inf);
	    ctrl->fileName = "./GEB_HFC -z -p " + ctrl->fileName;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	}
      } else {
	cout << "Apologies -- multiple file analysis at present is not possible with compressed files. " << endl;
	return(-1);
      }
	
    } else if (ctrl->compressedFileB) {

      if (!ctrl->analyze2AND3) {
	if (ctrl->noHFC) {
	  if (!ctrl->fileName.EndsWith(".bz2")) { 
	    ctrl->fileName = ctrl->fileName + ".bz2";
	  }	 
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    fclose(*inf);
	    ctrl->fileName = "bzcat " + ctrl->fileName;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	  
	} else if (!ctrl->noHFC) {
	  
	  if (!ctrl->fileName.EndsWith(".bz2")) { 
	    ctrl->fileName = ctrl->fileName + ".bz2";
	  }
	  *inf = fopen(ctrl->fileName.Data(), "r");
	  if (*inf) {
	    ctrl->fileName = "./GEB_HFC -bz -p " + ctrl->fileName;
	    *inf = popen(ctrl->fileName.Data(), "r");
	  }
	}
      } else {
	cout << "Apologies -- multiple file analysis at present is not possible with compressed files. " << endl;
	return(-1);
      }	
    } else if (ctrl->noHFC) {

      if (!ctrl->analyze2AND3) {
	*inf = fopen(ctrl->fileName.Data(), "r");
      }
	
    } else {

      if (!ctrl->analyze2AND3) {
	*inf = fopen(ctrl->fileName.Data(), "r");
	if (*inf) {
	  ctrl->fileName = "./GEB_HFC -p " + ctrl->fileName;
	  *inf = popen(ctrl->fileName.Data(), "r");
	}
      }
      
    }
  
    if (ctrl->analyze2AND3) {

      if (ctrl->compressedFile) {
	cout << "Apologies -- multiple file analysis at present is not possible with compressed files. " << endl;
	return(-1);
      } else if (ctrl->compressedFileB) {
	cout << "Apologies -- multiple file analysis at present is not possible with compressed files. " << endl;
	return(-1);
      } else {
	ctrl->fileName = "./MergeArbFiles -f1 " + ctrl->fileName + " -f2 " + ctrl->fileName2 + " -fOut pipe";	
	cout << ctrl->fileName.Data() << endl;
	*inf = popen(ctrl->fileName.Data(), "r");
      }
    }

  }
  
  if (!*inf) {

    printf("Cannot open: %s \n", ctrl->fileName.Data());
    return(2);

  } else {
    
    printf("Opened: %s \n", ctrl->fileName.Data());  
    
    if (ctrl->fileType != "f" && ctrl->fileType != "f1" && ctrl->fileType != "f2") {
      ctrl->outfileName = (ctrl->directory + "Run" + runNumber + "/Run" + runNumber +
			   ctrl->outputSuffix + ".root");
    } else {
      
      if (ctrl->outfileName == "") {
	ctrl->outfileName = ("./ROOTFiles/test.root");
      } else { /* Do nothing, we have a filename. */ } 

    }
    
    return(0);
  }
  
  return(0);
}

int ProcessEvent(Float_t currTS, controlVariables* ctrl, counterVariables* cnt) {

  int badCrystal = 0;

  /* Go ahead and write out the event, reset for the next one... */
  if (ctrl->superPulse) { /* This is a superpulse analysis -- quite different. */
    gret->checkSPIntegrity();
    gret->sp.MakeSuperPulses();
  } else if (ctrl->xtalkAnalysis) { /* Cross-talk analysis, also different... */
    Int_t multTemp = gret->g3out.crystalMult();
    //FillFoldHist(histos, gVar, multTemp, ctrl->xtLowE, ctrl->xtHighE);
  }

  if (gret->g3Temp.size() > 0) { gret->analyzeMode3(ctrl); }

  if (gret->g2out.crystalMult() > 0) {
    if (ctrl->doTRACK) {
      Int_t error = gret->fillShell2Track();
      if (error != 8) {
	gret->track.findTargetPos();
	Int_t trackStatus;
	trackStatus = gret->track.trackEvent();
	gret->fillMode1(trackStatus);
      }
    }
  }

  if (badCrystal >= 0) {
    if (ctrl->analyze2AND3) {
      if (gret->g3out.crystalMult() > 0) {
	for (UInt_t i=0; i<gret->g2out.crystalMult(); i++) {
	  gret->g2out.xtals[i].waveAll.resize(40*gret->g3out.xtals[0].traceLength, 0);
	  for (UInt_t um=0; um<gret->g3out.xtals.size(); um++) {
	    if (gret->g3out.xtals[um].module == gret->g2out.xtals[i].crystalID) { /* CHECK THIS MODULE = CRYSTALID */
	      for (Int_t j=0; j<40; j++) {
		UInt_t waveSize = gret->g3out.xtals[um].chn[j].wf.raw.size();
		if (waveSize > 0) {
		  for (UInt_t uk=0; uk<waveSize; uk++) {
		    gret->g2out.xtals[i].waveAll[j*waveSize + uk] = gret->g3out.xtals[um].chn[j].wf.raw[uk];
		  }
		}
	      }
	    }
	  }
	}
	for (UInt_t i=0; i<gret->g2out.crystalMult(); i++) {
	  if (gret->g2out.xtals[i].waveAll.size() <= 0) {
	    cnt->nMode2NoTraces[gret->g2out.xtals[i].crystalNum-1]++;
	  }
	}
	for (UInt_t i=0; i<gret->g3out.crystalMult(); i++) {
	  Int_t INCLUDED = 0;
	  for (UInt_t j=0; j<gret->g2out.crystalMult(); j++) {
	    if (gret->g3out.xtals[i].crystalNum == gret->g2out.xtals[j].crystalNum) {
	      INCLUDED = 1;
	    }
	  }
	  if (!INCLUDED) { 
	    cnt->nMode3NoMode2[gret->g3out.xtals[i].crystalNum-1]++;
	    g2CrystalEvent g2;
	    g2.crystalNum = gret->g3out.xtals[i].crystalNum;
	    g2.quadNum = gret->g3out.xtals[i].quadNum;
	    g2.timestamp = gret->g3out.xtals[i].LEDHigh();
	    g2.t0 = -1.0;
	    g2.chiSq = 0.0;
	    if (gret->g3out.xtals[i].crystalBuild() == 40) {
	      g2.error = 0xBAD;
	    } else {
	      g2.error = 0xBADC;
	    }
	    g2.cc = gret->g3out.xtals[i].cc4();
	    g2.cc1 = -1.0;
	    g2.cc2 = -1.0;
	    g2.cc3 = -1.0;
	    g2.cc4 = -1.0;
	    g2.segSum = -1.0;
	    g2.doppler = 0.0;  g2.dopplerSeg = 0.0;  g2.dopplerCrystal = 0.0;
	    
	    g2.waveAll.resize(40*gret->g3out.xtals[0].traceLength, 0);
	    for (Int_t j=0; j<40; j++) {
	      UInt_t waveSize = gret->g3out.xtals[i].chn[j].wf.raw.size();
	      if (waveSize > 0) {
		for (UInt_t uk=0; uk<waveSize; uk++) {
		  g2.waveAll[j*waveSize + uk] = gret->g3out.xtals[i].chn[j].wf.raw[uk];
		}
	      }
	    }
	    gret->g2out.xtals.push_back(g2);
	  }
	}
	gret->g3out.Clear();
      }   
    } /* End of analyze2AND3 */
  }
  
  /* Do particle-gamma type things -- i.e. correct Doppler
     for particle trajectories, etc. */
  if (!ctrl->dopplerSimple) {
    if (gret->g2out.crystalMult() > 0) {
      for (UInt_t ui=0; ui<gret->g2out.crystalMult(); ui++) {
#ifdef WITH_S800
	if ((cnt->event & 0x10) || (cnt->event & 0x100)) {
	  
	  /* First scale the radius... this trick comes from Dirk.  We still
	     don't really know why we need to do this, but we do. It makes
	     things better.  Probably something funny in the decomp bases.  */
	  TVector3 xyz = gret->g2out.xtals[ui].maxIntPtXYZ();
	  //  Double_t pR = xyz.XYvector().Mod();
	  // Double_t pTheta = TMath::ATan(xyz.Y()/xyz.X());
	  //Float_t xPrime = gret->var.radiusCor[gret->g2out.xtals[ui].crystalNum]*pR*TMath::Cos(pTheta);
	  //Float_t yPrime = gret->var.radiusCor[gret->g2out.xtals[ui].crystalNum]*pR*TMath::Sin(pTheta);
	  
	  /* This is a check against flipping the sign of the x/y coordinates...*/
	  //if (xyz.X() < 0 && xPrime > 0) { xPrime = -xPrime; }
	  //if (xyz.X() > 0 && xPrime < 0) { xPrime = -xPrime; }
	  //if (xyz.Y() < 0 && yPrime > 0) { yPrime = -yPrime; }
	  //if (xyz.Y() > 0 && yPrime < 0) { yPrime = -yPrime; }
	  //xyz.SetX(xPrime);  xyz.SetY(yPrime);
	  
	  /* And translate to world coordinates... */
	  TVector3 xyzL = gret->rot.crys2Lab(gret->g2out.xtals[ui].crystalID, xyz);
	  //printf("Here %f,%f,%f  ---> ", xyzL.X(), xyzL.Y(), xyzL.Z());
	  xyzL -= gret->var.targetXYZ;      
	  //printf("%f %f %f\n", xyzL.X(), xyzL.Y(), xyzL.Z());
	  /* Calculate vector from target to interaction point, 
	   including shifts in position (in cm). */



	  gret->g2out.xtals[ui].doppler = s800->getDoppler(xyzL, gret->var.beta, &gret->var);
	  //printf("Doppler %f\n", gret->g2out.xtals[ui].doppler);
	  if (gret->g2out.xtals[ui].doppler == 0) { 
	    gret->g2out.xtals[ui].doppler = gret->getDopplerSimple(gret->g2out.xtals[ui].maxIntPtXYZLab(), gret->var.beta); 
	  }
	}
#else /* WITH_S800 */
#ifdef WITH_CHICO
	TVector3 xyzL = gret->g2out.xtals[ui].maxIntPtXYZLab();    

	/*printf("xyzL.Y(),Z() before %f and %f ",xyzL.Y(),xyzL.Z());*/

      	if (xyzL.Y()<0) xyzL.SetY(xyzL.Y() - GTPosOffsetY1) ;
	else            xyzL.SetY(xyzL.Y() - GTPosOffsetY2) ;

        xyzL.SetZ(xyzL.Z() - GTPosOffsetZ) ;

	/*printf("and after %f and %f \n",xyzL.Y(),xyzL.Z());*/

	Float_t chicoTheta = chico->particle.fThetaL/(TMath::Pi()) *180.;
	if (chico->gotParticle && (chico->idParticle > 0)) {
	  Float_t gTheta, gPhi;
	  gTheta = xyzL.Theta();
	  if (xyzL.Phi() < 0) { gPhi = xyzL.Phi() + 2*TMath::Pi(); } 
	  else { gPhi = xyzL.Phi(); }

	  /*
	    printf("pgCosL, chico->particle.fThetaL %f chico->particle.fPhiL %f gTheta %f gPhi %f\n",chico->particle.fThetaL, chico->particle.fPhiL, gTheta, gPhi);
	  */
	  /* chico->idParticle defined as:
	     1 : Target = Right, Projectile = Left
	     2 : Target = Left , Projectile = Right
	  */

	  if (chico->idParticle == 1) {
	    chico->particle.pgCosP = chico->calcCos(chico->particle.fThetaL, chico->particle.fPhiL, gTheta, gPhi);
	    chico->particle.pgCosT = chico->calcCos(chico->particle.fThetaR, chico->particle.fPhiR, gTheta, gPhi);
	  } else if (chico->idParticle == 2) {
	    chico->particle.pgCosT = chico->calcCos(chico->particle.fThetaL, chico->particle.fPhiL, gTheta, gPhi);
	    chico->particle.pgCosP = chico->calcCos(chico->particle.fThetaR, chico->particle.fPhiR, gTheta, gPhi);
	  } else {
	    printf("Unexpected chico->idParticle = %d\n", chico->idParticle);
	  }

	  /* not applicable for deep inelastic
	     if (chicoTheta >=96. && chicoTheta <=180.) {
	     lbeta = chico->betaB[(Int_t)chicoTheta];
	     } else {
	     lbeta = chico->betaP0 + chico->betaP1*chico->particle.fThetaL;
	     lbeta += chico->betaP2*pow(chico->particle.fThetaL,2);
	     lbeta += chico->betaP3*pow(chico->particle.fThetaL,3);
	     lbeta += chico->betaP4*pow(chico->particle.fThetaL,4);
	     lbeta += chico->betaP5*pow(chico->particle.fThetaL,5);
	     }

	     if (chicoTheta >=96. && chicoTheta <=180.) {
	     lbeta = chico->betaB[(Int_t)chicoTheta];
	     } else {
	     lbeta = chico->betaT0 + chico->betaT1*chico->particle.fThetaL;
	     lbeta += chico->betaT2*pow(chico->particle.fThetaL,2);
	     lbeta += chico->betaT3*pow(chico->particle.fThetaL,3);
	     lbeta += chico->betaT4*pow(chico->particle.fThetaL,4);
	     lbeta += chico->betaT5*pow(chico->particle.fThetaL,5);
	     }
	  */

	  
	  Float_t lbeta, lcosTheta, lgamma, ldoppler;

	  lbeta = chico->particle.betaT;
	  lcosTheta = chico->particle.pgCosT;

	  lgamma = 1/sqrt(1.-lbeta*lbeta);	  
	  ldoppler = (lgamma*(1-lbeta*lcosTheta));
          
	  gret->g2out.xtals[ui].dopplerT = ldoppler;

	  lbeta = chico->particle.betaP;
	  lcosTheta = chico->particle.pgCosP;

	  lgamma = 1/sqrt(1.-lbeta*lbeta);	  
	  ldoppler = (lgamma*(1-lbeta*lcosTheta));
          
	  gret->g2out.xtals[ui].doppler = ldoppler;


	  /*
	    printf("chico->particle.t %llu chico->particle.LEDts %llu chicoTheta %f chico->particle.fThetaL %f lbeta %f chico->particle.pgCosL %f dopplerL %f \n",chico->particle.t, chico->particle.LEDts, chicoTheta, chico->particle.fThetaL, lbeta, chico->particle.pgCosL, dopplerL);
	  */

	} else { /* No CHICO particle info -- simple doppler */
	  gret->g2out.xtals[ui].dopplerT = 0.0;
	  gret->g2out.xtals[ui].doppler = gret->getDopplerSimple(gret->g2out.xtals[ui].maxIntPtXYZLab(), gret->var.beta);
	}
#else
	/* Nothing to do, no particle information available. */
#endif
#endif
      }
    } /* End of Mode2 data? */
    if (gret->g1out.gammaMult() > 0) {
      for (UInt_t ui=0; ui<gret->g1out.gammaMult(); ui++) {
#ifdef WITH_S800
	gret->g1out.gammas[ui].doppler = s800->getDoppler(gret->g1out.gammas[ui].xyzLab1, gret->var.beta, &gret->var);
	if (gret->g1out.gammas[ui].doppler == 0) {
	  gret->g1out.gammas[ui].doppler = gret->getDopplerSimple(gret->g1out.gammas[ui].xyzLab1, gret->var.beta);
	}
#else
	/* Nothing to do, no particle information available. */
#endif
      }
    }
  } /* End of !ctrl->dopplerSimple */    

  /* Write the histograms/tree */
  if (currTS > 0 && badCrystal >= 0) {
    if (ctrl->withHISTOS && ctrl->calibration && !ctrl->xtalkAnalysis) { 
      gret->fillHistos(1);
    } else if (ctrl->withHISTOS) {
      gret->fillHistos(2);
    }
    if (ctrl->withTREE) { 
      teb->Fill(); cnt->treeWrites++;
#ifdef WITH_PWALL
      /* Reset Phoswall */
      phosWall->Reset();
      phosWall->aux.Reset();
#endif
#ifdef WITH_CHICO
      chico->Reset();
#endif
#ifdef WITH_LENDA
      for (UInt_t ui=0; ui<ddasEv->getData().size(); ui++) {
	delete ddasEv->getData()[ui];
      }
      ddasEv->getData().clear();
      lendaEv->Clear();
#endif
#ifdef WITH_GOD
      goddess->Clear();
#endif
    }
  } 
  
  return (badCrystal);
}

void ResetEvent(controlVariables* ctrl, counterVariables* cnt) {
    
  /* Clear event structures. */
  if ((cnt->getEventBit(RAW)) || (cnt->getEventBit(DECOMP)) ||
      (cnt->getEventBit(BANK88)) || (cnt->getEventBit(RAWHISTORY)) ||
      (cnt->getEventBit(TRACK)) || (cnt->getEventBit(GRETSCALER))) {
    gret->Reset();
  }
  /* Reset temporary crystal event structures. */
#ifdef WITH_S800
  if ((cnt->event & 0x10) || (cnt->event & 0x100)) { s800->Reset(); }
#endif
#ifdef WITH_CHICO
  if (cnt->event & 0x800) { chico->Reset(); }
#endif
  
  cnt->event = 0x0000;
}



