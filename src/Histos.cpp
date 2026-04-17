#include "Histos.h"

ClassImp(Histos);

/****************************************************/

void Histos::InitializeHist(controlVariables* ctrl) {

  Int_t i=0; 

  /* If this is a calibration run, for energy or for cross-talk, 
     we use different histograms. */
  if (!ctrl->calibration && !ctrl->xtalkAnalysis) {

    /* Individual crystal spectra. */
    for (i=0; i<MAXCRYSTALS; i++) {
      char str[300];
      sprintf(str, "cc1_xtal%d", i+1);
      cc1cal[i] = new TH1F(str, str, 2000, 0, 2000);
      
      sprintf(str, "cc2_xtal%d", i+1);
      cc2cal[i] = new TH1F(str, str, 2000, 0, 2000);
      
      sprintf(str, "cc3_xtal%d", i+1);
      cc3cal[i] = new TH1F(str, str, 2000, 0, 2000);
      
      sprintf(str, "cc4_xtal%d", i+1);
      cc4cal[i] = new TH1F(str, str, 2000, 0, 2000);
    }
    
    /* Total projections. */
    for (i=0; i<4; i++) {
      char str[300];
      sprintf(str, "cc%dtotal", i+1);
      ccTot[i] = new TH1F(str, str, 3000, 0, 3000);
      sprintf(str, "cc%ddoptotal", i+1);
      ccDop[i] = new TH1F(str, str, 3000, 0, 3000);
    }
    
    /* Summary histograms...i.e. crystal # vs. something. */
    ccSummary = new TH2F("ccsummary", "ccsummary", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 3000, 0, 3000);
    ccDopSummary = new TH2F("ccdopsummary", "ccdopsummary", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 
  			    3000, 0, 3000);
    
    /* Multiplicity, hit patterns, etc. */
    mult = new TH1F("mult", "mult", MAXCRYSTALS + 1, 0, MAXCRYSTALS + 1);
    ccHit = new TH1F("cchitpat", "cchitpat", MAXCRYSTALS, 1, MAXCRYSTALS + 1);
    
    /* Gamma-gamma matrices. */
    ggMat = new TH2F("ggmat", "ggmat", 2500, 0, 10000, 2500, 0, 10000);
    ggMatDop = new TH2F("ggmatdop", "ggmatdop", 3000, 0, 3000, 3000, 0, 3000);
    
    
    /*********** Mode3-specific GRETINA histograms. ***********/
    for (i=0; i<MAXCRYSTALS; i++) {
      char str[300];
      sprintf(str, "basevsE_xtal%d", i+1);
      baseVsE[i] = new TH2F(str, str, 1500, 0, 1500, 1000, -500, 500);
      
      sprintf(str, "bankTSdiff_xtal%d", i+1);
      bankTSdiff[i] = new TH1I(str, str, 20, -10, 10);
      
      sprintf(str, "bankTSlow_xtal%d", i+1);
      bankTSlow[i] = new TH1I(str, str, 40, 0, 40);
      
      sprintf(str, "bankTShigh_xtal%d", i+1);
      bankTShigh[i] = new TH1I(str, str, 40, 0, 40);
    }
    
    /*********** Mode2-specific GRETINA histograms. ***********/
    xyProj = new TH2F("xyproj", "xyproj", 1000, -40, 40, 1000, -40, 40);
    xzProj = new TH2F("xzproj", "xzproj", 1000, -40, 40, 1000, -40, 40);
    yzProj = new TH2F("yzproj", "yzproj", 1000, -40, 40, 1000, -40, 40);
    xyz = new TH3F("xyz", "xyz", 500, -40, 40, 500, -40, 40, 500, -40, 40);
    thetaPhi = new TH2F("thetaphi", "thetaphi", 360, 0, 180, 360, 0, 360);
    
    thetaE = new TH2F("thetaE", "thetaE", 180, 0, 180, 3000, 0, 3000);
    thetaEDop = new TH2F("thetaEdop", "thetaEdop", 180, 0, 180, 3000, 0, 3000);
    
    chi2VsXtal = new TH2F("chi2vsxtal", "chi2vsxtal", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 1000, 0, 100);
    
    for (i=0; i<MAXCRYSTALS; i++) {
      char str[300];
      sprintf(str, "phivsE_xtal%d", i+1);
      phiVsE[i] = new TH2F(str, str, 1500, 0, 3000, 1500, 0, 2*TMath::Pi());
    }
    
    /***************** Total event histograms *****************/
    deltaT = new TH1F("dt", "dt", 3000, -3000, 3000);
    
    /***************** Spare spectrum *************************/
    test = new TH1F ("test", "test", 3000, 0, 3000);
    
  } else if (ctrl->calibration && !ctrl->xtalkAnalysis) {
    for (i=0; i<MAXCHANNELS; i++) {
      char str[300];
      sprintf(str, "eraw%i", i);
      eraw[i] = new TH1F(str, str, 30000, 0, 30000);
      eraw[i]->SetXTitle("channels");
      eraw[i]->SetYTitle("counts");
      
      sprintf(str, "ecal%i", i);
      ecal[i] = new TH1F(str, str, 5000, 0, 5000);
      ecal[i]->SetXTitle("energy (keV)");
      ecal[i]->SetYTitle("counts");

      sprintf(str, "smoke%i", i);
      if (i%40 == 9) { /* 10 MeV CC */
  	smoke[(4*(Int_t)(i/40))] = new TH2F(str, str, 2000, 1000, 1500, 
  					    4500, -500, 1000);
      } else if (i%40 == 19) { /* 30 MeV CC */
  	smoke[(4*(Int_t)(i/40)) + 1] = new TH2F(str, str, 200, 1000, 2000,
  						1000, -1000, 1000);
      } else if (i%40 == 29) { /* 2.5 MeV CC */
  	smoke[(4*(Int_t)(i/40)) + 2] = new TH2F(str, str, 200, 1000, 2000, 
  						1000, -5000, 5000);
      } else if (i%40 == 39) { /* 5 MeV CC */
  	smoke[(4*(Int_t)(i/40)) + 3] = new TH2F(str, str, 2000, 1000, 2000,
  						3000, -500, 2500);
      }
    }
     
    for (i=0; i<MAXCRYSTALS; i++) {
      char str[300];
      sprintf(str, "xtalSummary%i", i);
      xtalSummary[i] = new TH2F(str, str, 40, 0, 40, 5000, 0, 5000);
      xtalSummary[i]->SetXTitle("segmentID");
      xtalSummary[i]->SetYTitle("energy (keV)");
    }

    cc1All = new TH1F("cc1Total", "cc1Total", 5000, 0, 5000);
    cc2All = new TH1F("cc2Total", "cc2Total", 5000, 0, 5000);
    cc3All = new TH1F("cc3Total", "cc3Total", 5000, 0, 5000);
    cc4All = new TH1F("cc4Total", "cc4Total", 5000, 0, 5000);
    cc1Summary = new TH2F("cc1Summary", "cc1Summary", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 
  		      5000, 0, 5000);
    cc2Summary = new TH2F("cc2Summary", "cc2Summary", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 
  		      5000, 0, 5000);
    cc3Summary = new TH2F("cc3Summary", "cc3Summary", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 
  		      5000, 0, 5000);
    cc4Summary = new TH2F("cc4Summary", "cc4Summary", MAXCRYSTALS, 1, MAXCRYSTALS + 1, 
  		      5000, 0, 5000);
    
  } else if (ctrl->xtalkAnalysis) {
    /* These are the histograms for a cross-talk determination,
       using the method of multiplicity 2 events. */
    Int_t xtalIDmin = (ctrl->xtalkID - 1)*40;
    
    char str[300];

    for (i=0; i<40; i++) {
      sprintf(str, "eraw%3.3i", i+xtalIDmin);
      erawXT[i] = new TH1F(str, str, 35000, 0, 35000);
      erawXT[i]->SetXTitle("channels");
      erawXT[i]->SetYTitle("counts");
      sprintf(str, "ecal%3.3i", i+xtalIDmin);
      ecalXT[i] = new TH1F(str, str, 10000, 0, 10000);
      ecalXT[i]->SetXTitle("keV");
      ecalXT[i]->SetYTitle("counts");
    }
  
    for (Int_t j=0; j<40; j++) {
      for (Int_t k=j+1; k<40; k++) {
  	sprintf(str, "xtal%i_%i_%i", ctrl->xtalkID, j, k);
  	mult2Pair[j][k] = new TH1F(str, str, 10000, 0, 10000);
  	mult2Pair[j][k]->SetXTitle("channels ~ keV");
  	mult2Pair[j][k]->SetYTitle("counts");
  	sprintf(str, "xtal%i_%i_%i_2D", ctrl->xtalkID, j, k);
  	mult2Pair_2D[j][k] = new TH2F(str, str, 700,0,1400,700,0,1400);
      }
      sprintf(str, "xtal%i_%imult1_singles", ctrl->xtalkID, j);
      mult1Singles[j] = new TH1F(str, str, 35000, 0, 35000);
      mult1Singles[j]->SetXTitle("channels");
      mult1Singles[j]->SetYTitle("counts");
      sprintf(str, "xtal%i_%imult1_cc2", ctrl->xtalkID, j);
      mult1CC[j] = new TH2F(str, str, 2000,0,1000,2000,0,1000);
      mult1CC[j]->SetXTitle("segment channel");
      mult1CC[j]->SetYTitle("cc2 energy (keV)");
    }
    sprintf(str, "xtal%i_mult1", ctrl->xtalkID);
    mult1 = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult2", ctrl->xtalkID);
    mult2 = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult3", ctrl->xtalkID);
    mult3 = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult4", ctrl->xtalkID);
    mult4 = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult5", ctrl->xtalkID);
    mult5 = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_segSum", ctrl->xtalkID);
    segSum = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult1cal", ctrl->xtalkID);
    mult1cal = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult2cal", ctrl->xtalkID);
    mult2cal = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult3cal", ctrl->xtalkID);
    mult3cal = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult4cal", ctrl->xtalkID);
    mult4cal = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_mult5cal", ctrl->xtalkID);
    mult5cal = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_segSumcal", ctrl->xtalkID);
    segSumCal = new TH1F(str, str, 10000, 0, 10000);
    sprintf(str, "xtal%i_segmult", ctrl->xtalkID);
    segMult = new TH1F(str, str, 40, 0, 40);
  }    
}

/****************************************************/

void Histos::WriteHistos(controlVariables* ctrl) {

  Int_t i=0;

  if (!ctrl->calibration && !ctrl->xtalkAnalysis) {
    for (i=0; i<MAXCRYSTALS; i++) {
      cc1cal[i]->Write(); cc1cal[i]->Delete();
      cc2cal[i]->Write(); cc2cal[i]->Delete();
      cc3cal[i]->Write(); cc3cal[i]->Delete();
      cc4cal[i]->Write(); cc4cal[i]->Delete();
    }
    for (i=0; i<4; i++) {
      ccTot[i]->Write(); ccTot[i]->Delete();
      ccDop[i]->Write(); ccDop[i]->Delete();
    }
    ccSummary->Write(); ccSummary->Delete();
    ccDopSummary->Write(); ccDopSummary->Delete();
    mult->Write(); mult->Delete();
    ccHit->Write(); ccHit->Delete();
    ggMat->Write(); ggMat->Delete();
    ggMatDop->Write(); ggMatDop->Delete();
    for (i=0; i<MAXCRYSTALS; i++) {
      baseVsE[i]->Write(); baseVsE[i]->Delete();
      bankTSdiff[i]->Write(); bankTSdiff[i]->Delete();
      bankTSlow[i]->Write(); bankTSlow[i]->Delete();
      bankTShigh[i]->Write(); bankTShigh[i]->Delete();
    }
    xyProj->Write(); xyProj->Delete();
    xzProj->Write(); xzProj->Delete();
    yzProj->Write(); yzProj->Delete();
    xyz->Write(); xyz->Delete();
    thetaPhi->Write(); thetaPhi->Delete();
    thetaE->Write(); thetaE->Delete();
    thetaEDop->Write(); thetaEDop->Delete();
    chi2VsXtal->Write(); chi2VsXtal->Delete();
    for (i=0; i<MAXCRYSTALS; i++) {
      phiVsE[i]->Write(); phiVsE[i]->Delete();
    }
    deltaT->Write(); deltaT->Delete();
    
    test->Write(); test->Delete();
  } else if (ctrl->calibration && !ctrl->xtalkAnalysis) { /* ctrl->calibration */
    for (i=0; i<MAXCHANNELS; i++) {
      eraw[i]->Write(); eraw[i]->Delete();
      ecal[i]->Write(); ecal[i]->Delete();
      if (i%40 == 9) {
  	smoke[(4*(Int_t)(i/40))]->Write();
  	smoke[(4*(Int_t)(i/40))]->Delete();
      } else if (i%40 == 19) {
  	smoke[(4*(Int_t)(i/40))+1]->Write();
  	smoke[(4*(Int_t)(i/40))+1]->Delete();
      } else if (i%40 == 29) {
  	smoke[(4*(Int_t)(i/40))+2]->Write();
  	smoke[(4*(Int_t)(i/40))+2]->Delete();
      } else if (i%40 == 39) {
  	smoke[(4*(Int_t)(i/40))+3]->Write();
  	smoke[(4*(Int_t)(i/40))+3]->Delete();
      }
    }
    for (i=0; i<MAXCRYSTALS; i++) {
      xtalSummary[i]->Write(); xtalSummary[i]->Delete();
    }
    cc1All->Write(); cc1All->Delete();
    cc2All->Write(); cc2All->Delete();
    cc3All->Write(); cc3All->Delete();
    cc4All->Write(); cc4All->Delete();
    cc1Summary->Write(); cc1Summary->Delete();
    cc2Summary->Write(); cc2Summary->Delete();
    cc3Summary->Write(); cc3Summary->Delete();
    cc4Summary->Write(); cc4Summary->Delete();
  } else if (ctrl->xtalkAnalysis) {
    for (Int_t i=0; i<40; i++) {
      erawXT[i]->Write(); erawXT[i]->Delete();
      ecalXT[i]->Write(); ecalXT[i]->Delete();
      mult1Singles[i]->Write(); mult1Singles[i]->Delete();
      mult1CC[i]->Write(); mult1CC[i]->Delete();
      for (Int_t j=i+1; j<40; j++) {
  	mult2Pair[i][j]->Write(); mult2Pair[i][j]->Delete();
  	mult2Pair_2D[i][j]->Write(); mult2Pair_2D[i][j]->Delete();
      }
    }
    mult1->Write(); mult1->Delete();
    mult2->Write(); mult2->Delete();
    mult3->Write(); mult3->Delete();
    mult4->Write(); mult4->Delete();
    mult5->Write(); mult5->Delete();
    mult1cal->Write(); mult1cal->Delete();
    mult2cal->Write(); mult2cal->Delete();
    mult3cal->Write(); mult3cal->Delete();
    mult4cal->Write(); mult4cal->Delete();
    mult5cal->Write(); mult5cal->Delete();
    segMult->Write(); segMult->Delete();
    segSum->Write(); segSum->Delete();
    segSumCal->Write(); segSumCal->Delete();
  }
}
