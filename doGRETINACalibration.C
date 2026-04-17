#include "TGraph.h"

/***************** Fit Function Stuff **********************/

Double_t back(Double_t *x, Double_t *par) {
  /* Linear background - par[0] = offset, par[1] = slope */
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[0] + xx*par[1] < 0) { 
    f = 1000000000000;
  } else {
    f = (par[0] + xx*par[1]);
  }
  return f;
}

/***************************************************************/

Double_t plainGaus(Double_t *x, Double_t *par) {
  /* Gaussian - par[2] = amplitude, par[3] = amplitude scaling,
     par[4] = centroid, par[5] = sigma */
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[5] != 0) {
    f = ((par[2]*((100-par[3])/100))*
	 TMath::Exp(-0.5*((xx-par[4])/(par[5]))*((xx-par[4])/(par[5]))));
  } else {
    f = 1000000000000;
  } 
  return f;
}

/***************************************************************/

Double_t skewGaus(Double_t *x, Double_t *par) {
  /* Skewed Gaussian peak -- a la Radford gf3 */
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[5] != 0) {
    f = ((par[2]*par[3]/100)*
	 TMath::Exp((xx-par[4])/par[6])*
	 TMath::Erfc(((xx-par[4])/(TMath::Sqrt(2)*par[5])) +
		     ((par[5])/(TMath::Sqrt(2)*par[6]))));
  } else {
    f = 1000000000000;
  } 
  return f;
}

/***************************************************************/

Double_t stepFnc(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = 0;
  if (par[5] != 0) {
    f = ((par[7])*
	 TMath::Erfc((xx-par[4])/(TMath::Sqrt(2)*par[5])));
    if (f < 0) {
      f =  1000000000000;
    }
  } else {
    f = 1000000000000;
  }
  return f;
}

/***************************************************************/

Double_t skewFit(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = (back(x, &par[0]) + plainGaus(x, &par[0]) +
		skewGaus(x, &par[0]) + stepFnc(x, &par[0]));
  if (f < 0) { f = 100000000000; }
  return f;
}

/***************************************************************/

Double_t skewLimitP(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = (back(x, &par[0]) + plainGaus(x, &par[0]) +
		skewGaus(x, &par[0]) + stepFnc(x, &par[0]));
  if (f < 0) { f = 100000000000; }
  return TMath::Sqrt(f);
}

Double_t skewLimitM(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  Double_t f = (back(x, &par[0]) + plainGaus(x, &par[0]) +
		skewGaus(x, &par[0]) + stepFnc(x, &par[0]));
  if (f < 0) { f = 100000000000; }
  return -1*TMath::Sqrt(f);
}

/***************************************************************/

Double_t gausFit(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  return (back(x, &par[0]) + plainGaus(x, &par[0]) +
	  stepFnc(x, &par[0]));
}

/***************************************************************/

Double_t gausLimitP(Double_t *x, Double_t *par) {
  return TMath::Sqrt((back(x, &par[0]) + plainGaus(x, &par[0]) +
		      stepFnc(x, &par[0])));
}


Double_t gausLimitM(Double_t *x, Double_t *par) {
  return -1*TMath::Sqrt((back(x, &par[0]) + plainGaus(x, &par[0]) +
			 stepFnc(x, &par[0])));
}

/***************************************************************/

Double_t peakOnly(Double_t *x, Double_t *par) {
  Float_t xx = x[0];
  return (plainGaus(x, &par[0]) + skewGaus(x, &par[0]));
}

/***************************************************************/

Float_t dofit(TH1F *h, Double_t ilo, Double_t ihi, Int_t skewYes,
	       Double_t params[], Double_t params_err[]) {
  
  Double_t lo = ilo; Double_t hi = ihi;
  
  /* Retrieve the width of a bin, use bin 1 just because. */
  /** NOTE: invalid for variable width bins. **/
  Double_t binwidth = h->GetBinWidth(1);
  
  /* Set up the linear fit. */
  TF1 *backlin = new TF1("backlin",back,lo,hi, 2);
  backlin->SetLineColor(2);
  
  /* Take range passed into fitting routine and convert into bin numbers. */
  Int_t blo = Int_t (h->FindBin(lo));
  Int_t bhi = Int_t (h->FindBin(hi));

  /* Assume centroid of the gaussian is the bin with the maximum counts
     that is in range.... */
  Double_t maximum = 0;
  Double_t binmax = 0;
  for(Int_t i=blo; i<bhi; i++) {
    Double_t val = h->GetBinContent(i);
    if(val > maximum) { maximum = val;  binmax = i; }
  }

  /* Convert bin value back to an xaxis value. */
  Int_t guess_centroid = h->GetBinCenter(binmax);

  //   printf("binmax, binwidth, maximum = %f %f %f\n", binmax, binwidth, maximum);
  //printf("Starting centroid/height = %d/%f\n", guess_centroid, maximum);
    
  /* Create Gaussian for display and intergration purposes. */
  TF1 *gausD = new TF1("gausD", plainGaus, lo, hi, 8);
  gausD->SetLineColor(4);

  TF1 *skewgausD;
  if (skewYes) {
    skewgausD = new TF1("skewgausD", skewGaus, lo, hi, 8);
    skewgausD->SetLineColor(5);
  }
  
  TF1 *stepD = new TF1("stepD", stepFnc, lo, hi, 8);
  stepD->SetLineColor(6);
  
  /* Fit the combined spectrum with a gaussian and linear curves. */
  TF1 *userfit, *peaks, *limitP, *limitM;
  if (skewYes) {
    userfit = new TF1("userfit", skewFit, lo, hi, 8);
    peaks = new TF1("peaks", peakOnly, lo, hi, 8);
    limitP = new TF1("limitP", skewLimitP, lo, hi, 8);
    limitM = new TF1("limitM", skewLimitM, lo, hi, 8);
  } else {
    userfit = new TF1("userfit", gausFit, lo, hi, 8);
    peaks = new TF1("peaks", plainGaus, lo, hi, 8);
    limitP = new TF1("limitP", gausLimitP, lo, hi, 8);
    limitM = new TF1("limitM", gausLimitM, lo, hi, 8);
  }
  
  /* Create arrays for storing values. */
  Int_t npar = 8;
     
  /* Set the total fit function parameters to start the fit. */
  userfit->SetParameters(params);
 
  /* Background offset and slope */
  userfit->SetParLimits(0, -1., 10000); 
  userfit->SetParameter(0, h->GetBinContent(blo));
  userfit->SetParLimits(1, -100., 1.);
  userfit->SetParameter(1, 0.0);
    
  /* Total height, H */
  userfit->SetParLimits(2, 0.5*maximum, 2*maximum);
  userfit->SetParameter(2, maximum);

  /* Ratio of heights, R */
  userfit->SetParLimits(3, 0, 100);
  userfit->SetParameter(3, 50.);
  if (!skewYes) { userfit->FixParameter(3, 0); }

  /* Peak Position */
  userfit->SetParLimits(4, 0.99*guess_centroid, 1.01*guess_centroid);
  userfit->SetParameter(4, guess_centroid);

  /* Peak width (sigma) */
  userfit->SetParLimits(5, guess_centroid*0.0001, guess_centroid*0.02);
  userfit->SetParameter(5, guess_centroid*0.001);
  
  /* Beta for skewness */
  userfit->SetParLimits(6, 0.00001, 12000);
  userfit->SetParameter(6, 50.);
  if (!skewYes) { userfit->FixParameter(6, 0); }
  
  /* Step height */
  userfit->SetParLimits(7, 0, h->GetBinContent(blo)*1.1);
  userfit->SetParameter(7, h->GetBinContent(blo));
  
  /* Perform the total fit within a range, all weights set to 1 
     with loglike method. */
  h->Fit("userfit","WLMRQ");
  
  /* Retrieve the fit parameters and associated errors. */
  userfit->GetParameters(&params[0]);
  params_err[0] = userfit->GetParError(0);
  params_err[1] = userfit->GetParError(1);
  params_err[2] = userfit->GetParError(2);
  params_err[3] = userfit->GetParError(3);
  params_err[4] = userfit->GetParError(4);
  params_err[5] = userfit->GetParError(5);
  params_err[6] = userfit->GetParError(6);
  params_err[7] = userfit->GetParError(7);
  //cout << "Background: offset " << params[0] << " +/- " << params_err[0] << endl;
  //cout << "Background: slope " << params[1] << " +/- " << params_err[1] << endl;
  //cout << "Peak height: " << params[2] << " +/- " << params_err[2] << endl;
  //cout << "Ratio of heights (R): " << params[3] << " +/- " << params_err[3] << endl;
  //cout << "Peak Centroid: " << params[4] << " +/- " << params_err[4] << endl;
  //cout << "Peak Sigma: " << params[5] << " +/- " << params_err[5] << endl;
  //cout << "Peak Skewness (beta): " << params[6] << " +/- " << params_err[6] << endl;
  //cout << "Step height: " << params[7] << " +/- " << params_err[7] << endl;
  
  peaks->SetParameters(&params[0]);

  Double_t maxValue = peaks->GetMaximum();
  //cout << "Maximum value: " << maxValue << endl;
  Double_t fwhm, fwtm;
  Double_t valueLow = peaks->GetX(maxValue*0.5, lo, params[4]);
  Double_t valueHigh = userfit->GetX(maxValue*0.5, params[4], hi);
  fwhm = valueHigh - valueLow;
  //cout << "FWHM: " << fwhm << endl;
  valueLow = peaks->GetX(maxValue*0.1, lo, params[4]);
  valueHigh = peaks->GetX(maxValue*0.1, params[4], hi);
  fwtm = valueHigh - valueLow;
  //cout << "FW(1/10)M: " << fwtm << endl;
  //cout << " Ratio: " << fwtm/fwhm << endl;

  return fwhm;
}

void calibrate(TString fileName = "Run0000.root", TString sourceID="60Co", Int_t Q=20, Int_t P=1, Int_t doSegments=1, Int_t doDetMaps=0) {

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");
  gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

  /* Open the input file containing the energy histograms */
  TFile* finput = new TFile(fileName.Data(), "READ");
  fileName.Remove(fileName.Length()-5,5); /* Remove ".root" */

  Double_t nominalGain = 0.09; /* This is for the 10-MeV only GRETINA digitizers for ALL channels */
  
  Double_t ccGain[4] = {0};
  Double_t ccOffset[4] = {0};
  Double_t segGain[36] = {0};
  Double_t segOffset[36] = {0};
  
  
  Int_t ccChannel[4] = {169, 179, 189, 199};
  Int_t segChannel[36] = {160, 161, 162, 163, 164, 165, 166, 167, 168,
			  170, 171, 172, 173, 174, 175, 176, 177, 178,
			  180, 181, 182, 183, 184, 185, 186, 187, 188,
			  190, 191, 192, 193, 194, 195, 196, 197, 198};
  
  vector<pair<Double_t, Double_t>> ccPeaks[4];
  vector<pair<Double_t, Double_t>> segPeaks[36];

  Double_t minE = 500;
  Double_t maxE = 2000;
  Double_t stepE = 50;
  Double_t lastPeak = 0;
  Double_t lastMax = 0;
  
  Double_t energies[2] = {0., 0.};
  
  if (sourceID == "207Bi") { // 207Bi
    minE = 300;
    maxE = 2000;
    energies[0] = 569.;
    energies[1] = 1063.;
  } else if (sourceID == "60Co") { // 60Co
    minE = 800;
    maxE = 2000;
    energies[0] = 1173.2;
    energies[1] = 1332.5;
  }
  
  TH1F *ccSpec[4];
  TH1F *segSpec[36];

  FILE *calFile = fopen(Form("calFile_Q%dpos%d.txt", Q, P), "w");
  
  for (Int_t i=0; i<4; i++) {
    
    finput->cd();
    if (finput->Get(Form("eRaw%d", ccChannel[i])) != NULL) {
      ccSpec[i] = (TH1F*)finput->Get(Form("eRaw%d", ccChannel[i]));
    }

    lastPeak = 0; lastMax = 0;
        
    Float_t upperEdge = (minE+stepE)/nominalGain;
    ccSpec[i]->GetXaxis()->SetRangeUser(minE/nominalGain, upperEdge);
    Int_t binsPerStep = (Int_t)((stepE/nominalGain) / (ccSpec[i]->GetXaxis()->GetBinWidth(100)));
    while (upperEdge < maxE/nominalGain) {
      Float_t maxVal = ccSpec[i]->GetBinContent(ccSpec[i]->GetMaximumBin());
      //printf("maxVal %f\n", maxVal);
      if (maxVal > 5*(ccSpec[i]->Integral()/binsPerStep)) {
	if ((TMath::Abs(lastPeak - ccSpec[i]->GetBinCenter(ccSpec[i]->GetMaximumBin())) > 5/nominalGain)) {
	  //printf("Pushing; lastPeak = %f\n", lastPeak);
	  ccPeaks[i].push_back(make_pair(maxVal, ccSpec[i]->GetBinCenter(ccSpec[i]->GetMaximumBin())));
	  lastPeak = ccSpec[i]->GetBinCenter(ccSpec[i]->GetMaximumBin());
	  lastMax = maxVal;
	} else if (maxVal > lastMax) {
	  Int_t lastEntry = ccPeaks[i].size()-1;
	  ccPeaks[i][lastEntry].first = maxVal;
	  ccPeaks[i][lastEntry].second = ccSpec[i]->GetBinCenter(ccSpec[i]->GetMaximumBin());
	}
      }
      upperEdge += (stepE/nominalGain/2);
      ccSpec[i]->GetXaxis()->SetRangeUser(upperEdge-(stepE/nominalGain), upperEdge);
      //printf("New range %f %f", upperEdge-(stepE/nominalGain), upperEdge);
    }
    if (ccPeaks[i].size() > 0) {
      sort(ccPeaks[i].begin(), ccPeaks[i].end(), greater<>());
      if (ccPeaks[i][0].second > maxE/nominalGain) {
	ccPeaks[i].erase(ccPeaks[i].begin());
      }
      //printf("Found %zu candidates for CC %d\n", ccPeaks[i].size(), i);
      for (Int_t j=0; j<ccPeaks[i].size(); j++) {
	//printf("%d -- %f, %f\n", j, ccPeaks[i][j].first, ccPeaks[i][j].second);
      }
      
      Int_t closest = -1;
      Int_t distance = 100000;
      Int_t peakIndex[2] = {-1, -1};
      
      for (Int_t j=0; j<2; j++) {
	distance = 100000;
	closest = -1;
	for (Int_t m=0; m<ccPeaks[i].size(); m++) {
	  Double_t mdist = TMath::Abs(ccPeaks[i][m].second - (energies[j]/nominalGain));
	  if (mdist < distance) { distance = mdist; closest = m; }
	  //printf(" %d %f\n", m, mdist);
	}
	//printf("Closest for %f: %d\n", energies[j], closest);
	peakIndex[j] = closest;
      }
      
      //for (auto num : peaks[i]) {
      //  cout << num.first << " " << num.second << endl;
      //}
      
      Float_t centroids[2] = {0};
      Float_t sigma[2] = {0};
      
      peakIndex[0] = 0;  peakIndex[1] = 1;
      
      for (Int_t j=0; j<2; j++) {
	Double_t params[14];
	Double_t paramsErr[14];
	Double_t integrals[2];
	Double_t integralsErr[2];
	dofit(ccSpec[i], 0.95*ccPeaks[i][peakIndex[j]].second, 1.05*ccPeaks[i][peakIndex[j]].second, 1,
	      params, paramsErr);
	centroids[j] = params[4];
	sigma[j] = params[5];      
	//cout << " Fit peak[" << j << "] - " << centroids[j] << " " << sigma[j] << endl;
      }
      
      Double_t res = 0;
      
      if (centroids[1] > centroids[0]) {
	ccGain[i] = (energies[1] - energies[0])/(centroids[1]-centroids[0]);
	ccOffset[i] = energies[0] - (ccGain[i]*centroids[0]);
	res = sigma[1]*ccGain[i]*2.35;
      } else {
	ccGain[i] = (energies[1] - energies[0])/(centroids[0]-centroids[1]);
	ccOffset[i] = energies[0] - (ccGain[i]*centroids[1]);
	res = sigma[0]*ccGain[i]*2.35;
      }
      
      ccSpec[i]->GetXaxis()->SetRangeUser(minE/nominalGain, maxE/nominalGain);
      printf("CC %d\t%f\t%f\n", i, ccGain[i], ccOffset[i]);
      fprintf(calFile, "%d\t%f\t%f\n", ccChannel[i], ccGain[i], ccOffset[i]);
      //cout << " " << res << endl;
    } else {
      cout << "No peaks, no calibration." << endl;
    }
  }


  if (doSegments == 1) {
    
    for (Int_t i=0; i<36; i++) {
      finput->cd();
      if (finput->Get(Form("eRaw%d", segChannel[i])) != NULL) {
	segSpec[i] = (TH1F*)finput->Get(Form("eRaw%d", segChannel[i]));
      }
      
      Float_t upperEdge = (minE+stepE)/nominalGain;
      segSpec[i]->GetXaxis()->SetRangeUser(minE/nominalGain, upperEdge/nominalGain);
      Int_t binsPerStep = (Int_t)((stepE/nominalGain) / (segSpec[i]->GetXaxis()->GetBinWidth(100)));
      while (upperEdge < maxE/nominalGain) {
	Float_t maxVal = segSpec[i]->GetBinContent(segSpec[i]->GetMaximumBin());
	if (maxVal > 5*(segSpec[i]->Integral()/binsPerStep)) {
	  if ((TMath::Abs(lastPeak - segSpec[i]->GetBinCenter(segSpec[i]->GetMaximumBin())) > 5/nominalGain)) {
	    segPeaks[i].push_back(make_pair(maxVal, segSpec[i]->GetBinCenter(segSpec[i]->GetMaximumBin())));
	    lastPeak = segSpec[i]->GetBinCenter(segSpec[i]->GetMaximumBin());
	    lastMax = maxVal;
	  } else if (maxVal > lastMax) {
	    Int_t lastEntry = segPeaks[i].size()-1;
	    segPeaks[i][lastEntry].first = maxVal;
	    segPeaks[i][lastEntry].second = segSpec[i]->GetBinCenter(segSpec[i]->GetMaximumBin());
	  }
	}
	upperEdge += (stepE/nominalGain/2);
	segSpec[i]->GetXaxis()->SetRangeUser(upperEdge-(stepE/nominalGain), upperEdge);
      }
      if (segPeaks[i].size() > 0) {
	sort(segPeaks[i].begin(), segPeaks[i].end(), greater<>());
	if (segPeaks[i][0].second > maxE/nominalGain) {
	  segPeaks[i].erase(segPeaks[i].begin());
	}
	//printf("Found %zu candidates for CC %d\n", peaks[i].size(), i+36);
	for (Int_t j=0; j<segPeaks[i].size(); j++) {
	  //printf("%d -- %f, %f\n", j, peaks[i][j].first, peaks[i][j].second);
	}
	
	Int_t closest = -1;
	Int_t distance = 100000;
	Int_t peakIndex[2] = {-1, -1};
	
	for (Int_t j=0; j<2; j++) {
	  distance = 100000;
	  closest = -1;
	  for (Int_t m=0; m<segPeaks[i].size(); m++) {
	    Double_t mdist = TMath::Abs(segPeaks[i][m].second - (energies[j]/nominalGain));
	    if (mdist < distance) { distance = mdist; closest = m; }
	    //printf(" %d %f\n", m, mdist);
	  }
	  //printf("Closest for %f: %d\n", energies[j], closest);
	  peakIndex[j] = closest;
	}
	
	//for (auto num : peaks[i]) {
	//  cout << num.first << " " << num.second << endl;
	//}
	
	Float_t centroids[2] = {0};
	Float_t sigma[2] = {0};
	
	peakIndex[0] = 0;  peakIndex[1] = 1;
	
	for (Int_t j=0; j<2; j++) {
	  Double_t params[14];
	  Double_t paramsErr[14];
	  Double_t integrals[2];
	  Double_t integralsErr[2];
	  dofit(segSpec[i], 0.95*segPeaks[i][peakIndex[j]].second, 1.05*segPeaks[i][peakIndex[j]].second, 1,
		params, paramsErr);
	  centroids[j] = params[4];
	  sigma[j] = params[5];      
	  //cout << " Fit peak[" << j << "] - " << centroids[j] << " " << sigma[j] << endl;
	}
	Double_t res = 0;
	
	if (centroids[1] > centroids[0]) {
	  segGain[i] = (energies[1] - energies[0])/(centroids[1]-centroids[0]);
	  segOffset[i] = energies[0] - (segGain[i]*centroids[0]);
	  res = sigma[1]*segGain[i]*2.35;
	} else {
	  segGain[i] = (energies[1] - energies[0])/(centroids[0]-centroids[1]);
	  segOffset[i] = energies[0] - (segGain[i]*centroids[1]);
	  res = sigma[0]*segGain[i]*2.35;
	}
	
	segSpec[i]->GetXaxis()->SetRangeUser(minE/nominalGain, maxE/nominalGain);
	// cout << xtal << " " << i << " "  << ccGain[xtal][i] << " " << ccOffset[xtal][i];
	// cout << " " << res << endl;
	printf("Seg %d\t%f\t%f\n", i, segGain[i], segOffset[i]);
	fprintf(calFile, "%d\t%f\t%f\n", segChannel[i], segGain[i], segOffset[i]);
      } else {
	cout << "No peaks, no calibration." << endl;
      }
    }
  }

  fclose(calFile);
  
  printf("GOING INTO DETMAP\n");
  
  if (doSegments && doDetMaps) {
    gSystem->mkdir("DetMaps");
    FILE *detMapFile9 = fopen(Form("DetMaps/detmap_Q%dpos%d_CC9.txt", Q, P), "w");
    FILE *trGainFile9 = fopen(Form("DetMaps/tr_gain_Q%dpos%d_CC9.txt", Q, P), "w");
    FILE *detMapFile19 = fopen(Form("DetMaps/detmap_Q%dpos%d_CC19.txt", Q, P), "w");
    FILE *trGainFile19 = fopen(Form("DetMaps/tr_gain_Q%dpos%d_CC19.txt", Q, P), "w");
    FILE *detMapFile29 = fopen(Form("DetMaps/detmap_Q%dpos%d_CC29.txt", Q, P), "w");
    FILE *trGainFile29 = fopen(Form("DetMaps/tr_gain_Q%dpos%d_CC29.txt", Q, P), "w");
    FILE *detMapFile39 = fopen(Form("DetMaps/detmap_Q%dpos%d_CC39.txt", Q, P), "w");
    FILE *trGainFile39 = fopen(Form("DetMaps/tr_gain_Q%dpos%d_CC39.txt", Q, P), "w");
    Int_t ch = 0;
    for (Int_t i=0; i<9; i++) {
      fprintf(detMapFile9, "0x3 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile19, "0x3 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile29, "0x3 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile39, "0x3 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(trGainFile9, "%f\n", (ccGain[0]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile19, "%f\n", (ccGain[1]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile29, "%f\n", (ccGain[2]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile39, "%f\n", (ccGain[3]/segGain[i])/1.12964119028508558);
      ch++;
    }
    fprintf(detMapFile9, "0x3 \t 9 \t 0 \t 36 \t %f \t %f \n", ccOffset[0], ccGain[0]);
    ch = 0;
    for (Int_t i=9; i<18; i++) {
      fprintf(detMapFile9, "0x4 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile19, "0x4 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile29, "0x4 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile39, "0x4 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(trGainFile9, "%f\n", (ccGain[0]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile19, "%f\n", (ccGain[1]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile29, "%f\n", (ccGain[2]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile39, "%f\n", (ccGain[3]/segGain[i])/1.12964119028508558);
      ch++;
    }
    fprintf(detMapFile19, "0x4 \t 9 \t 0 \t 36 \t %f \t %f \n", ccOffset[1], ccGain[1]);
    ch = 0;
    for (Int_t i=18; i<27; i++) {
      fprintf(detMapFile9, "0x5 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile19, "0x5 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile29, "0x5 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile39, "0x5 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(trGainFile9, "%f\n", (ccGain[0]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile19, "%f\n", (ccGain[1]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile29, "%f\n", (ccGain[2]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile39, "%f\n", (ccGain[3]/segGain[i])/1.12964119028508558);
      ch++;
    }
    fprintf(detMapFile29, "0x5 \t 9 \t 0 \t 36 \t %f \t %f \n", ccOffset[2], ccGain[2]);
    ch = 0;
    for (Int_t i=27; i<36; i++) {
      fprintf(detMapFile9, "0x6 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile19, "0x6 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile29, "0x6 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(detMapFile39, "0x6 \t %d \t 0 \t %d \t %f \t %f \n", ch, i, segOffset[i], segGain[i]);
      fprintf(trGainFile9, "%f\n", (ccGain[0]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile19, "%f\n", (ccGain[1]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile29, "%f\n", (ccGain[2]/segGain[i])/1.12964119028508558);
      fprintf(trGainFile39, "%f\n", (ccGain[3]/segGain[i])/1.12964119028508558);
      ch++;
    }
    fprintf(detMapFile39, "0x6 \t 9 \t 0 \t 36 \t %f \t %f \n", ccOffset[3], ccGain[3]);
    fprintf(trGainFile9, "%f\n", (ccGain[0]/ccGain[0])/1.12964119028508558);
    fprintf(trGainFile19, "%f\n", (ccGain[1]/ccGain[1])/1.12964119028508558);
    fprintf(trGainFile29, "%f\n", (ccGain[2]/ccGain[2])/1.12964119028508558);
    fprintf(trGainFile39, "%f\n", (ccGain[3]/ccGain[3])/1.12964119028508558);
    fclose(detMapFile9);
    fclose(detMapFile19);
    fclose(detMapFile29);
    fclose(detMapFile39);
    fclose(trGainFile9);
    fclose(trGainFile19);
    fclose(trGainFile29);
    fclose(trGainFile39);
  }
  
}


void checkCalibration(TString fileName = "Run0000.root", TString sourceID="60Co", Int_t Q=20, Int_t P=1, Int_t doSeg=1) {

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Migrad");
  gROOT->ProcessLine("gErrorIgnoreLevel = 6001;");

  /* Open the input file containing the energy histograms */
  TFile* finput = new TFile(fileName.Data(), "READ");
  
  Int_t ccChannel[4] = {169, 179, 189, 199};
  Int_t segChannel[36] = {160, 161, 162, 163, 164, 165, 166, 167, 168,
			  170, 171, 172, 173, 174, 175, 176, 177, 178,
			  180, 181, 182, 183, 184, 185, 186, 187, 188,
			  190, 191, 192, 193, 194, 195, 196, 197, 198};
  
  TH1F *ccSpec[4];
  TH1F *segSpec[36];

  vector<Double_t> energies;

  FILE *fitFile = fopen(Form("resolutionFile_Q%dpos%d.txt", Q, P), "w");
  
  if (sourceID == "207Bi") { // 207Bi
    energies.push_back(569.);
    energies.push_back(1063.);
  } else if (sourceID == "60Co") { // 60Co
    energies.push_back(1173.2);
    energies.push_back(1332.5);
  } else if (sourceID == "152Eu") { // 152Eu
    energies.push_back(121.7817);
    energies.push_back(244.6974);
    energies.push_back(344.2785);
    energies.push_back(411.1165);
    energies.push_back(443.965);
    energies.push_back(778.9045);
    energies.push_back(867.380);
    energies.push_back(964.072);
    energies.push_back(1112.076);
    energies.push_back(1408.013);
  } else if (sourceID == "56Co") { // 56Co
    energies.push_back(846.7638);
    energies.push_back(1037.8333);
    energies.push_back(1771.327);
    energies.push_back(2598.438);
    energies.push_back(3253.402);
  }


  
  // CC first 
  for (Int_t i=0; i<4; i++) {

    finput->cd();
    if (finput->Get(Form("eCal%d", ccChannel[i])) != NULL) {
      ccSpec[i] = (TH1F*)finput->Get(Form("eCal%d", ccChannel[i]));
    }
          
    for (Int_t j=0; j<energies.size(); j++) {
      Double_t params[14];
      Double_t paramsErr[14];
      Double_t integrals[2];
      Double_t integralsErr[2];
      Double_t centroids[2];
      Double_t sigma[2];
      Float_t fwhm =  dofit(ccSpec[i], 0.98*energies[j], 1.02*energies[j], 1,
			    params, paramsErr);
      centroids[j] = params[4];
      sigma[j] = params[5];      
      
      // printf("%d %f %f %f %f\n", ccChannel[i], energies[j], centroids[j], centroids[j]-energies[j], fwhm);
      fprintf(fitFile, "%d,%f,%f,%f,%f\n", ccChannel[i], energies[j], centroids[j], centroids[j]-energies[j], fwhm);
            
    }
  }

  // Segments
  if (doSeg) {
    for (Int_t i=0; i<36; i++) {
      
      finput->cd();
      if (finput->Get(Form("eCal%d", segChannel[i])) != NULL) {
	segSpec[i] = (TH1F*)finput->Get(Form("eCal%d", segChannel[i]));
      }
      
      for (Int_t j=0; j<energies.size(); j++) {
	Double_t params[14];
	Double_t paramsErr[14];
	Double_t integrals[2];
	Double_t integralsErr[2];
	Double_t centroids[2];
	Double_t sigma[2];
	Float_t fwhm =  dofit(segSpec[i], 0.98*energies[j], 1.02*energies[j], 1,
			      params, paramsErr);
	centroids[j] = params[4];
	sigma[j] = params[5];      
	
	//printf("%d %f %f %f %f\n",  segChannel[i], energies[j], centroids[j], centroids[j]-energies[j], fwhm);
	fprintf(fitFile, "%d,%f,%f,%f,%f\n", segChannel[i], energies[j], centroids[j], centroids[j]-energies[j], fwhm);
	
      }
    }
  }

  fclose(fitFile);
}
  
  
