/**********************************************************************/
/* File: GRETINAWavefunction.C                                        */
/* Description: Functions for GRETINA waveform analysis               */
/* Author: H. Crawford                                                */
/* Date: January 2013                                                 */
/**********************************************************************/

#include "GRETINAWavefunction.h"
#include "TMath.h"

ClassImp(GRETINAWF);

/****************************************************/

GRETINAWF::GRETINAWF() {
  minCrossTime = 0; maxCrossTime = 2048;
  eFilterK = 100; eFilterM = 350;
  puFilterK = 10; puFilterM = 20; puFilterThresh = 150;
  
  segmentThresh = 25;
  
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<40; j++) {
      LEDcrossing[i][j] = 0.;
    }
  }
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    cc1BL[i] = 0; cc2BL[i] = 0;
    cc3BL[i] = 0; cc4BL[i] = 0;
  }
  for (Int_t i=0; i<MAXCHANNELS; i++) {
    restingBase[i] = 0;
    baseSamples[i] = 0;
    tau[i] = 0;
    baseline[i] = 0;
    sa[i] = 0; sb[i] = 0; sc[i] = 0;
    sd[i] = 0; se[i] = 0;
    for (Int_t j=0; j<3; j++) {
      for (Int_t k=0; k<4096; k++) {
	his[i][j][k] = 0;
      }
    }
  }
  
  minT = 140; maxT = 180; /* For short (186 sample) waveforms */
}

/****************************************************/

void GRETINAWF::Initialize() {
  minCrossTime = 0; maxCrossTime = 2048;
  
  segmentThresh = 25;
  
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<40; j++) {
      LEDcrossing[i][j] = 0.;
    }
  }
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    cc1BL[i] = 0; cc2BL[i] = 0;
    cc3BL[i] = 0; cc4BL[i] = 0;
  }
  for (Int_t i=0; i<MAXCHANNELS; i++) {
    restingBase[i] = 0;
    baseSamples[i] = 0;
    tau[i] = 0;
    baseline[i] = 0;
    sa[i] = 0; sb[i] = 0; sc[i] = 0;
    sd[i] = 0; se[i] = 0;
    for (Int_t j=0; j<3; j++) {
      for (Int_t k=0; k<4096; k++) {
	his[i][j][k] = 0;
      }
    }
  }
  
  minT = 120; maxT = 185;
}

/****************************************************/

void GRETINAWF::Reset() {
  minCrossTime = 0; maxCrossTime = 2048;
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<40; j++) {
      LEDcrossing[i][j] = 0.;
    }
  }
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    cc1BL[i] = 0; cc2BL[i] = 0;
    cc3BL[i] = 0; cc4BL[i] = 0;
    for (Int_t j=0; j<40; j++) {
      waveform2Out[i][j].clear();
    }
  }
  for (Int_t i=0; i<MAXCHANNELS; i++) {
    restingBase[i] = 0;
    baseSamples[i] = 0;
    tau[i] = 0;
    baseline[i] = 0;
  }
  waveformConCat.clear();
  waveE.clear();
  waveCE.clear();
}

/*****************************************************************/
/* ARGONNE ALGORITHMS FOR TRACE ANALYSIS - FROM S. ZHU           */
/*****************************************************************/

Int_t GRETINAWF::ANLAmplitude(Int_t T, Int_t cNum, Int_t chNum) {
  Int_t s1=0, s2=0;
  Int_t i=0;
  
  for (i=T; i<T-1; i++) {
    s1 += (Int_t)waveform2Out[cNum][chNum][i];
    s2 += (Int_t)waveform2Out[cNum][chNum][i + eFilterK];
  }
  s1/=4;
  s2/=4;
  return (s2-s1);
}

/****************************************************/

Int_t GRETINAWF::ANLBaseCal(Int_t T, Int_t cNum, Int_t chNum) {
  Int_t s=0;
  Int_t i=0;

  if (T<17) return 0;
  for (i=T-17; i<T-1; i++) { s += (Int_t)waveform2Out[cNum][chNum][i]; }
  s >>= 4;
  return s;
}

/****************************************************/

Double_t GRETINAWF::ANLEnergy1(Int_t T, Int_t id, Int_t cNum, Int_t chNum) {
  Double_t v[1024] = {0.}, s1=0., s2=0., e;
  Int_t i=0;
  
  v[0] = (Double_t)waveform2Out[cNum][chNum][0];
  for (i=1; i<1024; i++) {
    v[i] = ((Double_t)v[i-1] + (Double_t)waveform2Out[cNum][chNum][i] - 
	    (Double_t)waveform2Out[cNum][chNum][i-1]);
    v[i] += (((Double_t)waveform2Out[cNum][chNum][i] - baseline[id])/
	     (Double_t)tau[id]);
  }
  for (i = T + eFilterK; i < T + eFilterM + eFilterK; i++) {
    s1 += v[i - eFilterM - eFilterK];
    s2 += v[i];
  }
  e = (s2 - s1) / (Double_t)eFilterM * 3.0;
  return e;
}

/****************************************************/

Double_t GRETINAWF::ANLEnergy2(Int_t T, Int_t id, Int_t cNum, Int_t chNum) {
  Double_t s1=0., s2=0., s3=0.;
  Double_t e, c;
  Int_t i=0;
  
  c = (Double_t)exp(-(Double_t)(eFilterM + eFilterK)/(Double_t)tau[id]);
  for (i = T; i < T + eFilterM + eFilterK; i++) {
    s1 += (Double_t)waveform2Out[cNum][chNum][i - eFilterM - eFilterK] - baseline[id];
    s2 += (Double_t)waveform2Out[cNum][chNum][i] - baseline[id];
  }
  e = ((s2 - s1*c + (Double_t)eFilterM*s3/(Double_t)tau[id]) / 
       (Double_t)eFilterM*3.0);
  return e;
}

/****************************************************/

Float_t GRETINAWF::ANLEnergy(Int_t T, Int_t id, Int_t cNum, Int_t chNum) {
  Int_t b, i=0;
  Float_t s1=0, s2=0;
  Float_t e, c, C;

  c = (Float_t)exp(-(Float_t)(eFilterM+eFilterK)/(Float_t)tau[id]);
  for (i = T + eFilterK; i < T + eFilterM + eFilterK; i++) {
    s1 += (Float_t)waveform2Out[cNum][chNum][i - eFilterM - eFilterK] - baseline[id];
    s2 += (Float_t)waveform2Out[cNum][chNum][i] - baseline[id];
  } 
  e = (s2 - s1*c)/(Float_t)eFilterM*3.0;
#if(1)
  C = s1*c/(Float_t)eFilterM;
  b = ANLBaseCal(T, cNum, chNum);
#endif
  return e;
}

/****************************************************/

Float_t GRETINAWF::ANLRiseIntegral(Int_t T, Int_t id, Int_t cNum, Int_t chNum) {
  Float_t s=0;
  Float_t e;
  Int_t i=0;
  
  for (i = T; i < T + eFilterK; i++) {
    s += (Float_t)waveform2Out[cNum][chNum][i] - baseline[id] - (Float_t)waveform2Out[cNum][chNum][T];
  }
  e = s/(Float_t)eFilterK;
  return e;
}

/****************************************************/

Int_t GRETINAWF::LEDLevel(Int_t index, Int_t thresh, Int_t cNum, Int_t chNum) {
  Short_t d1, d2, d3, d4, d5, dd1, dd2, dd3;
  Short_t dd;
  d1 = waveform2Out[cNum][chNum][index-2] - waveform2Out[cNum][chNum][index - eFilterK - 2];
  d2 = waveform2Out[cNum][chNum][index-1] - waveform2Out[cNum][chNum][index - eFilterK - 1];
  d3 = waveform2Out[cNum][chNum][index] - waveform2Out[cNum][chNum][index - eFilterK];
  d4 = waveform2Out[cNum][chNum][index + 1] - waveform2Out[cNum][chNum][index - eFilterK + 1];
  d5 = waveform2Out[cNum][chNum][index + 2] - waveform2Out[cNum][chNum][index - eFilterK + 2];
  dd1 = (d1 + 2*d2 + d3)/4;
  dd2 = (d2 + 2*d3 + d4)/4;
  dd3 = (d3 + 2*d4 + d5)/4;
  dd = (dd1 + 2*dd2 + dd3)/4;
  if (thresh > 0) {
    if (dd > thresh) { return 1; }
    else { return 0; }
  } else {
    if (dd < thresh) { return 1; }
    else { return 0; }
  }
}

/****************************************************/

Int_t GRETINAWF::LED(Int_t index, Int_t thresh, Int_t cNum, Int_t chNum) {
  Int_t lled[1024] = {0};
  Int_t i=0;
  
  for (i=index+2; i<tracelength; i++) {
    lled[i] = LEDLevel(i, thresh, cNum, chNum);
    if ((lled[i] - lled[i-1])==1) { 
      return i; 
    }
  }
  return tracelength;
}

/*****************************************************************/
/* My fitting routine to extract baselines for our pulses.       */
/*****************************************************************/

Double_t GRETINAWF::ComputeChiSquare(Int_t startIndex, Int_t nPoints, 
				     Double_t tau, Double_t Nstart, 
				     Double_t baseline, Int_t cNum, Int_t chNum) {
  Double_t chiSq = 0;
  
  for (Int_t i=startIndex; i<(startIndex+nPoints); i++) {
    chiSq += ((waveform2Out[cNum][chNum][i] - (Nstart*TMath::Exp(-((i-startIndex)/tau))
				   + baseline))*
	      (waveform2Out[cNum][chNum][i] - (Nstart*TMath::Exp(-((i-startIndex)/tau))
				   + baseline)));
  }
  return chiSq;
}

/****************************************************/

Int_t GRETINAWF::LSFitExponential(Int_t startIndex, Int_t nPoints, 
				  Double_t Nstart, Double_t tau, 
				  Double_t &NstartFit, Double_t &baselineFit,
				  Double_t &NstartFitError, 
				  Double_t &baselineFitError, Int_t cNum, Int_t chNum) {
  
  /* Do a Levenberg-Marquardt minimization of a exponential decay + offset
     with a fixed tau, and get N0 and baseline(offset) value.
     
     Equation: N(t) = Nstart*exp(-(t-t0)/tau) + baseline

     We will define all matrices as:
        (  Nstart,Nstart     Nstart,baseline )
        ( baseline,Nstart   baseline,baseline)             */
  
  /* Set-up starting point for minimization...
     starting parameter, lambda, etc. */
  Double_t baselineStart = Nstart - 10; // Assume Nstart-10 is a decent guess...
  Double_t NstartStart = Nstart;
  Double_t baselineNow, baselineLast;
  Double_t NstartNow, NstartLast;

  Double_t lambda = 0.001;
  Double_t chiSqNow;
  Double_t chiSqLast;

  Double_t alpha[2][2] = {{0,0},{0,0}};
  Double_t beta[2] = {0};
  Double_t error[2][2] = {{0,0},{0,0}};
  Double_t stepping[2] = {0};

  Int_t notConverged = 1;
  Int_t iterations = 0;

  baselineNow = baselineStart;
  NstartNow = NstartStart;

  chiSqNow = ComputeChiSquare(startIndex, nPoints, tau, NstartNow, 
			      baselineNow, cNum, chNum);

  // cout << "Conditions: iter " << iterations << endl;
  // cout << "  Baseline now " << baselineNow << " Nstart now " 
  //      << NstartNow << endl;
  
  while (notConverged && iterations<10) {
    /* Clear the alpha, beta and error arrays... */
    alpha[0][0] = 0;  alpha[0][1] = 0;  alpha[1][0] = 0;  alpha[1][1] = 0;
    beta[0] = 0;  beta[1] = 0;
    error[0][0] = 0;  error[0][1] = 0;  error[1][0] = 0;  error[1][1] = 0;
    stepping[0] = 0;  stepping[1] = 0;

    /* Compute the curvature matrix and extremum vector first. */
    for (Int_t i=startIndex; i<(startIndex+nPoints); i++) {
      alpha[0][0] += TMath::Exp(-2*(i-startIndex)/tau);
      alpha[0][1] += TMath::Exp(-(i-startIndex)/tau);
      alpha[1][0] += TMath::Exp(-(i-startIndex)/tau);
      alpha[1][1] += 1;
      
      beta[0] += ((waveform2Out[cNum][chNum][i] -
		   (NstartNow*TMath::Exp(-(i-startIndex)/tau)) -
		   baselineNow)*
		  (TMath::Exp(-(i-startIndex)/tau)));
      beta[1] += ((waveform2Out[cNum][chNum][i] -
		   (NstartNow*TMath::Exp(-(i-startIndex)/tau)) -
		   baselineNow));
    }
    
    alpha[0][0] = alpha[0][0]*(1+lambda);
    alpha[1][1] = alpha[1][1]*(1+lambda);
    
    /* Now compute the error matrix and stepping vector. */
    error[0][0] = (1/((alpha[0][0]*alpha[1][1]) - 
		      (alpha[0][1]*alpha[0][1])))*alpha[1][1];
    error[0][1] = (1/((alpha[0][0]*alpha[1][1]) - 
		      (alpha[0][1]*alpha[0][1])))*(-alpha[0][1]);
    error[1][0] = (1/((alpha[0][0]*alpha[1][1]) - 
		      (alpha[0][1]*alpha[0][1])))*(-alpha[1][0]);
    error[1][1] = (1/((alpha[0][0]*alpha[1][1]) - 
		      (alpha[0][1]*alpha[0][1])))*alpha[0][0];

    stepping[0] = error[0][0]*beta[0] + error[0][1]*beta[1];
    stepping[1] = error[1][0]*beta[0] + error[1][1]*beta[1];

    /* Compute new values for variables. */
    NstartLast = NstartNow;
    NstartNow = NstartLast + stepping[0];
    baselineLast = baselineNow;
    baselineNow = baselineLast + stepping[1];

    chiSqLast = chiSqNow;
    chiSqNow = ComputeChiSquare(startIndex, nPoints, tau, 
				NstartNow, baselineNow, cNum, chNum);

    /* Evaluate our progress... */
    if (chiSqNow < chiSqLast) { // Success, decrease lambda
      lambda /= 10;
    } else {
      lambda *= 10;
      NstartNow = NstartLast;
      baselineNow = baselineLast;
    }

    iterations++;

    // cout << "Conditions: iter " << iterations << endl;
    // cout << "  Baseline now " << baselineNow << " Nstart now " << NstartNow << endl;
    // cout << "  Curvature:   " << alpha[0][0] << "     " << alpha[0][1] << endl;
    // cout << "               " << alpha[1][0] << "     " << alpha[1][1] << endl;
    // cout << "  Extremum:    " << beta[0] << "     " << beta[1] << endl;
    // cout << "  Error:       " << error[0][0] << "     " << error[0][1] << endl;
    // cout << "               " << error[1][0] << "     " << error[1][1] << endl;
    // cout << "  Stepping:    " << stepping[0] << "     " << stepping[1] << endl;
    // cout << "  Lambda       " << lambda      << " Chi2 " << chiSqNow << endl;

    if ( (TMath::Abs(stepping[0]) < (0.00001*TMath::Sqrt(error[0][0]))) &&
	 (TMath::Abs(stepping[1]) < (0.00001*TMath::Sqrt(error[1][1]))) ) {
      notConverged = 0;
    }
  }

  /* We're done, get the final answer... */
  alpha[0][0] = 0;  alpha[0][1] = 0;  alpha[1][0] = 0;  alpha[1][1] = 0;
  beta[0] = 0;  beta[1] = 0;
  error[0][0] = 0;  error[0][1] = 0;  error[1][0] = 0;  error[1][1] = 0;
  stepping[0] = 0;  stepping[1] = 0;
  
  /* Compute the curvature matrix and extremum vector first. */
  for (Int_t i=startIndex; i<(startIndex+nPoints); i++) {
    alpha[0][0] += TMath::Exp(-2*(i-startIndex)/tau);
    alpha[0][1] += TMath::Exp(-(i-startIndex)/tau);
    alpha[1][0] += TMath::Exp(-(i-startIndex)/tau);
    alpha[1][1] += 1;
    
    beta[0] += ((waveform2Out[cNum][chNum][i] -
		 (NstartNow*TMath::Exp(-(i-startIndex)/tau)) -
		 baselineNow)*
		(TMath::Exp(-(i-startIndex)/tau)));
    beta[1] += ((waveform2Out[cNum][chNum][i] -
		 (NstartNow*TMath::Exp(-(i-startIndex)/tau)) -
		 baselineNow));
  }
  
  /* Now compute the error matrix and stepping vector. */
  error[0][0] = (1/((alpha[0][0]*alpha[1][1]) - 
		    (alpha[0][1]*alpha[0][1])))*alpha[1][1];
  error[0][1] = (1/((alpha[0][0]*alpha[1][1]) - 
		    (alpha[0][1]*alpha[0][1])))*(-alpha[0][1]);
  error[1][0] = (1/((alpha[0][0]*alpha[1][1]) - 
		    (alpha[0][1]*alpha[0][1])))*(-alpha[1][0]);
  error[1][1] = (1/((alpha[0][0]*alpha[1][1]) - 
		    (alpha[0][1]*alpha[0][1])))*alpha[0][0];
  
  NstartFit = NstartNow;
  NstartFitError = TMath::Sqrt(error[0][0]);
  baselineFit = baselineNow;
  baselineFitError = TMath::Sqrt(error[1][1]);

  return iterations;
}

/****************************************************/

Double_t GRETINAWF::LSFitLinear(Int_t startIndex, Int_t nPoints,
				Double_t &slopeFit, Double_t &offsetFit) {

  /* Do a linear regression from startIndex to startIndex + nPoints */
  Double_t sumx = 0;
  Double_t sumy = 0;
  Double_t sumxx = 0;
  Double_t sumxy = 0;
  Double_t stddev = 0;
  Int_t n=0;

  for (n = startIndex; n <= (startIndex + nPoints); n++) {
    sumx += (Double_t)n;
    sumy += pzWaveform[n];
    sumxx += (Double_t)n*(Double_t)n;
    sumxy += (Double_t)n*pzWaveform[n];
  }

  slopeFit = (nPoints*sumxy - sumx*sumy)/(nPoints*sumxx - sumx*sumx);
  offsetFit = (sumy - slopeFit*sumx)/nPoints;
  for (n = startIndex; n<= (startIndex + nPoints); n++) {
    stddev += pow((pzWaveform[n] - (slopeFit*(Double_t)n + offsetFit)), 2.);
  }
  stddev = sqrt(stddev/(nPoints-1));

  return stddev;
}

/*****************************************************************/
/* Simple search routine to find baseline, basically look for    */
/* the slope to be zero on the baseline and the flat-top.        */
/*****************************************************************/

Float_t GRETINAWF::Baseline2Flatten(Int_t id, Float_t tau, 
				    Int_t startIndexBase, Int_t lengthBase,
				    Int_t startIndexTop, Int_t lengthTop,
				    Float_t smallStep, Int_t cNum, Int_t chNum) {

  Float_t pz = 0.0;  Float_t sum = 0.0;
  Int_t num = 0;  Int_t stddevnum = 0;
  
  Float_t mean = 0.0;  Float_t stddev = 0.0;
  Float_t stddevLow = 50000000.;  Float_t stddevLast = 50000000.;

  Float_t base = 0.0;  Float_t baseLow = 3000.0;

  Float_t step = 1.0;

  Int_t j=0; UInt_t uj=0;
  Int_t iter = 0;
  Int_t keepGoing = 1;

  for (j=startIndexBase; j<(startIndexBase + lengthBase); j++) {
    sum += (Float_t)waveform2Out[cNum][chNum][j]; 
    num++; 
  }
  base = sum/num;

  while (keepGoing && iter < 1000) {
    pzWaveform.clear();
    pz = 0.0;
    for (uj=0; uj<waveform2Out[cNum][chNum].size(); uj++) {
      pzWaveform.push_back(waveform2Out[cNum][chNum][uj] - base + pz*tau);
      pz += (waveform2Out[cNum][chNum][uj] - base);
    }

    num = 0;
    mean = 0.0;
    for (j=startIndexBase; j<(startIndexBase + lengthBase); j++) {
      mean += pzWaveform[j];
      num++;
    }
    mean = mean/num;
    
    stddevnum = 0;
    stddev = 0.0;
    for (j=startIndexBase; j<(startIndexBase + lengthBase); j++) {
      stddev += (mean - pzWaveform[j])*(mean - pzWaveform[j]);
      stddevnum++;
    }

    num = 0;
    mean = 0.0;
    for (j=startIndexTop; j<(startIndexTop + lengthTop); j++) {
      mean += pzWaveform[j];
      num++;
    }
    mean = mean/num;
    
    for (j=startIndexTop; j<(startIndexTop + lengthTop); j++) {
      stddev += (mean - pzWaveform[j])*(mean - pzWaveform[j]);
      stddevnum++;
    }
    
    stddev = stddev/stddevnum;

    if (stddev < stddevLow) {
      stddevLow = stddev;
      baseLow = base;
    }

    if (stddev > stddevLast) { keepGoing = 0; }
    stddevLast = stddev;

    if (keepGoing == 0 && step > smallStep) {
      step /= 10;
      base = baseLow + 20*step;
      keepGoing = 1;
    } else {
      base -= step;
    }

    iter++;
  }
  
  return baseLow;

}

/*****************************************************************/
/* This algorithm flattens the baseline to zero, at zero.        */
/*****************************************************************/

Int_t GRETINAWF::FlattenBaselineToZero(Float_t baseStart, Float_t ampStart, 
				       Float_t itau, Int_t startIndex, 
				       Int_t nPoints, Float_t &baseOut,
				       Float_t &ampOut, Int_t cNum, Int_t chNum) {
  Float_t pz = 0.0;  Int_t attempts = 0;
  Double_t stddev = 0.;  Double_t lineFit[2] = {0};

  Float_t ibase = baseStart;  Float_t iamp = ampStart;

  Int_t j=0;

  pzWaveform.clear();
  for (j=0; j<tracelength; j++) {
    pzWaveform.push_back(waveform2Out[cNum][chNum][j] - baseStart - ampStart + 
			 pz*itau);
    pz += waveform2Out[cNum][chNum][j] - baseStart;
  }

  stddev = LSFitLinear(startIndex, nPoints, lineFit[0], lineFit[1]);

  while ( ((fabs(lineFit[0]) > 0.005) || (fabs(lineFit[1]) > 0.01)) &&
	  (attempts < 100) ) {

    ibase += lineFit[1];
    pzWaveform.clear();
    pz = 0.0;
    for (j=0; j<tracelength; j++) {
      pzWaveform.push_back(waveform2Out[cNum][chNum][j] - ibase - iamp + 
			   pz*itau);
      pz += waveform2Out[cNum][chNum][j] - baseStart;
    }
    stddev = LSFitLinear(startIndex, nPoints, lineFit[0], lineFit[1]);

    attempts++;
  }
  
  baseOut = ibase;
  ampOut = iamp;

  return 1;

}

/*****************************************************************/
/* Algorithm from David for extracting energy from pulse... this */
/* requires knowledge of baseline...fit using above routines.    */
/*****************************************************************/

void GRETINAWF::RadfordBaselineFill(Int_t id, Int_t cNum, Int_t chNum) {
  Int_t s1=0, s2=0, s3=0, s4=0;
  Int_t dd[MAX_TRACE_LENGTH] = {0};
  Int_t pu = 0;
  
  Int_t j=0, k=0;
  
  /* Pile-up check first. Uses PU_* trapezoid to calculate
     signal derivative. */
  for (j=minT; j<minT+puFilterM; j++) {
    s1 += waveform2Out[cNum][chNum][j+puFilterM+puFilterK] - waveform2Out[cNum][chNum][j];
  }
  dd[minT] = s1;
  s2 += s1;
  s3++;
  for (j=minT; j<maxT-(2*puFilterM+puFilterK); j++) {
    s1 += (waveform2Out[cNum][chNum][j+2*puFilterM+puFilterK] -
	   waveform2Out[cNum][chNum][j+puFilterM+puFilterK] -
	   waveform2Out[cNum][chNum][j+puFilterM] + waveform2Out[cNum][chNum][j]);
    dd[j+1] = s1;
    s2 += s1;
    s3++;
  }

  /* Check for pile-up. */
  s2 /= s3;
  pu = 0;
  j = (minT + maxT - (2*puFilterM+puFilterK))/2; /* Midpoint of 
						  interesting trace 
						  segment. */
  for (k=minT; k<=maxT-(2*puFilterM+puFilterK); k++) {
    s4 = abs(s2 + s2*(j-k)/5000 - dd[k]); /* Rough PZ correction included. */
    if (pu < s4) pu = s4;
  }

  if (id%40 == 19) { /* Nominally 30MeV channel */
    if (pu > puFilterThresh/3) {
      return;
    }
  } else if (id%40 == 29) { /* Nominally 2.5MeV channel */
    if (pu > puFilterThresh*4) {
      return;
    }
  } else if (id%40 == 39) { /* Nominally 5MeV channel */
    if (pu > puFilterThresh*2) {
      return;
    }
  } else {
    if (pu > puFilterThresh) {
      return;
    }
  }

  /* No pileup!  Calculate difference in trace height over 30 
     samples (from preamp decay). */
  for (j=minT; j+45 < maxT; j+=30) {
    s1 = s2 = 0;
    for (k=j; k<j+15; k++) {
      s1 += waveform2Out[cNum][chNum][k];
      s2 += waveform2Out[cNum][chNum][k+30];
    }
    
    s3 = (10000 + s1)/15; /* Initial height + 500 */
    s4 = s1 - s2; /* 20 x change in height */
    
    if (s3 >= 0 && s3 < 4096) {
      his[id][1][s3] += s4/2;
      his[id][2][s3]++;

      /* Calculate sums for least-squares fit of preamp baseline. */
      sa[id] += s4;
      sb[id] += s3*s4;
      sc[id] += s3*s3;
      sd[id] += s3;
      se[id] += 1.0;
    }
    
  }

}

/*****************************************************************/

void GRETINAWF::RadfordBaselineFit(const char *file) {
  FILE *f;
  Double_t a, b, base, lo_ht, hi_ht;

  Int_t i=0, j=0;

  f = fopen(file, "w");
  fprintf(f, "# Input file = %s\n", file); 
  fprintf(f, "# minT, maxT = %d, %d\n", minT, maxT);
  fprintf(f, "# Pileup rise, flat, thresh = %d, %d, %d\n", 
	  puFilterK, puFilterM, puFilterThresh);
  fprintf(f, "#par_num     a        b     baseline  (+500)     decay    min_ht   max_ht\n");

  /* Do a least-squares fit of preamp baseline and decay time. */
  for (i=0; i<MAXCHANNELS; i++) {
    if (sa[i] > 0) {
      a = (sb[i]*se[i] - sa[i]*sd[i]) / (sc[i]*se[i] - sd[i]*sd[i]);
      if (a >= 0.5 || a <= 1.5) {
	b = (sb[i] - a*sc[i]) / sd[i];
	base = -b/a;
	lo_ht = (60.0/a) * 30.0/(maxT-minT); /* ~3 counts drop from minT to maxT */
	for (j=4094; j>500 && his[i][2][j]<5; j++) { ; }
	hi_ht = j - (Int_t)(base+0.5);
	/* Print results to file. */
	fprintf(f, " %6d %8.5f %9.2f %8.2f %8.2f %11.7f %6.1f %6.1f\n",
		i, a, b, base-500.0, base, 1.0 - pow((1.0 - a/20.0), 0.005), 
		lo_ht, hi_ht);
      }
    }
  }
  fclose(f);
}

/*****************************************************************/

Float_t GRETINAWF::RadfordEnergy(Int_t id, Int_t cNum, Int_t chNum) {
  
  UInt_t ui=0, uj=0;

  vector<Double_t> pzCorr;
  vector<Double_t> trap;
  Float_t pz = 0.0;

  Float_t sum = 0.0;

  /* Pole zero correction, assumes preamp baseline measured. */
  for (ui=0; ui<waveform2Out[cNum][chNum].size(); ui++) {
    pzCorr.push_back((waveform2Out[cNum][chNum][ui] - baseline[id] + 
		      pz*tau[id]));
    pz += waveform2Out[cNum][chNum][ui] - baseline[id];
  }

  for (ui=eFilterK; ui<pzCorr.size(); ui++) {
    trap.push_back(pzCorr[ui] - pzCorr[ui-eFilterK]); 
  }

  sum = 0;
  for (ui=0; ui<trap.size(); ui++) {
    sum += trap[ui];
    sum /= 2;
    if ((trap[ui] - sum) > 10) {
      sum = 0;
      
      for (uj=ui+25; uj<ui+45; uj++) {
	sum += trap[uj];
      }
      return sum;
    }
  }

  /* We should never get here... */
  return 0;
}

/*****************************************************************/
/* Algorithm from Ken for pulse shape analysis.                  */
/*****************************************************************/

GRETINAWF::linfit GRETINAWF::LinReg(Double_t *trace, Int_t startReg, Int_t stopReg) {
  /* Do a linear regression from startReg to stopReg */
  Double_t n = (Double_t)(stopReg - startReg);
  Double_t sumx = 0;
  Double_t sumy = 0;
  Double_t sumxx = 0;
  Double_t sumxy = 0;
  Double_t stddev = 0;
  Int_t nn=0;
  struct linfit mbs;

  for (nn = startReg; nn <= stopReg; nn++) {
    sumx += (Double_t)nn;
    sumy += trace[nn];
    sumxx += (Double_t)nn*(Double_t)nn;
    sumxy += (Double_t)nn*trace[nn];
  }

  mbs.m = (n*sumxy - sumx*sumy)/(n*sumxx - sumx*sumx);
  mbs.b = (sumy - mbs.m*sumx)/n;
  for (nn = startReg; nn<= stopReg; nn++) {
    stddev += pow((trace[nn] - (mbs.m*(Double_t)nn + mbs.b)), 2.);
  }
  stddev = sqrt(stddev/(n-1));
  mbs.s = stddev;

  return mbs;
}

GRETINAWF::peak GRETINAWF::GregorichTrapFilter(Int_t id, Int_t cNum, Int_t chNum) {
  /* Processes the trace and fills a structure with info:
     peak.amp contains amplitude of the largest peak
     peak.cen contains the centroid of the largest peak
     peak.base contains the calculated baseline of the trace
     Note: if peak.cen is 0., no peak was found                 */
  
  Int_t risen = 10; // typical pulse rise time
  Int_t basestart = 0; // first sample for baseline search (no pulse!)
  Int_t basestop = 300;  // lase sample for baseline search

  Double_t tracetrap[tracelength];
  Double_t peak[tracelength];
  
  Double_t expcorr[tracelength];
  Double_t tracecor1[tracelength];
  Double_t tracecor2[tracelength];

  Int_t startn;
  Int_t stopn;

  struct linfit baseline = {(Double_t)tracelength-1, 0., 0.};
  struct linfit topline = {(Double_t)tracelength-1, 0., 0.};
  struct peak result = {0, 0., 0., 0.};

  Double_t aa=4.; // multiply center channel counts by this (peak should be high)
  Double_t bb=1.; // multiply edge channel counts by this (edges should be low)
  Double_t cc=2.; // multiply abs(difference in edge channels)
  Double_t dd=0.25; // How far down trapezoid edges (0.5 for halfway, etc.)
  Int_t edgeoff;
  Int_t peakx = -1;
  Double_t peaky = 0;
  Double_t peakavg;
  Int_t halfgap = (eFilterK+2)/2;
  Int_t halffil = (halfgap+eFilterM);
  Int_t nn=0;

  for (nn=0; nn<tracelength; nn++) {
    tracetrap[nn] = 0;
    peak[nn] = 0;
  }

  /* Correct trace for exponential decay. */
  for (nn=0; nn<tracelength; nn++) {
    expcorr[nn] = exp((Double_t)nn/(Double_t)tau[id]);
    tracecor1[nn] = (Double_t)waveform2Out[cNum][chNum][nn]*expcorr[nn];
  }
  
  /* Baseline determination */
  baseline = LinReg(tracecor1, basestart, basestop);  
  result.base = (Int_t)(baseline.m*(Double_t)tau[id]*exp(((Double_t)basestart+
						  (Double_t)basestop)/-2./
						 (Double_t)tau[id]));
  
  /* Subtract baseline */
  for (nn=0; nn<tracelength; nn++) {
    tracecor2[nn] = ((Double_t)waveform2Out[cNum][chNum][nn]-result.base)*expcorr[nn]; 
  }

  startn = eFilterM + halfgap - 1;
  stopn = tracelength - eFilterM - halfgap + 1;
  
  /* Initialize trapezoidal filter at first point */
  for (nn=startn-halfgap-eFilterM+1; nn<=startn-halfgap; nn++) {
    tracetrap[startn] -= tracecor2[nn];
  }
  for (nn=startn+halfgap; nn<=startn+halfgap+eFilterM-1; nn++) {
    tracetrap[startn] += tracecor2[nn];
  }

  /* Run the filter over the rest of the trace */
  for (nn=startn+1; nn<=stopn; nn++) {
    tracetrap[nn] = tracetrap[nn-1] + tracecor2[nn-halffil] - tracecor2[nn-halfgap] - tracecor2[nn-1+halfgap] + tracecor2[nn+halffil-1];
  }

  /* Normalize the trapezoidal filter to give a gain of 1 */
  for (nn=startn; nn<=stopn; nn++) { tracetrap[nn] /= (Double_t)eFilterM; }
 
  /* Part way down the upper trapezoid side is this far from center. */
  edgeoff = (Int_t)((Double_t)(halfgap)+dd*(Double_t)eFilterM);
  
  /* Make a peak searcher run through the trapezoid */
  for (nn=startn+edgeoff; nn<stopn-edgeoff; nn++) {
    peak[nn] = aa*tracetrap[nn] - bb*tracetrap[nn-edgeoff] - bb*tracetrap[nn+edgeoff] - cc*fabs(tracetrap[nn-edgeoff] - tracetrap[nn+edgeoff]);
  }

  for (nn=1; nn<tracelength-1; nn++) {
    if (peak[nn]>1. && peak[nn-1]<peak[nn] && peak[nn+1]<peak[nn] && peak[nn]>peaky) {
      peaky = peak[nn];
      peakx = nn;
    }
  }

  /* Average the center of the trapezoid and divide out 
     exponential correction to get proper amplitudes.  Do a 
     linear regression on the flat top to see if it is a pileup. */
  if (peakx > -1) {
    peakavg = 0;
    for (nn = peakx-halfgap+1+3*risen; nn <= peakx+halfgap-1-3*risen; nn++) {
      peakavg += tracetrap[nn];
         }
    result.amp = peakavg/(Double_t)(eFilterK-6*risen)/expcorr[peakx];
    result.cen = (Double_t)peakx;
    topline = LinReg(tracetrap, peakx-halfgap+1+3*risen, peakx+halfgap-1-3*risen);
    if (fabs(topline.m)>0.03) { result.pileup = 1; }
    
  }

  return result;
}

/*****************************************************************/

void GRETINAWF::PoleZeroTrace(Int_t id, Int_t cNum, Int_t chNum) {
  Float_t pz = 0.0;
  pzWaveform.clear();
  for (Int_t i=0; i<tracelength; i++) {
    pzWaveform.push_back(waveform2Out[cNum][chNum][i] - restingBase[id] + pz*tau[id]);
    pz += waveform2Out[cNum][chNum][i] - restingBase[id];
  }
}

/*****************************************************************/

Float_t GRETINAWF::FPGAFilter(Int_t id, Int_t cNum, Int_t chNum) {
  Double_t e = 0;
  Double_t baseline[20] = {0};
  Double_t averagebaseline = 0;
  Double_t pzSum = 0;
  
  vector<Double_t> wave;
  vector<Double_t> pz;
  vector<Double_t> trap;
  
  for (UInt_t uj=0; uj<20; uj++) {
    baseline[uj] = waveform2Out[cNum][chNum][uj];
    averagebaseline += baseline[uj];
  }
  averagebaseline /= 20;
  
  averagebaseline = 0;
  
  for (Int_t j=0; j<tracelength; j++) {
    wave.push_back(waveform2Out[cNum][chNum][j] - averagebaseline);
  }
  
  pz.push_back((tau[id] + 1)*(wave[eFilterM - 1]));
  trap.push_back((wave[eFilterM - 1]));
  
  for (Int_t index = eFilterM; index < (Int_t)wave.size(); index++) {
    double scratch = ( (wave[index]) - 
		       ((index >= eFilterM) ? (wave[index - eFilterM]):0.0) -
		       ((index >= eFilterM + eFilterK) ? 
			(wave[index - eFilterM - eFilterK]):0.0) +
		       ((index >= 2*eFilterM + eFilterK) ? 
			(wave[index - 2*eFilterM - eFilterK]):0.0) );
    trap.push_back(trap[index - eFilterM - 1] + scratch);
    pzSum += trap[trap.size() - 1];
    pz.push_back(trap[index - eFilterM] + (pzSum)/tau[id]);
  }
  
  e = (pz[500] - pz[275])/32;
  
  pz.clear();
  trap.clear();
  wave.clear();
	    
  return e;	
}

/*****************************************************************/

Float_t GRETINAWF::SimpleEnergy(Int_t cNum, Int_t chNum) {
  Double_t sumA=0, sumB=0;

  Float_t e = 0;
  
  if (tracelength > 2000) {
    if (pzWaveform.size() > 0) {
      for (Int_t i=0; i<300; i++) { sumA += pzWaveform[i]; }
      for (Int_t i=500; i<1000; i++) { sumB += pzWaveform[i]; }
      sumA *= (1.666666666666667);
    } else {
      for (Int_t i=0; i<300; i++) { sumA += waveform2Out[cNum][chNum][i]; }
      for (Int_t i=500; i<1000; i++) { sumB += waveform2Out[cNum][chNum][i]; }
      sumA *= (1.666666666666667);
    }
  } else if (tracelength > 1000) {
    if (pzWaveform.size() > 0) {	 
      for (Int_t i=10; i<310; i++) { sumA += pzWaveform[i]; }
      for (Int_t i=500; i<1000; i++) { sumB += pzWaveform[i]; }
    } else {
      for (Int_t i=10; i<310; i++) { sumA += waveform2Out[cNum][chNum][i]; }
      for (Int_t i=500; i<1000; i++) { sumB += waveform2Out[cNum][chNum][i]; }
    }
  } else {
    if (pzWaveform.size() > 0) {	 
      for (Int_t i=10; i<60; i++) { sumA += pzWaveform[i]; }
      for (Int_t i=130; i<180; i++) { sumB += pzWaveform[i]; }
    } else {
      for (Int_t i=10; i<60; i++) { sumA += waveform2Out[cNum][chNum][i]; }
      for (Int_t i=130; i<180; i++) { sumB += waveform2Out[cNum][chNum][i]; }
    }
  }
  
  e = ((sumB) - (sumA));
  e /= 32.;
   
  return e;
}
 
/*****************************************************************/
/* CFD algorithm from Paul for getting T0 from a trace.          */
/*****************************************************************/

Float_t GRETINAWF::CalcT0(Int_t startSample=0, Int_t cNum=0, Int_t chNum=0) {
  Int_t maxj =0;
  Float_t sint[MAX_TRACE_LENGTH] = {0.};
  Float_t sdiff[MAX_TRACE_LENGTH] = {0.};
  Float_t s[MAX_TRACE_LENGTH] = {0.};

  Int_t i=0;
  UInt_t ui=0, uj=0;

  scfd.clear(); // Clear the global vector used to do the cfd

  Float_t start_cfd=0., start1=0., maxPulse=0., offset=0.;
  Float_t y1=0, x1=0, y2=0, x2=0;
  Float_t tcross=-1, tau_diff=6, delay=5, fraction=0.25;
  
  Float_t wave[MAX_TRACE_LENGTH] = {0.};
  
  /*----- Remove DC level -----*/
  offset = 0.;
  for (i=5; i<45; i++) { offset += (Float_t)waveform2Out[cNum][chNum][i]; }
  
  offset = offset/40;

  for (i=0; i<tracelength-1; i++) { 
    wave[i] = (Float_t)waveform2Out[cNum][chNum][i] - offset; 
  }
  
  /*----- Integrate and differentiate -----*/
  //  for (i=1; i<tracelength; i++) { sd[i] = wave[i]; }
  sint[0] = 0.333*(wave[0]*2 + wave[1]);
  for (i=1; i<tracelength-1; i++) { sint[i] = 0.25*(wave[i-1] + 2*wave[i] + wave[i+1] ); }
  s[0] = 0.333*(2*sint[0] + sint[1]);
  for (i=1; i<tracelength-1; i++) { s[i] = 0.25*(sint[i-1] + 2*sint[i] + sint[i+1]); }
  sdiff[0] = s[0];
  for (i=1; i<tracelength-1; i++) { sdiff[i] = ((s[i] - s[i-1]) + 
						       (sdiff[i-1] - (sdiff[i-1]/tau_diff)));}
  
  /*--- CFD ---*/
  for (i=(Int_t)delay; i<tracelength-1; i++) {
    if (i+delay < tracelength-1) {
      scfd.push_back((Short_t)(sdiff[i-(Int_t)delay] - fraction*sdiff[i]));
    } else {
      scfd.push_back((Short_t)(sdiff[i] - fraction*sdiff[i]));
    }
  }

  /*--- Find max pulse height ---*/
  maxPulse = 0; maxj = 0;
  for (uj=startSample; uj<scfd.size(); uj++) { 
    if (scfd[uj] > maxPulse) {
      maxPulse = scfd[uj];
      maxj = uj;
    }
  }
  
  /*--- Crossing point ---*/
  if (scfd.size() > (UInt_t)maxj) {
    if (scfd[maxj] < 0) { 
      for (uj=0; uj<scfd.size(); uj++) {
	scfd[uj] = -scfd[uj];
      }
    }
  
    for (i=maxj; i>0; i--) {
      if (scfd[i] < 0) { tcross = i; break; }
    }
  }

  ui = (UInt_t)tcross;
  if (scfd.size() > ui) { 
    x1 = ui;
    y1 = scfd[ui];
    x2 = ui+1;
    y2 = scfd[ui+1];
  
    start1 = x2 - (y2/((y2-y1)/(x2-x1)));
  } else {
    start1 = -10;
  }
  start_cfd = start1;
  return (start_cfd);
}

/*****************************************************************/
/* Pile-up check in the waveforms.  Assumes positive pulses (CC) */
/* for now.  Can extend later.                                   */ 
/*****************************************************************/

Int_t GRETINAWF::Look4Pileup(Int_t cNum, Int_t chNum) {
  
  Float_t wave[MAX_TRACE_LENGTH] = {0.};
  Float_t smoothed[MAX_TRACE_LENGTH] = {0.};
  Float_t smoothed2[MAX_TRACE_LENGTH] = {0.};
  
  Float_t sum1=0., sum2=0.;
  Int_t numPileUp = 0;

  Int_t i=0, j=0; 
  
  /*----- Remove DC level -----*/
  Float_t offset = 0;
  for (i=5; i<45; i++) { offset += waveform2Out[cNum][chNum][i]; }
    
  offset = offset/45.;
  
  for (i=0; i<tracelength; i++) { 
    wave[i] = (Float_t)waveform2Out[cNum][chNum][i] - offset; 
  }
    
  /* Smooth the pulses to get rid of the really bad noise that
     messes up the filter. */
  for (i=0; i<tracelength-1; i++) {
    smoothed[i+1] = 0.25*(wave[i] + 2*wave[i+1] + wave[i+2]);
  }
  for (i=0; i<tracelength-1; i++) {
    smoothed2[i+1] = 0.25*(smoothed[i] + 2*smoothed[i+1] + smoothed[i+2]);
  }
  
  for (i=0; i<tracelength-100; i++) {
    for (j=i; j<i+50; j++) {
      sum1 += smoothed2[j];
    }
    for (j=i+50; j<i+100; j++) {
      sum2 += smoothed2[j];
    }
    if ((sum2 - sum1)/50 > 25) {
      numPileUp++;
      i+=150; /* Skip over the rise time. */
    }
    sum1 = 0;
    sum2 = 0;
  }
  
  return numPileUp;

}

/*****************************************************************/
/* Mario timing algorithm for killing segments with bad timing   */ 
/* in preprocessor for Mode3->2 conversion.                      */
/*****************************************************************/

Int_t GRETINAWF::TimeSeg(Int_t len, Int_t cNum, Int_t chNum) {
  /* Returns index at which short energy filter has max value
     This is from David's pileup detector.                       */

  Int_t s1, s2, s3, s4;
  Int_t max_i;
  Float_t max;
  Int_t dd[1024];

  Int_t i=0;

  /* Sums for derivative / pile-up check */
  s1 = s2 = 0;
  s3 = 0;
  
  /* Use PU_* trapezoid to calculate derivative of signal. */
  for (i=minT; i<minT+puFilterM; i++) {
    /* First value of derivative */
    s1 += (waveform2Out[cNum][chNum][i+puFilterM+puFilterK] - 
	   waveform2Out[cNum][chNum][i]); 
  }
  dd[minT] = s1;
  s2 += s1;
  s3++;

  for (i=minT; i<maxT - (2*puFilterM+puFilterK); i++) {
    /* Remaining values of derivative. */
    s1 += (waveform2Out[cNum][chNum][i+2*puFilterM+puFilterK] - 
	   waveform2Out[cNum][chNum][i+puFilterM+puFilterK] -
	   waveform2Out[cNum][chNum][i+puFilterM] + 
	   waveform2Out[cNum][chNum][i]);
    dd[i+1] = s1;
    s2 += s1;
    s3++;
  }

  /* Check for maximum */
  s2 /= s3;
  max = -1.0;
  max_i = -1;
  for (i = minT; i<=maxT - (2*puFilterM+puFilterK); i++) {
    s4 = abs(s2 - dd[i]); /* Variation of derivative from average */
    if (max < s4) {
      max = s4;
      max_i = i;
    }
  }
  return max_i;

}
