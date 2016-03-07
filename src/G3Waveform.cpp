#include "GRETINA.h"

void g3Waveform::Clear() {
  raw.clear();
  pz.clear();
}

Int_t g3Waveform::LEDLevel(Int_t index, Int_t thresh, Int_t filterK) {
  Short_t d1, d2, d3, d4, d5, dd, dd1, dd2, dd3;
  d1 = raw[index-2] - raw[index-filterK-2];
  d2 = raw[index-1] - raw[index-filterK-1];
  d3 = raw[index] - raw[index-filterK];
  d4 = raw[index+1] - raw[index-filterK+1];
  d5 = raw[index+2] - raw[index-filterK+2];
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

Int_t g3Waveform::LED(Int_t index, Int_t thresh, Int_t filterK) {
  Int_t lled[1024] = {0};

  for (Int_t i=index+2; i<raw.size(); i++) {
    lled[i] = LEDLevel(i, thresh, filterK);
    if ((lled[i] - lled[i-1]) == 1) {
      return i;
    }
  }
  return raw.size();
}

/* CFD algorithm from Paul Fallon for getting T0 from the trace. */

Float_t g3Waveform::CFD(Int_t startSample = 0) {
  Int_t tracelength = raw.size();
  Int_t maxJ = 0;
  Float_t sint[MAX_TRACE_LENGTH] = {0.};
  Float_t sdiff[MAX_TRACE_LENGTH] = {0.};
  Float_t s[MAX_TRACE_LENGTH] = {0.};

  Float_t startCFD = 0., start1 = 0., maxPulse = 0., offset = 0.;
  Float_t y1 = 0., x1 = 0., y2 = 0., x2 = 0.;
  Float_t tcross = -1., tauDiff = 6., delay = 5., fraction = 0.25;

  Float_t wave[MAX_TRACE_LENGTH] = {0.};
  vector<Short_t> cfd;
  
  // cout << tracelength << endl;

  /* Remove DC level */
  Float_t baseline = BL(5, 45);
  for (Int_t i=0; i<tracelength; i++) {
    wave[i] = raw[i] - baseline;
  }

  /* Integrate and differentiate */
  sint[0] = (wave[0]*2 + wave[1])/3.;
  for (Int_t i=1; i<tracelength-1; i++) {
    sint[i] = (wave[i-1] + 2*wave[i] + wave[i+1])/4.;
  }
  s[0] = (2*sint[0] + sint[1])/3.;
  for (Int_t i=1; i<tracelength-1; i++) {
    s[i] = (sint[i-1] + 2*sint[i] + sint[i+1])/4.;
  }
  sdiff[0] = s[0];
  for (Int_t i=1; i<tracelength-1; i++) {
    sdiff[i] = ((s[i] - s[i-1]) +
		(sdiff[i-1] - (sdiff[i-1]/tauDiff)));
  }

  /* CFD */
  for (Int_t i = (Int_t)delay; i<tracelength-1; i++) {
    if (i+delay < tracelength-1) {
      cfd.push_back((Short_t)(sdiff[i - (Int_t)delay] - fraction*sdiff[i]));
    } else {
      cfd.push_back((Short_t)(sdiff[i] - fraction*sdiff[i]));
    }
  }

  /* Find max pulse height */
  maxPulse = 0; maxJ = 0;
  for (UInt_t uj = startSample; uj < cfd.size(); uj++) {
    if (cfd[uj] > maxPulse) {
      maxPulse = cfd[uj];
      maxJ = uj;
    }
  }

  /* Crossing point */
  if (cfd.size() > (UInt_t)maxJ) {
    if (cfd[maxJ] < 0) {
      for (UInt_t uj=0; uj<cfd.size(); uj++) {
	cfd[uj] = -cfd[uj];
      }
    }
    for (Int_t i=maxJ; i>0; i--) {
      if (cfd[i] < 0) { tcross = i; break; }
    }
  }

  UInt_t ui = (UInt_t)tcross;
  if (cfd.size() > ui) {
    x1 = ui; y1 = cfd[ui]; x2 = ui+1; y2 = cfd[ui+1];
    start1 = x2 - (y2/((y2-y1)/(x2-x1)));
  } else {
    start1 = -10;
  }
  startCFD = start1;
  cfd.clear();
  
  return (startCFD);
}

Float_t g3Waveform::riseTime(Float_t fLow, Float_t fHigh) {
  Float_t wave[MAX_TRACE_LENGTH] = {0.};
  Int_t tracelength = raw.size();

  /* Remove DC level */
  Float_t baseline = BL(5, 50);
  for (Int_t i=0; i<tracelength; i++) {
    wave[i] = raw[i] - baseline;
  }

  Float_t flatTop = BL(200,250);
  Float_t low = fLow * flatTop;
  Float_t high = fHigh * flatTop;

  Float_t tL = -1., tH = -1.;
  for (Int_t ui = 5; ui<tracelength; ui++) {
    if (wave[ui] > low) { tL = (Float_t)ui; break; }
  }
  for (Int_t ui = tracelength-5; ui > 0; ui--) {
    if (wave[ui] < high) { tH = (Float_t)ui; break; }
  }

  if (tL > 0 && tH > 0 && (tH > tL)) { return (tH - tL); }
  else { return -1; }
  
}

Float_t g3Waveform::BL(Int_t start, Int_t end) {
  Int_t sum = 0;
  for (Int_t i=start; i<end; i++) {
    sum += raw[i];
  }

  return ((Float_t)sum/(Float_t)(end - start));
}

Int_t g3Waveform::Look4Pileup() {
  Float_t wave[MAX_TRACE_LENGTH] = {0.};
  Float_t smooth[MAX_TRACE_LENGTH] = {0.};
  Float_t smooth2[MAX_TRACE_LENGTH] = {0.};
  Int_t tracelength = raw.size();

  Float_t sum1 = 0., sum2 = 0.;
  Int_t numPU = 0;

  /* Remove DC level */
  Float_t baseline = BL(5, 45);
  for (Int_t i=0; i<tracelength; i++) {
    wave[i] = raw[i] - baseline;
  }

  /* Smooth pulses to get rid of really bad noise... */
  for (Int_t i=0; i<tracelength-2; i++) {
    smooth[i+1] = (wave[i] + 2*wave[i+1] + wave[i+2])/4.;
  }
  for (Int_t i=0; i<tracelength-3; i++) {
    smooth2[i+1] = (smooth[i] + 2*smooth[i+1] + smooth[i+2])/4.;
  }

  for (Int_t i=0; i<tracelength-100; i++) {
    for (Int_t j=i; j<i+46; j++) {
      sum1 += smooth2[j];
    }
    for (Int_t j=i+50; j<i+46; j++) {
      sum2 += smooth2[j];
    }
    if ((sum2 - sum1)/46 > 25) {
      numPU++;
      i += 150; /* Skip over the rise time. */
    }
    sum1 = 0;  sum2 = 0;
  }

  return numPU;
  
}

void g3Waveform::calcPZ(Float_t base, Float_t tau) {
  Float_t ipz = 0.;
  Int_t tracelength = raw.size();

  pz.clear();
  for (Int_t i=0; i<tracelength; i++) {
    pz.push_back(raw[i] - base + ipz*tau);
    ipz += raw[i] - base;
  }
}
  

Float_t g3Waveform::simpleE() {
  Double_t sumA = 0., sumB = 0.;
  Float_t e = 0.;
  Int_t tracelength = raw.size();

  if (tracelength > 2000) {
    if (pz.size() > 0) {
      for (Int_t i=10; i<510; i++) { sumA += pz[i]; }
      for (Int_t i=650; i<1150; i++) { sumB += pz[i]; }     
    } else {
      for (Int_t i=10; i<510; i++) { sumA += raw[i]; }
      for (Int_t i=650; i<1150; i++) { sumB += raw[i]; }
    }
  } else if (tracelength > 1000) {
    if (pz.size() > 0) {
      for (Int_t i=10; i<310; i++) { sumA += pz[i]; }
      for (Int_t i=500; i<1000; i++) { sumB += pz[i]; }
      sumA *= (5./3.);
    } else {
      for (Int_t i=10; i<310; i++) { sumA += raw[i]; }
      for (Int_t i=500; i<1000; i++) { sumB += raw[i]; }
      sumA *= (5./3.);
    }
  } else {
    if (pz.size() > 0) {
      for (Int_t i=10; i<60; i++) { sumA += pz[i]; }
      for (Int_t i=130; i<180; i++) { sumB += pz[i]; }
    } else {
      for (Int_t i=10; i<60; i++) { sumA += raw[i]; }
      for (Int_t i=130; i<180; i++) { sumB += raw[i]; }
    }
  }

  e = ((sumB) - (sumA))/32.;
  return e;
}

Double_t g3Waveform::computeChiSquare(Int_t startIndex, Int_t nPoints,
				      Double_t tau, Double_t nStart,
				      Double_t baseline) {
  Double_t chiSq = 0;
  for (Int_t i=startIndex; i<(startIndex + nPoints); i++) {
    chiSq += ((raw[i] - (nStart*TMath::Exp(-((i-startIndex)/tau)) + baseline))*
	      (raw[i] - (nStart*TMath::Exp(-((i-startIndex)/tau)) + baseline)));
  }
  return chiSq;
  
}

Double_t g3Waveform::LSFitExpo(Int_t startIndex, Int_t nPoints,
			       Double_t nStart, Double_t tau,
			       Double_t &nStartFit, Double_t &baselineFit,
			       Double_t &nStartFitError,
			       Double_t &baselineFitError) {

  /* Levenberg-Marquardt minimization of exponential decay + offset
     with a fixed tau, to get N0 and baseline (offset) value. 
     
     Eqn: N(t) = nStart * Exp(-(t-t0)/tau) + baseline

     We will define all matrices as:
       ( nStart, nStart          nStart, baseline )
       ( baseline, nStart      baseline, baseline )    */

  /* Set-up starting point for minimization -- starting parameter, lambda... */
  Double_t baselineStart = nStart - 10;
  Double_t nStartStart = nStart;
  Double_t baselineNow, baselineLast;
  Double_t nStartNow, nStartLast;

  Double_t lambda = 0.001;
  Double_t chiSqNow, chiSqLast;
  Double_t alpha[2][2] = {{0.,0.},{0.,0.}};
  Double_t beta[2] = {0.};
  Double_t error[2][2] = {{0.,0.},{0.,0.}};
  Double_t stepping[2] = {0.};

  Int_t notConverged = 1;
  Int_t iterations = 0;

  baselineNow = baselineStart;
  nStartNow = nStartStart;
  chiSqNow = computeChiSquare(startIndex, nPoints, tau,
			      nStartNow, baselineNow);

  while (notConverged && iterations < 10) {
    /* Clear alpha, beta and error arrays */
    alpha[0][0] = 0.;  alpha[0][1] = 0.;  alpha[1][0] = 0.; alpha[1][1] = 0.;
    beta[0] = 0.;  beta[1] = 0.;
    error[0][0] = 0.;  error[0][1] = 0.;  error[1][0] = 0.; error[1][1] = 0.;
    stepping[0] = 0.;  stepping[1] = 0.;

    /* Compute the curvature matrix and extremum vector first */
    for (Int_t i = startIndex; i<(startIndex + nPoints); i++) {
      alpha[0][0] += TMath::Exp(-2*(i-startIndex)/tau);
      alpha[0][1] += TMath::Exp(-(i-startIndex)/tau);
      alpha[1][0] += TMath::Exp(-(i-startIndex)/tau);
      alpha[1][1] += 1.;

      beta[0] += ((raw[i] - (nStartNow*TMath::Exp(-(i-startIndex)/tau)) -
		   baselineNow)*(TMath::Exp(-(i-startIndex)/tau)));
      beta[1] += ((raw[i] - (nStartNow*TMath::Exp(-(i-startIndex)/tau)) -
		   baselineNow));
    }

    alpha[0][0] = alpha[0][0]*(1+lambda);
    alpha[1][1] = alpha[1][1]*(1+lambda);

    /* Compute error matrix and stepping vector */
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

    /* Compute new values for variables */
    nStartLast = nStartNow;
    nStartNow = nStartLast + stepping[0];
    baselineLast = baselineNow;
    baselineNow = baselineLast + stepping[1];
    chiSqLast = chiSqNow;
    chiSqNow = computeChiSquare(startIndex, nPoints, tau,
				nStartNow, baselineNow);

    /* Evaluate progress... */
    if (chiSqNow < chiSqLast) { lambda /= 10.; }
    else { lambda *= 10.; nStartNow = nStartLast; baselineNow = baselineLast; }

    iterations++;

    if ((TMath::Abs(stepping[0]) < (0.00001*TMath::Sqrt(error[0][0]))) &&
	(TMath::Abs(stepping[1]) < (0.00001*TMath::Sqrt(error[1][1]))) ) {
      notConverged = 0;
    }
  }

  /* We're done, get the final answer */
  
  /* Clear alpha, beta and error arrays */
  alpha[0][0] = 0.;  alpha[0][1] = 0.;  alpha[1][0] = 0.; alpha[1][1] = 0.;
  beta[0] = 0.;  beta[1] = 0.;
  error[0][0] = 0.;  error[0][1] = 0.;  error[1][0] = 0.; error[1][1] = 0.;
  stepping[0] = 0.;  stepping[1] = 0.;
  
  /* Compute the curvature matrix and extremum vector first */
  for (Int_t i = startIndex; i<(startIndex + nPoints); i++) {
    alpha[0][0] += TMath::Exp(-2*(i-startIndex)/tau);
    alpha[0][1] += TMath::Exp(-(i-startIndex)/tau);
    alpha[1][0] += TMath::Exp(-(i-startIndex)/tau);
    alpha[1][1] += 1.;
    
    beta[0] += ((raw[i] - (nStartNow*TMath::Exp(-(i-startIndex)/tau)) -
		 baselineNow)*(TMath::Exp(-(i-startIndex)/tau)));
    beta[1] += ((raw[i] - (nStartNow*TMath::Exp(-(i-startIndex)/tau)) -
		 baselineNow));
  }

  /* Compute error matrix  */
  error[0][0] = (1/((alpha[0][0]*alpha[1][1]) -
		    (alpha[0][1]*alpha[0][1])))*alpha[1][1];
  error[0][1] = (1/((alpha[0][0]*alpha[1][1]) -
		    (alpha[0][1]*alpha[0][1])))*(-alpha[0][1]);
  error[1][0] = (1/((alpha[0][0]*alpha[1][1]) -
		    (alpha[0][1]*alpha[0][1])))*(-alpha[1][0]);
  error[1][1] = (1/((alpha[0][0]*alpha[1][1]) -
		    (alpha[0][1]*alpha[0][1])))*alpha[0][0];
 
  nStartFit = nStartNow;
  nStartFitError = TMath::Sqrt(error[0][0]);
  baselineFit = baselineNow;
  baselineFitError = TMath::Sqrt(error[1][1]);
  return iterations;
}

Double_t g3Waveform::LSFitLinear(Int_t startIndex, Int_t nPoints,
				 Double_t &slopeFit, Double_t &offsetFit) {

  if (pz.size() > 0) {
    
  } else {
    calcPZ(BL(5,35), 50.);
  }
  
  /* Do linear regression from startIndex to startIndex + nPoints */
  Double_t sumX = 0., sumY = 0., sumXX = 0., sumXY = 0.;
  Double_t stdDev = 0.;
  
  for (Int_t n=startIndex; n<=(startIndex + nPoints); n++) {
    sumX += (Double_t)n;
    sumY += pz[n];
    sumXX += (Double_t)n*(Double_t)n;
    sumXY += (Double_t)n*pz[n];
  }
  slopeFit = (nPoints*sumXY - sumX*sumY)/(nPoints*sumXX - sumX*sumX);
  offsetFit = (sumY - slopeFit*sumX)/nPoints;
  for (Int_t n=startIndex; n<=(startIndex + nPoints); n++) {
    stdDev += pow((pz[n] - (slopeFit*(Double_t)n + offsetFit)), 2.);
  }
  stdDev = sqrt(stdDev/(nPoints-1));

  return stdDev;
}

Float_t g3Waveform::baseline2Flatten(Float_t tau, Int_t startIndexBase,
				     Int_t lengthBase, Int_t startIndexTop,
				     Int_t lengthTop, Float_t smallStep) {

  Float_t sum = 0.;
  Int_t num = 0, stdDevNum = 0;
  Float_t mean = 0., stdDev = 0.;
  Float_t stdDevLow = 50000000., stdDevLast = 50000000.;
  Float_t base = 0., baseLow = 3000.;
  Float_t step = 1.;

  Int_t iter = 0, keepGoing = 1;

  base = BL(startIndexBase, startIndexBase + lengthBase);
  while(keepGoing && iter < 1000) {
    pz.clear();
    calcPZ(base, tau);

    num = 0;  mean = 0.0;
    for (Int_t j=startIndexBase; j<(startIndexBase + lengthBase); j++) {
      mean += pz[j];
      num++;
    }
    mean = mean/num;

    stdDevNum = 0;  stdDev = 0.0;
    for (Int_t j=startIndexBase; j<(startIndexBase + lengthBase); j++) {
      stdDev += (mean - pz[j])*(mean - pz[j]);
      stdDevNum++;
    }

    num = 0;  mean = 0.0;
    for (Int_t j=startIndexTop; j<(startIndexTop + lengthTop); j++) {
      mean += pz[j];
      num++;
    }
    mean = mean/num;

    for (Int_t j=startIndexTop; j<(startIndexTop + lengthTop); j++) {
      stdDev += (mean - pz[j])*(mean - pz[j]);
      stdDevNum++;
    }
    stdDev = stdDev/stdDevNum;

    if (stdDev < stdDevLow) { stdDevLow = stdDev; baseLow = base; }
    if (stdDev > stdDevLast) { keepGoing = 0; }
    stdDevLast = stdDev;

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
