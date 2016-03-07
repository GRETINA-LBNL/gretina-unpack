/*! \file PhosWall.cpp
    \brief Parameters and functions for phoswich wall analysis.

    Author: H.L. Crawford
    Date: May 23, 2014
*/

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <sys/ioctl.h>
#include "TClass.h"

#include "PhosWall.h"
#include "colors.h"

/**********************************************************************/
/*  Coordinates class                                                 */
/**********************************************************************/

Coordinates::Coordinates() { ; }

/**********************************************************************/

Coordinates::Coordinates(Int_t nTele0, Int_t nPix0, TString fName, Bool_t printErrors) {
  nTele = nTele0; nPix = nPix0;
  
  if (printErrors) {
    cout << "PhosWall: # telescopes: " << nTele << ", # pixels: " << nPix << ", from file: "
	 << fName.Data() << endl;
  } else {
    cout << "PhosWall: Applying coordinates: " << fName.Data() << endl;
  }

  pars = new Pars*[nTele];
  for (Int_t i=0; i<nTele; i++) { pars[i] = new Pars[nPix]; }

  ifstream file; file.open(fName.Data());
  string title;
  getline(file, title);
  if (printErrors) { cout << "PhosWall: " << title << endl; }

  Int_t itele, ipix;
  Float_t axL, bxL, axR, bxR, ayD, byD, ayU, byU;
  for ( ; ; ) {
    file >> itele >> ipix >> axL >> bxL >> axR >> bxR >> ayD >> byD >> ayU >> byU;
    if (file.eof()) break;
    if (file.bad()) break;
    pars[itele][ipix].axL = axL;  pars[itele][ipix].bxL = bxL;
    pars[itele][ipix].axR = axR;  pars[itele][ipix].bxR = bxR;
    pars[itele][ipix].ayD = ayD;  pars[itele][ipix].byD = byD;
    pars[itele][ipix].axR = ayU;  pars[itele][ipix].bxR = byU;
    
    if (printErrors) {
      cout << "PhosWall: " << itele << ", " << ipix << ": " << axL << ", " << bxL << ", " 
	   << axR << ", " << bxR << ", " << ayD << ", " << byD << ", "
	   << ayU << ", " << byU << endl;
      cout << " -------------------------" << endl;
    }
  }

  file.close();  file.clear();

}

/**********************************************************************/

Coordinates::~Coordinates() {
  for (Int_t i=0; i<nTele; i++) { delete [] pars[i]; }
  delete pars;
}

/**********************************************************************/

Float_t Coordinates::algo(Int_t iTele, Int_t pixel, Int_t ipat, Float_t f, Float_t h) {
  Float_t p = -1.0;
  Float_t a, b;
  if (ipat == 0) { a = pars[iTele][pixel].axL;  b = pars[iTele][pixel].bxL; }
  if (ipat == 1) { a = pars[iTele][pixel].axR;  b = pars[iTele][pixel].bxR; }
  if (ipat == 2) { a = pars[iTele][pixel].ayD;  b = pars[iTele][pixel].byD; }
  if (ipat == 3) { a = pars[iTele][pixel].ayU;  b = pars[iTele][pixel].byU; }
  Float_t cm = a*(2.0 - pow(f,p)) + b*(f-0.50) + 6.08;
  return cm;
}

/**********************************************************************/
/* Neighbours class                                                   */
/**********************************************************************/

Neighbours::Neighbours() { ; }

/**********************************************************************/

Neighbours::Neighbours(Int_t nTele0, Int_t nPix0, TString fileName, Coordinates* coor0) {
  coor = coor0;  nTele = nTele0;  nPix = nPix0;
  
  coeffs = new Coeffs*[nTele];
  for (Int_t i=0; i<nTele; i++) { coeffs[i] = new Coeffs[nPix]; }
  ifstream file(fileName.Data());
  Int_t itele, ipix;
  Int_t up, dn, lf, rt, ul, ur, dl, dr;
  for ( ; ; ) {
    file >> itele >> ipix >> up >> dn >> lf >> rt >> ul >> ur >> dl >> dr;
    if (file.eof()) break;
    if (file.bad()) break;
    
    coeffs[itele][ipix].up = up;  coeffs[itele][ipix].dn = dn;
    coeffs[itele][ipix].lf = lf;  coeffs[itele][ipix].rt = rt;
    coeffs[itele][ipix].uL = ul;  coeffs[itele][ipix].uR = ur;
    coeffs[itele][ipix].dL = dl;  coeffs[itele][ipix].dR = dr;
  }

  file.close();  file.clear();  
}

/**********************************************************************/

Neighbours::~Neighbours() {
  for (Int_t i=0; i<nTele; i++) { delete [] coeffs[i]; }
  delete coeffs;
}

/**********************************************************************/

Float_t Neighbours::getBigOne(Int_t iTele, Int_t *pPix, Float_t *phA, Int_t size,
			      Int_t *pBigA, Float_t *bigA, Int_t maxIndex,
			      Int_t abc, Int_t iHist, Int_t LRDU) {
  if (maxIndex == 0) {
    for (Int_t i=0; i<size; i++) {
      if (*(phA+i) > *(phA + maxIndex)) { maxIndex = i; } 
    }
  }

  if (size == 1) {
    pBigA[0] = *pPix; /* Pixel of biggest PH */
    bigA[0] = *(phA); /* Biggest PH */
  } else if (size > 1) {
    pBigA[0] = *(pPix + maxIndex); /* Pixel of biggest PH */
    bigA[0] = *(phA + maxIndex); /* Biggest PH */
  }

  Int_t includeCorners = 0; /* = 1 include, = 0 skip */
  Float_t norm;
  if (includeCorners == 1) { norm = 2.92; }
  if (includeCorners == 0) { norm = 1.0; }
  
  Int_t pc = pBigA[0];
  bigA[1] = bigA[0]; 
  
  Float_t uA = 0.0, dA = 0.0, lA = 0.0, rA = 0.0;
  Float_t uL = 0.0, uR = 0.0, dL = 0.0, dR = 0.0;
  Float_t uF, dF, rF, lF, uLf, uRf, dLf, dRf;
  if (bigA[0] > 0) {
    for (Int_t i=0; i<size; i++) { /* Add in the neighbours PH's */
      Int_t p = *(pPix + i);
      if (p == coeffs[iTele][pc].up) { uA = *(phA+i);  bigA[1] = bigA[1] + uA; }
      if (p == coeffs[iTele][pc].dn) { dA = *(phA+i);  bigA[1] = bigA[1] + dA; }
      if (p == coeffs[iTele][pc].lf) { lA = *(phA+i);  bigA[1] = bigA[1] + lA; }
      if (p == coeffs[iTele][pc].rt) { dA = *(phA+i);  bigA[1] = bigA[1] + rA; }
      if (includeCorners) {
	if (p == coeffs[iTele][pc].uL) { uL = *(phA+i);  bigA[1] = bigA[1] + uL; }
	if (p == coeffs[iTele][pc].uR) { uR = *(phA+i);  bigA[1] = bigA[1] + uR; }
	if (p == coeffs[iTele][pc].dL) { dL = *(phA+i);  bigA[1] = bigA[1] + dL; }
	if (p == coeffs[iTele][pc].dR) { dR = *(phA+i);  bigA[1] = bigA[1] + dR; }
      }
    }
  } /* End if (bigA[0] > 0) */

  bigA[1] = bigA[1]/norm;
  if (uA == 0 && bigA[0] == 0) { uF = 0; }
  else { uF = uA/(uA + bigA[0]); } /* Up face fraction */
  if (dA == 0 && bigA[0] == 0) { dF = 0; }
  else { dF = dA/(dA + bigA[0]); } /* Down face fraction */
  if (lA == 0 && bigA[0] == 0) { lF = 0; }
  else { lF = lA/(lA + bigA[0]); } /* Left face fraction */
  if (rA == 0 && bigA[0] == 0) { rF = 0; }
  else { rF = rA/(rA + bigA[0]); } /* Right face fraction */

  /* Correct the edge cases now */
  Float_t pConst = 0.09;
  if ((uF == 0) && (dF != 0) && (rF != 0) && (lF != 0)) {
    uF = pConst/dF;
    bigA[1] = bigA[1] + bigA[0]*uF/(1+uF);
  }
  if ((uF != 0) && (dF == 0) && (rF != 0) && (lF != 0)) {
    dF = pConst/uF;
    bigA[1] = bigA[1] + bigA[0]*dF/(1+dF);
  }
  if ((uF != 0) && (dF != 0) && (rF == 0) && (lF != 0)) {
    rF = pConst/lF;
    bigA[1] = bigA[1] + bigA[0]*rF/(1+rF);
  }
  if ((uF != 0) && (dF != 0) && (rF != 0) && (lF == 0)) {
    lF = pConst/rF;
    bigA[1] = bigA[1] + bigA[0]*lF/(1+lF);
  }

  uLf = uL/(uL + bigA[0]);
  uRf = uR/(uR + bigA[0]);
  dLf = dL/(dL + bigA[0]);
  dRf = dR/(dR + bigA[0]);

  Int_t nTel0, nPix;
  Float_t xA = 0.0, yA = 0.0, xC = 0.0, yC = 0.0;
  Float_t ymm = 0.0, xmm = 0.0;
  
  Int_t left, right, up, down;
  left = 0; 
  right = 0;
  up = 0;
  down = 0;
  
  Float_t rat = 1.2;
  if (lF > rF) {
    xmm = 6.08 - coor->algo(iTele, pc, 0, lF, 9.0);
    if (right == 1) { xmm = 7.; }
    if (rF > 0.0) {
      Float_t xmm1 = coor->algo(iTele, pc, 1, rF, 9.0);
      xmm = (xmm + 0.5*xmm1)/1.5;
    }
  } /* End if (lF > rF) */
  if (rF > lF) {
    xmm = coor->algo(iTele, pc, 1, rF, 9.0);
    if (left == 1) { xmm = 7.; }
    if (lF > 0.0) {
      Float_t xmm2 = 6.08 - coor->algo(iTele, pc, 0, lF, 9.0);
      xmm = (xmm + 0.5*xmm2)/1.5;
    }
  } /* End if (rF > lF) */
  if (dF > uF) {
    ymm = 6.08 - coor->algo(iTele, pc, 2, dF, 9.0);
    if (up == 1) { ymm = 7.; }
    if (uF > 0.0) {
      Float_t ymm1 = coor->algo(iTele, pc, 3, uF, 9.0);
      ymm = (ymm + 0.5*ymm1)/1.5;
    }
  } /* End if (dF > uf) */
  if (uF > dF) {
    ymm = coor->algo(iTele, pc, 3, uF, 9.0);
    if (down == 1) { ymm = 7.; }
    if (dF > 0.0) {
      Float_t ymm2 = 6.08 - coor->algo(iTele, pc, 2, dF, 9.0);
      ymm = (ymm + 0.5*ymm2)/1.5;
    }
  } /* End if (uF > dF) */
  
  if (dF > uF && lF > rF) { yA = dF; xA = lF; }
  if (dF > uF && rF > lF) { yA = dF; xA = 0.5 + rF; }
  if (uF > dF && rF > lF) { yA = 0.5 + uF;  xA = 0.5 + rF; }
  if (uF > dF && lF > rF) { yA = 0.5 + uF;  xA = lF; }

  if (bigA[0] > 0) {
    if (LRDU == 1 || LRDU == 2) { /* 1 for plastic, 2 for CsI */
      bigA[2] = lF*100; /* Left cross-talk fraction on arrays A */
      bigA[3] = rF*100; /* Right cross-talk fraction on arrays B */
    } 
    if (LRDU == 3 || LRDU == 4) { /* 3 for plastic, 4 for CsI */
      bigA[2] = dF*100; /* Down cross-talk fraction on arrays A */
      bigA[3] = uF*100; /* Up cross-talk fraction on arrays B */
    }
    if (LRDU == 5) {
      bigA[2] = xA*100;
      bigA[3] = yA*100;
    }
  } else {
    bigA[2] = 0.0;
    bigA[3] = 0.0;
  }

  /* Correction for solid angle (x,y) dependence to be implemented */
  Float_t dX, d0 = 0.95;
  
  /* Left edge */
  if ((lF == 0.0) && (ymm <= 3.04)) {
    dX = d0/6 + d0*(1-xmm/6.08)*(ymm/3.04-1);
    xmm = xmm + dX;
  } else if ((lF == 0.0) && (ymm > 3.04) && (ymm <= 6.08)) {
    dX = d0/6 + d0*(1-xmm/6.08)*(1.0-ymm/3.04); 
    xmm = xmm + dX;
  }
  /* Right edge */
  if ((rF == 0.0) && (ymm <= 3.04)) {
    dX = d0/6 + d0*(xmm/6.08)*(1-ymm/3.04); 
    xmm = xmm + dX;
  } else if ((rF == 0.0) && ymm > 3.04 && ymm <= 6.08) {
    dX = d0/6 + d0*(xmm/6.08)*(ymm/3.04 - 1.0);
    xmm = xmm + dX;
  }
  d0 = d0/2.0;
  /* Center ones, note reduction in offset d0 */
  if ((lF > 0.0) && (rF > 0.0) && (xmm <= 3.04) && (ymm <= 3.04)) {
    dX = d0/4 + d0*(xmm/6.08)*(ymm/3.04-1); 
    xmm = xmm + dX;
  } else if ((lF > 0.0) && (rF > 0.0) && (xmm <= 3.04) && (ymm > 3.04)) {
    dX = d0/4 + d0*(xmm/6.08)*(1-ymm/3.04);
    xmm = xmm + dX;
  }
  if ((lF > 0.0) && (rF > 0.0) && (xmm > 3.04) && (ymm <= 3.04)) {
    dX = -d0/4 + d0*(1-xmm/6.08)*(1-ymm/3.04);
    xmm = xmm + dX;
  } else if ((lF > 0.0) && (rF > 0.0) && (xmm > 3.04) && (ymm > 3.04) && (ymm <= 6.08)) {
    dX = -d0/4 + d0*(1-xmm/6.08)*(ymm/3.04-1);
    xmm = xmm + dX;
  }
  
  bigA[4] = 10*xmm;
  bigA[4] = 10*ymm;

  /* Correct the edge pulse heights */
  Float_t k1 = 0.0751;  Float_t k2 = 0.00160;
  
  if (bigA[4] <= 30.4 && (pc <= 14 && pc >= 0 && pc%2 == 0)) {
    if (rA == 0 && rF == 0) {
      bigA[1] = bigA[1] + 0.0;
    } else {
      bigA[1] = bigA[1] + rA*(k2 + k1/rF);
    }
  }

  if (bigA[4] > 30.4 && (pc <= 31 && pc >= 17 && pc%2 == 1)) {
    if (lA == 0 && lF == 0) {
      bigA[1] = bigA[1] + 0.0;
    } else {
      bigA[1] = bigA[1] + lA*(k2 + k1/lF);
    }
  }
  
  if (bigA[5] <= 30.4 && (pc == 14 || pc == 15 || pc == 30 || pc == 31)) {
    if (uA == 0 && uF == 0) {
      bigA[1] = bigA[1] + 0.0;
    } else {
      bigA[1] = bigA[1] + uA*(k2 + k1/uF);
    }
  }

  if (bigA[5] > 30.4 && (pc == 0 || pc == 1 || pc == 16 || pc == 17)) {
    if (dA == 0 && dF == 0) {
      bigA[1] = bigA[1] + 0.0;
    } else {
      bigA[1] = bigA[1] + dA*(k2 + k1/dF);
    }
  }

  return 0;

}

/**********************************************************************/
/* Calibrate class                                                    */
/**********************************************************************/

Calibrate::Calibrate() { ; }

/**********************************************************************/

Calibrate::Calibrate(int nTele0, Int_t nPix0, TString fName, Bool_t printErrors) {
  nTele = nTele0;
  nPix = nPix0;

  if (printErrors) {
    cout << "# telescopes: " << nTele << ", # pixels: " << nPix << ", from file: " 
	 << fName.Data() << endl;
  }
  
  coeff = new Coeff* [nTele];
  for (Int_t i=0; i<nTele; i++) { coeff[i] = new Coeff[nPix]; }

  ifstream file(fName.Data());
  Int_t itele, ipix;
  Float_t x0, x1; /* Intercept and slope */
  for (;;) {
    file >> itele >> ipix >> x0 >> x1;
    if (file.eof()) break;
    if (file.bad()) break;
    coeff[itele][ipix].x1 = x1;
    coeff[itele][ipix].x0 = x0;

    if (printErrors) {
      cout << itele << ", " << ipix << ", " << coeff[itele][ipix].x0 << ", " 
	   << coeff[itele][ipix].x1 << endl;
    }
  }
  
  file.close();
  file.clear();
}

/**********************************************************************/

Calibrate::~Calibrate() {
  for (Int_t i=0; i<nTele; i++) { delete [] coeff[i]; }
  delete coeff;
}

/**********************************************************************/

Float_t Calibrate::getEnergy(Int_t iTele, Int_t iPix, Float_t channel, Int_t type, Int_t check) {
  Float_t x0 = coeff[iTele][iPix].x0; /* first point is pedestal */
  Float_t x1 = coeff[iTele][iPix].x1; /* second point is peak centroid */
  Float_t a, b, put, put0, put1, E;  Float_t shift = 0.0;
  Float_t eAlpha = 6096.0;
  Float_t r = ((Float_t)rand() / (Float_t)RAND_MAX) - 0.5; /* Random between -0.5 and 0.5 */
  
  /* Use in calibrations of the raw data for initial gain matching.
     For gain matching through plastic, x1 - x0 should be corrected for
     each pixel by dividing by the energy incident on the CsI: x0 or x1
     should be multiplied by 6.096/eAlphaPlastic */
  
  E = x1*(channel + 4*r - x0) + shift;
  put = 2000;
  put0 = 2000;
  put1 = 1000;
  
  b = (put1 - put0)/(x1 - x0);
  a = put1 - b*x1;

  Int_t gainSetting = 433;
 
  if (gainSetting == 433) {
    if (type == 1) { E = E*1000/2.5; } /* A: 433: 2.5, 322:  3.5 (all in keV/channel) */
    if (type == 2) { E = E*1000/10.0; } /* B: 433: 10.0, 322: 20.0 */
    if (type == 3) { E = E*1000/15.0; } /* C: 433: 15.0, 322: 26.0 */
  } else if (gainSetting == 322) {
    if (type == 1) { E = E*5000/3.5; }
    if (type == 2) { E = E*10000/20.0; }
    if (type == 3) { E = E*10000/26.0; }
  }

  if (type == 4) {
    channel = a + b*channel;
    return channel;
  } /* T */

  return E;
}
  
/**********************************************************************/
/* GainMatch class                                                    */
/**********************************************************************/

GainMatch::GainMatch() { ; }

/**********************************************************************/

GainMatch::GainMatch(Int_t nTele0, Int_t nPix0, TString fName, Bool_t printErrors) {
  nTele = nTele0;
  nPix = nPix0;

  if (printErrors) {
    cout << "# telescopes: " << nTele << ", # pixels: " << nPix << ", from file: " 
	 << fName.Data() << endl;
  }
  
  coeff = new Coeff* [nTele];
  for (Int_t i=0; i<nTele; i++) { coeff[i] = new Coeff[nPix]; }

  ifstream file(fName.Data());
  Int_t itele, ipix;
  Float_t x0, x1; /* Intercept and slope */
  for (;;) {
    file >> itele >> ipix >> x0 >> x1;
    if (file.eof()) break;
    if (file.bad()) break;
    coeff[itele][ipix].x1 = x1;
    coeff[itele][ipix].x0 = x0;

    if (printErrors) {
      cout << itele << ", " << ipix << ", " << coeff[itele][ipix].x0 << ", " 
	   << coeff[itele][ipix].x1 << endl;
    }
  }
  
  file.close();
  file.clear();
}

/**********************************************************************/

GainMatch::~GainMatch() {
  for (Int_t i=0; i<nTele; i++) { delete [] coeff[i]; }
  delete coeff;
}

/**********************************************************************/

Float_t GainMatch::GetGainEnergy(Int_t iTele, Int_t iPix, Float_t channel, Int_t type, Int_t check) {
  Float_t x0 = coeff[iTele][iPix].x0;
  Float_t x1 = coeff[iTele][iPix].x1;

  Float_t a, b, put, E;
  Float_t shiftG = 0.0;
  Float_t eAlpha = 8785.;
  Float_t r = ((Float_t)rand()/(Float_t)RAND_MAX)-0.5;
  
  if (type == 1) { put = eAlpha/20.0; } /* A: 20 9 = eAlpha/9.0 keV/channel for recalibration */
  if (type == 2) { put = eAlpha/10.0; } /* B: 10 1.5 = eAlpha/1.5 */
  if (type == 3) { put = eAlpha/15.0; } /* C: 15 2.0 = eAlpha/2.0 */
  if (type == 4) { E = channel;  return E; } /* Simple shifting in T */

  b = (put - x0)/x1;
  a = x0;
  
  E = a + b*(channel + r); 
  return E;
}

/**********************************************************************/
/* phosWallRaw class                                                  */
/**********************************************************************/

phosWallRaw::phosWallRaw() { ; }

/**********************************************************************/

phosWallRaw::~phosWallRaw() { ; }

/**********************************************************************/

void phosWallRaw::Initialize() {
  Reset();
}

/**********************************************************************/

void phosWallRaw::Reset() {
  hit.clear();
}

/**********************************************************************/

Int_t phosWallRaw::nHits() {
  return (Int_t)hit.size();
}

/**********************************************************************/
/* phosWallCalc class                                                 */
/**********************************************************************/

phosWallCalc::phosWallCalc() { ; }

/**********************************************************************/

phosWallCalc::~phosWallCalc() { ; }

/**********************************************************************/

void phosWallCalc::Initialize() {
  Reset();
}

/**********************************************************************/

void phosWallCalc::Reset() {
  hugeA = -1.0; hugeB = -1.0; hugeC = -1.0;
  uA = -1.0; dA = -1.0; lA = -1.0; rA = -1.0;
  uF = -1.0; dF = -1.0; lF = -1.0; rF = -1.0;
  ulA = -1.0; urA = -1.0; dlA = -1.0; drA = -1.0;
  ulF = -1.0; urF = -1.0; dlF = -1.0; drF = -1.0;
  hit.clear();
}

/**********************************************************************/

Float_t phosWallCalc::biggestA(TString var = "B") {
  Float_t biggestVar = -100.0;
  Float_t biggestA = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { biggestVar = hit[i].A; biggestA = hit[i].A; }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestA = hit[i].A; }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { biggestVar = hit[i].C; biggestA = hit[i].A; }
    } else {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestA = hit[i].A; }
    }
  }
  return biggestA;
}

/**********************************************************************/

Float_t phosWallCalc::biggestB(TString var = "B") {
  Float_t biggestVar = -100.0;
  Float_t biggestB = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { biggestVar = hit[i].A; biggestB = hit[i].B; }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestB = hit[i].B; }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { biggestVar = hit[i].C; biggestB = hit[i].B; }
    } else {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestB = hit[i].B; }
    }
  }
  return biggestB;
}

/**********************************************************************/

Float_t phosWallCalc::biggestC(TString var = "B") {
  Float_t biggestVar = -100.0;
  Float_t biggestC = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { biggestVar = hit[i].A; biggestC = hit[i].C; }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestC = hit[i].C; }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { biggestVar = hit[i].C; biggestC = hit[i].C; }
    } else {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestC = hit[i].C; }
    }
  }
  return biggestC;
}

/**********************************************************************/

Float_t phosWallCalc::biggestT(TString var = "B") {
  Float_t biggestVar = -100.0;
  Float_t biggestT = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { biggestVar = hit[i].A; biggestT = hit[i].T; }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestT = hit[i].T; }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { biggestVar = hit[i].C; biggestT = hit[i].T; }
    } else {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestT = hit[i].T; }
    }
  }
  return biggestT;
}

/**********************************************************************/

Int_t phosWallCalc::biggestPixel(TString var = "B") {
  Float_t biggestVar = -100.0;
  Int_t biggestPixel = -1;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { biggestVar = hit[i].A; biggestPixel = hit[i].pixel; }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestPixel = hit[i].pixel; }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { biggestVar = hit[i].C; biggestPixel = hit[i].pixel; }
    } else {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestPixel = hit[i].pixel; }
    }
  }
  return biggestPixel;
}

/**********************************************************************/

Int_t phosWallCalc::biggestHit(TString var = "B") {
  Float_t biggestVar = -100.0;
  Int_t biggestHit = -1;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { biggestVar = hit[i].A; biggestHit = i; }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestHit = i; }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { biggestVar = hit[i].C; biggestHit = i; }
    } else {
      if (hit[i].B > biggestVar) { biggestVar = hit[i].B; biggestHit = i; }
    }
  }
  return biggestHit;
}

/**********************************************************************/

Float_t phosWallCalc::secondBiggestA(TString var = "B") {
  Float_t biggestVar = -100.0;  Float_t secondBiggestVar = -100.0;
  Float_t biggestA = -100.0;  Float_t secondBiggestA = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestA = biggestA; 
	biggestVar = hit[i].A;  biggestA = hit[i].A; 
      } else if (hit[i].A > secondBiggestVar) {
	secondBiggestVar = hit[i].A;  secondBiggestA = hit[i].A;
      }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestA = biggestA; 
	biggestVar = hit[i].B;  biggestA = hit[i].A; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestA = hit[i].A;
      }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestA = biggestA; 
	biggestVar = hit[i].C;  biggestA = hit[i].A; 
      } else if (hit[i].C > secondBiggestVar) {
	secondBiggestVar = hit[i].C;  secondBiggestA = hit[i].A;
      }
    } else {
       if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestA = biggestA; 
	biggestVar = hit[i].B;  biggestA = hit[i].A; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestA = hit[i].A;
      }
    }
  }
  return secondBiggestA;
}

/**********************************************************************/

Float_t phosWallCalc::secondBiggestB(TString var = "B") {
  Float_t biggestVar = -100.0;  Float_t secondBiggestVar = -100.0;
  Float_t biggestB = -100.0;  Float_t secondBiggestB = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestB = biggestB; 
	biggestVar = hit[i].A;  biggestB = hit[i].B; 
      } else if (hit[i].A > secondBiggestVar) {
	secondBiggestVar = hit[i].A;  secondBiggestB = hit[i].B;
      }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestB = biggestB; 
	biggestVar = hit[i].B;  biggestB = hit[i].B; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestB = hit[i].B;
      }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestB = biggestB; 
	biggestVar = hit[i].C;  biggestB = hit[i].B; 
      } else if (hit[i].C > secondBiggestVar) {
	secondBiggestVar = hit[i].C;  secondBiggestB = hit[i].B;
      }
    } else {
       if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestB = biggestB; 
	biggestVar = hit[i].B;  biggestB = hit[i].B; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestB = hit[i].B;
      }
    }
  }
  return secondBiggestB;
}

/**********************************************************************/

Float_t phosWallCalc::secondBiggestC(TString var = "B") {
  Float_t biggestVar = -100.0;  Float_t secondBiggestVar = -100.0;
  Float_t biggestC = -100.0;  Float_t secondBiggestC = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestC = biggestC; 
	biggestVar = hit[i].A;  biggestC = hit[i].C; 
      } else if (hit[i].A > secondBiggestVar) {
	secondBiggestVar = hit[i].A;  secondBiggestC = hit[i].C;
      }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestC = biggestC; 
	biggestVar = hit[i].B;  biggestC = hit[i].C; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestC = hit[i].C;
      }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestC = biggestC; 
	biggestVar = hit[i].C;  biggestC = hit[i].C; 
      } else if (hit[i].C > secondBiggestVar) {
	secondBiggestVar = hit[i].C;  secondBiggestC = hit[i].C;
      }
    } else {
       if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestC = biggestC; 
	biggestVar = hit[i].B;  biggestC = hit[i].C; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestC = hit[i].C;
      }
    }
  }
  return secondBiggestC;
}

/**********************************************************************/

Float_t phosWallCalc::secondBiggestT(TString var = "B") {
  Float_t biggestVar = -100.0;  Float_t secondBiggestVar = -100.0;
  Float_t biggestT = -100.0;  Float_t secondBiggestT = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestT = biggestT; 
	biggestVar = hit[i].A;  biggestT = hit[i].T; 
      } else if (hit[i].A > secondBiggestVar) {
	secondBiggestVar = hit[i].A;  secondBiggestT = hit[i].T;
      }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestT = biggestT; 
	biggestVar = hit[i].B;  biggestT = hit[i].T; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestT = hit[i].T;
      }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestT = biggestT; 
	biggestVar = hit[i].C;  biggestT = hit[i].T; 
      } else if (hit[i].C > secondBiggestVar) {
	secondBiggestVar = hit[i].C;  secondBiggestT = hit[i].T;
      }
    } else {
       if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestT = biggestT; 
	biggestVar = hit[i].B;  biggestT = hit[i].T; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestT = hit[i].T;
      }
    }
  }
  return secondBiggestT;
}

/**********************************************************************/

Int_t phosWallCalc::secondBiggestPixel(TString var = "B") {
  Float_t biggestVar = -100.0;  Float_t secondBiggestVar = -100.0;
  Int_t biggestPixel = -100.0;  Int_t secondBiggestPixel = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestPixel = biggestPixel; 
	biggestVar = hit[i].A;  biggestPixel = hit[i].pixel; 
      } else if (hit[i].A > secondBiggestVar) {
	secondBiggestVar = hit[i].A;  secondBiggestPixel = hit[i].pixel;
      }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestPixel = biggestPixel; 
	biggestVar = hit[i].B;  biggestPixel = hit[i].pixel; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestPixel = hit[i].pixel;
      }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestPixel = biggestPixel; 
	biggestVar = hit[i].C;  biggestPixel = hit[i].pixel; 
      } else if (hit[i].C > secondBiggestVar) {
	secondBiggestVar = hit[i].C;  secondBiggestPixel = hit[i].pixel;
      }
    } else {
       if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestPixel = biggestPixel; 
	biggestVar = hit[i].B;  biggestPixel = hit[i].pixel; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestPixel = hit[i].pixel;
      }
    }
  }
  return secondBiggestPixel;
}

/**********************************************************************/

Int_t phosWallCalc::secondBiggestHit(TString var = "B") {
  Float_t biggestVar = -100.0;  Float_t secondBiggestVar = -100.0;
  Int_t biggestHit = -100.0;  Int_t secondBiggestHit = -100.0;
  for (int i=0; i<hit.size(); i++) {
    if (var.EqualTo('A')) {
      if (hit[i].A > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestHit = biggestHit; 
	biggestVar = hit[i].A;  biggestHit = i; 
      } else if (hit[i].A > secondBiggestVar) {
	secondBiggestVar = hit[i].A;  secondBiggestHit = i;
      }
    } else if (var.EqualTo('B')) {
      if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestHit = biggestHit; 
	biggestVar = hit[i].B;  biggestHit = i; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestHit = i;
      }
    } else if (var.EqualTo('C')) {
      if (hit[i].C > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestHit = biggestHit; 
	biggestVar = hit[i].C;  biggestHit = i; 
      } else if (hit[i].C > secondBiggestVar) {
	secondBiggestVar = hit[i].C;  secondBiggestHit = i;
      }
    } else {
       if (hit[i].B > biggestVar) { 
	secondBiggestVar = biggestVar;  secondBiggestHit = biggestHit; 
	biggestVar = hit[i].B;  biggestHit = i; 
      } else if (hit[i].B > secondBiggestVar) {
	secondBiggestVar = hit[i].B;  secondBiggestHit = i;
      }
    }
  }
  return secondBiggestHit;
}

/**********************************************************************/

Double_t phosWallCalc::distanceBetweenParticles(TString var = "B") {
  if ( (biggestHit(var) >= 0) && (secondBiggestHit(var) >= 0) ) {
    return ((hit[biggestHit(var)].Vec - hit[secondBiggestHit(var)].Vec).Mag());
  }
  return 0;
}

/**********************************************************************/
/* phosWallAux class                                                 */
/**********************************************************************/

phosWallAux::phosWallAux() { Reset(); }

/**********************************************************************/

phosWallAux::~phosWallAux() { ; }

/**********************************************************************/

void phosWallAux::Reset() {
  qdc.clear();
  tdc.clear();
}

/**********************************************************************/
/* phosWallFull class                                                 */
/**********************************************************************/

phosWallFull::phosWallFull() { Reset(); }

/**********************************************************************/

phosWallFull::~phosWallFull() { ; }

/**********************************************************************/

void phosWallFull::Initialize() {
  raw.Initialize();
  calc.Initialize();

  timestamp = 0;
  
  Int_t Ntel0 = 1;  Int_t Npix = 256;

  if (neighbourFile.EqualTo("")) {
    if (pixOrder == 1) { neighbourFile = "all_neighbors256-8x8.xy"; }
    else { neighbourFile = "all_neighbors256.xy"; }
  }
    
  ifExists(coefFileA);  ifExists(coefFileB);  
  ifExists(coefFileC);  ifExists(coefFileT);
  ifExists(coefFileAg);  ifExists(coefFileBg);  
  ifExists(coefFileCg);
  ifExists(coordFile);  ifExists(neighbourFile);  

  if (failCount > 0) { 
    cerr << endl << "PhosWall ----> Oops! " << failCount 
	 << " config file(s) missing. Aborting. " << endl; 
  } else if (failCount == 0) {
    cout << endl << "PhosWall ----> All config files present and accounted for. "
	 << "Continuing. " << endl;
  }
  cout << endl;
  cout << "PhosWall ----> Setting up calibration and gain matching. " << endl;

  Am = new Calibrate(Ntel0, Npix, coefFileA, showSetup);
  Bm = new Calibrate(Ntel0, Npix, coefFileB, showSetup);
  Cm = new Calibrate(Ntel0, Npix, coefFileC, showSetup);
  Tm = new Calibrate(Ntel0, Npix, coefFileT, showSetup);
  Ag = new GainMatch(Ntel0, Npix, coefFileAg, showSetup);
  Bg = new GainMatch(Ntel0, Npix, coefFileBg, showSetup);
  Cg = new GainMatch(Ntel0, Npix, coefFileCg, showSetup);
  myCoor = new Coordinates(Ntel0, Npix, coordFile, showSetup);

  ReadNeighbours(neighbourFile);
  ReadTimeGates(coefFileT);
}

/**********************************************************************/

Int_t phosWallFull::InitializeParameters(TString fileName) {
  
  FILE *phosIn;
  if ((phosIn = fopen(fileName.Data(), "r")) == NULL) {
    cerr << DRED << "PhosWall ---->  Control file " << fileName.Data() 
	 << " could not be opened." << RESET_COLOR << endl;
  }
  cout << "PhosWall: Setting controls based on " << fileName.Data() << endl;

  char line[300];  char junk[300];  char filename[300];
  Int_t value, value2;

  while ( !feof(phosIn) ) {
    fgets(line, 300, phosIn);
    if (strlen(line) == 1) { continue; }
    if (strncmp(line, "#", 1) == 0) { continue; }
    if (strncmp(line, ";", 1) == 0) { continue; }
    if (strncmp(line, "showSetup", 9) == 0) { 
      sscanf(line, "%s %d", junk, &value);
      showSetup = value;
    }
    if (strncmp(line, "eventPrinting", 13) == 0) {
      sscanf(line, "%s %d", junk, &value);
      doPrinting = value;
    }
    if (strncmp(line, "deBug", 5) == 0) {
      sscanf(line, "%s %d", junk, &value);
      deBug = value;
    }
    if (strncmp(line, "LRDU", 4) == 0) {
      sscanf(line, "%s %d", junk, &value);
      LRDU = value;
    } 
    if (strncmp(line, "pixOrder", 8) == 0) {
      sscanf(line, "%s %d", junk, &value);
      pixOrder = value;
    }
    if (strncmp(line, "minNeighbours", 13) == 0) {
      sscanf(line, "%s %d", junk, &value);
      minNeighbours = value;
    } 
    if (strncmp(line, "tLimits", 7) == 0) {
      sscanf(line, "%s %d %d", junk, &value, &value2);
      if (value < value2) { tLow = value; tHigh = value2; }
      if (value2 < value) { tLow = value2; tHigh = value; }
    }
    if (strncmp(line, "noTimes", 7) == 0) {
      sscanf(line, "%s %d", junk, &value);
      noTimes = value;
    }
    if (strncmp(line, "maxCheck", 8) == 0) {
      sscanf(line, "%s %d", junk, &value);
      maxCheck = value;
    }
    if (strncmp(line, "ACalibrationFile", 16) == 0) {
      printf(RED "line = %s" RESET_COLOR  "\n",line);
      sscanf(line, "%s %s", junk, filename);
      coefFileA = filename;
    }
    if (strncmp(line, "BCalibrationFile", 16) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coefFileB = filename;
    }
    if (strncmp(line, "CCalibrationFile", 16) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coefFileC = filename;
    }
    if (strncmp(line, "TCalibrationFile", 16) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coefFileT = filename;
    }
    if (strncmp(line, "AGainFile", 9) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coefFileAg = filename;
    }
    if (strncmp(line, "BGainFile", 9) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coefFileBg = filename;
    }
    if (strncmp(line, "CGainFile", 9) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coefFileCg = filename;
    }
    if (strncmp(line, "CoordinatesFile", 15) == 0) {
      sscanf(line, "%s %s", junk, filename);
      coordFile = filename;
    }
    if (strncmp(line, "NeighbourFile", 13) == 0) {
      sscanf(line, "%s %s", junk, filename);
      neighbourFile = filename;
    }
    if (strncmp(line, "dataType", 8) == 0) {
      sscanf(line, "%s %d", junk, &value);
      dataType = value;
    }
  }
  fclose(phosIn);
  return 0;
}

/**********************************************************************/

void phosWallFull::Reset() {
  raw.Reset();
  calc.Reset();
  timestamp = 0;
}

/**********************************************************************/

Int_t phosWallFull::getAndUnpackPhosWall(FILE *inf, Int_t length) {
  Reset();
  
  Int_t siz = fread(buf, 1, length, inf);
  
  UShort_t *p = buf;
  
  Int_t nHits = length/10; /* Each hit = 10 bytes */

   hitEvent ihit;
  
  for (Int_t i=0; i<nHits; i++) {
    ihit.pixel = *p; p++;
    ihit.det = ihit.pixel/64 + 1;
    ihit.A = *p; p++;
    ihit.B = *p; p++;
    ihit.C = *p; p++;
    ihit.T = *p; p++;

    raw.hit.push_back(ihit);    
  }
  return 0;
}

/**********************************************************************/

Int_t phosWallFull::getAndUnpackPhosWallAux(FILE *inf, Int_t length) {
  aux.Reset();
  
  Int_t siz = fread(buf, 1, length, inf);
  
  UShort_t *p = buf;  

  auxHit ihit;

  Int_t n = *p; p++;

  for (Int_t i=0; i<n; i++) {
    ihit.channel = *p; p++;
    ihit.value = *p; p++;
    aux.qdc.push_back(ihit);
  }

  n = *p; p++;

  for (Int_t i=0; i<n; i++) {
    ihit.channel = *p; p++;
    ihit.value = *p; p++;
    aux.tdc.push_back(ihit);
  }

  return 0;

}

/**********************************************************************/

void phosWallFull::ProcessPhosWall() {
  calEvent icalc;

  if (dataType == 1) { /* All hits in data file */
   
    for (Int_t hitNum=0; hitNum<raw.nHits(); hitNum++) {
      
      if ( ((raw.hit[hitNum].A < 0) || (raw.hit[hitNum].B < 0) ||
	    (raw.hit[hitNum].C < 0) || (raw.hit[hitNum].T < 0)) ||
	   ((raw.hit[hitNum].A > 8192) || (raw.hit[hitNum].B > 8192) ||
	    (raw.hit[hitNum].C > 8192) || (raw.hit[hitNum].T > 8192)) ||
	   ((raw.hit[hitNum].pixel > 255) || (raw.hit[hitNum].det > 4)) ||
	   (raw.nHits() > 256) ) {
	badEvents++;
	oobEvents++;
	
	if (deBug) {
	  cout << "PhosWall: " << eventNum << " had out-of-bounds values. " 
	       << "Skipping. " << endl;
	}
      } /* Protection against crazy values */
      
      if (raw.nHits() == 1) {
	oneEvents++;
	if (deBug) {
	  cout << "PhosWall: " <<  eventNum << " was a one-hit wonder. " 
	       << "Skipping. " << endl;
	}
      }
      
      Int_t timeStart = 0; Int_t timeEnd = 0;
      
      /* Set-up gain matching. */
      icalc.pixel = raw.hit[hitNum].pixel;
      icalc.A = Am->getEnergy(0, raw.hit[hitNum].pixel, raw.hit[hitNum].A, 1, 0);
      icalc.B = Bm->getEnergy(0, raw.hit[hitNum].pixel, raw.hit[hitNum].B, 2, 0);
      icalc.C = Cm->getEnergy(0, raw.hit[hitNum].pixel, raw.hit[hitNum].C, 3, 0);
      icalc.T = Tm->getEnergy(0, raw.hit[hitNum].pixel, raw.hit[hitNum].T, 4, 0);
      icalc.Vec = *pWallPositions.at(raw.hit[hitNum].pixel);
      calc.hit.push_back(icalc);

      if (0) {
	if (hitNum == raw.nHits()-1) {
	  cout << "GainMatched hits: " << endl;
	  for (int v=0; v<hitNum; v++) {
	    cout << "  " << v << ": " << calc.hit[v].pixel << ": (" << calc.hit[v].A << ", " << calc.hit[v].B << ", " << calc.hit[v].C << ") " << endl;
	  }
	  cout << "Biggest Hit " << calc.biggestHit() << " " << calc.biggestB() << endl;
	}
      }

      if (deBug) {
	if (hitNum == raw.nHits()-1) {
	  cout << endl;
	  cout << "PhosWall: " << " The biggest B pulse height of event " << eventNum 
	       << " was " << setw(5) << calc.biggestB() << " on pixel " << calc.biggestPixel() << ".";
	  cout << endl;
	  
	  for (Int_t v=0; v<raw.nHits(); v++) {
	    cout << "    Hit " << setw(3) << v << " : Pixel " << setw(3) << raw.hit[v].pixel
		 << " A: " << setw(4) << raw.hit[v].A << " B: " << setw(4) << raw.hit[v].B 
		 << " C: " << setw(4) << raw.hit[v].C << " T: " << setw(4) << raw.hit[v].T 
		 << " Am: " << setw(4) << calc.hit[v].A << " Bm: " << setw(4) << calc.hit[v].B
		 << " Cm: " << setw(4) << calc.hit[v].C << " Tm: " << setw(4) << calc.hit[v].T;
	    if (calc.biggestPixel("B") == raw.hit[v].pixel) { cout << " (B)"; }
	    if (isNeighbour(calc.biggestPixel("B"), raw.hit[v].pixel)) { cout << " (N)"; }
	    cout << endl;
	  }
	  
	  cout << " Events marked (N) are neighbours of the event's biggest hit (B)." << endl;
	}
      } /* End if (deBug) */
      
      if (hitNum == raw.nHits()-1 && raw.nHits() > 1) {
	Int_t nP = 0;
	nP = DetermineNParticles(20.0, calc.biggestHit());
	
	if (nP == 2) {
	  // cout << "Got two - event " << eventNum << endl;
	}

	Int_t neighbourCount = 0;
	
	for (Int_t check = 0; check < raw.nHits(); check++) {
	  if (isNeighbour(calc.biggestPixel("B"), raw.hit[check].pixel)) {
	    neighbourCount++;
	  }
	}
	if (neighbourCount < minNeighbours) {
	  if (deBug) {
	    cout << endl;
	    cout << "PhosWall: " << eventNum << " had too few neighbours (" 
		 << neighbourCount << "). Skipping. " << endl;
	  }
	} /* Check on the minimum number of neighbours */
	
	Int_t itele = 0;
	
	if ( (calc.biggestA("B") > 0) && (calc.biggestB("B") > 0) && (calc.biggestC("B") > 0) && (calc.biggestT("B") > 0) ) {
	  calc.hugeA = calc.biggestA("B"); calc.hugeB = calc.biggestB("B"); calc.hugeC = calc.biggestC("B");
	  calc.hugeA = correctNeighbours(calc.biggestPixel("B"), 1, 3, calc.biggestA("B"));
	  calc.hugeB = correctNeighbours(calc.biggestPixel("B"), 1, 3, calc.biggestB("B"));
	  calc.hugeC = correctNeighbours(calc.biggestPixel("B"), 1, 3, calc.biggestC("B"));
	  calc.hugeA = Ag->GetGainEnergy(itele, calc.biggestPixel("B"), calc.hugeA, 1, 1);
	  calc.hugeB = Bg->GetGainEnergy(itele, calc.biggestPixel("B"), calc.hugeB, 2, 0);
	  calc.hugeC = Cg->GetGainEnergy(itele, calc.biggestPixel("B"), calc.hugeC, 3, 0);
	}	
      } /* End if (hitNum == raw.nHits()) */
      
    } /* End loop over hits */
    
  } else if (dataType == 2) {
    if (raw.nHits() == 1) { 
      cout << "PhosWall: More than expected hits -- set for biggest hit report only." << endl;
    } 
  }

  eventNum++;
}

/**********************************************************************/

void phosWallFull::ProcessPhosWallAux() {

}

/**********************************************************************/

Bool_t phosWallFull::ifExists(TString fileName) {
  struct winsize w; struct stat st;
  ioctl(0, TIOCGWINSZ, &w);
  Int_t termColumns = w.ws_col;
  cout << "Checking " << fileName.Data();
  Int_t length = fileName.Length();
  for (Int_t i=0; i< termColumns - length - 20; i++) {
    cout << ".";
  }
  cout << ": ";
  if (stat(fileName.Data(), &st) == -1) {
    cout << "[FAILED]" << endl; failCount++;
    return false;
  } else {
    cout << "[OK]" << endl;
    return true;
  }
}

/**********************************************************************/

void phosWallFull::ReadNeighbours(TString file) {
  string line;
  Int_t lineCount = 0;
  Int_t throwAway, pixel, up, down, left, right;
  Int_t upleft, upright, downleft, downright;
  
  ifstream readFile;
  readFile.open(file.Data());
  if (!readFile.is_open()) {
    cout << "Could not open neighbours file " << file.Data() << " for parsing. " << endl;
    exit(1);
  }
  
  while ( lineCount <= 256 )  {
    getline(readFile, line);
    istringstream iss; iss.str(line);
    if (line[0] == '0') {
      iss >> throwAway >> pixel >> up >> down >> left >> right
	  >> upleft >> upright >> downleft >> downright;
    } /* End reading neighbour data -- valid lines begin with '0' */
    
    /* Populate array */
    neighbourhood[lineCount][0] = pixel;
    neighbourhood[lineCount][1] = up;        neighbourhood[lineCount][2] = down;
    neighbourhood[lineCount][3] = left;      neighbourhood[lineCount][4] = right;
    neighbourhood[lineCount][5] = upleft;    neighbourhood[lineCount][6] = upright;
    neighbourhood[lineCount][7] = downleft;  neighbourhood[lineCount][8] = downright;
    
    lineCount++;
  } /* End while (lineCount <= 256) */
  
  readFile.close();
  readFile.clear();
  
  if (showSetup) {
    cout << "Neighbour coordinates: " << endl;
    for (Int_t i=0; i<lineCount-1; i++) {
      for (Int_t e=0; e<9; e++) {
	cout << setw(5) << neighbourhood[i][e];
      } /* End e loop */
      cout << endl;
    } /* End i loop */
  }

}

/**********************************************************************/

void phosWallFull::ReadTimeGates(TString file) {
  string line;
  int lineCount = 0;
  int throwAway, pixel, start, end;
  
  ifstream readFile;
  readFile.open(file.Data());
  if (!readFile.is_open()) {
    cout << "Couldn't open the time gate file " << file.Data() << " for parsing." << endl;
    cout << "Opening all time gates from 0 to 8192. " << endl;
    for (Int_t i=0; i<256; i++) {
      timeGates[i][0] = i; /* Pixel */
      timeGates[i][1] = 0; /* Low end */
      timeGates[i][2] = 8192; /* High end */
    }
  }
  
  while(lineCount < 256) {
    getline(readFile, line);
    istringstream iss(line);
    if (line[0] == '0') {
      iss >> throwAway >> pixel >> start >> end;
    } /* Valid lines start with '0' */
    
    /* Populate array. */
    timeGates[lineCount][0] = pixel;
    timeGates[lineCount][1] = start;
    timeGates[lineCount][2] = end;
    
    lineCount++;
  }
  
  readFile.close();
  readFile.clear();
  
  if (showSetup) {
    cout << endl << "Time gates (" << file.Data() << ")" << endl << endl;
    for (int i=0; i<lineCount-1; i++) {
      cout << "Pixel " << setw(3) << timeGates[i][0] << ": ";
      cout << setw(4) << timeGates[i][1] << "-";
      cout << setw(4) << timeGates[i][2] << endl;
    } /* End i loop */
  }
  
}

/**********************************************************************/

Bool_t phosWallFull::isNeighbour(Int_t first, Int_t second) {
  Bool_t checkCorners = false;
  Bool_t result = false;
  if (first < 256) {
    if (!checkCorners) { maxCheck = 4; } else { maxCheck = 8; }
    for (Int_t n=1; n<= maxCheck; n++) {
      if (neighbourhood[first][n] == second) {
	result = true;
      } 
    }
  } else {
    cerr << "PhosWall ----> " << eventNum << ": Bad pixel passed to isNeighbour: " << first << endl;
  }
  return result;
}

/**********************************************************************/

Float_t phosWallFull::correctNeighbours(Int_t pixel, Int_t type, Int_t iHist, Float_t input) {
  /* Change normalization depending on if corners are checked */
  Float_t norm;
  if (maxCheck == 4) { norm = 1.0; } else if (maxCheck == 8) { norm = 2.92; }
  Float_t bigZero = input; 
  Float_t bigOne = bigZero; 

  /* Correct for "edge cases" (missing neighbours */
  for (Int_t b=0; b<calc.hit.size(); b++) {
    if (calc.hit[b].pixel == neighbourhood[pixel][1]) {
      if (type == 1) { calc.uA = calc.hit[b].A;  bigOne = bigOne + calc.uA; }
      if (type == 2) { calc.uA = calc.hit[b].B;  bigOne = bigOne + calc.uA; }
      if (type == 3) { calc.uA = calc.hit[b].C;  bigOne = bigOne + calc.uA; }
    }
    if (calc.hit[b].pixel == neighbourhood[pixel][2]) {
      if (type == 1) { calc.dA = calc.hit[b].A;  bigOne = bigOne + calc.dA; }
      if (type == 2) { calc.dA = calc.hit[b].B;  bigOne = bigOne + calc.dA; }
      if (type == 3) { calc.dA = calc.hit[b].C;  bigOne = bigOne + calc.dA; }
    }
    if (calc.hit[b].pixel == neighbourhood[pixel][3]) {
      if (type == 1) { calc.lA = calc.hit[b].A;  bigOne = bigOne + calc.lA; }
      if (type == 2) { calc.lA = calc.hit[b].B;  bigOne = bigOne + calc.lA; }
      if (type == 3) { calc.lA = calc.hit[b].C;  bigOne = bigOne + calc.lA; }
    }
    if (calc.hit[b].pixel == neighbourhood[pixel][4]) {
      if (type == 1) { calc.rA = calc.hit[b].A;  bigOne = bigOne + calc.rA; }
      if (type == 2) { calc.rA = calc.hit[b].B;  bigOne = bigOne + calc.rA; }
      if (type == 3) { calc.rA = calc.hit[b].C;  bigOne = bigOne + calc.rA; }
    }
    
    if (maxCheck == 8) {
      if (calc.hit[b].pixel == neighbourhood[pixel][5]) {
	if (type == 1) { calc.ulA = calc.hit[b].A;  bigOne = bigOne + calc.ulA; }
	if (type == 2) { calc.ulA = calc.hit[b].B;  bigOne = bigOne + calc.ulA; }
	if (type == 3) { calc.ulA = calc.hit[b].C;  bigOne = bigOne + calc.ulA; }
      }
      if (calc.hit[b].pixel == neighbourhood[pixel][6]) {
	if (type == 1) { calc.urA = calc.hit[b].A;  bigOne = bigOne + calc.urA; }
	if (type == 2) { calc.urA = calc.hit[b].B;  bigOne = bigOne + calc.urA; }
	if (type == 3) { calc.urA = calc.hit[b].C;  bigOne = bigOne + calc.urA; }
      }
      if (calc.hit[b].pixel == neighbourhood[pixel][7]) {
	if (type == 1) { calc.dlA = calc.hit[b].A;  bigOne = bigOne + calc.dlA; }
	if (type == 2) { calc.dlA = calc.hit[b].B;  bigOne = bigOne + calc.dlA; }
	if (type == 3) { calc.dlA = calc.hit[b].C;  bigOne = bigOne + calc.dlA; }
      }
      if (calc.hit[b].pixel == neighbourhood[pixel][8]) {
	if (type == 1) { calc.drA = calc.hit[b].A;  bigOne = bigOne + calc.drA; }
	if (type == 2) { calc.drA = calc.hit[b].B;  bigOne = bigOne + calc.drA; }
	if (type == 3) { calc.drA = calc.hit[b].C;  bigOne = bigOne + calc.drA; }
      }
    } /* maxCheck = 4 means don't check corners, 8 checks them */
  } /* End nHits loop */

  bigOne = bigOne/norm;  
  
  if (calc.uA == 0 && bigZero == 0) { calc.uF = 0; } else { calc.uF = calc.uA/(calc.uA + bigZero); }
  if (calc.dA == 0 && bigZero == 0) { calc.dF = 0; } else { calc.dF = calc.dA/(calc.dA + bigZero); }
  if (calc.lA == 0 && bigZero == 0) { calc.lF = 0; } else { calc.lF = calc.lA/(calc.lA + bigZero); }
  if (calc.rA == 0 && bigZero == 0) { calc.rF = 0; } else { calc.rF = calc.rA/(calc.rA + bigZero); }
  
  if (showSetup) {
    cout << endl;
    cout << "PhosWall:  Type: " << type << ". Neighbour check performed on pixel "
	 << pixel << "." << " Up: " << setw(3) << neighbourhood[pixel][1]
	 << " Down: " << setw(3) << neighbourhood[pixel][2] << " Left: "
	 << setw(3) << neighbourhood[pixel][3] << " Right: " << setw(3)
	 << neighbourhood[pixel][4] << endl;
    cout << " BigZero: " << bigZero << " BigOne: " << bigOne << endl;
  }

  Float_t pConst = 0.09;  Float_t ratio = 0.;
  if ((calc.uF == 0) && (calc.dF != 0)) {
    ratio = pConst/calc.dF;
    calc.uF = ratio/(ratio + 1);
    bigOne = bigOne + bigZero + calc.uF;
  }
  if ((calc.uF != 0) && (calc.dF == 0)) {
    ratio = pConst/calc.uF;
    calc.dF = ratio/(ratio + 1);
    bigOne = bigOne + bigZero * calc.dF;
  }
  if ((calc.rF == 0) && (calc.lF != 0)) {
    ratio = pConst/calc.lF;
    calc.rF = ratio/(ratio + 1);
    bigOne = bigOne + bigZero * calc.rF;
  }
  if ((calc.rF != 0) && (calc.lF == 0)) {
    ratio = pConst/calc.rF;
    calc.lF = ratio/(ratio + 1);
    bigOne = bigOne + bigZero * calc.lF;
  }

  if (doPrinting) {
    cout << "uA: " << setw(8) << calc.uA << " dA: " << setw(8) << calc.dA 
	 << " lA: " << setw(8) << calc.lA << " rA: " << setw(8) << calc.rA;
    cout << " (Sum: " << calc.uA + calc.dA + calc.lA + calc.rA << ")" << endl;
    cout << "uF: " << setw(8) << calc.uF << " dF: " << setw(8) << calc.dF 
	 << " lF: " << setw(8) << calc.lF << " rF: " << setw(8) << calc.rF;
    cout << endl;
  }

  calc.ulF = calc.ulA / (calc.ulA + bigZero);
  calc.urF = calc.urA / (calc.urA + bigZero);
  calc.dlF = calc.dlA / (calc.dlA + bigZero);
  calc.drF = calc.drA / (calc.drA + bigZero);

  return input;
}

/**********************************************************************/

/* Average distance = 5.5 cm;  phi1 = phi;  phi2 = phi + 270; 
                               phi3 = phi + 180;  phi4 = phi +  90 */

void phosWallFull::SetPwallPositions() {
   double detoffset, phi;
   for(int i=0; i<4; i++) {
     
      if(i==0)      {detoffset = 0;    phi = TMath::DegToRad()*0;}
      else if(i==1) {detoffset = 64;   phi = TMath::DegToRad()*270;}
      else if(i==2) {detoffset = 64*2; phi = TMath::DegToRad()*180;}
      else if(i==3) {detoffset = 64*3; phi = TMath::DegToRad()*90;}
    
      pWallPositions[0 + detoffset] = new TVector3;  
      pWallPositions[0 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*37.141,TMath::DegToRad()*20.966 +  phi);
      pWallPositions[1 + detoffset] = new TVector3;  
      pWallPositions[1 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*32.895,TMath::DegToRad()*24.772 +  phi);
      pWallPositions[2 + detoffset] = new TVector3;  
      pWallPositions[2 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*28.445,TMath::DegToRad()*30.019 +  phi);
      pWallPositions[3 + detoffset] = new TVector3;  
      pWallPositions[3 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*23.908,TMath::DegToRad()*37.687 +  phi);
      pWallPositions[4 + detoffset] = new TVector3;  
      pWallPositions[4 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*19.653,TMath::DegToRad()*49.364 +  phi);
      pWallPositions[5 + detoffset] = new TVector3;  
      pWallPositions[5 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*16.391,TMath::DegToRad()*67.124 +  phi);   
      pWallPositions[6 + detoffset] = new TVector3;  
      pWallPositions[6 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*15.165,TMath::DegToRad()*-89.178 + phi);  
      pWallPositions[7 + detoffset] = new TVector3;  
      pWallPositions[7 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*16.580,TMath::DegToRad()*-65.538 + phi);  
                                                                
      pWallPositions[8 + detoffset] = new TVector3;  
      pWallPositions[8 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*40.971,TMath::DegToRad()*26.144 +  phi);   
      pWallPositions[9 + detoffset] = new TVector3;  
      pWallPositions[9 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*36.941,TMath::DegToRad()*30.589 +  phi);   
      pWallPositions[10 + detoffset] = new TVector3;  
      pWallPositions[10 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*32.748,TMath::DegToRad()*36.507 +  phi);
      pWallPositions[11 + detoffset] = new TVector3;  
      pWallPositions[11 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*28.546,TMath::DegToRad()*44.701 +  phi);
      pWallPositions[12 + detoffset] = new TVector3;  
      pWallPositions[12 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*24.730,TMath::DegToRad()*56.18  +  phi);
      pWallPositions[13 + detoffset] = new TVector3;  
      pWallPositions[13 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*21.943,TMath::DegToRad()*71.77 +   phi);
      pWallPositions[14 + detoffset] = new TVector3;  
      pWallPositions[14 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*20.940,TMath::DegToRad()*-89.358 + phi);
      pWallPositions[15 + detoffset] = new TVector3;  
      pWallPositions[15 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*22.100,TMath::DegToRad()*-70.448 + phi);
                                                                 
      pWallPositions[16 + detoffset] = new TVector3;  
      pWallPositions[16 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*45.272,TMath::DegToRad()*30.836  + phi);
      pWallPositions[17 + detoffset] = new TVector3;  
      pWallPositions[17 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*41.556,TMath::DegToRad()*35.713  + phi);
      pWallPositions[18 + detoffset] = new TVector3;  
      pWallPositions[18 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*37.721,TMath::DegToRad()*41.992  + phi);
      pWallPositions[19 + detoffset] = new TVector3;  
      pWallPositions[19 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*33.932,TMath::DegToRad()*50.277  + phi);
      pWallPositions[20 + detoffset] = new TVector3;  
      pWallPositions[20 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*30.574,TMath::DegToRad()*61.151  + phi);
      pWallPositions[21 + detoffset] = new TVector3;  
      pWallPositions[21 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*28.196,TMath::DegToRad()*74.846  + phi);
      pWallPositions[22 + detoffset] = new TVector3;  
      pWallPositions[22 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*27.361,TMath::DegToRad()*-89.472 + phi);
      pWallPositions[23 + detoffset] = new TVector3;  
      pWallPositions[23 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*28.328,TMath::DegToRad()*-73.722 + phi);
      
      pWallPositions[24 + detoffset] = new TVector3;  
      pWallPositions[24 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*50.035,TMath::DegToRad()*35.110  + phi);
      pWallPositions[25 + detoffset] = new TVector3;  
      pWallPositions[25 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*46.724,TMath::DegToRad()*40.255  + phi);
      pWallPositions[26 + detoffset] = new TVector3;  
      pWallPositions[26 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*43.330,TMath::DegToRad()*46.672  + phi);
      pWallPositions[27 + detoffset] = new TVector3;  
      pWallPositions[27 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*40.027,TMath::DegToRad()*54.797  + phi);
      pWallPositions[28 + detoffset] = new TVector3;  
      pWallPositions[28 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*37.149,TMath::DegToRad()*64.933  + phi);
      pWallPositions[29 + detoffset] = new TVector3;  
      pWallPositions[29 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*35.153,TMath::DegToRad()*77.050  + phi);
      pWallPositions[30 + detoffset] = new TVector3;  
      pWallPositions[30 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*34.462,TMath::DegToRad()*-89.552 + phi);
      pWallPositions[31 + detoffset] = new TVector3;  
      pWallPositions[31 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*35.263,TMath::DegToRad()*-76.075 + phi);
      
      pWallPositions[32 + detoffset] = new TVector3;  
      pWallPositions[32 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*55.164,TMath::DegToRad()*38.979  + phi);
      pWallPositions[33 + detoffset] = new TVector3;  
      pWallPositions[33 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*52.326,TMath::DegToRad()*44.260  + phi);
      pWallPositions[34 + detoffset] = new TVector3;  
      pWallPositions[34 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*49.446,TMath::DegToRad()*50.663  + phi);
      pWallPositions[35 + detoffset] = new TVector3;  
      pWallPositions[35 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*46.672,TMath::DegToRad()*58.492  + phi);
      pWallPositions[36 + detoffset] = new TVector3;  
      pWallPositions[36 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*44.293,TMath::DegToRad()*67.884  + phi);
      pWallPositions[37 + detoffset] = new TVector3;  
      pWallPositions[37 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*42.667,TMath::DegToRad()*78.701  + phi);
      pWallPositions[38 + detoffset] = new TVector3;  
      pWallPositions[38 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*42.109,TMath::DegToRad()*-89.611 + phi);
      pWallPositions[39 + detoffset] = new TVector3;  
      pWallPositions[39 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*42.756,TMath::DegToRad()*-77.843 + phi);
      
      pWallPositions[40 + detoffset] = new TVector3;  
      pWallPositions[40 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*60.528,TMath::DegToRad()*42.468  + phi);
      pWallPositions[41 + detoffset] = new TVector3;  
      pWallPositions[41 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*58.205,TMath::DegToRad()*47.786  + phi);
      pWallPositions[42 + detoffset] = new TVector3;  
      pWallPositions[42 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*55.871,TMath::DegToRad()*54.074  + phi);
      pWallPositions[43 + detoffset] = new TVector3;  
      pWallPositions[43 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*53.651,TMath::DegToRad()*61.546  + phi);
      pWallPositions[44 + detoffset] = new TVector3;  
      pWallPositions[44 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*51.772,TMath::DegToRad()*70.238  + phi);
      pWallPositions[45 + detoffset] = new TVector3;  
      pWallPositions[45 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*50.503,TMath::DegToRad()*79.983  + phi);
      pWallPositions[46 + detoffset] = new TVector3;  
      pWallPositions[46 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*50.071,TMath::DegToRad()*-89.656 + phi);
      pWallPositions[47 + detoffset] = new TVector3;  
      pWallPositions[47 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*50.572,TMath::DegToRad()*-79.217 + phi);
      
      pWallPositions[48 + detoffset] = new TVector3;  
      pWallPositions[48 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*65.983,TMath::DegToRad()*45.606  + phi);
      pWallPositions[49 + detoffset] = new TVector3;  
      pWallPositions[49 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*64.180,TMath::DegToRad()*50.890  + phi);
      pWallPositions[50 + detoffset] = new TVector3;  
      pWallPositions[50 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*62.390,TMath::DegToRad()*57.005  + phi);
      pWallPositions[51 + detoffset] = new TVector3;  
      pWallPositions[51 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*60.711,TMath::DegToRad()*64.098  + phi);
      pWallPositions[52 + detoffset] = new TVector3;  
      pWallPositions[52 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*59.308,TMath::DegToRad()*72.154  + phi);
      pWallPositions[53 + detoffset] = new TVector3;  
      pWallPositions[53 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*58.371,TMath::DegToRad()*81.006  + phi);
      pWallPositions[54 + detoffset] = new TVector3;  
      pWallPositions[54 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*58.054,TMath::DegToRad()*-89.692 + phi);
      pWallPositions[55 + detoffset] = new TVector3;  
      pWallPositions[55 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*58.422,TMath::DegToRad()*-80.315 + phi);
      
      pWallPositions[56 + detoffset] = new TVector3;  
      pWallPositions[56 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*71.343,TMath::DegToRad()*48.456  + phi);
      pWallPositions[57 + detoffset] = new TVector3;  
      pWallPositions[57 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*70.127,TMath::DegToRad()*53.656  + phi);
      pWallPositions[58 + detoffset] = new TVector3;  
      pWallPositions[58 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*68.848,TMath::DegToRad()*59.560  + phi);
      pWallPositions[59 + detoffset] = new TVector3;  
      pWallPositions[59 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*67.665,TMath::DegToRad()*66.274  + phi);
      pWallPositions[60 + detoffset] = new TVector3;  
      pWallPositions[60 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*66.691,TMath::DegToRad()*73.755  + phi);
      pWallPositions[61 + detoffset] = new TVector3;  
      pWallPositions[61 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*66.047,TMath::DegToRad()*81.484  + phi);
      pWallPositions[62 + detoffset] = new TVector3;  
      pWallPositions[62 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*65.831,TMath::DegToRad()*-89.721 + phi);
      pWallPositions[63 + detoffset] = new TVector3;  
      pWallPositions[63 + detoffset]->SetMagThetaPhi(43.409,TMath::DegToRad()*66.082,TMath::DegToRad()*-81.219 + phi);
   }
}

/**********************************************************************/

Int_t phosWallFull::DetermineNParticles(Float_t radius, Int_t biggestHit) {

  Int_t nPixels = raw.hit.size();

  Bool_t used[256] = {0}; Bool_t checked[256] = {0};
  for (Int_t i=0; i<nPixels; i++) { 
    if (0) {
      cout << i << ": " << ( ( (calc.hit[i].Vec.X() - calc.hit[biggestHit].Vec.X())*
			       (calc.hit[i].Vec.X() - calc.hit[biggestHit].Vec.X())) + 
			     ( (calc.hit[i].Vec.Y() - calc.hit[biggestHit].Vec.Y())*
			       (calc.hit[i].Vec.Y() - calc.hit[biggestHit].Vec.Y())) ); 
    }
    
    if ( (calc.hit[i].Vec - calc.hit[biggestHit].Vec).Mag() < radius ) {
      /* We look within a sphere of some radius to decide if pixels are
	 related to the biggest hit we know of to this point. */
      used[i] = 1;
      if (0) { cout << " -- Clumped with biggest! (" << calc.hit[i].B << ")" <<  endl; }
    } else {
      if (0) { cout << " -- Possible second highest (" << calc.hit[i].B << ")" << endl; }
    }
  }

  /* As a check use another method...
     calculate the number of regions of pixels, based on the assumption that
     events will not have neighbouring pixels, and will not have pixels light up 
     with a gap in between, which is kind of an over-simplification. */
  for (Int_t i=0; i<256; i++) { used[i] = 0; }
  for (Int_t i=0; i<nPixels; i++) {
    if (i != biggestHit) {
      if (isNeighbour(calc.hit[biggestHit].pixel, calc.hit[i].pixel)) {
	used[i] = 1;
      }
    }
    if (i == biggestHit) { used[i] = 1; checked[i] = 1; }
  }
  
  Int_t done = 0; Int_t loop = 0; 
  while (!done) {
    Int_t i;
    for (i=0; i<nPixels; i++) {
      if (used[i] == 1 && checked[i] == 0) {
	for (Int_t j=0; j<nPixels; j++) {
	  if (isNeighbour(calc.hit[i].pixel, calc.hit[j].pixel)) {
	    used[i] = 1;
	  }
	}
	checked[i] = 1;
        i = nPixels;
      }
    }
    if (i == nPixels) { done = 1; }
  }

  Float_t secondBiggestB = 0;
  Int_t secondBiggestHit = -1;

  for (Int_t i=0; i<nPixels; i++) {
    if ((calc.hit[i].B > calc.secondBiggestB("B")) && used[i] == 0) {
      secondBiggestB = calc.hit[i].B;
      secondBiggestHit = i;
    }
  }

  if (0) { 
    if (secondBiggestB > 100  && calc.secondBiggestB("B") > 100) {
      cout << "PhosWall ----> " << eventNum << ": Two hits (both methods)" << endl;
    } else if (secondBiggestB > 100 && calc.secondBiggestB("B") < 100) {
      cout << "PhosWall ----> " << eventNum << ": Two hits (from pixel analysis only)" << endl;
    } else if (secondBiggestB < 100 && calc.secondBiggestB("B") > 100) {
      cout << "PhosWall ----> " << eventNum << ": Two hits (from radius analysis only)" << endl;
    }
  }

  if (calc.secondBiggestB("B") > 100) { return 2; } else { return 1; }

}


