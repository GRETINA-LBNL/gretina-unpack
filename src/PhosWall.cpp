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

TVector3 phosWallCalc::biggestPixelPosition(TString var = "B") {
  return hit[biggestHit(var)].Vec;
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
    char * str = fgets(line, 300, phosIn);
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
  for (Int_t i=0; i<256; i++) {
    pWallPositions[i] = new TVector3();
  }

pWallPositions[0]->SetXYZ(21.3698773964654, 36.6405948781801, 53.6457281229654);
pWallPositions[1]->SetXYZ(15.2003923142873, 36.6409081567024, 53.6464069528786);
pWallPositions[2]->SetXYZ(9.11941242843985, 36.6386155574124, 53.642393797128);
pWallPositions[3]->SetXYZ(3.04018226933476, 36.6416913438391, 53.6454849026649);
pWallPositions[4]->SetXYZ(-3.04018226933475, 36.6416913438391, 53.6454849026649);
pWallPositions[5]->SetXYZ(-9.11941242843984, 36.6386155574124, 53.642393797128);
pWallPositions[6]->SetXYZ(-15.2003923142872, 36.6409081567024, 53.6464069528786);
pWallPositions[7]->SetXYZ(-21.3698773964653, 36.6405948781801, 53.6457281229654);
pWallPositions[8]->SetXYZ(21.3538393676626, 40.564634135286, 48.8764187294901);
pWallPositions[9]->SetXYZ(15.1988659374766, 40.5562406385103, 48.8764025121107);
pWallPositions[10]->SetXYZ(9.1186828545775, 40.5584944322787, 48.8785028619377);
pWallPositions[11]->SetXYZ(3.04018734931569, 40.5581676516405, 48.8784090025697);
pWallPositions[12]->SetXYZ(-3.04018734931568, 40.5581676516405, 48.8784090025697);
pWallPositions[13]->SetXYZ(-9.1186828545775, 40.5584944322787, 48.8785028619377);
pWallPositions[14]->SetXYZ(-15.1988659374766, 40.5562406385103, 48.8764025121107);
pWallPositions[15]->SetXYZ(-21.3538393676626, 40.564634135286, 48.8764187294901);
pWallPositions[16]->SetXYZ(21.3704754068713, 44.4164705920435, 44.1789696008231);
pWallPositions[17]->SetXYZ(15.1998917095267, 44.4154863611401, 44.1778275514091);
pWallPositions[18]->SetXYZ(9.11949937670316, 44.4149300533203, 44.1774424623809);
pWallPositions[19]->SetXYZ(3.03970183185058, 44.4165902367194, 44.1798948587116);
pWallPositions[20]->SetXYZ(-3.03970183185057, 44.4165902367194, 44.1798948587116);
pWallPositions[21]->SetXYZ(-9.11949937670315, 44.4149300533203, 44.1774424623809);
pWallPositions[22]->SetXYZ(-15.1998917095267, 44.4154863611401, 44.1778275514091);
pWallPositions[23]->SetXYZ(-21.3704754068713, 44.4164705920435, 44.1789696008231);
pWallPositions[24]->SetXYZ(21.3694349957309, 48.2751801909771, 39.4801591929658);
pWallPositions[25]->SetXYZ(15.1992501799208, 48.2734240247271, 39.4794920926939);
pWallPositions[26]->SetXYZ(9.1196663089296, 48.2737064645332, 39.4798275785224);
pWallPositions[27]->SetXYZ(3.03985446842056, 48.2767663174297, 39.4824862279861);
pWallPositions[28]->SetXYZ(-3.03985446842056, 48.2767663174297, 39.4824862279861);
pWallPositions[29]->SetXYZ(-9.11966630892959, 48.2737064645332, 39.4798275785224);
pWallPositions[30]->SetXYZ(-15.1992501799208, 48.2734240247271, 39.4794920926939);
pWallPositions[31]->SetXYZ(-21.3694349957308, 48.2751801909771, 39.4801591929658);
pWallPositions[32]->SetXYZ(21.3696169704978, 52.1330168056971, 34.7817277098351);
pWallPositions[33]->SetXYZ(15.2001729055072, 52.1358820724433, 34.7838960096619);
pWallPositions[34]->SetXYZ(9.11948876816495, 52.1341146246288, 34.7828871746349);
pWallPositions[35]->SetXYZ(3.03993360088104, 52.1361646311879, 34.7834823623499);
pWallPositions[36]->SetXYZ(-3.03993360088103, 52.1361646311879, 34.7834823623499);
pWallPositions[37]->SetXYZ(-9.11948876816493, 52.1341146246288, 34.7828871746349);
pWallPositions[38]->SetXYZ(-15.2001729055072, 52.1358820724433, 34.7838960096619);
pWallPositions[39]->SetXYZ(-21.3696169704978, 52.1330168056971, 34.7817277098351);
pWallPositions[40]->SetXYZ(21.3685287140451, 55.9908977735959, 30.0826549370857);
pWallPositions[41]->SetXYZ(15.1994293623042, 55.9934028019465, 30.0848264638387);
pWallPositions[42]->SetXYZ(9.1201734296604, 55.9938050572881, 30.0848018080448);
pWallPositions[43]->SetXYZ(3.04030231979183, 55.9927927943668, 30.0836807929377);
pWallPositions[44]->SetXYZ(-3.04030231979183, 55.9927927943668, 30.0836807929377);
pWallPositions[45]->SetXYZ(-9.12017342966039, 55.9938050572881, 30.0848018080448);
pWallPositions[46]->SetXYZ(-15.1994293623042, 55.9934028019466, 30.0848264638387);
pWallPositions[47]->SetXYZ(-21.3685287140451, 55.9908977735959, 30.0826549370857);
pWallPositions[48]->SetXYZ(21.3706198957972, 59.8536453966802, 25.3865002394544);
pWallPositions[49]->SetXYZ(15.1992607055909, 59.8516106128116, 25.3852816493311);
pWallPositions[50]->SetXYZ(9.12012553947569, 59.8517738392358, 25.3851519412663);
pWallPositions[51]->SetXYZ(3.04048249165567, 59.854609344751, 25.3864676677396);
pWallPositions[52]->SetXYZ(-3.04048249165567, 59.854609344751, 25.3864676677396);
pWallPositions[53]->SetXYZ(-9.12012553947569, 59.8517738392358, 25.3851519412663);
pWallPositions[54]->SetXYZ(-15.1992607055909, 59.8516106128116, 25.3852816493311);
pWallPositions[55]->SetXYZ(-21.3706198957972, 59.8536453966802, 25.3865002394544);
pWallPositions[56]->SetXYZ(21.3709136930911, 63.7710434253813, 20.6189231544857);
pWallPositions[57]->SetXYZ(15.1994851588037, 63.7662025462, 20.6171181862136);
pWallPositions[58]->SetXYZ(9.12009319251732, 63.7699148125208, 20.6181323444593);
pWallPositions[59]->SetXYZ(3.03955969028272, 63.7677631081098, 20.6185806717435);
pWallPositions[60]->SetXYZ(-3.03955969028273, 63.7677631081098, 20.6185806717435);
pWallPositions[61]->SetXYZ(-9.12009319251731, 63.7699148125208, 20.6181323444593);
pWallPositions[62]->SetXYZ(-15.1994851588037, 63.7662025462, 20.6171181862136);
pWallPositions[63]->SetXYZ(-21.3709136930911, 63.7710434253813, 20.6189231544857);
pWallPositions[64]->SetXYZ(38.2723612979321, -18.2874605576452, 53.6457281229654);
pWallPositions[65]->SetXYZ(37.7658244073703, -12.1388048788921, 53.6464069528786);
pWallPositions[66]->SetXYZ(37.2639615214091, -6.07856929999032, 53.642393797128);
pWallPositions[67]->SetXYZ(36.7675926214365, -0.0196364976748579, 53.6454849026649);
pWallPositions[68]->SetXYZ(36.2680651382871, 6.04017416081012, 53.6454849026649);
pWallPositions[69]->SetXYZ(35.765649923353, 12.0983518945351, 53.642393797128);
pWallPositions[70]->SetXYZ(35.268399047373, 18.1589676383708, 53.6464069528786);
pWallPositions[71]->SetXYZ(34.7612806431286, 24.3075753833354, 53.6457281229654);
pWallPositions[72]->SetXYZ(42.1818182601232, -17.9491004381691, 48.8764187294901);
pWallPositions[73]->SetXYZ(41.6677962195402, -11.8156226548166, 48.8764025121107);
pWallPositions[74]->SetXYZ(41.1705298188565, -5.75580768009093, 48.8785028619377);
pWallPositions[75]->SetXYZ(40.6708302086523, 0.302113416608436, 48.8784090025697);
pWallPositions[76]->SetXYZ(40.1713018908193, 6.36193420071081, 48.8784090025697);
pWallPositions[77]->SetXYZ(39.6723403389678, 12.4196320253322, 48.8785028619377);
pWallPositions[78]->SetXYZ(39.1706238909761, 18.4790801856, 48.8764025121107);
pWallPositions[79]->SetXYZ(38.6733749270927, 24.6139405625768, 48.8764187294901);
pWallPositions[80]->SetXYZ(46.0220008060376, -17.6492356995602, 44.1789696008231);
pWallPositions[81]->SetXYZ(45.5140805362618, -11.4995917152147, 44.1778275514091);
pWallPositions[82]->SetXYZ(45.0139963424122, -5.43979905955398, 44.1774424623809);
pWallPositions[83]->SetXYZ(44.5161700115515, 0.619582912801444, 44.1798948587116);
pWallPositions[84]->SetXYZ(44.0167214683162, 6.6784359444307, 44.1798948587116);
pWallPositions[85]->SetXYZ(43.5156749185047, 12.7372413382978, 44.1774424623809);
pWallPositions[86]->SetXYZ(43.0167418852264, 18.7971288833634, 44.1778275514091);
pWallPositions[87]->SetXYZ(42.5108263568841, 24.94693811722, 44.1789696008231);
pWallPositions[88]->SetXYZ(49.8675810657328, -17.331189606044, 39.4801591929658);
pWallPositions[89]->SetXYZ(49.3589242398305, -11.1820065725967, 39.4794920926939);
pWallPositions[90]->SetXYZ(48.8597423768646, -5.1229507392284, 39.4798275785224);
pWallPositions[91]->SetXYZ(48.3633098096527, 0.936560469186382, 39.4824862279861);
pWallPositions[92]->SetXYZ(47.8638361869468, 6.99571774201955, 39.4824862279861);
pWallPositions[93]->SetXYZ(47.361395738838, 13.0543955463025, 39.4798275785224);
pWallPositions[94]->SetXYZ(46.8616932056396, 19.1134084611677, 39.4794920926939);
pWallPositions[95]->SetXYZ(46.3565797706778, 25.2628835739879, 39.4801591929658);
pWallPositions[96]->SetXYZ(53.7123917163598, -17.0144334856385, 34.7817277098351);
pWallPositions[97]->SetXYZ(53.2084015570923, -10.8656090288693, 34.7838960096619);
pWallPositions[98]->SetXYZ(52.7070863443822, -4.80562505533417, 34.7828871746349);
pWallPositions[99]->SetXYZ(52.2096684312703, 1.25354738438824, 34.7834823623499);
pWallPositions[100]->SetXYZ(51.7101818064366, 7.31286238714833, 34.7834823623499);
pWallPositions[101]->SetXYZ(51.2087710905813, 13.3713404892351, 34.7828871746349);
pWallPositions[102]->SetXYZ(50.7110211329549, 19.4316183443899, 34.7838960096619);
pWallPositions[103]->SetXYZ(50.2013627351355, 25.5799755720024, 34.7817277098351);
pWallPositions[104]->SetXYZ(57.5571422154059, -16.6964077840646, 30.0826549370857);
pWallPositions[105]->SetXYZ(57.0528213551232, -10.5479564703839, 30.0848264638387);
pWallPositions[106]->SetXYZ(56.5537858441844, -4.48921762346064, 30.0848018080448);
pWallPositions[107]->SetXYZ(56.053290057185, 1.57001811254481, 30.0836807929377);
pWallPositions[108]->SetXYZ(55.5537428487371, 7.63006806031409, 30.0836807929377);
pWallPositions[109]->SetXYZ(55.0553603131232, 13.6890857610304, 30.0848018080448);
pWallPositions[110]->SetXYZ(54.5555653085722, 19.7477620033972, 30.0848264638387);
pWallPositions[111]->SetXYZ(54.0462942488587, 25.8958052754111, 30.0826549370857);
pWallPositions[112]->SetXYZ(61.4070041227454, -16.3811509568563, 25.3865002394544);
pWallPositions[113]->SetXYZ(60.8979731412662, -10.2308204084399, 25.3852816493311);
pWallPositions[114]->SetXYZ(60.3987093308938, -4.17222155703923, 25.3851519412663);
pWallPositions[115]->SetXYZ(59.902067041706, 1.8871029981887, 25.3864676677396);
pWallPositions[116]->SetXYZ(59.4024902295075, 7.94751207159068, 25.3864676677396);
pWallPositions[117]->SetXYZ(58.9002938808821, 14.0059595281632, 25.3851519412663);
pWallPositions[118]->SetXYZ(58.4007470180167, 20.064535047734, 25.3852816493311);
pWallPositions[119]->SetXYZ(57.895814787385, 26.2152034534732, 25.3865002394544);
pWallPositions[120]->SetXYZ(65.3111840344433, -16.0596130606717, 20.6189231544857);
pWallPositions[121]->SetXYZ(64.7993507463837, -9.9094439350979, 20.6171181862136);
pWallPositions[122]->SetXYZ(64.3036028816282, -3.85029758307284, 20.6181323444593);
pWallPositions[123]->SetXYZ(63.8019170865578, 2.20950469626389, 20.6185806717435);
pWallPositions[124]->SetXYZ(63.302491898376, 8.26807440573948, 20.6185806717435);
pWallPositions[125]->SetXYZ(62.8051949934015, 14.3277917655575, 20.6181323444593);
pWallPositions[126]->SetXYZ(62.3020899902453, 20.3863316735277, 20.6171181862136);
pWallPositions[127]->SetXYZ(61.7999486744709, 26.5372997019701, 20.6189231544857);
pWallPositions[128]->SetXYZ(-19.3297179396293, -37.7566796387242, 53.6457281229654);
pWallPositions[129]->SetXYZ(-13.1694875434478, -37.4188812619928, 53.6464069528786);
pWallPositions[130]->SetXYZ(-7.09777214491732, -37.0833313441859, 53.642393797128);
pWallPositions[131]->SetXYZ(-1.02750963563929, -36.7532376358413, 53.6454849026649);
pWallPositions[132]->SetXYZ(5.04371698362952, -36.4200105952346, 53.6454849026649);
pWallPositions[133]->SetXYZ(11.1138942753794, -36.083696970232, 53.642393797128);
pWallPositions[134]->SetXYZ(17.1858586494188, -35.7526829697904, 53.6464069528786);
pWallPositions[135]->SetXYZ(23.3460523451024, -35.4142159713669, 53.6457281229654);
pWallPositions[136]->SetXYZ(-19.0986517778354, -41.6739426805112, 48.8764187294901);
pWallPositions[137]->SetXYZ(-12.9533883890185, -41.3282459070958, 48.8764025121107);
pWallPositions[138]->SetXYZ(-6.8822194363705, -40.9972792175919, 48.8785028619377);
pWallPositions[139]->SetXYZ(-0.812876950437627, -40.6638283175719, 48.8784090025697);
pWallPositions[140]->SetXYZ(5.25835981352404, -40.3306007201608, 48.8784090025697);
pWallPositions[141]->SetXYZ(11.328017354156, -39.9977233157025, 48.8785028619377);
pWallPositions[142]->SetXYZ(17.3989369317561, -39.6622134298997, 48.8764025121107);
pWallPositions[143]->SetXYZ(23.5451180149463, -39.3332355131357, 48.8764187294901);
pWallPositions[144]->SetXYZ(-18.9041675669329, -45.5209020953108, 44.1789696008231);
pWallPositions[145]->SetXYZ(-12.7429113150848, -45.1817479508468, 44.1778275514091);
pWallPositions[146]->SetXYZ(-6.67168743121123, -44.8479639152464, 44.1774424623809);
pWallPositions[147]->SetXYZ(-0.600935969300934, -44.5164256364321, 44.1798948587116);
pWallPositions[148]->SetXYZ(5.4693312190573, -44.1832512554185, 44.1798948587116);
pWallPositions[149]->SetXYZ(11.540206831555, -43.8483170350694, 44.1774424623809);
pWallPositions[150]->SetXYZ(17.6114893683416, -43.5156015574501, 44.1778275514091);
pWallPositions[151]->SetXYZ(23.7728511451926, -43.1783699020478, 44.1789696008231);
pWallPositions[152]->SetXYZ(-18.6916567965989, -49.3737555864068, 39.4801591929658);
pWallPositions[153]->SetXYZ(-12.5308411318083, -49.0338525269258, 39.4794920926939);
pWallPositions[154]->SetXYZ(-6.46037852821314, -48.7009502852085, 39.4798275785224);
pWallPositions[155]->SetXYZ(-0.389536084933133, -48.370808788932, 39.4824862279861);
pWallPositions[156]->SetXYZ(5.68103591778321, -48.0376176777923, 39.4824862279861);
pWallPositions[157]->SetXYZ(11.7518759960418, -47.7012836303574, 39.4798275785224);
pWallPositions[158]->SetXYZ(17.8223053139098, -47.3677749783617, 39.4794920926939);
pWallPositions[159]->SetXYZ(23.983311119207, -47.0313359609038, 39.4801591929658);
pWallPositions[160]->SetXYZ(-18.4804144179472, -53.2258043966441, 34.7817277098351);
pWallPositions[161]->SetXYZ(-12.3200851186752, -52.8905564207603, 34.7838960096619);
pWallPositions[162]->SetXYZ(-6.24863624388938, -52.5555470733411, 34.7828871746349);
pWallPositions[163]->SetXYZ(-0.178105431501133, -52.2244113150485, 34.7834823623499);
pWallPositions[164]->SetXYZ(5.89262459828664, -51.8912115303908, 34.7834823623499);
pWallPositions[165]->SetXYZ(11.9632906427658, -51.555898402446, 34.7828871746349);
pWallPositions[166]->SetXYZ(18.0349309289274, -51.224376250159, 34.7838960096619);
pWallPositions[167]->SetXYZ(24.1949437926579, -50.8833633479452, 34.7817277098351);
pWallPositions[168]->SetXYZ(-18.2679012863295, -57.0778278799146, 30.0826549370857);
pWallPositions[169]->SetXYZ(-12.107935924683, -56.7422390986158, 30.0848264638387);
pWallPositions[170]->SetXYZ(-6.03779420024486, -56.4094744646337, 30.0848018080448);
pWallPositions[171]->SetXYZ(0.0328842559282171, -56.0752637241622, 30.0836807929377);
pWallPositions[172]->SetXYZ(6.10435061527376, -55.7420235251147, 30.0836807929377);
pWallPositions[173]->SetXYZ(12.1755268564642, -55.4097492680502, 30.0848018080448);
pWallPositions[174]->SetXYZ(18.2456221617738, -55.0761389551819, 30.0848264638387);
pWallPositions[175]->SetXYZ(24.4053105753829, -54.7355046438355, 30.0826549370857);
pWallPositions[176]->SetXYZ(-18.0582961035108, -60.9348849497569, 25.3865002394544);
pWallPositions[177]->SetXYZ(-11.8963230985432, -60.5946393312082, 25.3852816493311);
pWallPositions[178]->SetXYZ(-5.82631505883025, -60.2616426459989, 25.3851519412663);
pWallPositions[179]->SetXYZ(0.244346550344572, -59.9312863899829, 25.3864676677396);
pWallPositions[180]->SetXYZ(6.31617271187264, -59.5980264427322, 25.3864676677396);
pWallPositions[181]->SetXYZ(12.3869372546923, -59.2619212227138, 25.3851519412663);
pWallPositions[182]->SetXYZ(18.456925076372, -58.9285561987508, 25.3852816493311);
pWallPositions[183]->SetXYZ(24.6191187619436, -58.5923310117245, 25.3865002394544);
pWallPositions[184]->SetXYZ(-17.8439011879196, -64.8464117899763, 20.6189231544857);
pWallPositions[185]->SetXYZ(-11.6820127277525, -64.5033604929567, 20.6171181862136);
pWallPositions[186]->SetXYZ(-5.61155377295372, -64.1738934403205, 20.6181323444593);
pWallPositions[187]->SetXYZ(0.459723634388469, -63.8385086691629, 20.6185806717435);
pWallPositions[188]->SetXYZ(6.52970696684797, -63.5053498679434, 20.6185806717435);
pWallPositions[189]->SetXYZ(12.6016612565602, -63.1741740635465, 20.6181323444593);
pWallPositions[190]->SetXYZ(18.6717109667907, -62.8372512593313, 20.6171181862136);
pWallPositions[191]->SetXYZ(24.8341276964305, -62.5038241485815, 20.6189231544857);
pWallPositions[192]->SetXYZ(-37.2126211995102, 20.3574476894488, 53.6457281229654);
pWallPositions[193]->SetXYZ(-37.0438152156662, 14.1902724189396, 53.6464069528786);
pWallPositions[194]->SetXYZ(-36.8748304551968, 8.11164051868269, 53.642393797128);
pWallPositions[195]->SetXYZ(-36.7112600269412, 2.03461052839669, 53.6454849026649);
pWallPositions[196]->SetXYZ(-36.544583872256, -4.04346910110546, 53.6454849026649);
pWallPositions[197]->SetXYZ(-36.3748288559759, -10.1204574354277, 53.642393797128);
pWallPositions[198]->SetXYZ(-36.2104063521723, -16.1992144521487, 53.6464069528786);
pWallPositions[199]->SetXYZ(-36.0409525192926, -22.3663719558844, 53.6457281229654);
pWallPositions[200]->SetXYZ(-41.1347462243086, 20.2338491476851, 48.8764187294901);
pWallPositions[201]->SetXYZ(-40.9576345333144, 14.081418747754, 48.8764025121107);
pWallPositions[202]->SetXYZ(-40.7932162995574, 8.00345872439357, 48.8785028619377);
pWallPositions[203]->SetXYZ(-40.6262647213465, 1.92725638372004, 48.8784090025697);
pWallPositions[204]->SetXYZ(-40.4595882881545, -4.15083340192598, 48.8784090025697);
pWallPositions[205]->SetXYZ(-40.2932543262243, -10.2271943082036, 48.8785028619377);
pWallPositions[206]->SetXYZ(-40.1243089831634, -16.3050301789899, 48.8764025121107);
pWallPositions[207]->SetXYZ(-39.9639564997831, -22.4579201876084, 48.8764187294901);
pWallPositions[208]->SetXYZ(-44.9855912494918, 20.1448916329476, 44.1789696008231);
pWallPositions[209]->SetXYZ(-44.8154581283505, 13.9766537276998, 44.1778275514091);
pWallPositions[210]->SetXYZ(-44.6482251129988, 7.89856156410809, 44.1774424623809);
pWallPositions[211]->SetXYZ(-44.4832240603875, 1.82100320607964, 44.1798948587116);
pWallPositions[212]->SetXYZ(-44.316574245397, -4.25611590953644, 44.1798948587116);
pWallPositions[213]->SetXYZ(-44.148218002452, -10.3337373554294, 44.1774424623809);
pWallPositions[214]->SetXYZ(-43.9820759679862, -16.4118594383529, 44.1778275514091);
pWallPositions[215]->SetXYZ(-43.813889037534, -22.5801507164193, 44.1789696008231);
pWallPositions[216]->SetXYZ(-48.842822283677, 20.038075902472, 39.4801591929658);
pWallPositions[217]->SetXYZ(-48.6719284514901, 13.8702578893093, 39.4794920926939);
pWallPositions[218]->SetXYZ(-48.5055560302912, 7.79295089182663, 39.4798275785224);
pWallPositions[219]->SetXYZ(-48.3419537293273, 1.7153398754968, 39.4824862279861);
pWallPositions[220]->SetXYZ(-48.1752955461291, -4.36208439854221, 39.4824862279861);
pWallPositions[221]->SetXYZ(-48.0055393979468, -10.439695231295, 39.4798275785224);
pWallPositions[222]->SetXYZ(-47.8385810956606, -16.5169861614287, 39.4794920926939);
pWallPositions[223]->SetXYZ(-47.6711767460108, -22.6848998710303, 39.4801591929658);
pWallPositions[224]->SetXYZ(-52.6992141699869, 19.9325060289655, 34.7817277098351);
pWallPositions[225]->SetXYZ(-52.5329603397394, 13.7653018046244, 34.7838960096619);
pWallPositions[226]->SetXYZ(-52.3645086405179, 7.68695114617879, 34.7828871746349);
pWallPositions[227]->SetXYZ(-52.1999039086789, 1.60962438872765, 34.7834823623499);
pWallPositions[228]->SetXYZ(-52.0332413870913, -4.46795809075879, 34.7834823623499);
pWallPositions[229]->SetXYZ(-51.8645013729451, -10.5453534991668, 34.7828871746349);
pWallPositions[230]->SetXYZ(-51.6995620232087, -16.6238004811603, 34.7838960096619);
pWallPositions[231]->SetXYZ(-51.5275582858624, -22.7908470186158, 34.7817277098351);
pWallPositions[232]->SetXYZ(-56.555615572979, 19.8256651857549, 30.0826549370857);
pWallPositions[233]->SetXYZ(-56.3890110890559, 13.6588154200516, 30.0848264638387);
pWallPositions[234]->SetXYZ(-56.2227674278784, 7.58183295327487, 30.0848018080448);
pWallPositions[235]->SetXYZ(-56.0550929166329, 1.50427431550184, 30.0836807929377);
pWallPositions[236]->SetXYZ(-55.8884101802541, -4.57404532468814, 30.0836807929377);
pWallPositions[237]->SetXYZ(-55.7227222524198, -10.6518539682229, 30.0848018080448);
pWallPositions[238]->SetXYZ(-55.5556531702896, -16.7288137984498, 30.0848264638387);
pWallPositions[239]->SetXYZ(-55.3840189864827, -22.8955256284095, 30.0826549370857);
pWallPositions[240]->SetXYZ(-60.4169689577297, 19.7218691802503, 25.3865002394544);
pWallPositions[241]->SetXYZ(-60.2457644204415, 13.5528848715248, 25.3852816493311);
pWallPositions[242]->SetXYZ(-60.079285130646, 7.47602967820232, 25.3851519412663);
pWallPositions[243]->SetXYZ(-59.9154631935857, 1.39859354098656, 25.3864676677396);
pWallPositions[244]->SetXYZ(-59.7487705793934, -4.68008630751941, 25.3864676677396);
pWallPositions[245]->SetXYZ(-59.5792422117255, -10.7575749607138, 25.3851519412663);
pWallPositions[246]->SetXYZ(-59.4124153795636, -16.8344206229477, 25.3852816493311);
pWallPositions[247]->SetXYZ(-59.2452573465246, -23.0035159040145, 25.3865002394544);
pWallPositions[248]->SetXYZ(-64.3329029410421, 19.6147783779954, 20.6189231544857);
pWallPositions[249]->SetXYZ(-64.1588914620558, 13.4458016724868, 20.6171181862136);
pWallPositions[250]->SetXYZ(-63.9959528390243, 7.36859248851148, 20.6181323444593);
pWallPositions[251]->SetXYZ(-63.8271211568435, 1.29040294188486, 20.6185806717435);
pWallPositions[252]->SetXYZ(-63.6604791346785, -4.78643199742483, 20.6185806717435);
pWallPositions[253]->SetXYZ(-63.4959492441637, -10.8647395099447, 20.6181323444593);
pWallPositions[254]->SetXYZ(-63.3255297402851, -16.9419662190657, 20.6171181862136);
pWallPositions[255]->SetXYZ(-63.1611748466864, -23.1112077491793, 20.6189231544857);
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


