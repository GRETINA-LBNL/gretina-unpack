/*! \file PhosWall.h
    \brief Parameter and function definitions for phoswich wall analysis.

    This file provides the basic data structures and function prototypes
    for Phoswich Wall analysis, based on analysis code from D. Sarantites
    and W. Reviol (Washington University).

    Author: H.L. Crawford
    Date: May 23, 2014
*/

#ifndef __PHOSWALL_H
#define __PHOSWALL_H

using namespace std;

#include <map>
#include <cstring>
#include <vector>
#include <iostream>

#include "TObject.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

/* Definitions */
// #define CHIPS 256

class phosWallFull;

struct Coeffs {
  Int_t up, dn;
  Int_t lf, rt;
  Int_t uL, uR, dL, dR;
};

struct Pars {
  Float_t axL, bxL, axR, bxR, ayD, byD, ayU, byU;
};

class Coordinates : public TObject {
 public:
  Float_t CM, cm, a, b, c;
  Pars **pars;

  Coordinates();
  Coordinates(Int_t nTele0, Int_t nPix0, TString fName, Bool_t printErrors);
  ~Coordinates();
  Float_t algo(Int_t iTele, Int_t pixel, Int_t ipat, Float_t f, Float_t h);
  
 private:
  Int_t nTele;
  Int_t nPix;

  ClassDef(Coordinates, 1);
};

class Neighbours : public TObject {
 public:
  Int_t nPix;
  
  Neighbours();
  Neighbours(Int_t nTele0, Int_t nPix0, TString fileName, Coordinates* corr0);
  ~Neighbours();
  Float_t getBigOne(Int_t iTele, Int_t *pPix, Float_t *phA, Int_t size,
		    Int_t *pBigA, Float_t *bigA, Int_t maxIndex,
		    Int_t abc, Int_t iHist, Int_t LRDU);

 private:
  Int_t nTele;
  Coeffs **coeffs;
  Coordinates *coor;

  ClassDef(Neighbours, 0);
};

struct Coeff {
  Float_t x1; /* Slope */
  Float_t x0; /* Offset */
};

class Calibrate : public TObject {
 public:
  Calibrate();
  Calibrate(Int_t nTele0, Int_t nPix0, TString fname, Bool_t printErrors);
  ~Calibrate();
  Float_t getEnergy(Int_t iTele, Int_t iPix, Float_t channel, 
		    Int_t type, Int_t check);

 private:
  Int_t nTele;
  Int_t nPix;
  Coeff **coeff;

  ClassDef(Calibrate, 1);
};

class GainMatch : public TObject {
 public:
  GainMatch();
  GainMatch(Int_t nTele0, Int_t nPix0, TString fname, Bool_t printErrors);
  ~GainMatch();
  Float_t GetGainEnergy(Int_t iTele, Int_t iPix, Float_t channel,
			Int_t type, Int_t check);
  
 private:
  Int_t nTele;
  Int_t nPix;
  Coeff **coeff;

  ClassDef(GainMatch, 1);
};

struct hitEvent {
  Int_t pixel;
  Int_t det;
  Int_t A, B, C, T;
};

class phosWallRaw: public TObject {
 public:
  std::vector<hitEvent> hit;

 public:
  phosWallRaw();
  ~phosWallRaw();
  void Initialize();
  void Reset();
  Int_t nHits();

  ClassDef(phosWallRaw, 1);
};

struct calEvent {
  Int_t pixel;
  Float_t A, B, C, T;
  TVector3 Vec;
};

class phosWallCalc: public TObject {
 public:
  std::vector<calEvent> hit;
  Float_t hugeA, hugeB, hugeC;
  Float_t uA, dA, lA, rA; /* Neighbour amplitudes */
  Float_t uF, dF, lF, rF; /* Neighbour fractions */
  Float_t ulA, urA, dlA, drA; /* Corner amplitudes */
  Float_t ulF, urF, dlF, drF; /* Corner fractions */
 
 public:
  phosWallCalc();
  ~phosWallCalc();
  void Initialize();
  void Reset();
  Float_t biggestA(TString var);
  Float_t biggestB(TString var);
  Float_t biggestC(TString var);
  Float_t biggestT(TString var);
  Int_t biggestPixel(TString var);
  TVector3 biggestPixelPosition(TString var);
  Int_t biggestHit(TString var);
  Float_t secondBiggestA(TString var);
  Float_t secondBiggestB(TString var);
  Float_t secondBiggestC(TString var);
  Float_t secondBiggestT(TString var);
  Int_t secondBiggestPixel(TString var);
  Int_t secondBiggestHit(TString var);
  Double_t distanceBetweenParticles(TString var);
  
  ClassDef(phosWallCalc, 1);
};

struct auxHit {
  Int_t channel;
  Int_t value;
};

class phosWallAux: public TObject {
 public:
  std::vector<auxHit> qdc;
  std::vector<auxHit> tdc;

 public:
  phosWallAux();
  ~phosWallAux();
  void Reset();

  ClassDef(phosWallAux, 1);
};

class phosWallFull: public TObject {
 public:
  UShort_t buf[2048];
  phosWallRaw raw;
  phosWallCalc calc;
  phosWallAux aux;
  
  /* Timestamp */
  long long timestamp;

  /* Control parameters initialized from file. */
  Int_t dataType;
  Bool_t showSetup; /* Show setup details with printing flag */
  Bool_t doPrinting; /* Print details about events */
  Bool_t deBug;
  Int_t LRDU; /* For CsI, 2: LR, 4: DU; for plastic, 1: LR, 3: DU */
  Int_t pixOrder; /* Should be 1 (horizontal banks of 8 from LR, DU) */
  Int_t minNeighbours; /* Don't include events with fewer than this
			  many neighbours. */
  Int_t tLow, tHigh; /* Time gates applied to bigA, bigB, bigC */
  Bool_t noTimes; /* Flag to enable or ignore time gates. */
  Int_t maxCheck; /* 4 for no corners, 8 with corners for neighbour sums */
  
  /* Calibration and mapping files */
  TString coefFileA, coefFileB, coefFileC, coefFileT;
  TString coefFileAg, coefFileBg, coefFileCg;
  TString coordFile, neighbourFile;
  Int_t failCount;

  /* Calibration parameters */
  Calibrate *Am, *Bm, *Cm, *Tm;
  GainMatch *Ag, *Bg, *Cg;
  Coordinates *myCoor;

  /* Counter variables */
  Int_t eventNum;
  Int_t oobEvents, oneEvents, badEvents;

  Int_t neighbourhood[256][9];
  Int_t timeGates[256][3];

  std::map<int, TVector3*> pWallPositions;

 public:
  phosWallFull();
  ~phosWallFull();
  
 public:
  void Initialize();
  Int_t InitializeParameters(TString fileName);
  void Reset();
  Int_t getAndUnpackPhosWall(FILE* inf, Int_t length);
  Int_t getAndUnpackPhosWallAux(FILE* inf, Int_t length);
  void ProcessPhosWall();
  void ProcessPhosWallAux();
  Bool_t ifExists(TString fileName);
  void ReadNeighbours(TString file);
  void ReadTimeGates(TString file);
  Bool_t isNeighbour(Int_t first, Int_t second);
  Float_t correctNeighbours(Int_t pixel, Int_t type, Int_t iHist, Float_t input);
  void SetPwallPositions();
  Int_t DetermineNParticles(Float_t radius, Int_t biggestHit);
    
 private:
  ClassDef(phosWallFull, 1);
};

#endif
