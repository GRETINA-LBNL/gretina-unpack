#ifndef __BGSStructures_H
#define __BGSStructures_H

#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"
#include "TString.h"
#include <vector>
#include <stdint.h>

/* Definitions */

#define maxBGSwords 480
#define nFPstrips 32

/* BGS Data Classes */

struct dataValue {
  Short_t param; /* Parameter number */
  UShort_t value; /* 16-bit data value */
};

struct bgsStruct {
  dataValue data[maxBGSwords];
};

struct rectPoint {
  Double_t x, y, z;
};

struct gGamma {
  Double_t segmentSum;
  Double_t energy;
  Float_t fom;
  long int time;
};

struct EVR {
  Bool_t implanted;
  ULong64_t timestamp;
  Double_t time; 
  Double_t energy;
  Double_t tof;
  Double_t eGamma[100];
  Double_t segEGamma[100];
  Float_t fomGamma[100];
  long int bgsTGamma[100];
};

struct DECAY { /* Covers all types: CE, alpha, fission, etc. */
  ULong64_t timestamp;
  Double_t time; 
  Double_t energy;
};

class BGSGamma : public TObject {
  
 public:
  vector<gGamma> vgam;

 public:
  BGSGamma() { ; }
  ~BGSGamma() { ; }

 private:
  ClassDef(BGSGamma, 1);
};

class BGSEvent : public TObject {

 public:
  ULong64_t timestamp;
  Double_t time;
  Double_t time10MHz;
  Double_t lastTime10MHz;
  Int_t overflowsTime10MHz;
  Double_t tof; /* TOF from MWPC to the FP */
  
  Double_t fCal[3][nFPstrips];
  Double_t bCal[nFPstrips];
  Float_t eGam[16];
  Float_t ecalGam[16];
  Double_t eMWPC;
  Double_t tMWPC;
  
  int EVR, CONV, ALPHA, ESCAPE, FISSION;
  int EVENTTYPE;
  
  struct EVR ievr;
  BGSGamma *ievrGamma;

  struct DECAY dparentCE, ddaughCE, dgdaughCE;
  struct DECAY dparentAL, ddaughAL, dgdaughAL;
  struct DECAY dparentX, ddaughX, dgdaughX;
  struct DECAY dfission;

  Double_t fEnergy1; Int_t fChip1, fStrip1; /* Largest front FP energy info */
  Double_t fEnergy2; Int_t fChip2, fStrip2; /* 2nd largest front FP energy info */
  Double_t bEnergy1; Int_t bStrip1; /* Largest back FP energy info. */
  Double_t bEnergy2; Int_t bStrip2; /* 2nd largest back FP energy info. */
  Double_t uEnergy1; Int_t uStrip1; /* Largest upstream energy info. */
  Double_t pEnergy1; Int_t pStrip1; /* Largest punchthru energy info. */
  Double_t wEnergy1; Int_t wStrip1; /* Largest wing detector energy info. */

  Double_t orEnergy; /* Energy in FP @ origin of event. */
  Double_t teEnergy; /* Energy in FP @ termination of event. */
  Double_t siEnergy; /* Total Si energy (reconstructed). */
  Double_t ceEnergy; /* C.E. (beta) energy */
  Double_t alEnergy; /* Recoil-corrected alpha energy */
  Double_t fiEnergy; /* Fission energy */
  
  struct rectPoint orPos; /* Position (x, y, z) in cm of origin of event */
  Int_t orChip, orFront, orBack; /* Chip and strip # of origin of event. */
  struct rectPoint tePos; /* Position (x, y, z) in cm of terminus of event */
  Int_t teChip, teFront, teBack; /* Chip and strip # of terminus of event. */
  
  Double_t evrEnergy;

 public:
  BGSEvent() { ClearBGS(); }
  ~BGSEvent() { ; }

  void ClearBGS() {
    timestamp = 0;
    time = 0.0;
    time10MHz = 0.0; lastTime10MHz = 0.0;
    overflowsTime10MHz = 0;

    tof = 0.0;
    
    for (Int_t s=0; s<nFPstrips; s++) {
      for (Int_t c=0; c<3; c++) {
	fCal[c][s] = -1.0;
      }
      bCal[s] = -1.0;
    }
    for (Int_t ng=0; ng<16; ng++) {
      eGam[ng] = -1.0;
      ecalGam[ng] = -1.0;
    }
    
    eMWPC = 0.0; tMWPC = 0.0;
    EVR = 0; CONV = 0; ALPHA = 0; ESCAPE = 0; FISSION = 0; 
    EVENTTYPE = 0;
    
    ievr.timestamp = 0; ievr.time = 0; ievr.tof = 0; ievr.energy = 0;
    ievrGamma.vgam.clear();

    dparentCE.timestamp = 0; dparentCE.time = 0; dparentCE.energy = 0;
    ddaughCE.timestamp = 0; ddaughCE.time = 0; ddaughCE.energy = 0;
    dgdaughCE.timestamp = 0; dgdaughCE.time = 0; dgdaughCE.energy = 0;
    dparentAL.timestamp = 0; dparentAL.time = 0; dparentAL.energy = 0;
    ddaughAL.timestamp = 0; ddaughAL.time = 0; ddaughAL.energy = 0;
    dgdaughAL.timestamp = 0; dgdaughAL.time = 0; dgdaughAL.energy = 0;
    dparentX.timestamp = 0; dparentX.time = 0; dparentX.energy = 0;
    ddaughX.timestamp = 0; ddaughX.time = 0; ddaughX.energy = 0;
    dgdaughX.timestamp = 0; dgdaughX.time = 0; dgdaughX.energy = 0;
    dfission.timestamp = 0; dfission.time = 0; dfission.energy = 0;

    fEnergy1 = 0; fChip1 = -1; fStrip1 = -1;
    fEnergy2 = 0; fChip2 = -1; fStrip2 = -1;
    bEnergy1 = 0; bStrip1 = -1;
    bEnergy2 = 0; bStrip2 = -1;
    uEnergy1 = 0; uStrip1 = -1;
    pEnergy1 = 0; pStrip1 = -1;
    wEnergy1 = 0; wStrip1 = -1;

    orEnergy = 0; teEnergy = 0; siEnergy = 0; 
    ceEnergy = 0; alEnergy = 0; fiEnergy = 0;

    orPos.x = 0; orPos.y = 0; orPos.z = 0;
    orChip = -1; orFront = -1; orBack = -1;
    tePos.x = 0; tePos.y = 0; tePos.z = 0;
    teChip = -1; teFront = -1; teBack = -1;
 
    evrEnergy = 0.;
  }
  
 private:
  ClassDef(BGSEvent, 1);
};

/* BGS Calibration Class - all the calibration parameters */

class BGSCalibration : public TObject {

 public:
  /* Low energy calibrations for the DSSD, punchthrus and upstreams
     E = (channel - La)*Lm; */
  Double_t FPLm[4][nfp], FPLa[4][nfp];
  Double_t PTLm[npt], PTLa[npt];
  Double_t USLm[nus], USLa[nus];

  /* High energy calibrations are a two-part function:
     Below ~ 15MeV it's linear, above it's exponential.
     The plan is to calculate an energy from both the linear
     and exponential curves, and take the higher of the two.
     The linear portion is: E = (channel - Ha)*Hm;
     The exponential part is: E = Hc*exp(Hd*channel); */
  Double_t FPHm[4][nfp], FPHa[4][nfp];
  Double_t FPHc[4][nfp], FPHd[4][nfp];
  Double_t WHm[nws], WHa[nws];
  Double_t WHc[nws], WHd[nws];
  Double_t USHm[nus], USHa[nus];
  Double_t USHc[nus], USHd[nus];

  /* Germanium calibration (not used in GRETINA + BGS) */
  Double_t GEa[nge], GEb[nge], GEc[nge], GEd[nge];
  Double_t GTm[nge], GTa[nge];

  /* TDC calibrations 
     0 = top, 1 = east, 2 = west, 3 = back, 4 = upstream, 5 = mwpc, 6 = punchthru */
  Double_t FPTm[7], FPTa[7];

  /* Thresholds for FP detectors*/
  Double_t FPLMinChan[4][nfp]; Double_t FPLMaxChan[4][nfp];
  Double_t PTMinChan[npt]; 
  Double_t USMinChan[nus]; Double_t USMaxChan[nus];
  Double_t WMinChan[nws];

  /* Thresholds for germaniums */
  Double_t GEMinChan[nge], GEMaxChan[nge], GTMinChan[nge], GTMaxChan[nge];

 public:
  BGSCalibration() { ; }
  ~BGSCalibration() { ; }

 private:
  ClassDef(BGSCalibration, 1);
};

/* BGS Conditions Class - all the analysis conditions (i.e. time windows) */

class BGSConditions : public TObject {

 public:
  Float_t rMinE, rMaxE; /* Energy window for recoils */
  
  Float_t pCEMinE, pCEMaxE; /* Energy window for parent CE (keV) */
  Float_t pCEMinT, pCEMaxT; /* Lifetime window for parent CE (sec) */
  Float_t dCEMinE, dCEMaxE; /* Energy window for daughter CE (keV) */
  Float_t dCEMinT, dCEMaxT; /* Lifetime window for daughter CE (sec) */
  Float_t gCEMinE, gCEMaxE; /* Energy window for granddaughter CE (keV) */
  Float_t gCEMinT, gCEMaxT; /* Lifetime window for granddaughter CE (sec) */

  Float_t pALMinE, pALMaxE; /* Energy window for parent alpha (keV) */
  Float_t pALMinT, pALMaxT; /* Lifetime window for parent alpha (sec) */
  Float_t dALMinE, dALMaxE; /* Energy window for daughter alpha (keV) */
  Float_t dALMinT, dALMaxT; /* Lifetime window for daughter alpha (sec) */
  Float_t gALMinE, gALMaxE; /* Energy window for granddaughter alpha (keV) */
  Float_t gALMinT, gALMaxT; /* Lifetime window for granddaughter alpha (sec) */

  Float_t pXMinE, pXMaxE; /* Energy window for parent escape (keV) */
  Float_t pXMinT, pXMaxT; /* Lifetime window for parent escape (sec) */
  Float_t dXMinE, dXMaxE; /* Energy window for daughter escape (keV) */
  Float_t dXMinT, dXMaxT; /* Lifetime window for daughter escape (sec) */
  Float_t gXMinE, gXMaxE; /* Energy window for granddaughter escape (keV) */
  Float_t gXMinT, gXMaxT; /* Lifetime window for granddaughter escape (sec) */

  Float_t fMinE, fMaxE; /* Energy window for fissions (keV) */
  Float_t fMinT, fMaxT; /* Lifetime window for fissions (sec) */
  
  Float_t mwaMin, mwaMax; /* MWPC anode channel window */
  Float_t mwaFPTDCMin, mwaFPTDCMax; /* MWPC anode - DSSD TDC channel window */
  Float_t mwcFPTDCMin, mwcFPTDCMax; /* MWPC cathode - DSSD TDC channel window */
  
  Float_t geTDCMin[nge], geTDCMax[nge]; /* Ge TDC channel windows */
  
  Int_t bufferDepth; /* Ring buffer depth */
  Float_t ADRC; /* Correction for energy calibration with external source */
  
  /* Decay conditions: reconstruct? (if > 1, with same backstrip OK)
                       minimum total electron E, maximum total electron E
	 	       front and back must agree within # keV             */
  Int_t reconBeta; Double_t betaMinE, betaMaxE, betaMatchE; 
  Int_t reconAlpha; Double_t alphaMinE, alphaMaxE, alphaMatchE;
  Float_t alphaMinOrE, alphaMinTeE, alphaMinTeERecon; /* Minimum terminus/origin energies */
  Int_t reconFiss; Double_t fissMinE, fissMaxE, fissMatchE;
  Float_t fissMinOrE;

  Float_t ptTDCMin, ptTDCMax; /* Punchthru TDC window */
  Float_t usTDCMin, usTDCMax; /* Upstream TDC window */
  
  Float_t fptTDCMin, fptTDCMax; /* FP top chip TDC window */
  Float_t fpwTDCMin, fpwTDCMax; /* FP west chip TDC window */
  Float_t fpeTDCMin, fpeTDCMax; /* FP east chip TDC window */
  Float_t fpbTDCMin, fpbTDCMax; /* FP back chip TDC window */

 public:
  BGSConditions() { ; }
  ~BGSConditions() { ; } 

 private:
  ClassDef(BGSConditions, 1);
};

#endif
