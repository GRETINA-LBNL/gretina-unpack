/*! \file GODDESS.h
    \brief Parameter and function definitions for GODDESS (ORRUBA) analysis.

    This file provides the data structures and function prototypes
    for GODDESS analysis, based on analysis code from Rutgers/ORNL.

    Author: H.L. Crawford
    Date: January 29, 2019
*/

#ifndef __GODDESS_H
#define __GODDESS_H

using namespace std;

#include <map>
#include <cstring>
#include <vector>
#include <iostream>

#include "TObject.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

/* Raw data holder, I think */
struct goddessEvent {
  std::vector<UShort_t> channels;
  std::vector<UShort_t> values;
  unsigned long long timestamp;
};

/* ORRUBA detector base class -- taken from Detector.h and siDet.h in ORRUBA code */
class orrubaDet : public TObject {
 public:
  typedef std::map<Short_t, Float_t> valueMap;
  typedef std::map<Short_t, unsigned long long> timeMap;
  Bool_t nType = kTRUE;
  Bool_t pType = kFALSE;

 private:
  UShort_t numPtype, numNtype;
  valueMap eRawPmap, eRawNmap, eCalPmap, eCalNmap;
  timeMap timePmap, timeNmap;
  
  // Calibration parameters and thresholds
  std::vector<std::vector<Float_t>> parECalP, parECalN;
  std::vector<Int_t> threshP, threshN;
  
 public:
  orrubaDet() { ; }
  virtual ~orrubaDet() { ; }
  void Clear();
  void SetNumContacts(Int_t pType, Int_t nType = 0);
  Bool_t ValidContact(UInt_t contact, Bool_t nType);

  virtual void SetRawValue(UInt_t detChannel, UInt_t rawValue, Int_t ignoreThresh) = 0;
  virtual void SetRawValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, 
			   Int_t ignoreThresh) = 0;
  virtual void SetTimestamp(UInt_t detChannel, unsigned long long timestamp) = 0;

  virtual Float_t GetESum(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE) = 0;
  virtual void GetMaxHitInfo(Int_t* stripMaxP, unsigned long long int* timestampMaxP, 
			     Int_t* stripMaxN, unsigned long long int* timestampMaxN,
			     Bool_t calibrated = kTRUE) = 0;
			     
  virtual void SortAndCalibrate(Bool_t doCalibrate = kTRUE) = 0;
  
  virtual std::vector<Float_t> GetHitsInfo(std::string info, std::vector<Float_t> dest = nullptr) = 0;
  virtual std::vector<Int_t> GetHitsInfo(std::string info, std::vector<Int_t> dest = nullptr) = 0;
  virutal std::vector<long long unsigned int> GetHitsInfo(std::string info, std::vector<long long unsigned int>* dest = nullptr) = 0;
  virtual Int_t GetMultiplicity(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE) = 0;

  Float_t GetCalEnergy(Int_t contact, Bool_t nType = kFALSE);
 
  virtual Int_t GetContactMult(Bool_t calibrated = kTRUE) = 0;
  virtual Int_t GetContactMult(Bool_t contactType, Bool_t calibrated) = 0;
  
  unsigned long long GetTimestamp();


virtual Int_t getNumChannels() = 0;

 private:
  ClassDef (orrubaDet, 1);
};

class superX3 : public orrubaDet {
 public:
  Double_t activeLength, activeWidth;
  std::vector<Int_t> stripsP;
  std::vector<Float_t> eNearRaw, eFarRaw;
  std::vector<Float_t> eNearCal, eFarCal;
  std::vector<long long unsigned int> timeNear, timeFar;
  
  std::vector<Int_t> stripsN;
  std::vector<Float_t> eRawN, eCalN;
  std::vector<long long unsigned int> timeN;

 private:
  TVector3 pStripEdgePos[5], nStripEdgePos[5]; // vector to mid point of edge, in mm
  TVector3 pStripCenterPos[4], nStripCenterPos[4]; // absolute center, in mm
  TVector3 eventPos;

  Float_t binsP[5], binsN[5]; 
  Float_t binsPcenter[4], binsNcenter[4];
  Float_t binsZ[5], binsZCenter[4];
  Float_t binsAzimuthal[5], binsAzimuthalCenter[4], binsPolar[5];

  orrubaDet::valueMap stripPosRaw, stripPosCal;
  Int_t stripContactMult[4];
  std::vector<Float_t> parPosCal[4], parStripECal[4];

 public:
  superX3();
  superX3( std::string serialNum, UShort_t sector, UShort_t depth, Bool_t upstream);
  virtual ~superX3();

  void Clear();
  void UpdatePosition(Int_t strip);
  
  TVector3 GetPStripCenterPos(Int_t strip) { return pStripCenterPos[strip]; }
  TVector3 GetNStripCenterPos(Int_t strip) { return nStripCenterPos[strip]; }
  Int_t GetNumBins() { return 4; }
  Float_t* GetPtypeBins() { return binsP; }
  Float_t* GetNtypeBins() { return binsN; }
  Float_t* GetPtypeCenterBins() { return binsPcenter; }
  Float_t* GetNtypeCenterBins() { return binsNcenter; }
  Float_t* GetZbins() { return binsZ; }
  Float_t* GetZCenterBins() { return binsZcenter; }
  Float_t* GetAzimuthalBins() { return binsAzimuthal; }
  Float_t* GetAzimuthalCenterBins() { return binsAzimuthalCenter; }
  Float_t* GetPolarBins() { return binsPolar; }
  
  orrubaDet::valueMap GetStripPosRaw() { return stripPosRaw; }
  orrubaDet::valueMap GetStripPosCal() { return stripPosCal; }
  std::vector<Float_t>* GetResStripParCal() { return parStripECal; }

  UShort_t GetNearContact(UShort_t strip);
  UShort_t GetFarContact(UShort_t strip);
  

 private:
  void ContructBins();

  ClassDef(superX3, 1);

};



class goddessFull: public TObject {
 public:
  UInt_t buf[2048];

  goddessEvent *rawGoddess;

 public:
  goddessFull() { ; }
  ~goddessFull() { ; }
  void getAnalogGoddess(FILE *inf, Int_t evtLength);
  void printAnalogRawEvent();

 public:
    
 private:
  ClassDef(goddessFull, 1);
};

#endif
