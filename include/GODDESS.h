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
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdint.h>

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
  typedef std::map<Short_t, uint64_t> timeMap; 
  typedef std::map<Short_t, Short_t> chMap;
  Bool_t nTYPE; 
  Bool_t pTYPE; 

 private: 
  UShort_t numPtype, numNtype; 
  valueMap eRawPmap, eRawNmap, eCalPmap, eCalNmap; 
  timeMap timePmap, timeNmap; 

  chMap chMapP, chMapN;
  
  // Calibration parameters and thresholds 
  std::vector< std::vector<Float_t> > parECalP, parECalN; 
  std::vector<Int_t> threshP, threshN; 
  
 protected:
  std::string serialNum; 
  std::string posID;
  UShort_t sector; 
  UShort_t depth; // # of detectors between this one and the target
  Bool_t upStream;
  
 public: 
  orrubaDet();
  orrubaDet(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_);
  virtual ~orrubaDet() { ; } 
  void Clear(); 

  std::string GetSerialNum() { return serialNum; }
  std::string GetPosID() { return posID; }
  UShort_t GetSector() { return sector; }
  UShort_t GetDepth() { return depth; }
  Bool_t GetUpStream() { return upStream; }

  void SetDetector(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_); // from orrubaDet

  void SetNumChannels(Int_t pType, Int_t nType = 0); 
  Int_t GetNumChannels(Bool_t nType); 
  Bool_t ValidChannel(UInt_t channel, Bool_t nType); 
  Bool_t ChannelHit(UInt_t detChannel, Bool_t nType);
  void SetChannelMap(std::map<Short_t, Short_t> channelMap, Bool_t nType);
  std::map<Short_t, Short_t> GetChannelMap(Bool_t nType);

  virtual void SetRawEValue(UInt_t detChannel, UInt_t rawValue, Int_t ignoreThresh); 
  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue,  
   			   Int_t ignoreThresh); 
  virtual orrubaDet::valueMap GetRawE(Bool_t nType);

  virtual orrubaDet::valueMap GetCalE(Bool_t nType);

  virtual void SetTimestamp(UInt_t detChannel, Bool_t nType, uint64_t timestamp); 
  virtual uint64_t GetTimestamp(UInt_t detChannel, Bool_t nType);
  virtual orrubaDet::timeMap GetTSmap(Bool_t nType);

  virtual Bool_t SetECalParameters(std::vector<Float_t> par, Int_t contact, Bool_t nType);
  virtual std::vector< std::vector<Float_t> > GetECalParameters(Bool_t nType);
  virtual Bool_t SetThresholds(std::vector<Int_t> thresholds, Bool_t nType, Int_t thrSize);
  virtual std::vector<Int_t> GetThresholds(Bool_t nType);
  
  //  virtual TVector3 GetEventPosition(Bool_t calibrated = kTRUE); // from orrubaDet

  //virtual Float_t GetESum(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE) = 0; 
  //virtual void GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
  //			     Int_t* stripMaxN, uint64_t* timestampMaxN, 
  //			     Bool_t calibrated = kTRUE) = 0; 
			     
  //virtual void SortAndCalibrate(Bool_t doCalibrate = kTRUE) = 0;  // from orrubaDet
  
  //virtual std::vector<Float_t> GetHitsInfo(std::string info, std::vector<Float_t> dest = nullptr) = 0; 
  //virtual std::vector<Int_t> GetHitsInfo(std::string info, std::vector<Int_t> dest = nullptr) = 0; 
  //virtual std::vector<long long unsigned int> GetHitsInfo(std::string info, std::vector<long long unsigned int>* dest = nullptr) = 0; 
  //virtual Int_t GetMultiplicity(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE) = 0; 

  Float_t GetCalEnergy(Int_t contact, Bool_t nType = kFALSE); 
 
  //virtual Int_t GetContactMult(Bool_t calibrated = kTRUE) = 0; 
  //virtual Int_t GetContactMult(Bool_t contactType, Bool_t calibrated) = 0; 
  
/*   unsigned long long GetTimestamp(); */

 protected:
  virtual void SetPosID();
   
 private: 
  ClassDef (orrubaDet, 1); 
}; 

class superX3 : public orrubaDet { 
 public: 
  Double_t activeLength, activeWidth; 
  std::vector<Int_t> stripsP; 
  std::vector<Float_t> eNearRaw, eFarRaw; 
  std::vector<Float_t> eNearCal, eFarCal; 
  std::vector<uint64_t> timeNear, timeFar; 
  
  std::vector<Int_t> stripsN; 
  std::vector<Float_t> eRawN, eCalN; 
  std::vector<uint64_t> timeN; 

 private: 
  /*   TVector3 pStripEdgePos[5], nStripEdgePos[5]; // vector to mid point of edge, in mm */
  /*   TVector3 pStripCenterPos[4], nStripCenterPos[4]; // absolute center, in mm */
  TVector3 eventPos; 
  
  /*   Float_t binsP[5], binsN[5];  */
  /*   Float_t binsPcenter[4], binsNcenter[4]; */
  /*   Float_t binsZ[5], binsZCenter[4]; */
  /*   Float_t binsAzimuthal[5], binsAzimuthalCenter[4], binsPolar[5]; */
  
  orrubaDet::valueMap stripPosRaw, stripPosCal; 
  Int_t stripContactMult[4]; 
  /*   std::vector<Float_t> parPosCal[4], parStripECal[4]; */
  
 public: 
  superX3(); 
  superX3( std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_); 
  virtual ~superX3() { ; } 
  
  void Clear(); 
  
  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, Int_t ignoreThresh); 
  
  /*   void UpdatePosition(Int_t strip); */
  
  /*   TVector3 GetPStripCenterPos(Int_t strip) { return pStripCenterPos[strip]; } */
  /*   TVector3 GetNStripCenterPos(Int_t strip) { return nStripCenterPos[strip]; } */
  /*   Int_t GetNumBins() { return 4; } */
  /*   Float_t* GetPtypeBins() { return binsP; } */
  /*   Float_t* GetNtypeBins() { return binsN; } */
  /*   Float_t* GetPtypeCenterBins() { return binsPcenter; } */
  /*   Float_t* GetNtypeCenterBins() { return binsNcenter; } */
  /*   Float_t* GetZbins() { return binsZ; } */
  /*   Float_t* GetZCenterBins() { return binsZcenter; } */
  /*   Float_t* GetAzimuthalBins() { return binsAzimuthal; } */
  /*   Float_t* GetAzimuthalCenterBins() { return binsAzimuthalCenter; } */
  /*   Float_t* GetPolarBins() { return binsPolar; } */
  
  /*   orrubaDet::valueMap GetStripPosRaw() { return stripPosRaw; } */
  /*   orrubaDet::valueMap GetStripPosCal() { return stripPosCal; } */
  /*   std::vector<Float_t>* GetResStripParCal() { return parStripECal; } */
  
  /*   UShort_t GetNearContact(UShort_t strip); */
  /*   UShort_t GetFarContact(UShort_t strip); */
  
 private: 
  void ContructBins(); 

  ClassDef(superX3, 1); 

}; 

class QQQ5 : public orrubaDet {
 public:
  Double_t firstStripWidth;
  Double_t deltaPitch;

  std::vector<Int_t> stripsP;
  std::vector<Float_t> eRawP, eCalP;
  std::vector<uint64_t> timeP;
  
  std::vector<Int_t> stripsN;
  std::vector<Float_t> eRawN, eCalN;
  std::vector<uint64_t> timeN;
  
 private:
  TVector3 eventPos; 
  
 public:
  QQQ5();
  QQQ5(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_);
  virtual ~QQQ5() { ; }

  void Clear();

  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue,  
   			   Int_t ignoreThresh); 


 protected:
  void SetPosID();

 private:
  ClassDef(QQQ5, 1);

};

class BB10 : public orrubaDet {
 public:
  Double_t activeWidth;
  std::vector<Int_t> stripsP;
  std::vector<Float_t> eRawP, eCalP;
  std::vector<uint64_t> timeP;
  
 private:
  TVector3 eventPos;

 public:
  BB10();
  BB10(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_);
  virtual ~BB10() { ; }

  void Clear();

  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue,  
			    Int_t ignoreThresh); 
  

 private:
  ClassDef(BB10, 1);

};

class goddessFull: public TObject {
 public:
  UInt_t buf[2048];
  goddessEvent *rawGoddess;

  QQQ5 *qqqU[4];
  superX3 *s3U[12], *s3D[12];
  BB10 *bb[8];
  QQQ5 *qqqD[4];

  uint64_t ts;

 public:
  goddessFull() { ; }
  ~goddessFull() { ; }
  void Initialize();
  

  void getAnalogGoddess(FILE *inf, Int_t evtLength);
  void printAnalogRawEvent();

  void InitializeDetectors();
  void ReportDetectors();

 public:
    
 private:
  ClassDef(goddessFull, 1);
};

#endif
