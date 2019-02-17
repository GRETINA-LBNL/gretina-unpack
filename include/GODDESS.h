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
#include <fstream>
#include <stdint.h>

#include "TObject.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

/* Raw data holder */
struct goddessEvent {
  std::vector<UShort_t> channels;
  std::vector<UShort_t> values;
  ULong64_t timestamp;
};

/**************************************************************************/
/***** SolidVector to define detector positions                           */
/**************************************************************************/

class solidVector : public TVector3 {
 private:
  Double_t rotZ;
  Double_t rotPhi;

 public:
  solidVector(Double_t x_=0, Double_t y_=0, Double_t z_=0, Double_t rot_Z=0, Double_t rot_Phi = 0);
  
  Double_t GetRotZ() { return rotZ; }
  Double_t GetRotPhi() { return rotPhi; }
  void SetRotationZ(Double_t angle) { rotZ = angle; }
  void SetRotationPhi(Double_t angle) { rotPhi = angle; }
  TVector3 GetTVector3() { return TVector3(X(), Y(), Z()); }

 private:
  ClassDef(solidVector, 1);
};

/**************************************************************************/
/***** ORRUBA Detector base class (from Detector.h, siDet.h, orrubaDet.h) */
/**************************************************************************/

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
  solidVector detPos;
  
 public: 
  /* Initializing and resetting... */
  orrubaDet();
  orrubaDet(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position);
  virtual ~orrubaDet() { ; } 
  void Clear(); 

  /* Detector identification stuff */
  std::string GetSerialNum() { return serialNum; }
  std::string GetPosID() { return posID; }
  UShort_t GetSector() { return sector; }
  UShort_t GetDepth() { return depth; }
  Bool_t GetUpStream() { return upStream; }
  void SetDetector(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position); 

  /* Electronics channel stuff */
  void SetChannelMap(std::map<Short_t, Short_t> channelMap, Bool_t nType);
  std::map<Short_t, Short_t> GetChannelMap(Bool_t nType);

  void SetNumChannels(Int_t pType, Int_t nType = 0); 
  Int_t GetNumChannels(Bool_t nType); 
  Bool_t ValidChannel(UInt_t channel, Bool_t nType); 
  Bool_t ChannelHit(UInt_t detChannel, Bool_t nType);

  /* Loading an event into the detector, setting/getting values. */
  virtual void LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds) = 0;
  
  virtual void SetRawEValue(UInt_t detChannel, UInt_t rawValue, Int_t ignoreThresh); 
  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue,  
			    Int_t ignoreThresh); 
  virtual Float_t GetRawEValue(UInt_t detChannel, Bool_t nType);
  virtual Float_t GetCalEValue(UInt_t detChannel, Bool_t nType);
  virtual orrubaDet::valueMap GetRawE(Bool_t nType);
  virtual orrubaDet::valueMap GetCalE(Bool_t nType);

  virtual void SetTimestamp(UInt_t detChannel, Bool_t nType, uint64_t timestamp); 
  virtual uint64_t GetTimestamp(UInt_t detChannel, Bool_t nType);
  virtual orrubaDet::timeMap GetTSmap(Bool_t nType);

  /* Calibration parameters and thresholds. */
  virtual Bool_t SetECalParameters(std::vector<Float_t> par, Int_t contact, Bool_t nType);
  virtual std::vector< std::vector<Float_t> > GetECalParameters(Bool_t nType);
  virtual Bool_t SetThresholds(std::vector<Int_t> thresholds, Bool_t nType, Int_t thrSize);
  virtual std::vector<Int_t> GetThresholds(Bool_t nType);

  /* Calibrate, and do analysis-ish stuff */  
  virtual void SortAndCalibrate(Bool_t doCalibrate = kTRUE) = 0;
  virtual TVector3 GetEventPosition(Bool_t calibrated = kTRUE) = 0;
  virtual void GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
  			     Int_t* stripMaxN, uint64_t* timestampMaxN, 
  			     Bool_t calibrated = kTRUE) = 0; 

  /* Multiplicity stuff and total energy */
  virtual Int_t GetChannelMult(Bool_t calibrated = kTRUE) = 0;
  virtual Int_t GetChannelMult(Bool_t nType, Bool_t calibrated) = 0;
  virtual Float_t GetESum(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE) = 0;

 protected:
  virtual void SetPosID();
   
 private: 
  ClassDef (orrubaDet, 1); 
}; 

/**************************************************************************/
/***** SuperX3 detector                                                   */
/**************************************************************************/

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
  TVector3 pStripEdgePos[5], nStripEdgePos[5]; // vector to mid point of edge, in mm 
  TVector3 pStripCenterPos[4], nStripCenterPos[4]; // absolute center, in mm 
  TVector3 eventPos; 
  
  Float_t binsP[5], binsPCenter[4];   
  Float_t binsN[5], binsNCenter[4]; 
  Float_t binsZ[5], binsZCenter[4]; 
  Float_t binsAzimuthal[5], binsAzimuthalCenter[4];
  Float_t binsPolar[5]; 
  
  orrubaDet::valueMap stripPosRaw, stripPosCal; 
  Int_t stripChannelMult[4]; 
  std::vector<Float_t> parStripPosCal[4], parStripECal[4]; 

 public: 
  superX3(); 
  superX3( std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position); 
  virtual ~superX3() { ; } 
  
  void Clear(); 

  /* Geometry functions */
  void SetGeomParameters(map<string, Double_t> geoInfo) { 
    activeWidth = geoInfo["SuperX3 Active Width"]; 
    activeLength = geoInfo["SuperX3 Active Length"];
  }
  TVector3 GetPStripCenterPos(Int_t strip) { return pStripCenterPos[strip]; } 
  TVector3 GetNStripCenterPos(Int_t strip) { return nStripCenterPos[strip]; } 
  Int_t GetNumBins() { return 4; } 
  Float_t* GetPtypeBins() { return binsP; } 
  Float_t* GetNtypeBins() { return binsN; } 
  Float_t* GetPtypeCenterBins() { return binsPCenter; } 
  Float_t* GetNtypeCenterBins() { return binsNCenter; } 
  Float_t* GetZbins() { return binsZ; } 
  Float_t* GetZCenterBins() { return binsZCenter; } 
  Float_t* GetAzimuthalBins() { return binsAzimuthal; } 
  Float_t* GetAzimuthalCenterBins() { return binsAzimuthalCenter; } 
  Float_t* GetPolarBins() { return binsPolar; } 
  void ConstructBins();

  /* Resistive strip stuff */
  Bool_t ValidStrip(UShort_t strip) { if (strip>=0 && strip<=3) { return kTRUE; } return kFALSE; }
  Int_t GetStrip(Int_t channel);
  UShort_t GetNearChannel(UShort_t strip);
  UShort_t GetFarChannel(UShort_t strip);
  Float_t GetNearE(Bool_t calibrated = kTRUE);
  Float_t GetFarE(Bool_t calibrated = kTRUE);
  void UpdatePosition(UShort_t strip);
  orrubaDet::valueMap GetStripPosRaw() { return stripPosRaw; } 
  orrubaDet::valueMap GetStripPosCal() { return stripPosCal; } 

  /* Resistive strip calibration.... */
  void SetStripPosCalParameters(Int_t strip, std::vector<Float_t> pars);
  void SetStripECalParameters(Int_t strip, std::vector<Float_t> pars);
  std::vector<Float_t>* GetResStripParCal() { return parStripECal; } 
  
  /* Multiplicity stuff and total energy */
  virtual Int_t GetChannelMult(Bool_t calibrated = kTRUE);
  virtual Int_t GetChannelMult(Bool_t nType, Bool_t calibrated);
  std::vector<Float_t> GetResE(Bool_t calibrated = kTRUE);

  /* Actual main analysis-ish functions */
  virtual void LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds);
  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, Int_t ignoreThresh); 
  
  virtual void SortAndCalibrate(Bool_t doCalibrate);
  virtual Float_t GetESum(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE);
  virtual void GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
  			     Int_t* stripMaxN, uint64_t* timestampMaxN, 
  			     Bool_t calibrated = kTRUE); 

  virtual TVector3 GetEventPosition(Bool_t calibrated = kTRUE);
  
 private: 
  ClassDef(superX3, 1); 

}; 

/**************************************************************************/
/***** QQQ5 detector                                                      */
/**************************************************************************/

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
  TVector3 pStripEdgePos[33], nStripEdgePos[5];
  TVector3 pStripCenterPos[32], nStripCenterPos[4];
  Float_t binsP[33], binsPCenter[32];
  Float_t binsN[5], binsNCenter[4];
  Float_t binsRho[33];
  Float_t binsAzimuthal[5];
  Float_t binsPolar[33], binsPolarCenter[32];

  TVector3 eventPos; 
  
 public:
  QQQ5();
  QQQ5(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position);
  virtual ~QQQ5() { ; }

  void Clear();

  /* Geometry functions */
  void SetGeomParameters(map<string, Double_t> geoInfo) { 
    firstStripWidth = geoInfo["QQQ5 First Strip Width"]; 
    deltaPitch = geoInfo["QQQ5 Delta Pitch"];
  }
  TVector3 GetPStripCenterPos(Int_t strip) { return pStripCenterPos[strip]; }
  TVector3* GetNStripCenterPos() { return nStripCenterPos; }
  Int_t GetNumNTypeBins() { return 4; }
  Int_t GetNumPTypeBins() { return 32; }
  Float_t* GetPTypeBins() { return binsP; }
  Float_t* GetPTypeBinsCenter() { return binsPCenter; }
  Float_t* GetNTypeBins() { return binsN; }
  Float_t* GetNTypeBinsCenter() { return binsNCenter; }
  Float_t* GetRhoBins() { return binsRho; }
  Float_t* GetAzimuthalBins() { return binsAzimuthal; }
  Float_t* GetPolarBins() { return binsPolar; }
  Float_t* GetPolarBinsCenter() { return binsPolarCenter; }
  void ConstructBins();

  /* Multiplicity stuff and total energy */
  virtual Int_t GetChannelMult(Bool_t calibrated = kTRUE);
  virtual Int_t GetChannelMult(Bool_t nType, Bool_t calibrated);
  virtual Float_t GetESum(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE);

  virtual void LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds);
  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue,  
   			   Int_t ignoreThresh); 
  virtual void SortAndCalibrate(Bool_t doCalibrate);

  virtual void GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
  			     Int_t* stripMaxN, uint64_t* timestampMaxN, 
  			     Bool_t calibrated = kTRUE); 

  virtual TVector3 GetEventPosition(Bool_t calibrated = kTRUE);

 protected:
  void SetPosID();

 private:
  ClassDef(QQQ5, 1);

};

/**************************************************************************/
/***** BB10 (not BB8...)  detector                                        */
/**************************************************************************/

class BB10 : public orrubaDet {
 public:
  Double_t activeWidth;
  std::vector<Int_t> stripsP;
  std::vector<Float_t> eRawP, eCalP;
  std::vector<uint64_t> timeP;
  
 private:
  /* Geometry information for BB10 detectors */
  /* Positions of edges and centers in mm/radians */
  TVector3 pStripEdgePos[9];
  TVector3 pStripCenterPos[8];
  Float_t binsP[9];
  Float_t binsAzimuthal[9];
  Float_t binsPCenter[8];
  Float_t binsAzimuthalCenter[8];

  TVector3 eventPos;

 public:
  BB10();
  BB10(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position);
  virtual ~BB10() { ; }

  void Clear();

  /* Geometry functions */
  void SetGeomParameters(map<string, Double_t> geoInfo) { activeWidth = geoInfo["BB10 Active Width"]; }
  TVector3 GetPStripCenterPos(Int_t strip) { return pStripCenterPos[strip]; }
  Int_t GetNumBins() { return 8; }
  Float_t* GetPTypeBins() { return binsP; }
  Float_t* GetAzimuthalBins() { return binsAzimuthal; }
  Float_t* GetPTypeBinCenters() { return binsPCenter; }
  Float_t* GetAzimuthalBinCenters() { return binsAzimuthalCenter; }
  void ConstructBins();
  
  /* Multiplicity stuff and total energy */
  virtual Int_t GetChannelMult(Bool_t calibrated = kTRUE);
  virtual Int_t GetChannelMult(Bool_t nType, Bool_t calibrated);
  virtual Float_t GetESum(Bool_t nType = kFALSE, Bool_t calibrated = kTRUE);

  virtual void LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds);
  virtual void SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue,  
			    Int_t ignoreThresh); 
  virtual void SortAndCalibrate(Bool_t doCalibrate = kTRUE);  

  virtual void GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
  			     Int_t* stripMaxN, uint64_t* timestampMaxN, 
  			     Bool_t calibrated = kTRUE); 

  virtual TVector3 GetEventPosition(Bool_t calibrated = kTRUE);

 private:
  ClassDef(BB10, 1);

};

/**************************************************************************/
/***** Ion chamber                                                        */
/**************************************************************************/

class ionChamber: public TObject {
 public:
  typedef std::map<Short_t, Float_t> valueMap; 
  typedef std::map<Short_t, Short_t> chMap;

 private:
  Int_t numAnode, numWireX, numWireY, numDE, numEres;
  valueMap anodeRaw, anodeCal;
  valueMap wireX, wireY;
  Float_t dE, resE, E; // All in MeV
  std::vector< std::vector<Float_t> > parAnodeECal;
  std::vector<Int_t> threshAnode, threshX, threshY; 

  chMap mapAnode, mapX, mapY;

 public:
  ionChamber() { ; }
  ionChamber(Int_t numAnode, Int_t numWireX, Int_t numWireY, Int_t numDE, Int_t numEres);
  ~ionChamber() { ; }

  void Clear();
  
  void SetChannelMap(std::map<Short_t, Short_t> channelMap, Int_t chType);
  std::map<Short_t, Short_t> GetChannelMap(Int_t chType);

  Bool_t ValidChannel(UInt_t channel, Int_t chType); 

  void SetRawValue(UInt_t channel, Int_t chType, UInt_t rawValue, Int_t ignoreThresh);
  
  Bool_t SetThresholds(std::vector<Int_t> thresholds, Int_t chType);
  std::vector<Int_t> GetThresholds(Int_t chType);
  // void SetAnodeECalPars(Int_t ch, std::vector<Float_t> pars);

  void LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds);

  Int_t GetNumChannels(Int_t chType);
  Float_t GetAnodeDE() { return dE; }
  Float_t GetAnodeResE() { return resE; }
  Float_t GetAnodeE() { return E; }
  Int_t GetMaxX();
  Int_t GetMaxY();

 private:
  ClassDef(ionChamber, 1);
};

struct detectorOUT {
  UShort_t depth;
  UShort_t upstream;
  Int_t pStrip, nStrip;
  Float_t pECal, nECal;
  TVector3 evPos;
  // sx3
  Float_t eNearCal, eFarCal;
  // Int_t nearStrip, farStrip;

};

class goddessOut: public TObject {
 public:
  std::vector<string> siID;
  std::vector<string> siType;
  std::vector<detectorOUT> si;
  Float_t icDE;
  Float_t icE;
  Float_t icX;
  Float_t icY;

 public:
  goddessOut() { ; }
  ~goddessOut() { ; }
  void Clear();

 private:
  ClassDef(goddessOut, 1);
};

/**************************************************************************/
/***** Overall GODDESS structure                                          */
/**************************************************************************/

class goddessFull: public TObject {
 public:
  UInt_t buf[2048];
  goddessEvent *rawGoddess;

  std::vector<QQQ5> qqs;
  std::vector<superX3> s3s;
  std::vector<BB10> bbs;
  std::vector<ionChamber> ics;

  goddessOut *eventOut;

  uint64_t ts;

  std::map<string, Double_t> geomInfo;
  
 public:
  goddessFull() { ; }
  ~goddessFull() { ; }
  void Initialize();
  void Clear();

  Bool_t ParseID(std::string id, Short_t &sector, Short_t &depth, Bool_t &upStream);
  void ReadPositions(TString filename);
  solidVector GetPosVector(std::string type, Short_t sector, Short_t depth, Bool_t upStream);
  void ReadConfiguration(TString filename);

  void InitializeDetectors();

  void getAnalogGoddess(FILE *inf, Int_t evtLength);
  void printAnalogRawEvent();
  void ProcessEvent();
  void ResetOutput(detectorOUT* det);

 public:
    
 private:
  ClassDef(goddessFull, 1);
};

#endif
