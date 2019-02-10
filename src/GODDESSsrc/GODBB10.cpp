#include "GODDESS.h"

/************************************************************/
/* BB10 detector functions                                  */
/* Si detector                                              */
/************************************************************/

BB10::BB10() {
  orrubaDet::SetNumChannels(8, 0);
  Clear();
}

BB10::BB10(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position) :
  orrubaDet(serial_Num, sector_, depth_, upstream_, position) 
{
  orrubaDet::SetNumChannels(8, 0);
  Clear();
}  

void BB10::Clear() {
  orrubaDet::Clear();
  
  stripsP.clear();
  eRawP.clear(); eCalP.clear();
  timeP.clear();
  
  eventPos.SetXYZ(0,0,0);
}

void BB10::ConstructBins() {
  activeWidth = 40.3;
  for (Int_t i = 0; i<8; i++) {
    TVector3 pStPosRefDetCenter( ((7./16.)*activeWidth) - (i*activeWidth/8.), 0., 0.);
    pStPosRefDetCenter.SetPhi(pStPosRefDetCenter.Phi() + detPos.GetRotZ());
    pStripCenterPos[i] = detPos.GetTVector3() + pStPosRefDetCenter;
  }
}

Int_t BB10::GetChannelMult(Bool_t calibrated) { 
  if (calibrated) { return eCalP.size(); }
  return eRawP.size();
}

Int_t BB10::GetChannelMult(Bool_t nType, Bool_t calibrated) {
  if (calibrated) { return eCalP.size(); }
  return eRawP.size();
}

void BB10::LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds) {
  std::map<Short_t, Short_t>::iterator itr;
  std::map<Short_t, Short_t> mapP;
  mapP = GetChannelMap(0);

  Short_t beginP, endP;
  itr = mapP.begin();  beginP = itr->second;  
  itr = mapP.end();  itr--;  endP = itr->second;
  
  for (Int_t i=0; i<ev->channels.size(); i++) {
    if (ev->channels[i] >= beginP && ev->channels[i] <= endP) {
      for (itr=mapP.begin(); itr!=mapP.end(); ++itr) {
	if (ev->channels[i] == itr->second) {
	  SetRawEValue(itr->first-1, 0, ev->values[i], ignoreThresholds);
	  break;
	}
      }
    }
  }
}

void BB10::SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, Int_t ignoreThresh) {
  orrubaDet::SetRawEValue(detChannel, nType, rawValue, ignoreThresh);
}

void BB10::SortAndCalibrate(Bool_t doCalibrate) {
  orrubaDet::valueMap ePMap = GetRawE(kFALSE);
  orrubaDet::timeMap tsPMap = GetTSmap(kFALSE);
  
  for (std::map<Short_t, Float_t>::iterator chItr = ePMap.begin(); chItr != ePMap.end(); ++chItr) {
    Int_t st_ = chItr->first;
    Float_t e_ = chItr->second;

    stripsP.push_back(st_);
    eRawP.push_back(e_);
    timeP.push_back(tsPMap[st_]);
    
    if(doCalibrate) {
      std::vector< std::vector<Float_t> > calPar = GetECalParameters(kFALSE);
      eCalP.push_back(e_*calPar[st_][1] + calPar[st_][0]);
    }
  }
}

void BB10::GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
			 Int_t* stripMaxN, uint64_t* timestampMaxN, 
			 Bool_t calibrated) {
  std::vector<Float_t>* energiesP_;
  if (calibrated) { energiesP_ = &eCalP; }
  else { energiesP_ = &eRawP; }

  if (stripMaxP != nullptr) { *stripMaxP = -1; }

  Float_t eMax = 0.0;
  for(UInt_t i=0; i<energiesP_->size(); i++) {
    if (energiesP_->at(i) > eMax) {
      if (stripMaxP != nullptr) *stripMaxP = stripsP.at(i);
      eMax = energiesP_->at(i);
      if (timestampMaxP != nullptr) *timestampMaxP = timeP.at(i);
    }
  }

}

Float_t BB10::GetESum(Bool_t nType, Bool_t calibrated) {
  Float_t eSum = 0.0;
  std::vector<Float_t>* toSum;
  if (!nType && calibrated) { toSum = &eCalP; }
  else if (!nType && !calibrated) { toSum = &eRawP; }
  else { return 0.0; }
  for (UInt_t i=0; i<toSum->size(); i++) {
    eSum += toSum->at(i);
  }
  return eSum;
}

TVector3 BB10::GetEventPosition(Bool_t calibrated) {
  Int_t pStripHit = 0;

  TVector3 interactionPos;
  GetMaxHitInfo(&pStripHit, nullptr, nullptr, nullptr, calibrated);

  if (pStripHit >= 0 && pStripHit < GetNumChannels(kFALSE)) {
    interactionPos = pStripCenterPos[pStripHit];
  } else {
    interactionPos.SetXYZ(0, 0, 0);
  }


  return interactionPos;
}
