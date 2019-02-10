#include "GODDESS.h"

/************************************************************/
/* QQQ5 detector functions                                  */
/* Annular detector                                         */
/************************************************************/

QQQ5::QQQ5() {
  orrubaDet::SetNumChannels(32, 4);
  Clear();
}

QQQ5::QQQ5(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position) :
  orrubaDet(serial_Num, sector_, depth_, upstream_, position) 
{
  orrubaDet::SetNumChannels(32, 4);
  SetPosID();
  Clear();
}  

void QQQ5::Clear() {
  orrubaDet::Clear();

  stripsP.clear(); 
  eRawP.clear(); eCalP.clear();
  timeP.clear();

  stripsN.clear(); 
  eRawN.clear(); eCalN.clear();
  timeN.clear();

  eventPos.SetXYZ(0, 0, 0);

}

void QQQ5::ConstructBins() {
  cout << "FirstStripWidth: " << firstStripWidth << endl;
  cout << "deltaPitch: " << deltaPitch << endl; 

  // Parameters from the geom file, set by hand for now as there is a problem getting the map back out ?!
  deltaPitch = 0.05;
  firstStripWidth = 2.55;

  TVector3 firstStripOffset(0., firstStripWidth/2., 0.);

  TVector3 prevStripRefDetCenter = firstStripOffset;
  pStripCenterPos[0] = detPos.GetTVector3() + firstStripOffset;

  for (Int_t i=1; i<32; i++) {
    Float_t prevStripWidth = firstStripWidth - (i-1)*deltaPitch;
    Float_t currStripWidth = firstStripWidth - i*deltaPitch;
    TVector3 pStPosRefDetCenter = prevStripRefDetCenter + TVector3(0., (prevStripWidth+currStripWidth)/2., 0.);
    prevStripRefDetCenter = pStPosRefDetCenter;
    pStPosRefDetCenter.SetPhi(pStPosRefDetCenter.Phi() + detPos.GetRotZ());
    pStripCenterPos[i] = detPos.GetTVector3() + pStPosRefDetCenter;
  }
}

void QQQ5::SetPosID() {
  posID.clear();
  if (upStream) { posID.append("U"); }
  else { posID.append("D"); }
  static const char *sectorCode[4] = {"A", "B", "C", "D"};
  posID.append(sectorCode[sector]);
  posID.append("_");
  if (depth == 0) { posID.append("dE"); }
  else if (depth == 1) { posID.append("E1"); }
  else if (depth == 2) { posID.append("E2"); }
}

Int_t QQQ5::GetChannelMult(Bool_t calibrated) {
  if(calibrated) { return eCalP.size() + eCalN.size(); }
  return eRawP.size() + eRawN.size(); 
}

Int_t QQQ5::GetChannelMult(Bool_t nType, Bool_t calibrated) {
  if (nType) { 
    if (calibrated) { return eCalN.size(); }
    return eRawN.size();
  } else {
    if (calibrated) { return eCalP.size(); }
    return eRawP.size();
  }
}
    
Float_t QQQ5::GetESum(Bool_t nType, Bool_t calibrated) {
  Float_t eSum = 0.0;
  vector<Float_t>* toSum;

  if (nType && calibrated) { toSum = &eCalN; }
  else if (nType && !calibrated) { toSum = &eRawN; }
  else if (!nType && calibrated) { toSum = &eCalP; }
  else if (!nType && !calibrated) { toSum = &eRawP; }
  else { return 0.0; }
  
  for (UInt_t i=0; i<toSum->size(); i++) {
    eSum += toSum->at(i);
  }
  return eSum;
}

void QQQ5::SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, Int_t ignoreThresh) {
  if (!ValidChannel(detChannel, nType)) {
    cout <<"ERROR: Unable to set raw value for QQQ5 " << serialNum << " " << (nType ? 'n':'p') <<  "-type channel." << endl;
  }
  orrubaDet::SetRawEValue(detChannel, nType, rawValue, ignoreThresh);
}

void QQQ5::SortAndCalibrate(Bool_t doCalibrate) {
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

  orrubaDet::valueMap eNMap = GetRawE(kTRUE);
  orrubaDet::timeMap tsNMap = GetTSmap(kTRUE);
  
  for (std::map<Short_t, Float_t>::iterator chItr = eNMap.begin(); chItr != eNMap.end(); ++chItr) {
    Int_t st_ = chItr->first;
    Float_t e_ = chItr->second;

    stripsN.push_back(st_);
    eRawN.push_back(e_);
    timeN.push_back(tsNMap[st_]);
    
    if(doCalibrate) {
      std::vector< std::vector<Float_t> > calPar = GetECalParameters(kTRUE);
      eCalN.push_back(e_*calPar[st_][1] + calPar[st_][0]);
    }
  }

}

void QQQ5::GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,  
			 Int_t* stripMaxN, uint64_t* timestampMaxN, 
			 Bool_t calibrated) {
  std::vector<Float_t>* energiesN_;
  std::vector<Float_t>* energiesP_;

  if (calibrated) { energiesN_ = &eCalN;  energiesP_ = &eCalP; }
  else { energiesN_ = &eRawN;  energiesP_ = &eRawP; }
  
  if (stripMaxP != nullptr) { *stripMaxP = -1; }
  if (stripMaxN != nullptr) { *stripMaxN = -1; }

  Float_t eMax = 0.0;
  for (UInt_t i=0; i<energiesP_->size(); i++) {
    if (energiesP_->at(i) > eMax) {
      if (stripMaxP != nullptr) *stripMaxP = stripsP.at(i);
      eMax = energiesP_->at(i);
      if (timestampMaxP != nullptr) *timestampMaxP = timeP.at(i);
    }
  }
  eMax = 0.0;
  for (UInt_t i=0; i<energiesN_->size(); i++) {
    if (energiesN_->at(i) > eMax) {
      if (stripMaxN != nullptr) *stripMaxN = stripsN.at(i);
      eMax = energiesN_->at(i);
      if (timestampMaxN != nullptr) *timestampMaxN = timeN.at(i);
    }
  }
} 

void QQQ5::LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds) {
  std::map<Short_t, Short_t>::iterator itr;
  std::map<Short_t, Short_t> mapP, mapN;
  mapP = GetChannelMap(0);
  mapN = GetChannelMap(1);

  Short_t beginP, endP, beginN, endN;
  itr = mapP.begin();  beginP = itr->second;  
  itr = mapP.end();  itr--;  
  endP = itr->second;
  itr = mapN.begin();  beginN = itr->second;
  itr = mapN.end();  itr--;  endN = itr->second;

  for (Int_t i=0; i<ev->channels.size(); i++) {
    if (ev->channels[i] >= beginP && ev->channels[i] <= endP) {
      for (itr=mapP.begin(); itr!=mapP.end(); ++itr) {
	if (ev->channels[i] == itr->second) {
	  SetRawEValue(itr->first-1, 0, ev->values[i], ignoreThresholds);
	  break;
	}
      }
    } else if (ev->channels[i] >= beginN && ev->channels[i] <= endN) {
      for (itr=mapN.begin(); itr!=mapN.end(); ++itr) {
	if (ev->channels[i] == itr->second) {
	  SetRawEValue(itr->first-1, 1, ev->values[i], ignoreThresholds);
	  break;
	}
      }
    }
  }
}

TVector3 QQQ5::GetEventPosition(Bool_t calibrated) {
  Int_t pStripHit, nStripHit;

  TVector3 interactionPos;
  GetMaxHitInfo(&pStripHit, nullptr, &nStripHit, nullptr, calibrated);

  if (pStripHit >= 0 && pStripHit < GetNumChannels(kFALSE)) {
    interactionPos = pStripCenterPos[pStripHit];
  } else { 
    interactionPos.SetXYZ(0.,0.,0.);
  }
  if (nStripHit >= 0 && nStripHit < GetNumChannels(kTRUE)) { 
    interactionPos.SetPhi(interactionPos.Phi() - 3./16.*TMath::Pi() + nStripHit/8.*TMath::Pi());
  }
  return interactionPos;
}

