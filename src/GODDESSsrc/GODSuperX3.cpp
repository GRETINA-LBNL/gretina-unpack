#include "GODDESS.h"

/************************************************************/
/* Super-X3 detector functions                              */
/* Double-sided Si strip with four strips each side         */
/* p-Type is position sensitive (resistive readout)         */
/************************************************************/

superX3::superX3() {
  orrubaDet::SetNumChannels(8, 4);
  Clear();
}

superX3::superX3(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position) :
  orrubaDet(serial_Num, sector_, depth_, upstream_, position) 
{
  orrubaDet::SetNumChannels(8, 4);
  Clear();
}  

void superX3::Clear() {
  orrubaDet::Clear();
  stripPosRaw.clear();  stripPosCal.clear();
  for (Int_t i=0; i<4; i++) { stripChannelMult[i] = 0; }

  stripsP.clear(); 
  eNearRaw.clear(); eFarRaw.clear();
  eNearCal.clear(); eFarCal.clear();
  timeNear.clear(); timeFar.clear();
  
  stripsN.clear();
  eRawN.clear(); eCalN.clear();
  timeN.clear();

  eventPos.SetXYZ(0, 0, 0);

}

void superX3::ConstructBins() {
  for (Int_t i=0; i<4; i++) {
    TVector3 pStPosRefDetCenter( ((3./8.)*activeWidth) - (i*activeWidth/4.), 0., 0.);
    pStPosRefDetCenter.SetPhi(pStPosRefDetCenter.Phi() + detPos.GetRotZ());
    pStripCenterPos[i] = detPos.GetTVector3() + pStPosRefDetCenter;
    TVector3 nStPosRefDetCenter(0., 0., ((3./8.)*activeLength) - (i*activeLength/4.));
    nStripCenterPos[i] = detPos.GetTVector3() + nStPosRefDetCenter;
  }
}

/****************************************************************************/

Int_t superX3::GetStrip(Int_t channel) {
  return (channel/2);
}

/* Strip 0:  Near 1  Far 0
   Strip 1:  Near 3  Far 2
   Strip 2:  Near 5  Far 4
   Strip 3:  Near 7  Far 6  */

UShort_t superX3::GetNearChannel(UShort_t strip) { // UPDATE WITH REAL CHANNEL INFORMATION
  return 2*strip+1;
}

UShort_t superX3::GetFarChannel(UShort_t strip) { // UPDATE WITH REAL CHANNEL INFORMATION
  return 2*strip;
}

Float_t superX3::GetNearE(Bool_t calibrated) {
  std::vector<Float_t> *energies;
  if (calibrated) { energies = &eNearCal; } 
  else { energies = &eNearRaw; }
  Float_t nearE = 0.0;
  for (UInt_t i=0; i<energies->size(); i++) {
    nearE += energies->at(i);
  }
  return nearE;
}

Float_t superX3::GetFarE(Bool_t calibrated) {
  std::vector<Float_t> *energies;
  if (calibrated) { energies = &eFarCal; } 
  else { energies = &eFarRaw; }
  Float_t farE = 0.0;
  for (UInt_t i=0; i<energies->size(); i++) {
    farE += energies->at(i);
  }
  return farE;
}

void superX3::UpdatePosition(UShort_t strip) {
  // Compute the near and far channels
  UShort_t nearChannel = GetNearChannel(strip);
  UShort_t farChannel = GetFarChannel(strip);
  if (!ChannelHit(nearChannel, kFALSE) || !ChannelHit(farChannel, kFALSE)) { return; }
  Float_t nearEnergy = GetCalEValue(nearChannel, kFALSE);
  Float_t farEnergy = GetCalEValue(farChannel, kFALSE);
  if ((nearEnergy == 0 || farEnergy == 0) || (nearEnergy < 0 || farEnergy < 0)) { return; }
  Float_t stripEnergy = nearEnergy + farEnergy;

  // Compute resistive strip position by Pos = (N-F)/E
  Float_t stripPosRaw_ = (nearEnergy - farEnergy) / stripEnergy;
  Float_t stripPosCal_ = stripPosRaw_; // USE CALIBRATION PARAMETERS!
  stripPosRaw[strip] += stripPosRaw_;
  stripPosCal[strip] += stripPosCal_;
  
  // Compute the event position
  Float_t xValue = (pStripEdgePos[strip].X() + pStripEdgePos[strip+1].X())/2.;
  Float_t yValue = (pStripEdgePos[strip].Y() + pStripEdgePos[strip+1].Y())/2.;
  eventPos.SetXYZ(xValue, yValue, stripPosCal_);
}

/****************************************************************************/

void superX3::SetStripPosCalParameters(Int_t strip, std::vector<Float_t> pars) {
  if (!ValidStrip(strip)) { 
    cout << "ERROR (superX3::SetStripPosCalParameters): Invalid strip " << strip << "!\n";
    return; 
  }
  if (pars.size() != 2) { 
    cout << "ERROR (superX3::SetStripPosCalParameters): Expected 2 position calibration"
	 << " parameters, got " << pars.size() << "!\n";
    return;
  }
  parStripPosCal[strip] = pars;
}

void superX3::SetStripECalParameters(Int_t strip, std::vector<Float_t> pars) {
  if (!ValidStrip(strip)) { 
    cout << "ERROR (superX3::SetStripECalParameters): Invalid strip " << strip << "!\n";
    return; 
  }
  parStripECal[strip] = pars;
}

/****************************************************************************/

Int_t superX3::GetChannelMult(Bool_t calibrated) {
  if (calibrated) { return eNearCal.size() + eCalN.size(); }
  return eNearRaw.size() + eRawN.size();
}

Int_t superX3::GetChannelMult(Bool_t nType, Bool_t calibrated) {
  if (nType) {
    if (calibrated) { return eCalN.size(); }
    return eRawN.size();
  } else {
    if (calibrated) { return eNearCal.size(); }
    return eNearRaw.size(); 
  }
}

std::vector<Float_t> superX3::GetResE(Bool_t calibrated) {
  std::vector<Float_t> resE;
  resE.clear();
  
  std::vector<Float_t>* toSumNear = 0;
  std::vector<Float_t>* toSumFar = 0;
  
  if (calibrated) { 
    toSumNear = &eNearCal;
    toSumFar = &eFarCal;
  } else {
    toSumNear = &eNearRaw;
    toSumFar = &eFarRaw;
  }

  if (toSumNear->size() == toSumFar->size()) {
    for (UInt_t i=0; i<toSumNear->size(); i++) {
      resE.push_back(toSumNear->at(i)+toSumFar->at(i));
    }
  }
  return resE;

}

/****************************************************************************/

void superX3::LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds) {
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

void superX3::SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, Int_t ignoreThresh) {
  if (!ValidChannel(detChannel, nType)) {
    char type = 'p';
    if (nType) { type = 'n'; }
    cout <<"ERROR: Unable to set raw value for SuperX3 " << serialNum << " " << (nType ? 'n':'p') <<  "-type channel." << endl;
  }
  orrubaDet::SetRawEValue(detChannel, nType, rawValue, ignoreThresh);
}


void superX3::SortAndCalibrate(Bool_t doCalibrate) {
  //cout << "------------------------------------- ENTERING SUPER X3 SORT AND CALIBRATE" << endl;
  orrubaDet::valueMap ePMap = GetRawE(kFALSE);
  orrubaDet::valueMap eNMap = GetRawE(kTRUE);

  orrubaDet::timeMap tsPMap = GetTSmap(kFALSE);
  orrubaDet::timeMap tsNMap = GetTSmap(kTRUE);
  
  std::vector<Int_t> alreadyDone;
  alreadyDone.clear();

  for (std::map<Short_t, Float_t>::iterator chItr = ePMap.begin(); chItr != ePMap.end(); ++chItr) {
    Int_t st_ = superX3::GetStrip(chItr->first);

    //cout << "GetStrip st_: " << st_ << endl;

    if (std::find(alreadyDone.begin(), alreadyDone.end(), st_) != alreadyDone.end()) continue;
    alreadyDone.push_back(st_);
    
    Int_t nearStrip = GetNearChannel(st_);  Int_t farStrip = GetFarChannel(st_);
    Float_t eNear = 0.0, eFar = 0.0;
    std::map<Short_t, Float_t>::iterator nearItr = ePMap.find(nearStrip);
    std::map<Short_t, Float_t>::iterator farItr = ePMap.find(farStrip);

    eNear = (nearItr != ePMap.end()) ? nearItr->second : 0.0;
    eFar = (farItr != ePMap.end()) ? farItr->second : 0.0;

    if (eNear > 0.0 && eFar > 0.0) {
      stripsP.push_back(st_);
      timeNear.push_back(tsPMap[nearStrip]);
      timeFar.push_back(tsPMap[farStrip]);
      eNearRaw.push_back(eNear);
      eFarRaw.push_back(eFar);

      if(doCalibrate) { //cout << "if(doCalibrate)" << endl;
	std::vector< std::vector<Float_t> > calPar = GetECalParameters(kFALSE);
	eNear = eNear*calPar[nearStrip][1] + calPar[nearStrip][0];
	eFar = eFar*calPar[farStrip][1] + calPar[farStrip][0];
	std::vector<Float_t> *rCalPar = GetResStripParCal();
	eNear = eNear*rCalPar[st_][1] + rCalPar[st_][0];	
	eFar = eFar*rCalPar[st_][1] + rCalPar[st_][0];
	
	if (eNear > 0.0 && eFar > 0.0) {
	  eNearCal.push_back(eNear);
	  eFarCal.push_back(eFar);
	  //cout << "eNearCal: " << eNearCal.back() << endl;
	  //cout << "eFarCal:  " << eFarCal.back() << endl;
	}
      }
    }
  }
  
  for (std::map<Short_t, Float_t>::iterator chItr = eNMap.begin(); chItr != eNMap.end(); ++chItr) {
    Int_t st_ = chItr->first; //cout << "Int_t st_: " << st_ << endl;
    Float_t e_ = chItr->second;
    
    stripsN.push_back(st_);
    eRawN.push_back(e_);
    timeN.push_back(tsNMap[st_]);
    
    if(doCalibrate) {
      std::vector< std::vector<Float_t> > calPar = GetECalParameters(kTRUE);
      eCalN.push_back(e_*calPar[st_][1] + calPar[st_][0]);
      // cout << "eCalN: " << eCalN.back() << endl;
    }
  }
  // cout << "-------------------------------------- EXITING SUPER X3 SORT AND CALIBRATE" << endl;
}

Float_t superX3::GetESum(Bool_t nType, Bool_t calibrated) {
  Float_t eSum = 0;
  vector<Float_t>* toSum;
  Bool_t deletePtr = kFALSE;

  if (nType && calibrated) { toSum = &eCalN; }
  else if (nType && !calibrated) { toSum = &eRawN; }
  else if (!nType) { 
    toSum = new std::vector<Float_t>;
    *toSum = GetResE(calibrated);
    deletePtr = kTRUE;
  } else { return 0.0; }

  for (UInt_t i=0; i<toSum->size(); i++) {
    eSum += toSum->at(i);
  }
  
  if (deletePtr) delete toSum;
  return eSum;
      
}

void superX3::GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,
			    Int_t* stripMaxN, uint64_t* timestampMaxN,
			    Bool_t calibrated) {
  std::vector<Float_t> *energiesN_;
  std::vector<Float_t> *energiesNear_;
  std::vector<Float_t> *energiesFar_;

  // cout << "superX3::GetMaxHitInfo" << endl;


  if (calibrated) {
    energiesN_ = &eCalN;  energiesNear_ = &eNearCal;  energiesFar_ = &eFarCal;
  } else {
    energiesN_ = &eRawN;  energiesNear_ = &eNearRaw;  energiesFar_ = &eFarRaw;
  }
  // cout << "energiesN_" << energiesN_ << endl;

  if (stripMaxP != nullptr) { *stripMaxP = -1; }
  if (stripMaxN != nullptr) { *stripMaxN = -1; }

  Float_t eMax = 0;
  for (UInt_t i=0; i<energiesNear_->size(); i++) {
    if (energiesNear_->at(i) + energiesFar_->at(i) > eMax) {
      if (stripMaxP != nullptr) *stripMaxP = stripsP.at(i);
      eMax = energiesNear_->at(i) + energiesFar_->at(i);
      if (timestampMaxP != nullptr) *timestampMaxP = (timeNear.at(i) + timeFar.at(i))/2.;
    }
  }

  eMax = 0;
  for (UInt_t i=0; i<energiesN_->size(); i++) {
    if (energiesN_->at(i) > eMax) {
      if (stripMaxN != nullptr) *stripMaxN = stripsN.at(i);
      eMax = energiesN_->at(i);
      if (timestampMaxN != nullptr) *timestampMaxN = timeN.at(i);
    }
  }

}

TVector3 superX3::GetEventPosition(Bool_t calibrated) {

  Int_t pStripHit, nStripHit;
  GetMaxHitInfo(&pStripHit, nullptr, &nStripHit, nullptr, calibrated);
  
  Float_t eNear, eFar;
  eNear = GetNearE(calibrated);
  eFar = GetFarE(calibrated);

  TVector3 interactionPos;
  
  if (pStripHit >= 0 && pStripHit < GetNumChannels(kFALSE) && eNear >= 0. && eFar >= 0.) {
    Float_t reCenter = (parStripPosCal[pStripHit].at(1) + parStripPosCal[pStripHit].at(0))/2.;
    Float_t normalize = parStripPosCal[pStripHit].at(1) - parStripPosCal[pStripHit].at(0);
    normalize = (normalize < 0.01) ? 1 : normalize; 
    Float_t zRes = (((eNear - eFar)/(eNear + eFar)) - reCenter)/normalize;
    if (!upStream) zRes *= -1;
    TVector3 zResPos(0., 0., zRes*activeLength);
    interactionPos = pStripCenterPos[pStripHit] + zResPos;
  } else {
    interactionPos.SetXYZ(0,0,0);
  }
  return interactionPos;
  
}
