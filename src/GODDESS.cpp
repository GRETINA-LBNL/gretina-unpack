#include "GODDESS.h"

orrubaDet::orrubaDet() :
  numNtype (0), numPtype(0) 
{ 
  Clear();
}

orrubaDet::orrubaDet(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_) :
  numNtype (0), numPtype(0), serialNum(serial_Num), sector(sector_), depth(depth_), upStream(upstream_)
{
  Clear();
  SetPosID();
}

void orrubaDet::Clear() {
  eRawPmap.clear();  eRawNmap.clear();
  eCalPmap.clear();  eCalNmap.clear();
  timePmap.clear();  timeNmap.clear();
}

/************************************************************/

void orrubaDet::SetDetector(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_) {
  serialNum = serial_Num;
  sector = sector_;
  depth = depth_;
  upStream = upstream_;

  SetPosID();
}

void orrubaDet::SetPosID() {
  posID.clear();
  if (upStream) { posID.append("U"); }
  else { posID.append("D"); }
  std::stringstream sect;  sect << sector;
  posID.append(sect.str());
  posID.append("_");
  if (depth == 0) { posID.append("dE"); }
  else if (depth == 1) { posID.append("E1"); }
  else if (depth == 2) { posID.append("E2"); }
}

/************************************************************/

void orrubaDet::SetNumChannels(Int_t pType, Int_t nType) {
  numPtype = pType;  numNtype = nType;
  parECalP.resize(pType);
  parECalN.resize(nType);
}

Int_t orrubaDet::GetNumChannels(Bool_t nType) {
  if (nType) { return numNtype; } 
  return numPtype;
}

Bool_t orrubaDet::ValidChannel(UInt_t channel, Bool_t nType) {
  size_t size = numPtype;
  if (nType) { size = numNtype; }
  if (channel >= size) {
    printf("ERROR: Contact specified, %u, is invalid!\n", channel);
    return kFALSE; 
  }
  return kTRUE;
}

Bool_t orrubaDet::ChannelHit(UInt_t channel, Bool_t nType) {
  valueMap *map = &eRawPmap;
  if (nType) { map = &eRawNmap; }
  if (map->find(channel) == map->end()) { return kFALSE; }
  return kTRUE;
}

void orrubaDet::SetChannelMap(std::map<Short_t, Short_t> channelMap, Bool_t nType) {
  if (nType) { chMapN = channelMap; } 
  else { chMapP = channelMap; }
}

std::map<Short_t, Short_t> orrubaDet::GetChannelMap(Bool_t nType) {
  if (nType) { return chMapN; }
  return chMapP;
}

/************************************************************/

void orrubaDet::SetRawEValue(UInt_t detChannel, UInt_t rawValue, Int_t ignoreThresh) {
  if (detChannel < numPtype) { SetRawEValue(detChannel, orrubaDet::pTYPE, rawValue, ignoreThresh); }
  else if (detChannel < numPtype+numNtype) {
    SetRawEValue(detChannel, orrubaDet::nTYPE, rawValue, ignoreThresh); 
  } else { printf("ERROR: Cannot SetRawValue for invalid channel %d\n", detChannel); }
}

void orrubaDet::SetRawEValue(UInt_t detChannel, Bool_t nType, UInt_t rawValue, Int_t ignoreThresh) {
  if (!ValidChannel(detChannel, nType)) { return; }
  UInt_t threshold = 0;
  if (nType) { 
    if (detChannel < threshN.size()) { threshold = threshN.at(detChannel); }
    if (ignoreThresh!=0 || (ignoreThresh==0 && rawValue > threshold)) {
      eRawNmap[detChannel] = rawValue;
    }
  } else {
    if (detChannel < threshP.size()) { threshold = threshP.at(detChannel); }
    if (ignoreThresh!=0 || (ignoreThresh==0 && rawValue > threshold)) {
      eRawPmap[detChannel] = rawValue;
    }
  }
}

orrubaDet::valueMap orrubaDet::GetRawE(Bool_t nType) {
  if (nType) { return eRawNmap; }
  return eRawPmap;
}

orrubaDet::valueMap orrubaDet::GetCalE(Bool_t nType) {
  if (nType) { return eCalNmap; }
  return eCalPmap;
}

/************************************************************/

void orrubaDet::SetTimestamp(UInt_t detChannel, Bool_t nType, uint64_t timestamp) {
  if (!ValidChannel(detChannel,nType)) { return; }
  if (nType == orrubaDet::nTYPE) { timeNmap[detChannel] = timestamp; }
  else { timePmap[detChannel] = timestamp; }
}

uint64_t orrubaDet::GetTimestamp(UInt_t detChannel, Bool_t nType) {
  uint64_t timestamp = 0xffffffffffffffff;
  for (map<Short_t, uint64_t>::iterator itr = timePmap.begin(); itr != timePmap.end(); ++itr) {
    if (timestamp > itr->second) { timestamp = itr->second; }
  }
  for (map<Short_t, uint64_t>::iterator itr = timeNmap.begin(); itr != timeNmap.end(); ++itr) {
    if (timestamp > itr->second) { timestamp = itr->second; }
  }
  return timestamp;
}

orrubaDet::timeMap orrubaDet::GetTSmap(Bool_t nType) {
  if (nType) { return timeNmap; }
  return timePmap;
}

/************************************************************/

Bool_t orrubaDet::SetECalParameters(std::vector<Float_t> par, Int_t contact, Bool_t nType) {
  if (!ValidChannel(contact,nType)) { return kFALSE; }
  if (nType) { parECalN.at(contact) = par; }
  else { parECalP.at(contact) = par; }
  return kTRUE;
}

std::vector< std::vector<Float_t> > orrubaDet::GetECalParameters(Bool_t nType) {
  if (nType) { return parECalN; }
  return parECalP;
}

Bool_t orrubaDet::SetThresholds(std::vector<Int_t> thresholds, Bool_t nType, Int_t thrSize) {
  if (thrSize == 0) { thrSize = (UInt_t)GetNumChannels(nType); }
  if (thresholds.size() != (UInt_t)thrSize) {
    printf("ERROR: Size of vector for thresholds (%d) was net as expected (%d)\n",
	   thresholds.size(), thrSize);
    return kFALSE;
  }
  if (nType == orrubaDet::nTYPE) { threshN = thresholds; }
  else { threshP = thresholds; }
  return kTRUE;
}

std::vector<Int_t> orrubaDet::GetThresholds(Bool_t nType) {
  if (nType) { return threshN; }
  return threshP;
}

Float_t orrubaDet::GetCalEnergy(Int_t contact, Bool_t nType) {
  if (!ValidChannel(contact, nType)) { return 0; }
  valueMap *eCal;
  if (nType) { eCal = &eCalNmap; } else { eCal = &eCalPmap; }
  map<Short_t, Float_t>::iterator itr = eCal->find(contact);
  if (itr != eCal->end()) { return itr->second; }
  return 0;
}

/************************************************************/
/* Super-X3 detector functions                              */
/* Double-sided Si strip with four strips each side         */
/* p-Type is position sensitive (resistive readout)         */
/************************************************************/

superX3::superX3() {
  orrubaDet::SetNumChannels(8, 4);
  Clear();
}

superX3::superX3(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_) :
  orrubaDet(serial_Num, sector_, depth_, upstream_) 
{
  orrubaDet::SetNumChannels(8, 4);
  Clear();
}  

void superX3::Clear() {
  orrubaDet::Clear();
  stripPosRaw.clear();  stripPosCal.clear();
  for (Int_t i=0; i<4; i++) { stripContactMult[i] = 0; }

  stripsP.clear(); 
  eNearRaw.clear(); eFarRaw.clear();
  eNearCal.clear(); eFarCal.clear();
  timeNear.clear(); timeFar.clear();
  
  stripsN.clear();
  eRawN.clear(); eCalN.clear();
  timeN.clear();

  eventPos.SetXYZ(0, 0, 0);

}

Int_t superX3::GetStrip(Int_t channel) {
  return (channel/2);
}

UShort_t superX3::GetNearChannel(UShort_t strip) { // UPDATE WITH REAL CHANNEL INFORMATION
  return strip > 1 ? 2*strip+1 : 2*strip;
}

UShort_t superX3::GetFarChannel(UShort_t strip) { // UPDATE WITH REAL CHANNEL INFORMATION
  return strip > 1 ? 2*strip : 2*strip + 1;
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
  orrubaDet::valueMap ePMap = GetRawE(kFALSE);
  orrubaDet::valueMap eNMap = GetRawE(kTRUE);

  orrubaDet::timeMap tsPMap = GetTSmap(kFALSE);
  orrubaDet::timeMap tsNMap = GetTSmap(kTRUE);
  
  std::vector<Int_t> alreadyDone;
  alreadyDone.clear();

  for (std::map<Short_t, Float_t>::iterator chItr = ePMap.begin(); chItr != ePMap.end(); ++chItr) {
    Int_t st_ = superX3::GetStrip(chItr->first);

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
      
      if(doCalibrate) {
	std::vector< std::vector<Float_t> > calPar = GetECalParameters(kFALSE);
	eNear = eNear*calPar[nearStrip][1] + calPar[nearStrip][0];
	eFar = eFar*calPar[farStrip][1] + calPar[farStrip][0];
	
	if (eNear > 0.0 && eFar > 0.0) {
	  eNearCal.push_back(eNear);
	  eFarCal.push_back(eFar);
	}
      }
    }
  }
  
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

void superX3::GetMaxHitInfo(Int_t* stripMaxP, uint64_t* timestampMaxP,
			    Int_t* stripMaxN, uint64_t* timestampMaxN,
			    Bool_t calibrated) {
  std::vector<Float_t> *energiesN_;
  std::vector<Float_t> *energiesNear_;
  std::vector<Float_t> *energiesFar_;
  if (calibrated) {
    energiesN_ = &eCalN;  energiesNear_ = &eNearCal;  energiesFar_ = &eFarCal;
  } else {
    energiesN_ = &eRawN;  energiesNear_ = &eNearRaw;  energiesFar_ = &eFarRaw;
  }

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
  Float_t reCenter = (parPosCal[pStripHit].at(1) + parPosCal[pStripHit].at(0))/2.;
  Float_t normalize = parPosCal[pStripHit].at(1) - parPosCal[pStripHit].at(0);
  normalize = (normalize < 0.01) ? 1 : normalize;
  Float_t zRes = (((eNear - eFar)/(eNear + eFar)) - reCenter)/normalize;
  if (!upStream) zRes *= -1;
  TVector3 zResPos(0., 0., zRes*activeLength);
  TVector3 interactionPos = pStripCenterPos[pStripHit] + zResPos;
  return interactionPos;

}

/************************************************************/
/* QQQ5 detector functions                                  */
/* Annular detector                                         */
/************************************************************/

QQQ5::QQQ5() {
  orrubaDet::SetNumChannels(32, 4);
  Clear();
}

QQQ5::QQQ5(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_) :
  orrubaDet(serial_Num, sector_, depth_, upstream_) 
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

TVector3 QQQ5::GetEventPosition(Bool_t calibrated) {

}

/************************************************************/
/* BB10 detector functions                                  */
/* Si detector                                              */
/************************************************************/

BB10::BB10() {
  orrubaDet::SetNumChannels(8, 0);
  Clear();
}

BB10::BB10(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_) :
  orrubaDet(serial_Num, sector_, depth_, upstream_) 
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

  Float_t eMax = 0.0;
  for(UInt_t i=0; i<energiesP_->size(); i++) {
    if (energiesP_->at(i) > eMax) {
      if (stripMaxP != nullptr) *stripMaxP = stripsP.at(i);
      eMax = energiesP_->at(i);
      if (timestampMaxP != nullptr) *timestampMaxP = timeP.at(i);
    }
  }

}

TVector3 BB10::GetEventPosition(Bool_t calibrated) {

}






void goddessFull::Initialize() {
  rawGoddess = new goddessEvent;
  InitializeDetectors();
}

void goddessFull::getAnalogGoddess(FILE *inf, Int_t evtLength) {
  
  Int_t siz = fread(buf, 1, evtLength, inf);
  UInt_t *p = buf;

  rawGoddess->channels.clear();
  rawGoddess->values.clear();
  rawGoddess->timestamp = 0;
  
  for (Int_t i=0; i<(evtLength/sizeof(UInt_t)); i++) { 
    UShort_t channel = *p & 0xffff;
    UShort_t value = (*p >> 16) & 0xffff;
    if (channel >= 1000 && channel <= 1003) {
      rawGoddess->timestamp |= (unsigned long long) value << (16*(246-channel));
    } else {
      rawGoddess->channels.push_back(channel);
      rawGoddess->values.push_back(value);
    }
    p++;
  }
}

void goddessFull::printAnalogRawEvent() {
  printf("-----------------------------\n");
  printf(" Raw GODDESS event: \n");
  printf("   TS = %ld\n", rawGoddess->timestamp);
  printf("   Channel\tValue\n");
  for (UInt_t i=0; i<rawGoddess->channels.size(); i++) {
    printf("    %d\t%d\n", rawGoddess->channels[i], rawGoddess->values[i]);
  }
  printf("-----------------------------\n");
}

void goddessFull::InitializeDetectors() {
  /* This is hardcoded right now -- needs to be changed to read config file at some point! */
  for (Int_t i=0; i<4; i++) {
    qqqU[i] = new QQQ5("abc", 0, 0, 1);
    qqqD[i] = new QQQ5("abc", 0, 0, 0);
  }
  for (Int_t i=0; i<8; i++) {
    bb[i] = new BB10("0x2000", 0, 0, 0);
  }
  for (Int_t i=0; i<12; i++) {
    s3U[i] = new superX3("0x3000", 0, 0, 1);
    s3D[i] = new superX3("0x4000", 0, 0, 0);
  }

  std::map<Short_t, Short_t> tempMap;

  for (Int_t i=0; i<4; i++) {
    tempMap = std::map<Short_t, Short_t>(); 
    for (Int_t j=1; j<=32; j++) {
      tempMap[j] = j+(i*32);
    }
    qqqU[i]->SetChannelMap(tempMap, 0); // p-type
  }
  for (Int_t i=0; i<4; i++) {
    tempMap.clear(); 
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+128+(i*4);
    }
    qqqU[i]->SetChannelMap(tempMap, 1); // n-type
  }

  for (Int_t i=0; i<4; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+144+(i*8);
    }
    s3U[i]->SetChannelMap(tempMap, 1); // n-type
  }
  for (Int_t i=6; i<10; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+144+((i-2)*8);
    }
    s3U[i]->SetChannelMap(tempMap, 1); // n-type
  }
  for (Int_t i=4; i<6; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+144+((i+4)*8);
    }
    s3U[i]->SetChannelMap(tempMap, 1); // n-type
  }
  for (Int_t i=10; i<12; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+144+(i*8);
    }
    s3U[i]->SetChannelMap(tempMap, 1); // n-type
  }

  for (Int_t i=0; i<12; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=8; j++) {
      tempMap[j] = j+192+(i*8);
    }
    s3U[i]->SetChannelMap(tempMap, 0); // p-type
  }
  for (Int_t i=0; i<12; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=8; j++) {
      tempMap[j] = j+288+(i*8);
    }
    s3D[i]->SetChannelMap(tempMap, 0); // p-type
  }
  for (Int_t i=0; i<12; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+384+(i*4);
    }
    s3D[i]->SetChannelMap(tempMap, 1); // n-type
  }
  for (Int_t i=0; i<8; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=8; j++) {
      tempMap[j] = j+432+(i*8);
    }
    bb[i]->SetChannelMap(tempMap, 0); // p-type
  }
  for (Int_t i=0; i<4; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=32; j++) {
      tempMap[j] = j+496+(i*32);
    }
    qqqD[i]->SetChannelMap(tempMap, 0); // p-type
  }
  for (Int_t i=0; i<4; i++) {
    tempMap.clear();
    for (Int_t j=1; j<=4; j++) {
      tempMap[j] = j+625+(i*4);
    }
    qqqD[i]->SetChannelMap(tempMap, 1); // n-type
  }
}

void goddessFull::ReportDetectors() {
  std::map<Short_t, Short_t>::iterator itr;
  std::map<Short_t, Short_t> tempMap;

  cout << "Reporting on QQQ5 upstream detectors: " << endl;
  for (Int_t i=0; i<4; i++) {
    tempMap = qqqU[i]->GetChannelMap(0);
    cout << " QQQ5 Det. " << i+1 << ",  PosID() = " << qqqU[i]->GetPosID() << endl;
    itr = tempMap.begin();
    cout << "   " << qqqU[i]->GetNumChannels(0) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl;
    tempMap = qqqU[i]->GetChannelMap(1);
    itr = tempMap.begin();
    cout << "   " << qqqU[i]->GetNumChannels(1) << " n-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl << endl;
  }
  cout << endl;

  cout << "Reporting on SuperX3 upstream detectors: " << endl;
  for (Int_t i=0; i<12; i++) {
    tempMap = s3U[i]->GetChannelMap(0);
    cout << " SuperX3 Det. " << i+1 << ",  PosID() = " << s3U[i]->GetPosID() << endl;
    itr = tempMap.begin();
    cout << "   " << s3U[i]->GetNumChannels(0) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl;
    tempMap = s3U[i]->GetChannelMap(1);
    itr = tempMap.begin();
    cout << "   " << s3U[i]->GetNumChannels(1) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl << endl;
  }
  cout << endl;

  cout << "Reporting on BB10 detectors:" << endl;
  for (Int_t i=0; i<8; i++) {
    tempMap = bb[i]->GetChannelMap(0);
    cout << " BB10 Det. " << i+1 << ",  PosID() = " << bb[i]->GetPosID() << endl;
    itr = tempMap.begin();
    cout << "   " << bb[i]->GetNumChannels(0) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl << endl;
  }
  cout << endl;

  cout << "Reporting on SuperX3 downstream detectors: " << endl;
  for (Int_t i=0; i<12; i++) {
    tempMap = s3D[i]->GetChannelMap(0);
    cout << " SuperX3 Det. " << i+1 << ",  PosID() = " << s3D[i]->GetPosID() << endl;
    itr = tempMap.begin();
    cout << "   " << s3D[i]->GetNumChannels(0) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl;
    tempMap = s3D[i]->GetChannelMap(1);
    itr = tempMap.begin();
    cout << "   " << s3D[i]->GetNumChannels(1) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl << endl;
  }
  cout << endl;

  cout << "Reporting on QQQ5 downstream detectors: " << endl;
  for (Int_t i=0; i<4; i++) {
    tempMap = qqqD[i]->GetChannelMap(0);
    cout << " QQQ5 Det. " << i+1 << ",  PosID() = " << qqqD[i]->GetPosID() << endl;
    itr = tempMap.begin();
    cout << "   " << qqqD[i]->GetNumChannels(0) << " p-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl;
    tempMap = qqqD[i]->GetChannelMap(1);
    itr = tempMap.begin(); 
    cout << "   " << qqqD[i]->GetNumChannels(1) << " n-type channels: DAQ " << itr->second << " to ";
    itr = tempMap.end(); itr--;
    cout << itr->second << endl << endl;
  }


  

}
