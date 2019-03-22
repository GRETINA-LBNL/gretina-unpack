#include "GODDESS.h"
//
//
solidVector::solidVector(Double_t x_, Double_t y_, Double_t z_, Double_t rot_Z, Double_t rot_Phi):
  TVector3(x_, y_, z_), rotZ(rot_Z), rotPhi(rot_Phi) 
{ ; }

orrubaDet::orrubaDet() :
  numNtype (0), numPtype(0) 
{ 
  Clear();
}

orrubaDet::orrubaDet(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position) :
  numNtype (0), numPtype(0), serialNum(serial_Num), sector(sector_), depth(depth_), upStream(upstream_), detPos(position)
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

void orrubaDet::SetDetector(std::string serial_Num, UShort_t sector_, UShort_t depth_, Bool_t upstream_, solidVector position) {
  serialNum = serial_Num;
  sector = sector_;
  depth = depth_;
  upStream = upstream_;
  detPos = position;

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

void orrubaDet::SetChannelMap(std::map<Short_t, Short_t> channelMap, Bool_t nType) {
  if (nType) { chMapN = channelMap; } 
  else { chMapP = channelMap; }
}

std::map<Short_t, Short_t> orrubaDet::GetChannelMap(Bool_t nType) {
  if (nType) { return chMapN; }
  return chMapP;
}

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

Float_t orrubaDet::GetRawEValue(UInt_t detChannel, Bool_t nType) {
  if (nType) { return eRawNmap[detChannel]; }
  return eRawPmap[detChannel];
}

Float_t orrubaDet::GetCalEValue(UInt_t detChannel, Bool_t nType) {
  if (nType) { return eCalNmap[detChannel]; }
  return eCalPmap[detChannel];
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
    printf("ERROR: Size of vector for thresholds (%lu) was not as expected (%d)\n",
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

/************************************************************/
/* GODDESS Tree Output                                      */
/************************************************************/

void goddessOut::Clear() {
  siID.clear();
  siType.clear();
  si.clear();
}


/************************************************************/
/* Overall GODDESS                                          */
/************************************************************/

void goddessFull::Initialize() {
  rawGoddess = new goddessEvent;
  eventOut = new goddessOut();
  eventOut->Clear();
  InitializeDetectors();
}

void goddessFull::Clear() {
  for (Int_t i=0; i<qqs.size(); i++) {
    qqs[i].Clear();
  }
  for (Int_t i=0; i<s3s.size(); i++) {
    s3s[i].Clear();
  }
  for (Int_t i=0; i<bbs.size(); i++) {
    bbs[i].Clear();
  }
  for (Int_t i=0; i<ics.size(); i++) {
    ics[i].Clear();
  }
  eventOut->Clear();
  rawGoddess->timestamp = 0;
  rawGoddess->channels.clear();
  rawGoddess->values.clear();
  ts = 0; /*2019-02-15 CMC test*/ 
}

void goddessFull::ResetOutput(detectorOUT *det) { 
  det->depth = 0;
  det->pStrip = -1;
  det->nStrip = -1;
  det->pECal = -1.0;
  det->nECal = -1.0;
  det->evPos.SetXYZ(-65000,-65000,-65000);
  //sx3 
  det->eNearCal = -1.0;
  det->eFarCal = -1.0;
  // det->nearStrip = -1;
  //det->farStrip = -1;

}

Bool_t goddessFull::ParseID(std::string id, Short_t &sector, Short_t &depth, Bool_t &upStream) {
  std::string subStr = id.substr(0, 1);
  if (subStr == "U") { upStream = kTRUE; }
  else if (subStr == "D") { upStream = kFALSE; }
  else {
    cout << "ERROR (goddessFull::ParseID): Unexpected character " << subStr 
	 << " in Up/Down stream position of id!\n";
    return kFALSE;
  }
  subStr = id.substr(1, id.length()-4);
  if (subStr == "A") { sector = 0; }
  else if (subStr == "B") { sector = 1; }
  else if (subStr == "C") { sector = 2; }
  else if (subStr == "D") { sector = 3; }
  else if (subStr.find_first_of("0123456789") != std::string::npos) {
    sector = atoi(subStr.c_str());
    if (sector > 11) { 
      cout << "ERROR (goddessFull::ParseID): Unexpected sector " << sector 
	   << " in id!\n";
      return kFALSE;
    }
  } else {
    cout << "ERROR (goddessFull::ParseID): Unexpected character " << subStr 
	 << " in sector position of id!\n";
    return kFALSE;
  }
  subStr = id.substr(id.length()-2, 2);
  if (subStr == "dE") { depth = 0; }
  else if (subStr == "E1") { depth = 1; }
  else if (subStr == "E2") { depth = 2; }
  else {
    cout << "ERROR (goddessFull::ParseID): Unexpected character " << subStr 
	 << " in depth position of id!\n";
    return kFALSE;
  }
  return kTRUE;
}

void goddessFull::ReadPositions(TString filename) {
  std::ifstream geomFile(filename.Data(), std::ifstream::in);
  if (!geomFile.good()) {
    cout << "ERROR (goddessFull::ReadPositions): Unable to read geometry file: "
	 << filename.Data() << "!\n";
  }
  cout << "Initializing GODDESS geometry from "  << filename.Data() << ".\n";
  
  string readLine;
  while (getline(geomFile, readLine)) {
    if (readLine.find('#') != std::string::npos) continue;
    //if (readLine.find('-') != std::string::npos) continue;
    if (readLine.empty()) continue;
    size_t firstCharPos = readLine.find_first_not_of(" \b\n\t\v");
    size_t equalPos = readLine.find('=');
    string geomKey = readLine.substr(firstCharPos, equalPos-firstCharPos);
    size_t lastCharPos = geomKey.find_last_not_of(" \b\n\t\v");
    geomKey = geomKey.substr(0, lastCharPos+1);
    if(geomKey.find("Offset") != string::npos) {
      string valStr = readLine.substr(equalPos + 1);
      size_t openPar = valStr.find("(");
      size_t closePar = valStr.find(")");
      size_t comma = valStr.find(",");
      if (openPar == string::npos || closePar == string::npos || comma == string::npos) {
	cout << "ERROR (goddessFull::ReadPositions): Geometry file is corrupted.\n";
	return;
      }
      geomInfo[geomKey+ " X"] = atof(valStr.substr(openPar+1, comma-(openPar+1)).c_str());
      cout << "Found parameter \"" << geomKey + " X" << "\" : " 
	   << geomInfo[geomKey + " X"] << endl;
      size_t nextComma = valStr.find(",", comma+1);
      if (nextComma == string::npos) {
	cout << "ERROR (goddessFull::ReadPositions): Geometry file is corrupted.\n";
	return;
      }
      geomInfo[geomKey+ " Y"] = atof(valStr.substr(comma+1, nextComma-(comma+1)).c_str());
      if (1) { cout << "Found parameter \"" << geomKey + " Y" << "\" : " 
		    << geomInfo[geomKey + " Y"] << endl; }
      geomInfo[geomKey+ " Z"] = atof(valStr.substr(nextComma+1, closePar-(nextComma+1)).c_str());
      if (1) { cout << "Found parameter \"" << geomKey + " Z" << "\" : " 
		    << geomInfo[geomKey + " Z"] << endl; }
    } else {
      geomInfo[geomKey] = atof(readLine.substr(equalPos+1).c_str());
      if (1) { cout << "Found parameter \"" << geomKey << "\" : "
		    << geomInfo[geomKey] << endl; }
    }
  }
 
  return;
}

solidVector goddessFull::GetPosVector(std::string type, Short_t sector, Short_t depth, Bool_t upStream) {
  Float_t barrelRadius = geomInfo["Barrel Radius"];
  Float_t barrelLength = geomInfo["Barrel Length"];
  Float_t sX3ActiveLength = geomInfo["SuperX3 Active Length"];
  Float_t qqq5FirstStripWidth = geomInfo["QQQ5 First Strip Width"];
  Float_t qqq5DeltaPitch = geomInfo["QQQ5 Delta Pitch"];


  for(int j=0; j<qqs.size(); j++){
    qqs[j].firstStripWidth = qqq5FirstStripWidth;
    qqs[j].deltaPitch = qqq5DeltaPitch;
  }

  solidVector pos(0.0, 0.0, 0.0);
  TVector3 zAxis(0, 0, 1);
  
  /* Compute positions for barrel detectors */
  if (type == "superX3" || type == "BB10") {
    TVector3 layerOffset(0.0, 0.0 + (depth-1)*geomInfo["Barrel Layer Spacing"], 0);
    TVector3 superX3Offset(geomInfo["SuperX3 Offset X"], geomInfo["SuperX3 Offset Y"],
			   geomInfo["SuperX3 Offset Z"]);
    TVector3 refSX3D_sect0 = TVector3(0, barrelRadius, sX3ActiveLength/2.) + superX3Offset + layerOffset;
    
    pos.SetXYZ(0, 0, 1);
    pos.SetTheta(upStream ? (TMath::Pi() - refSX3D_sect0.Angle(zAxis)) : refSX3D_sect0.Angle(zAxis));
    pos.SetPhi(geomInfo["SuperX3 Ref Phi"] + sector*geomInfo["SuperX3 Delta Phi"]);
    pos.SetMag(refSX3D_sect0.Mag());
    pos.SetRotationZ(geomInfo["SuperX3 Ref Zrot"] + sector*geomInfo["SuperX3 Delta Zrot"]);
  } else if (type == "QQQ5") {
    TVector3 layerOffset(0.0, 0.0, (depth-1)*geomInfo["Endcap Layer Spacing"]);
    TVector3 qqq5Offset(geomInfo["QQQ5 Offset X"], geomInfo["QQQ5 Offset Y"],
			geomInfo["QQQ5 Offset Z"]);
    TVector3 refQQQ5D_sectA = TVector3(0, 0, barrelLength/2.) + qqq5Offset + layerOffset;

    pos.SetXYZ(0, 0, 1);
    pos.SetTheta(upStream ? (TMath::Pi() - refQQQ5D_sectA.Angle(zAxis)) : refQQQ5D_sectA.Angle(zAxis));
    pos.SetPhi(geomInfo["QQQ5 Ref Phi"] + sector*geomInfo["QQQ5 Delta Phi"]);
    pos.SetMag(refQQQ5D_sectA.Mag());
    pos.SetRotationZ(geomInfo["QQQ5 Ref Zrot"] + sector*geomInfo["QQQ5 Delta Zrot"]);
  }
  return pos; 
}

void goddessFull::getAnalogGoddess(FILE *inf, Int_t evtLength) {
  
  Int_t siz = fread(buf, 1, evtLength, inf);
  UInt_t *p = buf;

  rawGoddess->channels.clear();
  rawGoddess->values.clear();
  rawGoddess->timestamp = 0;
  
  for (Int_t i=0; i<(evtLength/sizeof(UInt_t))-1; i++) { 
    UShort_t channel = (*p - 0x8000) & 0xffff;
    UShort_t value = (*p >> 16) & 0xffff;
    if (channel >= 1000 && channel <= 1003) {
      rawGoddess->timestamp |= (unsigned long long) value << (16*(channel-1000));
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
  printf("   TS = %llu\n", rawGoddess->timestamp);
  printf("   Channel\tValue\n");
  for (UInt_t i=0; i<rawGoddess->channels.size(); i++) {
    printf("    %d\t%d\n", rawGoddess->channels[i], rawGoddess->values[i]);
  }
  printf("-----------------------------\n");
  cin.get();
}

void goddessFull::ReadConfiguration(TString filename) {
  superX3 *s3 = 0;
  BB10 *bb = 0;
  QQQ5 *qq = 0;
  ionChamber *ic = 0;

  std::ifstream inFile(filename.Data(), std::ifstream::in);

  if (!inFile.good()) {
    cout << "Error! (goddessFull::ReadConfiguration) - Unable to read configuration file: " 
	 << filename.Data() << endl;
    return;
  }

  std::string line;
  UInt_t nn = 0;
  
  while(std::getline(inFile, line)) {
    nn++;
    if (line.find('#') != std::string::npos) continue;
    if (line.empty()) continue;

    size_t wordStart = line.find_first_not_of(" \b\n\t\v");
    if (wordStart == std::string::npos) continue;
    line.erase(0, wordStart);
    std::istringstream lineStream(line);

    /* Read first line for a specific detector block */
    std::string detType, serialNum = "", id;
    if (!(lineStream >> detType)) { break; }
   
    cout << "\nRegistering " << detType << " detector:\n";
    
    if (detType == "ion") {
      std::map<Short_t, Short_t> tempMap = std::map<Short_t, Short_t>(); 
      
      Int_t aDaqCh, xDaqCh, yDaqCh = 0;
      if (!(lineStream >> aDaqCh >> xDaqCh >> yDaqCh)) { break; }
      ic = new ionChamber(2, 32, 32, 1, 1);
      for (Int_t j=1; j<=2; j++) {
	tempMap[j] = j + aDaqCh;
      }
      ic->SetChannelMap(tempMap, 1);
      tempMap = std::map<Short_t, Short_t>(); 
      for (Int_t j=1; j<=32; j++) {
	tempMap[j] = j + xDaqCh;
      }
      ic->SetChannelMap(tempMap, 2);
      tempMap = std::map<Short_t, Short_t>(); 
      for (Int_t j=1; j<=32; j++) {
	tempMap[j] = j + yDaqCh;
      }
      ic->SetChannelMap(tempMap, 3);

    } else if (detType == "liquidScint") {


    } else if (detType == "QQQ5" || detType == "BB10" || detType == "superX3") {
      Int_t pTypeDaqCh, nTypeDaqCh = 0;
      
      if (detType == "BB10") {
	if (!(lineStream >> serialNum >> id >> pTypeDaqCh)) { break; }
      } else {
	if (!(lineStream >> serialNum >> id >> pTypeDaqCh >> nTypeDaqCh)) { break; }
      }

      cout << serialNum << ", LocID: " << id << ", pDAQ: " << pTypeDaqCh;
      if (detType == "BB10") { cout << "\n"; }
      else { cout << ", nDAQ: " << nTypeDaqCh << "\n"; }

      Short_t sector, depth;  Bool_t upStream;
      if (!ParseID(id, sector, depth, upStream)) {
      	cout << "Error! (goddessFull::ReadConfiguration) - Unable to parse position ID: " 
	     << id << " of detector " << serialNum << "!\n";
	break;
      }

      cout << " Upstream: " << (upStream ? "U" : "D");
      cout << ", Sector: " << sector << ", Depth: " << depth << "\n";

      solidVector pos_ = GetPosVector(detType, sector, depth, upStream);
      cout << " Position: " << pos_.x() << ", " << pos_.y() << ", " << pos_.z()
	   << ", RotZ: " << pos_.GetRotZ()*TMath::RadToDeg()
	   << ", RotPhi: " << pos_.GetRotPhi()*TMath::RadToDeg() << "\n";
      
      std::map<Short_t, Short_t> tempMap = std::map<Short_t, Short_t>(); 

      /* Construct the detector! */
      if (detType == "superX3") {
	s3 = new superX3();
	s3->SetDetector(serialNum, sector, depth, upStream, pos_);
	for (Int_t j=1; j<=8; j++) {
	  tempMap[j] = j+pTypeDaqCh;
	}
	s3->SetChannelMap(tempMap, 0); // p-type
	tempMap = std::map<Short_t, Short_t>(); 
	for (Int_t j=1; j<=4; j++) {
	  tempMap[j] = j+nTypeDaqCh;
	}
	s3->SetChannelMap(tempMap, 1); // n-type
	s3->SetGeomParameters(geomInfo);
	s3->ConstructBins();
      } else if (detType == "BB10") {
	bb = new BB10();
	bb->SetDetector(serialNum, sector, depth, upStream, pos_);
	for (Int_t j=1; j<=8; j++) {
	  tempMap[j] = j+pTypeDaqCh;
	}
	bb->SetChannelMap(tempMap, 0); // p-type
	bb->SetGeomParameters(geomInfo);
	bb->ConstructBins();
      } else if (detType == "QQQ5") {
	qq = new QQQ5();
	qq->SetDetector(serialNum, sector, depth, upStream, pos_);
	for (Int_t j=1; j<=32; j++) {
	  tempMap[j] = j+pTypeDaqCh;
	}
	qq->SetChannelMap(tempMap, 0); // p-type
	tempMap = std::map<Short_t, Short_t>(); 
	for (Int_t j=1; j<=4; j++) {
	  tempMap[j] = j+nTypeDaqCh;
	}
	qq->SetChannelMap(tempMap, 1); // n-type
	qq->SetGeomParameters(geomInfo);
	qq->ConstructBins();
      } else {
	cout << "Error! (goddessFull::ReadConfiguration) - Unknown detector type: " 
	     << detType << "!\n";
	break;
      }
    } else { continue; }
    
    /* Read calibration information */
    std::istream::streampos prevPos = inFile.tellg();
    while (std::getline(inFile, line)) {
      nn++;
      size_t posT = line.find('#');
      if (posT != std::string::npos) { line.erase(posT); }
      size_t wordStart = line.find_first_not_of(" \t");
      if (wordStart == std::string::npos) { 
	prevPos = inFile.tellg();
	continue;
      }
      line.erase(0, wordStart);
      std::istringstream lineStream(line);
      std::string calType, subType;
      Int_t detChannel;
      if (!(lineStream >> calType >> subType)) { inFile.seekg(prevPos); nn--; break; }
      
      Bool_t posCal = kFALSE;
      Bool_t timeCal = kFALSE;
      Bool_t thresh = kFALSE;
      if (calType == "posCal") {
	if (!(detType == "superX3" && subType == "resStrip") && !(detType=="ion" && subType=="scint")) {
	  cout << "Error! (goddessFull::ReadConfiguration) - Ignoring position calibration for " 
	       << subType << "-type." << "\n";
	  prevPos = inFile.tellg();
	  continue;
	} 
	posCal = kTRUE;
      } else if (calType == "timeCal") {
	if (!(detType=="ion" && subType == "scint")) {
	  cout << "Error! (goddessFull::ReadConfiguration) - Ignoring time calibration for " 
	       << subType << "-type." << "\n";
	  prevPos = inFile.tellg();
	  continue;
	}
	timeCal = kTRUE;
      } else if (calType == "thresh") { thresh = kTRUE; }
      else if (calType != "eCal") { 
	inFile.seekg(prevPos);
	nn--;
	break;
      }
      
      Bool_t secondaryType = kFALSE;
      
      /* Check subtype is valid */
      if (detType == "ion") { 
	if (subType != "x" && subType != "anode" && subType != "y") {
	  cout << "Error! (goddessFull::ReadConfiguration) - Unknown calibration subtype: " 
	       << subType << "!\n";
	  continue;
	}
      } else {
	if (subType == "n") { secondaryType = kTRUE; }
	else if (subType != "p" && subType != "resStrip") {
	  cout << "Error! (goddessFull::ReadConfiguration) - Unknown calibration subtype: " 
	       << subType << "!\n";
	  continue;
	}
      }
      
      if (thresh) {
	Int_t threshold;
	std::vector<Int_t> thresholds;
	while (lineStream >> threshold) { thresholds.push_back(threshold); }
	if (detType == "ion") {
	  if (subType == "anode") { 
	    cout << "Setting anode thresholds (" << thresholds.size() << ").\n";
	    ic->SetThresholds(thresholds, 1);
	  } else if (subType == "x") { 
	    cout << "Setting x thresholds (" << thresholds.size() << ").\n";
	    ic->SetThresholds(thresholds, 2);
	  } else if (subType == "y") {
	    cout << "Setting y thresholds (" << thresholds.size() << ").\n";
	    ic->SetThresholds(thresholds, 3);
	  }
	} else {
	  if (secondaryType) {
	    cout << "Setting n-type thresholds (" << thresholds.size() << ").\n";
	  } else {
	    cout << "Setting p-type thresholds (" << thresholds.size() << ").\n";
	  }
	  if (detType == "BB10") { bb->SetThresholds(thresholds, secondaryType, thresholds.size()); }
	  if (detType == "QQQ5") { qq->SetThresholds(thresholds, secondaryType, thresholds.size()); }
	  if (detType == "superX3") { s3->SetThresholds(thresholds, secondaryType, thresholds.size()); }
	}
	continue;
      } else {
	if (!(lineStream >> detChannel)) {
	  cout << "Error! (goddessFull::ReadConfiguration) - No detector channel specified for calibration.\n";
	  continue;
	}
      }
      
      /* Get the calibration parameters */
      Float_t param;
      std::vector<Float_t> params;
      while (lineStream >> param) { params.push_back(param); }
      if (detType == "superX3" && subType == "resStrip") {
	if (posCal) { s3->SetStripPosCalParameters(detChannel, params); }
	else { s3->SetStripECalParameters(detChannel, params); }
      } else {
	if (detType == "superX3" && subType != "resStrip") { s3->SetECalParameters(params, detChannel, secondaryType); }
	else if (detType == "BB10") { bb->SetECalParameters(params, detChannel, secondaryType); }
	else if (detType == "QQQ5") { qq->SetECalParameters(params, detChannel, secondaryType); }
      }
      prevPos = inFile.tellg();
    }
    
    if (detType == "QQQ5") { 
      cout << "Read " << qq->GetECalParameters(0).size() << " energy calibrations (p-type).\n";
      cout << "Read " << qq->GetECalParameters(1).size() << " energy calibrations (n-type).\n";
      qqs.push_back(*qq); 
    } else if (detType == "superX3") { 
      cout << "Read " << s3->GetECalParameters(0).size() << " energy calibrations (p-type).\n";
      cout << "Read " << s3->GetECalParameters(1).size() << " energy calibrations (n-type).\n";
      s3s.push_back(*s3); 
    } else if (detType == "BB10") { 
      cout << "Read " << bb->GetECalParameters(0).size() << " energy calibrations (p-type).\n";
      bbs.push_back(*bb); 
    } else if (detType == "ion") {
      ics.push_back(*ic);
    }
  }
  
  if (inFile.good()) { 
    cout << "Error! (goddessFull::ReadConfiguration) - Failed on line " << nn << "\n";
  }
}

void goddessFull::InitializeDetectors() {
  ReadPositions("goddessConfiguration/DetGeom.dat");
  ReadConfiguration("goddessConfiguration/DetConfig.dat");
}

void goddessFull::ProcessEvent() {
  detectorOUT detOut;

  for (UInt_t i=0; i<qqs.size(); i++) {
    qqs[i].LoadEvent(rawGoddess, 0);
    qqs[i].SortAndCalibrate(kTRUE);
    ResetOutput(&detOut);
    qqs[i].GetMaxHitInfo(&detOut.pStrip, nullptr, &detOut.nStrip, nullptr);
    //  cout << "QQS stripsP.size(): " << qqs[i].stripsP.size() << endl;
    //  cout << "QQS stripsN.size(): " << qqs[i].stripsN.size() << endl;
    if (detOut.pStrip >= 0) {
      detOut.depth = qqs[i].GetDepth();
      detOut.upstream = qqs[i].GetUpStream();
      detOut.pECal = -1.0;  detOut.nECal = -1.0;
      for (Int_t j=0; j<qqs[i].stripsP.size(); j++) {
	if (qqs[i].stripsP[j] == detOut.pStrip) {
	  detOut.pECal = qqs[i].eCalP[j];
	}
      }
      for (Int_t j=0; j<qqs[i].stripsN.size(); j++) {
	if (qqs[i].stripsN[j] == detOut.nStrip) {
	  detOut.nECal = qqs[i].eCalN[j];
	}
      }
      detOut.evPos = qqs[i].GetEventPosition();
      eventOut->siID.push_back(qqs[i].GetPosID());
      eventOut->siType.push_back("QQQ5");
      eventOut->si.push_back(detOut);
    }
  }
  // cout << "i<s3s.size(): " << s3s.size() << endl;
  for (UInt_t i=0; i<s3s.size(); i++) {
    s3s[i].LoadEvent(rawGoddess, 0);
    s3s[i].SortAndCalibrate(kTRUE);//test
    ResetOutput(&detOut);
    s3s[i].GetMaxHitInfo(&detOut.pStrip, nullptr, &detOut.nStrip, nullptr, kTRUE);
    // cout << "stripsP.size(): " << s3s[i].stripsP.size() << endl;
    // cout << "stripsN.size(): " << s3s[i].stripsN.size() << endl; 
    // cout << "detOut.pStrip: " << detOut.pStrip << endl;
    // cout << "detOut.nStrip: " << detOut.nStrip << endl;
     if (detOut.pStrip >= 0) {
       detOut.depth = s3s[i].GetDepth();
       detOut.upstream = s3s[i].GetUpStream();
       detOut.pECal    = -1.0; detOut.nECal   = -1.0;
       detOut.eNearCal = -1.0; detOut.eFarCal = -1.0;
      
      for (Int_t j=0; j<s3s[i].stripsP.size(); j++) { 
	if (s3s[i].stripsP[j] == detOut.pStrip) {
	  detOut.eNearCal = s3s[i].eNearCal[j];
	  detOut.eFarCal  = s3s[i].eFarCal[j];
	  detOut.pECal = s3s[i].eNearCal[j] + s3s[i].eFarCal[j];
	}
      }
      for (Int_t j=0; j<s3s[i].stripsN.size(); j++) {
	if (s3s[i].stripsN[j] == detOut.nStrip) {
	  detOut.nECal = s3s[i].eCalN[j];
	}
      }
   
      detOut.evPos = s3s[i].GetEventPosition();// NEED POSITION CALIB PARAMETERS IN .DAT FILE
      eventOut->siID.push_back(s3s[i].GetPosID());
      eventOut->siType.push_back("SuperX3");
      eventOut->si.push_back(detOut);
     }
  } //cout << "leaving s3s loop --- " << endl;
  
  for (UInt_t i=0; i<bbs.size(); i++) {
    bbs[i].LoadEvent(rawGoddess, 0);
    bbs[i].SortAndCalibrate(kTRUE);
    ResetOutput(&detOut);
    bbs[i].GetMaxHitInfo(&detOut.pStrip, nullptr, nullptr, nullptr);
    if (detOut.pStrip >= 0) {
      detOut.nStrip = -1;
      detOut.depth = bbs[i].GetDepth();
      detOut.upstream = bbs[i].GetUpStream();
      detOut.pECal = -1.0;  detOut.nECal = -1.0;
      for (Int_t j=0; j<bbs[i].stripsP.size(); j++) {
	if (bbs[i].stripsP[j] == detOut.pStrip) {
	  detOut.pECal = bbs[i].eCalP[j];
	}
      }
      detOut.evPos = bbs[i].GetEventPosition();
      eventOut->siID.push_back(bbs[i].GetPosID());
      eventOut->siType.push_back("BB10");
      eventOut->si.push_back(detOut);
    }
  }
  for (UInt_t i=0; i<ics.size(); i++) {
    ics[i].LoadEvent(rawGoddess, 0);
    eventOut->icDE = ics[i].GetAnodeDE();
    eventOut->icE = ics[i].GetAnodeE();
    eventOut->icX = ics[i].GetMaxX();
    eventOut->icY = ics[i].GetMaxY();
  }

  if (0) {
    std::map<Short_t, Float_t> eMap;
    std::map<Short_t, Float_t>::iterator itr;
    for (Int_t i=0; i<qqs.size(); i++) {
      cout << "QQQ5 ID: " << qqs[i].GetPosID() << endl;
      eMap = qqs[i].GetRawE(0);
      if (eMap.size() > 0) {
	for (itr=eMap.begin(); itr!= eMap.end(); ++itr) {
	  cout << "QQQ5-" << i+1 << " (P) Ch. " << itr->first << ", E = " << itr->second << endl;
	}
      }
      if (eMap.size() > 0) {
	eMap = qqs[i].GetRawE(1);
	for (itr=eMap.begin(); itr!= eMap.end(); ++itr) {
	  cout << "QQQ5-" << i+1 << " (N) Ch. " << itr->first << ", E = " << itr->second << endl;
	}
      }
    }
    for (Int_t i=0; i<s3s.size(); i++) {
      eMap = s3s[i].GetRawE(0);
      if (eMap.size() > 0) {
	for (itr=eMap.begin(); itr!= eMap.end(); ++itr) {
	  cout << "SuperX3-" << i+1 << " (P) Ch. " << itr->first << ", E = " << itr->second << endl;
	}
      }
      if (eMap.size() > 0) {
	eMap = s3s[i].GetRawE(1);
	for (itr=eMap.begin(); itr!= eMap.end(); ++itr) {
	  cout << "SuperX3-" << i+1 << " (N) Ch. " << itr->first << ", E = " << itr->second << endl;
	}
      }
    }
    for (Int_t i=0; i<bbs.size(); i++) {
      eMap = bbs[i].GetRawE(0);
      if (eMap.size() > 0) {
	for (itr=eMap.begin(); itr!= eMap.end(); ++itr) {
	  cout << "BB10-" << i+1 << " (P) Ch. " << itr->first << ", E = " << itr->second << endl;
	}
      }
    }
    cin.get();
  }

  if (0) {
    cout << "-------------------------" << endl;
    for (Int_t i=0; i<qqs.size(); i++) {
      if (qqs[i].eCalP.size() > 0) {
	for (Int_t j=0; j<qqs[i].eCalP.size(); j++) {
	  cout << "QQQ5-" << i+1 << " (P) Ch. " << qqs[i].stripsP[j] << ", E = " << qqs[i].eCalP[j] << endl;
	}
      }
      if (qqs[i].eCalN.size() > 0) {
	for (Int_t j=0; j<qqs[i].eCalN.size(); j++) {
	  cout << "QQQ5-" << i+1 << " (N) Ch. " << qqs[i].stripsN[j] << ", E = " << qqs[i].eCalN[j] << endl;
	}
      }
    }
    for (Int_t i=0; i<s3s.size(); i++) {
      // if (s3s[i].eCalP.size() > 0) {
      //	for (Int_t j=0; j<s3s[i].eCalP.size(); j++) {
      //	  cout << "SuperX3-" << i+1 << " (P) Ch. " << s3s[i].stripsP[j] << ", E = " << s3s[i].eCalP[j] << endl;
      //	}
      // }
      if (s3s[i].eCalN.size() > 0) {
	for (Int_t j=0; j<s3s[i].eCalN.size(); j++) {
	  cout << "SuperX3-" << i+1 << " (N) Ch. " << s3s[i].stripsN[j] << ", E = " << s3s[i].eCalN[j] << endl;
	}
      }
    }
    for (Int_t i=0; i<bbs.size(); i++) {
      if (bbs[i].eCalP.size() > 0) {
	for (Int_t j=0; j<bbs[i].eCalP.size(); j++) {
	  cout << "BB10-" << i+1 << " (P) Ch. " << bbs[i].stripsP[j] << ", E = " << bbs[i].eCalP[j] << endl;
	}
      }
    }

  }



}
