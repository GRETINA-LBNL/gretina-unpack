#include "GODDESS.h"

ionChamber::ionChamber(Int_t numAnode_, Int_t numWireX_, Int_t numWireY_, Int_t numDE_, Int_t numEres_) :
  numAnode(numAnode_), numWireX(numWireX_), numWireY(numWireY_), numDE(numDE_), numEres(numEres_) 
{ ; }

void ionChamber::Clear() {
  anodeRaw.clear(); anodeCal.clear();
  wireX.clear(); wireY.clear();
  dE = 0.0;  resE = 0.0;  E = 0.0;
}

/************************************************************/

Int_t ionChamber::GetNumChannels(Int_t chType) {
  if (chType==1) { // Anode channels
    return numAnode; 
  } else if (chType == 2) { // X position channels
    return numWireX;
  } 
  return numWireY;
}

Bool_t ionChamber::ValidChannel(UInt_t channel, Int_t chType) { 
  size_t size = numAnode;
  if (chType == 2) { size = numWireX; }
  else if (chType == 3) { size = numWireY; }
  if (channel >= size || channel < 0) { 
    cout << "ERROR (ionChamber::ValidChannel): Channel specified "
	 << channel << " is not valid\n";
    return kFALSE;
  }
  return kTRUE;
}

/************************************************************/

void ionChamber::SetRawValue(UInt_t detChannel, Int_t chType, UInt_t rawValue, Int_t ignoreThresh) {

  if (!ValidChannel(detChannel, chType)) { return; }
  UInt_t threshold = 0;
  if (chType == 1) { // Anode 
    if (detChannel < threshAnode.size()) { 
      threshold = threshAnode.at(detChannel); 
    }
    if (ignoreThresh!=0 || (ignoreThresh == 0 && rawValue > threshold)) {
      anodeRaw[detChannel] = rawValue;
    }
  } else if (chType == 2) { // X position signals
    if (detChannel < threshX.size()) { 
      threshold = threshX.at(detChannel); 
    }
    if (ignoreThresh!=0 || (ignoreThresh == 0 && rawValue > threshold)) {
      wireX[detChannel] = rawValue;
    }
  } else if (chType == 3) { // Y position signals
    if (detChannel < threshY.size()) { 
      threshold = threshY.at(detChannel); 
    }
    if (ignoreThresh!=0 || (ignoreThresh == 0 && rawValue > threshold)) {
      wireY[detChannel] = rawValue;
    }
  }

}

/************************************************************/

void ionChamber::SetChannelMap(std::map<Short_t, Short_t> channelMap, Int_t chType) {
  if (chType == 1) { mapAnode = channelMap; }
  else if (chType == 2) { mapX = channelMap; }
  else if (chType == 3) { mapY = channelMap; }
}

std::map<Short_t, Short_t> ionChamber::GetChannelMap(Int_t chType) {
  if (chType == 1) { return mapAnode; }
  else if (chType == 2) { return mapX; }
  return mapY; 
}

/************************************************************/

Bool_t ionChamber::SetThresholds(std::vector<Int_t> thresholds, Int_t chType) {
  if (chType == 1) { 
    if (thresholds.size() != numAnode) { 
      printf("ERROR: Size of vector for thresholds (%lu) not as expected (%d).\n", thresholds.size(), numAnode);
      return kFALSE;
    }
    threshAnode = thresholds;
  } else if (chType == 2) {
    if (thresholds.size() != numWireX) { 
      printf("ERROR: Size of vector for thresholds (%lu) not as expected (%d).\n", thresholds.size(), numWireX);
      return kFALSE;
    }    
    threshX = thresholds; 
  } else if (chType == 3) {
    if (thresholds.size() != numWireY) { 
      printf("ERROR: Size of vector for thresholds (%lu) not as expected (%d).\n", thresholds.size(), numWireY);
      return kFALSE;
    }
    threshY = thresholds;
  }
  return kTRUE;
}

std::vector<Int_t> ionChamber::GetThresholds(Int_t chType) {
  if (chType == 1) { return threshAnode; }
  else if (chType == 2) { return threshX; }
  return threshY; 
}

/************************************************************/

void ionChamber::LoadEvent(goddessEvent *ev, Bool_t ignoreThresholds) {
  std::map<Short_t, Short_t>::iterator itr;
  std::map<Short_t, Short_t> mapA, mapX, mapY;
  mapA = GetChannelMap(1);
  mapX = GetChannelMap(2);
  mapY = GetChannelMap(3);
  
  Short_t beginA, beginX, beginY;
  Short_t endA, endX, endY;
  itr = mapA.begin(); beginA = itr->second;
  itr = mapA.end(); itr--; endA = itr->second;
  itr = mapX.begin(); beginX = itr->second;
  itr = mapX.end(); itr--; endX = itr->second;
  itr = mapY.begin(); beginY = itr->second;
  itr = mapY.end(); itr--; endY = itr->second;
  
  for (Int_t i=0; i<ev->channels.size(); i++) {
    if (ev->channels[i] >= beginA && ev->channels[i] <= endA) {
      for (itr = mapA.begin(); itr!=mapA.end(); ++itr) {
	if (ev->channels[i] == itr->second) {
	  SetRawValue(itr->first-1, 1, ev->values[i], ignoreThresholds); 
	  break;
	}
      }
    } else if (ev->channels[i] >= beginX && ev->channels[i] <= endX) {
      for (itr = mapX.begin(); itr!=mapX.end(); ++itr) {
	if (ev->channels[i] == itr->second) {
	  SetRawValue(itr->first-1, 2, ev->values[i], ignoreThresholds); 
	  break;
	}
      }
    } else if (ev->channels[i] >= beginY && ev->channels[i] <= endY) {
      for (itr = mapY.begin(); itr!=mapY.end(); ++itr) {
	if (ev->channels[i] == itr->second) {
	  SetRawValue(itr->first-1, 3, ev->values[i], ignoreThresholds); 
	  break;
	}
      }
    }
  }
}
 
Int_t ionChamber::GetMaxX() {
  std::map<Short_t, Float_t>::iterator itr;
  
  Float_t maxSignal = 0;
  Int_t maxCh = 0;
  for (itr = wireX.begin(); itr!=wireX.end(); itr++) { 
    if (itr->second > maxSignal) {
      maxCh = itr->first; 
    }
  }
  return maxCh;
}

Int_t ionChamber::GetMaxY() {
  std::map<Short_t, Float_t>::iterator itr;
  
  Float_t maxSignal = 0;
  Int_t maxCh = 0;
  for (itr = wireY.begin(); itr!=wireY.end(); itr++) { 
    if (itr->second > maxSignal) {
      maxCh = itr->first; 
    }
  }
  return maxCh;
}
