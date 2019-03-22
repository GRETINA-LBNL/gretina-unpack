#include "LENDA-DDAS.h"

ClassImp(ddasEvent);

ClassImp(lendaChannel);
ClassImp(lendaBar);
ClassImp(lendaEvent);
ClassImp(lendaFilter);

/************************************************************/
/* ddasEvent class functions ********************************/
/************************************************************/

ddasEvent::ddasEvent() : TObject(), cdata() { ; }

ddasEvent::ddasEvent(const ddasEvent& obj) : TObject(obj), cdata() {
  /* Create new copies of the ddasChannel events. */
  for (UInt_t ui=0; ui<cdata.size(); ++ui) {
    cdata[ui] = new ddasChannel(*obj.cdata[ui]);
  } 
}

ddasEvent& ddasEvent::operator=(const ddasEvent& obj) {
  if (this!=&obj) {
    /* Create new copies of the ddasChannel events. */
    for (UInt_t ui=0; ui < cdata.size(); ++ui) {
      delete cdata[ui];
    }
    cdata.resize(obj.cdata.size());
    for (UInt_t ui=0; ui < cdata.size(); ++ui) {
      cdata[ui] = new ddasChannel(*obj.cdata[ui]);
    }  
  }
  return *this;
}

ddasEvent::~ddasEvent() {
  Reset();
}

void ddasEvent::addChannelData(ddasChannel* channel) {
  cdata.push_back(channel);
}

Double_t ddasEvent::getFirstTime() const {
  Double_t time = 0;
  if (cdata.size()>0) {  time = cdata.front()->GetTime(); }
  return time;
}

Double_t ddasEvent::getLastTime() const {
  Double_t time = 0;
  if (cdata.size()>0) { time = cdata.back()->GetTime(); }    
  return time;
}

Double_t ddasEvent::getTimeWidth() const {
  Double_t end_time   = getLastTime();
  Double_t start_time = getFirstTime();
  return (end_time-start_time);
}

void ddasEvent::Reset() {
  for (UInt_t ui=0; ui < cdata.size(); ++ui) {
    cdata[ui]->Reset();
  }
  /* Clear the array and resize it to zero */
  cdata.clear();
}

void ddasEvent::getEvent(FILE *inf, Int_t length) {
  Reset();

  Int_t deBugFn = 0;

  uint32_t sizeOfWholeThing = 0;
  uint32_t typeOfWholeThing = 0;   

  uint64_t timestamp;
  uint32_t sourceID, fragPayloadSize, barrier, sizeOfRingItem, typeOfRingItem;

  uint64_t timestampBH;
  uint32_t sizeOfBH, sourceIDBH, barrierBH;

  UInt_t bytesReadTotal = 0;
  UInt_t bytesReadInBody = 0, totalSizeOfBody = 0;

  Int_t count = 0;

  size_t siz = fread(&sizeOfWholeThing, sizeof(uint32_t), 1, inf);
  bytesReadTotal += sizeof(uint32_t);
  // fread(&typeOfWholeThing, sizeof(uint32_t), 1, inf);
  // bytesReadTotal += sizeof(uint32_t);
  typeOfWholeThing = 30;

  if(deBugFn) {
    printf("sizeofwholething = %d\n", sizeOfWholeThing);
    printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
  }

  switch(typeOfWholeThing) {
  case 30: /* Physics event */
    {  
      totalSizeOfBody = sizeOfWholeThing - 4;
      while (1) {
	size_t siz = fread(&timestamp, sizeof(uint64_t), 1, inf);
	bytesReadTotal += sizeof(uint64_t);
	bytesReadInBody += sizeof(uint64_t);

	if(deBugFn) {
	  printf("timestamp = %lu\n", timestamp);
	  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
	}

	siz = fread(&sourceID, sizeof(uint32_t), 1, inf);
	bytesReadTotal += sizeof(uint32_t);
	bytesReadInBody += sizeof(uint32_t);

	if(deBugFn) {
	  printf("sourceID = %u\n", sourceID);
	  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
	}

	siz = fread(&fragPayloadSize, sizeof(uint32_t), 1, inf);
	bytesReadTotal += sizeof(uint32_t);
	bytesReadInBody += sizeof(uint32_t);

	if(deBugFn) {
	  printf("fragPayloadSize = %u\n", fragPayloadSize);
	  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
	}

	siz = fread(&barrier, sizeof(uint32_t), 1, inf);
	bytesReadTotal += sizeof(uint32_t);
	bytesReadInBody += sizeof(uint32_t);

	if(deBugFn) {
	  printf("barrier = %u\n", barrier);
	  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
	}

	siz = fread(&sizeOfRingItem, sizeof(uint32_t), 1, inf);
	bytesReadTotal += sizeof(uint32_t);
	bytesReadInBody += sizeof(uint32_t);

	if(deBugFn) {
	  printf("sizeOfRingItem = %u\n", sizeOfRingItem);
	  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
	}

	siz = fread(&typeOfRingItem, sizeof(uint32_t), 1, inf);
	bytesReadTotal += sizeof(uint32_t);
	bytesReadInBody += sizeof(uint32_t);

	if(deBugFn) {
	  printf("typeOfRingItem = %u\n", typeOfRingItem);
	  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
	}

	Int_t NSCLDAQ11 = 1;

	switch(sourceID) {
	case 3: /* DDAS SourceID */
	  {
	    switch(typeOfRingItem) {
	    case 30: /* Physics */
	      {
		if (NSCLDAQ11) {
		  size_t siz = fread(&sizeOfBH, sizeof(uint32_t), 1, inf);
		  bytesReadTotal += sizeof(uint32_t);
		  bytesReadInBody += sizeof(uint32_t);
	
		  if(deBugFn) {
		    printf("sizeOfBH = %u\n", sizeOfBH);
		    
		    printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
		  }

		  siz = fread(&timestampBH, sizeof(uint64_t), 1, inf);
		  bytesReadTotal += sizeof(uint64_t);
		  bytesReadInBody += sizeof(uint64_t);

		  if(deBugFn) {
		    printf("timestampBH = %lu\n", timestampBH); 
		    printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
		  }

		  siz = fread(&sourceIDBH, sizeof(uint32_t), 1, inf);
		  bytesReadTotal += sizeof(uint32_t);
		  bytesReadInBody += sizeof(uint32_t);
		  
		  if(deBugFn) {
		    printf("sourceIDBH = %u\n", sourceIDBH);
		    printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
		  }

		  siz = fread(&barrierBH, sizeof(uint32_t), 1, inf);
		  bytesReadTotal += sizeof(uint32_t);
		  bytesReadInBody += sizeof(uint32_t);

		  if(deBugFn) {		  
		    printf("barrierBH = %u\n", barrierBH);
		    printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
		  }

		  count = sizeOfRingItem - 7*sizeof(uint32_t);
		} else {
		  count = sizeOfRingItem - 2*sizeof(uint32_t);
		}
		Int_t countRead = fread(buffer32, sizeof(int8_t), count, inf);
		bytesReadTotal += countRead * sizeof(int8_t);
		bytesReadInBody += countRead * sizeof(int8_t);
		
		if(deBugFn) {		
		  printf("bytesReadTotal: %u, bytesReadInBody: %u\n", bytesReadTotal, bytesReadInBody);
		}
		
		uint32_t *pt32 = buffer32;
		tempdchan = new ddasChannel();
		tempdchan->UnpackChannelData(pt32);
		addChannelData(tempdchan);
  	      }
	      break;
	    default:
	      {
		count = sizeOfRingItem - 2*sizeof(uint32_t);
		Int_t countRead = fread(buffer32, sizeof(int8_t), count, inf);
		bytesReadTotal += countRead * sizeof(int8_t);
		bytesReadInBody += countRead * sizeof(int8_t);
	      }
	      break;
	    }
	  }
	  break;
	default:
	  {
	    size_t size = fread(buffer08, sizeof(int8_t), sizeOfRingItem - 2*sizeof(uint32_t), inf);
	    bytesReadTotal += (sizeOfRingItem - 2*sizeof(uint32_t))*sizeof(int8_t);
	    bytesReadInBody += (sizeOfRingItem - 2*sizeof(uint32_t))*sizeof(int8_t);
	  }
	  break;
	}
	if (totalSizeOfBody == bytesReadInBody) {
	  break;
	}
      } /* while loop */
    }
    break;
  default:
    {
      count = sizeOfWholeThing - 2*sizeof(uint32_t);
      size_t siz = fread(buffer08, sizeof(uint8_t), count / sizeof(uint8_t), inf);
    }
    break;
  }

}

/************************************************************/
/* lendaChannel class functions *****************************/
/************************************************************/

lendaChannel::lendaChannel() {
  Clear();
}

lendaChannel::~lendaChannel() {
  Clear();
}

Bool_t lendaChannel::operator==(const lendaChannel & RHS){
  cout << "<lendaChannel::operator==()> Needs to be updated" << endl;
  return false;
}

void lendaChannel::Clear() {
  /* The vectors that hold the many times and the many energies
     needed for when this channel hold info for high rate event. */
  rSoftwareTimes.clear();
  rEnergies.clear();
  rCubicTimes.clear();
  rPulseHeights.clear();
  rSoftwareReferenceTimes.clear();
  rCubicReferenceTimes.clear();
  
  rSoftwareTimes.push_back(BADNUM);
  rEnergies.push_back(BADNUM);
  rCubicTimes.push_back(BADNUM);
  rPulseHeights.push_back(BADNUM);
  rSoftwareReferenceTimes.push_back(BADNUM);
  rCubicReferenceTimes.push_back(BADNUM);
  
  rChannel = BADNUM;
  rSlot = BADNUM;
  rGlobalID = BADNUM;
  rRefGlobalID = BADNUM;

  rOverflow = 0;
  rNumZeroCrossings = BADNUM;
  
  rCorrectedEnergy = BADNUM;
  
  rInternalEnergy = BADNUM;
  
  rFilterHeight = BADNUM;

  rTime = BADNUM;
  rCubicFitTime = BADNUM;

  rCorrectedTime = BADNUM;
  rCorrectedCubicFitTime = BADNUM;
  rCorrectedCubicTime = BADNUM;
  rCorrectedSoftTime = BADNUM;

  rTimeLow = 0;  rTimeHigh = 0;
  
  rEclockLow = 0;  rEclockHigh = 0;
  rEclockHigher = 0;  rEclock = 0;
  
  rSoftwareCFD = BADNUM;
  rCubicCFD = BADNUM;
  rCubicFitCFD = BADNUM;
  rInternalCFD = BADNUM;
  rCFDTrigBit = 0;
  
  rShortGate = BADNUM;
  rLongGate = BADNUM;

  rJentry = BADNUM;

  rRefTime = BADNUM;

  rTrace.clear();
  rFilter.clear();
  rCFD.clear();
  rUnderShoots.clear();

  rCFDResidual = BADNUM;

  rChannelName = "";
  rReferenceChannelName = "";
  rCoarseTime = BADNUM;
  rOtherTime = BADNUM;
}

void lendaChannel::SetSoftwareReferenceTimes(const vector<Double_t> & v) { 
  if (rSoftwareReferenceTimes.size() == 1 &&
      rSoftwareReferenceTimes[0] == BADNUM) {
    rSoftwareReferenceTimes = v; 
  } else {
    for (UInt_t ui=0; ui<v.size(); ui++) {
      rSoftwareReferenceTimes.push_back(v[ui]);
    }
  }
}

void lendaChannel::SetCubicReferenceTimes(const vector<Double_t> & v) { 
  if (rCubicReferenceTimes.size() == 1 &&
      rCubicReferenceTimes[0] == BADNUM) {
    rCubicReferenceTimes = v; 
  } else {
    for (UInt_t ui=0; ui<v.size(); ui++) {
      rCubicReferenceTimes.push_back(v[ui]);
    }
  }
}

/************************************************************/
/* lendaBar class functions *********************************/
/************************************************************/

lendaBar::lendaBar() {
  Clear();
}

lendaBar::~lendaBar() {
  Clear();
}

void lendaBar::Clear() {
  Tops.clear();
  Bottoms.clear();
  NumTops = 0;
  NumBottoms = 0;
  BarMultiplicity = 0;

  Dt = BADNUM;
  CubicDt = BADNUM;

  SyncDt = BADNUM;
  CubicSyncDt = BADNUM;
  
  COG = BADNUM;
  AvgEnergy = BADNUM;
  AvgPulseHeight = BADNUM;

  CorrectedCOG = BADNUM;
  CorrectedDt = BADNUM;
  CorrectedAvgEnergy = BADNUM;

  AvgT = BADNUM;
  CorrectedAvgT = BADNUM;
  SimpleEventBit = 0;

  Name = "_No_NAME_";
  BarId = BADNUM;

  TopTOF = BADNUM;
  BottomTOF = BADNUM;
  
  CorrectedTopTOF = BADNUM;
  CorrectedBottomTOF = BADNUM;

  CorrectedCubicTopTOF = BADNUM;
  CorrectedCubicBottomTOF = BADNUM;

  CorrectedCubicFitTopTOF = BADNUM;
  CorrectedCubicFitBottomTOF = BADNUM;

  CorrectedSoftTopTOF = BADNUM;
  CorrectedSoftBottomTOF = BADNUM;
}

void lendaBar::Finalize() {
  NumTops = Tops.size();
  NumBottoms = Bottoms.size();
  BarMultiplicity=NumTops + NumBottoms;
  
  if (NumTops == 1) {
    TopTOF = (Tops[0].GetTime() - Tops[0].GetReferenceTime());
    CubicTopTOF = (Tops[0].GetCubicTime() - Tops[0].GetCubicReferenceTime());
    CorrectedTopTOF = (Tops[0].GetCorrectedTime() - Tops[0].GetReferenceTime());
    CorrectedCubicTopTOF = (Tops[0].GetCorrectedCubicTime() -
			    Tops[0].GetCubicReferenceTime());
    CorrectedSoftTopTOF = (Tops[0].GetCorrectedSoftTime() -
			   Tops[0].GetSoftwareReferenceTime());
  }
  
  if (NumBottoms == 1) {
    BottomTOF = (Bottoms[0].GetTime() - Bottoms[0].GetReferenceTime());
    CubicBottomTOF = (Bottoms[0].GetCubicTime() -
		      Bottoms[0].GetCubicReferenceTime());
    CorrectedBottomTOF = (Bottoms[0].GetCorrectedTime() -
			  Bottoms[0].GetReferenceTime());
    CorrectedCubicBottomTOF = (Bottoms[0].GetCorrectedCubicTime() -
			       Bottoms[0].GetCubicReferenceTime());
    CorrectedSoftBottomTOF = (Bottoms[0].GetCorrectedSoftTime() -
			      Bottoms[0].GetSoftwareReferenceTime());
  }

  /* For simple events -- ones with only 1 Top and 1 Bottom --
     calculate some nice quantities */
  if (NumTops==1 && NumBottoms==1) {
    SimpleEventBit = 1;
    
    Dt = (Tops[0].GetTime() - Bottoms[0].GetTime());
    CubicDt = (Tops[0].GetCubicTime() - Bottoms[0].GetCubicTime());
    
    SyncDt = (CorrectedTopTOF - CorrectedBottomTOF);
    CubicSyncDt = (CorrectedCubicTopTOF - CorrectedCubicBottomTOF);
    
    COG = ((Tops[0].GetEnergy() - Bottoms[0].GetEnergy())/
	   (Tops[0].GetEnergy() + Bottoms[0].GetEnergy()));

    if (Tops[0].GetEnergy() > 0 && Bottoms[0].GetEnergy() > 0) {
      AvgEnergy = TMath::Sqrt(Tops[0].GetEnergy()*Bottoms[0].GetEnergy());
    }

    AvgPulseHeight = TMath::Sqrt(Tops[0].GetPulseHeight()*Bottoms[0].GetPulseHeight());
    
    CorrectedCOG=((Tops[0].GetCorrectedEnergy() - Bottoms[0].GetCorrectedEnergy())/
		  (Tops[0].GetCorrectedEnergy() + Bottoms[0].GetCorrectedEnergy()));
    CorrectedDt = (Tops[0].GetCorrectedTime() - Bottoms[0].GetCorrectedTime());
    
    if (Tops[0].GetCorrectedEnergy() > 0 &&
	Bottoms[0].GetCorrectedEnergy() > 0) {
      CorrectedAvgEnergy = TMath::Sqrt(Tops[0].GetCorrectedEnergy()*
				       Bottoms[0].GetCorrectedEnergy());
    }
    
    AvgT = 0.5*(Tops[0].GetTime() + Bottoms[0].GetTime());
    CorrectedAvgT = 0.5*(Tops[0].GetCorrectedTime() +
			 Bottoms[0].GetCorrectedTime());
  }
}

Bool_t lendaBar::operator==(const lendaBar & RHS) {
  Bool_t CalcValuesEqual = true;
  if (this->Dt == RHS.Dt &&
      this->CubicDt == RHS.CubicDt &&
      this->AvgPulseHeight == RHS.AvgPulseHeight &&
      this->CorrectedDt == RHS.CorrectedDt &&
      this->AvgT == RHS.AvgT&&
      this->CorrectedAvgT == RHS.CorrectedAvgT&&
      this->BarId == RHS.BarId&&
      this->TopTOF == RHS.TopTOF&&
      this->BottomTOF == RHS.BottomTOF){
    CalcValuesEqual=true;
  } else {
    printf("lendaBar::operator== failed\n");
    this->Dump();
    RHS.Dump();
    exit(1);
    return false;
  }
  for (Int_t i=0; i<Tops.size(); i++) {
    if (!(this->Tops[i] == RHS.Tops[i])){
      return false;
    }
  }
  for (Int_t i=0; i<Bottoms.size(); i++) {
    if (!(this->Bottoms[i] == RHS.Bottoms[i])){
      return false;
    }
  }
  return true;
}

/************************************************************/
/* lendaEvent class functions *******************************/
/************************************************************/

lendaEvent::lendaEvent() {  
  Clear();
}

void lendaEvent::Clear() {
  for (Int_t i=0; i<bars.size(); i++) {
    bars[i].Clear();
  }
  
  bars.clear(); /* Will call deconstructor on sub objects. */
  N = 0;
  numBars = 0;
  NumOfChannelsInEvent = 0;
  refScint.clear();

  NumUnMappedChannels = 0;
  numObjectScintillators = 0;
  etime = 0;
  internaltime = 0;
  UnMappedChannels.clear();

  ScintIntLow = 0; ScintIntHigh = 0; ScintIntAll = 0;
}

void lendaEvent::Finalize() {
  Int_t tot=0;
  for (Int_t i=0; i<bars.size(); i++){
    bars[i].Finalize();
    tot += bars[i].BarMultiplicity;
  }
  numObjectScintillators = refScint.size();
  
  N = tot + NumUnMappedChannels + numObjectScintillators;
  NumOfChannelsInEvent = N;
}

Bool_t lendaEvent::operator==(const lendaEvent & RHS) {

  bool CalcValuesEqual = true;
  
  if (this->N == RHS.N &&
      this->numBars == RHS.numBars &&
      this->NumOfChannelsInEvent == RHS.NumOfChannelsInEvent &&
      this->NumUnMappedChannels == RHS.NumUnMappedChannels &&
      this->numObjectScintillators == RHS.numObjectScintillators){
    CalcValuesEqual=true;
  } else {
    cout << "lendaEvent::operator== --> event failed." << endl;
    return false;
  }
  
  for (UInt_t ui=0; ui<bars.size(); ui++) {
    if (!(this->bars[ui] == RHS.bars[ui])) {
      return false;
    }
  }
  
  for (UInt_t ui=0; ui<refScint.size(); ui++) {
    if (!(this->refScint[ui] == RHS.refScint[ui])) {
      return false;
    }
  }
  
  for (UInt_t ui=0; ui<UnMappedChannels.size(); ui++) {
    if (!(this->UnMappedChannels[ui] == RHS.UnMappedChannels[ui])) {
      return false;
    }
  }

  return true;
}

lendaBar * lendaEvent::FindBar(Int_t BarId) {
  for (UInt_t ui=0; ui<bars.size(); ui++) {
    if (bars[ui].BarId == BarId){
      return &(bars[ui]);
    }
  }
  return NULL;
}


lendaBar * lendaEvent::FindBar(string Name) {
  for (UInt_t ui=0; ui<bars.size(); ui++){
    if (bars[ui].Name == Name){
      return &(bars[ui]);
    }
  }
  return NULL;
}

lendaBar lendaEvent::FindBarSafe(string Name) {
   for (UInt_t ui=0; ui<bars.size(); ui++) {
     if (bars[ui].Name == Name) {
       return (bars[ui]);
     }
   }
   return lendaBar();
}

lendaBar lendaEvent::FindBarSafe(Int_t BarId) {
  for (UInt_t ui=0; ui<bars.size(); ui++) {
    if (bars[ui].BarId == BarId){
      return (bars[ui]);
    }
  }
  return lendaBar();
}

lendaChannel* lendaEvent::FindChannel(string Name) {
  /* Perform linear search over the tops, bottoms, and reference channels.
     Does not search in the UnMapped channels because those don't 
     have a name */

  for (UInt_t ui=0; ui<bars.size(); ui++) {
    for (UInt_t uj=0; uj<bars[ui].Tops.size(); uj++) {
      if (Name == bars[ui].Tops[uj].GetName()) {
	return &bars[ui].Tops[uj];
      }
    }
    for (UInt_t uj=0; uj<bars[ui].Bottoms.size(); uj++) {
      if (Name == bars[ui].Bottoms[uj].GetName()) {
	return &bars[ui].Bottoms[uj];
      }
    }
  }
  
  for (UInt_t ui=0; ui<refScint.size(); ui++) {
    if (Name == refScint[ui].GetName()) {
      return &refScint[ui];
    }
  }
  return NULL;  
}

lendaChannel* lendaEvent::FindReferenceChannel(string Name) {
  for (UInt_t ui=0; ui<refScint.size(); ui++) {
    if (Name == refScint[ui].GetName()) {
      return &refScint[ui];
    }
  }
  return NULL;
}

/************************************************************/
/* lendaFilter class functions ******************************/
/************************************************************/

lendaFilter::lendaFilter() {
  numOfBadFits=0;
}

lendaFilter::~lendaFilter() { ; }

void lendaFilter::FastFilter(std::vector<UShort_t>& trace, std::vector<Double_t>& thisEventsFF,
			     Double_t FL, Double_t FG) {
  Double_t sumNum1 = 0, sumNum2 = 0;  
  Int_t start = 2*FL + FG - 1;
  
  for (Int_t i=0; i<(Int_t)trace.size(); i++) {
      if (i>=start) {
	for (Int_t j= i-(FL-1); j<=i; j++) {
	  if (j>=0) { sumNum1 += trace[j]; }
	  else { cout << "Oh no..." << endl; }
	}
	for (Int_t j=i-(2*FL+FG-1); j<=i-(FL+FG); j++) {
	  if (j>=0) { sumNum2 += trace[j]; }
	  else { cout << "Oh no..." << endl; }
	}
      }
      thisEventsFF.push_back(sumNum1 - sumNum2);
      
      sumNum1 = 0;  sumNum2 = 0;
  } /* End loop over trace */ 
}

void lendaFilter::FastFilterOp(std::vector<UShort_t>& trace, std::vector<Double_t>& thisEventsFF,
			       Double_t FL, Double_t FG) {
  Double_t sumNum1 = 0;
  Double_t sumNum2 = 0;
  
  Int_t size = trace.size();
  Int_t half = size/2;
  thisEventsFF.resize(size,0);

  Int_t start =half - 0.15*size;
  Int_t end = half + 0.15*size;

  for (Int_t i=start; i<end; i++) {
    
    for (Int_t j= i-(FL-1); j<=i; j++) {
      if (j>=0) { sumNum1 += trace[j]; }
      else { cout << "Oh no..." << endl; }
    }
    for (Int_t j=i-(2*FL+FG-1); j<=i-(FL+FG); j++) {
      if (j>=0) { sumNum2 += trace[j]; }
      else { cout << "Oh no..." << endl; }
    }
    thisEventsFF[i] = sumNum1 - sumNum2;
    sumNum1=0;  sumNum2=0;
  }
}

void lendaFilter::FastFilterFull(std::vector<UShort_t> &trace,
				 std::vector<Double_t> &thisEventsFF,
				 Double_t FL, Double_t FG, Double_t decayTime) {
  printf("lendaFilter::FastFilterFull --> not implemented.\n");
}

std::vector<Double_t> lendaFilter::CFD(std::vector<Double_t>& thisEventsFF,
				       Double_t CFD_delay,
				       Double_t CFD_scale_factor) {
  std::vector<Double_t> thisEventsCFD;
  thisEventsCFD.resize(thisEventsFF.size(), 0);

  for (Int_t j=0; j<(Int_t)thisEventsFF.size() - CFD_delay; j++) {
    thisEventsCFD[j + CFD_delay] = (thisEventsFF[j + CFD_delay] - 
				    thisEventsFF[j]/ ( TMath::Power(2,CFD_scale_factor+1) ));
  }
  return thisEventsCFD;
}

std::vector<Double_t> lendaFilter::CFDOp(std::vector<Double_t>& thisEventsFF,
					 Double_t CFD_delay,
					 Double_t CFD_scale_factor) {
  
  std::vector<Double_t> thisEventsCFD;
  Int_t size = thisEventsFF.size();
  thisEventsCFD.resize(thisEventsFF.size(),0);
  
  Int_t start = size/2 - 0.15*size;
  Int_t end = size/2 + 0.15*size;
  
  for (Int_t j=start; j<end; j++) {
    thisEventsCFD[j + CFD_delay] = (thisEventsFF[j+CFD_delay] - 
				    thisEventsFF[j]/ ( TMath::Power(2,CFD_scale_factor+1) ));
  }
  return thisEventsCFD;
}

Double_t lendaFilter::GetZeroCrossing(std::vector<Double_t>& CFD, Int_t& NumZeroCrossings,
				      Double_t& residual) {

  Double_t softwareCFD;
  std::vector<Double_t> thisEventsZeroCrossings(0);
  Double_t MaxValue = 0;
  Int_t MaxIndex = -1;
  Int_t Window = 40;
  for (Int_t j=(CFD.size()/2.0)-Window; j<(Int_t)(CFD.size()/2.0)+Window; j++) { 
    if (CFD[j] >= 0 && CFD[j+1] < 0 && TMath::Abs(CFD[j] - CFD[j+1]) > 40) {
      softwareCFD =j + CFD[j] / ( CFD[j] + TMath::Abs(CFD[j+1]) );
      thisEventsZeroCrossings.push_back(softwareCFD);
      if (TMath::Abs(CFD[j] - CFD[j+1]) > MaxValue) {
	MaxValue=TMath::Abs(CFD[j] - CFD[j+1]);
	MaxIndex =thisEventsZeroCrossings.size()-1;
	residual=CFD[j];
      }
    }
  }
  NumZeroCrossings=thisEventsZeroCrossings.size();
  if (thisEventsZeroCrossings.size() == 0) {
    return BADNUM;
  } else {
    return thisEventsZeroCrossings[MaxIndex];
  }

  return BADNUM;
}

Double_t lendaFilter::GetZeroCrossingImproved(std::vector<Double_t>& CFD, Int_t& MaxInTrace,
					      Double_t& residual) {
 
  Double_t softwareCFD = BADNUM;
  Double_t thisEventsZeroCrossing;
  Int_t Window=40;

  for (Int_t j=MaxInTrace; j<(Int_t)(CFD.size()/2)+Window; j++) { 
    if (CFD[j] >= 0 && CFD[j+1] < 0) {
      softwareCFD = j + CFD[j] / ( CFD[j] + TMath::Abs(CFD[j+1]) );
      residual = CFD[j];
      break;
    }
  }
  return softwareCFD;
}

vector<Double_t> lendaFilter::GetAllZeroCrossings(std::vector<Double_t> & CFD) {
  Double_t softwareCFD;
  std::vector<Double_t> thisEventsZeroCrossings(0);

  for (Int_t j=0; j<(Int_t)CFD.size()-1; j++) {
    if (CFD[j] >= 0 && CFD[j+1] < 0 ) {
      softwareCFD = j + CFD[j] / ( CFD[j] + TMath::Abs(CFD[j+1]) );
      thisEventsZeroCrossings.push_back(softwareCFD);
    }
  }
  return thisEventsZeroCrossings;
}

Double_t lendaFilter::GetZeroCrossingOp(std::vector<Double_t> & CFD,Int_t & NumZeroCrossings) {

  Double_t softwareCFD = -10;
  std::vector<Double_t> thisEventsZeroCrossings(10,0);
  Double_t MaxValue = 0;
  Int_t MaxIndex = -1;

  for (Int_t j=(CFD.size()/2.0)-20; j<(Int_t)(CFD.size()/2.0)+20; j++) { 
    if (CFD[j] >= 0 && CFD[j+1] < 0 && TMath::Abs(CFD[j] - CFD[j+1]) > 5) {
      softwareCFD =j + CFD[j] / ( CFD[j] + TMath::Abs(CFD[j+1]) );
      break;
    }
  }
  NumZeroCrossings = -1;
  if (softwareCFD == -10) {
    return BADNUM;
  } else  {
    return softwareCFD;
  }
}


vector<Double_t> lendaFilter::GetZeroCrossingHighRate(const std::vector<Double_t>& CFD,
						      const std::vector<Int_t>& PeakSpots) {
  std::vector<Double_t> thisEventsZeroCrossings(0);  
  for (UInt_t ui = 0; ui<PeakSpots.size(); ui++) {
    for (Int_t j=PeakSpots[ui]; j<(Int_t)CFD.size()-1; j++) { 
      if (CFD.at(j)>= 0 && CFD.at(j+1) < 0) {
	Double_t softwareCFD = j + CFD[j] / ( CFD[j] + TMath::Abs(CFD[j+1]));
	thisEventsZeroCrossings.push_back(softwareCFD);
	j = CFD.size()+1000;
	break;
      }
    }
  }
  return thisEventsZeroCrossings;
}

vector<Double_t> lendaFilter::GetZeroCrossingCubicHighRate(const std::vector<Double_t>& CFD,
							   const std::vector<Int_t>& PeakSpots) {
  std::vector<Double_t> thisEventsZeroCrossings(0);  
  for (UInt_t ui = 0; ui<PeakSpots.size(); ui++) {
    for (Int_t j=PeakSpots[ui]; j<(Int_t)CFD.size()-1; j++) { 
      if (CFD.at(j)>= 0 && CFD.at(j+1) < 0) {
	Double_t cubicCFD = DoMatrixInversionAlgorithm(CFD,j);
	thisEventsZeroCrossings.push_back(cubicCFD);
	j = CFD.size()+1000;
	break;
      }
    }
  }
  return thisEventsZeroCrossings;
}

vector<Int_t> lendaFilter::GetPulseHeightHighRate(const std::vector<UShort_t>& trace,
						  const std::vector<Int_t>& PeakSpots) {
  vector<Int_t> result;
  for (UInt_t ui = 0; ui<PeakSpots.size(); ui++) {
    result.push_back(trace[PeakSpots[ui]]);
  }
  return result;
}

Double_t lendaFilter::DoMatrixInversionAlgorithm(const std::vector<Double_t>& CFD,
						 Int_t theSpotAbove) {
  
  Double_t x[4];
  TMatrixT<Double_t> Y(4,1);
  
  for (Int_t i=0;i<4;i++) {
    x[i] = theSpotAbove - 1 + i; 
    Y[i][0]=CFD[theSpotAbove - 1 + i];
  }
  TMatrixT<Double_t> A(4,4);

  for (Int_t row=0; row<4; row++) {
    for (Int_t col=0; col<4; col++) {
      A[row][col]= pow(x[row], 3 - col);
    }
  }

  TMatrixT<Double_t> invertA = A.Invert();
  TMatrixT<Double_t> Coeffs(4,1);
  Coeffs = invertA*Y;
  
  bool notDone =true;
  Double_t left = x[1];
  Double_t right =x[2];
  Double_t valUp = getFunc(Coeffs, left);
  Double_t valDown = getFunc(Coeffs, right);

  Int_t loopCount=0;
  while (notDone) {
    loopCount++;
    if (TMath::Abs(TMath::Abs(valUp)-TMath::Abs(valDown) ) <0.001) {
      notDone = false;
    }
    
    Double_t mid = (left+right)/2.0;
    Double_t midVal = getFunc(Coeffs,mid);
    
    if (midVal > 0) { left=mid;}
    else { right=mid; }
    
    valUp = getFunc(Coeffs, left);
    valDown = getFunc(Coeffs, right);
    if (loopCount >30) {
      notDone = false;
      left = BADNUM;
    }
  }
  return left;
}

vector<Double_t> lendaFilter::GetMatrixInversionAlgorithmCoeffients(const std::vector<Double_t>& CFD,
								    Int_t& ReturnSpotAbove) {
  std::map <Double_t,Int_t> zeroCrossings;
  Double_t max=0;

  Int_t begin = (CFD.size()/2)-40;
  Int_t end = (CFD.size()/2)+40;
  
  for (Int_t i=begin; i<end; i++) {
    if (CFD[i]>=0 && CFD[i+1]<0) {
      Double_t val = CFD[i] - CFD[i+1];
      if ( val > max) { max = val; }
      zeroCrossings[val] = i;
    }
  }
  
  Int_t theSpotAbove = zeroCrossings[max];
  ReturnSpotAbove = theSpotAbove;

  Double_t x[4];
  TMatrixT<Double_t> Y(4,1);
  for (Int_t i=0; i<4; i++) {
    x[i]= theSpotAbove - 1 + i; 
    Y[i][0]=CFD[ theSpotAbove - 1 + i];
  }

  TMatrixT<Double_t> A(4,4);
  for (Int_t row=0; row<4; row++) {
    for (Int_t col=0; col<4; col++) {
      A[row][col]= pow(x[row], 3 - col);
    }
  }

  TMatrixT<Double_t> invertA = A.Invert();
  TMatrixT<Double_t> Coeffs(4,1);
  Coeffs = invertA*Y;
  
  vector<Double_t> retVec;
  for (Int_t i=0; i<4; i++) {
    retVec.push_back(Coeffs[i][0]);
  }
  return retVec;
}

Double_t lendaFilter::GetZeroCubic(std::vector<Double_t>& CFD, Int_t &MaxInTrace) {
  Int_t begin = MaxInTrace;
  Int_t end = (CFD.size()/2)+40;

  for (Int_t i=begin; i<end; i++) {
    if (CFD[i]>=0 && CFD[i+1]<0) {
      return DoMatrixInversionAlgorithm(CFD,i);
      break;
    }
  }

  return -1.;
}

Double_t lendaFilter::GetZeroFitCubic(std::vector<Double_t>& CFD) {  

  std::map <Double_t,Int_t> zeroCrossings;
  Double_t max = 0;
  
  /* Restrict search to center part of trace */
  Int_t begin = (CFD.size()/2) - 40;
  Int_t end = (CFD.size()/2) + 40;

  for (Int_t i=begin; i<end; i++) {
    if (CFD[i]>=0 && CFD[i+1]<0) {
      Double_t val = CFD[i] - CFD[i+1];
      if (val > max) { max = val; }
      zeroCrossings[val] = i;
    }
  }
  Int_t theSpotAbove = zeroCrossings[max];

  /* Points to consider are: 
     {spotAbove-2, spotAbove-1, spotAbove, spotBelow, spotBelow+1, spotBelow+2}
     
     Use fitting to find 3rd order polynomial from these 6 points. */
  
  Double_t x[6];
  Double_t y[6];
  for (Int_t i=0; i<6; i++) {
    x[i]= theSpotAbove - 2 + i; //first point is the one before zerocrossing
    y[i]=CFD[theSpotAbove - 2 + i];
  }
  TGraph theGraphForFitting(6,x,y);
  
  TFitResultPtr fitPointer = theGraphForFitting.Fit("pol3", "QSN");
  Int_t status = fitPointer;
 
  if ( status == 0) { 
    vector<Double_t> Coeffs(4);
    Coeffs[0]=fitPointer->Value(0);
    Coeffs[1]=fitPointer->Value(1);
    Coeffs[2]=fitPointer->Value(2);
    Coeffs[3]=fitPointer->Value(3);
  
    bool notDone =true;
    Double_t left = x[2];
    Double_t right =x[3];
    Double_t valUp = getFunc(Coeffs, left);
    Double_t valDown = getFunc(Coeffs, right);
    Int_t loopCount = 0;
    while (notDone) {
      loopCount++;
      if (TMath::Abs(TMath::Abs(valUp)-TMath::Abs(valDown) ) <0.001) {
	notDone = false;
      }
 
      Double_t mid = (left+right)/2.0;
      Double_t midVal = getFunc(Coeffs,mid);
   
      if (midVal > 0) { left = mid; }
      else { right = mid; }
      valUp = getFunc(Coeffs, left);
      valDown = getFunc(Coeffs, right);
      if (loopCount > 30) {
	notDone = false;
	left = BADNUM;
      }
    }
    return left;
  } else {
    cout << "FIT FAILED" << endl;
    return BADNUM;
  }
}

Double_t lendaFilter::getFunc(TMatrixT<Double_t>Coeffs, Double_t x) {
  Double_t total = 0;
  for (Int_t i=0; i<4; i++) {
    total += Coeffs[i][0]*TMath::Power(x,3-i);
  }
  return total;
}

Double_t lendaFilter::getFunc(vector<Double_t>& Coeffs,Double_t x) {
  Double_t total =0;
  for (Int_t i=0; i<4; i++) {
    total += Coeffs[i]*TMath::Power(x,i);
  }
  return total;
}

Double_t lendaFilter::fitTrace(std::vector<UShort_t> & trace,Double_t sigma,Double_t num) {

  Int_t size = (Int_t) trace.size();
  std::vector<Double_t> y_values, x_values;
  for ( Int_t i=0; i<size; ++i) {
    y_values.push_back((Double_t) trace[i]);
    x_values.push_back((Double_t) i);
  }

  Double_t max = -1;
  Int_t maxBin = -1;
  for (Int_t i=0; i<size; ++i) {
    if (y_values[i] > max) {
      max = y_values[i];
      maxBin = i;
    }
  }

  /* Setup fit function */
  Double_t base = 300;
  Double_t A = 0;
  Int_t fitWindowWidth = 10;
  Double_t mu = -1000;
  
  /* 1.7711 the sigma determined from looking at traces... */
  std::stringstream stream;
  stream << "[2]+[0]*exp(-0.5*( ((x-[1])/" << sigma << ")^2) )";
  
  TF1 *myfit = new TF1("myfit", stream.str().c_str(), 0, 200); 

  myfit->SetParameter(0, 1);
  myfit->SetParameter(1, 100);
  myfit->SetParameter(2, 300);
  myfit->SetParLimits(0, 0, 1000);

  TGraph * gr = new TGraph(size, x_values.data(), y_values.data());

  TFitResultPtr fitPointer = gr->Fit("myfit", "S0Q", "",
				     maxBin-fitWindowWidth, maxBin+fitWindowWidth);
  Int_t fitStatus = fitPointer;
  
  if ( fitStatus == 0 ) { 
    mu = fitPointer->Value(1);
    A = fitPointer->Value(0);
    base = fitPointer->Value(2);
  } else {
    numOfBadFits++;
  }
  
  gr->Delete();
  myfit->Delete();
  return mu;
}

Double_t lendaFilter::GetEnergyOld(std::vector<UShort_t>& trace) {
  
  /////////////////////////////////////////////////////////////////////
  // Basic Algorithm for extracting pulse integral from trace	     //
  // uses first 10 points as background level. Sums up the	     //
  // whole trace then subtracts background to background subtracted  //
  // energy 							     //
  /////////////////////////////////////////////////////////////////////

  Double_t thisEventsIntegral;
  Double_t sum = 0;
  Double_t signalTotalIntegral = 0;
  for ( Int_t i=0; i<10; i++) { sum += trace[i]; }
  sum = sum/10.0; 
  for (Int_t i=0; i<(Int_t)trace.size(); ++i) {
    signalTotalIntegral += trace[i];
  }
  if ((signalTotalIntegral - sum*trace.size()) > 0) {
    thisEventsIntegral = signalTotalIntegral - sum*trace.size();
  } else {
    thisEventsIntegral = BADNUM;
  }
  return thisEventsIntegral;
}

Double_t lendaFilter::GetEnergy(std::vector<UShort_t>& trace, Int_t MaxSpot) {
  
  ////////////////////////////////////////////////////////////////////////////
  // Slightly more complicated algorithm for getting energy of pulse.	    //
  // Takes the first 20% of the trace as a background.  Then it takes only  //
  // a window around the maximum spot as the pulse part of the trace.  This //
  // helps with strange pulse where there might be additional structures in //
  // the trace that would distort things (afterpulse, crazy noise..)	    //
  ////////////////////////////////////////////////////////////////////////////

  Double_t thisEventsIntegral;
  Double_t sumBegin = 0;
  Double_t sumEnd = 0;
  Double_t signalIntegral = 0;

  Int_t traceLength = trace.size();
  Int_t LengthForBackGround = 0.2*traceLength;

  Int_t windowForEnergy = 20;

  if (MaxSpot - windowForEnergy < 0) {
    return BADNUM;
  } else if (MaxSpot + windowForEnergy > (traceLength -1) ) {
    return BADNUM;
  }

  for (Int_t i=0; i<LengthForBackGround; i++) {
    sumBegin = sumBegin + trace[i];
    sumEnd = sumEnd + trace[traceLength - 1 - i];
  }
  
  Double_t BackGround = BADNUM;

  if (MaxSpot > LengthForBackGround && MaxSpot < (traceLength-LengthForBackGround) ) {
    /* The peak is in the middle use the beginning of trace as background */
    BackGround = sumBegin/LengthForBackGround;    
  } else if (MaxSpot < LengthForBackGround) {
    /* The peak is in the beginning of trace use end as background */
    BackGround = sumEnd/LengthForBackGround;
  } else if (MaxSpot > (traceLength-LengthForBackGround)) {
    /* Peak is at end use beginning for background */
    BackGround = sumBegin/LengthForBackGround;
  } else{
    return BADNUM;
  }
  
  for (Int_t i=MaxSpot-windowForEnergy;i< MaxSpot+windowForEnergy;++i) {
    signalIntegral += trace[i];
  }

  if (signalIntegral - BackGround *(2*windowForEnergy)>0) {
    thisEventsIntegral = signalIntegral - BackGround *(2*windowForEnergy);
  }  else{
    thisEventsIntegral = BADNUM;
  }
  return thisEventsIntegral;
}

vector<Double_t> lendaFilter::GetEnergyHighRate(const std::vector<UShort_t>& trace,
						std::vector<Int_t>& PeakSpots,
						vector<Double_t>& theUnderShoots,
						Double_t& MaxValueOut, Int_t& MaxIndexOut) {
  Int_t NumberOfSamples = 8;
  Int_t TraceLength = trace.size();
 
  /* Each chunk of sample size trace points is used to calculate the point in the trace with
     the most stable baseline */
  Int_t SampleSize = TMath::Floor(TraceLength/Float_t(NumberOfSamples));
  Int_t Remainder = TraceLength-SampleSize*NumberOfSamples;
  
  vector<Double_t> StandardDeviations(NumberOfSamples,0);
  vector<Double_t> Averages(NumberOfSamples,0);

  /* First for each sample window find the average of those points and find the 
     maximum point (will exclude the reamining points at the end of the trace */
  Double_t Max = 0;
  Int_t MaxIndex = -1;
  for (Int_t i=0; i<NumberOfSamples; i++) {
    for (Int_t j=0; j<SampleSize; j++) {
      Double_t val = trace[j + i*SampleSize];
      Averages[i] += val;
      if (val > Max) {
	Max = val;
	MaxIndex = j + i*SampleSize;
      }
    }
    Averages[i] /= SampleSize;
  }

  /* Now look at the end of the trace to see if the maximum is there */
  for (Int_t i=SampleSize*NumberOfSamples; i<TraceLength; i++) {
    Double_t val = trace[i];
    if (val > Max) {
      Max = val;
      MaxIndex = i;
    }
  }
  
  /* Copy values to output variables */
  MaxValueOut = Max;
  MaxIndexOut = MaxIndex;
  
  //Find standard deviation for each base line chunk */
  map<Double_t,pair<Int_t,Int_t> > baseLineRanges;
  map<Double_t,Double_t> StandardDeviation2Average;
  for (Int_t i=0; i<NumberOfSamples; i++) {
    for (Int_t j=0; j<SampleSize; j++) {
      StandardDeviations[i] += TMath::Power((trace[j +i*SampleSize]-Averages[i]),2);
    }
    baseLineRanges.insert( make_pair(StandardDeviations[i],make_pair(i*SampleSize,(i+1)*SampleSize)));
    StandardDeviation2Average.insert(make_pair(StandardDeviations[i],Averages[i]));
  }

  /* Take beginning of map which will give the smallest stardard deviation. */
  map<Double_t, Double_t>::iterator it = StandardDeviation2Average.begin();
  /* The best baseline will be the average of the chunk with the smallest standard deviation */
  Double_t BaseLine = it->second;
  Double_t stanDev =it->first;
  
  /* Threshold should be the amount above baseline. */
  Int_t ThresholdForOtherPulseInTrace= MaxValueOut*0.15;


  Int_t windowForEnergy=2;
  Int_t EndTraceCut =2;//2*windowForEnergy;
  Double_t currentMax=0;
  bool HasCrossedThreshold=false;

  /* Loop through the trace and try and find the peak values.  Exclude the first point.  
     It cannot be a maximum. */
  for (Int_t i=1; i<TraceLength-EndTraceCut; i++) {
    Double_t val = trace[i] - BaseLine;
    if (val > ThresholdForOtherPulseInTrace) {
      HasCrossedThreshold = true;
    }
    if (val > currentMax) {
      currentMax = val;
    } else if (val < currentMax && HasCrossedThreshold) {
      /* Possible peak spot. */
      Int_t PossibleSpot = i-1;
      if (trace[PossibleSpot-1] < trace[PossibleSpot]&&
	  trace[PossibleSpot+1] < trace[PossibleSpot]) {
	PeakSpots.push_back(PossibleSpot);

	currentMax = 0;
	HasCrossedThreshold = false;

	/* Now must skip foward in time around the length of one RF bucket, the minimum 
	   amount of time before a second real pulse could be there. */
	Int_t index = i + 8;
	if (index > TraceLength - windowForEnergy) {
	  i = TraceLength+1000;
	  break;
	} else {
	  i = index;
	}
      }
    }
  }

  /* The Maximums have been found, now finally we calculate pulse integrals. */
  vector<Double_t> theEnergies;
  for (UInt_t ui=0; ui<PeakSpots.size(); ui++) {
    Double_t temp=0;
    for (Int_t j=0; j<2; j++) {
      temp += (trace[PeakSpots[ui]+j] - BaseLine);
    }
    theEnergies.push_back(temp);
    
    Double_t min = 100000000.0;
    for (Int_t j=0; j<6; j++) {
      if (trace[PeakSpots[ui]-j] < min) {
	min = trace[PeakSpots[ui]-j];
      }
    }
    Double_t cut = 0.94*BaseLine;
    if (min < cut) { 
      theUnderShoots.push_back(BaseLine-min);
    } else {
      theUnderShoots.push_back(BADNUM);
    }
  }
  return theEnergies;
}

Double_t lendaFilter::GetGate(std::vector<UShort_t>& trace, Int_t start, Int_t L) {
  Int_t range = L;
  Int_t window = floor(trace.size()/5.0);
  Double_t bgk = 0;

  for (Int_t i=0; i<window; i++) {
    bgk = trace[i] + bgk;
  }
  bgk = bgk/(window);

  Double_t total=0;
  for (Int_t i=start; i<start+L; i++) {
    total += trace[i];
  }
  return (total - (bgk*range));
}

Int_t lendaFilter::GetMaxPulseHeight(vector<UShort_t>& trace, Int_t& MaxSpot) {
  Int_t imaxSpot = -1;
  Double_t max = 0;
  for (Int_t i=0; i<trace.size(); i++) {
    if (trace[i] > max) {
      max = trace[i];
      imaxSpot = i;
    }
  }
  MaxSpot = imaxSpot;
  return max;
}

Int_t lendaFilter::GetMaxPulseHeight(vector<Double_t>& trace, Int_t& MaxSpot) {
  Int_t imaxSpot = -1;
  Double_t max = 0;
  for (Int_t i=0; i<trace.size(); i++) {
    if (trace[i] > max) {
      max = trace[i];
      imaxSpot = i;
    }
  }
  MaxSpot = imaxSpot;
  return max;
}

Int_t lendaFilter::GetStartForPulseShape(Int_t MaxSpot) {
  if (MaxSpot-4 > 0) {
    return (MaxSpot-4);
  } else {
    return 0;
  }
}

Int_t lendaFilter::CalculateCFD(vector<UShort_t> trace,vector<Double_t> &CFDOut) {
  CFDOut.clear();
  Int_t TFLength = 6, TFGap = 0, CFDSF = 4, CFDDelay = 6;
  Double_t CFDWeight = 0.5;
  Int_t CFDThreshold = 50;

  Int_t CFDTime = 0, CFDTriggerPoint;

  Double_t cfdvalue = -1;
  Int_t cfdtrigpoint = -1;
  Int_t A = 0, m = 0;
 
  Double_t S0, S1, S2, S3;
  
  if(TFLength < 0) {
    cout << "Set the Trigger Filter Length!" << endl;
    return -107;
  }
  if(TFGap < 0) {
    cout << "Set the Trigger Filter Gap!" << endl;
    return -108;
  }  
  if(CFDSF < 0) {
    cout << "Set CFD Scale Factor!" << endl;
    return -115;
  }
  if(CFDDelay < 0) {
    cout << "Set CFD Delay!" << endl;
    return -116;
  } 
  
  if(trace.size() > 0) {    
    for(m=0; m<(CFDDelay+TFLength+TFGap+1); m++) {
      CFDOut.push_back(0);
    }
    m = CFDDelay + 2*TFLength + TFGap;
    
    for(A=(m-CFDDelay-2*TFLength-TFGap); A<(m-CFDDelay-TFGap-TFLength); A++) {
      S3 += trace[A];
    }
    for(A=(m-CFDDelay-TFLength); A<(m-CFDDelay); A++) {
     S2 += trace[A];
    }  
    for(A=(m-2*TFLength-TFGap); A<(m-TFGap-TFLength); A++) {
      S1 += trace[A];
    }
    for(A=m-TFLength; A<m; A++) {
      S0 += trace[A];
    }
    
    for(m=(CFDDelay+2*TFLength+TFGap+1); m<(Int_t)(trace.size()); m++) {
      S3 = S3-trace[m-CFDDelay-2*TFLength-TFGap-1]+trace[m-CFDDelay-TFGap-TFLength-1];
      S2 = S2-trace[m-CFDDelay-TFLength-1]+trace[m-CFDDelay-1];
      S1 = S1-trace[m-2*TFLength-TFGap-1]+trace[m-TFGap-TFLength-1];
      S0 = S0-trace[m-TFLength-1]+trace[m-1];
      
      Double_t cfdvalue2 = CFDWeight*(S0-S1) - (S2-S3);
      CFDOut.push_back(cfdvalue2);     
    }
    
  } else {
    cout << "Trace is length 0!" << endl;
    return -105;
  }
  return 0;  
}

vector<Double_t> lendaFilter::GetNewFirmwareCFD(const vector<UShort_t>& trace, Int_t FL,
						Int_t FG, Int_t d, Int_t w) {
  vector<Double_t> CFDOut;
  if(trace.size() == 0) {
    CFDOut.resize(1,BADNUM);
    return CFDOut;
  } else { ; }

  Int_t TFLength = FL;
  Int_t TFGap = FG;
  Int_t CFDSF = 0;
  Int_t CFDDelay = d;

  Double_t CFDWeight = GetNewFirmwareCFDWeight(w);
  Int_t CFDThreshold = 50, CFDTime = 0, CFDTriggerPoint;

  Int_t cfdtrigpoint = -1, A = 0, m = 0;

  ///////////////////////////////////////////
  // There are 4 sums in the CFD algorithm //
  ///////////////////////////////////////////
  Double_t S0 = 0., S1 = 0., S2 = 0., S3 = 0.;
  
  for(m=0; m<(CFDDelay+2*TFLength+TFGap+1); m++) {
    CFDOut.push_back(0);
  }
    
  m = CFDDelay + 2*TFLength + TFGap;

  for(A=(m-CFDDelay-2*TFLength-TFGap); A<(m-CFDDelay-TFGap-TFLength); A++) {
    S3 += trace[A];
  }
  for(A=(m-CFDDelay-TFLength); A<(m-CFDDelay); A++) {
    S2 += trace[A];
  }  
  for(A=(m-2*TFLength-TFGap); A<(m-TFGap-TFLength); A++) {
    S1 += trace[A];
  }
  for(A=m-TFLength; A<m; A++) {
    S0 += trace[A];
  }

  for(m=(CFDDelay+2*TFLength+TFGap+1); m<(Int_t)(trace.size()); m++) {
    S3 = S3-trace[m-CFDDelay-2*TFLength-TFGap-1]+trace[m-CFDDelay-TFGap-TFLength-1];
    S2 = S2-trace[m-CFDDelay-TFLength-1]+trace[m-CFDDelay-1];
    S1 = S1-trace[m-2*TFLength-TFGap-1]+trace[m-TFGap-TFLength-1];
    S0 = S0-trace[m-TFLength-1]+trace[m-1];
    Double_t cfdvalue = CFDWeight*(S0-S1) - (S2-S3);

    CFDOut.push_back(cfdvalue);     
  }
  return CFDOut;
}

Double_t lendaFilter::GetNewFirmwareCFDWeight(Int_t CFDScaleFactor) {
  Double_t temp = (Double_t) CFDScaleFactor;
  return (8.0-temp)/8.0;  
}
