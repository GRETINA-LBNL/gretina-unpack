#include "LENDA-DDAS.h"
#include "LENDA-Controls.h"

/*************************************************************/

lendaSettings::lendaSettings(){
  this->fName = "TheSettings";
  this->fTitle = "LENDA Settings";
}

lendaSettings::~lendaSettings() { ; }

void lendaSettings::PrintAll() {
  std::map<string, Double_t>::iterator it;
  for (it = TheTimingOffsets.begin(); it != TheTimingOffsets.end(); it++) {
    printf("/n");
    PrintChannelCorrections(it->first);
    PrintChannelMapInfo(Name2GlobalID[it->first]);
    PrintFilterInfo(it->first);
  }
}

void lendaSettings::PrintChannelCorrections(string Name) {
  printf("Channel Name %6s has slope %10.4lf",
	 Name.c_str(), TheEnergySlopes[Name]);
  printf(" has intercept %10.4lf and timing offset %10.4lf\n",
	 TheEnergyIntercepts[Name], TheTimingOffsets[Name]);
}

void lendaSettings::PrintChannelMapInfo(int GlobalID) {
  printf("Channel Name %6s has DDAS ID %4d",
	 GlobalID2FullName[GlobalID].c_str(), GlobalID);
  printf(" it's reference channel is %6s whith DDAS ID %4d\n",
	 GlobalID2RefName[GlobalID].c_str(),GlobalID2RefGlobalID[GlobalID]);
}

void lendaSettings::PrintChannelMapInfo(string Name) {
  PrintChannelMapInfo(Name2GlobalID[Name]);
}

void lendaSettings::PrintFilterInfo(string Name) {
  printf("Channel Name %6s has FL:%5d FG:%5d", Name.c_str(),
	 TheFLs[Name], TheFGs[Name]);
  printf(" d:%5d w:%5d NoTraceAnalysis:%5d \n",
	 Theds[Name],Thews[Name],TheDontTraceAnalyzeFlags[Name]);
}

void lendaSettings::AddCorrectionSettings(string Name, double slope,
					  double inter, double toff){
  TheTimingOffsets[Name] = toff;
  TheEnergySlopes[Name]=slope;
  TheEnergyIntercepts[Name]=inter;
}

void lendaSettings::AddMapSettings(string Name, int GlobalID,
				   string RefName, int refGlobalID){
  GlobalID2FullName[GlobalID] = Name;
  GlobalID2RefGlobalID[GlobalID] = refGlobalID;
  GlobalID2RefName[GlobalID] = RefName;
  Name2GlobalID[Name] = GlobalID;
}

void lendaSettings::AddFilterSettings(string Name, int FL, int FG,
				      int d, int w, bool flag) {
  TheFLs[Name] = FL;
  TheFGs[Name] = FG;
  Theds[Name] = d;
  Thews[Name] = w;
  TheDontTraceAnalyzeFlags[Name] = flag;
}

void lendaSettings::SetBarIds(map<string,int> v) {
  BarIds = v;
  BuildReverseMap();
}

void lendaSettings::BuildReverseMap() {
  std::map<string, Int_t>::iterator it;
  for (it = BarIds.begin(); it != BarIds.end(); it++) {  
    BarId2Name[it->second] = it->first;
  }
}

/*************************************************************/

lendaPacker::lendaPacker(lendaSettings* v) {
  /* Will use bad defaults to ensure that the calling program set 
     all filter parameters. */
   lg = 15;
   sg = 8;
   lg2 = -1;
   sg2 = -1;
   traceDelay = 0;
   jentry = -1;
   theSettings = v;
   
   referenceChannelPattern = "OBJ";

   saveTraces=true;
   
   Reset();
}

lendaPacker::~lendaPacker() {
  Reset();
}

void lendaPacker::SetSettingFileNames(string MapFileName,
				      string CorrectionsFileName){
  mapFileName = MapFileName;
  correctionsFileName = CorrectionsFileName;

  Info("lendaPacker","Using MapFile %s and CorrectionsFile %s",
       mapFileName.c_str(), correctionsFileName.c_str());
  
  BuildMaps();
}

void lendaPacker::Reset(){

  /* Clear the temporary filter vectors. */
  thisEventsFF.clear();
  thisEventsCFD.clear();
  
  thisChannelsEnergies.clear();
  thisChannelsPulseHeights.clear();
  
  thisChannelsSoftwareCFDs.clear();
  thisChannelsCubicCFDs.clear();
  
  thisChannelsPeakSpots.clear();
  
  thisChannelsUnderShoots.clear();
  
  /* Reset packer variables. */
  thisEventsFilterHeight = 0;
  longGate = 0;
  shortGate = 0;
  cubicFitCFD = 0;
  start = 0;
  thisChannelsNumZeroCrossings = 0;
  thisChannelsCFDResidual = 0;
  jentry = -1;
}

void lendaPacker::CalcTimeFilters(vector<UShort_t>& theTrace, mapInfo info) {
  printf(DMAGENTA "<lendaPacker::CalcTimeFilters>: DO NOT CALL!!!\n" RESET_COLOR);
  throw -1;
}

void lendaPacker::CalcEnergyGates(vector<UShort_t>& theTrace, mapInfo info) {
  printf(DMAGENTA "<lendaPacker::CalcEnergyGates>: DO NOT CALL!!!\n" RESET_COLOR);
  throw -1;
}

void lendaPacker::CalcAll(vector<UShort_t>& theTrace, mapInfo info) {

  //  CalcTimeFilters(theChannel->trace,info);
  //  CalcEnergyGates(theChannel->trace,info);
  //  vector <UShort_t> theTrace = theChannel->trace;
  Int_t FL = info.FL;
  Int_t FG = info.FG;
  Int_t d  = info.d;
  Int_t w  = info.w;
  //////////////////////////////////////////////////////////////////////////
  // If there are trace present in the data preform timing related trace  //
  // anaylsis routines			        			  //
  //////////////////////////////////////////////////////////////////////////
  if (theTrace.size() != 0 && info.DontTraceAnalyze == false) {
    theFilter.FastFilter(theTrace, thisEventsFF, FL, FG);
    // thisEventsCFD = theFilter.CFD(thisEventsFF,d,w);
    thisEventsCFD = theFilter.GetNewFirmwareCFD(theTrace, FL, FG, d, w);
    
    /* If this channel is NOT a reference channel run the basic algorithms. */ 
    if (! info.IsAReferenceChannel ) {
      Int_t imaxInTrace = -1;
      Int_t imaxInFilter = -1;
      /* Since this not a reference channel push the maximum point on 
	 to the vector holding the pulseheights.  There will be just the 
	 one entry in the vector. */
      thisChannelsPulseHeights.push_back( theFilter.GetMaxPulseHeight(theTrace, imaxInTrace) );
      
      /* For the filter things are not implemented as a vector. */
      thisEventsFilterHeight = theFilter.GetMaxPulseHeight(thisEventsFF, imaxInFilter);
      thisChannelsEnergies.push_back( theFilter.GetEnergy(theTrace, imaxInTrace) );
      shortGate = theFilter.GetGate(theTrace, imaxInTrace-4, 10);

      /* Run the basic timing algorithms for this channel. */
      thisChannelsSoftwareCFDs.push_back( theFilter.GetZeroCrossingImproved(thisEventsCFD, imaxInTrace,
									    thisChannelsCFDResidual) ); 
      thisChannelsCubicCFDs.push_back( theFilter.GetZeroCubic(thisEventsCFD, imaxInTrace) );
      // cubicFitCFD=theFilter.GetZeroFitCubic(thisEventsCFD);
    
    } else { /* This channel is a reference channel.  Perform the 
		highrate related algorithms. */

      Double_t trash;
      Int_t t;
      thisChannelsEnergies = theFilter.GetEnergyHighRate(theTrace, thisChannelsPeakSpots, thisChannelsUnderShoots, trash, t);
      if (thisChannelsPeakSpots.size() != 0) {
	thisChannelsPulseHeights = theFilter.GetPulseHeightHighRate(theTrace, thisChannelsPeakSpots);
	thisChannelsSoftwareCFDs = theFilter.GetZeroCrossingHighRate(thisEventsCFD, thisChannelsPeakSpots);
	thisChannelsCubicCFDs = theFilter.GetZeroCrossingCubicHighRate(thisEventsCFD, thisChannelsPeakSpots);
	shortGate = theFilter.GetGate(theTrace,thisChannelsPeakSpots[0]-4, 10);
	longGate = theFilter.GetGate(theTrace,thisChannelsPeakSpots[0]-4, 30);
      }
    }
  } /* End trace is not 0 and do trace analysis */
  if (thisChannelsPulseHeights.size() == 0) {
    thisChannelsPulseHeights.push_back(BADNUM);
  }
  if(thisChannelsSoftwareCFDs.size() == 0) {
    thisChannelsSoftwareCFDs.push_back(BADNUM);
  }
  if (thisChannelsCubicCFDs.size() == 0) {
    thisChannelsCubicCFDs.push_back(BADNUM);
  }
  if(thisChannelsEnergies.size() == 0) {
    thisChannelsEnergies.push_back(BADNUM);
  }
  if(thisChannelsUnderShoots.size() == 0) {
    thisChannelsUnderShoots.push_back(BADNUM);
  }
}

void lendaPacker::CalcAll(ddasChannel* c, mapInfo info) {
  printf(DMAGENTA "<lendaPacker::CalcAll(ddasChannel*)> not yet implemented\n" RESET_COLOR);
}

void lendaPacker::CalcAll(lendaChannel* c, mapInfo info) {
  printf(DMAGENTA "<lendaPacker::CalcAll(lendaChannel*)> not yet implemented\n" RESET_COLOR);
}

void lendaPacker::ForceAllBarFilters(Int_t FL, Int_t FG, Int_t d, Int_t w){
  std::map<Int_t, mapInfo>::iterator it;
  for (it = GlobalIDToMapInfo.begin(); it != GlobalIDToMapInfo.end(); it++) {
    if (it->second.IsAReferenceChannel == false){
      it->second.FL = FL;
      it->second.FG = FG;
      it->second.d = d;
      it->second.w = w;
    }
  }
  UpdateSettings();
}

void lendaPacker::ForceAllReferenceFilters(Int_t FL, Int_t FG, Int_t d, Int_t w) {
  std::map<Int_t, mapInfo>::iterator it;
  for (it = GlobalIDToMapInfo.begin(); it != GlobalIDToMapInfo.end(); it++) { 
    if (it->second.IsAReferenceChannel == true) {
      it->second.FL = FL;
      it->second.FG = FG;
      it->second.d = d;
      it->second.w = w;
    }
  }
  UpdateSettings();
}

void lendaPacker::ForceAllFilters(Int_t FL, Int_t FG, Int_t d, Int_t w) {
  std::map<Int_t, mapInfo>::iterator it;
  for (it = GlobalIDToMapInfo.begin(); it != GlobalIDToMapInfo.end(); it++) {
    it->second.FL = FL;
    it->second.FG = FG;
    it->second.d = d;
    it->second.w = w;
  }
  UpdateSettings();
}

void lendaPacker::SetGates(Double_t ilg, Double_t isg,
			   Double_t ilg2, Double_t isg2){
   lg = ilg;  sg = isg;  lg2 = ilg2;  sg2 = isg2;
} 

void lendaPacker::MakeLendaEvent(lendaEvent *Event, ddasEvent *theDDASEvent,
				  Long64_t jentry) {

  vector<ddasChannel*> theDDASChannels = theDDASEvent->getData();
  multimap<Int_t, refTimeContainer> gIDToRefTimes; 

  /* Container to hold the Bars that have fired in the event. 
     Is copied into the lendaEvent. */
  map<string, lendaBar> ThisEventsBars;
  /* Container to hold the Reference channels for the event.  
     Is a multimap to allow for more than one of the same 
     channel to be in one event. */
  multimap<string,lendaChannel> ThisEventsObjectScintillators;

  for (int i=0; i<(int)theDDASChannels.size(); i++) {
    Event->SeteTime(theDDASChannels[0]->eclock);
    Event->SetInternalTime(theDDASChannels[0]->time);

    int id = (theDDASChannels[i]->chanid +
	      CHANPERMOD*(theDDASChannels[i]->slotid-2));
    /* Look for this channel in the MapInfo Map */
    map<int,mapInfo>::iterator it = GlobalIDToMapInfo.find(id);
    string fullName;
    SetJEntry(jentry);  
    /* Check if the channel is in the map. */
    if (it != GlobalIDToMapInfo.end()) { 
      fullName = it->second.FullName;
      /* Special check for the reference channels. */
      if (fullName.find(referenceChannelPattern) != string::npos) { 
	/* If the channel is one of the Object Scintillators */
	lendaChannel Temp = DDASChannel2LendaChannel(theDDASChannels[i],it->second);
	gIDToRefTimes.insert(make_pair(id, refTimeContainer(&Temp)));
	ThisEventsObjectScintillators.insert(make_pair(fullName,Temp));	 
      
      } else if (fullName == "IGNORE") {
	// Do nothing
      } else { 
	/* It is a LENDA BAR */
	string nameOfBar=it->second.BarName;
	/* Check to see if this bar has been found in this event yet. */
	if (ThisEventsBars.count(nameOfBar) == 0) { 
	  /* The bar object will be put in a map to keep track of things */
	  lendaBar tempBar(nameOfBar,it->second.BarAngle);
	  int UniqueBarNum = BarNameToUniqueBarNumber[nameOfBar];
	  tempBar.SetBarId(UniqueBarNum);
	  PutDDASChannelInBar(it->second,tempBar,theDDASChannels[i]);
	  ThisEventsBars[nameOfBar] = tempBar; 
	} else {
	  /* The bar has already had a channel in it from a previous 
	     iteration of this loop over DDAS Channels. Instead of 
	     making a new bar take the already allocated one from the 
	     map and push this channel on to it. */
	  PutDDASChannelInBar(it->second, ThisEventsBars[nameOfBar],
			      theDDASChannels[i]);
	}
      }
    } else { /* This is when the channel was not in the map */
      //Event->PushUnMappedChannel(DDASChannel2LendaChannel(theDDASChannels[i],MapInfo()));
    }
  } /* End loop over DDASChannels */
  
  /////////////////////////////////////////////////////////////////////////////
  // All channels have now been extracted from the DDASEvent.  The channels  //
  // are also now sorted by channelName/barName.	     		     //
  /////////////////////////////////////////////////////////////////////////////
  
  /* Push the reference channels into the lendaEvent */
  
  std::multimap<string,lendaChannel>::iterator it; 
  for (it=ThisEventsObjectScintillators.begin(); it!=ThisEventsObjectScintillators.end(); it++) {
    Event->refScint.push_back((*it).second);
  }
  
  /* Now put the bars into the lenda Event by iterating over the 
     temporary bar map and set the reference times for each lendaChannel */
  FillReferenceTimesInEvent(Event,ThisEventsBars,gIDToRefTimes);

  ThisEventsBars.clear();
  gIDToRefTimes.clear();
  ThisEventsObjectScintillators.clear();
}

void lendaPacker::ReMakeLendaEvent(lendaEvent* inEvent, lendaEvent* outEvent) {

  map<string,lendaBar> ThisEventsBars;
  
  for (int i=0;i<inEvent->numBars;i++) {
    int numTops = inEvent->bars[i].NumTops;
    int numBottoms= inEvent->bars[i].NumBottoms;
  
    for (int t=0;t<numTops;t++) {
      RePackChannel(&inEvent->bars[i].Tops[t]);
    }
    for (int b=0;b<numBottoms;b++) {
      RePackChannel(&inEvent->bars[i].Bottoms[b]);
    }
    ThisEventsBars[inEvent->bars[i].Name] = inEvent->bars[i];
  }

  multimap <int, refTimeContainer >  gIDToRefTimes;
  
  for (int i=0; i<inEvent->numObjectScintillators; i++) {
    outEvent->refScint.push_back(inEvent->refScint[i]);
    lendaChannel* thisObjectScintillator = &(outEvent->refScint[i]);
    RePackChannel(thisObjectScintillator, true);
    gIDToRefTimes.insert(make_pair(thisObjectScintillator->GetGlobalID(),refTimeContainer(thisObjectScintillator)));
  }
  
  FillReferenceTimesInEvent(outEvent,ThisEventsBars,gIDToRefTimes);
  ThisEventsBars.clear();
}

void lendaPacker::RePackChannel(lendaChannel *channel, bool isAnObject){
 
  int GlobalID = channel->GetGlobalID();
  map<int,mapInfo>::iterator it = GlobalIDToMapInfo.find(GlobalID); 
  
  if (it != GlobalIDToMapInfo.end()) { 
    mapInfo info = it->second;    
    Reset(); 
    if (channel->GetTrace().size() !=0) { 
      vector <UShort_t> theTrace = channel->GetTrace();
      CalcAll(theTrace,info);
    }
    /* Set the waveform analysis results to the channel.
       including things like energy cubic times and
       the corrected times and energies. */
    if (info.DontTraceAnalyze==false) {
      PackCalculatedValues(channel,info);
    }
  } else { /* This Global Id has no map information. */ }
  Reset();
}

lendaChannel lendaPacker::DDASChannel2LendaChannel(ddasChannel* c, mapInfo info){
  vector <UShort_t> theTrace = c->GetTrace();
  CalcAll(theTrace,info);
   
  ///////////////////////////////////////////////////////////////////////////
  // CalcAll has done all the waveform analysis and put the results into   //
  // temporary  member variables of lendaPacker.  Now all things (basic    //
  // info and calculated info) will be copied to the lendachannel object.  //
  ///////////////////////////////////////////////////////////////////////////

  lendaChannel templenda;
  
  templenda.SetChannel(c->chanid);
  templenda.SetSlot(c->slotid);
  templenda.SetGlobalID(info.GlobalID);
  templenda.SetReferenceGlobalID(info.ReferenceGlobalID);
  templenda.SetChannelName(info.FullName);
  templenda.SetReferenceChannelName(info.ReferenceName);
  templenda.SetInternalEnergy(c->energy);
  templenda.SetTime(c->time);
  templenda.SetTimeLow(c->timelow);
  templenda.SetTimeHigh(c->timehigh);
  templenda.SetCFDTrigBit(c->GetCFDTriggerSourceBit());
  templenda.SetInternalCFD(c->timecfd/16384.0);
  templenda.SetEclock(c->eclock);
  templenda.SetEclockLow(c->eclocklow);
  templenda.SetEclockHigh(c->eclockhigh);
  templenda.SetEclockHigher(c->eclockhigher);
  if (saveTraces){
    templenda.SetTrace(c->trace);
  }
  if (info.DontTraceAnalyze == false){
    /* If there is not trace analysis for this channel
       don't call PackCalculatedValues. */
    PackCalculatedValues(&templenda,info);
  }
  
  Reset();
  return templenda;
}

void lendaPacker::RePackSoftwareTimes(lendaEvent *Event){
  printf(DMAGENTA "<lendaPacker::RePackSoftwareTimes>: This should not be called.\n" RESET_COLOR);
}

void lendaPacker::PackEvent(lendaEvent * Event){
  printf(DMAGENTA "<lendaPacker::PackEvent> is OBSOLETE!!\n" RESET_COLOR);
}

void lendaPacker::RePackEvent(lendaEvent * Event){
  stringstream stream;
  ifstream MapFile;
  ifstream CorrectionsFile;
  
  stringstream ss,ss1;
  ss << string(getenv("R00TLe_PRM"));
  ss1 << ss.str();
  ss << "/" << mapFileName;
  ss1 << "/" << correctionsFileName;
}

void lendaPacker::FindAndSetMapAndCorrectionsFileNames(int RunNumber) {
  /////////////////////////////////////////////////////////////////////////////
  // Look for MapFiles named MapFileRunNuber.txt and correction files	     //
  // name CorrectionsRunNumber.txt		      			     //
  /////////////////////////////////////////////////////////////////////////////

  stringstream filePath, mapFile, correctionsFile;
  ifstream testFile;

  filePath << string(getenv("LendaSettings"));
  filePath << "/";

  mapFile << filePath.str() << "MapFile" << RunNumber << ".txt";  
  correctionsFile << filePath.str() << "Corrections" << RunNumber << ".txt";  
  
  string corName = "Corrections.txt";
  string mapName = "MapFile.txt";

  testFile.open(mapFile.str().c_str());
  if (testFile.good()) {
    stringstream s;
    s << "MapFile" << RunNumber << ".txt";
    mapName = s.str();
  }
  testFile.close();
  
  testFile.open(correctionsFile.str().c_str());
  if (testFile.good()) {
    stringstream s;
    s << "Corrections" << RunNumber << ".txt";
    corName = s.str();
  }
  testFile.close();
  
  SetSettingFileNames(mapName,corName);
}

void lendaPacker::SetSettingsandCorrections(string MapFileName,
					    string CorrectionsFileName){
  already_called = 100;
  stringstream filePath, mapfilename, correctionfilename;
  filePath << "./";
  mapfilename << filePath.str() << MapFileName;
  correctionfilename << filePath.str() << CorrectionsFileName;  
  ifstream testFile;
  testFile.open(mapfilename.str().c_str());
  if (!testFile.good()) {
    printf(DMAGENTA "ERROR: Could not find %s\n" RESET_COLOR, mapfilename.str().c_str());
  }
  testFile.close();
  testFile.open(correctionfilename.str().c_str());
  if(!testFile.good()) {
    printf(DMAGENTA "ERROR: Could not find %s\n" RESET_COLOR, correctionfilename.str().c_str());
  }
  mapFileName = mapfilename.str();
  correctionsFileName = correctionfilename.str();
  Info("lendaPacker","Using MapFile %s and CorrectionsFile %s",
       mapFileName.c_str(), correctionsFileName.c_str());
  BuildMaps();
}

mapInfo lendaPacker::GetMapInfo(string FullName) {
  return GlobalIDToMapInfo[FullLocalToGlobalID[FullName]];
}

void lendaPacker::BuildMaps() {
  
  stringstream stream;
  ifstream MapFile, CorrectionsFile;

  stringstream ss,ss1;
  ss << mapFileName;
  ss1 << correctionsFileName;

  MapFile.open(ss.str().c_str());
  if (!MapFile.is_open()) {
    Error("lendaPacker", "Cannot find Map file.  Looking in %s",
	  ss.str().c_str());
    throw -1;
  }

  CorrectionsFile.open(ss1.str().c_str());
  if ( ! CorrectionsFile.is_open()) {
    Error("lendaPacker", "Cannot find Corrections file. Looking in %s",
	  ss1.str().c_str());
    throw -1;
  }

  /* Read in the map information -- should be 5 columns. 
           slot channel Name Angle ReferenceName           */
  int slot,channel;
  string name,ReferenceName;
  double Angle;
  int UniqueBarNumber=0;
   
  std::string line;
  while (std::getline(MapFile>>std::ws, line)) {
    int FL = 0, FG = 0, d = 0, w = 0;
    int dontTrace;
    if (line[0] != '#') {
      std::istringstream iss(line);
      iss >> slot >> channel >> name >> Angle 
	  >> ReferenceName >> FL >> FG >> d >> w >> dontTrace;
  
      /* Make Global ID */
      int spot = CHANPERMOD*(slot-2) + channel;
       
      mapInfo tempInfo;
      if(name.find(referenceChannelPattern) != string::npos) {
	/* It is a reference channel.  Set flag in mapinfo for later. */
	tempInfo.IsAReferenceChannel = true;
      } else{
	tempInfo.IsAReferenceChannel = false;
      }
      GlobalIDToFullLocal[spot] = name;
      FullLocalToGlobalID[name] = spot;
      string BarName = name.substr(0,name.size()-1);
      GlobalIDToBar[spot] = BarName;
      
      tempInfo.FL = FL;
      tempInfo.FG = FG;
      tempInfo.d = d;
      tempInfo.w = w;
      tempInfo.DontTraceAnalyze = dontTrace;
      tempInfo.GlobalID = spot;
      tempInfo.FullName = name;
      tempInfo.BarName = BarName;
      tempInfo.ReferenceName = ReferenceName;
      tempInfo.BarAngle = Angle;
      if (Angle !=0) {
	BarNameToBarAngle[BarName] = Angle;
      }
      if (BarNameToUniqueBarNumber.count(BarName)==0){ /* Isn't already there */
	BarNameToUniqueBarNumber[BarName]=UniqueBarNumber;
	UniqueBarNumber++;
      } 
      GlobalIDToMapInfo[spot] = tempInfo;
    } /* end if not a comment */
  } /* end while for Map file */

  Double_t slope,intercept,timeOffSet;
  while (std::getline(CorrectionsFile>>std::ws, line)) {
    if (line[0] != '#') {
      std::istringstream iss(line);
      iss >> name >> slope >> intercept >> timeOffSet;
      if (FullLocalToGlobalID.count(name) == 0) {
	/* There is a name in the corrections file that 
	   isn't in the Cable Map file. */  
	cout << "Found a name in the corrections file that" << endl;
	cout << "wasn't in the cable map file" << endl;
	cout << "Name is " << name << endl;
	throw -99;
      }
      int GlobalID = FullLocalToGlobalID[name];
      mapInfo * temp = &GlobalIDToMapInfo[GlobalID];
      temp->EnergySlope = slope;
      temp->EnergyIntercept = intercept;
      temp->TOFOffset = timeOffSet;
      temp->HasCorrections = true;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////
  // Now that the maps have been built.  Loop over the main conatainer    //
  // and determine the global ID of the reference name.  If it is not in  //
  // the map throw error	                                          //
  //////////////////////////////////////////////////////////////////////////

  std::map<Int_t, mapInfo>::iterator it;
  for (it=GlobalIDToMapInfo.begin(); it != GlobalIDToMapInfo.end(); it++) {
    int GlobalID = it->first;
    if (FullLocalToGlobalID.count(it->second.ReferenceName) !=0) {
      int RefGlobalID = FullLocalToGlobalID[it->second.ReferenceName];
      it->second.GlobalID=GlobalID;
      it->second.ReferenceGlobalID=RefGlobalID;
    } else {
      printf(DMAGENTA "Found a reference name in the map file that does not map to a channel.\n");
      printf("The name was %s -- it is from map info of %s\n" RESET_COLOR, it->second.ReferenceName.c_str(), it->second.FullName.c_str());
      throw -12;
    }

    if (it->second.BarAngle == 0) { /* This channel doesn't have a Bar angle.
				      Get angle from the BarNameToBarAngle Map
				      This rigamaroll is to let you type the 
				      angle of the bar once in the bar map.  
				      Instead of doing it both for the top 
				      and bottom channel.  One can be 0. */
      it->second.BarAngle = BarNameToBarAngle[it->second.BarName];
    }
  }

  UpdateSettings();
  for (it=GlobalIDToMapInfo.begin(); it != GlobalIDToMapInfo.end(); it++) {
    it->second.Print(); 
  }
  it = GlobalIDToMapInfo.find(1);
  std::cout<< "ID MAP INFO:" << endl;
  it->second.Print();
  std::cout << " SIZE: " << GlobalIDToMapInfo.size();
}

void lendaPacker::PutDDASChannelInBar(mapInfo info, lendaBar &theBar,
				      ddasChannel*theChannel) {
  string fullName = info.FullName;
  string lastLetter = fullName.substr(fullName.size()-1,1);

  /* All we need to do is determine whether this fullName coresponds to
     the TOP or bottom PMT of the bar */
  if (lastLetter == "T" ) { /* Top PMT */
    lendaChannel temp = DDASChannel2LendaChannel(theChannel,info);
    theBar.Tops.push_back(temp);
  } else if (lastLetter == "B") { /* Bottom PMT */
    lendaChannel temp =DDASChannel2LendaChannel(theChannel,info);
    theBar.Bottoms.push_back(temp);
  } else {
    printf(DMAGENTA "*************************************************************\n");
    printf("Found a ddasChannel name without a T or B as the last letter!\n");
    printf("***** Cable map file must have names that end in T or B *****\n");
    printf("***** Name was: %s\n", fullName.c_str());
    printf("*************************************************************\n" RESET_COLOR);
    throw -99;
  } 
}

void lendaPacker::PackCalculatedValues(lendaChannel* theChannel,
				       mapInfo & info){
  Double_t CoarseTime = (2*(theChannel->GetTimeLow() +
			    theChannel->GetTimeHigh() * 4294967296.0));
  
  theChannel->SetSoftwareCFD(thisChannelsSoftwareCFDs[0]);
  theChannel->SetCubicCFD(thisChannelsCubicCFDs[0]);
  theChannel->SetCubicFitCFD(cubicFitCFD);
  theChannel->SetFilterHeight(thisEventsFilterHeight);

  /* Now set the vectors for when this channel was a reference channel */
  vector <Double_t> TempTimes;
  vector <Double_t> TempCubicTimes;
  for (int i=0;i<(int)thisChannelsSoftwareCFDs.size();i++) {
    TempTimes.push_back( CoarseTime + thisChannelsSoftwareCFDs[i]);
    TempCubicTimes.push_back( CoarseTime + thisChannelsCubicCFDs[i]);
  }
  
  /* Setting the vectors for the Times/energies will allow one to
     call GetEnergy() or GetEnergies() from the lendachannel.
     GetEnergy() will just return the first one. */
  theChannel->SetEnergies(thisChannelsEnergies);
  theChannel->SetSoftwareTimes(TempTimes);
  theChannel->SetCubicTimes(TempCubicTimes);
  theChannel->SetPulseHeights(thisChannelsPulseHeights);
  theChannel->SetUnderShoots(thisChannelsUnderShoots);
  
  if (saveTraces) {
    theChannel->SetFilter(thisEventsFF);
    theChannel->SetCFD(thisEventsCFD);
  }
  theChannel->SetLongGate(longGate);
  theChannel->SetShortGate(shortGate);
  
  theChannel->SetJentry(jentry);
  theChannel->SetNumZeroCrossings(thisChannelsNumZeroCrossings);
  theChannel->SetCFDResidual(thisChannelsCFDResidual);
 
  if (info.HasCorrections) {
    theChannel->SetCorrectedEnergy( theChannel->GetEnergy()*info.EnergySlope + info.EnergyIntercept);
    theChannel->SetCorrectedTime( theChannel->GetTime() - info.TOFOffset);
    theChannel->SetCorrectedSoftTime(theChannel->GetSoftTime()-info.TOFOffset);
    theChannel->SetCorrectedCubicFitTime( theChannel->GetCubicFitTime() - info.TOFOffset);
    theChannel->SetCorrectedCubicTime( theChannel->GetCubicTime() - info.TOFOffset);
  } 
}

void lendaPacker::FillReferenceTimesInEvent(lendaEvent* Event,
					    map<string,lendaBar>& ThisEventsBars,
					    multimap <int,refTimeContainer>& GlobalIDToReferenceTimes) {


  std::pair<std::multimap<int, refTimeContainer>::iterator, std::multimap<int, refTimeContainer>::iterator> MultiMapRange;
  multimap<int, refTimeContainer>::iterator MultiMapIt;

  /* Loop over the ThisEventsBars container and set all the reference times */
  
  for (map<string,lendaBar>::iterator ii=ThisEventsBars.begin();
       ii!=ThisEventsBars.end();ii++) {
    
    /* If this event had TWO copies of OBJ1T for example both should 
       be in the multimap and we can loop over all timestamps in BOTH 
       and push those into the Bar's list of reference times */
    /* Top signals */
    for (int t=0; t<ii->second.Tops.size(); t++) {
      MultiMapRange = GlobalIDToReferenceTimes.equal_range(ii->second.Tops[t].GetReferenceGlobalID());  
      for (MultiMapIt = MultiMapRange.first; MultiMapIt!=MultiMapRange.second;  ++MultiMapIt) {
	refTimeContainer * temp = &MultiMapIt->second;
	ii->second.Tops[t].SetReferenceTime(temp->RefTime);
	ii->second.Tops[t].SetSoftwareReferenceTimes(temp->RefSoftTime);
	ii->second.Tops[t].SetCubicReferenceTimes(temp->RefCubicTime);
      }
    }

    /* Bottom signals */
    for (int b=0; b<ii->second.Bottoms.size(); b++) {
      MultiMapRange = GlobalIDToReferenceTimes.equal_range(ii->second.Bottoms[b].GetReferenceGlobalID());
      for (MultiMapIt = MultiMapRange.first; MultiMapIt!=MultiMapRange.second;  ++MultiMapIt) {
	refTimeContainer * temp = &MultiMapIt->second;
	ii->second.Bottoms[b].SetReferenceTime(temp->RefTime);
	ii->second.Bottoms[b].SetSoftwareReferenceTimes(temp->RefSoftTime);
	ii->second.Bottoms[b].SetCubicReferenceTimes(temp->RefCubicTime);
      }
    }
    Event->PushABar(ii->second);
  }
}

void lendaPacker::UpdateSettings() {
  std::map<Int_t, mapInfo>::iterator it;
  for (it = GlobalIDToMapInfo.begin(); it != GlobalIDToMapInfo.end(); it++) {
    int GlobalID = it->first;
    theSettings->AddCorrectionSettings(it->second.FullName,
				       it->second.EnergySlope,
				       it->second.EnergyIntercept,
				       it->second.TOFOffset);
    theSettings->AddMapSettings(it->second.FullName, GlobalID,
				it->second.ReferenceName,
				it->second.ReferenceGlobalID);
    theSettings->AddFilterSettings(it->second.FullName, it->second.FL,
				   it->second.FG, it->second.d, it->second.w,
				   it->second.DontTraceAnalyze);
  }
  theSettings->SetBarIds(BarNameToUniqueBarNumber);
}
