#ifndef __LENDACONTROL_H
#define __LENDACONTROL_H

/* Taken from "LendaPacker.hh" by Sam Lipschutz -- Copyright 2014 */

#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include "TFile.h"

#include "colors.h"

#include "LENDA-DDAS.h"

#define CHANPERMOD 16  /* Number of DDAS chans per DDAS module */

class lendaSettings : public TNamed {
 private:
  map<string, Double_t> TheTimingOffsets;
  map<string, Double_t> TheEnergySlopes;
  map<string, Double_t> TheEnergyIntercepts;

  map<string, Int_t> TheFLs;
  map<string,Int_t> TheFGs;
  map<string,Int_t> Theds;
  map<string,Int_t> Thews;
  map<string,bool> TheDontTraceAnalyzeFlags;
  
  map<Int_t, string> GlobalID2FullName;
  map<Int_t, Int_t> GlobalID2RefGlobalID;
  map<Int_t,string> GlobalID2RefName;

  map<string,Int_t> Name2GlobalID;

  map<string,Int_t> BarIds;
  map<Int_t,string> BarId2Name;

 public:
  lendaSettings();
  ~lendaSettings();
  
  void PrintAll();
  void PrintChannelCorrections(string Name);
  void PrintChannelMapInfo(Int_t GlobalID);
  void PrintChannelMapInfo(string Name);
  void PrintFilterInfo(string Name);

  Double_t GetChannelsSlope(string Name) { return TheEnergySlopes[Name]; }
  Double_t GetChannelsIntercept(string Name) { return TheEnergyIntercepts[Name]; }
  Double_t GetChannelsTimeOffset(string Name) { return TheTimingOffsets[Name]; }
  
  Int_t GetBarId(string Name) { return BarIds[Name]; }
  string GetBarName(Int_t BarId) { return BarId2Name[BarId]; }

  void AddCorrectionSettings(string Name,Double_t slope,Double_t inter,Double_t toff);
  void AddMapSettings(string Name,Int_t GlobalID,string RefName, Int_t refGlobalID);
  void AddFilterSettings(string Name,Int_t FL,Int_t FG,Int_t d,Int_t w,bool flag);
 
  void SetBarIds(map<string,Int_t> v);
  
  inline Int_t GetNumBars() { return BarIds.size(); }
  void BuildReverseMap();

 public:
  ClassDef(lendaSettings, 1);
};

/* MapInfo is a small Container class to hold the information 
   that is stored in the Mapfile/Corrections file. */

class mapInfo {
 public:
  Double_t EnergySlope; 
  Double_t EnergyIntercept;
  Double_t TOFOffset; 
  
  string FullName; 
  string BarName; 
  string ReferenceName; 

  Double_t BarAngle;
  bool HasCorrections;
  bool IsAReferenceChannel;

  Int_t GlobalID;
  Int_t ReferenceGlobalID;

  Int_t FL; /* Timing filter Rise length */
  Int_t FG; /* Timing filter Gap Length */
  Int_t d; /* Timing filter CFD delay */
  Int_t w; /* Timing filter CFD scale factor */
  
  bool DontTraceAnalyze;
  
 public:
 mapInfo(): EnergySlope(BADNUM), EnergyIntercept(BADNUM),
    TOFOffset(BADNUM), FullName(""), BarName(""), ReferenceName(""),
    BarAngle(BADNUM), HasCorrections(false), IsAReferenceChannel(false),
    ReferenceGlobalID(BADNUM), GlobalID(BADNUM), FL(BADNUM),
    FG(BADNUM), d(BADNUM), w(BADNUM){ ; }
  
  void Print() {
    cout << "**** NEW BAR: ******" << endl;
    cout << "	** Informations of MapFile:" << endl; 
    cout<<  "	FullName "<<FullName<<" \t Bar name: "<<BarName<<"\t Global Id: "<<GlobalID << endl;
    cout << "	Reference name: "<<ReferenceName<<"\t RefGlobal: "<<ReferenceGlobalID<<endl;
    cout << "	** Informations of CorrectionFile:" << endl;
    cout<<	"	Pulse Height Slope: "<<EnergySlope<<" \t Pulse Height Intercept: "<<EnergyIntercept<<" \t Timing Offset: "<<TOFOffset<<" \t Bar Angle: "<<BarAngle<<endl;
    cout<<	"	FL: "<<FL<<" \t FG: "<<FG<<"\t d: "<<d<<"\t w: "<<w<<endl;
  }

};

/* refTimeContainer is a small container class.  It will hold the 
   different types of times that will be set as the reference times 
   for a particular channel. */

class refTimeContainer {
 public:
  Double_t RefTime; 
  vector<Double_t> RefSoftTime;
  vector<Double_t> RefCubicTime;
  
 public:
  refTimeContainer(lendaChannel* c) {
    RefTime = c->GetTime();
    RefSoftTime = c->GetSoftwareTimes();
    RefCubicTime = c->GetCubicTimes();
    if (RefSoftTime.size() == 0) {
      throw "SOFT";
    }
    if (RefCubicTime.size() == 0) {
      throw "CUBIC";
    }
  }
};

/*  The main purpose of LendaPacker class is to take DDASEvents and 
    make LendaEvents.  The Packer is responsible for performing the 
    various trace analysis routines (defined in the lendaFilter class) 
    and saving them in the lendaEvent. 
    
    The Packer needs a cable map file which tells it which   
    DDAS channels to associate together into lendaBars.  Will Throw 
    exception if it cannot find this file.  Currently only supports single 
    crate setups.            
                                                                                   
    The map file also provides the following information to packer:  
       The scattering angle that the bar is at at, which other channel 
       provides the time of flight reference time.  The format of the 
       MapFile is:                                    
                                                                                  
       DDAS_slot_number DDAS_channel_number Bar_name_ending_in_T/B \\
       angle Name_of_refTime 
                                                                                  
    The reference time name must also be somewhere else in the mapfile 
    so that the actual DDAS slot and channel can be determined for it.  
    
    Written by Sam Lipschutz Copyright 2014  			          */

class lendaPacker {

 public:
  lendaSettings * theSettings;
  lendaFilter theFilter;

  Int_t already_called;
  string mapFileName, correctionsFileName;

  Bool_t saveTraces;

  Double_t lg, sg, lg2, sg2;
  Int_t traceDelay;
  Long64_t jentry;

  Int_t thisChannelsNumZeroCrossings;
  
  vector <Double_t> thisChannelsEnergies;
  vector <Int_t> thisChannelsPulseHeights;
  
  vector <Double_t> thisChannelsSoftwareCFDs;
  vector <Double_t> thisChannelsCubicCFDs;

  vector <Int_t> thisChannelsPeakSpots;
  vector <Double_t> thisChannelsUnderShoots;

  vector <Double_t> thisEventsFF; 
  vector <Double_t> thisEventsCFD;
  
  map<Int_t, string> GlobalIDToFullLocal;
  map<string, Int_t> FullLocalToGlobalID;
  map<Int_t, string> GlobalIDToBar;
  map<string, Int_t> BarNameToUniqueBarNumber;
  map<string, Double_t> BarNameToBarAngle;  
  map<Int_t, mapInfo > GlobalIDToMapInfo;
  
  Int_t thisEventsFilterHeight;

  Double_t longGate;
  Double_t shortGate;

  Double_t cubicFitCFD;

  Double_t start;
  Double_t thisChannelsCFDResidual;

  string referenceChannelPattern;
  
 public:
  lendaPacker() { ; }
  lendaPacker(lendaSettings*);
  ~lendaPacker();

  /* Set the filenames for Corrections/Map file.  Calls BuildMaps. */
  void SetSettingFileNames(string MapFileName, string CorrectionsFileName);  
  void Reset();

  void CalcTimeFilters(vector<UShort_t>& theTrace, mapInfo info);
  void CalcEnergyGates(vector<UShort_t>& theTrace, mapInfo info);
  void CalcAll(vector<UShort_t>& theTrace, mapInfo info);
  void CalcAll(ddasChannel* c, mapInfo info);
  void CalcAll(lendaChannel* c, mapInfo info);

  void ForceAllBarFilters(Int_t FL, Int_t FG, Int_t d, Int_t w);
  void ForceAllReferenceFilters(Int_t FL, Int_t FG, Int_t d, Int_t w);
  void ForceAllFilters(Int_t FL, Int_t FG, Int_t d, Int_t w);

  void SetGates(Double_t, Double_t, Double_t, Double_t);

  inline void SetTraceDelay(Int_t x) { traceDelay=x; }
  inline void SetJEntry(Long64_t n) { jentry=n; }
  
  /* Main method used in building lendaEvents.
     Takes the ddasChannels in the DDASEvent and maps them to the appropriate
     detectors.  Performs all waveform analysis and saves infomration in 
     a lendaEvent. */
  void MakeLendaEvent(lendaEvent *Event, ddasEvent *theDDASEvent,
		      Long64_t jentry);
  
  /** Take a lendaEvent and rebuild the event.
      Used to recalculate timing algorithms with new parameters.
      Does NOT re-map the channels.  */
  void ReMakeLendaEvent(lendaEvent*inEvent, lendaEvent*outEvent);
  void RePackChannel(lendaChannel *, bool isAnObject = false);

  /* Convert a DDASChannel to a LendaChannel.  Will call CalcAllMethod
     and save waveform analysis results into the LendaChannel. */
  lendaChannel DDASChannel2LendaChannel(ddasChannel* c, mapInfo info);
  
  /* OBSOLETE!!! */
  void RePackSoftwareTimes(lendaEvent *Event);
  void PackEvent(lendaEvent* Event);
  void RePackEvent(lendaEvent* Event);

  void FindAndSetMapAndCorrectionsFileNames(Int_t RunNumber);
  void SetSettingsandCorrections(string MapFileName,string CorrectionsFileName);

  mapInfo GetMapInfo(string FullName);

private:
  void BuildMaps();  
  void PutDDASChannelInBar(mapInfo info, lendaBar &theBar,
			   ddasChannel *theChannel);
  void PackCalculatedValues(lendaChannel*, mapInfo &);

  void FillReferenceTimesInEvent(lendaEvent* Event,
				 map<string, lendaBar>& ThisEventsBars,
				 multimap<Int_t, refTimeContainer>& GlobalIDToReferenceTimes);

  void UpdateSettings();

};

#endif
