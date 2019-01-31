#ifndef __LENDA_H
#define __LENDA_H

#include <iostream>
#include <map>
#include <sstream>
#include <stdint.h>
#include <string>
#include <vector>

#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TMath.h"
#include "TMatrixT.h"
#include "TObject.h"

#include "ddasChannel.h"

#define BADNUM -10008

/* Class for correlated ddasChannel objects */
class ddasEvent : public TObject {
 private:
  std::vector<ddasChannel*> cdata;
  ddasChannel *tempdchan;
  
  uint8_t buffer08[4096];
  uint32_t buffer32[4096];

 public:
  ddasEvent();
  ddasEvent(const ddasEvent& obj);
  ddasEvent& operator=(const ddasEvent& obj);
  ~ddasEvent();

  std::vector<ddasChannel*>& getData() { return cdata; }
  void addChannelData(ddasChannel* channel);
  UInt_t getNEvents() const { return cdata.size(); }
  Double_t getFirstTime() const;
  Double_t getLastTime() const;
  Double_t getTimeWidth() const;
  void Reset();

  void getEvent(FILE *inf, Int_t length);

  ClassDef(ddasEvent, 1);
};

/* lendaChannel -- lowest level container for DDAS Lenda information.  
   Each ddasChannel will have a corresponding lendaChannel.  Channel 
   by channel information for the module and the results of waveform 
   analysis are saved in lendaChannels. */

class lendaChannel : public TObject {
 private:
  Int_t rChannel, rSlot, rGlobalID, rRefGlobalID;
  Bool_t rOverflow;
  Int_t rNumZeroCrossings;

  Double_t rCorrectedEnergy, rInternalEnergy;
  Int_t rFilterHeight;

  Double_t rTime, rCubicFitTime;
  Double_t rCorrectedTime, rCorrectedCubicFitTime;
  Double_t rCorrectedSoftTime, rCorrectedCubicTime;

  UInt_t rTimeLow, rTimeHigh;
  
  UInt_t rEclockLow, rEclockHigh, rEclockHigher;
  Double_t rEclock;

  Double_t rCoarseTime;
  
  /* CFD Values for sub clock tick timing.
     These cfds refer to the cfd from calculated by
     LendaFilter.GetZeroCrossingImproved which 
     will give the value for the largest zero-crossing. */
  Double_t rSoftwareCFD, rCubicCFD, rCubicFitCFD;
  /* Internal CFD is for wherever the module triggered on. */
  Double_t rInternalCFD;
  UInt_t rCFDTrigBit;

  Double_t rShortGate, rLongGate;
  Long64_t rJentry;
  
  vector <UShort_t> rTrace;
  vector <Double_t> rFilter;
  vector <Double_t> rCFD;

  Double_t rCFDResidual;

  string rChannelName;
  string rReferenceChannelName;

  Double_t rRefTime;
  // Double_t rSoftRefTime;
  // Double_t rCubicRefTime;

  Double_t rOtherTime;

  vector <Double_t> rEnergies;
  vector <Double_t> rSoftwareTimes;
  vector <Double_t> rCubicTimes;
  vector <Int_t> rPulseHeights;
  vector <Double_t> rSoftwareReferenceTimes;
  vector <Double_t> rCubicReferenceTimes;
  vector <Double_t> rUnderShoots;

 public:
  lendaChannel(); 
  ~lendaChannel();
  Bool_t operator==(const lendaChannel & RHS);
  
  /* Resets all the values in the lendaChannel.  Should be called after every write to a tree. */
  void Clear(); 
  void Clear(Option_t *){ Clear(); } 

  void SetEnergies(const vector<Double_t> & v) { rEnergies=v; } 
  void SetSoftwareTimes(const vector<Double_t> & v) { rSoftwareTimes=v; } 
  void SetCubicTimes(const vector<Double_t> & v) { rCubicTimes=v; } 
  void SetPulseHeights(const vector<Int_t> & v) { rPulseHeights=v; }

  void SetSoftwareReferenceTimes(const vector<Double_t> & v); 
  void SetCubicReferenceTimes(const vector<Double_t> & v); 

  vector<Double_t> GetEnergies() const { return rEnergies; }
  vector<Double_t> GetSoftwareTimes() const { return rSoftwareTimes; }
  vector<Double_t> GetCubicTimes() const { return rCubicTimes; }
  vector<Int_t> GetPulseHeights() const { return rPulseHeights; }
  vector<Double_t> GetSoftwareReferenceTimes() const { return rSoftwareReferenceTimes; }
  vector<Double_t> GetCubicReferenceTimes() const { return rCubicReferenceTimes; }

  /* Gets the undershoot value for trace.  If an undershoot was present it will give the difference between
     the peak of the under shoot and the calculated baseLine.  Will store information in a vector to accomidate
     many different pulses in the trace. */
  vector <Double_t> GetUnderShoots() const { return rUnderShoots; }
  void SetUnderShoots(vector <Double_t> v){ rUnderShoots=v; }
  
  void SetEclockLow(UInt_t v) { rEclockLow=v; }
  void SetEclockHigh(UInt_t v) { rEclockHigh=v; }
  void SetEclockHigher(UInt_t v) { rEclockHigher=v; }
  void SetEclock(Double_t v) { rEclock=v; }
  inline Int_t GetEclockLow() { return rEclockLow; }
  inline Int_t  GetEclockHigh() { return rEclockHigh; }
  inline Int_t GetEclockHigher() { return rEclockHigher; }
  inline Double_t  GetEclock() { return rEclock; }
  
  inline Int_t GetChannel() { return rChannel; }
  inline Int_t GetSlot() { return rSlot; }
  inline Int_t GetGlobalID() { return rGlobalID; }
  inline Int_t GetReferenceGlobalID() { return rRefGlobalID; }

  inline Bool_t GetOverFlow() { return rOverflow; }
  inline Int_t GetNumZeroCrossings() { return rNumZeroCrossings; }

  inline Double_t GetCorrectedEnergy() { return rCorrectedEnergy; }
  inline Double_t GetEnergy() { return rEnergies[0]; }
  inline Double_t GetInternalEnergy() { return rInternalEnergy; }
  inline Int_t GetPulseHeight() { return rPulseHeights[0]; }
  inline Int_t GetFilterHeight() { return rFilterHeight; }

  inline Double_t GetTime() { return rTime; }
  inline Double_t GetSoftTime() { return rSoftwareTimes[0]; }
  inline Double_t GetCubicTime() { return rCubicTimes[0]; }
  inline Double_t GetCubicFitTime() { return rCubicFitTime; }

  inline Double_t GetCorrectedTime() { return rCorrectedTime; }
  inline Double_t GetCorrectedCubicTime() { return rCorrectedCubicTime; }
  inline Double_t GetCorrectedCubicFitTime() { return rCorrectedCubicFitTime; }
  inline Double_t GetCorrectedSoftTime() { return rCorrectedSoftTime; }

  inline UInt_t GetCFDTrigBit() { return rCFDTrigBit; }
  
  inline Double_t GetSoftwareCFD() { return rSoftwareCFD; }
  inline Double_t GetCubicCFD() { return rCubicCFD; }
  inline Double_t GetCubicFitCFD() { return rCubicFitCFD; }
  inline Double_t GetInternalCFD() { return rInternalCFD; }

  inline Double_t GetShortGate() { return rShortGate; }
  inline Double_t GetLongGate() { return rLongGate; }

  inline Long64_t GetJentry() { return rJentry; }

  inline UInt_t GetTimeLow() { return rTimeLow; }
  inline UInt_t GetTimeHigh() { return rTimeHigh; }

  vector <UShort_t> GetTrace() { return rTrace; }
  vector <Double_t> GetFilter() { return rFilter; }
  vector <Double_t> GetCFD() { return rCFD; }

  inline Double_t GetCFDResidual() { return rCFDResidual; }

  inline Double_t GetReferenceTime() { return rRefTime; } 
  inline Double_t GetCubicReferenceTime() { return rCubicReferenceTimes[0]; }
  inline Double_t GetSoftwareReferenceTime() { return rSoftwareReferenceTimes[0]; }

  string GetChannelName() { return rChannelName; }
  string GetName() { return rChannelName; }
  string GetReferenceChannelName() { return rReferenceChannelName; }

  Double_t GetCoarseTime() { return 2*(GetTimeLow() + GetTimeHigh() * 4294967296.0); }

  inline Double_t GetOtherTime() { return rOtherTime; }

  void SetChannel(Int_t v) { rChannel=v; }
  void SetSlot(Int_t v) { rSlot=v; }
  void SetGlobalID(Int_t v) { rGlobalID=v; }
  void SetReferenceGlobalID(Int_t v) { rRefGlobalID=v; }

  void SetOverFlow(Bool_t v) { rOverflow=v; }
  void SetNumZeroCrossings(Int_t v) { rNumZeroCrossings=v; }
  
  void SetCorrectedEnergy(Double_t v) { rCorrectedEnergy=v; }

  void SetInternalEnergy(Double_t v) { rInternalEnergy=v; }

  void SetFilterHeight(Int_t v) { rFilterHeight=v; }

  void SetTime(Double_t v) { rTime=v; }
  void SetCubicFitTime(Double_t v) { rCubicFitTime=v; }

  void SetCorrectedTime(Double_t v) { rCorrectedTime=v; }
  void SetCorrectedCubicFitTime(Double_t v) { rCorrectedCubicFitTime=v; }
  void SetCorrectedSoftTime(Double_t v) { rCorrectedSoftTime=v; }
  void SetCorrectedCubicTime(Double_t v) { rCorrectedCubicTime=v; }
  
  void SetTimeLow(UInt_t v) { rTimeLow=v; }
  void SetTimeHigh(UInt_t v) { rTimeHigh=v; }

  void SetCFDTrigBit(UInt_t v) { rCFDTrigBit=v; }
   
  void SetSoftwareCFD(Double_t v) { rSoftwareCFD=v; }
  void SetCubicCFD(Double_t v) { rCubicCFD=v; }
  void SetCubicFitCFD(Double_t v) { rCubicFitCFD=v; }
  void SetInternalCFD(Double_t v) { rInternalCFD=v; }

  void SetShortGate(Double_t v) { rShortGate=v; }
  void SetLongGate(Double_t v) { rLongGate=v; }
  
  void SetJentry(Long64_t v) { rJentry=v; }

  void SetTrace(vector<UShort_t> v) { rTrace=v; }
  void SetFilter(vector <Double_t> v) { rFilter=v; }
  void SetCFD(vector <Double_t> v) { rCFD=v; }

  void SetCFDResidual(Double_t v) { rCFDResidual=v; }

  void SetChannelName(string s) { rChannelName=s; }
  void SetReferenceChannelName(string s) { rReferenceChannelName=s; }
  
  void SetReferenceTime(Double_t v) { rRefTime=v; }
  void SetOtherTime(Double_t v) { rOtherTime=v; }
    
 public:
  ClassDef(lendaChannel, 1);
};

/* lendaBar -- holds information from individual Lenda bars.  Stored
   within lendaEvents, contains lendaChannel objects.  Each bar has
   a name and unique bar ID. */

class lendaBar : public TObject {
 public:
  string Name;
  Int_t BarId;
  Double_t BarAngle;
  
  Int_t NumTops;
  Int_t NumBottoms;
  Int_t BarMultiplicity;
  
  Bool_t SimpleEventBit;

  /* Storing a vector of Tops and bottoms to allow for 0 or more than one of each channel to be stored */
  vector <lendaChannel> Tops;
  vector <lendaChannel> Bottoms;
  
 private:
  Double_t Dt; /* Top time - Bottom time */
  Double_t CubicDt;

  /* The top-bottom time difference through the reference signal. */
  Double_t SyncDt;
  Double_t CubicSyncDt;

  Double_t COG; /* ET - EB/(ET+EB) */
  Double_t AvgEnergy; /* Sqrt(ET*EB) */
  Double_t AvgPulseHeight; /* Sqrt(PH_T*PH_B) */

  Double_t AvgT, CorrectedAvgT;

  Double_t CorrectedCOG; /* Gain corrected COG */
  Double_t CorrectedDt; //Time Offset Corrected Dt
  Double_t CorrectedAvgEnergy; //Gain Corrected AvgE

  //Time of flight calculated from the Top/Bottom channel
  Double_t TopTOF, BottomTOF; 

  Double_t CubicTopTOF, CubicBottomTOF;
  Double_t CorrectedCubicTopTOF, CorrectedCubicBottomTOF;
  Double_t CorrectedCubicFitTopTOF, CorrectedCubicFitBottomTOF;
  Double_t CorrectedSoftTopTOF, CorrectedSoftBottomTOF;
  Double_t CorrectedTopTOF, CorrectedBottomTOF;

 public:
  lendaBar();
  lendaBar(string s, Double_t a) { Clear(); Name=s; BarAngle=a;} /* Constrcutor specifying name and angle of bar. */
  ~lendaBar();
  
  void PushTopChannel(lendaChannel c) { Tops.push_back(c); }
  void PushBottomChannel(lendaChannel c) { Bottoms.push_back(c); }

  inline string GetBarName() { return Name; }
  inline Int_t GetNumTops() { return NumTops; }
  inline Int_t GetNumBottoms() { return NumBottoms; }
  inline Int_t GetBarMultiplicity() { return BarMultiplicity; }

  inline Double_t GetCubicDt() { return CubicDt; }
  inline Double_t GetDt() { return Dt; }

  inline Double_t GetSyncDt() { return SyncDt; }
  inline Double_t GetCubicSyncDt() { return CubicSyncDt; }

  inline Double_t GetCOG() { return COG; }
  inline Double_t GetAvgEnergy() { return AvgEnergy; } 
  inline Double_t GetAvgPulseHeight() { return AvgPulseHeight; }

  inline Double_t GetCorrectedDt() { return CorrectedDt; }
  inline Double_t GetCorrectedCOG() { return CorrectedCOG; }
  inline Double_t GetCorrectedAvgEnergy() { return CorrectedAvgEnergy; }

  inline Double_t GetAvgT() { return AvgT; }
  inline Double_t GetCorrectedAvgT() { return AvgT; }
  
  inline Double_t GetTopTOF() { return TopTOF; }
  inline Double_t GetBottomTOF() { return BottomTOF; }

  inline Double_t GetAvgTOF() { return 0.5*(TopTOF+BottomTOF); }
  inline Double_t GetCorrectedAvgTOF() { return 0.5*(CorrectedTopTOF + CorrectedBottomTOF); }

  inline Double_t GetCubicTopTOF() { return CubicTopTOF; }
  inline Double_t GetCubicBottomTOF() { return CubicBottomTOF; }

  inline Double_t GetCorrectedTopTOF() { return CorrectedTopTOF; }
  inline Double_t GetCorrectedBottomTOF() { return CorrectedBottomTOF; }

  inline Double_t GetCorrectedCubicTopTOF() { return CorrectedCubicTopTOF; }
  inline Double_t GetCorrectedCubicBottomTOF() { return CorrectedCubicBottomTOF; }

  inline Double_t GetCorrectedSoftTopTOF() { return CorrectedSoftTopTOF; }
  inline Double_t GetCorrectedSoftBottomTOF() { return CorrectedSoftBottomTOF; }

  inline Double_t GetCorrectedAvgSoftTOF() { return 0.5*(CorrectedSoftTopTOF+CorrectedSoftBottomTOF); }

  void SetBarId(Int_t v) { BarId=v; }
  void SetBarName(string v) { Name=v; }

  void Clear();
  void Clear(Option_t*) { Clear(); }

  void Finalize(); /* Calcuates convient quantities. */
  Bool_t operator==(const lendaBar & RHS); /*Probably not up to date. */

 public:
  ClassDef(lendaBar, 1);
};

/* Main container for Lenda events in the ROOT output. */

class lendaEvent : public TObject {
 public:
  vector <lendaBar> bars; 
  vector <lendaChannel> UnMappedChannels;
  vector <lendaChannel> refScint;
  
  Int_t NumOfChannelsInEvent; 
  Int_t N;
  Int_t numBars;
  Int_t NumUnMappedChannels;
  Int_t numObjectScintillators;

  Double_t ScintIntLow, ScintIntHigh, ScintIntAll;
  
  uint64_t etime;
  Double_t internaltime;
  
 public:
  lendaEvent();
  void Clear();
  void Clear(Option_t*){Clear();}
  
  /* Finalize should be called before writing to a tree.  Finalize calculates 
     convient quantities that will appear as branches and leaves in the tree.  
     Calls the finalize methods the Bars and ObjectScinitillators. */
  void Finalize(); 
  
  void PushABar(lendaBar aBar) { bars.push_back(aBar); numBars++; }
  void PushUnMappedChannel(lendaChannel c) {
    UnMappedChannels.push_back(c);
    NumUnMappedChannels++;
  }

  void SeteTime(uint64_t thetime) { etime = thetime; }
  void SetInternalTime(Double_t thetime) { internaltime = thetime;}
  Bool_t operator==(const lendaEvent& RHS);
  lendaBar * FindBar(int BarId);
  lendaBar * FindBar(string Name);
  lendaBar FindBarSafe(string Name);
  lendaBar FindBarSafe(int BarId);

  lendaChannel * FindChannel(string Name);
  lendaChannel * FindReferenceChannel(string Name);
  
 public:
  ClassDef(lendaEvent, 22);
};

class lendaFilter : public TObject {
 public:
  Double_t numOfBadFits;

 public:
  lendaFilter();
  ~lendaFilter();

  /* Software implementation of Pixie-16 fast trigger filter algorithm. */
  void FastFilter(std::vector <UShort_t> &trace,
		  std::vector <Double_t> &thisEventsFF, Double_t FL, Double_t FG);
  void FastFilterOp(std::vector <UShort_t> &trace,
		    std::vector <Double_t> &thisEventsFF, Double_t FL, Double_t FG);
  void FastFilterFull(std::vector <UShort_t> &trace,
		      std::vector <Double_t> &thisEventsFF,
		      Double_t FL, Double_t FG, Double_t decayTime);

  std::vector<Double_t> CFD (std::vector<Double_t>&, Double_t, Double_t);
  std::vector<Double_t> CFDOp (std::vector<Double_t>&, Double_t, Double_t);

  Double_t GetZeroCrossing(std::vector<Double_t>&, Int_t&, Double_t&);
  vector<Double_t> GetAllZeroCrossings(std::vector<Double_t>&);

  Double_t GetZeroCrossingImproved(std::vector<Double_t>&, Int_t&, Double_t&);
  Double_t GetZeroCrossingOp(std::vector<Double_t>&, Int_t&);

  Double_t GetZeroCubic(std::vector<Double_t>&, Int_t&);
  Double_t GetZeroFitCubic(std::vector<Double_t>&);   
  
  Double_t DoMatrixInversionAlgorithm(const std::vector<Double_t>& CFD, Int_t Spot);
  vector<Double_t> GetMatrixInversionAlgorithmCoeffients(const std::vector<Double_t>& CFD,
							 Int_t& Spot);

  Double_t fitTrace(std::vector<UShort_t>&, Double_t, Double_t);
  
  Double_t GetEnergyOld(std::vector<UShort_t>& trace);
  Double_t GetEnergy(std::vector<UShort_t>& trace, Int_t MaxSpot);

  vector<Double_t> GetEnergyHighRate(const std::vector<UShort_t>& trace,
				     std::vector<Int_t>& PeakSpots,
				     vector<Double_t>& theUnderShoots,
				     Double_t& MaxValueOut, Int_t& MaxIndexOut);

  vector<Int_t> GetPulseHeightHighRate(const std::vector<UShort_t>& trace,
				       const std::vector<Int_t>& PeakSpots);
  vector<Double_t> GetZeroCrossingHighRate(const std::vector<Double_t>& CFD,
					   const std::vector<Int_t>& PeakSpots);

  vector<Double_t> GetZeroCrossingCubicHighRate(const std::vector<Double_t>& CFD,
						const std::vector<Int_t>& PeakSpots);

  Double_t GetGate(std::vector<UShort_t>& trace, int start, int L);
  Int_t GetMaxPulseHeight(std::vector<UShort_t>& trace, Int_t& MaxSpot);
  Int_t GetMaxPulseHeight(std::vector<Double_t>& trace, Int_t& MaxSpot);

  double getFunc(TMatrixT<Double_t>, double);
  double getFunc(vector<double>&, double);
  Int_t GetStartForPulseShape(Int_t MaxSpot);

  Int_t CalculateCFD(vector<UShort_t> trace, vector<Double_t>&);
  vector<Double_t> GetNewFirmwareCFD(const vector<UShort_t>& trace, Int_t FL,
				     Int_t FG, Int_t d, Int_t w);
  Double_t GetNewFirmwareCFDWeight(Int_t CFDScaleFactor);
  
 public:
  ClassDef(lendaFilter, 1);
};

#endif
