/*! \file ddasChannel.h
    \brief DdasChannel class definition header file 
*/

#ifndef __DDASCHANNEL_H
#define __DDASCHANNEL_H

#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>

#include "TObject.h"

/* Define bit masks to extract data from specific locations
   of pixie16 data */

#define CHANNELIDMASK             0xF  // Bits 0-3 inclusive
#define SLOTIDMASK               0xF0  // Bits 4-7 inclusive
#define CRATEIDMASK             0xF00  // Bits 8-11 inclusive
#define HEADERLENGTHMASK      0x1F000  // Bits 12-16 inclusive
#define CHANNELLENGTHMASK  0x3FFE0000  // Bits 17-29 inclusive
#define OVERFLOWMASK       0x40000000  // Bit 30 has overflow information (1 - overflow)
#define FINISHCODEMASK     0x80000000  // Bit 31 has pileup information (1 - pileup)
#define LOWER16BITMASK         0xFFFF  // Lower 16 bits
#define UPPER16BITMASK     0xFFFF0000  // Upper 16 bits
#define LOWER12BITMASK          0xFFF  // Lower 12 bits
#define BIT31MASK          0x80000000  // Bit 31 
#define BIT30to16MASK      0x7FFF0000  // Bits 30 through 16
#define BIT29to16MASK      0x3FFF0000
#define BIT30MASK          0X40000000
#define MODIDENTIFIER      0x0F0E00FA  // Module idenfier (250 MHz, 14-bit, revision F) 
/* Number of words added to pixie16 channel event header when energy sums
   and baselines are recorded */
#define SIZEOFESUMS                 4
/* Number of words added to pixie16 channel event header when QDC sums
   are recorded */
#define SIZEOFQDCSUMS               8

using namespace std;

/** Holds basic information for individual Pixie-16 channel. 
    Written by Sean Liddick  */

class ddasChannel : public TObject {

public:
  /* Channel events always have the following info. */
  Double_t time;              ///< Assembled time including cfd from module
  Double_t cfd;               ///< Cfd time only from module

  UInt_t energy;              ///< Energy of event from energy filter
  UInt_t timehigh;            ///< portion of timestamp
  UInt_t timelow;             ///< portion of timestamp
  UInt_t timecfd;             ///< portion of timestamp

  Int_t channelnum;           ///< obsolete
  Int_t finishcode;           ///< indicate whether event is 
  Int_t channellength;        
  Int_t channelheaderlength;
  Int_t overflowcode;
  Int_t chanid;              ///<Channel in module
  Int_t slotid;              ///<Slot in crate
  Int_t crateid;             ///<Crate ID
  Int_t id;
  UInt_t cfdtrigsourcebit;     ///< cfd trigger source bit, 0- even, 1 - odd

  Int_t tracelength;
  /* External clock additional information */
  UInt_t eclocklow;            ///< portion of ex timestamp
  UInt_t eclockhigh;           ///< portion of ex timestamp
  UInt_t eclockhigher;         ///< portion of ex timestamp
  uint64_t eclock;	       ///< external time

  /* A channel may have extra information... */
  vector<UInt_t> energySums;
  vector<UShort_t> qdcSums;
  
  /* A waveform (trace) may be stored too. */
  vector<UShort_t> trace;

public:
  ddasChannel();
  ddasChannel(const ddasChannel& obj);
  ddasChannel& operator=(const ddasChannel& obj);
  ~ddasChannel();
  void Reset();

  UInt_t GetEnergy() { return energy; }
  UInt_t GetTimeHigh() { return timehigh; }
  UInt_t GetTimeLow() { return timelow; }
  UInt_t GetCFDTime() { return timecfd; }
  Double_t GetTime() { return time; }
  Double_t GetCFD() { return cfd; }
  UInt_t GetEnergySums(Int_t i) { return energySums[i]; }
  Int_t GetChannelNum() { return channelnum; }
  Int_t GetFinishCode() { return finishcode; }
  Int_t GetChannelLength() { return channellength; }
  Int_t GetChannelLengthHeader() { return channelheaderlength; }
  Int_t GetOverflowCode() { return overflowcode; }
  Int_t GetSlotID() { return slotid; }
  Int_t GetCrateID() { return crateid; }
  Int_t GetChannelID() { return chanid; }
  Int_t GetID() { return id;}
  UInt_t GetCFDTriggerSourceBit() { return cfdtrigsourcebit; }
  void SetEclockLow(UInt_t data);
  void SetEclockHigh(UInt_t data);
  void SetEclockHigher(UInt_t data);
  void SetEclock();
  
  vector<UShort_t> GetTrace() { return trace; }

  void SetChannelID(UInt_t data);
  void SetSlotID(UInt_t data);
  void SetCrateID(UInt_t data);
  void SetChannelHeaderLength(UInt_t data);
  void SetChannelLength(UInt_t data);
  void SetOverflowCode(UInt_t data);
  void SetFinishCode(UInt_t data);
  void SetID(UInt_t data);
  void SetCFDTriggerSourceBit(UInt_t data);

  void SetTimeLow(UInt_t data);
  void SetTimeHigh(UInt_t data);
  void SetTimeCFD(UInt_t data);
  void SetTime();
  void SetEnergy(UInt_t data);
  void SetTraceLength(UInt_t data);

  void SetEnergySums(UInt_t data);
  void SetQDCSums(UInt_t data);
  void SetTraceValues(UInt_t data);

  void UnpackChannelData(const uint32_t *data);

  ClassDef(ddasChannel, 1);
};

#endif
