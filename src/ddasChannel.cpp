/*! \file ddasChannel.cpp
    \ brief Function code for ddas channel event operations 
*/

#include "ddasChannel.h"

ClassImp(ddasChannel);

ddasChannel::ddasChannel() : TObject() {
  channelnum = 0;
  chanid = 0;
  slotid = 0;
  crateid = 0;
  channelheaderlength = 0;
  channellength = 0;
  finishcode = 0;
  overflowcode = 0;
  tracelength = 0;
  timelow = 0;
  timehigh = 0;
  timecfd = 0;
  cfdtrigsourcebit=0;
  energy = 0;
  
  time = 0;
  cfd = 0;
  
  energySums.reserve(4);
  qdcSums.reserve(8);
  trace.reserve(200);
  
  eclocklow=0;            ///< portion of ex timestamp
  eclockhigh=0;           ///< portion of ex timestamp
  eclockhigher=0;         ///< portion of ex timestamp
  eclock=0;			   ///< external time
}

ddasChannel::ddasChannel(const ddasChannel& obj)
  : TObject(obj), 
    channelnum(obj.channelnum),
    chanid(obj.chanid),
    slotid(obj.slotid),
    crateid(obj.crateid),
    channelheaderlength(obj.channelheaderlength),
    channellength(obj.channellength),
    finishcode(obj.finishcode),
    overflowcode(obj.overflowcode),
    tracelength(obj.tracelength),
    timelow(obj.timelow),
    timehigh(obj.timehigh),
    timecfd(obj.timecfd),
    cfdtrigsourcebit(obj.cfdtrigsourcebit),  
    energy(obj.energy),
    time(obj.time),
    cfd(obj.cfd),
    eclocklow(obj.eclocklow),
    eclockhigh(obj.eclockhigh),
    eclockhigher(obj.eclockhigher),
    eclock(obj.eclock),
    energySums(obj.energySums),
    qdcSums(obj.qdcSums),
    trace(obj.trace) { ; } 

ddasChannel& ddasChannel::operator=(const ddasChannel& obj) {
  if (this!=&obj) {
    channelnum = obj.channelnum;
    chanid = obj.chanid;
    slotid = obj.slotid;
    crateid = obj.crateid;
    channelheaderlength = obj.channelheaderlength;
    channellength = obj.channellength;
    finishcode = obj.finishcode;
    overflowcode = obj.overflowcode;
    tracelength = obj.tracelength;
    timelow = obj.timelow;
    timehigh = obj.timehigh;
    timecfd = obj.timecfd;
    cfdtrigsourcebit = obj.cfdtrigsourcebit;  //added CP
    energy = obj.energy;
    time = obj.time;
    cfd = obj.cfd;
    energySums = obj.energySums;
    qdcSums = obj.qdcSums;
    trace = obj.trace;
    eclocklow = obj.eclocklow;
    eclockhigh = obj.eclockhigh;
    eclockhigher = obj.eclockhigher;
    eclock = obj.eclock;
  }
  return *this;
}

void ddasChannel::Reset() {
  channelnum = 0;
  chanid = 0;
  slotid = 0;
  crateid = 0;
  channelheaderlength = 0;
  channellength = 0;
  finishcode = 0;
  tracelength = 0;
  cfdtrigsourcebit = 0;  //Added CP
  
  timelow = 0;
  timehigh = 0;
  timecfd = 0;
  energy = 0;
  
  time = 0;
  cfd = 0;
  
  energySums.clear();
  qdcSums.clear();
  trace.clear();
  
  eclocklow=0;            ///< portion of ex timestamp
  eclockhigh=0;           ///< portion of ex timestamp
  eclockhigher=0;         ///< portion of ex timestamp
  eclock=0;		  ///< external time
}

ddasChannel::~ddasChannel() { ; }

void ddasChannel::SetChannelID(UInt_t data) {
  chanid = ((Int_t) ((data & CHANNELIDMASK)));
}

void ddasChannel::SetSlotID(UInt_t data) {
  slotid = ((Int_t) ((data & SLOTIDMASK) >> 4));
}

void ddasChannel::SetCrateID(UInt_t data) {
  crateid = ((Int_t) ((data & CRATEIDMASK) >> 8));
}

void ddasChannel::SetChannelHeaderLength(UInt_t data) {
  channelheaderlength = ((Int_t) ((data & HEADERLENGTHMASK) >> 12 ));
}

void ddasChannel::SetChannelLength(UInt_t data) {
  channellength = ((Int_t) ((data & CHANNELLENGTHMASK) >> 17));
}

void ddasChannel::SetOverflowCode(UInt_t data) {
  overflowcode = ((Int_t) ((data & OVERFLOWMASK) >> 30));
}

void ddasChannel::SetFinishCode(UInt_t data) {
  finishcode = ((Int_t) ((data & FINISHCODEMASK) >> 31 ));
}

void ddasChannel::SetID(UInt_t data) {
  /* Broken function ... */
}

void ddasChannel::SetTimeLow(UInt_t data) {
  timelow = data;
}

void ddasChannel::SetTimeHigh(UInt_t data) {
  timehigh = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddasChannel::SetTimeCFD(UInt_t data) {
  timecfd = ((Int_t) ((data & BIT29to16MASK) >> 16));
}

void ddasChannel::SetCFDTriggerSourceBit(UInt_t data) {
  cfdtrigsourcebit = ((Int_t) ((data & BIT30MASK) >> 30 ));
}

void ddasChannel::SetTime() {
  time = timecfd/16384.0 + 2*(timelow + timehigh * 4294967296.0) - cfdtrigsourcebit;
}

void ddasChannel::SetEnergy(UInt_t data) {
  energy = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddasChannel::SetTraceLength(UInt_t data) {
  tracelength = ((Int_t) ((data & UPPER16BITMASK) >> 16));
}

// CWO: Adding EClock
void ddasChannel::SetEclockLow(UInt_t data) {
	eclocklow = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddasChannel::SetEclockHigh(UInt_t data) {
	eclockhigh = ((Int_t) ((data & UPPER16BITMASK)>>16));
}

void ddasChannel::SetEclockHigher(UInt_t data) {
	eclockhigher = ((Int_t) ((data & LOWER16BITMASK)));
}

void ddasChannel::SetEclock() {
  eclock = eclocklow + double(eclockhigh*65536.) + double(eclockhigher*4294967296.);
}

void ddasChannel::SetEnergySums(UInt_t data) {
  energySums.push_back(data);
}

void ddasChannel::SetQDCSums(UInt_t data) {
  qdcSums.push_back(data);
}

void ddasChannel::SetTraceValues(UInt_t data) {
  trace.push_back((data & LOWER16BITMASK));
  trace.push_back((data & UPPER16BITMASK)>>16);
}

void ddasChannel::UnpackChannelData(const uint32_t *data) {
  /* Put all variables in a known state before unpacking */
  Reset();

  int show_me=1;
  /* First four words of data identifiers, times, and energies. */
  /* data[0] is the number of half words in the event, skip it... */
  data++;

  /* A module-identifying 32-bit word is included in current ddas version (March 2016) */
  uint32_t module_identifier = *data;
  if (module_identifier == MODIDENTIFIER) {
    data++; // First word is Module-identifying word -- skip it. JP, March 2016
            // Otherwise, first word is the channel header
  } 

  /* Using the first word of DDAS information extract channel identifiers */
  SetChannelID(*data);
  SetSlotID(*data);
  SetCrateID(*data);
  SetChannelHeaderLength(*data);
  SetChannelLength(*data);
  SetOverflowCode(*data);
  SetFinishCode(*data++);
  /* Second word -- set most significant bits of timestamp. */
  SetTimeLow(*data++);

  /* Third word -- set least significant bits of timestamp. */
  SetTimeHigh(*data);

  SetTimeCFD(*data);
  SetCFDTriggerSourceBit(*data++);
  SetTime();

  /* Fourth word -- set energy. */
  SetEnergy(*data);
  SetTraceLength(*data++);

  /* Finished upacking the minimum set of data.  Determine if more data should be 
     unpacked based on a comparison between channel header length and channel length */

  /* If channel header length is 6 then the extra 2 words are external clock info. */
  if(channelheaderlength == 6) {
    SetEclockLow(*data);
    SetEclockHigh(*data++);
    SetEclockHigher(*data++);
    SetEclock();
  }
  /* If channel header length is 8 then the extra 4 words are energy sums and baselines. */
  if(channelheaderlength == 8) {
    for(int z=0; z<4; z++) { SetEnergySums(*data++); }
  }
  /* If channel header length ist 10 then extra 6 words are energy and external time. */
  if(channelheaderlength == 10) {
    for(int z=0; z<4; z++) { SetEnergySums(*data++); }
    SetEclockLow(*data);
    SetEclockHigh(*data++);
    SetEclockHigher(*data++);
    SetEclock();
  }
	
  /* If channel header length is 12 then the extra 8 words are QDC sums. */
  if(channelheaderlength == 12) {
      for(int z=0; z<8; z++) { SetQDCSums(*data++); }
  }
  
  /* If channel header length is 14 then the extra 10 words are QDC sums and external clock. */
  if(channelheaderlength == 14) {
    for(int z=0; z<8; z++) { SetQDCSums(*data++); }
    SetEclockLow(*data);
    SetEclockHigh(*data++);
    SetEclockHigher(*data++);
    SetEclock();
  }

  /* If channel header length is 16 then the extra 12 words are energy and QDC sums. */
  if(channelheaderlength == 16) { 
      for(int z=0; z<4; z++) { SetEnergySums(*data++); }
      for(int z=0; z<8; z++) { SetQDCSums(*data++); }
  }

  /* If channel header length is 18 then the extra 14 words are energy, qdc and external clock. */
  if(channelheaderlength == 18) {
    for(int z=0; z<4; z++) { SetEnergySums(*data++); }
    for(int z=0; z<8; z++) { SetQDCSums(*data++); }
    SetEclockLow(*data);
    SetEclockHigh(*data++);
    SetEclockHigher(*data++);
    SetEclock();
  }
    
  if(channelheaderlength != channellength) {
    if(channellength != (channelheaderlength + tracelength/2)) {
      cout << "Inconsistent lengths between channel length - " << channellength
           << " , header length - " << channelheaderlength 
           << " , and trace length - " << tracelength << endl;
      cout << "This is a fatal error - data is corrupted, please investigate " << endl;
      exit(1);
    }
    
    /* Starting traces.  If tracelength is non-zero, retrieve the trace. */
    if(tracelength !=0) {
      for(Int_t z = 0; z < tracelength/2; z++) {
	SetTraceValues(*data++);
      }
    } 
  } /* Finished unpacking extra data. */

}
