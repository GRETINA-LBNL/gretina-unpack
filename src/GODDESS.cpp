#include "GODDESS.h"

void orrubaDet::Clear() {
  enRawPmap.clear();  enRawNmap.clear();
  enCalPmap.clear();  enCalNmap.clear();
  timePmap.clear();   timeNmap.clear();
}

void orrubaDet::SetNumContacts(Int_t pType, Int_t nType) {
  numPtype = pType;  numNtype = nType;
  parECalP.resize(pType);
  parECalN.resize(nType);
}

Bool_t orrubaDet::ValidContact(UInt_t contact, Bool_t nType) {
  size_t size = numPtype;
  if (nType) { size = numNtype; }
  if (contact >= size) {
    printf("ERROR: Contact specified, %u, is invalid!\n", contact);
    return kFALSE; 
  }
  return kTRUE;
}

Float_t orrubaDet::GetCalEnergy(Int_t contact, Bool_t nType) {
  if (!ValidContact(contact, nType)) { return 0; }
  valueMap *eCal;
  if (nType) { eCal = &eCalNmap; } else { eCal = &eCalPmap; }
  
  auto itr = eCal->find(contact);
  if (itr != eCal->end()) { return itr->second; }
  return 0;
}

unsigned long

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
  printf("-----------------------------");
  printf(" Raw GODDESS event: \n");
  printf("   TS = %ld\n", rawGoddess->timestamp);
  printf("   Channel\tValue\n");
  for (UInt_t i=0; i<rawGoddess->channels.size(); i++) {
    printf("    %d\t%d\n", rawGoddess->channels[i], rawGoddess->values[i]);
  }
  printf("-----------------------------");
}
