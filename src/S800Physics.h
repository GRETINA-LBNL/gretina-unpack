#ifndef __S800PHYSICS_H
#define __S800PHYSICS_H

#include "Riostream.h"
#include "Rtypes.h"

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TClass.h"

#include <vector>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include "colors.h"

struct s800Phys {
  int32_t type;
  float crdc1_x;
  float crdc1_y;
  float crdc2_x;
  float crdc2_y;
  float ic_sum;
  float tof_xfp;
  float tof_obj;
  float rf;
  int32_t trigger;
  float ic_de;
  float tof_xfpe1;
  float tof_obje1;
  float tof_rfe1;
  float ata;
  float bta;
  float dta;
  float yta;
};

class s800Physics : public TObject {
 public:
  s800Phys s800Ph;
  
 public:
  s800Physics() { ; }
  ~s800Physics() { ; }
  Reset();
  getData();
  
 private:
  ClassDef(s800Physics, 1);
};

#endif
