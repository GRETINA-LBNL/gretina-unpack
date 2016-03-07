#include "S800Physics.h"

S800Physics::S800Physics() { ; }

S800Physics::~S800Physics() { ; }

S800Physics::Reset() {
  s800Ph.type = 0;
  s800Ph.crdc1_x = 0.0;
  s800Ph.crdc1_y = 0.0;
  s800Ph.crdc2_x = 0.0;
  s800Ph.crdc2_y = 0.0;
  s800Ph.ic_sum = 0.0;
  s800Ph.tof_xfp = 0.0;
  s800Ph.tof_obj = 0.0;
  s800Ph.rf = 0.0;
  s800Ph.trigger = 0;
  s800Ph.ic_de = 0.0;
  s800Ph.tof_xfpe1 = 0.0;
  s800Ph.tof_obje1 = 0.0;
  s800Ph.tof_rfe1 = 0.0;
  s800Ph.ata = 0.0;
  s800Ph.bta = 0.0;
  s800Ph.dta = 0.0;
  s800Ph.yta = 0.0;
}

S800Physics::getData(FILE *inf) {

}
