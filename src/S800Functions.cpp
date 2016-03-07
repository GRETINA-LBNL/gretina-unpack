/*! \file S800Functions.cpp
    \brief Implementation of functions for utilities required, but not included 
    in S800 classes.
    
    Implementation of two functions related to S800 analysis, including 
    checking gate conditions, and writing distilled 'S800Physics' outputs
    to filtered files.  Implementation of functions declared in 
    S800Functions.h.

    Author: H. Crawford
*/


#include "S800Functions.h"

/****************************************************/

#ifdef WITH_S800

Int_t CheckS800PIDGates(TCutG* incomingBeam, TCutG* outgoingBeam) {

  Float_t inX = 0., inY = 0., outX = 0., outY = 0.;
  TString inXlabel = incomingBeam->GetVarX(); 
  TString inYlabel = incomingBeam->GetVarY();
  TString outXlabel = outgoingBeam->GetVarX();
  TString outYlabel = outgoingBeam->GetVarY();

  if (inXlabel == "tof.xfp") { inX = s800->tof.xfp; }
  else if (inXlabel == "tof.obj") { inX = s800->tof.obj; }
  else if (inXlabel == "tof.tac_xfp") { inX = s800->tof.tac_xfp; }
  else if (inXlabel == "tof.tac_obj") { inX = s800->tof.tac_obj; }
  else { cout << "Unexpected X coordinate on incoming beam PID gate = "
	      << inXlabel << endl; }
  if (inYlabel == "tof.xfp") { inY = s800->tof.xfp; }
  else if (inYlabel == "tof.obj") { inY = s800->tof.obj; }
  else if (inYlabel == "tof.tac_xfp") { inY = s800->tof.tac_xfp; }
  else if (inYlabel == "tof.tac_obj") { inY = s800->tof.tac_obj; }
  else { cout << "Unexpected Y coordinate on incoming beam PID gate = "
	      << inYlabel << endl; } 

  if (outXlabel == "tof.xfp") { outX = s800->tof.xfp; }
  else if (outXlabel == "tof.obj") { outX = s800->tof.obj; }
  else if (outXlabel == "tof.tac_xfp") { outX = s800->tof.tac_xfp; }
  else if (outXlabel == "tof.tac_obj") { outX = s800->tof.tac_obj; }
  else if (outXlabel == "tof.rf") { outX = s800->tof.rf; }
  else if (outXlabel == "tof.xfpe1") { outX = s800->tof.xfpe1; }
  else if (outXlabel == "tof.obje1") { outX = s800->tof.obje1; }
  else if (outXlabel == "tof.tac_xfpe1") { outX = s800->tof.tac_xfpe1; }
  else if (outXlabel == "tof.tac_obje1") { outX = s800->tof.tac_obje1; }
  else if (outXlabel == "tof.rfe1") { outX = s800->tof.rfe1; }
  else { cout << "Unexpected X coordinate on outgoing beam PID gate = "
	      << outXlabel << endl; }

  if (outYlabel == "fp.ic.de") { outY = s800->fp.ic.de; }
  else { cout << "Unexpected Y coordinate on outgoing beam PID gate = "
	      << outYlabel << endl; }

  // cout << "Got gates: " << endl;
  // cout << "  Incoming: " << inXlabel << " " << inYlabel << endl;
  // cout << "  Outgoing: " << outXlabel << " " << outYlabel << endl;

  // cout << "And values: " << endl;
  // cout << "  Incoming: " << inX << " " << inY << endl;
  // cout << "  Outgoing: " << outX << " " << outY << endl;

  Int_t inGates = (incomingBeam->IsInside(inX, inY) && 
		   outgoingBeam->IsInside(outX, outY));

  return inGates;
}

#endif
