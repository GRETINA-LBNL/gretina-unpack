/*! \file CHICO.cpp
    \brief Parameter and functions for CHICO analysis.

    This file provides the analysis functions for CHICO analysis, based on 
    an analysis code CheckChico.cc from S. Zhu (ANL).

    Author: H.L. Crawford
    Date: April 29, 2014
*/

#include "DFMA.h"

void DFMARaw::Initialize() {

}

void DFMARaw::Reset() {
  Initialize();
}

DFMAFull::DFMAFull() {
 
}

void DFMAFull::Initialize() {
  raw.Initialize();
}

void DFMAFull::Reset() {
  ts = -1;
  raw.Reset();
}

