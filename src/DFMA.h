/*! \file DFMA.h
    \brief Parameter and function definitions for CHICO analysis.

    This file provides the data structures and function prototypes for CHICO
    analysis, based on an analysis code CheckChico.cc from S. Zhu (ANL).

    Author: H.L. Crawford
    Date: April 29, 2014
*/

#ifndef __DFMA_H
#define __DFMA_H

#include "colors.h"

#include "TObject.h"
#include "TMath.h"
#include "TString.h"
#include <cstring>
#include <vector>
#include <iostream>
#include <RConfigure.h>

class DFMAFull;

class DFMARaw : public TObject {
 private:
  DFMAFull* m_top;
  
 public:

 public:
  void Initialize();
  void Reset();
  
  ClassDef(DFMARaw, 1);
};

class DFMAFull : public TObject {

 public:
  DFMARaw raw;
 
  long long int ts;

 public:
  DFMAFull();
  
 public:
  void Initialize();
  void Reset();
  
 private:
  
  ClassDef(DFMAFull, 1);
};

#endif
