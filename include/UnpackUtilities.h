/*! \file UnpackUtilities.h
    \brief Function definitions meant to remove some stuff from the main
           body of Unpack and simplify the main loop
    
    This provides a few functions which can stand along, and allow an overall
    simplification of the main unpacking loop contained in Unpack.cpp.  
    Utilities include opening input files, adding data to an event, processing
    finished events and resetting event structures.

    Author: H. Crawford

    Date: March 2013
*/

#ifndef __UNPACKUTILITIES_H
#define __UNPACKUTILITIES_H

#include <stdio.h>
#include <stdlib.h>

#include "TVector3.h"

#include "Globals.h"

#include "Unpack.h"
#include "GRETINA.h"
#include "Track.h"

/* These are utilites to simplify the main unpacking code. */

Int_t OpenInputFile(FILE** inf, controlVariables* ctrl, TString runNumber); 
/*! \fn Int_t OpenInputFile(FILE** inf, controlVariables* ctrl, TString runNumber)
    \brief Determines the full path/filename for the input file, and opens it.
    \param inf FILE* pointer for the input file to be opened.
    \param ctrl An instance of the controlVariables class dictating file types, etc.
    \param runNumber TString containing the run number to build the filename from.
    \return Int_t -- returns 0 if successful in opening file.

    According to the flags provided via the controlVariables* object, builds 
    up the filename including the directory path for the input file.  Attempts 
    to open the file, iterating through likely directories.  Returns an 
    integer value to indicate success or failure.
*/

int ProcessEvent(Float_t currTS, controlVariables* ctrl, 
		  counterVariables* cnt);
/*! \fn void ProcessEvent(Float_t currTS, controlVariables* ctrl, counterVariables* cnt, GRETINAVariables* gVar, SuperPulse* sp, Histos* histos)
    \brief Analyzes time-grouped physics events, making calls to do segment 
           energy analysis, tracking etc.
    \param currTS Float_t value of the last timestamp in the event.
    \param ctrl An instance of the controlVariables class, providing the flags 
           to identify cross-talk or superpulse analysis
    \param cnt An instance of the counterVariables class, which is required to 
           pass to analysis functions for tallying event characteristics
    \param histos An instance of the Histos class, if histogram filling is 
           active in the analysis.
    \return Returns int = 0 if no problems, -1 if weird crystal number observed

    This function accesses the globally stored event data, and based on control 
    flags determines the required analysis.  From this function, calls to 
    superpulse analysis, cross-talk analysis, segment analysis and tracking 
    for GRETINA data are possible.  Also in this function, following analysis 
    of the physics event, ROOT histograms and/or trees are filled.
*/

void ResetEvent(controlVariables* ctrl, counterVariables* cnt);
/*! \fn void ResetEvent(controlVariables* ctrl)
    \brief Resets event structures, called after an event has been analyzed and written out.
    \param ctrl An instance of the controlVariables class, to indicate if waveform information needs to be reset.
    \return No return -- void.

    This function resets the global event-built data structures following 
    analysis of time-grouped event.  It clears the globally-defined vectors 
    for all GRETINA data types, resets waveform vectors if they are being 
    used, and clears any auxiliary detector data structures.
*/

#endif
