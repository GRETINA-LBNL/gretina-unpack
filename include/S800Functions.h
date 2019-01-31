/*! \file S800Functions.h
    \brief Function definitions for utilities required, but not included 
    in S800 classes.
    
    Provides two functions related to S800 analysis, including checking
    gate conditions, and writing distilled 'S800Physics' outputs to 
    filtered files.

    Author: H. Crawford
*/

#ifndef __S800FUNCTIONS_H
#define __S800FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>

#include "TCutG.h"

#include "Unpack.h"
#include "Globals.h"
#include "S800Parameters.h"

/****************************************************/

#ifdef WITH_S800

Int_t CheckS800PIDGates(TCutG* incomingBeam, TCutG* outgoingBeam);
/*! \fn Int_t CheckS800PIDGates(TCutG* incomingBeam, TCutG* outgoingBeam)
    \brief
    \param incomingBeam Pointer to ROOT TCutG for gating the incoming beam
    in an S800 analysis
    \param outgoingBeam Pointer to ROOT TCutG for gating on the outgoing beam
    in the S800
    \return Int_t value = 1 if gates are satisfied in the event, 0 if not

    Determines the x and y parameters of the passed gates, and then checks
    whether gate conditions are met for the event.  Returns 1 if both incoming
    and outgoing gate conditions are satisfied, and 0 if not.
*/

#endif

#endif
