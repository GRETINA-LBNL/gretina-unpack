/*! \file INLCorrection.h
    \brief Parameter and function definitions for non-linearity correction.
    
    This provides the functionality to perform offline non-linearity 
    corrections to Mode3 GRETINA data, based on the waveforms.  The necessary
    raw INL data is stored within the class, as well as calculated 
    corrections, and baseline information.

    Author: H. Crawford

    Date: March 2013
*/

#ifndef __INLCORRECTION_H
#define __INLCORRECTION_H

#include <stdio.h>
#include <stdlib.h>

#include "SortingStructures.h"
#include "GRETINA.h"
#include "Defines.h"

/*! \class INLCorrection
    \brief Class containing the raw data and calculated corrections for 
    offline non-linearity correction of raw (Mode3) data.

    This class contains the raw non-linearity data, and the calculated
    parameters required for offline implementation of the non-linearity 
    correction to raw (Mode3) GRETINA data.  Data stored within the class
    includes the raw measured INL, and calculated corrections based on
    resting baseline values.  Functions include those required to read in
    raw INL data and digitizer SN information, and calculate corrections.

*/

class INLCorrection : public TObject {
 public:
  Int_t crystalSN[MAXCRYSTALS][4];
  Float_t inl[MAXCRYSTALS][40][8192];
  Float_t enl[MAXCRYSTALS][40][3][8192];

  Int_t lastCalcBase[MAXCHANNELS];
  Float_t restingBase[MAXCHANNELS];
  Float_t Q[MAXCHANNELS];
  Float_t sigmaBase[MAXCHANNELS];
  Int_t baseSamples[MAXCHANNELS];

 public:
  INLCorrection() { ; }
  ~INLCorrection() { ; }

  void Initialize(controlVariables* ctrl, GRETINAVariables* gVar);
  /*! \fn void Initialize(controlVariables* ctrl, GRETINAVariables* gVar)
      \brief Initialization of non-linearity correction variables.
      \param ctrl An instance of the controlVariables class.
      \param gVar An instance of the GRETINAVariables class.
      \return No return -- void.

      Initializes calculated corrections, etc. to 0, and if INL correction
      is being carried out, reads the digitizer map information using the 
      member ReadDigitizerMap(TString filename) function, and the raw INL 
      data using the ReadRawINL() function.
  */

  void calculateENL(Int_t xtal, Int_t chn, Int_t s1_width, Int_t s2_width,
		    Float_t s1, Float_t s2, Float_t base, Float_t tau);
  /*! \fn void calculateENL(Int_t xtal, Int_t chn, Int_t s1_width, Int_t s2_width, Float_t s1, Float_t s2, Float_t base, Float_t tau)
      \brief Using the measured raw INL, baseline and tau information, 
      calculates the required correction values as a function of ADC level 
      in the waveform.
      \param xtal Int_t crystal number, within the scheme of the Unpacking code
      \param chn Int_t channel number within the crystal
      \param s1_width # of channels of baseline to be summed over in 
      determining INL correction
      \param s2_width # of channels of waveform flattop to be summed over for 
      determining INL correction
      \param s1 Value of the sum over baseline samples for a given trace 
      -- not really used in the function
      \param s2 Float_t values of the sum over flattop samples for a given 
      trace -- not really used
      \param base Float_t baseline value to be used in calculating the
      correction values
      \param tau Float_t decay constant for the waveforms
      \return No return -- void

      Based on raw INL values, and the baseline and decay constant of pulses,
      calculates the required correction as a function of ADC samples summed
      over the widths given as arguments for the baseline and flattop of the
      pulses.  Fills the member variable enl arrays with the results.
  */

  void ReadDigitizerMap(TString filename);
  /*! \fn ReadDigitizerMap(TString filename)
      \brief Reads the digitizer map file to obtain the serial number for 
      each digitizer.
      \param filename TString of the digitizer map filename
      \return No return -- void

      Reads and interprets the digitizer map file for the array, and populates
      the member crystalSN variable array.
  */

  void ReadRawINL();
  /*! \fn ReadRawINL()
      \brief Reads and stores the raw INL correction data.
      \param None -- accesses class variables.
      \return No return -- void

      Reads the raw INL values from files according to digitizer SN.  
      Populates the inl class member arrays with the values indexed 
      according to crystal and channel numbers.
  */

 private:
  ClassDef(INLCorrection, 1);
};
  
#endif
