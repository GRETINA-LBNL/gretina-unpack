/*! \file SuperPulse.h
    \brief Parameter and function defintitions for superpulse analysis.
    
    This provides the variables, parameters and functions required for
    superpulse analysis, based on the original C codes from David Radford 
    and Mario Cromaz.  Superpulse analysis as implemented produces .spn
    files for viewing in Radware, and use in all fitting codes, etc. for 
    basis production, as the original C codes did.

    Author: H. Crawford

    Date: January 2013
*/

#ifndef __SUPERPULSE_H
#define __SUPERPULSE_H

using namespace std;

#include <stdio.h>
#include <stdlib.h>

#include "SortingStructures.h"
#include "GRETINA.h"
#include "Defines.h"

#define DEBUG 0

/*! \class SuperPulse
    \brief Class containing the parameters required to peform a superpulse 
    analysis.

    This class contains the calibration parameters (and functions required
    to get them) for the superpulse analysis, as well as crystal-by-crystal
    values for important variables including segment multiplicity, hit segment
    energies, etc.  The averaged traces are also stored within this class.  
    Member functions include reading in calibration information, making 
    superpulses in terms of averaging traces, and writing out .spn files.

*/

class SuperPulse : public TObject {
 public:
  Int_t CFD_INT_LEN;
  Int_t CFD_DELAY;
  Int_t CFD_FRACTION;
  Float_t TR_SCALE;

  /* Calibration stuff */
  Float_t ehiGeOffset[MAXCHANNELS];
  Float_t ehiGeGain[MAXCHANNELS];
  Float_t trGain[MAXCHANNELS];
  Int_t map[MAXCRYSTALS][40];
  Float_t delay1[MAXCRYSTALS][40];

  Float_t lowE, highE;
  Int_t trLength;
  Int_t mult[MAXCRYSTALS];
  Int_t crystalBuild[MAXCRYSTALS];
  Float_t ccE[MAXCRYSTALS];
  Float_t segE[MAXCRYSTALS];
  Int_t netSeg[MAXCRYSTALS];
  Int_t segEventIndex[MAXCRYSTALS][40];
  Int_t segsNet[MAXCRYSTALS][40];
  Int_t segs[MAXCRYSTALS][40];
  Int_t data4net[MAXCRYSTALS][36];

  Float_t gain[MAXCRYSTALS][40];
  Int_t waves[MAXCRYSTALS][40][2048];
  Float_t averageTrace[MAXCRYSTALS][40][4096];
  Int_t averageTraceINT[MAXCRYSTALS][40][4096];

 public:
  SuperPulse() { ; }
  ~SuperPulse() { ; }

  void Initialize(controlVariables* ctrl, GRETINAVariables* gVar);
  /*! \fn void Initialize(controlVariables* ctrl, GRETINAVariables* gVar)
      \brief Initialization for superpulse analysis variables. 
      \param ctrl An instance of the controlVariables class.
      \param gVar An instance of the GRETINAVariables class.
      \return No return -- void.

      Initializes all counters, etc. to 0, and gains/offsets to 1.0/0.0.
      If superpulse analysis is called for in the unpacking, then calls
      the ReadDetMaps and ReadParams functions.     
  */

  Int_t ReadDetMaps(char *fn, GRETINAVariables* gVar);
  /*! \fn Int_t ReadDetMaps(char *fn, GRETINAVariables* gVar)
      \brief Reads in calibrations from detector map files. 
      \param fn String for the directory path to the detector map files.
      \param gVar An instance of the GRETINAVariables class.
      \return Int_t -- indicates success if 0. 
      
      Loops over all crystals/CC detector maps within the given directory, 
      and fills the gain and offset parameters in the SuperPulse class.
  */

  Int_t ReadParams(TString filename, const char *label, 
		   Float_t x[][40], Int_t len, GRETINAVariables* gVar);
  /*! \fn Int_t ReadParams(TString filename, const char *label, 
          Float_t x[][40], Int_t len, GRETINAVariables* gVar)
      \brief Reads data from the cross-talk GRETINA files.
      \param filename TString for the file containing the list of cross-talk files.
      \param label String for the parameter in the cross-talk file to be obtained.
      \param x[][40] Array of floats, to be filled with the parameters from the file.
      \param len Int_t number of parameters to be filled, per crystal (usually 40 channels worth).
      \param gVar An instance of the GRETINAVariables class. 
      \return Int_t -- indicates success if 0.

      Reads and stores the specified parameter from the cross-talk files for
      the GRETINA crystals.  In superpulse analysis, the parameter of interest
      is the delay1 value for time alignment considerations. 
  */

  void MakeSuperPulses();
  /*! \fn void MakeSuperPulses()
      \brief Adds the waveforms for an event to the average, if all requirements are met.
      \param None -- accesses global variables for data, all others are class variables.
      \return No return -- void.

      Checks requirements of segment multiplicity and energy cuts, calls AlignCFD routine 
      to ensure time alignment, and if all is good, averages event trace data with previous 
      traces.
  */

  Int_t AlignCFD(Int_t crystalNum);
  /*! \fn Int_t AlignCFD(Int_t crystalNum)
      \brief Shifts waveforms as appropriate based on CFD timing and delays from cross-talk file.
      \param crystalNum Int_t crystal number, within the numbering scheme of the unpacking code
      \return Int_t -- indicates success if 0.

      For a given crystal, determines the CFD time of the net segment and CC (calls cfdTime
      function), and calculates the shift required for time alignment.  Then literally 
      shifts all of the waveforms.
  */

  Float_t cfdTime(Int_t crystalNum, Int_t segNum);
  /*! \fn Float_t cfdTime(Int_t crystalNum, Int_t segNum)
      \brief Calculates the sample number of the CFD crossing for a pulse.
      \param crystalNum Int_t crystal number, within the numbering scheme of the unpacking code.
      \param segNum, Int_t segment number 
      \return Float_t values of the CFD crossing time (in samples).

      Calculates a Ge-style CFD to get the signal timing.  This algorithm comes from
      David Radford initially. 
  */

  void FinishSuperPulses();
  /*! \fn void FinishSuperPulses()
      \brief Finishes superpulse analysis -- normalizes traces, etc. 
      \param None -- all data required are stored in class variables.
      \return No return -- void.

      Calculates trace gains, and renormalizes traces accordingly.
  */

  void WriteSuperPulses();
  /*! \fn void WriteSuperPulses()
      \brief Writes superpulses to .spn files. 
      \param None -- all data required are stored in class variables.
      \return No return -- void.

      Writes superpulse data to .spn files with the naming convention of
      "SPCrystal<CrystalNum>_<Energy>.spn".
  */
  
 private:
  ClassDef(SuperPulse, 1);
};
  
#endif
