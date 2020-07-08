#ifndef __GRETINA_H
#define __GRETINA_H

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

#include "Track.h"
#include "SortingStructures.h"
#include "Defines.h"
#include "colors.h"

using std::vector;
using std::cout;  using std::endl;
using std::swap;

class GRETINA;

/******** Rotation matrix for crystal to lab transformation ********/
class rotationMatrix : public TObject {
 
 public:
  Float_t crmat[MAXDETPOS+1][MAXCRYSTALNUM+1][4][4];
  
 public:
  rotationMatrix() { ; }
  rotationMatrix(TString file) { ReadMatrix(file); }
  ~rotationMatrix() { ; }
  Int_t ReadMatrix(TString file);
  TVector3 crys2Lab(Int_t crystalID, TVector3 xyz);
  
  ClassDef(rotationMatrix, 1);
};

/*--------------------------------------------------------*/
/*---          GRETINA ANALYSIS VARIABLE CLASS         ---*/
/*--------------------------------------------------------*/

class GRETINAVariables : public TObject {

 public:
  /* Geometry/set-up stuff */
  Int_t hole[MAXQUADS];
  Int_t electronicsOrder[MAXQUADS];

  Double_t segCenter[2][3][36]; /* Segment centers */

  /* Doppler correction stuff */
  Float_t beta;
  TVector3 targetXYZ;
  Float_t radiusCor[MAXCRYSTALS];
  Float_t ataCor, btaCor;

  /* Calibration */
  Float_t ehiGeOffset[MAXCHANNELS];
  Float_t ehiGeGain[MAXCHANNELS];

  /* Non-linearity look-up table correction */
  Float_t dnlLU[MAXCHANNELS][5000];

  /* Old xtalk correction parameters */
  Float_t dinoFactor[MAXCRYSTALS][40][40];

  /* Variables related to waveform analysis */
  Int_t minCrossTime, maxCrossTime;

  /* Segment mapping... */
  Int_t Q1Special[4][36];
  Int_t Q2Special[4][36];
  Int_t QSegmentOrder[4][36];
  
 public:
  GRETINAVariables();
  ~GRETINAVariables() { ; }

 public:
  void Initialize();
  void InitializeGRETINAVariables(TString inputFilename);
  Double_t GetNextValue(FILE *file);
  void Reset();

  void ReadGeCalFile(TString filename);
  void ReadSegmentCenters(TString filename);

  void PrintSegCenters();
  


 private:
  
  ClassDef(GRETINAVariables, 1);
};

/******** Structures to store raw data ********/

struct globalHeader {
  Int_t type;
  Int_t length;
  long long timestamp;
};

/*--------------------------------------------------------*/
/*---       GRETINA MODE 3 - RAW DATA STRUCTURES       ---*/
/*--------------------------------------------------------*/

struct mode3DataPacket {
  unsigned short aahdr[2];
  unsigned short hdr[14];
  /* with traces... */
  unsigned short waveform[MAX_TRACE_LENGTH];
};

class g3Waveform : public TObject {

 public:
  vector<Short_t> raw;
  vector<Float_t> pz;
  
 public:
  void Clear();
  Int_t LEDLevel(Int_t index, Int_t thresh, Int_t filterK);
  Int_t LED(Int_t index, Int_t thresh, Int_t filterK);
  Float_t CFD(Int_t startSample);
  Float_t BL(Int_t start, Int_t end);
  Int_t Look4Pileup();
  void calcPZ(Float_t base, Float_t tau);
  Float_t simpleE();
  Float_t riseTime(Float_t fLow, Float_t fHigh);
  Double_t computeChiSquare(Int_t startIndex, Int_t nPoints,
			    Double_t tau, Double_t nStart, Double_t baseline);
  Double_t LSFitExpo(Int_t startIndex, Int_t nPoints,
		     Double_t nStart, Double_t tau,
		     Double_t &nStartFit, Double_t &baselineFit,
		     Double_t &nStartFitError,
		     Double_t &baselineFitError);
  Double_t LSFitLinear(Int_t startIndex, Int_t nPoints,
		       Double_t &slopeFit, Double_t &offsetFit);
  Float_t baseline2Flatten(Float_t tau, Int_t startIndexBase,
			   Int_t lengthBase, Int_t startIndexTop,
			   Int_t lengthTop, Float_t smallStep);
  
  ClassDef(g3Waveform, 1);
};

class g3ChannelEvent : public TObject {

 public:
  UShort_t hdr0, hdr1, hdr7;
  Int_t ID;
  Int_t segNum;
  Float_t eRaw, eCal;
  /* These are for event history, available in firmware 2.00_006f and beyond. */
  Float_t eCalPO, prevE1, prevE2;
  UShort_t deltaT1, deltaT2;
  UShort_t PZrollover;

  Float_t eCalc;
  Float_t calcTime;
  long long int timestamp;
  long long int CFDtimestamp;
  Float_t baseline;
  Float_t riseTime;
  g3Waveform wf;

 public:
  Int_t module() { return (hdr1 >> 4); }
  Int_t boardID() { return (hdr0 >> 11); }
  Int_t chanID() { return (hdr1 & 0xf); }
  Int_t tracelength() { return (((hdr0 & 0x7ff)*2) - 14); }
  Int_t chnNum() { return (ID%40); }
  Bool_t sign() { return (hdr7 & 0x0100); }
  Int_t pileUp() { return ((hdr7 & 0x8000) >> 15); }    

 public:

  ClassDef(g3ChannelEvent, 3);
};

class g3CrystalEvent : public TObject {
 
 public:
  Int_t crystalNum, quadNum, module;
  Int_t runNumber;
  Float_t dopplerSeg, dopplerCrystal;
  Int_t traceLength;
  Float_t segSum, segSumTimed, segSumTimed2;
  vector<g3ChannelEvent> chn;

 public:
  void OrderChannels();
  Float_t cc1();
  Float_t cc2();
  Float_t cc3();
  Float_t cc4();
  Float_t cc1Raw();
  Float_t cc2Raw();
  Float_t cc3Raw();
  Float_t cc4Raw();
  Float_t cc1Calc();
  Float_t cc2Calc();
  Float_t cc3Calc();
  Float_t cc4Calc();
  Float_t segSumRaw();
  Float_t segSumRawThresh(Float_t thresh);
  Float_t segSumRawTimed(Float_t minT, Float_t maxT);
  Float_t segSumRawThreshTimed(Float_t thresh, Float_t minT, Float_t maxT);
  Int_t segsHit(Float_t thresh);
  Float_t segSumCalc();
  Float_t segSumCalcThresh(Float_t thresh);
  Int_t segsHitCalc(Float_t thresh);
  Float_t maxSegE();
  Float_t maxSegECalc();
  Int_t maxSegNum();
  Int_t secondSegNum();
  Int_t maxSegNumCalc();
  Int_t deepRingHit(Float_t thresh);
  long long int LEDLow();
  long long int LEDHigh();
  UInt_t crystalBuild();
  long long int LEDRange();
  void mode3ToFile(FILE* outFile);

 private:
  ClassDef(g3CrystalEvent, 2);
};

class g3OUT : public TObject {

 public:
  vector <g3CrystalEvent> xtals;

 public:
  g3OUT() { ; }
  ~g3OUT() { ; }
  void Reset();
  UInt_t crystalMult();
  Float_t calorimeterE();

 private:
  ClassDef(g3OUT, 1);
};

struct mode3HistoryPacket {
  unsigned short aahdr[2];
  unsigned short hdr[8];
  unsigned short data[MAX_TRACE_LENGTH];
};

struct historyEvent{
  Float_t energy;
  long long int TS;
  Int_t BLpreSum;
  // 2016-07-23 CMC added module to gH to differentiate between digitizers
  Int_t module;
};

class g3HistoryEvent : public TObject {
  
 public:
  std::vector<historyEvent> past;
  
 public:
  g3HistoryEvent() { ; }
  ~g3HistoryEvent() { ; }
  UInt_t lengthOfPast() { return past.size(); }

 private: 
  ClassDef(g3HistoryEvent, 1);
};

class Bank88 : public TObject {
 
 public:
  vector<g3ChannelEvent> chn;
  long long int timestamp;
  long long int wfCFD;

 public:
  Bank88() { ; }
  ~Bank88() { ; }
 
  ClassDef(Bank88, 2);
};

/*--------------------------------------------------------*/
/*---    GRETINA MODE 2 - DECOMPOSED DATA STRUCTURE    ---*/
/*--------------------------------------------------------*/

/* Single interaction point structure */
struct ip { 
  Float_t x, y, z, e; /* Here, e refers to the fraction. */
  Int_t seg; /* Segment hit. */
  Float_t seg_energy; /* Energy of hit segment. */
};

struct ipOLD {
  Float_t x, y, z, e;
};

struct mode2Old {
  Int_t num; 
  Int_t crystal_id;
  Float_t tot_e;
  Float_t t0;
  Float_t chisq;
  Float_t norm_chisq;
  long long int timestamp;
  ipOLD intpts[MAX_INTPTS];
};

struct mode2ABCD1234 {
  Int_t type;
  Int_t crystal_id;
  Int_t num;
  Float_t tot_e;
  long long int timestamp;
  long long trig_time; 
  Float_t t0; 
  Float_t cfd;
  Float_t chisq;
  Float_t norm_chisq;
  Float_t baseline;
  Int_t pad; /* To ensure 8-byte alignment of struct. */
  ip intpts[MAX_INTPTS];
};

struct mode2ABCD5678 {
  Int_t type;
  Int_t crystal_id;
  Int_t num;
  Float_t tot_e;
  Int_t core_e[4];
  long long int timestamp;
  long long trig_time;
  Float_t t0;
  Float_t cfd;
  Float_t chisq;
  Float_t norm_chisq;
  Float_t baseline;
  Float_t prestep;
  Float_t poststep;
  Int_t pad;
  ip intpts[MAX_INTPTS];
};

struct mode2ABCD6789 {
  Int_t type;
  Int_t crystal_id;
  Int_t num;
  Float_t tot_e;
  Int_t core_e[4];
  long long int timestamp;
  /* Used for history-based energy correction at rate */
  Float_t totE_fixedPickOff_prior1;
  Float_t totE_fixedPickOff_prior2;
  Float_t t0;
  /* Used for history-based energy correction at rate */
  UShort_t deltaT_prior1;
  UShort_t deltaT_prior2;
  Float_t chisq;
  Float_t norm_chisq;
  /* Used for history-based energy correction at rate */
  Float_t totE_fixedPickOff_current;
  Float_t prestep;
  Float_t poststep;
  Int_t pad;
  ip intpts[MAX_INTPTS];

};

/* \class g2IntPt
   \brief Class containing the information relevant to a single gamma-ray interaction point resulting from
   signal decomposition (mode2 data).

   This class defines a single interaction point in a given GRETINA crystal.  It includes all relevant 
   parameters --> energy, position and segment number.
*/

class g2IntPt : public TObject {
  
 public:
  Int_t segNum;
  TVector3 xyz;
  TVector3 xyzLab;
  TVector3 xyzLabSeg;
  TVector3 xyzLabCrys;
  Float_t e, segE;
  
 public:
  g2IntPt() { ; }
  ~g2IntPt() { ; }

 public:
  ClassDef(g2IntPt, 1);
};

/* \class g2CrystalEvent
   \brief Class containing the information relevant for mode2 (decomposed data) crystal-level event.

   This class defines a single crystal event for mode2 (decomposed) data from GRETINA.  The
   class includes crystal-level information (such as the fit information from signal decomp and
   the crystal ID etc.), as well as a vector of interaction point information 
   (vector of g2IntPt class objects).
*/

class g2CrystalEvent : public TObject {

 public:
  Int_t crystalID;
  Int_t crystalNum;
  Int_t quadNum;
  long long int timestamp;
  Float_t t0;
  Float_t chiSq;
  Float_t normChiSq;
  Float_t bl;
  Int_t error;
  Float_t cc;
  Float_t ccCurrent, ccPrior1, ccPrior2;
  UShort_t deltaT1, deltaT2;
  Float_t cc1, cc2, cc3, cc4;
  Float_t segSum;
  Float_t doppler;
  Float_t dopplerT; /*20180517 CMC added for CHICO target-like gammas*/
  Float_t dopplerSeg;
  Float_t dopplerCrystal;
  vector<g2IntPt> intpts;
  vector<Short_t> waveAll;
  
 public:
  g2CrystalEvent() { ; }
  ~g2CrystalEvent() { ; }
  void Reset();
  Int_t numIntPts();
  Int_t numSegHit();
  Int_t segsHitAboveThresh(Float_t thresh);
  Float_t averageIntPtsPerSeg();
  Float_t segSumRaw();
  Float_t segSumRawThresh(Float_t thresh);
  TVector3 maxIntPtXYZ();
  Double_t maxIntPtX() { return maxIntPtXYZ().X(); }
  Double_t maxIntPtY() { return maxIntPtXYZ().Y(); }
  Double_t maxIntPtZ() { return maxIntPtXYZ().Z(); }
  Double_t maxIntPtR() { return maxIntPtXYZ().XYvector().Mod(); }
  TVector3 maxIntPtXYZLab();
  TVector3 maxIntPtXYZLabSeg();
  TVector3 maxIntPtXYZLabCrystal();
  Float_t maxIntPtE();
  Int_t maxIntPtSegNum();
  Float_t maxIntPtSegE();
  Float_t gTheta();
  Float_t gPhi();
  Float_t CosDop();
  Float_t cDoppler(Float_t beta, Float_t x, Float_t y, Float_t z);
  Float_t cDopplerSeg(Float_t beta);
  Float_t cDopplerCrystal(Float_t beta);

 public:
  ClassDef(g2CrystalEvent, 1);
};

class g2OUT : public TObject {

 public:
  Int_t runNumber;
  vector<g2CrystalEvent> xtals;
  
 public:
  g2OUT() { ; }
  ~g2OUT() { ; }
  void Reset();
  Float_t calorimeterDoppler(Float_t beta);
  Float_t calorimeterE(Float_t thresh);
  Float_t calorimeterE();
  UInt_t crystalMult();

 public:
  ClassDef(g2OUT, 1);
};

/*--------------------------------------------------------*/
/*---     GRETINA MODE 1 - TRACKED DATA STRUCTURES     ---*/
/*--------------------------------------------------------*/

struct trackedGamma {
  Float_t esum; /* gamma ray energy */
  Int_t ndet; /* number of interactions */
  Float_t fom; /* figure of merit */
  Int_t tracked; /* tracked == 1 if successful */
  long long int timestamp; /* timestamp of first interaction point */
  Float_t x0, y0, z0, e0; /* first interaction point */
  Float_t x1, y1, z1, e1; /* second interaction point */
  //Short_t fCryHit; /* First crystal hit ID */
};

class GTrackEvent : public TObject {
 public:
  vector<trackedGamma> obj; // [mult]

  GTrackEvent() { ; }
  ~GTrackEvent() { ; }
  
  Int_t GetVecSize() { return obj.size(); }

  ClassDef(GTrackEvent, 1);
};

class g1GammaEvent : public TObject {

 public:
  Int_t tracked;
  Float_t FOM;
  Int_t nIntPts;
  Float_t cc;
  Float_t doppler;
  TVector3 xyzLab1, xyzLab2;
  Float_t e1, e2;
  long long int timestamp;
  Float_t t0;

 public:
  g1GammaEvent() { ; }
  ~g1GammaEvent() { ; }
  Float_t gTheta();
  Float_t gPhi();
  /* Polarization functions */
  TVector3 normalVector1();
  TVector3 normalVector2();
  Float_t phiCompton();
  Float_t intPtDistance();
  Float_t scatterAngle();
  Float_t cDoppler(Float_t beta, Float_t x, Float_t y, Float_t z);
  
  ClassDef(g1GammaEvent, 1);

};

class g1OUT : public TObject {
 
 public:
  vector<g1GammaEvent> gammas;

 public:
  g1OUT() { ; }
  ~g1OUT() { ; }
  void Reset();
  Float_t calorimeterE();
  UInt_t gammaMult();

  ClassDef(g1OUT, 1);
};

/*--------------------------------------------------------*/
/*---        GRETINA SIMULATED DATA STRUCTURES         ---*/
/*--------------------------------------------------------*/

struct g4Sim_emittedGamma {
  Float_t e;
  Float_t x, y, z;
  Float_t phi, theta;
  Float_t beta;
};

struct g4Sim_abcd1234 {
  Int_t type;
  Int_t num;
  Int_t full;
  g4Sim_emittedGamma gammas[MAX_SIM_GAMMAS];
};

class g4SimOUT : public TObject {
 public:
  vector<g4Sim_abcd1234> sim;
  
 public:
  g4SimOUT() { ; }
  ~g4SimOUT() { ; }
  UInt_t GetSimMult() { return sim.size(); }
  void Reset() { sim.clear(); }
  
  ClassDef(g4SimOUT, 2);
};  

/* \class gHistos
   \brief Class containing GRETINA-related histograms for filling during analysis.

   This class holds the TH* histogram objects that may (or may not) be filled
   during analysis.  These histograms are only those related to GRETINA, and 
   are primarily for diagnostic purposes.

   The class contains the histograms themselves, and a function for writing out
   the histograms to the ROOT file populated in analysis.
*/

class gHistos : public TObject {
 public:
  Float_t xtalkThresh;
  TH1F *eRaw[MAXCHANNELS], *eCal[MAXCHANNELS];
  TH1F *mult2Pair[40][40];

 public:
  gHistos();
  ~gHistos();
  void writeHistos(Int_t ctrl);
  
  ClassDef(gHistos, 1);
};

/* \class SuperPulse
   \brief Class containing the parameters required to perform a superpulse
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


class GRETINA : public TObject {

 public:
  rotationMatrix rot;
  GRETINAVariables var;

  SuperPulse sp;

  /* Temporary structures for getting GRETINA data */
  g1GammaEvent g1X;
  g2CrystalEvent g2X; g2IntPt pt;
  g3CrystalEvent g3X; g3ChannelEvent g3ch;
  vector<g3ChannelEvent> g3Temp;
  historyEvent gH; 

  unsigned char gBuf[32*32*1024];

  g3HistoryEvent g3H;
  g3OUT g3out;
  g2OUT g2out;
  g1OUT g1out;
  g4SimOUT gSimOut;
  
  Bank88 b88;

  Track track;

  gHistos gHist;

 public:
  /* Basic geometry-related variables.  Maybe there is a better
     place for these? */

  Int_t cloverModuleNumber;
  Float_t cloverSlope[40];
  Float_t cloverOffset[40];

  Float_t wfMinCrossTime, wfMaxCrossTime;
  
  /* Segment ordering from channel mapping */
  Double_t QNormal[4][40], Q1Special[4][40], Q2Special[4][40];

  UShort_t scalerBuf[8192];

 public:
  GRETINA() { ; }
  ~GRETINA() { ; }
  void Initialize();
  void Reset();
  Float_t getDopplerSimple(TVector3 xyz, Float_t beta);

  Int_t getMode3(FILE *inf, Int_t evtLength, counterVariables *cnt,
		 controlVariables *ctrl);
  Int_t getMode3History(FILE *inf, Int_t evtLength, long long int hTS, counterVariables *cnt);
  Int_t getMode2(FILE *inf, Int_t evtLength, counterVariables *cnt);
  Int_t getMode1(FILE *inf, counterVariables *cnt);
  Int_t getSimulated(FILE *inf);
  void getScaler(FILE *inf, Int_t evtLength);
  Int_t getBank88(FILE *inf, Int_t evtLength, counterVariables *cnt);

  Int_t analyzeMode2(g2CrystalEvent *g2);
  void calibrateMode2CC(Int_t crystal, mode2ABCD5678 *g2, 
			g2CrystalEvent *g2crystal);
  void calibrateMode2CC(Int_t crystal, mode2ABCD6789 *g2, 
			g2CrystalEvent *g2crystal);
  void analyzeMode3(controlVariables *ctrl);
  void calibrateMode3(g3ChannelEvent *g3);
  void calibrateMode3SP(g3ChannelEvent *g3);

  void fillHistos(Int_t ctrl);

  void checkSPIntegrity();

  Int_t fillShell2Track();
  void fillMode1(Int_t trackStatus);

 public:
  ClassDef(GRETINA, 1);

};

#endif
