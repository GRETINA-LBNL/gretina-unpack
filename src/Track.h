#ifndef __TRACK_H
#define __TRACK_H

#include "Riostream.h"
#include "Rtypes.h"
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>

#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

/* Tracking define statements */

#define MAXNUMDET 30
#define TOOMANYHITS 98

#define MAXSHELLHITS 500
#define MAXCLUSTERHITS 200
#define MAXNUMSEG 9

#define NRETURNCODES 30

#define MAXPERM MAXNUMSEG*8*7*6*5*4*3*2
#define MAXPNUM 10000
#define MAXFACULTY 20

#define MAXFOM 101

//#define MAXFLOAT 3.40282347e+38
//#define MAXINT 2147483647
#define SHORTLEN 1024
#define LONGLEN 16384

#define MAXDETNUM 200

struct PM {
  Int_t np1, np2;
  Int_t p1[20], p2[20];
  Float_t s1, s2;
  Float_t likely1, likely2;
  Float_t dist;
  Float_t FOM1, FOM2;
};

class TrackClusterIP : public TObject {
 public:
  TVector3 xyz;
  Float_t eDet;
  Int_t order;
  long long int timestamp;
  Float_t t0;
  Int_t shellPos;
  Int_t detNum;
  /* For simulated data only */
  Float_t eSumOrig;
  Int_t origPos;

 public:
  Double_t xx();
  Double_t yy();
  Double_t zz();

  ClassDef(TrackClusterIP, 1);
};

class TrackCluster : public TObject {

 public:
  Int_t valid;
  Int_t nDet;
  Int_t tracked;
  Float_t FOM;
  Float_t eSum;
  Int_t trackNum;
  Int_t bestPermutation;
  Int_t processed;
  
  TrackClusterIP intpts[MAXNUMDET];

  Float_t badTheoAngPenalty;
  
 public:
  Int_t findFOM(TVector3 targetPos);
  Float_t findComptonAngle(Float_t eg, Float_t ee, Float_t *thc);
  Int_t sortSegEnergies();
  void Reset();
  void Clear();
  void PrintCluster();

  ClassDef(TrackCluster, 1);
};

class TrackShell : public TObject {

 public:
  Int_t nHit;
  Float_t t0[MAXSHELLHITS];
  Float_t chiSq[MAXSHELLHITS];
  Float_t normChiSq[MAXSHELLHITS];
  long long int timestamp[MAXSHELLHITS];
  Int_t crystalID[MAXSHELLHITS];
  TVector3 xyz[MAXSHELLHITS];
  Float_t eDet[MAXSHELLHITS];
  Float_t eSum[MAXSHELLHITS];
  Int_t origPos[MAXSHELLHITS];
  Int_t detNum[MAXSHELLHITS];
  
  Int_t module[MAXSHELLHITS];
  Int_t crystalType[MAXSHELLHITS];
  
  /* Deduced cluster information */

  Int_t numClusters;
  Int_t clusterNum[MAXSHELLHITS];
  
  /* For simulated data only... */
  Int_t knownNumClusters;
  Int_t knownClusterNum[MAXSHELLHITS];

 public:
  Double_t xx(Int_t i);
  Double_t yy(Int_t i);
  Double_t zz(Int_t i);
  Double_t nxx(Int_t i);
  Double_t nyy(Int_t i);
  Double_t nzz(Int_t i);
  Double_t relAngle(Int_t i, Int_t j);
  Double_t polAngle(Int_t i);
  Double_t aziAngle(Int_t i);
  void Print();
  void PrintClusters();

  ClassDef(TrackShell, 1);
};

class TrackCtrl : public TObject {

 public:
  Float_t targetX, targetY, targetZ;

  Int_t maxEvents;
  Int_t nPrint;
  Int_t nTracks;

  Int_t reTrack;

  Float_t alpha[MAXSHELLHITS], gamma, data;
  Float_t FOMJump, FOMGoodEnough;

  Int_t iterations;
  Int_t removeInvalids;

  Int_t fixDoubleHitsCode;
  Int_t fixDoubleHits;

  Int_t enabled[MAXDETNUM];

  TString trackDataFileName;

  /* Energy - interaction point stuff */
  Float_t nDetELim_lo[MAXNUMSEG], nDetELim_hi[MAXNUMSEG];
  Float_t nDetELim_FOM[MAXNUMSEG];
  Int_t useCCEnergy, useSegEnergy;

  /* Split cluster function controls */
  Int_t splitClusters[2];
  Float_t splitClusters_kickoutFOM[2];
  Float_t splitClusters_thresh[2];
  Int_t splitClusters_minN[2];
  Int_t splitClusters_maxN[2];
  Int_t splitClusters_maxTry[2];
  Float_t goodImproveFraction;

  /* Recluster function controls */
  Int_t reCluster[2];
  Float_t reCluster_kickoutFOM[2];
  Float_t reCluster_thresh[2];
  Int_t reCluster_minN[2];
  Int_t reCluster_maxTry[2];
  Float_t reCluster_reduxFactor[2];

  Float_t untracked_kickoutFOM;
  
  /* Combine cluster function controls */
  Int_t combineClusters;
  Float_t combineClusters_kickoutFOM;
  Float_t combineClusters_thresh;
  Int_t combineClusters_maxN;
  Float_t combineMaxDist;

  /* Matchmaker function controls */
  Int_t matchMaker;
  Float_t matchMakerMaxDist;
  Float_t matchMaker_kickoutFOM;

  /* Single hit contols */
  Float_t singlesFOM;
  Int_t markSingleHitRange;
  Float_t singleHitRange[LONGLEN];
  Float_t singleHitRangeFOM;

  Int_t trackOps[MAXNUMSEG];
  Int_t jumpGroupLength[MAXNUMSEG];

  Float_t badTheoAngPenalty;

 public:
  Int_t ReadChatFile(TString fileName);
  void CheckNumArgs(Int_t required, Int_t actual, char *str);
  Int_t strDecomp(char *str, Int_t dim, Int_t *yy);
  void PrintPartial();
  Int_t trackOpt(Int_t i); 

  ClassDef(TrackCtrl, 1);
};

class TrackStat : public TObject {

 public:
  UInt_t nEvents;
  long trackGetCalls;
  Int_t nHit[164];
  long trackPassCnt;
  UInt_t trackingCalls;
  UInt_t trackError[6][NRETURNCODES];
  UInt_t firstClusterTrackCalls;
  long trackFOMokay;
  Int_t permHit[MAXNUMSEG][MAXPERM];
  UInt_t reClusterReturns[2][NRETURNCODES];
  UInt_t reClusterTrackCalls[2][NRETURNCODES];
  UInt_t splitClusterReturns[2][NRETURNCODES];
  UInt_t splitClusterTrackCalls[2][NRETURNCODES];
  UInt_t combineClusterReturns[NRETURNCODES];
  UInt_t combineClusterTrackCalls[NRETURNCODES];
  UInt_t matchMakerReturns[NRETURNCODES];
  UInt_t matchMakerTrackCalls[NRETURNCODES];
  long long int singleHitOutOfRange;
  long long int badNDetELim;
  Double_t sumFOM;
  Double_t sumFOM_NT;
  Int_t nSumFOM_NT;
  UInt_t nClusters;
  long notTracked;
  long nPerm;
  long doubleSegHits;
  long segHits;
  long numInputDataLines;
  long goodTrack;
  long badTrack;
  UInt_t noTrueClusters;
  long long int nDetTooLarge;
  long long int fullEnergyEv;
  long long int comptonLossEnergyEv;
  long long int gammasForOutput;
  long long int gammasForOutputOK1;
  long long int gammasForOutputOK2;
  long long int gammasWritten;
  long long int trackedEventsWritten;

 public: 
  void Zero();

  ClassDef(TrackStat, 1);
};

class TrackPerm : public TObject {

 public:
  Int_t nPerm[MAXNUMSEG];
  Int_t lookup[MAXNUMSEG][MAXPERM][MAXNUMSEG];

 public:
  void Initialize();

  ClassDef(TrackPerm, 1);
};

class Track : public TObject {
 public:
  TrackShell shell;
  TrackCluster clust[MAXCLUSTERHITS];
  
  TrackStat stat;
  TrackCtrl ctrl;
  TrackPerm perm;

  TVector3 targetPos;
  Int_t nClusters;

  long long faculty[MAXFACULTY];

 public:
  Track() { ; }
  ~Track() { ; }
  void Initialize();
  Int_t findTargetPos();
  Int_t findClusters();
  Int_t transferClusters();
  Int_t trackEvent();
  Int_t doTrack(Int_t mode, Int_t iClust);
  Int_t reCluster(Int_t ii, Int_t m);
  Int_t splitCluster(Int_t ii, Int_t m);
  Int_t combineCluster (Int_t focusC);
  Int_t matchMaker(Int_t ii, Int_t jj);
  void PrintValidClusters();

  ClassDef(Track, 1);
};

#endif
