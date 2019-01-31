// @(#)root/test:$Id: guiviewer.h 20882 2007-11-19 11:31:26Z rdm $
// Author: Brett Viren   04/15/2001

#ifndef ROOT_TGFrame
#include "TGFrame.h"
#endif

#include "GRETINA.h"

class TList;
class TCanvas;
class TRootEmbeddedCanvas;
class TGaxis;
class TGDoubleSlider;
class TFile;
class TTree;
class TH1F;
class TGTextEntry;
class TGMainFrame;
class TGNumberEntryField;
class TGCheckButton;

class Viewer : public TGMainFrame {

private:
  TGMainFrame *fMain;
  TList               *fCleanup;
  TCanvas             *segmentCanvas, *ccCanvas, *xyCanvas, *rzCanvas;
  TRootEmbeddedCanvas *segmentCanvasE, *ccCanvasE, *xyCanvasE, *rzCanvasE;
  TFile               *dataIN;
  TTree               *treeIN;

  TGTextEntry *evNumInfo, *cryNumInfo, *ccInfo, *errorInfo;
  TGNumberEntryField *qNumInfo, *posNumInfo, *eLowLimit, *eHighLimit;
  TGCheckButton *fixCrystal;

  Int_t selectSpecificCrystal;

  TGTextEntry *intptsETitle, *intptsXTitle, *intptsYTitle, *intptsZTitle, *intptsSegTitle, *intptsRTitle;
  TGTextEntry *intptsSeg[16];
  TGTextEntry *intptsE[16];
  TGTextEntry *intptsX[16];
  TGTextEntry *intptsY[16];
  TGTextEntry *intptsZ[16];
  TGTextEntry *intptsR[16];
    
  g2OUT *g2;
   
   Int_t eventNum;
   Int_t crystalNum; Int_t crystalsInEvent;
   TH1F *totalPlot;
   TH1F *segmentPlot[36];
   TH1F *ccPlot[4];

   Int_t Q1Special[4][40];
   Int_t Q2Special[4][40];
   Int_t QNormal[4][40];
   
public:
   Viewer(const TGWindow *win, TString filename);
   virtual ~Viewer();
   void DoButton();
   void UpdateEventDisplay();
   void DoExpand(Int_t event, Int_t px, Int_t py, TObject *sel);
   ClassDef(Viewer,0) //GUI example
};
