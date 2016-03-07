// @(#)root/test:$Id: guiviewer.cxx 20882 2007-11-19 11:31:26Z rdm $
// Author: Brett Viren   04/15/2001

// guiviewer.cxx: GUI test program showing TGTableLayout widget manager,
// embedded ROOT canvases, axis and sliders.
// To run it do: make guiviewer; guiviewer

#include "decompViewer.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TGTableLayout.h"
#include "TGFrame.h"
#include "TList.h"
#include "TGWindow.h"
#include "TRootEmbeddedCanvas.h"
#include "TEllipse.h"
#include "TBox.h"
#include "TArrow.h"
#include "TPolyLine.h"
#include "TArc.h"
#include "TGButton.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TGDoubleSlider.h"
#include "TGaxis.h"
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TGTextEntry.h"
#include "TGraph.h"
#include "TMath.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"

Viewer::Viewer(const TGWindow *win, TString filename) : TGMainFrame(win,1000,800)
{
  /* Initialize detector segment maps... */
   for (Int_t i=0; i<9; i++) {
     for (Int_t j=0; j<4; j++) {
       Q1Special[j][i] = i;
       Q2Special[j][i] = i;
       QNormal[j][i] = i;
       Q1Special[j][9] = 36;
       Q2Special[j][9] = 36;
       QNormal[j][9] = 36;
     }
   }
   for (Int_t i=10; i<19; i++) {
     for (Int_t j=0; j<4; j++) {
       Q1Special[j][i] = i-1;
       Q2Special[j][i] = i-1;
       QNormal[j][i] = i-1;
       Q1Special[j][19] = 37;
       Q2Special[j][19] = 37;
       QNormal[j][19] = 37;
     }
   }
   for (Int_t i=20; i<29; i++) {
     for (Int_t j=0; j<4; j++) {
       Q1Special[j][i] = i-2;
       Q2Special[j][i] = i-2;
       QNormal[j][i] = i-2;
       Q1Special[j][29] = 38;
       Q2Special[j][29] = 38;
       QNormal[j][29] = 38;
     }
   }
   for (Int_t i=30; i<39; i++) {
     for (Int_t j=0; j<4; j++) {
       Q1Special[j][i] = i-3;
       Q2Special[j][i] = i-3;
       QNormal[j][i] = i-3;
       Q1Special[j][39] = 39;
       Q2Special[j][39] = 39;
       QNormal[j][39] = 39;
     }
   }
   
   for (Int_t i=0; i<3; i+=2) {
     Q1Special[i][ 0] =  1; Q1Special[i][ 1] =  8; Q1Special[i][ 2] =  3; Q1Special[i][ 3] =  4;
     Q1Special[i][ 4] = 35; Q1Special[i][ 5] =  0; Q1Special[i][ 6] =  7; Q1Special[i][ 7] = 14;
     Q1Special[i][ 8] =  9; Q1Special[i][10] = 10; Q1Special[i][11] = 29; Q1Special[i][12] = 30;
     Q1Special[i][13] = 13; Q1Special[i][14] = 20; Q1Special[i][15] = 15; Q1Special[i][16] = 16;
     Q1Special[i][17] = 23; Q1Special[i][18] = 24;
   }
   
   Q1Special[0][20] = 19; Q1Special[0][21] = 26; Q1Special[0][22] = 21; Q1Special[0][23] = 22;
   Q1Special[0][24] = 17; Q1Special[0][25] = 18; Q1Special[0][26] = 25; Q1Special[0][27] = 32;
   Q1Special[0][28] = 27; Q1Special[0][30] = 28; Q1Special[0][31] = 11; Q1Special[0][32] = 12;
   Q1Special[0][33] = 31; Q1Special[0][34] =  2; Q1Special[0][35] = 33; Q1Special[0][36] = 34;
   Q1Special[0][37] =  5; Q1Special[0][38] =  6;
   
   Q1Special[2][20] = 19; Q1Special[2][21] = 26; Q1Special[2][22] = 21; Q1Special[2][23] = 22;
   Q1Special[2][24] = 17; Q1Special[2][25] = 18; Q1Special[2][26] = 25; Q1Special[2][27] = 32;
   Q1Special[2][28] = 27; Q1Special[2][30] = 28; Q1Special[2][31] = 11; Q1Special[2][32] = 12;
   Q1Special[2][33] = 31; Q1Special[2][34] =  2; Q1Special[2][35] = 33; Q1Special[2][36] = 34;
   Q1Special[2][37] =  5; Q1Special[2][38] =  6;
   
   Q2Special[0][17] = 22; Q2Special[0][24] = 16;
   Q2Special[2][17] = 22; Q2Special[2][24] = 16;

   /* Initialize eventNum */
   eventNum = -1;
   crystalNum = 0; crystalsInEvent = 0;
   totalPlot = new TH1F();
   for (Int_t i=0; i<36; i++) {
     if (i<4) { ccPlot[i] = new TH1F(); }
     segmentPlot[i] = new TH1F();
   }

   dataIN = new TFile(filename.Data(), "READ");
   treeIN = (TTree*)dataIN->Get("teb");
   g2 = new g2OUT();
   treeIN->SetBranchAddress("g2", &g2);

   selectSpecificCrystal = 0;

   fCleanup = new TList;

   SetLayoutManager(new TGHorizontalLayout(this));

   /* Vertical frame -- on the left hand side */
   TGVerticalFrame *leftMenu = new TGVerticalFrame(this, 224, 576, kVerticalFrame);
   leftMenu->SetName("DecompViewerMenu");
   
   fCleanup->Add(leftMenu);

   TGTextButton *buttonNext = new TGTextButton(leftMenu, "Next Event", 0, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
   buttonNext->SetTextJustify(36);
   buttonNext->SetMargins(0,0,0,0);
   buttonNext->SetWrapLength(-1);
   buttonNext->Resize(224,24);
   leftMenu->AddFrame(buttonNext, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
   buttonNext->Connect("Clicked()","Viewer",this,"DoButton()");  

   TGTextButton *buttonPrevious = new TGTextButton(leftMenu, "Previous Event", 1, TGTextButton::GetDefaultGC()(), TGTextButton::GetDefaultFontStruct(), kRaisedFrame);
   buttonPrevious->SetTextJustify(36);
   buttonPrevious->SetMargins(0,0,0,0);
   buttonPrevious->SetWrapLength(-1);
   buttonPrevious->Resize(224,24);
   leftMenu->AddFrame(buttonPrevious, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
   buttonPrevious->Connect("Clicked()","Viewer",this,"DoButton()");   
   
   TGFont *ufont;
   ufont = gClient->GetFont("-*-helvetica-medium-r-*-*-12-*-*-*-*-*-iso8859-1");
   TGGC *uGC;
   GCValues_t ipEntry;
   ipEntry.fMask = kGCForeground | kGCBackground | kGCFillStyle | kGCFont | kGCGraphicsExposures;
   gClient->GetColorByName("#000000", ipEntry.fForeground);
   gClient->GetColorByName("#e0e0e0", ipEntry.fBackground);
   ipEntry.fFillStyle = kFillSolid;
   ipEntry.fFont = ufont->GetFontHandle();
   ipEntry.fGraphicsExposures = kFALSE;
   uGC = gClient->GetGC(&ipEntry, kTRUE);
   ULong_t ucolor;
   gClient->GetColorByName("#cccccc", ucolor);
 
   /* Event details */
   TGGroupFrame *evFrame = new TGGroupFrame(leftMenu, "Event Detail", kVerticalFrame);

   TGHorizontalFrame *evNumbers = new TGHorizontalFrame(evFrame);
   TGLabel *evNumLabel = new TGLabel(evNumbers, "Event#");
   evNumInfo = new TGTextEntry(evNumbers, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   evNumInfo->SetMaxLength(4096);
   evNumInfo->SetAlignment(kTextCenterX);
   evNumInfo->SetText("");
   evNumInfo->Resize(60, 20);
   TGLabel *cryNumLabel = new TGLabel(evNumbers, "CrystalEvent#");
   cryNumInfo = new TGTextEntry(evNumbers, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   cryNumInfo->SetMaxLength(4096);
   cryNumInfo->SetAlignment(kTextCenterX);
   cryNumInfo->SetText("");
   cryNumInfo->Resize(60, 20);

   evNumbers->AddFrame(evNumLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evNumbers->AddFrame(evNumInfo, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evNumbers->AddFrame(cryNumLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evNumbers->AddFrame(cryNumInfo, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evFrame->AddFrame(evNumbers, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 1, 1, 1, 0));

   TGHorizontalFrame *evCrystalDetails = new TGHorizontalFrame(evFrame);
   TGLabel *qNumLabel = new TGLabel(evCrystalDetails, "Quad#");
   qNumInfo = new TGNumberEntryField(evCrystalDetails, -1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
   qNumInfo->SetBackgroundColor(ucolor);
   qNumInfo->Resize(60, 20);
   TGLabel *posNumLabel = new TGLabel(evCrystalDetails, "Pos#");
   posNumInfo = new TGNumberEntryField(evCrystalDetails, -1, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
   posNumInfo->SetBackgroundColor(ucolor);
   posNumInfo->Resize(60, 20);
   fixCrystal = new TGCheckButton(evCrystalDetails, "Fix?", 2);
   fixCrystal->Connect("Clicked()","Viewer",this,"DoButton()");  

   evCrystalDetails->AddFrame(qNumLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evCrystalDetails->AddFrame(qNumInfo, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evCrystalDetails->AddFrame(posNumLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evCrystalDetails->AddFrame(posNumInfo, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evCrystalDetails->AddFrame(fixCrystal, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evFrame->AddFrame(evCrystalDetails, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 1, 1, 1, 0));
   
   TGHorizontalFrame *evDetails = new TGHorizontalFrame(evFrame);
   TGLabel *ccLabel = new TGLabel(evDetails, "CC E");
   ccInfo = new TGTextEntry(evDetails, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   ccInfo->SetMaxLength(4096);
   ccInfo->SetAlignment(kTextCenterX);
   ccInfo->SetText("");
   ccInfo->Resize(60, 20);
   TGLabel *errorLabel = new TGLabel(evDetails, "Decomp Error");
   errorInfo = new TGTextEntry(evDetails, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   errorInfo->SetMaxLength(4096);
   errorInfo->SetAlignment(kTextCenterX);
   errorInfo->SetText("");
   errorInfo->Resize(60, 20);

   evDetails->AddFrame(ccLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evDetails->AddFrame(ccInfo, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evDetails->AddFrame(errorLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evDetails->AddFrame(errorInfo, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evFrame->AddFrame(evDetails, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 1, 1, 1, 0));

   TGHorizontalFrame *eCut = new TGHorizontalFrame(evFrame);
   TGLabel *eCutLabel = new TGLabel(eCut, "CC E Limits");
   eLowLimit = new TGNumberEntryField(eCut, -1, 0, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
   eLowLimit->SetBackgroundColor(ucolor);
   eLowLimit->Resize(60, 20);
   TGLabel *ccELabel = new TGLabel(eCut, "<= E <=");
   eHighLimit = new TGNumberEntryField(eCut, -1, 20000, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive);
   eHighLimit->SetBackgroundColor(ucolor);
   eHighLimit->Resize(60, 20);

   eCut->AddFrame(eCutLabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   eCut->AddFrame(eLowLimit, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   eCut->AddFrame(ccELabel, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   eCut->AddFrame(eHighLimit, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   evFrame->AddFrame(eCut, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 1, 1, 1, 0));

   evFrame->MapWindow();
   leftMenu->AddFrame(evFrame, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX, 2, 2, 2, 2));

   /* Group menu */
   TGGroupFrame *ipFrame = new TGGroupFrame(leftMenu, "Interaction Points", kVerticalFrame);
   
   TGHorizontalFrame *titles = new TGHorizontalFrame(ipFrame);
   TGHorizontalFrame *ipLines[16];

   intptsSegTitle = new TGTextEntry(titles, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   intptsSegTitle->SetMaxLength(4096);
   intptsSegTitle->SetAlignment(kTextCenterX);
   intptsSegTitle->SetText("Seg#");
   intptsSegTitle->Resize(40, 20);
   titles->AddFrame(intptsSegTitle, new TGLayoutHints(kLHintsLeft | kLHintsTop, 1, 1, 1, 1));
   
   intptsXTitle = new TGTextEntry(titles, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   intptsXTitle->SetMaxLength(4096);
   intptsXTitle->SetAlignment(kTextCenterX);
   intptsXTitle->SetText("X");
   intptsXTitle->Resize(40, 20);
   titles->AddFrame(intptsXTitle, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
   
   intptsYTitle = new TGTextEntry(titles, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   intptsYTitle->SetMaxLength(4096);
   intptsYTitle->SetAlignment(kTextCenterX);
   intptsYTitle->SetText("Y");
   intptsYTitle->Resize(40, 20);
   titles->AddFrame(intptsYTitle, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
   
   intptsZTitle = new TGTextEntry(titles, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   intptsZTitle->SetMaxLength(4096);
   intptsZTitle->SetAlignment(kTextCenterX);
   intptsZTitle->SetText("Z");
   intptsZTitle->Resize(40, 20);
   titles->AddFrame(intptsZTitle, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));

   intptsRTitle = new TGTextEntry(titles, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   intptsRTitle->SetMaxLength(4096);
   intptsRTitle->SetAlignment(kTextCenterX);
   intptsRTitle->SetText("R");
   intptsRTitle->Resize(40, 20);
   titles->AddFrame(intptsRTitle, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
   
   intptsETitle = new TGTextEntry(titles, new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
   intptsETitle->SetMaxLength(4096);
   intptsETitle->SetAlignment(kTextCenterX);
   intptsETitle->SetText("E");
   intptsETitle->Resize(50, 20);
   titles->AddFrame(intptsETitle, new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
   
   ipFrame->AddFrame(titles, new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 1, 1, 1, 0));

   for (Int_t i=0; i<16; i++) {
     ipLines[i] = new TGHorizontalFrame(ipFrame);
     intptsSeg[i] = new TGTextEntry(ipLines[i], new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
     intptsSeg[i]->SetMaxLength(4096);
     intptsSeg[i]->SetAlignment(kTextCenterX);
     intptsSeg[i]->SetText("");
     intptsSeg[i]->Resize(40, 20);
     ipLines[i]->AddFrame(intptsSeg[i], new TGLayoutHints(kLHintsLeft | kLHintsTop, 1, 1, 1, 1));
     
     intptsX[i] = new TGTextEntry(ipLines[i], new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
     intptsX[i]->SetMaxLength(4096);
     intptsX[i]->SetAlignment(kTextCenterX);
     intptsX[i]->SetText("");
     intptsX[i]->Resize(40, 20);
     ipLines[i]->AddFrame(intptsX[i], new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
     
     intptsY[i] = new TGTextEntry(ipLines[i], new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
     intptsY[i]->SetMaxLength(4096);
     intptsY[i]->SetAlignment(kTextCenterX);
     intptsY[i]->SetText("");
     intptsY[i]->Resize(40, 20);
     ipLines[i]->AddFrame(intptsY[i], new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
     
     intptsZ[i] = new TGTextEntry(ipLines[i], new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
     intptsZ[i]->SetMaxLength(4096);
     intptsZ[i]->SetAlignment(kTextCenterX);
     intptsZ[i]->SetText("");
     intptsZ[i]->Resize(40, 20);
     ipLines[i]->AddFrame(intptsZ[i], new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));

     intptsR[i] = new TGTextEntry(ipLines[i], new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
     intptsR[i]->SetMaxLength(4096);
     intptsR[i]->SetAlignment(kTextCenterX);
     intptsR[i]->SetText("");
     intptsR[i]->Resize(40, 20);
     ipLines[i]->AddFrame(intptsR[i], new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
         
     intptsE[i] = new TGTextEntry(ipLines[i], new TGTextBuffer(14), -1, uGC->GetGC(), ufont->GetFontStruct(), kOwnBackground, ucolor);
     intptsE[i]->SetMaxLength(4096);
     intptsE[i]->SetAlignment(kTextCenterX);
     intptsE[i]->SetText("");
     intptsE[i]->Resize(50, 20);
     ipLines[i]->AddFrame(intptsE[i], new TGLayoutHints(kLHintsCenterX | kLHintsTop, 1, 1, 1, 1));
    
     ipFrame->AddFrame(ipLines[i], new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 1, 1, 0, 0));
   }

   ipFrame->MapWindow();
   leftMenu->AddFrame(ipFrame, new TGLayoutHints(kLHintsCenterX | kLHintsExpandX, 2, 2, 2, 2));
   ipFrame->MoveResize(0,224,224,576);
   
   AddFrame(leftMenu, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
   leftMenu->MoveResize(0,0,224,576);

 /* Embedded canvas --  main for waveforms */
   TGVerticalFrame *centerCanvases = new TGVerticalFrame(this, 224, 352, kVerticalFrame);
   segmentCanvasE = new TRootEmbeddedCanvas(0, centerCanvases, 616, 440);
   segmentCanvasE->SetName("SegmentWaveformCanvas");
   Int_t segmentCanvasID = segmentCanvasE->GetCanvasWindowId();
   segmentCanvas = new TCanvas("SegmentCanvas", 10, 10, segmentCanvasID);
   segmentCanvasE->AdoptCanvas(segmentCanvas);
   centerCanvases->AddFrame(segmentCanvasE, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
   segmentCanvasE->MoveResize(224, 0, 616, 440);
   segmentCanvas->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)","Viewer",this,"DoExpand(Int_t, Int_t, Int_t, TObject*)");

   /* Embedded canvas -- CC waveforms */
   ccCanvasE = new TRootEmbeddedCanvas(0, centerCanvases, 616, 136);
   ccCanvasE->SetName("CCWaveformCanvas");
   Int_t ccCanvasID = ccCanvasE->GetCanvasWindowId();
   ccCanvas = new TCanvas("ccCanvas", 10, 10, ccCanvasID);
   ccCanvasE->AdoptCanvas(ccCanvas);
   centerCanvases->AddFrame(ccCanvasE, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
   ccCanvasE->MoveResize(224, 440, 616, 145);
   ccCanvas->Connect("ProcessedEvent(Int_t, Int_t, Int_t, TObject*)","Viewer",this,"DoExpand(Int_t, Int_t, Int_t, TObject*)");
 
   AddFrame(centerCanvases, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

   /* Interaction point location canvases */
   TGVerticalFrame *rightCanvases = new TGVerticalFrame(this, 224, 352, kVerticalFrame);
   xyCanvasE = new TRootEmbeddedCanvas(0,rightCanvases,291,291);
   xyCanvasE->SetName("xyCanvas");
   Int_t xyCanvasEID = xyCanvasE->GetCanvasWindowId();
   xyCanvas = new TCanvas("xy", 10, 10, xyCanvasEID);
   xyCanvasE->AdoptCanvas(xyCanvas);
   rightCanvases->AddFrame(xyCanvasE, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   xyCanvasE->MoveResize(840,0,291,292);
  
   rzCanvasE = new TRootEmbeddedCanvas(0,rightCanvases,291,291);
   rzCanvasE->SetName("rzCanvas");
   Int_t rzCanvasEID = rzCanvasE->GetCanvasWindowId();
   rzCanvas = new TCanvas("rz", 10, 10, rzCanvasEID);
   rzCanvasE->AdoptCanvas(rzCanvas);
   rightCanvases->AddFrame(rzCanvasE, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));
   rzCanvasE->MoveResize(840,291,291,292);
  
   AddFrame(rightCanvases, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

   // exit on close window
   Connect("CloseWindow()","TApplication",gApplication,"Terminate(=0)");

   MapSubwindows();
   Resize(GetDefaultSize());
   MapWindow();
   Resize(1230, 600);
}

Viewer::~Viewer()
{
   segmentCanvas->Clear();
   delete segmentCanvas;
   fCleanup->Delete();
   delete fCleanup;
}

void Viewer::DoButton()
{
   TGButton* button = (TGButton*)gTQSender;
   UInt_t id = button->WidgetId();

   segmentCanvas->Clear();
   segmentCanvas->Divide(6, 6, 0.01, 0.01);
   segmentCanvas->cd(1);

   ccCanvas->Clear();
   ccCanvas->Divide(4, 1, 0.01, 0.01);
   ccCanvas->cd(1);

   xyCanvas->Clear();
   rzCanvas->Clear();

   switch (id) {
   case 0: { // Next Event
     /* Clear old histograms... */
     totalPlot->Delete();
     for (Int_t i=0; i<36; i++) {
       if (i<4) { ccPlot[i]->Delete(); }
       segmentPlot[i]->Delete();
     }
     
     Int_t gotTheNextOne = 0;

     /* Are there still more crystals in the current event? */
     if (crystalNum >= crystalsInEvent) {
       eventNum++;
       treeIN->GetEvent(eventNum);
       while (!gotTheNextOne) {
	 while (g2->xtals.size() <= 0) { eventNum++;  treeIN->GetEvent(eventNum); }
	 crystalsInEvent = g2->xtals.size();
	 crystalNum=0;
	 if (!selectSpecificCrystal && 
	     g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
	     g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) { gotTheNextOne = 1; }
	 else {
	   while (!gotTheNextOne && crystalNum < crystalsInEvent) {
	     if (g2->xtals[crystalNum].quadNum == qNumInfo->GetIntNumber() &&
		 (((g2->xtals[crystalNum].crystalNum-1)%4)+1) == posNumInfo->GetIntNumber() && 
		 g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
		 g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	       gotTheNextOne = 1;
	     } else {
	       crystalNum++;
	       gotTheNextOne = 0;
	     }
	   }
	   if (crystalNum >= crystalsInEvent && !gotTheNextOne) {
	     eventNum++;  treeIN->GetEvent(eventNum);
	   }
	 }
       }
     } else {
       while (!gotTheNextOne && crystalNum < crystalsInEvent) {
	 if (!selectSpecificCrystal &&
	     g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
	     g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	   gotTheNextOne = 1;
	 } else if (g2->xtals[crystalNum].quadNum == qNumInfo->GetIntNumber() &&
		    (((g2->xtals[crystalNum].crystalNum-1)%4)+1) == posNumInfo->GetIntNumber() &&
		    g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
		    g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	   gotTheNextOne = 1;
	 } else {
	   crystalNum++;
	   gotTheNextOne = 0;
	 }
       }
       while (!gotTheNextOne) {
	 eventNum++;  treeIN->GetEvent(eventNum);
	 while (g2->xtals.size() <= 0) { eventNum++;  treeIN->GetEvent(eventNum); }
	 crystalsInEvent = g2->xtals.size();
	 crystalNum=0;
	 if (!selectSpecificCrystal &&
	     g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
	     g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) { gotTheNextOne = 1; }
	 else {
	   while (!gotTheNextOne && crystalNum < crystalsInEvent) {
	     if (g2->xtals[crystalNum].quadNum == qNumInfo->GetIntNumber() &&
		 (((g2->xtals[crystalNum].crystalNum-1)%4)+1) == posNumInfo->GetIntNumber() &&
		 g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
		 g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	       gotTheNextOne = 1;
	     } else {
	       crystalNum++;
	       gotTheNextOne = 0;
	     }
	   }
	 }
       }
     }

     UpdateEventDisplay();
     break;
   }
   case 1: { // Previous Event
     /* Clear old histograms... */
     totalPlot->Delete();
     for (Int_t i=0; i<36; i++) {
       if (i<4) { ccPlot[i]->Delete(); }
       segmentPlot[i]->Delete();
     }
     
     Int_t gotTheNextOne = 0;

     /* Are there previous crystals in the current event? */
     if (crystalNum>1 && crystalsInEvent>1) {
       crystalNum--; crystalNum--;
       while (!gotTheNextOne && crystalNum>=0) {
	 if (!selectSpecificCrystal && 
	     g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
	     g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	   gotTheNextOne = 1;
	 } else if (g2->xtals[crystalNum].quadNum == qNumInfo->GetIntNumber() &&
	     (((g2->xtals[crystalNum].crystalNum-1)%4)+1) == posNumInfo->GetIntNumber() &&
		    g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
		    g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	   gotTheNextOne = 1;
	 } else {
	   crystalNum--;
	   gotTheNextOne = 0;
	 }
       }
       if (crystalNum < 0 && !gotTheNextOne) {
	 eventNum--;  treeIN->GetEvent(eventNum);
       }
       while (!gotTheNextOne) {
	 while (g2->xtals.size() <= 0) { eventNum--;  treeIN->GetEvent(eventNum); }
	 crystalsInEvent = g2->xtals.size();
	 crystalNum = (crystalsInEvent-1);
	 if (!selectSpecificCrystal &&
	     g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
	     g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) { gotTheNextOne = 1; }
	 else {
	   while (!gotTheNextOne && crystalNum >= 0) {
	     if (g2->xtals[crystalNum].quadNum == qNumInfo->GetIntNumber() &&
		 (((g2->xtals[crystalNum].crystalNum-1)%4)+1) == posNumInfo->GetIntNumber() &&
		 g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
		 g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	       gotTheNextOne = 1;
	     } else {
	       crystalNum--;
	       gotTheNextOne = 0;
	     }
	   }
	   if (crystalNum < 0 && !gotTheNextOne) {
	     eventNum--;  treeIN->GetEvent(eventNum);
	   }
	 }
       }
     } else {
       eventNum--;  treeIN->GetEvent(eventNum);
       while (!gotTheNextOne) {
	 while (g2->xtals.size() <= 0) { eventNum--;  treeIN->GetEvent(eventNum); }
	 crystalsInEvent = g2->xtals.size();
	 crystalNum = (crystalsInEvent-1);
	 if (!selectSpecificCrystal &&
	     g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
	     g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) { gotTheNextOne = 1; }
	 else {
	   while (!gotTheNextOne && crystalNum >= 0) {
	     if (g2->xtals[crystalNum].quadNum == qNumInfo->GetIntNumber() &&
		 (((g2->xtals[crystalNum].crystalNum-1)%4)+1) == posNumInfo->GetIntNumber() &&
		 g2->xtals[crystalNum].cc >= eLowLimit->GetIntNumber() && 
		 g2->xtals[crystalNum].cc <= eHighLimit->GetIntNumber()) {
	       gotTheNextOne = 1;
	     } else {
	       crystalNum--;
	       gotTheNextOne = 0;
	     }
	   }
	   if (crystalNum < 0 && !gotTheNextOne) {
	     eventNum--;  treeIN->GetEvent(eventNum);
	   }
	 }
       }
     }

     UpdateEventDisplay();
     break;
   } 
   case 2: {
     if (fixCrystal->IsOn()) {
       selectSpecificCrystal = 1;
     } else {
       selectSpecificCrystal = 0;
     }
     break;
   }
   } // end switch(id)
}

void Viewer::UpdateEventDisplay() {
  evNumInfo->SetText(Form("%d", eventNum));
  cryNumInfo->SetText(Form("%d/%d", crystalNum+1, crystalsInEvent));

  qNumInfo->SetIntNumber(g2->xtals[crystalNum].quadNum);
  posNumInfo->SetIntNumber(((g2->xtals[crystalNum].crystalNum-1)%4) + 1);
  
  ccInfo->SetText(Form("%0.2f", g2->xtals[crystalNum].cc));
  errorInfo->SetText(Form("0x%X", g2->xtals[crystalNum].error));

  /* Define histograms, we have an event! */
  totalPlot = new TH1F("totalPlot", "totalPlot", 
		       g2->xtals[crystalNum].waveAll.size(), 
		       0, g2->xtals[crystalNum].waveAll.size());
  for (Int_t i=0; i<36; i++) {
    if (i<4) {
      ccPlot[i] = new TH1F(Form("CC%02d", i), 
			   Form("CC%02d", i), 
			   g2->xtals[crystalNum].waveAll.size()/40, 
			   0, g2->xtals[crystalNum].waveAll.size()/40);
    }
    segmentPlot[i] = new TH1F(Form("Segment%02d", i), 
			      Form("Segment%02d", i), 
			      g2->xtals[crystalNum].waveAll.size()/40, 
			      0, g2->xtals[crystalNum].waveAll.size()/40);
    segmentPlot[i]->SetLineColor(kBlue-2);
  }
  
  for (UInt_t i=0; i<g2->xtals[crystalNum].waveAll.size(); i++) {
    for (UInt_t j=0; j<40; j++) {
      
      Int_t nCrys = (g2->xtals[crystalNum].crystalNum-1)%4;
      Int_t segmentNumber = 0;	
      
      if (g2->xtals[crystalNum].quadNum==1) {
	segmentNumber = Q1Special[nCrys][j];
      } else if (g2->xtals[crystalNum].quadNum==2) {
	segmentNumber = Q2Special[nCrys][j];
      } else {
	segmentNumber = QNormal[nCrys][j];
      }
      if (segmentNumber >= 36 && (i >= ((j)*(g2->xtals[crystalNum].waveAll.size()/40))&&
				  i<((j+1)*(g2->xtals[crystalNum].waveAll.size()/40)))) {
	ccPlot[segmentNumber-36]->Fill(i-((j)*(g2->xtals[crystalNum].waveAll.size()/40)),
				       g2->xtals[crystalNum].waveAll[i]);
      } else if (i >= ((j)*(g2->xtals[crystalNum].waveAll.size()/40))&&
		 i<((j+1)*(g2->xtals[crystalNum].waveAll.size()/40))) {
	segmentPlot[segmentNumber]->Fill(i-((j)*(g2->xtals[crystalNum].waveAll.size()/40)),
					 g2->xtals[crystalNum].waveAll[i]);
      }
    }
    totalPlot->Fill(i, g2->xtals[crystalNum].waveAll[i]);
  }
  
  if (g2->xtals[crystalNum].crystalNum%2==0) {
    xyCanvas->cd();
    xyCanvas->Range(-40,-40,40,40);
    Int_t n=4;
    Float_t xg[4] = {0, 26.48, 34.59, 25.22};
    Float_t yg[4] = {0, 14.31, 0.67, -14.43};
    TPolyLine *segment1 = new TPolyLine(n, xg, yg, "F");
    segment1->SetLineColor(2);
    segment1->SetFillColor(2);
    segment1->SetFillStyle(0);
    segment1->Draw();
    Float_t xg2[4] = {0, 25.22, 15.85, -1.24};
    Float_t yg2[4] = {0, -14.43, -29.52, -27.52};
    TPolyLine *segment2 = new TPolyLine(n, xg2, yg2, "F");
    segment2->SetLineColor(3);
    segment2->SetFillColor(3);
    segment2->SetFillStyle(0);
    segment2->Draw();
    Float_t xg3[4] = {0, -1.24, -18.33, -24.24};
    Float_t yg3[4] = {0, -27.52, -25.53, -11.83};
    TPolyLine *segment3 = new TPolyLine(n, xg3, yg3, "F");
    segment3->SetLineColor(4);
    segment3->SetFillColor(4);
    segment3->SetFillStyle(0);
    segment3->Draw();
    Float_t xg4[4] = {0, -24.24, -30.15, -22.6};
    Float_t yg4[4] = {0, -11.83, 1.87, 14.03};
    TPolyLine *segment4 = new TPolyLine(n, xg4, yg4, "F");
    segment4->SetLineColor(5);
    segment4->SetFillColor(5);
    segment4->SetFillStyle(0);
    segment4->Draw();
    Float_t xg5[4] = {0, -22.6, -15.05, 1.67};
    Float_t yg5[4] = {0, 14.03, 26.19, 27.07};
    TPolyLine *segment5 = new TPolyLine(n, xg5, yg5, "F");
    segment5->SetLineColor(6);
    segment5->SetFillColor(6);
    segment5->SetFillStyle(0);
    segment5->Draw();
    Float_t xg6[4] = {0, 1.67, 18.38, 26.48};
    Float_t yg6[4] = {0, 27.07, 27.95, 14.31};
    TPolyLine *segment6 = new TPolyLine(n, xg6, yg6, "F");
    segment6->SetLineColor(7);
    segment6->SetFillColor(7);
    segment6->SetFillStyle(0);
    segment6->Draw();
    TArc *backArc = new TArc(0,0,40,0,360);
    backArc->SetLineColor(2);
    backArc->SetFillStyle(0);
    backArc->Draw();
    
    rzCanvas->cd();
    rzCanvas->Range(-3,-3,93,43);
    Int_t n2=5;
    Float_t r[5] = {0, 12, 18, 10, 0};
    Float_t z[5] = {0, 0, 6, 36.15, 34.59};
    TPolyLine *segment7 = new TPolyLine(n2, r, z, "F");
    segment7->SetLineColor(2);
    segment7->SetFillColor(2);
    segment7->SetFillStyle(0);
    segment7->Draw();
    Float_t r2[4] = {18, 25, 24, 10};
    Float_t z2[4] = {6, 6, 38.2, 36.15};
    TPolyLine *segment8 = new TPolyLine(n, r2, z2, "F");
    segment8->SetLineColor(3);
    segment8->SetFillColor(3);
    segment8->SetFillStyle(0);
    segment8->Draw();
    Float_t r3[4] = {25, 39, 38, 24};
    Float_t z3[4] = {6, 6, 40, 38.2};
    TPolyLine *segment9 = new TPolyLine(n, r3, z3, "F");
    segment9->SetLineColor(4);
    segment9->SetFillColor(4);
    segment9->SetFillStyle(0);
    segment9->Draw();
    Float_t r4[4] = {39, 38, 56, 56};
    Float_t z4[4] = {6, 40, 40, 6};
    TPolyLine *segment10 = new TPolyLine(n, r4, z4, "F");
    segment10->SetLineColor(5);
    segment10->SetFillColor(5);
    segment10->SetFillStyle(0);
    segment10->Draw(); 
    Float_t r5[4] = {56, 56, 76, 76};
    Float_t z5[4] = {40, 6, 6, 40};
    TPolyLine *segment11 = new TPolyLine(n, r5, z5, "F");
    segment11->SetLineColor(6);
    segment11->SetFillColor(6);
    segment11->SetFillStyle(0);
    segment11->Draw();      
    Float_t r6[4] = {76, 76, 90, 90};
    Float_t z6[4] = {6, 40, 40, 6};
    TPolyLine *segment12 = new TPolyLine(n, r6, z6, "F");
    segment12->SetLineColor(7);
    segment12->SetFillColor(7);
    segment12->SetFillStyle(0);
    segment12->Draw();
    
  } else if (g2->xtals[crystalNum].crystalNum%2==1) {
    xyCanvas->cd();
    xyCanvas->Range(-40,-40,40,40);
    Int_t n=4;
    Float_t xg[4] = {0, 30.3325, 30.3325, 16.7878};
    Float_t yg[4] = {0, -0.0242, -16.6962, -24.984};
    TPolyLine *segment1 = new TPolyLine(n, xg, yg, "F");
    segment1->SetLineColor(2);
    segment1->SetFillColor(2);
    segment1->SetFillStyle(0);
    segment1->Draw();
    Float_t xg2[4] = {0, 16.7878, 3.243, -12.0297};
    Float_t yg2[4] = {0, -24.984, -33.272, -24.76};
    TPolyLine *segment2 = new TPolyLine(n, xg2, yg2, "F");
    segment2->SetLineColor(3);
    segment2->SetFillColor(3);
    segment2->SetFillStyle(0);
    segment2->Draw();
    Float_t xg3[4] = {0, -12.0297, -27.302, -28.4815};
    Float_t yg3[4] = {0, -24.76, -16.2481, -0.8425};
    TPolyLine *segment3 = new TPolyLine(n, xg3, yg3, "F");
    segment3->SetLineColor(4);
    segment3->SetFillColor(4);
    segment3->SetFillStyle(0);
    segment3->Draw();
    Float_t xg4[4] = {0, -28.4815, -29.661, -15.2858};
    Float_t yg4[4] = {0, -0.8425, 14.5631, 24.0134};
    TPolyLine *segment4 = new TPolyLine(n, xg4, yg4, "F");
    segment4->SetLineColor(5);
    segment4->SetFillColor(5);
    segment4->SetFillStyle(0);
    segment4->Draw();
    Float_t xg5[4] = {0, -15.2858, -0.9109, 14.7108};
    Float_t yg5[4] = {0, 24.0134, 33.4637, 25.0558};
    TPolyLine *segment5 = new TPolyLine(n, xg5, yg5, "F");
    segment5->SetLineColor(6);
    segment5->SetFillColor(6);
    segment5->SetFillStyle(0);
    segment5->Draw();
    Float_t xg6[4] = {0, 14.7108, 30.3325, 30.3325};
    Float_t yg6[4] = {0, 25.0558, 16.6479, -0.0242};
    TPolyLine *segment6 = new TPolyLine(n, xg6, yg6, "F");
    segment6->SetLineColor(7);
    segment6->SetFillColor(7);
    segment6->SetFillStyle(0);
    segment6->Draw();
    TArc *backArc = new TArc(0,0,40,0,360);
    backArc->SetLineColor(2);
    backArc->SetFillStyle(0);
    backArc->Draw();
    
    rzCanvas->cd();
    rzCanvas->Range(-3,-3,93,43);
    Int_t n2=5;
    Float_t r[5] = {0, 12, 18, 8, 0};
    Float_t z[5] = {0, 0, 6, 36.184, 34.62};
    TPolyLine *segment7 = new TPolyLine(n2, r, z, "F");
    segment7->SetLineColor(2);
    segment7->SetFillColor(2);
    segment7->SetFillStyle(0);
    segment7->Draw();
    Float_t r2[4] = {18, 25, 22, 8};
    Float_t z2[4] = {6, 6, 38.9215, 36.184};
    TPolyLine *segment8 = new TPolyLine(n, r2, z2, "F");
    segment8->SetLineColor(3);
    segment8->SetFillColor(3);
    segment8->SetFillStyle(0);
    segment8->Draw();
    Float_t r3[5] = {25, 22, 28, 38, 39};
    Float_t z3[5] = {6, 38.9215, 40, 40, 6};
    TPolyLine *segment9 = new TPolyLine(n2, r3, z3, "F");
    segment9->SetLineColor(4);
    segment9->SetFillColor(4);
    segment9->SetFillStyle(0);
    segment9->Draw();
    Float_t r4[4] = {38, 39, 55, 55};
    Float_t z4[4] = {40, 6, 6, 40};
    TPolyLine *segment10 = new TPolyLine(n, r4, z4, "F");
    segment10->SetLineColor(5);
    segment10->SetFillColor(5);
    segment10->SetFillStyle(0);
    segment10->Draw(); 
    Float_t r5[4] = {55, 55, 76, 76};
    Float_t z5[4] = {6, 40, 40, 6};
    TPolyLine *segment11 = new TPolyLine(n, r5, z5, "F");
    segment11->SetLineColor(6);
    segment11->SetFillColor(6);
    segment11->SetFillStyle(0);
    segment11->Draw();      
    Float_t r6[4] = {76, 76, 90, 90};
    Float_t z6[4] = {40, 6, 6, 40};
    TPolyLine *segment12 = new TPolyLine(n, r6, z6, "F");
    segment12->SetLineColor(7);
    segment12->SetFillColor(7);
    segment12->SetFillStyle(0);
    segment12->Draw();
  }
  
  Float_t xPlot[16] = {0};
  Float_t yPlot[16] = {0};
  Float_t rPlot[16] = {0};
  Float_t zPlot[16] = {0};
  Int_t nPoints = 0;
  
  for (Int_t i=0; i<MAX_INTPTS; i++) { /* Make 'hit' segments red... */
    if (g2->xtals[crystalNum].intpts[i].e != 0 && (i<15 || (i==15 && g2->xtals[crystalNum].intpts[i].e>0)) ) {
      segmentPlot[(Int_t)g2->xtals[crystalNum].intpts[i].segNum]->SetLineColor(kRed);
      xPlot[nPoints] = g2->xtals[crystalNum].intpts[i].xyz.X();
      yPlot[nPoints] = g2->xtals[crystalNum].intpts[i].xyz.Y();
      zPlot[nPoints] = g2->xtals[crystalNum].intpts[i].xyz.Z();
      rPlot[nPoints] = TMath::Sqrt((g2->xtals[crystalNum].intpts[i].xyz.X()*g2->xtals[crystalNum].intpts[i].xyz.X()) + (g2->xtals[crystalNum].intpts[i].xyz.Y()*g2->xtals[crystalNum].intpts[i].xyz.Y()));
      nPoints++;	
      intptsSeg[i]->SetText(Form("%d", (Int_t)g2->xtals[crystalNum].intpts[i].segNum));
      intptsX[i]->SetText(Form("%0.1f", g2->xtals[crystalNum].intpts[i].xyz.X()));
      intptsY[i]->SetText(Form("%0.1f", g2->xtals[crystalNum].intpts[i].xyz.Y()));
      intptsZ[i]->SetText(Form("%0.1f", g2->xtals[crystalNum].intpts[i].xyz.Z()));
      intptsR[i]->SetText(Form("%0.1f", rPlot[nPoints-1]));
      intptsE[i]->SetText(Form("%0.1f", g2->xtals[crystalNum].intpts[i].e));
    } else {
      intptsSeg[i]->SetText("");
      intptsX[i]->SetText("");
      intptsY[i]->SetText("");
      intptsZ[i]->SetText("");
      intptsR[i]->SetText("");
      intptsE[i]->SetText("");
    }
  }
  for (Int_t i=1; i<=36; i++) {
    segmentCanvas->cd(i);
    segmentPlot[i-1]->Draw();
  }
  for (Int_t i=1; i<=4; i++) {
    ccCanvas->cd(i);
    ccPlot[i-1]->Draw();
  }
  xyCanvas->cd();
  if (nPoints > 0) {
    TGraph *xyGraph = new TGraph(nPoints, xPlot, yPlot);
    xyGraph->Draw("* SAME");
  }
  rzCanvas->cd();
  if (nPoints > 0) {
    TGraph *rzGraph = new TGraph(nPoints, zPlot, rPlot);
    rzGraph->Draw("* SAME"); 
  }
  crystalNum++;
  segmentCanvas->Modified();
  segmentCanvas->Update();
  ccCanvas->Modified();
  ccCanvas->Update();
  xyCanvas->Modified();
  xyCanvas->Update();
  rzCanvas->Modified();
  rzCanvas->Update();
}

void Viewer::DoExpand(Int_t event, Int_t px, Int_t py, TObject *sel) {
  TCanvas *c = (TCanvas *) gTQSender;
  TPad *pad = (TPad *) c->GetSelectedPad();

  if (!pad) return;
  if (event == kButton1Double) {
    pad->Pop();
    pad->cd();
    TCanvas *blownUp;
    TIter next(pad->GetListOfPrimitives());
    if ((TCanvas*)gROOT->GetListOfCanvases()->FindObject("blownUp")) {
      blownUp = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("blownUp");
      blownUp->Show();
      blownUp->Clear();
    } else {
      blownUp = new TCanvas("blownUp", "Expanded", 750, 0, 700, 500);
      blownUp->SetFillColor(10);
      blownUp->ToggleEventStatus();
      blownUp->Draw();
    }
    blownUp->cd();

    gROOT->SetSelectedPad((TPad*)blownUp);
    TPad *clone = (TPad*)pad->Clone();
    clone->SetPad(0,0,1,1);
    clone->Draw();
    clone->Modified();
    clone->Update();
  }
}

//---- Main program ------------------------------------------------------------

int main(int argc, char **argv)
{
  TString fileName;
  if (argc<2) {
    cout << "You need to specify the ROOT file too!" << endl;
    cout << "Usage: " << argv[0] << " <ROOT file with path> " << endl;
    exit(0);
  } else {
    fileName = argv[1];
    cout << "Using ROOT file " << fileName << endl;
  }
  
  TApplication theApp("App", &argc, argv);

   if (gROOT->IsBatch()) {
      fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
      return 1;
   }

   Viewer viewer(gClient->GetRoot(), fileName);

   gStyle->SetOptStat(0);
   
   theApp.Run();

   return 0;
}
