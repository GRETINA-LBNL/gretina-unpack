#include "TFile.h"
#include "TTree.h"
#include "TMath.h"


#include "include/GRETINA.h"

void plotS2(TTree *dataT, Int_t nEvents) {

  Int_t segmentID[16] = {160, 161, 162, 163, 164, 165, 166, 167,
			 170, 171, 172, 173, 174, 175, 176, 177};
  Int_t ringID[48] = {220, 230, 221, 231, 222, 232, 223, 233,
		      224, 234, 225, 235, 226, 236, 227, 237,
		      180, 200, 181, 201, 182, 202, 183, 203,
		      184, 204, 185, 205, 186, 206, 187, 207,
		      190, 210, 191, 211, 192, 212, 193, 213,
		      194, 214, 195, 215, 196, 216, 197, 217};

  g3OUT *g3 = new g3OUT();
  dataT->SetBranchAddress("g3", &g3);

  Int_t nEntries = dataT->GetEntries();
  if (nEvents < nEntries) { nEntries = nEvents; }

  TH2D *mSegE = new TH2D("maxSegE", "maxSegE", 16, 0, 16, 6000, 0, 6000);
  TH2D *mRingE = new TH2D("maxRingE", "maxRingE", 48, 0, 48, 6000, 0, 6000);
  TH1D *mSeg = new TH1D("maxSeg", "maxSeg", 16, 0, 16);
  TH1D *mRing = new TH1D("maxRing", "maxRing", 48, 0, 48);
  
  //  TH2F *dist = new TH2F("map", "map", 360, 0, 6.28, 600, 0, 60);
  TCanvas *c1 = new TCanvas("c1", "S2 Map", 600, 600);
  TGraphPolar *gr = new TGraphPolar();

  
  for (Long64_t i=0; i<nEntries; i++) {
    dataT->GetEntry(i);

    Int_t maxSeg = -1; Int_t maxRing = -1; 
    Float_t maxSegE = -1; Float_t maxRingE = -1;
    
    for (Int_t i=0; i<g3->xtals.size(); i++) {
      for (Int_t j=0; j<g3->xtals[i].chn.size(); j++) {
	for (Int_t seg = 0; seg<16; seg++) {
	  if (g3->xtals[i].chn[j].ID == segmentID[seg]) {
	    if (g3->xtals[i].chn[j].eRaw > maxSegE) {
	      maxSeg = seg;
	      maxSegE = g3->xtals[i].chn[j].eRaw;
	    }
	    //printf("Segment %d - eRaw %f\n", seg, g3->xtals[i].chn[j].eRaw);
	  }
	}
	for (Int_t ring = 0; ring<48; ring++) {
	  if (g3->xtals[i].chn[j].ID == ringID[ring]) {
	    if (g3->xtals[i].chn[j].eRaw > maxRingE) {
	      maxRing = ring;
	      maxRingE = g3->xtals[i].chn[j].eRaw;
	    }
	    //printf("Ring %d - eRaw %f\n", ring, g3->xtals[i].chn[j].eRaw);
	  }
	}
      }
    }

    printf("\n");
    if (maxSeg >= 0 && maxRing >= 0) {
      // dist->Fill(gRandom->Uniform(22.5*maxSeg*0.017453, 0.017453*22.5*(maxSeg+1)), gRandom->Uniform(maxRing+2, maxRing+3));
      gr->AddPoint(1.5708+gRandom->Uniform(22.5*maxSeg*0.017453, 0.017453*22.5*(maxSeg+1)), gRandom->Uniform(maxRing+2, maxRing+3));
      mRingE->Fill(maxRing, maxRingE);  mRing->Fill(maxRing);
      mSegE->Fill(maxSeg, maxSegE);  mSeg->Fill(maxSeg);
    }
  }

  // dist->Draw("pol col");
  gr->Draw("AP");
  c1->Update();
  gr->GetPolargram()->SetNdivPolar(0);
  gr->GetPolargram()->SetNdivRadial(0);
  gr->GetPolargram()->SetLineColor(kWhite);
  gr->GetPolargram()->SetToDegree();
  gr->GetPolargram()->SetAxisAngle(3.14/2.);
  c1->Modified();
  c1->Update();

  TCanvas *c2 = new TCanvas("max", "max", 600, 600);
  c2->Divide(2,2);
  c2->cd(1); mSegE->Draw("colz");
  c2->cd(2); mRingE->Draw("colz");
  c2->cd(3); mSeg->Draw();
  c2->cd(4); mRing->Draw();
  
}
