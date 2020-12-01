#include "TCanvas.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"

#include <string>

#include "include/GRETINA.h"

void doAmScan1D(Int_t runNumStart, Int_t runNumStop, TString directory = "./", Int_t calibrated = 0, Int_t threshold = 100,  Int_t layer = 0) {

  Int_t files = runNumStop - runNumStart + 1;
  printf("Analysis includes %d files.  Starting now.\n", files);

  vector<Int_t> index;
  vector<Int_t> counts1;
  vector<Int_t> counts2;
  vector<Int_t> counts3; 
  vector<Int_t> counts4;
  vector<Int_t> counts5;
  vector<Int_t> counts6;
  vector<Int_t> xPos;
  vector<Int_t> yPos;
  vector<Int_t> zPos;
  
  for (Int_t fileNum = 0; fileNum < files; fileNum++) {
    TString filename = (directory + "Run0" + std::to_string(runNumStart+fileNum) + ".root");
    TFile *f = new TFile(filename.Data(), "READ");
    TTree *data = (TTree*)f->Get("teb");

    g3OUT *g3 = new g3OUT();
    data->SetBranchAddress("g3", &g3);
    Double_t collX, collY, collZ;
    data->SetBranchAddress("collX", &collX);
    data->SetBranchAddress("collY", &collY);
    data->SetBranchAddress("collZ", &collZ);

    Int_t counts[6] = {0};

    Long64_t nEntries = data->GetEntries();
    for (Int_t i=0; i<nEntries; i++) {
      data->GetEvent(i);
      if (i==0) { 
	index.push_back(fileNum);
	xPos.push_back(collX);
	yPos.push_back(collY);
	zPos.push_back(collZ);
      }
      
      for (Int_t xtal = 0; xtal<g3->xtals.size(); xtal++) {
	for (Int_t chn = 0; chn<g3->xtals[xtal].chn.size(); chn++) {
	  if (g3->xtals[xtal].chn[chn].segNum/6 == layer) {
	    if (!calibrated && g3->xtals[xtal].chn[chn].eRaw > threshold) {
	      counts[g3->xtals[xtal].chn[chn].segNum%6]++;
	    } else if (calibrated && g3->xtals[xtal].chn[chn].eCal > 58. && g3->xtals[xtal].chn[chn].eCal < 61.) {
	      counts[g3->xtals[xtal].chn[chn].segNum%6]++;
	    }
	  }
	}
      }
    } // loop over entries
    
    counts1.push_back(counts[0]);     counts2.push_back(counts[1]);     counts3.push_back(counts[2]);
    counts4.push_back(counts[3]);     counts5.push_back(counts[4]);     counts6.push_back(counts[5]);

    f->Close();
    f->Delete();
  
  }

  TGraph *gr1 = new TGraph(files, &(index[0]), &(counts1[0]));
  TGraph *gr2 = new TGraph(files, &(index[0]), &(counts2[0]));
  TGraph *gr3 = new TGraph(files, &(index[0]), &(counts3[0]));
  TGraph *gr4 = new TGraph(files, &(index[0]), &(counts4[0]));
  TGraph *gr5 = new TGraph(files, &(index[0]), &(counts5[0]));
  TGraph *gr6 = new TGraph(files, &(index[0]), &(counts6[0]));
  gr1->SetLineColor(1);  gr2->SetLineColor(2);  gr3->SetLineColor(3);
  gr4->SetLineColor(4);  gr5->SetLineColor(5);  gr6->SetLineColor(6);

  // Determine Y-scale for plot...
  vector<Int_t> max;
  max.push_back(*max_element(counts1.begin(), counts1.end()));
  max.push_back(*max_element(counts2.begin(), counts2.end()));
  max.push_back(*max_element(counts3.begin(), counts3.end()));
  max.push_back(*max_element(counts4.begin(), counts4.end()));
  max.push_back(*max_element(counts5.begin(), counts5.end()));
  max.push_back(*max_element(counts6.begin(), counts6.end()));
  gr1->SetMaximum(*max_element(max.begin(), max.end()));

  gr1->Draw("A*L");
  gr2->Draw("*Lsame");
  gr3->Draw("*Lsame");
  gr4->Draw("*Lsame");
  gr5->Draw("*Lsame");
  gr6->Draw("*Lsame");


}



void doAmScan2D(Int_t runNumStart, Int_t runNumStop, TString directory = "./", Int_t calibrated = 0, Int_t threshold = 100,  Int_t layer = 0) {

  Int_t files = runNumStop - runNumStart + 1;
  printf("Analysis includes %d files.  Starting now.\n", files);

  vector<Int_t> index;
  vector<Int_t> counts1;
  vector<Int_t> counts2;
  vector<Int_t> counts3; 
  vector<Int_t> counts4;
  vector<Int_t> counts5;
  vector<Int_t> counts6;
  vector<Int_t> xPos;
  vector<Int_t> yPos;
  vector<Int_t> zPos;
  
  for (Int_t fileNum = 0; fileNum < files; fileNum++) {
    TString filename = (directory + "Run0" + std::to_string(runNumStart+fileNum) + ".root");
    TFile *f = new TFile(filename.Data(), "READ");
    TTree *data = (TTree*)f->Get("teb");

    g3OUT *g3 = new g3OUT();
    data->SetBranchAddress("g3", &g3);
    Double_t collX, collY, collZ;
    data->SetBranchAddress("collX", &collX);
    data->SetBranchAddress("collY", &collY);
    data->SetBranchAddress("collZ", &collZ);

    Int_t counts[6] = {0};

    Long64_t nEntries = data->GetEntries();
    for (Int_t i=0; i<nEntries; i++) {
      data->GetEvent(i);
      if (i==0) { 
	index.push_back(fileNum);
	xPos.push_back(collX);
	yPos.push_back(collY);
	zPos.push_back(collZ);
	printf("File: %s; Index = %d; Collimator (x,y,z) = (%lf, %lf, %lf)\n", filename.Data(), fileNum, collX, collY, collZ);
      }
      
      for (Int_t xtal = 0; xtal<g3->xtals.size(); xtal++) {
	for (Int_t chn = 0; chn<g3->xtals[xtal].chn.size(); chn++) {
	  if (g3->xtals[xtal].chn[chn].segNum/6 == layer) {
	    if (!calibrated && g3->xtals[xtal].chn[chn].eRaw > threshold) {
	      counts[g3->xtals[xtal].chn[chn].segNum%6]++;
	    } else if (calibrated && g3->xtals[xtal].chn[chn].eCal > 58. && g3->xtals[xtal].chn[chn].eCal < 61.) {
	      counts[g3->xtals[xtal].chn[chn].segNum%6]++;
	    }
	  }
	}
      }
    } // loop over entries
    
    counts1.push_back(counts[0]);     counts2.push_back(counts[1]);     counts3.push_back(counts[2]);
    counts4.push_back(counts[3]);     counts5.push_back(counts[4]);     counts6.push_back(counts[5]);

    f->Close();
    f->Delete();
  
  }

  Double_t xmin = *min_element(xPos.begin(), xPos.end());
  Double_t xmax = *max_element(xPos.begin(), xPos.end());
  Double_t ymin = *min_element(yPos.begin(), yPos.end());
  Double_t ymax = *max_element(yPos.begin(), yPos.end());
  
  TH2D *hist[6];
  for (Int_t i=0; i<6; i++) {
    hist[i] = new TH2D(Form("hist%d", i+1), Form("hist%d", i+1), (Int_t) (xmax - xmin+1), xmin-1, xmax+1, 
		       (Int_t) (ymax-ymin+1), ymin-1, ymax+1);
  } 

  for (Int_t i=0; i<files; i++) {
    hist[0]->Fill(xPos[i], yPos[i], counts1[i]);
    hist[1]->Fill(xPos[i], yPos[i], counts2[i]);
    hist[2]->Fill(xPos[i], yPos[i], counts3[i]);
    hist[3]->Fill(xPos[i], yPos[i], counts4[i]);
    hist[4]->Fill(xPos[i], yPos[i], counts5[i]);
    hist[5]->Fill(xPos[i], yPos[i], counts6[i]);
  }

  for (Int_t i=0; i<6; i++) {
    hist[i]->SetLineColor(i+1);
  }

  hist[0]->Draw("BOX");
  for (Int_t i=1; i<6; i++) {
    hist[i]->Draw("BOX SAME");
  }


}
