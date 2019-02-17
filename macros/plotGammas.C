////////////////////////////////////////////////////////////////////////////////
//*** Created 2019/02/13 for GODDESS campaign - 56Fe (p,p') 1730 Jones exp ***//
//**************** H. Crawford, M. Jones, C. Santamaria **********************//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//***      Run with the command:         > root macros/plotGammas.C        ***//
////////////////////////////////////////////////////////////////////////////////

void plotGammas(){

  //*** Parameters ***//
  Int_t iRunBegin = 16; //First good production run with 10 MeV range in GRETINA
  Int_t iRunEnd = 93;
  Double_t pMinGate = 0.0, pMaxGate = 12.0;
  
  TH2F *hcc_847gate = new TH2F("hcc_847gate","p-Gamma with 842-852 keV gate", 4000,0,10000, 1000,0,14);
  TH2F *hcc_1810gate = new TH2F("hcc_1810gate","p-Gamma with 1800-1820 keV gate", 4000,0,10000, 1000,0,14);
  TH1F *hcc_Pgate = new TH1F("hcc_Pgate",Form("Gamma with %i-%i MeV p gate",int(pMinGate),int(pMaxGate)), 4000,0,10000);
  TH2F *hCalorimeter_Si = new TH2F("hCalorimeter_Si","hCalorimeter_Si", 4000,0,12000, 1000,0,14);
  
 
  //*** RootFile Readout ***//
  TString dataFolder ="/global/data1b/gretina/1730_data1b/";
  TChain *teb = new TChain("teb");

  for(Int_t iFile=iRunBegin; iFile<=iRunEnd; iFile++){
    if((iFile>=22 && iFile<=31) || iFile==43 || iFile==44 || iFile==47 || (iFile>=53 && iFile<=55) || iFile==63 || iFile==72 || (iFile>=82 && iFile>=84) || (iFile>=88 && iFile<=91) || iFile==92) continue; // runs removed for now for analysis -- Run88=Empty Target -- Runs 89-91=Si target
  TString dataFile = dataFolder + Form("Run%04d/Run%04d.root",iFile,iFile);
  std::cout << "Add File ::: " << dataFile << endl; 
  teb->Add(dataFile);
  }

  g2OUT *g2 = new g2OUT();
  goddessOut *eventOut = new goddessOut();
  teb->SetBranchAddress("g2",&g2);

  teb->SetBranchAddress("godEv",&eventOut);
  Long64_t nentries = teb->GetEntries();

  //*** Loop over entries ***//
  for(Long64_t i=0; i<nentries; i++){
    
    if( i%1000000==0) cout << i*100./nentries << " % events processed" << endl;
    teb->GetEntry(i);
    UInt_t gMult = g2->crystalMult();
    Bool_t InG1Gate = false, InG2Gate = false, InPGate = false;
    Int_t Gamma1Gate= 0, Gamma2Gate = 0;
    Int_t SiMult = eventOut->siType.size();

    for(int iSi = 0; iSi<SiMult; iSi++){
      if(eventOut->siType[iSi]=="SuperX3" && eventOut->si[iSi].upstream==1){
        if((12 - eventOut->si[iSi].pECal)>pMinGate && (12 - eventOut->si[iSi].pECal)<pMaxGate){
          InPGate = true;
          for(Int_t kk=0; kk<gMult; kk++) hcc_Pgate->Fill(g2->xtals[kk].cc);
	  }
      }
    }

    Double_t ECalorimeter = 0.0;

    for(UInt_t j=0; j<gMult; j++) {
      if(Gamma1Gate== false && g2->xtals[j].cc>842. && g2->xtals[j].cc<852.){
        InG1Gate = true;
        Gamma1Gate = j;
      }
      if(Gamma2Gate==false && g2->xtals[j].cc>1800. && g2->xtals[j].cc<1820.){
        InG2Gate = true;
        Gamma2Gate = j;
      }
      ECalorimeter += g2->xtals[j].cc;

    }

    for(Int_t iiSi = 0; iiSi<SiMult; iiSi++){
        if(eventOut->siType[iiSi]=="SuperX3" && eventOut->si[iiSi].upstream==1) hCalorimeter_Si->Fill(ECalorimeter, (12 - eventOut->si[iiSi].pECal));

      if(InG1Gate==true || InG2Gate==true){
	for(Int_t k=0; k<gMult; k++) {
	  if(k!=Gamma1Gate && k!=Gamma2Gate){
	    hcc_847gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	    hcc_1810gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k!=Gamma1Gate && k==Gamma2Gate){
	    hcc_847gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k==Gamma1Gate && k!=Gamma2Gate){
	    hcc_1810gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	}
      }
    }
  } //End loop over entries


  //*** Draw ***// 
  TCanvas *cG = new TCanvas("cG", "Gamma Plots", 1800,600);
  cG->Divide(2,2);
  cG->cd(1);
  hcc_Pgate->Draw();
  cG->cd(2);
  gPad->SetLogz();
  hcc_847gate->Draw("colz");
  cG->cd(3);
  gPad->SetLogz();
  hCalorimeter_Si->Draw("colz");
  cG->cd(4);
  gPad->SetLogz();
  hcc_1810gate->Draw("colz");


  //*** Save Histograms in root file ***//
  //TFile *fileOut = new TFile(Form("macros/Hists_Runs%i-%i_pGate%d-%d.root", iRunBegin, iRunEnd, pMinGate, pMaxGate), "RECREATE");
  //fileOut->Open();
  //hcc_Pgate->Write();
  //hcc_847gate->Write();
  //hCalorimeter_Si->Write();
  //hcc_1810gate->Write();
  //fileOut->Close();
}
