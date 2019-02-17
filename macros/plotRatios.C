////////////////////////////////////////////////////////////////////////////////
//*** Created 2019/02/13 for GODDESS campaign - 56Fe (p,p') 1730 Jones exp ***//
//**************** H. Crawford, M. Jones, C. Santamaria **********************//
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//***      Run with the command:         > root macros/plotGammas.C        ***//
////////////////////////////////////////////////////////////////////////////////


double integrate(TH1F* h, double xmin, double xmax){
  TAxis *axis = h->GetXaxis();
  int bmin = axis->FindBin(xmin); //in your case xmin=-1.5
  int bmax = axis->FindBin(xmax); //in your case xmax=0.8
  double integral = h->Integral(bmin,bmax);
  integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
  integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/axis->GetBinWidth(bmax);

  return integral;

}

vector<double> projectQ(TH2F* hst, double Q, double Elevel){
  // Q units in [MeV]

  double Qlo = Q - 0.5;
  double Qhi = Q + 0.5;

  int qlo = hst->GetYaxis()->FindBin(Qlo);
  int qhi = hst->GetYaxis()->FindBin(Qhi);
  

  TH1F* hst_proj = (TH1F*)hst->ProjectionX("_px", qlo, qhi, "");

  // primary gamma energy
  double Egam_lo = Qlo - Elevel;
  double Egam_hi = Qhi - Elevel;

  Egam_lo *= 1e3;
  Egam_hi *= 1e3;

  hst_proj->Draw();
  TLine *line0 = new TLine(Egam_lo, 1, Egam_lo, 100);
  TLine *line1 = new TLine(Egam_hi, 1, Egam_hi, 100);
  line0->SetLineColor(2);
  line1->SetLineColor(2);
  line0->Draw("SAME");
  line1->Draw("SAME");

  double counts = integrate(hst_proj, Egam_lo, Egam_hi);

  //cout << "Integrating from: " << Egam_lo << " to: " << Egam_hi << endl;
  //cout << "Feeding to " << Elevel << " is: " << counts << " counts." << endl;
  
  double Eavg = 0.5*(Egam_hi + Egam_lo);// really this is the efficiency-corrected average of the primary spectrum

  vector<double> vec;

  vec.push_back(counts);
  vec.push_back(Eavg);

  return vec;
  

}




void plotRatios(){

  //*** Parameters ***//
  Int_t iRunBegin = 16; //First good production run with 10 MeV range in GRETINA
  Int_t iRunEnd = 18;//93
  Double_t pMinGate = 0.0, pMaxGate = 12.0;
  
  TH2F *hcc_847gate = new TH2F("hcc_847gate","p-Gamma with 842-852 keV gate", 4000,0,10000, 1000,0,14);
  TH2F *hcc_1810gate = new TH2F("hcc_1810gate","p-Gamma with 1800-1820 keV gate", 4000,0,10000, 1000,0,14);
  TH2F *hcc_2657gate = new TH2F("hcc_2657gate","p-Gamma with 2600-2700 keV gate", 4000,0,10000, 1000,0,14);

  TH2F *hcc_847bgate = new TH2F("hcc_847bgate","p-Gamma with 842-852 keV Bgate", 4000,0,10000, 1000,0,14);
  TH2F *hcc_1810bgate = new TH2F("hcc_1810bgate","p-Gamma with 1800-1820 keV Bgate", 4000,0,10000, 1000,0,14);
  TH2F *hcc_2657bgate = new TH2F("hcc_2657bgate","p-Gamma with 2600-2700 keV Bgate", 4000,0,10000, 1000,0,14);

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
    Bool_t InG1Gate = false, InG2Gate = false, InG3Gate = false, InPGate = false;
    Int_t Gamma1Gate= 0, Gamma2Gate = 0, Gamma3Gate = 0;

    // background contributions
    Bool_t InB1Gate = false, InB2Gate = false, InB3Gate = false;
    Int_t BGamma1Gate= 0, BGamma2Gate = 0, BGamma3Gate = 0;


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
      if(Gamma1Gate== false && g2->xtals[j].cc>2110. && g2->xtals[j].cc<2120.){// 842 -- 852
        InG1Gate = true;
        Gamma1Gate = j;
      }
      if(Gamma2Gate==false && g2->xtals[j].cc>1800. && g2->xtals[j].cc<1820.){
        InG2Gate = true;
        Gamma2Gate = j;
      }
      if(Gamma3Gate==false && g2->xtals[j].cc>2600. && g2->xtals[j].cc<2700.){
        InG3Gate = true;
        Gamma3Gate = j;
      }
      // background gates
      if(BGamma1Gate== false && g2->xtals[j].cc>860. && g2->xtals[j].cc<870.){
        InB1Gate = true;
        BGamma1Gate = j;
      }
      if(BGamma2Gate==false && g2->xtals[j].cc>1820. && g2->xtals[j].cc<1840.){
        InB2Gate = true;
        BGamma2Gate = j;
      }
      if(BGamma3Gate==false && g2->xtals[j].cc>2700. && g2->xtals[j].cc<2800.){
        InB3Gate = true;
        BGamma3Gate = j;
      }
      ECalorimeter += g2->xtals[j].cc;

    }

    for(Int_t iiSi = 0; iiSi<SiMult; iiSi++){
        if(eventOut->siType[iiSi]=="SuperX3" && eventOut->si[iiSi].upstream==1) hCalorimeter_Si->Fill(ECalorimeter, (12 - eventOut->si[iiSi].pECal));

	if(InG1Gate==true || InG2Gate==true || InG3Gate==true){
	for(Int_t k=0; k<gMult; k++) {
	  if(k!=Gamma1Gate && k!=Gamma2Gate && k!=Gamma3Gate){
	    hcc_847gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	    hcc_1810gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	    hcc_2657gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k!=Gamma1Gate && k==Gamma2Gate && k==Gamma3Gate){
	    hcc_847gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k==Gamma1Gate && k!=Gamma2Gate && k==Gamma3Gate){
	    hcc_1810gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k==Gamma1Gate && k==Gamma2Gate && k!=Gamma3Gate){
	    hcc_2657gate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }

	  ////
	  if(InB1Gate==true || InB2Gate==true || InB3Gate==true){
	  if(k!=BGamma1Gate && k!=BGamma2Gate && k!=BGamma3Gate){
	    hcc_847bgate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	    hcc_1810bgate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	    hcc_2657bgate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k!=BGamma1Gate && k==BGamma2Gate && k==BGamma3Gate){
	    hcc_847bgate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k==BGamma1Gate && k!=BGamma2Gate && k==BGamma3Gate){
	    hcc_1810bgate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	  }
	  else if(k==BGamma1Gate && k==BGamma2Gate && k!=BGamma3Gate){
	    hcc_2657bgate->Fill(g2->xtals[k].cc,(12 - eventOut->si[iiSi].pECal));
	    }
	  }

	}
      }
    }
  } //End loop over entries


  //*** Draw ***// 
  TCanvas *cG = new TCanvas("cG", "Gamma Plots", 1800,600);
  cG->Divide(3,2);
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
  cG->cd(5);
  hcc_2657gate->Draw("colz");

  
  // Form ratios
  // background is estimated by gating above + below the photopeak, then subtracting the superposition of the primary spectra in proportion 
  // to the compton background underneath the photopeak
  double Qarr[8] = {4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5};
  
  double N847[8] = {0};
  double Nb847[8] = {0};
  double E847[8] = {0};
  
  double N1810[8] = {0};// branching to 2657 needs to be included
  double Nb1810[8] = {0};
  double E1810[8] = {0};
  
  double N2657[8] = {0};
  double Nb2657[8] = {0};
  double E2657[8] = {0};
 

  double R1[8] = {0};

  vector<double> v847, vb847;
  vector<double> v1810, vb1810;
  vector<double> v2657, vb2657;

  for(int qi=0;qi<8;qi++){
    // v847 = projectQ(hcc_847gate, Qarr[qi], 0.847);
    v847 = projectQ(hcc_847gate, Qarr[qi], 2.959);
    N847[qi] = v847[0];
    E847[qi] = v847[1];
    v847.clear();

    v1810 = projectQ(hcc_1810gate, Qarr[qi], 2.657);
    N1810[qi] = v1810[0];
    E1810[qi] = v1810[1];
    v1810.clear();

    v2657 = projectQ(hcc_2657gate, Qarr[qi], 2.657);
    N2657[qi] = v2657[0];
    E2657[qi] = v2657[1];
    v2657.clear();

    //_________________________
    vb847 = projectQ(hcc_847bgate, Qarr[qi], 0.847);
    Nb847[qi] = vb847[0];
    vb847.clear();

    vb1810 = projectQ(hcc_1810bgate, Qarr[qi], 2.657);
    Nb1810[qi] = vb1810[0];
    vb1810.clear();

    vb2657 = projectQ(hcc_2657bgate, Qarr[qi], 2.657);
    Nb2657[qi] = vb2657[0];
    vb2657.clear();

    // we need efficiencies
    //R1[qi] = (N847[qi] - Nb847[qi]) / (N1810[qi] - Nb1810[qi] + N2657[qi] - Nb2657[qi]);
    R1[qi] = N847[qi]/N1810[qi];
    R1[qi] = R1[qi]*pow( E2657[qi] / E847[qi], 3);


  }

  cG->cd(6);
  TGraph *graph1 = new TGraph(6, Qarr, R1);
  graph1->Draw("ALP");

  //*** Save Histograms in root file ***//
  //TFile *fileOut = new TFile(Form("macros/Hists_Runs%i-%i_pGate%d-%d.root", iRunBegin, iRunEnd, pMinGate, pMaxGate), "RECREATE");
  //fileOut->Open();
  //hcc_Pgate->Write();
  //hcc_847gate->Write();
  //hCalorimeter_Si->Write();
  //hcc_1810gate->Write();
  //fileOut->Close();
}
