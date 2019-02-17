void plotSi(){

  //*** Root File Readout ***//
  TString dataFolder ="/global/data1b/gretina/1730_data1b/";
  TChain *teb = new TChain("teb");
  Int_t iRunBegin = 89;
  Int_t iRunEnd = 91;


  for(Int_t iFile=iRunBegin; iFile<=iRunEnd; iFile++){
    TString dataFile = dataFolder + Form("Run00%i/Run00%i.root",iFile,iFile);
    std::cout << "Add File ::: " << dataFile << endl; 
    teb->Add(dataFile);
  }

  //*** Draw ***//
  TCanvas *cSiTarget = new TCanvas("cSiTarget", "Plots for Si target runs", 1200,600);
  cSiTarget->Divide(2,1);
  cSiTarget->cd(1);  
  teb->Draw("g2.cc>>hst(5000,0,10000)");
  cSiTarget->cd(2);
  teb->Draw("si.pECal>>hSi(3500,0,14)");
 
}
