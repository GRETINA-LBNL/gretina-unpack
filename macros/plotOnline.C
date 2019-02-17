void plotOnline(Int_t iRun=16){

  //*** Root File Readout ***//
  TString dataFile = Form("/global/data1b/gretina/1730_data1b/Run%04d/Run%04d.root",iRun, iRun);
  std::cout << "Read File ::: " << dataFile << endl; 
  TChain *teb = new TChain("teb");
  teb->Add(dataFile);

  //*** Draw ***//
  TCanvas *cOnline = new TCanvas("cOnline", "Online Plots", 1200,600);
  cOnline->Divide(2,1);
  cOnline->cd(1);  
  teb->Draw("g2.cc>>hst(4000,0,10000)","");
  cOnline->cd(2);
  gPad->SetLogz();
  teb->Draw("god.eventRaw.values:god.eventRaw.channels>>sichan_vals(600,0,600,1000,0,4000)","","colz");
 
}


void plotOnline(){


}
