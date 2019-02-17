void plotCoincMatrix(){

  //*** Root File Readout ***//
  TString dataFolder ="/global/data1b/gretina/1730_data1b/";
  TChain *teb = new TChain("teb");
  
  //for(Int_t i=22; i<=34; i++){
  int i = 33;
    //if(i>=22 && i<=31) continue;
    TString dataFile = dataFolder + Form("Run00%i/Run00%i.root",i,i);
    std::cout << "Add File ::: " << dataFile << endl; 
    teb->Add(dataFile);
  //}

  //*** Draw ***//
  teb->SetAlias("Qval","12-si.pECal");
  
  TCanvas *cOnline = new TCanvas("cOnline", "Online Plots", 1800,1000);
  cOnline->Divide(2,2);
  cOnline->cd(1);  
  teb->Draw("g2.cc>>hst(4000,0,10000)","");
  cOnline->cd(2);
  gPad->SetLogz();
  teb->Draw("god.eventRaw.values:god.eventRaw.channels>>sichan_vals(600,0,600,1000,0,4000)","","colz");
  cOnline->cd(3);
  teb->Draw("12-si.pECal:g2.cc>>(7000,0,14000,500,0,14)", "(siType==\"SuperX3\") && si.upstream==1", "colz");
  cOnline->cd(4);
  teb->Draw("g2.cc>>hst_gated(4000,0,10000)","g2.cc[0] > 842 && g2.cc[0] < 852 && Qval > 5 && Qval < 6");
}
