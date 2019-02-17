void draw1D(TTree *data, Int_t upstream) {

  Long64_t nEvents = data->GetEntries();
  printf("Sorting %ld events\n", nEvents);
  
  goddessOut *evOut = new goddessOut();
  data->SetBranchAddress("godEv", &evOut);

  TH1F *strip[48];
  for (Int_t st=0; st<48; st++) {
    strip[st] = new TH1F(Form("strip%d", st), Form("strip%d", st), 1000, 0, 6000);
  }
  
  for (Long64_t i=0; i<nEvents; i++) {
    data->GetEvent(i);
    if (i%10000==0) { printf("Done %f\r", 100.*(Float_t)i/(Float_t)nEvents); }
    Int_t det = 0;
    for (UInt_t j=0; j<evOut->si.size(); j++) {
      if (upstream) {
	if (evOut->siID[j]=="U0_E1") { det = 0; }
	else if (evOut->siID[j]=="U1_E1") { det = 1; }
	else if (evOut->siID[j]=="U2_E1") { det = 2; }
	else if (evOut->siID[j]=="U3_E1") { det = 3; }
	else if (evOut->siID[j]=="U4_E1") { det = 4; }
	else if (evOut->siID[j]=="U5_E1") { det = 5; }
	else if (evOut->siID[j]=="U6_E1") { det = 6; }
	else if (evOut->siID[j]=="U7_E1") { det = 7; }
	else if (evOut->siID[j]=="U8_E1") { det = 8; }
	else if (evOut->siID[j]=="U9_E1") { det = 9; }
	else if (evOut->siID[j]=="U10_E1") { det = 10; }
	else if (evOut->siID[j]=="U11_E1") { det = 11; }
      } else if (!upstream) {
	if (evOut->siID[j]=="D0_E1") { det = 0; }
	else if (evOut->siID[j]=="D1_E1") { det = 1; }
	else if (evOut->siID[j]=="D2_E1") { det = 2; }
	else if (evOut->siID[j]=="D3_E1") { det = 3; }
	else if (evOut->siID[j]=="D4_E1") { det = 4; }
	else if (evOut->siID[j]=="D5_E1") { det = 5; }
	else if (evOut->siID[j]=="D6_E1") { det = 6; }
	else if (evOut->siID[j]=="D7_E1") { det = 7; }
	else if (evOut->siID[j]=="D8_E1") { det = 8; }
	else if (evOut->siID[j]=="D9_E1") { det = 9; }
	else if (evOut->siID[j]=="D10_E1") { det = 10; }
	else if (evOut->siID[j]=="D11_E1") { det = 11; }
      }
      for (Int_t k=0; k<4; k++) {
	if (evOut->si[j].pStrip == k) {
	  strip[det*4 + k]->Fill(evOut->si[j].pECal);
	}
      }
    }
  }
  
  printf("Making canvases\n");
  TCanvas* m[12];
  for (Int_t p=0; p<12; p++) {
    m[p] = new TCanvas();
    m[p]->Divide(4,1);
    printf("Formed %d\n", p);
  }
  printf("Made the canvases\n");
  for (Int_t q=0; q<12; q++) {
    for (Int_t p=0; p<4; p++) {
      printf("%d %d %d\n", q, p+1, q*4+p);
      m[q]->cd(p+1);  strip[q*4 + p]->Draw();

    }
  }
}

void draw2D(TTree *data, Int_t upstream) {

  Long64_t nEvents = data->GetEntries();
  printf("Sorting %ld events\n", nEvents);
  
  goddessOut *evOut = new goddessOut();
  data->SetBranchAddress("godEv", &evOut);

  TH2F *strip[48];
  for (Int_t st=0; st<48; st++) {
    strip[st] = new TH2F(Form("strip%d", st), Form("strip%d", st), 1000, 0, 6000, 1000, 0, 6000);
  }
  
  for (Long64_t i=0; i<nEvents; i++) {
    data->GetEvent(i);
    if (i%10000==0) { printf("Done %f\r", 100.*(Float_t)i/(Float_t)nEvents); }
    Int_t det = -1;
    for (UInt_t j=0; j<evOut->si.size(); j++) {
      if (upstream) {
	if (evOut->siID[j]=="U0_E1") { det = 0; }
	else if (evOut->siID[j]=="U1_E1") { det = 1; }
	else if (evOut->siID[j]=="U2_E1") { det = 2; }
	else if (evOut->siID[j]=="U3_E1") { det = 3; }
	else if (evOut->siID[j]=="U4_E1") { det = 4; }
	else if (evOut->siID[j]=="U5_E1") { det = 5; }
	else if (evOut->siID[j]=="U6_E1") { det = 6; }
	else if (evOut->siID[j]=="U7_E1") { det = 7; }
	else if (evOut->siID[j]=="U8_E1") { det = 8; }
	else if (evOut->siID[j]=="U9_E1") { det = 9; }
	else if (evOut->siID[j]=="U10_E1") { det = 10; }
	else if (evOut->siID[j]=="U11_E1") { det = 11; }
      } else if (!upstream) {
	if (evOut->siID[j]=="D0_E1") { det = 0; }
	else if (evOut->siID[j]=="D1_E1") { det = 1; }
	else if (evOut->siID[j]=="D2_E1") { det = 2; }
	else if (evOut->siID[j]=="D3_E1") { det = 3; }
	else if (evOut->siID[j]=="D4_E1") { det = 4; }
	else if (evOut->siID[j]=="D5_E1") { det = 5; }
	else if (evOut->siID[j]=="D6_E1") { det = 6; }
	else if (evOut->siID[j]=="D7_E1") { det = 7; }
	else if (evOut->siID[j]=="D8_E1") { det = 8; }
	else if (evOut->siID[j]=="D9_E1") { det = 9; }
	else if (evOut->siID[j]=="D10_E1") { det = 10; }
	else if (evOut->siID[j]=="D11_E1") { det = 11; }
      }
      if (det >= 0) {
	for (Int_t k=0; k<4; k++) {
	  if (evOut->si[j].pStrip == k) {
	    strip[det*4 + k]->Fill(evOut->si[j].eNearCal, evOut->si[j].eFarCal);
	  }
	}
      }
    }
  }
  
  printf("Making canvases\n");
  TCanvas* m[12];
  for (Int_t p=0; p<12; p++) {
    m[p] = new TCanvas();
    m[p]->Divide(4,1);
    printf("Formed %d\n", p);
  }
  printf("Made the canvases\n");
  for (Int_t q=0; q<12; q++) {
    for (Int_t p=0; p<4; p++) {
      printf("%d %d %d\n", q, p+1, q*4+p);
      m[q]->cd(p+1);  strip[q*4 + p]->Draw("colz");
    }
  }
}
