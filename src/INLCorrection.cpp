/**********************************************************************/
/* File: INLCorrection.cpp                                            */
/* Description: Functions for GRETINA non-linearity correctiom        */
/* Author: H. Crawford                                                */
/* Date: March 2013                                                   */
/**********************************************************************/

#include "INLCorrection.h"

ClassImp(INLCorrection);

void INLCorrection::Initialize(controlVariables* ctrl, 
			       GRETINAVariables* gVar) { 
  
  for (Int_t crys=0; crys<MAXCRYSTALS; crys++) {
    for (Int_t bd=0; bd<4; bd++) {
      crystalSN[crys][bd] = 0;
    }
  }

  for (Int_t crys=0; crys<MAXCRYSTALS; crys++) {
    for (Int_t chn=0; chn<4; chn++) {
      for (Int_t j=0; j<8192; j++) {
	for (Int_t i=0; i<3; i++) {
	  enl[crys][chn][i][j] = 0;
	}
	inl[crys][chn][j] = 0;
      }
    }
  }

  for (Int_t chn=0; chn<MAXCHANNELS; chn++) {
    lastCalcBase[chn] = -10000;
    restingBase[chn] = 0.;
    Q[chn] = 0.;
    sigmaBase[chn] = 0.;
    baseSamples[chn] = 0;
  }

  if (ctrl->digMapFileName.Length() > 2) {
    ReadDigitizerMap(ctrl->digMapFileName);
    ReadRawINL();
  } else {
    /* No digitizer map filename provided; assume
       we aren't actually doing INL corrections, and let
       the user know... */
    cout << "No digitizer map filename provided; no INL correction. " << endl;
    ctrl->INLcorrection = 0;
  }
}

void INLCorrection::calculateENL(Int_t xtal, Int_t chn, Int_t s1_width,
				 Int_t s2_width, Float_t s1, Float_t s2,
				 Float_t base, Float_t tau) {
  
  Float_t ss1 = 0., ss2 = 0., ss3 = 0., ss4 = 0.;
  Float_t sum1 = 0., sum2 = 0., sum3 = 0., sum4 = 0.;

  Float_t e1 = 0.;
  Int_t i, j;

  /* Clear the ENL correction calculated previously. */
  for (i=0; i<3; i++) {
    for (j=0; j<8192; j++) {
      enl[xtal][chn][i][j] = 0.;
    }
  }

  for (i=0; i<8192; i++) {
    e1 = i*0.5 - 500;
    if ((chn%10) != 9) { e1 = -e1; }
    sum1 = sum2 = sum3 = sum4 = 0.;
    
    for (j=0; j<500; j++) {
      sum1 += e1 + base;
      sum2 += inl[xtal][chn][(Int_t)((8192.0 + e1 + base)/2.0)];
      
      if (j < s2_width) {
	sum3 += (e1 + base - inl[xtal][chn][(Int_t)((8192.0 + e1 + base)/2.0)]);
      }
      if (j >= (500 - s1_width)) {
	sum4 += (e1 + base - inl[xtal][chn][(Int_t)((8192.0 + e1 + base)/2.0)]);
      }
      e1 *= (1.0 - (tau));
    }
    
    if ((chn%10) != 9) {
      sum3 /= -s2_width;
      sum4 /= -s1_width;
    } else {
      sum3 /= s2_width;
      sum4 /= s1_width;
    }

    sum1 += 2000.0 - sum2;
    sum3 += 2000.0;
    sum4 += 2000.0;

    if ( i>0 && ss1>0 && ss1<sum1 && sum1<8190 ) {
      for (j = (1 + (Int_t)ss1); (e1 = j) <= sum1; j++) {
	enl[xtal][chn][0][j] = ss2 + (sum2 - ss2)*(e1 - ss1)/(sum1 - ss1);
      }
    }
    if ( i>0 && ss3>0 && ss3<sum3 && sum3<8190 ) { /* Index by beginning of region sum -- flattop */
      for (j = (1 + (Int_t)ss3); (e1 = j) <= sum3; j++) {
	enl[xtal][chn][1][j] = ss2 + (sum2 - ss2)*(e1 - ss3)/(sum3 - ss3);
      }
    }
    if ( i>0 && ss4>0 && ss4<sum4 && sum4<8190 ) { /* Index by end of region sum -- baseline */
      for (j = (1 + (Int_t)ss4); (e1 = j) <= sum4; j++) {
	enl[xtal][chn][2][j] = ss2 + (sum2 - ss2)*(e1 - ss4)/(sum4 - ss4);
      }
    }
    ss1 = sum1;
    ss2 = sum2;
    ss3 = sum3;
    ss4 = sum4;
  }

  for (i=0; i<8192; i++) {
    printf("%d %f %f\n", i, enl[xtal][chn][1][i], enl[xtal][chn][2][i]);
  }

}


void INLCorrection::ReadDigitizerMap(TString filename) {
  
  FILE *digMap;
  
  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<4; j++) {
      crystalSN[i][j] = 0;
    }
  }

  if ((digMap = fopen(filename.Data(), "r")) == NULL) {
    cerr << "Failed to open digitizer SN map: " << filename.Data() << endl;
  } else {
    cout << "Reading digitizer map \"" << filename.Data() << "\"" << endl;
  }

  Int_t dorder, dsn, i1, dhole;
  char *st, str[256];
  
  st = fgets(str, 200, digMap);
  while (st != NULL) {
    if (str[0] == 35) {
    } else if (str[0] == 59) {
    } else if (str[0] == 10) {
    } else {
      sscanf(str, "%i %i %X %X", &dorder, &i1, &dsn, &dhole);
      if (dorder%4 == 1) {
	if (i1 == 0) { crystalSN[(Int_t)(dorder/4)][0] = dsn; }
	if (i1 == 10) { crystalSN[(Int_t)(dorder/4)][0] = dsn | (2 << 10); }
	if (i1 == 01) { crystalSN[(Int_t)(dorder/4)][0] = dsn | (1 << 10); }
	if (i1 == 11) { crystalSN[(Int_t)(dorder/4)][0] = dsn | (3 << 10); }
      } else if (dorder%4 == 2) {
	if (i1 == 0) { crystalSN[(Int_t)(dorder/4)][1] = dsn; }
	if (i1 == 10) { crystalSN[(Int_t)(dorder/4)][1] = dsn | (2 << 10); }
	if (i1 == 11) { crystalSN[(Int_t)(dorder/4)][1] = dsn | (3 << 10); }
	if (i1 == 01) { crystalSN[(Int_t)(dorder/4)][1] = dsn | (1 << 10); }
      } else if (dorder%4 == 3) {
	if (i1 == 0) { crystalSN[(Int_t)(dorder/4)][2] = dsn; }
	if (i1 == 10) { crystalSN[(Int_t)(dorder/4)][2] = dsn | (2 << 10); }
	if (i1 == 11) { crystalSN[(Int_t)(dorder/4)][2] = dsn | (3 << 10); }
	if (i1 == 01) { crystalSN[(Int_t)(dorder/4)][2] = dsn | (1 << 10); }
      } else if (dorder%4 == 0) {
	if (i1 == 0) { crystalSN[(Int_t)(dorder/4)-1][3] = dsn; }
	if (i1 == 10) { crystalSN[(Int_t)(dorder/4)-1][3] = dsn | (2 << 10); }
	if (i1 == 11) { crystalSN[(Int_t)(dorder/4)-1][3] = dsn | (3 << 10); }
	if (i1 == 01) { crystalSN[(Int_t)(dorder/4)-1][3] = dsn | (1 << 10); }
      }
    }
    
    st = fgets(str, 200, digMap);
  }

  fclose(digMap);
  
}

void INLCorrection::ReadRawINL() {
  char serialN[50];
  TString serN;

  for (Int_t i=0; i<MAXCRYSTALS; i++) {
    for (Int_t j=0; j<4; j++) {
      sprintf(serialN, "%X", crystalSN[i][j]);
      serN = serialN;

      TString inlName;
      inlName = "/nfs/fs1/home/users/crawford/INL/INL-0x" + serN + ".dat";
  
      FILE *inlIn;
      if ((inlIn = fopen(inlName.Data(), "r")) == NULL) {
	cerr << "Failed to open INL file " << inlName.Data() << ". " << endl;
	cerr << "Skipping it and moving on the next one... " << endl;
      } else {
	cout << "Reading INL corrections from " << inlName.Data() << ". " << endl;
	size_t siz = fread(inl[i][j*10], sizeof(inl), 1, inlIn);
	fclose(inlIn);
      }

      /* Re-expand the inl corrections */
      for (Int_t k=0; k<10; k++) {
      	for (Int_t m=0; m<8192; m++) {
	  inl[i][((j*10) + k)][m] = -2.0f*(inl[i][((j*10) + k)][m] + 0.0f);
	}
      } 
      
    } /* Loop over 4 digitizers / crystal */
  } /* Loop over MAXCRYSTALS crystals */

}
