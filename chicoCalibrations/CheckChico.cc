#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"
#include "TCutG.h"
using namespace std;

#include "GEBSort.h"

#define CHECKCHICO 0
#define PPAC_NUM 20
#define microsec 100
#define OverLap 1000
#define PGWinL -5
#define PGWinH 5
#define PGTKWinL 55
#define PGTKWinH 81
#define GTPosOffsetZ 5.334
#define GTPosOffsetY1 -14.71
#define GTPosOffsetY2 27.46
#define ch2ns 0.100
#define Ebeam 291.0
#define Abeam 76.0
#define Atarget 208.0
#define Qval 1.0
//#define P0 6667.7
#define P0 5800.0

#include "recalculatetheta.cc"
FILE *fp1;
unsigned long long int tsFirst;

float offTheta[PPAC_NUM],gainTheta[PPAC_NUM],quadTheta[PPAC_NUM];
float offPhi[PPAC_NUM],gainPhi[PPAC_NUM];
float offGT[QNum*4],gainGT[QNum*4],offdt[QNum*4];
int anodemap[PPAC_NUM][PPAC_NUM];
int cathodemap[4*PPAC_NUM][4*PPAC_NUM];
float betaP0,betaP1,betaP2,betaP3,betaP4,betaP5;
float betaT0,betaT1,betaT2,betaT3,betaT4,betaT5;
float betaB[180];

#if(CHECKCHICO)
TH2F *AnodeQDC,*AnodeTDC,*CathodeTDC;
TH1D *AnodeQDCNum;
TH2F *ACTDCNum;
TH2F *AnodeCathodePair,*AnodePair,*CathodePair;
TH2F *THETARAW,*PHIRAW;
TH2F *ThetaUD[PPAC_NUM],*PhiUD[PPAC_NUM];
#endif

TH1D *cRate,*gRate,*cgRate;
TH1D *ChicoDataLength;
TH1D *TKDataLength;

TH2F *projcosP,*projxtalP;
TH2F *projcosB,*projxtalB;
TH2F *projcosT,*projxtalT;
TH2F *ThetaGamP;
TH2F *ThetaGamT;
TH2F *nEgamP,*nEgamT;
TH2F *n1EgamP,*n1EgamT;
TH2F *projcosTKP,*projfomP;
TH2F *projcosTKB,*projfomB;
TH2F *projcosTKT,*projfomT;
TH2F *ThetaGamTKP;
TH2F *ThetaGamTKT;
TH2F *nEgamTKP,*nEgamTKT,*nEgamTKB;
TH2F *NEgamTKP,*NEgamTKT,*NEgamTKB;
TH2F *nfomP,*NfomP;
TH2F *nfomT,*NfomT;
TH2F *nfomB,*NfomB;
TH2F *MASSL,*MASSR;
TH2F *MASSLsingle,*MASSRsingle;
TH2F *histdTL[PPAC_NUM/2],*histdTR[PPAC_NUM/2];
TH2F *ThetaLR;
TH2F *fThetaL[PPAC_NUM/2],*fThetaR[PPAC_NUM/2];
TH2F *fThetaB,*fPhiB;
TH2F *ThetaB,*PhiB;
TH2F *fThetaBsingle;
TH2F *PhiLR;
TH2F *ThetaDL;
TH2F *ThetaDR;
TH2F *fPhiLR;
TH1D *dtGT_ChicoRaw;
TH1D *dtTK_ChicoRaw;
TH2F *dtGT_Chico;
TH2F *dtTK_Chico;
TH2F *dtGT_ChicoXtal;
TH1D *dtGT;
TH1D *dtChico;
TH1D *GTx,*GTy,*GTz;
TH2F *GTe;
TH2F *gGTMap,*GTMap,*ChicoMapALL,*gChicoMap;
TH2F *TKe,*gTKMap;
TH2F *nGT_Chico,*nGT_ChicoRaw;
TH2F *AnodeMapDouble,*AnodeMapPair,*AnodeMapNoPair,*AnodeMapFinal;
TH1D *AnodeNumAll,*AnodeNumPair,*AnodeNumNoPair,*AnodeNumFinal;
TH1D *CathodeNumAll,*CathodeNumPair,*CathodeNumNoPair,*CathodeNumFinal;
TH1D *CathodeNumTLRonly,*CathodeNumPLRonly;
TH2F *mass;
TH2F *ggP,*ggT;
TH2F *ggTKP,*ggTKT;

TCutG *massLCut,*massRCut;

float crmat[MAXDETPOS+1][MAXCRYSTALNO+1][4][4];
float recalculatethetaP(float, float);
float recalculatethetaT(float, float);
//////////////////////////////////////////////////////////////////////
void RdCRMAT(const char *fn) {

  FILE *fp;

  float f1, f2, f3, f4;
  int pos, xtal;
  int nn = 0;
  char *st, str[256];

  fp = fopen(fn, "r");
  if (fp == NULL) {
    printf("Could not open \"%s\".\n", fn);
    exit(1);
  }
  printf("\"%s\" open....", fn);

  /* Read values. */
  nn = 0;
  st = fgets(str, 256, fp);
  while (st != NULL) {
    if (str[0] == 35) {
      /* '#' comment line, do nothing */
    } else if (str[0] == 59) {
      /* ';' comment line, do nothing */
    } else if (str[0] == 10) {
      /* Empty line, do nothing */
    } else {
      sscanf(str, "%i %i", &pos, &xtal);
      for (int i=0; i<4; i++) {
        st = fgets(str, 256, fp);
        sscanf(str, "%f %f %f %f", &f1, &f2, &f3, &f4);
        crmat[pos-1][xtal][i][0] = f1;
        crmat[pos-1][xtal][i][1] = f2;
        crmat[pos-1][xtal][i][2] = f3;
        crmat[pos-1][xtal][i][3] = f4;
      }
      nn++;
    }

    /* Attempt to read the next line. */
    st = fgets(str, 256, fp);
  }

  printf("Read %i rotation matrix coefficients.\n", nn);

  /* Done! */
  fclose(fp);
}
//////////////////////////////////////////////////////////////////////
void ChicoInit(const char *GTCalib){
  int i,j,k;
  float angle;
  string OneLine;

  ifstream THETACALFILE("ppacTheta.cal", ios::in);
  if(!THETACALFILE.is_open()) {
    cerr << "Error opening Calibration file ppacTheta.cal"<< endl;  
    exit(1);
  }
  getline(THETACALFILE,OneLine);
  //cout << OneLine << endl;
  getline(THETACALFILE,OneLine);
  //cout << OneLine << endl;
  for(i=0;i<PPAC_NUM;i++){
    THETACALFILE >> k >> offTheta[i] >> gainTheta[i];
    //THETACALFILE >> k >> offTheta[i] >> gainTheta[i] >> quadTheta[i];
    //cout << setw(3) << k << setw(15) << offTheta[i] << setw(15) << gainTheta[i] << endl;
  }
  THETACALFILE.close();

  ifstream PHICALFILE("ppacPhi.cal", ios::in);
  if(!PHICALFILE.is_open()) {
    cerr << "Error opening Calibration file ppacPhi.cal "<< endl;  
    exit(1);
  }
  getline(PHICALFILE,OneLine);
  //cout << OneLine << endl;
  getline(PHICALFILE,OneLine);
  //cout << OneLine << endl;
  for(i=0;i<PPAC_NUM;i++){
    PHICALFILE >> k >> offPhi[i] >> gainPhi[i];
    //cout << setw(3) << k << setw(15) << offPhi[i] << setw(15) << gainPhi[i] << endl;
  }
  PHICALFILE.close();

  ifstream BETAFILE("betaE.dat", ios::in);
  if(!BETAFILE.is_open()) {
    cerr << "Error opening Doppler factor Calibration file!"<< endl;  
    exit(1);
  }
  BETAFILE >> betaP0 >> betaP1 >> betaP2 >> betaP3 >> betaP4 >> betaP5;
  BETAFILE >> betaT0 >> betaT1 >> betaT2 >> betaT3 >> betaT4 >> betaT5;
  cout << setw(15) << betaP0 << setw(15) << betaP1 << setw(15) << betaP2 << setw(15) << betaP3 << endl;
  cout << setw(15) << betaP4 << setw(15) << betaP5 << endl;
  cout << setw(15) << betaT0 << setw(15) << betaT1 << setw(15) << betaT2 << setw(15) << betaT3 << endl;
  cout << setw(15) << betaT4 << setw(15) << betaT5 << endl;
  for(i=0;i<96;i++)betaB[i]=0.0;
  for(i=96;i<=168;i++){
    BETAFILE >> angle >> betaB[i];
    //cout << setw(5) << angle << setw(15) << betaB[i] << endl;
  }
  BETAFILE.close();
  for(i=169;i<180;i++)betaB[i]=0.0;

  ifstream GTCALFILE(GTCalib, ios::in);
  if(!GTCALFILE.is_open()) {
    cerr << "Error opening Gretina Calibration file!"<< endl;  
    exit(1);
  }
  for(i=0;i<QNum*4;i++){
    GTCALFILE >> offGT[i] >> gainGT[i] >> offdt[i];
    //cout << setw(3) << i+1 << setw(15) << offGT[i] << setw(15) << gainGT[i] << setw(15) << offdt[i] << endl;
  }
  GTCALFILE.close();

  //cout << "anodemap:" << endl;
  for(i=0;i<10;i++){
    for(j=0;j<10;j++){
      anodemap[i][j]=0;
    }
  }
  for(i=0;i<5;i++){
    for(j=0;j<10;j++){
      if(j==(i+5)){
        anodemap[i][j]=1;
        anodemap[j][i]=-1;
      }
    }
  }

  /*for(i=0;i<5;i++){
    for(j=0;j<10;j++){
      if(anodemap[i][j]!=0){
        cout << setw(3) << i << setw(3) << j << setw(3) << anodemap[i][j] << endl;
        cout << setw(3) << j << setw(3) << i << setw(3) << anodemap[j][i] << endl;
      }
    }
  }*/

  //cout << "cathodemap:" << endl;
  for(i=0;i<PPAC_NUM*4;i++){
    for(j=0;j<PPAC_NUM*4;j++){
      cathodemap[i][j]=0;
    }
  }
  for(i=0;i<PPAC_NUM*4;i++){
    for(j=0;j<PPAC_NUM*4;j++){
      if((i%4)==0 && (j==i+1)){
        cathodemap[i][j]=1;
        cathodemap[j][i]=-1;
      }
      if((i%4)==2 && (j==i+1)){
        cathodemap[i][j]=2;
        cathodemap[j][i]=-2;
      }
    }
  }
  cathodemap[56][57]=0;
  cathodemap[57][56]=0;
  cathodemap[58][59]=0;
  cathodemap[59][58]=0;

  cathodemap[56][59]=1;
  cathodemap[59][56]=-1;
  cathodemap[58][57]=2;
  cathodemap[57][58]=-2;

  /*for(i=0;i<PPAC_NUM*4;i++){
    for(j=0;j<PPAC_NUM*4;j++){
      if(cathodemap[i][j]!=0){
        cout << setw(3) << i << setw(3) << j << setw(3) << cathodemap[i][j] << endl;
        cout << setw(3) << j << setw(3) << i << setw(3) << cathodemap[j][i] << endl;
      }
    }
  }*/

}
//////////////////////////////////////////////////////////////////////
void setuproot(TFile *fRoot){
  char str[127];
  char fn[127];
  int i;

  sprintf(str,"GTx");
  sprintf(fn,"GTx");
  GTx=(TH1D*)fRoot->Get(str);
  if(GTx == 0)
    GTx = new TH1D(str,fn,400,-200,200);
  
  sprintf(str,"GTy");
  sprintf(fn,"GTy");
  GTy=(TH1D*)fRoot->Get(str);
  if(GTy == 0)
    GTy = new TH1D(str,fn,400,-200,200);
  
  sprintf(str,"GTz");
  sprintf(fn,"GTz");
  GTz=(TH1D*)fRoot->Get(str);
  if(GTz == 0)
    GTz = new TH1D(str,fn,400,-200,200);
  
  sprintf(str,"GTe");
  sprintf(fn,"GTe");
  GTe=(TH2F*)fRoot->Get(str);
  if(GTe == 0)
    GTe = new TH2F(str,fn,4096,0,4096,QNum*4,0,QNum*4);
  
  sprintf(str,"GTMap");
  sprintf(fn,"GTMap");
  GTMap=(TH2F*)fRoot->Get(str);
  if(GTMap == 0)
    GTMap = new TH2F(str,fn,180,0,180,360,0,360);
  
  sprintf(str,"gGTMap");
  sprintf(fn,"gGTMap");
  gGTMap=(TH2F*)fRoot->Get(str);
  if(gGTMap == 0)
    gGTMap = new TH2F(str,fn,180,0,180,360,0,360);
  
  sprintf(str,"TKe");
  sprintf(fn,"TKe");
  TKe=(TH2F*)fRoot->Get(str);
  if(TKe == 0)
    TKe = new TH2F(str,fn,3072,0,3072,200,0,2.);
  
  sprintf(str,"gTKMap");
  sprintf(fn,"gTKMap");
  gTKMap=(TH2F*)fRoot->Get(str);
  if(gTKMap == 0)
    gTKMap = new TH2F(str,fn,180,0,180,360,0,360);
  
  sprintf(str,"ChicoMapALL");
  sprintf(fn,"ChicoMapALL");
  ChicoMapALL=(TH2F*)fRoot->Get(str);
  if(ChicoMapALL == 0)
    ChicoMapALL = new TH2F(str,fn,180,0,180,360,0,360);
 
  sprintf(str,"gChicoMap");
  sprintf(fn,"gChicoMap");
  gChicoMap=(TH2F*)fRoot->Get(str);
  if(gChicoMap == 0)
    gChicoMap = new TH2F(str,fn,180,0,180,360,0,360);
 
  sprintf(str,"cgRate");
  sprintf(fn,"cgRate");
  cgRate=(TH1D*)fRoot->Get(str);
  if(cgRate == 0)
    cgRate = new TH1D(str,fn,90000,0,90000);
  
  sprintf(str,"cRate");
  sprintf(fn,"cRate");
  cRate=(TH1D*)fRoot->Get(str);
  if(cRate == 0)
    cRate = new TH1D(str,fn,90000,0,90000);
  
  sprintf(str,"gRate");
  sprintf(fn,"gRate");
  gRate=(TH1D*)fRoot->Get(str);
  if(gRate == 0)
    gRate = new TH1D(str,fn,90000,0,90000);
 
  sprintf(str,"ChicoDataLength");
  sprintf(fn,"ChicoDataLength");
  ChicoDataLength=(TH1D*)fRoot->Get(str);
  if(ChicoDataLength == 0)
    ChicoDataLength = new TH1D(str,fn,9000,0,9000);
 
  sprintf(str,"TKDataLength");
  sprintf(fn,"TKDataLength");
  TKDataLength=(TH1D*)fRoot->Get(str);
  if(TKDataLength == 0)
    TKDataLength = new TH1D(str,fn,9000,0,9000);
 
#if(CHECKCHICO) 
  sprintf(str,"AnodeQDC");
  sprintf(fn,"AnodeQDC");
  AnodeQDC=(TH2F*)fRoot->Get(str);
  if(AnodeQDC == 0)
    AnodeQDC = new TH2F(str,fn,32,0,32,4096,0,4096);
 
  sprintf(str,"AnodeTDC");
  sprintf(fn,"AnodeTDC");
  AnodeTDC=(TH2F*)fRoot->Get(str);
  if(AnodeTDC == 0)
    AnodeTDC = new TH2F(str,fn,128,0,128,8192,-4096,4096);
 
  sprintf(str,"CathodeTDC");
  sprintf(fn,"CathodeTDC");
  CathodeTDC=(TH2F*)fRoot->Get(str);
  if(CathodeTDC == 0)
    CathodeTDC = new TH2F(str,fn,128,0,128,4096,0,4096);
 
  sprintf(str,"AnodeQDCNum");
  sprintf(fn,"AnodeQDCNum");
  AnodeQDCNum=(TH1D*)fRoot->Get(str);
  if(AnodeQDCNum == 0)
    AnodeQDCNum = new TH1D(str,fn,32,0,32);
 
  sprintf(str,"ACTDCNum");
  sprintf(fn,"ACTDCNum");
  ACTDCNum=(TH2F*)fRoot->Get(str);
  if(ACTDCNum == 0)
    ACTDCNum = new TH2F(str,fn,64,0,64,64,0,64);
  
  sprintf(str,"AnodeCathodePair");
  sprintf(fn,"AnodeCathodePair");
  AnodeCathodePair = (TH2F*)fRoot->Get(str);
  if(AnodeCathodePair == 0)
    AnodeCathodePair = new TH2F(str,fn,200,0,200,200,0,200);
  
  sprintf(str,"AnodePair");
  sprintf(fn,"AnodePair");
  AnodePair = (TH2F*)fRoot->Get(str);
  if(AnodePair == 0)
    AnodePair = new TH2F(str,fn,200,0,200,200,0,200);
  
  sprintf(str,"CathodePair");
  sprintf(fn,"CathodePair");
  CathodePair=(TH2F*)fRoot->Get(str);
  if(CathodePair== 0)
    CathodePair= new TH2F(str,fn,200,0,200,200,0,200);
  
  sprintf(str,"THETARAW");
  sprintf(fn,"THETARAW");
  THETARAW=(TH2F*)fRoot->Get(str);
  if(THETARAW == 0)
    THETARAW = new TH2F(str,fn,PPAC_NUM,0,PPAC_NUM,4096,-4096,4096);
  
  sprintf(str,"PHIRAW");
  sprintf(fn,"PHIRAW");
  PHIRAW=(TH2F*)fRoot->Get(str);
  if(PHIRAW == 0)
    PHIRAW = new TH2F(str,fn,PPAC_NUM,0,PPAC_NUM,4096,-4096,4096);
  
  for(i=0;i<PPAC_NUM;i++){  
    sprintf(str,"ThetaUD%2.2i",i+1);
    sprintf(fn,"ThetaUD%2.2i",i+1);
    ThetaUD[i]=(TH2F*)fRoot->Get(str);
    if(ThetaUD[i] == 0)
      ThetaUD[i] = new TH2F(str,fn,1024,0,1024,1024,0,1024);

    sprintf(str,"PhiUD%2.2i",i+1);
    sprintf(fn,"PhiUD%2.2i",i+1);
    PhiUD[i]=(TH2F*)fRoot->Get(str);
    if(PhiUD[i] == 0)
      PhiUD[i] = new TH2F(str,fn,1024,0,1024,1024,0,1024);
  }
#endif

  sprintf(str,"CathodeNumAll");
  sprintf(fn,"CathodeNumAll");
  CathodeNumAll=(TH1D*)fRoot->Get(str);
  if(CathodeNumAll == 0)
    CathodeNumAll = new TH1D(str,fn,64,0,64);
 
  sprintf(str,"AnodeNumAll");
  sprintf(fn,"AnodeNumAll");
  AnodeNumAll=(TH1D*)fRoot->Get(str);
  if(AnodeNumAll == 0)
    AnodeNumAll = new TH1D(str,fn,32,0,32);
 
  sprintf(str,"AnodeMapDouble");
  sprintf(fn,"AnodeMapDouble");
  AnodeMapDouble = (TH2F*)fRoot->Get(str);
  if(AnodeMapDouble == 0)
    AnodeMapDouble = new TH2F(str,fn,100,0,100,100,0,100);
  
  sprintf(str,"AnodeMapPair");
  sprintf(fn,"AnodeMapPair");
  AnodeMapPair = (TH2F*)fRoot->Get(str);
  if(AnodeMapPair == 0)
    AnodeMapPair = new TH2F(str,fn,100,0,100,100,0,100);
  
  sprintf(str,"CathodeNumPair");
  sprintf(fn,"CathodeNumPair");
  CathodeNumPair=(TH1D*)fRoot->Get(str);
  if(CathodeNumPair == 0)
    CathodeNumPair = new TH1D(str,fn,64,0,64);
 
  sprintf(str,"AnodeNumPair");
  sprintf(fn,"AnodeNumPair");
  AnodeNumPair=(TH1D*)fRoot->Get(str);
  if(AnodeNumPair == 0)
    AnodeNumPair = new TH1D(str,fn,32,0,32);
 
  sprintf(str,"AnodeMapNoPair");
  sprintf(fn,"AnodeMapNoPair");
  AnodeMapNoPair = (TH2F*)fRoot->Get(str);
  if(AnodeMapNoPair == 0)
    AnodeMapNoPair = new TH2F(str,fn,100,0,100,100,0,100);
  
  sprintf(str,"CathodeNumNoPair");
  sprintf(fn,"CathodeNumNoPair");
  CathodeNumNoPair=(TH1D*)fRoot->Get(str);
  if(CathodeNumNoPair == 0)
    CathodeNumNoPair = new TH1D(str,fn,64,0,64);
 
  sprintf(str,"AnodeNumNoPair");
  sprintf(fn,"AnodeNumNoPair");
  AnodeNumNoPair=(TH1D*)fRoot->Get(str);
  if(AnodeNumNoPair == 0)
    AnodeNumNoPair = new TH1D(str,fn,32,0,32);
 
  sprintf(str,"AnodeMapFinal");
  sprintf(fn,"AnodeMapFinal");
  AnodeMapFinal = (TH2F*)fRoot->Get(str);
  if(AnodeMapFinal == 0)
    AnodeMapFinal = new TH2F(str,fn,100,0,100,100,0,100);
  
  sprintf(str,"CathodeNumFinal");
  sprintf(fn,"CathodeNumFinal");
  CathodeNumFinal=(TH1D*)fRoot->Get(str);
  if(CathodeNumFinal == 0)
    CathodeNumFinal = new TH1D(str,fn,64,0,64);
 
  sprintf(str,"AnodeNumFinal");
  sprintf(fn,"AnodeNumFinal");
  AnodeNumFinal=(TH1D*)fRoot->Get(str);
  if(AnodeNumFinal == 0)
    AnodeNumFinal = new TH1D(str,fn,32,0,32);
 
  sprintf(str,"CathodeNumTLRonly");
  sprintf(fn,"CathodeNumTLRonly");
  CathodeNumTLRonly=(TH1D*)fRoot->Get(str);
  if(CathodeNumTLRonly == 0)
    CathodeNumTLRonly = new TH1D(str,fn,64,0,64);
 
  sprintf(str,"CathodeNumPLRonly");
  sprintf(fn,"CathodeNumPLRonly");
  CathodeNumPLRonly=(TH1D*)fRoot->Get(str);
  if(CathodeNumPLRonly == 0)
    CathodeNumPLRonly = new TH1D(str,fn,64,0,64);
 
  sprintf(str,"projcosP");
  sprintf(fn,"projcosP");
  projcosP=(TH2F*)fRoot->Get(str);
  if(projcosP == 0)
    projcosP = new TH2F(str,fn,3072,0,3072,200,-1,1);

  sprintf(str,"projxtalP");
  sprintf(fn,"projxtalP");
  projxtalP=(TH2F*)fRoot->Get(str);
  if(projxtalP == 0)
    projxtalP = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"projcosB");
  sprintf(fn,"projcosB");
  projcosB=(TH2F*)fRoot->Get(str);
  if(projcosB == 0)
    projcosB = new TH2F(str,fn,3072,0,3072,200,-1,1);

  sprintf(str,"projxtalB");
  sprintf(fn,"projxtalB");
  projxtalB=(TH2F*)fRoot->Get(str);
  if(projxtalB == 0)
    projxtalB = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"projcosT");
  sprintf(fn,"projcosT");
  projcosT=(TH2F*)fRoot->Get(str);
  if(projcosT == 0)
    projcosT = new TH2F(str,fn,3072,0,3072,200,-1,1);

  sprintf(str,"projxtalT");
  sprintf(fn,"projxtalT");
  projxtalT=(TH2F*)fRoot->Get(str);
  if(projxtalT == 0)
    projxtalT = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"nEgamP");
  sprintf(fn,"nEgamP");
  nEgamP=(TH2F*)fRoot->Get(str);
  if(nEgamP == 0)
    nEgamP = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"nEgamT");
  sprintf(fn,"nEgamT");
  nEgamT=(TH2F*)fRoot->Get(str);
  if(nEgamT == 0)
    nEgamT = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"n1EgamP");
  sprintf(fn,"n1EgamP");
  n1EgamP=(TH2F*)fRoot->Get(str);
  if(n1EgamP == 0)
    n1EgamP = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"n1EgamT");
  sprintf(fn,"n1EgamT");
  n1EgamT=(TH2F*)fRoot->Get(str);
  if(n1EgamT == 0)
    n1EgamT = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"projcosTKP");
  sprintf(fn,"projcosTKP");
  projcosTKP=(TH2F*)fRoot->Get(str);
  if(projcosTKP == 0)
    projcosTKP = new TH2F(str,fn,3072,0,3072,200,-1,1);

  sprintf(str,"projfomP");
  sprintf(fn,"projfomP");
  projfomP=(TH2F*)fRoot->Get(str);
  if(projfomP == 0)
    projfomP = new TH2F(str,fn,3072,0,3072,200,0,2);

  sprintf(str,"projcosTKB");
  sprintf(fn,"projcosTKB");
  projcosTKB=(TH2F*)fRoot->Get(str);
  if(projcosTKB == 0)
    projcosTKB = new TH2F(str,fn,3072,0,3072,200,-1,1);

  sprintf(str,"projfomB");
  sprintf(fn,"projfomB");
  projfomB=(TH2F*)fRoot->Get(str);
  if(projfomB == 0)
    projfomB = new TH2F(str,fn,3072,0,3072,200,0,2);

  sprintf(str,"projcosTKT");
  sprintf(fn,"projcosTKT");
  projcosTKT=(TH2F*)fRoot->Get(str);
  if(projcosTKT == 0)
    projcosTKT = new TH2F(str,fn,3072,0,3072,200,-1,1);

  sprintf(str,"projfomT");
  sprintf(fn,"projfomT");
  projfomT=(TH2F*)fRoot->Get(str);
  if(projfomT == 0)
    projfomT = new TH2F(str,fn,3072,0,3072,200,0,2);

  sprintf(str,"nEgamTKP");
  sprintf(fn,"nEgamTKP");
  nEgamTKP=(TH2F*)fRoot->Get(str);
  if(nEgamTKP == 0)
    nEgamTKP = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"nEgamTKT");
  sprintf(fn,"nEgamTKT");
  nEgamTKT=(TH2F*)fRoot->Get(str);
  if(nEgamTKT == 0)
    nEgamTKT = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"nEgamTKB");
  sprintf(fn,"nEgamTKB");
  nEgamTKB=(TH2F*)fRoot->Get(str);
  if(nEgamTKB == 0)
    nEgamTKB = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"NEgamTKP");
  sprintf(fn,"NEgamTKP");
  NEgamTKP=(TH2F*)fRoot->Get(str);
  if(NEgamTKP == 0)
    NEgamTKP = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"NEgamTKT");
  sprintf(fn,"NEgamTKT");
  NEgamTKT=(TH2F*)fRoot->Get(str);
  if(NEgamTKT == 0)
    NEgamTKT = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"NEgamTKB");
  sprintf(fn,"NEgamTKB");
  NEgamTKB=(TH2F*)fRoot->Get(str);
  if(NEgamTKB == 0)
    NEgamTKB = new TH2F(str,fn,3072,0,3072,QNum*4,0,QNum*4);

  sprintf(str,"nfomP");
  sprintf(fn,"nfomP");
  nfomP=(TH2F*)fRoot->Get(str);
  if(nfomP == 0)
    nfomP = new TH2F(str,fn,QNum*4,0,QNum*4,200,0,2);

  sprintf(str,"NfomP");
  sprintf(fn,"NfomP");
  NfomP=(TH2F*)fRoot->Get(str);
  if(NfomP == 0)
    NfomP = new TH2F(str,fn,QNum*4,0,QNum*4,200,0,2);

  sprintf(str,"nfomT");
  sprintf(fn,"nfomT");
  nfomT=(TH2F*)fRoot->Get(str);
  if(nfomT == 0)
    nfomT = new TH2F(str,fn,QNum*4,0,QNum*4,200,0,2);

  sprintf(str,"NfomT");
  sprintf(fn,"NfomT");
  NfomT=(TH2F*)fRoot->Get(str);
  if(NfomT == 0)
    NfomT = new TH2F(str,fn,QNum*4,0,QNum*4,200,0,2);

  sprintf(str,"nfomB");
  sprintf(fn,"nfomB");
  nfomB=(TH2F*)fRoot->Get(str);
  if(nfomB == 0)
    nfomB = new TH2F(str,fn,QNum*4,0,QNum*4,200,0,2);

  sprintf(str,"NfomB");
  sprintf(fn,"NfomB");
  NfomB=(TH2F*)fRoot->Get(str);
  if(NfomB == 0)
    NfomB = new TH2F(str,fn,QNum*4,0,QNum*4,200,0,2);

  sprintf(str,"ggP");
  sprintf(fn,"ggP");
  ggP=(TH2F*)fRoot->Get(str);
  if(ggP == 0)
    ggP = new TH2F(str,fn,3072,0,3072,3072,0,3072);

  sprintf(str,"ggT");
  sprintf(fn,"ggT");
  ggT=(TH2F*)fRoot->Get(str);
  if(ggT == 0)
    ggT = new TH2F(str,fn,3072,0,3072,3072,0,3072);

  sprintf(str,"ggTKP");
  sprintf(fn,"ggTKP");
  ggTKP=(TH2F*)fRoot->Get(str);
  if(ggTKP == 0)
    ggTKP = new TH2F(str,fn,3072,0,3072,3072,0,3072);

  sprintf(str,"ggTKT");
  sprintf(fn,"ggTKT");
  ggTKT=(TH2F*)fRoot->Get(str);
  if(ggTKT == 0)
    ggTKT = new TH2F(str,fn,3072,0,3072,3072,0,3072);

  sprintf(str,"dtGT_ChicoRaw");
  sprintf(fn,"dtGT_ChicoRaw");
  dtGT_ChicoRaw=(TH1D*)fRoot->Get(str);
  if(dtGT_ChicoRaw == 0)
    dtGT_ChicoRaw = new TH1D(str,fn,20000,-10000,10000);
  
  sprintf(str,"dtTK_ChicoRaw");
  sprintf(fn,"dtTK_ChicoRaw");
  dtTK_ChicoRaw=(TH1D*)fRoot->Get(str);
  if(dtTK_ChicoRaw == 0)
    dtTK_ChicoRaw = new TH1D(str,fn,2000,-1000,1000);
  
  sprintf(str,"dtGT_Chico");
  sprintf(fn,"dtGT_Chico");
  dtGT_Chico=(TH2F*)fRoot->Get(str);
  if(dtGT_Chico == 0)
    dtGT_Chico = new TH2F(str,fn,2000,-1000,1000,2048,0,2048);
  
  sprintf(str,"dtTK_Chico");
  sprintf(fn,"dtTK_Chico");
  dtTK_Chico=(TH2F*)fRoot->Get(str);
  if(dtTK_Chico == 0)
    dtTK_Chico = new TH2F(str,fn,2000,-1000,1000,2048,0,2048);
  
  sprintf(str,"dtGT_ChicoXtal");
  sprintf(fn,"dtGT_ChicoXtal");
  dtGT_ChicoXtal=(TH2F*)fRoot->Get(str);
  if(dtGT_ChicoXtal == 0)
    dtGT_ChicoXtal = new TH2F(str,fn,2000,-1000,1000,QNum*4,0,QNum*4);
  
  sprintf(str,"dtChico");
  sprintf(fn,"dtChico");
  dtChico=(TH1D*)fRoot->Get(str);
  if(dtChico == 0)
    dtChico = new TH1D(str,fn,8000,0,8000);
  
  sprintf(str,"dtGT");
  sprintf(fn,"dtGT");
  dtGT=(TH1D*)fRoot->Get(str);
  if(dtGT == 0)
    dtGT = new TH1D(str,fn,6000,-1000,5000);
  
  sprintf(str,"nGT_ChicoRaw");
  sprintf(fn,"nGT_ChicoRaw");
  nGT_ChicoRaw=(TH2F*)fRoot->Get(str);
  if(nGT_ChicoRaw == 0)
    nGT_ChicoRaw = new TH2F(str,fn,MaxGTNum,0,MaxGTNum,MaxChicoNum,0,MaxChicoNum);
  
  sprintf(str,"nGT_Chico");
  sprintf(fn,"nGT_Chico");
  nGT_Chico=(TH2F*)fRoot->Get(str);
  if(nGT_Chico == 0)
    nGT_Chico = new TH2F(str,fn,MaxGTNum,0,MaxGTNum,MaxChicoNum,0,MaxChicoNum);
 
  sprintf(str,"ThetaGamP");
  sprintf(fn,"ThetaGamP");
  ThetaGamP=(TH2F*)fRoot->Get(str);
  if(ThetaGamP == 0)
    ThetaGamP = new TH2F(str,fn,3072,0,3072,180,1,181);

  sprintf(str,"ThetaGamT");
  sprintf(fn,"ThetaGamT");
  ThetaGamT=(TH2F*)fRoot->Get(str);
  if(ThetaGamT == 0)
    ThetaGamT = new TH2F(str,fn,3072,0,3072,90,1,91);

  sprintf(str,"ThetaGamTKP");
  sprintf(fn,"ThetaGamTKP");
  ThetaGamTKP=(TH2F*)fRoot->Get(str);
  if(ThetaGamTKP == 0)
    ThetaGamTKP = new TH2F(str,fn,3072,0,3072,180,1,181);

  sprintf(str,"ThetaGamTKT");
  sprintf(fn,"ThetaGamTKT");
  ThetaGamTKT=(TH2F*)fRoot->Get(str);
  if(ThetaGamTKT == 0)
    ThetaGamTKT = new TH2F(str,fn,3072,0,3072,90,1,91);

  sprintf(str,"ThetaLR");
  sprintf(fn,"ThetaLR");
  ThetaLR=(TH2F*)fRoot->Get(str);
  if(ThetaLR == 0)
    ThetaLR = new TH2F(str,fn,2048,-1024,1024,2048,-1024,1024);

  for(i=0;i<PPAC_NUM/2;i++){
    sprintf(str,"fThetaL%2.2i",i+1);
    sprintf(fn,"fThetaL%2.2i",i+1);
    fThetaL[i]=(TH2F*)fRoot->Get(str);
    if(fThetaL[i] == 0)
      fThetaL[i] = new TH2F(str,fn,900,0,90,900,0,90);
  
    sprintf(str,"fThetaR%2.2i",i+1);
    sprintf(fn,"fThetaR%2.2i",i+1);
    fThetaR[i]=(TH2F*)fRoot->Get(str);
    if(fThetaR[i] == 0)
      fThetaR[i] = new TH2F(str,fn,900,0,90,900,0,90);
  
    sprintf(str,"histdTL%2.2i",i+1);
    sprintf(fn,"histdTL%2.2i",i+1);
    histdTL[i]=(TH2F*)fRoot->Get(str);
    if(histdTL[i] == 0)
      histdTL[i] = new TH2F(str,fn,4000,-2000,2000,2048,-1024,1024);
  
    sprintf(str,"histdTR%2.2i",i+1);
    sprintf(fn,"histdTR%2.2i",i+1);
    histdTR[i]=(TH2F*)fRoot->Get(str);
    if(histdTR[i]== 0)
      histdTR[i]= new TH2F(str,fn,4000,-2000,2000,2048,-1024,1024);
  }

  sprintf(str,"ThetaDL");
  sprintf(fn,"ThetaDL");
  ThetaDL=(TH2F*)fRoot->Get(str);
  if(ThetaDL == 0)
    ThetaDL = new TH2F(str,fn,900,0,180,900,0,18);

  sprintf(str,"ThetaDR");
  sprintf(fn,"ThetaDR");
  ThetaDR=(TH2F*)fRoot->Get(str);
  if(ThetaDR == 0)
    ThetaDR = new TH2F(str,fn,900,0,180,900,0,18);

  sprintf(str,"PhiLR");
  sprintf(fn,"PhiLR");
  PhiLR=(TH2F*)fRoot->Get(str);
  if(PhiLR == 0)
    PhiLR = new TH2F(str,fn,2048,-1024,1024,2048,-1024,1024);

  sprintf(str,"fPhiLR");
  sprintf(fn,"fPhiLR");
  fPhiLR=(TH2F*)fRoot->Get(str);
  if(fPhiLR == 0)
    fPhiLR = new TH2F(str,fn,3600,0,360,3600,0,360);

  sprintf(str,"MASSL");
  sprintf(fn,"MASSL");
  MASSL=(TH2F*)fRoot->Get(str);
  if(MASSL == 0)
    MASSL = new TH2F(str,fn,900,0,90,1000,-500,500);
  
  sprintf(str,"MASSR");
  sprintf(fn,"MASSR");
  MASSR=(TH2F*)fRoot->Get(str);
  if(MASSR == 0)
    MASSR = new TH2F(str,fn,900,0,90,1000,-500,500);
  
  sprintf(str,"MASSLsingle");
  sprintf(fn,"MASSLsingle");
  MASSLsingle=(TH2F*)fRoot->Get(str);
  if(MASSLsingle == 0)
    MASSLsingle = new TH2F(str,fn,900,0,90,1000,-500,500);
  
  sprintf(str,"MASSRsingle");
  sprintf(fn,"MASSRsingle");
  MASSRsingle=(TH2F*)fRoot->Get(str);
  if(MASSRsingle == 0)
    MASSRsingle = new TH2F(str,fn,900,0,90,1000,-500,500);
  
  sprintf(str,"mass");
  sprintf(fn,"mass");
  mass=(TH2F*)fRoot->Get(str);
  if(mass == 0)
    mass = new TH2F(str,fn,900,0,180,500,0,500);
  
  sprintf(str,"fThetaBsingle");
  sprintf(fn,"fThetaBsingle");
  fThetaBsingle=(TH2F*)fRoot->Get(str);
  if(fThetaBsingle == 0)
    fThetaBsingle = new TH2F(str,fn,900,90,180,10,10,20);
  
  sprintf(str,"fThetaB");
  sprintf(fn,"fThetaB");
  fThetaB=(TH2F*)fRoot->Get(str);
  if(fThetaB == 0)
    fThetaB = new TH2F(str,fn,900,90,180,10,10,20);
  
  sprintf(str,"fPhiB");
  sprintf(fn,"fPhiB");
  fPhiB=(TH2F*)fRoot->Get(str);
  if(fPhiB == 0)
    fPhiB = new TH2F(str,fn,900,0,360,10,10,20);
  
  sprintf(str,"ThetaB");
  sprintf(fn,"ThetaB");
  ThetaB=(TH2F*)fRoot->Get(str);
  if(ThetaB == 0)
    ThetaB = new TH2F(str,fn,2048,-1024,1024,10,10,20);
  
  sprintf(str,"PhiB");
  sprintf(fn,"PhiB");
  PhiB=(TH2F*)fRoot->Get(str);
  if(PhiB == 0)
    PhiB = new TH2F(str,fn,2048,-1024,1024,10,10,20);
  
  //fRoot->ls();
}
//////////////////////////////////////////////////////////////////////
void setupcut(){
  TFile *fcutg;
  char WinName[64];

  fcutg = new TFile("ChicoCut", "read");

  sprintf(WinName, "massLCut");
  massLCut = (TCutG*)fcutg->Get(WinName);
  if(massLCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  sprintf(WinName, "massRCut");
  massRCut = (TCutG*)fcutg->Get(WinName);
  if(massRCut == NULL){
    printf("Could Not Read 2d CutG File %s\n",WinName);
    exit(-1);
  }

  //fcutg->ls();
  fcutg->Close();
}
//////////////////////////////////////////////////////////////////////
float GetMass(float thetaL, float thetaR, float dL, float dR, float dT){
  float pL,pR;
  float Mass;
  pL = P0*sinf(thetaR)/sinf(thetaL+thetaR);
  pR = P0*sinf(thetaL)/sinf(thetaL+thetaR);
  Mass = (-0.032206*dT*ch2ns + dL/pL*(Atarget+Abeam))/(dL/pL+dR/pR);
  return Mass;
}
//////////////////////////////////////////////////////////////////////
static unsigned int*
GetGEBEvBuf(GEBheader *Header, const char *FileName) { 
  /* function will read the gretina data file and */
  /* extract the event buffer following each GEB headr, but not interpret them. */

  unsigned int      i;
  unsigned int      t1, t2, t3, t4;
  unsigned int	    *TEMP;

  if(!(TEMP = (unsigned int*) malloc(GEBHDRLENBYTES))) {
    printf("\007  ERROR: Could not malloc data buffer %i bytes.\n",GEBHDRLENBYTES);
    exit(-1);
  }
  if(fread(TEMP,GEBHDRLENBYTES,1,fp1) != 1){
    if (feof(fp1)){
      printf("End of file %s\n",FileName);
      return NULL;
    }
    printf("file read error %s\n",FileName);
    return NULL;
  }
  memcpy(Header,TEMP,GEBHDRLENBYTES);
  free(TEMP);
  if(!(TEMP = (unsigned int*) malloc(Header->length))) {
    printf("\007  ERROR: Could not malloc data buffer %i bytes (length).\n",Header->length);
    exit(-1);
  }
  if(fread(TEMP,Header->length,1,fp1) != 1){
    if (feof(fp1)){
      printf("End of file %s\n",FileName);
      return NULL;
    }
    printf("file read error %s\n",FileName);
    return NULL;
  }
  if(Header->type == MODETHREE){
    for(i=0;i<Header->length/4;i++){
      t1 = (TEMP[i] & 0x000000ff) << 24;
      t2 = (TEMP[i] & 0x0000ff00) << 8;
      t3 = (TEMP[i] & 0x00ff0000) >> 8;
      t4 = (TEMP[i] & 0xff000000) >> 24;
      TEMP[i] = t1 + t2 + t3 + t4;
    }
  }
  return TEMP;
}
//////////////////////////////////////////////////////////////////////
int GetChicoEvent(unsigned int *DataRecord, CHICOEVENT *ChicoEvent){
  unsigned short int EvSize;
  unsigned short int chan=0;
  int val=0,refval=0;
  unsigned int NextInt;
  int seenTrailer=0;
  int i,j,k=0;
  int CH_Counter=0;
  static int multiAnodeTDCNum=0;
  static int multiCathodeTDCNum=0;

  EvSize = DataRecord[0]/4;
  ChicoEvent->status = (unsigned short int)((DataRecord[0] & 0xffff0000)>>16);
  ChicoEvent->LEDts = (unsigned long long int) (DataRecord[1] & 0xffffffff);
  ChicoEvent->LEDts += ((unsigned long long int) (DataRecord[2] & 0xffff)<<32);
  j = 2;
  NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
  j++;
  if(NextInt != 0xffffffff){					//Anode_QDC
    //assert(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) == QDCHEADER);
    if(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT)!=QDCHEADER)return 0;
    //assert(((NextInt & QDCGEOMASK) >> QDCGEOSHIFT) == ANODE_E_VSN);
    if(((NextInt & QDCGEOMASK) >> QDCGEOSHIFT) != ANODE_E_VSN)return 0;
    CH_Counter = (NextInt & COUNTMASK) >> COUNTSHIFT;
    k=0;
    for(i=0;i<CH_Counter;i++){
      NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
      j++;
      //assert(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) == DATA);
      if(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != DATA)return 0;
      chan = (unsigned short int)((NextInt & QDCCHANMASK) >> QDCCHANSHIFT);
      val = (NextInt & QDCDATAMASK);
      if(chan<PPAC_NUM && val >0 && k<32){
        ChicoEvent->anode_qdc_ch[k] = chan;
        ChicoEvent->anode_qdc_val[k] = val;
        k++;
      }
    }
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    //assert(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) == QDCTRAILER);
    if(((NextInt & QDCTYPEMASK) >> QDCTYPESHIFT) != QDCTRAILER)return 0;
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    if(NextInt != 0xffffffff){
      while (1){
        NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
        j++;
        if(NextInt == 0xffffffff)break;
      }
    }
  }
  ChicoEvent->anode_qdc_num = k;
  CH_Counter=0;
  chan=0;val=0;
  ChicoEvent->SINGLE = false;
  NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
  j++;
  if(NextInt != 0xffffffff){					//Anode_TDC
    //assert((NextInt & TDCTYPEMASK) == TDCHEADER);
    if((NextInt & TDCTYPEMASK) != TDCHEADER)return 0;
    //assert((NextInt & TDCGEOMASK) == ANODE_T_VSN);
    if((NextInt & TDCGEOMASK) != ANODE_T_VSN)return 0;
    while(!seenTrailer){
      NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
      j++;
      switch(NextInt & TDCTYPEMASK){
        case DATA:
          chan = (unsigned short int)((NextInt & TDCCHANMASK) >> TDCCHANSHIFT); 
          val =(NextInt & TDCDATAMASK);
          if(chan != ANODE_REFCH && chan != RFCH && chan != SingleFlag){
            if(chan<PPAC_NUM && CH_Counter<128){
              ChicoEvent->anode_tdc_ch[CH_Counter] = chan;
              ChicoEvent->anode_tdc_val[CH_Counter] =val;
              CH_Counter++;
            }
          }
          else if (chan == RFCH){
            ChicoEvent->RF = val;
          }
          else if (chan == SingleFlag){
            ChicoEvent->SINGLE = true;
          }
          else if(chan == ANODE_REFCH){
            refval = (NextInt & TDCDATAMASK);
          }
          break;
        case TDCTRAILER:
          seenTrailer = 1;
          break;
        default:
          break;
      }
    }
    //if(refval > 0){
      for(i=0;i<CH_Counter;i++){
        ChicoEvent->anode_tdc_val[i] -= refval;
      }
      //ChicoEvent->RF -= refval;
    //}
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    if(NextInt != 0xffffffff) {
      multiAnodeTDCNum++;
      if((multiAnodeTDCNum%100000)==0)printf("*Warning* anode TDC package with multiple events: %i\n",multiAnodeTDCNum);
      while (1){
        NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
        j++;
        if(NextInt == 0xffffffff)break;
      }
    }
  }
  ChicoEvent->anode_tdc_num = CH_Counter;
  CH_Counter=0; seenTrailer=0;
  chan=0;val=0;
  NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
  j++;
  if(NextInt != 0xffffffff){					//Cathode_TDC
    //assert((NextInt & TDCTYPEMASK) == TDCHEADER);
    if((NextInt & TDCTYPEMASK) != TDCHEADER)return 0;
    //assert((NextInt & TDCGEOMASK) == CATHODE_T_VSN);
    if((NextInt & TDCGEOMASK) != CATHODE_T_VSN)return 0;
    while(!seenTrailer){
      NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
      j++;
      switch(NextInt & TDCTYPEMASK){
        case DATA:
          chan = (unsigned short int)((NextInt & TDCCHANMASK) >> TDCCHANSHIFT); 
          val =(NextInt & TDCDATAMASK);
          if(chan != CATHODE_REFCH){
            if(chan < PPAC_NUM*4 && CH_Counter<128){
              ChicoEvent->cathode_tdc_ch[CH_Counter] = chan;
              ChicoEvent->cathode_tdc_val[CH_Counter] = val;
              CH_Counter++;
            }
          }
          else if(chan == CATHODE_REFCH){
            refval = (NextInt & TDCDATAMASK);
          }
          break;
        case TDCTRAILER:
          seenTrailer = 1;
          break;
        default:
          break;
      }
    }
    //if(refval > 0){
      //for(i=0;i<CH_Counter;i++){
      //  ChicoEvent->cathode_tdc_val[i] -= refval;
      //}
    //}
    NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
    j++;
    if(NextInt != 0xffffffff) {
      multiCathodeTDCNum++;
      if((multiCathodeTDCNum%100000)==0)printf("*Warning* cathode TDC package with multiple events: %i\n",multiCathodeTDCNum);
      while (1){
        NextInt = ((DataRecord[j] & 0xffff0000) >> 16) +  ((DataRecord[j+1] & 0xffff) << 16);
        j++;
        if(NextInt == 0xffffffff)break;
      }
    }
  }
  ChicoEvent->cathode_tdc_num = CH_Counter;
  return 1;
}
//////////////////////////////////////////////////////////////////////
int GetParticle(CHICOEVENT *ChicoEvent, PARTICLE *particle){
  static float d = (float) (RAND_MAX) + 1.0;
  int i,j,l,m;
  int dT=0;
  int anode[2]={160};
  //int valid = 0;
  int validTheta=0,validPhi=0;
  int validT[2] = {0};
  int validP[2] = {0};
  int theta[2]={0};
  int phi[2]={0};
  float ftheta[2]={0.};
  float fphi[2]={0.};
  float dL,dR;
  unsigned short int temp1=0,temp2=0;
  int vtemp1=0,vtemp2=0,doubleOK;
  int vanodemap=0,vcathodemap=0;
  int cathodesign;

#if(CHECKCHICO)
      HistChico(ChicoEvent);
#endif
  
  doubleOK = 0;
  AnodeNumAll->Fill(ChicoEvent->anode_tdc_num);
  CathodeNumAll->Fill(ChicoEvent->cathode_tdc_num);
  for(i=0;i<ChicoEvent->anode_tdc_num;i++){
    if(ChicoEvent->anode_tdc_ch[i]>=10){
      for(l=0;l<ChicoEvent->cathode_tdc_num-1;l++){
        for(m=l+1;m<ChicoEvent->cathode_tdc_num;m++){
          vcathodemap=cathodemap[ChicoEvent->cathode_tdc_ch[l]][ChicoEvent->cathode_tdc_ch[m]];
          if(vcathodemap!=0){
            cathodesign=vcathodemap/abs(vcathodemap);
            if(ChicoEvent->cathode_tdc_ch[l]/4==ChicoEvent->anode_tdc_ch[i]){
              if(abs(vcathodemap)==1){
                theta[0]=(ChicoEvent->cathode_tdc_val[l]-ChicoEvent->cathode_tdc_val[m])*cathodesign;
                ftheta[0] = gainTheta[ChicoEvent->anode_tdc_ch[i]]*(float)theta[0];
                //ftheta[0] += quadTheta[ChicoEvent->anode_tdc_ch[i]]*pow((float)theta[0],2);
                ftheta[0] += offTheta[ChicoEvent->anode_tdc_ch[i]];
                ftheta[0] += (float)rand()/d-0.5;
                validT[0]=1;
              }
              if(abs(vcathodemap)==2){
                phi[0]=(ChicoEvent->cathode_tdc_val[l]-ChicoEvent->cathode_tdc_val[m])*cathodesign;
                fphi[0] = gainPhi[ChicoEvent->anode_tdc_ch[i]]*(float)phi[0];
                fphi[0] += offPhi[ChicoEvent->anode_tdc_ch[i]];
                //if(fphi[0]>=16.0 && fphi[l]<=20){
                  fphi[0] += 36.0*(float)(ChicoEvent->anode_tdc_ch[i]%10);       //10: back and front number of ppac
                  fphi[0] += (float)rand()/d-0.5;
                  validP[0]=1;
                //}
              }
            }
          }
        }
      }
      if(validT[0] && validP[0]){
        particle->t= (double)ChicoEvent->LEDts;
        particle->id=ChicoEvent->anode_tdc_ch[i];
        particle->thetaL=theta[0];
        particle->fthetaL=ftheta[0]*M_PI/180.;
        particle->phiL=phi[0];
        particle->fphiL=fphi[0]*M_PI/180.;
        particle->single = ChicoEvent->SINGLE;
        particle->back = true;
        return 1;
      }
    }
    if(ChicoEvent->anode_tdc_num >=2){
      for(j=i+1;j<ChicoEvent->anode_tdc_num;j++){
        vanodemap = anodemap[ChicoEvent->anode_tdc_ch[i]][ChicoEvent->anode_tdc_ch[j]];
        if(vanodemap!=0){
          temp1 = ChicoEvent->anode_tdc_ch[i];
          temp2 = ChicoEvent->anode_tdc_ch[j];
          vtemp1 = ChicoEvent->anode_tdc_val[i];
          vtemp2 = ChicoEvent->anode_tdc_val[j];
          doubleOK = 1;
        }
        AnodeMapDouble->Fill(ChicoEvent->anode_tdc_ch[i],ChicoEvent->anode_tdc_ch[j]);
      }
    }
  }
  validT[0]=0;validP[0]=0;

  if(doubleOK == 1){
    ChicoEvent->anode_tdc_num = 2;
    ChicoEvent->anode_tdc_ch[0] = temp1;
    ChicoEvent->anode_tdc_ch[1] = temp2;
    ChicoEvent->anode_tdc_val[0] = vtemp1;
    ChicoEvent->anode_tdc_val[1] = vtemp2;
  }

  if(ChicoEvent->anode_tdc_num == 2){
    vanodemap = anodemap[ChicoEvent->anode_tdc_ch[0]][ChicoEvent->anode_tdc_ch[1]];
    if(vanodemap!=0){
      for(i=0;i<2;i++){
        //anode[i] = ChicoEvent->anode_tdc_ch[(((i-vanodemap)+1)/2+i)%2];
        anode[i] = ChicoEvent->anode_tdc_ch[((vanodemap+1)/2==i)];      //0: left; 1: right
      }
      AnodeMapPair->Fill(ChicoEvent->anode_tdc_ch[0],ChicoEvent->anode_tdc_ch[1]);
      AnodeNumPair->Fill(ChicoEvent->anode_tdc_num);
      CathodeNumPair->Fill(ChicoEvent->cathode_tdc_num);
      dT = (float)(ChicoEvent->anode_tdc_val[0]-ChicoEvent->anode_tdc_val[1])*vanodemap;
      for(i=0;i<ChicoEvent->cathode_tdc_num-1;i++){
        for(j=i+1;j<ChicoEvent->cathode_tdc_num;j++){
          vcathodemap=cathodemap[ChicoEvent->cathode_tdc_ch[i]][ChicoEvent->cathode_tdc_ch[j]];
          if(vcathodemap!=0){
            cathodesign=vcathodemap/abs(vcathodemap);
            for(l=0;l<2;l++){
              if(ChicoEvent->cathode_tdc_ch[i]/4==anode[l]){
                if(abs(vcathodemap)==1){
                  theta[l]=(ChicoEvent->cathode_tdc_val[i]-ChicoEvent->cathode_tdc_val[j])*cathodesign;
                  ftheta[l] = gainTheta[anode[l]]*(float)theta[l];
                  //ftheta[l] += quadTheta[anode[l]]*pow((float)theta[l],2);
                  ftheta[l] += offTheta[anode[l]];
                  ftheta[l] += (float)rand()/d-0.5;
                  validT[l]=1;
                }
                if(abs(vcathodemap)==2){
                  phi[l]=(ChicoEvent->cathode_tdc_val[i]-ChicoEvent->cathode_tdc_val[j])*cathodesign;
                  fphi[l] = gainPhi[anode[l]]*(float)phi[l];
                  fphi[l] += offPhi[anode[l]];
                  //if(fphi[l]>=16.0 && fphi[l]<=20){
                    fphi[l] += 36.0*(float)(anode[l]%10);       //10: back and front number of ppac
                    fphi[l] += (float)rand()/d-0.5;
                    validP[l]=1;
                  //}
                }
              }
            }
          }
        }
      }
    }
    else{
      AnodeMapNoPair->Fill(ChicoEvent->anode_tdc_ch[0],ChicoEvent->anode_tdc_ch[1]);
      AnodeNumNoPair->Fill(ChicoEvent->anode_tdc_num);
      CathodeNumNoPair->Fill(ChicoEvent->cathode_tdc_num);
    } 
  }

  //if((validT[0]*validP[0]*validT[1]*validP[1])!=0)valid = 1;
  if(validT[0] == 1 && validT[1] ==1)validTheta = 1;
  if((validP[0]+validP[1]) >= 1)validPhi = 1;
  //if((validP[0]*validP[1]) == 1)validPhi = 1;
  if(((validT[0]+validT[1]) ==1)&&validPhi==1)CathodeNumTLRonly->Fill(ChicoEvent->cathode_tdc_num);
  if(((validP[0]+validP[1]) ==1)&&validTheta==1)CathodeNumPLRonly->Fill(ChicoEvent->cathode_tdc_num);
  if(validTheta==1 && validPhi==1){
    AnodeMapFinal->Fill(ChicoEvent->anode_tdc_ch[0],ChicoEvent->anode_tdc_ch[1]);
    AnodeNumFinal->Fill(ChicoEvent->anode_tdc_num);
    CathodeNumFinal->Fill(ChicoEvent->cathode_tdc_num);
    particle->t= (double)ChicoEvent->LEDts;
    particle->rf = ((double)ChicoEvent->RF * ch2ns)*0.1 + (double)rand()/d-0.5;
    particle->id=anode[0];
    particle->dT=dT;
    particle->thetaL=theta[0];
    particle->fthetaL=ftheta[0]*M_PI/180.;
    particle->thetaR=theta[1];
    particle->fthetaR=ftheta[1]*M_PI/180.;
    particle->phiL=phi[0];
    particle->phiR=phi[1];
    if(validP[0]==1 && validP[1] == 0){
      particle->fphiL = fphi[0]*M_PI/180.;
      particle->fphiR = fphi[0]*M_PI/180. + M_PI;
    }
    else if(validP[0]==0 && validP[1] == 1){
      particle->fphiR = fphi[1]*M_PI/180.;
      particle->fphiL = fphi[1]*M_PI/180. - M_PI;
    }
    else if(validP[0]==1 && validP[1] == 1){
      if(fphi[0]<180. && fphi[1]<360.){
        particle->fphiL=(fphi[0]+fphi[1]-180.)*0.5*M_PI/180.;
        particle->fphiR=(fphi[0]+fphi[1]+180.)*0.5*M_PI/180.;
      }
    }
    dL = 12.8/(0.75471*sinf(particle->fthetaL)*cosf(particle->fphiL-(18.0+(float)(anode[0]%10)*36.0)*M_PI/180.)
        +0.65606*cosf(particle->fthetaL));
    dR = 12.8/(0.75471*sinf(particle->fthetaR)*cosf(particle->fphiR-(18.0+(float)(anode[1]%10)*36.0)*M_PI/180.)
        +0.65606*cosf(particle->fthetaR));
    particle->dL = dL;
    particle->dR = dR;
    particle->single = ChicoEvent->SINGLE;
    particle->back = false;
    particle->Mass = GetMass(particle->fthetaL,particle->fthetaR,dL,dR,(float)dT);
    return 1;
  }
  else{
    return 0;
  }
}

//////////////////////////////////////////////////////////////////////
void PrintParticle(PARTICLE *particle){
    printf("particle->id=%i\n",particle->id);
    printf("particle->dT=%i\n",particle->dT);
    printf("particle->thetaL=%i\n",particle->thetaL);
    printf("particle->thetaR=%i\n",particle->thetaR);
    printf("particle->fthetaL=%f\n",particle->fthetaL);
    printf("particle->fthetaR=%f\n",particle->fthetaR);
    printf("particle->phiL=%i\n",particle->phiL);
    printf("particle->phiR=%i\n",particle->phiR);
    printf("particle->fphiL=%f\n",particle->fphiL);
    printf("particle->fphiR=%f\n",particle->fphiR);
    printf("particle->dL=%f\n",particle->dL);
    printf("particle->dR=%f\n",particle->dR);
    printf("particle->Mass=%f\n",particle->Mass);
}
//////////////////////////////////////////////////////////////////////
int GetGamma(GTEVENT2 *GTEvent2, GAMMA *gamma){
  int i,nmax=0;
  float emax = 0.0;
  int detectorPosition,crystalNumber;
  float R,r;

  if(GTEvent2->pad != 0) return 0;
  for(i=0;i<GTEvent2->num;i++){
    if(GTEvent2->intpts[i].e > emax)nmax = i;
  }
  detectorPosition = ((GTEvent2->crystal_id & 0xfffc)>>2)-1;
  crystalNumber = (GTEvent2->crystal_id & 0x0003);

  gamma->nintpts = GTEvent2->num;
  gamma->id = crystalNumber;
  switch(detectorPosition+1){
    case Q1:
      gamma->cc_id = 1;
      break;
    case Q2:
      gamma->cc_id = 2;
      break;
    case Q3:
      gamma->cc_id = 3;
      break;
    case Q4:
      gamma->cc_id = 4;
      break;
    case Q5:
      gamma->cc_id = 5;
      break;
    case Q6:
      gamma->cc_id = 6;
      break;
    case Q7:
      gamma->cc_id = 7;
      break;
    case Q8:
      gamma->cc_id = 8;
      break;
    case Q9:
      gamma->cc_id = 9;
      break;
    case Q10:
      gamma->cc_id = 10;
      break;
    case Q11:
      gamma->cc_id = 11;
      break;
    default:
      break;
  }

  gamma->x = ( (crmat[detectorPosition][crystalNumber][0][0] * GTEvent2->intpts[nmax].x) +
               (crmat[detectorPosition][crystalNumber][0][1] * GTEvent2->intpts[nmax].y) +
               (crmat[detectorPosition][crystalNumber][0][2] * GTEvent2->intpts[nmax].z) +
               (crmat[detectorPosition][crystalNumber][0][3]) );

  gamma->y = ( (crmat[detectorPosition][crystalNumber][1][0] * GTEvent2->intpts[nmax].x) +
               (crmat[detectorPosition][crystalNumber][1][1] * GTEvent2->intpts[nmax].y) +
               (crmat[detectorPosition][crystalNumber][1][2] * GTEvent2->intpts[nmax].z) +
               (crmat[detectorPosition][crystalNumber][1][3]) );

  gamma->z = ( (crmat[detectorPosition][crystalNumber][2][0] * GTEvent2->intpts[nmax].x) +
               (crmat[detectorPosition][crystalNumber][2][1] * GTEvent2->intpts[nmax].y) +
               (crmat[detectorPosition][crystalNumber][2][2] * GTEvent2->intpts[nmax].z) +
               (crmat[detectorPosition][crystalNumber][2][3]) );

  if(gamma->y<0)gamma->y+=GTPosOffsetY1;
  else if(gamma->y>=0)gamma->y+=GTPosOffsetY2;
  gamma->z+=GTPosOffsetZ;

  R = sqrtf( powf(gamma->x,2.0) + powf(gamma->y,2.0) + powf(gamma->z,2.0) );
  r = sqrtf( powf(gamma->x,2.0) + powf(gamma->y,2.0) );
  gamma->theta = acosf(gamma->z/R);
  gamma->phi = acosf(gamma->x/r);
  if(gamma->y < 0) gamma->phi = 2*M_PI - gamma->phi;

  gamma->e = GTEvent2->tot_e;

  //gamma->t = (double)GTEvent2->timestamp + (double)GTEvent2->t0/10.;
  gamma->t = (double)GTEvent2->timestamp;
  gamma->t0 = (double)GTEvent2->t0;
  
  return 1;
}
//////////////////////////////////////////////////////////////////////
void ProcessChico(PARTICLE *particle){
  if(!particle->back){
    ThetaLR->Fill(particle->thetaL,particle->thetaR);
    ThetaLR->Fill(particle->thetaL,particle->thetaR);
    PhiLR->Fill(particle->phiL,particle->phiR);
    if(massLCut->IsInside(particle->fthetaL/M_PI*180.,particle->dT)){
      fThetaL[particle->id]->Fill(particle->fthetaL/M_PI*180.,particle->fthetaR/M_PI*180.);
    }
    else if(massRCut->IsInside(particle->fthetaR/M_PI*180.,particle->dT)){
      fThetaR[particle->id]->Fill(particle->fthetaL/M_PI*180.,particle->fthetaR/M_PI*180.);
    }
    fPhiLR->Fill(particle->fphiL/M_PI*180.,particle->fphiR/M_PI*180.);
    ThetaDL->Fill(particle->fthetaL/M_PI*180.,particle->dL);
    ThetaDR->Fill(particle->fthetaR/M_PI*180.,particle->dR);
    histdTL[particle->id]->Fill(particle->thetaL,particle->dT);
    histdTR[particle->id]->Fill(particle->thetaR,particle->dT);
    ChicoMapALL->Fill(particle->fthetaL/M_PI*180.,particle->fphiL/M_PI*180.);
    ChicoMapALL->Fill(particle->fthetaR/M_PI*180.,particle->fphiR/M_PI*180.);
    MASSL->Fill(particle->fthetaL/M_PI*180.,particle->dT);
    MASSR->Fill(particle->fthetaR/M_PI*180.,particle->dT);
    if(particle->single){
      MASSLsingle->Fill(particle->fthetaL/M_PI*180.,particle->dT);
      MASSRsingle->Fill(particle->fthetaR/M_PI*180.,particle->dT);
    }
    mass->Fill(particle->fthetaL/M_PI*180.,particle->Mass);
  }
  else if(particle->back){
    ChicoMapALL->Fill(particle->fthetaL/M_PI*180.,particle->fphiL/M_PI*180.);
    ThetaB->Fill(particle->thetaL,particle->id);
    PhiB->Fill(particle->phiL,particle->id);
    fThetaB->Fill(particle->fthetaL/M_PI*180.,particle->id);
    fPhiB->Fill(particle->fphiL/M_PI*180.,particle->id);
    if(particle->single){
      fThetaBsingle->Fill(particle->fthetaL/M_PI*180.,particle->id);
    }
  }
}
//////////////////////////////////////////////////////////////////////
void ProcessGamma(GAMMA *gamma){
  float e;
  int id;
  GTMap->Fill(gamma->theta/M_PI*180.,gamma->phi/M_PI*180.);
  GTx->Fill(gamma->x);
  GTy->Fill(gamma->y);
  GTz->Fill(gamma->z);
  id = gamma->id+(gamma->cc_id-1)*4;
  e = gamma->e*gainGT[id]+ offGT[id];
  GTe->Fill(e,id);
}
//////////////////////////////////////////////////////////////////////
void ProcessGammaTK(TRACKED_GAMMA_HIT *GTEvent1){
  int i;
  for(i=0;i<GTEvent1->ngam;i++){
    TKe->Fill(GTEvent1->gr[i].esum,GTEvent1->gr[i].fom);
  }
}
//////////////////////////////////////////////////////////////////////
float GCcos(float ptheta, float pphi, float gtheta, float gphi){
  float gccos;
  gccos=sinf(ptheta)*sinf(gtheta)*cosf(pphi-gphi)+cosf(ptheta)*cosf(gtheta);
  return gccos;
}
//////////////////////////////////////////////////////////////////////
void HistChico(CHICOEVENT *ChicoEvent){
#if(CHECKCHICO)
  int i,j;
  int phiRaw,thetaRaw,id;

  ACTDCNum->Fill(ChicoEvent->anode_tdc_num,ChicoEvent->cathode_tdc_num);

  if(ChicoEvent->anode_tdc_num == 2)

  for(i=0;i<ChicoEvent->anode_tdc_num-1;i++){
    for(j=i+1;j<ChicoEvent->anode_tdc_num;j++){
      AnodePair->Fill(ChicoEvent->anode_tdc_ch[i],ChicoEvent->anode_tdc_ch[j]);
    }
  }
  for(i=0;i<ChicoEvent->cathode_tdc_num-1;i++){
    for(j=i+1;j<ChicoEvent->cathode_tdc_num;j++){
      CathodePair->Fill(ChicoEvent->cathode_tdc_ch[i],ChicoEvent->cathode_tdc_ch[j]);
    }
  }
  for(i=0;i<ChicoEvent->anode_tdc_num;i++){
    for(j=0;j<ChicoEvent->cathode_tdc_num;j++){
      AnodeCathodePair->Fill(ChicoEvent->anode_tdc_ch[i],ChicoEvent->cathode_tdc_ch[j]);
    }
  }

  for(i=0;i<ChicoEvent->cathode_tdc_num;i++){
    for(j=i+1;j<ChicoEvent->cathode_tdc_num;j++){
      CathodeTDC->Fill(ChicoEvent->cathode_tdc_ch[i],ChicoEvent->cathode_tdc_val[j]);
      if((ChicoEvent->cathode_tdc_ch[i]/4)==(ChicoEvent->cathode_tdc_ch[j]/4)){
        id = ChicoEvent->cathode_tdc_ch[i]/4;
        if((ChicoEvent->cathode_tdc_ch[i]%4)==0 && (ChicoEvent->cathode_tdc_ch[j]%4)==1){
          thetaRaw=ChicoEvent->cathode_tdc_val[i]-ChicoEvent->cathode_tdc_val[j];
          THETARAW->Fill(id,thetaRaw);
          ThetaUD[id]->Fill(ChicoEvent->cathode_tdc_val[i]/4,ChicoEvent->cathode_tdc_val[j]/4);
        }
        if((ChicoEvent->cathode_tdc_ch[j]%4)==0 && (ChicoEvent->cathode_tdc_ch[i]%4)==1){
          thetaRaw=ChicoEvent->cathode_tdc_val[j]-ChicoEvent->cathode_tdc_val[i];
          THETARAW->Fill(id,thetaRaw);
          ThetaUD[id]->Fill(ChicoEvent->cathode_tdc_val[j]/4,ChicoEvent->cathode_tdc_val[i]/4);
        }
        if((ChicoEvent->cathode_tdc_ch[i]%4)==2 && (ChicoEvent->cathode_tdc_ch[j]%4)==3){
          phiRaw=ChicoEvent->cathode_tdc_val[i]-ChicoEvent->cathode_tdc_val[j];
          PHIRAW->Fill(id,phiRaw);
          PhiUD[id]->Fill(ChicoEvent->cathode_tdc_val[i]/4,ChicoEvent->cathode_tdc_val[j]/4);
        }
        if((ChicoEvent->cathode_tdc_ch[j]%4)==3 && (ChicoEvent->cathode_tdc_ch[i]%4)==2){
          phiRaw=ChicoEvent->cathode_tdc_val[j]-ChicoEvent->cathode_tdc_val[i];
          PHIRAW->Fill(id,phiRaw);
          PhiUD[id]->Fill(ChicoEvent->cathode_tdc_val[j]/4,ChicoEvent->cathode_tdc_val[i]/4);
        }
      }
    }
  }
#endif
}
//////////////////////////////////////////////////////////////////////////////
void InitCoinEvent(COINCIDENCE *CoinEvent){
  CoinEvent->nCoinGT = 0;
  CoinEvent->nCoinTK = 0;
  CoinEvent->nCoinChico = 0;
}
//////////////////////////////////////////////////////////////////////////////
void IncrementCoinEvent(COINCIDENCE *CoinEvent, GEBheader *Header, unsigned int *EventBuf){
  GTEVENT2   GTEvent2;
  TRACKED_GAMMA_HIT GTEvent1;
  CHICOEVENT ChicoEvent;
  GAMMA      gamma;
  PARTICLE   particle;
  static unsigned long long int tsLastChico=0;
  static unsigned long long int tsLastGT=0;

  if(Header->type == MODETWO && CoinEvent->nCoinGT < MaxGTNum){
    assert(Header->length == sizeof(GTEVENT2));
    memcpy(&GTEvent2,EventBuf,Header->length);
    dtGT->Fill(GTEvent2.timestamp-tsLastGT);
    tsLastGT = GTEvent2.timestamp;
    gRate->Fill((double)(GTEvent2.timestamp-tsFirst)*1E-8);
    if(GetGamma(&GTEvent2,&gamma)!=0){
      ProcessGamma(&gamma);
      CoinEvent->gamma[CoinEvent->nCoinGT] = gamma;
      CoinEvent->nCoinGT++;
    }
  }
  else if(Header->type == MODEONE){
    memcpy(&GTEvent1,EventBuf,Header->length);
    ProcessGammaTK(&GTEvent1);
    CoinEvent->gammatk = GTEvent1;
    CoinEvent->nCoinTK++;
  }
  else if(Header->type == CHICOTYPE && CoinEvent->nCoinChico < MaxChicoNum && Header->length < 1024){
    if(GetChicoEvent(EventBuf,&ChicoEvent)==1){
      dtChico->Fill(ChicoEvent.LEDts-tsLastChico);
      tsLastChico = ChicoEvent.LEDts;
      cRate->Fill((double)(ChicoEvent.LEDts-tsFirst)*1E-8);
      if(GetParticle(&ChicoEvent,&particle)==1){
        ProcessChico(&particle);
        //PrintParticle(&particle);
        CoinEvent->particle[CoinEvent->nCoinChico] = particle;
        CoinEvent->nCoinChico++;
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
void ProcessEventMode2(COINCIDENCE *CoinEvent){
  int i,j;
  double deltaT;
  float gccos;
  float BETA,eg;
  float theta;
  int   xtal_id;
  int ngam = 0;
  int ngamB = 0;
  float egamP[MaxGTNum];
  float egamB[MaxGTNum];
  float egamT[MaxGTNum];
  int nintpts[MaxGTNum];

  nGT_Chico->Fill(CoinEvent->nCoinGT,CoinEvent->nCoinChico);
  //for(j=0;j<CoinEvent->nCoinChico;j++){
  j=0;
    for(i=0;i<CoinEvent->nCoinGT;i++){
      xtal_id = (CoinEvent->gamma[i].cc_id-1)*4 + CoinEvent->gamma[i].id;
      deltaT=CoinEvent->particle[j].t - (CoinEvent->gamma[i].t + CoinEvent->gamma[i].t0);
      deltaT-=(double)offdt[xtal_id];
      dtGT_Chico->Fill(deltaT,CoinEvent->gamma[i].e);
      if(CoinEvent->gamma[i].e>=200.0)dtGT_ChicoXtal->Fill(deltaT,xtal_id);
      if(deltaT >= PGWinL && deltaT <= PGWinH ){
        gGTMap->Fill(CoinEvent->gamma[i].theta/M_PI*180.,CoinEvent->gamma[i].phi/M_PI*180.);
        cgRate->Fill((CoinEvent->particle[j].t-(double)tsFirst)*1e-8);
        gChicoMap->Fill(CoinEvent->particle[j].fthetaL/M_PI*180.,CoinEvent->particle[j].fphiL/M_PI*180.);
        gChicoMap->Fill(CoinEvent->particle[j].fthetaR/M_PI*180.,CoinEvent->particle[j].fphiR/M_PI*180.);
        if(CoinEvent->particle[j].back){
          theta=CoinEvent->particle[j].fthetaL;
          if((theta/M_PI*180.)>=96. && (theta/M_PI*180.)<=168.){
            gccos = GCcos(theta,CoinEvent->particle[j].fphiL,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
            BETA = betaB[(int)(theta/M_PI*180.)];
            eg = CoinEvent->gamma[i].e*gainGT[xtal_id];
            eg+= offGT[xtal_id];
            eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
            ThetaGamP->Fill(eg,theta/M_PI*180.);
            projcosB->Fill(eg,gccos);
            projxtalB->Fill(eg,xtal_id);
            egamB[ngam]=eg;
            ngamB++;
          }
        }
        else if(massLCut->IsInside(CoinEvent->particle[j].fthetaL/M_PI*180.,CoinEvent->particle[j].dT)){
        //if(CoinEvent->particle[j].Mass>=65 && CoinEvent->particle[j].Mass <= 85){
          //theta=recalculatethetaP(CoinEvent->particle[j].fthetaL,CoinEvent->particle[j].fthetaR);
          theta=CoinEvent->particle[j].fthetaL;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiL,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaL/M_PI*180.)];
          BETA = betaP0 + betaP1*theta;
          BETA+= betaP2*pow(theta,2);
          BETA+= betaP3*pow(theta,3);
          BETA+= betaP4*pow(theta,4);
          BETA+= betaP5*pow(theta,5);
          eg = CoinEvent->gamma[i].e*gainGT[xtal_id];
          eg+= offGT[xtal_id];
          if(gccos>=-0.05 && gccos<=0.05)projxtalP->Fill(eg,xtal_id);
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          ThetaGamP->Fill(eg,theta/M_PI*180.);
          //eg = CoinEvent->gamma[i].e*(1-BETA*gccos)/sqrt(1-pow(BETA,2)); 
          projcosP->Fill(eg,gccos);
          egamP[ngam]=eg;
          //theta=recalculatethetaT(CoinEvent->particle[j].fthetaL,CoinEvent->particle[j].fthetaR);
          theta=CoinEvent->particle[j].fthetaR;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiR,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaL/M_PI*180.)];
          BETA = betaT0 + betaT1*theta;
          BETA+= betaT2*pow(theta,2);
          BETA+= betaT3*pow(theta,3);
          BETA+= betaT4*pow(theta,4);
          BETA+= betaT5*pow(theta,5);
          eg = CoinEvent->gamma[i].e*gainGT[xtal_id];
          eg+= offGT[xtal_id];
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          ThetaGamT->Fill(eg,theta/M_PI*180.);
          //eg = CoinEvent->gamma[i].e*(1-BETA*gccos)/sqrt(1-pow(BETA,2)); 
          projcosT->Fill(eg,gccos);
          projxtalT->Fill(eg,xtal_id);
          nintpts[ngam] = CoinEvent->gamma[i].nintpts;
          egamT[ngam]=eg;
          ngam++;
        }
        else if(massRCut->IsInside(CoinEvent->particle[j].fthetaR/M_PI*180.,CoinEvent->particle[j].dT)){	
        //else if(CoinEvent->particle[j].Mass>=198 && CoinEvent->particle[j].Mass <= 218){
          //theta=recalculatethetaP(CoinEvent->particle[j].fthetaR,CoinEvent->particle[j].fthetaL);
          theta=CoinEvent->particle[j].fthetaR;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiR,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaR/M_PI*180.)];
          BETA = betaP0 + betaP1*theta;
          BETA+= betaP2*pow(theta,2);
          BETA+= betaP3*pow(theta,3);
          BETA+= betaP4*pow(theta,4);
          BETA+= betaP5*pow(theta,5);
          eg = CoinEvent->gamma[i].e*gainGT[CoinEvent->gamma[i].id+(CoinEvent->gamma[i].cc_id-1)*4];
          eg+= offGT[CoinEvent->gamma[i].id+(CoinEvent->gamma[i].cc_id-1)*4];
          if(gccos>=-0.05 && gccos<=0.05)projxtalP->Fill(eg,xtal_id);
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          ThetaGamP->Fill(eg,theta/M_PI*180.);
          //eg = CoinEvent->gamma[i].e*(1-BETA*gccos)/sqrt(1-pow(BETA,2)); 
	  projcosP->Fill(eg,gccos);
          egamP[ngam]=eg;

          //theta=recalculatethetaT(CoinEvent->particle[j].fthetaR,CoinEvent->particle[j].fthetaL);
          theta=CoinEvent->particle[j].fthetaL;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiL,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaR/M_PI*180.)];
          BETA = betaT0 + betaT1*theta;
          BETA+= betaT2*pow(theta,2);
          BETA+= betaT3*pow(theta,3);
          BETA+= betaT4*pow(theta,4);
          BETA+= betaT5*pow(theta,5);
          eg = CoinEvent->gamma[i].e*gainGT[CoinEvent->gamma[i].id+(CoinEvent->gamma[i].cc_id-1)*4];
          eg+= offGT[CoinEvent->gamma[i].id+(CoinEvent->gamma[i].cc_id-1)*4];
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          ThetaGamT->Fill(eg,theta/M_PI*180.);
          //eg = CoinEvent->gamma[i].e*(1-BETA*gccos)/sqrt(1-pow(BETA,2)); 
	  projcosT->Fill(eg,gccos);
          projxtalT->Fill(eg,CoinEvent->gamma[i].id+(CoinEvent->gamma[i].cc_id-1)*4);
          nintpts[ngam] = CoinEvent->gamma[i].nintpts;
          egamT[ngam]=eg;
          ngam++;
        }
      }
    }
    if(ngam == 1){
      n1EgamP->Fill(egamP[0],nintpts[0]);
      n1EgamT->Fill(egamT[0],nintpts[0]);
    }
    if(ngam>0){
      for(i=0;i<ngam;i++){
        nEgamP->Fill(egamP[i],ngam);
        nEgamT->Fill(egamT[i],ngam);
        for(j=i+1;j<ngam;j++){
          ggP->Fill(egamP[i],egamP[j]);
          ggP->Fill(egamP[j],egamP[i]);
          ggT->Fill(egamT[i],egamT[j]);
          ggT->Fill(egamT[j],egamT[i]);
        }
      }
    }
  //}
}
//////////////////////////////////////////////////////////////////////////////
void ProcessEventMode1(COINCIDENCE *CoinEvent){
  int i,j;
  double deltaT;
  float R,r;
  float gtheta,gphi;
  float theta;
  float gccos;
  float BETA,eg;
  int ngam = 0;
  float egamP[MaxGTNum];
  float egamT[MaxGTNum];

  //for(j=0;j<CoinEvent->nCoinChico;j++){
  j=0;
    for(i=0;i<CoinEvent->gammatk.ngam;i++){
      deltaT=CoinEvent->particle[j].t - (double)CoinEvent->gammatk.gr[i].timestamp;
      dtTK_Chico->Fill(deltaT,CoinEvent->gammatk.gr[i].esum);
      if(deltaT >= PGTKWinL && deltaT <= PGTKWinH ){
        R = sqrtf( powf(CoinEvent->gammatk.gr[i].x0,2.0) + powf((CoinEvent->gammatk.gr[i].y0),2.0) + powf(CoinEvent->gammatk.gr[i].z0,2.0) );
        r = sqrtf( powf(CoinEvent->gammatk.gr[i].x0,2.0) + powf((CoinEvent->gammatk.gr[i].y0),2.0) );
        gtheta = acosf(CoinEvent->gammatk.gr[i].z0/R);
        gphi = acosf(CoinEvent->gammatk.gr[i].x0/r);
        if(CoinEvent->gammatk.gr[i].y0 < 0) gphi = 2*M_PI - gphi;
        gTKMap->Fill(gtheta/M_PI*180.,gphi/M_PI*180.);
        if(CoinEvent->particle[j].back){
          theta=CoinEvent->particle[j].fthetaL;
          if((theta/M_PI*180.)>=96. && (theta/M_PI*180.)<=168.){
            gccos = GCcos(theta,CoinEvent->particle[j].fphiL,gtheta,gphi);
            BETA = betaB[(int)(theta/M_PI*180.)];
            eg = CoinEvent->gammatk.gr[i].esum;
            eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
            projfomB->Fill(eg,CoinEvent->gammatk.gr[i].fom);
            nEgamTKB->Fill(eg,CoinEvent->gammatk.gr[i].ndet);
            NEgamTKB->Fill(eg,CoinEvent->gammatk.ngam);
            nfomB->Fill(CoinEvent->gammatk.gr[i].ndet,CoinEvent->gammatk.gr[i].fom);
            NfomB->Fill(CoinEvent->gammatk.ngam,CoinEvent->gammatk.gr[i].fom);
            if(CoinEvent->gammatk.gr[i].fom <= 0.8){
              ThetaGamTKP->Fill(eg,theta/M_PI*180.);
              projcosTKB->Fill(eg,gccos);
              egamP[ngam]=eg;
              ngam++;
            }
          }
        }
        else if(massLCut->IsInside(CoinEvent->particle[j].fthetaL/M_PI*180.,CoinEvent->particle[j].dT)){
        //if(CoinEvent->particle[j].Mass>=65 && CoinEvent->particle[j].Mass <= 85){
          //theta=recalculatethetaP(CoinEvent->particle[j].fthetaL,CoinEvent->particle[j].fthetaR);
          theta=CoinEvent->particle[j].fthetaL;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiL,gtheta,gphi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaL/M_PI*180.)];
          BETA = betaP0 + betaP1*theta;
          BETA+= betaP2*pow(theta,2);
          BETA+= betaP3*pow(theta,3);
          BETA+= betaP4*pow(theta,4);
          BETA+= betaP5*pow(theta,5);
          eg = CoinEvent->gammatk.gr[i].esum;
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          projfomP->Fill(eg,CoinEvent->gammatk.gr[i].fom);
          nEgamTKP->Fill(eg,CoinEvent->gammatk.gr[i].ndet);
          NEgamTKP->Fill(eg,CoinEvent->gammatk.ngam);
          nfomP->Fill(CoinEvent->gammatk.gr[i].ndet,CoinEvent->gammatk.gr[i].fom);
          NfomP->Fill(CoinEvent->gammatk.ngam,CoinEvent->gammatk.gr[i].fom);
          if(CoinEvent->gammatk.gr[i].fom <= 0.8){
            ThetaGamTKP->Fill(eg,theta/M_PI*180.);
            projcosTKP->Fill(eg,gccos);
            egamP[ngam]=eg;
          }
          //theta=recalculatethetaT(CoinEvent->particle[j].fthetaL,CoinEvent->particle[j].fthetaR);
          theta=CoinEvent->particle[j].fthetaR;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiR,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          BETA = betaT0 + betaT1*theta;
          BETA+= betaT2*pow(theta,2);
          BETA+= betaT3*pow(theta,3);
          BETA+= betaT4*pow(theta,4);
          BETA+= betaT5*pow(theta,5);
          eg = CoinEvent->gammatk.gr[i].esum;
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          projfomT->Fill(eg,CoinEvent->gammatk.gr[i].fom);
          nEgamTKT->Fill(eg,CoinEvent->gammatk.gr[i].ndet);
          NEgamTKT->Fill(eg,CoinEvent->gammatk.ngam);
          nfomT->Fill(CoinEvent->gammatk.gr[i].ndet,CoinEvent->gammatk.gr[i].fom);
          NfomT->Fill(CoinEvent->gammatk.ngam,CoinEvent->gammatk.gr[i].fom);
          if(CoinEvent->gammatk.gr[i].fom <= 0.8){
            ThetaGamTKT->Fill(eg,theta/M_PI*180.);
            projcosTKT->Fill(eg,gccos);
            egamT[ngam]=eg;
            ngam++;
          }
        }
        else if(massRCut->IsInside(CoinEvent->particle[j].fthetaR/M_PI*180.,CoinEvent->particle[j].dT)){	
        //else if(CoinEvent->particle[j].Mass>=198 && CoinEvent->particle[j].Mass <= 218){
          //theta=recalculatethetaP(CoinEvent->particle[j].fthetaR,CoinEvent->particle[j].fthetaL);
          theta=CoinEvent->particle[j].fthetaR;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiR,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaR/M_PI*180.)];
          BETA = betaP0 + betaP1*theta;
          BETA+= betaP2*pow(theta,2);
          BETA+= betaP3*pow(theta,3);
          BETA+= betaP4*pow(theta,4);
          BETA+= betaP5*pow(theta,5);
          eg = CoinEvent->gammatk.gr[i].esum;
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          projfomP->Fill(eg,CoinEvent->gammatk.gr[i].fom);
          nEgamTKP->Fill(eg,CoinEvent->gammatk.gr[i].ndet);
          NEgamTKP->Fill(eg,CoinEvent->gammatk.ngam);
          nfomP->Fill(CoinEvent->gammatk.gr[i].ndet,CoinEvent->gammatk.gr[i].fom);
          NfomP->Fill(CoinEvent->gammatk.ngam,CoinEvent->gammatk.gr[i].fom);
          if(CoinEvent->gammatk.gr[i].fom <= 0.8){
            ThetaGamTKP->Fill(eg,theta/M_PI*180.);
	    projcosTKP->Fill(eg,gccos);
            egamP[ngam]=eg;
          }
          //theta=recalculatethetaT(CoinEvent->particle[j].fthetaR,CoinEvent->particle[j].fthetaL);
          theta=CoinEvent->particle[j].fthetaL;
          gccos = GCcos(theta,CoinEvent->particle[j].fphiL,CoinEvent->gamma[i].theta,CoinEvent->gamma[i].phi);
          //BETA = betaP[(int)(CoinEvent->particle[j].fthetaR/M_PI*180.)];
          BETA = betaT0 + betaT1*theta;
          BETA+= betaT2*pow(theta,2);
          BETA+= betaT3*pow(theta,3);
          BETA+= betaT4*pow(theta,4);
          BETA+= betaT5*pow(theta,5);
          eg = CoinEvent->gammatk.gr[i].esum;
          eg*= (1.-BETA*gccos)/sqrt(1.-pow(BETA,2)); 
          projfomT->Fill(eg,CoinEvent->gammatk.gr[i].fom);
          nEgamTKT->Fill(eg,CoinEvent->gammatk.gr[i].ndet);
          NEgamTKT->Fill(eg,CoinEvent->gammatk.ngam);
          nfomT->Fill(CoinEvent->gammatk.gr[i].ndet,CoinEvent->gammatk.gr[i].fom);
          NfomT->Fill(CoinEvent->gammatk.ngam,CoinEvent->gammatk.gr[i].fom);
          if(CoinEvent->gammatk.gr[i].fom <= 0.8){
            ThetaGamTKT->Fill(eg,theta/M_PI*180.);
	    projcosTKT->Fill(eg,gccos);
            egamT[ngam]=eg;
            ngam++;
          }
        }
      }
    }
    if(ngam>1){
      for(i=0;i<ngam;i++){
        for(j=i+1;j<ngam;j++){
          ggTKP->Fill(egamP[i],egamP[j]);
          ggTKP->Fill(egamP[j],egamP[i]);
          ggTKT->Fill(egamT[i],egamT[j]);
          ggTKT->Fill(egamT[j],egamT[i]);
        }
      }
    }
  //}
}
//////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv)
{
  TFile *fRoot;
  const char *FileName,*GTCalib;
  const char *RootName,*OPTION;
  unsigned int *EventBuf=NULL;
  unsigned long long int tsEarly=0;
  unsigned long long int tsLate=0;

  unsigned long long int tsGT[MaxGTNum]={0};
  unsigned long long int tsChico[MaxChicoNum]={0};
  unsigned long long int tsTK=0;
  int j;
  int nGTRaw,nChicoRaw;

  GEBheader  Header;
  COINCIDENCE CoinEvent;

  int status;
  unsigned long EventNum=0;
  unsigned long ChicoEventNum=0;
  unsigned long GTEventNum=0;
  unsigned long TKEventNum=0;

  if (argc!=5) {
    printf("USAGE: FileName RootFilleName OPTION GTCalib\n");
    exit(EXIT_FAILURE);
  }

  printf("Data FileName: %s\n",argv[1]);
  FileName = argv[1];
  fp1 = fopen64(FileName,"rb");
  if(!fp1){
    status = errno;
    printf(" Unable to open file %s: %s\n", FileName,strerror(status));
    exit(EXIT_FAILURE);
  }

  printf("RootName: %s\n",argv[2]);
  RootName = argv[2];
  printf("OPTION: %s\n",argv[3]);
  OPTION = argv[3];
  fRoot = new TFile(RootName,OPTION);     // if adding to file
  //printf("setting up root\n");
  setuproot(fRoot);
  //printf("setting up cut\n");
  setupcut();
  printf("GTCalib File: %s\n",argv[4]);
  GTCalib = argv[4];
  //printf("ChicoInit\n");
  ChicoInit(GTCalib);
  //printf("Read CRmat\n");
  RdCRMAT("crmat.dat"); 

  nGTRaw = 0;
  nChicoRaw = 0;
  
  InitCoinEvent(&CoinEvent);
  while((EventBuf=GetGEBEvBuf(&Header,FileName))!=NULL){
    if(EventNum==0)tsFirst = Header.timestamp;
    if(Header.type == CHICOTYPE){
      ChicoEventNum++;
      ChicoDataLength->Fill(Header.length);
    }
    else if(Header.type == MODETWO)GTEventNum++;
    else if(Header.type == MODEONE){
      TKEventNum++;
      TKDataLength->Fill(Header.length);
    }

    if(Header.timestamp<tsEarly)tsEarly=Header.timestamp;
    if(Header.timestamp>tsLate)tsLate=Header.timestamp;
    if((tsLate-tsEarly)>OverLap){
      if(CoinEvent.nCoinChico>0 && CoinEvent.nCoinGT>0)ProcessEventMode2(&CoinEvent);
      if(CoinEvent.nCoinChico>0 && CoinEvent.nCoinTK>0)ProcessEventMode1(&CoinEvent);
      tsEarly=Header.timestamp;
      tsLate=Header.timestamp;
      InitCoinEvent(&CoinEvent);
      nGT_ChicoRaw->Fill(nGTRaw,nChicoRaw);
      for(j=0;j<nGTRaw;j++){
        dtGT_ChicoRaw->Fill((double)(tsChico[0]-tsGT[j]));
      }
      dtTK_ChicoRaw->Fill((double)(tsChico[0]-tsTK));
      nGTRaw = 0;
      nChicoRaw = 0;
    }
    if(Header.type == CHICOTYPE && nChicoRaw < MaxChicoNum){tsChico[nChicoRaw++]=Header.timestamp;}
    else if(Header.type == MODETWO && nGTRaw < MaxGTNum){tsGT[nGTRaw++]=Header.timestamp;}
    else if(Header.type == MODEONE){tsTK=Header.timestamp;}
    IncrementCoinEvent(&CoinEvent,&Header,EventBuf);
    free(EventBuf);
    EventNum++;
  }
  fclose(fp1);
  printf("writing histogram....\n");
  printf("Done!\n Total Event Number: %ld; GT Event Number: %ld; TK Event Number: %ld; Chico Event Number: %ld\n",
          EventNum,GTEventNum,TKEventNum,ChicoEventNum);
  fRoot->Write(NULL, TObject::kOverwrite);
  //fRoot->ls();
  fRoot->Close();
  exit(0);
}
