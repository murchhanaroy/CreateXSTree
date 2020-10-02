//this code is used to plot RCS physics
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include <iomanip>
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

const int kNZbins=1;
const double kZLength=40.0;
const double kCosThBinWidth=0.0005;
const double kPhiBinWidth=0.01;
const double kPBinWidth=0.005;
//conclusion: The bin width here will affect the integrated geometry acceptance by
//several percent. For example, if all bin width are doubled, the acceptance will increase
//by ~1%.
//const double kCosThBinWidth=0.001;
//const double kPhiBinWidth=0.02;
//const double kPBinWidth=0.01;

static double gBeam=0.0;

TH3F *h3ETP0[kNZbins];       //store thrown counts
TH3F *h3ETP_XS[kNZbins][3];  //store detected counts, weighted by XS, 1 for H2, 2 for N2, 3 for 3He 

void GetTH1BoundaryValues(TH1 *h1,double &pXStart, double &pXEnd, double entriescut=-1)
{
  // Get the start bin and end bin 
  int pNXBin=h1->GetNbinsX();
  int pXBinStart=0, pXBinEnd=-1;

  if(entriescut<=0) entriescut=0.1*h1->GetMaximum(); 
  for(int i=1;i<=pNXBin;i++) {
      //continuous 4 bins none zero
    if(pXBinStart<=0 && i+3<=pNXBin) {
      if (h1->GetBinContent(i)>entriescut   && h1->GetBinContent(i+1)>0 &&
        h1->GetBinContent(i+2)>0 && h1->GetBinContent(i+3)>0 )
      pXBinStart=i;
    }  

    if(pXBinEnd<=0 && pNXBin-i-2>=1) {
      if (h1->GetBinContent(pNXBin-i+1)>entriescut && h1->GetBinContent(pNXBin-i)>0 &&
        h1->GetBinContent(pNXBin-i-1)>0 && h1->GetBinContent(pNXBin-i-2)>0)
      pXBinEnd=pNXBin-i+1;
    }
      if(pXBinStart>0 && pXBinEnd>0) break;
  }
  if(pXBinEnd<=0) pXBinEnd=pNXBin;

  pXStart=h1->GetBinLowEdge(pXBinStart);
  pXEnd=h1->GetBinLowEdge(pXBinEnd)+h1->GetBinWidth(pXBinEnd);
  cout<<setw(15)<<h1->GetName()
      <<":  Xmin="<<setw(10)<<pXStart
      <<"  Xmax="<<setw(10)<<pXEnd
      <<"  HalfWidth="<<setw(10)<<(pXEnd-pXStart)/2<<endl;

  double pYmax=h1->GetMaximum();
  TLine *L1=new TLine(pXStart,0,pXStart,pYmax);
  L1->SetLineColor(8);
  L1->SetLineWidth(2);
  h1->GetListOfFunctions()->Add(L1); 
  TLine *L2=new TLine(pXEnd,0,pXEnd,pYmax);
  h1->GetListOfFunctions()->Add(L2); 
  L2->SetLineColor(8);
  L2->SetLineWidth(2);
}

//declare histo
void FillHisto(const char *extracut="1"){
  TTree * T=(TTree*) gROOT->FindObject("T");
  //find the limit of costheta and phi
  int NxP=50,NxCosTh=50,NxPhi=50;
  double Pmin,Pmax,CosThmin,CosThmax,Phimin,Phimax;
  TH1 *h1;
  
  T->Draw("p0>>h1p",Form("p0<%.4f && istop==0",gBeam));  
  h1=(TH1F*) gROOT->FindObject("h1p");
  GetTH1BoundaryValues(h1,Pmin,Pmax,5);
  NxP = (int) ceil((Pmax-Pmin)/kPBinWidth);
  
  T->Draw("cos(theta0)>>h1costh",Form("p0<%.4f && istop==0",gBeam));  
  h1=(TH1F*) gROOT->FindObject("h1costh");
  GetTH1BoundaryValues(h1,CosThmin,CosThmax,5);
  NxCosTh = (int) ceil((CosThmax-CosThmin)/kCosThBinWidth);
  
  T->Draw("phi0>>h1phi",Form("p0<%.4f && istop==0",gBeam));  
  h1=(TH1F*) gROOT->FindObject("h1phi");
  GetTH1BoundaryValues(h1,Phimin,Phimax,5);
  NxPhi = (int) ceil((Phimax-Phimin)/kPhiBinWidth);

  cout<<"\n ================== binning: "<<NxP<<" x "<<NxCosTh<<" x "<<NxPhi<<" = "<<NxP*NxCosTh*NxPhi<<" ==================\n"; 

  //create histo
  for(int i=0;i<kNZbins;i++) {
    h3ETP0[i]=new TH3F(Form("h3ETP0_Z%02d",i),"Thrown Conut",
       NxP,Pmin,Pmax,NxCosTh,CosThmin,CosThmax,NxPhi,Phimin,Phimax);
    h3ETP_XS[i][0]=new TH3F(Form("h3ETP_1H_Z%02d",i),"Yield, Weighted by XS",
       NxP,Pmin,Pmax,NxCosTh,CosThmin,CosThmax,NxPhi,Phimin,Phimax);
    h3ETP_XS[i][1]=new TH3F(Form("h3ETP_14N_Z%02d",i),"Yield, Weighted by XS",
       NxP,Pmin,Pmax,NxCosTh,CosThmin,CosThmax,NxPhi,Phimin,Phimax);
    h3ETP_XS[i][2]=new TH3F(Form("h3ETP_3He_Z%02d",i),"Yield, Weighted by XS",
       NxP,Pmin,Pmax,NxCosTh,CosThmin,CosThmax,NxPhi,Phimin,Phimax);
  }

  //fill histo
  char cut0[255];
  char cut[3][255];
  double dZ = kZLength/kNZbins;
  for(int i=0;i<kNZbins;i++) {
    double zz=-kZLength/2.0 + (i+0.5)*dZ;
    sprintf(cut0,"abs(vz0-(%.3f)) <= %.3f && p0>=%.4f && p0<=%.4f && cos(theta0)>=%.4f && cos(theta0)<=%.4f && phi0>=%.4f && phi0<=%.4f",
      zz,dZ/2,Pmin,Pmax,CosThmin,CosThmax,Phimin,Phimax);
    sprintf(cut[0],"( (%s) && (%s) && (istop==0) && xs_1h>0) * xs_1h",cut0,extracut);
    //sprintf(cut[0],"( (%s) && (%s) && (istop==0)) * 1.0",cut0,extracut);  //debug only
    sprintf(cut[1],"( (%s) && (%s) && (istop==0) && xs_3he>0) * xs_14n",cut0,extracut);
    sprintf(cut[2],"( (%s) && (%s) && (istop==0) && xs_14n>0) * xs_3he",cut0,extracut);
    
    cout<<" char* cut = \""<<cut0<<"\";"<<endl;
    cout<<" char* cut0 = \""<<cut[0]<<"\";"<<endl;
    cout<<" char* cut1 = \""<<cut[1]<<"\";"<<endl;
    cout<<" char* cut2 = \""<<cut[2]<<"\";"<<endl;

    T->Project(Form("h3ETP0_Z%02d",i),"phi0:cos(theta0):p0",cut0);
    T->Project(Form("h3ETP_1H_Z%02d",i),"phi0:cos(theta0):p0",cut[0]);
    T->Project(Form("h3ETP_14N_Z%02d",i),"phi0:cos(theta0):p0",cut[1]);
    T->Project(Form("h3ETP_3He_Z%02d",i),"phi0:cos(theta0):p0",cut[2]);
    
    //cout<<h3ETP0[i]->GetName()<<": enntires = "<<h3ETP0[i]->GetEntries()<<",   cut0 = \""<<cut0<<"\""<<endl;
    //cout<<h3ETP_XS[i][0]->GetName()<<": enntires = "<<h3ETP_XS[i][0]->GetEntries()<<endl;
    //cout<<h3ETP_XS[i][1]->GetName()<<": enntires = "<<h3ETP_XS[i][1]->GetEntries()<<endl;
    //cout<<h3ETP_XS[i][2]->GetName()<<": enntires = "<<h3ETP_XS[i][2]->GetEntries()<<endl;
  }
}


//get integated XS for 40cm long target
void GetInteXS(double* pInteXS)
{ 
  TH3F *h3=0;
  for(int j=0;j<3;j++) {
    pInteXS[j]=0.0;
    for(int i=0;i<kNZbins;i++) {
      h3 = (TH3F*) h3ETP_XS[i][j]->Clone();
      h3->Divide(h3ETP0[i]);
      double pInteXS_Z = h3->Integral("width");
      pInteXS[j] += pInteXS_Z;
      cout<<setw(14)<<h3ETP_XS[i][j]->GetName()<<": N0= "<<setw(12)<<h3ETP0[i]->GetEntries()<<", N="
          <<setw(12)<<h3ETP_XS[i][j]->GetEntries()<<", this_InteXS_Z= "<<setw(12)<<pInteXS_Z*1.0E6<<" (pb),  <InteXS>= "
          <<setw(12)<<pInteXS[j]*1.0E6/(i+1)<<" (pb)"<<endl;
    }
     pInteXS[j] /= kNZbins;  //take the average
     pInteXS[j] *= 1000.0;  //turn into nb
  }
}


//get rate for method 2
void GetRate(double pBeamCurrent_uA, double pBeamE_GeV, double pDetectorAngle_deg, double pDetectorMomentum_GeV, int Det, int pFullAcc)
{
  if(pFullAcc) cout<<"\n================ Method2: Full Acceptance ================\n";
  else         cout<<"\n================ Method2: 2-SC-bar only ==================\n";
  
  gBeam=pBeamE_GeV;
  const char* DetName[]= {"HMS","SHMS"};
  const char* TgName[]= {"H2","N2","3He"};
  double I_uA = pBeamCurrent_uA;
  double Lumi_1uA[3]={124.55, 124.55, 62.27};  //in unit of 10^33
  double pRate[3];

  double pInteXS[3]={0.0,0.0,0.0};  //in unit of nb
  
  if(pFullAcc) FillHisto("1");
  else FillHisto("abs(xfp)<8");
  GetInteXS(pInteXS); 
  
  //00000"00000008|0000007|0000007|0000007|0000007|00000008|000007|000000009|00000000001|00000000901|0000000090|"
  printf("Detector  Target   BeamE Current   DetP0 DetAngle VZ(cm) Thick(cm) Lumi(10^33)  InteXS(pb)   Rate(Hz)\n");
  for(int j=0;j<3;j++) {
    pRate[j] = I_uA * Lumi_1uA[j] * pInteXS[j];  //in Hz

    printf("%8s %7s %7.4f %7.2f %7.4f %7.2f %7.2f %9.2f %11.2f %11.2f %10.2f\n",
      DetName[Det-1],TgName[j],pBeamE_GeV,I_uA,pDetectorMomentum_GeV,pDetectorAngle_deg,0.0,40.0,
      I_uA*Lumi_1uA[j],pInteXS[j]*1000,pRate[j]);
  }
}


void GetRate(int det=1, int fullAcc=1)
{
  if(det==1) GetRate(1.0,2.1,11.7,2.068,1,fullAcc);
  if(det==2) GetRate(1.0,2.1,8.5,2.083,2,fullAcc);
}

void GetRate(const char *filename,int det=1, int fullAcc=1)
{
  TFile *file = TFile::Open(filename);
  if(det==1) GetRate(1.0,2.1,11.7,2.068,1,fullAcc);
  if(det==2) GetRate(1.0,2.1,8.5,2.083,2,fullAcc);
  file->Close();
}


