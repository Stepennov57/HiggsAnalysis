#include <TError.h>
#include <TMath.h>


#include <TCanvas.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TF1.h>
#include <TStyle.h>
#include <TVector.h>
#include <TGraph.h>
#include "TMatrixD.h"
#include "TUnfoldDensity.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>

#include <TF1.h>

#include <TFile.h>

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"

// #define VERBOSE_LCURVE_SCAN

using namespace std;

void test()
{

 TCanvas *c1 = new TCanvas("c1","comparison",100,10,1100,800);
    TPaveText *pt = new TPaveText(.05,.1,.95,.8);
    c1->SetWindowSize(700, 800);

 gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

TFile *f2 = new TFile("unfolding_v58E_tightCtag_xs_genj.root");
//TFile *f2 = new TFile("/Users/Anton/closure_v62.root");


double Ptbinning[5] =  {30, 45, 60,  90, 250};
//double Ptbinning[5] =  {0, 30, 50, 95, 300};

double binning2[14] = {40, 42, 45, 50, 56, 66, 81, 104, 139, 190, 267, 382, 555, 640};
int N_gen = 4;

TH2D* response = (TH2D*)f2->Get("response_w_0"); // responce matrix (MC)
TH1D* Z_pt_gen_matched = (TH1D*)f2->Get("acceptance2_w_0"); // pt distribution of gen Z, matched with measured reco Z (MC)
TH1D* Z_pt_reco_matched = (TH1D*)f2->Get("background2_w_0"); // pt distribution of measured reco Z, matched with gen Z (MC)
TH1D* Z_pt_gen = (TH1D*)f2->Get("acceptance1_w_0"); // pt distribution of gen Z (MC)
TH1D* Z_pt_reco = (TH1D*)f2->Get("background1_w_0"); // pt distribution of measured reco Z (MC)


TH1D* hTrue = (TH1D*)f2->Get("acceptance1_w_0");
TH1D* hTrue2 = (TH1D*)f2->Get("acceptance1_w_0");
TH1D* hMeas = (TH1D*)f2->Get("dataToUnfold_w_0");
TH1D* hClos  = (TH1D*)f2->Get("background1_w_0");



 TH1D *bkgN = new TH1D(*Z_pt_reco_matched);
bkgN->Divide(bkgN, Z_pt_reco , 1.0, 1.0, "B");

//hMeas->Multiply(bkgN);

// multiplying distribution, which is about to be unfolded, by fiducial purity
    
hClos->Multiply(bkgN);

cout<<"Here1"<<endl;
TUnfoldDensity unfold(response,TUnfold::kHistMapOutputVert);
//TUnfold unfold(response,TUnfold::kHistMapOutputVert,                 TUnfold::kRegModeNone);
cout<<"Here2"<<endl;

TUnfold::ERegMode regMode=TUnfold::kRegModeSize ;
//TUnfold::ERegMode regMode=TUnfold::kRegModeCurvature    ;
//unfold.RegularizeBins(1,1,7,regMode);

// if(unfold.SetInput(hMeas)>=10000) {
    
//unfolding MC "measured" distrubution. The result should be the same as original distribution at generator level
if(unfold.SetInput(hClos)>=10000) {
  std::cout<<"Unfolding result may be wrong\n";
  }


cout<<"Here3"<<endl;
//========================================================================
  // the unfolding is done here
  //
  // scan L curve and find best point
  Int_t nScan=60;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;

#ifdef VERBOSE_LCURVE_SCAN
  Int_t oldinfo=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kInfo;
#endif
//  iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
//cout<<"Here4"<<endl;
unfoldingReturn=unfold.DoUnfold(0.0);

#ifdef VERBOSE_LCURVE_SCAN
  gErrorIgnoreLevel=oldinfo;
#endif

double binning[6] = {30, 50,  70, 100,  130, 300};
TH1 *histResult = unfold.GetOutput("Unfolded");

TH2 *histEmatTotal=unfold.GetEmatrixTotal("EmatTotal");

double covError = 0;
for(int i = 1; i < 5; i++)
{
for(int j = 1; j < 5; j++)
{
covError += histEmatTotal->GetBinContent(i, j);
}
}
cout<<"!!!!!!!!!1  "<<covError<<endl;

    
// dividing unfolded result by acceptance
TH1D *accN = new TH1D(*Z_pt_gen_matched);
accN->Divide(accN, Z_pt_gen , 1.0, 1.0, "B");
  histResult->Divide(accN);

    
    hTrue->SetLineColor(kRed);
hTrue->Draw("HIST");
histResult->Draw("PEsame");



}



