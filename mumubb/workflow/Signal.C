

#define Signal_cxx

#include "Signal.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "cuts.h"
#include "sphericity.h"
//#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

double dr(double eta1, double phi1, double eta2, double phi2)
{
//Double_t M_PI = 3.14159265358979323846;
    Double_t M_PI2 = 3.14159265358979323846*2.;
double deta = fabs(eta1 - eta2);
double dphi = fabs(phi1 - phi2);
 if( dphi > M_PI ) dphi = M_PI2 - dphi;
return sqrt(deta*deta + dphi*dphi);
}

void Signal::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Signal.C
//      root> Signal t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


#include "hists.icc"

bool isMC = true;
int matchedSignal = 0;
int matchedSignalAK4 = 0;
int genSignal = 0;
int counter = 0;
int ak4Events = 0;
int ak8Events = 0;
int ak4ak8Events = 0;

double centralMass = 15.0;






TFile *MyFile = new TFile("/afs/cern.ch/work/a/astepenn/haa/CMSSW_8_0_26/src/Signal_baseline_v5_trk_all.root","recreate");
RoccoR  rc("RoccoR2016.txt");


   if (fChain == 0) return;
//Double_t M_PI = 3.14159265358979323846;
    Double_t M_PI2 = 3.14159265358979323846*2.;






TLorentzVector VMu;
TLorentzVector Vb1;
TLorentzVector Vb2;
TLorentzVector VA1;
TLorentzVector VA2;
 TLorentzVector VJP;
 TLorentzVector VJPtmp;
TLorentzVector VJJ;
TLorentzVector VJJtmp;
 TLorentzVector VMuP;
  TLorentzVector VMuM;
  TLorentzVector VMuMu;
TLorentzVector VJM; 
TLorentzVector VSV;
TLorentzVector VSVs;
TLorentzVector VSVstmp;
TLorentzVector VH;
TLorentzVector VMuMuJJ;
TLorentzVector VMuMuJ;
TLorentzVector VZJ;
TLorentzVector VZSV;

TLorentzVector VSV1;
TLorentzVector VSV2;
TLorentzVector VGENP;
TLorentzVector VM1;
TLorentzVector VM2;
TLorentzVector VS1;
TLorentzVector VS2;
TLorentzVector VGenMuP;
TLorentzVector VGenMuM;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//for (Long64_t jentry=0; jentry<300000;jentry++) {  
    Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
double w = 1.0;
if(weight < 0) w = -1.0;
int jetsCounter = 0;
int fJetsCounter = 0;
//if(event != 23) continue;

if(isMC){
if(int(lumi) < 50)
{
w = w*wPu[int(lumi)];
}
}


Int_t MUMINUS=-1;
     Int_t nminus=0;
     Int_t MUPLUS=-1;
     Int_t nplus=0;
     Int_t nZgood = 0;
Double_t ptm =0;
Double_t ptp =0;

double matrix[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double matrixNorm = 0;

double ptMax = 0;
double ptMax2 = 0;

bool mu1matched = false;
bool mu2matched	= false;

//cout<<"new events:  "<<endl;
//======================================================= muons block ==========================================================


#include "muonsBlock2.icc"

trig = false;

//cout<<"triges:  "<<triggers->size()<<endl;
for(int i = 0; i < triggers->size(); i++)
{
if( ((int)((*triggers)[i]).find("HLT_IsoMu24_v") > -1 || (int)((*triggers)[i]).find("HLT_IsoTkMu24_v") > -1) &&  (VMuP.Pt() > 26 || VMuM.Pt() > 26) ) trig = true;
}
if(trig == false) continue;

//cout<<"+ and -:  "<<VMuP.Pt()<<" , "<<VMuM.Pt()<<endl;

if(VMuP.Pt() > VMuM.Pt())
{
mu1Pt_nojet->Fill(VMuP.Pt(),w);
mu2Pt_nojet->Fill(VMuM.Pt(), w);
mu1Eta_nojet->Fill(VMuP.Eta(),w);
mu2Eta_nojet->Fill(VMuM.Eta(), w);
}
else
{
mu1Pt_nojet->Fill(VMuM.Pt(),w);
mu2Pt_nojet->Fill(VMuP.Pt(), w);
mu1Eta_nojet->Fill(VMuM.Eta(),w);
mu2Eta_nojet->Fill(VMuP.Eta(), w);
}

if(M_Z > 14.8 && M_Z < 15.2)
{
if(VMuP.Pt() > VMuM.Pt())
{
mu1Pt_nojet_M->Fill(VMuP.Pt(),w);
mu2Pt_nojet_M->Fill(VMuM.Pt(), w);
mu1Eta_nojet_M->Fill(VMuP.Eta(),w);
mu2Eta_nojet_M->Fill(VMuM.Eta(), w);
}
else
{
mu1Pt_nojet_M->Fill(VMuM.Pt(),w);
mu2Pt_nojet_M->Fill(VMuP.Pt(), w);
mu1Eta_nojet_M->Fill(VMuM.Eta(),w);
mu2Eta_nojet_M->Fill(VMuP.Eta(), w);
}
}


//if(fabs(VMuMu.M() - 15.0) > 1.0) continue;
//if(pfMET > 55.0) continue;
//========================================================================================================================================
//==============================================================================================================
VA1.SetPxPyPzE(GenPx[3], GenPy[3], GenPz[3], GenE[3]);
Vb1.SetPxPyPzE(GenPx[5], GenPy[5], GenPz[5], GenE[5]);
Vb2.SetPxPyPzE(GenPx[6], GenPy[6], GenPz[6], GenE[6]);
VA2.SetPxPyPzE(GenPx[4], GenPy[4], GenPz[4], GenE[4]);
VGenMuP.SetPxPyPzE(GenPx[7], GenPy[7], GenPz[7], GenE[7]);
VGenMuM.SetPxPyPzE(GenPx[8], GenPy[8], GenPz[8], GenE[8]);

//cout<<"new event"<<endl;
//cout<<"id 7, 8:  "<<GenId[7]<<" , "<<GenId[8]<<endl;
//if(dr(VMuP.Eta(), VMuP.Phi(), VGenMuP.Eta(), VGenMuP.Phi()) < 0.3) mu1matched = true;
//if(dr(VMuM.Eta(), VMuM.Phi(), VGenMuM.Eta(), VGenMuM.Phi()) < 0.3) mu2matched =	true;

//cout<<"dr:  "<<dr(VMuM.Eta(), VMuM.Phi(), VMuP.Eta(), VMuP.Phi())<<endl;

//====================================================== SVs block ============================================
int isig1 = -1;
int isig2 = -1;

int SVcounter = 0;
int SVcounter2 = 0;
/*
for(int w_i=0; w_i<NSV; w_i++)
                        {
                        if((*sv_mass)[w_i] == 0 || (*sv_significance)[w_i] == 0) continue;
                        double sv_e = sqrt(pow((*sv_pt_x)[w_i],2) + pow((*sv_pt_y)[w_i],2) + pow((*sv_pt_z)[w_i],2) + pow((*sv_mass)[w_i],2) );
                        VSV.SetPxPyPzE((*sv_pt_x)[w_i], (*sv_pt_y)[w_i], (*sv_pt_z)[w_i], sv_e);
//                if(dr(VJPtmp.Eta(), VJPtmp.Phi(), VSV.Eta(), VSV.Phi()) > 0.8) continue;

SVcounter2++;

}
*/

//==============================================================================================================


int indMax = -1;
double vsvMass = 0;
double vsvPt = 0;
double dR = 999;
for(Int_t jj=0;jj<NumRecoJetsAK8PFCorrected; jj++) {
//============================ jet id
if(ak8PFJetsCHSValueMappt[jj] < 0) continue;
Float_t Pt1 =  JetRecoPtAK8PFCorrected[jj];
        Float_t Et1 =  JetRecoEtAK8PFCorrected[jj];
        Float_t Eta1 = JetRecoEtaAK8PFCorrected[jj];
        Float_t Phi1 = JetRecoPhiAK8PFCorrected[jj];





//VJPtmp.SetPtEtaPhiE(Pt1, Eta1, Phi1, Et1*cosh(Eta1));
VJPtmp.SetPtEtaPhiM(ak8PFJetsCHSValueMappt[jj], ak8PFJetsCHSValueMapeta[jj], ak8PFJetsCHSValueMapphi[jj], ak8PFJetsCHSValueMapmass[jj]);
//VJPtmp.SetPtEtaPhiM(ak8PFJetsCHSValueMappt[jj], ak8PFJetsCHSValueMapeta[jj], ak8PFJetsCHSValueMapphi[jj], 15.0);

if(VJPtmp.Pt() < 60.0) continue;
if(fabs(VJPtmp.Eta()) > 2.4) fJetsCounter++;
if(fabs(VJPtmp.Eta()) > 2.4) continue;
if(dr(VMuP.Eta(), VMuP.Phi(), VJPtmp.Eta(), VJPtmp.Phi()) < 0.8) continue;
if(dr(VMuM.Eta(), VMuM.Phi(), VJPtmp.Eta(), VJPtmp.Phi()) < 0.8) continue;



SVcounter = 0;
double sig1Tmp = 0;
double sig2Tmp = 0;
int isig1Tmp = -1;
int isig2Tmp = -1;
 for(int w_i=0; w_i<NSV; w_i++)
                        {
                //        if((*sv_mass)[w_i] == 0 || (*sv_significance)[w_i] == 0) continue;
        if((*sv_mass)[w_i] == 0 ) continue;             
      double sv_e = sqrt(pow((*sv_pt_x)[w_i],2) + pow((*sv_pt_y)[w_i],2) + pow((*sv_pt_z)[w_i],2) + pow((*sv_mass)[w_i],2) );
                        VSV.SetPxPyPzE((*sv_pt_x)[w_i], (*sv_pt_y)[w_i], (*sv_pt_z)[w_i], sv_e);
                if(dr(VJPtmp.Eta(), VJPtmp.Phi(), VSV.Eta(), VSV.Phi()) > 0.8) continue;
//GlobalVector flightDir((*sv_pt_x)[w_i], (*sv_pt_y)[w_i],(*sv_pt_z)[w_i]);
  //    GlobalVector jetDir(VJPtmp.Px(),VJPtmp.Py(),VJPtmp.Pz());
//if( Geom::deltaR2( flightDir, jetDir ) > 0.8 ) continue;

if(SVcounter == 0)
{
VSVs = VSV;
}
else
{
VSVs = VSVs + VSV;
}

SVcounter++;
if((*sv_mass)[w_i] > sig1Tmp)
{
sig2Tmp = sig1Tmp;
isig2Tmp = isig1Tmp;
VSV2 = VSV1;
sig1Tmp = (*sv_mass)[w_i];
isig1Tmp = w_i;
VSV1 = VSV;
}
else
{
if((*sv_mass)[w_i] > sig2Tmp)
{
VSV2 = VSV;
isig2Tmp = w_i;
sig2Tmp = (*sv_mass)[w_i];
}
}

}

if(SVcounter < 1) continue;

if(ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass[jj] < 0) ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass[jj] = -1;
if(ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass[jj] < 0) ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass[jj] = -1;  


int nSVs = 0;
double maxS = 0.0;
double maxM = 0.0;
double ntmaxS = 0.0;
double ntmaxM = 0.0;

jetsCounter++;

	if(ptMax < VJPtmp.Pt())
	{
	vsvMass = VSVs.M();
	vsvPt = VSVs.Pt();
	isig1 = isig1Tmp;
        isig2 = isig2Tmp;
	ptMax2 = ptMax;
	VJM = VJP;
	SVcounter2 = SVcounter;
	ptMax = VJPtmp.Pt();
	VJP = VJPtmp;
indMax = jj;
	}
	else
	{
	if(ptMax2 < VJPtmp.Pt())
        {
	ptMax2 = VJPtmp.Pt();
        VJM = VJPtmp;
        }
	}


}

hNJets->Fill(jetsCounter,w);
hNFJets->Fill(fJetsCounter,w);



//if(bTaggedCandidatesAK4 > 0) continue;

//if(bTaggedCandidatesAK4 > 0 && candidatesAK4 > 1)  ak4Events++;
//if(ptMax > 0)  ak8Events++;
//if(ptMax > 0 && bTaggedCandidatesAK4 > 0 && candidatesAK4 > 1 ) ak4ak8Events++;


if(!(ptMax > 0)) continue;



if(VMuP.Pt() > VMuM.Pt())
{
mu1Pt_1jet->Fill(VMuP.Pt(),w);
mu2Pt_1jet->Fill(VMuM.Pt(), w);
mu1Eta_1jet->Fill(VMuP.Eta(),w);
mu2Eta_1jet->Fill(VMuM.Eta(), w);
}
else
{
mu1Pt_1jet->Fill(VMuM.Pt(),w);
mu2Pt_1jet->Fill(VMuP.Pt(), w);
mu1Eta_1jet->Fill(VMuM.Eta(),w);
mu2Eta_1jet->Fill(VMuP.Eta(), w);
}

if(M_Z > 14.8 && M_Z < 15.2)
{
ZPt_1jet_M->Fill(VMuMu.Pt(),w);
if(VMuP.Pt() > VMuM.Pt())
{
mu1Pt_1jet_M->Fill(VMuP.Pt(),w);
mu2Pt_1jet_M->Fill(VMuM.Pt(), w);
mu1Eta_1jet_M->Fill(VMuP.Eta(),w);
mu2Eta_1jet_M->Fill(VMuM.Eta(), w);
}
else
{
mu1Pt_1jet_M->Fill(VMuM.Pt(),w);
mu2Pt_1jet_M->Fill(VMuP.Pt(), w);
mu1Eta_1jet_M->Fill(VMuM.Eta(),w);
mu2Eta_1jet_M->Fill(VMuP.Eta(), w);
}
}


//#include "firstJet.icc"

hNSV->Fill(SVcounter2,w);

//continue;
htau1->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax],w);
//htau2->Fill((*sv_mass)[isig1]/(*sv_ntrks)[isig1],w);
htau3->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau3[indMax],w);
htau4->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau4[indMax],w);

//htau2->Fill(VMuMu.() - VJP.Pt(),w);

double sv_e = sqrt(pow((*sv_pt_x)[isig1],2) + pow((*sv_pt_y)[isig1],2) + pow((*sv_pt_z)[isig1],2) + pow((*sv_mass)[isig1],2) );
VSV1.SetPxPyPzE((*sv_pt_x)[isig1], (*sv_pt_y)[isig1], (*sv_pt_z)[isig1], sv_e);


dimuon_Mass_1jet->Fill(M_Z,w);
ZPt_1jet->Fill(VMuMu.Pt(),w);
HPt_1jet->Fill(fabs(VMuMu.Pt() - VJP.Pt()),w);
ZPhi_1jet->Fill(VMuMu.Phi(), w);
hMET->Fill(pfMET,w);
hPtJ1->Fill(VJP.Pt(),w);
hPhiJ1->Fill(VJP.Phi(),w);
hEtaJ1->Fill(VJP.Eta(),w);

hdrJ1Z->Fill(dr(VJP.Eta(), VJP.Phi(), VMuMu.Eta(), VMuMu.Phi()),w);

htau24->Fill(vsvPt,w);
htau34->Fill(vsvPt/VJP.Pt(),w);

//int matchedSVcounter = 0;
VZJ = VJP + VMuMu;
VZSV = VMuMu + VSV1;
hMJ1Z->Fill(VZJ.M(),w);
hMJ1->Fill(VJP.M(), w);
hMJ1_1->Fill(ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass[indMax],w); //indMax
hMJ1_2->Fill(ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass[indMax],w);
hMSVZ->Fill(VZSV.M(), w);

if(VZJ.M() > 118 && VZJ.M() < 145)
// && pfMET < 65.0 && (*sv_mass)[isig1] > 1.2 && ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax] > 0.2 && (*sv_significance)[isig1] > 6.0)
{
hSVMN->Fill((*sv_mass)[isig1]/(*sv_ntrks)[isig1],w);
dimuon_Mass_1jet_precut_1->Fill(M_Z,w);
if((*sv_mass)[isig1] > 1.2)
{
dimuon_Mass_1jet_precut_2->Fill(M_Z,w);
if(pfMET < 65.0)
{
dimuon_Mass_1jet_precut_3->Fill(M_Z,w);
if(ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax] > 0.2)
{
dimuon_Mass_1jet_precut_4->Fill(M_Z,w);

if((*sv_significance)[isig1] > 6.0)
{
dimuon_Mass_1jet_precut_5->Fill(M_Z,w);
}

}

}

}


}

double dphiHEt = fabs(VZJ.Phi() - pfMETphi);
 if( dphiHEt > M_PI ) dphiHEt = M_PI2 - dphiHEt;

double dphiJEt = fabs(VJP.Phi() - pfMETphi);
 if( dphiJEt > M_PI ) dphiJEt = M_PI2 - dphiJEt;

double dphiZEt = fabs(VMuMu.Phi() - pfMETphi);
 if( dphiZEt > M_PI ) dphiZEt = M_PI2 - dphiZEt;

double mT = sqrt(2*pfMET*VZJ.Et()*(1 - cos(dphiHEt)));

if(M_Z > 14.8 && M_Z < 15.2 && VZJ.M() > 118 && VZJ.M() < 145 && (*sv_mass)[isig1] > 1.2) hSVNT->Fill((*sv_ntrks)[isig1],w);


hMT->Fill(mT,w);
htau14->Fill(dphiHEt,w);
//if(M_Z > 14.8 && M_Z < 15.2 && VZJ.M() > 115.0 && VZJ.M() < 165.0)
//if(M_Z > 19.75 && M_Z < 20.25 && VZJ.M() > 115.0 && VZJ.M() < 165.0)
//if(M_Z > 10.0 && M_Z < 20.0 && VZJ.M() > 110.0 && VZJ.M() < 170.0)
//if(M_Z > 14.8 && M_Z < 15.2 )
//if(M_Z > 14.8 && M_Z < 15.2 && VZJ.M() > 110.0 && VZJ.M() < 170.0)
//if(M_Z > 13 && M_Z < 17 &  VZJ.M() > 115.0 && VZJ.M() < 170     ) // <--- 15 GeV train
//if(M_Z > 14.8 && M_Z < 15.2 &  VZJ.M() > 110.0 && VZJ.M() < 170   )
//if(M_Z < 30)
//if( M_Z > 0.0)
if(M_Z > 14.0 && M_Z < 16.0   )
//if(M_Z > 14.8 && M_Z < 15.2 && VZJ.M() > 110.0 && VZJ.M() < 170   ) // <--- 15 GeV main
//if(M_Z > 19.75 && M_Z < 20.25 && VZJ.M() > 110.0 && VZJ.M() < 170   ) // <--- 20 GeV main
//if(M_Z > 24.6 && M_Z < 25.35 && VZJ.M() > 110.0 && VZJ.M() < 170   ) // <--- 25 GeV main
//if(M_Z > 29.5 && M_Z < 30.4 && VZJ.M() > 110.0 && VZJ.M() < 170   ) // <--- 30 GeV main

//if(M_Z > 19.75 && M_Z < 20.25 && VZJ.M() > 117.0  &&  VZJ.M() < 155 ) //<--- 20 test
//if(M_Z > 19.75 && M_Z < 20.25 && VZJ.M() > 100.0  &&  VZJ.M() < 200.0 ) //<--- 20 train
//if(M_Z > 24.6 && M_Z < 25.35 && VZJ.M() > 117.0  &&  VZJ.M() < 155 ) //<--- 25 test
//if(M_Z > 29.5 && M_Z < 30.35 && VZJ.M() > 117.0  &&  VZJ.M() < 155 ) //<--- 30 test
{

//hNAK4J_BDT->Fill(candidatesAK4,w);
 W = w;
SVM1 = VSV1.M();
 MZ = VMuMu.M();



DRAP = VMuMu.Rapidity() - VJP.Rapidity();
ZY_1jet->Fill(VMuMu.Rapidity() - VJP.Rapidity(), w);
if(VMuP.Pt() > VMuM.Pt())
{
PTMU1 = VMuP.Pt();
PTMu2 = VMuM.Pt(); 
mu1Pt_BDT->Fill(VMuP.Pt(),w);
mu2Pt_BDT->Fill(VMuM.Pt(),w);
}
else
{
PTMU2 =	VMuP.Pt();
PTMu1 =	VMuM.Pt();
mu2Pt_BDT->Fill(VMuP.Pt(),w);
mu1Pt_BDT->Fill(VMuM.Pt(),w);
}
ZPt_BDT->Fill(VMuMu.Pt(), w);
PTZ = VMuMu.Pt();
PTH = VZJ.E() + pfMET;
MET = pfMET;
MZSV = VZSV.M();
 MZJ = VZJ.M();
TAU1 = ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax];
PTJ1 = VJP.Pt();
SG1 = (*sv_significance)[isig1];
NTK = (*sv_ntrks)[isig1];

DPHIHET = dphiHEt;
hDPhiHEt_BDT->Fill(dphiHEt,w);

hSVM_BDT->Fill(SVM1,w);
dimuon_Mass_BDT->Fill(fabs(VMuMu.M() - 15.0),w);
hMET_BDT->Fill(pfMET,w);
hMJ1Z_BDT->Fill(fabs(VZJ.M() - 136.278),w);
htau1_BDT->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax],w);
hPtJ_BDT->Fill(VJP.Pt(),w);
hSVS_BDT->Fill((*sv_significance)[isig1],w);
hSVNT_BDT->Fill((*sv_ntrks)[isig1],w);

tVars->Fill();
}

if(isig2 > -1) hSVS2->Fill((*sv_significance)[isig2],w);
hSVM->Fill((*sv_mass)[isig1],w);
hMF->Fill(VSV1.E()/(VJP.E()),w);
//hMF->Fill(VMuMu.E()/(VMuMu.E() + VJP.E()),w);
if(isig2 > -1)
{
 sv_e = sqrt(pow((*sv_pt_x)[isig2],2) + pow((*sv_pt_y)[isig2],2) + pow((*sv_pt_z)[isig2],2) + pow((*sv_mass)[isig2],2) );
VSV2.SetPxPyPzE((*sv_pt_x)[isig2], (*sv_pt_y)[isig2], (*sv_pt_z)[isig2], sv_e);
//hMF->Fill(dr(VSV2.Eta(), VSV2.Phi(), VSV1.Eta(), VSV1.Phi()),w);
}
else
{
//hMF->Fill(-1,w);
}

if(M_Z > 14.8 && M_Z < 15.2)
//if(M_Z > 19.75 && M_Z < 20.25)
//if(M_Z > 24.6 && M_Z < 25.35)
//if(M_Z > 29.5 && M_Z < 30.35)
{
htau23->Fill(VSV1.M()/VJP.M(),w);
hDPhiZEt->Fill(dphiZEt,w);
hDPhiJEt->Fill(dphiJEt,w);
hMJ1Z_precut->Fill(VZJ.M(),w);
if(VZJ.M() > 118 && VZJ.M() < 148 )
//if(VZJ.M() > 117 && VZJ.M() < 143)
//if(VZJ.M() > 117 && VZJ.M() < 143)
//if(VZJ.M() > 115  && VZJ.M() < 136)
{
//hSVNT->Fill((*sv_ntrks)[isig1],w);
hSVM_precut->Fill((*sv_mass)[isig1],w);
if((*sv_mass)[isig1] > 1.2)
//if((*sv_mass)[isig1] > 1.3)
//if((*sv_mass)[isig1] > 1.0)
//if((*sv_mass)[isig1] > 0.8)
{
//if(pfMET < 50.0)
if(pfMET < 65.0)
//if(pfMET < 50.0)
{
if(ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax] > 0.2)
//if(ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax] > 0.3)
{
hSVS_precut->Fill((*sv_significance)[isig1],w);
}
htau1_precut->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[indMax],w);
htau2_precut->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau2[indMax],w);
htau3_precut->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau3[indMax],w);
htau4_precut->Fill(ak8PFJetsCHSValueMapNjettinessAK8CHSTau4[indMax],w);
hMF_precut->Fill(VMuMu.E()/(VMuMu.E() + VJP.E()),w);
hMT_precut->Fill(mT,w);
}
htau14_precut->Fill(dphiHEt,w);
hMET_precut->Fill(pfMET,w);
}
}
}

}



mu1Pt_BDT->Write();
mu2Pt_BDT->Write();
ZPt_BDT->Write();
hSVM_BDT->Write();
dimuon_Mass_BDT->Write();
hMET_BDT->Write();
hMJ1Z_BDT->Write();
htau1_BDT->Write();
hPtJ_BDT->Write();
hSVS_BDT->Write();
hSVNT_BDT->Write();



//cout<<"counter:  "<<counter<<endl;
//cout<<"integral:  "<<hNSV->Integral(3,10)<<endl;
dimuon_Mass->Write();
hNJets->Write();
hPtJ1->Write();
hPhiJ1->Write();
hEtaJ1->Write();
dimuon_Mass_1jet->Write();
dimuon_Mass_1jet_matched->Write();
ZPt_1jet->Write();
ZPt_1jet_M->Write();
ZY_1jet->Write();
ZPhi_1jet->Write();
hNSV->Write();
hMET->Write();
hMET_precut->Write();
hMJ1Z->Write();
hSVS->Write();
hSVS2->Write();
hMF->Write();
hSVMs->Write();
hSVM_precut->Write();
hSVM->Write();;
hSVM2->Write();


sv1Eta_nojet->Write();
sv2Eta_nojet->Write();
sv1Pt_nojet->Write();
sv2Pt_nojet->Write();
svsEta_nojet->Write();
svsPt_nojet->Write();
svsPhi_nojet->Write();
hDRJ1VV->Write();
dimuon_Mass_2v->Write();

hMJ1Z_precut->Write();
hMF_precut->Write();
hMJ1->Write();
dxy1->Write();
dz1->Write();
dxy2->Write();
dz2->Write();
muPId->Write();
matchingIds->Write();
matchingIds2->Write();


mu1Pt_nojet_M->Write();
mu2Pt_nojet_M->Write();
mu1Eta_nojet_M->Write();
mu2Eta_nojet_M->Write();
mu1Pt_nojet->Write();
mu2Pt_nojet->Write();
mu1Eta_nojet->Write();
mu2Eta_nojet->Write();


mu1Pt_1jet_M->Write();
mu2Pt_1jet_M->Write();
mu1Eta_1jet_M->Write();
mu2Eta_1jet_M->Write();
mu1Pt_1jet->Write();
mu2Pt_1jet->Write();
mu1Eta_1jet->Write();
mu2Eta_1jet->Write();



mu2Pt_mismatch->Write();
ZPt_postcut->Write();
matchingIds_postcut->Write();

hdr1->Write();
hdr2->Write();
hdr1_lhe->Write();
hdr2_lhe->Write();
hMult->Write();

hEnergies->Write();

htau1->Write();
htau2->Write();
htau3->Write();
htau4->Write();
htau1_precut->Write();
htau2_precut->Write();
htau3_precut->Write();
htau4_precut->Write();

hDRJ1V2->Write();
hDRJ1V1->Write();
hdrJ1Z->Write();

hMJ1_1->Write();
hMJ1_2->Write();

hDPhiHEt_BDT->Write();

htau12->Write();
htau13->Write();
htau14->Write();
htau23->Write();
htau24->Write();
htau34->Write();
hNFJets->Write();
hDPhiZEt->Write();
hDPhiJEt->Write();
htau14_precut->Write();
hMT->Write();
hMT_precut->Write();
tVars->Write();
//cout<<"ak4, ak8, both:  "<<ak4Events<<" , "<<ak8Events<<" , "<<ak4ak8Events<<endl;
hSVS_precut->Write();
//dimuon_Mass_1jet_precut->Write();
dimuon_Mass_1jet_precut_1->Write();
dimuon_Mass_1jet_precut_2->Write();
dimuon_Mass_1jet_precut_3->Write();
dimuon_Mass_1jet_precut_4->Write();
dimuon_Mass_1jet_precut_5->Write();
dimuon_Mass_1jet_precut_6->Write();
hSVNT->Write();
hSVMN->Write();
hNAK4J->Write();
hNAK4J_BDT->Write();
ZY_1jet->Write();
HPt_1jet->Write();
hMSVZ->Write();
cout<<"number of entries for Signal:  "<<tVars->GetEntriesFast()<<endl;
if(tVars->GetEntriesFast() < 1) cout<<"Signal has no entries"<<endl;
}



