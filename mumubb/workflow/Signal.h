//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  1 21:50:43 2013 by ROOT version 5.32/00
// from TChain /
//////////////////////////////////////////////////////////

#ifndef Signal_h
#define Signal_h

#include <iostream>
using namespace std;

#include "TLorentzVector.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "RoccoR.cc"


class Signal {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types

Float_t weight;

Float_t             lumi;
   Int_t           event;
   Int_t           run;
bool trig = false;

Int_t NSV = 0;
std::vector<double> * sv_x = 0;
std::vector<double> * sv_y = 0;
std::vector<double> * sv_z = 0;
std::vector<double> * sv_pt_x = 0;
std::vector<double> * sv_pt_y = 0;
std::vector<double> * sv_pt_z = 0;
std::vector<double> * sv_significance = 0;
std::vector<double> *  sv_value = 0;
std::vector<double> * sv_error = 0;
std::vector<double> * sv_mass = 0;
std::vector<double> * sv_mass_corr = 0;
std::vector<double> * sv_pt = 0;
std::vector<int> * sv_ntrks = 0;

Int_t           NumGenP;
 Float_t        GenM[20];
Float_t        GenPx[20];
Float_t        GenPy[20];
Float_t        GenPz[20];
Float_t        GenE[20];
Int_t        GenId[20];
Int_t        GenStatus[20];

std::vector<std::string> * triggers = 0;

double pfMET;
        double pfMETphi;

  Int_t           NumRecoJetsAK8PFCorrected;
   Float_t         JetRecoPtAK8PFCorrected[100];
   Float_t         JetRecoEtAK8PFCorrected[100];
   Float_t         JetRecoEtaAK8PFCorrected[100];
   Float_t         JetRecoMassAK8PFCorrected[100];
   Float_t         JetRecoPhiAK8PFCorrected[100];
Float_t         NjettinessAK8Puppitau1[15];
 Float_t         NjettinessAK8Puppitau2[15];
 Float_t         NjettinessAK8Puppitau3[15];
 Float_t         NjettinessAK8Puppitau4[15];
 Float_t         ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[15];
 Float_t         ak8PFJetsCHSValueMapNjettinessAK8CHSTau2[15];
 Float_t         ak8PFJetsCHSValueMapNjettinessAK8CHSTau3[15];
 Float_t         ak8PFJetsCHSValueMapNjettinessAK8CHSTau4[15];
 Float_t         ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass[15];
 Float_t         ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass[15];
 Float_t         ak8PFJetsCHSValueMapeta[15];
 Float_t         ak8PFJetsCHSValueMapmass[15];
 Float_t         ak8PFJetsCHSValueMapphi[15];
 Float_t         ak8PFJetsCHSValueMappt[15];
 Float_t         ak8PFJetsPuppiSoftDropMass[15];
 Float_t         ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2[15];
 Float_t         ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3[15];
 Float_t         ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2[15];
 Float_t         ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3[15];

  Int_t           NumRecoJetsPFCorrected;
   Float_t         JetRecoPtPFCorrected[100];
   Float_t         JetRecoEtPFCorrected[100];
   Float_t         JetRecoEtaPFCorrected[100];
   Float_t         JetRecoRapidityPFCorrected[100];
   Float_t         JetRecoPhiPFCorrected[100];
   Float_t        JetRecoUNCPFCorrected[100];
   Float_t          JetRecoPUIDPFCorrected[100];

Float_t         DiscMVA[100];
Float_t         DiscCSV[100];
Float_t         DiscCvB[100];
Float_t         DiscCvL[100];

Float_t         DiscDeepCSVb[100];
Float_t         DiscDeepCSVbb[100];
Float_t         DiscDeepCSVc[100];
Float_t         DiscDeepCSVudsg[100];


Float_t JetRecoSVPt[100];
Float_t JetRecoSVMassCorrectedPFCorrected[100];
Float_t         JetRecoNHFPFCorrected[100];
Float_t         JetRecoNEMFPFCorrected[100];
Float_t         JetRecoCHFPFCorrected[100] ;
Float_t         JetRecoMUFPFCorrected[100] ;
Float_t         JetRecoCEMFPFCorrected[100];
Int_t           JetRecoNumConstPFCorrected[100] ;
Int_t           JetRecoNumNeutralParticlesPFCorrected[100] ;
Float_t         JetRecoCHMPFCorrected[100] ;
Float_t         JetRecoDiscLPFCorrected[100];
Float_t         JetRecoDiscBPFCorrected[100];
Int_t         JetRecoPFlavPFCorrected[100];
Int_t         JetRecoHFlavPFCorrected[100];

// TBranch        *run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_NvtxEv;   //!

 TBranch *b_pfMET;
TBranch *b_pfMET_corr;
  TBranch *b_pfMETphi; 



Float_t         RecoPtMuons[25];
   Float_t         RecoEtaMuons[25];
   Float_t         RecoPhiMuons[25];
   Int_t         RecoChargeMuons[25];
Int_t           NumRecoMuons;
Int_t           NumLayers[25];
   Float_t         RecoIsoMuonsDxyvtx[25];
   Float_t         RecoIsoMuonsDzvtx[25];

float trkIsolation[25];
float pfIsolation[25];

bool isLoose[25];
bool isMedium[25];
bool isTight[25];


   TBranch        *b_NumRecoJetsPFCorrected;   //!
   TBranch        *b_JetRecoPtPFCorrected;   //!
   TBranch        *b_JetRecoEtPFCorrected;   //!
   TBranch        *b_JetRecoEtaPFCorrected;   //!
   TBranch        *b_JetRecoRapidityPFCorrected;   //!
   TBranch        *b_JetRecoPhiPFCorrected;   //!
 TBranch *b_lumi;
 TBranch *b_JetRecoUNCPFCorrected; 

TBranch         *b_genpid;
TBranch         *b_genppt;
TBranch        *b_genpeta;
TBranch         *b_genpphi;



   TBranch        *b_NumRecoMuons;   //!
   TBranch        *b_RecoPtMuons;   //!
   TBranch        *b_RecoEtaMuons;   //!
   TBranch        *b_RecoPhiMuons;   //!
   TBranch        *b_RecoChargeMuons;   //!
   TBranch        *b_RecoIsoMuonsSumPt;   //!
   //   TBranch        *b_RecoImpactDxyMuons;   //!
   //   TBranch        *b_RecoImpactDzMuons;   //!
   //   TBranch        *b_RecoImpactDxyErrMuons;   //!
   //   TBranch        *b_RecoImpactDzErrMuons;   //!


/*
TBranch         *b_SoftPtMuons;
  TBranch         *b_SoftEtaMuons;
  TBranch         *b_SoftPhiMuons;
  TBranch        *b_SoftChargeMuons;
TBranch         *b_SoftIsoMuonsSumPt;
*/






   Signal(TTree *tree=0);
   virtual ~Signal();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Signal_cxx
Signal::Signal(TTree *tree)
{
    cout<<"test1"<<endl;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f) {
         f = new TFile("Memory Directory");
         f->cd("Rint:/");
      }
      tree = (TTree*)gDirectory->Get("MuEG");
       cout <<"single tree mode"<<endl;

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("MuEG","");
 cout <<"chain  mode"<<endl;

////////////////////////////      mc:
//afs/cern.ch/work/a/astepenn/haa4b/CMSSW_9_4_13/src/test/Select/crab_SUSYgg20
//chain->Add("/afs/cern.ch/work/a/astepenn/haa4b/CMSSW_9_4_13/src/test/Select/crab_SUSYggM15_2/results/default-filename*.root/wztree");

//chain->Add("/mnt/sdd1/stepennov/haa/crab_Signal15/results/default-filename*.root/wztree");
//chain->Add("/eos/user/a/astepenn/haa/crab_Signal15/results/default-filename*.root/wztree");
//chain->Add("/afs/cern.ch/work/a/astepenn/haa4b/CMSSW_9_4_13/src/test/Select/crab_Signal15/results/default-filename*.root/wztree");

chain->Add("/afs/cern.ch/work/a/astepenn/haa/CMSSW_8_0_26/src/Signal_preselection_v5_trk.root/wztree");

//chain->Add("/mnt/sdb1/stepennov/haa/crab_DY/results/default-filename*.root/wztree;1");
//chain->Add("/mnt/sdb1/stepennov/haa/crab_DY_2/results/default-filename*.root/wztree;1");

//chain->Add("/afs/cern.ch/work/a/astepenn/haa4b/CMSSW_9_4_13/src/test/Select/default-filename.root/wztree");

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Signal::~Signal()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Signal::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Signal::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Signal::Init(TTree *tree)
{
    cout<<"test2"<<endl;

   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   // General 






fChain->SetBranchAddress("lumi", &lumi);
fChain->SetBranchAddress("weight",       &weight);

fChain->SetBranchAddress("pfMET",       &pfMET);
fChain->SetBranchAddress("pfMETphi",       &pfMETphi);

fChain->SetBranchAddress("trig",       &trig);
fChain->SetBranchAddress("run",       &run);
   fChain->SetBranchAddress("event",     &event);

fChain->SetBranchAddress("NSV", &NSV);


fChain->SetBranchAddress("sv_x", &sv_x);
fChain->SetBranchAddress("sv_y", &sv_y);
fChain->SetBranchAddress("sv_z", &sv_z);
fChain->SetBranchAddress("sv_pt_x", &sv_pt_x);
fChain->SetBranchAddress("sv_pt_y", &sv_pt_y);
fChain->SetBranchAddress("sv_pt_z", &sv_pt_z);
fChain->SetBranchAddress("sv_significance", &sv_significance);
fChain->SetBranchAddress("sv_value", &sv_value);
fChain->SetBranchAddress("sv_error", &sv_error);
fChain->SetBranchAddress("sv_mass", &sv_mass);
fChain->SetBranchAddress("sv_mass_corr", &sv_mass_corr);
fChain->SetBranchAddress("sv_pt", &sv_pt);
fChain->SetBranchAddress("sv_ntrks", &sv_ntrks);
/*
fChain->SetBranchAddress("NumGenP", &NumGenP);
fChain->SetBranchAddress("GenM", GenM);
fChain->SetBranchAddress("GenPx", GenPx);
fChain->SetBranchAddress("GenPy", GenPy);
fChain->SetBranchAddress("GenPz", GenPz);
fChain->SetBranchAddress("GenE", GenE);
fChain->SetBranchAddress("GenId", GenId);
fChain->SetBranchAddress("GenStatus", GenStatus);
*/

fChain->SetBranchAddress("NumRecoMuons", &NumRecoMuons);
fChain->SetBranchAddress("RecoPtMuons", RecoPtMuons);
fChain->SetBranchAddress("RecoEtaMuons", RecoEtaMuons);
fChain->SetBranchAddress("RecoPhiMuons", RecoPhiMuons);
fChain->SetBranchAddress("RecoChargeMuons", RecoChargeMuons);
fChain->SetBranchAddress("RecoIsoMuonsDxyvtx", RecoIsoMuonsDxyvtx);
fChain->SetBranchAddress("RecoIsoMuonsDzvtx", RecoIsoMuonsDzvtx);
fChain->SetBranchAddress("isLoose", isLoose);
fChain->SetBranchAddress("isMedium", isMedium);
fChain->SetBranchAddress("isTight", isTight);
fChain->SetBranchAddress("NumLayers", NumLayers);
fChain->SetBranchAddress("trkIsolation", trkIsolation);
fChain->SetBranchAddress("pfIsolation", pfIsolation);

fChain->SetBranchAddress("triggers", &triggers);

fChain->SetBranchAddress("NumRecoJetsAK8PFCorrected", &NumRecoJetsAK8PFCorrected);
    fChain->SetBranchAddress("JetRecoPtAK8PFCorrected", JetRecoPtAK8PFCorrected);
    fChain->SetBranchAddress("JetRecoEtAK8PFCorrected", JetRecoEtAK8PFCorrected);
    fChain->SetBranchAddress("JetRecoEtaAK8PFCorrected", JetRecoEtaAK8PFCorrected);
    fChain->SetBranchAddress("JetRecoMassAK8PFCorrected", JetRecoMassAK8PFCorrected);
    fChain->SetBranchAddress("JetRecoPhiAK8PFCorrected", JetRecoPhiAK8PFCorrected);

fChain->SetBranchAddress("NjettinessAK8Puppitau1", NjettinessAK8Puppitau1);
fChain->SetBranchAddress("NjettinessAK8Puppitau2", NjettinessAK8Puppitau2);
fChain->SetBranchAddress("NjettinessAK8Puppitau3", NjettinessAK8Puppitau3);
fChain->SetBranchAddress("NjettinessAK8Puppitau4", NjettinessAK8Puppitau4);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapNjettinessAK8CHSTau1", ak8PFJetsCHSValueMapNjettinessAK8CHSTau1);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapNjettinessAK8CHSTau2", ak8PFJetsCHSValueMapNjettinessAK8CHSTau2);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapNjettinessAK8CHSTau3", ak8PFJetsCHSValueMapNjettinessAK8CHSTau3);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapNjettinessAK8CHSTau4", ak8PFJetsCHSValueMapNjettinessAK8CHSTau4);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass", ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass", ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapeta", ak8PFJetsCHSValueMapeta);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapmass", ak8PFJetsCHSValueMapmass);
fChain->SetBranchAddress("ak8PFJetsCHSValueMapphi", ak8PFJetsCHSValueMapphi);
fChain->SetBranchAddress("ak8PFJetsCHSValueMappt", ak8PFJetsCHSValueMappt);
fChain->SetBranchAddress("ak8PFJetsPuppiSoftDropMass", ak8PFJetsPuppiSoftDropMass);
fChain->SetBranchAddress("ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2", ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2);
fChain->SetBranchAddress("ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3", ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3);
fChain->SetBranchAddress("ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2", ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2);
fChain->SetBranchAddress("ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3", ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3);




fChain->SetBranchAddress("NumRecoJetsPFCorrected", &NumRecoJetsPFCorrected);
    fChain->SetBranchAddress("JetRecoPtPFCorrected", JetRecoPtPFCorrected);
    fChain->SetBranchAddress("JetRecoEtPFCorrected", JetRecoEtPFCorrected);
    fChain->SetBranchAddress("JetRecoEtaPFCorrected", JetRecoEtaPFCorrected);
    fChain->SetBranchAddress("JetRecoRapidityPFCorrected", JetRecoRapidityPFCorrected);
    fChain->SetBranchAddress("JetRecoPhiPFCorrected", JetRecoPhiPFCorrected);
fChain->SetBranchAddress("JetRecoUNCPFCorrected", JetRecoUNCPFCorrected);
fChain->SetBranchAddress("JetRecoPUIDPFCorrected", JetRecoPUIDPFCorrected);
fChain->SetBranchAddress("JetRecoHFlavPFCorrected", JetRecoHFlavPFCorrected);

fChain->SetBranchAddress("JetRecoNHFPFCorrected", JetRecoNHFPFCorrected);
fChain->SetBranchAddress("JetRecoNEMFPFCorrected", JetRecoNEMFPFCorrected);
fChain->SetBranchAddress("JetRecoCHFPFCorrected", JetRecoCHFPFCorrected);
fChain->SetBranchAddress("JetRecoMUFPFCorrected", JetRecoMUFPFCorrected);
fChain->SetBranchAddress("JetRecoCEMFPFCorrected", JetRecoCEMFPFCorrected);
fChain->SetBranchAddress("JetRecoNumConstPFCorrected", JetRecoNumConstPFCorrected);
fChain->SetBranchAddress("JetRecoNumNeutralParticlesPFCorrected", JetRecoNumNeutralParticlesPFCorrected);
fChain->SetBranchAddress("JetRecoCHMPFCorrected", JetRecoCHMPFCorrected);


fChain->SetBranchAddress("DiscMVA", DiscMVA);
fChain->SetBranchAddress("DiscCSV", DiscCSV);
fChain->SetBranchAddress("DiscCvB", DiscCvB);
fChain->SetBranchAddress("DiscCvL", DiscCvL);
fChain->SetBranchAddress("DiscDeepCSVb", DiscDeepCSVb);
fChain->SetBranchAddress("DiscDeepCSVbb", DiscDeepCSVbb);
fChain->SetBranchAddress("DiscDeepCSVc", DiscDeepCSVc);
fChain->SetBranchAddress("DiscDeepCSVudsg", DiscDeepCSVudsg);


   Notify();

}

Bool_t Signal::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Signal::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Signal::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Signalw_cxx



