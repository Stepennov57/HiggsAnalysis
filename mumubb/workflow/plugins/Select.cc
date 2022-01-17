// -*- C++ -*-

//





// Package:    test/Select
// Class:      Select
// 
/**\class Select Select.cc test/Select/plugins/Select.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Anton Stepennov
//         Created:  Fri, 22 Feb 2019 16:19:11 GMT
//
//

#include <TFile.h>
#include "TProfile.h"
#include "TTree.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>


#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTSuperLayer.h"
#include "DataFormats/DTRecHit/interface/DTSLRecSegment2D.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentUpdator.h"
#include "RecoLocalMuon/DTSegment/src/DTSegmentCleaner.h"
#include "RecoLocalMuon/DTSegment/src/DTHitPairForFit.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSecondaryVertexTagInfo.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/JetReco/interface/Jet.h"

using namespace std;
 using namespace edm;

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TLorentzVector.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Select : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Select(const edm::ParameterSet&);
      ~Select();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
edm::EDGetTokenT<pat::MuonCollection> tok_muon_;
edm::EDGetTokenT<reco::VertexCollection> tok_vrtx_;
   edm::EDGetTokenT<TriggerResults> tok_trig_;
edm::EDGetTokenT<pat::ElectronCollection> tok_ele_;

edm::EDGetTokenT<reco::GenJetCollection> genjet_tag_;
edm::EDGetTokenT<reco::GenParticleCollection> genpTok;

edm::EDGetTokenT<reco::GenJetCollection> genpTok2;
edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;

edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
edm::EDGetTokenT<GenEventInfoProduct> genI_tok_;
edm::EDGetTokenT<std::vector<PileupSummaryInfo > >  _puSummaryInfo;

edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> tok_sv_;
edm::EDGetTokenT<edm::View<reco::Vertex>> tok_pv_;


//edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> GenJetFlavourInfosToken_;

edm::EDGetTokenT<pat::JetCollection>  tok_ak8jet_;

edm::EDGetTokenT<pat::JetCollection>  tok_jet_;
//edm::EDGetTokenT<pat::JetCollection>
edm::EDGetTokenT<edm::View<pat::Jet > > tok_jet_updated_;

edm::EDGetTokenT<double>            tok_Rho_;

edm::EDGetTokenT<pat::METCollection> tok_met_;

edm::EDGetTokenT< double > prefweight_token;
edm::EDGetTokenT< double > prefweightup_token;
edm::EDGetTokenT< double > prefweightdown_token;


TLorentzVector   VJP ;
TLorentzVector   VMuP ;
TLorentzVector  VMuM  ;
TLorentzVector  VJJ  ;

TLorentzVector   VMuP2 ;
TLorentzVector  VMuM2  ;
TLorentzVector  VMuMu ;

TLorentzVector  VMuMuLHE ;
TLorentzVector   VMuPLHE ;
TLorentzVector  VMuMLHE  ;

TLorentzVector  VMuMuGEN ;
TLorentzVector   VMuPGEN ;
TLorentzVector  VMuMGEN  ;

TLorentzVector  VMuMuGEN2 ;
TLorentzVector   VMuPGEN2 ;
TLorentzVector  VMuMGEN2  ;

bool tightCtag;
bool mediumBtag;

bool trig;
double weight2;
   Int_t           event;

UInt_t           event32u;
Long64_t           event64;
ULong64_t           event64u;

Float_t weightPS;
Float_t weightPSUp;
Float_t weightPSDown;

Float_t weight;
Float_t         genppt;
Float_t         genpptP;
Float_t         genpptM;
Int_t           run;
Int_t truecounter;
Int_t mode1;
Int_t mode2;

//   Float_t         weight;
   Float_t         VtxZ;
   Float_t         VtxRho;
   Int_t           NvtxEv;
Float_t             lumi;
Int_t lumisection;

Double_t Rho;

Int_t         NumRecoJetsAK8PFCorrected;
 Float_t         JetRecoPtAK8PFCorrected[15];
   Float_t         JetRecoEtAK8PFCorrected[15];
   Float_t         JetRecoEtaAK8PFCorrected[15];
   Float_t         JetRecoPhiAK8PFCorrected[15];
	Float_t         JetRecoMassAK8PFCorrected[15];

 Float_t         pfBoostedDoubleSecondaryVertexAK8BJetTags[15];
 Float_t         pfMassIndependentDeepDoubleBvLJetTags[15];

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
   Float_t         JetRecoPtPFCorrected[30];
   Float_t         JetRecoEtPFCorrected[30];
   Float_t         JetRecoEtaPFCorrected[30];
   Float_t         JetRecoRapidityPFCorrected[30];
   Float_t         JetRecoPhiPFCorrected[30];
   Float_t         JetRecoPUMVAPFCorrected[30];
   Float_t         JetRecoPUIDPFCorrected[30];
Float_t         JetRecoPUIDNEWPFCorrected[30];

Float_t         JetRecoNHFPFCorrected[30];
Float_t         JetRecoNEMFPFCorrected[30];
Float_t         JetRecoCHFPFCorrected[30] ;
Float_t         JetRecoMUFPFCorrected[30] ;
Float_t         JetRecoCEMFPFCorrected[30];
Int_t           JetRecoNumConstPFCorrected[30] ;
Int_t           JetRecoNumNeutralParticlesPFCorrected[30] ;
Float_t         JetRecoCHMPFCorrected[30] ;

std::vector<float> oldpfpt;
std::vector<float> oldpfeta;
std::vector<float> oldpfphi;
std::vector<int> oldpfid;
Int_t oldpfnum;
std::vector<std::string> triggers;

std::vector<float> generatorpt;
std::vector<float> generatoreta;
std::vector<float> generatorphi;
std::vector<float> generatore;
std::vector<int> generatorid;
std::vector<int> generatormotherid;
std::vector<int> generatorstatus;
Int_t NumGeneratorP;


Float_t         JetRecoPUIDCBPFCorrected[30];
Float_t         JetRecoUNCPFCorrected[30];

Float_t         DiscMVA[30];
Float_t         DiscCSV[30];
Float_t         DiscCvB[30];
Float_t         DiscCvL[30];

Float_t         DiscDeepCSVb[30];
Float_t         DiscDeepCSVbb[30];
Float_t         DiscDeepCSVc[30];
Float_t         DiscDeepCSVudsg[30];

Float_t         DiscDeepJetb[30];
Float_t         DiscDeepJetbb[30];
Float_t         DiscDeepJetlepb[30];

Float_t JetRecoSVMassAVRPFCorrected[30];
Float_t JetRecoSVMassIVFPFCorrected[30];
Float_t JetRecoSVMassCorrectedPFCorrected[30];
Float_t JetRecoSVPt[30];
Int_t JetRecoSVN[30];

Int_t         JetRecoPFlavPFCorrected[30];
Int_t         JetRecoHFlavPFCorrected[30];


int muplus;
int muminus;

Int_t           NumRecoElectrons;
Float_t ElectronPt[25];
Float_t ElectronEta[25];
Float_t ElectronPhi[25];
Float_t ElectronEnergy[25];
Int_t   ElectronCharge[25];

Float_t ElectronIdLoose[25];
Float_t ElectronIdMedium[25];
Float_t ElectronIdTight[25];

Float_t  trackMomentumAtVtx[25];
Float_t  ecalEnergy[25];


Float_t full5x5_sigmaIetaIeta[25];
Float_t ooEmooP[25];
Float_t dEtaIn[25];
Float_t dPhiIn[25];
Float_t hOverE[25];
Float_t sigmaIetaIeta[25];
Float_t dxy[25];
Float_t dz[25];
Float_t relIsoWithEA[25];
Float_t eMissingInnerHits[25];
bool passCVeto[25];
 Float_t mva[25];
Float_t Esc[25];
Float_t Etasc[25];


Int_t NumTrigObj;
Float_t TrigEta[20];
Float_t TrigPhi[20];

Float_t         RecoPtMuons[25];
   Float_t         RecoEtaMuons[25];
   Float_t         RecoPhiMuons[25];
   Int_t         RecoChargeMuons[25];
Int_t           NumRecoMuons;
Int_t           NumLayers[25];
   Float_t         RecoIsoMuonsDxyvtx[25];
   Float_t         RecoIsoMuonsDzvtx[25];

 Float_t         RecoIsoMuonsDxyvtx2[25];
   Float_t         RecoIsoMuonsDzvtx2[25];

Float_t         RecoIsoMuonsIsol[25];
   Float_t         RecoIsoMuonsSumPt[25];

Float_t         RecoIsoMuonsSumChargedHpt[25];
Float_t         RecoIsoMuonsSumNeutralHpt[25];
Float_t         RecoIsoMuonsPhoton[25];
Float_t         RecoIsoMuonsSumPUPt[25];

   Float_t         RecoIsoMuonsSumPtPV[25];

float ValidFraction[25];
bool isLoose[25];
bool isMedium[25];
bool isTight[25];
bool isGlobal[25];
float NChi[25];
int numberOfValidMuonHits[25];
int numberOfMatchedStations[25];
int numberOfValidPixelHits[25];
int trackerLayersWithMeasurement[25];
float ChiLocPos[25];
float trkKink[25];
float segmentCompatibility[25];
float trkIsolation[25];
float pfIsolation[25];

Int_t           NumGenJets;
   Float_t         PtGenJets[30];
   Float_t         EtGenJets[30];
   Float_t         EtaGenJets[30];
   Float_t         PhiGenJets[30];
Float_t         PFlavGenJets[30];
Float_t         HFlavGenJets[30];

Int_t           NumGenPFinal;
Int_t           NumGenP;
 Float_t        GenM[20];
Float_t        GenPx[20];
Float_t        GenPy[20];
Float_t        GenPz[20];
Float_t        GenE[20];
Int_t        GenId[20];
Int_t        GenStatus[20];

std::vector<double> evtweight;
std::vector<string> evtwid;

Int_t        NSV;
std::vector<double> sv_x;
std::vector<double> sv_y;
std::vector<double> sv_z;
std::vector<double> sv_pt_x;
std::vector<double> sv_pt_y;
std::vector<double> sv_pt_z;
std::vector<double> sv_significance;
std::vector<double> sv_value;
std::vector<double> sv_error;
std::vector<double> sv_mass;
std::vector<double> sv_mass_corr;
std::vector<double> sv_pt;
std::vector<int> sv_ntrks;



Int_t NumGenMu;
Float_t GenMuPt[25];
Float_t GenMuEta[25];
Float_t GenMuPhi[25];
Float_t GenMuE[25];
Int_t GenMuStatus[25];
Int_t GenMuId[25];

Int_t NumGenMu2;
Float_t GenMuPt2[25];
Float_t GenMuEta2[25];
Float_t GenMuPhi2[25];
Float_t GenMuE2[25];
Int_t GenMuStatus2[25];
Int_t GenMuId2[25];

double pfMET;
	double pfMETphi;
double pfMET_corr;
        double pfMETphi_corr;

int numberOfSVertices;
float sVertexX[25];
float sVertexY[25]; 
float sVertexZ[25]; 
float sVertexEta[25]; 
float sVertexPhi[25]; 
float sVertexMass[25]; 
float sVertexMassCorrected[25];

TTree *wwtree;
        TTree *wwtree2;
      TFile *f;
std::string fileName;
   virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Select::Select(const edm::ParameterSet& iConfig)

{
fileName = iConfig.getUntrackedParameter<std::string>("fileName","default-filename.root");
   //now do what ever initialization is needed
usesResource("TFileService");
tok_muon_ = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"));
tok_vrtx_ = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
tok_trig_  = consumes<TriggerResults>(edm::InputTag("TriggerResults","","HLT"));
tok_ele_ = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));

prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));

genjet_tag_                    = consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"));
genpTok = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
genpTok2 = consumes<reco::GenJetCollection>(edm::InputTag("particleLevel", "leptons"));


triggerObjects_  = consumes<std::vector<pat::TriggerObjectStandAlone>>(edm::InputTag("slimmedPatTrigger","","PAT"));

generatorlheToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer","")) ;
genI_tok_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
_puSummaryInfo                  = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));


//jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("jetFlavourInfos") );
//GenJetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( iConfig.getParameter<edm::InputTag>("flavourMap") );


tok_Rho_ = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));

tok_sv_ = consumes<edm::View<reco::VertexCompositePtrCandidate>>(edm::InputTag("slimmedSecondaryVertices"));
tok_pv_ = consumes<edm::View<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));


tok_ak8jet_            = consumes<pat::JetCollection>(edm::InputTag("slimmedJetsAK8"));
tok_jet_            = consumes<pat::JetCollection>(edm::InputTag("slimmedJets"));
//tok_jet_updated_                = consumes<edm::View<pat::Jet > >(iConfig.getParameter<edm::InputTag>("updatedPatJetsUpdatedJEC"));

tok_met_ = consumes<pat::METCollection>(edm::InputTag("slimmedMETs"));

   //now do what ever initialization is needed
   usesResource("TFileService");

}


Select::~Select()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

double testF()
{
return 1.0;
}



double SVM(const pat::Jet* jet)
{
double mass = -1;
 double en = 0.;
      double px = 0.;
      double py = 0.;
      double pz = 0.;
const reco::CandSecondaryVertexTagInfo *candSVTagInfo = jet->tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
if( candSVTagInfo!=nullptr) {
        if ( candSVTagInfo->nVertices() >= 1 ) {

//TrackKinematics kin(candSVTagInfo->secondaryVertex(0));
//JetRecoSVMassPFCorrected[NumRecoJetsPFCorrected]  =  kin.vectorSum().M();	
//cout<<"mass:   "<<kin.vectorSum().M()<<endl;
//cout<<"sv eta:  "<<candSVTagInfo->secondaryVertex(0).eta();
//return candSVTagInfo->secondaryVertex(0).eta();
double ntracks = candSVTagInfo->nVertexTracks(0);
for (unsigned int ktrack =0 ;ktrack < ntracks; ++ktrack){
	double m_hadron = 0.13957;
	m_hadron = 0.13957;
	  double ptsigned = candSVTagInfo->vertexTracks(0)[ktrack]->charge() * candSVTagInfo->vertexTracks(0)[ktrack]->pt();
	  double eta_track = candSVTagInfo->vertexTracks(0)[ktrack]->eta();
	  double phi_track = candSVTagInfo->vertexTracks(0)[ktrack]->phi();
	  double pt_track = fabs(ptsigned);
	  double px_track = pt_track*cos(phi_track);
	  double py_track = pt_track*sin(phi_track);
	  double pz_track = pt_track*sinh(eta_track);
	  en += sqrt(m_hadron*m_hadron+pt_track*pt_track+pz_track*pz_track);
	  px += px_track;
	  py += py_track;
	  pz += pz_track;

	}
mass = sqrt(en*en - px*px - py*py - pz*pz);

}
}
return mass;
}

//
// member functions
//

// ------------ method called for each event  ------------
void
Select::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

//cout<<"test:  "<<testF()<<endl;
Rho = 0;
truecounter = 0; 
weight = 1;
trig = false;
tightCtag = false;
mediumBtag = false;
run = iEvent.id().run();
//run = 1;
event = iEvent.id().event();
event32u = iEvent.id().event();
event64 = iEvent.id().event();
event64u = iEvent.id().event();
lumisection = iEvent.luminosityBlock();


edm::Handle< double > theprefweight;
iEvent.getByToken(prefweight_token, theprefweight ) ;
double _prefiringweight =(*theprefweight);

edm::Handle< double > theprefweightup;
iEvent.getByToken(prefweightup_token, theprefweightup ) ;
double _prefiringweightup =(*theprefweightup);

edm::Handle< double > theprefweightdown;
iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
double _prefiringweightdown =(*theprefweightdown);

//cout<<"weight:  "<<_prefiringweightup<<endl;
weightPS = _prefiringweight;
weightPSUp = _prefiringweightup;
weightPSDown = _prefiringweightdown;

edm::Handle<GenEventInfoProduct> genEvt;
      iEvent.getByToken(genI_tok_,genEvt);
      int procid = (int) genEvt->signalProcessID();
    weight = genEvt->weight();


NSV = 0;
sv_ntrks = {};
sv_x = {};
sv_y = {};
sv_z = {};
sv_pt_x = {};
 sv_pt_y = {};
 sv_pt_z = {};
 sv_significance = {};
 sv_value = {};
 sv_error = {};
 sv_mass = {};
 sv_mass_corr = {};
 sv_pt = {};

generatorpt = {};
generatoreta = {};
generatorphi = {};
generatore = {};
generatorid = {};
generatormotherid = {};
generatorstatus = {};



evtweight = {};
evtwid = {};
//cout<<"main weight:  "<<weight<<endl;
evtweight = genEvt->weights();
cout<<"size:   "<<evtweight.size()<<endl;
for(int t = 0; t< 107; t++) //111
{
//cout<<evtweight[t]<<endl;
}

//cout<<"genweight   "<<weight<<endl;
 edm::Handle<std::vector< PileupSummaryInfo > > infoPU;
      iEvent.getByToken(_puSummaryInfo,infoPU);

      for(std::vector<PileupSummaryInfo>::const_iterator it = infoPU->begin(); it != infoPU->end(); it++)
        {
          if(it->getBunchCrossing() == 0) lumi = it->getTrueNumInteractions();
        }



Double_t mupluspt=0;
Double_t muminuspt=0;
muplus=0;
muminus=0;
int muplus25 = 0;
int muminus25 = 0;
trig = false;

edm::Handle<double> rh;
  iEvent.getByToken(tok_Rho_,rh);
const double rho_val = *(rh.product());
//cout<<"rho:   "<<rho_val<<endl;
Rho = rho_val;


edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);
JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);


//==================================== triggers ===================================
edm::Handle<edm::TriggerResults> triggerResults;
//edm::InputTag trigResultsTag("TriggerResults","","HLT");
//iEvent.getByLabel(trigResultsTag,triggerResults);

 Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

triggers = {};

iEvent.getByToken(tok_trig_, triggerResults);
const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);

 for(unsigned i=0; i<triggerNames.size(); i++) {
if(triggerResults->accept(i)) {
//std::cout<<", "<<triggerNames.triggerName(i)<<std::endl;
//if((int)(triggerNames.triggerName(i)).find("HLT_IsoMu24_v")>-1) trig = true;
//if((int)(triggerNames.triggerName(i)).find("HLT_IsoTkMu24")>-1) trig = true;
//if((int)(triggerNames.triggerName(i)).find("HLT_Ele27_eta2p1_WPTight_Gsf")>-1) trigE = true;
// if((int)(triggerNames.triggerName(i)).find("HLT_")>-1) std::cout<<"shoot!!  "<<triggerNames.triggerName(i)<<std::endl;;
//if((int)(triggerNames.triggerName(i)).find("HLT_Mu")>-1) std::cout<<"here1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<triggerNames.triggerName(i)<<std::endl;
//std::cout<<"contains:  "<<(int)(triggerNames.triggerName(i)).find("HLT")<<std::endl;
triggers.push_back(triggerNames.triggerName(i));
//std::string::size_type pos = s.find("text");

//if(contains(triggerNames.triggerName(i), s2) == true) std::cout<<"here1  "<<triggerNames.triggerName(i)<<std::endl;
           // if((int)(triggerNames.triggerName(i)).find("HLT_IsoMu")>-1) std::cout<<"shoot!!"<<std::endl;
               // if((int)(triggerNames.triggerName(i)).find("HLT_IsoMu24_eta2p1")>-1) {HLT_IsoMu24_eta2p1=1; }
}
}

NumTrigObj = 0;
/*
std::vector<std::string> filterLabels_;
for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(triggerNames);
        obj.unpackFilterLabels(iEvent, *triggerResults);
//	std::cout << "\t   Filters:    ";w
//        std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.hasFilterLabel("hltEle27WPTightGsfTrackIsoFilter") << std::endl;
//	cout<<"trigger name: "
if(obj.hasFilterLabel("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") || obj.hasFilterLabel("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"))
{
cout<<"name0:  "<<obj.filterLabels()[0]<<endl;
}

if(obj.hasFilterLabel("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") || obj.hasFilterLabel("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09"))
{
TrigEta[NumTrigObj] = obj.eta();
TrigPhi[NumTrigObj] = obj.phi();
NumTrigObj++;
}

}
*/

//===================================== genp collection ============================


Handle<reco::GenParticleCollection> particles;
   iEvent.getByToken(genpTok, particles);

Handle<reco::GenJetCollection> particles2;
   iEvent.getByToken(genpTok2, particles2);

NumGenMu = 0;
NumGenMu2 = 0;
double genptp = 0;
double genptm = 0;
double genptp2 = 0;
double genptm2 = 0;

double ptm = 0;
double ptp = 0;

double MZGEN = 0;
double MZGEN2 = 0;

 NumGenP = 0;
 NumGenPFinal = 0;
NumGeneratorP = 0;


edm::Handle<LHEEventProduct> lheEventProduct;
iEvent.getByToken(generatorlheToken_, lheEventProduct);

// =============================================== gen particles ====================================

for ( reco::GenParticleCollection::const_iterator genParticle = particles->begin();
        genParticle != particles->end(); ++genParticle ) {
    int absPdgId = TMath::Abs(genParticle->pdgId());
    int status = genParticle->status();
//if(!(absPdgId == 25)) continue;
//if(absPdgId == 5) cout<<"status:  "<<status<<" id: "<<genParticle->pdgId()<<" pt: "<<genParticle->pt()<<" ,parent:  "<<genParticle->mother()->pdgId()<<endl;
NumGeneratorP++;
int motherid = 0;
if(genParticle->mother() != NULL) motherid = genParticle->mother()->pdgId();
generatorpt.push_back(genParticle->pt());
generatoreta.push_back(genParticle->eta());
generatorphi.push_back(genParticle->phi());
generatore.push_back(genParticle->energy());
generatorid.push_back(genParticle->pdgId());
generatormotherid.push_back(motherid);
generatorstatus.push_back(genParticle->status());
}

// =============================================== lhe particles ====================================================================


const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
                        std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
                        for ( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle ) {
                        //	int id = std::abs(lheEvent.IDUP[idxParticle]);
                        //	int status = lheEvent.ISTUP[idxParticle];
                        //	double mass = lheEvent.PUP[idxParticle][4];

GenId[NumGenP] = lheEvent.IDUP[idxParticle];
GenStatus[NumGenP] = lheEvent.ISTUP[idxParticle];
GenPx[NumGenP] = lheEvent.PUP[idxParticle][0];
GenPy[NumGenP] = lheEvent.PUP[idxParticle][1];
GenPz[NumGenP] = lheEvent.PUP[idxParticle][2];
GenE[NumGenP] = lheEvent.PUP[idxParticle][3];
GenM[NumGenP] = lheEvent.PUP[idxParticle][4];

//cout<<"status, id, num, idx:  "<<lheEvent.ISTUP[idxParticle]<<" , "<<lheEvent.IDUP[idxParticle]<<" , "<<NumGenP<<" , "<<idxParticle<<" , "<<lheEvent.PUP[idxParticle][0]<<" , "<<lheEvent.PUP[idxParticle][1]<<" , "<<lheEvent.PUP[idxParticle][2]<<endl;
if(lheEvent.ISTUP[idxParticle] == 1) NumGenPFinal++;
//cout<<"mothers:  "<<lheEvent.MOTHUP[idxParticle].first<<" , "<<lheEvent.MOTHUP[idxParticle].second<<endl;


VMuP.SetPxPyPzE(GenPx[NumGenP], GenPy[NumGenP], GenPz[NumGenP], GenE[NumGenP]);
//if(fabs(VMuP.M() - GenM[NumGenP]) > 2.0)
//{
//cout<<"A pt, eta, phi, M1, M2:  "<<VMuP.Pt()<<" , "<<VMuP.Eta()<<" , "<<VMuP.Phi()<<" , "<<VMuP.M()<<" , "<<GenM[NumGenP]<<endl;
//cout<<"px, py, pz, E, M:   "<<GenPx[NumGenP]<<" , "<<GenPy[NumGenP]<<" , "<<GenPz[NumGenP]<<" , "<<GenE[NumGenP]<<" , "<<GenM[NumGenP]<<endl;
//}

                                NumGenP++;

                }


evtweight = {};
evtwid = {};

  std::vector<gen::WeightsInfo> weights = lheEventProduct->weights();
  //if(weights.size() != 0) {evtweight->clear();}
  for( unsigned int i = 0; i < weights.size(); i++) {
    gen::WeightsInfo lhew = weights[i];
   evtweight.push_back(lhew.wgt);
        evtwid.push_back(lhew.id);
  //   cout <<"  weights -> i = " << lhew.id <<"  weight = " << lhew.wgt  << endl;
  }




//===================================================================================== genjets ========================================================


//================================================================================================================

//=================================================================== missing Et ===================================================================================================$
//=================================================================== missing Et ===================================================================================================$
edm::Handle<pat::METCollection> pfMEThandle;
iEvent.getByToken(tok_met_, pfMEThandle);
pfMET	   = (pfMEThandle->front() ).et();
  pfMETphi = (pfMEThandle->front() ).phi();
//cout<<"missing et:   "<<pfMET<<endl;

/*
edm::Handle< edm::View<reco::PFMET> > pfMEThandle2;
iEvent.getByToken(tok_met_corr_, pfMEThandle2);
pfMET_corr      = (pfMEThandle2->front() ).et();
  pfMETphi_corr = (pfMEThandle2->front() ).phi();
//cout<<"missing et, corr:   "<<pfMET<<" , "<<pfMET_corr<<endl;
*/

//========================================= muons ===============================


//========================================= muons ===============================

edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByToken(tok_vrtx_,pvHandle);
    const reco::VertexCollection & vertices = *pvHandle.product();
edm::Handle<pat::MuonCollection> muonColl;
iEvent.getByToken(tok_muon_, muonColl);
reco::Vertex vertx = *(vertices.begin());

NvtxEv = 0;
reco::Vertex::Point vert;
for(reco::VertexCollection::const_iterator it = vertices.begin() ; it != vertices.end() ; ++it)  {
      if(!(*it).isFake() && (*it).ndof() > 4  ) {
        if(it->tracksSize() > 0 &&  ( fabs(it->z()) <= 24. ) &&  ( fabs(it->position().rho()) <= 2. ) ) {
//        result = true;
          NvtxEv++;
          if(it == vertices.begin()) {
        vert = (*it).position();
        }
        }
      } // nonfake vertex
    }



NumRecoMuons=0;
int numMuP = 0;
int numMuM = 0;
for (pat::MuonCollection::const_iterator muon=muonColl->begin(), muonCollEnd=muonColl->end(); muon!=muonCollEnd;
++muon)
{

if(muon->pt() > 5.0)
{
        //reco::TrackRef trkTrack = muon->track();
        RecoPtMuons[NumRecoMuons] = muon->pt();
        RecoEtaMuons[NumRecoMuons] = muon->eta();
        RecoPhiMuons[NumRecoMuons] = muon->phi();
isLoose[NumRecoMuons] = muon::isLooseMuon(*muon);
isMedium[NumRecoMuons] = muon::isMediumMuon(*muon);
isGlobal[NumRecoMuons] = muon->isGlobalMuon();
isTight[NumRecoMuons] = muon::isTightMuon(*muon, vertx);

        RecoChargeMuons[NumRecoMuons] = muon->charge();

if(muon->charge() > 0)
{
numMuP++;
}
else
{
numMuM++;
}

 if( !(muon->muonBestTrack().isNull()))
{
        RecoIsoMuonsDxyvtx[NumRecoMuons] =  fabs(muon->muonBestTrack()->dxy(vertx.position()));
 RecoIsoMuonsDzvtx[NumRecoMuons] =  fabs(muon->muonBestTrack()->dz(vertx.position()));
}
else
{
RecoIsoMuonsDxyvtx[NumRecoMuons] = -1;
RecoIsoMuonsDzvtx[NumRecoMuons] = -1;
}

if(muon->isTrackerMuon() || muon->isGlobalMuon())
{
NumLayers[NumRecoMuons] = muon->innerTrack()->hitPattern().trackerLayersWithMeasurement();
}
else
{
NumLayers[NumRecoMuons] = -1;
}

trkIsolation[NumRecoMuons] = muon->isolationR03().sumPt;
pfIsolation[NumRecoMuons] = (muon->pfIsolationR04().sumChargedHadronPt + fmax(0, muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - 0.5*muon->pfIsolationR04().sumPUPt));


  NumRecoMuons++;

}
}

//============================================= jet ========================================================
NumRecoJetsPFCorrected = 0;
NumRecoJetsAK8PFCorrected = 0;

//edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
  //iEvent.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );



edm::Handle<pat::JetCollection> ak8pfjetH;
iEvent.getByToken(tok_ak8jet_, ak8pfjetH);

edm::Handle<pat::JetCollection> pfjetH;
iEvent.getByToken(tok_jet_, pfjetH);

edm::Handle<edm::View<reco::VertexCompositePtrCandidate>> svs_;
iEvent.getByToken(tok_sv_, svs_);

edm::Handle<edm::View<reco::Vertex>> pvs_;
iEvent.getByToken(tok_pv_, pvs_);
const auto & pv = (*pvs_)[0];

//======================================================== sv block ==========================================

double ptSV = -1;
double massIVF = -1;
double m_vertex_ptcorr = -1;

int nSV= 0;

for(const auto &sv: *svs_){
nSV++;

VertexDistance3D vdist;
double en = 0.;
      double px = 0.;
      double py = 0.;
      double pz = 0.;
double m_hadron = 0.13957;
reco::TrackKinematics kin(sv);

int nTracks = 0;

for(size_t i=0; i < sv.numberOfSourceCandidatePtrs(); ++i)
{

  en += sqrt(m_hadron*m_hadron+(sv.daughterPtr(i)->pt()*sv.daughterPtr(i)->pt())+(sv.daughterPtr(i)->pz()*sv.daughterPtr(i)->pz()));
          px += sv.daughterPtr(i)->px();
          py += sv.daughterPtr(i)->py();
          pz += sv.daughterPtr(i)->pz();
        nTracks++;
}

GlobalVector flightDir(sv.vertex().x() - pv.x(), sv.vertex().y() - pv.y(),sv.vertex().z() - pv.z());
Measurement1D dl= vdist.distance(pv,VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));

sv_ntrks.push_back(nTracks);
sv_x.push_back(sv.vertex().x() - pv.x());
sv_y.push_back(sv.vertex().y() - pv.y());
sv_z.push_back(sv.vertex().z() - pv.z());
sv_pt_x.push_back(sv.vertex().x() - pv.x());
sv_pt_y.push_back(sv.vertex().y() - pv.y());
sv_pt_z.push_back(sv.vertex().z() - pv.z());
sv_significance.push_back(dl.significance());
sv_value.push_back(dl.value());
sv_error.push_back(dl.error());

double massIVF = sv.p4().M();
double m_vertex = sqrt(en*en - px*px - py*py - pz*pz);
double m_vertex2 = (en*en - px*px - py*py - pz*pz);
double ptSV = sqrt(px*px + py*py);
double alpha = (px*flightDir.x() + py*flightDir.y())/(sqrt(px*px+py*py)*sqrt(flightDir.x()*flightDir.x()+flightDir.y()*flightDir.y()));
 alpha = acos(alpha);
double ptcorr = sqrt(px*px+py*py)*sin(alpha);
m_vertex_ptcorr = sqrt(m_vertex2+ptcorr*ptcorr)+ptcorr;

sv_mass.push_back(massIVF);
sv_mass_corr.push_back(m_vertex_ptcorr);
sv_pt.push_back(px*px + py*py);

}
NSV = nSV;

//================================================================================================================================================================================================================
//======================================================== jets block =======================================================

for (unsigned ijet = 0; ijet < ak8pfjetH->size(); ++ijet) {
const pat::Jet* jet = dynamic_cast<const pat::Jet*>(&(*ak8pfjetH)[ijet]);
//cout<<"ak8 jet "<<jet->pt()<<endl;
//cout<<"ak8 jet id:  "<<jet->userFloat("pileupJetId:fullId")<<endl;
//cout<<"other mass"<<jet->userFloat("asdfa")<<endl;
//cout<<"ak8 pt, eta, phi, mass:  "<<jet->pt()<<" , "<<jet->eta()<<" , "<<jet->phi()<<" , "<<jet->mass()<<" , "<<jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")<<" , "<<jet->userFloat("ak8PFJetsCHSValueMap:mass")<<endl;

//cout<<jet->userInt("asdfs")<<endl;


//cout<<"mass1, mass2, mass3:  "<<jet->mass()<<" , "<<jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")<<" , "<<jet->userFloat("ak8PFJetsCHSValueMap:mass")<<endl;

         JetRecoPtAK8PFCorrected[NumRecoJetsAK8PFCorrected] = jet->pt();
         JetRecoEtAK8PFCorrected[NumRecoJetsAK8PFCorrected] = jet->et();
         JetRecoEtaAK8PFCorrected[NumRecoJetsAK8PFCorrected] = jet->eta();
         JetRecoPhiAK8PFCorrected[NumRecoJetsAK8PFCorrected] = jet->phi();
	JetRecoMassAK8PFCorrected[NumRecoJetsAK8PFCorrected] = jet->mass();


pfBoostedDoubleSecondaryVertexAK8BJetTags[NumRecoJetsAK8PFCorrected] = jet->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
//pfMassIndependentDeepDoubleBvLJetTags[NumRecoJetsAK8PFCorrected] = jet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb");


//cout<<"discr:  "<<jet->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags")<<" , "<<jet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb")<<endl;
NjettinessAK8Puppitau1[NumRecoJetsAK8PFCorrected] = jet->userFloat("NjettinessAK8Puppi:tau1");
NjettinessAK8Puppitau2[NumRecoJetsAK8PFCorrected] = jet->userFloat("NjettinessAK8Puppi:tau2");
NjettinessAK8Puppitau3[NumRecoJetsAK8PFCorrected] = jet->userFloat("NjettinessAK8Puppi:tau3");
NjettinessAK8Puppitau4[NumRecoJetsAK8PFCorrected] = jet->userFloat("NjettinessAK8Puppi:tau4");
ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1");
ak8PFJetsCHSValueMapNjettinessAK8CHSTau2[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2");
ak8PFJetsCHSValueMapNjettinessAK8CHSTau3[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3");
ak8PFJetsCHSValueMapNjettinessAK8CHSTau4[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau4");
ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass");
ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass");
ak8PFJetsCHSValueMapeta[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:eta");
ak8PFJetsCHSValueMapmass[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:mass");
ak8PFJetsCHSValueMapphi[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:phi");
ak8PFJetsCHSValueMappt[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsCHSValueMap:pt");
ak8PFJetsPuppiSoftDropMass[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsPuppiSoftDropMass");
ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2");
ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3");
ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2");
ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3[NumRecoJetsAK8PFCorrected] = jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3");

NumRecoJetsAK8PFCorrected++;
}

reco::JetRefBaseProd jetref(pfjetH);

for (unsigned ijet = 0; ijet < pfjetH->size(); ++ijet) {

const pat::Jet* jet = dynamic_cast<const pat::Jet*>(&(*pfjetH)[ijet]);

float mva   = jet->userFloat("pileupJetId:fullDiscriminant");
int    idflag = jet->userInt("pileupJetId:fullId");


if(jet->pt() < 15.0 || fabs(jet->eta()) > 5) continue;
// == old old old 
//cout<<endl;
//cout<<"n sv:  "<<nSV<<endl;
//cout<<nTracks<<endl;
//cout<<"mass vs corrected:  "<<massIVF<<" , "<<m_vertex_ptcorr<<endl;
//JetRecoSVMassAVRPFCorrected[NumRecoJetsPFCorrected] = jet->userFloat("vtxMass");
//JetRecoSVMassIVFPFCorrected[NumRecoJetsPFCorrected] = massIVF;
//JetRecoSVMassCorrectedPFCorrected[NumRecoJetsPFCorrected] = m_vertex_ptcorr;
//JetRecoSVPt[NumRecoJetsPFCorrected] = ptSV;
//JetRecoSVN[NumRecoJetsPFCorrected] = nTracks;
JetRecoPUMVAPFCorrected[NumRecoJetsPFCorrected] = mva;
JetRecoPUIDPFCorrected[NumRecoJetsPFCorrected] = idflag ;

float NHF = jet->neutralHadronEnergyFraction();
float NEMF = jet->neutralEmEnergyFraction();
float CHF = jet->chargedHadronEnergyFraction();
float MUF = jet->muonEnergyFraction();
float CEMF = jet->chargedEmEnergyFraction();
int NumConst = jet->chargedMultiplicity()+jet->neutralMultiplicity();
int NumNeutralParticles =jet->neutralMultiplicity();
int CHM = jet->chargedMultiplicity();

//cout<<jet->mass()<<endl;
VJP.SetPtEtaPhiE(jet->pt(), jet->eta(), jet->phi(), jet->et()*cosh(jet->eta()));

JetRecoPtPFCorrected[NumRecoJetsPFCorrected] = jet->pt();
        JetRecoEtPFCorrected[NumRecoJetsPFCorrected] = jet->et();
        JetRecoEtaPFCorrected[NumRecoJetsPFCorrected] = jet->eta();
        JetRecoPhiPFCorrected[NumRecoJetsPFCorrected] = jet->phi();
        JetRecoRapidityPFCorrected[NumRecoJetsPFCorrected] = jet->rapidity();

//========= jetid things ===========
JetRecoNHFPFCorrected[NumRecoJetsPFCorrected] = NHF;
JetRecoNEMFPFCorrected[NumRecoJetsPFCorrected] = NEMF;
JetRecoCHFPFCorrected[NumRecoJetsPFCorrected] = CHF;
JetRecoMUFPFCorrected[NumRecoJetsPFCorrected] = MUF;
JetRecoCEMFPFCorrected[NumRecoJetsPFCorrected] = CEMF;
JetRecoNumConstPFCorrected[NumRecoJetsPFCorrected] = NumConst;
JetRecoNumNeutralParticlesPFCorrected[NumRecoJetsPFCorrected] = NumNeutralParticles;
JetRecoCHMPFCorrected[NumRecoJetsPFCorrected] = CHM;

jecUnc->setJetEta(jet->eta());
    jecUnc->setJetPt(jet->pt());
     double unc_up = jecUnc->getUncertainty(true);
        JetRecoUNCPFCorrected[NumRecoJetsPFCorrected] = unc_up;

//cout<<"unc:   "<<unc_up<<endl;



DiscMVA[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfCombinedMVAV2BJetTags");
DiscCSV[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
DiscCvB[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfCombinedCvsBJetTags");
DiscCvL[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfCombinedCvsLJetTags");
DiscDeepCSVb[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepCSVJetTags:probb");
DiscDeepCSVbb[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepCSVJetTags:probbb");
DiscDeepCSVc[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepCSVJetTags:probc");
DiscDeepCSVudsg[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepCSVJetTags:probudsg");

DiscDeepJetb[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepFlavourJetTags:probb");
DiscDeepJetbb[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepFlavourJetTags:probbb");
DiscDeepJetlepb[NumRecoJetsPFCorrected] = jet->bDiscriminator("pfDeepFlavourJetTags:problepb");



double CvsL = DiscDeepCSVc[NumRecoJetsPFCorrected]/(DiscDeepCSVc[NumRecoJetsPFCorrected] + DiscDeepCSVudsg[NumRecoJetsPFCorrected]);
double CvsB = DiscDeepCSVc[NumRecoJetsPFCorrected]/(DiscDeepCSVb[NumRecoJetsPFCorrected] + DiscDeepCSVbb[NumRecoJetsPFCorrected]);

if(CvsL > 0.59 && CvsB > 0.05) tightCtag = true;
if(DiscDeepCSVb[NumRecoJetsPFCorrected] +   DiscDeepCSVbb[NumRecoJetsPFCorrected] > 0.6321) mediumBtag = true;

 JetRecoPFlavPFCorrected[NumRecoJetsPFCorrected] = jet->partonFlavour() ;
        JetRecoHFlavPFCorrected[NumRecoJetsPFCorrected] = jet->hadronFlavour() ;

//cout<<"jet pt:  "<<jet->pt()<<" , "<<jet->eta()<<" , "<<jet->phi()<<endl;
//cout<<"jet pt, eta, flav:  "<<jet->pt()<<" , "<<jet->eta()<<" , "<<(*theJetFlavourInfos)[jetref->refAt(ijet)].getHadronFlavour()<<endl;
//cout<<"jet pt, eta, flav:  "<<jet->hadronFlavour()<<endl;
NumRecoJetsPFCorrected++;

}

//==========================================================================================================
//===================================== vertices ===========================================================
//========================================= electrons ===============================
/*
edm::Handle<pat::ElectronCollection> electronColl;
iEvent.getByToken(tok_ele_, electronColl);
NumRecoElectrons=0;

for (pat::ElectronCollection::const_iterator el=electronColl->begin(), electronCollEnd=electronColl->end(); el!=electronCollEnd;
++el)
{
const int nEtaBins = 7;
    const float etaBinLimits[nEtaBins+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    const float effectiveAreaValues[nEtaBins] = {0.1440, 0.1562, 0.1032, 0.0859, 0.1116, 0.1321, 0.1654};

if(el->pt() < 5.0) continue;
const reco::GsfElectron::PflowIsolationVariables& pfIso = el->pfIsolationVariables();
float etaSC = el->superCluster()->eta();

Esc[NumRecoElectrons]  = el->superCluster()->energy();
Etasc[NumRecoElectrons]  = el->superCluster()->eta();

 int etaBin = 0;
    while(etaBin < nEtaBins-1 && abs(etaSC) > etaBinLimits[etaBin+1])  ++etaBin;

float area = effectiveAreaValues[etaBin];
    relIsoWithEA[NumRecoElectrons] = ( pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - (Rho) * area ) )/el->pt();

trackMomentumAtVtx[NumRecoElectrons]               = sqrt(el->trackMomentumAtVtx().mag2());
ecalEnergy[NumRecoElectrons]                       = el->ecalEnergy();

ElectronIdLoose[NumRecoElectrons] = el->electronID("cutBasedElectronID-Fall17-94X-V2-loose");
ElectronIdMedium[NumRecoElectrons] = el->electronID("cutBasedElectronID-Fall17-94X-V2-medium");
ElectronIdTight[NumRecoElectrons] = el->electronID("cutBasedElectronID-Fall17-94X-V2-tight");

ElectronPt[NumRecoElectrons] = el->pt();
ElectronEta[NumRecoElectrons] = el->eta();
ElectronPhi[NumRecoElectrons] = el->phi();
ElectronEnergy[NumRecoElectrons] = el->energy();
ElectronCharge[NumRecoElectrons] = el->charge();

dxy[NumRecoElectrons] = el->gsfTrack()->dxy(pv.position());
dz[NumRecoElectrons] = el->gsfTrack()->dz(pv.position());
NumRecoElectrons++;
}
*/
//==========================================================================================================
genppt = -1;
genpptP = -1;
genpptM = -1;
cout<<"fill"<<endl;
//if(numMuP > 0 && numMuM > 0  && NumGenPFinal == 2)
if(numMuP > 0 && numMuM > 0)
{
wwtree->Fill();
}
truecounter++;
wwtree2->Fill();
//==============================================================================================================================================


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Select::beginJob()
{
using namespace std;
f = new TFile(fileName.c_str(), "RECREATE");	  //// the name of the Root file should be set in Jet.cfg!
wwtree  = new TTree("wztree", "tracks tree");
wwtree2  = new TTree("wztree2", "tracks tree");
wwtree->Branch("genppt", &genppt, "genppt/F");
wwtree->Branch("genpptP", &genpptP, "genpptP/F");
wwtree->Branch("genpptM", &genpptM, "genpptM/F");
wwtree->Branch("run", &run, "run/I");
wwtree->Branch("lumisection", &lumisection, "lumisection/I");
wwtree->Branch("event", &event, "event/I");

wwtree->Branch("tightCtag", &tightCtag, "tightCtag/B");
wwtree->Branch("mediumBtag", &mediumBtag, "mediumBtag/B");

wwtree->Branch("trig", &trig, "trig/B");
wwtree->Branch("muplus", &muplus, "muplus/I");
wwtree->Branch("muminus", &muminus, "muminus/I");


wwtree->Branch("pfMET", &pfMET, "pfMET/D");
wwtree->Branch("pfMETphi", &pfMETphi, "pfMETphi/D");

wwtree->Branch("evtweight", &evtweight);
wwtree->Branch("evtwid", &evtwid);
wwtree->Branch("triggers", &triggers);

/*
wwtree->Branch("NumGenP", &NumGenP, "NumGenP/I");
wwtree->Branch("GenM", &GenM, "GenM[20]/F");
wwtree->Branch("GenPx", &GenPx, "GenPx[20]/F");
wwtree->Branch("GenPy", &GenPy, "GenPy[20]/F");
wwtree->Branch("GenPz", &GenPz, "GenPz[20]/F");
wwtree->Branch("GenE", &GenE, "GenE[20]/F");
wwtree->Branch("GenId", &GenId, "GenId[20]/I");
wwtree->Branch("GenStatus", &GenStatus, "GenStatus[20]/I");


wwtree->Branch("NumGeneratorP", &NumGeneratorP, "NumGeneratorP/I");
wwtree->Branch("generatorpt", &generatorpt);
wwtree->Branch("generatoreta", &generatoreta);
wwtree->Branch("generatorphi", &generatorphi);
wwtree->Branch("generatore", &generatore);
wwtree->Branch("generatorid", &generatorid);
wwtree->Branch("generatormotherid", &generatormotherid);
wwtree->Branch("generatorstatus", &generatorstatus);
*/

wwtree->Branch("NSV", &NSV, "NSV/I");
wwtree->Branch("sv_x", &sv_x);
wwtree->Branch("sv_y", &sv_y);
wwtree->Branch("sv_z", &sv_z);
wwtree->Branch("sv_pt_x", &sv_pt_x);
wwtree->Branch("sv_pt_y", &sv_pt_y);
wwtree->Branch("sv_pt_z", &sv_pt_z);
wwtree->Branch("sv_significance", &sv_significance);
wwtree->Branch("sv_value", &sv_value);
wwtree->Branch("sv_error", &sv_error);
wwtree->Branch("sv_mass", &sv_mass);
wwtree->Branch("sv_mass_corr", &sv_mass_corr);
wwtree->Branch("sv_pt", &sv_pt);
wwtree->Branch("sv_ntrks", &sv_ntrks);

wwtree->Branch("event32u", &event32u, "event32u/i"); // 32 unsigned
wwtree->Branch("event64", &event64, "event64/L"); //64
wwtree->Branch("event64u", &event64u, "event64u/l"); //65 unsigned

wwtree->Branch("NvtxEv", &NvtxEv, "NvtxEv/I");
wwtree->Branch("lumi", &lumi, "lumi/F");
wwtree2->Branch("lumi", &lumi, "lumi/F");
wwtree->Branch("Rho", &Rho, "Rho/D");

wwtree2->Branch("evtwid", &evtwid);
wwtree2->Branch("evtweight", &evtweight);
wwtree2->Branch("truecounter", &truecounter, "truecounter/I");
wwtree2->Branch("weight", &weight, "weight/F");
wwtree->Branch("weight", &weight, "weight/F");

wwtree->Branch("weightPS", &weightPS, "weightPS/F");
wwtree->Branch("weightPSUp", &weightPSUp, "weightPSUp/F");
wwtree->Branch("weightPSDown", &weightPSDown, "weightPSDown/F");

wwtree->Branch("mode1", &mode1, "mode1/I");
wwtree->Branch("mode2", &mode2, "mode2/I");


/*
wwtree->Branch("NumTrigObj", &NumTrigObj, "NumTrigObj/I");
wwtree->Branch("TrigEta", &TrigEta, "TrigEta[20]/F");
wwtree->Branch("TrigPhi", &TrigPhi, "TrigPhi[20]/F");
*/


wwtree->Branch("DiscMVA", &DiscMVA, "DiscMVA[30]/F");
wwtree->Branch("DiscCSV", &DiscCSV, "DiscCSV[30]/F");
wwtree->Branch("DiscCvB", &DiscCvB, "DiscCvB[30]/F");
wwtree->Branch("DiscCvL", &DiscCvL, "DiscCvL[30]/F");
wwtree->Branch("DiscDeepCSVb", &DiscDeepCSVb, "DiscDeepCSVb[30]/F");
wwtree->Branch("DiscDeepCSVbb", &DiscDeepCSVbb, "DiscDeepCSVbb[30]/F");
wwtree->Branch("DiscDeepCSVc", &DiscDeepCSVc, "DiscDeepCSVc[30]/F");
wwtree->Branch("DiscDeepCSVudsg", &DiscDeepCSVudsg, "DiscDeepCSVudsg[30]/F");

wwtree->Branch("DiscDeepJetb", &DiscDeepJetb, "DiscDeepJetb[30]/F");
wwtree->Branch("DiscDeepJetbb", &DiscDeepJetbb, "DiscDeepJetbb[30]/F");
wwtree->Branch("DiscDeepJetlepb", &DiscDeepJetlepb, "DiscDeepJetlepb[30]/F");

wwtree->Branch("NumRecoJetsAK8PFCorrected", &NumRecoJetsAK8PFCorrected, "NumRecoJetsAK8PFCorrected/I");
wwtree->Branch("JetRecoPtAK8PFCorrected", &JetRecoPtAK8PFCorrected, "JetRecoPtAK8PFCorrected[15]/F");
wwtree->Branch("JetRecoEtAK8PFCorrected", &JetRecoEtAK8PFCorrected, "JetRecoEtAK8PFCorrected[15]/F");
wwtree->Branch("JetRecoEtaAK8PFCorrected", &JetRecoEtaAK8PFCorrected, "JetRecoEtaAK8PFCorrected[15]/F");
wwtree->Branch("JetRecoPhiAK8PFCorrected", &JetRecoPhiAK8PFCorrected, "JetRecoPhiAK8PFCorrected[15]/F");
wwtree->Branch("JetRecoMassAK8PFCorrected", &JetRecoMassAK8PFCorrected, "JetRecoMassAK8PFCorrected[15]/F");



wwtree->Branch("pfBoostedDoubleSecondaryVertexAK8BJetTags", &pfBoostedDoubleSecondaryVertexAK8BJetTags, "pfBoostedDoubleSecondaryVertexAK8BJetTags[15]/F");
//wwtree->Branch("pfMassIndependentDeepDoubleBvLJetTags", &pfMassIndependentDeepDoubleBvLJetTags, "pfMassIndependentDeepDoubleBvLJetTags[15]/F");

wwtree->Branch("NjettinessAK8Puppitau1", &NjettinessAK8Puppitau1, "NjettinessAK8Puppitau1[15]/F");
wwtree->Branch("NjettinessAK8Puppitau2", &NjettinessAK8Puppitau2, "NjettinessAK8Puppitau2[15]/F");
wwtree->Branch("NjettinessAK8Puppitau3", &NjettinessAK8Puppitau3, "NjettinessAK8Puppitau3[15]/F");
wwtree->Branch("NjettinessAK8Puppitau4", &NjettinessAK8Puppitau4, "NjettinessAK8Puppitau4[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapNjettinessAK8CHSTau1", &ak8PFJetsCHSValueMapNjettinessAK8CHSTau1, "ak8PFJetsCHSValueMapNjettinessAK8CHSTau1[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapNjettinessAK8CHSTau2", &ak8PFJetsCHSValueMapNjettinessAK8CHSTau2, "ak8PFJetsCHSValueMapNjettinessAK8CHSTau2[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapNjettinessAK8CHSTau3", &ak8PFJetsCHSValueMapNjettinessAK8CHSTau3, "ak8PFJetsCHSValueMapNjettinessAK8CHSTau3[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapNjettinessAK8CHSTau4", &ak8PFJetsCHSValueMapNjettinessAK8CHSTau4, "ak8PFJetsCHSValueMapNjettinessAK8CHSTau4[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass", &ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass, "ak8PFJetsCHSValueMapak8PFJetsCHSPrunedMass[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass", &ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass, "ak8PFJetsCHSValueMapak8PFJetsCHSSoftDropMass[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapeta", &ak8PFJetsCHSValueMapeta, "ak8PFJetsCHSValueMapeta[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapmass", &ak8PFJetsCHSValueMapmass, "ak8PFJetsCHSValueMapmass[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMapphi", &ak8PFJetsCHSValueMapphi, "ak8PFJetsCHSValueMapphi[15]/F");
wwtree->Branch("ak8PFJetsCHSValueMappt", &ak8PFJetsCHSValueMappt, "ak8PFJetsCHSValueMappt[15]/F");
wwtree->Branch("ak8PFJetsPuppiSoftDropMass", &ak8PFJetsPuppiSoftDropMass, "ak8PFJetsPuppiSoftDropMass[15]/F");
wwtree->Branch("ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2", &ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2, "ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN2[15]/F");
wwtree->Branch("ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3", &ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3, "ak8PFJetsPuppiSoftDropValueMapnb1AK8PuppiSoftDropN3[15]/F");
wwtree->Branch("ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2", &ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2, "ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN2[15]/F");
wwtree->Branch("ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3", &ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3, "ak8PFJetsPuppiSoftDropValueMapnb2AK8PuppiSoftDropN3[15]/F");






wwtree->Branch("NumRecoJetsPFCorrected", &NumRecoJetsPFCorrected, "NumRecoJetsPFCorrected/I");
wwtree->Branch("JetRecoPtPFCorrected", &JetRecoPtPFCorrected, "JetRecoPtPFCorrected[30]/F");
wwtree->Branch("JetRecoEtPFCorrected", &JetRecoEtPFCorrected, "JetRecoEtPFCorrected[30]/F");
wwtree->Branch("JetRecoEtaPFCorrected", &JetRecoEtaPFCorrected, "JetRecoEtaPFCorrected[30]/F");
wwtree->Branch("JetRecoRapidityPFCorrected", &JetRecoRapidityPFCorrected, "JetRecoRapidityPFCorrected[30]/F");
wwtree->Branch("JetRecoPhiPFCorrected", &JetRecoPhiPFCorrected, "JetRecoPhiPFCorrected[30]/F");
wwtree->Branch("JetRecoUNCPFCorrected", &JetRecoUNCPFCorrected, "JetRecoUNCPFCorrected[30]/F");

wwtree->Branch("JetRecoNHFPFCorrected", &JetRecoNHFPFCorrected, "JetRecoNHFPFCorrected[30]/F");
wwtree->Branch("JetRecoNEMFPFCorrected", &JetRecoNEMFPFCorrected, "JetRecoNEMFPFCorrected[30]/F");
wwtree->Branch("JetRecoCHFPFCorrected", &JetRecoCHFPFCorrected, "JetRecoCHFPFCorrected[30]/F");
wwtree->Branch("JetRecoMUFPFCorrected", &JetRecoMUFPFCorrected, "JetRecoMUFPFCorrected[30]/F");
wwtree->Branch("JetRecoCEMFPFCorrected", &JetRecoCEMFPFCorrected, "JetRecoCEMFPFCorrected[30]/F");
wwtree->Branch("JetRecoNumConstPFCorrected", &JetRecoNumConstPFCorrected, "JetRecoNumConstPFCorrected[30]/I");
wwtree->Branch("JetRecoNumNeutralParticlesPFCorrected", &JetRecoNumNeutralParticlesPFCorrected, "JetRecoNumNeutralParticlesPFCorrected[30]/I");
wwtree->Branch("JetRecoCHMPFCorrected", &JetRecoCHMPFCorrected, "JetRecoCHMPFCorrected[30]/F");


wwtree->Branch("JetRecoPUIDPFCorrected", &JetRecoPUIDPFCorrected, "JetRecoPUIDPFCorrected[30]/F");
wwtree->Branch("JetRecoPUMVAPFCorrected", &JetRecoPUMVAPFCorrected, "JetRecoPUMVAPFCorrected[30]/F");
wwtree->Branch("JetRecoSVN", &JetRecoSVN, "JetRecoSVN[30]/I");
  wwtree->Branch("JetRecoPFlavPFCorrected", &JetRecoPFlavPFCorrected, "JetRecoPFlavPFCorrected[30]/I");
    wwtree->Branch("JetRecoHFlavPFCorrected", &JetRecoHFlavPFCorrected, "JetRecoHFlavPFCorrected[30]/I");
wwtree->Branch("JetRecoPUIDNEWPFCorrected", &JetRecoPUIDNEWPFCorrected, "JetRecoPUIDNEWPFCorrected[30]/F");



wwtree->Branch("NumRecoMuons", &NumRecoMuons, "NumRecoMuons/I");
wwtree->Branch("RecoPtMuons", &RecoPtMuons, "RecoPtMuons[25]/F");
wwtree->Branch("RecoEtaMuons", &RecoEtaMuons, "RecoEtaMuons[25]/F");
wwtree->Branch("RecoPhiMuons", &RecoPhiMuons, "RecoPhiMuons[25]/F");
wwtree->Branch("RecoChargeMuons", &RecoChargeMuons, "RecoChargeMuons[25]/I");
wwtree->Branch("RecoIsoMuonsDxyvtx", &RecoIsoMuonsDxyvtx, "RecoIsoMuonsDxyvtx[25]/F");
wwtree->Branch("RecoIsoMuonsDzvtx", &RecoIsoMuonsDzvtx, "RecoIsoMuonsDzvtx[25]/F");
wwtree->Branch("isLoose", &isLoose, "isLoose[25]/B");
wwtree->Branch("isMedium", &isMedium, "isMedium[25]/B");
wwtree->Branch("isTight", &isTight, "isTight[25]/B");
wwtree->Branch("NumLayers", &NumLayers, "NumLayers[25]/I");

wwtree->Branch("trkIsolation", &trkIsolation, "trkIsolation[25]/F");
wwtree->Branch("pfIsolation", &pfIsolation, "pfIsolation[25]/F");

/*
wwtree->Branch("NumRecoElectrons", &NumRecoElectrons, "NumRecoElectrons/I");
wwtree->Branch("mva", &mva, "mva[25]/F");
wwtree->Branch("ElectronPt", &ElectronPt, "ElectronPt[25]/F");
wwtree->Branch("ElectronEta", &ElectronEta, "ElectronEta[25]/F");
wwtree->Branch("ElectronPhi", &ElectronPhi, "ElectronPhi[25]/F");
wwtree->Branch("ElectronEnergy", &ElectronEnergy, "ElectronEnergy[25]/F");
wwtree->Branch("ElectronCharge", &ElectronCharge, "ElectronCharge[25]/I");
wwtree->Branch("trackMomentumAtVtx", &trackMomentumAtVtx, "trackMomentumAtVtx[25]/F");
wwtree->Branch("ecalEnergy", &ecalEnergy, "ecalEnergy[25]/F");
wwtree->Branch("dxy", &dxy, "dxy[25]/F");
wwtree->Branch("dz", &dz, "dz[25]/F");
wwtree->Branch("Esc", &Esc, "Esc[25]/F");
wwtree->Branch("Etasc", &Etasc, "Eetasc[25]/F");
wwtree->Branch("ElectronIdLoose", &ElectronIdLoose, "ElectronIdLoose[25]/F");
wwtree->Branch("ElectronIdMedium", &ElectronIdMedium, "ElectronIdMedium[25]/F");
wwtree->Branch("ElectronIdTight", &ElectronIdTight, "ElectronIdTight[25]/F");
*/


cout << "here ------------- here ------------ here ----------- here " <<endl;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Select::endJob() 
{
std::cout<<"Select::endjob"<<std::endl;

//MET_phi->Write();
//MET_phi_corr->Write();

using namespace std;
   cout << "------------------------------->>>>> End Job" << endl;
f->WriteTObject(wwtree);
f->WriteTObject(wwtree2);
delete wwtree;
delete wwtree2;
f->Close();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Select::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Select);

