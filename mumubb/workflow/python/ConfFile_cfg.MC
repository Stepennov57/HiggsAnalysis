import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

from Configuration.AlCa.GlobalTag import GlobalTag

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')


process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')



process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")  #Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softLepton_cff import *
from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoBTag.Combined.combinedMVA_cff import *
from RecoBTag.CTagging.RecoCTagging_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3')


from PhysicsTools.PatAlgos.patTemplate_cfg import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.options.numberOfThreads = cms.untracked.uint32(4)

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#	'/store/mc/RunIISummer16MiniAODv3/SUSYGluGluToHToAA_AToMuMu_AToBB_M-15_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/20BFE2AE-813B-E911-8EC9-0242AC130002.root',
#	'/store/mc/RunIISummer16MiniAODv3/SUSYGluGluToHToAA_AToMuMu_AToBB_M-15_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/601BF2D4-8C3B-E911-911D-0242AC130002.root',
#	'/store/mc/RunIISummer16MiniAODv3/SUSYGluGluToHToAA_AToMuMu_AToBB_M-15_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/EA6F5CFF-733B-E911-B579-0242AC130002.root'
#	'/store/mc/RunIISummer16MiniAODv3/SUSYGluGluToHToAA_AToBB_AToTauTau_M-45_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/20000/10DC3143-6700-EA11-8089-0CC47AFCC376.root'
#	'/store/mc/RunIIFall17MiniAODv2/GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/08F8F8C0-41B5-E811-9A00-0CC47AFC3D32.root'
#	'/store/mc/RunIIFall17MiniAODv2/GluGluToHToTauTauMaxmixDecay_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/0012BD3B-3EB3-E811-B9E6-FA163EF96BC3.root'
#	'/store/mc/RunIISummer16MiniAODv2/SUSYGluGluToHToAA_AToMuMu_AToBB_M-20_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/70000/8480F340-26DB-E611-A91B-0CC47A13CDB0.root'
	'/store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver1-v1/80000/3A6CBD55-208C-E811-BE87-0242AC1C0502.root'
#	'/store/mc/RunIISummer16MiniAODv3/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_backup_94X_mcRun2_asymptotic_v3-v2/90000/EE035E2C-3338-E911-8CE3-AC1F6BAC7D1A.root '
#	'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/120000/007A2300-43DF-E811-A43F-842B2B688F18.root'	
#       	'/store/mc/RunIISummer16MiniAODv3/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/270000/94686A2A-DFE0-E811-9B88-00266CFCC861.root'
#	'/store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/270000/DEE4F680-14E9-E811-9CE7-48FD8EE73A8D.root'
#	'/store/mc/RunIISummer16MiniAODv3/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/270000/FCA5F186-66E4-E811-9912-509A4C748A3D.root'
#	'/store/mc/RunIISummer16MiniAODv3/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1/60000/3EEF0067-4ABD-E811-A178-0CC47AB35D36.root'
#	'/store/mc/RunIISummer16MiniAODv3/DY01234jets_13TeV-sherpa/MINIAODSIM/PUMoriond17_pilot_94X_mcRun2_asymptotic_v3-v2/260000/C650C59D-1827-EA11-98D0-003048F596AE.root'
#	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-12_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/130000/101D20D8-C0CE-E911-A7DA-7CD30AB053DC.root',
# 	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-12_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/130000/1A82D2F8-ABCE-E911-9726-24BE05CEADA1.root',#
#	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-12_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/130000/26252F9B-92CB-E911-A424-FA163EC5EAEC.root'
#	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-20_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/260000/0644D2AB-F0D3-E911-B414-EC0D9A82264E.root',
#	'/store/mc/RunIISummer16MiniAODv3/SUSYGluGluToHToAA_AToMuMu_AToBB_M-20_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/110000/DC0DDEC3-8440-E911-A7DD-0025905C3D40.root'
#	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-20_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/260000/1231F28D-80D0-E911-87EC-0025901C0610.root',
#	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-20_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/260000/022C7391-7CD3-E911-80B3-0CC47A4D7640.root',
#	'/store/mc/RunIIFall17MiniAODv2/SUSYWHToAA_AATo4B_M-20_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/260000/24FB15B6-F0D3-E911-99FF-B083FED42ECF.root'
 )

)



from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                      connect = cms.string('sqlite_file:Summer16_07Aug2017_V11_MC.db'),
                          toGet =  cms.VPSet(
        cms.PSet(record = cms.string("JetCorrectionsRecord"),
                         tag = cms.string("JetCorrectorParametersCollection_Summer16_07Aug2017_V11_MC_AK4PFchs"),
                         label=cms.untracked.string("AK4PFchs")
                           )
                                         )
                           )
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

process.demo = cms.EDAnalyzer('Select',
                                        isMC = cms.bool(False),
                                        updatedPatJetsUpdatedJEC = cms.InputTag("updatedPatJets"),
					flavourMap  = cms.InputTag("genJetFlavourInfos")
#					jetFlavourInfos = cms.InputTag("jetFlavourInfosAK4PFJets")
)

#process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
 #pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
  # svSource = cms.InputTag('slimmedSecondaryVertices'),
#   btagInfos = bTagInfos,
   btagDiscriminators = ['pfCombinedInclusiveSecondaryVertexV2BJetTags',
                              'pfJetProbabilityBJetTags',
                              'pfCombinedMVAV2BJetTags',

                              'pfCombinedCvsLJetTags',
                              'pfCombinedCvsBJetTags',

                              'pfDeepCSVJetTags:probudsg',
                              'pfDeepCSVJetTags:probbb',
                              'pfDeepCSVJetTags:probb',
                              'pfDeepCSVJetTags:probc'] 
)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

#process.updatedPatJets.addBTagInfo     = cms.bool(True)
process.updatedPatJets.addTagInfos = cms.bool(True)
#process.updatedPatJetsTransientCorrected.addBTagInfo     = cms.bool(True)
process.updatedPatJetsTransientCorrected.addTagInfos = cms.bool(True)
process.pfInclusiveSecondaryVertexFinderCvsLTagInfos.extSVCollection = cms.InputTag("slimmedSecondaryVertices")


from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
#     jets=cms.InputTag("slimmedGenJets"),
#    jets = cms.InputTag("ak4GenJets"),
     particles=cms.InputTag("prunedGenParticles")
)

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
  #  jets = cms.InputTag("ak4GenJets")
   jets=cms.InputTag("slimmedGenJets")
    )

from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
#    genParticles = cms.InputTag("ak4GenJets"),
    genParticles=cms.InputTag("prunedGenParticles"),
    jetFlavourInfos = cms.InputTag("genJetFlavourInfos"),
#    jetFlavourInfos = "genJetFlavourInfos"
    )

from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
process.matchGenCHadron = matchGenCHadron.clone(
#    genParticles = cms.InputTag("ak4GenJets"),
    genParticles=cms.InputTag("prunedGenParticles"),
    jetFlavourInfos = cms.InputTag("genJetFlavourInfos"),
    #jetFlavourInfos = "genJetFlavourInfos"
    )

#from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
#process.jetFlavourInfosAK4PFJets = ak4JetFlavourInfos.clone()
#process.jetFlavourInfosAK4PFJets.jets = cms.InputTag("slimmedJets")

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mcRun2_asymptotic_v3', '')

#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                       era='2016-Legacy')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.genParticles2HepMC.genParticles = cms.InputTag("mergedGenParticles")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
process.particleLevel.HepMCCollection = cms.InputTag("genParticles2HepMC:unsmeared")
process.particleLevel.lepMaxEta = cms.double(10.0)
process.particleLevel.lepMinPt  = cms.double(0.0)
process.particleLevel.particleMaxEta = cms.double(10.0)

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)

process.p = cms.Path(
process.egammaPostRecoSeq *
process.patJetCorrFactors *
process.updatedPatJets *
process.mergedGenParticles*process.genParticles2HepMC*process.particleLevel *
process.selectedHadronsAndPartons*
process.genJetFlavourInfos *
#process.dump *
#process.jetFlavourInfosAK4PFJets*
process.prefiringweight *
process.demo)


