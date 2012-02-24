import FWCore.ParameterSet.Config as cms
#from TopAnalysis.TopUtils.EventWeightPU_cfi import *

writeNTuple = cms.EDAnalyzer('NTupleWriter',
                       PUSource = cms.InputTag("addPileupInfo"),
                       hypoKey = cms.InputTag("ttFullLepHypKinSolution","Key"),
                       genParticles = cms.InputTag("genParticles"),
                       FullLepEvent = cms.InputTag("ttFullLepEvent"),
                       src = cms.InputTag("genEvt"),
                       muons = cms.InputTag("muons"),
                       elecs = cms.InputTag("elecs"),
                       jets      = cms.InputTag("jets"),
                       met       = cms.InputTag("met"),
                       vertices  = cms.InputTag('goodOfflinePrimaryVertices'),
#                       vertices  = cms.InputTag('offlinePrimaryVertices'),
                       genJets       = cms.InputTag("ak5GenJets", "", "HLT"),
                       BJetIndex     = cms.InputTag("produceHadronLevelBJets", "BJetIndex"),
                       AntiBJetIndex = cms.InputTag("produceHadronLevelBJets", "AntiBJetIndex"),
                       BHadJetIndex = cms.InputTag("produceGenLevelBJets", "BHadJetIndex"),
                       AntiBHadJetIndex = cms.InputTag("produceGenLevelBJets", "AntiBHadJetIndex"),
                       weightPU    =  cms.InputTag('eventWeightPU', 'eventWeightPU'),
                       weightPU3D  =  cms.InputTag('eventWeightPU', 'eventWeightPU3D'),
                       weightPU_Up    =  cms.InputTag('eventWeightPU', 'eventWeightPUUp'),
                       weightPU_Down    =  cms.InputTag('eventWeightPU', 'eventWeightPUDown'),
                       weightLepSF =  cms.InputTag('eventWeightDileptonSF', 'eventWeightDileptonSF'),
                       weightKinFit = cms.InputTag("eventWeightKinEffForNtuple", "eventWeightDileptonKinEffSF"),
                       datatype  =  cms.InputTag("data"),
                       includeTrigger = cms.bool(False),
                       triggerResults = cms.InputTag('TriggerResults','','HLT'),
                       decayMode = cms.InputTag("generatorTopFilter", "decayMode"),
                       isTtBarSample = cms.bool(True),
                       #directory = cms.string(""), #../NTupDir
                       directory = cms.string("../NTupDir"), #../NTupDir
)

