import FWCore.ParameterSet.Config as cms
from DileptonEventWeight_cfi import *

analyzeFullLepGenEvent = cms.EDAnalyzer("FullLepGenAnalyzer",
    src = cms.InputTag("genEvt"),      
    weightPU     = eventWeightPUTag,
    weightLepSF  = cms.InputTag(""),
    genJets             = cms.InputTag("ak5GenJets", "", "SIM"), #REPLACED HLT FOR SIM
    BJetIndex           = cms.InputTag("produceHadronLevelBJets", "BJetIndex"),
    AntiBJetIndex       = cms.InputTag("produceHadronLevelBJets", "AntiBJetIndex"),
)
