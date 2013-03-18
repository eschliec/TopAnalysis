import FWCore.ParameterSet.Config as cms

generatorTtDileptonFilter = cms.EDFilter("GeneratorTtDileptonFilter",
    src                 = cms.InputTag("genEvt"),
    genJets             = cms.InputTag("ak5GenJets", '', 'SIM'), #REPLACED HLT FOR SIM
    BJetIndex           = cms.InputTag("produceHadronLevelBJets", "BJetIndex"),
    AntiBJetIndex       = cms.InputTag("produceHadronLevelBJets", "AntiBJetIndex"),
    leptonPt	        = cms.double(0. ),
    centralLeptonEta    = cms.double(100),
    leptonEta           = cms.double(100),
    bPt		        = cms.double(0. ),
    bEta	        = cms.double(100),
    bHadronPt           = cms.double(0),
    bHadronEta          = cms.double(100),
    met	                = cms.double(0. ),
    invariantLepLepMass = cms.double(0. ),
    zVeto               = cms.vdouble(0.,0.)
)
