import FWCore.ParameterSet.Config as cms

filterFullLepObjects = cms.EDAnalyzer("FullLepObjectCountFilter",
    weight  = cms.InputTag("eventWeight"),
    name    = cms.string(''),
    objects = cms.InputTag(''),
    n       = cms.uint32(0)    
)
