import FWCore.ParameterSet.Config as cms

countObjects = cms.EDAnalyzer("ObjectCountFilter",
     weight  = cms.InputTag("eventWeight"),
     name    = cms.string(''),
     objects = cms.InputTag(''),
     n       = cms.uint32(0)    
)
