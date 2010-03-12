import FWCore.ParameterSet.Config as cms

analyzeSingleObjectJetTrigger = cms.EDAnalyzer("SingleObjectJetTrigger",
           patTriggerEvent    = cms.InputTag('patTriggerEvent'),
           patTrigger         = cms.InputTag('patTrigger'),
           jets               = cms.InputTag('selectedLayer1Jets'),
           triggerMatchedJets = cms.InputTag('triggerMatchedSelectedLayer1Jets')
)



