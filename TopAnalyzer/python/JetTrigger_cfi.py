import FWCore.ParameterSet.Config as cms

analyzeJetTrigger = cms.EDAnalyzer("JetTrigger",
                                   TriggerResults     = cms.InputTag('TriggerResults','','HLT'),
                                   TriggerSummary     = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
                                   patTriggerEvent    = cms.InputTag('patTriggerEvent'),
                                   patTrigger         = cms.InputTag('patTrigger'),
                                   jets               = cms.InputTag('selectedLayer1Jets'),
                                   triggerMatchedJets = cms.InputTag('triggerMatchedSelectedLayer1Jets'),
                                   analyzedTrigger    = cms.InputTag('hlt4jet30', '', 'HLT')
)



