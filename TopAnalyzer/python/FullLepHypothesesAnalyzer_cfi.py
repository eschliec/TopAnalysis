import FWCore.ParameterSet.Config as cms

analyzeFullLepHypotheses = cms.EDAnalyzer("FullLepHypothesesAnalyzer",
    hist = cms.string('analyzeFullLepHypotheses.hist'),
    hypoKey = cms.InputTag("ttFullLepHypKinSolution","Key"),
    weight = cms.InputTag("eventWeight"),
    FullLepEvent = cms.InputTag("ttFullLepEvent"),
    alsoWrongCharge = cms.bool(True)
)



