import FWCore.ParameterSet.Config as cms

analyzeSemiLepHypotheses = cms.EDAnalyzer("SemiLepHypothesesAnalyzer",
    hist = cms.string('analyzeSemiLepHypotheses.hist'),
    hypoKey = cms.InputTag("ttSemiLepHypWMassMaxSumPt","Key"),
    weight = cms.InputTag("eventWeight"),
    semiLepEvent = cms.InputTag("ttSemiLepEvent"),
    maxSumDRGenMatch = cms.double(0.8),
    nJetsMax = cms.uint32(8),
    minMVADisc = cms.double(0.6),
    minProbKinFit = cms.double(0.05)
)



