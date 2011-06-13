import FWCore.ParameterSet.Config as cms

analyzePartonMatch = cms.EDAnalyzer("PartonMatchAnalyzer",
    matches = cms.VInputTag(cms.InputTag("ttSemiEvtJetPartonMatch")),
    sumDR = cms.InputTag("ttSemiEvtJetPartonMatch","SumDR"),
    quarks = cms.vint32(),
    hist = cms.string('analyzePartonMatch.hist'),
    sumDPt = cms.InputTag("ttSemiEvtJetPartonMatch","SumPt"),
    jets = cms.vint32(),
    ref = cms.uint32(0)
)



