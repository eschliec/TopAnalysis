import FWCore.ParameterSet.Config as cms

analyzeEventShapeMuon = cms.EDFilter("EventShapeAnalyzerMuon",
    hist = cms.string('analyzeEventShapeMuon.hist'),
    maxDeltaRJet = cms.double(0.3),
    minPt = cms.double(40.0),
    nJets = cms.uint32(4),
    maxEta = cms.double(1.4),
    lepton = cms.InputTag("selectedLayer1Muons"),
    doMatching = cms.bool(True),
    maxDeltaRLepton = cms.double(0.3),
    jets = cms.InputTag("selectedLayer1Jets")
)



