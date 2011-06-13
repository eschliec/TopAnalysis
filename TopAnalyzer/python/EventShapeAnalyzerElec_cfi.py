import FWCore.ParameterSet.Config as cms

analyzeEventShapeElec = cms.EDFilter("EventShapeAnalyzerElec",
    hist = cms.string('analyzeEventShapeElec.hist'),
    maxDeltaRJet = cms.double(0.3),
    minPt = cms.double(40.0),
    nJets = cms.uint32(4),
    maxEta = cms.double(1.4),
    lepton = cms.InputTag("selectedLayer1Electrons"),
    doMatching = cms.bool(True),
    maxDeltaRLepton = cms.double(0.3),
    jets = cms.InputTag("selectedLayer1Jets")
)


