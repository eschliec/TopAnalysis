import FWCore.ParameterSet.Config as cms

countLeptons = cms.EDFilter("LeptonCounter",
    minNelec = cms.double(0.0),
    maxRatio = cms.double(10.0),
    numberOfMuonBins = cms.int32(30),
    minRatio = cms.double(0.0),
    muons = cms.InputTag("allLayer1Muons"),
    numberOfElecBins = cms.int32(30),
    maxNelec = cms.double(30.0),
    maxNmuon = cms.double(30.0),
    numberOfRatioBins = cms.int32(200),
    minNmuon = cms.double(0.0)
)



