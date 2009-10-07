import FWCore.ParameterSet.Config as cms

filterDiMuonMass = cms.EDAnalyzer("DiMuonMassFilter",

    ## sources
    muons = cms.InputTag("selectedLayer1Muons"),
    
    ## cuts on muon muon mass
    qcdCut = cms.double(20.),
    zCut   = cms.vdouble(80.,105.)
)


