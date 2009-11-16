import FWCore.ParameterSet.Config as cms

analyzeMuonQuality = cms.EDAnalyzer("MuonQualityAnalyzer",
    ## input collection                             
    src = cms.InputTag("selectedLayer1Muons"),
    ## event weight
    weight = cms.InputTag("eventWeight"),
    ## use the weight or not                             
    useWeight = cms.bool(False),
    ## special parameters for muon quality analysis
    analyze   = cms.PSet(
      ## fill histograms for 1.,2.,3.,... leading
      ## muon? -1 corresponds to 'all'
      ## counting starts with 0=leading Muon! 
      index = cms.int32(-1)
    )
)



