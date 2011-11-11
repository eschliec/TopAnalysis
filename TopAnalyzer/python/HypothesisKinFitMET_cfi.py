import FWCore.ParameterSet.Config as cms

analyzeHypothesisKinFitMET = cms.EDAnalyzer("HypothesisKinFitMETAnalyzer",
    ## input collectionA
    srcA = cms.InputTag("ttSemiLepEvent"),
    ## input collectionB                            
    srcB = cms.InputTag("patMETs"),                                         
    ## analyzer specific configurables
    analyze   = cms.PSet(
      ## no additional configurables needed in this case
    )    
)



