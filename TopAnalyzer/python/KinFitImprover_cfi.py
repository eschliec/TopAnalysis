import FWCore.ParameterSet.Config as cms

analyzeKinFitImprover = cms.EDAnalyzer("KinFitImproveAnalyzer",
    ## input collectionA
    srcA = cms.InputTag("ttFullHadEvent"),
    ## input collectionB                            
    srcB = cms.InputTag("selectedLayer1Jets"),                                         
    ## event weight
    weight = cms.InputTag("eventWeight"),
    ## use the weight or not                             
    useWeight = cms.bool(False),
    ## analyzer specific configurables
    analyze   = cms.PSet(
      ## no additional configurables needed in this case
    )    
)



