import FWCore.ParameterSet.Config as cms

analyzeKinFitQuality = cms.EDAnalyzer("KinFitQualityAnalyzer",
    ## input collectionA
    srcA = cms.InputTag("ttFullHadEvent"),
    ## input collectionB                            
    srcB = cms.InputTag("selectedPatJets"),                                         
    ## event weight
    weight = cms.InputTag("eventWeight"),
    ## use the weight or not                             
    useWeight = cms.bool(False),
    ## analyzer specific configurables
    analyze   = cms.PSet(
      ## no additional configurables needed in this case
    )    
)



