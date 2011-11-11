import FWCore.ParameterSet.Config as cms

analyzeHypothesisKinFitJets = cms.EDAnalyzer("HypothesisKinFitJetsAnalyzer",
    ## input collectionA
    srcA = cms.InputTag("ttSemiLepEvent"),
    ## input collectionB                            
    srcB = cms.InputTag("selectedPatJets"),                                         
    ## analyzer specific configurables
    analyze   = cms.PSet(
      ## no additional configurables needed in this case
    )    
)



