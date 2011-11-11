import FWCore.ParameterSet.Config as cms

analyzeHypothesisKinFitLepton = cms.EDAnalyzer("HypothesisKinFitLeptonAnalyzer",
    ## input collectionA
    srcA = cms.InputTag("ttSemiLepEvent"),
    ## input collectionB                            
    srcB = cms.InputTag("selectedPatMuons"),                                         
    ## analyzer specific configurables
    analyze   = cms.PSet(
      ## no additional configurables needed in this case
    )    
)



