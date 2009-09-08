import FWCore.ParameterSet.Config as cms

from TopAnalysis.TopAnalyzer.FullLepHypothesesAnalyzer_cfi import analyzeFullLepHypotheses
## clone module
analyzeKinSolution = analyzeFullLepHypotheses.clone()
analyzeGenMatch    = analyzeFullLepHypotheses.clone()

## do replacements
analyzeKinSolution.hypoKey = 'ttFullLepHypKinSolution:Key'
analyzeGenMatch.hypoKey    = 'ttFullLepHypGenMatch:Key'

analyzeAllFullLepHypotheses = cms.Sequence(analyzeGenMatch +
                                           analyzeKinSolution     
                                          )



