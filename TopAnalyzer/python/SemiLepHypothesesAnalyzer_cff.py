import FWCore.ParameterSet.Config as cms

from TopAnalysis.TopAnalyzer.SemiLepHypothesesAnalyzer_cfi import analyzeSemiLepHypotheses
## clone module
analyzeGeom          = analyzeSemiLepHypotheses.clone()
analyzeWMassMaxSumPt = analyzeSemiLepHypotheses.clone()
analyzeMaxSumPtWMass = analyzeSemiLepHypotheses.clone()
analyzeKinFit        = analyzeSemiLepHypotheses.clone()
analyzeGenMatch      = analyzeSemiLepHypotheses.clone()
analyzeMVADisc       = analyzeSemiLepHypotheses.clone()

## do replacements
analyzeGeom.hypoKey          = 'ttSemiLepHypGeom:Key'
analyzeWMassMaxSumPt.hypoKey = 'ttSemiLepHypWMassMaxSumPt:Key'
analyzeMaxSumPtWMass.hypoKey = 'ttSemiLepHypMaxSumPtWMass:Key'
analyzeKinFit.hypoKey        = 'ttSemiLepHypKinFit:Key'
analyzeGenMatch.hypoKey      = 'ttSemiLepHypGenMatch:Key'
analyzeMVADisc.hypoKey       = 'ttSemiLepHypMVADisc:Key'

analyzeAllSemiLepHypotheses = cms.Sequence(analyzeGeom          +
                                           analyzeWMassMaxSumPt +
                                           analyzeMaxSumPtWMass +
                                           analyzeKinFit        +
                                           analyzeGenMatch      +
                                           analyzeMVADisc
                                           )



