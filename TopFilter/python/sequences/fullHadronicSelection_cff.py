import FWCore.ParameterSet.Config as cms

## all leptons kept for possible rejection of events later on (not yet implemented)

## jet selector
from TopAnalysis.TopFilter.sequences.jetSelection_cff import *
## muon selector
#from TopAnalysis.TopFilter.sequences.muonSelection_cff import *

## jet selector
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
## muon selector
#from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
## electron selector
#from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *

## jet count filter
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
## muon count filter
#from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
## electron count filter
#from PhysicsTools.PatAlgos.selectionLayer1.electronCountFilter_cfi import *

## ---
##    setup the collections for the fully-hadronic event selection
## ---

## setup the jet selection collection
tightLeadingJets = selectedPatJets.clone(src = 'goodJets',
                                         cut = 'abs(eta) < 2.4 & pt > 40'
                                         )
tightBottomJets  = selectedPatJets.clone(src = 'trackCountingHighPurBJets',
                                         cut = 'abs(eta) < 2.4 & pt > 50'
                                         )

## setting up the collections for the fully-hadronic
## event selection; on these collection monitoring
## can still be performed
fullHadronicSelection = cms.Sequence(goodJets *
                                     trackCountingHighPurBJets *
                                     tightLeadingJets *
                                     tightBottomJets
                                     )

## ---
##    configure the cutflow scenario
## ---

## jet quality analyzer
from TopAnalysis.TopAnalyzer.JetQuality_cfi import *
## muon quality analyzer
#from TopAnalysis.TopAnalyzer.MuonQuality_cfi import *
## jet kinematics analyzer
from TopAnalysis.TopAnalyzer.JetKinematics_cfi import *
## muon kinematics analyzer
#from TopAnalysis.TopAnalyzer.MuonKinematics_cfi import *
## event shape analyzer
from TopAnalysis.TopAnalyzer.EventShapes_cfi import *
## genParticle analyzer
from TopAnalysis.TopAnalyzer.GenParticle_cfi import *
## analyzer for special variables in full hadronic channel
from TopAnalysis.TopAnalyzer.FullHadSpecial_cfi import *
## kinfit analyzer
from TopAnalysis.TopAnalyzer.KinFitQuality_cfi import *
from TopAnalysis.TopAnalyzer.KinFitImprover_cfi import *
## analyzer for fully hadronic event reco
from TopAnalysis.TopAnalyzer.FullHadTopReco_cfi import *
## high level trigger filter
from TopAnalysis.TopFilter.sequences.triggerFilter_cff import *
## semileptonic selection
from TopAnalysis.TopFilter.sequences.fullHadronicSelection_cff import *
## generator matching
from TopAnalysis.TopFilter.sequences.generatorMatching_cff import *
## kinFit producer
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullHadEvtBuilder_cff import *
from TopAnalysis.TopAnalyzer.FullHadHypothesisAnalyzer_cff import *

## to switch verbosity modes of the kinFit
#ttFullHadEvent.verbosity = 1

## configuration of kinematic fit
kinFitTtFullHadEventHypothesis.maxNComb = -1

kinFitTtFullHadEventHypothesis.bTags = 2
kinFitTtFullHadEventHypothesis.bTagAlgo = 'trackCountingHighPurBJetTags'
kinFitTtFullHadEventHypothesis.minBTagValueBJet    = 2.17 #2.02
kinFitTtFullHadEventHypothesis.maxBTagValueNonBJet = 4.31 #3.4

#setForAllTtFullHadHypotheses(process, 'maxNJets', -1)
kinFitTtFullHadEventHypothesis.maxNJets = -1
ttFullHadJetPartonMatch.maxNJets        = -1

#setForAllTtFullHadHypotheses(process, 'jetCorrectionLevel', 'had')
kinFitTtFullHadEventHypothesis.jetCorrectionLevel = 'had'
ttFullHadHypGenMatch.jetCorrectionLevel           = 'had'

#setForAllTtFullHadHypotheses(process, 'jets', 'tightLeadingJets')
kinFitTtFullHadEventHypothesis.jets = 'tightLeadingJets'
ttFullHadJetPartonMatch.jets        = 'tightLeadingJets'
ttFullHadHypGenMatch.jets           = 'tightLeadingJets'
ttFullHadHypKinFit.jets             = 'tightLeadingJets'

## define ordered jets
uds0    = cms.PSet(index = cms.int32(0), correctionLevel = cms.string('abs'), flavor = cms.string("uds") )
uds1    = cms.PSet(index = cms.int32(1), correctionLevel = cms.string('abs'), flavor = cms.string("uds") )
uds2    = cms.PSet(index = cms.int32(2), correctionLevel = cms.string('abs'), flavor = cms.string("uds") )
uds3    = cms.PSet(index = cms.int32(3), correctionLevel = cms.string('abs'), flavor = cms.string("uds") )
uds4    = cms.PSet(index = cms.int32(4), correctionLevel = cms.string('abs'), flavor = cms.string("uds") )
uds5    = cms.PSet(index = cms.int32(5), correctionLevel = cms.string('abs'), flavor = cms.string("uds") )
bottom0 = cms.PSet(index = cms.int32(0), correctionLevel = cms.string('abs'), flavor = cms.string("b")   )
bottom1 = cms.PSet(index = cms.int32(1), correctionLevel = cms.string('abs'), flavor = cms.string("b")   )

## ---
##    MONITOR STEP 0
## ---

## JET KINEMATICS

## kinematics analyzers
tightBottomJetKinematics_0  = analyzeJetKinematics.clone (src = 'trackCountingHighPurBJets' )
tightLeadingJetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets')
tightLead_0_JetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets', analyze = uds0 )
tightLead_1_JetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets', analyze = uds1 )
tightLead_2_JetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets', analyze = uds2 )
tightLead_3_JetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets', analyze = uds3 )
tightLead_4_JetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets', analyze = uds4 )
tightLead_5_JetKinematics_0 = analyzeJetKinematics.clone (src = 'goodJets', analyze = uds5 )
tightBJet_0_JetKinematics_0 = analyzeJetKinematics.clone (src = 'trackCountingHighPurBJets' , analyze = bottom0)
tightBJet_1_JetKinematics_0 = analyzeJetKinematics.clone (src = 'trackCountingHighPurBJets' , analyze = bottom1)

## collect kinematics analyzers
monitorJetsKinematics_0 = cms.Sequence(tightBottomJetKinematics_0  +
                                       tightBJet_0_JetKinematics_0 +
                                       tightBJet_1_JetKinematics_0 +
                                       tightLeadingJetKinematics_0 +
                                       tightLead_0_JetKinematics_0 +
                                       tightLead_1_JetKinematics_0 +
                                       tightLead_2_JetKinematics_0 +
                                       tightLead_3_JetKinematics_0 +
                                       tightLead_4_JetKinematics_0 +
                                       tightLead_5_JetKinematics_0  
                                       )

## JET QUALITY

## quality analyzers
tightBottomJetQuality_0  = analyzeJetQuality.clone (src = 'trackCountingHighPurBJets' )
tightLeadingJetQuality_0 = analyzeJetQuality.clone (src = 'goodJets')
tightLead_0_JetQuality_0 = analyzeJetQuality.clone (src = 'goodJets', analyze = uds0 )
tightLead_1_JetQuality_0 = analyzeJetQuality.clone (src = 'goodJets', analyze = uds1 )
tightLead_2_JetQuality_0 = analyzeJetQuality.clone (src = 'goodJets', analyze = uds2 )
tightLead_3_JetQuality_0 = analyzeJetQuality.clone (src = 'goodJets', analyze = uds3 )
tightLead_4_JetQuality_0 = analyzeJetQuality.clone (src = 'goodJets', analyze = uds4 )
tightLead_5_JetQuality_0 = analyzeJetQuality.clone (src = 'goodJets', analyze = uds5 )
tightBJet_0_JetQuality_0 = analyzeJetQuality.clone (src = 'trackCountingHighPurBJets' , analyze = bottom0)
tightBJet_1_JetQuality_0 = analyzeJetQuality.clone (src = 'trackCountingHighPurBJets' , analyze = bottom1)

## collect quality analyzers
monitorJetsQuality_0 = cms.Sequence(tightBottomJetQuality_0  +
                                    tightBJet_0_JetQuality_0 +
                                    tightBJet_1_JetQuality_0 +
                                    tightLeadingJetQuality_0 +
                                    tightLead_0_JetQuality_0 +
                                    tightLead_1_JetQuality_0 +
                                    tightLead_2_JetQuality_0 +
                                    tightLead_3_JetQuality_0 +
                                    tightLead_4_JetQuality_0 +
                                    tightLead_5_JetQuality_0  
                                    )

## EVENT SHAPES

## collect event shape analyzers
eventShapes_0 = analyzeEventShapes.clone( src = 'goodJets' )

## monitor sequence for event shape analyzers
monitorEventShapes_0 = cms.Sequence( eventShapes_0 )

## FULL HAD SPECIAL

## collect analyzers specially for full hadronic analysis
fullHadSpecial_0 = analyzeFullHadSpecials.clone( src = 'goodJets' )

## monitor sequence for specially for full hadronic analyzers
monitorFullHadSpecials_0 = cms.Sequence( fullHadSpecial_0 )

## GEN PARTICLE

## collect analyzers for genParticles
genParticles_0 = analyzeGenParticles.clone()

## monitor sequence for genParticles
monitorGenParticles_0 = cms.Sequence( genParticles_0 )

## ---
##    FILTER STEP 1
## ---

## select all events with at least 6 jets
leadingJetSelection = countPatJets.clone( src = 'tightLeadingJets',
                                          minNumber = 6
                                         )

## ---
##    MONITOR STEP 1
## ---

## JET KINEMATICS

## kinematics analyzers
tightBottomJetKinematics_1  = analyzeJetKinematics.clone (src = 'tightBottomJets' )
tightLeadingJetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets')
tightLead_0_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds0 )
tightLead_1_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds1 )
tightLead_2_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds2 )
tightLead_3_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds3 )
tightLead_4_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds4 )
tightLead_5_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds5 )
tightBJet_0_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightBottomJets' , analyze = bottom0)
tightBJet_1_JetKinematics_1 = analyzeJetKinematics.clone (src = 'tightBottomJets' , analyze = bottom1)

## collect kinematics analyzers
monitorJetsKinematics_1 = cms.Sequence(tightBottomJetKinematics_1  +
                                       tightBJet_0_JetKinematics_1 +
                                       tightBJet_1_JetKinematics_1 +
                                       tightLeadingJetKinematics_1 +
                                       tightLead_0_JetKinematics_1 +
                                       tightLead_1_JetKinematics_1 +
                                       tightLead_2_JetKinematics_1 +
                                       tightLead_3_JetKinematics_1 +
                                       tightLead_4_JetKinematics_1 +
                                       tightLead_5_JetKinematics_1  
                                       )

## JET QUALITY

## quality analyzers
tightBottomJetQuality_1  = analyzeJetQuality.clone (src = 'tightBottomJets' )
tightLeadingJetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets')
tightLead_0_JetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds0 )
tightLead_1_JetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds1 )
tightLead_2_JetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds2 )
tightLead_3_JetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds3 )
tightLead_4_JetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds4 )
tightLead_5_JetQuality_1 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds5 )
tightBJet_0_JetQuality_1 = analyzeJetQuality.clone (src = 'tightBottomJets' , analyze = bottom0)
tightBJet_1_JetQuality_1 = analyzeJetQuality.clone (src = 'tightBottomJets' , analyze = bottom1)

## collect quality analyzers
monitorJetsQuality_1 = cms.Sequence(tightBottomJetQuality_1  +
                                    tightBJet_0_JetQuality_1 +
                                    tightBJet_1_JetQuality_1 +
                                    tightLeadingJetQuality_1 +
                                    tightLead_0_JetQuality_1 +
                                    tightLead_1_JetQuality_1 +
                                    tightLead_2_JetQuality_1 +
                                    tightLead_3_JetQuality_1 +
                                    tightLead_4_JetQuality_1 +
                                    tightLead_5_JetQuality_1  
                                    )

## EVENT SHAPES

## collect event shape analyzers
eventShapes_1 = analyzeEventShapes.clone( src = 'tightLeadingJets' )

## monitor sequence for event shape analyzers
monitorEventShapes_1 = cms.Sequence( eventShapes_1 )

## FULL HAD SPECIAL

## collect analyzers specially for full hadronic analysis
fullHadSpecial_1 = analyzeFullHadSpecials.clone( src = 'tightLeadingJets' )

## monitor sequence for specially for full hadronic analyzers
monitorFullHadSpecials_1 = cms.Sequence( fullHadSpecial_1 )

## GEN PARTICLE

## collect analyzers for genParticles
genParticles_1 = analyzeGenParticles.clone()

## monitor sequence for genParticles
monitorGenParticles_1 = cms.Sequence( genParticles_1 )

## ---
##    FILTER STEP 2
## ---

## select events with at least 2 b jets
bottomJetSelection  = countPatJets.clone( src = 'tightBottomJets',
                                          minNumber = 2
                                         )

## ---
##    MONITOR STEP 2
## ---

## JET KINEMATICS

## collect kinematics analyzers
tightBottomJetKinematics_2  = analyzeJetKinematics.clone (src = 'tightBottomJets' )
tightLeadingJetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets')
tightLead_0_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds0 )
tightLead_1_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds1 )
tightLead_2_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds2 )
tightLead_3_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds3 )
tightLead_4_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds4 )
tightLead_5_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds5 )
tightBJet_0_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightBottomJets' , analyze = bottom0)
tightBJet_1_JetKinematics_2 = analyzeJetKinematics.clone (src = 'tightBottomJets' , analyze = bottom1)

## to be called with fullHadronicSelection
monitorJetsKinematics_2 = cms.Sequence(tightBottomJetKinematics_2  +
                                       tightBJet_0_JetKinematics_2 +
                                       tightBJet_1_JetKinematics_2 +
                                       tightLeadingJetKinematics_2 +
                                       tightLead_0_JetKinematics_2 +
                                       tightLead_1_JetKinematics_2 +
                                       tightLead_2_JetKinematics_2 +
                                       tightLead_3_JetKinematics_2 +
                                       tightLead_4_JetKinematics_2 +
                                       tightLead_5_JetKinematics_2  
                                       )

## JET QUALITY

## quality analyzers
tightBottomJetQuality_2  = analyzeJetQuality.clone (src = 'tightBottomJets' )
tightLeadingJetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets')
tightLead_0_JetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds0 )
tightLead_1_JetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds1 )
tightLead_2_JetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds2 )
tightLead_3_JetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds3 )
tightLead_4_JetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds4 )
tightLead_5_JetQuality_2 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds5 )
tightBJet_0_JetQuality_2 = analyzeJetQuality.clone (src = 'tightBottomJets' , analyze = bottom0)
tightBJet_1_JetQuality_2 = analyzeJetQuality.clone (src = 'tightBottomJets' , analyze = bottom1)

## collect quality analyzers
monitorJetsQuality_2 = cms.Sequence(tightBottomJetQuality_2  +
                                    tightBJet_0_JetQuality_2 +
                                    tightBJet_1_JetQuality_2 +
                                    tightLeadingJetQuality_2 +
                                    tightLead_0_JetQuality_2 +
                                    tightLead_1_JetQuality_2 +
                                    tightLead_2_JetQuality_2 +
                                    tightLead_3_JetQuality_2 +
                                    tightLead_4_JetQuality_2 +
                                    tightLead_5_JetQuality_2  
                                    )

## EVENT SHAPES

## collect event shape analyzers
eventShapes_2 = analyzeEventShapes.clone( src = 'tightLeadingJets' )

## monitor sequence for event shape analyzers
monitorEventShapes_2 = cms.Sequence( eventShapes_2 )

## KINFIT QUALITY AND FULL HAD TOP RECO

## kinfit quality analyzer
## collect kinfit quality analyzers
kinFitQuality_2  = analyzeKinFitQuality.clone  ( srcB = 'tightLeadingJets' )
kinFitImprover0_2 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(0) ) )
kinFitImprover1_2 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(1) ) )
kinFitImprover2_2 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(2) ) )
kinFitImprover3_2 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(3) ) )
kinFitImprover4_2 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(4) ) )

## collect fully hadronic top reco analyzers
fullHadTopReco_2 = analyzeFullHadTopReco.clone( srcB = 'tightLeadingJets' )

## monitor sequence for kinfit quality analyzers
monitorKinFit_2 = cms.Sequence( kinFitQuality_2   *
                                kinFitImprover0_2 *
                                kinFitImprover1_2 *
                                kinFitImprover2_2 *
                                kinFitImprover3_2 *
                                kinFitImprover4_2 *
                                fullHadTopReco_2
                               )

## FULL HAD SPECIAL

## collect analyzers specially for full hadronic analysis
fullHadSpecial_2 = analyzeFullHadSpecials.clone( src = 'tightLeadingJets' )

## monitor sequence for specially for full hadronic analyzers
monitorFullHadSpecials_2 = cms.Sequence( fullHadSpecial_2 )

## GEN PARTICLE

## collect analyzers for genParticles
genParticles_2 = analyzeGenParticles.clone()

## monitor sequence for genParticles
monitorGenParticles_2 = cms.Sequence( genParticles_2 )

## ---
##    FILTER STEP 3
## ---

## kinfit quality filter
from TopAnalysis.TopFilter.filters.KinFitQualityFilter_cfi import *
filterKinFitQuality = filterKinFitQuality.clone( srcB = 'tightLeadingJets', minProb = 0.01 )

## ---
##    MONITOR STEP 3
## ---

## JET KINEMATICS

## collect kinematics analyzers
tightBottomJetKinematics_3  = analyzeJetKinematics.clone (src = 'tightBottomJets' )
tightLeadingJetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets')
tightLead_0_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds0 )
tightLead_1_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds1 )
tightLead_2_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds2 )
tightLead_3_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds3 )
tightLead_4_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds4 )
tightLead_5_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightLeadingJets', analyze = uds5 )
tightBJet_0_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightBottomJets' , analyze = bottom0)
tightBJet_1_JetKinematics_3 = analyzeJetKinematics.clone (src = 'tightBottomJets' , analyze = bottom1)

## to be called with fullHadronicSelection
monitorJetsKinematics_3 = cms.Sequence(tightBottomJetKinematics_3  +
                                       tightBJet_0_JetKinematics_3 +
                                       tightBJet_1_JetKinematics_3 +
                                       tightLeadingJetKinematics_3 +
                                       tightLead_0_JetKinematics_3 +
                                       tightLead_1_JetKinematics_3 +
                                       tightLead_2_JetKinematics_3 +
                                       tightLead_3_JetKinematics_3 +
                                       tightLead_4_JetKinematics_3 +
                                       tightLead_5_JetKinematics_3  
                                       )

## JET QUALITY

## quality analyzers
tightBottomJetQuality_3  = analyzeJetQuality.clone (src = 'tightBottomJets' )
tightLeadingJetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets')
tightLead_0_JetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds0 )
tightLead_1_JetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds1 )
tightLead_2_JetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds2 )
tightLead_3_JetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds3 )
tightLead_4_JetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds4 )
tightLead_5_JetQuality_3 = analyzeJetQuality.clone (src = 'tightLeadingJets', analyze = uds5 )
tightBJet_0_JetQuality_3 = analyzeJetQuality.clone (src = 'tightBottomJets' , analyze = bottom0)
tightBJet_1_JetQuality_3 = analyzeJetQuality.clone (src = 'tightBottomJets' , analyze = bottom1)

## collect quality analyzers
monitorJetsQuality_3 = cms.Sequence(tightBottomJetQuality_3  +
                                    tightBJet_0_JetQuality_3 +
                                    tightBJet_1_JetQuality_3 +
                                    tightLeadingJetQuality_3 +
                                    tightLead_0_JetQuality_3 +
                                    tightLead_1_JetQuality_3 +
                                    tightLead_2_JetQuality_3 +
                                    tightLead_3_JetQuality_3 +
                                    tightLead_4_JetQuality_3 +
                                    tightLead_5_JetQuality_3  
                                    )

## EVENT SHAPES

## collect event shape analyzers
eventShapes_3 = analyzeEventShapes.clone( src = 'tightLeadingJets' )

## monitor sequence for event shape analyzers
monitorEventShapes_3 = cms.Sequence( eventShapes_3 )

## KINFIT QUALITY AND FULL HAD TOP RECO

## kinfit quality analyzer
## collect kinfit quality analyzers
kinFitQuality_3  = analyzeKinFitQuality.clone  ( srcB = 'tightLeadingJets' )
kinFitImprover0_3 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(0) ) )
kinFitImprover1_3 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(1) ) )
kinFitImprover2_3 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(2) ) )
kinFitImprover3_3 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(3) ) )
kinFitImprover4_3 = analyzeKinFitImprover.clone( srcB = 'tightLeadingJets' , analyze = cms.PSet( comboType = cms.uint32(4) ) )

## collect fully hadronic top reco analyzers
fullHadTopReco_3 = analyzeFullHadTopReco.clone( srcB = 'tightLeadingJets' )

## monitor sequence for kinfit quality analyzers
monitorKinFit_3 = cms.Sequence( kinFitQuality_3  *
                                kinFitImprover0_3 *
                                kinFitImprover1_3 *
                                kinFitImprover2_3 *
                                kinFitImprover3_3 *
                                kinFitImprover4_3 *
                                fullHadTopReco_3
                               )

## FULL HAD SPECIAL

## collect analyzers specially for full hadronic analysis
fullHadSpecial_3 = analyzeFullHadSpecials.clone( src = 'tightLeadingJets' )

## monitor sequence for specially for full hadronic analyzers
monitorFullHadSpecials_3 = cms.Sequence( fullHadSpecial_3 )

## GEN PARTICLE

## collect analyzers for genParticles
genParticles_3 = analyzeGenParticles.clone()

## monitor sequence for genParticles
monitorGenParticles_3 = cms.Sequence( genParticles_3 )

## ---
##    FILTER STEP 4 (not used)
## ---

## event shape filter
from TopAnalysis.TopFilter.filters.EventShapeFilter_cfi import *
filterEventShapes = filterEventShape.clone( minC = 0.75 )

## ---
##    MONITOR STEP 4
## ---

## To be added in time


## ---
##    run the final sequence
## ---
analyseFullHadronicSelection = cms.Sequence(## do the hlt triggering
                                            hltQuadJet30          *
                                            #hltHt200              *
                                            ## do the selections
                                            fullHadronicSelection *
                                            ## do the matching
                                            matchJetsToPartons    *
                                            ## do the monitoring
                                            monitorJetsKinematics_0  *
                                            monitorJetsQuality_0     *
                                            monitorEventShapes_0     *
                                            monitorFullHadSpecials_0 *
                                            monitorGenParticles_0    *
                                            ## do the 1. event selection
                                            leadingJetSelection      *
                                            monitorJetsKinematics_1  *
                                            monitorJetsQuality_1     *
                                            monitorEventShapes_1     *
                                            monitorFullHadSpecials_1 *
                                            monitorGenParticles_1    *
                                            ## do the 2. event selection
                                            bottomJetSelection       *
                                            makeTtFullHadEvent       *
                                            monitorKinFit_2          *
                                            monitorJetsKinematics_2  *
                                            monitorJetsQuality_2     *
                                            monitorEventShapes_2     *
                                            monitorFullHadSpecials_2 *
                                            monitorGenParticles_2    *
                                            ## do the 3. event selection
                                            filterKinFitQuality      *
                                            monitorKinFit_3          *
                                            monitorJetsKinematics_3  *
                                            monitorJetsQuality_3     *
                                            monitorEventShapes_3     *
                                            monitorFullHadSpecials_3 *
                                            monitorGenParticles_3
                                            )

## ---
##    provide a function to disable parts of the selection
## ---
def disableCountFilter(whichCountFilter):
    whichCountFilter.minNumber = 0
    whichCountFilter.maxNumber = 999999

## ---
##    run on real data
## ---
def runOnRealData(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'removing all elements from the sequence '
    print '*analyseFullHadronicSelection* that rely '
    print 'on generator information to run properly '
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.analyseFullHadronicSelection.remove(process.matchJetsToPartons)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_0)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_1)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_2)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_3)

## ---
##    remove modules that produce monitoring plots during the cutflow
## ---
def removeMonitoringOfCutflow(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'removing all monitoring elements from the '
    print 'sequence *analyseFullHadronicSelection* so'
    print 'only a pure selection of events is done '
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.analyseFullHadronicSelection.remove(process.monitorJetsKinematics_0)
    process.analyseFullHadronicSelection.remove(process.monitorJetsKinematics_1)
    process.analyseFullHadronicSelection.remove(process.monitorJetsKinematics_2)
    process.analyseFullHadronicSelection.remove(process.monitorJetsKinematics_3)
    process.analyseFullHadronicSelection.remove(process.monitorJetsQuality_0)
    process.analyseFullHadronicSelection.remove(process.monitorJetsQuality_1)
    process.analyseFullHadronicSelection.remove(process.monitorJetsQuality_2)
    process.analyseFullHadronicSelection.remove(process.monitorJetsQuality_3)
    process.analyseFullHadronicSelection.remove(process.monitorEventShapes_0)
    process.analyseFullHadronicSelection.remove(process.monitorEventShapes_1)
    process.analyseFullHadronicSelection.remove(process.monitorEventShapes_2)
    process.analyseFullHadronicSelection.remove(process.monitorEventShapes_3)
    process.analyseFullHadronicSelection.remove(process.monitorFullHadSpecials_0)
    process.analyseFullHadronicSelection.remove(process.monitorFullHadSpecials_1)
    process.analyseFullHadronicSelection.remove(process.monitorFullHadSpecials_2)
    process.analyseFullHadronicSelection.remove(process.monitorFullHadSpecials_3)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_0)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_1)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_2)
    process.analyseFullHadronicSelection.remove(process.monitorGenParticles_3)
    process.analyseFullHadronicSelection.remove(process.monitorKinFit_2)
    process.analyseFullHadronicSelection.remove(process.monitorKinFit_3)

## ---
##    remove default trigger
## ---
def removeDefaultTrigger(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'removing the default trigger from '
    print 'standard fully hadronic event selection '
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.analyseFullHadronicSelection.remove(process.hltHt200)
    process.analyseFullHadronicSelection.remove(process.hltQuadJet30)
    
## ---
##    switch all necessary filters to run this sequence for background estimation
## ---
def runAsBackgroundEstimation(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'switching *analyseFullHadronicSelection* to'
    print 'background estimation'
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.bottomJetSelection.minNumber = 0
    process.bottomJetSelection.maxNumber = 0
    process.kinFitTtFullHadEventHypothesis.bTags = 0

## ---
##    run analysis on PFJets instead of caloJets
## ---
def runOnPF(process):
    print '++++++++++++++++++++++++++++++++++++++++++++'
    print 'switching all inputs in to run on PFJets'
    print 'instead of caloJets: goodJets -> goodJetsPF'
    print '++++++++++++++++++++++++++++++++++++++++++++'
    process.analyseFullHadronicSelection.replace(process.goodJets, process.goodJetsPF)
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.analyseFullHadronicSelection, 'goodJets', 'goodJetsPF')

## ---
##    switch to simpleSecondaryVertex bTagger
## ---
def switchToSSV(process):
    process.analyseFullHadronicSelection.replace(process.trackCountingHighPurBJets, process.simpleSecondaryVertexBJets)
    process.kinFitTtFullHadEventHypothesis.bTagAlgo            = 'simpleSecondaryVertexBJetTags'
    process.kinFitTtFullHadEventHypothesis.minBTagValueBJet    = 2.02 #2.17
    process.kinFitTtFullHadEventHypothesis.maxBTagValueNonBJet = 3.4  #4.31
    from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
    massSearchReplaceAnyInputTag(process.analyseFullHadronicSelection, 'trackCountingHighPurBJets', 'simpleSecondaryVertexBJets')

