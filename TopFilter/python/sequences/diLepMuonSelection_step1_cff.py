import FWCore.ParameterSet.Config as cms

from TopAnalysis.TopAnalyzer.MuonAnalyzer_cfi import analyzeMuon
## -------------------------------------Muon eta  |  Muon pt  |   Jet eta |  Jet pt   | muon iso  |  muon dist
diLepMuon00 = analyzeMuon.clone()    ##     -     |     -     |     -     |     -     |     -     |     -
diLepMuon10 = analyzeMuon.clone()    ##     x     |     -     |     -     |     -     |     -     |     -
diLepMuon20 = analyzeMuon.clone()    ##     x     |     x     |     -     |     -     |     -     |     -
diLepMuon01 = analyzeMuon.clone()    ##     -     |     -     |     x     |     -     |     -     |     -
diLepMuon02 = analyzeMuon.clone()    ##     -     |     -     |     x     |     x     |     -     |     -
diLepMuon11 = analyzeMuon.clone()    ##     x     |     -     |     x     |     -     |     -     |     -
diLepMuon12 = analyzeMuon.clone()    ##     x     |     -     |     x     |     x     |     -     |     -
diLepMuon21 = analyzeMuon.clone()    ##     x     |     x     |     x     |     -     |     -     |     -
diLepMuon22 = analyzeMuon.clone()    ##     x     |     x     |     x     |     x     |     -     |     -
diLepMuon32 = analyzeMuon.clone()    ##     x     |     x     |     x     |     x     |     x     |     -
diLepMuon42 = analyzeMuon.clone()    ##     x     |     x     |     x     |     x     |     x     |     x

## put any extra configuration here:

from TopAnalysis.TopAnalyzer.JetAnalyzer_cfi  import analyzeJets
## -------------------------------------Muon eta  |  Muon pt  | Jets eta  |  Jets pt
diLepJets00 = analyzeJets.clone()    ##     -     |     -     |     -     |     -     |     -     |     -
diLepJets10 = analyzeJets.clone()    ##     x     |     -     |     -     |     -     |     -     |     -
diLepJets20 = analyzeJets.clone()    ##     x     |     x     |     -     |     -     |     -     |     -
diLepJets01 = analyzeJets.clone()    ##     -     |     -     |     x     |     -     |     -     |     -
diLepJets02 = analyzeJets.clone()    ##     -     |     -     |     x     |     x     |     -     |     -
diLepJets11 = analyzeJets.clone()    ##     x     |     -     |     x     |     -     |     -     |     -
diLepJets12 = analyzeJets.clone()    ##     x     |     -     |     x     |     x     |     -     |     -
diLepJets21 = analyzeJets.clone()    ##     x     |     x     |     x     |     -     |     -     |     -
diLepJets22 = analyzeJets.clone()    ##     x     |     x     |     x     |     x     |     -     |     -
diLepJets32 = analyzeJets.clone()    ##     x     |     x     |     x     |     x     |     x     |     -
diLepJets42 = analyzeJets.clone()    ##     x     |     x     |     x     |     x     |     x     |     x

## put any extra configuration here:


## import selection cuts here
from TopAnalysis.TopFilter.selections.diLepMuonSelection_step1_cff import *

from TopAnalysis.TopFilter.filters.EtaEventFilter_cfi import filterMuonEta
filterDiLep10 = filterMuonEta.clone()
from TopAnalysis.TopFilter.filters.EtaEventFilter_cfi import filterJetsEta
filterDiLep01 = filterJetsEta.clone()
from TopAnalysis.TopFilter.filters.PtEventFilter_cfi  import filterMuonPt
filterDiLep20 = filterMuonPt.clone()
from TopAnalysis.TopFilter.filters.PtEventFilter_cfi  import filterJetsPt
filterDiLep02 = filterJetsPt.clone()

from TopAnalysis.TopFilter.filters.IsolationEventFilter_cfi  import filterMuonIsolation
filterDiLep30 = filterMuonIsolation.clone()
from TopAnalysis.TopFilter.filters.DistanceEventFilter_cfi  import filterMuonDistance
filterDiLep40 = filterMuonDistance.clone()

## put any extra configuration here:
filterDiLep10.cuts = diLepMuonEta
filterDiLep01.cuts = diLepJetsEta
filterDiLep20.cuts = diLepMuonPt
filterDiLep02.cuts = diLepJetsPt
filterDiLep30.cuts = diLepCombIsolation
filterDiLep40.cuts = diLepMuonJetDistance

## define sequences
makeDiMuonMuonKinematics = cms.Sequence(diLepMuon00 + diLepJets00 *
                                        filterDiLep10 *
                                        diLepMuon10 + diLepJets10 *
                                        filterDiLep20 *
                                        diLepMuon20 + diLepJets20 
				       )

makeDiMuonJetsKinematics = cms.Sequence(filterDiLep01 *
                                        diLepMuon21 + diLepJets21 *
                                        filterDiLep02 *
                                        diLepMuon22 + diLepJets22 
				       )

makeDiMuonMuonIsolation  = cms.Sequence(filterDiLep30 *
				        diLepMuon32 + diLepJets32 *
				        filterDiLep40 *
				        diLepMuon42 + diLepJets42 
				       )

selectDiLepMuon = cms.Sequence(makeDiMuonMuonKinematics *
                               makeDiMuonJetsKinematics * 
                               makeDiMuonMuonIsolation 				 				 
                              )
