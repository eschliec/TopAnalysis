import FWCore.ParameterSet.Config as cms

from TopAnalysis.TopAnalyzer.MuonAnalyzer_cfi import analyzeMuon
## -------------------------------------Muon eta  |  Muon pt  | Jets eta  |  Jets pt
fullLepMuon00 = analyzeMuon.clone()    ##     -     |     -     |     -     |     -
fullLepMuon10 = analyzeMuon.clone()    ##     x     |     -     |     -     |     -
fullLepMuon20 = analyzeMuon.clone()    ##     x     |     x     |     -     |     -
fullLepMuon01 = analyzeMuon.clone()    ##     -     |     -     |     x     |     -
fullLepMuon02 = analyzeMuon.clone()    ##     -     |     -     |     x     |     x
fullLepMuon11 = analyzeMuon.clone()    ##     x     |     -     |     x     |     -
fullLepMuon12 = analyzeMuon.clone()    ##     x     |     -     |     x     |     x
fullLepMuon21 = analyzeMuon.clone()    ##     x     |     x     |     x     |     -
fullLepMuon22 = analyzeMuon.clone()    ##     x     |     x     |     x     |     x

## put any extra configuration here:

from TopAnalysis.TopAnalyzer.JetAnalyzer_cfi  import analyzeJets
## -------------------------------------Muon eta  |  Muon pt  | Jets eta  |  Jets pt
fullLepJets00 = analyzeJets.clone()    ##     -     |     -     |     -     |     -
fullLepJets10 = analyzeJets.clone()    ##     x     |     -     |     -     |     -
fullLepJets20 = analyzeJets.clone()    ##     x     |     x     |     -     |     -
fullLepJets01 = analyzeJets.clone()    ##     -     |     -     |     x     |     -
fullLepJets02 = analyzeJets.clone()    ##     -     |     -     |     x     |     x
fullLepJets11 = analyzeJets.clone()    ##     x     |     -     |     x     |     -
fullLepJets12 = analyzeJets.clone()    ##     x     |     -     |     x     |     x
fullLepJets21 = analyzeJets.clone()    ##     x     |     x     |     x     |     -
fullLepJets22 = analyzeJets.clone()    ##     x     |     x     |     x     |     x

## put any extra configuration here:


## import selection cuts here
from TopAnalysis.TopFilter.selections.fullLepMuonSelection_step0_cff import *

from TopAnalysis.TopFilter.filters.EtaEventFilter_cfi import filterMuonEta
filterFullLep10 = filterMuonEta.clone()
from TopAnalysis.TopFilter.filters.EtaEventFilter_cfi import filterJetsEta
filterFullLep01 = filterJetsEta.clone()
from TopAnalysis.TopFilter.filters.PtEventFilter_cfi  import filterMuonPt
filterFullLep20 = filterMuonPt.clone()
from TopAnalysis.TopFilter.filters.PtEventFilter_cfi  import filterJetsPt
filterFullLep02 = filterJetsPt.clone()

## put any extra configuration here:
filterFullLep10.cuts = fullLepMuonEta
filterFullLep01.cuts = fullLepJetsEta
filterFullLep20.cuts = fullLepMuonPt
filterFullLep02.cuts = fullLepJetsPt

## define sequences
selectFullLepMuon = cms.Sequence(fullLepMuon00 + fullLepJets00 *
                                 filterFullLep10 *
                                 fullLepMuon10 + fullLepJets10 *
                                 filterFullLep20 *
                                 fullLepMuon20 + fullLepJets20 *
                                 filterFullLep01 *
                                 fullLepMuon21 + fullLepJets21 *
                                 filterFullLep02 *
                                 fullLepMuon22 + fullLepJets22
                                 )


