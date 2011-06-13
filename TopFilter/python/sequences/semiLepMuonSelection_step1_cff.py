import FWCore.ParameterSet.Config as cms

from TopAnalysis.TopAnalyzer.MuonAnalyzer_cfi import analyzeMuon
## -------------------------------------Muon eta  |  Muon pt  | Jets eta  |  Jets pt
semiLepMuon00 = analyzeMuon.clone()  ##     -     |     -     |     -     |     -
semiLepMuon10 = analyzeMuon.clone()  ##     x     |     -     |     -     |     -
semiLepMuon20 = analyzeMuon.clone()  ##     x     |     x     |     -     |     -
semiLepMuon01 = analyzeMuon.clone()  ##     -     |     -     |     x     |     -
semiLepMuon02 = analyzeMuon.clone()  ##     -     |     -     |     x     |     x
semiLepMuon11 = analyzeMuon.clone()  ##     x     |     -     |     x     |     -
semiLepMuon12 = analyzeMuon.clone()  ##     x     |     -     |     x     |     x
semiLepMuon21 = analyzeMuon.clone()  ##     x     |     x     |     x     |     -
semiLepMuon22 = analyzeMuon.clone()  ##     x     |     x     |     x     |     x

## put any extra configuration here:

from TopAnalysis.TopAnalyzer.JetAnalyzer_cfi  import analyzeJets
## -------------------------------------Muon eta  |  Muon pt  | Jets eta  |  Jets pt
semiLepJets00 = analyzeJets.clone()  ##     -     |     -     |     -     |     -
semiLepJets10 = analyzeJets.clone()  ##     x     |     -     |     -     |     -
semiLepJets20 = analyzeJets.clone()  ##     x     |     x     |     -     |     -
semiLepJets01 = analyzeJets.clone()  ##     -     |     -     |     x     |     -
semiLepJets02 = analyzeJets.clone()  ##     -     |     -     |     x     |     x
semiLepJets11 = analyzeJets.clone()  ##     x     |     -     |     x     |     -
semiLepJets12 = analyzeJets.clone()  ##     x     |     -     |     x     |     x
semiLepJets21 = analyzeJets.clone()  ##     x     |     x     |     x     |     -
semiLepJets22 = analyzeJets.clone()  ##     x     |     x     |     x     |     x

## put any extra configuration here:


## import selection cuts here
from TopAnalysis.TopFilter.selections.semiLepMuonSelection_step1_cff import *

from TopAnalysis.TopFilter.filters.EtaEventFilter_cfi import filterMuonEta
filterSemiLep10 = filterMuonEta.clone()
from TopAnalysis.TopFilter.filters.EtaEventFilter_cfi import filterJetsEta
filterSemiLep01 = filterJetsEta.clone()
from TopAnalysis.TopFilter.filters.PtEventFilter_cfi  import filterMuonPt
filterSemiLep20 = filterMuonPt.clone()
from TopAnalysis.TopFilter.filters.PtEventFilter_cfi  import filterJetsPt
filterSemiLep02 = filterJetsPt.clone()

## put any extra configuration here:
filterSemiLep10.cuts = semiLepMuonEta
filterSemiLep01.cuts = semiLepJetsEta
filterSemiLep20.cuts = semiLepMuonPt
filterSemiLep02.cuts = semiLepJetsPt

## define sequences
selectSemiLepMuon = cms.Sequence(semiLepMuon00 + semiLepJets00 *
                                 filterSemiLep10 *
                                 semiLepMuon10 + semiLepJets10 *
                                 filterSemiLep20 *
                                 semiLepMuon20 + semiLepJets20 *
                                 filterSemiLep01 *
                                 semiLepMuon21 + semiLepJets21 *
                                 filterSemiLep02 *
                                 semiLepMuon22 + semiLepJets22
                                 )
