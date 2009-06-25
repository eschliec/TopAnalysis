import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# cfg file for the analysis of selected muon
# quantities
#-------------------------------------------------
process = cms.Process("Muon")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
# ignoreTotal=cms.untracked.int32(1),
# oncePerEventMode=cms.untracked.bool(False)
#)

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    ## add your favourite file here
    )
)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

#-------------------------------------------------
# muon analysis
#-------------------------------------------------

## analyze muon isolation
process.load("TopAnalysis.TopAnalyzer.IsolationAnalyzer_cfi")
from TopAnalysis.TopAnalyzer.IsolationAnalyzer_cfi import analyzeisolationMET
process.before = analyzeisolationMET.clone()
process.after = analyzeisolationMET.clone()
process.before.modulename = "findQCDBkgMVA"
process.after.modulename = "findQCDBkgMVA"
process.load("TopAnalysis.TopUtils.QCDBkgMVAComputer_cff")
## import selection cuts here
from TopAnalysis.TopFilter.selections.semiLepMuonSelection_step1_cff import *
from TopAnalysis.TopFilter.filters.SemiLepMuonEventFilter_cfi import filterSemiLepMuonEvent
filterSemiLepMuonEvent.lepEtaFilter = semiLepMuonEta
filterSemiLepMuonEvent.lepPtFilter = semiLepMuonPt
filterSemiLepMuonEvent.lepPtFilter.min = [20.0]
filterSemiLepMuonEvent.jetEtaFilter = semiLepJetsEta
filterSemiLepMuonEvent.jetPtFilter = semiLepJetsPt
#use lower cuts for higher QCD statistics
filterSemiLepMuonEvent.jetPtFilter.min = [20.0, 20.0, 20.0, 20.0]
filterSemiLepMuonEvent.trkIsoFilter = semiLepMuonIsolationTrk
filterSemiLepMuonEvent.calIsoFilter = semiLepMuonIsolationCal
#filterSemiLepMuonEvent.jetDistFilter = semiLepMuonJetDistance
filterSemiLepMuonEvent.combIsoFilter = semiLepMuonIsolationRelComb
filterSemiLepMuonEvent.combIsoFilter.max = [0.3]
#process.load("TopAnalysis.TopFilter.sequences.semiLepMuonSelection_step1_cff")
process.tight = filterSemiLepMuonEvent.clone()
process.loose = filterSemiLepMuonEvent.clone()
process.loose.relComb = False
## std sequence to produce the ttGenEvt
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
## configure ttGenEventFilters
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEventFilters_cff")
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.electron = False
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.muon     = True
process.ttSemiLeptonicFilter.allowedTopDecays.decayBranchA.tau      = False
process.load("TopAnalysis.TopAnalyzer.MuonAnalyzer_cfi")
process.analyzeAllMuons = process.analyzeMuon.clone()
process.analyzeIsolatedMuons = process.analyzeMuon.clone()
## register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeMuonIso.root')
)

process.p1 = cms.Path(process.makeGenEvt *
                      process.ttSemiLeptonicFilter * 
                      process.findQCDBkgMVA * process.analyzeAllMuons *process.loose * process.before * process.tight * process.analyzeIsolatedMuons * process.after)
