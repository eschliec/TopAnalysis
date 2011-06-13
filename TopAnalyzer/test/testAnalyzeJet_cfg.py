import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for tqaflayer1 & 2 production from
# fullsim for semi-leptonic ttbar events 
#-------------------------------------------------
process = cms.Process("TEST")

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.FwkReport.reportEvery = 100

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #PAT test sample
    'file:/afs/cern.ch/cms/PRS/top/cmssw-data/relval200-for-pat-testing/FullSimTTBar-2_1_X_2008-07-08_STARTUP_V4-AODSIM.100.root'
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

## configure geometry
process.load("Configuration.StandardSequences.Geometry_cff")

## configure conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V4::All')

## Magnetic field now needs to be in the high-level py
process.load("Configuration.StandardSequences.MagneticField_cff")

#-------------------------------------------------
# tqaf configuration; if you want just to produce 
# tqafLayer2 on top of an already existing
# tqafLayer1 just comment the standard tqafLayer1
# production sequence
#-------------------------------------------------

## std sequence for tqaf layer1
process.load("TopQuarkAnalysis.TopObjectProducers.tqafLayer1_full_cff")
process.p0 = cms.Path(process.tqafLayer1)
#-------------------------------------------------
# private uni Hamburg analysis code
#-------------------------------------------------

## add event weight information
from TopAnalysis.TopUtils.EventWeightPlain_cfi import *
process.eventWeight = eventWeight

# analyze jets
from TopAnalysis.TopAnalyzer.JetAnalyzer_cfi import analyzeJets
process.analyzeJets = analyzeJets

process.analyzeAllJets = process.analyzeJets.clone()
process.analyzeSelJets = process.analyzeJets.clone()

## reconfigure
process.analyzeAllJets.input = 'allLayer1Jets'
process.analyzeSelJets.input = 'selectedLayer1Jets'

## register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('analyzeJets.root')
)

process.p1 = cms.Path(process.eventWeight    *
                      process.analyzeAllJets +
                      process.analyzeSelJets
                      )
