import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for patTuple production
#-------------------------------------------------

process = cms.Process("PatTuple")

#-------------------------------------------------
# process configuration
#-------------------------------------------------

## configure message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('topFilter')
process.MessageLogger.cerr.topFilter = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## define input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/0A71AE7F-4DE2-DE11-8B2F-001D09F251CC.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/14D22A92-62E2-DE11-9B14-000423D990CC.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/263E80C6-41E2-DE11-A194-001617C3B66C.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/26EC6965-4CE2-DE11-ABCF-003048D373AE.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/2A6903C0-68E2-DE11-B6BA-001D09F28D54.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/2EC732CC-5CE2-DE11-A781-001D09F290BF.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/4693BF9A-40E2-DE11-BDBD-000423D944F8.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/54A3A373-4CE2-DE11-8658-000423D99AAA.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/5C4B3A8E-63E2-DE11-A02A-000423D99E46.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/68ED95B5-50E2-DE11-B4C8-001D09F27003.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/76CCDA0D-8AE2-DE11-AF65-0030487A3C9A.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/7E34865F-45E2-DE11-896A-000423D98F98.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/7EFA67BE-66E2-DE11-AE17-001617C3B79A.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/903D2066-61E2-DE11-9F6E-0019B9F704D6.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/A0FC9BDF-65E2-DE11-A6A1-000423D174FE.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/A4D9D21B-58E2-DE11-8F7A-000423D986A8.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/B60CC58F-5CE2-DE11-9FC7-001D09F24DDF.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/D0586F90-5FE2-DE11-8976-001D09F24691.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/E432BCD7-55E2-DE11-B670-001617C3B6CC.root',
    '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/596/F67BCF17-48E2-DE11-98B1-000423D94534.root'
    )
)
process.source.firstRun = cms.untracked.uint32(123596)
process.source.lastRun  = cms.untracked.uint32(123596)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False)
)

## configure geometry & conditions
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('GR09_P_V7::All')

#-------------------------------------------------
# trigger configuration
#-------------------------------------------------

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)

## select the "physics bits"
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('40 OR 41')

#-------------------------------------------------
# pat configuration
#-------------------------------------------------

## std sequence for pat
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## use the correct jet energy corrections
process.jetCorrFactors.corrSample = "900GeV"

## switch off MC matching
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, 'All')

## add jet collections
from PhysicsTools.PatAlgos.tools.jetTools import *

#addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetIcone5'),
#                 'JPTc',
#                 doJTA        = True,
#                 doBTagging   = True,
#                 jetCorrLabel = ('IC5','JPT'),
#                 doType1MET   = True,
#                 doL1Counters = False,
#                 genJetCollection = cms.InputTag("iterativeCone5GenJets")
#                 )

#addJetCollection(process,cms.InputTag('sisCone5CaloJets'),
#                 'SC5',
#                 doJTA        = True,
#                 doBTagging   = True,
#                 jetCorrLabel = ('SC5','Calo'),
#                 doType1MET   = True,
#                 doL1Counters = False,
#                 genJetCollection=cms.InputTag("sisCone5GenJets")
#                 )

#-------------------------------------------------
# private Uni Hamburg analysis code
#-------------------------------------------------

## add event weight information
process.load("TopAnalysis.TopUtils.EventWeightPlain_cfi")
process.wght = cms.Path(process.eventWeight)   

## basic event selection
#process.load("TopAnalysis.TopFilter.sequences.semiLepMuonSelection_step0_cff")
#process.semiLepMuonSelection = cms.Path(process.slmFilterFullMon )   

#process.load("TopAnalysis.TopFilter.sequences.fullLepElMuSelection_step0_cff")
#process.fullLepElMuSelection = cms.Path(process.flemFilterFullMon)

#process.load("TopAnalysis.TopFilter.sequences.fullLepMuonSelection_step0_cff")
#process.fullLepMuonSelection = cms.Path(process.flmmFilterFullMon)

## register TFileService
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('commonTemplate.root')
)

#-------------------------------------------------
# process output:
# - first the event selection is defined: only
# those events that have passed the full
# production path are selected and written to file
# - for the event content, ALL objects
# are dropped before the patTuple content is added
#-------------------------------------------------

## process path
process.p1 = cms.Path(process.hltLevel1GTSeed *
                      process.patDefaultSequence)

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    )
)

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    outputCommands = cms.untracked.vstring('drop *'),
    dropMetaDataForDroppedData = cms.untracked.bool(True),                                     
    fileName = cms.untracked.string('patTuple_Run123596.root')
)
# save pat output
from PhysicsTools.PatAlgos.patEventContent_cff import *
process.out.outputCommands += patTriggerEventContent
process.out.outputCommands += patExtraAodEventContent
process.out.outputCommands += patEventContentTriggerMatch
process.out.outputCommands += patEventContentNoLayer1Cleaning
process.out.outputCommands += ["keep *_selectedLayer1Jets*_*_*",
                               "keep *_layer1METs*_*_*"]
process.out.outputCommands += ["keep *_eventWeight_*_*"]
# drop stuff which is not needed
process.out.outputCommands += ["drop *_towerMaker_*_*"]

process.outpath = cms.EndPath(process.out)
