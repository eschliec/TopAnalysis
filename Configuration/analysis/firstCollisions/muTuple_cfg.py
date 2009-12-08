import FWCore.ParameterSet.Config as cms

#-------------------------------------------------
# test cfg file for muTuple production
#-------------------------------------------------

process = cms.Process("MuTuple")

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
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_1.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_10.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_11.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_12.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_13.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_14.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_15.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_16.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_17.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_18.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_19.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_2.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_20.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_21.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_22.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_23.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_24.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_25.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_26.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_27.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_28.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_29.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_3.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_30.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_31.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_32.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_33.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_34.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_35.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_36.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_37.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_38.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_39.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_4.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_40.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_41.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_42.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_43.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_44.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_45.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_46.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_47.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_48.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_49.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_5.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_50.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_51.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_52.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_53.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_54.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_55.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_56.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_57.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_58.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_59.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_6.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_60.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_61.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_62.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_63.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_64.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_65.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_66.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_67.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_68.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_69.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_7.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_70.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_71.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_72.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_73.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_74.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_75.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_76.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_77.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_78.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_79.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_8.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_80.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_81.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_82.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_83.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_84.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_85.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_86.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_87.root',
    '/store/user/snaumann/firstCollisions/patTuple_Run123596_9.root'
    )
)
process.source.firstRun = cms.untracked.uint32(123596)
process.source.lastRun  = cms.untracked.uint32(123596)

## define maximal number of events to loop over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

## configure process options
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

#-------------------------------------------------
# jet selection
#-------------------------------------------------

## select muons
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
process.goodMuons = selectedLayer1Muons.clone(src = 'selectedLayer1Muons',
                                              cut = 'combinedMuon.isNull = 0'
                                              )

## select events with good muons
from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *
process.countMuons = countLayer1Muons.clone(src = 'goodMuons',
                                            minNumber = 1
                                            )

#-------------------------------------------------
# process output
#-------------------------------------------------

## process path
process.p1 = cms.Path(process.goodMuons *
                      process.countMuons)

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    )
)

## configure output module
process.out = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    outputCommands = cms.untracked.vstring('keep *'),
    dropMetaDataForDroppedData = cms.untracked.bool(True),                                     
    fileName = cms.untracked.string('muTuple_Run123596.root')
)

process.outpath = cms.EndPath(process.out)
