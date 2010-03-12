import FWCore.ParameterSet.Config as cms

filterTrigger        = cms.EDFilter("NewTriggerTestFilter",
      whichTrigger   = cms.string("QuadJet80303030"),
      useEventWeight = cms.bool(False),
      weight         = cms.InputTag("eventWeight")
)
