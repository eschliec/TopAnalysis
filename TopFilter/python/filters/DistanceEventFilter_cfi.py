import FWCore.ParameterSet.Config as cms

## define basic event filters
filterMuonDistance = cms.EDFilter("MuonDistanceEventFilter",
    weight = cms.InputTag("eventWeight"),
    input  = cms.VInputTag(cms.InputTag("selectedLayer1Muons")),
    cuts    = cms.PSet(
      refs  = cms.VInputTag(cms.InputTag("selectedLayer1Jets")),
      thresh= cms.uint32 ( 0),
      mode  = cms.uint32 ( 0),
      name  = cms.string (''),
      min   = cms.vdouble(  ),
      max   = cms.vdouble(  )
    )    
)
