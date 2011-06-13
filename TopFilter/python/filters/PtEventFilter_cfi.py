import FWCore.ParameterSet.Config as cms

## define basic event filters
filterMuonPt = cms.EDFilter("PtEventFilter",
    weight = cms.InputTag("eventWeight"),
    input  = cms.VInputTag(cms.InputTag("selectedLayer1Muons")),
    cuts   = cms.PSet(
      mode = cms.uint32 ( 0),
      name = cms.string (''),
      min  = cms.vdouble(  ),
      max  = cms.vdouble(  )
    )
)
filterJetsPt = cms.EDFilter("PtEventFilter",
    weight = cms.InputTag("eventWeight"),
    input  = cms.VInputTag(cms.InputTag("selectedLayer1Jets")),
    cuts   = cms.PSet(
      mode = cms.uint32 ( 0),
      name = cms.string (''),
      min  = cms.vdouble(  ),
      max  = cms.vdouble(  )
    ) 
)
