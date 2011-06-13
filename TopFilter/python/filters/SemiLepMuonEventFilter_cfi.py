import FWCore.ParameterSet.Config as cms

## define basic event filter
filterSemiLepMuonEvent = cms.EDFilter("SemiLepMuonEventFilter",
    weight  = cms.InputTag("eventWeight"),
    leptons = cms.VInputTag(cms.InputTag("selectedLayer1Muons")),
    jets    = cms.VInputTag(cms.InputTag("selectedLayer1Jets")),

    ## selection switches                                      
    lepEta = cms.bool(True),
    lepPt  = cms.bool(True),
    jetEta = cms.bool(True),
    jetPt  = cms.bool(True),                                      
    trkIso = cms.bool(True),
    calIso = cms.bool(True),
    jetDist= cms.bool(True),

    ## predefine selection structure
    lepEtaFilter = cms.PSet(
        mode = cms.uint32 ( 0),
        name = cms.string (''),
        min  = cms.vdouble(  ),
        max  = cms.vdouble(  )
    ),
    lepPtFilter = cms.PSet(
        mode = cms.uint32 ( 0),
        name = cms.string (''),
        min  = cms.vdouble(  ),
        max  = cms.vdouble(  )
    ),
    jetEtaFilter = cms.PSet(
        mode = cms.uint32 ( 0),
        name = cms.string (''),
        min  = cms.vdouble(  ),
        max  = cms.vdouble(  )
    ),
    jetPtFilter = cms.PSet(
        mode = cms.uint32 ( 0),
        name = cms.string (''),
        min  = cms.vdouble(  ),
        max  = cms.vdouble(  )
    ),
    trkIsoFilter = cms.PSet(
        type = cms.uint32 ( 0),
        mode = cms.uint32 ( 0),
        name = cms.string (''),
        min  = cms.vdouble(  ),
        max  = cms.vdouble(  )
    ),
    calIsoFilter = cms.PSet(
        type = cms.uint32 ( 0),
        mode = cms.uint32 ( 0),
        name = cms.string (''),
        min  = cms.vdouble(  ),
        max  = cms.vdouble(  )
    ),
    jetDistFilter = cms.PSet(
        refs  = cms.InputTag("selectedLayer1Jets"),
        thresh= cms.vdouble( 0),
        mode  = cms.uint32( 0),
        name  = cms.string (''),
        min   = cms.vdouble(  ),
        max   = cms.vdouble(  )
    )
)


