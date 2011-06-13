import FWCore.ParameterSet.Config as cms

## muon preselection
semiLepMuonEta = cms.PSet(
    mode = cms.uint32 ( 0),
    name = cms.string ('semiLepMuonEta'),
    min  = cms.vdouble(-2.5),
    max  = cms.vdouble( 2.5)
)
semiLepMuonPt  = cms.PSet(
    mode = cms.uint32 ( 0),
    name = cms.string ('semiLepMuonPt' ),
    min  = cms.vdouble(   20.0),
    max  = cms.vdouble(99999.0)
)

## jets preselection
semiLepJetsEta = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('semiLepJetsEta'),
    min  = cms.vdouble(-3.0, -3.0, -3.0, -3.0),
    max  = cms.vdouble( 3.0,  3.0,  3.0,  3.0)
)
semiLepJetsPt  = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('semiLepJetsPt' ),
    min  = cms.vdouble(    5.0,     5.0,     5.0,     5.0),
    max  = cms.vdouble(99999.0, 99999.0, 99999.0, 99999.0)
)

## isolation preselection
semiLepMuonIsolationTrk = cms.PSet(
    type = cms.uint32 (  0  ),
    mode = cms.uint32 (  0  ),
    name = cms.string('semiLepMuonIsolationTrk'),
    min  = cms.vdouble(  0.0),
    max  = cms.vdouble(999.0)
)
semiLepMuonIsolationCal = cms.PSet(
    type = cms.uint32 (  0  ),
    mode = cms.uint32 (  0  ),
    name = cms.string('semiLepMuonIsolationCal'),
    min  = cms.vdouble(  0.0),
    max  = cms.vdouble(999.0)
)
semiLepMuonJetDistance = cms.PSet(
    refs  = cms.VInputTag(cms.InputTag("selectedLayer1Jets")),
    thresh= cms.double(  0  ),
    mode  = cms.uint32(  0  ),
    name = cms.string('semiLepMuonJetDistance'),
    min  = cms.vdouble(  0.0),
    max  = cms.vdouble(999.0)
)

