import FWCore.ParameterSet.Config as cms

## muon preselection
semiLepMuonEta = cms.PSet(
    mode = cms.uint32 ( 0),
    name = cms.string ('semiLepMuonEta'),
    min  = cms.vdouble(-2.1),
    max  = cms.vdouble( 2.1)
)
semiLepMuonPt  = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('semiLepMuonPt' ),
    min  = cms.vdouble(   30.0),
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
    min  = cms.vdouble(   40.0,    40.0,    40.0,    40.0),
    max  = cms.vdouble(99999.0, 99999.0, 99999.0, 99999.0)
)

## isolation preselection
semiLepMuonIsolationTrk = cms.PSet(
    type = cms.uint32 (  1  ),
    mode = cms.uint32 (  0  ),    
    name = cms.string('semiLepMuonIsolationTrk'),
    min  = cms.vdouble(  0.0),
    max  = cms.vdouble(  3.0)
)
semiLepMuonIsolationCal = cms.PSet(
    type = cms.uint32 (  3  ),
    mode = cms.uint32 (  0  ),    
    name = cms.string('semiLepMuonIsolationCal'),
    min  = cms.vdouble(  0.0),
    max  = cms.vdouble(  5.0)
)
semiLepMuonJetDistance = cms.PSet(
    refs  = cms.VInputTag(cms.InputTag("selectedLayer1Jets")),
    thresh= cms.double(  20.),
    mode  = cms.uint32(   0 ),    
    name = cms.string('semiLepMuonIsolationJet'),
    min  = cms.vdouble(  0.0),
    max  = cms.vdouble(  0.5)
)

