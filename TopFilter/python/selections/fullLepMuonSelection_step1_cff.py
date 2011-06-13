import FWCore.ParameterSet.Config as cms

## muon selection
fullLepMuonEta = cms.PSet(
    mode = cms.uint32 (0),
    name = cms.string ('fullLepMuonEta'),
    min  = cms.vdouble(-2.1,-2.1),
    max  = cms.vdouble( 2.1, 2.1)
)
fullLepMuonPt = cms.PSet(
    mode = cms.uint32 (0),    
    name = cms.string ('fullLepMuonPt'),
    min  = cms.vdouble(   20.0,   20.0),
    max  = cms.vdouble( 9999.0, 9999.0)
)

## jets selection
fullLepJetsEta = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('fullLepJetsEta'),
    min  = cms.vdouble(-2.4, -2.4),
    max  = cms.vdouble( 2.4,  2.4)
)
fullLepJetsPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('fullLepJetsPt'),
    min  = cms.vdouble(  40.0,   40.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## isolation selection
fullLepCombIsolation = cms.PSet(
    type = cms.uint32 (  4  ),
    mode = cms.uint32 (  0  ),    
    name = cms.string('fullLepMuonCombIso'),
    min  = cms.vdouble(  0.9, 0.9),
    max  = cms.vdouble(  1.0, 1.0)   
)    
fullLepMuonJetDistance = cms.PSet(
    refs  = cms.InputTag("selectedLayer1Jets"),
    thresh= cms.vdouble(  20.),
    mode  = cms.uint32(   0 ),    
    name = cms.string('semiLepMuonIsoJet'),
    min  = cms.vdouble(  0.0, 0.0),
    max  = cms.vdouble(  0.3, 0.3)
)



