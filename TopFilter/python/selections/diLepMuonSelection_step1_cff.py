import FWCore.ParameterSet.Config as cms

## muon preselection
diLepMuonEta = cms.PSet(
    mode = cms.uint32 ( 0),
    name = cms.string ('diLepMuonEta'),
    min  = cms.vdouble(-2.1,-2.1),
    max  = cms.vdouble( 2.1, 2.1)
)
diLepMuonPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('diLepMuonPt' ),
    min  = cms.vdouble(  20.0,   20.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## jets preselection
diLepJetsEta = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('diLepJetsEta'),
    min  = cms.vdouble(-2.4, -2.4),
    max  = cms.vdouble( 2.4,  2.4)
)
diLepJetsPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('diLepJetsPt' ),
    min  = cms.vdouble(  40.0,   40.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## isolation preselection
diLepCombIsolation = cms.PSet(    
    type = cms.uint32 (  4  ),
    mode = cms.uint32 (  0  ),    
    name = cms.string('diLepMuonCombIso'),
    min  = cms.vdouble(  0.9, 0.9),
    max  = cms.vdouble(  1.0, 1.0)   
)

diLepMuonJetDistance = cms.PSet(
    refs  = cms.InputTag("selectedLayer1Jets"),
    thresh= cms.vdouble(  20.),
    mode  = cms.uint32(   0 ),    
    name = cms.string('semiLepMuonIsoJet'),
    min  = cms.vdouble(  0.0, 0.0),
    max  = cms.vdouble(  0.3, 0.3)
)

