import FWCore.ParameterSet.Config as cms

## muon preselection
diLepMuonEta = cms.PSet(
    mode = cms.uint32 ( 0),
    name = cms.string ('diLepMuonEta'),
    min  = cms.vdouble(-2.5,-2.5),
    max  = cms.vdouble( 2.5, 2.5)
)
diLepMuonPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('diLepMuonPt' ),
    min  = cms.vdouble(  15.0,   15.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## jets preselection
diLepJetsEta = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('diLepJetsEta'),
    min  = cms.vdouble(-3.25, -3.25),
    max  = cms.vdouble( 3.25,  3.25)
)
diLepJetsPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('diLepJetsPt' ),
    min  = cms.vdouble(   5.0,    5.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## isolation preselection
# not needed for preselection since dimuon is seldom enough
