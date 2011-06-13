import FWCore.ParameterSet.Config as cms

## muon preselection
fullLepMuonEta = cms.PSet(
    mode = cms.uint32 ( 0),
    name = cms.string ('fullLepMuonEta'),
    min  = cms.vdouble(-2.5,-2.5),
    max  = cms.vdouble( 2.5, 2.5)
)
fullLepMuonPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('fullLepMuonPt' ),
    min  = cms.vdouble(  15.0,   15.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## jets preselection
fullLepJetsEta = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('fullLepJetsEta'),
    min  = cms.vdouble(-3.25, -3.25),
    max  = cms.vdouble( 3.25,  3.25)
)
fullLepJetsPt = cms.PSet(
    mode = cms.uint32 ( 0),    
    name = cms.string ('fullLepJetsPt' ),
    min  = cms.vdouble(   5.0,    5.0),
    max  = cms.vdouble(9999.0, 9999.0)
)

## isolation preselection
# not needed for preselection since dimuon is seldom enough


