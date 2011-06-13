import FWCore.ParameterSet.Config as cms

analyzeMuon = cms.EDAnalyzer("MuonAnalyzer",
    weight  = cms.InputTag("eventWeight"),
    input   = cms.InputTag("allLayer1Muons"),
    hist    = cms.string('analyzeMuon.hist'),
                             
    ## fill id  histograms
    id = cms.PSet(
    ),    
    doId = cms.bool(True),

    ## fill kin histograms
    kin = cms.PSet(
        jets = cms.InputTag("allLayer1Jets")
    ),
    doKin = cms.bool(True),

    ## fill res histograms
    res = cms.PSet(
        binsPt  = cms.vdouble( 0.0, 20.0, 30.0, 50.0, 80.0, 120.0, 180.0),
        binsEta = cms.vdouble(-2.5, -1.3, -0.5,  0.0,  0.5,   1.3,   2.5),
        matchDR = cms.double (0.1),
        binsPhi = cms.vdouble(-3.14, -1.57, 0.0, 1.57, 3.14)
    ),                             
    doRes = cms.bool(True)
)



