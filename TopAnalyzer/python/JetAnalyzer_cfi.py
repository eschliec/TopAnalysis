import FWCore.ParameterSet.Config as cms

analyzeJets = cms.EDAnalyzer("JetAnalyzer",
    weight  = cms.InputTag("eventWeight"),
    input   = cms.InputTag("allLayer1Jets"),
    hist    = cms.string('analyzeJets.hist'),

    ## fill id  histograms
    id = cms.PSet(
        nJets = cms.int32(4)
    ),                          
    doId = cms.bool(True),

    ## fill kin histograms
    kin = cms.PSet(
        nJets = cms.int32(4)
    ),                          
    doKin = cms.bool(True),

    ## fill res histograms
    res = cms.PSet(
        binsPt  = cms.vdouble(  0.0,  10.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0),
        binsEta = cms.vdouble( -2.5,  -1.3, -0.5,  0.0,  0.5,   1.3,  2.5),
        binsPhi = cms.vdouble(-3.14, -1.57,  0.0, 1.57, 3.14),
        matchDR = cms.double (0.5),
        nJets   = cms.int32  (4)
    ),
    doRes = cms.bool(True)
)



