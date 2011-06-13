import FWCore.ParameterSet.Config as cms

matrixmodule = cms.EDAnalyzer("MatrixAnalyzer",
    hist_directory = cms.string('.'),
    make = cms.bool(False),
    varBins = cms.vdouble(0.0, 20.0, 40.0, 60.0, 100.0, 
        150.0, 300.0),
    hist = cms.string('QCDMatrix.hist'),
    var = cms.InputTag("selectedLayer1METs"),
    muons = cms.InputTag("allLayer1Muons"),
    pmn = cms.string(''),
    before = cms.bool(False)
)



