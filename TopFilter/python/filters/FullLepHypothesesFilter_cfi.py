import FWCore.ParameterSet.Config as cms

filterFullLepHypothesis = cms.EDAnalyzer("FullLepHypothesesFilter",
    hypoKey = cms.InputTag("ttFullLepHypKinSolution","Key"),
    FullLepEvent = cms.InputTag("ttFullLepEvent"),
    weight = cms.InputTag("eventWeight"),
    jets = cms.InputTag("selectedLayer1Jets"),
        
    ## cuts on muon muon mass
    qcdCut = cms.double(10.),
    zCut   = cms.vdouble(82., 100.),
        
    ## b-tagging (if bDisc has two arguments both jets have 
    ## to pass the cuts, if only one number is given only one jet
    ## has to be tagged.)
    bAlgorithm     = cms.string('simpleSecondaryVertexBJetTags'),
    bDiscriminator = cms.vdouble(1.)
)
