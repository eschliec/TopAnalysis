import FWCore.ParameterSet.Config as cms

analyzeMuonPair = cms.EDAnalyzer("DimuonAnalyzer",

    # Sources
    muons = cms.InputTag("selectedLayer1Muons"),
    ref = cms.InputTag("decaySubset"), 
    
    # Are the events weighted?
    useEventWeight = cms.bool(True), 
    
    # HistFile          
    hist = cms.string('muonPair.hist'), 
    
    # match with mc        
    doMatch = cms.bool(False)
)


