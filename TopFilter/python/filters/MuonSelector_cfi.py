import FWCore.ParameterSet.Config as cms

# module to make quality and isolation cuts on muon tracks see:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string
# for comments about the cuts see: 
#  https://twiki.cern.ch/twiki/bin/view/CMS/VplusJets

goodMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedLayer1Muons"),
    cut = cms.string('combinedMuon.isNull = 0'
		     '& pt > 20.'
		     '& abs(eta) < 2.1'    
                     '& track.numberOfValidHits >= 11' 
		     '& abs(track.d0) < 0.2'
		     '& combinedMuon.normalizedChi2 < 10.0'
		     '& ecalIsoDeposit.candEnergy < 4'
		     '& hcalIsoDeposit.candEnergy < 6'
		     '& (trackIso+caloIso)/pt < 0.1' 
		    )
)
