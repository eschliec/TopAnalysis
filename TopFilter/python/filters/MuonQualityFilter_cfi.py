import FWCore.ParameterSet.Config as cms

# module to make quality cuts on muon tracks
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string
#
# these are the cuts suggested at the lepton+jets workshop 09/18/2008
# see http://indico.cern.ch/getFile.py/access?contribId=7&sessionId=
# 0&resId=1&materialId=slides&confId=40028

qualityMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("selectedLayer1Muons"),
    cut = cms.string('combinedMuon.isNull = 0'
                     '& track.numberOfValidHits >= 7' 
		     '& abs(track.d0) < 0.25'
		     '& combinedMuon.normalizedChi2 < 5.0'
		     #'& pt/(trackIso+caloIso+pt) > 0.9' 
		    )
)
