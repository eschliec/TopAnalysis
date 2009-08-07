import FWCore.ParameterSet.Config as cms
from TopAnalysis.TopFilter.filters.FullLepObjectCountFilter_cfi import *

countMuons         = filterFullLepObjects.clone()
countMuons.name    = cms.string('MuonCounter')
countMuons.objects = cms.InputTag('goodMuons')
countMuons.n       = cms.uint32(2)

countJets         = filterFullLepObjects.clone()
countJets.name    = cms.string('JetCounter')
countJets.objects = cms.InputTag('goodJets')
countJets.n       = cms.uint32(2)

countObjects = cms.Sequence(countMuons *
                            countJets
			   )
