import FWCore.ParameterSet.Config as cms

# module to make cuts on jets
# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string
#
# these cuts assuere that the jets are in a kinematic range to make
# b-tagging usable

goodJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("selectedLayer1Jets"),
    cut = cms.string('pt > 40.'
		     '& abs(eta) < 2.3'
		    )
)
