import FWCore.ParameterSet.Config as cms

#
# module to make kinematic distributions for b-jets on  gen and rec level
# gen level b-jet identification done with TopAnalysis/TopUtils/GenLevelBJetProducer
#

analyzeSemiLepBJets = cms.EDAnalyzer("SemiLepBjetAnalyzer",
                                     # get kinematic fit output
                                     semiLepEvent = cms.InputTag("ttSemiLepEvent"),
                                     # specify kinematic fit hypothesis
                                     hypoKey = cms.string("kKinFit"),
                                     # choose genjet collection
                                     # NOTE: as the output of TopAnalysis/TopUtils/GenLevelBJetProducer
                                     #       is used to identify the b-jets, the SAME jet collection
                                     #       has to be chosen
                                     genJets = cms.InputTag('ak5GenJets','','HLT'),
                                     # specify if 'genJets' is a pure b gen jet collection:
                                     # for true, BHadJetIndex and AntiBHadJetIndex will be ignored
                                     # and the leading two jets are chosen!
                                     # NOTE: the assumed ordering within the collection is:
                                     #       1st entry=b-jet, 2nd entry=anti-bjet
                                     bJetCollection = cms.bool(False),
                                     # output manager:
                                     # 0: no output, 1: info, >=2: debug
                                     output = cms.int32(0),
                                     # specify event weight
                                     weight = cms.InputTag(""),
                                     # produce gen level plots?
                                     genPlots = cms.bool(False), 
                                     # produce rec level plots?
                                     recPlots = cms.bool(True),
                                     # input tags to get the (anti)-b-jet indices as created by
                                     # TopAnalysis/TopUtils/GenLevelBJetProducer
                                     # please adapt if cloning and renaming produceGenLevelBJets!
                                     BHadJetIndex     = cms.InputTag("produceGenLevelBJets", "BHadJetIndex"    ),
                                     AntiBHadJetIndex = cms.InputTag("produceGenLevelBJets", "AntiBHadJetIndex"),
                                     # take minimum dR as criteria to assign rec and gen level b-jets to each other?
                                     # this is useful for simultaneous measurements of b-jet and anti-bjet 
                                     # where only the resolution and not the association matters
                                     useClosestDrBs = cms.bool(True),
                                     # create tree?
                                     useTree = cms.bool(False)
                                     )
