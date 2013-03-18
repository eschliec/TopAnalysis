import FWCore.ParameterSet.Config as cms

writeFullHadTree = cms.EDAnalyzer("FullHadTreeWriter",
                                  ## jets
                                  JetSrc = cms.InputTag('selectedPatJets'),
                                  ## MET
                                  METSrc = cms.InputTag('patMETs'),
                                  ## muons
                                  MuonSrc = cms.InputTag('selectedPatMuons'),
                                  ## electrons
                                  ElectronSrc = cms.InputTag('selectedPatElectrons'),
                                  ## kinematic fit result
                                  FitSrc = cms.InputTag("ttFullHadEvent"),
                                  ## GenEventInfoProduct for PDF uncertainties
                                  GenSrc = cms.InputTag("generator"),
                                  ## multi-jet MVA discriminator
                                  MultiJetMVADiscSrc = cms.InputTag("mvaDisc","DiscSel"),
                                  ## PU info src for MC
                                  PUSrc = cms.InputTag("addPileupInfo"),
                                  ## primary vertex info
                                  VertexSrc = cms.InputTag("goodOfflinePrimaryVertices"),
                                  ## weight from PU reweighting
                                  PUweightSrc = cms.InputTag("eventWeightPU","eventWeightPU"),
                                  ## weight of MC sample
                                  MCweight = cms.double(-1.),
                                  ## do the PDF uncertainty
                                  DoPDFUncertainty = cms.bool(False),
                                  ## information on b-tagging
                                  bTagName = cms.vstring("BTAGTCHEM"   ,"BTAGTCHEM"   ,
                                                         "MISTAGTCHEM" ,"MISTAGTCHEM" ,
                                                         "BTAGTCHPT"   ,"BTAGTCHPT"   ,
                                                         "MISTAGTCHPT" ,"MISTAGTCHPT" ,
                                                         "BTAGSSVHEM"  ,"BTAGSSVHEM"  ,
                                                         "MISTAGSSVHEM","MISTAGSSVHEM",
                                                         "BTAGSSVHPT"  ,"BTAGSSVHPT"  ,
                                                         "MISTAGSSVHPT","MISTAGSSVHPT",
                                                         "BTAGCSVT"    ,"BTAGCSVT"    ,
                                                         "MISTAGCSVT"  ,"MISTAGCSVT"  ),
                                  bTagVal  = cms.vstring("BTAGBEFFCORR","BTAGBERRCORR",
                                                         "BTAGLEFFCORR","BTAGLERRCORR",
                                                         "BTAGBEFFCORR","BTAGBERRCORR",
                                                         "BTAGLEFFCORR","BTAGLERRCORR",
                                                         "BTAGBEFFCORR","BTAGBERRCORR",
                                                         "BTAGLEFFCORR","BTAGLERRCORR",
                                                         "BTAGBEFFCORR","BTAGBERRCORR",
                                                         "BTAGLEFFCORR","BTAGLERRCORR",
                                                         "BTAGBEFFCORR","BTAGBERRCORR",
                                                         "BTAGLEFFCORR","BTAGLERRCORR"),
                                  )



