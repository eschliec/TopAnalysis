import FWCore.ParameterSet.Config as cms

from TopQuarkAnalysis.TopObjectResolutions.stringResolutions_etEtaPhi_cff import *

analyzeFullHadQCDEstimation = cms.EDAnalyzer("FullHadQCDEstimation",
    ## input collection                             
        
    ## src of jets
    JetSrc             = cms.InputTag("selectedPatJets"),
    ## src of multijet mva discriminator
    MultiJetMVADiscSrc = cms.InputTag("mvaDisc","DiscSel"),

    ## weight of mc sample
    MCweight = cms.double(1.0),

    ## special parameters for 
    ## QCD background estimation for
    ## fully hadronic ttbar events

    ## choose TTree for output instead of histograms, if applicable
    useTree    = cms.bool(False),
    ## choose b-tagging algorithm and working point to be used
    ## available are the following so far:
    ## 1. TCHEM (TCHE algo @ medium WP)
    ## 2. TCHPT (TCHP algo @ tight WP)
    bTagAlgoWP = cms.string("TCHEM"),
    # ------------------------------------------------
    # resolutions used for the kinematic fit
    # ------------------------------------------------
    udscResolutions = udscResolution.functions,
    bResolutions    = bjetResolution.functions,
    ## resolution scale factors
    jetEnergyResolutionScaleFactors = cms.vdouble(1.0),
    jetEnergyResolutionEtaBinning = cms.vdouble(0.0,-1.0)                                          
)



