import FWCore.ParameterSet.Config as cms

analyzeJets = cms.EDAnalyzer("BJetAnalyzer",
    
    # sources    
    jets = cms.InputTag("selectedLayer1Jets"),    
    
    # histfile        
    hist = cms.string('jets.hist'),
        
    # domains    
    show_kinematics  = cms.bool(True), 
    show_energyscale = cms.bool(True),    
    show_btagging    = cms.bool(True),
     
    # take only leading jets    
    only_leading = cms.int32(999),
            
    # b-tagging discriminator cuts        
    trk_count_high_eff_cut = cms.double(0.3),
    trk_count_high_pur_cut = cms.double(0.3),
    jet_prob_cut           = cms.double(0.3),
    jet_b_prob_cut         = cms.double(0.3),
    impact_par_mva_cut     = cms.double(0.2),
    simple_sv_cut          = cms.double(0.4),
    comb_sv_cut            = cms.double(0.8),
    comb_sv_mva_cut        = cms.double(0.8),
    soft_mu_cut            = cms.double(0.2),
    soft_mu_no_ip_cut      = cms.double(0.2),
    soft_elec_cut          = cms.double(0.2),
        
    # bins for b-tagging
    eta_bins = cms.vdouble(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 
                           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                           2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0
			  ),
    
    pt_bins  = cms.vdouble(15., 20., 25., 30., 35., 40., 50., 60., 70.,
                           85., 100., 120., 150., 200., 250.
			  )
) 

