#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFitJets.h"

/// default constructor for fw lite
HypothesisKinFitJets::HypothesisKinFitJets()
{
}

/// default constructor for full fw
HypothesisKinFitJets::HypothesisKinFitJets(const edm::ParameterSet& cfg)
{
}

/// histogramm booking for fwlite 
void HypothesisKinFitJets::book()
{
  /** 
      Pull Distributions (KinFit Relative to the Reco Input)
  **/
  // hadronic top b pt
  hists_["hadBQuarkPullPt"   ] = new TH1F( "hadBQuarkPullPt"    , "hadBQuarkPullPt"     ,  1000,  -5.,   5. );
  // hadronic top b eta
  hists_["hadBQuarkPullEta"  ] = new TH1F( "hadBQuarkPullEta"   , "hadBQuarkPullEta"    ,  1000,  -5.,   5. );
  // hadronic top b phi
  hists_["hadBQuarkPullPhi"  ] = new TH1F( "hadBQuarkPullPhi"   , "hadBQuarkPullPhi"    ,  1000,  -5.,   5. );
  // leptonic top b pt
  hists_["lepBQuarkPullPt"   ] = new TH1F( "lepBQuarkPullPt"    , "lepBQuarkPullPt"     ,  1000,  -5.,   5. );
  // leptonic top b eta
  hists_["lepBQuarkPullEta"  ] = new TH1F( "lepBQuarkPullEta"   , "lepBQuarkPullEta"    ,  1000,  -5.,   5. );
  // leptonic top b phi
  hists_["lepBQuarkPullPhi"  ] = new TH1F( "lepBQuarkPullPhi"   , "lepBQuarkPullPhi"    ,  1000,  -5.,   5. );
  // lightQuarkPull pt
  hists_["lightQuarkPullPt"  ] = new TH1F( "lightQuarkPullPt"   , "lightQuarkPullPt"    ,  1000,  -5.,   5. );
  // lightQuarkPull b eta
  hists_["lightQuarkPullEta" ] = new TH1F( "lightQuarkPullEta"  , "lightQuarkPullEta"   ,  1000,  -5.,   5. );
  // lightQuarkPull b phi
  hists_["lightQuarkPullPhi" ] = new TH1F( "lightQuarkPullPhi"  , "lightQuarkPullPhi"   ,  1000,  -5.,   5. );
}

/// histogramm booking for fw
void HypothesisKinFitJets::book(edm::Service<TFileService>& fs)
{
  /** 
      Pull Distributions (KinFit Relative to the Reco Input)
  **/
  // hadronic top b pt
  hists_["hadBQuarkPullPtKinFitRec"  ] = fs->make<TH1F>( "hadBQuarkPullPtKinFitRec"  , "hadBQuarkPullPtKinFitRec"    ,  1000,  -5.,   5. );
  // hadronic top b eta								   
  hists_["hadBQuarkPullEtaKinFitRec" ] = fs->make<TH1F>( "hadBQuarkPullEtaKinFitRec" , "hadBQuarkPullEtaKinFitRec"   ,  1000,  -5.,   5. );
  // hadronic top b phi								   
  hists_["hadBQuarkPullPhiKinFitRec" ] = fs->make<TH1F>( "hadBQuarkPullPhiKinFitRec" , "hadBQuarkPullPhiKinFitRec"   ,  1000,  -5.,   5. );
  // leptonic top b pt								   
  hists_["lepBQuarkPullPtKinFitRec"  ] = fs->make<TH1F>( "lepBQuarkPullPtKinFitRec"  , "lepBQuarkPullPtKinFitRec"    ,  1000,  -5.,   5. );
  // leptonic top b eta								   
  hists_["lepBQuarkPullEtaKinFitRec" ] = fs->make<TH1F>( "lepBQuarkPullEtaKinFitRec" , "lepBQuarkPullEtaKinFitRec"   ,  1000,  -5.,   5. );
  // leptonic top b phi								   
  hists_["lepBQuarkPullPhiKinFitRec" ] = fs->make<TH1F>( "lepBQuarkPullPhiKinFitRec" , "lepBQuarkPullPhiKinFitRec"   ,  1000,  -5.,   5. );
  // lightQuarkPull pt								   
  hists_["lightQuarkPullPtKinFitRec" ] = fs->make<TH1F>( "lightQuarkPullPtKinFitRec" , "lightQuarkPullPtKinFitRec"   ,  1000,  -5.,   5. );
  // lightQuarkPull b eta							   
  hists_["lightQuarkPullEtaKinFitRec"] = fs->make<TH1F>( "lightQuarkPullEtaKinFitRec", "lightQuarkPullEtaKinFitRec"  ,  1000,  -5.,   5. );
  // lightQuarkPull b phi							   
  hists_["lightQuarkPullPhiKinFitRec"] = fs->make<TH1F>( "lightQuarkPullPhiKinFitRec", "lightQuarkPullPhiKinFitRec"  ,  1000,  -5.,   5. );
  /** 
      Pull Distributions (Reco Relative to parton truth)
  **/
  // hadronic top b pt
  hists_["hadBQuarkPullPtRecPartonTruth"  ] = fs->make<TH1F>( "hadBQuarkPullPtRecPartonTruth"  , "hadBQuarkPullPtRecPartonTruth"    ,  1000,  -5.,   5. );
  // hadronic top b eta								   
  hists_["hadBQuarkPullEtaRecPartonTruth" ] = fs->make<TH1F>( "hadBQuarkPullEtaRecPartonTruth" , "hadBQuarkPullEtaRecPartonTruth"   ,  1000,  -5.,   5. );
  // hadronic top b phi								   
  hists_["hadBQuarkPullPhiRecPartonTruth" ] = fs->make<TH1F>( "hadBQuarkPullPhiRecPartonTruth" , "hadBQuarkPullPhiRecPartonTruth"   ,  1000,  -5.,   5. );
  // leptonic top b pt								   
  hists_["lepBQuarkPullPtRecPartonTruth"  ] = fs->make<TH1F>( "lepBQuarkPullPtRecPartonTruth"  , "lepBQuarkPullPtRecPartonTruth"    ,  1000,  -5.,   5. );
  // leptonic top b eta								   
  hists_["lepBQuarkPullEtaRecPartonTruth" ] = fs->make<TH1F>( "lepBQuarkPullEtaRecPartonTruth" , "lepBQuarkPullEtaRecPartonTruth"   ,  1000,  -5.,   5. );
  // leptonic top b phi								   
  hists_["lepBQuarkPullPhiRecPartonTruth" ] = fs->make<TH1F>( "lepBQuarkPullPhiRecPartonTruth" , "lepBQuarkPullPhiRecPartonTruth"   ,  1000,  -5.,   5. );
  // lightQuarkPull pt								   
  hists_["lightQuarkPullPtRecPartonTruth" ] = fs->make<TH1F>( "lightQuarkPullPtRecPartonTruth" , "lightQuarkPullPtRecPartonTruth"   ,  1000,  -5.,   5. );
  // lightQuarkPull b eta							   
  hists_["lightQuarkPullEtaRecPartonTruth"] = fs->make<TH1F>( "lightQuarkPullEtaRecPartonTruth", "lightQuarkPullEtaRecPartonTruth"  ,  1000,  -5.,   5. );
  // lightQuarkPull b phi							   
  hists_["lightQuarkPullPhiRecPartonTruth"] = fs->make<TH1F>( "lightQuarkPullPhiRecPartonTruth", "lightQuarkPullPhiRecPartonTruth"  ,  1000,  -5.,   5. );

  /** 
      Pull Distributions (KinFit Relative to parton truth)
  **/
  // hadronic top b pt
  hists_["hadBQuarkPullPtKinFitPartonTruth"  ] = fs->make<TH1F>( "hadBQuarkPullPtKinFitPartonTruth"  , "hadBQuarkPullPtKinFitPartonTruth"    ,  1000,  -5.,   5. );
  // hadronic top b eta								   
  hists_["hadBQuarkPullEtaKinFitPartonTruth" ] = fs->make<TH1F>( "hadBQuarkPullEtaKinFitPartonTruth" , "hadBQuarkPullEtaKinFitPartonTruth"   ,  1000,  -5.,   5. );
  // hadronic top b phi								   
  hists_["hadBQuarkPullPhiKinFitPartonTruth" ] = fs->make<TH1F>( "hadBQuarkPullPhiKinFitPartonTruth" , "hadBQuarkPullPhiKinFitPartonTruth"   ,  1000,  -5.,   5. );
  // leptonic top b pt								   
  hists_["lepBQuarkPullPtKinFitPartonTruth"  ] = fs->make<TH1F>( "lepBQuarkPullPtKinFitPartonTruth"  , "lepBQuarkPullPtKinFitPartonTruth"    ,  1000,  -5.,   5. );
  // leptonic top b eta								   
  hists_["lepBQuarkPullEtaKinFitPartonTruth" ] = fs->make<TH1F>( "lepBQuarkPullEtaKinFitPartonTruth" , "lepBQuarkPullEtaKinFitPartonTruth"   ,  1000,  -5.,   5. );
  // leptonic top b phi								   
  hists_["lepBQuarkPullPhiKinFitPartonTruth" ] = fs->make<TH1F>( "lepBQuarkPullPhiKinFitPartonTruth" , "lepBQuarkPullPhiKinFitPartonTruth"   ,  1000,  -5.,   5. );
  // lightQuarkPull pt								   
  hists_["lightQuarkPullPtKinFitPartonTruth" ] = fs->make<TH1F>( "lightQuarkPullPtKinFitPartonTruth" , "lightQuarkPullPtKinFitPartonTruth"   ,  1000,  -5.,   5. );
  // lightQuarkPull b eta							   
  hists_["lightQuarkPullEtaKinFitPartonTruth"] = fs->make<TH1F>( "lightQuarkPullEtaKinFitPartonTruth", "lightQuarkPullEtaKinFitPartonTruth"  ,  1000,  -5.,   5. );
  // lightQuarkPull b phi							   
  hists_["lightQuarkPullPhiKinFitPartonTruth"] = fs->make<TH1F>( "lightQuarkPullPhiKinFitPartonTruth", "lightQuarkPullPhiKinFitPartonTruth"  ,  1000,  -5.,   5. );
}

/// histogram filling interface for reconstruction level for access with fwlite or full framework
void
HypothesisKinFitJets::fill(const TtSemiLeptonicEvent& tops, const edm::View<pat::Jet>& jets, const double& weight)
{
  // make sure to have a valid hypothesis on reconstruction level
  if( tops.isHypoValid("kKinFit") ){
    // define indices in pat::Jet collection
    int lightQ    = tops.jetLeptonCombination("kKinFit")[TtSemiLepEvtPartons::LightQ   ];
    int lightQBar = tops.jetLeptonCombination("kKinFit")[TtSemiLepEvtPartons::LightQBar];
    int hadB      = tops.jetLeptonCombination("kKinFit")[TtSemiLepEvtPartons::HadB     ];
    int lepB      = tops.jetLeptonCombination("kKinFit")[TtSemiLepEvtPartons::LepB     ];
    /** 
	Fill the Pull Distributions
    **/    
    // a) hadronic b-jet
    // make sure the hadronic b-jet index is in the range of the jet collection
    if( hadB < (int)jets.size() ){
      // get values (parton truth, reconstruction, kinematic fit)
      double hadBPtPartonTruth = tops.hadronicDecayB()->pt();
      double hadBPtRec    = jets[hadB].pt();
      double hadBPtKinFit    = tops.hadronicDecayB("kKinFit")->pt();
      double hadBEtaPartonTruth = tops.hadronicDecayB()->eta();
      double hadBEtaRec    = jets[hadB].eta();
      double hadBEtaKinFit    = tops.hadronicDecayB("kKinFit")->eta();
      double hadBPhiPartonTruth = tops.hadronicDecayB()->phi();
      double hadBPhiRec    = jets[hadB].phi();
      double hadBPhiKinFit    = tops.hadronicDecayB("kKinFit")->phi();
      // fill plots
      // (i) reconstructed vs kin.fitted
      // hadronic top b pt
      hists_.find("hadBQuarkPullPtKinFitRec"   )->second->Fill( (hadBPtKinFit -hadBPtRec )/hadBPtRec , weight);
      // hadronic top b eta					
      hists_.find("hadBQuarkPullEtaKinFitRec"  )->second->Fill( (hadBEtaKinFit-hadBEtaRec)/hadBEtaRec, weight);
      // hadronic top b phi					
      hists_.find("hadBQuarkPullPhiKinFitRec"  )->second->Fill( (hadBPhiKinFit-hadBPhiRec)/hadBPhiRec, weight);
      // (ii) reconstructed vs parton truth
      // hadronic top b pt
      hists_.find("hadBQuarkPullPtRecPartonTruth"   )->second->Fill( (hadBPtRec -hadBPtPartonTruth )/hadBPtPartonTruth , weight);
      // hadronic top b eta					
      hists_.find("hadBQuarkPullEtaRecPartonTruth"  )->second->Fill( (hadBEtaRec-hadBEtaPartonTruth)/hadBEtaPartonTruth, weight);
      // hadronic top b phi					
      hists_.find("hadBQuarkPullPhiRecPartonTruth"  )->second->Fill( (hadBPhiRec-hadBPhiPartonTruth)/hadBPhiPartonTruth, weight);
      // (iii) kin.fitted vs parton truth
      // hadronic top b pt
      hists_.find("hadBQuarkPullPtKinFitPartonTruth"   )->second->Fill( (hadBPtKinFit -hadBPtPartonTruth )/hadBPtPartonTruth , weight);
      // hadronic top b eta					
      hists_.find("hadBQuarkPullEtaKinFitPartonTruth"  )->second->Fill( (hadBEtaKinFit-hadBEtaPartonTruth)/hadBEtaPartonTruth, weight);
      // hadronic top b phi					
      hists_.find("hadBQuarkPullPhiKinFitPartonTruth"  )->second->Fill( (hadBPhiKinFit-hadBPhiPartonTruth)/hadBPhiPartonTruth, weight);
    }

    // b) leptonic b-jet
    // make sure the leptonic b-jet index is in the range of the jet collection
    if( lepB < (int)jets.size() ){
      // get values (parton truth, reconstruction, kinematic fit)
      double lepBPtPartonTruth = tops.leptonicDecayB()->pt();
      double lepBPtRec    = jets[lepB].pt();
      double lepBPtKinFit    = tops.leptonicDecayB("kKinFit")->pt();
      double lepBEtaPartonTruth = tops.leptonicDecayB()->eta();
      double lepBEtaRec    = jets[lepB].eta();
      double lepBEtaKinFit    = tops.leptonicDecayB("kKinFit")->eta();
      double lepBPhiPartonTruth = tops.leptonicDecayB()->phi();
      double lepBPhiRec    = jets[lepB].phi();
      double lepBPhiKinFit    = tops.leptonicDecayB("kKinFit")->phi();
      // fill plots
      // (i) reconstructed vs kin.fitted
      // leptonic top b pt
      hists_.find("lepBQuarkPullPtKinFitRec"   )->second->Fill( (lepBPtKinFit -lepBPtRec )/lepBPtRec , weight);
      // leptonic top b eta					
      hists_.find("lepBQuarkPullEtaKinFitRec"  )->second->Fill( (lepBEtaKinFit-lepBEtaRec)/lepBEtaRec, weight);
      // leptonic top b phi					
      hists_.find("lepBQuarkPullPhiKinFitRec"  )->second->Fill( (lepBPhiKinFit-lepBPhiRec)/lepBPhiRec, weight);
      // (ii) reconstructed vs parton truth
      // leptonic top b pt
      hists_.find("lepBQuarkPullPtRecPartonTruth"   )->second->Fill( (lepBPtRec -lepBPtPartonTruth )/lepBPtPartonTruth , weight);
      // leptonic top b eta					
      hists_.find("lepBQuarkPullEtaRecPartonTruth"  )->second->Fill( (lepBEtaRec-lepBEtaPartonTruth)/lepBEtaPartonTruth, weight);
      // leptonic top b phi					
      hists_.find("lepBQuarkPullPhiRecPartonTruth"  )->second->Fill( (lepBPhiRec-lepBPhiPartonTruth)/lepBPhiPartonTruth, weight);
      // (iii) kin.fitted vs parton truth
      // leptonic top b pt
      hists_.find("lepBQuarkPullPtKinFitPartonTruth"   )->second->Fill( (lepBPtKinFit -lepBPtPartonTruth )/lepBPtPartonTruth , weight);
      // leptonic top b eta					
      hists_.find("lepBQuarkPullEtaKinFitPartonTruth"  )->second->Fill( (lepBEtaKinFit-lepBEtaPartonTruth)/lepBEtaPartonTruth, weight);
      // leptonic top b phi					
      hists_.find("lepBQuarkPullPhiKinFitPartonTruth"  )->second->Fill( (lepBPhiKinFit-lepBPhiPartonTruth)/lepBPhiPartonTruth, weight);
    }

   // make sure the light quark indices is in the range of the jet collection
    if( lightQ < (int)jets.size() && lightQBar < (int)jets.size() ){
      // get values (parton truth, reconstruction, kinematic fit)
      double lightQBarPtPartonTruth = tops.hadronicDecayQuarkBar()->pt();
      double lightQBarPtRec    = jets[lightQBar].pt();
      double lightQBarPtKinFit    = tops.hadronicDecayQuarkBar("kKinFit")->pt();
      double lightQBarEtaPartonTruth = tops.hadronicDecayQuarkBar()->eta();
      double lightQBarEtaRec    = jets[lightQBar].eta();
      double lightQBarEtaKinFit    = tops.hadronicDecayQuarkBar("kKinFit")->eta();
      double lightQBarPhiPartonTruth = tops.hadronicDecayQuarkBar()->phi();
      double lightQBarPhiRec    = jets[lightQBar].phi();
      double lightQBarPhiKinFit    = tops.hadronicDecayQuarkBar("kKinFit")->phi();
      double lightQPtPartonTruth = tops.hadronicDecayQuark()->pt();
      double lightQPtRec    = jets[lightQ].pt();
      double lightQPtKinFit    = tops.hadronicDecayQuark("kKinFit")->pt();
      double lightQEtaPartonTruth = tops.hadronicDecayQuark()->eta();
      double lightQEtaRec    = jets[lightQ].eta();
      double lightQEtaKinFit    = tops.hadronicDecayQuark("kKinFit")->eta();
      double lightQPhiPartonTruth = tops.hadronicDecayQuark()->phi();
      double lightQPhiRec    = jets[lightQ].phi();
      double lightQPhiKinFit    = tops.hadronicDecayQuark("kKinFit")->phi();
      // choose quark - antiquark assignment that fits best to truth values
      double pullQQbar=std::abs((lightQPtKinFit -lightQPtPartonTruth )/lightQPtPartonTruth );
      pullQQbar +=     std::abs((lightQEtaKinFit-lightQEtaPartonTruth)/lightQEtaPartonTruth);
      pullQQbar +=     std::abs((lightQPhiKinFit-lightQPhiPartonTruth)/lightQPhiPartonTruth);
      pullQQbar +=     std::abs((lightQBarPtKinFit -lightQBarPtPartonTruth )/lightQBarPtPartonTruth );
      pullQQbar +=     std::abs((lightQBarEtaKinFit-lightQBarEtaPartonTruth)/lightQBarEtaPartonTruth);
      pullQQbar +=     std::abs((lightQBarPhiKinFit-lightQBarPhiPartonTruth)/lightQBarPhiPartonTruth);
      double pullQQbarSwitched=std::abs((lightQPtKinFit -lightQBarPtPartonTruth )/lightQBarPtPartonTruth );
      pullQQbarSwitched +=     std::abs((lightQEtaKinFit-lightQBarEtaPartonTruth)/lightQBarEtaPartonTruth);
      pullQQbarSwitched +=     std::abs((lightQPhiKinFit-lightQBarPhiPartonTruth)/lightQBarPhiPartonTruth);
      pullQQbarSwitched +=     std::abs((lightQBarPtKinFit -lightQPtPartonTruth )/lightQPtPartonTruth );
      pullQQbarSwitched +=     std::abs((lightQBarEtaKinFit-lightQEtaPartonTruth)/lightQEtaPartonTruth);
      pullQQbarSwitched +=     std::abs((lightQBarPhiKinFit-lightQPhiPartonTruth)/lightQPhiPartonTruth);
      if(pullQQbarSwitched<pullQQbar){
	lightQBarPtPartonTruth  = tops.hadronicDecayQuark()->pt();
	lightQBarEtaPartonTruth = tops.hadronicDecayQuark()->eta();
	lightQBarPhiPartonTruth = tops.hadronicDecayQuark()->phi();
	lightQPtPartonTruth  = tops.hadronicDecayQuarkBar()->pt();
	lightQEtaPartonTruth = tops.hadronicDecayQuarkBar()->eta();
	lightQPhiPartonTruth = tops.hadronicDecayQuarkBar()->phi();
      }
      // fill plots
      // (i) reconstructed vs kin.fitted
      // light q pt
      hists_.find("lightQuarkPullPtKinFitRec" )->second->Fill( (lightQBarPtKinFit-lightQBarPtRec )/lightQBarPtRec , weight);
      hists_.find("lightQuarkPullPtKinFitRec" )->second->Fill( (lightQPtKinFit   -lightQPtRec    )/lightQPtRec    , weight);
      // light q eta			      	
      hists_.find("lightQuarkPullEtaKinFitRec")->second->Fill( (lightQBarEtaKinFit-lightQBarEtaRec)/lightQBarEtaRec, weight);
      hists_.find("lightQuarkPullEtaKinFitRec")->second->Fill( (lightQEtaKinFit   -lightQEtaRec  )/lightQEtaRec   , weight);
      // light q phi			      	
      hists_.find("lightQuarkPullPhiKinFitRec")->second->Fill( (lightQBarPhiKinFit-lightQBarPhiRec)/lightQBarPhiRec, weight);
      hists_.find("lightQuarkPullPhiKinFitRec")->second->Fill( (lightQPhiKinFit   -lightQPhiRec   )/lightQPhiRec   , weight);
      // (ii) reconstructed vs parton truth
      // light q pt
      hists_.find("lightQuarkPullPtRecPartonTruth" )->second->Fill( (lightQBarPtRec-lightQBarPtPartonTruth)/lightQBarPtPartonTruth, weight);
      hists_.find("lightQuarkPullPtRecPartonTruth" )->second->Fill( (lightQPtRec   -lightQPtPartonTruth   )/lightQPtPartonTruth   , weight);
      // light q eta				   	
      hists_.find("lightQuarkPullEtaRecPartonTruth")->second->Fill( (lightQBarEtaRec-lightQBarEtaPartonTruth)/lightQBarEtaPartonTruth, weight);
      hists_.find("lightQuarkPullEtaRecPartonTruth")->second->Fill( (lightQEtaRec   -lightQEtaPartonTruth   )/lightQEtaPartonTruth   , weight);
      // light q phi				   	
      hists_.find("lightQuarkPullPhiRecPartonTruth")->second->Fill( (lightQBarPhiRec-lightQBarPhiPartonTruth)/lightQBarPhiPartonTruth, weight);
      hists_.find("lightQuarkPullPhiRecPartonTruth")->second->Fill( (lightQPhiRec   -lightQPhiPartonTruth   )/lightQPhiPartonTruth   , weight);
      // (iii) kin.fitted vs parton truth
      // light q pt
      hists_.find("lightQuarkPullPtKinFitPartonTruth" )->second->Fill( (lightQBarPtKinFit-lightQBarPtPartonTruth )/lightQBarPtPartonTruth, weight);
      hists_.find("lightQuarkPullPtKinFitPartonTruth" )->second->Fill( (lightQPtKinFit   -lightQPtPartonTruth    )/lightQPtPartonTruth   , weight);
      // light q eta				      
      hists_.find("lightQuarkPullEtaKinFitPartonTruth")->second->Fill( (lightQBarEtaKinFit-lightQBarEtaPartonTruth)/lightQBarEtaPartonTruth, weight);
      hists_.find("lightQuarkPullEtaKinFitPartonTruth")->second->Fill( (lightQEtaKinFit   -lightQEtaPartonTruth   )/lightQEtaPartonTruth   , weight);
      // light q phi				      
      hists_.find("lightQuarkPullPhiKinFitPartonTruth")->second->Fill( (lightQBarPhiKinFit-lightQBarPhiPartonTruth)/lightQBarPhiPartonTruth, weight);
      hists_.find("lightQuarkPullPhiKinFitPartonTruth")->second->Fill( (lightQPhiKinFit   -lightQPhiPartonTruth   )/lightQPhiPartonTruth   , weight);
    }
  }
}
