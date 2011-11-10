#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFitLepton.h"

/// default constructor for fw lite
HypothesisKinFitLepton::HypothesisKinFitLepton()
{
}

/// default constructor for full fw
HypothesisKinFitLepton::HypothesisKinFitLepton(const edm::ParameterSet& cfg) :
  corrPerm_( cfg.getParameter<bool>   ("corrPerm") ),
  maxChi2_ ( cfg.getParameter<double> ("maxChi2" ) )
{
}

/// histogramm booking for fwlite 
void HypothesisKinFitLepton::book()
{
  /** 
      Pull Distributions (Relative to the Reco Input)
  **/
  // lepton pt
  hists_["leptonPullPtKinFitRec"      ] = new TH1F( "leptonPullPtKinFitRec"       , "leptonPullPtKinFitRec"        ,  200,  -1.,   1. );
  // lepton eta
  hists_["leptonPullEtaKinFitRec"     ] = new TH1F( "leptonPullEtaKinFitRec"      , "leptonPullEtaKinFitRec"       ,  200,  -1.,   1. );
  // lepton phi
  hists_["leptonPullPhiKinFitRec"     ] = new TH1F( "leptonPullPhiKinFitRec"      , "leptonPullPhiKinFitRec"       ,  200,  -1.,   1. );
}

/// histogramm booking for fw
void HypothesisKinFitLepton::book(edm::Service<TFileService>& fs)
{
  /** 
      Pull Distributions (KinFit Relative to the Reco Input)
  **/
  // lepton pt
  hists_["leptonPullPtKinFitRec"      ] = fs->make<TH1F>( "leptonPullPtKinFitRec" , "leptonPullPtKinFitRec"        ,  200,  -1.,   1. );
  // lepton eta
  hists_["leptonPullEtaKinFitRec"     ] = fs->make<TH1F>( "leptonPullEtaKinFitRec", "leptonPullEtaKinFitRec"       ,  200,  -1.,   1. );
  // lepton phi
  hists_["leptonPullPhiKinFitRec"     ] = fs->make<TH1F>( "leptonPullPhiKinFitRec", "leptonPullPhiKinFitRec"       ,  200,  -1.,   1. );
  /** 
      Pull Distributions (Reco Relative to parton truth)
  **/
  // lepton pt
  hists_["leptonPullPtRecPartonTruth"      ] = fs->make<TH1F>( "leptonPullPtRecPartonTruth" , "leptonPullPtRecPartonTruth"        ,  200,  -1.,   1. );
  // lepton eta
  hists_["leptonPullEtaRecPartonTruth"     ] = fs->make<TH1F>( "leptonPullEtaRecPartonTruth", "leptonPullEtaRecPartonTruth"       ,  200,  -1.,   1. );
  // lepton phi
  hists_["leptonPullPhiRecPartonTruth"     ] = fs->make<TH1F>( "leptonPullPhiRecPartonTruth", "leptonPullPhiRecPartonTruth"       ,  200,  -1.,   1. );
  /** 
      Pull Distributions (KinFit Relative to parton truth)
  **/
  // lepton pt
  hists_["leptonPullPtKinFitPartonTruth"      ] = fs->make<TH1F>( "leptonPullPtKinFitPartonTruth" , "leptonPullPtKinFitPartonTruth"        ,  200,  -1.,   1. );
  // lepton eta
  hists_["leptonPullEtaKinFitPartonTruth"     ] = fs->make<TH1F>( "leptonPullEtaKinFitPartonTruth", "leptonPullEtaKinFitPartonTruth"       ,  200,  -1.,   1. );
  // lepton phi
  hists_["leptonPullPhiKinFitPartonTruth"     ] = fs->make<TH1F>( "leptonPullPhiKinFitPartonTruth", "leptonPullPhiKinFitPartonTruth"       ,  200,  -1.,   1. );
}

/// histogram filling interface for reconstruction level for access with fwlite or full framework
void
HypothesisKinFitLepton::fill(const TtSemiLeptonicEvent& tops, const edm::View<reco::Candidate>& leptons, const double& weight)
{
  // make sure to have a valid hypothesis on reconstruction level
  if( tops.isHypoValid("kKinFit") ){
    // a) define index in pat::Lepton collection
    int lepton = tops.jetLeptonCombination("kKinFit")[TtSemiLepEvtPartons::Lepton];
    // b) get requested permutation 
    bool permutation=true;
    if(corrPerm_){
      permutation=false;
      // if jet parton match exists:
      if(tops.isHypoValid("kGenMatch")){
	// indices for all quarks from Kinfit Hypothesis and genmatch
	int lepBIndex         = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LepB     ];
	int hadBIndex         = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::HadB     ];
	int lightQIndex       = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LightQ   ];
	int lightQBarIndex    = tops.jetLeptonCombination("kKinFit"  )[TtSemiLepEvtPartons::LightQBar];
	int lepBIndexGen      = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LepB     ];
	int hadBIndexGen      = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::HadB     ];
	int lightQIndexGen    = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LightQ   ];
	int lightQBarIndexGen = tops.jetLeptonCombination("kGenMatch")[TtSemiLepEvtPartons::LightQBar];
	// check for correct permutation
	if((lepBIndex==lepBIndexGen)&&(hadBIndex==hadBIndexGen)&&
	   (((lightQIndex==lightQIndexGen   )&&(lightQBarIndex==lightQBarIndexGen))||
	    ((lightQIndex==lightQBarIndexGen)&&(lightQBarIndex==lightQIndexGen   )))) permutation=true;
      }
    }
    /** 
	Fill the Pull Distributions
    **/
    // check chi2 and permutation requirement
    if(permutation&&tops.fitChi2()<maxChi2_){
      // make sure the lepton index is in the range of the input collection
      if( lepton < (int)leptons.size() ){
	// c) get values (parton truth, reconstruction, kinematic fit)
	// lepton pt
	double lepPtPartonTruth = tops.singleLepton()->pt();
	double lepPtRec    = leptons[lepton].pt();
	double lepPtKinFit    = tops.singleLepton("kKinFit")->pt();
	// lepton eta
	double lepEtaPartonTruth = tops.singleLepton()->eta();
	double lepEtaRec    = leptons[lepton].eta();
	double lepEtaKinFit    = tops.singleLepton("kKinFit")->eta();
	// lepton phi
	double lepPhiPartonTruth = tops.singleLepton()->phi();
	double lepPhiRec    = leptons[lepton].phi();
	double lepPhiKinFit    = tops.singleLepton("kKinFit")->phi();
	// d) fill plots
	// (i) reconstructed vs kin.fitted
	hists_.find("leptonPullPtKinFitRec"   )->second->Fill( (lepPtKinFit- lepPtRec )/lepPtRec, weight);
	hists_.find("leptonPullEtaKinFitRec"  )->second->Fill( (lepEtaKinFit-lepEtaRec)/lepPtRec, weight);
	hists_.find("leptonPullPhiKinFitRec"  )->second->Fill( (lepPhiKinFit-lepEtaRec)/lepPtRec, weight);
	// (ii) reconstructed vs parton truth
	hists_.find("leptonPullPtRecPartonTruth" )->second->Fill( (lepPtRec -lepPtPartonTruth )/lepPtPartonTruth , weight);
	hists_.find("leptonPullEtaRecPartonTruth")->second->Fill( (lepEtaRec-lepEtaPartonTruth)/lepEtaPartonTruth, weight);
	hists_.find("leptonPullPhiRecPartonTruth")->second->Fill( (lepPhiRec-lepPhiPartonTruth)/lepPhiPartonTruth, weight);
	// (iii) kin.fitted vs parton truth
	hists_.find("leptonPullPtKinFitPartonTruth" )->second->Fill( (lepPtKinFit- lepPtPartonTruth )/lepPtPartonTruth , weight);
	hists_.find("leptonPullEtaKinFitPartonTruth")->second->Fill( (lepEtaKinFit-lepEtaPartonTruth)/lepEtaPartonTruth, weight);
	hists_.find("leptonPullPhiKinFitPartonTruth")->second->Fill( (lepPhiKinFit-lepPhiPartonTruth)/lepPhiPartonTruth, weight);
      }
    }
  }
}
