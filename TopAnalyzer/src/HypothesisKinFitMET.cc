#include "TopAnalysis/TopAnalyzer/interface/HypothesisKinFitMET.h"

/// default constructor for fw lite
HypothesisKinFitMET::HypothesisKinFitMET()
{
}

/// default constructor for full fw
HypothesisKinFitMET::HypothesisKinFitMET(const edm::ParameterSet& cfg)
{
}

/// histogramm booking for fwlite 
void HypothesisKinFitMET::book()
{
  /** 
      Pull Distributions (KinFit Relative to the Reco Input)
  **/
  // neutrino pt
  hists_["neutrinoPullPtKinFitRec" ] = new TH1F( "neutrinoPullPtKinFitRec" , "neutrinoPullPtKinFitRec" , 1000,  -5.,   5. );
  // neutrino eta	       					   			   
  hists_["neutrinoPullEtaKinFitRec"] = new TH1F( "neutrinoPullEtaKinFitRec", "neutrinoPullEtaKinFitRec", 1000,  -5.,   5. );
  // neutrino phi	       					   			   
  hists_["neutrinoPullPhiKinFitRec"] = new TH1F( "neutrinoPullPhiKinFitRec", "neutrinoPullPhiKinFitRec", 1000,  -5.,   5. );
}

/// histogramm booking for fw
void HypothesisKinFitMET::book(edm::Service<TFileService>& fs)
{
  /** 
     Pull Distributions (KinFit Relative to the Reco Input)
  **/
  // neutrino pt
  hists_["neutrinoPullPtKinFitRec" ] = fs->make<TH1F>( "neutrinoPullPtKinFitRec" , "neutrinoPullPtKinFitRec" , 1000, -5., 5. );
  // neutrino eta	       								      
  hists_["neutrinoPullEtaKinFitRec"] = fs->make<TH1F>( "neutrinoPullEtaKinFitRec", "neutrinoPullEtaKinFitRec", 1000, -5., 5. );
  // neutrino phi	       								      
  hists_["neutrinoPullPhiKinFitRec"] = fs->make<TH1F>( "neutrinoPullPhiKinFitRec", "neutrinoPullPhiKinFitRec", 1000, -5., 5. );
  /** 
      Pull Distributions (Reco Relative to parton truth)
  **/
  hists_["neutrinoPullPtRecPartonTruth" ] = fs->make<TH1F>( "neutrinoPullPtRecPartonTruth" , "neutrinoPullPtRecPartonTruth" , 1000, -5., 5. );
  // neutrino eta			  
  hists_["neutrinoPullEtaRecPartonTruth"] = fs->make<TH1F>( "neutrinoPullEtaRecPartonTruth", "neutrinoPullEtaRecPartonTruth", 1000, -5., 5. );
  // neutrino phi			  
  hists_["neutrinoPullPhiRecPartonTruth"] = fs->make<TH1F>( "neutrinoPullPhiRecPartonTruth", "neutrinoPullPhiRecPartonTruth", 1000, -5., 5. );
  /** 
      Pull Distributions (KinFit Relative to parton truth)
  **/
  hists_["neutrinoPullPtKinFitPartonTruth" ] = fs->make<TH1F>( "neutrinoPullPtKinFitPartonTruth" , "neutrinoPullPtKinFitPartonTruth" , 1000, -5., 5. );
  // neutrino eta			  
  hists_["neutrinoPullEtaKinFitPartonTruth"] = fs->make<TH1F>( "neutrinoPullEtaKinFitPartonTruth", "neutrinoPullEtaKinFitPartonTruth", 1000, -5., 5. );
  // neutrino phi			  
  hists_["neutrinoPullPhiKinFitPartonTruth"] = fs->make<TH1F>( "neutrinoPullPhiKinFitPartonTruth", "neutrinoPullPhiKinFitPartonTruth", 1000, -5., 5. );
}

/// histogram filling interface for reconstruction level for access with fwlite or full framework
void
HypothesisKinFitMET::fill(const TtSemiLeptonicEvent& tops, const edm::View<pat::MET>& mets, const double& weight)
{
  // make sure to have a valid hypothesis on reconstruction level
  if( tops.isHypoValid("kKinFit") ){
    /** 
	Fill the Pull Distributions (Relative to the Reco Input)
    **/
    // make sure the met index is in the range of the met collection
    if( 0 < (int)mets.size() ){
      // a) get values (parton truth, reconstruction, kinematic fit)
      // neutrino pt
      double neutrinoPtPartonTruth = tops.singleNeutrino()->pt();
      double neutrinoPtRec    = mets[0].pt();
      double neutrinoPtKinFit    = tops.singleNeutrino("kKinFit")->pt();
      // neutrino eta
      double neutrinoEtaPartonTruth = tops.singleNeutrino()->eta();
      double neutrinoEtaRec    = mets[0].eta();
      double neutrinoEtaKinFit    = tops.singleNeutrino("kKinFit")->eta();
      // neutrino phi
      double neutrinoPhiPartonTruth = tops.singleNeutrino()->phi();
      double neutrinoPhiRec    = mets[0].phi();
      double neutrinoPhiKinFit    = tops.singleNeutrino("kKinFit")->phi();
      // b) fill plots
      // (i) reconstructed vs kin.fitted
      hists_.find("neutrinoPullPtKinFitRec"   )->second->Fill( (neutrinoPtKinFit- neutrinoPtRec )/neutrinoPtRec, weight);
      hists_.find("neutrinoPullEtaKinFitRec"  )->second->Fill( (neutrinoEtaKinFit-neutrinoEtaRec)/neutrinoPtRec, weight);
      hists_.find("neutrinoPullPhiKinFitRec"  )->second->Fill( (neutrinoPhiKinFit-neutrinoEtaRec)/neutrinoPtRec, weight);
      // (ii) reconstructed vs parton truth
      hists_.find("neutrinoPullPtRecPartonTruth" )->second->Fill( (neutrinoPtRec -neutrinoPtPartonTruth )/neutrinoPtPartonTruth , weight);
      hists_.find("neutrinoPullEtaRecPartonTruth")->second->Fill( (neutrinoEtaRec-neutrinoEtaPartonTruth)/neutrinoEtaPartonTruth, weight);
      hists_.find("neutrinoPullPhiRecPartonTruth")->second->Fill( (neutrinoPhiRec-neutrinoPhiPartonTruth)/neutrinoPhiPartonTruth, weight);
      // (iii) kin.fitted vs parton truth
      hists_.find("neutrinoPullPtKinFitPartonTruth" )->second->Fill( (neutrinoPtKinFit- neutrinoPtPartonTruth )/neutrinoPtPartonTruth , weight);
      hists_.find("neutrinoPullEtaKinFitPartonTruth")->second->Fill( (neutrinoEtaKinFit-neutrinoEtaPartonTruth)/neutrinoEtaPartonTruth, weight);
      hists_.find("neutrinoPullPhiKinFitPartonTruth")->second->Fill( (neutrinoPhiKinFit-neutrinoPhiPartonTruth)/neutrinoPhiPartonTruth, weight);
    }
  }
} 
