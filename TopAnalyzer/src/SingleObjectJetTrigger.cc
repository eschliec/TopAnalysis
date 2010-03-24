#include "TopAnalysis/TopAnalyzer/interface/SingleObjectJetTrigger.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

#include "TMath.h"

//#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

/// default constructor
SingleObjectJetTrigger::SingleObjectJetTrigger(const edm::ParameterSet& cfg)
{
  patTriggerEvent_   = cfg.getParameter<edm::InputTag>("patTriggerEvent");
  patTrigger_        = cfg.getParameter<edm::InputTag>("patTrigger");
  jets_              = cfg.getParameter<edm::InputTag>("jets");
  triggerMatchedJets_= cfg.getParameter<edm::InputTag>("triggerMatchedJets");
}

/// calculate the trigger efficiency for QuadJet (or FiveJet) triggers based on single object trigger efficiency
double 
SingleObjectJetTrigger::triggerEfficiency(const edm::Handle<pat::JetCollection> jets, const bool fiveJet)
{
  double effy_ = 1.; //   X0,       X1,       X2,       X3,       X4,       X5,       X6,       X7,       X8,       X9,
  double effies_[100] = { 0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        // 0X
			  0,        0,        0,        0,        0,        0,        0,        0,        0,        0,	      // 1X
			  0.103389, 0.13367,  0.158888, 0.176171, 0.181272, 0.19171,  0.194168, 0.206794, 0.22368,  0.269992, // 2X
			  0.609654, 0.768382, 0.837768, 0.903173, 0.957051, 0.980878, 0.985011, 0.988416, 0.989534, 0.991462, // 3X
			  0.992501, 0.994188, 0.994625, 0.995198, 0.995896, 0.99635,  0.997574, 0.997853, 0.99841,  0.998711, // 4X
			  0.99913,  0.999176, 0.999175, 0.999411, 0.999525, 0.999681, 0.99987,  0.999767, 0.999799, 0.999829, // 5X
			  0.999792, 0.99993,  0.999892, 0.999964, 1,        0.999924, 1,        0.99996,  1,        1,	      // 6X
			  1,        0.999957, 1,        0.999955, 0.999955, 1,        1,        0.999951, 1,        0.999948, // 7X
			  1,        0.999946, 1,        1,        1,        1,        1,        0.999875, 1,        1,	      // 8X
			  0.999863, 1,        1,        1,        1,        1,        1,        1,        1,        1 };      // 9X
  int ptBin;
  
  for(unsigned int idx = 0; idx <= 4; ++idx){

    if(!fiveJet && idx == 4) break;
    ptBin = floor(jets->at(idx).pt());
    ptBin > 99 ? effy_ *= 1. : effy_ *= effies_[ptBin];
  }
  /*
  effy_ *= 0.5*TMath::Erf((jets->at(0).pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets->at(1).pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets->at(2).pt()-29.9592)/2.64442961111111)+0.5;
  effy_ *= 0.5*TMath::Erf((jets->at(3).pt()-29.9592)/2.64442961111111)+0.5;
  if (fiveJet) effy_ *= 0.5*TMath::Erf((jets->at(4).pt()-29.9592)/2.64442961111111)+0.5;
  */
  return effy_;
}

/// initiate histograms
void
SingleObjectJetTrigger::beginJob()
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  // pt jets
  hists_["pt" ] = fs->make<TH1F>( "pt"  , "pt"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1"] = fs->make<TH1F>( "pt1" , "pt1" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2"] = fs->make<TH1F>( "pt2" , "pt2" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3"] = fs->make<TH1F>( "pt3" , "pt3" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4"] = fs->make<TH1F>( "pt4" , "pt4" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5"] = fs->make<TH1F>( "pt5" , "pt5" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6"] = fs->make<TH1F>( "pt6" , "pt6" , 100,  0. , 100. );
  // jet multiplicity
  hists_["n"] = fs->make<TH1F>( "n" , "n" , 10,  -0.5 , 9.5 );
  // pt jets
  hists_["pt_match" ] = fs->make<TH1F>( "pt_match"  , "pt_match"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1_match"] = fs->make<TH1F>( "pt1_match" , "pt1_match" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2_match"] = fs->make<TH1F>( "pt2_match" , "pt2_match" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3_match"] = fs->make<TH1F>( "pt3_match" , "pt3_match" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4_match"] = fs->make<TH1F>( "pt4_match" , "pt4_match" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5_match"] = fs->make<TH1F>( "pt5_match" , "pt5_match" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6_match"] = fs->make<TH1F>( "pt6_match" , "pt6_match" , 100,  0. , 100. );
  // pt trigger jets
  hists_["pt_trigger" ] = fs->make<TH1F>( "pt_trigger"  , "pt_trigger"  , 100,  0. , 100. );
  // pt 1. trigger jet
  hists_["pt1_trigger"] = fs->make<TH1F>( "pt1_trigger" , "pt1_trigger" , 100,  0. , 100. );
  // pt 2. trigger jet
  hists_["pt2_trigger"] = fs->make<TH1F>( "pt2_trigger" , "pt2_trigger" , 100,  0. , 100. );
  // pt 3. trigger jet
  hists_["pt3_trigger"] = fs->make<TH1F>( "pt3_trigger" , "pt3_trigger" , 100,  0. , 100. );
  // pt 4. trigger jet
  hists_["pt4_trigger"] = fs->make<TH1F>( "pt4_trigger" , "pt4_trigger" , 100,  0. , 100. );
  // pt 5. trigger jet
  hists_["pt5_trigger"] = fs->make<TH1F>( "pt5_trigger" , "pt5_trigger" , 100,  0. , 100. );
  // pt 6. trigger jet
  hists_["pt6_trigger"] = fs->make<TH1F>( "pt6_trigger" , "pt6_trigger" , 100,  0. , 100. );
  // trigger jet multiplicity
  hists_["n_trigger"] = fs->make<TH1F>( "n_trigger" , "n_trigger" , 10,  -0.5 , 9.5 );
  // pt jets
  hists_["pt_effy" ] = fs->make<TH1F>( "pt_effy"  , "pt_effy"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1_effy"] = fs->make<TH1F>( "pt1_effy" , "pt1_effy" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2_effy"] = fs->make<TH1F>( "pt2_effy" , "pt2_effy" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3_effy"] = fs->make<TH1F>( "pt3_effy" , "pt3_effy" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4_effy"] = fs->make<TH1F>( "pt4_effy" , "pt4_effy" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5_effy"] = fs->make<TH1F>( "pt5_effy" , "pt5_effy" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6_effy"] = fs->make<TH1F>( "pt6_effy" , "pt6_effy" , 100,  0. , 100. );
  // Trigger vs. offline p_{T} (GeV)
  //hists2D_["ptTrigOff" ] = fs->make<TH2F>( "ptTrigOff" , "Trigger vs. offline p_{T} (GeV)", 60,  0., 300., 60,  0., 300. );
  // Trigger vs. offline eta
  //hists2D_["etaTrigOff"] = fs->make<TH2F>( "etaTrigOff", "Trigger vs. offline #eta"       , 50, -2.5, 2.5, 50, -2.5, 2.5 );
  // Trigger vs. offline phi
  //hists2D_["phiTrigOff"] = fs->make<TH2F>( "phiTrigOff", "Trigger vs. offline #phi"       , 64, -3.2, 3.2, 64, -3.2, 3.2 );

  /**
     Efficiencies vs. reco pt
  **/

  // QuadJet30 efficiency for 1. jet and normalization histogram
  hists_["QJ30_1"]   = fs->make<TH1F>( "QJ30_1"   , "QuadJet30_1"      , 100,  0. , 100. );
  hists_["QJ30_1_n"] = fs->make<TH1F>( "QJ30_1_n" , "QuadJet30_1_norm" , 100,  0. , 100. );
  hists_["QJ30_1_e"] = fs->make<TH1F>( "QJ30_1_e" , "QuadJet30_1_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 2. jet and normalization histogram
  hists_["QJ30_2"]   = fs->make<TH1F>( "QJ30_2"   , "QuadJet30_2"      , 100,  0. , 100. );
  hists_["QJ30_2_n"] = fs->make<TH1F>( "QJ30_2_n" , "QuadJet30_2_norm" , 100,  0. , 100. );
  hists_["QJ30_2_e"] = fs->make<TH1F>( "QJ30_2_e" , "QuadJet30_2_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 3. jet and normalization histogram
  hists_["QJ30_3"]   = fs->make<TH1F>( "QJ30_3"   , "QuadJet30_3"      , 100,  0. , 100. );
  hists_["QJ30_3_n"] = fs->make<TH1F>( "QJ30_3_n" , "QuadJet30_3_norm" , 100,  0. , 100. );
  hists_["QJ30_3_e"] = fs->make<TH1F>( "QJ30_3_e" , "QuadJet30_3_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 4. jet and normalization histogram
  hists_["QJ30_4"]   = fs->make<TH1F>( "QJ30_4"   , "QuadJet30_4"      , 100,  0. , 100. );
  hists_["QJ30_4_n"] = fs->make<TH1F>( "QJ30_4_n" , "QuadJet30_5_norm" , 100,  0. , 100. );
  hists_["QJ30_4_e"] = fs->make<TH1F>( "QJ30_4_e" , "QuadJet30_4_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 1. jet and normalization histogram
  hists_["QJ40_1"]   = fs->make<TH1F>( "QJ40_1"   , "QuadJet40_1"      , 100,  0. , 100. );
  hists_["QJ40_1_n"] = fs->make<TH1F>( "QJ40_1_n" , "QuadJet40_1_norm" , 100,  0. , 100. );
  hists_["QJ40_1_e"] = fs->make<TH1F>( "QJ40_1_e" , "QuadJet40_1_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 2. jet and normalization histogram
  hists_["QJ40_2"]   = fs->make<TH1F>( "QJ40_2"   , "QuadJet40_2"      , 100,  0. , 100. );
  hists_["QJ40_2_n"] = fs->make<TH1F>( "QJ40_2_n" , "QuadJet40_2_norm" , 100,  0. , 100. );
  hists_["QJ40_2_e"] = fs->make<TH1F>( "QJ40_2_e" , "QuadJet40_2_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 3. jet and normalization histogram
  hists_["QJ40_3"]   = fs->make<TH1F>( "QJ40_3"   , "QuadJet40_3"      , 100,  0. , 100. );
  hists_["QJ40_3_n"] = fs->make<TH1F>( "QJ40_3_n" , "QuadJet40_3_norm" , 100,  0. , 100. );
  hists_["QJ40_3_e"] = fs->make<TH1F>( "QJ40_3_e" , "QuadJet40_3_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 4. jet and normalization histogram
  hists_["QJ40_4"]   = fs->make<TH1F>( "QJ40_4"   , "QuadJet40_4"      , 100,  0. , 100. );
  hists_["QJ40_4_n"] = fs->make<TH1F>( "QJ40_4_n" , "QuadJet40_5_norm" , 100,  0. , 100. );
  hists_["QJ40_4_e"] = fs->make<TH1F>( "QJ40_4_e" , "QuadJet40_4_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 1. jet and normalization histogram
  hists_["QJ8333_1"]   = fs->make<TH1F>( "QJ8333_1"   , "QuadJet80303030_1"      , 100,  0. , 100. );
  hists_["QJ8333_1_n"] = fs->make<TH1F>( "QJ8333_1_n" , "QuadJet80303030_1_norm" , 100,  0. , 100. );
  hists_["QJ8333_1_e"] = fs->make<TH1F>( "QJ8333_1_e" , "QuadJet80303030_1_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 2. jet and normalization histogram
  hists_["QJ8333_2"]   = fs->make<TH1F>( "QJ8333_2"   , "QuadJet80303030_2"      , 100,  0. , 100. );
  hists_["QJ8333_2_n"] = fs->make<TH1F>( "QJ8333_2_n" , "QuadJet80303030_2_norm" , 100,  0. , 100. );
  hists_["QJ8333_2_e"] = fs->make<TH1F>( "QJ8333_2_e" , "QuadJet80303030_2_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 3. jet and normalization histogram
  hists_["QJ8333_3"]   = fs->make<TH1F>( "QJ8333_3"   , "QuadJet80303030_3"      , 100,  0. , 100. );
  hists_["QJ8333_3_n"] = fs->make<TH1F>( "QJ8333_3_n" , "QuadJet80303030_3_norm" , 100,  0. , 100. );
  hists_["QJ8333_3_e"] = fs->make<TH1F>( "QJ8333_3_e" , "QuadJet80303030_3_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 4. jet and normalization histogram
  hists_["QJ8333_4"]   = fs->make<TH1F>( "QJ8333_4"   , "QuadJet80303030_4"      , 100,  0. , 100. );
  hists_["QJ8333_4_n"] = fs->make<TH1F>( "QJ8333_4_n" , "QuadJet80303030_5_norm" , 100,  0. , 100. );
  hists_["QJ8333_4_e"] = fs->make<TH1F>( "QJ8333_4_e" , "QuadJet80303030_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 1. jet and normalization histogram
  hists_["FJ30_1"]   = fs->make<TH1F>( "FJ30_1"   , "FiveJet30_1"      , 100,  0. , 100. );
  hists_["FJ30_1_n"] = fs->make<TH1F>( "FJ30_1_n" , "FiveJet30_1_norm" , 100,  0. , 100. );
  hists_["FJ30_1_e"] = fs->make<TH1F>( "FJ30_1_e" , "FiveJet30_1_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 2. jet and normalization histogram
  hists_["FJ30_2"]   = fs->make<TH1F>( "FJ30_2"   , "FiveJet30_2"      , 100,  0. , 100. );
  hists_["FJ30_2_n"] = fs->make<TH1F>( "FJ30_2_n" , "FiveJet30_2_norm" , 100,  0. , 100. );
  hists_["FJ30_2_e"] = fs->make<TH1F>( "FJ30_2_e" , "FiveJet30_2_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 3. jet and normalization histogram
  hists_["FJ30_3"]   = fs->make<TH1F>( "FJ30_3"   , "FiveJet30_3"      , 100,  0. , 100. );
  hists_["FJ30_3_n"] = fs->make<TH1F>( "FJ30_3_n" , "FiveJet30_3_norm" , 100,  0. , 100. );
  hists_["FJ30_3_e"] = fs->make<TH1F>( "FJ30_3_e" , "FiveJet30_3_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 4. jet and normalization histogram
  hists_["FJ30_4"]   = fs->make<TH1F>( "FJ30_4"   , "FiveJet30_4"      , 100,  0. , 100. );
  hists_["FJ30_4_n"] = fs->make<TH1F>( "FJ30_4_n" , "FiveJet30_4_norm" , 100,  0. , 100. );
  hists_["FJ30_4_e"] = fs->make<TH1F>( "FJ30_4_e" , "FiveJet30_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 5. jet and normalization histogram
  hists_["FJ30_5"]   = fs->make<TH1F>( "FJ30_5"   , "FiveJet30_5"      , 100,  0. , 100. );
  hists_["FJ30_5_n"] = fs->make<TH1F>( "FJ30_5_n" , "FiveJet30_5_norm" , 100,  0. , 100. );
  hists_["FJ30_5_e"] = fs->make<TH1F>( "FJ30_5_e" , "FiveJet30_5_effy" , 100,  0. , 100. );

  /**
     Efficiencies vs. trigger pt
  **/

  // QuadJet30 efficiency for 1. jet and normalization histogram
  hists_["QJ30_trigger_1"]   = fs->make<TH1F>( "QJ30_trigger_1"   , "QuadJet30_trigger_1"      , 100,  0. , 100. );
  hists_["QJ30_trigger_1_n"] = fs->make<TH1F>( "QJ30_trigger_1_n" , "QuadJet30_trigger_1_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_1_e"] = fs->make<TH1F>( "QJ30_trigger_1_e" , "QuadJet30_trigger_1_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 2. jet and normalization histogram
  hists_["QJ30_trigger_2"]   = fs->make<TH1F>( "QJ30_trigger_2"   , "QuadJet30_trigger_2"      , 100,  0. , 100. );
  hists_["QJ30_trigger_2_n"] = fs->make<TH1F>( "QJ30_trigger_2_n" , "QuadJet30_trigger_2_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_2_e"] = fs->make<TH1F>( "QJ30_trigger_2_e" , "QuadJet30_trigger_2_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 3. jet and normalization histogram
  hists_["QJ30_trigger_3"]   = fs->make<TH1F>( "QJ30_trigger_3"   , "QuadJet30_trigger_3"      , 100,  0. , 100. );
  hists_["QJ30_trigger_3_n"] = fs->make<TH1F>( "QJ30_trigger_3_n" , "QuadJet30_trigger_3_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_3_e"] = fs->make<TH1F>( "QJ30_trigger_3_e" , "QuadJet30_trigger_3_effy" , 100,  0. , 100. );
  // QuadJet30 efficiency for 4. jet and normalization histogram
  hists_["QJ30_trigger_4"]   = fs->make<TH1F>( "QJ30_trigger_4"   , "QuadJet30_trigger_4"      , 100,  0. , 100. );
  hists_["QJ30_trigger_4_n"] = fs->make<TH1F>( "QJ30_trigger_4_n" , "QuadJet30_trigger_5_norm" , 100,  0. , 100. );
  hists_["QJ30_trigger_4_e"] = fs->make<TH1F>( "QJ30_trigger_4_e" , "QuadJet30_trigger_4_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 1. jet and normalization histogram
  hists_["QJ40_trigger_1"]   = fs->make<TH1F>( "QJ40_trigger_1"   , "QuadJet40_trigger_1"      , 100,  0. , 100. );
  hists_["QJ40_trigger_1_n"] = fs->make<TH1F>( "QJ40_trigger_1_n" , "QuadJet40_trigger_1_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_1_e"] = fs->make<TH1F>( "QJ40_trigger_1_e" , "QuadJet40_trigger_1_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 2. jet and normalization histogram
  hists_["QJ40_trigger_2"]   = fs->make<TH1F>( "QJ40_trigger_2"   , "QuadJet40_trigger_2"      , 100,  0. , 100. );
  hists_["QJ40_trigger_2_n"] = fs->make<TH1F>( "QJ40_trigger_2_n" , "QuadJet40_trigger_2_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_2_e"] = fs->make<TH1F>( "QJ40_trigger_2_e" , "QuadJet40_trigger_2_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 3. jet and normalization histogram
  hists_["QJ40_trigger_3"]   = fs->make<TH1F>( "QJ40_trigger_3"   , "QuadJet40_trigger_3"      , 100,  0. , 100. );
  hists_["QJ40_trigger_3_n"] = fs->make<TH1F>( "QJ40_trigger_3_n" , "QuadJet40_trigger_3_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_3_e"] = fs->make<TH1F>( "QJ40_trigger_3_e" , "QuadJet40_trigger_3_effy" , 100,  0. , 100. );
  // QuadJet40 efficiency for 4. jet and normalization histogram
  hists_["QJ40_trigger_4"]   = fs->make<TH1F>( "QJ40_trigger_4"   , "QuadJet40_trigger_4"      , 100,  0. , 100. );
  hists_["QJ40_trigger_4_n"] = fs->make<TH1F>( "QJ40_trigger_4_n" , "QuadJet40_trigger_5_norm" , 100,  0. , 100. );
  hists_["QJ40_trigger_4_e"] = fs->make<TH1F>( "QJ40_trigger_4_e" , "QuadJet40_trigger_4_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 1. jet and normalization histogram
  hists_["QJ8333_trigger_1"]   = fs->make<TH1F>( "QJ8333_trigger_1"   , "QuadJet80303030_trigger_1"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_1_n"] = fs->make<TH1F>( "QJ8333_trigger_1_n" , "QuadJet80303030_trigger_1_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_1_e"] = fs->make<TH1F>( "QJ8333_trigger_1_e" , "QuadJet80303030_trigger_1_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 2. jet and normalization histogram
  hists_["QJ8333_trigger_2"]   = fs->make<TH1F>( "QJ8333_trigger_2"   , "QuadJet80303030_trigger_2"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_2_n"] = fs->make<TH1F>( "QJ8333_trigger_2_n" , "QuadJet80303030_trigger_2_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_2_e"] = fs->make<TH1F>( "QJ8333_trigger_2_e" , "QuadJet80303030_trigger_2_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 3. jet and normalization histogram
  hists_["QJ8333_trigger_3"]   = fs->make<TH1F>( "QJ8333_trigger_3"   , "QuadJet80303030_trigger_3"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_3_n"] = fs->make<TH1F>( "QJ8333_trigger_3_n" , "QuadJet80303030_trigger_3_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_3_e"] = fs->make<TH1F>( "QJ8333_trigger_3_e" , "QuadJet80303030_trigger_3_effy" , 100,  0. , 100. );
  // QuadJet80303030 efficiency for 4. jet and normalization histogram
  hists_["QJ8333_trigger_4"]   = fs->make<TH1F>( "QJ8333_trigger_4"   , "QuadJet80303030_trigger_4"      , 100,  0. , 100. );
  hists_["QJ8333_trigger_4_n"] = fs->make<TH1F>( "QJ8333_trigger_4_n" , "QuadJet80303030_trigger_5_norm" , 100,  0. , 100. );
  hists_["QJ8333_trigger_4_e"] = fs->make<TH1F>( "QJ8333_trigger_4_e" , "QuadJet80303030_trigger_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 1. jet and normalization histogram
  hists_["FJ30_trigger_1"]   = fs->make<TH1F>( "FJ30_trigger_1"   , "FiveJet30_trigger_1"      , 100,  0. , 100. );
  hists_["FJ30_trigger_1_n"] = fs->make<TH1F>( "FJ30_trigger_1_n" , "FiveJet30_trigger_1_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_1_e"] = fs->make<TH1F>( "FJ30_trigger_1_e" , "FiveJet30_trigger_1_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 2. jet and normalization histogram
  hists_["FJ30_trigger_2"]   = fs->make<TH1F>( "FJ30_trigger_2"   , "FiveJet30_trigger_2"      , 100,  0. , 100. );
  hists_["FJ30_trigger_2_n"] = fs->make<TH1F>( "FJ30_trigger_2_n" , "FiveJet30_trigger_2_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_2_e"] = fs->make<TH1F>( "FJ30_trigger_2_e" , "FiveJet30_trigger_2_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 3. jet and normalization histogram
  hists_["FJ30_trigger_3"]   = fs->make<TH1F>( "FJ30_trigger_3"   , "FiveJet30_trigger_3"      , 100,  0. , 100. );
  hists_["FJ30_trigger_3_n"] = fs->make<TH1F>( "FJ30_trigger_3_n" , "FiveJet30_trigger_3_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_3_e"] = fs->make<TH1F>( "FJ30_trigger_3_e" , "FiveJet30_trigger_3_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 4. jet and normalization histogram
  hists_["FJ30_trigger_4"]   = fs->make<TH1F>( "FJ30_trigger_4"   , "FiveJet30_trigger_4"      , 100,  0. , 100. );
  hists_["FJ30_trigger_4_n"] = fs->make<TH1F>( "FJ30_trigger_4_n" , "FiveJet30_trigger_4_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_4_e"] = fs->make<TH1F>( "FJ30_trigger_4_e" , "FiveJet30_trigger_4_effy" , 100,  0. , 100. );
  // FiveJet30 efficiency for 5. jet and normalization histogram
  hists_["FJ30_trigger_5"]   = fs->make<TH1F>( "FJ30_trigger_5"   , "FiveJet30_trigger_5"      , 100,  0. , 100. );
  hists_["FJ30_trigger_5_n"] = fs->make<TH1F>( "FJ30_trigger_5_n" , "FiveJet30_trigger_5_norm" , 100,  0. , 100. );
  hists_["FJ30_trigger_5_e"] = fs->make<TH1F>( "FJ30_trigger_5_e" , "FiveJet30_trigger_5_effy" , 100,  0. , 100. );
 }

/// analyze triggers and fill histograms
void
SingleObjectJetTrigger::analyze(const edm::Event& event, const edm::EventSetup&)
{
  //edm::Handle<std::vector<pat::TriggerFilter> > triggerFilters;
  //event.getByLabel(patTrigger_, triggerFilters);
  //edm::Handle<std::vector<pat::TriggerPath> > triggerPaths;
  //event.getByLabel(patTrigger_, triggerPaths);
  //edm::Handle<std::vector<pat::TriggerObject> > triggerObjects;
  //event.getByLabel(patTrigger_, triggerObjects);
  edm::Handle<pat::TriggerEvent> triggerEvent;
  event.getByLabel(patTriggerEvent_, triggerEvent);
  //edm::Handle<edm::View<pat::Jet> > jets;
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jets_, jets);
  //edm::Handle<edm::View<pat::Jet> > triggerMatchedJets;
  //event.getByLabel(triggerMatchedJets_, triggerMatchedJets);

  if( !triggerEvent.failedToGet() ) {
    // PAT trigger helper for trigger matching information
    const pat::helper::TriggerMatchHelper matchHelper;

    const pat::TriggerObjectMatch * triggerMatch( triggerEvent->triggerObjectMatchResult( "jetTriggerMatchHLTJets" ) );
    unsigned int jetCounter = 1;
    hists_.find("n")->second->Fill( jets->size() );

    double effy = 1;
    double jet1_pt = 0;
    double jet2_pt = 0;
    double jet3_pt = 0;
    double jet4_pt = 0;
    double jet5_pt = 0;
    double triggerJet1_pt = 0;
    double triggerJet2_pt = 0;
    double triggerJet3_pt = 0;
    double triggerJet4_pt = 0;
    double triggerJet5_pt = 0;
    
    for (size_t iJet = 0; iJet < jets->size(); ++iJet, ++jetCounter) {
      double jet_pt = (*jets)[iJet].pt();
      hists_.find("pt")->second->Fill( jet_pt );
      if     (jetCounter == 1) { hists_.find("pt1")->second->Fill( jet_pt ); jet1_pt = jet_pt; }
      else if(jetCounter == 2) { hists_.find("pt2")->second->Fill( jet_pt ); jet2_pt = jet_pt; }
      else if(jetCounter == 3) { hists_.find("pt3")->second->Fill( jet_pt ); jet3_pt = jet_pt; }
      else if(jetCounter == 4) { hists_.find("pt4")->second->Fill( jet_pt ); jet4_pt = jet_pt; }
      else if(jetCounter == 5) { hists_.find("pt5")->second->Fill( jet_pt ); jet5_pt = jet_pt; }
      else if(jetCounter == 6)   hists_.find("pt6")->second->Fill( jet_pt );

      const reco::CandidateBaseRef candBaseRef( pat::JetRef( jets , iJet) );
      const pat::TriggerObjectRef trigRef( matchHelper.triggerMatchObject( candBaseRef, triggerMatch, event, *triggerEvent ) );

      if ( trigRef.isAvailable() ) {
	double triggerMatchedJet_pt = candBaseRef->pt();
	hists_.find("pt_match")->second->Fill( triggerMatchedJet_pt );
	if     (jetCounter == 1) { hists_.find("pt1_match")->second->Fill( triggerMatchedJet_pt ); /* triggerJet1_pt = trigRef->pt(); */ }
	else if(jetCounter == 2) { hists_.find("pt2_match")->second->Fill( triggerMatchedJet_pt ); /* triggerJet2_pt = trigRef->pt(); */ }
	else if(jetCounter == 3) { hists_.find("pt3_match")->second->Fill( triggerMatchedJet_pt ); /* triggerJet3_pt = trigRef->pt(); */ }
	else if(jetCounter == 4) { hists_.find("pt4_match")->second->Fill( triggerMatchedJet_pt ); /* triggerJet4_pt = trigRef->pt(); */ }
	else if(jetCounter == 5) { hists_.find("pt5_match")->second->Fill( triggerMatchedJet_pt ); /* triggerJet5_pt = trigRef->pt(); */ }
	else if(jetCounter == 6)   hists_.find("pt6_match")->second->Fill( triggerMatchedJet_pt );
      }
    }

    pat::TriggerFilterRefVector filters = triggerEvent->acceptedFilters();
    for(pat::TriggerFilterRefVector::const_iterator filter = filters.begin(); filter != filters.end(); ++filter){
      if((*filter)->label() == "hlt4jet30"){
	
	pat::TriggerObjectRefVector objects = triggerEvent->filterObjects("hlt4jet30");

	hists_.find("n_trigger")->second->Fill( objects.size() );
	unsigned int triggerJetCounter = 1;
	for (size_t iTrigJet = 0; iTrigJet < objects.size(); ++iTrigJet, ++triggerJetCounter) {
	  double triggerJet_pt = objects[iTrigJet]->pt();
	  hists_.find("pt_trigger")->second->Fill( triggerJet_pt );
	  if(triggerJetCounter == 1) { hists_.find("pt1_trigger")->second->Fill( triggerJet_pt ); triggerJet1_pt = triggerJet_pt; }
	  if(triggerJetCounter == 2) { hists_.find("pt2_trigger")->second->Fill( triggerJet_pt ); triggerJet2_pt = triggerJet_pt; }
	  if(triggerJetCounter == 3) { hists_.find("pt3_trigger")->second->Fill( triggerJet_pt ); triggerJet3_pt = triggerJet_pt; }
	  if(triggerJetCounter == 4) { hists_.find("pt4_trigger")->second->Fill( triggerJet_pt ); triggerJet4_pt = triggerJet_pt; }
	  if(triggerJetCounter == 5) { hists_.find("pt5_trigger")->second->Fill( triggerJet_pt ); triggerJet5_pt = triggerJet_pt; }
	  if(triggerJetCounter == 6) { hists_.find("pt6_trigger")->second->Fill( triggerJet_pt ); }
	}
	
	if(jets->size() >= 4 && objects[0]->pt() >= 30. && objects[1]->pt() >= 30. && objects[2]->pt() >= 30. && objects[3]->pt() >= 30.) {
	  
	  effy = triggerEfficiency(jets, false);
	  
	  hists_.find("QJ30_1")->second->Fill( jet1_pt, effy );
	  hists_.find("QJ30_2")->second->Fill( jet2_pt, effy );
	  hists_.find("QJ30_3")->second->Fill( jet3_pt, effy );
	  hists_.find("QJ30_4")->second->Fill( jet4_pt, effy );
	  hists_.find("QJ30_1_n")->second->Fill( jet1_pt );
	  hists_.find("QJ30_2_n")->second->Fill( jet2_pt );
	  hists_.find("QJ30_3_n")->second->Fill( jet3_pt );
	  hists_.find("QJ30_4_n")->second->Fill( jet4_pt );

	  hists_.find("QJ30_trigger_1")->second->Fill( triggerJet1_pt, effy );
	  hists_.find("QJ30_trigger_2")->second->Fill( triggerJet2_pt, effy );
	  hists_.find("QJ30_trigger_3")->second->Fill( triggerJet3_pt, effy );
	  hists_.find("QJ30_trigger_4")->second->Fill( triggerJet4_pt, effy );
	  hists_.find("QJ30_trigger_1_n")->second->Fill( triggerJet1_pt );
	  hists_.find("QJ30_trigger_2_n")->second->Fill( triggerJet2_pt );
	  hists_.find("QJ30_trigger_3_n")->second->Fill( triggerJet3_pt );
	  hists_.find("QJ30_trigger_4_n")->second->Fill( triggerJet4_pt );

	  if(jets->size() >= 5 && objects.size() >= 5){
	    if(objects[4]->pt() >= 30.){

	      effy = triggerEfficiency(jets, true);

	      hists_.find("FJ30_1")->second->Fill( jet1_pt, effy );
	      hists_.find("FJ30_2")->second->Fill( jet2_pt, effy );
	      hists_.find("FJ30_3")->second->Fill( jet3_pt, effy );
	      hists_.find("FJ30_4")->second->Fill( jet4_pt, effy );
	      hists_.find("FJ30_5")->second->Fill( jet5_pt, effy );
	      hists_.find("FJ30_1_n")->second->Fill( jet1_pt );
	      hists_.find("FJ30_2_n")->second->Fill( jet2_pt );
	      hists_.find("FJ30_3_n")->second->Fill( jet3_pt );
	      hists_.find("FJ30_4_n")->second->Fill( jet4_pt );
	      hists_.find("FJ30_5_n")->second->Fill( jet5_pt );


	      hists_.find("FJ30_trigger_1")->second->Fill( triggerJet1_pt, effy );
	      hists_.find("FJ30_trigger_2")->second->Fill( triggerJet2_pt, effy );
	      hists_.find("FJ30_trigger_3")->second->Fill( triggerJet3_pt, effy );
	      hists_.find("FJ30_trigger_4")->second->Fill( triggerJet4_pt, effy );
	      hists_.find("FJ30_trigger_5")->second->Fill( triggerJet5_pt, effy );
	      hists_.find("FJ30_trigger_1_n")->second->Fill( triggerJet1_pt );
	      hists_.find("FJ30_trigger_2_n")->second->Fill( triggerJet2_pt );
	      hists_.find("FJ30_trigger_3_n")->second->Fill( triggerJet3_pt );
	      hists_.find("FJ30_trigger_4_n")->second->Fill( triggerJet4_pt );
	      hists_.find("FJ30_trigger_5_n")->second->Fill( triggerJet5_pt );

	      effy = triggerEfficiency(jets, false);

	    }
	  }
	}
    
	if(jets->size() >= 4 && objects[0]->pt() >= 40. && objects[1]->pt() >= 40. && objects[2]->pt() >= 40. && objects[3]->pt() >= 40.) {
	  hists_.find("QJ40_1")->second->Fill( jet1_pt, effy );
	  hists_.find("QJ40_2")->second->Fill( jet2_pt, effy );
	  hists_.find("QJ40_3")->second->Fill( jet3_pt, effy );
	  hists_.find("QJ40_4")->second->Fill( jet4_pt, effy );
	  hists_.find("QJ40_1_n")->second->Fill( jet1_pt );
	  hists_.find("QJ40_2_n")->second->Fill( jet2_pt );
	  hists_.find("QJ40_3_n")->second->Fill( jet3_pt );
	  hists_.find("QJ40_4_n")->second->Fill( jet4_pt );

	  hists_.find("QJ40_trigger_1")->second->Fill( triggerJet1_pt, effy );
	  hists_.find("QJ40_trigger_2")->second->Fill( triggerJet2_pt, effy );
	  hists_.find("QJ40_trigger_3")->second->Fill( triggerJet3_pt, effy );
	  hists_.find("QJ40_trigger_4")->second->Fill( triggerJet4_pt, effy );
	  hists_.find("QJ40_trigger_1_n")->second->Fill( triggerJet1_pt );
	  hists_.find("QJ40_trigger_2_n")->second->Fill( triggerJet2_pt );
	  hists_.find("QJ40_trigger_3_n")->second->Fill( triggerJet3_pt );
	  hists_.find("QJ40_trigger_4_n")->second->Fill( triggerJet4_pt );
	}

	if(jets->size() >= 4 && objects[0]->pt() >= 80. && objects[1]->pt() >= 30. && objects[2]->pt() >= 30. && objects[3]->pt() >= 30.) {
	  hists_.find("QJ8333_1")->second->Fill( jet1_pt, effy );
	  hists_.find("QJ8333_2")->second->Fill( jet2_pt, effy );
	  hists_.find("QJ8333_3")->second->Fill( jet3_pt, effy );
	  hists_.find("QJ8333_4")->second->Fill( jet4_pt, effy );
	  hists_.find("QJ8333_1_n")->second->Fill( jet1_pt );
	  hists_.find("QJ8333_2_n")->second->Fill( jet2_pt );
	  hists_.find("QJ8333_3_n")->second->Fill( jet3_pt );
	  hists_.find("QJ8333_4_n")->second->Fill( jet4_pt );

	  hists_.find("QJ8333_trigger_1")->second->Fill( triggerJet1_pt, effy );
	  hists_.find("QJ8333_trigger_2")->second->Fill( triggerJet2_pt, effy );
	  hists_.find("QJ8333_trigger_3")->second->Fill( triggerJet3_pt, effy );
	  hists_.find("QJ8333_trigger_4")->second->Fill( triggerJet4_pt, effy );
	  hists_.find("QJ8333_trigger_1_n")->second->Fill( triggerJet1_pt );
	  hists_.find("QJ8333_trigger_2_n")->second->Fill( triggerJet2_pt );
	  hists_.find("QJ8333_trigger_3_n")->second->Fill( triggerJet3_pt );
	  hists_.find("QJ8333_trigger_4_n")->second->Fill( triggerJet4_pt );
	}
      }
    }
  }
}

void
SingleObjectJetTrigger::endJob()
{
  std::vector<std::string> hists2Norm;
  hists2Norm.push_back("pt");
  hists2Norm.push_back("pt1");
  hists2Norm.push_back("pt2");
  hists2Norm.push_back("pt3");
  hists2Norm.push_back("pt4");
  hists2Norm.push_back("pt5");
  hists2Norm.push_back("pt6");
  for( std::vector<std::string>::const_iterator hist = hists2Norm.begin(); hist != hists2Norm.end(); ++hist ){
    hists_.find(*hist + "_effy")->second->Divide(hists_.find(*hist + "_match")->second, hists_.find(*hist)->second);
  }

  hists2Norm.clear();
  hists2Norm.push_back("QJ30_1");
  hists2Norm.push_back("QJ40_1");
  hists2Norm.push_back("QJ8333_1");
  hists2Norm.push_back("FJ30_1");
  hists2Norm.push_back("QJ30_2");
  hists2Norm.push_back("QJ40_2");
  hists2Norm.push_back("QJ8333_2");
  hists2Norm.push_back("FJ30_2");
  hists2Norm.push_back("QJ30_3");
  hists2Norm.push_back("QJ40_3");
  hists2Norm.push_back("QJ8333_3");
  hists2Norm.push_back("FJ30_3");
  hists2Norm.push_back("QJ30_4");
  hists2Norm.push_back("QJ40_4");
  hists2Norm.push_back("QJ8333_4");
  hists2Norm.push_back("FJ30_4");
  hists2Norm.push_back("FJ30_5");
  hists2Norm.push_back("QJ30_trigger_1");
  hists2Norm.push_back("QJ40_trigger_1");
  hists2Norm.push_back("QJ8333_trigger_1");
  hists2Norm.push_back("FJ30_trigger_1");
  hists2Norm.push_back("QJ30_trigger_2");
  hists2Norm.push_back("QJ40_trigger_2");
  hists2Norm.push_back("QJ8333_trigger_2");
  hists2Norm.push_back("FJ30_trigger_2");
  hists2Norm.push_back("QJ30_trigger_3");
  hists2Norm.push_back("QJ40_trigger_3");
  hists2Norm.push_back("QJ8333_trigger_3");
  hists2Norm.push_back("FJ30_trigger_3");
  hists2Norm.push_back("QJ30_trigger_4");
  hists2Norm.push_back("QJ40_trigger_4");
  hists2Norm.push_back("QJ8333_trigger_4");
  hists2Norm.push_back("FJ30_trigger_4");
  hists2Norm.push_back("FJ30_trigger_5");
  for( std::vector<std::string>::const_iterator hist = hists2Norm.begin(); hist != hists2Norm.end(); ++hist ){
    hists_.find(*hist + "_e")->second->Divide(hists_.find(*hist)->second, hists_.find(*hist + "_n")->second);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( SingleObjectJetTrigger );
