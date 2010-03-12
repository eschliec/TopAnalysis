#include <TMath.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "TopAnalysis/TopAnalyzer/interface/KinFitImprover.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullHadEvtPartons.h"

/// default constructor for fw lite
KinFitImprover::KinFitImprover()
{
}

/// default constructor for full fw
KinFitImprover::KinFitImprover(const edm::ParameterSet& cfg)
{
}

/// histogramm booking for fwlite 
void KinFitImprover::book()
{
  /** 
      
  **/

  // show how many fits are right, only have swapped b-quarks or are totally wrong
  hists_["fitQuality"] = new TH1F( "fitQuality" , "fitQuality" ,    3 , 1. ,    4. );
  // top mass for right combinations
  hists_["topMass1"  ] = new TH1F( "topMass1"   , "topMass1"   ,  100 , 0. , 1000. );
  // top mass for swapped b-quark combinations
  hists_["topMass2"  ] = new TH1F( "topMass2"   , "topMass2"   ,  100 , 0. , 1000. );
  // top mass for wrong combinations
  hists_["topMass3"  ] = new TH1F( "topMass3"   , "topMass3"   ,  100 , 0. , 1000. );
  // min deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists_["minDRFit1" ] = new TH1F( "minDRFit1"  , "minDRFit1"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists_["minDRFit2" ] = new TH1F( "minDRFit2"  , "minDRFit2"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists_["minDRFit3" ] = new TH1F( "minDRFit3"  , "minDRFit3"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists_["minDRRec1" ] = new TH1F( "minDRRec1"  , "minDRRec1"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (reco) for b-swapped combinations
  hists_["minDRRec2" ] = new TH1F( "minDRRec2"  , "minDRRec2"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (reco) for wrong combinations
  hists_["minDRRec3" ] = new TH1F( "minDRRec3"  , "minDRRec3"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["minDRFitRec1"] = new TH2F( "minDRFitRec1" , "minDRFitRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists2D_["minDRFitRec2"] = new TH2F( "minDRFitRec2" , "minDRFitRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists2D_["minDRFitRec3"] = new TH2F( "minDRFitRec3" , "minDRFitRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists_["maxDRFit1" ] = new TH1F( "maxDRFit1"  , "maxDRFit1"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists_["maxDRFit2" ] = new TH1F( "maxDRFit2"  , "maxDRFit2"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists_["maxDRFit3" ] = new TH1F( "maxDRFit3"  , "maxDRFit3"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists_["maxDRRec1" ] = new TH1F( "maxDRRec1"  , "maxDRRec1"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (reco) for b-swapped combinations
  hists_["maxDRRec2" ] = new TH1F( "maxDRRec2"  , "maxDRRec2"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (reco) for wrong combinations
  hists_["maxDRRec3" ] = new TH1F( "maxDRRec3"  , "maxDRRec3"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["maxDRFitRec1"] = new TH2F( "maxDRFitRec1" , "maxDRFitRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // max deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists2D_["maxDRFitRec2"] = new TH2F( "maxDRFitRec2" , "maxDRFitRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // max deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists2D_["maxDRFitRec3"] = new TH2F( "maxDRFitRec3" , "maxDRFitRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // all deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists_["allDRFit1" ] = new TH1F( "allDRFit1"  , "allDRFit1"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists_["allDRFit2" ] = new TH1F( "allDRFit2"  , "allDRFit2"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists_["allDRFit3" ] = new TH1F( "allDRFit3"  , "allDRFit3"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists_["allDRRec1" ] = new TH1F( "allDRRec1"  , "allDRRec1"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (reco) for b-swapped combinations
  hists_["allDRRec2" ] = new TH1F( "allDRRec2"  , "allDRRec2"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (reco) for wrong combinations
  hists_["allDRRec3" ] = new TH1F( "allDRRec3"  , "allDRRec3"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["allDRFitRec1"] = new TH2F( "allDRFitRec1" , "allDRFitRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // all deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists2D_["allDRFitRec2"] = new TH2F( "allDRFitRec2" , "allDRFitRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // all deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists2D_["allDRFitRec3"] = new TH2F( "allDRFitRec3" , "allDRFitRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
}

/// histogramm booking for fw
void KinFitImprover::book(edm::Service<TFileService>& fs)
{
  /** 
      
  **/

  // show how many fits are right, only have swapped b-quarks or are totally wrong
  hists_["fitQuality"] = fs->make<TH1F>( "fitQuality" , "fitQuality" ,    3 , 1. ,    4. );
  // top mass for right combinations
  hists_["topMass1"  ] = fs->make<TH1F>( "topMass1"   , "topMass1"   ,  100 , 0. , 1000. );
  // top mass for swapped b-quark combinations
  hists_["topMass2"  ] = fs->make<TH1F>( "topMass2"   , "topMass2"   ,  100 , 0. , 1000. );
  // top mass for wrong combinations
  hists_["topMass3"  ] = fs->make<TH1F>( "topMass3"   , "topMass3"   ,  100 , 0. , 1000. );
  // min deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists_["minDRFit1" ] = fs->make<TH1F>( "minDRFit1"  , "minDRFit1"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists_["minDRFit2" ] = fs->make<TH1F>( "minDRFit2"  , "minDRFit2"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists_["minDRFit3" ] = fs->make<TH1F>( "minDRFit3"  , "minDRFit3"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists_["minDRRec1" ] = fs->make<TH1F>( "minDRRec1"  , "minDRRec1"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (reco) for b-swapped combinations
  hists_["minDRRec2" ] = fs->make<TH1F>( "minDRRec2"  , "minDRRec2"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (reco) for wrong combinations
  hists_["minDRRec3" ] = fs->make<TH1F>( "minDRRec3"  , "minDRRec3"  ,  100 , 0. ,   10. );
  // min deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["minDRFitRec1"] = fs->make<TH2F>( "minDRFitRec1" , "minDRFitRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists2D_["minDRFitRec2"] = fs->make<TH2F>( "minDRFitRec2" , "minDRFitRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists2D_["minDRFitRec3"] = fs->make<TH2F>( "minDRFitRec3" , "minDRFitRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists_["maxDRFit1" ] = fs->make<TH1F>( "maxDRFit1"  , "maxDRFit1"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists_["maxDRFit2" ] = fs->make<TH1F>( "maxDRFit2"  , "maxDRFit2"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists_["maxDRFit3" ] = fs->make<TH1F>( "maxDRFit3"  , "maxDRFit3"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists_["maxDRRec1" ] = fs->make<TH1F>( "maxDRRec1"  , "maxDRRec1"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (reco) for b-swapped combinations
  hists_["maxDRRec2" ] = fs->make<TH1F>( "maxDRRec2"  , "maxDRRec2"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (reco) for wrong combinations
  hists_["maxDRRec3" ] = fs->make<TH1F>( "maxDRRec3"  , "maxDRRec3"  ,  100 , 0. ,   10. );
  // max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["maxDRFitRec1"] = fs->make<TH2F>( "maxDRFitRec1" , "maxDRFitRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // max deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists2D_["maxDRFitRec2"] = fs->make<TH2F>( "maxDRFitRec2" , "maxDRFitRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // max deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists2D_["maxDRFitRec3"] = fs->make<TH2F>( "maxDRFitRec3" , "maxDRFitRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min vs. max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["minMaxDRFit1"] = fs->make<TH2F>( "minMaxDRFit1" , "minMaxDRFit1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min vs. max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["minMaxDRFit2"] = fs->make<TH2F>( "minMaxDRFit2" , "minMaxDRFit2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min vs. max deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["minMaxDRFit3"] = fs->make<TH2F>( "minMaxDRFit3" , "minMaxDRFit3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min vs. max deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists2D_["minMaxDRRec1"] = fs->make<TH2F>( "minMaxDRRec1" , "minMaxDRRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min vs. max deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists2D_["minMaxDRRec2"] = fs->make<TH2F>( "minMaxDRRec2" , "minMaxDRRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // min vs. max deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists2D_["minMaxDRRec3"] = fs->make<TH2F>( "minMaxDRRec3" , "minMaxDRRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // all deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists_["allDRFit1" ] = fs->make<TH1F>( "allDRFit1"  , "allDRFit1"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists_["allDRFit2" ] = fs->make<TH1F>( "allDRFit2"  , "allDRFit2"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists_["allDRFit3" ] = fs->make<TH1F>( "allDRFit3"  , "allDRFit3"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (reco) for right combinations
  hists_["allDRRec1" ] = fs->make<TH1F>( "allDRRec1"  , "allDRRec1"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (reco) for b-swapped combinations
  hists_["allDRRec2" ] = fs->make<TH1F>( "allDRRec2"  , "allDRRec2"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (reco) for wrong combinations
  hists_["allDRRec3" ] = fs->make<TH1F>( "allDRRec3"  , "allDRRec3"  ,  100 , 0. ,   10. );
  // all deltaR between b quark and corresponding light quarks (fit) for right combinations
  hists2D_["allDRFitRec1"] = fs->make<TH2F>( "allDRFitRec1" , "allDRFitRec1" , 100 , 0. , 10. , 100 , 0. , 10. );
  // all deltaR between b quark and corresponding light quarks (fit) for b-swapped combinations
  hists2D_["allDRFitRec2"] = fs->make<TH2F>( "allDRFitRec2" , "allDRFitRec2" , 100 , 0. , 10. , 100 , 0. , 10. );
  // all deltaR between b quark and corresponding light quarks (fit) for wrong combinations
  hists2D_["allDRFitRec3"] = fs->make<TH2F>( "allDRFitRec3" , "allDRFitRec3" , 100 , 0. , 10. , 100 , 0. , 10. );
  // energy of process for right combinations
  hists_["sqrt_s_fit1"  ] = fs->make<TH1F>( "sqrt_s_fit1"   , "sqrt_s_fit1"   ,  500 , 0. ,   5000. );
  // pt of final system for right combinations
  hists_["pt_final_fit1"] = fs->make<TH1F>( "pt_final_fit1" , "pt_final_fit1" ,  100 , 0. ,   1000. );
  // energy of process for b-swapped combinations
  hists_["sqrt_s_fit2"  ] = fs->make<TH1F>( "sqrt_s_fit2"   , "sqrt_s_fit2"   ,  500 , 0. ,   5000. );
  // pt of final system for b-swapped combinations
  hists_["pt_final_fit2"] = fs->make<TH1F>( "pt_final_fit2" , "pt_final_fit2" ,  100 , 0. ,   1000. );
  // energy of process for wrong combinations
  hists_["sqrt_s_fit3"  ] = fs->make<TH1F>( "sqrt_s_fit3"   , "sqrt_s_fit3"   ,  500 , 0. ,   5000. );
  // pt of final system for wrong combinations
  hists_["pt_final_fit3"] = fs->make<TH1F>( "pt_final_fit3" , "pt_final_fit3" ,  100 , 0. ,   1000. );
  // energy of process for right combinations
  hists_["sqrt_s_rec1"  ] = fs->make<TH1F>( "sqrt_s_rec1"   , "sqrt_s_rec1"   ,  500 , 0. ,   5000. );
  // pt of final system for right combinations
  hists_["pt_final_rec1"] = fs->make<TH1F>( "pt_final_rec1" , "pt_final_rec1" ,  100 , 0. ,   1000. );
  // energy of process for b-swapped combinations
  hists_["sqrt_s_rec2"  ] = fs->make<TH1F>( "sqrt_s_rec2"   , "sqrt_s_rec2"   ,  500 , 0. ,   5000. );
  // pt of final system for b-swapped combinations
  hists_["pt_final_rec2"] = fs->make<TH1F>( "pt_final_rec2" , "pt_final_rec2" ,  100 , 0. ,   1000. );
  // energy of process for wrong combinations
  hists_["sqrt_s_rec3"  ] = fs->make<TH1F>( "sqrt_s_rec3"   , "sqrt_s_rec3"   ,  500 , 0. ,   5000. );
  // pt of final system for wrong combinations
  hists_["pt_final_rec3"] = fs->make<TH1F>( "pt_final_rec3" , "pt_final_rec3" ,  100 , 0. ,   1000. );
}

/// function to find types of jet-combinations in KinFits (1 right, 2 swapped b-quarks, 3 wrong)
int
KinFitImprover::comboType()
{
  int lightQFit    = jetIndex[0];
  int lightQBarFit = jetIndex[1];
  int bFit         = jetIndex[2];
  int bBarFit      = jetIndex[3];
  int lightPFit    = jetIndex[4];
  int lightPBarFit = jetIndex[5];

  int lightQGen    = jetIndex[6];
  int lightQBarGen = jetIndex[7];
  int bGen         = jetIndex[8];
  int bBarGen      = jetIndex[9];
  int lightPGen    = jetIndex[10];
  int lightPBarGen = jetIndex[11];

  if( ( bFit == bGen && bBarFit == bBarGen ) || ( bFit == bBarGen && bBarFit == bGen ) ) {

	if( bFit == bGen &&
	    ( ( lightQFit == lightQGen    && lightQBarFit == lightQBarGen ) ||
	      ( lightQFit == lightQBarGen && lightQBarFit == lightQGen    ) ) &&
	    ( ( lightPFit == lightPGen    && lightPBarFit == lightPBarGen ) ||
	      ( lightPFit == lightPBarGen && lightPBarFit == lightPGen    ) ) )
	  return 1;
	
	if( bBarFit == bGen &&
	    ( ( lightPFit == lightQGen    && lightPBarFit == lightQBarGen ) ||
	      ( lightPFit == lightQBarGen && lightPBarFit == lightQGen    ) ) &&
	    ( ( lightQFit == lightPGen    && lightQBarFit == lightPBarGen ) ||
	      ( lightQFit == lightPBarGen && lightQBarFit == lightPGen    ) ) )
	  return 1;

	if( bFit == bBarGen &&
	    ( ( lightQFit == lightQGen    && lightQBarFit == lightQBarGen ) ||
	      ( lightQFit == lightQBarGen && lightQBarFit == lightQGen    ) ) &&
	    ( ( lightPFit == lightPGen    && lightPBarFit == lightPBarGen ) ||
	      ( lightPFit == lightPBarGen && lightPBarFit == lightPGen    ) ) )
	  return 2;
	
	if( bBarFit == bBarGen &&
	    ( ( lightPFit == lightQGen    && lightPBarFit == lightQBarGen ) ||
	      ( lightPFit == lightQBarGen && lightPBarFit == lightQGen    ) ) &&
	    ( ( lightQFit == lightPGen    && lightQBarFit == lightPBarGen ) ||
	      ( lightQFit == lightPBarGen && lightQBarFit == lightPGen    ) ) )
	  return 2;
      }

  return 3;
}
/// histogram filling interface for reconstruction level for access with fwlite or full framework
void
KinFitImprover::fill(const TtFullHadronicEvent& tops, const edm::View<pat::Jet>& jets, const double& weight)
{
  // make sure to have a valid hypothesis on reconstruction level
  if( tops.isHypoValid("kKinFit") ){
    if( tops.isHypoValid("kGenMatch") ){

    // define indices in pat::Jet collection

      jetIndex.clear();
      jetIndex.resize(12);

      jetIndex[0]  = tops.jetLeptonCombination("kKinFit")[TtFullHadEvtPartons::LightQ   ];
      jetIndex[1]  = tops.jetLeptonCombination("kKinFit")[TtFullHadEvtPartons::LightQBar];
      jetIndex[2]  = tops.jetLeptonCombination("kKinFit")[TtFullHadEvtPartons::B        ];
      jetIndex[3]  = tops.jetLeptonCombination("kKinFit")[TtFullHadEvtPartons::BBar     ];
      jetIndex[4]  = tops.jetLeptonCombination("kKinFit")[TtFullHadEvtPartons::LightP   ];
      jetIndex[5]  = tops.jetLeptonCombination("kKinFit")[TtFullHadEvtPartons::LightPBar];

      jetIndex[6]  = tops.jetLeptonCombination("kGenMatch")[TtFullHadEvtPartons::LightQ   ];
      jetIndex[7]  = tops.jetLeptonCombination("kGenMatch")[TtFullHadEvtPartons::LightQBar];
      jetIndex[8]  = tops.jetLeptonCombination("kGenMatch")[TtFullHadEvtPartons::B        ];
      jetIndex[9]  = tops.jetLeptonCombination("kGenMatch")[TtFullHadEvtPartons::BBar     ];
      jetIndex[10] = tops.jetLeptonCombination("kGenMatch")[TtFullHadEvtPartons::LightP   ];
      jetIndex[11] = tops.jetLeptonCombination("kGenMatch")[TtFullHadEvtPartons::LightPBar];

/*
	std::cout << "-------- New Event --------" << std::endl;
	std::cout << "LightQ-mass   : " << tops.lightQ("kKinFit")->mass()    << std::endl;
	std::cout << "LightQBar-mass: " << tops.lightQBar("kKinFit")->mass() << std::endl;
	std::cout << "LightP-mass   : " << tops.lightP("kKinFit")->mass()    << std::endl;
	std::cout << "LightPBar-mass: " << tops.lightPBar("kKinFit")->mass() << std::endl;
	std::cout << "B-mass        : " << tops.b("kKinFit")->mass()         << std::endl;
	std::cout << "BBar-mass     : " << tops.bBar("kKinFit")->mass()      << std::endl;
	std::cout << "WPlus-mass    : " << tops.wPlus("kKinFit")->mass()     << std::endl;
	std::cout << "WMinus-mass   : " << tops.wMinus("kKinFit")->mass()    << std::endl;
	std::cout << "Top-mass      : " << tops.top("kKinFit")->mass()       << std::endl;
	std::cout << "TopBar-mass   : " << tops.topBar("kKinFit")->mass()    << std::endl;


	if(tops.numberOfAvailableHypos("kKinFit")>1){
	}

*/

      float dR1 = deltaR(tops.b("kKinFit")->phi(), tops.b("kKinFit")->eta(), tops.lightQ("kKinFit")->phi(), tops.lightQ("kKinFit")->eta());
      float dR2 = deltaR(tops.b("kKinFit")->phi(), tops.b("kKinFit")->eta(), tops.lightQBar("kKinFit")->phi(), tops.lightQBar("kKinFit")->eta());
      float dR3 = deltaR(tops.bBar("kKinFit")->phi(), tops.bBar("kKinFit")->eta(), tops.lightP("kKinFit")->phi(), tops.lightP("kKinFit")->eta());
      float dR4 = deltaR(tops.bBar("kKinFit")->phi(), tops.bBar("kKinFit")->eta(), tops.lightPBar("kKinFit")->phi(), tops.lightPBar("kKinFit")->eta());

      float dR1Rec;
      float dR2Rec;
      float dR3Rec;
      float dR4Rec;
      dR1Rec = deltaR(jets[jetIndex[2]].phi(), jets[jetIndex[2]].eta(), jets[jetIndex[0]].phi(), jets[jetIndex[0]].eta());
      dR2Rec = deltaR(jets[jetIndex[2]].phi(), jets[jetIndex[2]].eta(), jets[jetIndex[1]].phi(), jets[jetIndex[1]].eta());
      dR3Rec = deltaR(jets[jetIndex[3]].phi(), jets[jetIndex[3]].eta(), jets[jetIndex[4]].phi(), jets[jetIndex[4]].eta());
      dR4Rec = deltaR(jets[jetIndex[3]].phi(), jets[jetIndex[3]].eta(), jets[jetIndex[5]].phi(), jets[jetIndex[5]].eta());

      std::vector<float> dRs;
      dRs.clear();
      
      dRs.push_back(dR1);
      
      if(dR2 < dR1) dRs.insert(dRs.begin(), dR2);
      else dRs.push_back(dR2);
      
      if(dR3 < dR1) dRs.insert(dRs.begin(), dR3);
      else if(dR3 < dR2) dRs.insert(++dRs.begin(), dR3);
      else dRs.push_back(dR3);
      
      if(dR4 < dR1) dRs.insert(dRs.begin(), dR4);
      else if(dR4 < dR2) dRs.insert(++dRs.begin(), dR4);
      else if(dR4 < dR3) dRs.insert(++(++dRs.begin()), dR4);
      else dRs.push_back(dR4);

      std::vector<float> dRecs;
      dRecs.clear();
      
      dRecs.push_back(dR1Rec);
      
      if(dR2Rec < dR1Rec) dRecs.insert(dRecs.begin(), dR2Rec);
      else dRecs.push_back(dR2Rec);
      
      if(dR3Rec < dR1Rec) dRecs.insert(dRecs.begin(), dR3Rec);
      else if(dR3Rec < dR2Rec) dRecs.insert(++dRecs.begin(), dR3Rec);
      else dRecs.push_back(dR3Rec);
      
      if(dR4Rec < dR1Rec) dRecs.insert(dRecs.begin(), dR4Rec);
      else if(dR4Rec < dR2Rec) dRecs.insert(++dRecs.begin(), dR4Rec);
      else if(dR4Rec < dR3Rec) dRecs.insert(++(++dRecs.begin()), dR4Rec);
      else dRecs.push_back(dR4Rec);
      
      if( comboType() == 1 ){
	hists_.find("fitQuality")->second->Fill( 1.5 );
	hists_.find("topMass1"  )->second->Fill( tops.top("kKinFit")->mass() );
	hists_.find("allDRFit1" )->second->Fill( dR1 );
	hists_.find("allDRFit1" )->second->Fill( dR2 );
	hists_.find("allDRFit1" )->second->Fill( dR3 );
	hists_.find("allDRFit1" )->second->Fill( dR4 );
	hists_.find("minDRFit1" )->second->Fill( dRs.front() );
	hists_.find("maxDRFit1" )->second->Fill( dRs.back() );
	hists2D_.find("minMaxDRFit1")->second->Fill( dRs.back() , dRs.front() );
	
	hists_.find("allDRRec1" )->second->Fill( dR1Rec );
	hists_.find("allDRRec1" )->second->Fill( dR2Rec );
	hists_.find("allDRRec1" )->second->Fill( dR3Rec );
	hists_.find("allDRRec1" )->second->Fill( dR4Rec );
	hists_.find("minDRRec1" )->second->Fill( dRecs.front() );
	hists_.find("maxDRRec1" )->second->Fill( dRecs.back() );
	hists2D_.find("minMaxDRRec1")->second->Fill( dRecs.back() , dRecs.front() );
	hists2D_.find("minDRFitRec1")->second->Fill( dRecs.front() , dRs.front() );
	hists2D_.find("maxDRFitRec1")->second->Fill( dRecs.back() , dRs.back() );

	hists_.find("sqrt_s_fit1"  )->second->Fill( (tops.b("kKinFit")->p4()    + tops.lightQ("kKinFit")->p4() + tops.lightQBar("kKinFit")->p4() +
						     tops.bBar("kKinFit")->p4() + tops.lightP("kKinFit")->p4() + tops.lightPBar("kKinFit")->p4()).mass() );
	hists_.find("pt_final_fit1")->second->Fill( (tops.b("kKinFit")->p4()    + tops.lightQ("kKinFit")->p4() + tops.lightQBar("kKinFit")->p4() +
						     tops.bBar("kKinFit")->p4() + tops.lightP("kKinFit")->p4() + tops.lightPBar("kKinFit")->p4()).pt() );
	hists_.find("sqrt_s_rec1"  )->second->Fill( (jets[jetIndex[0]].p4() + jets[jetIndex[1]].p4() + jets[jetIndex[2]].p4() +
						     jets[jetIndex[3]].p4() + jets[jetIndex[4]].p4() + jets[jetIndex[5]].p4()).mass() );
	hists_.find("pt_final_rec1")->second->Fill( (jets[jetIndex[0]].p4() + jets[jetIndex[1]].p4() + jets[jetIndex[2]].p4() +
						     jets[jetIndex[3]].p4() + jets[jetIndex[4]].p4() + jets[jetIndex[5]].p4()).pt() );
      }
      
      if( comboType() == 2 ){
	hists_.find("fitQuality")->second->Fill( 2.5 );
	hists_.find("topMass2"  )->second->Fill( tops.top("kKinFit")->mass() );
	hists_.find("allDRFit2" )->second->Fill( dR1 );
	hists_.find("allDRFit2" )->second->Fill( dR2 );
	hists_.find("allDRFit2" )->second->Fill( dR3 );
	hists_.find("allDRFit2" )->second->Fill( dR4 );
	hists_.find("minDRFit2" )->second->Fill( dRs.front() );
	hists_.find("maxDRFit2" )->second->Fill( dRs.back() );
	hists2D_.find("minMaxDRFit2")->second->Fill( dRs.back() , dRs.front() );

	hists_.find("allDRRec2" )->second->Fill( dR1Rec );
	hists_.find("allDRRec2" )->second->Fill( dR2Rec );
	hists_.find("allDRRec2" )->second->Fill( dR3Rec );
	hists_.find("allDRRec2" )->second->Fill( dR4Rec );
	hists_.find("minDRRec2" )->second->Fill( dRecs.front() );
	hists_.find("maxDRRec2" )->second->Fill( dRecs.back() );
	hists2D_.find("minMaxDRRec2")->second->Fill( dRecs.back() , dRecs.front() );
	hists2D_.find("minDRFitRec2")->second->Fill( dRecs.front() , dRs.front() );
	hists2D_.find("maxDRFitRec2")->second->Fill( dRecs.back() , dRs.back() );

	hists_.find("sqrt_s_fit2"  )->second->Fill( (tops.b("kKinFit")->p4()    + tops.lightQ("kKinFit")->p4() + tops.lightQBar("kKinFit")->p4() +
						     tops.bBar("kKinFit")->p4() + tops.lightP("kKinFit")->p4() + tops.lightPBar("kKinFit")->p4()).mass() );
	hists_.find("pt_final_fit2")->second->Fill( (tops.b("kKinFit")->p4()    + tops.lightQ("kKinFit")->p4() + tops.lightQBar("kKinFit")->p4() +
						     tops.bBar("kKinFit")->p4() + tops.lightP("kKinFit")->p4() + tops.lightPBar("kKinFit")->p4()).pt() );
	hists_.find("sqrt_s_rec2"  )->second->Fill( (jets[jetIndex[0]].p4() + jets[jetIndex[1]].p4() + jets[jetIndex[2]].p4() +
						     jets[jetIndex[3]].p4() + jets[jetIndex[4]].p4() + jets[jetIndex[5]].p4()).mass() );
	hists_.find("pt_final_rec2")->second->Fill( (jets[jetIndex[0]].p4() + jets[jetIndex[1]].p4() + jets[jetIndex[2]].p4() +
						     jets[jetIndex[3]].p4() + jets[jetIndex[4]].p4() + jets[jetIndex[5]].p4()).pt() );
      }

      if( comboType() == 3 ){
	hists_.find("fitQuality")->second->Fill( 3.5 );
	hists_.find("topMass3"  )->second->Fill( tops.top("kKinFit")->mass() );
	hists_.find("allDRFit3" )->second->Fill( dR1 );
	hists_.find("allDRFit3" )->second->Fill( dR2 );
	hists_.find("allDRFit3" )->second->Fill( dR3 );
	hists_.find("allDRFit3" )->second->Fill( dR4 );
	hists_.find("minDRFit3" )->second->Fill( dRs.front() );
	hists_.find("maxDRFit3" )->second->Fill( dRs.back() );
	hists2D_.find("minMaxDRFit3")->second->Fill( dRs.back() , dRs.front() );

	hists_.find("allDRRec3" )->second->Fill( dR1Rec );
	hists_.find("allDRRec3" )->second->Fill( dR2Rec );
	hists_.find("allDRRec3" )->second->Fill( dR3Rec );
	hists_.find("allDRRec3" )->second->Fill( dR4Rec );
	hists_.find("minDRRec3" )->second->Fill( dRecs.front() );
	hists_.find("maxDRRec3" )->second->Fill( dRecs.back() );
	hists2D_.find("minMaxDRRec3")->second->Fill( dRecs.back() , dRecs.front() );
	hists2D_.find("minDRFitRec3")->second->Fill( dRecs.front() , dRs.front() );
	hists2D_.find("maxDRFitRec3")->second->Fill( dRecs.back() , dRs.back() );

	hists_.find("sqrt_s_fit3"  )->second->Fill( (tops.b("kKinFit")->p4()    + tops.lightQ("kKinFit")->p4() + tops.lightQBar("kKinFit")->p4() +
						     tops.bBar("kKinFit")->p4() + tops.lightP("kKinFit")->p4() + tops.lightPBar("kKinFit")->p4()).mass() );
	hists_.find("pt_final_fit3")->second->Fill( (tops.b("kKinFit")->p4()    + tops.lightQ("kKinFit")->p4() + tops.lightQBar("kKinFit")->p4() +
						     tops.bBar("kKinFit")->p4() + tops.lightP("kKinFit")->p4() + tops.lightPBar("kKinFit")->p4()).pt() );
	hists_.find("sqrt_s_rec3"  )->second->Fill( (jets[jetIndex[0]].p4() + jets[jetIndex[1]].p4() + jets[jetIndex[2]].p4() +
						     jets[jetIndex[3]].p4() + jets[jetIndex[4]].p4() + jets[jetIndex[5]].p4()).mass() );
	hists_.find("pt_final_rec3")->second->Fill( (jets[jetIndex[0]].p4() + jets[jetIndex[1]].p4() + jets[jetIndex[2]].p4() +
						     jets[jetIndex[3]].p4() + jets[jetIndex[4]].p4() + jets[jetIndex[5]].p4()).pt() );
      }
      dRs.clear();
      dRecs.clear();
    }
  }
  jetIndex.clear();
}
