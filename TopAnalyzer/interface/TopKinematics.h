#ifndef TopKinematics_h
#define TopKinematics_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "TopAnalysis/TopAnalyzer/interface/SingleObject.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiLeptonicEvent.h"
#include "TMath.h"

/**
   \class   TopKinematics TopKinematics.h "TopAnalysis/TopAnalyzer/interface/TopKinematics.h"

   \brief   Derived class to analyze the kinematics of top quarks on reconstruction and generator level

   The structure keeps histograms for the kinematics of top quarks and ttbar pairs. These histograms 
   can be filled from the TtSemiLeptonicEvent class (with an extra option to require matching for 
   purity & stability studies) or from the TtGenEvent class. The class is derived from the 
   SingleObject<Collection> interface, which makes it usable in full framework or fwlite. 
*/
double pi = TMath::Pi();

namespace CrossSection{
 // binning for top & ttbar cross section histograms
 double topPt      [] = {             0. ,    10. ,    35. ,    65. ,   100. ,   135. ,   170. ,   215. ,   300. ,   400.  };
 double topY       [] = {   -4. ,    -3. ,  -1.92 ,   -1.2 ,  -0.48 ,     0. ,   0.48 ,    1.2 ,   1.92 ,     3. ,     4.  };
 double topPhi     [] = {                                       -pi ,   -2.0 ,   -1.0 ,     0. ,     1. ,     2. ,     pi  };
 double ttbarMass  [] = {                                      300. ,   345. ,   400. ,   490. ,   615. ,   790. ,  1000.  };
 double ttbarPt    [] = {                      0. ,    10. ,    20. ,    30. ,    40. ,    60. ,    80. ,   130. ,   200.  };
 double ttbarY     [] = {   -5. ,    -4. ,    -3. ,   -2.3 ,   -1.6 ,     0. ,    1.6 ,    2.3 ,     3. ,     4. ,     5.  };
 double ttbarPhi   [] = {                                       -pi ,   -2.0 ,   -1.0 ,     0. ,     1. ,     2. ,     pi  };
 double ttbarDelPhi[] = {                                       -pi ,   -2.0 ,   -1.0 ,     0. ,     1. ,     2. ,     pi  };
 double topWAngle  [] = {                                        0. , 1./6*pi, 1./3*pi, 1./2*pi, 2./3*pi, 5./6*pi,     pi  };
 double ttbarHT    [] = {                                       50. ,   150. ,   225. ,   300. ,   400. ,   600. ,  1000.  };
}

class TopKinematics : public SingleObject<TtSemiLeptonicEvent> {
  
 public:
  /// default constructor for generator level analysis in fw lite
  explicit TopKinematics();
  /// default constructor for reco level analysis in fw lite
  explicit TopKinematics(const std::string& hypoKey, const std::string& lepton, const bool& matchForStabilityAndPurity);
  /// default constructor for fwfull
  explicit TopKinematics(const edm::ParameterSet& configFile);
  /// default destructor
  ~TopKinematics(){};

  /**
     The following functions have to be implemented for any class
     derived from SingleObject<Collection>. The fill function is
     overloaded to serve both for the TtGenEvent class directoy or 
     the TtSemiLeptonicEvent class.
  **/
  /// histogramm booking for fwlite 
  void book();
  /// histogramm booking for fw
  void book(edm::Service<TFileService>& fileService);
  /// histogram filling interface for reconstruction level for access with fw or fwlite
  void fill(const TtSemiLeptonicEvent& tops, const double& weight=1.);
  /// histogram filling interface for generator level for access with fw or fwlite 
  void fill(const TtGenEvent& tops, const double& weight=1.);
  // histogram filling and in addition save information for event identification
  // double definition to avoid problems with different module inputs
  void fill2(const TtGenEvent& tops, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& weight=1.){
    fillValue("runNumber", runNumber, weight);
    fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
    fillValue("eventNumber", eventNumber, weight);
    fill(tops, weight);
  }
  void fill2(const TtSemiLeptonicEvent& tops, const double& runNumber, const double& luminosityBlockNumber, const double& eventNumber, const double& weight=1.){
    fillValue("runNumber", runNumber, weight);
    fillValue("luminosityBlockNumber", luminosityBlockNumber, weight);
    fillValue("eventNumber", eventNumber, weight);
    fill(tops, weight);
  }
  /// everything which needs to be done after the event loop
  void process(){};

 private:
  /**
     Helper functions for special histogram management
  **/
  /// histogram filling for candidates topA and topB
  void fill(const reco::Candidate* leptonicTop, const reco::Candidate* hadronicTop, const reco::Candidate* leptonicW, const reco::Candidate* hadronicW, double HT, const double& weight=1.);
  /// histogram filling for candidate topA and topB (for stability and purity calculation)
  void fill(const reco::Candidate* leptonicTopRec, const reco::Candidate* leptonicTopGen, const reco::Candidate* hadronicTopRec, const reco::Candidate* hadronicTopGen, const reco::Candidate* leptonicWRec, const reco::Candidate* leptonicWGen, const reco::Candidate* hadronicWRec, const reco::Candidate* hadronicWGen, double HTrec, double HTgen, const double& weight=1.);
  /// helper function for determining stability and purity
  void match(const std::string& histo, const double& genValue, const double& recValue, const double& weight);
  ///  helper functions to fill 1D angle histos
  void fillAngles(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& muon    ,
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
		  const double& weight=1.);
  ///  helper functions to fill event shape variables
  void fillEventShapes(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& muon    ,
		       const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
		       const double& weight, const bool useMu, const bool useNu);
  ///  helper functions to fill final state objects
  void fillFinalStateObjects(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& hadB    ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& q       ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& qbar    ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepB    ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& lepton  ,
			     const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& neutrino,
			     const double& weight);
  double getDecayChannel(const TtGenEvent& tops){
    double eventType=-4;
    // identify event by the following number coding:
    // -4 : undefined
    // -3 : non ttbar MC
    // -2 : dileptonic unknown    ttbar MC
    // -1 : semileptonic unknown  ttbar MC
    // 01 : semileptonic electron ttbar MC
    // 02 : semileptonic muon     ttbar MC
    // 03 : semileptonic tau      ttbar MC
    // 11 : dileptonic e   e      ttbar MC
    // 12 : dileptonic e   mu     ttbar MC
    // 13 : dileptonic e   tau    ttbar MC
    // 14 : dileptonic mu  mu     ttbar MC
    // 15 : dileptonic mu  tau    ttbar MC
    // 16 : dileptonic tau tau    ttbar MC
    // 20 : fully hadronic        ttbar MC
    if( !tops.isTtBar()        ) eventType = -3;
    else if(  tops.isFullHadronic() ) eventType = 20;
    else if(  tops.isSemiLeptonic() ) {
      switch( tops.semiLeptonicChannel() ) {
      case WDecay::kElec : eventType = 1; break;
      case WDecay::kMuon : eventType = 2; break;
      case WDecay::kTau  : eventType = 3; break;
      default            : eventType = -1; break;
      }
    }
    else if( tops.isFullLeptonic() ) {
      if     (  tops.fullLeptonicChannel().first  == WDecay::kElec  && 
		tops.fullLeptonicChannel().second == WDecay::kElec    )    eventType = 11;
      else if( ( tops.fullLeptonicChannel().first  == WDecay::kElec && 
		 tops.fullLeptonicChannel().second == WDecay::kMuon   ) ||
	       ( tops.fullLeptonicChannel().first  == WDecay::kMuon && 
		 tops.fullLeptonicChannel().second == WDecay::kElec   )  ) eventType = 12;
      else if( ( tops.fullLeptonicChannel().first  == WDecay::kElec && 
		 tops.fullLeptonicChannel().second == WDecay::kTau    ) ||
	       ( tops.fullLeptonicChannel().first  == WDecay::kTau && 
		 tops.fullLeptonicChannel().second == WDecay::kElec   )  ) eventType = 13;
      else if(  tops.fullLeptonicChannel().first  == WDecay::kMuon && 
		tops.fullLeptonicChannel().second == WDecay::kMuon  )      eventType = 14;
      else if( ( tops.fullLeptonicChannel().first  == WDecay::kMuon && 
		 tops.fullLeptonicChannel().second == WDecay::kTau    ) ||
	       ( tops.fullLeptonicChannel().first  == WDecay::kTau  &&
		 tops.fullLeptonicChannel().second  == WDecay::kMuon) )    eventType = 15;
      else if(  tops.fullLeptonicChannel().first  == WDecay::kTau &&
		tops.fullLeptonicChannel().second == WDecay::kTau   )      eventType = 16;
      else eventType = -2;
    }
    return eventType;
  }
  
 private:
  /// class key of hypothesis
  std::string hypoKey_;
  /// lepton flavour in semileptonic event: muon (default) or electron
  std::string lepton_;
  /// bool to decide whether to use a ttree 
  bool useTree_;
  /// apply matching for stability and purity or not
  bool matchForStabilityAndPurity_;
  /// choose whether you want to destinguish between top/antitop instead of leptonic/hadronic top
  bool ttbarInsteadOfLepHadTop_;
  /// number of jets considered in Kinematic fit (needed for the classification of permutations)
  int maxNJets;
  /// histogram container for correlation plots
  std::map<std::string, TH2*> corrs_;
};

#endif
