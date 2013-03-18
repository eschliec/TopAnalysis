#ifndef SemiLepBjetAnalyzer_h
#define SemiLepBjetAnalyzer_h

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"


/**
   \class   SemiLepBjetAnalyzer.cc SemiLepBjetAnalyzer.h "TopAnalysis/TopAnalyzer/interface/SemiLepBjetAnalyzer.h"

   \brief   Derived class to analyze the kinematics of b-jets on reconstruction level (kinfit output) and generator level (genjets from)
   
   The structure keeps histograms for the kinematics of bjets. These histograms are filled either from the TtSemiLeptonicEvent class (reconstruction level) or from the specified genjet collection. In the latter case the identification as done in the GenLevelBJetProducer (from TopAnalysis/TopUtils) is used. Furthermore, also the correlation plot between reconstructed and generated quantity is filled. The lepton charge is used to identify the b and bbar jets from the kinematic fit output. Both (b and bbar) are always drawn in the same histogram.   
*/

class SemiLepBjetAnalyzer : public edm::EDAnalyzer {

 public:

  explicit SemiLepBjetAnalyzer(const edm::ParameterSet&);
  ~SemiLepBjetAnalyzer(){};
  
 private:

  // ---
  //    functions
  // ---  
  // default (framework)
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  // helper functions to get b and bbar jet 
  // a) from kinfit hypothesis
  std::pair<const reco::Candidate*,const reco::Candidate*> getbJets(const edm::Handle<TtSemiLeptonicEvent> semiLepEvt, const std::string hypoKey_);
  // b) from genJet collection
  const reco::GenJet* getJetFromCollection(const reco::GenJetCollection& genJets, int IX);
  // ---
  //    inputs
  // ---
  // kinfit output
  edm::InputTag semiLepEvt_;
  // class key of kinfit hypothesis
  std::string hypoKey_;
  // genJet collection 
  edm::InputTag genJets_;
  // pure b gen jet collection?
  bool bJetCollection_;
  // output manager
  int verbose;
  // event weight
  edm::InputTag weight_;
  // create gen plots?
  bool genPlots_;
  // create rec plots?
  bool recPlots_;
  // from gen level b-jet identification  
  edm::InputTag bHadJetIdx_, antibHadJetIdx_;
  // take minimum dR as criteria to assign rec and gen level b-jets to each other?
  // this is useful for simultaneous measurements of b-jet and anti-bjet 
  // where only the resolution and not the association matters
  bool useClosestDrBs_;

  // ---
  //    define Tree for event content
  // ---
  TTree * tree;
  bool useTree_;
  float weight;
  unsigned int runNumber, luminosityBlockNumber, eventNumber;
  float valueBqPtRec;
  float valueBqPtGen;
  float valueLeadBqPtRec;
  float valueLeadBqPtGen;
  float valueBqEtaRec;
  float valueBqEtaGen;
  float valueBqYRec;
  float valueBqYGen;
  float valueBbarqPtRec;
  float valueBbarqPtGen;
  float valueBbarqEtaRec;
  float valueBbarqEtaGen;
  float valueBbarqYRec;
  float valueBbarqYGen;
  float valueBbbarPtRec;
  float valueBbbarPtGen;
  float valueBbbarEtaRec;
  float valueBbbarEtaGen;
  float valueBbbarYRec;
  float valueBbbarYGen;
  float valueBbbarMassRec;
  float valueBbbarMassGen;
  bool bbSwapBetter;

  // ---
  //    histos
  // ---
  TH1F* bqPtRec;
  TH1F* bqPtGen;

  TH1F* leadbqPtRec;
  TH1F* leadbqPtGen;

  TH1F* bqEtaRec;
  TH1F* bqEtaGen;

  TH1F* bqYRec;
  TH1F* bqYGen;

  TH1F* bbbarPtRec;
  TH1F* bbbarPtGen;

  TH1F* bbbarEtaRec;
  TH1F* bbbarEtaGen;

  TH1F* bbbarYRec;
  TH1F* bbbarYGen;

  TH1F* bbbarMassRec;
  TH1F* bbbarMassGen;

  TH2F* bqPt_;
  TH2F* leadbqPt_;
  TH2F* bqEta_;
  TH2F* bqY_;
  TH2F* bqPtClosestPt_;
  TH2F* bqEtaClosestPt_;
  TH2F* bqYClosestPt_;
  TH2F* bbbarPt_;
  TH2F* bbbarEta_;
  TH2F* bbbarY_;
  TH2F* bbbarMass_;
};

#endif
