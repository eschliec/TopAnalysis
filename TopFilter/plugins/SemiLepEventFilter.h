#ifndef SemiLepEventFilter_h
#define SemiLepEventFilter_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TopAnalysis/TopFilter/interface/PtFilter.h"
#include "TopAnalysis/TopFilter/interface/EtaFilter.h"
#include "TopAnalysis/TopFilter/interface/DistanceFilter.h"
#include "TopAnalysis/TopFilter/interface/IsolationFilter.h"


template <typename Lep> 
class SemiLepEventFilter : public edm::EDFilter {

 public:
  
  typedef DistanceFilter<Lep> LepDistanceFilter;
  typedef IsolationFilter<Lep> LepIsolationFilter;

  explicit SemiLepEventFilter(const edm::ParameterSet&);
  ~SemiLepEventFilter(){};
  
 private:
  
  virtual void beginJob(const edm::EventSetup&) ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
 private:

  edm::InputTag wgt_;
  std::vector<edm::InputTag> jets_;
  std::vector<edm::InputTag> leps_;

  /// switches 
  bool doJetEta_;
  bool doJetPt_;
  bool doLepEta_;
  bool doLepPt_;
  bool doTrkIso_;
  bool doCalIso_;
  bool doDist_;

  /// filters
  PtFilter  jetPt_;
  PtFilter  lepPt_;
  EtaFilter jetEta_; 
  EtaFilter lepEta_; 
  LepDistanceFilter jetDist_;
  LepIsolationFilter trkIso_;
  LepIsolationFilter calIso_;
};

template <typename Lep> 
SemiLepEventFilter<Lep>::SemiLepEventFilter(const edm::ParameterSet& cfg):
  wgt_     ( cfg.template getParameter<edm::InputTag>( "weight" ) ),
  jets_    ( cfg.template getParameter<std::vector<edm::InputTag> >("jets" ) ),
  leps_    ( cfg.template getParameter<std::vector<edm::InputTag> >("leptons" ) ),
  jetPt_   ( cfg.template getParameter<edm::ParameterSet>("jetPtFilter"  ) ), 
  lepPt_   ( cfg.template getParameter<edm::ParameterSet>("lepPtFilter"  ) ),
  jetEta_  ( cfg.template getParameter<edm::ParameterSet>("jetEtaFilter" ) ), 
  lepEta_  ( cfg.template getParameter<edm::ParameterSet>("lepEtaFilter" ) ), 
  jetDist_ ( cfg.template getParameter<edm::ParameterSet>("jetDistFilter") ),
  trkIso_  ( cfg.template getParameter<edm::ParameterSet>("trkIsoFilter" ) ), 
  calIso_  ( cfg.template getParameter<edm::ParameterSet>("calIsoFilter" ) ) 
{
  doJetEta_= cfg.template getParameter<bool>("jetEta" );
  doJetPt_ = cfg.template getParameter<bool>("jetPt"  );
  doLepEta_= cfg.template getParameter<bool>("lepEta" );
  doLepPt_ = cfg.template getParameter<bool>("lepPt"  );
  doTrkIso_= cfg.template getParameter<bool>("trkIso" );
  doCalIso_= cfg.template getParameter<bool>("calIso" );
  doDist_  = cfg.template getParameter<bool>("jetDist");
}

template <typename Lep> 
bool SemiLepEventFilter<Lep>::filter(edm::Event& evt, const edm::EventSetup& setup)
{
  // get event weight
  edm::Handle<double> wgt;
  evt.getByLabel(wgt_, wgt);
  double weight = *wgt;

  // receive input tags for lepton kinematic
  std::vector<edm::View<reco::Candidate> > kinLeps;
  for(std::vector<edm::InputTag>::const_iterator tag = leps_.begin(); 
      tag!=leps_.end(); ++tag){
    edm::Handle<edm::View<reco::Candidate> > src;
    evt.getByLabel(*tag, src);
    kinLeps.push_back(*src);
  }

  // receive input tags for lepton isolation
  std::vector<Lep> isoLeps;
  for(std::vector<edm::InputTag>::const_iterator tag = leps_.begin(); 
      tag!=leps_.end(); ++tag){
    edm::Handle<Lep> src;
    evt.getByLabel(*tag, src);
    isoLeps.push_back(*src);
  }

  // receive input tags for jet kinematic
  std::vector<edm::View<reco::Candidate> > kinJets;
  for(std::vector<edm::InputTag>::const_iterator tag = jets_.begin(); 
      tag!=jets_.end(); ++tag){
    edm::Handle<edm::View<reco::Candidate> > src;
    evt.getByLabel(*tag, src);
    kinJets.push_back(*src);
  }

  // do the event selection
  bool passed=true;
  if(passed && doLepEta_) passed=lepEta_ (evt, kinLeps, weight);
  if(passed && doLepPt_ ) passed=lepPt_  (evt, kinLeps, weight);
  if(passed && doTrkIso_) passed=trkIso_ (evt, isoLeps, weight);
  if(passed && doCalIso_) passed=calIso_ (evt, isoLeps, weight);
  if(passed && doDist_  ) passed=jetDist_(evt, isoLeps, weight);
  if(passed && doJetEta_) passed=jetEta_ (evt, kinJets, weight);
  if(passed && doJetPt_ ) passed=jetPt_  (evt, kinJets, weight);
  return passed;
}

template <typename Lep> 
void SemiLepEventFilter<Lep>::beginJob(const edm::EventSetup&)
{
}

template <typename Lep> 
void SemiLepEventFilter<Lep>::endJob()
{
  using std::cout;
  using std::endl;

  cout << "----------------------------------"  << endl;
  cout << " applyLepEta  [" << doLepEta_ << "]" << endl;
  cout << " applyLepPt   [" << doLepPt_  << "]" << endl;
  cout << " applyTrkIso  [" << doTrkIso_ << "]" << endl;
  cout << " applyCalIso  [" << doCalIso_ << "]" << endl;
  cout << " applyDist    [" << doDist_   << "]" << endl;
  cout << " applyJetEta  [" << doJetEta_ << "]" << endl;
  cout << " applyJetPt   [" << doJetPt_  << "]" << endl;
  cout << "----------------------------------"  << endl;

  if( !doLepEta_ && !doLepPt_ && !doTrkIso_ && !doCalIso_ && !doDist_
      && !doJetEta_ && !doJetPt_ )
    cout << "no cuts were applied..." << endl;
  else {
    if(doLepEta_ ) lepEta_ .summarize();
    if(doLepPt_  ) lepPt_  .summarize();
    if(doTrkIso_ ) trkIso_ .summarize();
    if(doCalIso_ ) calIso_ .summarize();
    if(doDist_   ) jetDist_.summarize();
    if(doJetEta_ ) jetEta_ .summarize();
    if(doJetPt_  ) jetPt_  .summarize();
  }
}

#endif
