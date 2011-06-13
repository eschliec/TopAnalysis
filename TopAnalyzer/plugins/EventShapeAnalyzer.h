#ifndef EventShapeAnalyzer_h
#define EventShapeAnalyzer_h

#include <memory>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TMath.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopQuarkAnalysis/TopTools/interface/EventShapeVariables.h"


template <typename LeptonType>
class EventShapeAnalyzer : public edm::EDAnalyzer {

 public:

  explicit EventShapeAnalyzer(const edm::ParameterSet&);
  ~EventShapeAnalyzer(){};
  
 private:

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

 private:

  edm::InputTag jets_;
  edm::InputTag leps_;
  std::string hist_;

  unsigned int nJets_;
  double minPt_, maxEta_;
  bool doMatching_;
  double maxDRLep_, maxDRJet_;

  TH1F *sph_;   
  TH1F *apl_;
  TH1F *cir_;
  TH1F *iso_;
};

template <typename LeptonType>
EventShapeAnalyzer<LeptonType>::EventShapeAnalyzer(const edm::ParameterSet& cfg):
  jets_( cfg.getParameter<edm::InputTag>( "jets" ) ),
  leps_( cfg.getParameter<edm::InputTag>( "lepton" ) ),
  hist_( cfg.getParameter<std::string>( "hist" ) ),
  nJets_ ( cfg.getParameter<unsigned int>( "nJets" ) ),
  minPt_ ( cfg.getParameter<double>( "minPt" ) ),
  maxEta_( cfg.getParameter<double>( "maxEta" ) ),
  doMatching_( cfg.getParameter<bool>( "doMatching" ) ),
  maxDRLep_( cfg.getParameter<double>( "maxDeltaRLepton" ) ),
  maxDRJet_( cfg.getParameter<double>( "maxDeltaRJet" ) )
{
}

template <typename LeptonType>
void EventShapeAnalyzer<LeptonType>::analyze(const edm::Event& evt, const edm::EventSetup& setup)
{
  edm::Handle<std::vector<pat::Jet> > jets; 
  evt.getByLabel(jets_, jets);
  
  edm::Handle<LeptonType> leps; 
  evt.getByLabel(leps_, leps);

  // prepare vectors for event shape variable calculation 
  std::vector<TVector3> vecs;
  vecs.clear();

  unsigned int idx=0;
  for(std::vector<pat::Jet>::const_iterator jet=jets->begin(); 
      jet!=jets->end(); ++jet) {
    if(jet->pt()<minPt_) 
      continue;
    if(fabs(jet->eta())>maxEta_)
      continue;
    ++idx;
    if( idx>nJets_ ) break;
    TVector3 vec(jet->px(),jet->py(),jet->pz());
    vecs.push_back( vec );
  }
  if(vecs.size()<nJets_){
    edm::LogWarning ( "TooViewJetsFound" ) << "too view jets found within requirement";
    return;
  }

  typename LeptonType::const_iterator lep=leps->begin();
  if(lep!=leps->end()){
    TVector3 vec(lep->px(),lep->py(),lep->pz());
    vecs.push_back( vec );
  }
  else{
    edm::LogWarning ( "NoLeptonFound" ) << "no required lepton found in event";
    return;
  }

  EventShapeVariables esv;
  sph_->Fill( esv.sphericity ( vecs ) );
  apl_->Fill( esv.aplanarity ( vecs ) );
  cir_->Fill( esv.circularity( vecs ) );
  iso_->Fill( esv.isotropy   ( vecs ) ); 
}

template <typename LeptonType>
void EventShapeAnalyzer<LeptonType>::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if( !fs )
    throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );
  
  NameScheme shp("shp");
  if(hist_.empty()){
    sph_ = fs->make<TH1F>(shp.name(      "sph"), shp.name("sph"), 100, 0., 1.);
    apl_ = fs->make<TH1F>(shp.name(      "apl"), shp.name("apl"), 100, 0., 1.);
    cir_ = fs->make<TH1F>(shp.name(      "cir"), shp.name("cir"), 100, 0., 1.);
    iso_ = fs->make<TH1F>(shp.name(      "iso"), shp.name("iso"), 100, 0., 1.);
  }
  else{
    ofstream hist(hist_.c_str(), std::ios::out);
    sph_ = fs->make<TH1F>(shp.name(hist, "sph"), shp.name("sph"), 100, 0., 1.);
    apl_ = fs->make<TH1F>(shp.name(hist, "apl"), shp.name("apl"), 100, 0., 1.);
    cir_ = fs->make<TH1F>(shp.name(hist, "cir"), shp.name("cir"), 100, 0., 1.);
    iso_ = fs->make<TH1F>(shp.name(hist, "iso"), shp.name("iso"), 100, 0., 1.);
  }
}

#endif
