#ifndef MuonKinematic_h
#define MuonKinematic_h

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"


class MuonKinematic{

 public:

  explicit MuonKinematic();
  explicit MuonKinematic(const edm::ParameterSet&);
  ~MuonKinematic(){
    // free allocated space
    if(fwLite_){
      delete en_;
      delete pt_;
      delete eta_;
      delete phi_;
      delete isoJet_;
      delete isoJet5_;
      delete isoJet10_;
      delete isoJet15_;
      delete isoJet20_;
      delete isoJet25_;
      delete isoJet30_;
      delete isoJet35_;
      delete isoJet40_;
      delete isoTrkPt_;
      delete isoCalPt_;
      delete isoRelPt_;
      delete isoRelComb_;
      delete isoTrkN_;
      delete isoEcalN_;
      delete isoHcalN_;
      delete dREcalPt_;
      delete dREcalN_;
      delete dRHcalPt_;
      delete dRHcalN_;
      delete dRTrkPt_;
      delete dRTrkN_;
      delete ptVsTrkIso_;
      delete ptVsCalIso_;
    }
  };

  void book();
  void book(edm::Service<TFileService>&);
  void book(edm::Service<TFileService>&, ofstream&);
  void fill(const edm::Event&, const std::vector<pat::Muon>&, const double&);
  void fill(const std::vector<pat::Jet>& jets,
	    const std::vector<pat::Muon>& muons, const double&);
  void norm(){
    double totalNumberOfMuons = en_->GetEntries();
    dREcalPt_->Scale( 1./totalNumberOfMuons );
    dREcalN_ ->Scale( 1./totalNumberOfMuons );
    dRHcalPt_->Scale( 1./totalNumberOfMuons );
    dRHcalN_ ->Scale( 1./totalNumberOfMuons );
    dRTrkPt_ ->Scale( 1./totalNumberOfMuons );
    dRTrkN_  ->Scale( 1./totalNumberOfMuons );
  };
  void write(TFile&, const char*);

 private:

  bool fwLite_;
  edm::InputTag jets_;

  TH1F *en_;
  TH1F *pt_;
  TH1F *eta_;
  TH1F *phi_;

  TH1F *isoJet_;
  TH1F *isoJet5_;
  TH1F *isoJet10_;
  TH1F *isoJet15_;
  TH1F *isoJet20_;
  TH1F *isoJet25_;
  TH1F *isoJet30_;
  TH1F *isoJet35_;
  TH1F *isoJet40_;
  TH1F *isoTrkPt_;
  TH1F *isoCalPt_;
  TH1F *isoRelPt_;
  TH1F *isoRelComb_;
  TH1F *isoTrkN_;
  TH1F *isoEcalN_;
  TH1F *isoHcalN_;
  TH1F *dREcalPt_;
  TH1F *dREcalN_;
  TH1F *dRHcalPt_;
  TH1F *dRHcalN_;
  TH1F *dRTrkPt_;
  TH1F *dRTrkN_;

  TH2F *ptVsTrkIso_;
  TH2F *ptVsCalIso_;
};

#endif
