#ifndef MuonResolution_h
#define MuonResolution_h

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TH1F.h"
#include "TFile.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"


class MuonResolution{

 public:

  explicit MuonResolution();
  explicit MuonResolution(const edm::ParameterSet&);
  ~MuonResolution(){
    // free allocated space
    if(fwLite_){
      delete calPt_;
      delete resPt_;
      delete calEta_;
      delete resEta_;
      delete calPhi_;
      delete resPhi_;
    }
  };
  
  void book();
  void book(edm::Service<TFileService>&);
  void book(edm::Service<TFileService>&, ofstream&);
  void fill(const edm::Event&, const std::vector<pat::Muon>&, const double&);
  void fill(const std::vector<pat::Muon>&, const double&);
  void norm(){};
  void write(TFile&, const char*);
  void setMatchDR(double match) {matchDR_=match;};
  void setBinsPt(std::vector<double> bins) {binsPt_=bins;};

 private:

  bool fwLite_;
  std::vector<double> binsPt_;
  std::vector<double> binsEta_;
  std::vector<double> binsPhi_;
  double matchDR_;

  TH1F *calPt_, *resPt_, *calEta_, *resEta_, *calPhi_, *resPhi_;
  std::vector<TH1F*> relPt_;
  std::vector<TH1F*> Eta_;
  std::vector<TH1F*> Phi_;
};

#endif
