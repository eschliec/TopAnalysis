#ifndef MuonId_h
#define MuonId_h

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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"


class MuonId{

 public:

  explicit MuonId();
  explicit MuonId(const edm::ParameterSet&);
  ~MuonId(){
    // free allocated space
    if(fwLite_){
      delete muComp_; 
      delete muEm_; 
      delete muEmS9_; 
      delete muHad_; 
      delete muHadS9_;
      delete muHo_; 
      delete muHoS9_; 
    }
  };

  void book();
  void book(edm::Service<TFileService>&);
  void book(edm::Service<TFileService>&, ofstream&);
  void fill(const edm::Event&, const std::vector<pat::Muon>&, const double&);
  void fill(const std::vector<pat::Muon>&, const double&);
  void norm(){};
  void write(TFile&, const char*);
  
 private:

  bool fwLite_;

  TH1F *muComp_;
  TH1F *muEm_;
  TH1F *muEmS9_;
  TH1F *muHad_;
  TH1F *muHadS9_;
  TH1F *muHo_;
  TH1F *muHoS9_;

  TH1F *emEt03_;
  TH1F *hadEt03_;
  TH1F *hoEt03_;   
  TH1F *nTracks03_; 
  TH1F *sumPt03_; 
  
  TH1F *emEt05_;
  TH1F *hadEt05_;
  TH1F *hoEt05_;  
  TH1F *nTracks05_; 
  TH1F *sumPt05_; 
};

#endif
