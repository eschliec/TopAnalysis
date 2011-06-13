#ifndef ElecId_h
#define ElecId_h

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

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"

#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


class ElecId{

 public:

  explicit ElecId();
  explicit ElecId(const edm::ParameterSet&);
  ~ElecId(){};

  void book();
  void book(edm::Service<TFileService>&);
  void book(edm::Service<TFileService>&, ofstream&);
  void fill(const edm::Event&, const std::vector<pat::Electron>&, const double&);
  void fill(const std::vector<pat::Electron>&, const double&);
  void norm(){};
  void write(const char*, const char*);
  
 private:

  // additional evt content/steerings  
  edm::InputTag bshp_, eshp_;

  TH1F *eops_;
  TH1F *eopb_;
  TH1F *nocl_;
  TH1F *hoem_;
  TH1F *deta_;
  TH1F *dphi_;  
  TH1F *drtk_;
  TH1F *tdpt_;
  TH1F *tdeta_;
  TH1F *tdphi_;
  TH1F *nohit_;
  TH1F *nvhit_;
  TH1F *chi2_;
  TH1F *s1os9_;
};

#endif
