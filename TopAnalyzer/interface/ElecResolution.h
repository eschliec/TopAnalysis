#ifndef ElecResolution_h
#define ElecResolution_h

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
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "DataFormats/PatCandidates/interface/Electron.h"


class ElecResolution{

 public:

  explicit ElecResolution(std::vector<double>, double);
  explicit ElecResolution(const edm::ParameterSet&);
  ~ElecResolution(){};

  void book();
  void book(edm::Service<TFileService>&);
  void book(edm::Service<TFileService>&, ofstream&);
  void fill(const edm::Event&, const std::vector<pat::Electron>&, const double&);
  void fill(const std::vector<pat::Electron>&, const double&);
  void norm(){};
  void write(const char*, const char*);

 private:

  // additional evt content/steerings
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
