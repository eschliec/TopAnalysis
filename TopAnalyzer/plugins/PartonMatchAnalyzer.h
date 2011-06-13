#ifndef PartonMatchAnalyzer_h
#define PartonMatchAnalyzer_h

#include <memory>
#include <string>
#include <vector>

#include "TH1F.h"
#include "TH2F.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TopAnalysis/TopUtils/interface/NameScheme.h"


class PartonMatchAnalyzer : public edm::EDAnalyzer {

 public:

  explicit PartonMatchAnalyzer(const edm::ParameterSet&);
  ~PartonMatchAnalyzer(){};
  
 private:

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:

  std::vector<edm::InputTag> matches_; // InputTag of matches
  std::vector<int> quarks_, jets_;     // quarks/jets in consideration (if empty any)
  edm::InputTag sumDR_, sumDPt_;
  std::string hist_;
  unsigned int ref_;
  unsigned int normQrks_;
  unsigned int normJets_;

  TH2F* corr_;
  TH1F* monDR_;
  TH1F* monDPt_;
  TH1F* algJets_;
  TH1F* algQrks_;
  TH1F* effJets_;
  TH1F* effQrks_;
  TH1F* purJets_;
  TH1F* purQrks_;
};

#endif
