#ifndef MuonAnalyzer_h
#define MuonAnalyzer_h

#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/Candidate/interface/Particle.h>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "TopAnalysis/TopAnalyzer/interface/BJetKin.h"
#include "TopAnalysis/TopAnalyzer/interface/BJetTag.h"

using namespace std;

class BJetAnalyzer : public edm::EDAnalyzer {

 public:

  explicit BJetAnalyzer(const edm::ParameterSet&);
  ~BJetAnalyzer();
  
 private:

  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  edm::InputTag jets_;
  
  std::string hist_;
  int onlyLeading_;
  bool doKin_, doES_, doBtag_;

  BJetKin *kin_uds_, *kin_c_, *kin_b_, *kin_g_;
  BJetTag *btg_uds_, *btg_c_, *btg_b_, *btg_g_;

  typedef std::vector<pat::Jet> PatJetCollection; 
  
  TH1F *nFlav_;
};

#endif
