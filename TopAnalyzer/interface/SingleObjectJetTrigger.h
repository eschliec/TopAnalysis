#ifndef SingleObjectJetTrigger_h
#define SingleObjectJetTrigger_h

#include <string>
#include <cstring>
#include <vector>
#include <map>

#include "TH1.h"
#include "TH2.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//#include "FWCore/Framework/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

/**
   \class   SingleObjectJetTrigger SingleObjectJetTrigger.h "TopAnalysis/TopAnalyzer/interface/SingleObjectJetTrigger.h"

   \brief   EDAnalyzer for jet trigger

   The structure keeps control histograms to analyze jet triggers. These histograms are to be
   filled from the trigger::TriggerEvent, edm:TriggerResults and from edm::View<pat::Jet>'s.
*/

class SingleObjectJetTrigger : public edm::EDAnalyzer {

 public:
  /// default constructor
  explicit SingleObjectJetTrigger(const edm::ParameterSet& configFile);
  /// default destructor
  ~SingleObjectJetTrigger(){};

 private:

  /// initiate histograms
  virtual void beginJob(const edm::EventSetup&);
  /// analyze triggers and fill histograms
  virtual void analyze(const edm::Event& event, const edm::EventSetup&);
  /// empty
  virtual void endJob();
        
  /// trigger results input collection	
  edm::InputTag triggerResults_;
  /// trigger summary input collection	
  edm::InputTag triggerSummary_;
  /// pat trigger event
  edm::InputTag patTriggerEvent_;
  /// pat trigger (filter, path, object)
  edm::InputTag patTrigger_;
  /// jets input collection	
  edm::InputTag jets_;
  /// trigger matched jets input collection	
  edm::InputTag triggerMatchedJets_;
  /// trigger to get analyzed
  edm::InputTag analyzedTrigger_;

  /// container to keep all histogramms
  std::map<std::string, TH1*> hists_;
  std::map<std::string, TH2*> hists2D_;
};

#endif
