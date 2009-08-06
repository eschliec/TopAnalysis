#ifndef FullLepHypothesesFilter_h
#define FullLepHypothesesFilter_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "AnalysisDataFormats/TopObjects/interface/TtFullLeptonicEvent.h"

class FullLepHypothesesFilter : public edm::EDFilter {

 public:

  explicit FullLepHypothesesFilter(const edm::ParameterSet&);
  ~FullLepHypothesesFilter(){};
  
 private:

  virtual void beginJob(const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:
 
  edm::InputTag hypoKey_;  
  edm::InputTag FullLepEvt_;
  edm::InputTag wgt_;
  edm::InputTag jets_;
  
  double qcdCut_;
  std::vector<double> zCut_;
  std::string bAlgo_;
  std::vector<double> bDisc_; 

  unsigned int beforeCuts_;
  unsigned int validCounter_;  
  unsigned int betweenCuts_;  
  unsigned int afterCuts_;
  double beforeCutsWeighted_;
  double validCounterWeighted_;  
  double betweenCutsWeighted_;
  double afterCutsWeighted_;
};

#endif
