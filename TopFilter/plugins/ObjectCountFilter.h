#ifndef ObjectCountFilter_h
#define ObjectCountFilter_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"

class ObjectCountFilter : public edm::EDFilter {

 public:

  explicit ObjectCountFilter(const edm::ParameterSet&);
  ~ObjectCountFilter(){};
  
 private:

  virtual void beginJob(const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 private:
 
  edm::InputTag wgt_;
  std::string name_;  
  edm::InputTag objects_;
  unsigned int size_;
  
  unsigned int beforeCuts_; 
  unsigned int afterCuts_;
  double beforeCutsWeighted_;
  double afterCutsWeighted_;
};

#endif
