#ifndef DiMuonMassFilter_h  
#define DiMuonMassFilter_h

#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class DiMuonMassFilter : public edm::EDFilter {

 public:
  
  explicit DiMuonMassFilter(const edm::ParameterSet&);
  ~DiMuonMassFilter();

 private:

  virtual void beginJob(const edm::EventSetup&);
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag muons_;
  };  

#endif  
