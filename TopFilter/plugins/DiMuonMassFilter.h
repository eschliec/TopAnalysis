#ifndef DiMuonMassFilter_h  
#define DiMuonMassFilter_h

#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/*
   \class   DiMuonMassFilter DiMuonMassFilter.h "TopAnalysis/TopFilter/plugins/DiMuonMassFilter.h"

   \brief   Plugin to veto events with two muons which give a certain invariant mass like Z, J/Psi, ...

   The class vetos events where the invariant mass of the two leading muons lies between two values given in the
   config-file or with an invariant mass lower than a certain value. Note that there is no selection of the muons 
   within the event filer. 
**/

class DiMuonMassFilter : public edm::EDFilter {

 public:
  /// default constructor
  explicit DiMuonMassFilter(const edm::ParameterSet& configFile);
  /// default destructor
  ~DiMuonMassFilter();
  
 private:
  /// sanity check 
  virtual void beginJob(const edm::EventSetup&);
  /// event veto
  virtual bool filter(edm::Event& event, const edm::EventSetup& setup);
  
 private:
  /// muon collection label
  edm::InputTag muons_;
  /// veto cut on low mass resonances like J/Psi, default is 20GeV 
  double qcdCut_;
  /// veto cuts on Z-mass, default values are 80GeV, 105GeV
  std::vector<double> zCut_;
  
};  

#endif  
