#ifndef DimuonAnalyzer_h
#define DimuonAnalyzer_h

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "TString.h"
#include "TFile.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"

class TH1F;
class TH2F;
class DimuonAnalyzer : public edm::EDAnalyzer {

  public:

    explicit DimuonAnalyzer(const edm::ParameterSet&);
    ~DimuonAnalyzer();
    
  private:
  
    virtual void beginJob(const edm::EventSetup&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    bool matches(const reco::GenParticle*, const reco::GenParticle*);
    
    edm::InputTag muons_, ref_;
    bool matchMuons_;     
    bool useEvtWgt_;  
    std::string hist_;     
    typedef std::vector<pat::Muon> PatMuonCollection; 
    
    TH1F * dimassRC_;
    TH1F * dimassWC_;
      
            
    TH2F * isoDimassCorrelation_;   
    
    TH1F * absCount_;
    TH1F * relCount_;
    TH1F * combCount_;
    TH1F * diCombCount_;
};

#endif
