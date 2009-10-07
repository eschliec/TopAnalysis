#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TopAnalysis/TopFilter/plugins/DiMuonMassFilter.h"

/// default constructor 
DiMuonMassFilter::DiMuonMassFilter(const edm::ParameterSet& cfg):
  muons_ (cfg.getParameter<edm::InputTag>(       "muons"  )),
  qcdCut_(cfg.getParameter<double>(              "qcdCut" )),
  zCut_  (cfg.getParameter<std::vector<double> >("zCut"   ))  
{
}

/// default destructor
DiMuonMassFilter::~DiMuonMassFilter()
{
}

/// sanity check
void 
DiMuonMassFilter::beginJob(const edm::EventSetup& setup)
{ 
  edm::LogError err("topFilter");
  if(zCut_.size()!=2)   err << "zCut has wrong size. Size has to be 2!\n";
  if(zCut_[0]>zCut_[1]) err << "Lower zCut is higher than upper zCut. All events will be skipped!\n";
}

/// event veto
bool
DiMuonMassFilter::filter(edm::Event& event, const edm::EventSetup& setup)
{   
  // fetch the input collection from the event content  
  edm::Handle<std::vector<pat::Muon> > muons;
  event.getByLabel(muons_, muons);
  
  bool pass=true;  
  if( muons->size()>1 ){ 
    reco::Particle::LorentzVector diMuon = (*muons)[0].p4() + (*muons)[1].p4();
    double diMuonMass = sqrt( diMuon.Dot(diMuon) );
    if( (diMuonMass>zCut_[0] && diMuonMass<=zCut_[1]) || diMuonMass < qcdCut_ ){ 
      pass=false;
    }    
  }    
  return pass;
}
