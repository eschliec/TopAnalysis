#include "TopAnalysis/TopFilter/plugins/NewTriggerTestFilter.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "FWCore/Framework/interface/Event.h"

NewTriggerTestFilter::NewTriggerTestFilter(const edm::ParameterSet& cfg):
  whichTrigger_  ( cfg.getParameter<std::string>  ("whichTrigger"  ) ),
  useEventWeight_( cfg.getParameter<bool>         ("useEventWeight") ),
  weight_        ( cfg.getParameter<edm::InputTag>("weight"        ) ),
  beforeCut_( 0 ), afterCut_( 0 ), beforeCutWeighted_( 0. ), afterCutWeighted_( 0. )
{
}

bool NewTriggerTestFilter::filter(edm::Event& event, const edm::EventSetup& setup)
{
  double weight = 1.;
  if(useEventWeight_) { 
    edm::Handle<double> weightHandle;
    event.getByLabel(weight_, weightHandle);
    weight = *weightHandle;
  }

  ++beforeCut_;
  beforeCutWeighted_ += weight;

  bool pass = false;

  edm::Handle<pat::TriggerEvent> triggerEvent;
  event.getByLabel("patTriggerEvent", triggerEvent);

  pat::TriggerObjectRefVector objects = triggerEvent->filterObjects("hlt4jet30");

  if(whichTrigger_ == "QuadJet80303030"){
    if(objects[0]->pt() >= 80){
      pass = true;
    }
  }
  else if(whichTrigger_ == "QuadJet40"){
    if(objects[0]->pt() >= 40){
      if(objects[1]->pt() >= 40){
	if(objects[2]->pt() >= 40){
	  if(objects[3]->pt() >= 40){
	    pass = true;
	  }
	}
      }
    }
  }
  else if(whichTrigger_ == "FiveJet30"){
    if(objects.size() >=5){
      if(objects[4]->pt() >= 30){
	pass = true;
      }
    }
  }
  if(pass){
    std::cout << "-------------------" << std::endl;
    std::cout << "pt 1: " << objects[0]->pt() << std::endl;
    std::cout << "pt 2: " << objects[1]->pt() << std::endl;
    std::cout << "pt 3: " << objects[2]->pt() << std::endl;
    std::cout << "pt 4: " << objects[3]->pt() << std::endl;
  }
  return pass;
}

void NewTriggerTestFilter::beginJob(const edm::EventSetup& setup)
{
}

void NewTriggerTestFilter::endJob()
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( NewTriggerTestFilter );
