#include "TopAnalysis/TopAnalyzer/interface/JetTrigger.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

//#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

/// default constructor
JetTrigger::JetTrigger(const edm::ParameterSet& cfg)
{
  triggerResults_    = cfg.getParameter<edm::InputTag>("TriggerResults");
  triggerSummary_    = cfg.getParameter<edm::InputTag>("TriggerSummary");
  patTriggerEvent_   = cfg.getParameter<edm::InputTag>("patTriggerEvent");
  patTrigger_        = cfg.getParameter<edm::InputTag>("patTrigger");
  jets_              = cfg.getParameter<edm::InputTag>("jets");
  triggerMatchedJets_= cfg.getParameter<edm::InputTag>("triggerMatchedJets");
  analyzedTrigger_   = cfg.getParameter<edm::InputTag>("analyzedTrigger");
}

/// initiate histograms
void
JetTrigger::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;
  if( !fs ) throw edm::Exception( edm::errors::Configuration, "TFile Service is not registered in cfg file" );

  // pt jets
  //hists_["pt" ] = fs->make<TH1F>( "pt"  , "pt"  , 100,  0. , 100. );
  // pt 1. jet
  //hists_["pt1"] = fs->make<TH1F>( "pt1" , "pt1" , 100,  0. , 100. );
  // pt 2. jet
  //hists_["pt2"] = fs->make<TH1F>( "pt2" , "pt2" , 100,  0. , 100. );
  // pt 3. jet
  //hists_["pt3"] = fs->make<TH1F>( "pt3" , "pt3" , 100,  0. , 100. );
  // pt 4. jet
  //hists_["pt4"] = fs->make<TH1F>( "pt4" , "pt4" , 100,  0. , 100. );
  // pt 5. jet
  //hists_["pt5"] = fs->make<TH1F>( "pt5" , "pt5" , 100,  0. , 100. );
  // pt 6. jet
  //hists_["pt6"] = fs->make<TH1F>( "pt6" , "pt6" , 100,  0. , 100. );
  // pt jets
  //hists_["pt_trig" ] = fs->make<TH1F>( "pt_trig"  , "pt_trig"  , 100,  0. , 100. );
  // pt 1. jet
  hists_["pt1_pat"] = fs->make<TH1F>( "pt1_pat" , "pt1_pat" , 100,  0. , 100. );
  // pt 2. jet
  hists_["pt2_pat"] = fs->make<TH1F>( "pt2_pat" , "pt2_pat" , 100,  0. , 100. );
  // pt 3. jet
  hists_["pt3_pat"] = fs->make<TH1F>( "pt3_pat" , "pt3_pat" , 100,  0. , 100. );
  // pt 4. jet
  hists_["pt4_pat"] = fs->make<TH1F>( "pt4_pat" , "pt4_pat" , 100,  0. , 100. );
  // pt 5. jet
  hists_["pt5_pat"] = fs->make<TH1F>( "pt5_pat" , "pt5_pat" , 100,  0. , 100. );
  // pt 6. jet
  hists_["pt6_pat"] = fs->make<TH1F>( "pt6_pat" , "pt6_pat" , 100,  0. , 100. );
  // jet multiplicity
  //hists_["n"] = fs->make<TH1F>( "n" , "n" , 10,  -0.5 , 9.5 );
  // jet multiplicity
  hists_["n_pat"] = fs->make<TH1F>( "n_pat" , "n_pat" , 10,  -0.5 , 9.5 );
  // Trigger vs. offline p_{T} (GeV)
  hists2D_["ptTrigOff" ] = fs->make<TH2F>( "ptTrigOff" , "Trigger vs. offline p_{T} (GeV)", 60,  0., 300., 60,  0., 300. );
  // Trigger vs. offline p_{T} (GeV)
  //hists2D_["ptTrigOff2"] = fs->make<TH2F>( "ptTrigOff2", "Trigger vs. offline p_{T} (GeV)", 60,  0., 300., 60,  0., 300. );
  // Trigger vs. offline eta
  hists2D_["etaTrigOff"] = fs->make<TH2F>( "etaTrigOff", "Trigger vs. offline #eta"       , 50, -2.5, 2.5, 50, -2.5, 2.5 );
  // Trigger vs. offline phi
  hists2D_["phiTrigOff"] = fs->make<TH2F>( "phiTrigOff", "Trigger vs. offline #phi"       , 64, -3.2, 3.2, 64, -3.2, 3.2 );
 }

/// analyze triggers and fill histograms
void
JetTrigger::analyze(const edm::Event& event, const edm::EventSetup&)
{
  //edm::Handle<edm::TriggerResults> triggerResults; 
  //event.getByLabel(triggerResults_, triggerResults); 
  //edm::Handle<trigger::TriggerEvent> triggerSummary; 
  //event.getByLabel(triggerSummary_, triggerSummary); 
  edm::Handle<pat::TriggerEvent> triggerEvent;
  event.getByLabel(patTriggerEvent_, triggerEvent);
  //edm::Handle<std::vector<pat::TriggerFilter> > triggerFilters;
  //event.getByLabel(patTrigger_, triggerFilters);
  //edm::Handle<std::vector<pat::TriggerPath> > triggerPaths;
  //event.getByLabel(patTrigger_, triggerPaths);
  //edm::Handle<std::vector<pat::TriggerObject> > triggerObjects;
  //event.getByLabel(patTrigger_, triggerObjects);
  //edm::Handle<edm::View<pat::Jet> > jets;
  edm::Handle<pat::JetCollection> jets;
  event.getByLabel(jets_, jets);
  //edm::Handle<edm::View<pat::Jet> > triggerMatchedJets;
  //event.getByLabel(triggerMatchedJets_, triggerMatchedJets);

//  if( triggerResults.failedToGet() ) {
//    std::cout << "---------------------------" << std::endl;
//    std::cout << "-- NO TRIGGER RESULTS !! --" << std::endl;
//    std::cout << "---------------------------" << std::endl;
//  }
//  if( triggerSummary.failedToGet() ) {
//    std::cout << "---------------------------" << std::endl;
//    std::cout << "-- NO TRIGGER SUMMARY !! --" << std::endl;
//    std::cout << "---------------------------" << std::endl;
//  }
//  if( triggerEvent.failedToGet() ) {
//    std::cout << "---------------------------" << std::endl;
//    std::cout << "--- NO TRIGGER EVENT !! ---" << std::endl;
//    std::cout << "---------------------------" << std::endl;
//  }
  if( !triggerEvent.failedToGet() ) {
    //std::cout << "---------------------------" << std::endl;
    //std::cout << "-- TRIGGER EVENT OPEN !! --" << std::endl;
    //std::cout << "---------------------------" << std::endl;
    pat::TriggerObjectRefVector objects = triggerEvent->filterObjects("hlt4jet30");
    unsigned int jetCounter = 1;
    hists_.find("n_pat")->second->Fill( objects.size() );
    for(pat::TriggerObjectRefVector::const_iterator obj = objects.begin(); obj != objects.end(); ++obj, ++jetCounter){
      if     (jetCounter == 1) hists_.find("pt1_pat")->second->Fill( (*obj)->pt() );
      else if(jetCounter == 2) hists_.find("pt2_pat")->second->Fill( (*obj)->pt() );
      else if(jetCounter == 3) hists_.find("pt3_pat")->second->Fill( (*obj)->pt() );
      else if(jetCounter == 4) hists_.find("pt4_pat")->second->Fill( (*obj)->pt() );
      else if(jetCounter == 5) hists_.find("pt5_pat")->second->Fill( (*obj)->pt() );
      else if(jetCounter == 6) hists_.find("pt6_pat")->second->Fill( (*obj)->pt() );
    }

    // PAT trigger helper for trigger matching information
    const pat::helper::TriggerMatchHelper matchHelper;

    const pat::TriggerObjectMatch * triggerMatch( triggerEvent->triggerObjectMatchResult( "jetTriggerMatchHLTJets" ) );
    for (size_t iJet = 0; iJet < jets->size(); ++iJet) {
      const reco::CandidateBaseRef candBaseRef( pat::JetRef( jets , iJet) );
      const pat::TriggerObjectRef trigRef( matchHelper.triggerMatchObject( candBaseRef, triggerMatch, event, *triggerEvent ) );
      // fill histograms
      if ( trigRef.isAvailable() ) {
	hists2D_.find("phiTrigOff")->second->Fill( candBaseRef->phi(), trigRef->phi() );
	hists2D_.find("etaTrigOff")->second->Fill( candBaseRef->eta(), trigRef->eta() );
	hists2D_.find("ptTrigOff" )->second->Fill( candBaseRef->pt() , trigRef->pt()  );
      }
    }
  }
//  if( triggerFilters.failedToGet() ) {
//    std::cout << "---------------------------" << std::endl;
//    std::cout << "-- NO TRIGGER FILTERS !! --" << std::endl;
//    std::cout << "---------------------------" << std::endl;
//  }
//  if( !triggerFilters.failedToGet() ) {
//    //std::cout << "---------------------------" << std::endl;
//    //std::cout << "- TRIGGER FILTERS OPEN !! -" << std::endl;
//    //std::cout << "---------------------------" << std::endl;
//  }
//  if( triggerPaths.failedToGet() ) {
//    std::cout << "---------------------------" << std::endl;
//    std::cout << "--- NO TRIGGER PATHS !! ---" << std::endl;
//    std::cout << "---------------------------" << std::endl;
//  }
//  if( !triggerPaths.failedToGet() ) {
//    //std::cout << "---------------------------" << std::endl;
//    //std::cout << "-- TRIGGER PATHS OPEN !! --" << std::endl;
//    //std::cout << "---------------------------" << std::endl;
//  }
//  if( triggerObjects.failedToGet() ) {
//    std::cout << "---------------------------" << std::endl;
//    std::cout << "-- NO TRIGGER OBJECTS !! --" << std::endl;
//    std::cout << "---------------------------" << std::endl;
//  }
//  if( !triggerObjects.failedToGet() ) {
//    //std::cout << "---------------------------" << std::endl;
//    //std::cout << "- TRIGGER OBJECTS OPEN !! -" << std::endl;
//    //std::cout << "---------------------------" << std::endl;
//  }
//  if( !triggerResults.failedToGet() ) {
//  }
//  /*
//  for(std::vector<pat::TriggerObject>::const_iterator obj = triggerObjects->begin(); obj != triggerObjects->end(); ++obj){
//    //std::cout << obj->collection() << std::endl;
//    for(std::vector<int>::const_iterator filt = obj->filterIds().begin(); filt != obj->filterIds().end(); ++filt){
//      //std::cout << *filt << std::endl;
//    }
//  }
//
//  int triggerMatchedJetCounter = 0;
//  for(edm::View<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet){
//    hists_.find("pt")->second->Fill( jet->pt());
//    if(jet == jets->begin())                               hists_.find("pt1")->second->Fill( jet->pt());
//    if(jets->size() > 1 && jet == ++jets->begin())         hists_.find("pt2")->second->Fill( jet->pt());
//    if(jets->size() > 2 && jet == ++(++jets->begin()))     hists_.find("pt3")->second->Fill( jet->pt());
//    if(jets->size() > 3 && jet == ++(++(++jets->begin()))) hists_.find("pt4")->second->Fill( jet->pt());
//    for(edm::View<pat::Jet>::const_iterator trigJet = triggerMatchedJets->begin(); trigJet != triggerMatchedJets->end(); ++trigJet){
//      if(jet == jets->begin()) {
//	hists_.find("pt_trig")->second->Fill( trigJet->pt());
//	if(trigJet == triggerMatchedJets->begin()) hists_.find("pt1_trig")->second->Fill( trigJet->pt());
//	if(triggerMatchedJets->size() > 1 && trigJet == ++triggerMatchedJets->begin())         hists_.find("pt2_trig")->second->Fill( trigJet->pt());
//	if(triggerMatchedJets->size() > 2 && trigJet == ++(++triggerMatchedJets->begin()))     hists_.find("pt3_trig")->second->Fill( trigJet->pt());
//	if(triggerMatchedJets->size() > 3 && trigJet == ++(++(++triggerMatchedJets->begin()))) hists_.find("pt4_trig")->second->Fill( trigJet->pt());
//      }
//      if(jet->pt() == trigJet->pt()){
//	if(jet->eta() == trigJet->eta()){
//	  if(jet->phi() == trigJet->phi()){
//	    ++triggerMatchedJetCounter;
//	  }
//	}
//      }
//    }
//  }
//  hists_.find("n"    )->second->Fill( jets->size() );
//  hists_.find("trign")->second->Fill( triggerMatchedJetCounter );
//  */
//  if( !triggerSummary.failedToGet() ) {
///*
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    for(unsigned int i=0; i<triggerSummary->collectionTags().size(); ++i){
//      std::cout << triggerSummary->collectionTag(i) << std::endl;
//    }
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//*/
///*
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    for(unsigned int i=0; i<triggerSummary->sizeFilters(); ++i){
//      std::cout << "Tag: " << triggerSummary->filterTag(i) << std::endl;
//      //for(unsigned int j=0; j<triggerSummary->filterKeys(i).size(); ++j){
//      //std::cout << "Key: " << triggerSummary->filterKeys(i)[j] << std::endl;
//      //}
//    }
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//*/
///*
//    std::cout << "Collection size : " << triggerSummary->collectionTags().size()                     << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "sizeCollections : " << triggerSummary->sizeCollections()                           << std::endl;
//    std::cout << "hltQuadJet30    : " << triggerSummary->collectionIndex((edm::InputTag)"hlt4jet30") << std::endl;
//    std::cout << "hltHT200        : " << triggerSummary->collectionIndex((edm::InputTag)"hltHT200")  << std::endl;
//
//    std::cout << "sizeFilters  : " << triggerSummary->sizeFilters()                           << std::endl;
//*/
///*
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "hltQuadJet30 filterIndex : " << triggerSummary->filterIndex(analyzedTrigger_) << std::endl;
//    std::cout << "hltQuadJet30 filterKeys  : " << triggerSummary->filterKeys(triggerSummary->filterIndex(analyzedTrigger_)).size() << std::endl;
//*/
//    int jetCounter = 1;
//    hists_.find("n")->second->Fill( triggerSummary->filterKeys(triggerSummary->filterIndex(analyzedTrigger_)).size() );
//    for(std::vector<short unsigned int>::const_iterator idx  = triggerSummary->filterKeys(triggerSummary->filterIndex(analyzedTrigger_)).begin();
//	                                                idx != triggerSummary->filterKeys(triggerSummary->filterIndex(analyzedTrigger_)).end(); ++idx, ++jetCounter){
//      /*
//      std::cout << "-------------------------" << std::endl;
//      std::cout << "hltQuadJet30 filterKey   : " << *idx << std::endl;
//      std::cout << "hltQuadJet30 filterObject: " << triggerSummary->getObjects().at(*idx).pt() << std::endl;
//      */
//      if     (jetCounter == 1) hists_.find("pt1")->second->Fill( triggerSummary->getObjects().at(*idx).pt() );
//      else if(jetCounter == 2) hists_.find("pt2")->second->Fill( triggerSummary->getObjects().at(*idx).pt() );
//      else if(jetCounter == 3) hists_.find("pt3")->second->Fill( triggerSummary->getObjects().at(*idx).pt() );
//      else if(jetCounter == 4) hists_.find("pt4")->second->Fill( triggerSummary->getObjects().at(*idx).pt() );
//      else if(jetCounter == 5) hists_.find("pt5")->second->Fill( triggerSummary->getObjects().at(*idx).pt() );
//      else if(jetCounter == 6) hists_.find("pt6")->second->Fill( triggerSummary->getObjects().at(*idx).pt() );
//    }
//    /*
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "hltQuadJet30 collectionIndex : " << triggerSummary->collectionIndex(analyzedTrigger_) << std::endl;
//    std::cout << "hltQuadJet30 collectionKeys  : " << triggerSummary->collectionKeys(triggerSummary->collectionIndex(analyzedTrigger_)).size() << std::endl;
//    for(std::vector<short unsigned int>::const_iterator idx = triggerSummary->collectionKeys(triggerSummary->collectionIndex(analyzedTrigger_)).begin(); idx != triggerSummary->collectionKeys(triggerSummary->collectionIndex(analyzedTrigger_)).end(); ++idx){
//      std::cout << "-------------------------" << std::endl;
//      std::cout << "hltQuadJet30 filterKey   : " << *idx << std::endl;
//      std::cout << "hltQuadJet30 filterObject: " << triggerSummary->getObjects().at(*idx).pt() << std::endl;
//    }
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "hltHT200     : " << triggerSummary->filterIndex((edm::InputTag)"hltHT200")  << std::endl;
//
//    std::cout << "sizeObjects     : " << triggerSummary->sizeObjects()                        << std::endl;
//*/
///*
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//    int jetCounter = 0;
//    for(unsigned int i=0; i<triggerSummary->getObjects().size(); ++i){
//      std::cout << "ID (" << i << ") : " << triggerSummary->getObjects()[i].id() << std::endl;
//      if(triggerSummary->getObjects()[i].id() == -84) ++jetCounter;
//      if     (jetCounter == 1) hists_.find("pt1")->second->Fill( triggerSummary->getObjects()[i].pt() );
//      else if(jetCounter == 2) hists_.find("pt2")->second->Fill( triggerSummary->getObjects()[i].pt() );
//      else if(jetCounter == 3) hists_.find("pt3")->second->Fill( triggerSummary->getObjects()[i].pt() );
//      else if(jetCounter == 4) hists_.find("pt4")->second->Fill( triggerSummary->getObjects()[i].pt() );
//    }
//    std::cout << "-------------------------" << std::endl;
//    std::cout << "-------------------------" << std::endl;
//*/
//  }

}

/// empty
void
JetTrigger::endJob()
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( JetTrigger );
