#ifndef TtSemiLepSignalSelectorMVAComputer_h
#define TtSemiLepSignalSelectorMVAComputer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/MVAComputer/interface/HelperMacros.h"
#include "PhysicsTools/MVAComputer/interface/MVAComputerCache.h"

#ifndef TtSemiLepSignalSelectorMVARcd_defined
#define TtSemiLepSignalSelectorMVARcd_defined
MVA_COMPUTER_CONTAINER_DEFINE(TtSemiLepSignalSelectorMVA);
#endif

class TtSemiLepSignalSelectorMVAComputer: public edm::EDProducer {

public:

	explicit TtSemiLepSignalSelectorMVAComputer(const edm::ParameterSet&);
	~TtSemiLepSignalSelectorMVAComputer();

private:

	virtual void beginJob(const edm::EventSetup&);
	virtual void produce(edm::Event& evt, const edm::EventSetup& setup);
	virtual void endJob();

	edm::InputTag leptons_;
	edm::InputTag jets_;
	edm::InputTag METs_;

	unsigned int nJetsMax_;

	PhysicsTools::MVAComputerCache mvaComputer;
};

#endif
