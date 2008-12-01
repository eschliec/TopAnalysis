#ifndef DiscAnalyzer_h
#define DiscAnalyzer_h

#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

class DiscAnalyzer : public edm::EDAnalyzer {
public:
	DiscAnalyzer(const edm::ParameterSet& set) {

	}
	~DiscAnalyzer() {
	}

	void beginJob(const edm::EventSetup& setup) {
		edm::Service<TFileService> fs;
		if (!fs) {
			throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
		}
		std::cout << "test" << std::endl;
		disc_ = fs->make<TH1F> ("MVA_disc", "MVA_disc", 25, 0., 1.);
	}
	void analyze(const edm::Event& evt, const edm::EventSetup& setup) {
		double weight, disc;
		edm::Handle<double> disc_handle;
		evt.getByLabel("findTtSemiLepSignalSelectorMVA", "DiscSel", disc_handle);
		//evt.getByLabel("DiscSel", disc_handle);
		edm::Handle<double> weightHandle;
		evt.getByLabel("eventWeight", weightHandle);

		disc = *disc_handle;
		weight = *weightHandle;
		disc_->Fill(disc, weight);
	}

	void endJob(){
		//normalize to 1
		if(disc_->Integral() != 0)
		{
			disc_->Scale(1/disc_->Integral());
		}
	}

private:
	TH1F *disc_;

};
#endif
