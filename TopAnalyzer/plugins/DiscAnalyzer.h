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
#include "TopAnalysis/TopAnalyzer/plugins/IsolationHelper.h"

class DiscAnalyzer: public edm::EDAnalyzer {
public:
	DiscAnalyzer(const edm::ParameterSet& cfg) :
		module_(cfg.getParameter<std::string> ("modulename")), discinput_(cfg.getParameter<std::string> (
				"discriminator")), hist_(cfg.getParameter<std::string> ("histfile")) {

	}
	~DiscAnalyzer() {
	}

	void beginJob(const edm::EventSetup& setup) {
		if (hist_.empty())
					return;
		edm::Service<TFileService> fs;
		if (!fs) {
			throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
		}

		helper_ = new IsolationHelper(fs, hist_);
		helper_->addHistogram("MVA_disc", 25, 0., 1.);
		//disc_ = fs->make<TH1F> ("MVA_disc", "MVA_disc", 25, 0., 1.);
		discNorm_ = fs->make<TH1F> ("norm_MVA_disc", "norm_MVA_disc", 25, 0., 1.);
	}
	void analyze(const edm::Event& evt, const edm::EventSetup& setup) {
		double weight, disc;
		edm::Handle<double> disc_handle;
		evt.getByLabel(module_, discinput_, disc_handle);
		//evt.getByLabel("DiscSel", disc_handle);
		edm::Handle<double> weightHandle;
		evt.getByLabel("eventWeight", weightHandle);

		disc = *disc_handle;
		weight = *weightHandle;
		disc_->Fill(disc, weight);
		discNorm_->Fill(disc, weight);
	}

	void endJob() {
		//normalize to 1
		if (discNorm_->Integral() != 0) {
			discNorm_->Sumw2();
			discNorm_->Scale(1 / discNorm_->Integral());
		}
	}

private:
	TH1F *disc_, *discNorm_;
	std::string module_;
	std::string discinput_;
	std::string hist_;
	IsolationsHelper helper_;

};
#endif
