#include "TopAnalysis/TopAnalyzer/plugins/CorrelationMonitor.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"

class IsolationHelper {
	//needs TFileService
public:
	IsolationHelper() {

	}
	IsolationHelper(edm::Service<TFileService> fs, std::string file) {
		mon = new CorrelationMonitor(fs, "Var2D.hist");
		file_ = file;
		fs_ = fs;
		caloMin_ = 0.;
		trackMin_ = 0.;
		caloMax_ = 10.;
		trackMax_ = 10.;
		caloIso_ = 1000.;
		trackIso_ = 1000.;

		caloBins_ = 10;
		trackBins_ = 10;
	}

	void makeSummaryPlots() {
		//for each histo, make 2 plots for correlation
		map<string, TH1F*>::iterator iter;
		unsigned int x = 0;

		//book histos for summary
		NameScheme nam("var");
		ofstream off(file_.c_str(), std::ios::out);
		summaryCalo_ = fs_->make<TH1F> (nam.name(off, "caloCorrSummary"), nam.name("caloCorrSummary"), histos_.size(),
				-1., 1.);
		summaryTrack_ = fs_->make<TH1F> (nam.name(off, "trackCorrSummary"), nam.name("trackCorrSummary"),
				histos_.size(), -1., 1.);

		for (iter = histos_.begin(); iter != histos_.end(); ++iter) {
			std::string varname = iter->first;
			double caloCorr = mon->getCorrelationFactor(varname + "VScalo");
			double trackCorr = mon->getCorrelationFactor(varname + "VStrack");
			x++;
			summaryCalo_->SetBinContent(x, caloCorr);
			summaryTrack_->SetBinContent(x, trackCorr);
			summaryCalo_->GetXaxis()->SetBinLabel(x, varname.c_str());
			summaryTrack_->GetXaxis()->SetBinLabel(x, varname.c_str());
		}
	}

	void addHistogram(std::string name, unsigned int numberOfBins, double min, double max) {
		NameScheme nam("var");
		ofstream off(file_.c_str(), std::ios::out);
		TH1F *hist = fs_->make<TH1F> (nam.name(off, name.c_str()), nam.name(name.c_str()), numberOfBins, min, max);
		histos_.insert(make_pair(name, hist));
		//adding correlation histogramms
		mon->addHist(name + "VScalo", numberOfBins, min, max, caloBins_, caloMin_, caloMax_);
		mon->addHist(name + "VStrack", numberOfBins, min, max, trackBins_, trackMin_, trackMax_);
	}

	void fill(std::string name, double var, double caloIso, double trackIso, double weight) {
		map<string, TH1F*>::iterator iter = histos_.find(name);
		if (iter != histos_.end()) {
			iter->second->Fill(var, weight);
			mon->fill(name + "VScalo", var, caloIso, weight);
			mon->fill(name + "VStrack", var, trackIso, weight);
		}
	}

	void fill(std::string name, double var) {
		map<string, TH1F*>::iterator iter = histos_.find(name);
		if (iter != histos_.end()) {
			iter->second->Fill(var, weight_);
			mon->fill(name + "VScalo", var, caloIso_, weight_);
			mon->fill(name + "VStrack", var, trackIso_, weight_);
		}
	}

	void setCaloIso(double caloIso) {
		caloIso_ = caloIso;
	}

	void setTrackIso(double trackIso) {
		trackIso_ = trackIso;
	}

	void setWeight(double weight) {
		weight_ = weight;
	}

	void normalize() {
		map<string, TH1F*>::iterator iter;
		for (iter = histos_.begin(); iter != histos_.end(); ++iter) {
			TH1F *hist = iter->second;
			double allEvents = hist->Integral();
			if(allEvents != 0) hist->Scale(1/allEvents);
		}
	}

private:
	//for 2D correlation plots
	CorrelationMonitor *mon;
	//1D histos for the distribution of the variables
	std::map<string, TH1F*> histos_;
	TH1F * summaryCalo_, *summaryTrack_;
	edm::Service<TFileService> fs_;
	unsigned int caloBins_, trackBins_;
	double caloMin_, caloMax_, trackMin_, trackMax_;
	double caloIso_, trackIso_, weight_;
	std::string file_;
};
