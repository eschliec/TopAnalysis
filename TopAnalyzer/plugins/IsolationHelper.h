#ifndef IsolationHelper_h
#define IsolationHelper_h

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
		jetMin_ = 0.;
		jetMax_ = 3.5;
		caloIso_ = 1000.;
		trackIso_ = 1000.;

		caloBins_ = 10;
		trackBins_ = 10;
		jetBins_ = 35;
	}

	void makeSummaryPlots(std::string prefix = "") {
		//for each histo, make 2 plots for correlation
		map<string, TH1F*>::iterator iter;
		unsigned int x = 0;
		string caloName = prefix + "caloCorrSummary";
		string trackName = prefix + "trackCorrSummary";

		//book histos for summary
		NameScheme nam("var");
		//open the file and don't delete the content (append)
		ofstream off(file_.c_str(), std::ios::app);
		summaryCalo_ = fs_->make<TH1F> (nam.name(off, caloName.c_str()), nam.name(caloName.c_str()), histos_.size(),
				-1., 1.);
		summaryTrack_ = fs_->make<TH1F> (nam.name(off, trackName.c_str()), nam.name(trackName.c_str()), histos_.size(),
				-1., 1.);

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
		ios_base::openmode mode;
		//rewrite the file when the first histogramm is written
		if (histos_.size() == 0) {
			mode = std::ios::out;
		} else {//otherwise append the histogramm to the hist file
			mode = std::ios::app;
			//cout << "hierjeztzt" << endl;
		}
		ofstream off(file_.c_str(), mode);
		std::string norm = "norm_" + name;
		TH1F *hist = fs_->make<TH1F> (nam.name(off, name.c_str()), nam.name(name.c_str()), numberOfBins, min, max);
		TH1F *histn = fs_->make<TH1F> (nam.name(off, norm.c_str()), nam.name(norm.c_str()), numberOfBins, min, max);
		histos_.insert(make_pair(name, hist));
		normhistos_.insert(make_pair(name, histn));
		//adding correlation histogramms
		mon->addHist(name + "VScalo", numberOfBins, min, max, caloBins_, caloMin_, caloMax_);
		mon->addHist(name + "VStrack", numberOfBins, min, max, trackBins_, trackMin_, trackMax_);
		mon->addHist(name + "VSjetiso", numberOfBins, min, max, jetBins_, jetMin_, jetMax_);
	}

	void fill(std::string name, double var, double caloIso, double trackIso, double jetIso, double weight) {
		map<string, TH1F*>::iterator iter = histos_.find(name);
		if (iter != histos_.end()) {
			iter->second->Fill(var, weight);
			mon->fill(name + "VScalo", var, caloIso, weight);
			mon->fill(name + "VStrack", var, trackIso, weight);
			mon->fill(name + "VSjetiso", var, jetIso, weight);
		}
		map<string, TH1F*>::iterator itern = normhistos_.find(name);
		if (itern != normhistos_.end()) {
			itern->second->Fill(var, weight);
		}
	}

	void fill(std::string name, double var) {
		map<string, TH1F*>::iterator iter = histos_.find(name);
		if (iter != histos_.end()) {
			iter->second->Fill(var, weight_);
			mon->fill(name + "VScalo", var, caloIso_, weight_);
			mon->fill(name + "VStrack", var, trackIso_, weight_);
			mon->fill(name + "VSjetiso", var, jetIso_, weight_);
		}
		map<string, TH1F*>::iterator itern = normhistos_.find(name);
		if (itern != normhistos_.end()) {
			itern->second->Fill(var, weight_);
		}
	}

	void setCaloIso(double caloIso) {
		caloIso_ = caloIso;
	}

	void setTrackIso(double trackIso) {
		trackIso_ = trackIso;
	}

	void setJetIso(double jetiso){
		jetIso_ = jetiso;
	}


	void setWeight(double weight) {
		weight_ = weight;
	}

	void normalize() {
		map<string, TH1F*>::iterator iter;
		for (iter = normhistos_.begin(); iter != normhistos_.end(); ++iter) {
			TH1F *hist = iter->second;
			double allEvents = hist->Integral();
			if (allEvents != 0)
				hist->Scale(1 / allEvents);
		}
	}

private:
	//for 2D correlation plots
	CorrelationMonitor *mon;
	//1D histos for the distribution of the variables
	std::map<string, TH1F*> histos_, normhistos_;
	TH1F * summaryCalo_, *summaryTrack_;
	edm::Service<TFileService> fs_;
	unsigned int caloBins_, trackBins_, jetBins_;
	double caloMin_, caloMax_, trackMin_, trackMax_, jetMin_, jetMax_;
	double caloIso_, trackIso_, jetIso_, weight_;
	std::string file_;
};
#endif

