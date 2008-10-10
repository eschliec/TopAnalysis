#ifndef CorrelationMonitor_h
#define CorrelationMonitor_h

#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

class CorrelationMonitor {
public:
	CorrelationMonitor(edm::Service<TFileService> fs, std::string file) {
		fs_ = fs;
		file_ = file;
	}

	void fill(std::string name, double v1, double v2, double weight) {
		map<string, TH2F*>::iterator iter = histos_.find(name);
		if (iter != histos_.end()) {
			iter->second->Fill(v1, v2, weight);
		}
	}
	/**
	 *  adds a histogram to the Isolation Monitor
	 */
	void addHist(std::string name, TH2F* hist) {
		hist->SetDrawOption("COLZ");
		histos_.insert(make_pair(name, hist));
	}

	void addHist(std::string name, unsigned int numberOfBins1, double min1, double max1, unsigned int numberOfBins2,
			double min2, double max2) {
		NameScheme nam("corr");
		ofstream out(file_.c_str(), std::ios::out);
		TH2F *hist = fs_->make<TH2F> (nam.name(out, name.c_str()), nam.name(name.c_str()), numberOfBins1, min1, max1,
				numberOfBins2, min2, max2);
		addHist(name, hist);
	}

	/**
	 * prints the CorrelationFactors of the histograms
	 */
	void printCorrelation() {
		map<string, TH2F*>::iterator iter;
		for (iter = histos_.begin(); iter != histos_.end(); ++iter) {
			cout << iter->second->GetTitle() << " " << iter->first << " Cor. " << iter->second->GetCorrelationFactor(1,
					2) << endl;
		}
	}
	/**
	 * Gives the correlation-factor of the histograms named 'string'
	 */
	Double_t getCorrelationFactor(std::string name) {
		map<string, TH2F*>::iterator iter = histos_.find(name);
		if (iter != histos_.end()) {
			return iter->second->GetCorrelationFactor(1, 2);
		} else {
			return -2.;// factor goes from -1 to 1
		}
	}
	/**
	 * Gives the correlation-factor error of the histograms named 'string'
	 */
	Double_t getCorrelationError(std::string name) {
		map<string, TH2F*>::iterator iter = histos_.find(name);
		Double_t corr = -2.;
		Double_t rms1 = 0.;
		Double_t rms2 = 0.;
		Double_t errrms1 = 0.;
		Double_t errrms2 = 0.;
		Double_t err2cov = 0.;
		Double_t stats[7];
		if (iter != histos_.end()) {
			corr = iter->second->GetCorrelationFactor(1, 2);
			rms1 = iter->second->GetRMS(1);
			rms2 = iter->second->GetRMS(2);
			errrms1 = iter->second->GetRMSError(1);
			errrms2 = iter->second->GetRMSError(2);
			iter->second->GetStats(stats);
		}
		if (rms1 == 0)
			return 0;
		if (rms2 == 0)
			return 0;
		//squared error
		err2cov = stats[1];
		Double_t rel1 = errrms1 / rms1;
		Double_t rel2 = errrms2 / rms2;
		//error calculation from http://.... or paper.....
		return sqrt((rel1 * rel1 + rel2 * rel2) * corr + err2cov / rms1 / rms1 / rms2 / rms2);
	}

	void printSomething() {
		cout << "test" << endl;
	}

private:
	std::map<string, TH2F*> histos_;
	edm::Service<TFileService> fs_;
	std::string file_;
};
#endif
