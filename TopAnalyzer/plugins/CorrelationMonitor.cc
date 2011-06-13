#include "TopAnalysis/TopAnalyzer/plugins/CorrelationMonitor.h"
//#include "FWCore/MessageLogger/interface/MessageLogger.h"

CorrelationMonitor::CorrelationMonitor() {
}
CorrelationMonitor::~CorrelationMonitor() {
}

/**
 * fills a stored histogramm with name = name with values v1 and v2 with a certain weight
 */
void CorrelationMonitor::fill(std::string name, double v1, double v2,
		double weight) {
	map<string, TH2F*>::iterator iter = histos_.find(name);
	if (iter != histos_.end()) {
		iter->second->Fill(v1, v2, weight);
	}
}

void CorrelationMonitor::book(ofstream& file) {

}

/**
 * prints the correlation of all histogramms to the std output
 */
void CorrelationMonitor::printCorrelation() {
	map<string, TH2F*>::iterator iter;
	for (iter = histos_.begin(); iter != histos_.end(); ++iter) {
		cout << iter->second->GetTitle() << " " << iter->first << " Cor. "
				<< iter->second->GetCorrelationFactor(1, 2) << endl;
	}
}

/**
 * returns the correlation error
 */
Double_t CorrelationMonitor::getCorrelationFactor(std::string name) {
	map<string, TH2F*>::iterator iter = histos_.find(name);
	if (iter != histos_.end()) {
		return iter->second->GetCorrelationFactor(1, 2);
	} else {
		return -2.;// factor goes from -1 to 1
	}
}

/**
 * returns the correlation error
 */
Double_t CorrelationMonitor::getCorrelationError(std::string name) {
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
	Double_t rel1 = errrms1/rms1;
	Double_t rel2 = errrms2/rms2;
	//error calculation from http://.... or paper.....
	return sqrt((rel1*rel1 + rel2*rel2)*corr + err2cov/rms1/rms1/rms2/rms2);
}

/**
 * adds a histogramm to the collection and makes its draw option to 'COLZ'
 */
void CorrelationMonitor::addHist(string name, TH2F* hist) {
	hist->SetDrawOption("COLZ");
	histos_.insert(make_pair(name, hist));
}
