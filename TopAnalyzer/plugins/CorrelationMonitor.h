#ifndef CorrelationMonitor_h
#define CorrelationMonitor_h

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

class CorrelationMonitor {
public:
	explicit CorrelationMonitor();
	//	IsolationMonitor(std::string);
	~CorrelationMonitor();
	void book(ofstream& file);

	void fill(std::string, double, double, double);
	/**
	 *  adds a histogram to the Isolation Monitor
	 */
	void addHist(std::string, TH2F*);
	/**
	 * prints the CorrelationFactors of the histograms
	 */
	void printCorrelation();
	/**
	 * Gives the correlation-factor of the histograms named 'string'
	 */
	Double_t getCorrelationFactor(std::string);
	/**
	 * Gives the correlation-factor error of the histograms named 'string'
	 */
	Double_t getCorrelationError(std::string);

private:
	std::map<string, TH2F*> histos_;
};
#endif
