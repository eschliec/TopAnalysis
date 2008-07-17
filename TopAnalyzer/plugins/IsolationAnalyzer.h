#ifndef IsolationAnalyzer_h
#define IsolationAnalyzer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TopAnalysis/TopUtils/interface/NameScheme.h"
#include "TopAnalysis/TopUtils/interface/RootSystem.h"
#include "TopAnalysis/TopUtils/interface/RootHistograms.h"
#include "TopAnalysis/TopAnalyzer/plugins/CorrelationMonitor.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class IsolationAnalyzer : public edm::EDAnalyzer {

public:
	explicit IsolationAnalyzer(const edm::ParameterSet&);
	~IsolationAnalyzer();

private:
	virtual void beginJob(const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	
	edm::InputTag muons_;
	edm::InputTag met_;
	edm::InputTag ttgen_;
	edm::InputTag jets_;
	
	CorrelationMonitor *isomon_;
	vector<double> ptBins_;
	vector<CorrelationMonitor*> hmonitors_;
	vector<CorrelationMonitor*> smonitors_;
	TH1F *thadpt_, *tleppt_;
	TH1F *hcaloCorr_, *htrackCorr_, *lcaloCorr_, *ltrackCorr_, *lptCorr_, *hptCorr_;
	TH1F *minDPhiMETJet_, *dPhiMETjet1_, *dPhiMETjet2_, *dPhiMETjet3_, *dPhiMETjet4_, *dPhiMETmuon_; 
	
	typedef std::vector<pat::Muon> TopMuonCollection;
	typedef std::vector<pat::Jet>  TopJetCollection;
};
#endif
