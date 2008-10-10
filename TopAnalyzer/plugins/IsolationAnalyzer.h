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
#include "TopAnalysis/TopAnalyzer/plugins/IsolationHelper.h"
#include <TGraph.h>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

class IsolationAnalyzer : public edm::EDAnalyzer {

public:
	explicit IsolationAnalyzer(const edm::ParameterSet&);
	~IsolationAnalyzer();

private:
	typedef std::vector<pat::Muon> TopMuonCollection;
	typedef std::vector<pat::Jet>  TopJetCollection;

	virtual void beginJob(const edm::EventSetup&);
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();
	virtual unsigned int getClosestJet(edm::Handle<TopJetCollection>&, const pat::MET* &);


	std::string hist_;
	edm::InputTag muons_;
	edm::InputTag met_;
	edm::InputTag ttgen_;
	edm::InputTag jets_;

	vector<double> ptBins_;
	bool ttbarMC_;

	TH1F *varCorrelationsCaloIso_;
	TH1F *varCorrelationsTrackIso_;
	IsolationHelper *helper_, *helperTbar_;
	//for ttbar Binning
	vector<IsolationHelper> ttBarHelper_;
};
#endif
