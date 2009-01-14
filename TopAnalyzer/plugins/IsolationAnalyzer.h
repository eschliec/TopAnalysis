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

#include "DataFormats/Math/interface/deltaR.h"

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
	template <class T>
	pat::Jet getClosestJet(edm::Handle<TopJetCollection>&, T);
	virtual pat::Jet getClosestJetInDeltaPhi(edm::Handle<TopJetCollection> &, double phi);

	double getJetIso(edm::Handle<TopJetCollection> jets, edm::Handle<TopMuonCollection> muons, double minJetPt) {
		TopMuonCollection::const_iterator muon = muons->begin();
		double minDR = -1.;
		if (muon != muons->end()) {
			for (TopJetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
				double dR = deltaR(muon->eta(), muon->phi(), jet->eta(), jet->phi());
				if ((minDR < 0 || dR < minDR) && jet->pt() > minJetPt)
					minDR = dR;
			}
		}
		return minDR;
	}


	std::string hist_;
	edm::InputTag muons_;
	edm::InputTag met_;
	edm::InputTag ttgen_;
	edm::InputTag jets_;

	vector<double> ptBins_;
	bool ttbarMC_, useMVA_;
	std::string module_;
	std::string discinput_;

	TH1F *recoMETUncorrectedMET_, *genMetRecoDiff_;

	IsolationHelper *helper_;
	//for ttbar Binning
	map<string, IsolationHelper*> ttBarHelper_;
};
#endif
