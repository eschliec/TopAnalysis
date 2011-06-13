#include "TopAnalysis/TopAnalyzer/plugins/LeptonNumberAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using std::cout;
using std::endl;
using reco::GenParticle;

LeptonNumberAnalyzer::LeptonNumberAnalyzer(const edm::ParameterSet& cfg) :
	muons_(cfg.getParameter<edm::InputTag>("muons")),
			numberOfRatioBins_(cfg.getParameter<int>("numberOfRatioBins")),
			numberOfMuonBins_(cfg.getParameter<int>("numberOfMuonBins")),
			numberOfElecBins_(cfg.getParameter<int>("numberOfElecBins")),
			minRatio_(cfg.getParameter<double>("minRatio")),
			maxRatio_(cfg.getParameter<double>("maxRatio")),
			minNmuon_(cfg.getParameter<double>("minNmuon")),
			maxNmuon_(cfg.getParameter<double>("maxNmuon")),
			minNelec_(cfg.getParameter<double>("minNelec")),
			maxNelec_(cfg.getParameter<double>("maxNelec")) {
	eleCounter_ = 0;
	muCounter_ = 0;
	nOfEvents_ = 0;
}

LeptonNumberAnalyzer::~LeptonNumberAnalyzer() {
}

void LeptonNumberAnalyzer::beginJob(const edm::EventSetup&) {
	edm::Service<TFileService> fs;
	if ( !fs) {
		throw edm::Exception( edm::errors::Configuration,
				"TFile Service is not registered in cfg file" );
	}

	muonElecRatio_ = fs->make<TH1F>("MuonToElectronRatio",
			"MuonToElectronRatio", numberOfRatioBins_, minRatio_, maxRatio_);

	numberOfMuons_ = fs->make<TH1F>("numberOfMuons", "numberOfMuons",
			numberOfMuonBins_, minNmuon_, maxNmuon_);

	numberOfElecs_ = fs->make<TH1F>("numberOfElecs", "numberOfElecs",
			numberOfElecBins_, minNelec_, maxNelec_);
}

void LeptonNumberAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup) {
	nOfEvents_++;
	edm::Handle<reco::GenParticleCollection> genParticles;
	evt.getByLabel("genParticles", genParticles);

//	edm::Handle<edm::View<pat::Muon> > recMuons;
//	evt.getByLabel("muons", recMuons);

	edm::Handle<double> weightHandle;
	evt.getByLabel("eventWeight", weightHandle);

	double weight = *weightHandle;

	int evtEl, evtMu;
	evtEl = 0;
	evtMu = 0;
//	evtRMu = recMuons->size();
	for (reco::GenParticleCollection::const_iterator part =
			genParticles->begin(); part!=genParticles->end(); ++part) {
		if (fabs(part->pdgId()) == 13) {
			muCounter_++;
			evtMu++;
		}
		if (fabs(part->pdgId()) == 11) {
			eleCounter_++;
			evtEl++;
		}
	}
	double ratio = 0.;
	if (evtEl > 0)
		ratio = evtMu/evtEl;

	muonElecRatio_->Fill(ratio, weight);
	numberOfMuons_->Fill(evtMu, weight);
	numberOfElecs_->Fill(evtEl, weight);

}

void LeptonNumberAnalyzer::endJob() {
	cout << "# of e :" << eleCounter_ << endl;
	cout << "average # of e per event :"
	<< eleCounter_/nOfEvents_ << endl;
	cout << "# of mu :" << muCounter_ << endl;
	cout << "average # of mu per event:"
	<< muCounter_/nOfEvents_ << endl;
}
