#include "TopAnalysis/TopAnalyzer/plugins/IsolationAnalyzer.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include <math.h>

using std::cout;
using std::endl;
using reco::GenParticle;
using edm::LogInfo;
/**
 * Tests several variables for correlation with lepton isolation
 */
IsolationAnalyzer::IsolationAnalyzer(const edm::ParameterSet& cfg) :
	hist_(cfg.getParameter<std::string> ("hist")), muons_(cfg.getParameter<edm::InputTag> ("muons")), met_(
			cfg.getParameter<edm::InputTag> ("missingEt")), ttgen_(cfg.getParameter<edm::InputTag> ("genEvent")),
			jets_(cfg.getParameter<edm::InputTag> ("jets")), ptBins_(cfg.getParameter<std::vector<double> > ("ptBins")),
			ttbarMC_(cfg.getParameter<bool> ("ttbarMC")) {
}

void IsolationAnalyzer::beginJob(const edm::EventSetup&) {
	if (hist_.empty())
		return;

	edm::Service<TFileService> fs;
	if (!fs) {
		throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
	}

	helper_ = new IsolationHelper(fs, hist_);

	if (ttbarMC_) {
		helperTbar_ = new IsolationHelper(fs, "ttbarMC_" + hist_);

		//for each ptBin new IsolationHelper
		for (unsigned int x = 0; x < ptBins_.size() - 1; x++) {
			//helperTbar
		}
	}
	helper_->addHistogram("invariantMassJ3andJ4", 50, 0., 200.);
	helper_->addHistogram("minDeltaPhiMETJets", 100, 0., 4.);
	helper_->addHistogram("deltaPhiMetJet1", 80, -4., 4.);
	helper_->addHistogram("deltaPhiMetJet2", 80, -4., 4.);
	helper_->addHistogram("deltaPhiMetJet3", 80, -4., 4.);
	helper_->addHistogram("deltaPhiMetJet4", 80, -4., 4.);
	helper_->addHistogram("deltaPhiMetleadingMuon", 80, -4., 4.);
	helper_->addHistogram("deltaPhiMetMuons", 80, -4., 4.);
	helper_->addHistogram("DeltaPhiTimesDeltaEta", 50, -10., 10.);
	helper_->addHistogram("METTimesleadingJetEt", 70, 0., 7000.);
	helper_->addHistogram("Jet3EtOverJet1EtJet3Et", 30, 0., 6.);
	helper_->addHistogram("Jet3EtOverJet2EtJet3Et", 30, 0., 6.);
	helper_->addHistogram("Jet4EtOverJet1EtJet4Et", 30, 0., 6.);
	helper_->addHistogram("Jet4EtOverJet2EtJet3Et", 30, 0., 6.);
	helper_->addHistogram("JetEtSum34", 150, 0., 300.);
	helper_->addHistogram("DeltaPhiMuonJet3", 80, -4., 4.);
	helper_->addHistogram("DeltaPhiMuonJet4", 80, -4., 4.);
	helper_->addHistogram("DeltaPhiMuonJet12", 80, -4., 4.);
	helper_->addHistogram("SumJet1ET2TimesMuEt", 80, 60., 500);
	helper_->addHistogram("Sum4JetsMuEt", 100, 100., 600);
	helper_->addHistogram("VectorSumJetsMu", 60, 0., 300.);
	helper_->addHistogram("MET", 100, 0., 500.);
	NameScheme nam("var");
	ofstream off(hist_.c_str(), std::ios::out);
	recoMETUncorrectedMET_ = fs->make<TH1F> (nam.name(off, "recoMETUncorrectedMET"),
			nam.name("recoMETUncorrectedMET"), 200, -100., 100.);
}

IsolationAnalyzer::~IsolationAnalyzer() {
}

void IsolationAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup) {

	const pat::MET *met;
	//	const pat::Muon *leadingMuon;
	reco::Particle::LorentzVector vec;

	edm::Handle<std::vector<pat::MET> > metH;
	evt.getByLabel(met_, metH);

	edm::Handle<TopMuonCollection> muons;
	evt.getByLabel(muons_, muons);

	//the event
	edm::Handle<TtGenEvent> genEvt;
	evt.getByLabel(ttgen_, genEvt);

	edm::Handle<double> weightHandle;
	evt.getByLabel("eventWeight", weightHandle);

	edm::Handle<TopJetCollection> jets;
	evt.getByLabel(jets_, jets);

	met = &(*metH)[0];

	double weight = *weightHandle;

	if (jets->size() >= 4) {
		//get the first for jets
		TopJetCollection::const_iterator jet = jets->begin();

		double jet1Phi, jet2Phi, jet3Phi, jet4Phi;
		double jet1Eta, jet2Eta, jet3Eta, jet4Eta;
		double jet1pt, jet2pt, jet3pt, jet4pt;
		double jet1Et, jet2Et, jet3Et, jet4Et;
		double dpTde = 10000.;
		unsigned int jetno = getClosestJet(jets, met);
		unsigned int x = 1;

		jet1pt = jet->pt();
		jet1Et = jet->et();
		jet1Eta = jet->eta();
		jet1Phi = jet->phi();

		vec = jet->p4();
		if (x == jetno)
			dpTde = fabs(deltaPhi(met->phi(), jet->phi())) * fabs(met->eta() - jet->eta());
		x++;
		++jet;
		jet2pt = jet->pt();
		jet2Et = jet->et();
		jet2Eta = jet->eta();
		jet2Phi = jet->phi();
		vec += jet->p4();

		if (x == jetno)
			dpTde = fabs(deltaPhi(met->phi(), jet->phi())) * fabs(met->eta() - jet->eta());
		x++;
		++jet;
		jet3pt = jet->pt();
		jet3Et = jet->et();
		jet3Eta = jet->eta();
		jet3Phi = jet->phi();
		vec += jet->p4();

		if (x == jetno)
			dpTde = fabs(deltaPhi(met->phi(), jet->phi())) * fabs(met->eta() - jet->eta());
		x++;
		++jet;
		jet4pt = jet->pt();
		jet4Et = jet->et();
		jet4Eta = jet->eta();
		jet4Phi = jet->phi();
		vec += jet->p4();

		if (x == jetno)
			dpTde = fabs(deltaPhi(met->phi(), jet->phi())) * fabs(met->eta() - jet->eta());
		x++;

		//how big is the minimal deltaPhi between MET and one of the 4 jets
		double mindp;
		mindp = min(fabs(deltaPhi(met->phi(), jet1Phi)), fabs(deltaPhi(met->phi(), jet2Phi)));
		mindp = min(mindp, fabs(deltaPhi(met->phi(), jet3Phi)));
		mindp = min(mindp, fabs(deltaPhi(met->phi(), jet4Phi)));
		//transversal mass of jet 3+4
		double mt34 = sqrt((jet3Et + jet4Et) * (jet3Et + jet4Et) - (jet3pt + jet4pt) * (jet3pt + jet4pt));
		//transversal mass of jet 1+3+4 (jet1 could be switched with jet2)
		double mt3Jet = sqrt((jet1Et + jet3Et + jet4Et) * (jet1Et + jet3Et + jet4Et) - (jet1pt + jet3pt + jet4pt)
				* (jet1pt + jet3pt + jet4pt));
		cout << mt34 << endl;
		cout << mt3Jet << endl;

		//set the event-weight
		helper_->setWeight(weight);
		//now the muons
		unsigned int size = muons->size();
		for (unsigned int i = 0; i < size; ++i) {
			pat::Muon mu = (pat::Muon) (*muons)[i];
			double caloIso = mu.caloIso();
			double trackIso = mu.trackIso();
			//set isolation
			helper_->setCaloIso(caloIso);
			helper_->setTrackIso(trackIso);
			if (i == 0) {
				//leading muon
				vec += mu.p4();
				helper_->fill("METTimesleadingJetEt", jet1Et * met->et());
				helper_->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
				helper_->fill("deltaPhiMetJet1", deltaPhi(met->phi(), jet1Phi));
				helper_->fill("deltaPhiMetJet2", deltaPhi(met->phi(), jet2Phi));
				helper_->fill("deltaPhiMetJet3", deltaPhi(met->phi(), jet3Phi));
				helper_->fill("deltaPhiMetJet4", deltaPhi(met->phi(), jet4Phi));
				helper_->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
				helper_->fill("DeltaPhiTimesDeltaEta", deltaPhi(met->phi(), jet1Phi) * (met->eta() * jet1Eta));
				helper_->fill("METTimesleadingJetEt", met->et() * jet1Et);
				helper_->fill("Jet3EtOverJet1EtJet3Et", jet3Et / (jet3Et + jet1Et));
				helper_->fill("Jet3EtOverJet2EtJet3Et", jet3Et / (jet3Et + jet2Et));
				helper_->fill("Jet4EtOverJet1EtJet4Et", jet4Et / (jet4Et + jet1Et));
				helper_->fill("Jet4EtOverJet2EtJet3Et", jet4Et / (jet4Et + jet2Et));
				helper_->fill("JetEtSum34", jet3Et + jet4Et);
				helper_->fill("DeltaPhiMuonJet3", deltaPhi(mu.phi(), jet3Phi));
				helper_->fill("DeltaPhiMuonJet4", deltaPhi(mu.phi(), jet4Phi));
				helper_->fill("DeltaPhiMuonJet12", deltaPhi(mu.phi(), jet1Phi) + deltaPhi(mu.phi(), jet2Phi));
				helper_->fill("invariantMassJ3andJ4", mt34);
				helper_->fill("Sum4JetsMuEt", jet1Et + jet2Et + jet3Et + jet4Et + mu.et());
				helper_->fill("SumJet1ET2TimesMuEt", jet1Et + 2 * mu.et());
				helper_->fill("VectorSumJetsMu", vec.pt());
				helper_->fill("MET", met->et());
				recoMETUncorrectedMET_->Fill(met->et() - met->uncorrectedPt(), weight);
			}
			if (i == 2) {
				//2nd leading muon
			}
			helper_->fill("deltaPhiMetMuons", deltaPhi(met->phi(), mu.phi()), caloIso, trackIso, weight);

			//TODO: get top-pt from reco lvl
			/*			if (ttbarMC_) {
			 //			cout << "ttbar iso muin" << endl;
			 if (thad)
			 hDeltaPhi_->SetPoint(event_, thad->pt(), fabs(deltaPhi(mu.phi(), met->phi())));
			 if (tlep)
			 lDeltaPhi_->SetPoint(event_, tlep->pt(), fabs(deltaPhi(mu.phi(), met->phi())));
			 for (unsigned int i = 0; i < ptBins_.size() - 1; i++) {
			 if (thad && tlep) {
			 if (thad->pt() >= ptBins_[i] && thad->pt() < ptBins_[i + 1]) {
			 hmonitors_[i]->fill("trackCorrelation", met->pt(), mu.trackIso(), weight);
			 hmonitors_[i]->fill("caloCorrelation", met->pt(), mu.caloIso(), weight);
			 hmonitors_[i]->fill("pt", met->pt(), mu.pt(), weight);
			 }

			 if (tlep->pt() >= ptBins_[i] && tlep->pt() < ptBins_[i + 1]) {
			 //					smonitors_[i]->fill(mu, *met, weight);
			 smonitors_[i]->fill("trackCorrelation", log(met->pt()), mu.trackIso(), weight);
			 smonitors_[i]->fill("caloCorrelation", met->pt(), mu.caloIso(), weight);
			 smonitors_[i]->fill("pt", met->pt(), mu.pt(), weight);
			 }
			 }
			 }
			 }*/
		}

	}
}

void IsolationAnalyzer::endJob() {
	helper_->makeSummaryPlots();
	helper_->normalize();
	if (ttbarMC_) {
		helperTbar_->makeSummaryPlots();
	}
}

unsigned int IsolationAnalyzer::getClosestJet(edm::Handle<TopJetCollection> & j, const pat::MET* &r) {
	//	cout << "closest" << endl;
	unsigned int i = 0;
	if (j->size() >= 4) {
		//		cout << "more than 4 jets" << endl;
		double dr2 = 999.;
		TopJetCollection::const_iterator jet = j->begin();
		for (unsigned x = 0; x < 4; x++) {
			if (jet == j->end())
				break;
			//			cout << "jet" << x << endl;
			double temp = deltaR2(*jet, *r);
			if (temp < dr2) {
				dr2 = temp;
				i = x + 1;
			}
			++jet;
		}
	}
	return i;
}
