#include "TopAnalysis/TopAnalyzer/plugins/IsolationAnalyzer.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "TopQuarkAnalysis/TopTools/interface/EventShapeVariables.h"
#include <math.h>

using std::cout;
using std::endl;
using reco::GenParticle;
using edm::LogInfo;

double newWeightPt(double var, bool signalLike = true) {
	double y1 = 0.1 * var - 1.8;
	double y2 = 1 / y1;
	double ret = 1.;
	if (var <= 80) {
		if (signalLike)
			ret = y1;
		else
			ret = y2;
	}
	return ret;
}
double newWeightDphiMuJ1J2(double var, bool signalLike = true) {
	var = fabs(var);
	double y1 = (1 / M_PI / M_PI) * var * var - 2 / M_PI * var + 1;
	double y2 = 1 - y1;
	double ret = 1.;
	if (signalLike)
		ret = y1;
	else
		ret = y2;
	return ret;
}

double newWeightDphiMETMu(double var, bool signalLike = true) {
	var = fabs(var);
	double y1 = 1.5 / M_PI * var + 0.5;
	double y2 = 1 / y1;
	double ret = 1;
	if (signalLike)
		ret = y1;
	else
		ret = y2;
	return ret;
}

double newWeightCirc(double var, bool signalLike = true) {
	double y1 = var * var - 2 * var + 1;
	double y2 = 1-y1;
	double ret = 1.;
	if (signalLike)
		ret = y2;
	else
		ret = y1;
	return ret;
}
/**
 * Tests several variables for correlation with lepton isolation
 * For more information see IsolationAnalyzer twiki page
 * Twiki page link:
 * http://
 */
IsolationAnalyzer::IsolationAnalyzer(const edm::ParameterSet& cfg) :
	hist_(cfg.getParameter<std::string> ("hist")), muons_(cfg.getParameter<edm::InputTag> ("muons")), met_(
			cfg.getParameter<edm::InputTag> ("missingEt")), ttgen_(cfg.getParameter<edm::InputTag> ("genEvent")),
			jets_(cfg.getParameter<edm::InputTag> ("jets")),
			ptBins_(cfg.getParameter<std::vector<double> > ("ptBins")), ttbarMC_(cfg.getParameter<bool> ("ttbarMC")),
			useMVA_(cfg.getParameter<bool> ("useMVA")), module_(cfg.getParameter<std::string> ("modulename")),
			discinput_(cfg.getParameter<std::string> ("discriminator")) {
	ttbarMC_ = false;
}

void IsolationAnalyzer::beginJob(const edm::EventSetup&) {
	if (hist_.empty())
		return;

	edm::Service<TFileService> fs;
	if (!fs) {
		throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
	}

	helper_ = new IsolationHelper(fs, hist_);

	//disable ttbarMC for now
	if (ttbarMC_) {
		//for each ptBin new IsolationHelper
		for (unsigned int x = 0; x < ptBins_.size() - 1; x++) {
			IsolationHelper *helper = new IsolationHelper(fs, "ttbarMC_" + hist_);
			double bin = ptBins_.at(x);
			double nextBin = ptBins_.at(x + 1);
			stringstream name;
			name << ptBins_.at(x);
			stringstream pres;
			pres << "Ttbar_ " << bin << "-" << nextBin << "_";
			string pre = pres.str();
			helper->addHistogram(pre + "invariantMassJ3andJ4", 40, 0., 80.);
			helper->addHistogram(pre + "minDeltaPhiMETJets", 100, 0., 4.);
			helper->addHistogram(pre + "deltaPhiMetJet1", 80, -4., 4.);
			helper->addHistogram(pre + "deltaPhiMetJet2", 80, -4., 4.);
			helper->addHistogram(pre + "deltaPhiMetJet3", 80, -4., 4.);
			helper->addHistogram(pre + "deltaPhiMetJet4", 80, -4., 4.);
			helper->addHistogram(pre + "deltaPhiMetleadingMuon", 80, -4., 4.);
			helper->addHistogram(pre + "deltaPhiMetMuons", 80, -4., 4.);
			helper->addHistogram(pre + "DeltaPhiTimesDeltaEta", 50, -10., 10.);
			helper->addHistogram(pre + "METTimesleadingJetEt", 70, 0., 7000.);
			helper->addHistogram(pre + "Jet3EtOverJet1EtJet3Et", 25, 0., 0.5);
			helper->addHistogram(pre + "Jet3EtOverJet2EtJet3Et", 25, 0., 0.5);
			helper->addHistogram(pre + "Jet4EtOverJet1EtJet4Et", 25, 0., 0.5);
			helper->addHistogram(pre + "Jet4EtOverJet2EtJet3Et", 25, 0., 0.5);
			helper->addHistogram(pre + "JetEtSum34", 150, 0., 300.);
			helper->addHistogram(pre + "DeltaPhiMuonJet3", 80, -4., 4.);
			helper->addHistogram(pre + "DeltaPhiMuonJet4", 80, -4., 4.);
			helper->addHistogram(pre + "DeltaPhiMuonJet12", 80, -4., 4.);
			helper->addHistogram(pre + "SumJet1ET2TimesMuEt", 80, 60., 500);
			helper->addHistogram(pre + "Sum4JetsMuEt", 100, 100., 600);
			helper->addHistogram(pre + "VectorSumJetsMu", 60, 0., 300.);
			helper->addHistogram(pre + "MET", 70, 0., 350.);
			ttBarHelper_.insert(make_pair(name.str(), helper));
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
	helper_->addHistogram("DeltaPhiMuonJet12", 80, -7., 7.);
	helper_->addHistogram("SumJet1ET2TimesMuEt", 80, 60., 500);
	helper_->addHistogram("Sum4JetsMuEt", 100, 100., 600);
	helper_->addHistogram("VectorSumJetsMu", 70, 0., 350.);
	helper_->addHistogram("MET", 70, 0., 350.);
	helper_->addHistogram("deltaPhiTtbar", 80, -4., 4.);
	helper_->addHistogram("fabsDeltaPhiTtbar", 80, 0, 4.);
	helper_->addHistogram("deltaPhiJet1Jet2", 80, -4, 4.);
	helper_->addHistogram("TriJetTMass", 70, 60., 350); //5GeV Schritte
	helper_->addHistogram("muon_pt", 24, 0., 120.);
	helper_->addHistogram("circularity", 100, 0., 1.);
	if (useMVA_)
		helper_->addHistogram("MVAdisc", 20, 0., 1.);
	NameScheme nam("var");
	ofstream off(hist_.c_str(), std::ios::app);
	sumDeltaPhiMuvsdeltaPhiJ1J2_ = fs->make<TH2F> (nam.name(off, "dPhiSumMuTwoJestvsdPhij1j2"), nam.name(
			"dPhiSumMuTwoJestvsdPhij1j2"), 80, -7., 7., 80, -4., 4.);
	recoMETUncorrectedMET_ = fs->make<TH1F> (nam.name(off, "recoMETUncorrectedMET"), nam.name("recoMETUncorrectedMET"),
			200, -100., 100.);
	genMetRecoDiff_ = fs->make<TH1F> (nam.name(off, "recoMETGenMET"), nam.name("recoMETGenMET"), 200, -100., 100.);
	norm_genMetRecoDiff_
			= fs->make<TH1F> (nam.name(off, "norm_recoMETGenMET"), nam.name("recoMETGenMET"), 51, -1., 50.);
}

IsolationAnalyzer::~IsolationAnalyzer() {
}

void IsolationAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& setup) {

	edm::Handle<double> disc_handle;
	double disc = 0.;
	if (useMVA_) {
		evt.getByLabel(module_, discinput_, disc_handle);
		disc = *disc_handle;
	}

	const pat::MET *met;
	reco::Particle::LorentzVector vec;

	edm::Handle<std::vector<pat::MET> > metH;
	evt.getByLabel(met_, metH);

	edm::Handle<TopMuonCollection> muons;
	evt.getByLabel(muons_, muons);

	//the genEvent
	edm::Handle<TtGenEvent> genEvt;
	evt.getByLabel(ttgen_, genEvt);

	edm::Handle<double> weightHandle;
	evt.getByLabel("eventWeight", weightHandle);

	edm::Handle<TopJetCollection> jets;
	evt.getByLabel(jets_, jets);

	met = &(*metH)[0];

	double weight = *weightHandle;
	std::vector<TVector3> p;

	if (jets->size() >= 4) {
		//get the first four jets
		TopJetCollection::const_iterator jet = jets->begin();

		double jet1Phi, jet2Phi, jet3Phi, jet4Phi;
		double jet1Eta, jet2Eta, jet3Eta, jet4Eta;
		double jet1pt, jet2pt, jet3pt, jet4pt;
		double jet1Et, jet2Et, jet3Et, jet4Et;
		double dpTde = 10000.;
		pat::Jet closMETJet = getClosestJet(jets, *met);

		jet1pt = jet->pt();
		jet1Et = jet->et();
		jet1Eta = jet->eta();
		jet1Phi = jet->phi();
		TVector3 jet1(jet->px(), jet->py(), jet->pz());
		p.push_back(jet1);

		vec = jet->p4();

		++jet;
		jet2pt = jet->pt();
		jet2Et = jet->et();
		jet2Eta = jet->eta();
		jet2Phi = jet->phi();
		vec += jet->p4();
		TVector3 jet2(jet->px(), jet->py(), jet->pz());
		p.push_back(jet2);

		++jet;
		jet3pt = jet->pt();
		jet3Et = jet->et();
		jet3Eta = jet->eta();
		jet3Phi = jet->phi();
		vec += jet->p4();
		TVector3 jet3(jet->px(), jet->py(), jet->pz());
		p.push_back(jet3);

		++jet;
		jet4pt = jet->pt();
		jet4Et = jet->et();
		jet4Eta = jet->eta();
		jet4Phi = jet->phi();
		TVector3 jet4(jet->px(), jet->py(), jet->pz());
		p.push_back(jet4);

		//how big is the minimal deltaPhi between MET and one of the 4 jets
		double mindp;
		mindp = min(fabs(deltaPhi(met->phi(), jet1Phi)), fabs(deltaPhi(met->phi(), jet2Phi)));
		mindp = min(mindp, fabs(deltaPhi(met->phi(), jet3Phi)));
		mindp = min(mindp, fabs(deltaPhi(met->phi(), jet4Phi)));
		//transversal mass of jet 3+4
		double mt34 = sqrt((jet3Et + jet4Et) * (jet3Et + jet4Et) - (jet3pt + jet4pt) * (jet3pt + jet4pt));
		//transversal mass of jet 1+3+4 (jet1 could be switched with jet2)
		double mt3Jet1 = sqrt((jet1Et + jet3Et + jet4Et) * (jet1Et + jet3Et + jet4Et) - (jet1pt + jet3pt + jet4pt)
				* (jet1pt + jet3pt + jet4pt));
		double mt3Jet2 = sqrt((jet2Et + jet3Et + jet4Et) * (jet2Et + jet3Et + jet4Et) - (jet2pt + jet3pt + jet4pt)
				* (jet2pt + jet3pt + jet4pt));
		double mt3Jet = 0.;
		//TODO: make W-mass as config parameter
		if (fabs(mt3Jet1 - 170) < fabs(mt3Jet2 - 170)) {
			//mt3Jet1 closer to W mass
			mt3Jet = mt3Jet1;
		} else
			mt3Jet = mt3Jet2;

		if (&closMETJet)
			dpTde = fabs(deltaPhi(met->phi(), closMETJet.phi())) * fabs(met->eta() - closMETJet.eta());

		//set the event-weight
		helper_->setWeight(weight);
		//now the muons
		double jetIso = getJetIso(jets, muons, 20);
		unsigned int size = muons->size();
		for (unsigned int i = 0; i < size; ++i) {
			pat::Muon mu = (pat::Muon)(*muons)[i];
			double caloIso = mu.caloIso();
			double trackIso = mu.trackIso();

			//set isolation
			helper_->setCaloIso(caloIso);
			helper_->setTrackIso(trackIso);
			helper_->setJetIso(jetIso);
			if (i == 0) {
				TVector3 lep(mu.px(), mu.py(), mu.pz());
				p.push_back(lep);
				EventShapeVariables eventshape;
				//leading muon
				vec += mu.p4();
				//complicated variable:
				double phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8;
				//get closest jet to mu (bjet from semileptonic t-decay)
				pat::Jet closestMuJet = getClosestJet(jets, mu);
				//get 4th jet + closest (hadronic W)
				phi1 = jet3Phi + (0.5 * deltaPhi(jet4Phi, jet3Phi));
				phi2 = jet4Phi + (0.5 * deltaPhi(jet3Phi, jet4Phi));
				//				if(phi1 == phi2) cout << "it's working!" << endl;

				pat::Jet closest34Jet = getClosestJetInDeltaPhi(jets, phi1);
				phi3 = phi1 + 0.5 * deltaPhi(closest34Jet.phi(), phi1);
				phi4 = closest34Jet.phi() + 0.5 * deltaPhi(phi1, closest34Jet.phi());
				//				if(phi3 == phi4) cout << "it's really working!" << endl;

				phi5 = mu.phi() + 0.5 * deltaPhi(met->phi(), mu.phi());
				phi6 = met->phi() + 0.5 * deltaPhi(mu.phi(), met->phi());
				//				if(phi5 == phi6) cout << "it's really, really working!" << endl;

				phi7 = phi5 + 0.5 * deltaPhi(closestMuJet.phi(), phi5);
				phi8 = closestMuJet.phi() + 0.5 * deltaPhi(phi5, closestMuJet.phi());
				//				if(phi7 == phi8) cout << "katsching!" << endl;
				//				cout << "<" << phi1 << "," << phi2 << ">" << endl;
				//				cout << "<" << phi3 << "," << phi4 << ">" << endl;
				//				cout << "<" << phi5 << "," << phi6 << ">" << endl;
				//				cout << "<" << phi7 << "," << phi8 << ">" << endl;
				//				cout << "here <" << deltaPhi(jet3Phi, jet4Phi) << "," << deltaPhi(jet4Phi, jet3Phi) << ">" << endl;
				//fill histogramms
				helper_->setWeight(weight * newWeightCirc(eventshape.circularity(p)));
				helper_->fill("METTimesleadingJetEt", jet1Et * met->et());
				helper_->fill("minDeltaPhiMETJets", mindp);
				helper_->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
				helper_->fill("deltaPhiMetJet1", deltaPhi(met->phi(), jet1Phi));
				helper_->fill("deltaPhiMetJet2", deltaPhi(met->phi(), jet2Phi));
				helper_->fill("deltaPhiMetJet3", deltaPhi(met->phi(), jet3Phi));
				helper_->fill("deltaPhiMetJet4", deltaPhi(met->phi(), jet4Phi));
				helper_->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
				helper_->fill("DeltaPhiTimesDeltaEta", deltaPhi(met->phi(), jet1Phi) * fabs(met->eta() - jet1Eta));
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
				helper_->fill("TriJetTMass", mt3Jet);
				helper_->fill("deltaPhiTtbar", deltaPhi(phi8, phi4));
				helper_->fill("fabsDeltaPhiTtbar", fabs(deltaPhi(phi8, phi4)));
				helper_->fill("deltaPhiJet1Jet2", deltaPhi(jet1Phi, jet2Phi));
				helper_->fill("muon_pt", mu.pt());
				helper_->fill("circularity", eventshape.circularity(p));

				sumDeltaPhiMuvsdeltaPhiJ1J2_->Fill(deltaPhi(mu.phi(), jet1Phi) + deltaPhi(mu.phi(), jet2Phi), deltaPhi(
						jet1Phi, jet2Phi), weight);
				recoMETUncorrectedMET_->Fill(met->et() - met->uncorrectedPt(), weight);
				double genMetreco = (met->et() - met->genMET()->et());
				genMetRecoDiff_->Fill(genMetreco, weight);
				norm_genMetRecoDiff_->Fill(genMetreco / met->genMET()->et(), weight);
				if (useMVA_)
					helper_->fill("MVAdisc", disc);
			}
			if (i == 1) {
				//2nd leading muon
			}
			helper_->fill("deltaPhiMetMuons", deltaPhi(met->phi(), mu.phi()));

			if (ttbarMC_) {
				TtGenEvent event = *genEvt;
				const reco::GenParticle* thad = event.hadronicDecayTop();
				const reco::GenParticle* tlep = event.leptonicDecayTop();

				if (i == 0) {
					for (unsigned int x = 0; x < ptBins_.size() - 1; x++) {
						if (thad && tlep) {
							//get the responsible helper here
							IsolationHelper *helper;
							stringstream name;
							name << ptBins_.at(x);
							map<string, IsolationHelper*>::iterator iter = ttBarHelper_.find(name.str());
							if (iter != ttBarHelper_.end()) {
								helper = iter->second;
								if (thad->pt() >= ptBins_[x] && thad->pt() < ptBins_[x + 1]) {
									helper->fill("METTimesleadingJetEt", jet1Et * met->et());
									helper->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
									helper->fill("deltaPhiMetJet1", deltaPhi(met->phi(), jet1Phi));
									helper->fill("deltaPhiMetJet2", deltaPhi(met->phi(), jet2Phi));
									helper->fill("deltaPhiMetJet3", deltaPhi(met->phi(), jet3Phi));
									helper->fill("deltaPhiMetJet4", deltaPhi(met->phi(), jet4Phi));
									helper->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
									helper->fill("DeltaPhiTimesDeltaEta", deltaPhi(met->phi(), jet1Phi) * (met->eta()
											* jet1Eta));
									helper->fill("METTimesleadingJetEt", met->et() * jet1Et);
									helper->fill("Jet3EtOverJet1EtJet3Et", jet3Et / (jet3Et + jet1Et));
									helper->fill("Jet3EtOverJet2EtJet3Et", jet3Et / (jet3Et + jet2Et));
									helper->fill("Jet4EtOverJet1EtJet4Et", jet4Et / (jet4Et + jet1Et));
									helper->fill("Jet4EtOverJet2EtJet3Et", jet4Et / (jet4Et + jet2Et));
									helper->fill("JetEtSum34", jet3Et + jet4Et);
									helper->fill("DeltaPhiMuonJet3", deltaPhi(mu.phi(), jet3Phi));
									helper->fill("DeltaPhiMuonJet4", deltaPhi(mu.phi(), jet4Phi));
									helper->fill("DeltaPhiMuonJet12", deltaPhi(mu.phi(), jet1Phi) + deltaPhi(mu.phi(),
											jet2Phi));
									helper->fill("invariantMassJ3andJ4", mt34);
									helper->fill("Sum4JetsMuEt", jet1Et + jet2Et + jet3Et + jet4Et + mu.et());
									helper->fill("SumJet1ET2TimesMuEt", jet1Et + 2 * mu.et());
									helper->fill("VectorSumJetsMu", vec.pt());
									helper->fill("MET", met->et());
								}

								if (tlep->pt() >= ptBins_[x] && tlep->pt() < ptBins_[x + 1]) {
									helper->fill("METTimesleadingJetEt", jet1Et * met->et());
									helper->fill("deltaPhiMetleadingMuon", deltaPhi(met->phi(), mu.phi()));
									helper->fill("deltaPhiMetJet1", deltaPhi(met->phi(), jet1Phi));
									helper->fill("deltaPhiMetJet2", deltaPhi(met->phi(), jet2Phi));
									helper->fill("deltaPhiMetJet3", deltaPhi(met->phi(), jet3Phi));
									helper->fill("deltaPhiMetJet4", deltaPhi(met->phi(), jet4Phi));
									helper->fill("DeltaPhiTimesDeltaEta", deltaPhi(met->phi(), jet1Phi) * (met->eta()
											* jet1Eta));
									helper->fill("METTimesleadingJetEt", met->et() * jet1Et);
									helper->fill("Jet3EtOverJet1EtJet3Et", jet3Et / (jet3Et + jet1Et));
									helper->fill("Jet3EtOverJet2EtJet3Et", jet3Et / (jet3Et + jet2Et));
									helper->fill("Jet4EtOverJet1EtJet4Et", jet4Et / (jet4Et + jet1Et));
									helper->fill("Jet4EtOverJet2EtJet3Et", jet4Et / (jet4Et + jet2Et));
									helper->fill("JetEtSum34", jet3Et + jet4Et);
									helper->fill("DeltaPhiMuonJet3", deltaPhi(mu.phi(), jet3Phi));
									helper->fill("DeltaPhiMuonJet4", deltaPhi(mu.phi(), jet4Phi));
									helper->fill("DeltaPhiMuonJet12", deltaPhi(mu.phi(), jet1Phi) + deltaPhi(mu.phi(),
											jet2Phi));
									helper->fill("invariantMassJ3andJ4", mt34);
									helper->fill("Sum4JetsMuEt", jet1Et + jet2Et + jet3Et + jet4Et + mu.et());
									helper->fill("SumJet1ET2TimesMuEt", jet1Et + 2 * mu.et());
									helper->fill("VectorSumJetsMu", vec.pt());
									helper->fill("MET", met->et());
								}
							}
						}
					}
				}
			}
		}

	}
}

void IsolationAnalyzer::endJob() {
	helper_->makeSummaryPlots();
	helper_->normalize();
	if (ttbarMC_) {
		//loop over all helpers and make summaryPlots and normalize them
		map<string, IsolationHelper*>::iterator iter;
		for (iter = ttBarHelper_.begin(); iter != ttBarHelper_.end(); ++iter) {
			IsolationHelper *helper = iter->second;
			helper->makeSummaryPlots("Ttbar");
			helper->normalize();
		}
	}
	if (recoMETUncorrectedMET_->Integral() != 0) {
		recoMETUncorrectedMET_->Scale(1 / recoMETUncorrectedMET_->Integral());
	}
}
template<class T>
pat::Jet IsolationAnalyzer::getClosestJet(edm::Handle<TopJetCollection> & j, T r) {
	pat::Jet closestJet;
	if (j->size() >= 4) {
		double dr2 = 999.;
		TopJetCollection::const_iterator jet = j->begin();
		//only the four highest jets
		for (unsigned x = 0; x < 4; x++) {
			if (jet == j->end())
				break;
			double temp = deltaR2(*jet, r);
			if (temp < dr2) {
				dr2 = temp;
				closestJet = *jet;
			}
			++jet;
		}
	}
	return closestJet;
}

pat::Jet IsolationAnalyzer::getClosestJetInDeltaPhi(edm::Handle<TopJetCollection> & j, double phi) {
	pat::Jet closestJet;
	if (j->size() >= 4) {
		double dr2 = 999.;
		TopJetCollection::const_iterator jet = j->begin();
		//only the four highest jets
		for (unsigned x = 0; x < 4; x++) {
			if (jet == j->end())
				break;
			double temp = fabs(deltaPhi(jet->phi(), phi));
			if (temp < dr2) {
				dr2 = temp;
				closestJet = *jet;
			}
			++jet;
		}
	}
	return closestJet;
}
