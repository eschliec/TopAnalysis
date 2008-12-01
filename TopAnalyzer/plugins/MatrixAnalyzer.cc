#include "TopAnalysis/TopAnalyzer/plugins/MatrixAnalyzer.h"

using std::cout;
using std::endl;
using reco::GenParticle;

MatrixAnalyzer::MatrixAnalyzer(const edm::ParameterSet& cfg) :
	hist_(cfg.getParameter<std::string> ("hist")),
	muons_(cfg.getParameter<edm::InputTag> ("muons")),
	var_(cfg.getParameter<edm::InputTag> ("var")),
	jets_(cfg.getParameter<edm::InputTag> ("jets")),
	bins_( cfg.getParameter<edm::ParameterSet> ("bins" ) ){
	varBins_ = bins_.getParameter<std::vector<double> >( "metBins" );
	mvaDiscBins_ = bins_.getParameter<std::vector<double> >( "mvaDiscBins" );
	debug_ = true;
	noBins_ = 5;
	sampleweight_ = 0.;
//	notNeededHists_ = "QCDnotNeeded.hist";
	mvamodule_ = "";


	Counters_ = new LeptonCounter();

	for (unsigned int x = 0; x < varBins_.size() - 1; x++) {
		std::stringstream tmp;
		tmp << "MET" << varBins_.at(x);
		Counters_->addCounter(tmp.str());
		Counters_->addCounter(tmp.str() + "simple");
	}

	for (unsigned int x = 0; x < mvaDiscBins_.size() - 1; x++) {
		std::stringstream tmp;
		tmp << "mva" << mvaDiscBins_.at(x);
		Counters_->addCounter(tmp.str());
		Counters_->addCounter(tmp.str() + "simple");
	}
	Counters_->addCounter("weighted");

	Counters_->addCounter("simple");
}

MatrixAnalyzer::~MatrixAnalyzer() {
}

void MatrixAnalyzer::beginJob(const edm::EventSetup&) {

	if (hist_.empty())
		return;

	edm::Service<TFileService> fs;
	if (!fs) {
		throw edm::Exception(edm::errors::Configuration,
				"TFile Service is not registered in cfg file");
	}
	//f_ = &fs->file();
	ofstream hist(hist_.c_str(), std::ios::out);
	//ofstream non(notNeededHists_.c_str(), std::ios::out);

	NameScheme nam("mbg");

	//make histograms
//	background_ = fs->make<TH1F> (nam.name(non, "N_background"), nam.name(
//			"background"), 5, 0, 5.);
//	lep_ = fs->make<TH1F> (nam.name(non, "N_semilepton"), nam.name("lepton"),
//			5, 0, 5.);
//	llep_ = fs->make<TH1F> (nam.name(non, "N_dilepton"), nam.name("ll"), 5, 0, 5.);
//	multilep_ = fs->make<TH1F> (nam.name(non, "N_multilepton"),
//			nam.name("ml"), 5, 0, 5.);
//	overall_ = fs->make<TH1F> (nam.name(non, "N_overall"),
//			nam.name("overall"), noBins_, 0, noBins_);
//
//	binnedBkg_ = fs->make<TH1F> (nam.name(non, "binnedBkg"), nam.name(
//			"binnedBkg"), varBins1_.size() -1, &varBins1_[0]);
//	binnedSemiLep_ = fs->make<TH1F> (nam.name(non, "binnedSemiLep"), nam.name(
//			"binnedSemiLep"), varBins1_.size() -1, &varBins1_[0]);
//	binnedDiLep_ = fs->make<TH1F> (nam.name(non, "binnedDiLep"), nam.name(
//			"binnedDiLep"), varBins1_.size() -1, &varBins1_[0]);
//	binnedMultiLep_ = fs->make<TH1F> (nam.name(non, "binnedMultiLep"),
//			nam.name("binnedMultiLep"), varBins1_.size() -1, &varBins1_[0]);
//
//	binnedOverall_ = fs->make<TH1F> (nam.name(non, "binnedOverall"), nam.name(
//			"binnedOverall"), varBins1_.size() -1, &varBins1_[0]);
//
//	binnedSimpleBkg_ = fs->make<TH1F> (nam.name(non, "binnedSimpleBkg"),
//			nam.name("binnedSimpleBkg"), varBins1_.size() -1, &varBins1_[0]);
//	binnedSimpleSemiLep_ = fs->make<TH1F> (
//			nam.name(non, "binnedSimpleSemiLep"), nam.name(
//					"binnedSimpleSemiLep"), varBins1_.size() -1, &varBins1_[0]);
//	binnedSimpleDiLep_ = fs->make<TH1F> (nam.name(non, "binnedSimpleDiLep"),
//			nam.name("binnedSimpleDiLep"), varBins1_.size() -1, &varBins1_[0]);
//	binnedSimpleMultiLep_ = fs->make<TH1F> (nam.name(non,
//			"binnedSimpleMultiLep"), nam.name("binnedSimpleMultiLep"),
//			varBins1_.size() -1, &varBins1_[0]);
//	binnedSimpleOverall_ = fs->make<TH1F> (
//			nam.name(non, "binnedSimpleOverall"), nam.name(
//					"binnedSimpleOverall"), varBins1_.size() -1, &varBins1_[0]);
	nVSmet_ = fs->make<TH1F> (nam.name(hist, "nVSmet"), nam.name("numberOfeventsVSmet"), varBins_.size() - 1,
			&varBins_[0]);
	nVSdisc_ = fs->make<TH1F> (nam.name(hist, "nVSdisc"), nam.name("numberOfeventsVSmvaDiscriminator"),
			mvaDiscBins_.size() - 1, &mvaDiscBins_[0]);

	nVSmetSimple_ = fs->make<TH1F> (nam.name(hist, "nVSmetSimple"), nam.name("numberOfeventsVSmetSimple"),
			varBins_.size() - 1, &varBins_[0]);
	nVSdiscSimple_ = fs->make<TH1F> (nam.name(hist, "nVSdiscSimple"),
			nam.name("numberOfeventsVSmvaDiscriminatorSimple"), mvaDiscBins_.size() - 1, &mvaDiscBins_[0]);
	varPlot_ = fs->make<TH1F> (nam.name(hist, var_.label().c_str()), nam.name(var_.label().c_str()), 250, 0, 500);
}

void MatrixAnalyzer::analyze(const edm::Event& evt,
		const edm::EventSetup& setup) {

	edm::Handle<double> disc_handle;
	evt.getByLabel("findTtSemiLepSignalSelectorMVA", "DiscSel", disc_handle);
	double disc = *disc_handle;

	edm::Handle<reco::GenParticleCollection> genParticles;
	evt.getByLabel("genParticles", genParticles);

	edm::Handle<TopMuonCollection> muons;
	evt.getByLabel(muons_, muons);

	edm::Handle<edm::View<pat::Muon> > recMuons;
	evt.getByLabel(muons_, recMuons);

	edm::Handle<edm::View<reco::RecoCandidate> > recVars;
	evt.getByLabel(var_, recVars);

//	edm::Handle<TopJetCollection> jets;
//	evt.getByLabel(jets_, jets);

	edm::Handle<edm::Association<reco::GenParticleCollection> > genMatch;
	evt.getByLabel("muonMatch", genMatch);

	edm::Handle<double> weightHandle;
	evt.getByLabel("eventWeight", weightHandle);

	sampleweight_ = *weightHandle;

	//TopJetCollection::const_iterator jet = jets->begin();
	Double_t pt = (*recVars)[0].pt();
//	Double_t pt1 = jet->pt();
//	jet++;
//	Double_t pt2 = jet->pt();
//	jet++;
//	Double_t pt3 = jet->pt();
//	jet++;
//	Double_t pt4 = jet->pt();
	varPlot_->Fill(pt, sampleweight_);
	for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
//		bool pass1 = pt1 >= varBins1_[i] && (pt1 < varBins1_[i + 1] || i == varBins1_.size() - 2);
//		bool pass2 = pt2 >= varBins2_[i];
//		bool pass3 = pt3 >= varBins2_[i];
//		bool pass4 = pt4 >= varBins2_[i];
		bool pass = pt >= varBins_[i] && (pt < varBins_[i + 1] || i == varBins_.size() - 2);
//		bool passAll = pass1 && pass2 && pass3 && pass4 && pass;
		//redefine passAll if MET instead of Jets should be used

		if (pass) {
			double t = varBins_.at(i);
			std::stringstream tmp;
			tmp << "MET" << t;
			if (muons->size() < 1) {
				Counters_->addPureHadronic(tmp.str() + "simple");
				Counters_->addPureHadronic(tmp.str(), sampleweight_);
			}
			else if (muons->size() == 1) {
				Counters_->addSemiLeptonic(tmp.str() + "simple");
				Counters_->addSemiLeptonic(tmp.str(), sampleweight_);
			}
			else if (muons->size() == 2) {
				Counters_->addDiLeptonic(tmp.str() + "simple");
				Counters_->addDiLeptonic(tmp.str(),sampleweight_);
			}
			else if (muons->size() > 2) {
				Counters_->addMultiLeptonic(tmp.str() + "simple");
				Counters_->addMultiLeptonic(tmp.str(), sampleweight_);
			}
		}
	}

	for (unsigned int i = 0; i < mvaDiscBins_.size() - 1; i++) {
			bool pass = disc >= mvaDiscBins_[i] && (disc < mvaDiscBins_[i + 1] || i == mvaDiscBins_.size() - 2);

			if (pass) {
				double t = mvaDiscBins_.at(i);
				std::stringstream tmp;
				tmp << "mva" << t;
				if (muons->size() < 1) {
					Counters_->addPureHadronic(tmp.str() + "simple");
					Counters_->addPureHadronic(tmp.str(), sampleweight_);
				}
				else if (muons->size() == 1) {
					Counters_->addSemiLeptonic(tmp.str() + "simple");
					Counters_->addSemiLeptonic(tmp.str(), sampleweight_);
				}
				else if (muons->size() == 2) {
					Counters_->addDiLeptonic(tmp.str() + "simple");
					Counters_->addDiLeptonic(tmp.str(),sampleweight_);
				}
				else if (muons->size() > 2) {
					Counters_->addMultiLeptonic(tmp.str() + "simple");
					Counters_->addMultiLeptonic(tmp.str(), sampleweight_);
				}
			}
		}

	if (muons->size() < 1) {
		Counters_->addPureHadronic("simple");
		Counters_->addPureHadronic("weighted", sampleweight_);
	}

	if (muons->size() == 1) {
		Counters_->addSemiLeptonic("simple");
		Counters_->addSemiLeptonic("weighted", sampleweight_);
	}
	if (muons->size() == 2) {
		Counters_->addDiLeptonic("simple");
		Counters_->addDiLeptonic("weighted", sampleweight_);
	}
	if (muons->size() > 2) {
		Counters_->addMultiLeptonic("simple");
		Counters_->addMultiLeptonic("weighted", sampleweight_);
	}
}

void MatrixAnalyzer::endJob() {
	setEnv();
}

void MatrixAnalyzer::setEnv() {
	for (unsigned int i = 0; i < varBins_.size() - 1; i++) {
		std::stringstream tmp;
		tmp << "MET" << varBins_.at(i);
		double s, _ss;//, b, d, m, o, _sb, _sd, _sm, _so;
		int bbin = i+1;

		//		b = Counters_->getPureHadronic(tmp.str());
		s = Counters_->getSemiLeptonic(tmp.str());
		//		d = Counters_->getDiLeptonic(tmp.str());
		//		m = Counters_->getMultiLeptonic(tmp.str());
		//		o = Counters_->getAllLeptonic(tmp.str());
		//
		//		_sb = Counters_->getPureHadronic(tmp.str() + "simple");
		_ss = Counters_->getSemiLeptonic(tmp.str() + "simple");
		//		_sd = Counters_->getDiLeptonic(tmp.str() + "simple");
		//		_sm = Counters_->getMultiLeptonic(tmp.str() + "simple");
		//		_so = Counters_->getAllLeptonic(tmp.str() + "simple");

		nVSmet_->SetBinContent(bbin, s);
		nVSmetSimple_->SetBinContent(bbin, _ss);

//		binnedBkg_->SetBinContent(bbin, b);
//		binnedSemiLep_->SetBinContent(bbin, s);
//		binnedDiLep_->SetBinContent(bbin, d);
//		binnedMultiLep_->SetBinContent(bbin, m);
//		binnedOverall_->SetBinContent(bbin, o);
//
//		binnedSimpleBkg_->SetBinContent(bbin, _sb);
//		binnedSimpleSemiLep_->SetBinContent(bbin, _ss);
//		binnedSimpleDiLep_->SetBinContent(bbin, _sd);
//		binnedSimpleMultiLep_->SetBinContent(bbin, _sm);
//		binnedSimpleOverall_->SetBinContent(bbin, _so);
	}

	for (unsigned int i = 0; i < mvaDiscBins_.size() - 1; i++) {
			std::stringstream tmp;
			tmp << "mva" << mvaDiscBins_.at(i);
			double s, _ss;//, b, d, m, o, _sb, _sd, _sm, _so;
			int bbin = i+1;

			s = Counters_->getSemiLeptonic(tmp.str());
			_ss = Counters_->getSemiLeptonic(tmp.str() + "simple");

			nVSdisc_->SetBinContent(bbin, s);
			nVSdiscSimple_->SetBinContent(bbin, _ss);
		}

//	background_->SetBinContent(2, Counters_->getPureHadronic("weighted"));
//	lep_->SetBinContent(2, Counters_->getSemiLeptonic("weighted"));
//	llep_->SetBinContent(2, Counters_->getDiLeptonic("weighted"));
//	multilep_->SetBinContent(2, Counters_->getMultiLeptonic("weighted"));
//	overall_->SetBinContent(2, Counters_->getAllLeptonic("weighted"));
}

int MatrixAnalyzer::numberOfmatchedMuons(
		edm::Handle<edm::View<pat::Muon> > recMuons,
		edm::Handle<edm::Association<reco::GenParticleCollection> > genMatch) {
	int matchedmu = 0;
	int size = recMuons->size();
	cout << "<pdgId, status, mother.pdgId>" << endl;
	cout << "ID: ";
	for (int i = 0; i < size; ++i) {
		edm::RefToBase<reco::Candidate> muonRef =
				recMuons->refAt(i)->originalObjectRef();
		reco::GenParticleRef genMuon = (*genMatch)[muonRef];

		int mid = 0;
		int st = 0;
		if (genMuon->numberOfMothers() > 0) {
			mid = genMuon->mother(0)->pdgId();
			st = genMuon->mother(0)->status();
			cout << "mothers " << genMuon->numberOfMothers() << endl;
			if (genMuon->mother(0)->numberOfMothers() > 0) {
				cout << genMuon->mother(0)->pdgId() << endl;
				cout << genMuon->mother(0)->status() << endl;
				cout << genMuon->mother(0)->numberOfMothers() << endl;
			}
		}

		cout << "<" << genMuon->pdgId() << "," << genMuon->status() << ","
				<< mid << "," << st << "> ";
		if (fabs(mid) == 13)
			matchedmu++;
	}
	cout << "\n";
	return matchedmu;
}

