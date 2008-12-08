#include "TopAnalysis/TopAnalyzer/plugins/MatrixAnalyzer.h"

using std::cout;
using std::endl;
using reco::GenParticle;

MatrixAnalyzer::MatrixAnalyzer(const edm::ParameterSet& cfg) :
	hist_(cfg.getParameter<std::string> ("hist")),
	muons_(cfg.getParameter<edm::InputTag> ("muons")),
	var_(cfg.getParameter<edm::InputTag> ("var")),
	jets_(cfg.getParameter<edm::InputTag> ("jets")),
	bins_(cfg.getParameter<edm::ParameterSet> ("bins" )),
	useMVA_(cfg.getParameter<bool> ("useMVA")),
	module_(cfg.getParameter<std::string> ("modulename")),
	discinput_(cfg.getParameter<std::string> ("discriminator")) {
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
	if (useMVA_) {
		for (unsigned int x = 0; x < mvaDiscBins_.size() - 1; x++) {
			std::stringstream tmp;
			tmp << "mva" << mvaDiscBins_.at(x);
			Counters_->addCounter(tmp.str());
			Counters_->addCounter(tmp.str() + "simple");
		}
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
	nVSmet_ = fs->make<TH1F> (nam.name(hist, "nVSmet"), nam.name("numberOfeventsVSmet"), varBins_.size() - 1,
			&varBins_[0]);
	if (useMVA_) {
		nVSdisc_ = fs->make<TH1F> (nam.name(hist, "nVSdisc"), nam.name("numberOfeventsVSmvaDiscriminator"),
				mvaDiscBins_.size() - 1, &mvaDiscBins_[0]);
		nVSdiscSimple_ = fs->make<TH1F> (nam.name(hist, "nVSdiscSimple"), nam.name(
				"numberOfeventsVSmvaDiscriminatorSimple"), mvaDiscBins_.size() - 1, &mvaDiscBins_[0]);
	}


	nVSmetSimple_ = fs->make<TH1F> (nam.name(hist, "nVSmetSimple"), nam.name("numberOfeventsVSmetSimple"),
			varBins_.size() - 1, &varBins_[0]);

	varPlot_ = fs->make<TH1F> (nam.name(hist, var_.label().c_str()), nam.name(var_.label().c_str()), 250, 0, 500);
}

void MatrixAnalyzer::analyze(const edm::Event& evt,
		const edm::EventSetup& setup) {

	edm::Handle<double> disc_handle;
	double disc = 0.;
	if (useMVA_){
		evt.getByLabel(module_, discinput_, disc_handle);
		disc = *disc_handle;
	}

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
		//bigger than bin, smaller than next bin if it is not the last bin
		bool pass = pt >= varBins_[i] && (pt < varBins_[i + 1] || i == varBins_.size() - 2);

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
	if (useMVA_) {
		for (unsigned int i = 0; i < mvaDiscBins_.size() - 1; i++) {
			bool pass = disc >= mvaDiscBins_[i] && (disc < mvaDiscBins_[i + 1] || i == mvaDiscBins_.size() - 2);

			if (pass) {
				double t = mvaDiscBins_.at(i);
				std::stringstream tmp;
				tmp << "mva" << t;
				if (muons->size() < 1) {
					Counters_->addPureHadronic(tmp.str() + "simple");
					Counters_->addPureHadronic(tmp.str(), sampleweight_);
				} else if (muons->size() == 1) {
					Counters_->addSemiLeptonic(tmp.str() + "simple");
					Counters_->addSemiLeptonic(tmp.str(), sampleweight_);
				} else if (muons->size() == 2) {
					Counters_->addDiLeptonic(tmp.str() + "simple");
					Counters_->addDiLeptonic(tmp.str(), sampleweight_);
				} else if (muons->size() > 2) {
					Counters_->addMultiLeptonic(tmp.str() + "simple");
					Counters_->addMultiLeptonic(tmp.str(), sampleweight_);
				}
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
		double s, _ss;
		int bbin = i+1;

		s = Counters_->getSemiLeptonic(tmp.str());

		_ss = Counters_->getSemiLeptonic(tmp.str() + "simple");

		nVSmet_->SetBinContent(bbin, s);
		nVSmetSimple_->SetBinContent(bbin, _ss);

	}
	if (useMVA_) {
		for (unsigned int i = 0; i < mvaDiscBins_.size() - 1; i++) {
			std::stringstream tmp;
			tmp << "mva" << mvaDiscBins_.at(i);
			double s, _ss;//, b, d, m, o, _sb, _sd, _sm, _so;
			int bbin = i + 1;

			s = Counters_->getSemiLeptonic(tmp.str());
			_ss = Counters_->getSemiLeptonic(tmp.str() + "simple");

			nVSdisc_->SetBinContent(bbin, s);
			nVSdiscSimple_->SetBinContent(bbin, _ss);
		}
	}
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

